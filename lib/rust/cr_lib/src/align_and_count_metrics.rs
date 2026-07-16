#![deny(missing_docs)]

use crate::align_metrics::{AlignAndCountVisitor, BarcodeMetrics, LibFeatThenBarcodeOrder};
use crate::aligner::BarcodeSummary;
use anyhow::Result;
use cr_types::types::LibraryType;
use itertools::Itertools;
use martian_filetypes::LazyWrite;
use metric::TxHashMap;
use rand::RngExt;
use rand::rngs::SmallRng;
use shardio::ShardSender;
use std::collections::HashSet;
use std::collections::hash_map::Entry::{Occupied, Vacant};
use std::fs::File;
use std::io::BufWriter;
use transcriptome::Gene;
use tx_annotation::read::ReadAnnotations;
use tx_annotation::visitor::AnnotatedReadVisitor;

pub(super) struct StageVisitor<W>
where
    W: LazyWrite<ReadAnnotations, BufWriter<File>>,
{
    visitors: TxHashMap<LibraryType, AlignAndCountVisitor>,
    ann_writer: Option<(W, f32, SmallRng)>,
    // Number of reads written to the `ann_writer`
    ann_writer_num_reads: usize,
    metrics_sender: ShardSender<BarcodeMetrics, LibFeatThenBarcodeOrder>,
    target_genes: Option<HashSet<Gene>>,
    is_spatial: bool,
}

impl<W> StageVisitor<W>
where
    W: LazyWrite<ReadAnnotations, BufWriter<File>>,
{
    pub(super) fn new(
        metrics_sender: ShardSender<BarcodeMetrics, LibFeatThenBarcodeOrder>,
        target_genes: Option<HashSet<Gene>>,
        is_spatial: bool,
    ) -> Self {
        StageVisitor {
            visitors: TxHashMap::default(),
            ann_writer: None,
            ann_writer_num_reads: 0,
            metrics_sender,
            target_genes,
            is_spatial,
        }
    }

    pub(super) fn with_ann_writer_sample(
        metrics_sender: ShardSender<BarcodeMetrics, LibFeatThenBarcodeOrder>,
        target_genes: Option<HashSet<Gene>>,
        ann_writer: W,
        sample_rate: f32,
        rng: SmallRng,
        is_spatial: bool,
    ) -> Self {
        StageVisitor {
            visitors: TxHashMap::default(),
            ann_writer: Some((ann_writer, sample_rate, rng)),
            ann_writer_num_reads: 0,
            metrics_sender,
            target_genes,
            is_spatial,
        }
    }

    pub(super) fn ann_writer_num_reads(&self) -> usize {
        self.ann_writer_num_reads
    }

    pub(super) fn finish(mut self) -> Result<()> {
        if let Some((writer, _, _)) = self.ann_writer.take() {
            writer.finish()?;
        }
        for (_, v) in self.visitors {
            v.finish()?; // Finish will be called on drop, but we will miss any IO errors
        }
        self.metrics_sender.finished()?;
        Ok(())
    }

    pub(super) fn has_barcode_summaries(&self) -> bool {
        self.visitors
            .values()
            .any(|v| v.barcode_summaries.is_some())
    }

    /// Iterate over all of the BarcodeSummary in totally-sorted order.
    pub(super) fn sorted_barcode_summaries(&self) -> impl Iterator<Item = &BarcodeSummary> {
        self.visitors
            .iter()
            .sorted_by_key(|(library_type, _)| **library_type)
            .filter_map(|(_, v)| v.barcode_summaries.as_ref())
            // NOTE: this produces sorted output because barcode summaries for each library
            // type are stored in different visitors, so even though BarcodeSummary
            // sorts by LibraryType before barcode, and the order of each collection
            // of BarcodeSummary is sorted by barcode first, this still allows us
            // to produce a totally-sorted collection by first sorting the visitors
            // by library type, then flattening all of their iterators together.
            // We could replace this with Itertools::kmerge() which would ensure
            // that we maintain total ordering even if the BarcodeSummary
            // order criteria change, at the cost of a lot of extra comparisons.
            .flatten()
    }
}

impl<W> AnnotatedReadVisitor for StageVisitor<W>
where
    W: LazyWrite<ReadAnnotations, BufWriter<File>>,
{
    fn visit_read_annotation(&mut self, annotation: &ReadAnnotations) {
        if let Some((writer, rate, rng)) = self.ann_writer.as_mut()
            && rng.random_bool(*rate as f64)
        {
            self.ann_writer_num_reads += 1;
            writer
                .write_item(annotation)
                .expect("Error while writing ReadAnnotations.");
        }

        match self.visitors.entry(annotation.read.library_type) {
            Occupied(mut occ) => {
                occ.get_mut().visit_read_annotation(annotation);
            }
            Vacant(vac) => {
                vac.insert(AlignAndCountVisitor::new(
                    annotation.read.library_type,
                    self.metrics_sender.clone(),
                    self.target_genes.clone(),
                    self.is_spatial,
                ))
                .visit_read_annotation(annotation);
            }
        }
    }
}
