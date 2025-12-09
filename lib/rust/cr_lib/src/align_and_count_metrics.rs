#![deny(missing_docs)]

use crate::align_metrics::{AlignAndCountVisitor, BarcodeMetrics, LibFeatThenBarcodeOrder};
use crate::aligner::BarcodeSummary;
use anyhow::Result;
use cr_types::types::LibraryType;
use martian_filetypes::LazyWrite;
use metric::TxHashMap;
use rand::Rng;
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
}

impl<W> StageVisitor<W>
where
    W: LazyWrite<ReadAnnotations, BufWriter<File>>,
{
    pub(super) fn new(
        metrics_sender: ShardSender<BarcodeMetrics, LibFeatThenBarcodeOrder>,
        target_genes: Option<HashSet<Gene>>,
    ) -> Self {
        StageVisitor {
            visitors: TxHashMap::default(),
            ann_writer: None,
            ann_writer_num_reads: 0,
            metrics_sender,
            target_genes,
        }
    }

    pub(super) fn with_ann_writer_sample(
        metrics_sender: ShardSender<BarcodeMetrics, LibFeatThenBarcodeOrder>,
        target_genes: Option<HashSet<Gene>>,
        ann_writer: W,
        sample_rate: f32,
        rng: SmallRng,
    ) -> Self {
        StageVisitor {
            visitors: TxHashMap::default(),
            ann_writer: Some((ann_writer, sample_rate, rng)),
            ann_writer_num_reads: 0,
            metrics_sender,
            target_genes,
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

    pub(super) fn barcode_summaries(&self) -> Vec<BarcodeSummary> {
        self.visitors
            .values()
            .flat_map(|v| v.barcode_summaries.iter().cloned())
            .collect()
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
                ))
                .visit_read_annotation(annotation);
            }
        }
    }
}
