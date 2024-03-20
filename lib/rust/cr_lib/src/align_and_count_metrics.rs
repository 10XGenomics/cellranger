use crate::align_metrics::{AlignAndCountVisitor, BarcodeMetrics, LibFeatThenBarcodeOrder};
use crate::aligner::BarcodeSummary;
use anyhow::Result;
use cr_types::types::LibraryType;
use fxhash::FxHashMap;
use martian_filetypes::LazyWrite;
use rand::Rng;
use rand_chacha::ChaCha20Rng;
use shardio::ShardSender;
use std::collections::hash_map::Entry::{Occupied, Vacant};
use std::collections::HashSet;
use std::fs::File;
use std::io::BufWriter;
use transcriptome::Gene;
use tx_annotation::read::ReadAnnotations;
use tx_annotation::visitor::AnnotatedReadVisitor;

pub(crate) struct StageVisitor<W>
where
    W: LazyWrite<ReadAnnotations, BufWriter<File>>,
{
    visitors: FxHashMap<LibraryType, AlignAndCountVisitor>,
    ann_writer: Option<(W, f32, ChaCha20Rng)>,
    // Number of reads written to the `ann_writer`
    ann_writer_num_reads: usize,
    metrics_sender: ShardSender<BarcodeMetrics, LibFeatThenBarcodeOrder>,
    target_genes: Option<HashSet<Gene>>,
}

impl<W> StageVisitor<W>
where
    W: LazyWrite<ReadAnnotations, BufWriter<File>>,
{
    pub(crate) fn new(
        metrics_sender: ShardSender<BarcodeMetrics, LibFeatThenBarcodeOrder>,
        target_genes: Option<HashSet<Gene>>,
    ) -> Self {
        StageVisitor {
            visitors: FxHashMap::default(),
            ann_writer: None,
            ann_writer_num_reads: 0,
            metrics_sender,
            target_genes,
        }
    }

    pub(crate) fn with_ann_writer_sample(
        metrics_sender: ShardSender<BarcodeMetrics, LibFeatThenBarcodeOrder>,
        target_genes: Option<HashSet<Gene>>,
        ann_writer: W,
        sample_rate: f32,
        rng: ChaCha20Rng,
    ) -> Self {
        StageVisitor {
            visitors: FxHashMap::default(),
            ann_writer: Some((ann_writer, sample_rate, rng)),
            ann_writer_num_reads: 0,
            metrics_sender,
            target_genes,
        }
    }

    pub(crate) fn ann_writer_num_reads(&self) -> usize {
        self.ann_writer_num_reads
    }

    pub(crate) fn finish(mut self) -> Result<()> {
        if let Some((writer, _, _)) = self.ann_writer.take() {
            writer.finish()?;
        }
        for (_, v) in self.visitors {
            v.finish()?; // Finish will be called on drop, but we will miss any IO errors
        }
        self.metrics_sender.finished()?;
        Ok(())
    }

    pub(crate) fn barcode_summaries(&self) -> Vec<BarcodeSummary> {
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
        if let Some((ref mut writer, rate, rng)) = self.ann_writer.as_mut() {
            if rng.gen_range(0.0..1.0) < *rate {
                self.ann_writer_num_reads += 1;
                writer
                    .write_item(annotation)
                    .expect("Error while writing ReadAnnotations.");
            }
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
