use crate::barcode_sort::ReadVisitor;
use crate::per_type_metric;
use anyhow::Result;
use barcode::{BarcodeConstructMetric, BarcodeSegmentState, HasBarcode};
use cr_types::rna_read::RnaRead;
use cr_types::types::LibraryFeatures;
use fastq_set::read_pair::ReadPair;
use json_report_derive::JsonReport;
use martian_derive::martian_filetype;
use martian_filetypes::bin_file::BinaryFormat;
use metric::{CountMetric, Metric, PercentMetric, TxHashMap};
use metric_derive::Metric;
use serde::{Deserialize, Serialize};

/// The metrics associated with BarcodeCorrection.
#[derive(Clone, Serialize, Deserialize, Metric, JsonReport)]
pub struct InnerBarcodeCorrectionMetrics {
    corrected_bc: PercentMetric,
    /// Fraction of corrected barcodes in each segment
    corrected_bc_in: BarcodeConstructMetric<PercentMetric>,
    good_bc: PercentMetric,
    good_bc_in: BarcodeConstructMetric<PercentMetric>,
}

impl InnerBarcodeCorrectionMetrics {
    /// Return a new BarcodeCorrectionMetrics initialized with the specified number of valid read pairs.
    pub fn with_valid(
        num_valid_read_pairs: i64,
        num_valid_segments: BarcodeConstructMetric<CountMetric>,
    ) -> Self {
        InnerBarcodeCorrectionMetrics {
            corrected_bc: (0, num_valid_read_pairs).into(),
            corrected_bc_in: num_valid_segments.map(|x| PercentMetric::from_parts(0.into(), x)),
            good_bc: (num_valid_read_pairs, num_valid_read_pairs).into(),
            good_bc_in: num_valid_segments.map(|x| PercentMetric::from_parts(x, x)),
        }
    }
}

per_type_metric!(
    BarcodeCorrectionMetrics,
    BarcodeCorrectionVisitor,
    LibraryFeatures,
    InnerBarcodeCorrectionMetrics
);

impl BarcodeCorrectionMetrics {
    pub fn with_valid(
        lib_feats: LibraryFeatures,
        num_valid_read_pairs: i64,
        num_valid_segments: BarcodeConstructMetric<CountMetric>,
    ) -> Self {
        let mut m = TxHashMap::default();
        m.insert(
            lib_feats,
            InnerBarcodeCorrectionMetrics::with_valid(num_valid_read_pairs, num_valid_segments),
        );
        BarcodeCorrectionMetrics(m)
    }
}

impl ReadVisitor for BarcodeCorrectionVisitor {
    type ReadType = RnaRead;

    // Visit reads and record barcode correction metrics.
    fn visit_processed_read(&mut self, read: &mut RnaRead) -> Result<()> {
        let metrics = InnerBarcodeCorrectionMetrics {
            // This visitor visits only uncorrected invalid barcodes.
            // The barcode is valid now if and only if it was corrected.
            corrected_bc: read.barcode().is_valid().into(),
            corrected_bc_in: read
                .segmented_barcode()
                .segments()
                .map(|seg| {
                    PercentMetric::from(seg.state == BarcodeSegmentState::ValidAfterCorrection)
                })
                .into(),
            good_bc: read.barcode().is_valid().into(),
            good_bc_in: read
                .segmented_barcode()
                .segments_valid()
                .map(PercentMetric::from)
                .into(),
        };

        self.metrics
            .0
            .entry(read.library_feats())
            .or_insert_with(Metric::new)
            .merge(metrics);

        Ok(())
    }

    fn visit_unprocessed_read(&mut self, _: ReadPair, _: String) -> Result<()> {
        unreachable!()
    }
}

martian_filetype! { BarcodeCorrectionMetricsFiletype, "bcm" }
pub type BarcodeCorrectionMetricsFormat =
    BinaryFormat<BarcodeCorrectionMetricsFiletype, BarcodeCorrectionMetrics>;
