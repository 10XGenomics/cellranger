#![deny(missing_docs)]
use barcode::{BarcodeConstructMetric, BarcodeSegmentState};
use cr_types::rna_read::RnaRead;
use cr_types::types::LibraryType;
use json_report_derive::JsonReport;
use martian_derive::martian_filetype;
use martian_filetypes::bin_file::BinaryFormat;
use metric::{CountMetric, Metric, PercentMetric, TxHashMap};
use metric_derive::Metric;
use serde::{Deserialize, Serialize};

/// The metrics associated with BarcodeCorrection.
#[derive(Default, Serialize, Deserialize, Metric, JsonReport)]
pub struct InnerBarcodeCorrectionMetrics {
    /// Fraction of corrected barcodes.
    corrected_bc: PercentMetric,
    /// Fraction of corrected barcodes in each segment.
    corrected_bc_in: BarcodeConstructMetric<PercentMetric>,
    /// Fraction of good barcodes after correction.
    good_bc: PercentMetric,
    /// Fraction of good barcodes after correction in each segment.
    good_bc_in: BarcodeConstructMetric<PercentMetric>,
    /// Fraction of barcodes that could not be corrected, because the correction was ambiguous.
    ambiguous_bc_in: BarcodeConstructMetric<PercentMetric>,
    /// Fraction of barcodes that could not be corrected, because the correction was an indel.
    indel_bc_in: BarcodeConstructMetric<PercentMetric>,
    /// Fraction of barcodes that could not be corrected, because the molecule has unusual structure.
    unusual_bc_in: BarcodeConstructMetric<PercentMetric>,
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
            ambiguous_bc_in: num_valid_segments.map(|x| PercentMetric::from_parts(0.into(), x)),
            indel_bc_in: num_valid_segments.map(|x| PercentMetric::from_parts(0.into(), x)),
            unusual_bc_in: num_valid_segments.map(|x| PercentMetric::from_parts(0.into(), x)),
        }
    }
}

#[derive(Default, Serialize, Deserialize, Metric)]
pub struct BarcodeCorrectionMetrics(pub TxHashMap<LibraryType, InnerBarcodeCorrectionMetrics>);

impl BarcodeCorrectionMetrics {
    pub fn with_valid(
        library_type: LibraryType,
        num_valid_read_pairs: i64,
        num_valid_segments: BarcodeConstructMetric<CountMetric>,
    ) -> Self {
        let mut m = TxHashMap::default();
        m.insert(
            library_type,
            InnerBarcodeCorrectionMetrics::with_valid(num_valid_read_pairs, num_valid_segments),
        );
        BarcodeCorrectionMetrics(m)
    }

    /// Visit reads and record barcode correction metrics.
    pub fn visit_processed_read(&mut self, read: &mut RnaRead) {
        let metrics = InnerBarcodeCorrectionMetrics {
            // This visitor visits only uncorrected invalid barcodes.
            // The barcode is valid now if and only if it was corrected.
            corrected_bc: read.barcode_is_valid().into(),
            corrected_bc_in: read
                .segmented_barcode
                .segments()
                .map(|seg| {
                    PercentMetric::from(seg.state == BarcodeSegmentState::ValidAfterCorrection)
                })
                .into(),
            good_bc: read.barcode_is_valid().into(),
            good_bc_in: read
                .segmented_barcode
                .segments_valid()
                .map(PercentMetric::from)
                .into(),
            ambiguous_bc_in: read
                .segmented_barcode
                .segments()
                .map(|seg| PercentMetric::from(seg.state == BarcodeSegmentState::InvalidAmbiguous))
                .into(),
            indel_bc_in: read
                .segmented_barcode
                .segments()
                .map(|seg| PercentMetric::from(seg.state == BarcodeSegmentState::InvalidIndel))
                .into(),
            unusual_bc_in: read
                .segmented_barcode
                .segments()
                .map(|seg| PercentMetric::from(seg.state == BarcodeSegmentState::InvalidStructure))
                .into(),
        };

        self.0.entry(read.library_type).or_default().merge(metrics);
    }
}

martian_filetype! { BarcodeCorrectionMetricsFiletype, "bcm" }
pub type BarcodeCorrectionMetricsFormat =
    BinaryFormat<BarcodeCorrectionMetricsFiletype, BarcodeCorrectionMetrics>;

/// Barcode diversity metrics
#[derive(JsonReport)]
pub(super) struct BarcodeDiversityMetrics {
    /// Number of barcodes detected
    pub barcodes_detected: usize,

    /// Effective barcode diversity
    pub effective_barcode_diversity: f64,
}
