use crate::count_matrix::CountMatrixFile;
use barcode::{barcode_string, Barcode};
use cr_types::reference::feature_reference::FeatureReference;
use cr_types::rna_read::{LegacyLibraryType, RnaRead};
use cr_types::types::FeatureBarcodeCount;
use cr_types::SampleAssignment;
use cr_websummary::multi::tables::SequencingMetricsTable;
use martian_derive::{martian_filetype, MartianStruct, MartianType};
use martian_filetypes::bin_file::BinaryFormat;
use martian_filetypes::json_file::JsonFormat;
use martian_filetypes::tabular_file::CsvFile;
use metric::TxHashMap;
use serde::{Deserialize, Serialize};
use std::borrow::Cow;
use std::marker::PhantomData;
use std::path::PathBuf;
use tx_annotation::read::ReadAnnotations;

#[derive(Serialize, Deserialize, Clone, MartianType)]
pub struct Primer {
    pub name: String,
    pub seq: String,
}

/// Aggregate barcode metrics for one barcode: aggregate_barcodes.csv
#[derive(Deserialize)]
pub struct AggregateBarcode {
    #[serde(with = "barcode_string")]
    pub(crate) barcode: Barcode,
    pub(crate) library_type: LegacyLibraryType,
    #[allow(dead_code)]
    umis: u64,
    #[allow(dead_code)]
    umi_corrected_reads: u64,
    #[allow(dead_code)]
    frac_corrected_reads: f64,
    #[allow(dead_code)]
    frac_total_reads: f64,
    pub(crate) frac_sample_reads: f64,
}

/// Produced by MULTI_WRITE_PER_SAMPLE_MATRICES, and consumed by STRUCTIFY_PER_SAMPLE_OUTS.
#[derive(Serialize, Deserialize, Clone, MartianStruct)]
pub struct SampleMatrices {
    pub sample: SampleAssignment,
    pub filtered_matrix_h5: CountMatrixFile,
    pub filtered_matrix_mex: PathBuf,
    pub raw_matrix_h5: CountMatrixFile,
    pub raw_matrix_mex: PathBuf,
    pub filtered_barcodes: CsvFile<()>,
    pub aggregate_barcodes: Option<CsvFile<AggregateBarcode>>,
}

// Shardio file containing BcUmiInfo records, sorted by barcode
martian_filetype!(BcUmiInfoShardFile, "bui");

/// Sort by barcode and then by feature.
pub struct BarcodeThenFeatureOrder<B>(PhantomData<B>);

impl<B: Copy + Ord> shardio::SortKey<FeatureBarcodeCount<B>> for BarcodeThenFeatureOrder<B> {
    type Key = (B, u32);

    fn sort_key(fbc: &FeatureBarcodeCount<B>) -> Cow<'_, Self::Key> {
        Cow::Owned((fbc.barcode, fbc.feature_idx))
    }
}

martian_filetype! { AnnSpillFile, "ann.spill" }
pub type AnnSpillFormat = BinaryFormat<AnnSpillFile, Vec<ReadAnnotations>>;

martian_filetype! { ReadSpillFile, "read.spill" }
pub type ReadSpillFormat = BinaryFormat<ReadSpillFile, Vec<RnaRead>>;

// Shardio file containing FeatureBarcodeCount records
martian_filetype!(CountShardFile, "csf");

// Shardio file containing BAM Record objects from rust_htslib,
// in position-sorted order
martian_filetype!(AlignShardFile, "asf");

martian_filetype!(ReadShardFile, "shard");
martian_filetype!(ReadPrefixCountFile, "rpc");
martian_filetype!(UmiCountFile, "umi");

// Shardio file storing per barcode metrics
martian_filetype!(BarcodeMetricsShardFile, "bmsf");

martian_filetype!(FastqFile, "fastq");
// Traditional BAM file
martian_filetype!(BamFile, "bam");
martian_filetype!(BcListFile, "bcl");

// Feature reference file, holding a bincode/serde serialized `FeatureReference`
martian_filetype!(_FeatureReferenceFile, "frf");
pub type FeatureReferenceFormat = BinaryFormat<_FeatureReferenceFile, FeatureReference>;

martian_filetype!(H5File, "h5");

martian_filetype!(_SequencingMetricsFile, "smf");
pub type SequencingMetricsFormat =
    JsonFormat<_SequencingMetricsFile, TxHashMap<LegacyLibraryType, SequencingMetricsTable>>;

martian_filetype!(SvgFile, "svg");
