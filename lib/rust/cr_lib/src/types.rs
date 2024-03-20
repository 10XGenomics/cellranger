use crate::probe_barcode_matrix::ProbeCounts;
use crate::utils::hard_link_martianfile;
use anyhow::Error;
use barcode::{barcode_string, Barcode};
use cr_h5::count_matrix::CountMatrixFile;
use cr_types::reference::feature_reference::FeatureReference;
use cr_types::rna_read::RnaRead;
use cr_types::{LibraryType, SampleAssignment};
use cr_websummary::multi::tables::SequencingMetricsTable;
use martian::MartianRover;
use martian_derive::{martian_filetype, MartianStruct, MartianType};
use martian_filetypes::bin_file::BinaryFormat;
use martian_filetypes::json_file::JsonFormat;
use martian_filetypes::tabular_file::CsvFile;
use metric::TxHashMap;
use serde::{Deserialize, Serialize};
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
    pub(crate) library_type: LibraryType,
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
#[derive(Serialize, Deserialize, Clone, MartianStruct, Debug)]
pub struct SampleMatrices {
    pub sample: SampleAssignment,
    pub filtered_matrix_h5: CountMatrixFile,
    pub filtered_matrix_mex: PathBuf,
    pub raw_matrix_h5: CountMatrixFile,
    pub raw_probe_bc_matrix: Option<CountMatrixFile>,
    pub raw_matrix_mex: PathBuf,
    pub filtered_barcodes: CsvFile<()>,
    pub aggregate_barcodes: Option<CsvFile<AggregateBarcode>>,
    pub per_probe_metrics: Option<CsvFile<ProbeCounts>>,
}

/// Produced by SETUP_VDJ_DEMUX, and consumed by VDJ_ANALYSIS.
#[derive(Serialize, Deserialize, Clone, MartianStruct, Debug, Default)]
pub struct GexMatrices {
    pub filtered_matrix_h5: Option<CountMatrixFile>,
    pub raw_matrix_h5: Option<CountMatrixFile>,
    pub filtered_barcodes: Option<CsvFile<()>>,
}

impl GexMatrices {
    pub fn hardlink(&self, rover: &MartianRover) -> Result<Self, Error> {
        let mut hardlinks = self.clone();
        if let Some(filtered_matrix_h5) = hardlinks.filtered_matrix_h5 {
            hardlinks.filtered_matrix_h5 = Some(hard_link_martianfile(filtered_matrix_h5, rover)?);
        }
        if let Some(raw_matrix_h5) = hardlinks.raw_matrix_h5 {
            hardlinks.raw_matrix_h5 = Some(hard_link_martianfile(raw_matrix_h5, rover)?);
        }
        if let Some(filtered_barcodes) = hardlinks.filtered_barcodes {
            hardlinks.filtered_barcodes = Some(hard_link_martianfile(filtered_barcodes, rover)?);
        }
        Ok(hardlinks)
    }
}

impl From<SampleMatrices> for GexMatrices {
    fn from(value: SampleMatrices) -> Self {
        GexMatrices {
            filtered_matrix_h5: Some(value.filtered_matrix_h5),
            raw_matrix_h5: Some(value.raw_matrix_h5),
            filtered_barcodes: Some(value.filtered_barcodes),
        }
    }
}

// Shardio file containing BcUmiInfo records, sorted by barcode
martian_filetype!(BcUmiInfoShardFile, "bui");

martian_filetype! { AnnSpillFile, "ann.spill" }
pub type AnnSpillFormat = BinaryFormat<AnnSpillFile, Vec<ReadAnnotations>>;

martian_filetype! { ReadSpillFile, "read.spill" }
pub type ReadSpillFormat = BinaryFormat<ReadSpillFile, Vec<RnaRead>>;

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

martian_filetype!(_SequencingMetricsFile, "smf");
pub type SequencingMetricsFormat =
    JsonFormat<_SequencingMetricsFile, TxHashMap<LibraryType, SequencingMetricsTable>>;

martian_filetype!(SvgFile, "svg");
