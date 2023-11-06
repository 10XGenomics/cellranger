use barcode::Barcode;
use cr_types::rna_read::RnaRead;
use martian_derive::{martian_filetype, MartianStruct};
use martian_filetypes::bin_file::{BinaryFormat, BincodeFile};
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::lz4_file::Lz4;
use metric::SimpleHistogram;
use serde::{Deserialize, Serialize};
use std::path::PathBuf;
use vdj_reference::VdjReceptor;

martian_filetype!(BamFile, "bam");
martian_filetype!(BamBaiFile, "bam.bai");
martian_filetype!(FastaFile, "fasta");
martian_filetype!(FastqFile, "fastq");
martian_filetype!(_AsmReadsPerBcFile, "arp");
pub type AsmReadsPerBcFormat = BinaryFormat<_AsmReadsPerBcFile, SimpleHistogram<String>>;

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct AssemblyStageInputs {
    pub bc_sorted_rna_reads: Vec<Lz4<BincodeFile<Vec<RnaRead>>>>,
    pub paired_end: bool,
    pub vdj_reference_path: Option<PathBuf>,
    pub receptor: Option<VdjReceptor>,
    pub n50_n50_rpu: u32,
    pub npairs: u64,
    pub denovo: bool,
    pub inner_enrichment_primers: Option<PathBuf>,
    pub total_read_pairs: i64,
    pub corrected_bc_counts: JsonFile<SimpleHistogram<Barcode>>,
    pub r2_revcomp: bool,
}
