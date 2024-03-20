use anyhow::Result;
use barcode::Barcode;
use cr_types::chemistry::ChemistryDefs;
use cr_types::rna_read::RnaRead;
use itertools::Itertools;
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

#[derive(Clone, Debug, Serialize, Deserialize, Default)]
#[serde(try_from = "UmiListString", into = "UmiListString")]
pub struct UmiList(pub Vec<i32>);

#[derive(Clone, Debug, Serialize, Deserialize, Default)]
struct UmiListString(String);

#[derive(Clone, Debug, Serialize, Deserialize, Default)]
pub struct ContigSummaryRow {
    pub barcode: String,
    pub contig_name: String,
    pub num_reads: usize,
    pub num_pairs: usize,
    pub num_umis: usize,
    pub umi_list: UmiList,
}

impl TryFrom<UmiListString> for UmiList {
    type Error = anyhow::Error;

    fn try_from(value: UmiListString) -> Result<Self> {
        let umi_list: Result<Vec<_>, _> = value.0.split(',').map(str::parse).collect();
        Ok(Self(umi_list?))
    }
}

impl From<UmiList> for UmiListString {
    fn from(value: UmiList) -> Self {
        Self(value.0.iter().join(","))
    }
}

#[derive(Clone, Debug, Serialize, Deserialize, Default)]
pub struct UmiSummaryRow {
    pub barcode: String,
    pub umi_id: i32,
    pub umi: String,
    pub reads: usize,
    pub min_umi_reads: i32,
    pub good_umi: bool,
    pub contigs: String,
}

#[derive(Clone, Debug, Serialize, Deserialize, Default)]
pub struct BarcodeSupport {
    pub barcode: String,
    // num non-solo surviving umis i.e. xucounts supporting this barcode
    pub count: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct AssemblyStageInputs {
    pub chemistry_defs: ChemistryDefs,
    pub bc_sorted_rna_reads: Vec<Lz4<BincodeFile<Vec<RnaRead>>>>,
    pub vdj_reference_path: Option<PathBuf>,
    pub receptor: Option<VdjReceptor>,
    pub n50_n50_rpu: u32,
    pub npairs: u64,
    pub denovo: bool,
    pub inner_enrichment_primers: Option<PathBuf>,
    pub total_read_pairs: i64,
    pub corrected_bc_counts: JsonFile<SimpleHistogram<Barcode>>,
    pub min_contig_length: Option<usize>,
}
