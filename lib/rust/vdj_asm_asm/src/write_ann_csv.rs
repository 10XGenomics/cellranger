//!
//! WriteAnnCsv stage produces the CS outputs `all_contig_annotations.csv` and
//! `filtered_contig_annotations.csv` from the contig annotations json file.
//!

use anyhow::Result;
use martian::prelude::*;
use martian_derive::{make_mro, MartianStruct};
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::tabular_file::CsvFile;
use martian_filetypes::{LazyFileTypeIO, LazyWrite};
use serde::{Deserialize, Serialize};
use vdj_ann::annotate::ContigAnnotation;
use vdj_types::VdjRegion;

// One row in contig annotations csv file.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ContigAnnotationCsvRow {
    pub barcode: String,
    pub is_cell: bool,
    pub contig_id: String,
    pub high_confidence: bool,
    pub length: usize,
    pub chain: String,
    pub v_gene: Option<String>,
    pub d_gene: Option<String>,
    pub j_gene: Option<String>,
    pub c_gene: Option<String>,
    pub full_length: bool,
    pub productive: bool,
    pub fwr1: Option<String>,
    pub fwr1_nt: Option<String>,
    pub cdr1: Option<String>,
    pub cdr1_nt: Option<String>,
    pub fwr2: Option<String>,
    pub fwr2_nt: Option<String>,
    pub cdr2: Option<String>,
    pub cdr2_nt: Option<String>,
    pub fwr3: Option<String>,
    pub fwr3_nt: Option<String>,
    pub cdr3: Option<String>,
    pub cdr3_nt: Option<String>,
    pub fwr4: Option<String>,
    pub fwr4_nt: Option<String>,
    pub reads: usize,
    pub umis: usize,
    pub raw_clonotype_id: Option<String>,
    pub raw_consensus_id: Option<String>,
    pub exact_subclonotype_id: Option<String>,
}

impl From<&ContigAnnotation> for ContigAnnotationCsvRow {
    fn from(ann: &ContigAnnotation) -> Self {
        ContigAnnotationCsvRow {
            barcode: ann.barcode.clone(),
            is_cell: ann.is_cell,
            contig_id: ann.contig_name.clone(),
            high_confidence: ann.high_confidence,
            length: ann.sequence.len(),
            chain: ann
                .chain_type()
                .map_or_else(|| (*"None").to_string(), |c| c.to_string()),
            v_gene: ann.get_gene_name(VdjRegion::V).cloned(),
            d_gene: ann.get_gene_name(VdjRegion::D).cloned(),
            j_gene: ann.get_gene_name(VdjRegion::J).cloned(),
            c_gene: ann.get_gene_name(VdjRegion::C).cloned(),
            full_length: ann.is_full_length(),
            productive: ann.is_productive(),
            fwr1: ann.fwr1.as_ref().map(|r| r.aa_seq.clone()),
            fwr1_nt: ann.fwr1.as_ref().map(|r| r.nt_seq.clone()),
            cdr1: ann.cdr1.as_ref().map(|r| r.aa_seq.clone()),
            cdr1_nt: ann.cdr1.as_ref().map(|r| r.nt_seq.clone()),
            fwr2: ann.fwr2.as_ref().map(|r| r.aa_seq.clone()),
            fwr2_nt: ann.fwr2.as_ref().map(|r| r.nt_seq.clone()),
            cdr2: ann.cdr2.as_ref().map(|r| r.aa_seq.clone()),
            cdr2_nt: ann.cdr2.as_ref().map(|r| r.nt_seq.clone()),
            fwr3: ann.fwr3.as_ref().map(|r| r.aa_seq.clone()),
            fwr3_nt: ann.fwr3.as_ref().map(|r| r.nt_seq.clone()),
            cdr3: ann.cdr3.clone(),
            cdr3_nt: ann.cdr3_seq.clone(),
            fwr4: ann.fwr4.as_ref().map(|r| r.aa_seq.clone()),
            fwr4_nt: ann.fwr4.as_ref().map(|r| r.nt_seq.clone()),
            reads: ann.read_count,
            umis: ann.umi_count,
            raw_clonotype_id: ann.info.raw_clonotype_id.clone(),
            raw_consensus_id: ann.info.raw_consensus_id.clone(),
            exact_subclonotype_id: ann.info.exact_subclonotype_id.clone(),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct WriteAnnCsvStageInputs {
    pub all_contig_annotations_json: JsonFile<Vec<ContigAnnotation>>,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct WriteAnnCsvStageOutputs {
    pub all_contig_annotations_csv: CsvFile<ContigAnnotationCsvRow>,
    pub filtered_contig_annotations_csv: CsvFile<ContigAnnotationCsvRow>,
}

pub struct WriteAnnCsv;

#[make_mro]
impl MartianMain for WriteAnnCsv {
    type StageInputs = WriteAnnCsvStageInputs;
    type StageOutputs = WriteAnnCsvStageOutputs;
    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        // Define contig_ann reader
        let contig_reader = args.all_contig_annotations_json.lazy_reader()?;

        let all_contig_annotations_csv: CsvFile<_> = rover.make_path("all_contig_annotations");
        let mut all_writer = all_contig_annotations_csv.lazy_writer()?;

        let filtered_contig_annotations_csv: CsvFile<_> =
            rover.make_path("filtered_contig_annotations");
        let mut filt_writer = filtered_contig_annotations_csv.lazy_writer()?;

        for ann in contig_reader {
            let ann: ContigAnnotation = ann?;
            let row = ContigAnnotationCsvRow::from(&ann);

            all_writer.write_item(&row)?;
            if row.productive && row.is_cell && row.high_confidence {
                filt_writer.write_item(&row)?;
            }
        }

        all_writer.finish()?;
        filt_writer.finish()?;

        Ok(WriteAnnCsvStageOutputs {
            all_contig_annotations_csv,
            filtered_contig_annotations_csv,
        })
    }
}
