//! MatchVdjOuts stage code

use crate::write_contig_proto::ProtoFile;
use anyhow::Result;
use martian::prelude::*;
use martian_derive::{make_mro, martian_filetype, MartianStruct};
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::tabular_file::{CsvFile, TsvFile};
use serde::{Deserialize, Serialize};
use std::path::PathBuf;
use vdj_reference::VdjReceptor;

martian_filetype!(FaFile, "fa");
martian_filetype!(FastaFile, "fasta");
martian_filetype!(VdjLoupeFile, "vloupe");
martian_filetype! {HtmlFile, "html"}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct VdjAggrResults {
    airr_rearrangement: TsvFile<()>,
    clonotypes: CsvFile<()>,
    donor_regions: FaFile,
    consensus_fasta: FastaFile,
    filtered_contig_annotations_csv: CsvFile<()>,
    consensus_annotations_csv: CsvFile<()>,
    web_summary_data: JsonFile<()>,
    vloupe: Option<VdjLoupeFile>,
    filter_summary: HtmlFile,
    enclone_output: ProtoFile,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct AntigenAggrResults {
    antigen_specificity_scores: Option<CsvFile<()>>,
    per_barcode_csv: Option<CsvFile<()>>,
}

// It would have been better if we received a single vec of a struct, but
// that does not play well with martian map call. So we will take each individual
// field as a vector. All the vectors in this struct should be the same length
#[derive(Debug, Clone, Deserialize, MartianStruct)]
pub struct MatchVdjOutsStageInputs {
    receptors: Vec<VdjReceptor>,
    clonotypes: Vec<CsvFile<()>>,
    donor_ref_fas: Vec<FaFile>,
    consensus_fastas: Vec<FastaFile>,
    vdj_reference_paths: Vec<PathBuf>, // Should be the same reference
    filtered_contig_annotations_csvs: Vec<CsvFile<()>>,
    consensus_annotations_csvs: Vec<CsvFile<()>>,
    web_summary_data: Vec<JsonFile<()>>,
    vloupes: Option<Vec<VdjLoupeFile>>, // TODO: Not an Option after vlconverter is in place
    antigen_analysis: Vec<Option<AntigenAggrResults>>,
    antigen_aggr_web_summary_data_in: Vec<Option<JsonFile<()>>>,
    airr_rearrangements: Vec<TsvFile<()>>,
    filter_summaries: Vec<HtmlFile>,
    enclone_outputs: Vec<ProtoFile>,
}

impl MatchVdjOutsStageInputs {
    fn pick_idx(&self, idx: Option<usize>) -> Option<VdjAggrResults> {
        idx.map(|i| VdjAggrResults {
            airr_rearrangement: self.airr_rearrangements[i].clone(),
            clonotypes: self.clonotypes[i].clone(),
            donor_regions: self.donor_ref_fas[i].clone(),
            consensus_fasta: self.consensus_fastas[i].clone(),
            filtered_contig_annotations_csv: self.filtered_contig_annotations_csvs[i].clone(),
            consensus_annotations_csv: self.consensus_annotations_csvs[i].clone(),
            web_summary_data: self.web_summary_data[i].clone(),
            vloupe: self.vloupes.as_ref().map(|v| v[i].clone()),
            filter_summary: self.filter_summaries[i].clone(),
            enclone_output: self.enclone_outputs[i].clone(),
        })
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct MatchVdjOutsStageOutputs {
    vdj_t_results: Option<VdjAggrResults>,
    vdj_t_gd_results: Option<VdjAggrResults>,
    vdj_b_results: Option<VdjAggrResults>,
    vdj_reference_path: Option<PathBuf>,
    antigen_results: Option<AntigenAggrResults>,
    antigen_aggr_web_summary_data: Option<JsonFile<()>>,
}

// This is our stage struct
pub struct MatchVdjOuts;

#[make_mro(stage_name = MATCH_VDJ_AGGR_OUTS)]
impl MartianMain for MatchVdjOuts {
    type StageInputs = MatchVdjOutsStageInputs;
    type StageOutputs = MatchVdjOutsStageOutputs;
    fn main(&self, args: Self::StageInputs, _rover: MartianRover) -> Result<Self::StageOutputs> {
        assert!(
            args.receptors.len() <= 2,
            "Cannot have more than 2 elements in aggr_results"
        );

        let tcr_idx = args.receptors.iter().position(|&r| r == VdjReceptor::TR);
        let tcrgd_idx = args.receptors.iter().position(|&r| r == VdjReceptor::TRGD);
        let ig_idx = args.receptors.iter().position(|&r| r == VdjReceptor::IG);

        let mut antigen_analysis: Vec<AntigenAggrResults> = args
            .antigen_analysis
            .iter()
            .filter_map(|aar| {
                aar.as_ref().and_then(|aar| {
                    if aar.antigen_specificity_scores.is_none() & aar.per_barcode_csv.is_none() {
                        None
                    } else {
                        Some(aar.clone())
                    }
                })
            })
            .collect();

        assert!(
            antigen_analysis.len() <= 1,
            "Cannot have more than 1 antigen analysis result"
        );

        let antigen_analysis = antigen_analysis.pop();

        let antigen_aggr_web_summary_data = args
            .antigen_aggr_web_summary_data_in
            .clone()
            .into_iter()
            .flatten()
            .next();

        Ok(MatchVdjOutsStageOutputs {
            vdj_t_results: args.pick_idx(tcr_idx),
            vdj_t_gd_results: args.pick_idx(tcrgd_idx),
            vdj_b_results: args.pick_idx(ig_idx),
            vdj_reference_path: args.vdj_reference_paths.first().cloned(),
            antigen_results: antigen_analysis,
            antigen_aggr_web_summary_data,
        })
    }
}
