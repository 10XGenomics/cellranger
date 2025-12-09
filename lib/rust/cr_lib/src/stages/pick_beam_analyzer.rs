//! Martian stage PICK_BEAM_ANALYZER
#![deny(missing_docs)]

use super::compute_antigen_vdj_metrics::AntigenVdjMetricsFormat;
use crate::utils::hard_link_martianfile;
use anyhow::Result;
use cr_types::SampleAssignment;
use cr_websummary::multi::antigen::AntigenSpecificityRow;
use martian::{MartianMain, MartianRover};
use martian_derive::{MartianStruct, make_mro};
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::tabular_file::CsvFile;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct BeamAnalyzerOutputs {
    pub antigen_specificity_scores: Option<CsvFile<AntigenSpecificityRow>>,
    pub antigen_assignment: Option<CsvFile<()>>,
    pub clonotype_concordance: Option<CsvFile<()>>,
    pub exact_subclonotype_concordance: Option<CsvFile<()>>,
    pub specificity_summary: Option<JsonFile<()>>,
    pub antigen_vdj_metrics_json: Option<JsonFile<()>>,
    pub antigen_vdj_metrics_bin: Option<AntigenVdjMetricsFormat>,
    pub per_barcode: Option<CsvFile<()>>,
}

impl BeamAnalyzerOutputs {
    fn hard_link(self, rover: &MartianRover) -> Result<Self> {
        let BeamAnalyzerOutputs {
            antigen_specificity_scores,
            antigen_assignment,
            clonotype_concordance,
            exact_subclonotype_concordance,
            specificity_summary,
            antigen_vdj_metrics_json,
            antigen_vdj_metrics_bin,
            per_barcode,
        } = self;
        Ok(BeamAnalyzerOutputs {
            antigen_specificity_scores: antigen_specificity_scores
                .map(|x| hard_link_martianfile(x, rover))
                .transpose()?,
            antigen_assignment: antigen_assignment
                .map(|x| hard_link_martianfile(x, rover))
                .transpose()?,
            clonotype_concordance: clonotype_concordance
                .map(|x| hard_link_martianfile(x, rover))
                .transpose()?,
            exact_subclonotype_concordance: exact_subclonotype_concordance
                .map(|x| hard_link_martianfile(x, rover))
                .transpose()?,
            specificity_summary: specificity_summary
                .map(|x| hard_link_martianfile(x, rover))
                .transpose()?,
            antigen_vdj_metrics_json: antigen_vdj_metrics_json
                .map(|x| hard_link_martianfile(x, rover))
                .transpose()?,
            antigen_vdj_metrics_bin: antigen_vdj_metrics_bin
                .map(|x| hard_link_martianfile(x, rover))
                .transpose()?,
            per_barcode: per_barcode
                .map(|x| hard_link_martianfile(x, rover))
                .transpose()?,
        })
    }

    fn is_empty(&self) -> bool {
        self.antigen_specificity_scores.is_none()
            && self.antigen_assignment.is_none()
            && self.clonotype_concordance.is_none()
            && self.exact_subclonotype_concordance.is_none()
            && self.specificity_summary.is_none()
            && self.antigen_vdj_metrics_json.is_none()
            && self.antigen_vdj_metrics_bin.is_none()
            && self.per_barcode.is_none()
    }
}

/// The Martian stage inputs.
#[derive(Clone, Deserialize, MartianStruct)]
pub struct PickBeamAnalyzerStageInputs {
    pub vdj_t: Option<HashMap<SampleAssignment, Option<BeamAnalyzerOutputs>>>,
    pub vdj_t_gd: Option<HashMap<SampleAssignment, Option<BeamAnalyzerOutputs>>>,
    pub vdj_b: Option<HashMap<SampleAssignment, Option<BeamAnalyzerOutputs>>>,
}

/// The Martian stage outputs.
#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct PickBeamAnalyzerStageOutputs {
    pub output: HashMap<SampleAssignment, Option<BeamAnalyzerOutputs>>,
}

/// Martian stage PICK_BEAM_ANALYZER
pub struct PickBeamAnalyzer;

#[make_mro(volatile = strict)]
impl MartianMain for PickBeamAnalyzer {
    type StageInputs = PickBeamAnalyzerStageInputs;
    type StageOutputs = PickBeamAnalyzerStageOutputs;

    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        let mut output: HashMap<SampleAssignment, Option<BeamAnalyzerOutputs>> = HashMap::new();
        let inputs: Vec<HashMap<SampleAssignment, Option<BeamAnalyzerOutputs>>> =
            vec![args.vdj_t, args.vdj_t_gd, args.vdj_b]
                .into_iter()
                .flatten()
                .collect();
        for sample in inputs[0].keys() {
            let mut options: Vec<_> = inputs
                .iter()
                .filter_map(|inputs| inputs[sample].clone())
                .filter(|beam_analyzer| !beam_analyzer.is_empty())
                .collect();
            assert!(options.len() <= 1);
            output.insert(
                sample.clone(),
                options
                    .pop()
                    .map(|beam_analyzer| beam_analyzer.hard_link(&rover))
                    .transpose()?,
            );
        }
        Ok(Self::StageOutputs { output })
    }
}
