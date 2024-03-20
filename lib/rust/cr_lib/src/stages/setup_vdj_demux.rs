//! SetupVDJDemux stage code

use crate::{GexMatrices, SampleMatrices};
use anyhow::{Context, Result};
use cr_types::{CrMultiGraph, FingerprintFile};
use martian::prelude::*;
use martian_derive::{make_mro, MartianStruct};
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::{FileTypeRead, FileTypeWrite};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct SetupVDJDemuxStageInputs {
    pub multi_matrices: Option<Vec<SampleMatrices>>,
    pub multi_graph: Option<JsonFile<CrMultiGraph>>,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct VdjDemuxSampleInfo {
    pub sample_id: String,
    pub sample_number: usize,
    pub gex_matrices: Option<GexMatrices>,
    pub fingerprint: Option<FingerprintFile>,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct SetupVDJDemuxStageOutputs {
    pub is_multi: bool,
    pub is_not_multi: bool,
    pub has_antigen: bool,
    pub per_sample_info: Option<HashMap<String, VdjDemuxSampleInfo>>,
}

pub struct SetupVDJDemux;

#[make_mro]
impl MartianMain for SetupVDJDemux {
    type StageInputs = SetupVDJDemuxStageInputs;
    type StageOutputs = SetupVDJDemuxStageOutputs;
    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        let multi_graph = match args.multi_graph {
            Some(multi_graph) => multi_graph.read()?,
            None => {
                return Ok(SetupVDJDemuxStageOutputs {
                    is_multi: false,
                    is_not_multi: true,
                    has_antigen: false,
                    per_sample_info: None,
                });
            }
        };

        let has_antigen = multi_graph.has_library_type(cr_types::LibraryType::FeatureBarcodes(
            cr_types::FeatureBarcodeType::Antigen,
        ));
        let mut per_sample_info = HashMap::new();
        for (sample_number, sample) in multi_graph.samples.into_iter().enumerate() {
            let fingerprint = sample
                .cell_multiplexing_type()
                .map(|_| {
                    rover
                        .make_path::<FingerprintFile>(format!("{}_fingerprint", sample.sample_id))
                        .with_content(&sample.fingerprints)
                })
                .transpose()?;
            per_sample_info.insert(
                sample.sample_id.clone(),
                VdjDemuxSampleInfo {
                    sample_id: sample.sample_id,
                    sample_number,
                    gex_matrices: None,
                    fingerprint,
                },
            );
        }

        if let Some(multi_matrices) = args.multi_matrices {
            for sample_matrices in multi_matrices {
                let sample_id = sample_matrices.sample.to_string();
                per_sample_info.get_mut(&sample_id).unwrap().gex_matrices = Some(
                    GexMatrices::from(sample_matrices)
                        .hardlink(&rover) // Hardlink to play nice with Martian
                        .context("Error hard linking SampleMatrices")?,
                );
            }
        }

        Ok(SetupVDJDemuxStageOutputs {
            is_multi: true,
            is_not_multi: false,
            has_antigen,
            per_sample_info: Some(per_sample_info),
        })
    }
}
