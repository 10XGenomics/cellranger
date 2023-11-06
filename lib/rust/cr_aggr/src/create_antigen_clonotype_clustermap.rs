//! CreateClonotypeClustermap stage code

use anyhow::Result;
use cr_websummary::multi::antigen::{clonotype_specificity_heatmap, AntigenSpecificityRow};
use cr_websummary::RawChartWithHelp;
use martian::prelude::*;
use martian_derive::{make_mro, MartianStruct};
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::tabular_file::CsvFile;
use martian_filetypes::FileTypeWrite;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Deserialize, MartianStruct)]
pub struct CreateClonotypeClustermapInputs {
    antigen_specificity: Option<CsvFile<AntigenSpecificityRow>>,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct CreateClonotypeClustermapOutputs {
    antigen_clonotype_clustermap: Option<JsonFile<RawChartWithHelp>>,
}

// This is our stage struct
pub struct CreateClonotypeClustermap;

#[make_mro(stage_name = CREATE_CLONOTYPE_CLUSTERMAP)]
impl MartianMain for CreateClonotypeClustermap {
    type StageInputs = CreateClonotypeClustermapInputs;
    type StageOutputs = CreateClonotypeClustermapOutputs;
    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        let clonotype_clustermap = args
            .antigen_specificity
            .map(clonotype_specificity_heatmap)
            .transpose()?
            .flatten();

        if let Some(clon) = clonotype_clustermap {
            let antigen_clonotype_clustermap: JsonFile<_> =
                rover.make_path("antigen_clonotype_clustermap.json");
            antigen_clonotype_clustermap.write(&clon)?;
            Ok(CreateClonotypeClustermapOutputs {
                antigen_clonotype_clustermap: Some(antigen_clonotype_clustermap),
            })
        } else {
            Ok(CreateClonotypeClustermapOutputs {
                antigen_clonotype_clustermap: None,
            })
        }
    }
}
