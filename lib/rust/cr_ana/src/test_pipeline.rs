//! Secondary analysis pipeline for testing

use crate::stage_testing::run_stage;
use crate::stages::pca::{PcaStage, PcaStageInputs};
use crate::stages::umap::{UmapImplementation, UmapStage, UmapStageInputs};
use crate::types::H5File;
use anyhow::Result;
use martian::prelude::*;
use martian_derive::{make_mro, MartianStruct};
use serde::{Deserialize, Serialize};
use std::path::{Path, PathBuf};

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
struct SecondaryTestInputs {
    matrix_h5: H5File,
    num_pca_genes: Option<usize>,
    num_principal_comps: Option<usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
struct SecondaryTestOutputs {
    umap_csv: PathBuf,
    umap: H5File,
    pca: H5File,
}

struct SecondaryTest;

#[make_mro]
impl MartianMain for SecondaryTest {
    type StageInputs = SecondaryTestInputs;
    type StageOutputs = SecondaryTestOutputs;
    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        // -----------------------------------------------------------------------------------
        // RUN_PCA_NG
        println!("run pca");
        let pca_args = PcaStageInputs {
            matrix_h5: args.matrix_h5.clone(),
            num_pca_genes: None,
            num_principal_comps: None,
            is_spatial: false,
            pca_map: None,
        };
        let pca_outs = run_stage(PcaStage, pca_args, &rover)?;

        println!("run umap");
        let umap_args = UmapStageInputs {
            matrix_h5: args.matrix_h5.clone(),
            pca_h5: pca_outs.pca_h5,
            random_seed: None,
            n_neighbors: None,
            input_pcs: None,
            max_dims: None,
            min_dist: None,
            metric: None,
            implementation: UmapImplementation::Parallel,
        };

        let umap_outs = run_stage(UmapStage, umap_args, &rover)?;

        Ok(SecondaryTestOutputs {
            umap_csv: umap_outs.umap_csv,
            umap: umap_outs.umap_h5,
            pca: args.matrix_h5,
        })
    }
}

#[cfg(test)]
fn run_pipeline_test(args: SecondaryTestInputs, path: &Path) -> Result<()> {
    // let dir = tempfile::tempdir()?;
    let _outs = SecondaryTest.test_run(path, args)?;
    Ok(())
}

#[test]
fn test_first() -> Result<()> {
    let matrix = "filtered_feature_bc_matrix.h5";

    let args = SecondaryTestInputs {
        matrix_h5: matrix.into(),
        num_pca_genes: None,
        num_principal_comps: None,
    };

    let path = "umap_test_pipeline_1419427_1460102";
    std::fs::create_dir(path)?;
    run_pipeline_test(args, Path::new(path))
}
