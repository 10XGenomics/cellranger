//! PCA stage code

use crate::io::{csv, h5};
use crate::pca::run_pca;
use crate::types::H5File;
use crate::EXCLUDED_FEATURE_TYPES;
use anyhow::Result;
use cr_types::reference::feature_reference::FeatureType;
use hdf5_io::matrix::read_adaptive_csr_matrix;
use martian::prelude::{MartianRover, MartianStage, Resource, StageDef};
use martian_derive::{make_mro, MartianStruct};
use scan_types::matrix::AdaptiveFeatureBarcodeMatrix as FBM;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::path::PathBuf;

const PCA_THRESHOLD: f64 = 3.0;
const PCA_COMPONENTS: usize = 10;

#[derive(Clone, Debug, Serialize, Deserialize, MartianStruct)]
pub struct PcaOutputs {
    pub pca_h5: H5File,
    pub pca_csv: PathBuf,
}

#[derive(Clone, Debug, Deserialize, MartianStruct)]
pub struct PcaStageInputs {
    pub matrix_h5: H5File,
    pub num_pca_genes: Option<usize>,
    pub num_principal_comps: Option<usize>,
    pub is_spatial: bool,
    pub pca_map: Option<HashMap<String, PcaOutputs>>,
}

#[derive(Clone, Debug, Serialize, Deserialize, MartianStruct)]
pub struct PcaStageOutputs {
    pub pca_h5: H5File,
    pub pca_csv: PathBuf,
}

#[derive(Clone, Debug, Serialize, Deserialize, MartianStruct)]
pub struct PcaChunkInputs {
    feature_type: FeatureType,
}

pub struct PcaStage;

#[make_mro(stage_name = RUN_PCA_NG, volatile = strict)]
impl MartianStage for PcaStage {
    type StageInputs = PcaStageInputs;
    type StageOutputs = PcaStageOutputs;
    type ChunkInputs = PcaChunkInputs;
    type ChunkOutputs = PcaStageOutputs;

    fn split(
        &self,
        args: Self::StageInputs,
        _rover: MartianRover,
    ) -> Result<StageDef<Self::ChunkInputs>> {
        let num_pcs = args.num_principal_comps.unwrap_or(PCA_COMPONENTS);
        let mem_gib = 3
            + (h5::estimate_mem_gib_from_nnz(&args.matrix_h5)?.ceil()
                + h5::estimate_mem_gib_for_pca_outs(&args.matrix_h5, num_pcs)?.ceil())
                as isize;
        let pca_map = args.pca_map.unwrap_or_default();

        Ok(h5::matrix_feature_types(&args.matrix_h5)?
            .into_iter()
            .filter(|&(feature_type, count)| {
                count >= 2
                    && !EXCLUDED_FEATURE_TYPES.contains(&feature_type)
                    && !pca_map.contains_key(&feature_type.to_string())
            })
            .map(|(feature_type, _)| {
                (
                    Self::ChunkInputs { feature_type },
                    Resource::with_mem_gb(mem_gib).threads(1),
                )
            })
            .collect())
    }

    fn main(
        &self,
        args: Self::StageInputs,
        chunk_args: Self::ChunkInputs,
        rover: MartianRover,
    ) -> Result<Self::ChunkOutputs> {
        rayon::ThreadPoolBuilder::new()
            .num_threads(rover.get_threads())
            .build_global()?;
        let retained = Some(chunk_args.feature_type.to_string());
        let FBM {
            barcodes,
            feature_ids,
            matrix,
            ..
        } = read_adaptive_csr_matrix(&args.matrix_h5, retained.as_deref(), Some(0))?.0;
        let max_features = args.num_pca_genes.unwrap_or_else(|| matrix.rows());
        let num_pcs = args.num_principal_comps.unwrap_or(PCA_COMPONENTS);
        let result = run_pca(
            &matrix,
            &feature_ids,
            chunk_args.feature_type,
            max_features,
            num_pcs,
            PCA_THRESHOLD,
            args.is_spatial,
        )
        .or_else(|_| {
            log::warn!(
                "PCA with threshold = {} failed, reattempting without threshold",
                PCA_THRESHOLD
            );
            run_pca(
                &matrix,
                &feature_ids,
                chunk_args.feature_type,
                max_features,
                num_pcs,
                0.0,
                args.is_spatial,
            )
        })?;

        let pca_h5 = rover.make_path("pca_h5");
        h5::save_pca(&pca_h5, &result)?;

        let pca_csv: PathBuf = rover.make_path("pca_csv");
        csv::save_pca(&pca_csv, &result, &barcodes, &feature_ids)?;

        Ok(Self::StageOutputs { pca_h5, pca_csv })
    }

    fn join(
        &self,
        args: Self::StageInputs,
        _chunk_defs: Vec<Self::ChunkInputs>,
        chunk_outs: Vec<Self::ChunkOutputs>,
        rover: MartianRover,
    ) -> Result<Self::StageOutputs> {
        let pca_h5 = rover.make_path("pca_h5");
        h5::combine_pcas(
            &pca_h5,
            chunk_outs.iter().map(|c| &c.pca_h5).chain(
                args.pca_map
                    .iter()
                    .flat_map(|m| m.values().map(|c| &c.pca_h5)),
            ),
        )?;

        let pca_csv: PathBuf = rover.make_path("pca_csv");
        csv::combine_pcas(
            &pca_csv,
            chunk_outs.iter().map(|c| &c.pca_csv).chain(
                args.pca_map
                    .iter()
                    .flat_map(|m| m.values().map(|c| &c.pca_csv)),
            ),
        )?;

        Ok(Self::StageOutputs { pca_h5, pca_csv })
    }
}
