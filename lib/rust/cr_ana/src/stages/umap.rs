//! Martian stage RUN_UMAP

use crate::io::{csv, h5};
use crate::types::{EmbeddingResult, EmbeddingType, H5File};
use crate::EXCLUDED_FEATURE_TYPES;
use anyhow::Result;
use cr_types::reference::feature_reference::FeatureType;
use hdf5_io::matrix::get_barcodes_between;
use martian::prelude::{MartianRover, MartianStage, Resource, StageDef};
use martian_derive::{make_mro, MartianStruct, MartianType};
use serde::{Deserialize, Serialize};
use std::path::PathBuf;
use umap_rs::dist::DistanceType;
use umap_rs::umap::Umap;

const N_NEIGHBORS: usize = 30;
const MIN_DIST: f64 = 0.3;
const N_COMPONENTS: usize = 2;

#[derive(Clone, Debug, Deserialize, MartianStruct)]
pub struct UmapStageInputs {
    pub matrix_h5: H5File,
    pub pca_h5: H5File,
    pub random_seed: Option<u64>,
    pub n_neighbors: Option<usize>,
    pub input_pcs: Option<usize>,
    pub max_dims: Option<usize>,
    pub min_dist: Option<f64>,
    pub metric: Option<String>,
    pub implementation: UmapImplementation,
}

#[derive(Serialize, Deserialize, MartianType, Clone, Debug)]
#[serde(rename_all = "snake_case")]
pub enum UmapImplementation {
    Original,
    Parallel,
}

use UmapImplementation::{Original, Parallel};

#[derive(Clone, Debug, Serialize, MartianStruct)]
pub struct UmapStageOutputs {
    pub umap_h5: H5File,
    pub umap_csv: PathBuf,
}

#[derive(Clone, Debug, Serialize, Deserialize, MartianStruct)]
pub struct UmapChunkInputs {
    umap_dims: usize,
    feature_type: FeatureType,
}

#[derive(Clone, Debug, Serialize, Deserialize, MartianStruct)]
pub struct UmapChunkOutputs {
    umap_h5: H5File,
    umap_csv: PathBuf,
}

pub struct UmapStage;

#[make_mro(stage_name = RUN_UMAP, volatile = strict)]
impl MartianStage for UmapStage {
    type StageInputs = UmapStageInputs;
    type StageOutputs = UmapStageOutputs;
    type ChunkInputs = UmapChunkInputs;
    type ChunkOutputs = UmapChunkOutputs;

    fn split(
        &self,
        args: Self::StageInputs,
        _rover: MartianRover,
    ) -> Result<StageDef<Self::ChunkInputs>> {
        let mut stage_def = StageDef::new();
        let (_nfeatures, ncells, _nnz) = h5::matrix_shape(&args.matrix_h5)?;
        let feature_types = h5::matrix_feature_types(&args.matrix_h5)?;

        let threads = match args.implementation {
            Original => 5,
            Parallel => 4,
        };

        for umap_dims in N_COMPONENTS..=args.max_dims.unwrap_or(N_COMPONENTS) {
            for (&feature_type, &count) in &feature_types {
                // if we have only 1 feature, skip!
                if count < 2 || EXCLUDED_FEATURE_TYPES.contains(&feature_type) {
                    continue;
                }
                let chunk_inputs = Self::ChunkInputs {
                    umap_dims,
                    feature_type,
                };
                let chunk_mem_gb = ((ncells as f64 / 125_000.0).ceil() as isize).max(4);
                let chunk_vmem_gb = chunk_mem_gb * 8;
                let chunk_resource = Resource::new()
                    .threads(threads)
                    .mem_gb(chunk_mem_gb)
                    .vmem_gb(chunk_vmem_gb);
                stage_def.add_chunk_with_resource(chunk_inputs, chunk_resource);
            }
        }
        Ok(stage_def)
    }

    fn main(
        &self,
        args: Self::StageInputs,
        chunk_args: Self::ChunkInputs,
        rover: MartianRover,
    ) -> Result<Self::ChunkOutputs> {
        let (proj, barcodes) = {
            let proj = h5::load_transformed_pca_matrix(
                &args.pca_h5,
                chunk_args.feature_type,
                args.input_pcs,
            )?;
            let matrix = hdf5::File::open(&args.matrix_h5)?.group("matrix")?;
            let barcodes = get_barcodes_between(0, None, &matrix)?;
            (proj, barcodes)
        };
        let (num_bcs, _) = proj.dim();

        let metric = match args.metric.as_deref() {
            Some("correlation") | Some("pearson") | None => DistanceType::pearson(),
            Some("cosine") => DistanceType::cosine(),
            Some("euclidean") => DistanceType::euclidean(),
            Some(metric) => unimplemented!("unsupported metric: {}", metric),
        };
        let spread = 1.0;
        let umap = Umap::new(
            Some(metric),
            chunk_args.umap_dims,
            args.min_dist.unwrap_or(MIN_DIST),
            spread,
            args.n_neighbors.unwrap_or(N_NEIGHBORS).min(num_bcs - 1),
            None,
        );

        let embedding = match args.implementation {
            Parallel => {
                let mut state =
                    umap.initialize_fit_parallelized(&proj, args.random_seed, rover.get_threads());
                state.optimize_multithreaded(rover.get_threads());
                state.embedding
            }
            Original => {
                let mut state = umap.initialize_fit(&proj, args.random_seed, rover.get_threads());
                state.optimize();
                state.embedding
            }
        };

        let result = EmbeddingResult::new(
            &barcodes,
            embedding,
            EmbeddingType::Umap,
            chunk_args.feature_type,
        );

        let umap_h5 = rover.make_path("umap_h5");
        h5::save_embedding(&umap_h5, &result)?;

        let umap_csv: PathBuf = rover.make_path("umap_csv");
        csv::save_embedding(&umap_csv, &result)?;

        Ok(Self::ChunkOutputs { umap_h5, umap_csv })
    }

    fn join(
        &self,
        _args: Self::StageInputs,
        _chunk_defs: Vec<Self::ChunkInputs>,
        chunk_outs: Vec<Self::ChunkOutputs>,
        rover: MartianRover,
    ) -> Result<Self::StageOutputs> {
        let umap_h5 = rover.make_path("umap_h5");
        h5::combine_embeddings(
            &umap_h5,
            chunk_outs.iter().map(|c| &c.umap_h5),
            EmbeddingType::Umap,
        )?;

        let umap_csv: PathBuf = rover.make_path("umap_csv");
        csv::combine_embeddings(&umap_csv, chunk_outs.iter().map(|c| &c.umap_csv))?;

        Ok(Self::StageOutputs { umap_h5, umap_csv })
    }
}
