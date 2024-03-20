//! Graph-clustering (Louvain) + merge-clusters stage code

use crate::io::{csv, h5};
use crate::louvain::{run_louvain, run_louvain_parallel};
use crate::types::{ClusteringResult, ClusteringType, H5File};
use anyhow::Result;
use cr_types::reference::feature_reference::FeatureType;
use hdf5_io::matrix::read_adaptive_csr_matrix;
use log::info;
use martian::prelude::{MartianRover, MartianStage, Resource, StageDef};
use martian_derive::{make_mro, MartianStruct, MartianType};
use scan_rs::merge_clusters::merge_clusters;
use scan_rs::nn::knn;
use serde::{Deserialize, Serialize};
use std::default::Default;
use std::path::PathBuf;

const DEFAULT_K: usize = 1;
const NEIGHBOR_A: f64 = -230.0;
const NEIGHBOR_B: f64 = 120.0;
const RESOLUTION: f64 = 1.0;
const RANDOM_SEED: usize = 0xBADC0FFEE0DDF00D;
const ACTIVE_FEATURE_TYPES: &[FeatureType] = &[
    FeatureType::Gene,
    FeatureType::Barcode(cr_types::FeatureBarcodeType::Antibody),
];

#[derive(
    Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize, MartianType, Default,
)]
pub enum SimilarityType {
    #[serde(rename = "nn")]
    #[default]
    NN,
    #[serde(rename = "snn")]
    SNN,
}

#[derive(Debug, Deserialize, MartianStruct)]
pub struct GraphClusteringStageInputs {
    matrix_h5: H5File,
    pca_h5: H5File,
    num_neighbors: Option<usize>,
    neighbor_a: Option<f64>,
    neighbor_b: Option<f64>,
    input_pcs: Option<usize>,
    resolution: Option<f64>,
    random_seed: Option<usize>,
    threads: usize,
    parallel_clustering: bool,
}

impl GraphClusteringStageInputs {
    fn compute_k(&self) -> Result<usize> {
        let feature_type = h5::matrix_feature_types(&self.matrix_h5)?
            .into_iter()
            .next()
            .unwrap()
            .0;
        let num_bcs = h5::load_transformed_pca(&self.pca_h5, feature_type, self.input_pcs)?
            .1
            .num_bcs;
        let neighbor_a = self.neighbor_a.unwrap_or(NEIGHBOR_A);
        let neighbor_b = self.neighbor_b.unwrap_or(NEIGHBOR_B);
        let k = (neighbor_a + neighbor_b * (num_bcs as f64).log10()).round() as usize;
        Ok(self
            .num_neighbors
            .unwrap_or(k)
            .clamp(DEFAULT_K, DEFAULT_K.max(num_bcs.saturating_sub(1))))
    }
}

#[derive(Debug, Serialize, Deserialize, MartianStruct)]
pub struct GraphClusteringStageOutputs {
    clusters_h5: H5File,
    clusters_csv: PathBuf,
}

#[derive(Debug, Serialize, Deserialize, MartianStruct)]
pub struct GraphClusteringChunkInputs {
    feature_type: FeatureType,
}

pub struct GraphClusteringStage;

#[make_mro(stage_name = RUN_GRAPH_CLUSTERING_NG, volatile = strict)]
impl MartianStage for GraphClusteringStage {
    type StageInputs = GraphClusteringStageInputs;
    type StageOutputs = GraphClusteringStageOutputs;
    type ChunkInputs = GraphClusteringChunkInputs;
    type ChunkOutputs = GraphClusteringStageOutputs;

    fn split(
        &self,
        args: Self::StageInputs,
        _rover: MartianRover,
    ) -> Result<StageDef<Self::ChunkInputs>> {
        let feature_types = h5::matrix_feature_types(&args.matrix_h5)?;
        let num_bcs = h5::load_transformed_pca(
            &args.pca_h5,
            *feature_types.iter().next().unwrap().0,
            args.input_pcs,
        )?
        .1
        .num_bcs;
        let k = args.compute_k()?;
        // (adjacency matrix + observed edge list + weighted edge list) * max no. edges
        // (u32 = 4 + u32 = 4 + usize = 8 + f64 = 8) = 24 * num_edges
        let mem_gib = (2.5
            + h5::estimate_mem_gib_from_nnz(&args.matrix_h5)?
            + ((4 + 4 + 8 + 8) * (num_bcs * k)) as f64 / 1e9)
            .ceil() as isize;

        Ok(ACTIVE_FEATURE_TYPES
            .iter()
            .filter(|feature_type| {
                feature_types
                    .get(feature_type)
                    .is_some_and(|&count| count >= 2)
            })
            .map(|&feature_type| {
                (
                    Self::ChunkInputs { feature_type },
                    Resource::with_mem_gb(mem_gib).threads(args.threads.try_into().unwrap()),
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
            .num_threads(args.threads)
            .build_global()?;
        let proj =
            h5::load_transformed_pca_matrix(&args.pca_h5, chunk_args.feature_type, args.input_pcs)?;
        let (matrix, _) = read_adaptive_csr_matrix(&args.matrix_h5, None, None)?;
        let k = args.compute_k()?;

        info!("computing k-nearest neighbors with k = {}", k);
        let neighbors = knn::<u32>(&proj.view(), k);

        info!("running louvain");
        let labels = if args.parallel_clustering {
            run_louvain_parallel(&neighbors, args.resolution.unwrap_or(RESOLUTION))
        } else {
            let seed = args.random_seed.or(Some(RANDOM_SEED));
            run_louvain(&neighbors, args.resolution.unwrap_or(RESOLUTION), seed)
        };

        info!("merging clusters");
        let labels = merge_clusters(&matrix, &proj, labels)
            .into_iter()
            .map(|v| v as i64 + 1)
            .collect::<Vec<_>>();
        let result =
            ClusteringResult::new(ClusteringType::Louvain, chunk_args.feature_type, labels);
        let clusters_h5 = rover.make_path("clusters_h5");
        h5::save_clustering(&clusters_h5, &result)?;

        let clusters_csv: PathBuf = rover.make_path("clusters_csv");
        csv::save_clustering(&clusters_csv, &result, &matrix.barcodes)?;

        Ok(Self::ChunkOutputs {
            clusters_h5,
            clusters_csv,
        })
    }

    fn join(
        &self,
        _args: Self::StageInputs,
        _chunk_defs: Vec<Self::ChunkInputs>,
        chunk_outs: Vec<Self::ChunkOutputs>,
        rover: MartianRover,
    ) -> Result<Self::StageOutputs> {
        let clusters_h5 = rover.make_path("clusters_h5");
        h5::combine_clusterings(&clusters_h5, chunk_outs.iter().map(|c| &c.clusters_h5))?;

        let clusters_csv: PathBuf = rover.make_path("clusters_csv");
        csv::combine_clusterings(&clusters_csv, chunk_outs.iter().map(|c| &c.clusters_csv))?;

        Ok(Self::StageOutputs {
            clusters_h5,
            clusters_csv,
        })
    }
}
