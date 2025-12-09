//! Graph-clustering (Louvain) + merge-clusters stage code
#![expect(missing_docs)]

use crate::io::{csv, h5};
use crate::louvain::{run_louvain, run_louvain_parallel};
use crate::types::{ClusteringResult, ClusteringType, H5File, ReductionType};
use anyhow::{Context, Result};
use cr_types::reference::feature_reference::FeatureType;
use hdf5_io::matrix::read_adaptive_csr_matrix;
use log::info;
use martian::prelude::{MartianRover, MartianStage, Resource, StageDef};
use martian_derive::{MartianStruct, make_mro};
use martian_filetypes::FileTypeRead;
use martian_filetypes::tabular_file::CsvFile;
use ndarray::Array;
use rand::SeedableRng;
use rand::rngs::SmallRng;
use scan_rs::merge_clusters::merge_clusters;
use scan_rs::nn::{find_nn, knn};
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
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
    FeatureType::ProteinExpression,
];

const SKETCH_SAMPLE_SIZE: usize = 500_000;

#[derive(Serialize, Deserialize)]
struct BarcodeAndCluster {
    #[serde(alias = "Barcode")]
    barcode: String,
    #[serde(alias = "Cluster")]
    cluster: usize,
}

#[derive(Debug, Deserialize, MartianStruct)]
pub struct GraphClusteringStageInputs {
    matrix_h5: H5File,
    pca_h5: H5File,
    graphclust_init: Option<CsvFile<BarcodeAndCluster>>,
    num_neighbors: Option<usize>,
    neighbor_a: Option<f64>,
    neighbor_b: Option<f64>,
    input_pcs: Option<usize>,
    resolution: Option<f64>,
    random_seed: Option<usize>,
    threads: usize,
    parallel_clustering: bool,
    sketch: bool,
}

impl GraphClusteringStageInputs {
    fn compute_k(&self) -> Result<usize> {
        let feature_type = h5::matrix_feature_types(&self.matrix_h5)?
            .into_iter()
            .next()
            .unwrap()
            .0;
        let num_bcs = h5::load_transformed_reduction_h5(
            &self.pca_h5,
            feature_type,
            self.input_pcs,
            ReductionType::pca,
        )?
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

#[make_mro(stage_name = RUN_GRAPH_CLUSTERING, volatile = strict)]
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
        let num_bcs = h5::load_transformed_reduction_h5(
            &args.pca_h5,
            *feature_types.iter().next().unwrap().0,
            args.input_pcs,
            ReductionType::pca,
        )?
        .1
        .num_bcs;
        let k: usize = args.compute_k()?;
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
                    Resource::with_mem_gb(mem_gib + 2).threads(args.threads.try_into().unwrap()),
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
        let proj = h5::load_transformed_matrix(
            &args.pca_h5,
            chunk_args.feature_type,
            args.input_pcs,
            ReductionType::pca,
        )?;
        let doing_sketching = args.sketch && SKETCH_SAMPLE_SIZE < proj.nrows();
        if doing_sketching {
            info!(
                "Decided to sketch {SKETCH_SAMPLE_SIZE} rows from proj which has {} rows",
                proj.nrows()
            );
        } else {
            info!(
                "Decided to not sketch because args.sketch={}, SKETCH_SAMPLE_SIZE={SKETCH_SAMPLE_SIZE}, number of rows in proj={}",
                args.sketch,
                proj.nrows()
            );
        }

        let (matrix, _) = read_adaptive_csr_matrix(&args.matrix_h5, None, None)?;
        let k = args.compute_k()?;
        let init_clusters = if let Some(graphclust_init_file) = args.graphclust_init {
            let bc_to_cluster_map: HashMap<_, _> = graphclust_init_file
                .read()?
                .into_iter()
                .map(|x| (x.barcode, x.cluster))
                .collect();
            let init_clusters: Vec<_> = matrix
                .barcodes
                .iter()
                .map(|x| {
                    bc_to_cluster_map
                        .get(x)
                        .copied()
                        .with_context(|| format!("barcode {x} not found in barcode to cluster map"))
                })
                .collect::<Result<_>>()?;
            Some(Array::from_vec(init_clusters))
        } else {
            None
        };

        info!("computing k-nearest neighbors with k = {k}");
        let sampled_indices = if doing_sketching {
            info!("sampling {SKETCH_SAMPLE_SIZE} rows from proj");
            Some(
                rand::seq::index::sample(
                    &mut SmallRng::seed_from_u64(0),
                    proj.nrows(),
                    SKETCH_SAMPLE_SIZE,
                )
                .into_vec(),
            )
        } else {
            None
        };

        let (selected_proj, selected_init_clusters) =
            if let Some(ref sampled_indices_vec) = sampled_indices {
                (
                    &proj.select(ndarray::Axis(0), sampled_indices_vec),
                    init_clusters.map(|x| x.select(ndarray::Axis(0), sampled_indices_vec)),
                )
            } else {
                (&proj, init_clusters)
            };

        let downsampled_indices_to_sampled_indices_map: Option<HashMap<_, _>> = sampled_indices
            .map(|sampled_indices_vec| {
                sampled_indices_vec
                    .into_iter()
                    .enumerate()
                    .map(|(ind, val)| (val, ind))
                    .collect()
            });

        let (kneighbors, ball_tree) = knn::<u32>(&selected_proj.view(), k);

        info!("running louvain");
        let labels = if args.parallel_clustering {
            run_louvain_parallel(
                &kneighbors,
                args.resolution.unwrap_or(RESOLUTION),
                selected_init_clusters.as_ref(),
            )
        } else {
            let seed = args.random_seed.or(Some(RANDOM_SEED));
            run_louvain(
                &kneighbors,
                args.resolution.unwrap_or(RESOLUTION),
                selected_init_clusters.as_ref(),
                seed,
            )
        };

        let labels = if let Some(downsampled_indices_to_sampled_indices_map) =
            downsampled_indices_to_sampled_indices_map
        {
            let all_neighbors = find_nn::<usize>(&proj.view(), 1, &ball_tree, true);
            all_neighbors
                .into_iter()
                .enumerate()
                .map(|(ind, nbr)| {
                    if downsampled_indices_to_sampled_indices_map.contains_key(&ind) {
                        labels[downsampled_indices_to_sampled_indices_map[&ind]]
                    } else {
                        labels[nbr]
                    }
                })
                .collect::<Vec<_>>()
        } else {
            labels
        };

        info!("merging clusters");
        let label_hashset: HashSet<_> = labels.iter().collect();
        info!("cluster labels: {label_hashset:?}");
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
