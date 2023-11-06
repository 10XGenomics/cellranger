//! Hierarchical clustering stage code

use crate::hclust_utils::get_cluster_representatives;
use crate::io::{csv, h5};
use crate::types::{feature_prefixed, ClusteringResult, ClusteringType, FeatureType, H5File};
use anyhow::Result;
use hclust::{ClusterDirection, DistanceMetric, HierarchicalCluster, LinkageMethod};
use hdf5_io::matrix::read_adaptive_csr_matrix;
use itertools::zip;
use martian::prelude::{MartianRover, MartianStage, Resource, StageDef};
use martian_derive::{make_mro, MartianStruct};
use serde::{Deserialize, Serialize};
use std::cmp::min;
use std::default::Default;
use std::path::PathBuf;

const ACTIVE_FEATURE_TYPES: &[FeatureType] = &[FeatureType::Gene, FeatureType::Antibody];
const HCLUST_STEP_SIZE: usize = 5;

#[derive(Debug, Deserialize, MartianStruct)]
pub struct HierarchicalClusteringStageInputs {
    matrix_h5: H5File,
    graph_clusters_h5: H5File,
    is_antibody_only: bool,
}

#[derive(Debug, Serialize, Deserialize, MartianStruct)]
pub struct HierarchicalClusteringStageOutputs {
    clusters_h5: Option<H5File>,
    clusters_csv: Option<PathBuf>,
}

#[derive(Debug, Serialize, Deserialize, MartianStruct)]
pub struct HierarchicalClusteringChunkInputs {
    feature_type: FeatureType,
}

pub struct HierarchicalClusteringStage;

#[make_mro(stage_name = RUN_HIERARCHICAL_CLUSTERING, mem_gb = 4, threads = 1)]
impl MartianStage for HierarchicalClusteringStage {
    type StageInputs = HierarchicalClusteringStageInputs;
    type StageOutputs = HierarchicalClusteringStageOutputs;
    type ChunkInputs = HierarchicalClusteringChunkInputs;
    type ChunkOutputs = HierarchicalClusteringStageOutputs;

    fn split(
        &self,
        args: Self::StageInputs,
        _rover: MartianRover,
    ) -> Result<StageDef<Self::ChunkInputs>> {
        let mem_gb = h5::est_mem_gb_from_nnz(&args.matrix_h5)?;
        let feature_types = h5::matrix_feature_types(&args.matrix_h5)?;

        let mut stage_def = StageDef::new();
        for &feature_type in ACTIVE_FEATURE_TYPES {
            // do not bother producing when the number of features is <2, or we don't have it
            if feature_types
                .get(&feature_type)
                .map(|&count| count < 2)
                .unwrap_or(true)
            {
                continue;
            }
            let chunk_resource = Resource::new()
                .mem_gb(mem_gb.ceil() as isize + 1)
                .vmem_gb(mem_gb.ceil() as isize + 4)
                .threads(1);
            stage_def.add_chunk_with_resource(Self::ChunkInputs { feature_type }, chunk_resource);
        }
        Ok(stage_def)
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
        let (matrix, _) = read_adaptive_csr_matrix(&args.matrix_h5, retained.as_deref(), None)?;
        let is_feature_prefixed = feature_prefixed(args.is_antibody_only, chunk_args.feature_type);

        let graphclust_results = h5::load_clustering(
            &args.graph_clusters_h5,
            ClusteringType::Louvain,
            chunk_args.feature_type,
            args.is_antibody_only,
        )?;
        let mat = matrix.matrix.to_csmat().to_csc();
        let cluster_reps = get_cluster_representatives(&graphclust_results.labels, &mat);

        let num_graphclust_clusters = cluster_reps.ncols();
        if num_graphclust_clusters.saturating_sub(HCLUST_STEP_SIZE) > 1 {
            let clusters = HierarchicalCluster::new(
                &cluster_reps,
                DistanceMetric::Euclidean,
                LinkageMethod::Ward,
                ClusterDirection::Columns,
            );

            let clusters_csv: PathBuf = rover.make_path("clusters_csv");
            let step_size = min(HCLUST_STEP_SIZE, num_graphclust_clusters);
            let num_clusters_to_output = (2..=num_graphclust_clusters
                .saturating_sub(HCLUST_STEP_SIZE))
                .rev()
                .step_by(step_size);
            let temp_h5_files: Vec<H5File> = num_clusters_to_output
                .clone()
                .map(|n| rover.make_path(format!("clusters_h5_{n}_clusters")))
                .collect();
            for (num_clusters, temp_h5_file) in zip(num_clusters_to_output, &temp_h5_files) {
                // Note that old_cluster_i is mapped to cluster_map[old_cluster_i - 1]
                // because clusters are labelled 1..max_cluster without zero indexing.
                let cluster_map = clusters.fcluster(&num_clusters);
                let result = ClusteringResult::new(
                    ClusteringType::Hierarchical(num_clusters),
                    chunk_args.feature_type,
                    graphclust_results
                        .labels
                        .iter()
                        .map(|x| cluster_map[(x - 1) as usize] as i64)
                        .collect(),
                    is_feature_prefixed,
                );
                csv::save_clustering(&clusters_csv, &result, &matrix.barcodes)?;
                h5::save_clustering(temp_h5_file, &result)?;
            }

            let clusters_h5 = rover.make_path("clusters_h5");
            h5::combine_clusterings(&clusters_h5, temp_h5_files.iter())?;
            Ok(Self::ChunkOutputs {
                clusters_csv: Some(clusters_csv),
                clusters_h5: Some(clusters_h5),
            })
        } else {
            Ok(Self::ChunkOutputs {
                clusters_csv: None,
                clusters_h5: None,
            })
        }
    }

    fn join(
        &self,
        _args: Self::StageInputs,
        _chunk_defs: Vec<Self::ChunkInputs>,
        chunk_outs: Vec<Self::ChunkOutputs>,
        rover: MartianRover,
    ) -> Result<Self::StageOutputs> {
        let (clusters_h5_raw, clusters_csv_raw): (Vec<_>, Vec<_>) = chunk_outs
            .into_iter()
            .map(|c| (c.clusters_h5, c.clusters_csv))
            .unzip();
        let clusters_h5_vec: Vec<H5File> = clusters_h5_raw.into_iter().flatten().collect();
        let clusters_csv_vec: Vec<PathBuf> = clusters_csv_raw.into_iter().flatten().collect();

        let clusters_h5 = if clusters_h5_vec.is_empty() {
            None
        } else {
            let clusters_h5 = rover.make_path("clusters_h5");
            h5::combine_clusterings(&clusters_h5, clusters_h5_vec.iter())?;
            Some(clusters_h5)
        };

        let clusters_csv = if clusters_csv_vec.is_empty() {
            None
        } else {
            let clusters_csv: PathBuf = rover.make_path("clusters_csv");
            csv::combine_clusterings(&clusters_csv, clusters_csv_vec.iter())?;
            Some(clusters_csv)
        };

        Ok(Self::StageOutputs {
            clusters_csv,
            clusters_h5,
        })
    }
}
