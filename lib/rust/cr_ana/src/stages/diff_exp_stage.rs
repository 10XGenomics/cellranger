//! Differential Expression stage code

use crate::io::{csv, h5};
use crate::types::{ClusteringType, H5File};
use anyhow::Result;
use cr_types::reference::feature_reference::FeatureType;
use cr_types::FeatureBarcodeType;
use diff_exp::{compute_sseq_params, sseq_differential_expression};
use hdf5_io::matrix::{read_adaptive_csr_matrix, read_matrix_metadata};
use itertools::Itertools;
use martian::prelude::{MartianRover, MartianStage, Resource, StageDef};
use martian_derive::{make_mro, MartianStruct, MartianType};
use martian_filetypes::bin_file::BincodeFile;
use martian_filetypes::lz4_file::Lz4;
use martian_filetypes::{LazyFileTypeIO, LazyWrite};
use ndarray::{s, Array1, Array2};
use scan_types::matrix::{AdaptiveFeatureBarcodeMatrix as FBM, MatrixMetadata};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::iter::zip;
use std::path::PathBuf;

const CHUNK_SIZE: usize = 8;

#[derive(Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize, MartianType)]
pub(crate) struct ClusterKey {
    pub(crate) feature_type: FeatureType,
    pub(crate) cluster_type: ClusteringType,
    pub(crate) cluster_id: i64,
}

#[derive(Serialize, Deserialize)]
pub struct DiffExpResult {
    norm_mean_a: Vec<f64>,
    log2_fold_change: Vec<f64>,
    adjusted_p_value: Vec<f64>,
}

impl DiffExpResult {
    fn new(diff_exp: diff_exp::diff_exp::DiffExpResult) -> Self {
        let diff_exp::diff_exp::DiffExpResult {
            normalized_mean_in,
            log2_fold_change,
            adjusted_p_values,
            ..
        } = diff_exp;
        DiffExpResult {
            norm_mean_a: normalized_mean_in.into_raw_vec(),
            log2_fold_change: log2_fold_change.into_raw_vec(),
            adjusted_p_value: adjusted_p_values.into_raw_vec(),
        }
    }
}

#[derive(Deserialize, MartianStruct)]
pub struct DiffExpStageInputs {
    matrix_h5: H5File,
    clustering_h5: H5File,
    is_antibody_only: bool,
}

#[derive(Serialize, MartianStruct)]
pub struct DiffExpStageOutputs {
    diffexp_h5: H5File,
    diffexp_csv: PathBuf,
}

#[derive(Serialize, Deserialize, MartianStruct)]
pub struct DiffExpChunkInputs {
    cluster_keys: Vec<ClusterKey>,
}

#[derive(Serialize, Deserialize, MartianStruct)]
pub struct DiffExpChunkOutputs {
    diffexp: Lz4<BincodeFile<Vec<DiffExpResult>>>,
}

pub struct DiffExpStage;

#[make_mro(stage_name = RUN_DIFFERENTIAL_EXPRESSION_NG, volatile = strict)]
impl MartianStage for DiffExpStage {
    type StageInputs = DiffExpStageInputs;
    type StageOutputs = DiffExpStageOutputs;
    type ChunkInputs = DiffExpChunkInputs;
    type ChunkOutputs = DiffExpChunkOutputs;

    fn split(
        &self,
        args: Self::StageInputs,
        _rover: MartianRover,
    ) -> Result<StageDef<Self::ChunkInputs>> {
        let mem_gib = 6 + h5::estimate_mem_gib_from_nnz(&args.matrix_h5)?.ceil() as isize;
        let stage_def = h5::load_clusterings(&args.clustering_h5, args.is_antibody_only)?
            .into_iter()
            .flat_map(|c| {
                (1..=c.num_clusters()).map(move |cluster_id| ClusterKey {
                    feature_type: c.feature_type,
                    cluster_type: c.clustering_type,
                    cluster_id,
                })
            })
            .sorted()
            .chunks(CHUNK_SIZE)
            .into_iter()
            .map(|chunk| {
                (
                    DiffExpChunkInputs {
                        cluster_keys: chunk.collect(),
                    },
                    Resource::with_mem_gb(mem_gib),
                )
            })
            .collect();
        Ok(stage_def)
    }

    fn main(
        &self,
        args: Self::StageInputs,
        chunk_args: Self::ChunkInputs,
        rover: MartianRover,
    ) -> Result<Self::ChunkOutputs> {
        rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .build_global()?;
        let clusterings: HashMap<_, _> =
            h5::load_clusterings(&args.clustering_h5, args.is_antibody_only)?
                .into_iter()
                .map(|c| ((c.feature_type, c.clustering_type), c))
                .collect();
        let retained = if args.is_antibody_only {
            Some(FeatureType::Barcode(FeatureBarcodeType::Antibody).to_string())
        } else {
            Some(FeatureType::Gene.to_string())
        };
        let FBM { matrix, .. } =
            read_adaptive_csr_matrix(&args.matrix_h5, retained.as_deref(), None)?.0;
        let sseq_params = compute_sseq_params(&matrix, None, None, None);
        let outputs = Self::ChunkOutputs {
            diffexp: rover.make_path("diffexp"),
        };
        let mut writer = outputs.diffexp.lazy_writer()?;
        for cluster_key in chunk_args.cluster_keys {
            let key = (cluster_key.feature_type, cluster_key.cluster_type);
            let clustering = &clusterings[&key];
            let mut cond_a = vec![];
            let mut cond_b = vec![];
            for (i, &l) in clustering.labels.iter().enumerate() {
                if l == cluster_key.cluster_id {
                    cond_a.push(i);
                } else {
                    cond_b.push(i);
                }
            }
            let result = DiffExpResult::new(sseq_differential_expression(
                &matrix,
                &cond_a,
                &cond_b,
                &sseq_params,
                None,
            ));
            writer.write_item(&result)?;
        }
        writer.finish()?;
        Ok(outputs)
    }

    fn join(
        &self,
        args: Self::StageInputs,
        chunk_defs: Vec<Self::ChunkInputs>,
        chunk_outs: Vec<Self::ChunkOutputs>,
        rover: MartianRover,
    ) -> Result<Self::StageOutputs> {
        rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .build_global()?;
        let retained = if args.is_antibody_only {
            Some(FeatureType::Barcode(FeatureBarcodeType::Antibody).to_string())
        } else {
            Some(FeatureType::Gene.to_string())
        };
        let MatrixMetadata {
            feature_ids,
            feature_names,
            feature_types,
            ..
        } = read_matrix_metadata(&args.matrix_h5, retained.as_deref())?;
        let clusterings: HashMap<_, _> =
            h5::load_clusterings(&args.clustering_h5, args.is_antibody_only)?
                .into_iter()
                .map(|c| ((c.feature_type, c.clustering_type.lc()), c))
                .collect();
        let n_genes = feature_ids.len();
        let diffexp_h5 = rover.make_path("diffexp_h5");
        let h5_file = hdf5::File::create(&diffexp_h5)?;
        let h5_group = h5_file.create_group("all_differential_expression")?;
        let diffexp_csv: PathBuf = rover.make_path("diffexp_csv");
        for ((feature_type, typ), group) in &zip(chunk_defs, chunk_outs)
            .flat_map(|(def, out)| {
                let reader = out
                    .diffexp
                    .lazy_reader()
                    .expect("failed to open DE result file");
                zip(def.cluster_keys, reader)
            })
            .group_by(|(key, _)| (key.feature_type, key.cluster_type))
        {
            let clusterings_key = (feature_type, typ.lc());
            let num_clusters = clusterings[&clusterings_key].num_clusters() as usize;
            let mut data = Array2::<f64>::zeros((n_genes, 3 * num_clusters));
            for (key, res) in group {
                // convert from 1-indexed to 0-indexed
                let i = key.cluster_id as usize - 1;
                let DiffExpResult {
                    norm_mean_a,
                    log2_fold_change,
                    adjusted_p_value,
                    ..
                } = res?;
                data.slice_mut(s![.., 3 * i])
                    .assign(&Array1::from(norm_mean_a));
                data.slice_mut(s![.., 3 * i + 1])
                    .assign(&Array1::from(log2_fold_change));
                data.slice_mut(s![.., 3 * i + 2])
                    .assign(&Array1::from(adjusted_p_value));
            }
            let group = format!("_{}_{}", feature_type.as_snake_case(), typ.lc());
            h5_group
                .create_group(&group)?
                .new_dataset::<f64>()
                .shape(data.dim())
                .create("data")?
                .write(&data)?;

            csv::save_differential_expression(
                &diffexp_csv,
                typ,
                feature_type,
                &feature_ids,
                &feature_names,
                &data,
            )?;
        }
        h5_file
            .new_dataset::<i64>()
            .shape((feature_types.indices.len(),))
            .create("diffexp_feature_indices")?
            .write(&feature_types.indices)?;

        Ok(Self::StageOutputs {
            diffexp_h5,
            diffexp_csv,
        })
    }
}
