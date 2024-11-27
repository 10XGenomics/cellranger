//! Replacement for RUN_FBPCA. Could be absorbed into RUN_PCA_NG stage after some refactoring

use crate::io::h5;
use crate::types::H5File;
use anyhow::Result;
use hdf5_io::matrix::read_adaptive_csr_matrix;
use martian::prelude::*;
use martian::MartianVoid;
use martian_derive::{make_mro, martian_filetype, MartianStruct};
use ndarray::Axis;
use scan_rs::dim_red::bk_svd::BkSvd;
use scan_rs::dim_red::Pca;
use scan_types::matrix::AdaptiveFeatureBarcodeMatrix as FBM;
use serde::{Deserialize, Serialize};
use sqz::ScaleAxis;

#[derive(Debug, Deserialize, MartianStruct, Clone)]
pub struct Pca2StageInputs {
    matrix_h5: H5File,
    num_pcs: Option<usize>,
}

martian_filetype!(NpyFile, "npy");

#[derive(Debug, Serialize, Deserialize, MartianStruct)]
pub struct Pca2StageOutputs {
    dimred_matrix: NpyFile,
}

pub struct Pca2Stage;

const NUM_PCS_DEFAULT: usize = 100;

#[make_mro(stage_name = RUN_PCA2, volatile = strict)]
impl MartianStage for Pca2Stage {
    type StageInputs = Pca2StageInputs;
    type StageOutputs = Pca2StageOutputs;
    type ChunkInputs = MartianVoid;
    type ChunkOutputs = MartianVoid;

    fn split(
        &self,
        args: Self::StageInputs,
        _rover: MartianRover,
    ) -> Result<StageDef<Self::ChunkInputs>> {
        let (num_features, num_cells, nnz) = h5::matrix_shape(&args.matrix_h5)?;
        let mem_gib = 10 + 72 * nnz / 1024 / 1024 / 1024;
        println!("num_features={num_features},num_cells={num_cells},nnz={nnz},mem_gib={mem_gib}");
        Ok(StageDef::with_join_resource(
            Resource::with_mem_gb(mem_gib).threads(4),
        ))
    }

    fn main(
        &self,
        _args: Self::StageInputs,
        _chunk_args: Self::ChunkInputs,
        _rover: MartianRover,
    ) -> Result<Self::ChunkOutputs> {
        unreachable!()
    }

    fn join(
        &self,
        args: Self::StageInputs,
        _chunk_defs: Vec<Self::ChunkInputs>,
        _chunk_outs: Vec<Self::ChunkOutputs>,
        rover: MartianRover,
    ) -> Result<Self::StageOutputs> {
        rayon::ThreadPoolBuilder::new()
            .num_threads(rover.get_threads())
            .build_global()?;
        let FBM {
            barcodes: _,
            feature_ids: _,
            matrix,
            ..
        } = read_adaptive_csr_matrix(args.matrix_h5, None, None)?.0;

        let sum_sq_umi_counts = &matrix
            .view()
            .apply(|x: u32| (x as f64).powi(2))
            .sum_axis(ndarray::Axis(0));
        let col_scales = sum_sq_umi_counts.mapv(|c: f64| 1.0 / c.sqrt());
        let scale_cols = ScaleAxis::new(Axis(1), col_scales);
        let norm_matrix = matrix.compose_map(scale_cols);

        let [num_features, num_cells] = norm_matrix.shape();
        let min_dim = std::cmp::min(num_features, num_cells);
        let num_pcs = args.num_pcs.unwrap_or(NUM_PCS_DEFAULT);

        let num_pcs = if min_dim < num_pcs {
            log::warn!(
                "matrix shape {:?} < requested PCs {}, reducing to {}",
                norm_matrix.shape(),
                num_pcs,
                min_dim
            );
            min_dim
        } else {
            num_pcs
        };
        let (u, s, _) = BkSvd {
            k_multiplier: 2.0,
            n_iter: 10,
        }
        .run_pca(&norm_matrix.t().center(Axis(0), None), num_pcs)?;
        // println!("s = {:?}", s);

        let dimred_matrix = u * s;

        let dimred_matrix_file: NpyFile = rover.make_path("dimred_matrix");

        ndarray_npy::write_npy(&dimred_matrix_file, &dimred_matrix)?;

        Ok(Pca2StageOutputs {
            dimred_matrix: dimred_matrix_file,
        })
    }
}
