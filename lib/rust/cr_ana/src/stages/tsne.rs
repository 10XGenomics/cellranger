//! Martian stage RUN_TSNE

use crate::io::{csv, h5};
use crate::types::{EmbeddingResult, EmbeddingType, H5File};
use crate::EXCLUDED_FEATURE_TYPES;
use anyhow::Result;
use bhtsne::BarnesHutTSNE;
use cr_types::reference::feature_reference::FeatureType;
use cr_types::FeatureBarcodeType;
use hdf5_io::matrix::{get_barcodes_between, read_adaptive_csr_matrix};
use martian::prelude::{MartianRover, MartianStage, Resource, StageDef};
use martian_derive::{make_mro, MartianStruct};
use scan_types::matrix::AdaptiveFeatureBarcodeMatrix as FBM;
use serde::{Deserialize, Serialize};
use std::path::PathBuf;

// this initial random state is derived from np.random.RandomState(0).randint(2**31-1)
const RANDOM_SEED: Option<u32> = Some(209652396);
const PERPLEXITY: f64 = 30.0;
const MAX_ITER: usize = 1000;
const STOP_LYING_ITER: usize = 250;
const MOM_SWITCH_ITER: usize = 250;
const THETA: f64 = 0.5;
const N_COMPONENTS: usize = 2;

#[derive(Debug, Deserialize, MartianStruct)]
pub struct TsneStageInputs {
    matrix_h5: H5File,
    pca_h5: H5File,
    random_seed: Option<u32>,
    perplexity: Option<f64>,
    input_pcs: Option<usize>,
    max_dims: Option<usize>,
    max_iter: Option<usize>,
    stop_lying_iter: Option<usize>,
    mom_switch_iter: Option<usize>,
    theta: Option<f64>,
}

#[derive(Debug, Serialize, MartianStruct)]
pub struct TsneStageOutputs {
    tsne_h5: H5File,
    tsne_csv: PathBuf,
}

#[derive(Debug, Serialize, Deserialize, MartianStruct)]
pub struct TsneChunkInputs {
    tsne_dims: usize,
    feature_type: FeatureType,
}

#[derive(Debug, Serialize, Deserialize, MartianStruct)]
pub struct TsneChunkOutputs {
    tsne_h5: H5File,
    tsne_csv: PathBuf,
}

pub struct RunTsne;

#[make_mro(volatile = strict)]
impl MartianStage for RunTsne {
    type StageInputs = TsneStageInputs;
    type StageOutputs = TsneStageOutputs;
    type ChunkInputs = TsneChunkInputs;
    type ChunkOutputs = TsneChunkOutputs;

    fn split(
        &self,
        args: Self::StageInputs,
        _rover: MartianRover,
    ) -> Result<StageDef<Self::ChunkInputs>> {
        let (_nfeatures, ncells, _nnz) = h5::matrix_shape(&args.matrix_h5)?;
        let ncells_mem_gib = (8e-6 * ncells as f64).ceil() as isize;
        let matrix_mem_gib =
            (0.5 + h5::estimate_mem_gib_from_nnz(&args.matrix_h5)?).ceil() as isize;

        let feature_types = h5::matrix_feature_types(&args.matrix_h5)?;
        let components = N_COMPONENTS..=args.max_dims.unwrap_or(N_COMPONENTS);
        let stage_def = components.flat_map(|tsne_dims| {
            feature_types
                .iter()
                .filter(|(&feature_type, &count)| {
                    count >= 2 && !EXCLUDED_FEATURE_TYPES.contains(&feature_type)
                })
                .map(move |(&feature_type, _count)| {
                    let chunk_mem_gib = 4.max(
                        if feature_type == FeatureType::Barcode(FeatureBarcodeType::Multiplexing) {
                            // For multiplexing the log2-transformed count matrix is used
                            matrix_mem_gib
                        } else {
                            ncells_mem_gib
                        },
                    );
                    (
                        Self::ChunkInputs {
                            tsne_dims,
                            feature_type,
                        },
                        Resource::with_mem_gb(chunk_mem_gib),
                    )
                })
        });
        Ok(stage_def.collect())
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
        let mut tsne = BarnesHutTSNE::default();
        tsne.seed = args.random_seed.or(RANDOM_SEED);
        tsne.perplexity = args.perplexity.unwrap_or(PERPLEXITY);
        tsne.max_iter = args.max_iter.unwrap_or(MAX_ITER);
        tsne.theta = args.theta.unwrap_or(THETA);
        tsne.n_dims = chunk_args.tsne_dims as u32;
        tsne.stop_lying_iter = args.stop_lying_iter.or(Some(STOP_LYING_ITER));
        tsne.mom_switch_iter = args.mom_switch_iter.or(Some(MOM_SWITCH_ITER));

        let (mut proj, barcodes) =
            if chunk_args.feature_type == FeatureType::Barcode(FeatureBarcodeType::Multiplexing) {
                // Use feature space for multiplexing due to much lower dimension than gene expression
                // Possible optimization: subsample the full matrix, this ends up much
                // more memory than the PCA case
                let FBM {
                    barcodes, matrix, ..
                } = read_adaptive_csr_matrix(
                    &args.matrix_h5,
                    Some(chunk_args.feature_type.to_string().as_str()),
                    None,
                )?
                .0;

                let proj = matrix
                    .to_dense()
                    .map_mut(|i| (*i as f64 + 1.).log2())
                    .reversed_axes()
                    .as_standard_layout()
                    .to_owned();

                (proj, barcodes)
            } else {
                let matrix = hdf5::File::open(&args.matrix_h5)?.group("matrix")?;
                let barcodes = get_barcodes_between(0, None, &matrix)?;

                (
                    h5::load_transformed_pca_matrix(
                        &args.pca_h5,
                        chunk_args.feature_type,
                        args.input_pcs,
                    )?,
                    barcodes,
                )
            };

        let (num_bcs, _) = proj.dim();
        tsne.perplexity = tsne
            .perplexity
            .min(1.0f64.max(-1.0 + (num_bcs - 1) as f64 / 3.0));

        let embedding = if ((num_bcs - 1) as f64) < 3. * tsne.perplexity {
            // At sufficiently low cell count, we cannot satisfy perplexity >= 1 and this condition below,
            // so the projection is defined as all zeros
            ndarray::Array2::<f64>::zeros((num_bcs, tsne.n_dims as usize))
        } else {
            tsne.init(&mut proj);
            while !tsne.run_n(tsne.max_iter) {}
            tsne.result()
        };

        let result = EmbeddingResult::new(
            &barcodes,
            embedding,
            EmbeddingType::Tsne,
            chunk_args.feature_type,
        );

        let tsne_h5: H5File = rover.make_path("tsne_h5");
        h5::save_embedding(&tsne_h5, &result)?;

        let tsne_csv: PathBuf = rover.make_path("tsne_csv");
        csv::save_embedding(&tsne_csv, &result)?;

        Ok(Self::ChunkOutputs { tsne_h5, tsne_csv })
    }

    fn join(
        &self,
        _args: Self::StageInputs,
        _chunk_defs: Vec<Self::ChunkInputs>,
        chunk_outs: Vec<Self::ChunkOutputs>,
        rover: MartianRover,
    ) -> Result<Self::StageOutputs> {
        let tsne_h5 = rover.make_path("tsne_h5");
        h5::combine_embeddings(
            &tsne_h5,
            chunk_outs.iter().map(|c| &c.tsne_h5),
            EmbeddingType::Tsne,
        )?;

        let tsne_csv: PathBuf = rover.make_path("tsne_csv");
        csv::combine_embeddings(&tsne_csv, chunk_outs.iter().map(|c| &c.tsne_csv))?;
        Ok(Self::StageOutputs { tsne_h5, tsne_csv })
    }
}
