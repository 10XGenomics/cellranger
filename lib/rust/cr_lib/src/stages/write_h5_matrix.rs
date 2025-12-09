//! Martian stage WRITE_H5_MATRIX.
#![deny(missing_docs)]

use crate::types::FeatureReferenceFormat;
use crate::utils::estimate_mem::barcode_mem_gib;
use anyhow::Result;
use cr_h5::count_matrix::write_matrix_h5;
use cr_types::chemistry::{ChemistryDefs, ChemistryDefsExt};
use cr_types::{BarcodeIndexOutput, CountShardFile, GemWell, H5File};
use martian::{MartianRover, MartianStage, MartianVoid, Resource, StageDef};
use martian_derive::{MartianStruct, make_mro};
use martian_filetypes::FileTypeRead;
use serde::{Deserialize, Serialize};

/// Martian stage WRITE_H5_MATRIX.
/// Output the feature-barcode matrix HDF5 file `raw_feature_bc_matrix.h5`.
pub struct WriteH5Matrix;

#[derive(Deserialize, Clone, MartianStruct)]
pub struct WriteH5MatrixStageInputs {
    pub gem_well: GemWell,
    pub counts: Vec<CountShardFile>,
    pub feature_reference: FeatureReferenceFormat,
    pub chemistry_defs: ChemistryDefs,
    pub sample_id: String,
    pub barcode_index_output: BarcodeIndexOutput,
}

#[derive(Serialize, Deserialize, Clone, MartianStruct)]
pub struct WriteH5MatrixStageOutputs {
    pub matrix: H5File,
}
#[make_mro(volatile = strict)]
impl MartianStage for WriteH5Matrix {
    type StageInputs = WriteH5MatrixStageInputs;
    type StageOutputs = WriteH5MatrixStageOutputs;
    type ChunkInputs = MartianVoid;
    type ChunkOutputs = MartianVoid;

    fn split(
        &self,
        args: Self::StageInputs,
        _rover: MartianRover,
    ) -> Result<StageDef<Self::ChunkInputs>> {
        // Multi uses more memory for sample_barcodes and other data.
        let barcodes_count = args.barcode_index_output.num_barcodes;
        // bytes_per_barcode and offset_gib are empirically determined.
        let mem_gib = barcode_mem_gib(barcodes_count, 220, 2);
        println!("barcode_count={barcodes_count},mem_gib={mem_gib}");
        Ok(StageDef::with_join_resource(Resource::with_mem_gb(mem_gib)))
    }

    fn main(
        &self,
        _args: Self::StageInputs,
        _chunk_args: MartianVoid,
        _rover: MartianRover,
    ) -> Result<Self::ChunkOutputs> {
        unreachable!()
    }

    fn join(
        &self,
        args: Self::StageInputs,
        _chunk_defs: Vec<MartianVoid>,
        _chunk_outs: Vec<MartianVoid>,
        rover: MartianRover,
    ) -> Result<Self::StageOutputs> {
        let raw_feature_bc_matrix = rover.make_path("raw_feature_bc_matrix");
        write_matrix_h5(
            &raw_feature_bc_matrix,
            &args.counts,
            &args.feature_reference.read()?,
            &args.sample_id,
            &args.chemistry_defs.description(),
            args.gem_well,
            &args.barcode_index_output.index.read()?,
            &rover.pipelines_version(),
        )?;
        Ok(Self::StageOutputs {
            matrix: raw_feature_bc_matrix,
        })
    }
}
