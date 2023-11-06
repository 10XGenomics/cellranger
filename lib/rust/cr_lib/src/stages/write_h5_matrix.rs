//! Martian stage WRITE_H5_MATRIX.

use crate::count_matrix::write_matrix_h5;
use crate::types::FeatureReferenceFormat;
use crate::{CountShardFile, H5File};
use anyhow::Result;
use barcode::Barcode;
use cr_types::chemistry::ChemistryDef;
use cr_types::types::{BarcodeIndexFormat, GemWell};
use martian::{MartianMain, MartianRover};
use martian_derive::{make_mro, MartianStruct};
use martian_filetypes::FileTypeRead;
use serde::{Deserialize, Serialize};
use std::hash::Hash;

/// Output the feature-barcode matrix HDF5 file `raw_feature_bc_matrix.h5`.
pub struct WriteH5Matrix;

#[derive(Deserialize, Clone, MartianStruct)]
pub struct StageInputs<B: Eq + Hash> {
    pub gem_well: GemWell,
    pub counts: Vec<CountShardFile>,
    pub feature_reference: FeatureReferenceFormat,
    pub chemistry_def: ChemistryDef,
    pub sample_id: String,
    pub barcode_index: BarcodeIndexFormat<B>,
}

#[derive(Serialize, Deserialize, Clone, MartianStruct)]
pub struct StageOutputs {
    pub matrix: H5File,
}

#[make_mro(mem_gb = 7, volatile = strict)]
impl MartianMain for WriteH5Matrix {
    type StageInputs = StageInputs<Barcode>;
    type StageOutputs = StageOutputs;

    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        let raw_feature_bc_matrix: H5File = rover.make_path("raw_feature_bc_matrix");
        let barcode_index = args.barcode_index.read()?;
        write_matrix_h5(
            &raw_feature_bc_matrix,
            &args.counts,
            &args.feature_reference.read()?,
            &args.sample_id,
            &args.chemistry_def,
            args.gem_well,
            &barcode_index,
            &rover.pipelines_version(),
        )?;
        Ok(StageOutputs {
            matrix: raw_feature_bc_matrix,
        })
    }
}
