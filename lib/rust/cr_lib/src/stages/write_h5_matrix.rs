//! Martian stage WRITE_H5_MATRIX.

use crate::types::FeatureReferenceFormat;
use anyhow::Result;
use cr_h5::count_matrix::write_matrix_h5;
use cr_types::chemistry::{ChemistryDefs, ChemistryDefsExt};
use cr_types::{BarcodeIndexFormat, CountShardFile, GemWell, H5File};
use martian::{MartianMain, MartianRover};
use martian_derive::{make_mro, MartianStruct};
use martian_filetypes::FileTypeRead;
use serde::{Deserialize, Serialize};

/// Output the feature-barcode matrix HDF5 file `raw_feature_bc_matrix.h5`.
pub struct WriteH5Matrix;

#[derive(Deserialize, Clone, MartianStruct)]
pub struct StageInputs {
    pub gem_well: GemWell,
    pub counts: Vec<CountShardFile>,
    pub feature_reference: FeatureReferenceFormat,
    pub chemistry_defs: ChemistryDefs,
    pub sample_id: String,
    pub barcode_index: BarcodeIndexFormat,
}

#[derive(Serialize, Deserialize, Clone, MartianStruct)]
pub struct StageOutputs {
    pub matrix: H5File,
}

#[make_mro(mem_gb = 7, volatile = strict)]
impl MartianMain for WriteH5Matrix {
    type StageInputs = StageInputs;
    type StageOutputs = StageOutputs;

    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        let raw_feature_bc_matrix = rover.make_path("raw_feature_bc_matrix");
        write_matrix_h5(
            &raw_feature_bc_matrix,
            &args.counts,
            &args.feature_reference.read()?,
            &args.sample_id,
            &args.chemistry_defs.description(),
            args.gem_well,
            &args.barcode_index.read()?,
            &rover.pipelines_version(),
        )?;
        Ok(StageOutputs {
            matrix: raw_feature_bc_matrix,
        })
    }
}
