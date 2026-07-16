//! Martian stage COMBINE_CHEMISTRY_DEFS
#![deny(missing_docs)]

use anyhow::{Result, ensure};
use cr_types::MetricsFile;
use cr_types::chemistry::{ChemistryDef, ChemistryDefs, ChemistryDefsExt};
use cr_types::types::{LibraryType, VdjChainType};
use itertools::{Itertools, zip_eq};
use martian::{MartianMain, MartianRover};
use martian_derive::{MartianStruct, make_mro};
use metric::JsonReporter;
use serde::{Deserialize, Serialize};
use std::iter::chain;
use vdj_reference::VdjReceptor;

/// Select the library type based on the receptor. The receptor is None for VDJ de-novo.
fn library_type_from_vdj_receptor(receptor: Option<VdjReceptor>) -> LibraryType {
    match receptor {
        None => LibraryType::Vdj(VdjChainType::Auto),
        Some(VdjReceptor::IG) => LibraryType::Vdj(VdjChainType::VdjB),
        Some(VdjReceptor::TR) => LibraryType::Vdj(VdjChainType::VdjT),
        Some(VdjReceptor::TRGD) => LibraryType::Vdj(VdjChainType::VdjTGD),
    }
}

#[derive(Clone, Deserialize, MartianStruct)]
pub struct CombineChemistryDefsStageInputs {
    pub count_chemistry_defs: Option<ChemistryDefs>,
    pub vdj_chemistry_defs: Option<Vec<ChemistryDef>>,
    pub vdj_receptors: Option<Vec<Option<VdjReceptor>>>,
}

#[derive(Deserialize, Serialize, MartianStruct)]
pub struct CombineChemistryDefsStageOutputs {
    pub chemistry_defs: ChemistryDefs,
    pub summary: MetricsFile,
}

/// Martian stage COMBINE_CHEMISTRY_DEFS
pub struct CombineChemistryDefs;

#[make_mro(volatile = strict)]
impl MartianMain for CombineChemistryDefs {
    type StageInputs = CombineChemistryDefsStageInputs;
    type StageOutputs = CombineChemistryDefsStageOutputs;

    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        assert_eq!(
            args.vdj_receptors.is_some(),
            args.vdj_chemistry_defs.is_some()
        );

        let vdj_receptors = args.vdj_receptors.unwrap_or_default();
        ensure!(vdj_receptors.iter().all_unique(), "{vdj_receptors:?}");

        let vdj_chemistry_defs: ChemistryDefs =
            zip_eq(vdj_receptors, args.vdj_chemistry_defs.unwrap_or_default())
                .map(|(receptor, chem_def)| (library_type_from_vdj_receptor(receptor), chem_def))
                .collect();

        let chemistry_defs: ChemistryDefs = chain(
            args.count_chemistry_defs.unwrap_or_default(),
            vdj_chemistry_defs,
        )
        .collect();

        let summary = MetricsFile::from_reporter(
            &rover,
            "chemistry_defs",
            &JsonReporter::from(("chemistry_defs", chemistry_defs.to_value())),
        )?;
        Ok(CombineChemistryDefsStageOutputs {
            chemistry_defs,
            summary,
        })
    }
}
