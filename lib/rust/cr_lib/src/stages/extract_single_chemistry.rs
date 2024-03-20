//! Martian stage EXTRACT_SINGLE_CHEMISTRY.
//! This is a utility stage to extract a single pre-declared library type from
//! the map of chemistries produced by DETECT_CHEMISTRY. If no library type is
//! provided, the stage asserts that only a single chemistry is present and
//! returns it.
//!
//! This is a shim to accomodate code written from when DETECT_CHEMISTRY only
//! produced a single output, and ideally should be eliminated from the pipeline
//! through long-term refactoring.

use anyhow::{anyhow, Result};
use cr_types::chemistry::{ChemistryDef, ChemistryDefs};
use cr_types::LibraryType;
use itertools::Itertools;
use martian::prelude::*;
use martian_derive::{make_mro, MartianStruct};
use serde::{Deserialize, Serialize};

#[derive(Clone, Deserialize, MartianStruct)]
pub struct Inputs {
    pub chemistry_defs: ChemistryDefs,
    pub library_to_extract: Option<LibraryType>,
}

#[derive(Serialize, Deserialize, MartianStruct)]
pub struct Outputs {
    pub chemistry_def: ChemistryDef,
}

pub struct ExtractSingleChemistry;

#[make_mro(volatile = strict)]
impl MartianMain for ExtractSingleChemistry {
    type StageInputs = Inputs;
    type StageOutputs = Outputs;

    fn main(&self, args: Self::StageInputs, _rover: MartianRover) -> Result<Self::StageOutputs> {
        Ok(Self::StageOutputs {
            chemistry_def: if let Some(library_to_extract) = args.library_to_extract {
                args.chemistry_defs
                    .get(&library_to_extract)
                    .ok_or_else(|| anyhow!("no chemistry for {} is present", library_to_extract))?
                    .clone()
            } else {
                args.chemistry_defs.into_values().dedup().exactly_one()?
            },
        })
    }
}
