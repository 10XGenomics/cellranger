//! Martian stage GET_CHEMISTRY_DEF
//! Load chemistry def from a chemistry name or a custom chemistry def

use anyhow::Result;
use cr_types::chemistry::{ChemistryDef, ChemistryName};
use martian::{MartianMain, MartianRover};
use martian_derive::{make_mro, MartianStruct};
use serde::{Deserialize, Serialize};

/// The Martian stage inputs.
#[derive(Clone, Deserialize, MartianStruct)]
pub struct StageInputs {
    chemistry_name: ChemistryName,
    custom_chemistry_def: Option<ChemistryDef>,
}

/// The Martian stage outputs.
#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct StageOutputs {
    chemistry_def: ChemistryDef,
}

pub struct GetChemistryDef;

#[make_mro(volatile = strict)]
impl MartianMain for GetChemistryDef {
    type StageInputs = StageInputs;
    type StageOutputs = StageOutputs;

    /// Run the Martian stage LOGIC_NOT.
    fn main(&self, args: Self::StageInputs, _rover: MartianRover) -> Result<Self::StageOutputs> {
        Ok(Self::StageOutputs {
            chemistry_def: match args.chemistry_name {
                ChemistryName::Custom => args
                    .custom_chemistry_def
                    .expect("expecting a custom chemistry def to be present"),
                name => ChemistryDef::named(name),
            },
        })
    }
}
