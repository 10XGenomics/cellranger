//! Martian stage LOGIC_NOT
//! Return the logical not of its input argument.
#![deny(missing_docs)]

use anyhow::Result;
use martian::{MartianMain, MartianRover};
use martian_derive::{MartianStruct, make_mro};
use serde::{Deserialize, Serialize};

/// The Martian stage inputs.
#[derive(Clone, Deserialize, MartianStruct)]
pub struct LogicNotStageInputs {
    pub boolean: Option<bool>,
}

/// The Martian stage outputs.
#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct LogicNotStageOutputs {
    pub not_boolean: Option<bool>,
}

/// Martian stage LOGIC_NOT.
/// Return the logical not of its input argument.
pub struct LogicNot;

#[make_mro(volatile = strict)]
impl MartianMain for LogicNot {
    type StageInputs = LogicNotStageInputs;
    type StageOutputs = LogicNotStageOutputs;

    /// Run the Martian stage LOGIC_NOT.
    fn main(&self, args: Self::StageInputs, _rover: MartianRover) -> Result<Self::StageOutputs> {
        Ok(Self::StageOutputs {
            not_boolean: args.boolean.map(|x| !x),
        })
    }
}
