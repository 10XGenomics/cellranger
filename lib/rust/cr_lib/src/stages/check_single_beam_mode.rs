//! Martian stage CHECK_SINGLE_BEAM_MODE

use anyhow::{bail, Result};
use cr_types::reference::feature_reference::BeamMode;
use itertools::Itertools;
use martian::prelude::*;
use martian_derive::{make_mro, MartianStruct};
use serde::{Deserialize, Serialize};

#[derive(Deserialize, MartianStruct)]
pub struct CheckSingleBeamModeStageInputs {
    pub beam_modes: Vec<Option<BeamMode>>,
}

#[derive(Serialize, Deserialize, MartianStruct)]
#[cfg_attr(test, derive(Debug))]
pub struct CheckSingleBeamModeStageOutputs {
    pub beam_mode: Option<BeamMode>,
}

pub struct CheckSingleBeamMode;

#[make_mro(volatile = strict)]
impl MartianMain for CheckSingleBeamMode {
    type StageInputs = CheckSingleBeamModeStageInputs;
    type StageOutputs = CheckSingleBeamModeStageOutputs;

    fn main(&self, args: Self::StageInputs, _rover: MartianRover) -> Result<Self::StageOutputs> {
        let beam_modes: Vec<_> = args.beam_modes.into_iter().unique().collect();
        match beam_modes.as_slice() {
            [] => bail!("no BEAM modes detected"),
            [one] => Ok(CheckSingleBeamModeStageOutputs { beam_mode: *one }),
            more => bail!("multiple BEAM modes detected: {more:?}"),
        }
    }
}
