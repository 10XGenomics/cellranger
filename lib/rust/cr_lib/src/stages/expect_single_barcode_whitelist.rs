//! Martian stage EXPECT_SINGLE_BARCODE_WHITELIST.

use anyhow::{ensure, Result};
use martian::prelude::*;
use martian_derive::{make_mro, MartianStruct};
use serde::{Deserialize, Serialize};

#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct Inputs {
    pub barcode_whitelists: Vec<String>,
}

#[derive(Clone, Serialize, Deserialize, MartianStruct, PartialEq)]
pub struct Outputs {
    pub barcode_whitelist: Option<String>,
}

pub struct ExpectSingleBarcodeWhitelist;

#[make_mro]
impl MartianMain for ExpectSingleBarcodeWhitelist {
    type StageInputs = Inputs;
    type StageOutputs = Outputs;
    fn main(
        &self,
        mut args: Self::StageInputs,
        _rover: MartianRover,
    ) -> Result<Self::StageOutputs> {
        ensure!(
            args.barcode_whitelists.len() < 2,
            "Expected at most one barcode whitelist but found {}:\n{}",
            args.barcode_whitelists.len(),
            args.barcode_whitelists.join("\n")
        );
        Ok(Outputs {
            barcode_whitelist: args.barcode_whitelists.pop(),
        })
    }
}
