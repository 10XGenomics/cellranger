#![allow(unused)]
use anyhow::Result;
use martian::AsMartianPrimaryType;
use serde::{Deserialize, Serialize};
use std::path::Path;

#[derive(Debug, Serialize, Deserialize, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[serde(rename_all = "snake_case")]
pub enum OligoPart {
    Bc1,
}

impl AsMartianPrimaryType for OligoPart {
    fn as_martian_primary_type() -> martian::MartianPrimaryType {
        martian::MartianPrimaryType::Str
    }
}
pub fn load_oligos(_: &Path, _: OligoPart) -> Result<Vec<String>> {
    unimplemented!()
}
