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
pub fn spot_pitch(_: &Path) -> Result<u32> {
    unimplemented!()
}
pub fn validate_slide_id_name(slide_id: &str) -> Result<String> {
    unimplemented!()
}
