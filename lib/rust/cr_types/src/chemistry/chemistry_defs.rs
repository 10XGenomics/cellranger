//! Standard chemistry definitions imported from a json file
#![deny(missing_docs)]

use crate::chemistry::{ChemistryDef, ChemistryName};
use anyhow::{Context, Result};
use itertools::Itertools;
use lazy_static::lazy_static;
use std::collections::HashMap;

const CHEMISTRY_DEFS_JSON_STR: &str = std::include_str!("chemistry_defs.json");

fn read_chemistry_defs_json() -> Result<String> {
    let lib_bin_exe = std::env::current_exe()?;
    let lib_bin_dir = lib_bin_exe.parent().context("parent failed")?;
    let libdir = lib_bin_dir.parent().context("parent failed")?;
    Ok(std::fs::read_to_string(
        libdir.join("python/cellranger/chemistry_defs.json"),
    )?)
}

lazy_static! {
    static ref CHEMISTRY_DEF_MAP: HashMap<ChemistryName, ChemistryDef> = {
        if let Ok(s) = read_chemistry_defs_json() {
            // Use the chemistry_defs.json found in the tarball.
            serde_json::from_str(&s).unwrap()
        } else {
            // Cannot find chemistry_defs.json. Use the chemistry defs stored in the executable.
            serde_json::from_str(CHEMISTRY_DEFS_JSON_STR).unwrap()
        }
    };
}

/// Return the mapping of all static chemistry definitions.
pub fn known_chemistry_defs() -> &'static HashMap<ChemistryName, ChemistryDef> {
    &CHEMISTRY_DEF_MAP
}

/// Get a chemistry def from the static mapping.
pub fn get_chemistry_def(name: ChemistryName) -> Option<&'static ChemistryDef> {
    let maybe_def = CHEMISTRY_DEF_MAP.get(&name);
    if let Some(def) = maybe_def {
        assert_eq!(def.name, name);
    }
    maybe_def
}

/// Normalize non-functional fields of a def if it otherwise matches a single known def.
/// This is used to collapse a custom chemistry into a known static chemistry
/// if it is effectively a duplicate.
/// Return None if no match is found, or if more than once match is found.
pub fn normalize_chemistry_def(mut def: ChemistryDef) -> Option<ChemistryDef> {
    CHEMISTRY_DEF_MAP
        .values()
        .filter(|known_def| {
            // This is a quick hack to avoid having to define a reduced equality method.
            // I think this is advantageous because it specifies what you ignore
            // instead of what you include and is thus resilient to adding fields.
            def.name = known_def.name;
            def.description.clone_from(&known_def.description);
            def == **known_def
        })
        .exactly_one()
        .ok()
        .cloned()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_chem_defs_json_parse() {
        let map: HashMap<ChemistryName, ChemistryDef> =
            serde_json::from_str(CHEMISTRY_DEFS_JSON_STR).unwrap();
        for (name, def) in &map {
            assert_eq!(*name, def.name);
        }
    }
}
