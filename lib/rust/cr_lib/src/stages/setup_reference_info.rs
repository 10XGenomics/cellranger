//! Martian stage SETUP_REFERENCE_INFO
//! Computes reference genomes and version from reference path or target set metadata so
//! that downstream stages do not to duplicate this logic (i.e. checking if we do not have
//! a reference path, then we should have a target set file).
#![deny(missing_docs)]

use anyhow::Result;
use cr_types::reference::probe_set_reference::TargetSetFile;
use cr_types::reference::reference_info::ReferenceInfo;
use martian::prelude::*;
use martian_derive::{MartianStruct, make_mro};
use serde::{Deserialize, Serialize};
use std::path::{Path, PathBuf};

/// Inputs for the SETUP_REFERENCE_INFO stage
#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct SetupReferenceInfoStageInputs {
    /// reference path used to extract list of reference genomes and version
    pub reference_path: Option<PathBuf>,
    /// probe target set CSV file used if the reference_path is not provided
    pub target_set: Option<TargetSetFile>,
}

/// Outputs for the SETUP_REFERENCE_INFO stage
#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct SetupReferenceInfoStageOutputs {
    /// ReferenceInfo built from reference_path or target_set
    pub reference_info: ReferenceInfo,
}

/// Martian stage SETUP_REFERENCE_INFO
pub struct SetupReferenceInfo;

#[make_mro(volatile = strict)]
impl MartianMain for SetupReferenceInfo {
    type StageInputs = SetupReferenceInfoStageInputs;
    type StageOutputs = SetupReferenceInfoStageOutputs;
    fn main(&self, args: Self::StageInputs, _rover: MartianRover) -> Result<Self::StageOutputs> {
        let reference_info =
            get_reference_info(args.reference_path.as_deref(), args.target_set.as_ref())?;
        Ok(SetupReferenceInfoStageOutputs { reference_info })
    }
}

/// Extract reference genomes and version from reference_path or target_set
pub fn get_reference_info(
    reference_path: Option<&Path>,
    target_set: Option<&TargetSetFile>,
) -> Result<ReferenceInfo> {
    println!("reference_path: {reference_path:?}");
    println!("target_set: {target_set:?}");
    let ref_info = if let Some(reference_path) = reference_path {
        ReferenceInfo::from_reference_path(reference_path)?
    } else if let Some(target_set) = target_set.as_ref() {
        ReferenceInfo::from_probe_set_csv(target_set)?
    } else {
        ReferenceInfo::default()
    };
    Ok(ref_info)
}

#[cfg(test)]
mod tests {
    use super::*;
    use metric::JsonReport;

    #[test]
    fn test_get_reference_info_probe_and_ref() -> Result<()> {
        insta::assert_snapshot!(
            get_reference_info(
                Some(Path::new("test/reference/GRCh38_ref_tiny")),
                Some(&TargetSetFile::from(
                    "test/probe_sets/GRCh38-fmt3-refv24.csv"
                )),
            )?
            .to_json_reporter()
        );
        Ok(())
    }

    #[test]
    fn test_get_reference_info_human_probe_only() -> Result<()> {
        insta::assert_snapshot!(
            get_reference_info(
                None,
                Some(&TargetSetFile::from(
                    "test/probe_sets/GRCh38-fmt3-refv24.csv"
                )),
            )?
            .to_json_reporter()
        );
        Ok(())
    }

    #[test]
    fn test_get_reference_info_mouse_probe_only() -> Result<()> {
        insta::assert_snapshot!(
            get_reference_info(
                None,
                Some(&TargetSetFile::from(
                    "test/probe_sets/GRCh38-fmt3-refv24.csv"
                )),
            )?
            .to_json_reporter()
        );
        Ok(())
    }

    #[test]
    fn test_get_reference_info_ref_only() -> Result<()> {
        insta::assert_snapshot!(
            get_reference_info(
                Some(Path::new("test/reference/GRCh38-and-mm10_ref_tiny")),
                None,
            )?
            .to_json_reporter()
        );
        Ok(())
    }

    #[test]
    fn test_get_reference_info_none() -> Result<()> {
        insta::assert_snapshot!(get_reference_info(None, None)?.to_json_reporter());
        Ok(())
    }
}
