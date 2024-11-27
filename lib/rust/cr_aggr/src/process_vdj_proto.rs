//! Martian stage PROCESS_VDJ_PROTO

use crate::parse_aggr_csv::VdjAggrCsvLibrary;
use anyhow::Result;
use itertools::Itertools;
use martian::prelude::*;
use martian_derive::{make_mro, MartianStruct};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap};
use std::fmt;
use vdj_proto::io::VdjProtoReader;
use vdj_proto::types::VdjMetadata;
use vdj_reference::VdjReceptor;

#[derive(Debug, thiserror::Error)]
enum AggrInconsistency {
    MultipleReferenceVersions {
        libs_of_ref: BTreeMap<String, Vec<String>>,
    },
    MultipleReceptors {
        libs_of_receptor: BTreeMap<VdjReceptor, Vec<String>>,
    },
}

impl fmt::Display for AggrInconsistency {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        use AggrInconsistency::{MultipleReceptors, MultipleReferenceVersions};
        let msg = match self {
            MultipleReferenceVersions { libs_of_ref } => {
                let refs_display = libs_of_ref
                    .iter()
                    .map(|(k, library_ids)| {
                        let libs = library_ids
                            .iter()
                            .map(|lib| format!(" - {lib}"))
                            .join("\n");
                        format!(
                            "The following library ID(s) use identical reference (with hash {k}):\n{libs}\n"
                        )
                    })
                    .join("\n");
                format!(
                    "The V(D)J datasets you are trying to aggregate were created with different \
                    references, but the 'aggr' command requires identical \
                    references in order to combine the datasets. Please re-run the pipeline with \
                    uniform reference in order to aggregate these data.\n\n{refs_display}"
                )
            }
            MultipleReceptors { libs_of_receptor } => {
                let receptor_display = libs_of_receptor
                    .iter()
                    .map(|(receptor, library_ids)| {
                        let libs = library_ids.iter().map(|lib| format!(" - {lib}")).join("\n");
                        format!(
                            "The following library ID(s) contain data from {receptor} chain:\n{libs}\n"
                        )
                    })
                    .join("\n");
                format!(
                    "The V(D)J datasets you are trying to aggregate were created with different \
                    chain types (i.e some are T cell based while some are B cell based), but the \
                    'aggr' command requires a single chain type in order to combine the datasets. \
                    Please check your input data.\n\n{receptor_display}"
                )
            }
        };
        write!(f, "{msg}")
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct ProcessVdjProtoStageInputs {
    pub libraries: Vec<VdjAggrCsvLibrary>,
    // Map from new_gem_well -> (library_id, old_gem_well)
    // as computed by the count aggr pipeline. We will reuse it
    // if available.
    #[mro_type = "map"]
    pub count_gem_well_map: Option<HashMap<u32, (String, u32)>>,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct ProcessVdjProtoStageOutputs {
    pub receptor: VdjReceptor,
    // Map from new_gem_well -> (library_id, old_gem_well)
    #[mro_type = "map"] // TODO: Fix this in martian
    pub gem_well_map: BTreeMap<u32, (String, u32)>,
}

// This is our stage struct
pub struct ProcessVdjProto;

#[make_mro]
impl MartianMain for ProcessVdjProto {
    type StageInputs = ProcessVdjProtoStageInputs;
    type StageOutputs = ProcessVdjProtoStageOutputs;
    fn main(&self, args: Self::StageInputs, _rover: MartianRover) -> Result<Self::StageOutputs> {
        assert!(!args.libraries.is_empty());

        let metadata_map: BTreeMap<String, VdjMetadata> = {
            let mut map = BTreeMap::new();
            for def in &args.libraries {
                map.insert(
                    def.library_id.clone(),
                    VdjProtoReader::read_metadata(&def.vdj_contig_info)?,
                );
            }
            map
        };

        // --------------------------------------------------------------------
        // Ensure that the same version of reference was used
        // We will populate a hashmap from reference_hash to library ids so that
        // we can report back a useful error message
        let libs_of_ref: BTreeMap<String, Vec<String>> = {
            let mut map = BTreeMap::new();
            for (lib, meta) in &metadata_map {
                if !map.contains_key(&meta.reference_fasta_hash) {
                    map.insert(meta.reference_fasta_hash.clone(), Vec::new());
                }
                map.get_mut(&meta.reference_fasta_hash)
                    .unwrap()
                    .push(lib.clone());
            }
            map
        };
        if libs_of_ref.len() > 1 {
            return Err(AggrInconsistency::MultipleReferenceVersions { libs_of_ref }.into());
        }
        drop(libs_of_ref);

        // --------------------------------------------------------------------
        // Ensure that it's the same receptor in all sample defs
        let libs_of_receptor: BTreeMap<VdjReceptor, Vec<String>> = {
            let mut map = BTreeMap::new();
            for (lib, meta) in &metadata_map {
                let receptor = meta.vdj_receptor();
                map.entry(receptor).or_insert_with(Vec::new);
                map.get_mut(&receptor).unwrap().push(lib.clone());
            }
            map
        };
        if libs_of_receptor.len() > 1 {
            return Err(AggrInconsistency::MultipleReceptors { libs_of_receptor }.into());
        }

        let receptor = *libs_of_receptor.keys().next().unwrap();

        let gem_well_map = match args.count_gem_well_map {
            Some(map) => map.into_iter().collect(),
            None => args
                .libraries
                .into_iter()
                .flat_map(|def| {
                    let lib_id = def.library_id;
                    metadata_map[&lib_id]
                        .gem_wells
                        .clone()
                        .into_iter()
                        .sorted()
                        .map(move |gw| (lib_id.clone(), gw))
                })
                .enumerate()
                .map(|(i, (id, gw))| ((i + 1) as u32, (id, gw)))
                .collect(),
        };

        Ok(ProcessVdjProtoStageOutputs {
            receptor,
            gem_well_map,
        })
    }
}

pub fn make_test_library(key: &'static str) -> VdjAggrCsvLibrary {
    VdjAggrCsvLibrary {
        library_id: key.into(),
        vdj_contig_info: format!("../cr_aggr/test_resources/vdj_contig_info/{key}.pb").into(),
        donor: "Donor".into(),
        origin: "Origin".into(),
        meta: HashMap::new(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use dui_tests::stage_test::StageFailTest;
    use dui_tests::{stage_fail_dui_test, DuiTest};
    use insta::assert_yaml_snapshot;

    #[test]
    fn test_conflicting_reference() -> Result<()> {
        stage_fail_dui_test!(
            VdjAggrConflictingReference,
            description: "User inputs proto files generated with different reference versions. \
            This is not supported in VDJ aggr.",
            stage: ProcessVdjProto,
            args: ProcessVdjProtoStageInputs {
                libraries: vec![
                    make_test_library("human_pbmc_ig_1"),
                    make_test_library("human_pbmc_ig_2"),
                    make_test_library("mouse_spleno_ig_1"),
                ],
                count_gem_well_map: None,
            },
        );
        Ok(())
    }

    #[test]
    fn test_conflicting_receptor() -> Result<()> {
        stage_fail_dui_test!(
            VdjAggrConflictingReceptor,
            description: "User inputs a mixture of TCR and IG proto files for VDJ aggr.",
            stage: ProcessVdjProto,
            args: ProcessVdjProtoStageInputs {
                libraries: vec![
                    make_test_library("human_pbmc_tcr_1"),
                    make_test_library("human_pbmc_tcr_2"),
                    make_test_library("human_pbmc_ig_1"),
                ],
                count_gem_well_map: None,
            },
        );
        Ok(())
    }

    #[test]
    fn test_with_count_gem_well_map() -> Result<()> {
        let count_gem_well_map = serde_json::from_str(
            r#"{
    "1": [
        "human_pbmc_tcr_1",
        1
    ],
    "2": [
        "human_pbmc_tcr_2",
        1
    ]
}"#,
        )?;
        let args = ProcessVdjProtoStageInputs {
            libraries: vec![
                make_test_library("human_pbmc_tcr_1"),
                make_test_library("human_pbmc_tcr_2"),
            ],
            count_gem_well_map: Some(count_gem_well_map),
        };
        let outs = ProcessVdjProto.test_run_tmpdir(args)?;
        assert_yaml_snapshot!(&outs);
        Ok(())
    }

    #[test]
    fn test_without_count_gem_well_map() -> Result<()> {
        let args = ProcessVdjProtoStageInputs {
            libraries: vec![
                make_test_library("mouse_spleno_ig_1"),
                make_test_library("mouse_spleno_ig_2"),
                make_test_library("mouse_spleno_ig_3"),
            ],
            count_gem_well_map: None,
        };
        let outs = ProcessVdjProto.test_run_tmpdir(args)?;
        assert_yaml_snapshot!(&outs);
        Ok(())
    }
}
