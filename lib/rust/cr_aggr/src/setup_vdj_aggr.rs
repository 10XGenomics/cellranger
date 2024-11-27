//! Martian stage SETUP_VDJ_AGGR

use crate::parse_aggr_csv::VdjAggrCsvLibrary;
use anyhow::Result;
use itertools::Itertools;
use martian::prelude::*;
use martian_derive::{make_mro, martian_filetype, MartianStruct, MartianType};
use martian_filetypes::json_file::{JsonFile, JsonFormat};
use martian_filetypes::tabular_file::CsvFile;
use martian_filetypes::{FileTypeWrite, LazyFileTypeIO, LazyWrite};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap};
use std::path::PathBuf;
use vdj_ann::annotate::ContigAnnotation;
use vdj_proto::io::VdjProtoReader;
use vdj_reference::VdjReceptor;

#[derive(Debug, Clone, Serialize, Deserialize, MartianType)]
pub struct EncloneMetaRow {
    #[serde(skip_serializing_if = "Option::is_none")]
    tcr: Option<PathBuf>,
    #[serde(skip_serializing_if = "Option::is_none")]
    tcrgd: Option<PathBuf>,
    #[serde(skip_serializing_if = "Option::is_none")]
    bcr: Option<PathBuf>,
    #[serde(skip_serializing_if = "Option::is_none")]
    gex: Option<PathBuf>,
    origin: String,
    donor: String,
}

martian_filetype! { _EncloneProtoMeta, "em" }
pub type EncloneProtoMetaFormat = JsonFormat<_EncloneProtoMeta, enclone_proto::types::Metadata>;

#[derive(Debug, Clone, Deserialize, MartianStruct)]
pub struct SetupVdjAggrStageInputs {
    pub libraries: Vec<VdjAggrCsvLibrary>,
    #[mro_type = "map"] // TODO: Fix this in martian
    pub gem_well_map: BTreeMap<u32, (String, u32)>,
    pub receptor: VdjReceptor,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct SetupVdjAggrStageOutputs {
    pub contig_ann_json_files: Vec<JsonFile<Vec<ContigAnnotation>>>, // One per VdjAggrSampleDef
    pub enclone_input_csv: CsvFile<EncloneMetaRow>,
    pub enclone_gem_well_meta: EncloneProtoMetaFormat,
    pub vdj_reference_path: PathBuf,
    pub combined_ann_json: JsonFile<Vec<ContigAnnotation>>,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct SetupVdjAggrChunkInputs {
    chunk_id: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct SetupVdjAggrChunkOutputs {
    chunk_ann_json: JsonFile<Vec<ContigAnnotation>>,
    enclone_meta_row: EncloneMetaRow,
    #[mro_type = "map"]
    enclone_gem_well_info: HashMap<u32, enclone_proto::types::GemWellInfo>,
}

pub struct SetupVdjAggr;

// This is hardcoded in enclone!
const ANNOTATIONS_JSON_FILENAME: &str = "contig_annotations";

#[make_mro]
impl MartianStage for SetupVdjAggr {
    type StageInputs = SetupVdjAggrStageInputs;
    type StageOutputs = SetupVdjAggrStageOutputs;
    type ChunkInputs = SetupVdjAggrChunkInputs;
    type ChunkOutputs = SetupVdjAggrChunkOutputs;

    fn split(
        &self,
        args: Self::StageInputs,
        _rover: MartianRover,
    ) -> Result<StageDef<Self::ChunkInputs>> {
        Ok((0..args.libraries.len())
            .map(|chunk_id| SetupVdjAggrChunkInputs { chunk_id })
            .collect())
    }

    fn main(
        &self,
        args: Self::StageInputs,
        chunk_args: Self::ChunkInputs,
        rover: MartianRover,
    ) -> Result<Self::ChunkOutputs> {
        let sdef = &args.libraries[chunk_args.chunk_id];

        let new_gem_well_map: HashMap<_, _> = args
            .gem_well_map
            .into_iter()
            .filter(|(_, (lib_id, _))| lib_id == &sdef.library_id)
            .map(|(new_gw, (_, old_gw))| (old_gw, new_gw))
            .collect();

        let chunk_ann_json: JsonFile<_> = rover.make_path(ANNOTATIONS_JSON_FILENAME);
        let mut ann_writer = chunk_ann_json.lazy_writer()?;

        for ann in VdjProtoReader::read_annotations(&sdef.vdj_contig_info)? {
            let mut ann = ann?;
            // Renumber gem well. This is a bit messy because of the string representations
            let old_barcode = ann.barcode.clone();
            let barcode_parts: Vec<_> = old_barcode.split('-').collect();
            let old_gw: u32 = barcode_parts[1].parse()?;
            let new_gw = new_gem_well_map[&old_gw];
            ann.barcode = format!("{}-{new_gw}", barcode_parts[0]);
            ann.contig_name = ann.contig_name.replace(&old_barcode, &ann.barcode);

            // TODO: Do we need to reset these?
            // ann.clonotype = None;
            // ann.info = Default::default();

            ann_writer.write_item(&ann)?;
        }

        ann_writer.finish()?;

        let mut enclone_gem_well_info = HashMap::new();
        for new_gw in new_gem_well_map.values() {
            enclone_gem_well_info.insert(
                *new_gw,
                enclone_proto::types::GemWellInfo {
                    library_id: sdef.library_id.clone(),
                    donor: sdef.donor.clone(),
                    origin: sdef.origin.clone(),
                    additional_data: sdef.meta.clone(),
                },
            );
        }

        let files_path = rover.files_path().to_path_buf();
        let (bcr, tcr, tcrgd) = match args.receptor {
            VdjReceptor::TR => (None, Some(files_path), None),
            VdjReceptor::IG => (Some(files_path), None, None),
            VdjReceptor::TRGD => (None, None, Some(files_path)),
        };

        Ok(SetupVdjAggrChunkOutputs {
            chunk_ann_json,
            enclone_meta_row: EncloneMetaRow {
                tcr,
                tcrgd,
                bcr,
                gex: None,
                origin: sdef.origin.clone(),
                donor: sdef.donor.clone(),
            },
            enclone_gem_well_info,
        })
    }

    fn join(
        &self,
        args: Self::StageInputs,
        _chunk_defs: Vec<Self::ChunkInputs>,
        chunk_outs: Vec<Self::ChunkOutputs>,
        rover: MartianRover,
    ) -> Result<Self::StageOutputs> {
        let vdj_reference_path: PathBuf = rover.make_path("reference");
        VdjProtoReader::read_reference(&args.libraries[0].vdj_contig_info)?
            .write_to_folder(&vdj_reference_path)?;

        let enclone_input_csv: CsvFile<_> = rover.make_path("enclone_meta");
        enclone_input_csv.write(
            &chunk_outs
                .iter()
                .map(|co| co.enclone_meta_row.clone())
                .collect(),
        )?;

        let enclone_gem_well_meta: EncloneProtoMetaFormat =
            rover.make_path("enclone_gem_well_meta");
        let per_gem_well_info: HashMap<u32, enclone_proto::types::GemWellInfo> = chunk_outs
            .iter()
            .flat_map(|co| co.enclone_gem_well_info.clone())
            .collect();
        enclone_gem_well_meta.write(&enclone_proto::types::Metadata {
            additional_columns: args.libraries[0].meta.keys().cloned().collect(),
            donors: args
                .libraries
                .into_iter()
                .map(|sdef| sdef.donor)
                .unique()
                .collect(),
            per_gem_well_info,
        })?;

        let contig_ann_json_files: Vec<_> =
            chunk_outs.into_iter().map(|co| co.chunk_ann_json).collect();

        // Create the merged annotation JSON. Would be better to use a binary format.
        let combined_ann_json: JsonFile<_> = rover.make_path("combined_contig_annotations");
        let mut contig_writer = combined_ann_json.lazy_writer()?;
        for ann_json in &contig_ann_json_files {
            for ann in ann_json.lazy_reader()? {
                let ann: ContigAnnotation = ann?;
                contig_writer.write_item(&ann)?;
            }
        }

        Ok(SetupVdjAggrStageOutputs {
            contig_ann_json_files,
            enclone_input_csv,
            enclone_gem_well_meta,
            vdj_reference_path,
            combined_ann_json,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::process_vdj_proto::make_test_library;
    use insta::assert_yaml_snapshot;
    use martian_filetypes::FileTypeRead;
    use pretty_assertions::assert_eq;

    #[test]
    fn test_gem_well_renumbering() -> Result<()> {
        let write_args = SetupVdjAggrStageInputs {
            libraries: vec![
                make_test_library("mouse_spleno_ig_1"),
                make_test_library("mouse_spleno_ig_2"),
                make_test_library("mouse_spleno_ig_3"),
            ],
            gem_well_map: vec![
                (1u32, ("mouse_spleno_ig_1".to_string(), 1u32)),
                (2u32, ("mouse_spleno_ig_2".to_string(), 1u32)),
                (3u32, ("mouse_spleno_ig_3".to_string(), 1u32)),
            ]
            .into_iter()
            .collect(),
            receptor: VdjReceptor::IG,
        };

        let tempdir = tempfile::tempdir()?;
        let write_outs = SetupVdjAggr.test_run(&tempdir, write_args)?;
        let mut expected_csv_rows = vec![String::from("bcr,origin,donor")];
        for (i, ann_json) in write_outs.contig_ann_json_files.into_iter().enumerate() {
            let anns: Vec<vdj_ann::annotate::ContigAnnotation> = ann_json.read()?;
            expected_csv_rows.push(format!(
                "{},Origin,Donor",
                ann_json.as_ref().parent().unwrap().display()
            ));
            let gem_well_suffix = format!("-{}", i + 1);
            for ann in anns {
                assert!(
                    ann.barcode.ends_with(&gem_well_suffix),
                    "Barcode {} does not end in {}",
                    ann.barcode,
                    gem_well_suffix
                );
                assert!(
                    ann.contig_name.contains(&ann.barcode),
                    "Barcode {} is not contained in {}",
                    ann.barcode,
                    ann.contig_name
                );
            }
        }

        assert_eq!(
            expected_csv_rows.join("\n"),
            std::fs::read_to_string(write_outs.enclone_input_csv)?.trim_end()
        );

        let donor = "Donor".to_string();
        let gw_info = |id: &str| enclone_proto::types::GemWellInfo {
            library_id: id.to_string(),
            donor: donor.clone(),
            origin: "Origin".into(),
            additional_data: Default::default(),
        };

        let expected_metadata = enclone_proto::types::Metadata {
            additional_columns: vec![],
            donors: vec![donor.clone()],
            per_gem_well_info: vec![
                (1, gw_info("mouse_spleno_ig_1")),
                (2, gw_info("mouse_spleno_ig_2")),
                (3, gw_info("mouse_spleno_ig_3")),
            ]
            .into_iter()
            .collect(),
        };

        let actual_metadata = write_outs.enclone_gem_well_meta.read()?;

        assert_eq!(expected_metadata, actual_metadata);

        Ok(())
    }

    #[test]
    fn test_meta_csv_tcr() -> Result<()> {
        let write_args = SetupVdjAggrStageInputs {
            libraries: vec![
                make_test_library("human_pbmc_tcr_1"),
                make_test_library("human_pbmc_tcr_2"),
            ],
            gem_well_map: vec![
                (1u32, ("human_pbmc_tcr_1".to_string(), 1u32)),
                (2u32, ("human_pbmc_tcr_2".to_string(), 1u32)),
            ]
            .into_iter()
            .collect(),
            receptor: VdjReceptor::TR,
        };

        let tempdir = tempfile::tempdir()?;
        let write_outs = SetupVdjAggr.test_run(&tempdir, write_args)?;
        let mut expected_csv_rows = vec![String::from("tcr,origin,donor")];
        for ann_json in write_outs.contig_ann_json_files {
            expected_csv_rows.push(format!(
                "{},Origin,Donor",
                ann_json.as_ref().parent().unwrap().display()
            ));
        }

        assert_eq!(
            expected_csv_rows.join("\n"),
            std::fs::read_to_string(write_outs.enclone_input_csv)?.trim_end()
        );
        Ok(())
    }

    #[test]
    fn test_gem_well_info() -> Result<()> {
        let write_args = SetupVdjAggrStageInputs {
            libraries: vec![
                VdjAggrCsvLibrary {
                    library_id: "human_pbmc_tcr_1".into(),
                    vdj_contig_info: "test_resources/vdj_contig_info/human_pbmc_tcr_1.pb".into(),
                    donor: "Donor1".into(),
                    origin: "Origin1".into(),
                    meta: vec![("Timepoint".to_string(), "T1".to_string())]
                        .into_iter()
                        .collect(),
                },
                VdjAggrCsvLibrary {
                    library_id: "human_pbmc_tcr_2".into(),
                    vdj_contig_info: "test_resources/vdj_contig_info/human_pbmc_tcr_2.pb".into(),
                    donor: "Donor2".into(),
                    origin: "Origin2".into(),
                    meta: vec![("Timepoint".to_string(), "T2".to_string())]
                        .into_iter()
                        .collect(),
                },
            ],
            gem_well_map: vec![
                (1u32, ("human_pbmc_tcr_1".to_string(), 1u32)),
                (2u32, ("human_pbmc_tcr_2".to_string(), 1u32)),
            ]
            .into_iter()
            .collect(),
            receptor: VdjReceptor::TR,
        };

        let tempdir = tempfile::tempdir()?;
        let outs = SetupVdjAggr.test_run(&tempdir, write_args)?;

        let metadata = outs.enclone_gem_well_meta.read()?;

        let mut settings = insta::Settings::clone_current();
        settings.set_sort_maps(true);
        settings.bind(|| {
            assert_yaml_snapshot!(&metadata);
        });

        Ok(())
    }
}
