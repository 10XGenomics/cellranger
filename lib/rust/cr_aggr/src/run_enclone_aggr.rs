//! RunEncloneAggr stage code

use crate::setup_vdj_aggr::{EncloneMetaRow, EncloneProtoMetaFormat};
use anyhow::Result;
use enclone_ranger::main_enclone::main_enclone_ranger;
use martian::prelude::*;
use martian_derive::{make_mro, martian_filetype, MartianStruct};
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::tabular_file::CsvFile;
use serde::{Deserialize, Serialize};
use std::path::PathBuf;
use vdj_ann::annotate::ContigAnnotation;

martian_filetype! {ProtoFile, "pb"}
martian_filetype! {FaFile, "fa"}

#[derive(Debug, Clone, Deserialize, MartianStruct)]
pub struct RunEncloneAggrStageInputs {
    // These files are read by enclone, via paths listed in the csv.
    #[allow(dead_code)]
    contig_ann_json_files: Vec<JsonFile<Vec<ContigAnnotation>>>, // One per VdjAggrSampleDef
    enclone_input_csv: CsvFile<EncloneMetaRow>,
    enclone_gem_well_meta: EncloneProtoMetaFormat,
    vdj_reference_path: PathBuf,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct RunEncloneAggrStageOutputs {
    enclone_output: ProtoFile,
    donor_ref_fa: FaFile,
}

// This is our stage struct
pub struct RunEncloneAggr;

#[make_mro(mem_gb = 9, threads = 4)]
impl MartianMain for RunEncloneAggr {
    type StageInputs = RunEncloneAggrStageInputs;
    type StageOutputs = RunEncloneAggrStageOutputs;
    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        // Call enclone.  We instruct it to generate a proto output file, and print nothing.
        // Need to bound number of threads used by enclone.

        let mut argsx = vec!["enclone".to_string()];
        let fasta_path = args.vdj_reference_path.join("fasta/regions.fa");

        let enclone_output: ProtoFile = rover.make_path("enclone_outputs");
        let donor_ref_fa: FaFile = rover.make_path("donor_ref");
        argsx.push("PRE=".to_string());
        argsx.push("CELLRANGER".to_string());
        argsx.push(format!(
            "META={}",
            args.enclone_input_csv.as_ref().to_str().unwrap()
        ));
        argsx.push(format!("REF={}", fasta_path.to_str().unwrap()));
        argsx.push("NOPRINT".to_string());
        argsx.push("NOPAGER".to_string());
        argsx.push("NOPRETTY".to_string());
        argsx.push("FORCE_EXTERNAL".to_string()); // do not test for internal run
        argsx.push(format!("MAX_CORES={}", rover.get_threads()));
        argsx.push(format!(
            "PROTO={}",
            enclone_output.as_ref().to_str().unwrap()
        ));
        argsx.push(format!(
            "PROTO_METADATA={}",
            args.enclone_gem_well_meta.as_ref().to_str().unwrap()
        ));
        argsx.push(format!(
            "DONOR_REF_FILE={}",
            donor_ref_fa.as_ref().to_str().unwrap()
        ));
        argsx.push(format!("MAX_CORES={}", rover.get_threads()));

        println!("{}", argsx.join(" "));
        main_enclone_ranger(&argsx).unwrap();

        println!("Done with enclone!");

        Ok(RunEncloneAggrStageOutputs {
            enclone_output,
            donor_ref_fa,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::process_vdj_proto::make_test_library;
    use crate::setup_vdj_aggr::{SetupVdjAggr, SetupVdjAggrStageInputs};
    use enclone_proto::proto_io::read_proto;
    use martian_filetypes::FileTypeRead;
    use vdj_reference::VdjReceptor;

    #[test]
    fn test_enclone_aggr() -> Result<()> {
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

        let tempdir1 = tempfile::tempdir()?;
        let write_outs = SetupVdjAggr.test_run(&tempdir1, write_args)?;

        let enclone_args = RunEncloneAggrStageInputs {
            contig_ann_json_files: write_outs.contig_ann_json_files,
            enclone_input_csv: write_outs.enclone_input_csv,
            enclone_gem_well_meta: write_outs.enclone_gem_well_meta.clone(),
            vdj_reference_path: write_outs.vdj_reference_path,
        };

        let tempdir2 = tempfile::tempdir()?;
        let outs = RunEncloneAggr.test_run(&tempdir2, enclone_args)?;

        let enclone_outputs = read_proto(outs.enclone_output)?;

        assert_eq!(
            enclone_outputs.metadata,
            write_outs.enclone_gem_well_meta.read()?
        );

        assert_eq!(
            enclone_outputs.universal_reference.items.len(),
            355 // grep "|IG|" regions.fa | wc -l
        );

        Ok(())
    }
}
