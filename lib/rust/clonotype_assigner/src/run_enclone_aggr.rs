//! RunEncloneAggr stage code

use anyhow::Result;
use cr_aggr::setup_vdj_aggr::{EncloneMetaRow, EncloneProtoMetaFormat};
use enclone_ranger::main_enclone::{main_enclone_ranger, CellrangerOpt, InputSpec};
use martian::prelude::*;
use martian_derive::{make_mro, martian_filetype, MartianStruct};
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::tabular_file::CsvFile;
use serde::{Deserialize, Serialize};
use std::path::PathBuf;
use vdj_ann::annotate::ContigAnnotation;
use vdj_filter_barcodes::filter_log::FilterSwitch;

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
    filter_switch: FilterSwitch,
    mix_donors: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct RunEncloneAggrStageOutputs {
    enclone_output: ProtoFile,
    donor_ref_fa: FaFile,
    barcode_fate: JsonFile<()>,
}

// This is our stage struct
pub struct RunEncloneAggr;

#[make_mro(mem_gb = 16, threads = 4)]
impl MartianMain for RunEncloneAggr {
    type StageInputs = RunEncloneAggrStageInputs;
    type StageOutputs = RunEncloneAggrStageOutputs;
    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        let enclone_output: ProtoFile = rover.make_path("enclone_outputs");
        let donor_ref_fa: FaFile = rover.make_path("donor_ref");
        let barcode_fate: JsonFile<()> = rover.make_path("barcode_fate");

        let mut opt = CellrangerOpt {
            input: Some(InputSpec::MetaFile(
                args.enclone_input_csv.to_str().unwrap().to_string(),
            )),
            refname: args
                .vdj_reference_path
                .join("fasta/regions.fa")
                .to_str()
                .unwrap()
                .to_string(),
            max_cores: Some(rover.get_threads()),
            proto: enclone_output.to_str().unwrap().to_string(),
            proto_metadata: args.enclone_gem_well_meta.to_str().unwrap().to_string(),
            dref_file: donor_ref_fa.to_str().unwrap().to_string(),
            fate_file: barcode_fate.to_str().unwrap().to_string(),
            // Option to split clonotypes that have 4 chains or more.
            // These are most likely false joins due to a shared chain
            split_max_chains: 4,
            mix_donors: args.mix_donors,
            ..Default::default()
        };

        args.filter_switch.update_enclone_args(&mut opt.filter);

        println!("Invoking enclone:\n{opt:?}");
        main_enclone_ranger(opt).unwrap();

        println!("Done with enclone!");

        Ok(RunEncloneAggrStageOutputs {
            enclone_output,
            donor_ref_fa,
            barcode_fate,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use cr_aggr::process_vdj_proto::make_test_library;
    use cr_aggr::setup_vdj_aggr::{SetupVdjAggr, SetupVdjAggrStageInputs};
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
            filter_switch: FilterSwitch::default(),
            mix_donors: false,
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
