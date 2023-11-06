//! HandleNoRef stage code

use anyhow::{Context, Result};
use martian::prelude::*;
use martian_derive::{make_mro, MartianStruct};
use martian_filetypes::json_file::JsonFile;
use serde::{Deserialize, Serialize};
use vdj_ann::annotate::ContigAnnotation;

#[derive(Debug, Clone, Deserialize, MartianStruct)]
pub struct HandleNoRefStageInputs {
    asm_contig_json: JsonFile<Vec<ContigAnnotation>>,
    clonotype_contig_json: Option<JsonFile<Vec<ContigAnnotation>>>,
    has_no_vdj_ref: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct HandleNoRefStageOutputs {
    final_contig_annotations: JsonFile<Vec<ContigAnnotation>>,
}

// This is our stage struct
pub struct HandleNoRef;

#[make_mro(stage_name = HANDLE_NO_VDJ_REF)]
impl MartianMain for HandleNoRef {
    type StageInputs = HandleNoRefStageInputs;
    type StageOutputs = HandleNoRefStageOutputs;
    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        assert!(args.has_no_vdj_ref ^ args.clonotype_contig_json.is_some());

        let dest: JsonFile<_> = rover.make_path("all_contig_annotations.json");

        let src = args.clonotype_contig_json.unwrap_or(args.asm_contig_json);

        std::fs::copy(&src, &dest).with_context(|| {
            format!(
                "Error: unable to copy {} to {}",
                &src.as_ref().display(),
                &dest.as_ref().display()
            )
        })?;

        Ok(HandleNoRefStageOutputs {
            final_contig_annotations: dest,
        })
    }
}
