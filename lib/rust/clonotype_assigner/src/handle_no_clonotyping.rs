//! HandleNoClonotyping stage code

use anyhow::{Context, Result};
use martian::prelude::*;
use martian_derive::{make_mro, MartianStruct};
use martian_filetypes::json_file::JsonFile;
use serde::{Deserialize, Serialize};
use vdj_ann::annotate::ContigAnnotation;

#[derive(Debug, Clone, Deserialize, MartianStruct)]
pub struct HandleNoClonotypingStageInputs {
    asm_contig_json: JsonFile<Vec<ContigAnnotation>>,
    clonotype_contig_json: Option<JsonFile<Vec<ContigAnnotation>>>,
    disable_clonotyping: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct HandleNoClonotypingStageOutputs {
    final_contig_annotations: JsonFile<Vec<ContigAnnotation>>,
}

// This is our stage struct
pub struct HandleNoClonotyping;

#[make_mro(stage_name = HANDLE_NO_CLONOTYPING)]
impl MartianMain for HandleNoClonotyping {
    type StageInputs = HandleNoClonotypingStageInputs;
    type StageOutputs = HandleNoClonotypingStageOutputs;
    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        assert!(args.disable_clonotyping ^ args.clonotype_contig_json.is_some());

        let dest: JsonFile<_> = rover.make_path("all_contig_annotations.json");

        let src = args.clonotype_contig_json.unwrap_or(args.asm_contig_json);

        std::fs::copy(&src, &dest).with_context(|| {
            format!(
                "Error: unable to copy {} to {}",
                &src.as_ref().display(),
                &dest.as_ref().display()
            )
        })?;

        Ok(HandleNoClonotypingStageOutputs {
            final_contig_annotations: dest,
        })
    }
}
