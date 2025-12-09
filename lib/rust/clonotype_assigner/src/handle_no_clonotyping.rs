//! HandleNoClonotyping stage code
#![expect(missing_docs)]

use anyhow::{Context, Result};
use martian::prelude::*;
use martian_derive::{MartianStruct, make_mro};
use martian_filetypes::json_file::JsonFile;
use serde::{Deserialize, Serialize};
use vdj_ann::annotate::ContigAnnotation;

#[derive(Debug, Clone, Deserialize, MartianStruct)]
pub struct ContigAnnotationSource {
    merged_liblevel: Option<JsonFile<Vec<ContigAnnotation>>>,
    post_cell_filtering: Option<JsonFile<Vec<ContigAnnotation>>>,
    post_clonotyping: Option<JsonFile<Vec<ContigAnnotation>>>,
}

impl ContigAnnotationSource {
    fn get_annotation(self, disable_clonotyping: bool) -> JsonFile<Vec<ContigAnnotation>> {
        match (
            self.merged_liblevel,
            self.post_cell_filtering,
            self.post_clonotyping,
            disable_clonotyping,
        ) {
            // Cellranger vdj OR Cellranger multi sample-level
            (None, Some(_), Some(annot), false) => annot,
            // Cellranger vdj denovo mode OR Cellranger multi sample-level with skip-clonotyping
            (None, Some(annot), None, true) => annot,
            // Cellranger multi library-level
            (Some(annot), None, None, true) => annot,
            (_, _, _, _) => unreachable!(
                "Incompatibe combination of contig annotations and disable_clonotyping boolean!"
            ),
        }
    }
}

#[derive(Debug, Clone, Deserialize, MartianStruct)]
pub struct HandleNoClonotypingStageInputs {
    contigs: ContigAnnotationSource,
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
        let src = args.contigs.get_annotation(args.disable_clonotyping);

        let dest: JsonFile<_> = rover.make_path("all_contig_annotations.json");
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
