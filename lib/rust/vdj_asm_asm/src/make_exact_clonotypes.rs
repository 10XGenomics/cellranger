//! MakeExactClonotypes stage code

use anyhow::Result;
use martian::prelude::*;
use martian_derive::{make_mro, MartianStruct};
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::FileTypeWrite;
use serde::{Deserialize, Serialize};
use vdj_ann::annotate::ContigAnnotation;
use vdj_asm_utils::exact_clonotyping::{generate_exact_clonotypes, ExactClonotype};

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct MakeExactClonotypesStageInputs {
    pub contig_annotations: JsonFile<Vec<ContigAnnotation>>,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct MakeExactClonotypesStageOutputs {
    pub exact_clonotypes: JsonFile<Vec<ExactClonotype>>,
}

pub struct MakeExactClonotypes;

#[make_mro(mem_gb = 4)]
impl MartianMain for MakeExactClonotypes {
    type StageInputs = MakeExactClonotypesStageInputs;
    type StageOutputs = MakeExactClonotypesStageOutputs;
    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        let exact_clonotypes = generate_exact_clonotypes(args.contig_annotations.clone())?;
        let exact_clonotypes_fn: JsonFile<_> = rover.make_path("exact_clonotypes");
        exact_clonotypes_fn.write(&exact_clonotypes)?;

        Ok(MakeExactClonotypesStageOutputs {
            exact_clonotypes: exact_clonotypes_fn,
        })
    }
}
