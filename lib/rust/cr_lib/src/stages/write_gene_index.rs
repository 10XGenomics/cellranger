//! Martian stage WRITE_GENE_INDEX

use anyhow::Result;
use martian::prelude::*;
use martian_derive::{make_mro, MartianStruct};
use martian_filetypes::json_file::JsonFile;
use serde::{Deserialize, Serialize};
use std::path::PathBuf;

pub struct WriteGeneIndex;

#[derive(Deserialize, Clone, MartianStruct)]
pub struct StageInputs {
    pub(crate) reference_path: PathBuf,
}

#[derive(Serialize, Deserialize, Clone, MartianStruct)]
pub struct StageOutputs {
    pub(crate) gene_index: JsonFile<()>,
}

#[make_mro(mem_gb = 6, volatile = strict)]
impl MartianMain for WriteGeneIndex {
    type StageInputs = StageInputs;
    type StageOutputs = StageOutputs;

    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        let gene_index: JsonFile<()> = rover.make_path("gene_index");
        transcriptome::python_gene_index::write_gene_index(&args.reference_path, &gene_index)?;
        Ok(StageOutputs { gene_index })
    }
}
