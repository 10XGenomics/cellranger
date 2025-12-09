//! Martian stage WRITE_GENE_INDEX
#![deny(missing_docs)]

use anyhow::Result;
use martian::prelude::*;
use martian_derive::{MartianStruct, make_mro};
use martian_filetypes::json_file::JsonFile;
use serde::{Deserialize, Serialize};
use std::path::PathBuf;

/// Martian stage WRITE_GENE_INDEX
pub struct WriteGeneIndex;

#[derive(Deserialize, Clone, MartianStruct)]
pub struct StageInputs {
    pub(super) reference_path: Option<PathBuf>,
}

#[derive(Serialize, Deserialize, Clone, MartianStruct)]
pub struct StageOutputs {
    pub(super) gene_index: Option<JsonFile<()>>,
}

#[make_mro(mem_gb = 6, volatile = strict)]
impl MartianMain for WriteGeneIndex {
    type StageInputs = StageInputs;
    type StageOutputs = StageOutputs;

    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        let Some(reference_path) = &args.reference_path else {
            return Ok(StageOutputs { gene_index: None });
        };

        let gene_index: JsonFile<()> = rover.make_path("gene_index");
        transcriptome::python_gene_index::write_gene_index(reference_path, &gene_index)?;
        Ok(StageOutputs {
            gene_index: Some(gene_index),
        })
    }
}
