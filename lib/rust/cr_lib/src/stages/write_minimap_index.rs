//! Martian stage WRITE_MINIMAP_INDEX
#![deny(missing_docs)]

use anyhow::{Result, bail};
use cr_types::AlignerParam;
use martian::prelude::*;
use martian_derive::{MartianStruct, make_mro, martian_filetype};
use martian_filetypes::{LazyFileTypeIO, LazyWrite};
use minimap2::{Aligner, Preset};
use serde::{Deserialize, Serialize};
use std::path::PathBuf;
use transcriptome::{Bed12Format, Transcriptome};

martian_filetype!(MinimapIndex, "mmi");

/// Write minimap2 genome index.
pub struct WriteMinimapIndex;

/// The Martian stage inputs.
#[derive(Deserialize, Clone, MartianStruct)]
pub struct StageInputs {
    pub reference_path: PathBuf,
    pub aligner: AlignerParam,
}

/// The Martian stage outputs.
#[derive(Serialize, Deserialize, Clone, MartianStruct)]
pub struct StageOutputs {
    pub minimap_index: MinimapIndex,
    pub tx_bed: Bed12Format,
}

#[make_mro(mem_gb = 24, threads = 6, volatile = strict)]
impl MartianMain for WriteMinimapIndex {
    type StageInputs = StageInputs;
    type StageOutputs = StageOutputs;

    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        assert!(args.aligner == AlignerParam::Minimap2);
        let fasta_path = args.reference_path.join("fasta/genome.fa");
        if !fasta_path.exists() {
            bail!("FASTA file doesn't exist: {}", fasta_path.display());
        }
        let minimap_index: MinimapIndex = rover.make_path("genome");
        Aligner::builder()
            .preset(Preset::Splice)
            .with_index_threads(rover.get_threads())
            .with_index(fasta_path, minimap_index.as_ref().to_str())
            .expect("Unable to build minimap2 index");

        // Generate BED12 file
        let bed_file: Bed12Format = rover.make_path("tx_bed");
        let mut bed_writer = bed_file.lazy_writer()?;
        for tx in Transcriptome::from_reference_path(&args.reference_path)
            .unwrap()
            .convert_to_bed12()
        {
            bed_writer.write_item(&tx)?;
        }
        bed_writer.finish()?;
        Ok(StageOutputs {
            minimap_index,
            tx_bed: bed_file,
        })
    }
}
