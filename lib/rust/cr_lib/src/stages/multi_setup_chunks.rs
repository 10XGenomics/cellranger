//! Martian stage MULTI_SETUP_CHUNKS

use anyhow::{ensure, Result};
use barcode::WhitelistSpec;
use cr_types::chemistry::ChemistryDefs;
use cr_types::rna_read::RnaChunk;
use cr_types::sample_def::SampleDef;
use cr_types::LibraryType;
use itertools::Itertools;
use martian::prelude::*;
use martian_derive::{make_mro, MartianStruct};
use metric::TxHashMap;
use serde::{Deserialize, Serialize};

#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct MultiSetupChunksStageInputs {
    pub sample_id: String,
    pub sample_def: Vec<SampleDef>,
    pub chemistry_defs: ChemistryDefs,
    pub default_library_type: Option<LibraryType>,
}

#[derive(Clone, Serialize, Deserialize, MartianStruct, PartialEq)]
pub struct MultiSetupChunksStageOutputs {
    pub chunks: Vec<RnaChunk>,
    /// Gel Bead barcode whitelist(s) or the spot barcode in Visium
    pub barcode_whitelists: Vec<String>,
    pub visium_hd_slide_name: Option<String>,
}

pub struct MultiSetupChunks;

#[make_mro]
impl MartianMain for MultiSetupChunks {
    type StageInputs = MultiSetupChunksStageInputs;
    type StageOutputs = MultiSetupChunksStageOutputs; // Use `MartianVoid` if empty
    fn main(
        &self,
        mut args: Self::StageInputs,
        _rover: MartianRover,
    ) -> Result<Self::StageOutputs> {
        // Make sure that the gem groups in sample def is either all null or all some
        ensure!(
            args.sample_def.iter().all(|def| def.gem_group.is_some())
                || args.sample_def.iter().all(|def| def.gem_group.is_none()),
            "Inconsistent gem_group tags in sample def. \
             Please specify all gem_group tags as null, or all gem_group tags with an integer."
        );

        // If all gem_groups are set to null, then set them all to 1
        for def in &mut args.sample_def {
            def.gem_group = def.gem_group.or(Some(1));
        }

        let library_id_map: TxHashMap<_, _> = args
            .sample_def
            .iter()
            .map(|x| (x.gem_group, x.library_type))
            .unique()
            .enumerate()
            .map(|(library_id, gem_group_library_type)| (gem_group_library_type, library_id as u16))
            .collect();

        let fastqs: Vec<_> = args
            .sample_def
            .iter()
            .map(|sample_def| {
                anyhow::Ok(
                    sample_def
                        .find_fastqs()?
                        .into_iter()
                        .map(move |fastqs| (sample_def, fastqs)),
                )
            })
            .flatten_ok()
            .try_collect()?;

        let default_library_type = args.default_library_type.unwrap_or_default();
        let chunks: Vec<RnaChunk> = fastqs
            .into_iter()
            .enumerate()
            .map(|(chunk_id, (sample_def, fastqs))| {
                let library_type = sample_def.library_type.unwrap_or(default_library_type);
                RnaChunk::new(
                    &args.chemistry_defs[&library_type],
                    sample_def,
                    default_library_type,
                    fastqs,
                    &args.sample_id,
                    library_id_map[&(sample_def.gem_group, sample_def.library_type)],
                    u16::try_from(chunk_id).unwrap(),
                )
            })
            .collect();

        let barcode_whitelists = args
            .chemistry_defs
            .values()
            .filter_map(|x| {
                x.barcode_whitelist()
                    .option_gel_bead()
                    .and_then(WhitelistSpec::whitelist_name)
            })
            .unique()
            .map(String::from)
            .collect();

        let visium_hd_slide_name = args
            .chemistry_defs
            .values()
            .filter_map(|x| x.barcode_whitelist().map_option(WhitelistSpec::slide_name))
            .flatten()
            .dedup()
            .at_most_one()
            .unwrap()
            .map(String::from);

        Ok(MultiSetupChunksStageOutputs {
            chunks,
            barcode_whitelists,
            visium_hd_slide_name,
        })
    }
}
