//! Martian stage MULTI_SETUP_CHUNKS

use anyhow::Result;
use barcode::WhitelistSpec;
use cr_types::chemistry::ChemistryDef;
use cr_types::rna_read::{LegacyLibraryType, RnaChunk};
use cr_types::sample_def::SampleDef;
use itertools::Itertools;
use martian::prelude::*;
use martian_derive::{make_mro, MartianStruct};
use metric::TxHashMap;
use serde::{Deserialize, Serialize};

#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct MultiSetupChunksStageInputs {
    pub sample_id: String,
    pub sample_def: Vec<SampleDef>,
    pub chemistry_def: ChemistryDef,
    pub default_library_type: Option<LegacyLibraryType>,
}

#[derive(Clone, Serialize, Deserialize, MartianStruct, PartialEq)]
pub struct MultiSetupChunksStageOutputs {
    pub chunks: Vec<RnaChunk>,
    /// Gel Bead barcode whitelist or the spot barcode in Visium
    pub barcode_whitelist: Option<String>,
    pub visium_hd_slide_name: Option<String>,
}

// This is our stage struct
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
        assert!(
            args.sample_def.iter().all(|def| def.gem_group.is_some())
                || args.sample_def.iter().all(|def| def.gem_group.is_none()),
            "Inconsistent gem_group tags in sample def. Please specify all gem_group tags as null, or all gem_group tags with an integer."
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
                Ok(sample_def
                    .find_fastqs()?
                    .into_iter()
                    .map(move |fastqs| (sample_def, fastqs)))
            })
            .flatten_ok()
            .collect::<Result<_>>()?;

        let chunks: Vec<RnaChunk> = fastqs
            .into_iter()
            .enumerate()
            .map(|(chunk_id, (sample_def, fastqs))| {
                RnaChunk::new(
                    &args.chemistry_def,
                    sample_def,
                    args.default_library_type.unwrap_or_default(),
                    fastqs,
                    &args.sample_id,
                    library_id_map[&(sample_def.gem_group, sample_def.library_type)],
                    u16::try_from(chunk_id).unwrap(),
                )
            })
            .collect();

        Ok(MultiSetupChunksStageOutputs {
            chunks,
            barcode_whitelist: args
                .chemistry_def
                .barcode_whitelist()
                .option_gel_bead()
                .and_then(|spec| match spec {
                    WhitelistSpec::TxtFile { name } => Some(name.to_string()),
                    _ => None,
                }),
            visium_hd_slide_name: {
                let slide_names: Vec<_> = args
                    .chemistry_def
                    .barcode_whitelist()
                    .iter()
                    .filter_map(|spec| match spec {
                        WhitelistSpec::TxtFile { .. }
                        | WhitelistSpec::DynamicTranslation { .. } => None,
                        WhitelistSpec::SlideFile { slide, .. } => Some(slide),
                    })
                    .unique()
                    .collect();
                assert!(slide_names.len() <= 1);
                slide_names.into_iter().next().cloned()
            },
        })
    }
}

#[cfg(test)]
#[cfg(feature = "slow_tests")]
mod tests {
    use super::*;
    use insta::assert_ron_snapshot;
    use martian_filetypes::json_file::JsonFile;
    use martian_filetypes::FileTypeRead;
    use std::collections::BTreeMap;

    #[test]
    fn test_setup_chunks_outs() -> Result<()> {
        let data: BTreeMap<String, MultiSetupChunksStageInputs> =
            JsonFile::from("test/setup_chunks_sample_defs.json").read()?;
        let total_items = data.len();
        for (i, (sample, args)) in data.into_iter().enumerate() {
            println!("({}/{}) Sample: {}", i, total_items, sample);
            let outs = MultiSetupChunks.test_run_tmpdir(args)?;
            assert_ron_snapshot!(sample, outs);
        }
        Ok(())
    }
}
