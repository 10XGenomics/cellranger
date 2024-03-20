//! Martian stage DETECT_CHEMISTRY_TEST

use crate::stages::detect_chemistry::{DetectChemistry, DetectChemistryStageInputs};
use anyhow::Result;
use cr_types::chemistry::ChemistryName;
use cr_types::LibraryType;
use itertools::Itertools;
use martian::prelude::*;
use martian_derive::{make_mro, MartianStruct};
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::{FileTypeRead, FileTypeWrite};
use metric::TxHashMap;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::path::PathBuf;

#[derive(Clone, Deserialize, MartianStruct)]
pub struct StageInputs {
    inputs: JsonFile<TxHashMap<String, DetectChemistryStageInputs>>,
    expected: JsonFile<TxHashMap<String, TxHashMap<LibraryType, ChemistryName>>>,
}

#[derive(Clone, Serialize, MartianStruct)]
pub struct StageOutputs {
    actual: JsonFile<BTreeMap<String, Option<TxHashMap<LibraryType, ChemistryName>>>>,
}

#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct ChunkInputs {
    chunk_inputs: TxHashMap<String, DetectChemistryStageInputs>,
}

#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct ChunkOutputs {
    chunk_outputs: TxHashMap<String, Option<TxHashMap<LibraryType, ChemistryName>>>,
}

pub struct DetectChemistryTest;

const ITEMS_PER_CHUNK: usize = 100;

#[make_mro(mem_gb = 32)]
impl MartianStage for DetectChemistryTest {
    type StageInputs = StageInputs;
    type StageOutputs = StageOutputs;
    type ChunkInputs = ChunkInputs;
    type ChunkOutputs = ChunkOutputs;

    fn split(
        &self,
        args: Self::StageInputs,
        _rover: MartianRover,
    ) -> Result<StageDef<Self::ChunkInputs>> {
        let inputs = args.inputs.read()?;
        let _expected = args.expected.read()?;
        Ok(inputs
            .into_iter()
            // .filter(|(k, _)| expected.contains_key(k))
            .chunks(ITEMS_PER_CHUNK)
            .into_iter()
            .map(|chunks| ChunkInputs {
                chunk_inputs: chunks.collect(),
            })
            .collect())
    }

    fn main(
        &self,
        _args: Self::StageInputs,
        chunk_args: Self::ChunkInputs,
        rover: MartianRover,
    ) -> Result<Self::ChunkOutputs> {
        let mut result = TxHashMap::default();
        for (k, v) in chunk_args.chunk_inputs {
            let run_dir: PathBuf = rover.make_path(&k);
            std::fs::create_dir(&run_dir)?;
            println!("RUNNING {k:?}");
            match DetectChemistry.test_run(&run_dir, v) {
                Ok(outs) => {
                    result.insert(
                        k,
                        Some(
                            outs.chemistry_defs
                                .into_iter()
                                .map(|(lib_type, def)| (lib_type, def.name))
                                .collect(),
                        ),
                    );
                }
                Err(e) => {
                    eprintln!("> FAILURE for {k}: {e}");
                    result.insert(k, None);
                }
            }
        }
        Ok(ChunkOutputs {
            chunk_outputs: result,
        })
    }
    fn join(
        &self,
        args: Self::StageInputs,
        _chunk_defs: Vec<Self::ChunkInputs>,
        chunk_outs: Vec<Self::ChunkOutputs>,
        rover: MartianRover,
    ) -> Result<Self::StageOutputs> {
        let outputs: BTreeMap<_, _> = chunk_outs
            .into_iter()
            .flat_map(|co| co.chunk_outputs)
            .collect();

        let actual: JsonFile<_> = rover.make_path("actual");
        actual.write(&outputs)?;

        let expected = args.expected.read()?;

        for (k, actual) in outputs {
            match (actual, expected.get(&k)) {
                (Some(a), Some(e)) => {
                    if a != *e {
                        println!("MISMATCH for k = {k} actual = {a:?} expected = {e:?}");
                    }
                }
                (Some(a), None) => {
                    println!("UNVETTED for k = {k} actual = {a:?}");
                }
                (None, Some(e)) => println!("ERROR for k = {k} expected = {e:?}"),
                (None, None) => {}
            }
        }

        Ok(StageOutputs { actual })
    }
}
