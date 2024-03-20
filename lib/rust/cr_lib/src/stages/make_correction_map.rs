//! Martian stage MAKE_CORRECTION_MAP

use anyhow::Result;
use barcode::BarcodeConstruct;
use barcode_extensions::{CorrectionMap, CorrectionMapBuilder};
use cr_types::chemistry::{BarcodeExtraction, ChemistryDef, ChemistryDefs, ChemistryDefsExt};
use cr_types::{BcSegmentCountFormat, LibraryType};
use itertools::Itertools;
use martian::{MartianRover, MartianStage, MartianVoid, Resource, StageDef};
use martian_derive::{make_mro, martian_filetype, MartianStruct};
use martian_filetypes::bin_file::BinaryFormat;
use martian_filetypes::{FileTypeRead, FileTypeWrite};
use metric::TxHashMap;
use serde::{Deserialize, Serialize};

martian_filetype!(CorrectionMapFile, "cmf");
pub type CorrectionMapFormat =
    BinaryFormat<CorrectionMapFile, TxHashMap<LibraryType, BarcodeConstruct<CorrectionMap>>>;

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct MakeCorrectionMapStageInputs {
    chemistry_defs: ChemistryDefs,

    /// Counts of uncorrected valid segments per library type
    pub barcode_segment_counts: BcSegmentCountFormat,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct MakeCorrectionMapStageOutputs {
    correction_map: Option<CorrectionMapFormat>,
}

pub struct MakeCorrectionMap;

fn is_two_part_joint_extraction(chemistry_def: &ChemistryDef) -> bool {
    if !matches!(
        chemistry_def.barcode_extraction(),
        Some(BarcodeExtraction::JointBc1Bc2 { .. })
    ) {
        return false;
    }

    match chemistry_def.barcode_construct() {
        BarcodeConstruct::Segmented(segments) if segments.is_two_part() => {
            let slide_names: Vec<_> = segments
                .into_iter()
                .filter_map(|s| s.whitelist().slide_name())
                .unique()
                .collect();
            assert!(
                slide_names.len() <= 1,
                "Barcode parts have different slides: {slide_names:?}"
            );
            !slide_names.is_empty()
        }
        _ => false,
    }
}

#[make_mro(volatile = strict)]
impl MartianStage for MakeCorrectionMap {
    type StageInputs = MakeCorrectionMapStageInputs;
    type StageOutputs = MakeCorrectionMapStageOutputs;
    type ChunkInputs = MartianVoid;
    type ChunkOutputs = MartianVoid;

    fn split(
        &self,
        args: Self::StageInputs,
        _rover: MartianRover,
    ) -> Result<StageDef<Self::ChunkInputs>> {
        let chemistry_def = args.chemistry_defs.primary();
        Ok(if is_two_part_joint_extraction(chemistry_def) {
            // ~6GB memory is needed to build all the sequences two edits away
            // from a whitelist of size 5k with 19 bases each in memory.
            // Once deduped, it reduces to ~2.5GB. Since we produce the sequences
            // one part at a time, the required memory is ~2.5+6 = 8.5GB.
            StageDef::with_join_resource(Resource::new().mem_gb(12).threads(4))
        } else {
            StageDef::new()
        })
    }

    fn main(
        &self,
        _args: Self::StageInputs,
        _chunk_args: Self::ChunkInputs,
        _rover: MartianRover,
    ) -> Result<Self::ChunkOutputs> {
        unreachable!()
    }

    fn join(
        &self,
        args: Self::StageInputs,
        _chunk_defs: Vec<Self::ChunkInputs>,
        _chunk_outs: Vec<Self::ChunkOutputs>,
        rover: MartianRover,
    ) -> Result<Self::StageOutputs> {
        let chemistry_def = &args.chemistry_defs.primary();
        if !is_two_part_joint_extraction(chemistry_def) {
            return Ok(MakeCorrectionMapStageOutputs {
                correction_map: None,
            });
        }

        rayon::ThreadPoolBuilder::new()
            .num_threads(rover.get_threads())
            .build_global()
            .unwrap();

        let correction_map_builders = chemistry_def.barcode_construct().map_result(|bc| {
            bc.whitelist()
                .as_source(false)
                .and_then(|src| src.as_whitelist().map(|wl| CorrectionMapBuilder::new(&wl)))
        })?;

        let correction_map: CorrectionMapFormat = rover.make_path("correction_map");

        let per_lib_segment_fst: TxHashMap<_, _> = args
            .barcode_segment_counts
            .read()?
            .into_iter()
            .map(|(lib_type, segment_counts)| {
                correction_map_builders
                    .as_ref()
                    .zip(segment_counts)
                    .map_result(|(builder, counts)| builder.build(&counts))
                    .map(|data| (lib_type, data))
            })
            .try_collect()?;
        #[allow(clippy::drop_non_drop)]
        drop(correction_map_builders);

        correction_map.write(&per_lib_segment_fst)?;

        Ok(MakeCorrectionMapStageOutputs {
            correction_map: Some(correction_map),
        })
    }
}
