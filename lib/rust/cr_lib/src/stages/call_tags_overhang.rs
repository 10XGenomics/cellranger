//! Martian stage CALL_TAGS_OH
//! Assign cells to samples using their overhang.
use crate::read_level_multiplexing::{
    get_barcodes_per_multiplexing_identifier, get_umi_per_multiplexing_identifier_for_feature_type,
};
use anyhow::Result;
use barcode::whitelist::BarcodeId;
use cr_h5::count_matrix::CountMatrixFile;
use cr_types::chemistry::ChemistryDefs;
use cr_types::reference::feature_reference::FeatureType;
use itertools::Itertools;
use martian::prelude::{MartianMain, MartianRover};
use martian_derive::{make_mro, MartianStruct};
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::FileTypeWrite;
use metric::TxHashMap;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// The Martian stage inputs.
#[derive(Clone, Deserialize, MartianStruct)]
pub struct CallTagsOHStageInputs {
    pub chemistry_defs: ChemistryDefs,
    pub raw_feature_bc_matrix: CountMatrixFile,
}

/// The Martian stage outputs.
#[derive(Serialize, Deserialize, MartianStruct)]
pub struct CallTagsOHStageOutputs {
    pub barcodes_per_tag: Option<JsonFile<TxHashMap<BarcodeId, Vec<String>>>>,
    pub summary: Option<JsonFile<CallTagsOHMetrics>>,
}

/// Martian stage CALL_TAGS_OH
/// Assign valid bcs to samples using their probe barcode sequence.
pub struct CallTagsOH;

/// Overhang metrics.
#[derive(Serialize)]
pub struct CallTagsOHMetrics {
    /// The number of UMI per overhang.
    umi_per_overhang: TxHashMap<BarcodeId, i64>,

    /// The number of valid barcodes per overhang
    valid_barcodes_per_overhang: TxHashMap<BarcodeId, i64>,
}

impl CallTagsOHMetrics {
    /// Calculate the overhang metrics.
    fn new(
        umi_per_overhang: TxHashMap<BarcodeId, i64>,
        barcodes_per_overhang: TxHashMap<BarcodeId, Vec<String>>,
    ) -> Self {
        let valid_barcodes_per_overhang = barcodes_per_overhang
            .iter()
            .map(|(overhang_id, gel_bead_barcodes)| (*overhang_id, gel_bead_barcodes.len() as i64))
            .collect();

        Self {
            umi_per_overhang,
            valid_barcodes_per_overhang,
        }
    }
}

#[make_mro(mem_gb = 8, volatile = strict)]
impl MartianMain for CallTagsOH {
    type StageInputs = CallTagsOHStageInputs;
    type StageOutputs = CallTagsOHStageOutputs;

    /// Run the Martian stage CALL_TAGS_OH
    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        let overhang_barcodes: HashMap<_, _> = args
            .chemistry_defs
            .iter()
            .filter_map(|(library_type, chemistry_def)| {
                Some((library_type, chemistry_def.overhang_read_barcode()?))
            })
            .collect();

        if overhang_barcodes.is_empty() {
            return Ok(Self::StageOutputs {
                barcodes_per_tag: None,
                summary: None,
            });
        }

        let (overhang_offset, overhang_length) = overhang_barcodes
            .values()
            .map(|x| (x.offset(), x.length()))
            .dedup()
            .exactly_one()
            .unwrap();

        let overhang_range = overhang_offset..(overhang_offset + overhang_length);

        let whitelist_sources: Vec<_> = overhang_barcodes
            .into_values()
            .map(|x| x.whitelist().as_source(true))
            .try_collect()?;

        // TODO(CELLRANGER-7847): Factor out duplicated seq_to_id code
        let overhang_seq_to_id: TxHashMap<_, _> = whitelist_sources
            .into_iter()
            .dedup()
            .exactly_one()
            .unwrap()
            .as_translation_seq_to_id()?;

        let matrix = args.raw_feature_bc_matrix.read()?;
        let barcodes_per_overhang = get_barcodes_per_multiplexing_identifier(
            &matrix,
            &overhang_seq_to_id,
            &overhang_range,
        )?;
        let barcodes_per_tag_file: JsonFile<_> = rover.make_path("barcodes_per_tag");
        barcodes_per_tag_file.write(&barcodes_per_overhang)?;

        let umi_per_overhang = get_umi_per_multiplexing_identifier_for_feature_type(
            &matrix,
            &overhang_seq_to_id,
            &overhang_range,
            FeatureType::Gene,
        );

        let summary: JsonFile<_> = rover.make_path("summary");
        summary.write(&CallTagsOHMetrics::new(
            umi_per_overhang,
            barcodes_per_overhang,
        ))?;

        Ok(Self::StageOutputs {
            barcodes_per_tag: Some(barcodes_per_tag_file),
            summary: Some(summary),
        })
    }
}
