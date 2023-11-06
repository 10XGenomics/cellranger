//! Martian stage CALL_TAGS_OH
//! Assign cells to samples using their overhang.
use crate::count_matrix::CountMatrixFile;
use crate::read_level_multiplexing::{
    get_barcodes_per_multiplexing_identifier, get_umi_per_multiplexing_identifier_for_feature_type,
};
use anyhow::Result;
use barcode::whitelist::BarcodeId;
use barcode::WhitelistSource;
use cr_types::chemistry::ChemistryDef;
use cr_types::{FeatureType, MetricsFile};
use json_report_derive::JsonReport;
use martian::prelude::{MartianMain, MartianRover};
use martian_derive::{make_mro, MartianStruct};
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::FileTypeWrite;
use metric::{Metric, TxHashMap};
use serde::{Deserialize, Serialize};

/// The Martian stage inputs.
#[derive(Clone, Deserialize, MartianStruct)]
pub struct CallTagsOHStageInputs {
    pub chemistry_def: ChemistryDef,
    pub raw_feature_bc_matrix: CountMatrixFile,
}

/// The Martian stage outputs.
#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct CallTagsOHStageOutputs {
    pub barcodes_per_tag: Option<JsonFile<TxHashMap<BarcodeId, Vec<String>>>>,
    pub summary: Option<MetricsFile>,
}

/// Martian stage CALL_TAGS_OH
/// Assign valid bcs to samples using their probe barcode sequence.
pub struct CallTagsOH;

/// Overhang metrics.
#[derive(JsonReport)]
struct CallTagsOHMetrics {
    /// The number of UMI per overhang.
    #[json_report(block)]
    umi_per_overhang: TxHashMap<BarcodeId, i64>,

    /// The number of valid barcodes per overhang
    #[json_report(block)]
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
        if !args.chemistry_def.name.is_overhang_multiplexed() {
            return Ok(Self::StageOutputs {
                barcodes_per_tag: None,
                summary: None,
            });
        }
        let overhang_read_component = args.chemistry_def.overhang_read_barcode().unwrap();
        let overhang_range = overhang_read_component.offset()
            ..overhang_read_component.offset() + overhang_read_component.length();

        let overhang_seq_to_id =
            WhitelistSource::from_spec(overhang_read_component.whitelist(), true, None)?
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

        let metrics = CallTagsOHMetrics::new(umi_per_overhang, barcodes_per_overhang);

        Ok(Self::StageOutputs {
            barcodes_per_tag: Some(barcodes_per_tag_file),
            summary: Some(MetricsFile::from_reporter(&rover, "summary", &metrics)?),
        })
    }
}
