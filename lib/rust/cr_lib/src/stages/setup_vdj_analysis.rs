//! Martians stage SETUP_VDJ_ANALYSIS
#![deny(missing_docs)]

use super::setup_vdj_demux::VdjDemuxSampleInfo;
use crate::GexMatrices;
use anyhow::Result;
use cr_h5::count_matrix::CountMatrixFile;
use cr_types::reference::feature_reference::BeamMode;
use cr_types::{BarcodeMultiplexingType, ReadLevel};
use itertools::Itertools;
use martian::prelude::*;
use martian_derive::{MartianStruct, make_mro};
use martian_filetypes::FileTypeRead;
use martian_filetypes::tabular_file::CsvFile;
use serde::{Deserialize, Serialize};
use vdj_reference::VdjReceptor;

#[derive(Serialize, Deserialize, Clone, MartianStruct, Debug)]
pub struct VdjAnalysisConfig {
    per_sample: bool,
    is_multi: bool,
    has_no_vdj_ref: bool,
    denovo: bool,
    has_antigen: bool,
    skip_clonotyping: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct SetupVdjAnalysisStageInputs {
    pub receptor: Option<VdjReceptor>,
    pub vdj_config: VdjAnalysisConfig,
    pub demux_sample_info: Option<VdjDemuxSampleInfo>,
    pub lib_level_gex: GexMatrices,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct SetupVdjAnalysisStageOutputs {
    pub receptor: VdjReceptor,
    pub disable_cell_calling: bool,
    pub disable_clonotyping: bool,
    pub disable_beam: bool,
    pub disable_asm_metrics: bool,
    pub multiplexing_method: Option<BarcodeMultiplexingType>,
    pub beam_mode: Option<BeamMode>,
    pub filtered_matrix_h5: Option<CountMatrixFile>,
    pub raw_matrix_h5: Option<CountMatrixFile>,
    pub filtered_barcodes: Option<CsvFile<()>>,
}

/// Martians stage SETUP_VDJ_ANALYSIS
pub struct SetupVdjAnalysis;

#[derive(Eq, PartialEq)]
pub enum VdjAnalysisType {
    /// Cellranger vdj
    LibraryLevelCount,

    /// Cellranger multi library-level VDJ analysis
    LibraryLevelMulti,

    /// Cellranger multi sample-level VDJ analysis
    SampleLevelMulti,
}

struct VdjBools {
    disable_cell_calling: bool,
    disable_clonotyping: bool,
}

impl VdjAnalysisType {
    fn default_bools(&self) -> VdjBools {
        #[allow(clippy::enum_glob_use)]
        use VdjAnalysisType::*;
        match self {
            LibraryLevelCount => VdjBools {
                disable_cell_calling: false,
                disable_clonotyping: false,
            },
            LibraryLevelMulti => VdjBools {
                disable_cell_calling: true,
                disable_clonotyping: true,
            },
            SampleLevelMulti => VdjBools {
                disable_cell_calling: false,
                disable_clonotyping: false,
            },
        }
    }
}

#[make_mro]
impl MartianMain for SetupVdjAnalysis {
    type StageInputs = SetupVdjAnalysisStageInputs;
    type StageOutputs = SetupVdjAnalysisStageOutputs;
    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        let analysis_type = match (args.vdj_config.per_sample, args.vdj_config.is_multi) {
            (false, false) => VdjAnalysisType::LibraryLevelCount,
            (false, true) => VdjAnalysisType::LibraryLevelMulti,
            (true, true) => VdjAnalysisType::SampleLevelMulti,
            (true, false) => unreachable!(),
        };

        let mut vdj_bools = analysis_type.default_bools();
        vdj_bools.disable_clonotyping |=
            args.vdj_config.has_no_vdj_ref || args.vdj_config.skip_clonotyping;

        let multiplexing_type = if let Some(fingerprint) = args
            .demux_sample_info
            .as_ref()
            .and_then(|info| info.fingerprint.as_ref())
        {
            fingerprint
                .read()?
                .into_iter()
                .map(|f| f.barcode_multiplexing_type())
                .dedup()
                .exactly_one()
                .unwrap()
        } else {
            None
        };

        let disable_asm_metrics = analysis_type == VdjAnalysisType::SampleLevelMulti
            && matches!(
                multiplexing_type,
                Some(BarcodeMultiplexingType::CellLevel(_))
            );

        let sample_level_files =
            if let Some(sample) = args.demux_sample_info.and_then(|info| info.gex_matrices) {
                // Hardlink to play nice with Martian
                sample.hardlink(&rover)?
            } else {
                GexMatrices::default()
            };

        let gex_matrices = match (analysis_type, multiplexing_type) {
            // celranger vdj
            (VdjAnalysisType::LibraryLevelCount, None) => GexMatrices::default(),
            // celranger vdj is incompatible with multiplexing
            (VdjAnalysisType::LibraryLevelCount, Some(_)) => unreachable!(),
            // cellranger multi library-level analysis
            (VdjAnalysisType::LibraryLevelMulti, _) => args.lib_level_gex.hardlink(&rover)?, // Hardlink to play nice with Martian
            // cellranger multi VDJ+GEX sample-level analysis
            (VdjAnalysisType::SampleLevelMulti, None) => sample_level_files,
            // cellranger multi CMO and Hashtag multiplexed sample-level analysis
            (VdjAnalysisType::SampleLevelMulti, Some(BarcodeMultiplexingType::CellLevel(_))) => {
                sample_level_files
            }
            (
                VdjAnalysisType::SampleLevelMulti,
                Some(BarcodeMultiplexingType::ReadLevel(ReadLevel::OH)),
            ) => sample_level_files,
            (
                VdjAnalysisType::SampleLevelMulti,
                Some(BarcodeMultiplexingType::ReadLevel(ReadLevel::RTL)),
            ) => {
                unreachable!()
            }
        };

        let receptor = match (args.receptor, args.vdj_config.denovo) {
            (None, true) => VdjReceptor::TR, // This is the default value of vdj receptor in denovo mode
            (None, false) => unreachable!(),
            (Some(r), _) => r,
        };

        let beam_mode = match (args.vdj_config.has_antigen, receptor) {
            (true, VdjReceptor::TR | VdjReceptor::TRGD) => Some(BeamMode::BeamT),
            (true, VdjReceptor::IG) => Some(BeamMode::BeamAB),
            (false, _) => None,
        };

        let disable_beam = !args.vdj_config.has_antigen
            || vdj_bools.disable_clonotyping
            || !args.vdj_config.per_sample;

        Ok(SetupVdjAnalysisStageOutputs {
            receptor,
            disable_cell_calling: vdj_bools.disable_cell_calling,
            disable_clonotyping: vdj_bools.disable_clonotyping,
            disable_beam,
            beam_mode,
            filtered_matrix_h5: gex_matrices.filtered_matrix_h5,
            raw_matrix_h5: gex_matrices.raw_matrix_h5,
            filtered_barcodes: gex_matrices.filtered_barcodes,
            disable_asm_metrics,
            multiplexing_method: multiplexing_type,
        })
    }
}
