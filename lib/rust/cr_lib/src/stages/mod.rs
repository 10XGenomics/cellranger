//! cr_lib::stages
#![deny(missing_docs)]

mod align_and_count;
mod barcode_correction;
mod build_per_sample_vdj_ws_contents;
mod call_tags_genetic;
mod call_tags_overhang;
mod call_tags_rtl;
mod check_barcodes_compatibility;
mod check_barcodes_compatibility_vdj;
mod check_single_beam_mode;
mod collate_metrics;
mod collate_probe_metrics;
mod compute_antigen_vdj_metrics;
mod copy_chemistry_spec;
mod count_alleles;
mod create_multi_graph;
mod demux_probe_bc_matrix;
mod detect_chemistry;
mod detect_chemistry_test;
mod detect_vdj_receptor;
mod expect_single_barcode_whitelist;
mod extract_single_chemistry;
mod generate_cas_websummary;
mod get_chemistry_def;
mod get_gdna_metrics;
mod logic_not;
mod make_correction_map;
mod make_shard;
mod merge_gem_well_files;
mod merge_metrics;
mod multi_preflight;
mod multi_setup_chunks;
mod parse_multi_config;
mod pick_beam_analyzer;
mod rust_bridge;
mod setup_reference_info;
mod setup_vdj_analysis;
mod setup_vdj_demux;
mod write_barcode_index;
mod write_barcode_summary;
mod write_gene_index;
mod write_h5_matrix;
mod write_matrix_market;
mod write_minimap_index;
mod write_molecule_info;
mod write_multi_web_summary;
mod write_pos_bam;

#[cfg(feature = "tenx_internal")]
mod internal;
#[cfg(feature = "tenx_source_available")]
mod stubs;

pub use align_and_count::{
    AlignAndCount, DEFAULT_TRANSCRIPTOME_MIN_SCORE, MAX_READ_ANN_SAMPLE, ReadShards,
};
pub use barcode_correction::{BarcodeCorrection, correct_barcode_in_read};
pub use build_per_sample_vdj_ws_contents::BuildPerSampleVdjWsContents;
pub use call_tags_genetic::CallTagsGenetic;
pub use call_tags_overhang::CallTagsOH;
pub use call_tags_rtl::CallTagsRTL;
pub use check_barcodes_compatibility::CheckBarcodesCompatibility;
pub(super) use check_barcodes_compatibility::sample_valid_barcodes;
pub use check_barcodes_compatibility_vdj::CheckBarcodesCompatibilityVdj;
pub use check_single_beam_mode::CheckSingleBeamMode;
pub use collate_metrics::CollateMetrics;
pub(super) use collate_metrics::SampleMetrics;
pub use collate_probe_metrics::CollateProbeMetrics;
pub use compute_antigen_vdj_metrics::ComputeAntigenVdjMetrics;
pub use copy_chemistry_spec::CopyChemistrySpec;
pub use count_alleles::CountAlleles;
pub use create_multi_graph::CreateMultiGraph;
pub use demux_probe_bc_matrix::DemuxProbeBcMatrix;
pub use detect_chemistry::DetectChemistry;
pub(super) use detect_chemistry::DetectedProbeBarcodePairingFile;
pub use detect_chemistry_test::DetectChemistryTest;
pub use detect_vdj_receptor::DetectVdjReceptor;
pub use expect_single_barcode_whitelist::ExpectSingleBarcodeWhitelist;
pub use extract_single_chemistry::ExtractSingleChemistry;
pub use generate_cas_websummary::GenerateCasWebsummary;
pub use get_chemistry_def::GetChemistryDef;
pub use get_gdna_metrics::GetGdnaMetrics;
pub use logic_not::LogicNot;
pub use make_correction_map::MakeCorrectionMap;
pub use make_shard::MakeShard;
pub use merge_gem_well_files::MergeGemWellFiles;
pub use merge_metrics::MergeMetrics;
pub use multi_preflight::MultiPreflight;
pub use multi_setup_chunks::MultiSetupChunks;
pub use parse_multi_config::{CellCalling, CountInputs, ParseMultiConfig};
pub use pick_beam_analyzer::PickBeamAnalyzer;
pub use rust_bridge::RustBridge;
pub use setup_reference_info::SetupReferenceInfo;
pub use setup_vdj_analysis::SetupVdjAnalysis;
pub use setup_vdj_demux::SetupVDJDemux;
pub use write_barcode_index::WriteBarcodeIndex;
pub use write_barcode_summary::WriteBarcodeSummary;
pub use write_gene_index::WriteGeneIndex;
pub use write_h5_matrix::WriteH5Matrix;
pub use write_matrix_market::WriteMatrixMarket;
pub use write_minimap_index::WriteMinimapIndex;
pub use write_molecule_info::WriteMoleculeInfo;
pub use write_multi_web_summary::{WriteMultiWebSummary, load_dist, write_web_summary_html};
pub use write_pos_bam::{WritePosBam, make_bam_header};

/// The StageInputs of a MartianStage.
pub type StageInputs<T> = <T as martian::MartianStage>::StageInputs;
