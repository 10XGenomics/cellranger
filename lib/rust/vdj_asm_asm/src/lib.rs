//! vdj_asm_asm
#![deny(missing_docs)]
pub mod airrfilter;
pub mod asm_call_cells;
pub mod asm_metrics;
mod assembly;
pub use assembly::{Assembly, BarcodeDataBriefFile};
mod adapter;
pub mod assembly_types;
mod contig_aligner;
pub mod filter_exact_clonotypes;
pub mod filter_sample_specific;
pub mod make_exact_clonotypes;
pub mod make_filter_switch;
pub mod merge_per_sample_annotations;
pub mod subset_assembly_outs;
mod translator;
pub mod write_ann_csv;
pub mod write_contig_outs;

// initialize insta test harness
#[cfg(test)]
#[ctor::ctor]
fn init() {
    // this ensures insta knows where to find its snap tests
    let cwd = std::env::current_dir().unwrap();
    let workspace_root = cwd.parent().unwrap();
    unsafe { std::env::set_var("INSTA_WORKSPACE_ROOT", workspace_root) }
}
