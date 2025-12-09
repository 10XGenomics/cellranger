//! clonotype_assigner
#![deny(missing_docs)]
pub mod assigner;
pub use assigner::Assigner;
pub mod fill_clonotype_info;
pub mod handle_no_clonotyping;
pub mod write_clonotype_outs;
pub mod write_concat_ref_outs;
pub mod write_consensus_bam;
pub mod write_consensus_txt;
