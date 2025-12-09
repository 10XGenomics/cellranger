//! cr_types::reference
#![deny(missing_docs)]

pub mod feature_checker;
pub mod feature_extraction;
pub mod feature_reference;
pub mod genome_of_chrom;
pub mod probe_set_reference;
/// Metadata of the reference from the `reference.json` (if reference_path
/// provided) or from the probe set reference CSV file.
pub mod reference_info;
