//! cr_ana
#![deny(missing_docs)]
use cr_types::FeatureBarcodeType;
use cr_types::reference::feature_reference::FeatureType;

mod hclust_utils;
mod io;
mod louvain;
mod pca;
#[cfg(test)]
mod stage_testing;
pub mod stages;
#[cfg(test)]
mod test_pipeline;
mod types;

pub(crate) const EXCLUDED_FEATURE_TYPES: &[FeatureType] =
    &[FeatureType::Barcode(FeatureBarcodeType::Antigen)];
