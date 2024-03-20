use crate::reference::feature_reference::FeatureReference;
use anyhow::{bail, Result};
use std::collections::HashMap;
use std::iter::zip;

/// Given an array of integer counts, one per feature.
/// Return the proportions, normalized within each feature type.
pub fn compute_feature_dist(raw_counts: Vec<i64>, feat_ref: &FeatureReference) -> Result<Vec<f64>> {
    let feature_types: Vec<_> = feat_ref
        .feature_defs
        .iter()
        .map(|fd| fd.feature_type)
        .collect();

    if raw_counts.len() != feature_types.len() {
        bail!(
            "Mismatch between number of features in counts JSON ({}) and in feature reference ({})",
            raw_counts.len(),
            feature_types.len()
        );
    }

    let mut feature_type_sums = HashMap::new();
    for (&count, &feature_type) in zip(&raw_counts, &feature_types) {
        *feature_type_sums.entry(feature_type).or_insert(0) += count;
    }

    let mut proportions = vec![0f64; raw_counts.len()];
    for (i, (count, feature_type)) in zip(raw_counts, feature_types).enumerate() {
        let sum = feature_type_sums[&feature_type];
        if sum > 0 {
            proportions[i] = count as f64 / sum as f64;
        }
    }

    // when all feature barcode features are zero, which can happen for e.g. cycle failure,
    //   initialize to one count equivalent each. TODO: implement laplacian smoothing
    let fbc_feature_indices = feat_ref
        .feature_defs
        .iter()
        .map(|fd| fd.index)
        .collect::<Vec<_>>();
    if fbc_feature_indices.iter().all(|&i| proportions[i] == 0.0) {
        for &i in &fbc_feature_indices {
            proportions[i] = 1.0 / fbc_feature_indices.len() as f64;
        }
    }

    Ok(proportions)
}
