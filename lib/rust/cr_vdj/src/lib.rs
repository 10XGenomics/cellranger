//! cr_vdj
#![deny(missing_docs)]

pub mod clonotype_hist;
pub mod clonotype_table;
/// Copy the VDJ reference
pub mod copy_vdj_reference;
pub mod create_barcode_csv;
pub mod make_vdj_plots;
mod matrix;
pub mod summarize_vdj_filters;

use serde::Serialize;

/// Assert that two JSON values are equal
pub fn check_eq_json(j1: &str, j2: &str) {
    pretty_assertions::assert_eq!(
        serde_json::from_str::<serde_json::value::Value>(j1).unwrap(),
        serde_json::from_str::<serde_json::value::Value>(j2).unwrap()
    );
}

/// Test JSON round trip
pub fn test_json_roundtrip<T: Serialize + serde::de::DeserializeOwned>(json: &str) -> T {
    let parsed: T = serde_json::from_str(json).unwrap();
    let parsed_str = serde_json::to_string(&parsed).unwrap();
    check_eq_json(&parsed_str, json);
    parsed
}
