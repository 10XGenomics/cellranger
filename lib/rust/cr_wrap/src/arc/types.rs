use anyhow::{bail, ensure, Result};
use clap::{self, Parser};
use ordered_float::NotNan;
use serde::Serialize;
use std::collections::HashMap;

#[derive(Parser, Debug, Clone)]
pub struct ForceCellsArgs {
    /// Cell caller override: define the minimum number of ATAC transposition events
    /// in peaks (ATAC counts) for a cell barcode. Note: this option must be specified in
    /// conjunction with `min-gex-count`. With `--min-atac-count=X` and
    /// `--min-gex-count=Y` a barcode is defined as a cell if it contains at least
    /// X ATAC counts AND at least Y GEX UMI counts.
    #[clap(long, value_name = "NUM")]
    pub min_atac_count: Option<usize>,

    /// Cell caller override: define the minimum number of GEX UMI counts for a cell
    /// barcode. Note: this option must be specified in conjunction with
    /// `min-atac-count`. With `--min-atac-count=X` and
    /// `--min-gex-count=Y` a barcode is defined as a cell if it contains at least
    /// X ATAC counts AND at least Y GEX UMI counts.
    #[clap(long, value_name = "NUM")]
    pub min_gex_count: Option<usize>,
}

impl ForceCellsArgs {
    pub fn to_mro_arg(&self) -> Result<Option<HashMap<String, MinCounts>>> {
        match (self.min_atac_count, self.min_gex_count) {
            (None, None) => Ok(None),
            (Some(a), Some(g)) => Ok(Some(
                [("default".to_string(), MinCounts { atac: a, gex: g })]
                    .iter()
                    .cloned()
                    .collect(),
            )),
            (Some(v), None) => bail!(
                "Invalid cell caller override: 'min_atac_count' = {v} but 'min_gex_count' \
                 unspecified. You must specify both to override the cell caller."
            ),
            (None, Some(v)) => bail!(
                "Invalid cell caller override: 'min_gex_count' = {v} but 'min_atac_count' \
                 unspecified. You must specify both to override the cell caller."
            ),
        }
    }
}

#[derive(Debug, Serialize, Clone)]
pub struct MinCounts {
    atac: usize,
    gex: usize,
}

/// Feature linkage max distance argument validation
///
/// Must be positive
pub fn validate_distance(d: &str) -> Result<NotNan<f64>> {
    let w: NotNan<f64> = d.parse()?;
    ensure!(w.into_inner() > 0.0, "value = {w} must be positive");
    Ok(w)
}

/// Peak Q-value argument validation
///
/// Must be strictly between 0 and 1
pub fn validate_strict_fraction(q: &str) -> Result<NotNan<f64>> {
    let w: NotNan<f64> = q.parse()?;
    let x = w.into_inner();
    ensure!(0.0 < x && x < 1.0, "value = {w} must satisfy 0 < value < 1");
    Ok(w)
}

/// K-means max clusters argument must satisfy 2 <= K <= 10.
pub const MAX_CLUSTERS_RANGE: std::ops::RangeInclusive<u64> = 2..=10;
