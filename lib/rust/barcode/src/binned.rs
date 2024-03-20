//!
//! Barcodes which are defined based on binning geometry
//!

use anyhow::ensure;
use pyo3::{pyclass, pymethods};
use serde::{Deserialize, Serialize};
use std::fmt::{Display, Formatter};
use std::str::FromStr;

/// A row or column index in a square binning scheme
#[derive(Serialize, Deserialize, Clone, Copy, PartialOrd, Ord, PartialEq, Eq, Hash, Debug)]
pub struct SquareBinRowOrColumnIndex {
    pub index: usize,
    pub size_um: u32,
}

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize, Debug)]
#[pyclass]
pub struct SquareBinIndex {
    #[pyo3(get)]
    pub row: usize,
    #[pyo3(get)]
    pub col: usize,
    #[pyo3(get)]
    pub size_um: u32,
}

impl SquareBinIndex {
    pub fn scale(self, pitch_um: u32) -> u32 {
        self.size_um / pitch_um
    }
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, anyhow::Error> {
        Self::from_str(std::str::from_utf8(bytes)?)
    }

    /// Extracts binned barcode from unbinned HD barcode
    pub fn extract_binned_barcode(self, bin_scale: i32) -> SquareBinIndex {
        SquareBinIndex {
            row: self.row / bin_scale as usize,
            col: self.col / bin_scale as usize,
            size_um: self.size_um * bin_scale as u32,
        }
    }
}

#[pymethods]
impl SquareBinIndex {
    #[new]
    pub fn new(
        barcode: Option<String>,
        row: Option<usize>,
        col: Option<usize>,
        size_um: Option<u32>,
    ) -> pyanyhow::Result<Self> {
        match (barcode, row, col, size_um) {
            (Some(barcode), None, None, None) => Ok(barcode.parse()?),
            (None, Some(row), Some(col), Some(size_um)) => Ok(SquareBinIndex { row, col, size_um }),
            bad_args => Err(anyhow::anyhow!("Invalid arguments: {:?}", bad_args).into()),
        }
    }
    fn __str__(&self) -> String {
        self.to_string()
    }
}

pub const SQUARE_BIN_PREFIX: &str = "s";
const MICROMETER: &str = "um";

impl Display for SquareBinIndex {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{SQUARE_BIN_PREFIX}_{:03}{MICROMETER}_{:05}_{:05}",
            self.size_um, self.row, self.col
        )
    }
}

impl FromStr for SquareBinIndex {
    type Err = anyhow::Error;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        fn error(s: &str) -> anyhow::Error {
            anyhow::anyhow!("Unable to parse {s} as a spatial barcode.")
        }
        let s = match s.split_once('-') {
            Some((barcode, _gem_group)) => barcode,
            None => s,
        };
        let mut parts = s.split('_');
        let prefix = parts.next().ok_or_else(|| error(s))?;
        ensure!(
            prefix == SQUARE_BIN_PREFIX,
            "Invalid prefix {} in {}",
            prefix,
            s
        );

        let size = parts
            .next()
            .and_then(|p| p.strip_suffix(MICROMETER))
            .map(str::parse)
            .transpose()?
            .ok_or_else(|| error(s))?;
        let row = parts
            .next()
            .map(str::parse)
            .transpose()?
            .ok_or_else(|| error(s))?;
        let col = parts
            .next()
            .map(str::parse)
            .transpose()?
            .ok_or_else(|| error(s))?;

        Ok(SquareBinIndex {
            row,
            col,
            size_um: size,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_square_bin() {
        let b = SquareBinIndex {
            row: 1,
            col: 2,
            size_um: 16,
        };
        assert_eq!(b.to_string(), "s_016um_00001_00002");
        assert_eq!(b, "s_016um_00001_00002".parse().unwrap());
        assert_eq!(b, "s_016um_00001_00002-1".parse().unwrap());
    }

    proptest::proptest! {
        #[test]
        fn test_roundtrip(row in 0..5000usize, col in 0..5000usize, size_um in 0..200u32) {
            let barcode = SquareBinIndex { row, col, size_um };
            assert_eq!(barcode, barcode.to_string().parse().unwrap());
        }
    }
}
