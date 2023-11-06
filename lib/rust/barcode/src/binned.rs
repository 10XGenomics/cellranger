//!
//! Barcodes which are defined based on binning geometry
//!

use serde::{Deserialize, Serialize};
use std::fmt::{Display, Formatter};

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct SpatialSquareBinBarcode {
    pub row: usize,
    pub col: usize,
    pub size_um: u32,
    pub gem_group: u16,
}

impl SpatialSquareBinBarcode {
    pub fn scale(self, pitch_um: u32) -> u32 {
        self.size_um / pitch_um
    }
}

impl Display for SpatialSquareBinBarcode {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        // TODO: Update based on the agreed convention, say "square_{size}um_{row}_{col}-{}"
        // We want to keep the length under crate::MAX_BARCODE_LENGTH
        write!(
            f,
            // ,
            "s_{:03}um_{:05}_{:05}-{}",
            self.size_um, self.row, self.col, self.gem_group
        )
    }
}
