// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
#![deny(missing_docs)]

// Miscellaneous functions.

use crate::core::defs::BarcodeContigs;
use std::cmp::Ordering;

pub(crate) fn sort_tig_bc(tig_bc: &mut [BarcodeContigs]) {
    tig_bc.sort_by(|x, y| -> Ordering {
        for i in 0..x.len() {
            // Order by number of chains.

            if i >= y.len() {
                return Ordering::Greater;
            }

            let ordering = x[i].cmp(&y[i]);
            if ordering != Ordering::Equal {
                return ordering;
            }
        }
        if x.len() < y.len() {
            return Ordering::Less;
        }
        Ordering::Equal
    });
}
