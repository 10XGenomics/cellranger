//! enclone
// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.
#![deny(missing_docs)]

use crate::core::barcode_fate::BarcodeFate;
use crate::core::enclone_structs::BarcodeFates;
use itertools::zip_eq;
use vector_utils::erase_if_iter;

pub(crate) mod allele;
pub(crate) mod graph_filter;
pub(crate) mod info;
pub(crate) mod innate;
pub(crate) mod join;
pub(crate) mod join2;
pub(crate) mod join_core;
pub(crate) mod misc1;
pub(crate) mod misc2;
pub(crate) mod misc3;

/// Filter barcodes and record their fates.
pub(crate) trait BarcodeFilter<T> {
    /// Filter the collection of items.
    /// For any item that the filter removes, return the reason the barcode was
    /// removed.
    // TODO: if all filters apply a single fate, we could extract this into a
    // trait constant.
    fn filter(&self, items: &[T]) -> Vec<Option<BarcodeFate>>;

    /// Get the dataset ID and barcodes from this item to record fate.
    fn fate_keys(&self, item: &T) -> impl Iterator<Item = (usize, String)>;

    fn filter_items(&self, items: &mut Vec<T>, fate: &mut [BarcodeFates]) {
        let filter_result = self.filter(items);
        for (item, filter_reason) in zip_eq(items.iter_mut(), &filter_result) {
            let Some(reason) = filter_reason else {
                continue;
            };
            for (dataset_id, barcode) in self.fate_keys(item) {
                fate[dataset_id].insert(barcode.to_string(), *reason);
            }
        }
        erase_if_iter(items, filter_result.iter().map(Option::is_some));
    }
}
