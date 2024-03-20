use anyhow::Result;
use barcode::whitelist::BarcodeId;
use barcode::BcSegSeq;
use cr_h5::count_matrix::{BarcodeWithGemGroup, CountMatrix};
use cr_types::reference::feature_reference::FeatureType;
use itertools::Itertools;
use metric::TxHashMap;
use std::collections::HashMap;
use std::ops::Range;

/// Return the ID of this read-level multiplexing sequence.
fn map_multiplexing_seq_to_id(
    barcode: &BarcodeWithGemGroup,
    seq_to_id_map: &TxHashMap<BcSegSeq, BarcodeId>,
    multiplexing_seq_range: &Range<usize>,
) -> BarcodeId {
    let seq = &BcSegSeq::from_bytes(&barcode.as_bytes()[multiplexing_seq_range.clone()]);
    seq_to_id_map[seq]
}

/// Return the valid barcodes for each read-level multiplexing identifier
pub fn get_barcodes_per_multiplexing_identifier(
    matrix: &CountMatrix,
    seq_to_id_map: &TxHashMap<BcSegSeq, BarcodeId>,
    multiplexing_seq_range: &Range<usize>,
) -> Result<TxHashMap<BarcodeId, Vec<String>>> {
    Ok(matrix
        .barcodes()
        .iter()
        .fold(TxHashMap::default(), |mut acc, barcode| {
            let multiplexing_id =
                map_multiplexing_seq_to_id(barcode, seq_to_id_map, multiplexing_seq_range);
            if let Some(barcodes) = acc.get_mut(&multiplexing_id) {
                barcodes.push(barcode.to_string());
            } else {
                acc.insert(multiplexing_id, vec![barcode.to_string()]);
            }
            acc
        }))
}

/// Return the number of UMIs for each read-level multiplexing identifier
/// for all feature types.
pub fn get_umi_per_multiplexing_identifier(
    matrix: &CountMatrix,
    seq_to_id_map: &TxHashMap<BcSegSeq, BarcodeId>,
    multiplexing_seq_range: &Range<usize>,
) -> HashMap<FeatureType, TxHashMap<BarcodeId, i64>> {
    matrix
        .counts()
        .group_by(|x| (x.barcode, x.feature.feature_type))
        .into_iter()
        .map(|((barcode, feature_type), counts)| {
            (
                (
                    feature_type,
                    map_multiplexing_seq_to_id(barcode, seq_to_id_map, multiplexing_seq_range),
                ),
                counts.map(|x| x.count as i64).sum(),
            )
        })
        .into_grouping_map()
        .sum()
        .into_iter()
        .map(|((feature_type, multiplexing_id), count)| (feature_type, (multiplexing_id, count)))
        .into_grouping_map()
        .collect()
}

/// Return the number of UMIs for each read-level multiplexing identifier
/// for the specified feature type.
pub fn get_umi_per_multiplexing_identifier_for_feature_type(
    matrix: &CountMatrix,
    seq_to_id_map: &TxHashMap<BcSegSeq, BarcodeId>,
    multiplexing_seq_range: &Range<usize>,
    feature_type: FeatureType,
) -> TxHashMap<BarcodeId, i64> {
    matrix.barcode_counts_for_feature_type(feature_type).fold(
        TxHashMap::default(),
        |mut acc, (barcode, count)| {
            let multiplexing_id =
                map_multiplexing_seq_to_id(barcode, seq_to_id_map, multiplexing_seq_range);
            *acc.entry(multiplexing_id).or_default() += count;
            acc
        },
    )
}
