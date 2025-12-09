//! Compute overlap coefficient matrix for a collection of barcodes.
#![deny(missing_docs)]
use barcode::whitelist::BarcodeId;
use barcode::{BcSegSeq, GelBeadAndProbeConstruct};
use itertools::Itertools;
use metric::TxHashMap;
use serde::{Deserialize, Serialize};
use std::cmp::min;

/// Associate a probe barcode ID to the set of all gel bead barcodes it was observed in.
/// Tracks the number of observations of the probe barcode in the gel bead barcode.
pub type GelBeadBarcodesPerProbeBarcode = TxHashMap<BarcodeId, TxHashMap<BcSegSeq, usize>>;

pub struct ProbeBarcodeGelBeadGrouper<'a> {
    probe_barcode_seq_to_id: &'a TxHashMap<BcSegSeq, BarcodeId>,
    barcodes: GelBeadBarcodesPerProbeBarcode,
}

impl<'a> ProbeBarcodeGelBeadGrouper<'a> {
    /// Group all gel bead barcodes by the probe barcodes they are associated with.
    /// Each construct is treated as a single observation.
    pub fn group_all(
        constructs: impl Iterator<Item = GelBeadAndProbeConstruct<BcSegSeq>>,
        probe_barcode_seq_to_id: &'a TxHashMap<BcSegSeq, BarcodeId>,
    ) -> GelBeadBarcodesPerProbeBarcode {
        let mut grouper = Self::new(probe_barcode_seq_to_id);
        for construct in constructs {
            grouper.group(construct, 1);
        }
        grouper.finish()
    }

    /// Initialize a new grouper.
    pub fn new(probe_barcode_seq_to_id: &'a TxHashMap<BcSegSeq, BarcodeId>) -> Self {
        Self {
            probe_barcode_seq_to_id,
            barcodes: GelBeadBarcodesPerProbeBarcode::default(),
        }
    }

    /// Group a single construct with the provided observation count.
    pub fn group(&mut self, construct: GelBeadAndProbeConstruct<BcSegSeq>, count: usize) {
        let probe_barcode_id = self
            .probe_barcode_seq_to_id
            .get(&construct.probe)
            .unwrap_or_else(|| panic!("probe barcode not found in whitelist: {}", construct.probe));
        *self
            .barcodes
            .entry(*probe_barcode_id)
            .or_default()
            .entry(construct.gel_bead)
            .or_default() += count;
    }

    /// Consume this grouper and return the barcode grouping.
    pub fn finish(self) -> GelBeadBarcodesPerProbeBarcode {
        self.barcodes
    }
}

/// Calculate the number of gel-bead barcodes in common for each pair of probe barcodes.
pub fn calculate_barcode_overlap_counts(
    gel_bead_barcodes_per_probe_barcode: &GelBeadBarcodesPerProbeBarcode,
) -> TxHashMap<(BarcodeId, BarcodeId), i64> {
    gel_bead_barcodes_per_probe_barcode
        .iter()
        .sorted_by_key(|&(probe_barcode_id, _)| probe_barcode_id)
        .tuple_combinations()
        .map(|(x, y)| {
            let overlap_count = x.1.keys().filter(|k| y.1.contains_key(*k)).count();
            ((*x.0, *y.0), overlap_count as i64)
        })
        .collect()
}

/// One row of frp_gem_barcode_overlap.csv.
#[derive(Deserialize, Serialize)]
pub struct FRPGemBarcodeOverlapRow {
    /// First probe barcode identifier
    pub barcode1_id: BarcodeId,
    /// Second probe barcode identifier
    pub barcode2_id: BarcodeId,
    /// Number of gel-bead barcodes of barcode1
    pub barcode1_gems: i64,
    /// Number of gel-bead barcodes of barcode2
    pub barcode2_gems: i64,
    /// Number of gel-bead barcodes in common
    pub common_gems: i64,
    /// The overlap coefficient of these two probe barcodes
    /// overlap(X, Y) = |X & Y| / min(|X|, |Y|)
    pub overlap: f64,
}

impl FRPGemBarcodeOverlapRow {
    /// Swap positions of barcode1 and barcode2.
    pub fn swap_order(&mut self) {
        (
            self.barcode1_id,
            self.barcode2_id,
            self.barcode1_gems,
            self.barcode2_gems,
        ) = (
            self.barcode2_id,
            self.barcode1_id,
            self.barcode2_gems,
            self.barcode1_gems,
        );
    }
}

/// Calculate overlap coefficients of all probe barcodes.
/// The results are sorted by barcode IDs.
pub fn calculate_frp_gem_barcode_overlap(
    gel_bead_barcodes_per_probe_barcode: &GelBeadBarcodesPerProbeBarcode,
) -> Vec<FRPGemBarcodeOverlapRow> {
    calculate_barcode_overlap_counts(gel_bead_barcodes_per_probe_barcode)
        .into_iter()
        .sorted_by_key(|&(barcode1_barcode2, _)| barcode1_barcode2)
        .map(|((barcode1_id, barcode2_id), common_gems)| {
            let barcode1_gems = gel_bead_barcodes_per_probe_barcode[&barcode1_id].len() as i64;
            let barcode2_gems = gel_bead_barcodes_per_probe_barcode[&barcode2_id].len() as i64;
            let overlap = common_gems as f64 / min(barcode1_gems, barcode2_gems) as f64;
            FRPGemBarcodeOverlapRow {
                barcode1_id,
                barcode2_id,
                barcode1_gems,
                barcode2_gems,
                common_gems,
                overlap,
            }
        })
        .collect()
}
