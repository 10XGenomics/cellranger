//! Detection of probe barcode pairing between gene expression and feature barcoding.
use crate::barcode_overlap::{
    calculate_frp_gem_barcode_overlap, FRPGemBarcodeOverlapRow, GelBeadBarcodesPerProbeBarcode,
    ProbeBarcodeGelBeadGrouper,
};
use anyhow::Result;
use barcode::whitelist::{categorize_multiplexing_barcode_id, BarcodeId, MultiplexingBarcodeType};
use barcode::{BarcodeConstruct, BcSegSeq, WhitelistSource};
use cr_types::chemistry::ChemistryDef;
use fastq_set::read_pair::{ReadPair, ReadPart};
use itertools::Itertools;
use multi::config::MultiConfigCsv;
use std::collections::HashSet;

/// Return true if we should detect probe barcode pairings.
pub fn should_detect_probe_barcode_pairing(multi_config: &MultiConfigCsv) -> bool {
    multi_config.is_rtl_multiplexed()
        && multi_config.libraries.has_gene_expression()
        && multi_config.libraries.has_feature_barcode()
}

/// Only consider probe barcodes that were seen in at least this fraction of all gel beads.
const GEL_BEAD_FRAC_THRESHOLD: f64 = 0.005;

pub type BestPairs = Vec<(BarcodeId, BarcodeId)>;

/// Identify the optimal pairing of Gex and Ab probe barcodes.
/// Ignore any probe barcode that was only observed in a trivial fraction of all GEMs.
/// Return the overlap matrix and the best pairings.
pub fn detect_probe_barcode_pairing(
    chemistry_def: &ChemistryDef,
    reads: &[Vec<ReadPair>],
) -> Result<(Vec<FRPGemBarcodeOverlapRow>, BestPairs)> {
    let whitelist_source =
        WhitelistSource::construct(chemistry_def.barcode_whitelist(), false, None)?;
    let probe_barcode_seq_to_id = whitelist_source.as_ref().probe().as_raw_seq_to_id()?;
    let whitelist = whitelist_source
        .as_ref()
        .map_result(WhitelistSource::as_whitelist)?;
    let bp_range = chemistry_def.barcode_range();
    let mut gel_bead_barcodes_per_probe_barcode = ProbeBarcodeGelBeadGrouper::group_all(
        reads
            .iter()
            .flatten()
            .filter_map(|read_pair| bp_range.map_option(|r| read_pair.get_range(r, ReadPart::Seq)))
            .map(|seqs| seqs.map(BcSegSeq::from_bytes))
            .filter_map(move |seqs| {
                seqs.zip(whitelist.as_ref())
                    .map_option(|(seq, wl)| wl.match_to_whitelist(seq))
            })
            .map(|barcode_components| match barcode_components {
                BarcodeConstruct::GelBeadAndProbe(comps) => comps,
                other => panic!("expected gel bead and probe construct but found {other:?}"),
            }),
        &probe_barcode_seq_to_id,
    );
    filter_gel_bead_barcodes_per_probe_barcode(&mut gel_bead_barcodes_per_probe_barcode);

    let overlaps = calculate_frp_gem_barcode_overlap(&gel_bead_barcodes_per_probe_barcode);
    let pairings = calculate_matching(
        overlaps
            .iter()
            .filter_map(|row| Some((get_rtl_and_ab_barcode_from_row(row)?, row.overlap))),
    );
    Ok((overlaps, pairings))
}

/// Filter the mapping from probe barcode to gel bead barcodes.
/// Ignore all probe barcodes seen in only a tiny fraction of GEMs.
fn filter_gel_bead_barcodes_per_probe_barcode(mapping: &mut GelBeadBarcodesPerProbeBarcode) {
    let total_gel_beads = mapping
        .values()
        .flatten()
        .map(|(gb_bc, _count)| *gb_bc)
        .unique()
        .count();
    mapping.retain(|_, gel_bead_bcs| {
        (gel_bead_bcs.len() as f64 / total_gel_beads as f64) > GEL_BEAD_FRAC_THRESHOLD
    });
}

/// Return Some((rtl_barcode, ab_barcode)) if one barcode is an RTL barcode
/// and the other is an AB barcode and None otherwise.
pub fn get_rtl_and_ab_barcode_from_row(
    row: &FRPGemBarcodeOverlapRow,
) -> Option<(BarcodeId, BarcodeId)> {
    match (
        categorize_multiplexing_barcode_id(&row.barcode1_id),
        categorize_multiplexing_barcode_id(&row.barcode2_id),
    ) {
        (MultiplexingBarcodeType::RTL, MultiplexingBarcodeType::Antibody) => {
            Some((row.barcode1_id, row.barcode2_id))
        }
        (MultiplexingBarcodeType::Antibody, MultiplexingBarcodeType::RTL) => {
            Some((row.barcode2_id, row.barcode1_id))
        }
        _ => None,
    }
}

/// Pair each RTL multiplexing barcode with an antibody multiplexing barcode.
/// Calculate a maximum weight matching of a bipartite graph using a greedy heuristic.
/// See https://en.wikipedia.org/wiki/Assignment_problem
/// and https://en.wikipedia.org/wiki/Maximum_weight_matching
fn calculate_matching(overlaps: impl Iterator<Item = ((BarcodeId, BarcodeId), f64)>) -> BestPairs {
    let mut matched = HashSet::new();
    overlaps
        .sorted_by(|a, b| {
            a.1.partial_cmp(&b.1)
                .unwrap()
                .reverse()
                .then_with(|| a.0.cmp(&b.0))
        })
        .filter_map(|(barcode_pair, _weight)| {
            if matched.contains(&barcode_pair.0) || matched.contains(&barcode_pair.1) {
                None
            } else {
                matched.insert(barcode_pair.0);
                matched.insert(barcode_pair.1);
                Some(barcode_pair)
            }
        })
        .sorted()
        .collect()
}

#[cfg(test)]
mod test {
    use super::calculate_matching;
    use barcode::whitelist::BarcodeId;

    #[test]
    fn test_calculate_matching() {
        let overlaps = [
            ("BC004", "AB004", 0.440),
            ("BC004", "AB005", 0.034),
            ("BC004", "AB006", 0.063),
            ("BC004", "AB007", 0.049),
            ("BC004", "AB008", 0.039),
            ("BC005", "AB004", 0.197),
            ("BC005", "AB005", 0.425),
            ("BC005", "AB006", 0.055),
            ("BC005", "AB007", 0.044),
            ("BC005", "AB008", 0.049),
            ("BC006", "AB004", 0.189),
            ("BC006", "AB005", 0.102),
            ("BC006", "AB006", 0.557),
            ("BC006", "AB007", 0.103),
            ("BC006", "AB008", 0.103),
            ("BC007", "AB004", 0.205),
            ("BC007", "AB005", 0.035),
            ("BC007", "AB006", 0.008),
            ("BC007", "AB007", 0.066),
            ("BC007", "AB008", 0.006),
            ("BC008", "AB004", 0.198),
            ("BC008", "AB005", 0.027),
            ("BC008", "AB006", 0.024),
            ("BC008", "AB007", 0.018),
            ("BC008", "AB008", 0.149),
        ]
        .map(|x| ((BarcodeId::pack(x.0), BarcodeId::pack(x.1)), x.2));

        let result = [
            ("BC004", "AB004"),
            ("BC005", "AB005"),
            ("BC006", "AB006"),
            ("BC007", "AB007"),
            ("BC008", "AB008"),
        ]
        .map(|x| (BarcodeId::pack(x.0), BarcodeId::pack(x.1)));

        assert_eq!(calculate_matching(overlaps.into_iter()), result);
    }
}
