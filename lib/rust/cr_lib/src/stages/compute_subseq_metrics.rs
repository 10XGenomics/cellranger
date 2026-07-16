//! Martian stage COMPUTE_SUBSEQ_METRICS
//!
//! Compute distinct-subsequence metrics over cell-associated gel-bead barcodes.
#![deny(missing_docs)]

use anyhow::{Result, bail};
use barcode::{BcSegSeq, WhitelistSpec};
use cr_types::LibraryType;
use cr_types::chemistry::{BarcodeReadComponent, ChemistryDef, ChemistryDefs, ChemistryDefsExt};
use cr_types::filtered_barcodes::{FilteredBarcodesCsv, read_filtered_barcodes_set};
use martian::{MartianMain, MartianRover};
use martian_derive::{MartianStruct, make_mro};
use martian_filetypes::FileTypeWrite;
use martian_filetypes::json_file::JsonFile;
use metric::{TxHashMap, TxHashSet};
use serde::{Deserialize, Serialize};
use std::ops::Range;

#[derive(Clone, Debug)]
struct SubseqRegion {
    /// Region index (0 or 1) used in metric names
    index: usize,
    /// Half-open byte range within the gel-bead barcode sequence.
    range: Range<usize>,
}

/// Subseq regions for the gel-bead barcode.
///
/// Splits the barcode into equal-length half subsequences when length >= 2.
/// Otherwise keeps a single region spanning the full-length barcode.
fn subseq_regions(gel_bead: &BarcodeReadComponent) -> Vec<SubseqRegion> {
    let length = gel_bead.length();
    if length < 2 {
        vec![SubseqRegion {
            index: 0,
            range: 0..length,
        }]
    } else {
        vec![
            SubseqRegion {
                index: 0,
                range: 0..length / 2,
            },
            SubseqRegion {
                index: 1,
                range: length / 2..length,
            },
        ]
    }
}

/// Stage inputs.
#[derive(Clone, Deserialize, MartianStruct)]
pub struct ComputeSubseqMetricsStageInputs {
    /// Cell-associated barcodes from FILTER_BARCODES.
    pub filtered_barcodes: FilteredBarcodesCsv,
    /// Chemistry definitions, used for whitelist + OCM-ID positions.
    pub chemistry_defs: ChemistryDefs,
}

/// Stage outputs.
#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct ComputeSubseqMetricsStageOutputs {
    /// Flat metrics JSON serialized from `SubseqMetrics`.
    pub summary: JsonFile<SubseqMetrics>,
}

/// Per-region distinct-subseq metrics emitted to the summary JSON.
/// Field names exactly match the metric keys.
/// Region 0 is always populated. Region 1 is None when length < 2.
#[derive(Clone, Serialize)]
pub struct SubseqMetrics {
    pub barcode_subseq_0_tot_count: u64,
    pub barcode_subseq_0_obs_count: u64,
    pub barcode_subseq_0_exp_count: f64,
    pub barcode_subseq_0_logp: f64,
    pub barcode_subseq_0_ge_0_count: u64,
    pub barcode_subseq_0_ge_1_count: u64,
    pub barcode_subseq_0_ge_2_count: u64,
    pub barcode_subseq_0_ge_3_count: u64,
    pub barcode_subseq_0_ge_4_count: u64,
    pub barcode_subseq_0_ge_5_count: u64,
    pub barcode_subseq_1_tot_count: Option<u64>,
    pub barcode_subseq_1_obs_count: Option<u64>,
    pub barcode_subseq_1_exp_count: Option<f64>,
    pub barcode_subseq_1_logp: Option<f64>,
    pub barcode_subseq_1_ge_0_count: Option<u64>,
    pub barcode_subseq_1_ge_1_count: Option<u64>,
    pub barcode_subseq_1_ge_2_count: Option<u64>,
    pub barcode_subseq_1_ge_3_count: Option<u64>,
    pub barcode_subseq_1_ge_4_count: Option<u64>,
    pub barcode_subseq_1_ge_5_count: Option<u64>,
    pub barcode_subseq_sum_logp: f64,
}

/// Martian stage COMPUTE_SUBSEQ_METRICS.
pub struct ComputeSubseqMetrics;

#[make_mro(mem_gb = 2, volatile = strict)]
impl MartianMain for ComputeSubseqMetrics {
    type StageInputs = ComputeSubseqMetricsStageInputs;
    type StageOutputs = ComputeSubseqMetricsStageOutputs;

    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        let summary: JsonFile<SubseqMetrics> = rover.make_path("summary");
        let chem = select_chemistry_def(&args.chemistry_defs)?;
        let metrics = match gel_bead_component(chem) {
            // If Chemistry lacks a gel-bead barcode, set `SubseqMetrics` values to zeroes.
            None => SubseqMetrics::zeros(),
            Some(gel_bead) => {
                let ocm_id_rel_range =
                    ocm_id_relative_range(gel_bead, args.chemistry_defs.overhang_read_barcode());
                let regions = subseq_regions(gel_bead);
                validate_regions(&regions, gel_bead.length())?;
                let stats = compute_metrics(
                    gel_bead,
                    ocm_id_rel_range.as_ref(),
                    &regions,
                    &args.filtered_barcodes,
                )?;
                stats.into_metrics(&regions)
            }
        };
        summary.write(&metrics)?;
        Ok(ComputeSubseqMetricsStageOutputs { summary })
    }
}

impl SubseqMetrics {
    /// All-zero `SubseqMetrics` (used when the chemistry has no unique
    /// gel-bead component, to pass a stable key set with zeroes to MERGE_METRICS).
    fn zeros() -> Self {
        SubseqMetrics {
            barcode_subseq_0_tot_count: 0,
            barcode_subseq_0_obs_count: 0,
            barcode_subseq_0_exp_count: 0.0,
            barcode_subseq_0_logp: 0.0,
            barcode_subseq_0_ge_0_count: 0,
            barcode_subseq_0_ge_1_count: 0,
            barcode_subseq_0_ge_2_count: 0,
            barcode_subseq_0_ge_3_count: 0,
            barcode_subseq_0_ge_4_count: 0,
            barcode_subseq_0_ge_5_count: 0,
            barcode_subseq_1_tot_count: Some(0),
            barcode_subseq_1_obs_count: Some(0),
            barcode_subseq_1_exp_count: Some(0.0),
            barcode_subseq_1_logp: Some(0.0),
            barcode_subseq_1_ge_0_count: Some(0),
            barcode_subseq_1_ge_1_count: Some(0),
            barcode_subseq_1_ge_2_count: Some(0),
            barcode_subseq_1_ge_3_count: Some(0),
            barcode_subseq_1_ge_4_count: Some(0),
            barcode_subseq_1_ge_5_count: Some(0),
            barcode_subseq_sum_logp: 0.0,
        }
    }
}

/// Select the chemistry definition to use for subseq metrics.
///
/// Returns the GEX chemistry definition when present. Otherwise, checks that
/// all non-null entries have gel-bead barcodes with identical whitelists,
/// and returns the first entry.
fn select_chemistry_def(defs: &ChemistryDefs) -> Result<&ChemistryDef> {
    if let Some(gex) = defs.get(&LibraryType::Gex) {
        return Ok(gex);
    }
    let mut reference_whitelist = None;
    for chem in defs.values() {
        let Some(gel_bead) = chem.barcode_construct().option_gel_bead() else {
            bail!("COMPUTE_SUBSEQ_METRICS: all chemistry definitions must have gel-bead barcodes");
        };
        let wl = gel_bead.whitelist();
        match &reference_whitelist {
            None => reference_whitelist = Some(wl.clone()),
            Some(ref_wl) => {
                if wl != ref_wl {
                    bail!(
                        "COMPUTE_SUBSEQ_METRICS: chemistry definitions have \
                         inconsistent gel-bead whitelists"
                    );
                }
            }
        }
    }
    defs.values()
        .next()
        .ok_or_else(|| anyhow::anyhow!("COMPUTE_SUBSEQ_METRICS: no chemistry definitions provided"))
}

/// Return the gel-bead barcode component for a chemistry that exposes one, else None.
fn gel_bead_component(chem: &ChemistryDef) -> Option<&BarcodeReadComponent> {
    chem.barcode_construct().option_gel_bead()
}

/// If the OCM-ID barcode is within the gel-bead barcode, return its
/// position relative to the start of the gel-bead sequence (else `None`).
fn ocm_id_relative_range(
    gel_bead: &BarcodeReadComponent,
    ocm_id: Option<&BarcodeReadComponent>,
) -> Option<Range<usize>> {
    let ocm_id = ocm_id?;
    if ocm_id.rprange().read() != gel_bead.rprange().read() {
        return None;
    }
    let gb_start = gel_bead.offset();
    let gb_end = gb_start + gel_bead.length();
    let ocm_start = ocm_id.offset();
    let ocm_end = ocm_start + ocm_id.length();
    if ocm_start < gb_start || ocm_end > gb_end {
        return None;
    }
    Some((ocm_start - gb_start)..(ocm_end - gb_start))
}

/// Validate that the configured regions have unique 0/1 indices, are
/// inside the gel-bead barcode, are pairwise non-overlapping, and
/// jointly span the full-length gel-bead barcode without gaps.
fn validate_regions(regions: &[SubseqRegion], gel_bead_len: usize) -> Result<()> {
    if regions.is_empty() {
        bail!("no subseq regions configured for gel-bead length {gel_bead_len}");
    }
    let mut seen = TxHashSet::default();
    for r in regions {
        if r.range.end > gel_bead_len || r.range.start >= r.range.end {
            bail!(
                "subseq region {} ({:?}) does not fit gel-bead length {}",
                r.index,
                r.range,
                gel_bead_len
            );
        }
        if r.index > 1 {
            bail!(
                "subseq region index {} not supported (must be 0 or 1)",
                r.index,
            );
        }
        if !seen.insert(r.index) {
            bail!("duplicate subseq region index {}", r.index);
        }
    }
    let mut sorted: Vec<_> = regions.iter().collect();
    sorted.sort_by_key(|r| r.range.start);
    if sorted[0].range.start != 0 {
        bail!(
            "subseq regions do not start at 0 (first start: {})",
            sorted[0].range.start
        );
    }
    if sorted.last().unwrap().range.end != gel_bead_len {
        bail!(
            "subseq regions do not span gel-bead length {} (last end: {})",
            gel_bead_len,
            sorted.last().unwrap().range.end
        );
    }
    for pair in sorted.windows(2) {
        if pair[0].range.end != pair[1].range.start {
            bail!(
                "subseq regions are not contiguous: {:?} -> {:?}",
                pair[0].range,
                pair[1].range
            );
        }
    }
    Ok(())
}

/// Maximum (inclusive) bin index for the per-subseq cell-count histogram.
const MAX_SUBSEQ_HIST_VAL: usize = 5;

/// Aggregated per-region statistics for the metric output.
struct PerRegionStats {
    tot_count: u64,
    obs_count: u64,
    exp_count: f64,
    logp: f64,
    /// `ge_count[x]` = number of whitelist subseqs observed in at least `x`
    /// cells (cell count `>= x`), for `x` in `0..=MAX_SUBSEQ_HIST_VAL`.
    ge_count: [u64; MAX_SUBSEQ_HIST_VAL + 1],
}

/// All region stats produced by a single stage invocation.
struct StageStats {
    per_region: Vec<PerRegionStats>,
    sum_logp: f64,
}

impl StageStats {
    /// Map per-region stats onto the explicit `SubseqMetrics` field names.
    fn into_metrics(self, regions: &[SubseqRegion]) -> SubseqMetrics {
        let mut by_index: [Option<PerRegionStats>; 2] = [None, None];
        for (region, stats) in regions.iter().zip(self.per_region) {
            by_index[region.index] = Some(stats);
        }
        let r0 = by_index[0]
            .take()
            .expect("validate_regions guarantees region index 0 is present");
        let r1 = by_index[1].take();
        SubseqMetrics {
            barcode_subseq_0_tot_count: r0.tot_count,
            barcode_subseq_0_obs_count: r0.obs_count,
            barcode_subseq_0_exp_count: r0.exp_count,
            barcode_subseq_0_logp: r0.logp,
            barcode_subseq_0_ge_0_count: r0.ge_count[0],
            barcode_subseq_0_ge_1_count: r0.ge_count[1],
            barcode_subseq_0_ge_2_count: r0.ge_count[2],
            barcode_subseq_0_ge_3_count: r0.ge_count[3],
            barcode_subseq_0_ge_4_count: r0.ge_count[4],
            barcode_subseq_0_ge_5_count: r0.ge_count[5],
            barcode_subseq_1_tot_count: r1.as_ref().map(|r| r.tot_count),
            barcode_subseq_1_obs_count: r1.as_ref().map(|r| r.obs_count),
            barcode_subseq_1_exp_count: r1.as_ref().map(|r| r.exp_count),
            barcode_subseq_1_logp: r1.as_ref().map(|r| r.logp),
            barcode_subseq_1_ge_0_count: r1.as_ref().map(|r| r.ge_count[0]),
            barcode_subseq_1_ge_1_count: r1.as_ref().map(|r| r.ge_count[1]),
            barcode_subseq_1_ge_2_count: r1.as_ref().map(|r| r.ge_count[2]),
            barcode_subseq_1_ge_3_count: r1.as_ref().map(|r| r.ge_count[3]),
            barcode_subseq_1_ge_4_count: r1.as_ref().map(|r| r.ge_count[4]),
            barcode_subseq_1_ge_5_count: r1.as_ref().map(|r| r.ge_count[5]),
            barcode_subseq_sum_logp: self.sum_logp,
        }
    }
}

/// Per-OCM-ID bucket key. An empty key represents the single-bucket fallback.
type OcmIdKey = BcSegSeq;
/// Per-OCM-ID count maps for one region: subseq -> count of whitelist barcodes.
type RegionWlCounts = TxHashMap<OcmIdKey, TxHashMap<BcSegSeq, u64>>;

/// Numerical floor on `log10(p)` to avoid `-inf` for vanishingly small p-values.
const LOGP_FLOOR: f64 = -300.0;

/// Stream the gel-bead whitelist and the cell-associated barcodes once each,
/// then compute the per-region metrics.
fn compute_metrics(
    gel_bead: &BarcodeReadComponent,
    ocm_id_rel_range: Option<&Range<usize>>,
    regions: &[SubseqRegion],
    filtered_barcodes: &FilteredBarcodesCsv,
) -> Result<StageStats> {
    // Whitelist pass: a single iteration over the whitelist file. We retain
    // only per-region/per-OCM-ID count maps and per-OCM-ID totals; the whitelist
    // itself is never materialised in memory.
    let wl_spec = gel_bead.whitelist();
    let use_translation = match wl_spec {
        WhitelistSpec::TxtFile { translation, .. } => *translation,
        WhitelistSpec::DynamicTranslation { .. } => true,
        WhitelistSpec::SlideFile { .. } => false,
    };
    let wl_source = wl_spec.as_source()?;
    let mut wl_total: TxHashMap<OcmIdKey, u64> = TxHashMap::default();
    let mut wl_counts: Vec<RegionWlCounts> = (0..regions.len())
        .map(|_| RegionWlCounts::default())
        .collect();

    for entry in wl_source.iter()? {
        let (seq, trans, _id) = entry?;
        let canonical = if use_translation {
            trans.ok_or_else(|| {
                anyhow::anyhow!(
                    "translation whitelist entry is missing the translated sequence column"
                )
            })?
        } else {
            seq
        };
        let bytes = canonical.seq();
        if bytes.len() != gel_bead.length() {
            bail!(
                "whitelist row has unexpected length {}, expected gel-bead \
                 length {} (whitelist + chemistry def out of sync)",
                bytes.len(),
                gel_bead.length(),
            );
        }
        let ocm_id = ocm_id_key(bytes, ocm_id_rel_range);
        *wl_total.entry(ocm_id).or_default() += 1;
        for (region, counts) in regions.iter().zip(wl_counts.iter_mut()) {
            let slice = BcSegSeq::from_bytes(&bytes[region.range.clone()]);
            *counts.entry(ocm_id).or_default().entry(slice).or_default() += 1;
        }
    }

    // Cell pass: read the filtered barcodes set, group by OCM-ID, and count
    // how many cells have each (region, OCM-ID) slice value.
    let cells = read_filtered_barcodes_set(filtered_barcodes)?;
    let mut cells_per_ocm_id: TxHashMap<OcmIdKey, u64> = TxHashMap::default();
    let mut observed: Vec<TxHashMap<OcmIdKey, TxHashMap<BcSegSeq, u64>>> =
        (0..regions.len()).map(|_| TxHashMap::default()).collect();

    for bc in &cells {
        let bytes = bc.sequence_bytes();
        if bytes.len() < gel_bead.length() {
            bail!(
                "cell barcode {} has unexpected length {}, expected at least \
                 gel-bead length {} (filtered barcodes + chemistry def out of sync)",
                bc,
                bytes.len(),
                gel_bead.length(),
            );
        }
        // Truncate to gel-bead portion of full barcode sequence
        let bytes = &bytes[..gel_bead.length()];
        let ocm_id = ocm_id_key(bytes, ocm_id_rel_range);
        *cells_per_ocm_id.entry(ocm_id).or_default() += 1;
        for (region, per_ocm) in regions.iter().zip(observed.iter_mut()) {
            let slice = BcSegSeq::from_bytes(&bytes[region.range.clone()]);
            *per_ocm.entry(ocm_id).or_default().entry(slice).or_default() += 1;
        }
    }

    // Per-region statistics: Poisson CDF (per OCM-ID, summed) plus the
    // observation-count histogram bins over the whitelist universe.
    let per_region: Vec<PerRegionStats> = regions
        .iter()
        .enumerate()
        .map(|(i, _)| {
            // Universe = distinct subseqs in the gel-bead whitelist for this
            // region, across all OCM-ID buckets.
            let universe: TxHashSet<BcSegSeq> = wl_counts[i]
                .values()
                .flat_map(|inner| inner.keys().copied())
                .collect();
            let tot_count = universe.len() as u64;

            // Observation count per slice = sum across OCM-ID buckets.
            let mut cell_count_per_slice: TxHashMap<BcSegSeq, u64> = TxHashMap::default();
            for inner in observed[i].values() {
                for (slice, &cnt) in inner {
                    *cell_count_per_slice.entry(*slice).or_default() += cnt;
                }
            }

            // Histogram bins over the possible whitelist subseqs.
            // `ge_count[x]` = number of subseqs with observed cell-count `>= x`.
            let mut ge_count = [0u64; MAX_SUBSEQ_HIST_VAL + 1];
            for slice in &universe {
                let cnt = cell_count_per_slice.get(slice).copied().unwrap_or(0);
                for (x, slot) in ge_count.iter_mut().enumerate() {
                    if cnt >= x as u64 {
                        *slot += 1;
                    }
                }
            }

            // Existing per-(region, OCM-ID) Poisson statistics.
            let mut obs_count = 0u64;
            let mut exp_count = 0.0;
            let mut logp = 0.0;
            for (ocm_id, &cells_j) in &cells_per_ocm_id {
                let w_j = wl_total.get(ocm_id).copied().unwrap_or(0);
                if cells_j == 0 || w_j == 0 {
                    continue;
                }
                let m_rj = expected_distinct_for_ocm_id(wl_counts[i].get(ocm_id), cells_j, w_j);
                let k_rj = observed[i].get(ocm_id).map_or(0, |m| m.len() as u64);
                obs_count += k_rj;
                exp_count += m_rj;
                logp += poisson_log10_cdf(m_rj, k_rj);
            }
            PerRegionStats {
                tot_count,
                obs_count,
                exp_count,
                logp: logp.max(LOGP_FLOOR),
                ge_count,
            }
        })
        .collect();

    let sum_logp = per_region.iter().map(|s| s.logp).sum::<f64>();
    Ok(StageStats {
        per_region,
        sum_logp,
    })
}

/// Extract the OCM-ID key (a slice of the gel-bead barcode) or the empty key
/// when the chemistry has no OCM-ID inside the gel-bead.
fn ocm_id_key(bytes: &[u8], ocm_id_rel_range: Option<&Range<usize>>) -> OcmIdKey {
    match ocm_id_rel_range {
        Some(r) => BcSegSeq::from_bytes(&bytes[r.clone()]),
        None => BcSegSeq::default(),
    }
}

/// Expected number of distinct subseqs observed for one (region, OCM-ID),
/// `m_{r,j} = sum over s in S_{r,j} of (1 - exp(-cells_j * c_{s,j,r} / W_j))`.
fn expected_distinct_for_ocm_id(
    wl_counts_for_ocm_id: Option<&TxHashMap<BcSegSeq, u64>>,
    cells_j: u64,
    w_j: u64,
) -> f64 {
    let Some(counts) = wl_counts_for_ocm_id else {
        return 0.0;
    };
    let scale = cells_j as f64 / w_j as f64;
    counts
        .values()
        .map(|&c| 1.0 - (-scale * c as f64).exp())
        .sum()
}

/// `log10` of the exact Poisson CDF, `P(X <= k | X ~ Poisson(m)) =
/// sum_{i=0}^{k} exp(-m) * m^i / i!`, computed with streaming log-sum-exp so
/// that values much smaller than `f64::MIN_POSITIVE` are still represented in
/// log space. Returns 0.0 (i.e. log10(1)) when `m == 0`.
fn poisson_log10_cdf(m: f64, k: u64) -> f64 {
    if m <= 0.0 {
        return 0.0;
    }
    let ln_m = m.ln();
    // Term i = 0: log(exp(-m)) = -m.
    let mut log_term = -m;
    let mut max_log = log_term;
    let mut sum_norm = 1.0; // = exp(log_term - max_log) for i = 0.
    for i in 1..=k {
        // log_term_i = log_term_{i-1} + ln(m) - ln(i).
        log_term += ln_m - (i as f64).ln();
        if log_term > max_log {
            sum_norm *= (max_log - log_term).exp();
            max_log = log_term;
        }
        sum_norm += (log_term - max_log).exp();
    }
    let log_cdf_natural = max_log + sum_norm.ln();
    (log_cdf_natural / std::f64::consts::LN_10).max(LOGP_FLOOR)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn r(start: usize, end: usize) -> Range<usize> {
        start..end
    }

    fn region(i: usize, range: Range<usize>) -> SubseqRegion {
        SubseqRegion { index: i, range }
    }

    #[test]
    fn validate_regions_accepts_partition() {
        assert!(validate_regions(&[region(0, r(0, 8)), region(1, r(8, 16))], 16).is_ok());
        assert!(validate_regions(&[region(0, r(0, 16))], 16).is_ok());
    }

    #[test]
    fn validate_regions_rejects_empty_set() {
        assert!(validate_regions(&[], 16).is_err());
    }

    #[test]
    fn validate_regions_rejects_out_of_bounds() {
        assert!(validate_regions(&[region(0, r(0, 17))], 16).is_err());
    }

    #[test]
    fn validate_regions_rejects_empty_or_inverted() {
        assert!(validate_regions(&[region(0, r(5, 5))], 16).is_err());
    }

    #[test]
    fn validate_regions_rejects_duplicate_index() {
        // Duplicate indices are caught even if the ranges happen to partition.
        assert!(validate_regions(&[region(0, r(0, 8)), region(0, r(8, 16))], 16).is_err());
    }

    #[test]
    fn validate_regions_rejects_gap() {
        assert!(validate_regions(&[region(0, r(0, 4)), region(1, r(8, 16))], 16).is_err());
    }

    #[test]
    fn validate_regions_rejects_overlap() {
        assert!(validate_regions(&[region(0, r(0, 10)), region(1, r(8, 16))], 16).is_err());
    }

    #[test]
    fn validate_regions_rejects_non_zero_start() {
        assert!(validate_regions(&[region(0, r(2, 16))], 16).is_err());
    }

    #[test]
    fn validate_regions_rejects_short_span() {
        assert!(validate_regions(&[region(0, r(0, 14))], 16).is_err());
    }

    #[test]
    fn validate_regions_rejects_index_above_cap() {
        assert!(validate_regions(&[region(0, r(0, 8)), region(2, r(8, 16))], 16).is_err());
    }

    #[test]
    fn ocm_id_relative_range_inside_gel_bead() {
        // We can't construct BarcodeReadComponent directly in this crate (its
        // fields are pub(super)); the function is exercised end-to-end in
        // the bazel pipeline test. Sanity-check the slicing helper here.
        let bytes = b"ACGTACGTACGTACGT";
        let r = r(7, 9);
        assert_eq!(ocm_id_key(bytes, Some(&r)), BcSegSeq::from_bytes(b"TA"));
        assert!(ocm_id_key(bytes, None).is_empty());
    }

    #[test]
    fn expected_distinct_for_uniform_whitelist() {
        // Whitelist of 4 distinct subseqs each appearing once, 8 cells drawn
        // uniformly. lambda_s = 8 * 1 / 4 = 2 for each subseq.
        let mut counts: TxHashMap<BcSegSeq, u64> = TxHashMap::default();
        for s in [b"AA", b"AC", b"AG", b"AT"] {
            counts.insert(BcSegSeq::from_bytes(s), 1);
        }
        let m = expected_distinct_for_ocm_id(Some(&counts), 8, 4);
        let expected = 4.0 * (1.0 - (-2.0_f64).exp());
        assert!((m - expected).abs() < 1e-12);
    }

    #[test]
    fn poisson_log10_cdf_m_zero_is_one() {
        assert_eq!(poisson_log10_cdf(0.0, 0), 0.0);
        assert_eq!(poisson_log10_cdf(0.0, 5), 0.0);
    }

    #[test]
    fn poisson_log10_cdf_k_zero_matches_minus_m_log10e() {
        // P(X = 0) = exp(-m) so log10(P) = -m * log10(e).
        let m = 3.5;
        let got = poisson_log10_cdf(m, 0);
        let expected = -m / std::f64::consts::LN_10;
        assert!((got - expected).abs() < 1e-12);
    }

    #[test]
    fn poisson_log10_cdf_known_values() {
        // P(Poisson(1) <= 1) = 2 / e.
        let got = poisson_log10_cdf(1.0, 1);
        let expected = (2.0 / std::f64::consts::E).log10();
        assert!((got - expected).abs() < 1e-12);

        // P(Poisson(5) <= 5) computed directly: sum_{i=0..=5} e^-5 5^i/i!.
        let mut p = 0.0_f64;
        let mut term = (-5.0_f64).exp();
        p += term;
        for i in 1..=5u64 {
            term *= 5.0 / i as f64;
            p += term;
        }
        let got = poisson_log10_cdf(5.0, 5);
        assert!((got - p.log10()).abs() < 1e-10);
    }

    #[test]
    fn poisson_log10_cdf_left_tail_underflow_is_floored() {
        // m = 100, k = 0 gives p = exp(-100), log10 ~ -43.4, well above floor.
        let got_close = poisson_log10_cdf(100.0, 0);
        assert!((got_close - (-100.0 / std::f64::consts::LN_10)).abs() < 1e-9);
        // m = 100000, k = 0 gives log10 ~ -43429, which would be -inf in a
        // naive linear-space computation; here it gets clamped to LOGP_FLOOR.
        let got_far = poisson_log10_cdf(100_000.0, 0);
        assert_eq!(got_far, LOGP_FLOOR);
    }

    #[test]
    fn poisson_log10_cdf_at_mean_is_about_half() {
        // P(Poisson(m) <= floor(m)) ~ 0.5 for moderate m; log10 ~ log10(0.5).
        let got = poisson_log10_cdf(20.0, 20);
        assert!(got < 0.0);
        assert!(got > -0.5); // log10(0.5) = -0.301
    }
}
