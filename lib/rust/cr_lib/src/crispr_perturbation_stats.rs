//! Streaming sufficient-statistics reducer for the CRISPR perturbation differential expression.
//!
//! `MEASURE_PERTURBATIONS`'s sSeq DE reduces to per-(gene, cell-group) sufficient statistics: the
//! kernel never needs cell-resolution data, only additive sums over each group's cells (the control
//! once + each perturbation). `fast_utils::diff_exp::reduce_group` computes those from an in-memory
//! matrix; this module computes the IDENTICAL statistics from a single streaming pass over the full
//! feature-barcode matrix, so `SPLIT_PERTURBATION_MATRIX` can emit them instead of materializing the
//! giant per-chunk matrices the old MEASURE chunk hstacked (the ~33 GB by-target peak).
//!
//! ## Bit-identity
//! The reducer is bit-for-bit (`rtol=0`) equal to `reduce_group`, not merely tolerant. `reduce_group`
//! derives `sum_w`/`sum_wsq` from `sqz`'s `mean_var_rows`, whose CSC branch accumulates per gene by
//! sweeping the group's columns in ascending column order, ascending row order within a column
//! (scan-rs `sqz/src/mat.rs::mean_var_rows`). A full-matrix stream (barcodes `0..n_bc`, features
//! ascending within a barcode) visits each group's cells in exactly that order, so the `f64` add
//! order — and therefore the result — is identical. `finalize` reproduces `reduce_group`'s
//! `sum_w = mean·n`, `sum_wsq = (var + mean²)·n` reconstruction verbatim. `feature_sums`, `sum_ct`,
//! `sum_inv`, `count_totals`, and `n` are exact by construction. This is proven by `tests` below and
//! is the basis for the bit-identical Part-2 QA gate.

use serde::{Deserialize, Serialize};

/// Per-group sufficient statistics, the streaming analogue of `fast_utils`'s `reduce_group` output.
/// Every field is an additive sum over the group's cells, in matrix storage order.
#[derive(Clone, Serialize, Deserialize)]
pub struct GroupStats {
    /// Σ raw count per gene over the group's cells (`n_genes`).
    pub feature_sums: Vec<u64>,
    /// Σ w per gene, `w = count/count_total`; stored as `mean·n` to match the matrix path (`n_genes`).
    pub sum_w: Vec<f64>,
    /// Σ w² per gene; stored as `(var + mean²)·n` to match the matrix path (`n_genes`).
    pub sum_wsq: Vec<f64>,
    /// Σ count_total over the group's cells.
    pub sum_ct: f64,
    /// Σ 1/count_total over the group's cells (cells with zero total skipped).
    pub sum_inv: f64,
    /// Number of cells in the group (includes cells with zero GEX counts).
    pub n: usize,
    /// Per-cell count_total in the group's storage order (for the `median_ab` over a ∪ control).
    pub count_totals: Vec<f64>,
}

/// Accumulates per-(gene, group) sufficient statistics over a streamed matrix for ONE mode.
///
/// Groups are indexed `0..n_groups`, where `0` is the shared non-targeting control and `1..n_groups`
/// are the perturbations (the same group ids `COMPUTE_PRESPLIT_ASSIGNMENT` assigns per cell). The
/// per-(group, gene) accumulators are laid out group-major (`group*n_genes + gene`) so each cell is a
/// sequential sweep of one gene block — the layout the Phase-0 wall-clock prototype confirmed.
pub struct PerturbationReducer {
    n_genes: usize,
    n_groups: usize,
    /// Σ w, group-major `[group*n_genes + gene]`.
    acc_w: Vec<f64>,
    /// Σ w², group-major `[group*n_genes + gene]`.
    acc_wsq: Vec<f64>,
    /// Σ count, group-major `[group*n_genes + gene]`.
    feature_sums: Vec<u64>,
    sum_ct: Vec<f64>,
    sum_inv: Vec<f64>,
    n: Vec<usize>,
    count_totals: Vec<Vec<f64>>,
}

impl PerturbationReducer {
    /// Create a reducer for `n_groups` groups (control + perturbations) over `n_genes` genes.
    pub fn new(n_genes: usize, n_groups: usize) -> Self {
        PerturbationReducer {
            n_genes,
            n_groups,
            acc_w: vec![0.0; n_genes * n_groups],
            acc_wsq: vec![0.0; n_genes * n_groups],
            feature_sums: vec![0; n_genes * n_groups],
            sum_ct: vec![0.0; n_groups],
            sum_inv: vec![0.0; n_groups],
            n: vec![0; n_groups],
            count_totals: vec![Vec::new(); n_groups],
        }
    }

    /// Add one cell (assigned to `group`) given its GEX `entries` as `(gene_idx, count)` pairs in
    /// ascending gene order (matrix storage order). The cell's `count_total` is the sum of its GEX
    /// counts — matching `compute_count_total` (a sum over the GEX-only matrix's rows) — so `w` and
    /// every aggregate equal the in-memory path. Must be called for EVERY assigned barcode, including
    /// cells with no GEX entries (`count_total == 0`), to keep `n`/`count_totals` aligned with
    /// `reduce_group` over the same cell set.
    pub fn add_cell(&mut self, group: usize, entries: &[(usize, u32)]) {
        debug_assert!(group < self.n_groups);
        let count_total: f64 = entries.iter().map(|&(_, v)| f64::from(v)).sum();
        let base = group * self.n_genes;
        // `CountNormalized` maps to 0 when count_total == 0; entries are empty in that case anyway.
        if count_total != 0.0 {
            for &(gene, v) in entries {
                let w = f64::from(v) / count_total;
                self.acc_w[base + gene] += w;
                self.acc_wsq[base + gene] += w * w;
                self.feature_sums[base + gene] += u64::from(v);
            }
            self.sum_inv[group] += 1.0 / count_total;
        }
        self.sum_ct[group] += count_total;
        self.n[group] += 1;
        self.count_totals[group].push(count_total);
    }

    /// Finalize one group's statistics, reproducing `reduce_group`'s reconstruction exactly:
    /// `mean = Σw/n`, `sum_w = mean·n`, `var = Σw²/n − mean²`, `sum_wsq = (var + mean²)·n`. With
    /// `n == 0` the group is empty and every aggregate is zero (matching `mean_var_rows` over no
    /// columns, which leaves the zero-initialised accumulators untouched).
    pub fn finalize_group(&self, group: usize) -> GroupStats {
        let n = self.n[group];
        let nf = n as f64;
        let base = group * self.n_genes;
        let mut sum_w = vec![0.0; self.n_genes];
        let mut sum_wsq = vec![0.0; self.n_genes];
        if n > 0 {
            for g in 0..self.n_genes {
                let mean = self.acc_w[base + g] / nf;
                let var = self.acc_wsq[base + g] / nf - mean * mean;
                sum_w[g] = mean * nf;
                sum_wsq[g] = (var + mean * mean) * nf;
            }
        }
        GroupStats {
            feature_sums: self.feature_sums[base..base + self.n_genes].to_vec(),
            sum_w,
            sum_wsq,
            sum_ct: self.sum_ct[group],
            sum_inv: self.sum_inv[group],
            n,
            count_totals: self.count_totals[group].clone(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Reference reduction copied verbatim from the production in-memory path, so the test pins the
    /// reducer to the exact arithmetic `MEASURE` runs today:
    ///   * `mean_var_rows` (scan-rs `sqz/src/mat.rs:328-342`, CSC branch): for each col (ascending),
    ///     for each (row, value) ascending, `mean[g] += w; var[g] += w²`; then `mean /= m`,
    ///     `var = var/m − mean²`.
    ///   * `reduce_group` (fast_utils `diff_exp.rs`): `sum_w = mean·n`, `sum_wsq = (var + mean²)·n`,
    ///     `count_total = Σ_g count` per cell, `sum_ct`/`sum_inv` folds.
    ///
    /// `cells` are passed in ascending order (as `np.flatnonzero` yields), matching the stream.
    fn reference_group(
        dense: &[Vec<u32>], // [cell][gene]
        n_genes: usize,
        cells: &[usize],
    ) -> GroupStats {
        let count_total: Vec<f64> = dense
            .iter()
            .map(|row| row.iter().map(|&v| f64::from(v)).sum())
            .collect();
        let m = cells.len() as f64;
        let mut mean = vec![0.0f64; n_genes];
        let mut var = vec![0.0f64; n_genes];
        // mean_var_rows CSC accumulation order: outer over cols (cells) ascending, inner over genes.
        for &c in cells {
            let ct = count_total[c];
            for g in 0..n_genes {
                let v = dense[c][g];
                if v == 0 {
                    continue; // sparse: implicit zeros are not visited
                }
                let w = if ct == 0.0 { 0.0 } else { f64::from(v) / ct };
                mean[g] += w;
                var[g] += w * w;
            }
        }
        let mut sum_w = vec![0.0f64; n_genes];
        let mut sum_wsq = vec![0.0f64; n_genes];
        if !cells.is_empty() {
            for g in 0..n_genes {
                mean[g] /= m;
                var[g] = var[g] / m - mean[g] * mean[g];
                sum_w[g] = mean[g] * m;
                sum_wsq[g] = (var[g] + mean[g] * mean[g]) * m;
            }
        }
        let mut feature_sums = vec![0u64; n_genes];
        for &c in cells {
            for g in 0..n_genes {
                feature_sums[g] += u64::from(dense[c][g]);
            }
        }
        let sum_ct: f64 = cells.iter().map(|&c| count_total[c]).sum();
        let sum_inv: f64 = cells
            .iter()
            .filter(|&&c| count_total[c] != 0.0)
            .map(|&c| 1.0 / count_total[c])
            .sum();
        let count_totals: Vec<f64> = cells.iter().map(|&c| count_total[c]).collect();
        GroupStats {
            feature_sums,
            sum_w,
            sum_wsq,
            sum_ct,
            sum_inv,
            n: cells.len(),
            count_totals,
        }
    }

    /// Deterministic dense genes×cells data with non-zero per-gene variance and some zero entries
    /// (to exercise the sparse-skip path). Mirrors the construction in the `diff_exp` batch tests.
    fn build_dense(n_genes: usize, n_cells: usize) -> Vec<Vec<u32>> {
        (0..n_cells)
            .map(|c| {
                (0..n_genes)
                    .map(|g| {
                        let v = (c * 5 + g * 7 + 3) % 11;
                        // knock out a deterministic ~1/4 of entries to create real sparsity
                        if (c + g) % 4 == 0 { 0 } else { v as u32 }
                    })
                    .collect()
            })
            .collect()
    }

    fn assert_bit_identical(got: &GroupStats, want: &GroupStats) {
        assert_eq!(got.n, want.n, "n");
        assert_eq!(got.feature_sums, want.feature_sums, "feature_sums");
        assert_eq!(got.count_totals, want.count_totals, "count_totals (rtol=0)");
        assert_eq!(got.sum_ct, want.sum_ct, "sum_ct (rtol=0)");
        assert_eq!(got.sum_inv, want.sum_inv, "sum_inv (rtol=0)");
        for g in 0..want.sum_w.len() {
            assert_eq!(got.sum_w[g], want.sum_w[g], "sum_w g={g} (rtol=0)");
            assert_eq!(got.sum_wsq[g], want.sum_wsq[g], "sum_wsq g={g} (rtol=0)");
        }
    }

    /// The streaming reducer must reproduce the in-memory `reduce_group` reference BIT-FOR-BIT
    /// (rtol=0) — including `sum_w`/`sum_wsq`, the only quantities the design feared would drift.
    /// Cells are fed in full-matrix storage order (barcode 0..n_bc), each routed to its group,
    /// exactly as `SPLIT_PERTURBATION_MATRIX` will drive it.
    #[test]
    fn test_streaming_matches_reduce_group_bit_identical() {
        let n_genes = 40usize;
        let n_cells = 200usize;
        let dense = build_dense(n_genes, n_cells);

        // group 0 = control (cells 0..120), group 1 = pert A (120..160), group 2 = pert B (160..200).
        // Interleave is irrelevant here; what matters is each group's cells stream in ascending order.
        let control: Vec<usize> = (0..120).collect();
        let pert_a: Vec<usize> = (120..160).collect();
        let pert_b: Vec<usize> = (160..200).collect();
        let group_of = |c: usize| -> usize {
            if c < 120 {
                0
            } else if c < 160 {
                1
            } else {
                2
            }
        };

        let mut reducer = PerturbationReducer::new(n_genes, 3);
        // Stream in storage order, building each cell's (gene, count) entry list (ascending gene).
        for (c, row) in dense.iter().enumerate() {
            let entries: Vec<(usize, u32)> = row
                .iter()
                .enumerate()
                .filter(|&(_, &v)| v != 0)
                .map(|(g, &v)| (g, v))
                .collect();
            reducer.add_cell(group_of(c), &entries);
        }

        assert_bit_identical(
            &reducer.finalize_group(0),
            &reference_group(&dense, n_genes, &control),
        );
        assert_bit_identical(
            &reducer.finalize_group(1),
            &reference_group(&dense, n_genes, &pert_a),
        );
        assert_bit_identical(
            &reducer.finalize_group(2),
            &reference_group(&dense, n_genes, &pert_b),
        );
    }

    /// A cell with zero GEX counts assigned to a group must still count toward `n`/`count_totals`
    /// (count_total = 0), matching `reduce_group` over a cell set that includes it.
    #[test]
    fn test_zero_count_cell_counted() {
        let n_genes = 5usize;
        // cell 0 all-zero, cells 1,2 nonzero.
        let dense = vec![
            vec![0, 0, 0, 0, 0],
            vec![1, 0, 2, 0, 3],
            vec![0, 4, 0, 5, 0],
        ];
        let cells: Vec<usize> = vec![0, 1, 2];
        let mut reducer = PerturbationReducer::new(n_genes, 1);
        for row in &dense {
            let entries: Vec<(usize, u32)> = row
                .iter()
                .enumerate()
                .filter(|&(_, &v)| v != 0)
                .map(|(g, &v)| (g, v))
                .collect();
            reducer.add_cell(0, &entries);
        }
        let got = reducer.finalize_group(0);
        let want = reference_group(&dense, n_genes, &cells);
        assert_bit_identical(&got, &want);
        assert_eq!(got.n, 3);
        assert_eq!(got.count_totals[0], 0.0);
    }
}
