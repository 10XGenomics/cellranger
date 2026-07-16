//! Martian stage SPLIT_PERTURBATION_MATRIX
//!
//! Stream the filtered feature-barcode matrix once and emit, for the two downstream
//! `MEASURE_PERTURBATIONS` instances (by-feature and by-target), everything they need WITHOUT ever
//! materializing a per-chunk feature×cell matrix:
//!
//!   - **Sufficient statistics** — per (GEX gene × cell-group) aggregates of the sSeq differential
//!     expression (`feature_sums`, `Σw`, `Σw²`, plus per-group `Σcount_total`, `Σ1/count_total`,
//!     cell count `n`, and the per-cell `count_total`s). The DE kernel needs only these, not
//!     cell-resolution data, so the chunk loads tens of MB of statistics instead of hstacking the
//!     full matrix (the old ~33 GB by-target peak). The shared non-targeting control is reduced
//!     ONCE (`control_stats`); each chunk gets its perturbations' statistics (`*_chunk_stats[k]`).
//!     The reducer (`crispr_perturbation_stats`) is bit-identical to the in-memory `reduce_group`.
//!   - **Bootstrap matrices** — the empirical fold-change CIs still resample per-cell target-gene
//!     counts, so a small sparse matrix of the target-gene rows is emitted: `control_bootstrap`
//!     (shared, control cells) + `*_bootstrap[k]` (chunk perturbed cells). The MEASURE chunk
//!     hstacks `[control_bootstrap ++ chunk_bootstrap]` exactly as it used to hstack
//!     `[control_matrix ++ chunk_slice]`, so `group_a`/`group_b` and the bootstrap draws stay
//!     bit-identical. Each carries the per-cell `count_total` (in cell storage order) for the
//!     size factors. The full feature reference rides in `control_stats` (for gene names + the
//!     CRISPR `target_info`).
//!
//! The biology — which barcode is control, dropped, or which perturbation/chunk it belongs to — is
//! decided upstream by `COMPUTE_PRESPLIT_ASSIGNMENT` (reusing the Python clustering) and handed in
//! as a per-mode `reduce_plan`. A single streaming pass scatters the per-group statistics AND
//! appends the target-gene bootstrap cells; peak memory is the two per-mode accumulators
//! (`n_genes × n_groups`), bounded and independent of the matrix's nnz.
#![deny(missing_docs)]

use crate::crispr_perturbation_stats::{GroupStats, PerturbationReducer};
use anyhow::{Context, Result, bail};
use cr_h5::count_matrix::{
    BarcodeCountOffset, BarcodeWithGemGroup, Count, CountMatrix, CountMatrixFile,
    FEATURE_REF_GROUP, MATRIX_GROUP,
};
use cr_h5::{
    extend_dataset, feature_reference_io, make_column_ds, scalar_attribute, write_column_ds,
};
use cr_types::H5File;
use cr_types::reference::feature_reference::{FeatureDef, FeatureReference, FeatureType};
use hdf5::Dataset;
use hdf5::types::{FixedAscii, VarLenAscii, VarLenUnicode};
use martian::{MartianRover, MartianStage, MartianVoid, Resource, StageDef};
use martian_derive::{MartianStruct, make_mro};
use martian_filetypes::FileTypeRead;
use martian_filetypes::json_file::JsonFile;
use ndarray::Array2;
use serde::{Deserialize, Serialize};
use std::path::Path;
use std::str::FromStr;

const STATS_GROUP: &str = "stats";
const FLUSH: usize = 1 << 18; // per-output append-buffer flush threshold

/// Per-mode reduction plan from `COMPUTE_PRESPLIT_ASSIGNMENT`.
#[derive(Clone, Deserialize)]
pub struct ModePlan {
    /// One entry per barcode in matrix storage order: -1 drop, 0 control, 1..=num_perturbations the
    /// perturbation group (in `get_perturbations` order).
    group_assignment: Vec<i32>,
    /// For each perturbation group `p` (1-based; index `p-1`), the 1-based MEASURE chunk that owns
    /// it (0 = not assigned to any chunk, e.g. the non-targeting "perturbation").
    group_to_chunk: Vec<i32>,
    /// Perturbation name for each group `p` (1-based; index `p-1`), in `get_perturbations` order.
    perturbation_names: Vec<String>,
    /// Number of perturbation groups (= max group id).
    num_perturbations: i64,
    /// Number of MEASURE chunks for this mode.
    num_chunks: i64,
}

/// Reduction plan (both modes + the shared bootstrap target-gene rows).
#[derive(Clone, Deserialize)]
pub struct ReducePlan {
    /// GEX rows (0-based among GEX features, storage order) of every targeted gene present — the
    /// feature rows of the bootstrap matrices.
    target_gene_rows: Vec<i64>,
    /// By-feature mode plan.
    by_feature: ModePlan,
    /// By-target mode plan.
    by_target: ModePlan,
}

/// SPLIT_PERTURBATION_MATRIX stage inputs.
#[derive(Clone, Deserialize, MartianStruct)]
pub struct SplitPerturbationMatrixStageInputs {
    /// The full filtered feature-barcode matrix (all feature types).
    pub filtered_feature_counts_matrix: CountMatrixFile,
    /// Per-mode reduction plan from `COMPUTE_PRESPLIT_ASSIGNMENT`.
    pub reduce_plan: JsonFile<ReducePlan>,
}

/// SPLIT_PERTURBATION_MATRIX stage outputs.
#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct SplitPerturbationMatrixStageOutputs {
    /// Shared non-targeting control sufficient statistics + the full feature reference.
    pub control_stats: H5File,
    /// Shared bootstrap matrix: target-gene rows × control cells (+ per-cell count_total).
    pub control_bootstrap: H5File,
    /// Per-chunk perturbation sufficient statistics for the by-feature MEASURE instance.
    pub by_feature_chunk_stats: Vec<H5File>,
    /// Per-chunk perturbation sufficient statistics for the by-target MEASURE instance.
    pub by_target_chunk_stats: Vec<H5File>,
    /// Per-chunk bootstrap matrices (target-gene rows × perturbed cells) for the by-feature instance.
    pub by_feature_bootstrap: Vec<H5File>,
    /// Per-chunk bootstrap matrices (target-gene rows × perturbed cells) for the by-target instance.
    pub by_target_bootstrap: Vec<H5File>,
}

/// Martian stage SPLIT_PERTURBATION_MATRIX.
pub struct SplitPerturbationMatrix;

/// One open bootstrap matrix file being assembled column-by-column (barcode = column), holding only
/// the target-gene rows plus each cell's `count_total` (the size-factor numerator the MEASURE
/// bootstrap divides by `median_ab`).
struct BootstrapOut {
    file: hdf5::File,
    group: hdf5::Group,
    data_ds: Dataset,
    idx_ds: Dataset,
    data_buf: Vec<Count>,
    idx_buf: Vec<BarcodeCountOffset>,
    bc_counts: Vec<i64>,
    bc_list: Vec<BarcodeWithGemGroup>,
    count_total: Vec<f64>,
}

impl BootstrapOut {
    fn create(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let file =
            hdf5::File::create(path).with_context(|| format!("creating {}", path.display()))?;
        let group = file.create_group(MATRIX_GROUP)?;
        let data_ds = make_column_ds::<Count>(&group, "data")?;
        let idx_ds = make_column_ds::<BarcodeCountOffset>(&group, "indices")?;
        Ok(BootstrapOut {
            file,
            group,
            data_ds,
            idx_ds,
            data_buf: Vec::with_capacity(FLUSH),
            idx_buf: Vec::with_capacity(FLUSH),
            bc_counts: Vec::new(),
            bc_list: Vec::new(),
            count_total: Vec::new(),
        })
    }

    /// Append one cell: its target-gene `entries` (`(bootstrap_row, count)`, ascending row) and its
    /// `count_total` (sum over ALL GEX genes — not just the target rows).
    fn push_cell(
        &mut self,
        bc: BarcodeWithGemGroup,
        entries: &[(BarcodeCountOffset, Count)],
        count_total: f64,
    ) -> Result<()> {
        for &(row, count) in entries {
            self.data_buf.push(count);
            self.idx_buf.push(row);
        }
        self.bc_counts.push(entries.len() as i64);
        self.bc_list.push(bc);
        self.count_total.push(count_total);
        if self.data_buf.len() >= FLUSH {
            extend_dataset(&self.data_ds, &self.data_buf)?;
            self.data_buf.clear();
            extend_dataset(&self.idx_ds, &self.idx_buf)?;
            self.idx_buf.clear();
        }
        Ok(())
    }

    fn finalize(&mut self, feature_ref: &FeatureReference, n_features: usize) -> Result<()> {
        if !self.data_buf.is_empty() {
            extend_dataset(&self.data_ds, &self.data_buf)?;
            self.data_buf.clear();
            extend_dataset(&self.idx_ds, &self.idx_buf)?;
            self.idx_buf.clear();
        }
        let mut indptr: Vec<i64> = Vec::with_capacity(self.bc_counts.len() + 1);
        indptr.push(0);
        let mut tot = 0i64;
        for &c in &self.bc_counts {
            tot += c;
            indptr.push(tot);
        }
        write_column_ds(&self.group, "indptr", &indptr)?;
        write_column_ds(
            &self.group,
            "shape",
            &[n_features as i32, self.bc_list.len() as i32],
        )?;
        let bc_ds = make_column_ds::<BarcodeWithGemGroup>(&self.group, "barcodes")?;
        extend_dataset(&bc_ds, &self.bc_list)?;
        // Per-cell count_total in the SAME cell order as the columns (control storage order, then
        // perturbed storage order once hstacked downstream) — the bootstrap size-factor numerators.
        write_column_ds(&self.group, "count_total", &self.count_total)?;
        // Feature reference (the target-gene subset) so `select_features_by_ids` works downstream.
        let mut fg = self.group.create_group(FEATURE_REF_GROUP)?;
        feature_reference_io::to_h5(feature_ref, &mut fg)?;
        write_matrix_file_attributes(&self.file)?;
        Ok(())
    }
}

/// File-level attributes mirroring `cr_h5::count_matrix::write_matrix_h5`, so the bootstrap files
/// load as standard `CountMatrix` H5s downstream.
fn write_matrix_file_attributes(file: &hdf5::File) -> Result<()> {
    scalar_attribute(file, "filetype", VarLenUnicode::from_str("matrix")?)?;
    scalar_attribute(file, "version", 2isize)?;
    scalar_attribute(
        file,
        "software_version",
        VarLenUnicode::from_str("split_perturbation_matrix")?,
    )?;
    scalar_attribute(
        file,
        "chemistry_description",
        VarLenAscii::from_ascii("custom")?,
    )?;
    file.new_attr_builder()
        .with_data(&[FixedAscii::<256>::from_ascii("0")?])
        .create("library_ids")?;
    file.new_attr_builder()
        .with_data(&[1isize])
        .create("original_gem_groups")?;
    Ok(())
}

/// Write a sufficient-statistics H5 holding `groups` (one per group, group-major), with their
/// `names`. When `feature_ref` is `Some`, also embeds the full feature reference under
/// `/matrix/features` (+ matrix file attrs) so `CountMatrix.load_feature_ref_from_h5_file` can read
/// it — used for the shared `control_stats`.
fn write_group_stats(
    path: impl AsRef<Path>,
    groups: &[GroupStats],
    names: &[&str],
    n_genes: usize,
    feature_ref: Option<&FeatureReference>,
) -> Result<()> {
    let path = path.as_ref();
    let file = hdf5::File::create(path).with_context(|| format!("creating {}", path.display()))?;
    let g = file.create_group(STATS_GROUP)?;
    let n = groups.len();

    let feature_sums =
        Array2::<u64>::from_shape_fn((n, n_genes), |(i, j)| groups[i].feature_sums[j]);
    let sum_w = Array2::<f64>::from_shape_fn((n, n_genes), |(i, j)| groups[i].sum_w[j]);
    let sum_wsq = Array2::<f64>::from_shape_fn((n, n_genes), |(i, j)| groups[i].sum_wsq[j]);
    g.new_dataset_builder()
        .with_data(&feature_sums)
        .create("feature_sums")?;
    g.new_dataset_builder().with_data(&sum_w).create("sum_w")?;
    g.new_dataset_builder()
        .with_data(&sum_wsq)
        .create("sum_wsq")?;

    let sum_ct: Vec<f64> = groups.iter().map(|s| s.sum_ct).collect();
    let sum_inv: Vec<f64> = groups.iter().map(|s| s.sum_inv).collect();
    let n_cells: Vec<i64> = groups.iter().map(|s| s.n as i64).collect();
    write_column_ds(&g, "sum_ct", &sum_ct)?;
    write_column_ds(&g, "sum_inv", &sum_inv)?;
    write_column_ds(&g, "n", &n_cells)?;

    // Per-cell count_total concatenated across groups with CSR-style offsets (groups are ragged).
    let mut count_totals: Vec<f64> = Vec::new();
    let mut offsets: Vec<i64> = Vec::with_capacity(n + 1);
    offsets.push(0);
    for s in groups {
        count_totals.extend_from_slice(&s.count_totals);
        offsets.push(count_totals.len() as i64);
    }
    write_column_ds(&g, "count_totals", &count_totals)?;
    write_column_ds(&g, "count_totals_offsets", &offsets)?;

    let names_h5: Vec<VarLenUnicode> = names
        .iter()
        .map(|s| VarLenUnicode::from_str(s))
        .collect::<std::result::Result<_, _>>()?;
    if names_h5.is_empty() {
        g.new_dataset::<VarLenUnicode>()
            .shape(0)
            .create("perturbation_names")?;
    } else {
        g.new_dataset_builder()
            .with_data(&names_h5)
            .create("perturbation_names")?;
    }
    scalar_attribute(&g, "n_genes", n_genes as i64)?;

    if let Some(fref) = feature_ref {
        let matrix_group = file.create_group(MATRIX_GROUP)?;
        let mut fg = matrix_group.create_group(FEATURE_REF_GROUP)?;
        feature_reference_io::to_h5(fref, &mut fg)?;
        write_matrix_file_attributes(&file)?;
    }
    Ok(())
}

/// Build the GEX feature remap (`storage pos -> gex row`, -1 for non-GEX) and the GEX feature defs
/// in storage order. Mirrors `select_features_by_type(GENE_EXPRESSION)`'s row selection.
fn gex_remap(feature_ref: &FeatureReference) -> (Vec<i64>, Vec<FeatureDef>) {
    let mut remap = vec![-1i64; feature_ref.feature_defs.len()];
    let mut gex_defs: Vec<FeatureDef> = Vec::new();
    for (pos, fdef) in feature_ref.feature_defs.iter().enumerate() {
        if fdef.feature_type == FeatureType::Gene {
            remap[pos] = gex_defs.len() as i64;
            gex_defs.push(fdef.clone());
        }
    }
    (remap, gex_defs)
}

/// Build the bootstrap feature reference (the target-gene GEX defs reindexed 0..n_target) and a
/// `gex row -> bootstrap row` map (-1 for non-target genes).
fn bootstrap_subset(
    gex_defs: &[FeatureDef],
    target_gene_rows: &[i64],
) -> (FeatureReference, Vec<i64>) {
    let n_gex = gex_defs.len();
    let mut boot_of_gex = vec![-1i64; n_gex];
    let mut new_defs: Vec<FeatureDef> = Vec::with_capacity(target_gene_rows.len());
    for (boot_row, &gex_row) in target_gene_rows.iter().enumerate() {
        boot_of_gex[gex_row as usize] = boot_row as i64;
        let mut d = gex_defs[gex_row as usize].clone();
        d.index = boot_row;
        new_defs.push(d);
    }
    let subset = FeatureReference {
        feature_defs: new_defs,
        ..Default::default()
    };
    (subset, boot_of_gex)
}

#[make_mro(volatile = strict)]
impl MartianStage for SplitPerturbationMatrix {
    type StageInputs = SplitPerturbationMatrixStageInputs;
    type StageOutputs = SplitPerturbationMatrixStageOutputs;
    type ChunkInputs = MartianVoid;
    type ChunkOutputs = MartianVoid;

    fn split(
        &self,
        args: Self::StageInputs,
        _rover: MartianRover,
    ) -> Result<StageDef<Self::ChunkInputs>> {
        // Peak memory is the two per-mode accumulators plus a bounded overhead; the matrix is never
        // resident (streamed) and the bootstrap matrices stream to H5 (FLUSH-bounded buffers + only
        // per-cell scalars retained). The accumulators are the only term that scales with the data:
        // each holds {Σw, Σw², feature_sums} = 24 B per (GEX gene × group), group-major, so their
        // combined ceiling is `24 · n_gex · (P_by_feature + P_by_target + 2 controls)`. Actual RSS is
        // ≤ that (sparse groups leave zero-init pages un-faulted), so it is a safe upper bound that
        // also scales to larger guide libraries / mouse references — unlike a fixed constant.
        //
        // On top of the accumulators sits a NON-constant overhead that grows with the group/chunk
        // count: async reader prefetch + hdf5 buffers, the (≤ num_chunks per mode) bootstrap writers'
        // FLUSH append buffers, per-cell bootstrap scalars, and the per-chunk `write_group_stats`
        // Array2 transients (≈ (P/num_chunks)·n_gex·24 B). Crucially, at large P the accumulators are
        // ~fully faulted (no sparse savings), so peak RSS ≈ bound + overhead. A FLAT overhead is
        // therefore unsafe at scale: a 1M-cell / 100K-guide synthetic run (P_bf=76,076, P_bt=18,447,
        // n_gex=18,446 → bound 38.97 GiB) actually peaked at 41.3 GiB, OOMing the old `bound + 1.5`
        // (= 41 GiB) request. Scale the margin WITH the bound instead: `1.25·bound + 2`. This covers
        // the growing overhead with headroom while barely moving small cases. Validated points:
        // case1 (peak 2.44, bound 1.86 → req 5), case2 (peak 3.67, bound 3.89 → req 7), synth-100k
        // (peak 41.3, bound 38.97 → req 51).
        const MARGIN: f64 = 1.25;
        const BASE_GIB: f64 = 2.0;
        const FLOOR_GIB: f64 = 2.0;
        let plan = args.reduce_plan.read()?;
        let n_gex = args
            .filtered_feature_counts_matrix
            .load_feature_type_counts()?
            .get(&FeatureType::Gene)
            .copied()
            .unwrap_or(0);
        let n_groups = plan.by_feature.num_perturbations + plan.by_target.num_perturbations + 2;
        let acc_gib = 24.0 * (n_gex as f64) * (n_groups as f64) / (1024.0 * 1024.0 * 1024.0);
        let mem_gib = (acc_gib * MARGIN + BASE_GIB).ceil().max(FLOOR_GIB);
        Ok(StageDef::with_join_resource(Resource::with_mem_gb(
            mem_gib as isize,
        )))
    }

    fn main(
        &self,
        _args: Self::StageInputs,
        _chunk_args: Self::ChunkInputs,
        _rover: MartianRover,
    ) -> Result<Self::ChunkOutputs> {
        unreachable!("SPLIT_PERTURBATION_MATRIX does all work in join")
    }

    fn join(
        &self,
        args: Self::StageInputs,
        _chunk_defs: Vec<Self::ChunkInputs>,
        _chunk_outs: Vec<Self::ChunkOutputs>,
        rover: MartianRover,
    ) -> Result<Self::StageOutputs> {
        let plan = args.reduce_plan.read()?;
        let bf = &plan.by_feature;
        let bt = &plan.by_target;
        let bf_n = bf.num_chunks as usize;
        let bt_n = bt.num_chunks as usize;

        let streaming = args.filtered_feature_counts_matrix.read_streaming_async()?;
        let bcs = streaming.barcodes();
        let feature_ref = streaming.feature_reference();
        let n_bc = bcs.len();
        if bf.group_assignment.len() != n_bc || bt.group_assignment.len() != n_bc {
            bail!(
                "assignment lengths (by_feature {}, by_target {}) != matrix barcode count {}",
                bf.group_assignment.len(),
                bt.group_assignment.len(),
                n_bc
            );
        }

        let (remap, gex_defs) = gex_remap(feature_ref);
        let n_gex = gex_defs.len();
        let (boot_ref, boot_of_gex) = bootstrap_subset(&gex_defs, &plan.target_gene_rows);
        let n_target = boot_ref.feature_defs.len();

        // Reducers: by-feature carries the shared control at group 0 + its perturbations; by-target
        // carries only its perturbations (group 0 unused — control is shared from by-feature).
        let mut bf_reducer = PerturbationReducer::new(n_gex, bf.num_perturbations as usize + 1);
        let mut bt_reducer = PerturbationReducer::new(n_gex, bt.num_perturbations as usize + 1);

        // Output paths.
        let control_stats: H5File = rover.make_path("control_stats");
        let control_bootstrap: H5File = rover.make_path("control_bootstrap");
        let by_feature_chunk_stats: Vec<H5File> = (0..bf_n)
            .map(|k| rover.make_path(format!("by_feature_chunk_stats_{k:03}")))
            .collect();
        let by_target_chunk_stats: Vec<H5File> = (0..bt_n)
            .map(|k| rover.make_path(format!("by_target_chunk_stats_{k:03}")))
            .collect();
        let by_feature_bootstrap: Vec<H5File> = (0..bf_n)
            .map(|k| rover.make_path(format!("by_feature_bootstrap_{k:03}")))
            .collect();
        let by_target_bootstrap: Vec<H5File> = (0..bt_n)
            .map(|k| rover.make_path(format!("by_target_bootstrap_{k:03}")))
            .collect();

        let mut control_boot_out = BootstrapOut::create(&control_bootstrap)?;
        let mut bf_boot_outs: Vec<BootstrapOut> = by_feature_bootstrap
            .iter()
            .map(BootstrapOut::create)
            .collect::<Result<_>>()?;
        let mut bt_boot_outs: Vec<BootstrapOut> = by_target_bootstrap
            .iter()
            .map(BootstrapOut::create)
            .collect::<Result<_>>()?;

        // One streaming pass: scatter per-group statistics AND append target-gene bootstrap cells.
        // `raw_counts()` skips zero-count barcodes, so walk 0..n_bc and pull the matching run for
        // each barcode (handling zero-count gaps), preserving storage order.
        let mut it = streaming.raw_counts()?.peekable();
        let mut gex_entries: Vec<(usize, u32)> = Vec::new();
        let mut target_entries: Vec<(BarcodeCountOffset, Count)> = Vec::new();
        for (b, (&bf_g, &bt_g)) in bf
            .group_assignment
            .iter()
            .zip(&bt.group_assignment)
            .enumerate()
        {
            gex_entries.clear();
            target_entries.clear();
            let mut count_total = 0.0f64;
            while let Some(peeked) = it.peek() {
                match peeked {
                    Ok(c) if c.barcode_idx == b => {
                        let rc = it.next().unwrap()?;
                        let gex_row = remap[rc.feature_idx];
                        if gex_row >= 0 {
                            // Counts are non-negative by matrix invariant (see `convert_matrix`'s
                            // transmute assert); the reducer's contract is u32.
                            gex_entries.push((gex_row as usize, rc.count as u32));
                            count_total += f64::from(rc.count);
                            let boot_row = boot_of_gex[gex_row as usize];
                            if boot_row >= 0 {
                                target_entries.push((boot_row as BarcodeCountOffset, rc.count));
                            }
                        }
                    }
                    Ok(_) => break, // belongs to a later barcode
                    Err(_) => {
                        return Err(it.next().unwrap().unwrap_err()).context("streaming counts");
                    }
                }
            }

            // A control cell is control in both modes; it can't also be a perturbation. The shared
            // control is reduced once (by-feature group 0) and written to control_bootstrap.
            if bf_g == 0 || bt_g == 0 {
                if bf_g > 0 || bt_g > 0 {
                    bail!(
                        "barcode {b} is control in one mode but perturbed in the other (by_feature={bf_g}, by_target={bt_g})"
                    );
                }
                bf_reducer.add_cell(0, &gex_entries);
                control_boot_out.push_cell(bcs[b], &target_entries, count_total)?;
                continue;
            }
            if bf_g >= 1 {
                bf_reducer.add_cell(bf_g as usize, &gex_entries);
                let k = bf.group_to_chunk[(bf_g - 1) as usize];
                if k >= 1 {
                    bf_boot_outs[(k - 1) as usize].push_cell(
                        bcs[b],
                        &target_entries,
                        count_total,
                    )?;
                }
            }
            if bt_g >= 1 {
                bt_reducer.add_cell(bt_g as usize, &gex_entries);
                let k = bt.group_to_chunk[(bt_g - 1) as usize];
                if k >= 1 {
                    bt_boot_outs[(k - 1) as usize].push_cell(
                        bcs[b],
                        &target_entries,
                        count_total,
                    )?;
                }
            }
        }

        // Finalize bootstrap matrices.
        control_boot_out.finalize(&boot_ref, n_target)?;
        for out in &mut bf_boot_outs {
            out.finalize(&boot_ref, n_target)?;
        }
        for out in &mut bt_boot_outs {
            out.finalize(&boot_ref, n_target)?;
        }

        // Shared control statistics (by-feature group 0) + the full feature reference.
        let control = bf_reducer.finalize_group(0);
        write_group_stats(
            &control_stats,
            &[control],
            &["Non-Targeting"],
            n_gex,
            Some(feature_ref),
        )?;

        // Per-chunk perturbation statistics, grouped by `group_to_chunk`.
        write_chunk_stats(&by_feature_chunk_stats, bf, &bf_reducer, n_gex)?;
        write_chunk_stats(&by_target_chunk_stats, bt, &bt_reducer, n_gex)?;

        Ok(SplitPerturbationMatrixStageOutputs {
            control_stats,
            control_bootstrap,
            by_feature_chunk_stats,
            by_target_chunk_stats,
            by_feature_bootstrap,
            by_target_bootstrap,
        })
    }
}

/// Write each chunk's perturbation statistics, bucketing perturbation groups by `group_to_chunk`.
fn write_chunk_stats(
    paths: &[H5File],
    mode: &ModePlan,
    reducer: &PerturbationReducer,
    n_gex: usize,
) -> Result<()> {
    for (k, path) in paths.iter().enumerate() {
        let chunk = (k + 1) as i32;
        let mut groups: Vec<GroupStats> = Vec::new();
        let mut names: Vec<&str> = Vec::new();
        for p in 1..=mode.num_perturbations as usize {
            if mode.group_to_chunk[p - 1] == chunk {
                groups.push(reducer.finalize_group(p));
                names.push(&mode.perturbation_names[p - 1]);
            }
        }
        write_group_stats(path, &groups, &names, n_gex, None)?;
    }
    Ok(())
}
