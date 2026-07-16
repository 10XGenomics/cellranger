//! Martian stage WRITE_CRISPR_ONLY_MATRIX
//!
//! Stream the full filtered feature-barcode matrix once and write a CRISPR-only
//! feature-barcode matrix: keep only the CRISPR Guide Capture feature rows (remapped to a
//! contiguous CRISPR-only feature reference) while preserving every barcode column. Every
//! downstream `CALL_PROTOSPACERS` chunk then loads this tiny matrix instead of re-inflating the
//! whole GEX+CRISPR matrix.
//!
//! This reproduces, byte-for-byte at the content level, what the prior Python stage did
//! (`CountMatrix.load_h5_file(full).select_features_by_type("CRISPR Guide Capture").save_h5_file()`):
//!   - feature axis: the CRISPR features in their original relative order, reindexed 0..n_crispr,
//!     with their feature reference (ids/names/tags/sequence) carried through verbatim;
//!   - barcode axis: ALL barcodes preserved in storage order (including ones that have no CRISPR
//!     counts), exactly as `select_features` keeps `self.bcs` unchanged;
//!   - file attributes: only `filetype`/`version`, matching `save_h5_file` called without a
//!     software version (it writes no `software_version`/`chemistry`/`library_ids` attributes).
//!
//! Memory is bounded by the streaming reader plus the (tiny) CRISPR append buffers, NOT by the
//! matrix size: the CRISPR submatrix carries a negligible fraction of the nonzeros, so the whole
//! GEX nnz is streamed past and discarded rather than held resident. The join memory request
//! scales with the feature count (see `split()`) to handle very large guide libraries.
#![deny(missing_docs)]

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
use cr_types::types::FeatureBarcodeType;
use hdf5::types::VarLenUnicode;
use martian::{MartianRover, MartianStage, MartianVoid, Resource, StageDef};
use martian_derive::{MartianStruct, make_mro};
use serde::{Deserialize, Serialize};
use std::str::FromStr;

const FLUSH: usize = 1 << 18; // append-buffer flush threshold

/// WRITE_CRISPR_ONLY_MATRIX stage inputs.
#[derive(Clone, Deserialize, MartianStruct)]
pub struct WriteCrisprOnlyMatrixStageInputs {
    /// The full filtered feature-barcode matrix (all feature types).
    pub filtered_feature_counts_matrix: CountMatrixFile,
}

/// WRITE_CRISPR_ONLY_MATRIX stage outputs.
#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct WriteCrisprOnlyMatrixStageOutputs {
    /// CRISPR-only feature-barcode matrix (all barcodes, CRISPR feature rows only).
    pub crispr_only_matrix: H5File,
}

/// Martian stage WRITE_CRISPR_ONLY_MATRIX.
pub struct WriteCrisprOnlyMatrix;

/// Build the CRISPR-only feature reference (the CRISPR feature defs reindexed 0..n_crispr) and a
/// per-feature remap: `remap[old_idx]` is the new contiguous CRISPR index, or -1 for non-CRISPR
/// features. Mirrors `FeatureReference.select_features(crispr_indices)`: CRISPR features are not
/// target genes, so any target set is dropped.
fn crispr_subset(feature_ref: &FeatureReference) -> (FeatureReference, Vec<i64>) {
    let mut remap = vec![-1i64; feature_ref.feature_defs.len()];
    let mut new_defs: Vec<FeatureDef> = Vec::new();
    // Key `remap` by storage position (`pos`), which is what the streaming reader yields as
    // `rc.feature_idx`; don't rely on `fdef.index` happening to equal it.
    for (pos, fdef) in feature_ref.feature_defs.iter().enumerate() {
        if fdef.feature_type == FeatureType::Barcode(FeatureBarcodeType::Crispr) {
            remap[pos] = new_defs.len() as i64;
            let mut d = fdef.clone();
            d.index = new_defs.len();
            new_defs.push(d);
        }
    }
    let subset = FeatureReference {
        feature_defs: new_defs,
        ..Default::default()
    };
    (subset, remap)
}

#[make_mro(volatile = strict)]
impl MartianStage for WriteCrisprOnlyMatrix {
    type StageInputs = WriteCrisprOnlyMatrixStageInputs;
    type StageOutputs = WriteCrisprOnlyMatrixStageOutputs;
    type ChunkInputs = MartianVoid;
    type ChunkOutputs = MartianVoid;

    fn split(
        &self,
        args: Self::StageInputs,
        _rover: MartianRover,
    ) -> Result<StageDef<Self::ChunkInputs>> {
        // Streaming: the counts are never resident (reader chunk + tiny CRISPR append buffers, a
        // rounding error of total nnz). The one term that scales is the in-memory feature reference
        // (held for the GEX/CRISPR remap + the CRISPR-only feature defs): ~order 1 KB per feature.
        // Negligible at the thousands-of-guides scale (peak ~0.16 GiB at ~1M cells, 7K guides), but a
        // ~1M-feature CRISPR library pushes it to ~1.2 GiB and a fixed 1 GiB request OOMs under
        // `--monitor` (the first stage to fail at extreme library size). Size from the feature count:
        // ~1.5 KB/feature + 0.5 GiB base, floor 1 GiB. (7K->1, 118K->1, 1M->2.)
        let n_features = args
            .filtered_feature_counts_matrix
            .load_dimensions()?
            .num_features;
        let mem_gib = (0.5 + (n_features as f64) * 1536.0 / (1u64 << 30) as f64)
            .ceil()
            .max(1.0);
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
        unreachable!("WRITE_CRISPR_ONLY_MATRIX does all work in join")
    }

    fn join(
        &self,
        args: Self::StageInputs,
        _chunk_defs: Vec<Self::ChunkInputs>,
        _chunk_outs: Vec<Self::ChunkOutputs>,
        rover: MartianRover,
    ) -> Result<Self::StageOutputs> {
        let streaming = args.filtered_feature_counts_matrix.read_streaming_async()?;
        let bcs = streaming.barcodes();
        let n_bc = bcs.len();
        let (crispr_ref, remap) = crispr_subset(streaming.feature_reference());
        let n_crispr = crispr_ref.feature_defs.len();
        if n_crispr == 0 {
            bail!("no CRISPR Guide Capture features found in the input matrix");
        }

        let crispr_only_matrix: H5File = rover.make_path("crispr_only_matrix");
        let file = hdf5::File::create(&crispr_only_matrix)
            .with_context(|| format!("creating {}", crispr_only_matrix.display()))?;
        let group = file.create_group(MATRIX_GROUP)?;
        let data_ds = make_column_ds::<Count>(&group, "data")?;
        let idx_ds = make_column_ds::<BarcodeCountOffset>(&group, "indices")?;

        // Stream counts in barcode order, keeping only CRISPR rows. `raw_counts()` skips
        // zero-count barcodes, so we walk every barcode 0..n_bc and pull its matching run (which
        // may be empty, or contain only discarded non-CRISPR entries), preserving the full barcode
        // axis exactly as `select_features` does.
        let mut it = streaming.raw_counts()?.peekable();
        let mut data_buf: Vec<Count> = Vec::with_capacity(FLUSH);
        let mut idx_buf: Vec<BarcodeCountOffset> = Vec::with_capacity(FLUSH);
        let mut bc_counts: Vec<i64> = Vec::with_capacity(n_bc);
        for b in 0..n_bc {
            let mut kept = 0i64;
            while let Some(peeked) = it.peek() {
                match peeked {
                    Ok(c) if c.barcode_idx == b => {
                        let rc = it.next().unwrap()?;
                        let new_idx = remap[rc.feature_idx];
                        if new_idx >= 0 {
                            data_buf.push(rc.count);
                            idx_buf.push(new_idx as BarcodeCountOffset);
                            kept += 1;
                        }
                    }
                    Ok(_) => break, // belongs to a later barcode
                    Err(_) => {
                        return Err(it.next().unwrap().unwrap_err()).context("streaming counts");
                    }
                }
            }
            bc_counts.push(kept);
            if data_buf.len() >= FLUSH {
                extend_dataset(&data_ds, &data_buf)?;
                data_buf.clear();
                extend_dataset(&idx_ds, &idx_buf)?;
                idx_buf.clear();
            }
        }
        if !data_buf.is_empty() {
            extend_dataset(&data_ds, &data_buf)?;
            extend_dataset(&idx_ds, &idx_buf)?;
        }

        // indptr = [0, cumsum(bc_counts)]
        let mut indptr: Vec<i64> = Vec::with_capacity(n_bc + 1);
        indptr.push(0);
        let mut tot = 0i64;
        for &c in &bc_counts {
            tot += c;
            indptr.push(tot);
        }
        write_column_ds(&group, "indptr", &indptr)?;
        write_column_ds(&group, "shape", &[n_crispr as i32, n_bc as i32])?;
        let bc_ds = make_column_ds::<BarcodeWithGemGroup>(&group, "barcodes")?;
        extend_dataset(&bc_ds, bcs)?;

        let mut fg = group.create_group(FEATURE_REF_GROUP)?;
        feature_reference_io::to_h5(&crispr_ref, &mut fg)?;

        // Match `save_h5_file` called without a software version: only filetype + version.
        scalar_attribute(&file, "filetype", VarLenUnicode::from_str("matrix")?)?;
        scalar_attribute(&file, "version", 2isize)?;
        file.close()?;

        Ok(WriteCrisprOnlyMatrixStageOutputs { crispr_only_matrix })
    }
}
