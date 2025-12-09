//! A feature-barcode count matrix.
#![expect(missing_docs)]

use crate::iter::{H5Iterator, ThreadedH5Iterator};
use crate::{
    ChunkedWriter, ColumnAction, extend_dataset, feature_reference_io, make_column_ds,
    scalar_attribute, write_column_ds,
};
use anyhow::{Context, Result, bail};
use barcode::{Barcode, MAX_BARCODE_AND_GEM_GROUP_LENGTH};
use cr_types::barcode_index::BarcodeIndex;
use cr_types::reference::feature_reference::{FeatureDef, FeatureReference, FeatureType};
use cr_types::{BarcodeThenFeatureOrder, CountShardFile, FeatureBarcodeCount, GemWell, H5File};
use hdf5::Group;
use hdf5::types::{FixedAscii, VarLenAscii, VarLenUnicode};
use itertools::{Itertools, zip_eq};
use martian_derive::martian_filetype;
use serde::{Deserialize, Serialize};
use shardio::ShardReader;
use std::io::Write;
use std::iter::zip;
use std::mem::size_of;
use std::ops::Range;
use std::str::FromStr;

// Group keys.
const MATRIX_GROUP: &str = "matrix";
const FEATURE_REF_GROUP: &str = "features";

const MATRIX_SHAPE_DATASET: &str = "shape";

// Dataset keys.
const DATASET_DATA: &str = "data";
const DATASET_BARCODES: &str = "barcodes";
const DATASET_BARCODE_COUNT_OFFSETS: &str = "indptr";
const DATASET_FEATURE_INDICES: &str = "indices";

pub const MAT_H5_BUF_SZ: usize = 1 << 20;

const H5_FILETYPE_KEY: &str = "filetype";
const MATRIX_H5_FILETYPE: &str = "matrix";
const MATRIX_H5_VERSION_KEY: &str = "version";
const MATRIX_H5_VERSION: isize = 2;
pub const SOFTWARE_H5_VERSION_KEY: &str = "software_version";
const H5_CHEMISTRY_DESC_KEY: &str = "chemistry_description";
const H5_LIBRARY_ID_MAPPING_KEY: &str = "library_ids";
const H5_ORIG_GEM_GROUP_MAPPING_KEY: &str = "original_gem_groups";

/// The number of times a particular (barcode, feature) pair was observed.
pub type Count = i32;
/// The index into the array of features for a particular count.
pub type FeatureIdx = i32;
/// A value representing the offset into the count array corresponding to
/// the beginning of the count range for that barcode.
pub type BarcodeCountOffset = i64;
/// A barcode, including the -N GEM well suffix.
pub type BarcodeWithGemGroup = FixedAscii<MAX_BARCODE_AND_GEM_GROUP_LENGTH>;

martian_filetype!(CountMatrixFile, "h5");

/// A feature-barcode count matrix HDF5 file.
impl CountMatrixFile {
    fn matrix_group(&self) -> Result<Group> {
        Ok(hdf5::File::open(self)
            .with_context(|| format!("While opening {self:?}"))?
            .group(MATRIX_GROUP)?)
    }

    pub fn load_dimensions(&self) -> Result<MatrixDimensions> {
        let matrix_group = self.matrix_group()?;
        let shape: Vec<i32> = matrix_group.dataset(MATRIX_SHAPE_DATASET)?.read_raw()?;
        let [num_features, num_barcodes] = shape.as_slice() else {
            unreachable!()
        };
        let nnz = matrix_group.dataset(DATASET_DATA)?.size();
        Ok(MatrixDimensions {
            num_features: *num_features as usize,
            num_barcodes: *num_barcodes as usize,
            num_non_zeros: nnz,
        })
    }

    // Return the amount of memory in GiB required for loading the matrix.
    pub fn estimate_mem_gib(&self) -> Result<f64> {
        Ok(self.load_dimensions()?.estimate_mem_gib())
    }

    // Return the amount of memory in GiB required for loading the matrix in streaming mode.
    pub fn estimate_mem_gib_without_counts(&self) -> Result<f64> {
        let bcs_and_feats_gib = self.load_dimensions()?.estimate_mem_gib_without_counts();
        // Add memory for a loaded chunk.
        let chunks_gib = ((size_of::<Count>() + size_of::<FeatureIdx>()) * MAT_H5_BUF_SZ) as f64
            / (1024. * 1024. * 1024.);
        Ok(bcs_and_feats_gib + chunks_gib)
    }

    /// Read all matrix data into memory.
    pub fn read(&self) -> Result<CountMatrixEager> {
        // Augment any read error with the path to the file.
        self.read_inner()
            .with_context(|| self.display().to_string())
    }

    /// Create a streaming reader.
    pub fn read_streaming(&self) -> Result<CountMatrixStreaming> {
        // Augment any read error with the path to the file.
        self.read_streaming_inner(false)
            .with_context(|| self.display().to_string())
    }

    /// Create a streaming reader that will use threads to read the matrix.
    pub fn read_streaming_async(&self) -> Result<CountMatrixStreaming> {
        // Augment any read error with the path to the file.
        self.read_streaming_inner(true)
            .with_context(|| self.display().to_string())
    }

    fn read_streaming_inner(&self, use_async: bool) -> Result<CountMatrixStreaming> {
        let bcs_and_features = Self::read_bcs_and_features(&self.matrix_group()?)?;
        Ok(CountMatrixStreaming {
            file: self.clone(),
            bcs_and_features,
            use_async,
        })
    }

    pub fn read_barcodes(&self) -> Result<Vec<BarcodeWithGemGroup>> {
        let matrix = self.matrix_group()?;
        matrix
            .dataset(DATASET_BARCODES)?
            .read_raw()
            .context("Ran into trouble reading barcodes from feature bc matrix.")
    }

    fn read_inner(&self) -> Result<CountMatrixEager> {
        let matrix = self.matrix_group()?;

        Ok(CountMatrixEager {
            counts: matrix.dataset(DATASET_DATA)?.read_raw()?,
            feature_indices: matrix.dataset(DATASET_FEATURE_INDICES)?.read_raw()?,
            bcs_and_features: Self::read_bcs_and_features(&matrix)?,
        })
    }

    fn read_bcs_and_features(matrix: &Group) -> Result<BarcodesAndFeatures> {
        Ok(BarcodesAndFeatures {
            barcodes: matrix.dataset(DATASET_BARCODES)?.read_raw()?,
            barcode_count_offsets: matrix.dataset(DATASET_BARCODE_COUNT_OFFSETS)?.read_raw()?,
            feature_reference: feature_reference_io::from_h5(&matrix.group(FEATURE_REF_GROUP)?)?,
        })
    }
}

pub struct MatrixDimensions {
    pub num_features: usize,
    pub num_barcodes: usize,
    pub num_non_zeros: usize,
}

impl MatrixDimensions {
    // Return the amount of memory in bytes required for the matrix.
    fn estimate_mem_bytes(&self) -> usize {
        size_of::<Count>() * self.num_non_zeros + // counts (i32)
        size_of::<FeatureIdx>() * self.num_non_zeros + // feature_indices (i32)
        self.estimate_mem_bytes_without_counts()
    }

    fn estimate_mem_bytes_without_counts(&self) -> usize {
        size_of::<BarcodeWithGemGroup>() * self.num_barcodes + // barcodes
        size_of::<BarcodeCountOffset>() * (1 + self.num_barcodes) + // barcode_count_offsets (i64)
        FeatureReference::estimate_mem_bytes_for_feature_count(self.num_features)
    }

    /// Return the amount of memory in GiB required for the matrix.
    pub fn estimate_mem_gib(&self) -> f64 {
        (self.estimate_mem_bytes() as f64) / (1024.0 * 1024.0 * 1024.0)
    }

    /// Return the amount of memory in GiB required for the matrix metadata.
    pub fn estimate_mem_gib_without_counts(&self) -> f64 {
        (self.estimate_mem_bytes_without_counts() as f64) / (1024.0 * 1024.0 * 1024.0)
    }
}

/// Count matrix that streams counts from disk rather than eagerly loading them.
///
/// The barcode and feature data is eagerly loaded.
///
/// If use_async is true, short-lived threads will be spawned for loading the
/// counts in the background. One chunk will be buffered at a time.
pub struct CountMatrixStreaming {
    #[allow(unused)]
    file: CountMatrixFile,
    pub bcs_and_features: BarcodesAndFeatures,
    use_async: bool,
}

impl CountMatrixStreaming {
    /// Return an iterator over raw counts.
    pub fn raw_counts(&self) -> Result<impl Iterator<Item = Result<RawCount>>> {
        self.raw_counts_inner()
            .with_context(|| self.file.display().to_string())
    }

    /// Return an iterator over full counts including metadata.
    pub fn counts(&self) -> Result<impl Iterator<Item = Result<AnnotatedCount<'_>>>> {
        Ok(self.raw_counts()?.map(|r| {
            let raw_count = r?;
            Ok(AnnotatedCount {
                count: raw_count.count,
                barcode: &self.bcs_and_features.barcodes[raw_count.barcode_idx],
                feature: &self.bcs_and_features.feature_reference.feature_defs
                    [raw_count.feature_idx],
            })
        }))
    }

    /// Return an iterator over all observed barcodes and their total counts
    /// for the specified feature type.
    /// Barcodes with zero counts for that feature type will still be included
    /// in the output, with a count of zero.
    pub fn barcode_counts_for_feature_type(
        &self,
        feature_type: FeatureType,
    ) -> Result<impl Iterator<Item = Result<(&BarcodeWithGemGroup, i64)>>> {
        self.filtered_barcode_counts(move |count| count.feature.feature_type == feature_type)
    }

    /// Return an iterator over all barcodes and their total observed counts.
    /// Apply the provided filter on each individual count to determine if they
    /// participate in the sum.
    /// Barcodes that have zero counts after filtering will still be included
    /// in the output, with a count of zero.
    /// Note that this will not include barcodes present in the barcodes array
    /// but which have no counts even before filtering.
    fn filtered_barcode_counts<F: Fn(&AnnotatedCount<'_>) -> bool>(
        &self,
        filter: F,
    ) -> Result<impl Iterator<Item = Result<(&BarcodeWithGemGroup, i64)>>> {
        // We can't use Itertools::chunk_by because we're returning an iterator.
        // This is basically a custom chunker iterator that folds each group.
        let mut bc_idx = Box::new(0usize);
        let mut sum = Box::new(0i64);

        let mut raw_counts_iter = self.raw_counts()?.peekable();

        Ok(std::iter::from_fn(move || {
            *sum = 0;

            // Fetch an initial value. If there are none left, we're done.
            let raw_count = match raw_counts_iter.next() {
                None => return None,
                Some(Ok(count)) => count,
                Some(Err(err)) => {
                    return Some(Err(err));
                }
            };
            *bc_idx = raw_count.barcode_idx;
            if filter(&self.annotate_count(&raw_count)) {
                *sum += raw_count.count as i64;
            }
            // Process values until we hit the next barcode or run out.
            loop {
                let Some(next_item) = raw_counts_iter
                    .next_if(|res| res.as_ref().is_ok_and(|count| count.barcode_idx == *bc_idx))
                else {
                    // Done with this barcode, or all barcodes, or we hit an error.
                    // Yield the current barcode and sum.
                    // If this was an error, it will be surfaced on the next iteration.
                    return Some(Ok((&self.bcs_and_features.barcodes[*bc_idx], *sum)));
                };
                let next_count = match next_item {
                    Err(err) => {
                        // This should be unreachable. Handle it anyway.
                        return Some(Err(err));
                    }
                    Ok(count) => count,
                };
                // Sum if passes filter.
                if filter(&self.annotate_count(&next_count)) {
                    *sum += next_count.count as i64;
                }
            }
        }))
    }

    fn annotate_count(&self, raw_count: &RawCount) -> AnnotatedCount<'_> {
        AnnotatedCount {
            count: raw_count.count,
            barcode: &self.bcs_and_features.barcodes[raw_count.barcode_idx],
            feature: &self.bcs_and_features.feature_reference.feature_defs[raw_count.feature_idx],
        }
    }

    fn raw_counts_inner(&self) -> Result<impl Iterator<Item = Result<RawCount>>> {
        let matrix = self.file.matrix_group()?;
        if self.use_async {
            let bare_counts_iter =
                ThreadedH5Iterator::new(matrix.dataset(DATASET_DATA)?, MAT_H5_BUF_SZ, 1);
            let feat_idx_iter =
                ThreadedH5Iterator::new(matrix.dataset(DATASET_FEATURE_INDICES)?, MAT_H5_BUF_SZ, 1);
            return self
                .raw_counts_from_iters(bare_counts_iter, feat_idx_iter)
                .map(|iter| Box::new(iter) as Box<dyn Iterator<Item = Result<RawCount>>>);
        }
        let bare_counts_iter = H5Iterator::new(matrix.dataset(DATASET_DATA)?, MAT_H5_BUF_SZ);
        let feat_idx_iter =
            H5Iterator::new(matrix.dataset(DATASET_FEATURE_INDICES)?, MAT_H5_BUF_SZ);
        self.raw_counts_from_iters(bare_counts_iter, feat_idx_iter)
            .map(|iter| Box::new(iter) as Box<dyn Iterator<Item = Result<RawCount>>>)
    }

    fn raw_counts_from_iters<C, F>(
        &self,
        bare_counts_iter: C,
        feat_idx_iter: F,
    ) -> Result<impl Iterator<Item = Result<RawCount>>>
    where
        C: Iterator<Item = Result<Count>>,
        F: Iterator<Item = Result<FeatureIdx>>,
    {
        // Assign barcode IDs by repeating the index as many times as expected
        // by each barcode range.
        let bc_idx_iter = self
            .bcs_and_features
            .barcode_count_offsets
            .iter()
            .copied()
            .tuple_windows()
            .enumerate()
            .flat_map(|(bc_idx, (start, end))| {
                let count = (end - start).max(0);
                std::iter::repeat_n(bc_idx, count as usize)
            });

        Ok(zip_eq(bare_counts_iter, feat_idx_iter)
            .zip_eq(bc_idx_iter)
            .map(|((count_result, feature_idx_result), barcode_idx)| {
                let count = count_result.with_context(|| self.file.display().to_string())?;
                let feature_idx =
                    feature_idx_result.with_context(|| self.file.display().to_string())?;
                Ok(RawCount {
                    count,
                    barcode_idx,
                    feature_idx: feature_idx as usize,
                })
            }))
    }
}

impl CountMatrix for CountMatrixStreaming {}

/// The subset of the matrix that isn't the actual counts themselves.
///
/// This makes sense to break out because this is usually the minority of
/// memory consumed by the matrix and thus can be eagerly loaded.
pub struct BarcodesAndFeatures {
    barcodes: Vec<BarcodeWithGemGroup>,
    barcode_count_offsets: Vec<BarcodeCountOffset>,
    pub feature_reference: FeatureReference,
}

impl BarcodesAndFeatures {
    /// Return an iterator over all barcodes and their range in the counts array.
    fn iter_barcode_ranges(&self) -> impl Iterator<Item = (&BarcodeWithGemGroup, Range<usize>)> {
        zip(
            &self.barcodes,
            self.barcode_count_offsets.iter().tuple_windows(),
        )
        .map(|(barcode, (&start, &end))| (barcode, (start as usize..end as usize)))
    }
}

/// In-memory representation of the count matrix.
/// All data is eagerly loaded.
pub struct CountMatrixEager {
    counts: Vec<Count>,
    feature_indices: Vec<FeatureIdx>,
    bcs_and_features: BarcodesAndFeatures,
}

impl CountMatrixEager {
    /// Return an iterator over all observed barcodes and their total counts
    /// for the specified feature type.
    /// Barcodes with zero counts for that feature type will still be included
    /// in the output, with a count of zero.
    pub fn barcode_counts_for_feature_type(
        &self,
        feature_type: FeatureType,
    ) -> impl Iterator<Item = (&BarcodeWithGemGroup, i64)> {
        self.filtered_barcode_counts(move |count| count.feature.feature_type == feature_type)
    }

    /// Return an iterator over all barcodes and their total observed counts.
    /// Apply the provided filter on each individual count to determine if they
    /// participate in the sum.
    /// Barcodes that have zero counts after filtering will still be included
    /// in the output, with a count of zero.
    /// Note that this will not include barcodes present in the barcodes array
    /// but which have no counts even before filtering.
    fn filtered_barcode_counts<F: Fn(&AnnotatedCount<'_>) -> bool>(
        &self,
        filter: F,
    ) -> impl Iterator<Item = (&BarcodeWithGemGroup, i64)> {
        self.bcs_and_features
            .iter_barcode_ranges()
            .filter_map(move |(barcode, range)| {
                if range.is_empty() {
                    return None;
                }
                let filtered_count = range
                    .filter_map(|i| {
                        let count = self.counts[i];
                        if filter(&self.annotated_count(i, barcode)) {
                            Some(count as i64)
                        } else {
                            None
                        }
                    })
                    .sum();
                Some((barcode, filtered_count))
            })
    }

    /// Return annotated count information for the provided count index.
    fn annotated_count<'a>(
        &'a self,
        index: usize,
        barcode: &'a BarcodeWithGemGroup,
    ) -> AnnotatedCount<'a> {
        AnnotatedCount {
            count: self.counts[index],
            barcode,
            feature: &self.bcs_and_features.feature_reference.feature_defs
                [self.feature_indices[index] as usize],
        }
    }

    /// Return an iterator over all counts, including their barcode and feature.
    pub fn counts(&self) -> impl Iterator<Item = AnnotatedCount<'_>> {
        self.bcs_and_features
            .iter_barcode_ranges()
            .flat_map(|(barcode, range)| range.map(|i| self.annotated_count(i, barcode)))
    }

    pub fn raw_counts(&self) -> impl Iterator<Item = RawCount> + '_ {
        self.bcs_and_features
            .barcode_count_offsets
            .iter()
            .tuple_windows()
            .enumerate()
            .flat_map(move |(barcode_idx, (&start, &end))| {
                (start as usize..end as usize).map(move |index| RawCount {
                    count: self.counts[index],
                    barcode_idx,
                    feature_idx: self.feature_indices[index] as usize,
                })
            })
    }
}

impl EagerFeatures for CountMatrixEager {
    fn bcs_and_features(&self) -> &BarcodesAndFeatures {
        &self.bcs_and_features
    }
}

impl EagerFeatures for CountMatrixStreaming {
    fn bcs_and_features(&self) -> &BarcodesAndFeatures {
        &self.bcs_and_features
    }
}

impl CountMatrix for CountMatrixEager {}

trait EagerFeatures {
    fn bcs_and_features(&self) -> &BarcodesAndFeatures;
}

/// Shared API provided by all count matrix implementations.
// TODO: we might want to pull all of the common methods in here (at least the
// fallible versions) so that eager vs streaming loading are drop-in replacements
// for each other.
#[allow(private_bounds)] // Nothing outside this crate will implement this trait.
pub trait CountMatrix: EagerFeatures {
    /// Return a slice of all barcodes, including GEM group.
    /// This will include all barcodes in the matrix, even if they have no counts.
    fn barcodes(&self) -> &[BarcodeWithGemGroup] {
        &self.bcs_and_features().barcodes
    }

    fn num_barcodes(&self) -> usize {
        self.barcodes().len()
    }

    fn num_features(&self) -> usize {
        self.feature_reference().num_features()
    }

    /// Return the feature reference used by this count matrix.
    fn feature_reference(&self) -> &FeatureReference {
        &self.bcs_and_features().feature_reference
    }
}

#[derive(PartialEq, Eq, Debug)]
pub struct RawCount {
    pub count: Count,
    pub barcode_idx: usize,
    pub feature_idx: usize,
}

/// All data for a feature-barcode count.
#[derive(PartialEq, Eq, Debug)]
pub struct AnnotatedCount<'a> {
    pub count: Count,
    pub barcode: &'a BarcodeWithGemGroup,
    pub feature: &'a FeatureDef,
}

/// Write a sorted iterator of barcodes (with GEM group) to an HDF5 file.
pub fn write_barcodes_column(
    hdf5_group: &Group,
    column_name: &str,
    sorted_barcodes: impl Iterator<Item = Barcode>,
) -> Result<()> {
    // Buffer for writing out the fixed-width string representation of each barcode.
    let mut buf = Vec::with_capacity(BarcodeWithGemGroup::capacity());
    let formatted_bc_iter = sorted_barcodes.map(|bc| {
        write!(&mut buf, "{bc}")?;
        let formatted = BarcodeWithGemGroup::from_ascii(&buf)?;
        buf.clear();
        anyhow::Ok(formatted)
    });
    formatted_bc_iter.process_results(|bc_iter| {
        ChunkedWriter::write_all(
            hdf5_group,
            column_name,
            MAT_H5_BUF_SZ,
            ColumnAction::CreateNew,
            bc_iter,
        )
    })?
}

#[allow(clippy::too_many_arguments)]
/// Write a feature-barcode matrix h5 file to disk
/// Requires a loaded FeatureReference and list of barcodes to include.
/// Needs a shardio of FeatureCount values in natural order.
pub fn write_matrix_h5(
    path: &H5File,
    chunks: &[CountShardFile],
    feature_ref: &FeatureReference,
    sample_id: &str,
    chemistry_description: &str,
    gem_well: GemWell,
    barcode_index: &BarcodeIndex,
    software_version: &str,
) -> Result<()> {
    let reader: ShardReader<FeatureBarcodeCount, BarcodeThenFeatureOrder> =
        ShardReader::open_set(chunks)?;

    let unique_gem_groups = [gem_well.inner() as usize];

    let f = hdf5::File::create(path)?;
    let mut group = f.create_group(MATRIX_GROUP)?;
    let mut feature_ref_group = group.create_group(FEATURE_REF_GROUP)?;
    write_matrix_h5_helper(&mut group, &reader, feature_ref, barcode_index)?;
    feature_reference_io::to_h5(feature_ref, &mut feature_ref_group)?;

    // Add h5 metadata attributes

    scalar_attribute(
        &f,
        H5_FILETYPE_KEY,
        VarLenUnicode::from_str(MATRIX_H5_FILETYPE)?,
    )?;

    scalar_attribute(&f, MATRIX_H5_VERSION_KEY, MATRIX_H5_VERSION)?;

    scalar_attribute(
        &f,
        SOFTWARE_H5_VERSION_KEY,
        VarLenUnicode::from_str(software_version)?,
    )?;

    scalar_attribute(
        &f,
        H5_CHEMISTRY_DESC_KEY,
        VarLenAscii::from_ascii(chemistry_description)?,
    )?;

    let library_id_mapping = ndarray::Array1::from_elem(
        unique_gem_groups.len(),
        FixedAscii::<256>::from_ascii(sample_id)?,
    );
    let attr = f
        .new_attr::<FixedAscii<256>>()
        .shape(library_id_mapping.shape())
        .create(H5_LIBRARY_ID_MAPPING_KEY)?;
    attr.as_writer().write(library_id_mapping.view())?;

    let orig_gem_group = ndarray::Array1::from(
        unique_gem_groups
            .iter()
            .map(|&x| x as isize)
            .collect::<Vec<_>>(),
    );
    let attr = f
        .new_attr::<isize>()
        .shape(orig_gem_group.shape())
        .create(H5_ORIG_GEM_GROUP_MAPPING_KEY)?;
    attr.as_writer().write(orig_gem_group.view())?;

    f.close()?;
    Ok(())
}

/// Save UMI count data to the Cell Ranger HDF5 sparse matrix format. `group` must be the h5 group named 'matrix' where the data will be deposited.
/// /matrix
///     /data (int32 count values)
///     /indices (int32 feature ids)
///     /indptr (int64 indexeinto data for each barcode)
fn write_matrix_h5_helper(
    group: &mut Group,
    reader: &ShardReader<FeatureBarcodeCount, BarcodeThenFeatureOrder>,
    feature_ref: &FeatureReference,
    barcode_index: &BarcodeIndex,
) -> Result<()> {
    // Avoid failure in writing a hdf5 dataset of length 0 with gzip chunk size 1
    if barcode_index.is_empty() {
        bail!(
            "No 10x barcodes were observed in the experiment. This is likely the consequence of a \
            sample mixup or very poor sequencing quality on the barcode bases. Further execution \
            is halted."
        );
    }
    write_barcodes_column(group, DATASET_BARCODES, barcode_index.sorted_barcodes())?;

    let data = make_column_ds::<Count>(group, DATASET_DATA)?;
    let indices = make_column_ds::<BarcodeCountOffset>(group, DATASET_FEATURE_INDICES)?;

    let mut data_buf = Vec::with_capacity(MAT_H5_BUF_SZ);
    let mut indices_buf = Vec::with_capacity(MAT_H5_BUF_SZ);
    let mut barcode_counts = vec![0; barcode_index.len()];

    reader.iter()?.process_results(|iter| {
        for (barcode, bc_counts) in &iter.chunk_by(|x| x.barcode) {
            let mut n = 0i64;
            for (feature_idx, feat_counts) in &bc_counts.chunk_by(|x| x.feature_idx) {
                let umi_count: u32 = feat_counts.into_iter().map(|x| x.umi_count).sum();
                data_buf.push(umi_count as Count);
                indices_buf.push(feature_idx as FeatureIdx);
                n += 1;
            }
            let barcode_index = barcode_index.must_get(&barcode);
            assert_eq!(barcode_counts[barcode_index], 0);
            barcode_counts[barcode_index] = n;

            if data_buf.len() > MAT_H5_BUF_SZ * 3 / 4 {
                extend_dataset::<Count>(&data, &data_buf)?;
                data_buf.clear();

                extend_dataset(&indices, &indices_buf)?;
                indices_buf.clear();
            }
        }
        extend_dataset(&data, &data_buf)?;
        extend_dataset(&indices, &indices_buf)?;

        let mut total = 0i64;
        let mut indptr: Vec<BarcodeCountOffset> = vec![0];
        indptr.extend(barcode_counts.iter().map(|x| {
            total += x;
            total
        }));
        write_column_ds(group, DATASET_BARCODE_COUNT_OFFSETS, &indptr)?;

        write_column_ds(
            group,
            MATRIX_SHAPE_DATASET,
            &[
                feature_ref.feature_defs.len() as i32,
                barcode_counts.len() as i32,
            ],
        )?;
        anyhow::Ok(())
    })??;
    Ok(())
}

#[cfg(test)]
mod test {
    use super::*;
    use martian::MartianFileType;
    use std::collections::HashSet;
    use std::path::PathBuf;

    fn matrix_path(name: &str) -> PathBuf {
        PathBuf::from("../cr_h5/test/h5").join(name)
    }

    fn nonzero_bcs(data: &BarcodesAndFeatures) -> Vec<BarcodeWithGemGroup> {
        data.barcode_count_offsets
            .iter()
            .tuple_windows()
            .enumerate()
            .filter_map(|(bc_idx, (start, end))| (end > start).then_some(data.barcodes[bc_idx]))
            .collect()
    }

    fn bcs_from_counts(counts: &[(&BarcodeWithGemGroup, i64)]) -> Vec<BarcodeWithGemGroup> {
        counts.iter().map(|(bc, _)| **bc).collect()
    }

    #[test]
    fn test_count_matrix() -> Result<()> {
        let mat_a = CountMatrixFile::from_path(&matrix_path("diff_test_fbm_a.h5"));
        let mat = mat_a.read()?;
        let mat_stream = mat_a.read_streaming()?;
        assert_eq!(
            mat.bcs_and_features.barcodes.len() + 1,
            mat.bcs_and_features.barcode_count_offsets.len()
        );

        // Ensure that our counts iterator skips barcodes with zero-size ranges.
        let all_counts: Vec<_> = mat.counts().collect();
        assert_eq!(mat.counts.len(), all_counts.len());
        let all_counts_stream = mat_stream.counts()?.map(|c| c.unwrap()).collect_vec();
        assert_eq!(all_counts, all_counts_stream,);

        let feature_types: HashSet<_> = mat.counts().map(|c| c.feature.feature_type).collect();
        assert_eq!(1, feature_types.len());
        assert_eq!(FeatureType::Gene, *feature_types.iter().next().unwrap());

        // Check that we iterate over barcodes in the expected order.

        let nonzero_bcs: Vec<_> = nonzero_bcs(&mat.bcs_and_features);

        let bc_counts: Vec<_> = mat.filtered_barcode_counts(|_| true).collect();
        assert_eq!(nonzero_bcs, bcs_from_counts(&bc_counts));

        let bc_counts_stream: Vec<_> = mat_stream
            .filtered_barcode_counts(|_| true)?
            .map(Result::unwrap)
            .collect();
        assert_eq!(bc_counts, bc_counts_stream);

        // Check that filtering by feature type works.
        let bcs_gene_feats: Vec<_> = bcs_from_counts(
            &mat.barcode_counts_for_feature_type(FeatureType::Gene)
                .collect_vec(),
        );
        assert!(!bcs_gene_feats.is_empty());

        let bcs_ab_feats: Vec<_> = bcs_from_counts(
            &mat.barcode_counts_for_feature_type(FeatureType::Barcode(
                cr_types::FeatureBarcodeType::Antibody,
            ))
            .collect_vec(),
        );
        assert!(!bcs_ab_feats.is_empty());

        // The gene and antibody bcs together should be all of the nonzero barcodes.
        let all_feat_bcs: Vec<_> = bcs_gene_feats
            .iter()
            .chain(&bcs_ab_feats)
            .copied()
            .unique()
            .collect();
        assert_eq!(nonzero_bcs, all_feat_bcs);

        let bcs_gene_feats_stream: Vec<_> = bcs_from_counts(
            &mat_stream
                .barcode_counts_for_feature_type(FeatureType::Gene)?
                .map(Result::unwrap)
                .collect_vec(),
        );
        assert_eq!(bcs_gene_feats, bcs_gene_feats_stream);

        // Test barcode count filtering. Pick an arbitrarty feature.
        let feature_index = mat.feature_indices[0] as usize;
        let feature = &mat.feature_reference().feature_defs[feature_index];
        let filtered_counts: Vec<_> = mat
            .filtered_barcode_counts(|c| c.feature == feature)
            .collect();
        // Ensure we include all barcodes even if zero counts after filtering.
        assert_eq!(bcs_gene_feats, bcs_from_counts(&filtered_counts));

        let filtered_counts_stream: Vec<_> = mat_stream
            .filtered_barcode_counts(|c| c.feature == feature)?
            .map(Result::unwrap)
            .collect();
        assert_eq!(filtered_counts, filtered_counts_stream);

        // Regression test - check filtering by non-zero counts gives a sensible
        // answer.
        let expected_barcode_count = mat
            .feature_indices
            .iter()
            .filter(|i| **i as usize == feature_index)
            .count();
        assert!(expected_barcode_count > 0);
        assert_eq!(
            expected_barcode_count,
            filtered_counts
                .into_iter()
                .filter(|(_bc, count)| *count > 0)
                .count()
        );
        Ok(())
    }

    #[test]
    fn test_load_dimensions() -> Result<()> {
        let mat = CountMatrixFile::from_path(&matrix_path("diff_test_fbm_a.h5"));
        let MatrixDimensions {
            num_features,
            num_barcodes,
            num_non_zeros,
        } = mat.load_dimensions()?;
        assert_eq!(num_features, 33538);
        assert_eq!(num_barcodes, 737280);
        assert_eq!(num_non_zeros, 108617);
        Ok(())
    }
}
