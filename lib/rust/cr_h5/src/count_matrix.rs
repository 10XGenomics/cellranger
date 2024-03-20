//! A feature-barcode count matrix.

use crate::{
    extend_dataset, feature_reference_io, make_column_ds, scalar_attribute, write_column_ds,
};
use anyhow::{bail, Context, Result};
use barcode::MAX_BARCODE_AND_GEM_GROUP_LENGTH;
use cr_types::barcode_index::BarcodeIndex;
use cr_types::reference::feature_reference::{FeatureDef, FeatureReference, FeatureType};
use cr_types::{BarcodeThenFeatureOrder, CountShardFile, FeatureBarcodeCount, GemWell, H5File};
use hdf5::types::{FixedAscii, VarLenAscii, VarLenUnicode};
use hdf5::Group;
use itertools::{process_results, Itertools};
use martian_derive::martian_filetype;
use serde::{Deserialize, Serialize};
use shardio::ShardReader;
use std::cell::OnceCell;
use std::fmt::Display;
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
type Count = i32;
/// The index into the array of features for a particular count.
type FeatureIdx = i32;
/// A value representing the offset into the count array corresponding to
/// the beginning of the count range for that barcode.
type BarcodeCountOffset = i64;
/// A barcode, including the -N GEM well suffix.
pub type BarcodeWithGemGroup = FixedAscii<MAX_BARCODE_AND_GEM_GROUP_LENGTH>;

martian_filetype!(CountMatrixFile, "h5");

/// A feature-barcode count matrix HDF5 file.
impl CountMatrixFile {
    fn matrix_group(&self) -> Result<hdf5::Group> {
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

    // Return the amount of memory in GiB required for the matrix.
    pub fn estimate_mem_gib(&self) -> Result<f64> {
        Ok(self.load_dimensions()?.estimate_mem_gib())
    }
    /// Read all matrix data into memory.
    pub fn read(&self) -> Result<CountMatrix> {
        // Augment any read error with the path to the file.
        self.read_inner()
            .with_context(|| self.display().to_string())
    }

    pub fn read_barcodes(&self) -> Result<Vec<BarcodeWithGemGroup>> {
        let matrix = self.matrix_group()?;
        matrix
            .dataset(DATASET_BARCODES)?
            .read_raw()
            .context("Ran into trouble reading barcodes from feature bc matrix.")
    }

    fn read_inner(&self) -> Result<CountMatrix> {
        let matrix = self.matrix_group()?;
        Ok(CountMatrix {
            counts: matrix.dataset(DATASET_DATA)?.read_raw()?,
            barcodes: matrix.dataset(DATASET_BARCODES)?.read_raw()?,
            feature_indices: matrix.dataset(DATASET_FEATURE_INDICES)?.read_raw()?,
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
        size_of::<BarcodeWithGemGroup>() * self.num_barcodes + // barcodes
        size_of::<BarcodeCountOffset>() * (1 + self.num_barcodes) + // barcode_count_offsets (i64)
        FeatureReference::estimate_mem_bytes_for_feature_count(self.num_features)
    }

    // Return the amount of memory in GiB required for the matrix.
    pub fn estimate_mem_gib(&self) -> f64 {
        (self.estimate_mem_bytes() as f64) / (1024.0 * 1024.0 * 1024.0)
    }
}

/// In-memory representation of the count matrix, lazily loaded once.
/// Use this type if you need to conditionally load the matrix in multiple places.
pub struct LazyCountMatrix {
    file: CountMatrixFile,
    matrix: OnceCell<CountMatrix>,
}

impl LazyCountMatrix {
    pub fn new(file: CountMatrixFile) -> Self {
        Self {
            file,
            matrix: Default::default(),
        }
    }

    pub fn loaded(&self) -> Result<&CountMatrix> {
        // TODO: refactor once once_cell_try is stabilized.
        if let Some(mat) = self.matrix.get() {
            return Ok(mat);
        }
        let mat = self.file.read()?;
        assert!(self.matrix.set(mat).is_ok());
        Ok(self.matrix.get().unwrap())
    }
}

/// In-memory representation of the count matrix.
/// All data is eagerly loaded.
pub struct CountMatrix {
    counts: Vec<Count>,
    barcodes: Vec<BarcodeWithGemGroup>,
    feature_indices: Vec<FeatureIdx>,
    barcode_count_offsets: Vec<BarcodeCountOffset>,
    feature_reference: FeatureReference,
}

impl CountMatrix {
    /// Return a slice of all barcodes, including GEM group.
    /// This will include all barcodes in the matrix, even if they have no counts.
    pub fn barcodes(&self) -> &[BarcodeWithGemGroup] {
        &self.barcodes
    }

    pub fn num_barcodes(&self) -> usize {
        self.barcodes.len()
    }

    pub fn num_features(&self) -> usize {
        self.feature_reference.num_features()
    }

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
    pub fn filtered_barcode_counts<F: Fn(&AnnotatedCount<'_>) -> bool>(
        &self,
        filter: F,
    ) -> impl Iterator<Item = (&BarcodeWithGemGroup, i64)> {
        self.iter_barcode_ranges().map(move |(barcode, range)| {
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
            (barcode, filtered_count)
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
            feature: &self.feature_reference.feature_defs[self.feature_indices[index] as usize],
        }
    }

    /// Return an iterator over all barcodes and their range in the counts array.
    fn iter_barcode_ranges(&self) -> impl Iterator<Item = (&BarcodeWithGemGroup, Range<usize>)> {
        zip(
            &self.barcodes,
            self.barcode_count_offsets.iter().tuple_windows(),
        )
        .map(|(barcode, (&start, &end))| (barcode, (start as usize..end as usize)))
    }

    /// Return an iterator over all counts, including their barcode and feature.
    pub fn counts(&self) -> impl Iterator<Item = AnnotatedCount<'_>> {
        self.iter_barcode_ranges()
            .flat_map(|(barcode, range)| range.map(|i| self.annotated_count(i, barcode)))
    }

    /// Return the feature reference used by this count matrix.
    pub fn feature_reference(&self) -> &FeatureReference {
        &self.feature_reference
    }

    pub fn raw_counts(&self) -> impl Iterator<Item = RawCount> + '_ {
        self.barcode_count_offsets
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

pub struct RawCount {
    pub count: Count,
    pub barcode_idx: usize,
    pub feature_idx: usize,
}

/// All data for a feature-barcode count.
pub struct AnnotatedCount<'a> {
    pub count: Count,
    pub barcode: &'a BarcodeWithGemGroup,
    pub feature: &'a FeatureDef,
}

/// Write a column of barcode strings to an HDF5 file.
pub fn write_barcodes_column<B>(
    hdf5_group: &Group,
    column_name: &str,
    sorted_barcodes: &[B],
) -> Result<()>
where
    B: Display,
{
    let strings = sorted_barcodes
        .iter()
        .map(|barcode| FixedAscii::from_ascii(&barcode.to_string()).unwrap())
        .collect::<Vec<BarcodeWithGemGroup>>();
    hdf5_group
        .new_dataset::<BarcodeWithGemGroup>()
        .deflate(1)
        .shape((strings.len(),))
        .create(column_name)?
        .write(&strings)?;
    Ok(())
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

    process_results(reader.iter()?, |iter| {
        for (barcode, bc_counts) in &iter.group_by(|x| x.barcode) {
            let mut n = 0i64;
            for (feature_idx, feat_counts) in &bc_counts.group_by(|x| x.feature_idx) {
                let umi_count: u32 = feat_counts.into_iter().map(|x| x.umi_count).sum();
                data_buf.push(umi_count as Count);
                indices_buf.push(feature_idx as FeatureIdx);
                n += 1;
            }
            let barcode_index = barcode_index.get_index(&barcode);
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

    #[test]
    fn test_count_matrix() -> Result<()> {
        let mat_a = CountMatrixFile::from_path(&matrix_path("diff_test_fbm_a.h5"));
        let mat = mat_a.read()?;
        assert_eq!(mat.barcodes.len() + 1, mat.barcode_count_offsets.len());
        // Ensure that our counts iterator skips barcodes with zero-size ranges.
        assert_eq!(mat.counts.len(), mat.counts().count());

        let feature_types: HashSet<_> = mat.counts().map(|c| c.feature.feature_type).collect();
        assert_eq!(1, feature_types.len());
        assert_eq!(FeatureType::Gene, *feature_types.iter().next().unwrap());

        // Check that we iterate over barcodes in the expected order.
        let bcs_nonzero_counts: Vec<_> = mat
            .barcode_counts_for_feature_type(FeatureType::Gene)
            .map(|(bc, _count)| *bc)
            .collect();
        assert!(!bcs_nonzero_counts.is_empty());
        assert_eq!(
            bcs_nonzero_counts.as_slice(),
            &mat.barcodes()[..bcs_nonzero_counts.len()],
        );

        // Test barcode count filtering. Pick an arbitrarty feature.
        let feature_index = mat.feature_indices[0] as usize;
        let feature = &mat.feature_reference().feature_defs[feature_index];
        let filtered_counts: Vec<_> = mat
            .filtered_barcode_counts(|c| c.feature == feature)
            .collect();
        // Ensure we include all barcodes even if zero counts after filtering.
        assert_eq!(bcs_nonzero_counts.len(), filtered_counts.len());

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
