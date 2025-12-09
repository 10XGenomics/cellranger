#![expect(missing_docs)]
///
/// Types for dealing with 10x barcodes
///
use crate::LibraryType;
use anyhow::{Context, Result, format_err};
use barcode::Barcode;
use itertools::Itertools;
use martian::{AsMartianBlanketType, AsMartianPrimaryType, MartianBlanketType};
use martian_derive::{MartianStruct, martian_filetype};
use martian_filetypes::bin_file::BinaryFormat;
use martian_filetypes::{LazyFileTypeIO, LazyWrite};
use metric::{Histogram, OrderedHistogram, TxHashMap, TxHashSet};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::fs::{hard_link, rename};
use std::iter::Peekable;
use std::path::{Path, PathBuf};

martian_filetype!(SortedBarcodeCount, "sbc");
/// A vector of barcode counts, written in barcode-sorted order.
/// The collection can be streamed item by item.
pub type SortedBarcodeCountFile = BinaryFormat<SortedBarcodeCount, Vec<(Barcode, usize)>>;

/// Write a sorted count file from an iterator of barcode, count pairs.
///
/// The caller must ensure that the iterator is already sorted.
pub fn write_sorted_count_file(
    items: impl Iterator<Item = (Barcode, usize)>,
    path: &Path,
) -> Result<SortedBarcodeCountFile> {
    let file: SortedBarcodeCountFile = path.into();
    let mut writer = file.lazy_writer()?;
    for item in items {
        writer.write_item(&item)?;
    }
    writer.finish()?;
    Ok(file)
}

/// Write a sorted count file from an OrderedHistogram.
pub fn write_sorted_count_file_from_histogram(
    hist: &OrderedHistogram<Barcode>,
    path: &Path,
) -> Result<SortedBarcodeCountFile> {
    write_sorted_count_file(
        hist.iter().map(|(bc, count)| (*bc, count.count() as usize)),
        path,
    )
}

/// A map of library type to sorted barcode count file for that library type.
///
/// The purpose of this type is to permit streaming all of the barcode counts
/// in sorted order in constant memory, as well as merging two or more collections
/// of counts in a constant-memory stream.
#[derive(Clone, Serialize, Deserialize)]
pub struct PerLibrarySortedBarcodeCounts(TxHashMap<LibraryType, SortedBarcodeCountFile>);

// Need to manually implement this since derive doesn't handle newtype structs.
impl AsMartianBlanketType for PerLibrarySortedBarcodeCounts {
    fn as_martian_blanket_type() -> MartianBlanketType {
        MartianBlanketType::TypedMap(Box::new(SortedBarcodeCountFile::as_martian_blanket_type()))
    }
}

impl PerLibrarySortedBarcodeCounts {
    /// Return an iterator of sorted barcode counts for the provided library type.
    ///
    /// Returns Err if the library type isn't in the collection.
    pub fn iter_for_library(
        &self,
        library_type: LibraryType,
    ) -> Result<impl Iterator<Item = Result<(Barcode, usize)>> + use<>> {
        self.0
            .get(&library_type)
            .ok_or_else(|| format_err!("sorted barcode counts not present for {library_type}"))?
            .lazy_reader()
    }

    /// Return iterators of sorted barcode counts for all library types.
    pub fn iter_each_lib(
        &self,
    ) -> Result<TxHashMap<LibraryType, impl Iterator<Item = Result<(Barcode, usize)>> + use<>>>
    {
        self.0
            .iter()
            .map(|(lib_type, f)| anyhow::Ok((*lib_type, f.lazy_reader()?)))
            .collect()
    }

    /// Return an iterator of sorted barcodes from all library types.
    /// The barcodes will be deduplicated.
    pub fn iter_barcodes(&self) -> Result<impl Iterator<Item = Result<Barcode>> + use<>> {
        Ok(merge_sorted_readers(self.readers()?)
            .map(|result| result.map(|(barcode, _count)| barcode))
            .dedup_by(dedup_results))
    }

    /// Return an iterator that merges all barcode counts from all library types.
    /// The counts will come in barcode-sorted order.
    pub fn iter_counts(&self) -> Result<MergeSortedCounts> {
        Ok(MergeSortedCounts::from_readers(self.readers()?))
    }

    /// Return an iterator over the library types.
    pub fn library_types(&self) -> impl Iterator<Item = LibraryType> + '_ {
        self.0.keys().copied()
    }

    /// Open readers for all files, in a consistent order.
    fn readers(&self) -> Result<Vec<impl Iterator<Item = Result<(Barcode, usize)>> + use<>>> {
        Ok(self
            .iter_each_lib()?
            .into_iter()
            .sorted_by_key(|(lib_type, _)| *lib_type)
            .map(|(_, r)| r)
            .collect())
    }

    /// Write a collection of ordered count histograms to disk.
    ///
    /// Call the provided make_path closure to create paths for each file.
    pub fn write_histograms<P>(
        histograms: &TxHashMap<LibraryType, OrderedHistogram<Barcode>>,
        make_path: P,
    ) -> Result<Self>
    where
        P: Fn(LibraryType) -> PathBuf,
    {
        let mut files = TxHashMap::default();
        for (library_type, hist) in histograms {
            let file = write_sorted_count_file_from_histogram(hist, &make_path(*library_type))?;
            files.insert(*library_type, file);
        }
        Ok(Self(files))
    }

    /// Merge a collection of per-lib input counts into a single output collection.
    ///
    /// This operation may move the input files directly into output paths
    /// for any library type that only has one input file, avoiding the need to
    /// stream anything.
    ///
    /// This is streamed directly from/to disk in constant memory.
    /// Call the provided make_path closure to create paths for each file.
    pub fn merge<P>(per_libs: impl Iterator<Item = Self>, make_path: P) -> Result<Self>
    where
        P: Fn(LibraryType) -> PathBuf,
    {
        Self::merge_from_iter(
            per_libs.flat_map(|per_lib| per_lib.0.into_iter()),
            make_path,
        )
    }

    /// Merge an iterator of per-lib input counts into a single output collection.
    ///
    /// This operation may move the input files directly into output paths
    /// for any library type that only has one input file, avoiding the need to
    /// stream anything.
    ///
    /// This is streamed directly from/to disk in constant memory.
    /// Call the provided make_path closure to create paths for each file.
    pub fn merge_from_iter<P>(
        collections: impl Iterator<Item = (LibraryType, SortedBarcodeCountFile)>,
        make_path: P,
    ) -> Result<Self>
    where
        P: Fn(LibraryType) -> PathBuf,
    {
        // Collect all of the input files for each library type.
        let files_by_lib_type = collections.into_group_map();

        let mut output_files = TxHashMap::default();
        for (library_type, input_files) in files_by_lib_type {
            let path: SortedBarcodeCountFile = make_path(library_type).into();
            // Optimization - if we only have one input file, there's no neeed to merge.
            // Move the file into the output path.
            let output_file = if input_files.len() == 1 {
                rename(&input_files[0], &path).with_context(|| {
                    format!(
                        "failed to move sorted counts file at path {} to path {}",
                        input_files[0].to_string_lossy(),
                        path.to_string_lossy()
                    )
                })?;
                path
            } else {
                // Merge-stream the input files.
                MergeSortedCounts::from_readers(
                    input_files
                        .into_iter()
                        .map(|f| f.lazy_reader())
                        .try_collect()?,
                )
                .process_results(|counts_iter| write_sorted_count_file(counts_iter, &path))??
            };
            output_files.insert(library_type, output_file);
        }
        Ok(Self(output_files))
    }

    pub fn files(&self) -> &TxHashMap<LibraryType, SortedBarcodeCountFile> {
        &self.0
    }
}

fn merge_sorted_readers<R>(readers: Vec<R>) -> impl Iterator<Item = Result<(Barcode, usize)>>
where
    R: Iterator<Item = Result<(Barcode, usize)>> + 'static,
{
    readers
        .into_iter()
        .kmerge_by(merge_results_by_key(|(bc, _count)| *bc))
}

/// Return a closure suitable to pass to Itertools::kmerge_by that handles merging
/// a stream of Result<T> by:
///     - always sorting Err to the front of the line
///     - uses the key projection function F to get a key value for each item
///         to use for sorting. The key function should be cheap to call.
pub fn merge_results_by_key<F, T, K, E>(
    get_key: F,
) -> impl FnMut(&Result<T, E>, &Result<T, E>) -> bool
where
    F: Fn(&T) -> K,
    K: Ord,
{
    move |a, b| match (a, b) {
        (Err(_), _) => true,
        (_, Err(_)) => false,
        (Ok(a), Ok(b)) => get_key(a) < get_key(b),
    }
}

/// A function to pass to Itertools::dedup_by that handles deduping
/// a stream of Result<T> by never deduping Err.
///
/// The items in each Ok are deduped by direct comparison.
pub fn dedup_results<T: Eq, E>(a: &Result<T, E>, b: &Result<T, E>) -> bool {
    match (a, b) {
        (Err(_), _) | (_, Err(_)) => false,
        (Ok(bc0), Ok(bc1)) => bc0 == bc1,
    }
}

/// An iterator that merges multiple fallible streams of sorted barcode counts.
pub struct MergeSortedCounts {
    #[allow(clippy::type_complexity)]
    merged_reader: Peekable<Box<dyn Iterator<Item = Result<(Barcode, usize)>>>>,
}

impl MergeSortedCounts {
    /// Create a count merger from a collection of individually sorted readers.
    pub fn from_readers<T>(readers: Vec<T>) -> Self
    where
        T: Iterator<Item = Result<(Barcode, usize)>> + 'static,
    {
        Self {
            merged_reader: (Box::new(merge_sorted_readers(readers))
                as Box<dyn Iterator<Item = Result<(Barcode, usize)>>>)
                .peekable(),
        }
    }

    /// Create a count merger from a single sorted reader.
    pub fn new<T>(sorted_reader: T) -> Self
    where
        T: Iterator<Item = Result<(Barcode, usize)>> + 'static,
    {
        Self {
            merged_reader: (Box::new(sorted_reader)
                as Box<dyn Iterator<Item = Result<(Barcode, usize)>>>)
                .peekable(),
        }
    }
}

// It would be nice to use Itertools::chunk_by but the fact that we get a stream
// of Result makes this rather unwieldly without consuming the entire iterator
// using process_results. This is essentially a hand-rolled version of that
// operation, but that handles Result in a straightforward fashion, assuming that
// we want to terminate promptly on any error.
impl Iterator for MergeSortedCounts {
    type Item = Result<(Barcode, usize)>;

    fn next(&mut self) -> Option<Self::Item> {
        // Get the first item in the next group.
        let (barcode, mut count) = match self.merged_reader.next() {
            None => {
                // End of iteration.
                return None;
            }
            Some(Err(err)) => {
                // Encountered an error.
                return Some(Err(err));
            }
            Some(Ok(item)) => item,
        };

        // Consume items as long as they match the current barcode, or are errors.
        while let Some(next_result) = self.merged_reader.next_if(|next_result| match next_result {
            Ok((next_barcode, _)) => *next_barcode == barcode,
            Err(_) => true,
        }) {
            match next_result {
                Err(err) => {
                    return Some(Err(err));
                }
                Ok((_, next_count)) => {
                    count += next_count;
                }
            }
        }

        Some(Ok((barcode, count)))
    }
}

martian_filetype!(BarcodeIndexFile, "bi");
type _BarcodeIndexFormat = BinaryFormat<BarcodeIndexFile, Vec<Barcode>>;

/// On-disk representation of the BarcodeIndex.
///
/// By storing as a streamable sorted vec, we can write the index in constant
/// memory from a pre-sorted collection, and we can also stream the index if
/// we don't need the by-barcode lookup.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct BarcodeIndexFormat(_BarcodeIndexFormat);

impl BarcodeIndexFormat {
    /// Write the barcode index from a sorted iterator of barcodes.
    ///
    /// Return the path, and the total number of barcodes.
    pub fn write(
        path: &Path,
        sorted_bc_iter: impl Iterator<Item = Barcode>,
    ) -> Result<BarcodeIndexOutput> {
        let file: BinaryFormat<BarcodeIndexFile, Vec<Barcode>> = path.into();
        let mut writer = file.lazy_writer()?;
        // FIXME: the writer already tracks count privately
        let mut count = 0;
        for bc in sorted_bc_iter {
            writer.write_item(&bc)?;
            count += 1;
        }
        writer.finish()?;
        Ok(BarcodeIndexOutput {
            index: Self(file),
            num_barcodes: count,
        })
    }

    /// Read the barcode index.
    pub fn read(&self) -> Result<BarcodeIndex> {
        self.0
            .lazy_reader()?
            .process_results(|bc_iter| BarcodeIndex::from_sorted(bc_iter))
    }

    /// Stream the sorted barcodes from disk item by item.
    pub fn lazy_reader(&self) -> Result<impl Iterator<Item = Result<Barcode>> + use<>> {
        self.0.lazy_reader()
    }
}

impl AsMartianPrimaryType for BarcodeIndexFormat {
    fn as_martian_primary_type() -> martian::MartianPrimaryType {
        _BarcodeIndexFormat::as_martian_primary_type()
    }
}

impl AsRef<Path> for BarcodeIndexFormat {
    fn as_ref(&self) -> &Path {
        &self.0
    }
}

/// A data structure to assist with looking up the numeric index for a barcode.
///
/// The index is created by enumerating all sorted barcodes.
#[derive(Debug, Clone)]
pub struct BarcodeIndex {
    index: BTreeMap<Barcode, usize>,
}

impl BarcodeIndex {
    /// Construct a barcode index from a unique, sorted iterator of barcodes.
    ///
    /// The caller must ensure that the iterator is producing sorted, deduped barcodes.
    pub fn from_sorted(bc_iter: impl IntoIterator<Item = Barcode>) -> Self {
        Self {
            index: bc_iter
                .into_iter()
                .enumerate()
                .map(|(index, bc)| (bc, index))
                .collect(),
        }
    }

    /// Return the numerical index of the specified barcode.
    pub fn get(&self, barcode: &Barcode) -> Option<usize> {
        self.index.get(barcode).copied()
    }

    /// Return the numerical index of the specified barcode.
    ///
    /// Panic if the barcode is not in the index.
    pub fn must_get(&self, barcode: &Barcode) -> usize {
        self.get(barcode)
            .ok_or_else(|| format!("barcode {barcode} not in index"))
            .unwrap()
    }

    /// Return an iterator of the sorted barcodes in index order.
    pub fn sorted_barcodes(&self) -> impl Iterator<Item = Barcode> + '_ {
        self.index.keys().copied()
    }

    /// Return an indicator vector corresponding to a set of barcodes
    /// A vector of bools returned where index i is true if
    /// sorted_barcodes[i] is in the set. Else false.
    /// If the set contains barcodes not in the BarcodeIndex, this is
    /// not caught
    pub fn into_indicator_vec(&self, filter_set: &TxHashSet<Barcode>) -> Vec<bool> {
        self.index.keys().map(|x| filter_set.contains(x)).collect()
    }

    /// Return whether this index contains no barcodes.
    pub fn is_empty(&self) -> bool {
        self.index.is_empty()
    }

    /// Return the number of barcodes.
    pub fn len(&self) -> usize {
        self.index.len()
    }
}

/// Data structure returned after construction of BarcodeIndex.
/// Includes the total number of barcodes in the index.
#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct BarcodeIndexOutput {
    pub index: BarcodeIndexFormat,
    pub num_barcodes: usize,
}

impl BarcodeIndexOutput {
    /// Hardlink or copy this output into a new path.
    pub fn hardlink_or_copy(&self, path: &Path) -> Result<Self> {
        let path: _BarcodeIndexFormat = path.into();
        // Try hard linking; if fails, fall back to copy.
        if hard_link(&self.index, &path).is_err() {
            std::fs::copy(&self.index, &path)?;
        }
        Ok(Self {
            index: BarcodeIndexFormat(path),
            num_barcodes: self.num_barcodes,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use barcode::BcSeq;

    fn bc(seq: &[u8]) -> Barcode {
        Barcode::with_seq(1, BcSeq::from_bytes(seq), true)
    }

    #[test]
    fn test_into_indicator_vec() {
        let vec_of_barcodes = vec![bc(b"A"), bc(b"G"), bc(b"T")];

        let query_vec_of_barcodes = vec![bc(b"G"), bc(b"A"), bc(b"C")];
        let b_index = BarcodeIndex::from_sorted(vec_of_barcodes);
        let set_of_barcodes = TxHashSet::from_iter(query_vec_of_barcodes);
        assert_eq!(
            b_index.into_indicator_vec(&set_of_barcodes),
            vec![true, true, false]
        );
    }

    fn make_hist(count_pairs: &[(Barcode, usize)]) -> OrderedHistogram<Barcode> {
        let mut hist = OrderedHistogram::default();
        for (bc, count) in count_pairs {
            hist.observe_by(bc, *count as i64);
        }
        hist
    }

    /// Test all key functions of sorted barcode counts.
    #[test]
    fn test_per_library_sorted_barcode_counts() -> Result<()> {
        let counts0 = [(bc(b"A"), 1), (bc(b"C"), 3), (bc(b"G"), 4)];
        let counts1 = [(bc(b"G"), 1), (bc(b"T"), 7)];
        let counts2 = [(bc(b"AA"), 11), (bc(b"CC"), 13)];

        let hists: TxHashMap<_, _> = [
            (LibraryType::Gex, make_hist(&counts0)),
            (LibraryType::Antibody, make_hist(&counts1)),
            (LibraryType::Crispr, make_hist(&counts2)),
        ]
        .into_iter()
        .collect();

        let dir = tempfile::tempdir()?;

        let per_lib = PerLibrarySortedBarcodeCounts::write_histograms(
            &hists,
            path_maker(dir.path(), "hists"),
        )?;

        let assert_lib_eq = |lib_type, counts| {
            let pulled_counts: Vec<_> = per_lib
                .iter_for_library(lib_type)
                .unwrap()
                .try_collect()
                .unwrap();
            assert_eq!(counts, pulled_counts.as_slice());
        };

        assert_lib_eq(LibraryType::GeneExpression, counts0.as_slice());
        assert_lib_eq(LibraryType::Antibody, counts1.as_slice());
        assert_lib_eq(LibraryType::Crispr, counts2.as_slice());

        // Ensure we merge correctly.
        let counts: Vec<_> = per_lib.iter_counts()?.try_collect()?;

        let expected_counts = vec![
            (bc(b"A"), 1),
            (bc(b"AA"), 11),
            (bc(b"C"), 3),
            (bc(b"CC"), 13),
            (bc(b"G"), 5),
            (bc(b"T"), 7),
        ];

        assert_eq!(counts, expected_counts);

        let barcodes: Vec<_> = per_lib.iter_barcodes()?.try_collect()?;
        let expected_barcodes = vec![bc(b"A"), bc(b"AA"), bc(b"C"), bc(b"CC"), bc(b"G"), bc(b"T")];
        assert_eq!(barcodes, expected_barcodes);

        // Merge the counts collection with itself, ensure we get double the counts.
        let merged = PerLibrarySortedBarcodeCounts::merge(
            [per_lib.clone(), per_lib.clone()].into_iter(),
            path_maker(dir.path(), "merged"),
        )?;

        let counts_merged: Vec<_> = merged.iter_counts()?.try_collect()?;
        let expected_counts_merged: Vec<_> = expected_counts
            .iter()
            .copied()
            .map(|(bc, count)| (bc, count * 2))
            .collect();
        assert_eq!(counts_merged, expected_counts_merged);

        // Test the "move file" optimization where a merge will move a file
        // instead of stream-merging it if a library type only has a single file.
        // Do this by merging an empty collection with our original.
        for file in per_lib.0.values() {
            assert!(std::fs::exists(file)?);
        }
        let moved = PerLibrarySortedBarcodeCounts::merge(
            [
                per_lib.clone(),
                PerLibrarySortedBarcodeCounts(Default::default()),
            ]
            .into_iter(),
            path_maker(dir.path(), "moved"),
        )?;
        for file in per_lib.0.values() {
            assert!(!std::fs::exists(file)?);
        }

        let counts_moved: Vec<_> = moved.iter_counts()?.try_collect()?;
        assert_eq!(counts_moved, expected_counts);
        Ok(())
    }

    fn path_maker(dir: &Path, ctx: &str) -> impl Fn(LibraryType) -> PathBuf + use<> {
        let dir = dir.to_path_buf();
        let ctx = ctx.to_string();
        move |lib_type: LibraryType| -> PathBuf {
            let mut path = PathBuf::new();
            path.push(&dir);
            path.push(format!("{lib_type}_{ctx}"));
            path
        }
    }
}
