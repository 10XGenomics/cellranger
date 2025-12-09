#![expect(missing_docs)]
use anyhow::Result;
use martian::MartianFileType;
use martian_derive::martian_filetype;
use parquet::basic::{Compression, ZstdLevel};
use parquet::file::metadata::ParquetMetaData;
use parquet::file::properties::WriterProperties;
use parquet::file::reader::FileReader;
use parquet::file::serialized_reader::SerializedFileReader;
use parquet::file::writer::SerializedFileWriter;
use parquet::record::{RecordReader, RecordWriter, Row};
use parquet_derive::ParquetRecordWriter;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::BufWriter;

martian_filetype! {ParquetFile, "parquet"}

impl ParquetFile {
    /// Instantiate a writer with a row group size.
    ///
    /// The row group size should be chosen based on a balance between the amount of chunking and
    /// compression/buffering.
    /// From https://parquet.apache.org/docs/file-format/configurations/#row-group-size
    /// Larger row groups allow for larger column chunks which makes it possible to do larger
    /// sequential IO. Larger groups also require more buffering in the write path
    /// (or a two pass write). We recommend large row groups (512MB - 1GB).
    pub fn writer<T>(&self, row_group_size: usize) -> Result<ParquetWriter<T>>
    where
        for<'a> &'a [T]: RecordWriter<T>,
    {
        assert!(row_group_size > 0);
        // T::schema() would have been a better design in ParquetRecordWriter
        let empty: &[T] = &[];
        Ok(ParquetWriter {
            writer: Some(SerializedFileWriter::new(
                self.buf_writer()?,
                empty.schema()?,
                WriterProperties::builder()
                    .set_compression(Compression::ZSTD(ZstdLevel::try_new(3)?))
                    .build()
                    .into(),
            )?),
            row_group_size,
            cache: Vec::with_capacity(row_group_size),
        })
    }

    /// Instantiate a reader. It is often better to use polars or similar libraries for reading
    /// parquet files.
    ///
    /// This struct exists in case it is more convenient to operate on a vec of struct rather than
    /// a dataframe.
    pub fn reader<T>(&self) -> Result<ParquetReader<T>>
    where
        Vec<T>: RecordReader<T>,
    {
        let reader = SerializedFileReader::try_from(File::open(self.as_ref())?)?;
        Ok(ParquetReader::new(reader))
    }

    pub fn metadata(&self) -> Result<ParquetMetaData> {
        let reader = SerializedFileReader::try_from(File::open(self.as_ref())?)?;
        Ok(reader.metadata().clone())
    }
}

pub struct ParquetWriter<T>
where
    for<'a> &'a [T]: RecordWriter<T>,
{
    writer: Option<SerializedFileWriter<BufWriter<File>>>,
    row_group_size: usize,
    cache: Vec<T>,
}

impl<T> ParquetWriter<T>
where
    for<'a> &'a [T]: RecordWriter<T>,
{
    pub fn push(&mut self, item: T) -> Result<()> {
        self.cache.push(item);
        if self.cache.len() >= self.row_group_size {
            self.write_and_clear_cache()?;
        }
        Ok(())
    }

    pub fn push_all<I>(&mut self, items: I) -> Result<()>
    where
        I: IntoIterator<Item = T>,
    {
        for item in items {
            self.push(item)?;
        }
        Ok(())
    }

    fn write_and_clear_cache(&mut self) -> Result<()> {
        if !self.cache.is_empty() {
            Self::write_to_new_row_group(self.writer.as_mut().unwrap(), &self.cache)?;
            self.cache.clear();
        }
        Ok(())
    }

    fn write_to_new_row_group(
        writer: &mut SerializedFileWriter<BufWriter<File>>,
        items: &[T],
    ) -> Result<()> {
        let mut row_group = writer.next_row_group()?;
        items.write_to_row_group(&mut row_group)?;
        row_group.close()?;
        Ok(())
    }

    pub fn write_all(&mut self, items: &[T]) -> Result<()> {
        for item_chunk in items.chunks(self.row_group_size) {
            Self::write_to_new_row_group(self.writer.as_mut().unwrap(), item_chunk)?;
        }
        Ok(())
    }
    fn finish(&mut self) -> Result<()> {
        self.write_and_clear_cache()?;
        if let Some(writer) = self.writer.take() {
            writer.close()?;
        }
        Ok(())
    }
    pub fn close(mut self) -> Result<()> {
        self.finish()?;
        Ok(())
    }
}

impl<T> Drop for ParquetWriter<T>
where
    for<'a> &'a [T]: RecordWriter<T>,
{
    fn drop(&mut self) {
        self.finish().unwrap();
    }
}

pub struct ParquetReader<T>
where
    Vec<T>: RecordReader<T>,
{
    reader: SerializedFileReader<File>,
    next_row_group_idx: usize,
    num_row_groups: usize,
    cache: std::vec::IntoIter<T>,
}

impl<T> ParquetReader<T>
where
    Vec<T>: RecordReader<T>,
{
    pub fn new(reader: SerializedFileReader<File>) -> Self {
        let num_row_groups = reader.num_row_groups();
        ParquetReader {
            reader,
            next_row_group_idx: 0,
            num_row_groups,
            cache: Vec::new().into_iter(),
        }
    }

    fn read_next_row_group(&mut self) -> Result<()> {
        if self.next_row_group_idx >= self.num_row_groups {
            return Ok(());
        }
        let mut row_group = self.reader.get_row_group(self.next_row_group_idx)?;
        self.next_row_group_idx += 1;
        let num_rows = row_group.metadata().num_rows() as usize;
        let mut cache: Vec<T> = Vec::with_capacity(num_rows);
        cache.read_from_row_group(&mut *row_group, num_rows)?;
        self.cache = cache.into_iter();

        Ok(())
    }
}

impl<T> Iterator for ParquetReader<T>
where
    Vec<T>: RecordReader<T>,
{
    type Item = Result<T>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(item) = self.cache.next() {
            return Some(Ok(item));
        }
        if self.next_row_group_idx >= self.num_row_groups {
            return None;
        }
        match self.read_next_row_group() {
            Ok(()) => self.next(),
            Err(e) => Some(Err(e)),
        }
    }
}

/// A row in the per-read gap alignment table.
///
/// This table contains the gap sequence between two gap align probes for reads
/// that have been aligned to gap align probes.
/// Since ParquetRecordReader does not support Option<>, we use Deserialize
#[derive(ParquetRecordWriter, Deserialize)]
pub struct PerReadGapAlignRow {
    pub barcode: String,
    pub is_valid_barcode: bool,
    pub umi: String,
    pub left_probe_id: String,
    pub right_probe_id: String,
    pub probe_type: String,
    pub read_gap_sequence: String,
    pub read_gap_len: usize,
    pub expected_gap_sequence: Option<String>,
    pub gap_levenshtein_distance: Option<u32>,
    pub gap_within_max_error: Option<bool>,
    pub gap_exactly_matches_expected: Option<bool>,
    // Alignment (gap and expected gap) stats
    pub num_matches: Option<usize>,
    pub num_mismatches: Option<usize>,
    pub num_insertions: Option<usize>,
    pub num_deletions: Option<usize>,
    pub ends_with_insertion: Option<bool>,
    pub ends_with_deletion: Option<bool>,
    pub starts_with_insertion: Option<bool>,
    pub starts_with_deletion: Option<bool>,
}

impl PerReadGapAlignRow {
    // About 200bytes per row, 1M rows would mean that each row group size is about 200MB
    pub const ROW_GROUP_SIZE: usize = 1_000_000;

    pub fn from_row(row: Row) -> Self {
        serde_json::from_value(row.to_json_value()).unwrap()
    }
}

/// A row in the per-umi gap alignment table
///
/// This table contains the consensus gap sequence from based on reads from this UMI
/// and includes alignment statistics.
/// Since ParquetRecordReader does not support Option<>, we use Deserialize
#[derive(ParquetRecordWriter, Deserialize)]
pub struct PerUmiGapAlignRow {
    pub barcode: String,
    pub is_valid_barcode: bool,
    pub sample_id: Option<String>,
    pub umi: String,
    pub left_probe_id: String,
    pub right_probe_id: String,
    pub probe_type: String,
    pub umi_gap_sequence: String,
    pub umi_gap_len: usize,
    pub expected_gap_sequence: Option<String>,
    pub gap_levenshtein_distance: Option<u32>,
    pub gap_within_max_error: Option<bool>,
    pub gap_exactly_matches_expected: Option<bool>,
    // Alignment (consensun gap and expected gap) stats
    pub num_matches: Option<usize>,
    pub num_mismatches: Option<usize>,
    pub num_insertions: Option<usize>,
    pub num_deletions: Option<usize>,
    pub ends_with_insertion: Option<bool>,
    pub ends_with_deletion: Option<bool>,
    pub starts_with_insertion: Option<bool>,
    pub starts_with_deletion: Option<bool>,
    // Aggregate read level info
    pub num_gap_reads: usize, // num reads within max error of most common gap seq
    pub num_wasted_reads: usize,
    // Sum of read level alignment info
    pub num_reads_with_aligned_gap: usize,
    pub num_reads_with_correct_gap: usize,
    pub num_read_gaps_with_mimatches: usize,
    pub num_read_gaps_starting_with_insertion: usize,
    pub num_read_gaps_starting_with_deletion: usize,
    pub num_read_gaps_ending_with_insertion: usize,
    pub num_read_gaps_ending_with_deletion: usize,
}

impl PerUmiGapAlignRow {
    pub const ROW_GROUP_SIZE: usize = PerReadGapAlignRow::ROW_GROUP_SIZE;

    pub fn from_row(row: Row) -> Self {
        serde_json::from_value(row.to_json_value()).unwrap()
    }
}

/// A row in the per-bc vdj cell filter table
///
#[derive(ParquetRecordWriter, Deserialize)]
pub struct PerBarcodeFilter {
    pub barcode: String,
    pub is_cell: bool,
    pub is_gex_cell: Option<bool>,
    pub is_asm_cell: Option<bool>,
    pub low_umi: bool,
    pub no_v_region: bool,
    pub low_junction_support: bool,
    pub no_conf_contig: bool,
    pub low_rpu: bool,
    pub non_dominant_junction: bool,
    pub weak_junction: bool,
    pub chimeric: bool,
    pub common_clone: bool,
    pub gel_bead_contamination: bool,
    pub gel_bead_indel: Option<bool>,
    pub enclone_fate: Option<String>,
    pub insert_priming: bool,
}

impl PerBarcodeFilter {
    // About 60bytes per row, 1M rows would mean that each row group size is 60MB
    pub const ROW_GROUP_SIZE: usize = 1_000_000;

    pub fn from_row(row: Row) -> Self {
        serde_json::from_value(row.to_json_value()).unwrap()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use itertools::Itertools;
    use martian::MartianTempFile;
    use parquet::file::reader::FileReader;
    use parquet::file::serialized_reader::SerializedFileReader;
    use parquet::record::Field;
    use parquet_derive::ParquetRecordReader;
    use pretty_assertions::assert_eq;
    use proptest::arbitrary::any;
    use proptest::collection::vec;
    use proptest::proptest;

    #[derive(ParquetRecordWriter, ParquetRecordReader, PartialEq, Debug)]
    struct BarcodeUmi {
        barcode: String,
        umi: String,
    }

    impl BarcodeUmi {
        fn test_data() -> Vec<Self> {
            vec![
                BarcodeUmi {
                    barcode: "ACGT".to_string(),
                    umi: "AA".to_string(),
                },
                BarcodeUmi {
                    barcode: "TAGC".to_string(),
                    umi: "AC".to_string(),
                },
                BarcodeUmi {
                    barcode: "TAGG".to_string(),
                    umi: "AT".to_string(),
                },
            ]
        }
    }

    #[test]
    fn test_parquet_file_writer() -> Result<()> {
        let file = ParquetFile::tempfile()?;
        let mut writer = file.writer::<BarcodeUmi>(2)?;
        writer.write_all(&BarcodeUmi::test_data())?;
        writer.close()?;

        let reader = SerializedFileReader::try_from(file.as_ref().as_ref())?;

        assert_eq!(reader.num_row_groups(), 2);

        let entries = reader
            .into_iter()
            .map(|row| {
                row.unwrap()
                    .get_column_iter()
                    .map(|col| (col.0.clone(), col.1.clone()))
                    .collect_vec()
            })
            .collect::<Vec<_>>();

        let f = |x: &str| Field::Str(x.to_string());
        assert_eq!(
            entries,
            vec![
                vec![
                    ("barcode".to_string(), f("ACGT")),
                    ("umi".to_string(), f("AA"))
                ],
                vec![
                    ("barcode".to_string(), f("TAGC")),
                    ("umi".to_string(), f("AC"))
                ],
                vec![
                    ("barcode".to_string(), f("TAGG")),
                    ("umi".to_string(), f("AT"))
                ]
            ]
        );

        Ok(())
    }

    #[test]
    fn test_parquet_file_reader() -> Result<()> {
        let data = BarcodeUmi::test_data();
        let file = ParquetFile::tempfile()?;
        {
            let mut writer = file.writer::<BarcodeUmi>(2)?;
            writer.write_all(&data)?;
            writer.close()?;
        }

        let reader = file.reader::<BarcodeUmi>()?;
        let entries = reader.collect::<Result<Vec<_>>>()?;

        assert_eq!(entries, data);
        Ok(())
    }

    proptest! {
        #[test]
        fn prop_test_parquet_roundtrip(
            items in vec(any::<(String, String)>(), 0usize..10_000usize),
            row_group_size in 1..10_000usize,
        ) {
            let data = items.into_iter().map(|(barcode, umi)| BarcodeUmi { barcode, umi }).collect::<Vec<_>>();
            let file = ParquetFile::tempfile().unwrap();
            let mut writer = file.writer::<BarcodeUmi>(row_group_size).unwrap();
            writer.write_all(&data).unwrap();
            writer.close().unwrap();

            let reader = file.reader::<BarcodeUmi>().unwrap();
            let entries = reader.collect::<Result<Vec<_>>>().unwrap();
            assert_eq!(entries, data);
        }
    }
}
