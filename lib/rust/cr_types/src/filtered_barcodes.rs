use crate::types::GenomeName;
use anyhow::Result;
use barcode::{barcode_string, Barcode};
use itertools::Itertools;
use martian_filetypes::tabular_file::CsvFileNoHeader;
use martian_filetypes::LazyFileTypeIO;
use metric::TxHashSet;
use serde::{Deserialize, Serialize};

/// Filtered barcodes CSV
pub type FilteredBarcodesCsv = CsvFileNoHeader<FilteredBarcodesCsvRow>;

/// Filtered barcodes CSV row
#[derive(Debug, Serialize, Deserialize, PartialEq)]
pub struct FilteredBarcodesCsvRow {
    pub genome: GenomeName,
    #[serde(with = "barcode_string")]
    pub barcode: Barcode,
}

/// Read the filtered_barcodes CSV file and return the set of barcodes
pub fn read_filtered_barcodes_set(
    filtered_barcodes_filename: &FilteredBarcodesCsv,
) -> Result<TxHashSet<Barcode>> {
    filtered_barcodes_filename
        .lazy_reader()?
        .map(|row| row.map(|r| r.barcode))
        .try_collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use martian_filetypes::FileTypeRead;

    #[test]
    fn test_load_barcode_csv() -> Result<()> {
        assert_eq!(
            FilteredBarcodesCsv::from("test/filtered_barcodes.csv").read()?,
            vec![
                FilteredBarcodesCsvRow {
                    genome: "GRCh38".into(),
                    barcode: "AAACCTGCATCCCATC-1".into()
                },
                FilteredBarcodesCsvRow {
                    genome: "GRCh38".into(),
                    barcode: "AAAGCAACACCGAAAG-1".into()
                },
                FilteredBarcodesCsvRow {
                    genome: "GRCh38".into(),
                    barcode: "AAATGCCCACATTTCT-1".into()
                }
            ]
        );
        Ok(())
    }
}
