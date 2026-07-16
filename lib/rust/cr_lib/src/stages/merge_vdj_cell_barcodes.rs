//! Martian stage MERGE_VDJ_CELL_BARCODES
//!
//! Merge the union of cell-associated barcodes from all VDJ libraries.
//! Output `vdj_cell_barcodes` as a `FilteredBarcodesCsv`.
#![deny(missing_docs)]

use anyhow::Result;
use barcode::Barcode;
use cr_types::filtered_barcodes::{FilteredBarcodesCsv, FilteredBarcodesCsvRow};
use martian::{MartianMain, MartianRover};
use martian_derive::{MartianStruct, make_mro};
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::{FileTypeRead, LazyFileTypeIO, LazyWrite};
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;

/// Stage inputs.
#[derive(Clone, Deserialize, MartianStruct)]
pub struct MergeVdjCellBarcodesStageInputs {
    /// Cell barcodes from the VDJ-T library-level analyzer (may be null).
    pub vdj_t_cell_barcodes: Option<JsonFile<Vec<String>>>,
    /// Cell barcodes from the VDJ-T-GD library-level analyzer (may be null).
    pub vdj_t_gd_cell_barcodes: Option<JsonFile<Vec<String>>>,
    /// Cell barcodes from the VDJ-B library-level analyzer (may be null).
    pub vdj_b_cell_barcodes: Option<JsonFile<Vec<String>>>,
}

/// Stage outputs.
#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct MergeVdjCellBarcodesStageOutputs {
    /// Union of all VDJ cell barcodes in `FilteredBarcodesCsv` format.
    pub vdj_cell_barcodes: FilteredBarcodesCsv,
}

/// Martian stage MERGE_VDJ_CELL_BARCODES.
pub struct MergeVdjCellBarcodes;

#[make_mro(mem_gb = 2, volatile = strict)]
impl MartianMain for MergeVdjCellBarcodes {
    type StageInputs = MergeVdjCellBarcodesStageInputs;
    type StageOutputs = MergeVdjCellBarcodesStageOutputs;

    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        // Read each non-null input file into a Vec<Barcode>.
        let vdj_t = args.vdj_t_cell_barcodes.map(read_barcodes).transpose()?;
        let vdj_t_gd = args.vdj_t_gd_cell_barcodes.map(read_barcodes).transpose()?;
        let vdj_b = args.vdj_b_cell_barcodes.map(read_barcodes).transpose()?;

        // Union of all three Vecs into a BTreeSet
        let union: BTreeSet<Barcode> = [vdj_t, vdj_t_gd, vdj_b]
            .into_iter()
            .flatten()
            .flatten()
            .collect();

        // Write cell barcodes as a FilteredBarcodesCsv with genome="VDJ"
        let vdj_cell_barcodes: FilteredBarcodesCsv = rover.make_path("vdj_cell_barcodes");
        {
            let mut writer = vdj_cell_barcodes.lazy_writer()?;
            for barcode in &union {
                writer.write_item(&FilteredBarcodesCsvRow {
                    genome: "VDJ".into(),
                    barcode: *barcode,
                })?;
            }
        }

        Ok(MergeVdjCellBarcodesStageOutputs { vdj_cell_barcodes })
    }
}

/// Read a VDJ cell-barcode JSON file into a `Vec<Barcode>`.
fn read_barcodes(file: JsonFile<Vec<String>>) -> Result<Vec<Barcode>> {
    let strings: Vec<String> = file.read()?;
    Ok(strings.iter().map(|s| Barcode::from(s.as_str())).collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use martian_filetypes::FileTypeWrite;
    use tempfile::tempdir;

    fn write_barcode_file(
        dir: &std::path::Path,
        name: &str,
        barcodes: &[&str],
    ) -> JsonFile<Vec<String>> {
        let file: JsonFile<Vec<String>> = JsonFile::from(dir.join(name));
        let vec: Vec<String> = barcodes.iter().map(ToString::to_string).collect();
        file.write(&vec).unwrap();
        file
    }

    #[test]
    fn test_union_deduplicates_across_inputs() {
        let dir = tempdir().unwrap();
        let f1 = write_barcode_file(dir.path(), "f1.json", &["AAAA-1", "CCCC-1"]);
        let f2 = write_barcode_file(dir.path(), "f2.json", &["CCCC-1", "TTTT-1"]);
        let f3 = write_barcode_file(dir.path(), "f3.json", &["GGGG-1"]);

        let v1 = read_barcodes(f1).unwrap();
        let v2 = read_barcodes(f2).unwrap();
        let v3 = read_barcodes(f3).unwrap();
        let union: BTreeSet<Barcode> = [v1, v2, v3].into_iter().flatten().collect();
        assert_eq!(union.len(), 4);
        assert!(union.contains(&Barcode::from("CCCC-1"))); // CCCC-1 appears in f1 and f2 but counted once
    }

    #[test]
    fn test_single_input_preserved() {
        let dir = tempdir().unwrap();
        let f = write_barcode_file(dir.path(), "b.json", &["AAAA-1", "CCCC-1", "GGGG-1"]);
        let barcodes = read_barcodes(f).unwrap();
        assert_eq!(barcodes.len(), 3);
    }

    #[test]
    fn test_empty_input_gives_empty_union() {
        let dir = tempdir().unwrap();
        let f = write_barcode_file(dir.path(), "empty.json", &[]);
        let barcodes = read_barcodes(f).unwrap();
        assert!(barcodes.is_empty());
    }
}
