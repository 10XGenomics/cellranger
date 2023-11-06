//! HandleGexCells stage code

use anyhow::Result;
use martian::prelude::*;
use martian_derive::{make_mro, MartianStruct};
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::tabular_file::CsvFileNoHeader;
use martian_filetypes::{FileTypeRead, LazyFileTypeIO, LazyWrite};
use serde::{Deserialize, Serialize};
use std::collections::HashSet;
use vdj_ann::annotate::ContigAnnotation;

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct FilteredBarcodesCsvRow {
    genome: String,
    barcode: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct HandleGexCellsStageInputs {
    pub asm_contig_annotations: JsonFile<Vec<ContigAnnotation>>,
    pub filtered_barcodes: Option<CsvFileNoHeader<FilteredBarcodesCsvRow>>,
    pub is_antibody_only: Option<bool>,
    pub is_non_targeted_gex: Option<bool>,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct HandleGexCellsStageOutputs {
    pub contig_annotations: JsonFile<Vec<ContigAnnotation>>,
}

pub struct HandleGexCells;

#[make_mro]
impl MartianMain for HandleGexCells {
    type StageInputs = HandleGexCellsStageInputs;
    type StageOutputs = HandleGexCellsStageOutputs;
    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        // We will update `is_gex_cell` only when we have gene expression data available. In the
        // case of antibody only or targeted gene expression, we do not use the cell calls made
        // in the count pipeline as of now.
        let gex_cells = match (
            args.filtered_barcodes,
            args.is_antibody_only,
            args.is_non_targeted_gex,
        ) {
            (Some(csv), Some(false), Some(true)) => {
                let rows = csv.read()?;
                Some(rows.into_iter().map(|r| r.barcode).collect::<HashSet<_>>())
            }
            _ => None,
        };

        let reader = args.asm_contig_annotations.lazy_reader()?;
        let contig_annotations: JsonFile<Vec<ContigAnnotation>> =
            rover.make_path("contig_annotations");
        let mut writer = contig_annotations.lazy_writer()?;
        for ann in reader {
            let mut ann: ContigAnnotation = ann?;
            ann.is_asm_cell = Some(ann.is_cell);
            ann.is_gex_cell = gex_cells.as_ref().map(|cells| cells.contains(&ann.barcode));
            writer.write_item(&ann)?;
        }
        writer.finish()?;

        Ok(HandleGexCellsStageOutputs { contig_annotations })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_load_barcode_csv() -> Result<()> {
        let rows: Vec<FilteredBarcodesCsvRow> =
            CsvFileNoHeader::from("test_resources/filtered_barcodes.csv").read()?;

        let genome = "GRCh38".to_string();
        assert_eq!(
            rows,
            vec![
                FilteredBarcodesCsvRow {
                    genome: genome.clone(),
                    barcode: "AAACCTGCATCCCATC-1".into()
                },
                FilteredBarcodesCsvRow {
                    genome: genome.clone(),
                    barcode: "AAAGCAACACCGAAAG-1".into()
                },
                FilteredBarcodesCsvRow {
                    genome,
                    barcode: "AAATGCCCACATTTCT-1".into()
                }
            ]
        );
        Ok(())
    }
}
