//! Martian stage HANDLE_GEX_CELLS

use anyhow::Result;
use cr_types::filtered_barcodes::FilteredBarcodesCsv;
use itertools::Itertools;
use martian::{MartianMain, MartianRover};
use martian_derive::{make_mro, MartianStruct};
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::{LazyFileTypeIO, LazyWrite};
use metric::TxHashSet;
use serde::{Deserialize, Serialize};
use vdj_ann::annotate::ContigAnnotation;

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct HandleGexCellsStageInputs {
    pub asm_contig_annotations: JsonFile<Vec<ContigAnnotation>>,
    pub filtered_barcodes: Option<FilteredBarcodesCsv>,
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
        let gex_cells: Option<TxHashSet<_>> =
            if let (Some(false), Some(true), Some(filtered_barcodes)) = (
                args.is_antibody_only,
                args.is_non_targeted_gex,
                args.filtered_barcodes,
            ) {
                Some(
                    filtered_barcodes
                        .lazy_reader()?
                        .map_ok(|x| x.barcode.to_string())
                        .try_collect()?,
                )
            } else {
                None
            };

        let contig_annotations: JsonFile<Vec<ContigAnnotation>> =
            rover.make_path("contig_annotations");
        let mut writer = contig_annotations.lazy_writer()?;
        for ann in args.asm_contig_annotations.lazy_reader()? {
            let ann: ContigAnnotation = ann?;
            writer.write_item(&ContigAnnotation {
                is_asm_cell: Some(ann.is_cell),
                is_gex_cell: gex_cells.as_ref().map(|cells| cells.contains(&ann.barcode)),
                ..ann
            })?;
        }
        writer.finish()?;

        Ok(HandleGexCellsStageOutputs { contig_annotations })
    }
}
