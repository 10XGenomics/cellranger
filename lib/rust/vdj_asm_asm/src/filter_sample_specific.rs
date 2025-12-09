//! Martian stage FILTER_SAMPLE_SPECIFIC
#![expect(missing_docs)]

use anyhow::Result;
use barcode::whitelist::BarcodeId;
use cr_types::chemistry::{ChemistryDefs, ChemistryDefsExt};
use cr_types::filtered_barcodes::FilteredBarcodesCsv;
use cr_types::{Fingerprint, FingerprintFile};
use itertools::Itertools;
use martian::{MartianMain, MartianRover};
use martian_derive::{MartianStruct, make_mro};
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::{FileTypeRead, LazyFileTypeIO, LazyWrite};
use metric::TxHashSet;
use serde::{Deserialize, Serialize};
use std::collections::HashSet;
use vdj_ann::annotate::ContigAnnotation;
use vdj_filter_barcodes::filter_barcode_level::overhang_demux_filter;
use vdj_filter_barcodes::filter_log::{VdjFilterLogFormat, extend_filter_log};

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct FilterSampleSpecificStageInputs {
    pub contig_annotations: JsonFile<Vec<ContigAnnotation>>,
    pub vdj_chemistry_def: ChemistryDefs,
    pub sample_fingerprint: Option<FingerprintFile>,
    pub filtered_barcodes: Option<FilteredBarcodesCsv>,
    pub is_antibody_only: Option<bool>,
    pub is_non_targeted_gex: Option<bool>,
    pub filter_diagnostics: VdjFilterLogFormat,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct FilterSampleSpecificStageOutputs {
    pub contig_annotations: JsonFile<Vec<ContigAnnotation>>,
    pub filter_diagnostics: VdjFilterLogFormat,
}

pub struct FilterSampleSpecific;

#[make_mro(mem_gb = 4)]
impl MartianMain for FilterSampleSpecific {
    type StageInputs = FilterSampleSpecificStageInputs;
    type StageOutputs = FilterSampleSpecificStageOutputs;

    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        // Set-up filter logging
        let filter_diagnostics_file: VdjFilterLogFormat = rover.make_path("filter_diagnostics");
        let mut filter_logger =
            extend_filter_log(args.filter_diagnostics, &filter_diagnostics_file);

        // OH demux filter
        let demuxed_bcs = if let (Some(overhang_read_barcode), Some(fingerprint)) = (
            args.vdj_chemistry_def.overhang_read_barcode(),
            args.sample_fingerprint,
        ) {
            let valid_overhang_ids: HashSet<BarcodeId> = fingerprint
                .read()?
                .iter()
                .flat_map(Fingerprint::tag_names)
                .map(|tag_name| BarcodeId::pack(tag_name))
                .collect();
            let barcodes: Vec<String> = args
                .contig_annotations
                .lazy_reader()?
                .map(|ann| ann.unwrap().barcode)
                .collect();
            let demuxed_bcs = overhang_demux_filter(
                overhang_read_barcode,
                valid_overhang_ids,
                barcodes,
                Some(&mut filter_logger),
            );
            Some(demuxed_bcs)
        } else {
            None
        };

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
        for ann in args.contig_annotations.lazy_reader()? {
            let ann: ContigAnnotation = ann?;
            writer.write_item(&ContigAnnotation {
                is_asm_cell: Some(
                    demuxed_bcs
                        .as_ref()
                        .map_or(ann.is_cell, |bcs| bcs.contains(&ann.barcode) && ann.is_cell),
                ),
                is_gex_cell: gex_cells.as_ref().map(|cells| cells.contains(&ann.barcode)),
                ..ann
            })?;
        }
        writer.finish()?;

        Ok(FilterSampleSpecificStageOutputs {
            contig_annotations,
            filter_diagnostics: filter_diagnostics_file,
        })
    }
}
