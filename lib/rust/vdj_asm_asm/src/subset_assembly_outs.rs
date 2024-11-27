//! SubsetAssemblyOuts stage code

use crate::assembly::BarcodeDataFile;
use crate::assembly_types::{AsmReadsPerBcFormat, BarcodeSupport, UmiSummaryRow};
use crate::BarcodeDataBriefFile;
use anyhow::Result;
use barcode::whitelist::BarcodeId;
use cr_types::chemistry::{ChemistryDefs, ChemistryDefsExt};
use cr_types::{BarcodeMultiplexingType, Fingerprint, FingerprintFile, ReadLevel};
use martian::prelude::*;
use martian_derive::{make_mro, MartianStruct};
use martian_filetypes::json_file::{Json, JsonFile, JsonFormat};
use martian_filetypes::tabular_file::{CsvFile, TsvFile};
use martian_filetypes::{FileTypeRead, FileTypeWrite, LazyFileTypeIO, LazyWrite};
use metric::{SimpleHistogram, TxHashSet};
use serde::{Deserialize, Serialize};
use std::collections::HashSet;
use vdj_ann::annotate::ContigAnnotation;
use vdj_filter_barcodes::filter_barcode_level::map_multiplexing_seq_to_id;

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct SubsetAssemblyOutsStageInputs {
    pub per_sample: bool,
    pub multiplexing_method: Option<BarcodeMultiplexingType>,
    pub vdj_chemistry_def: ChemistryDefs,
    pub sample_fingerprint: Option<FingerprintFile>,
    pub contig_annotations: Option<JsonFile<Vec<ContigAnnotation>>>,
    pub merged_annotations: Option<JsonFile<Vec<ContigAnnotation>>>,
    pub total_read_pairs: i64,
    pub corrected_barcode_counts: JsonFile<SimpleHistogram<String>>,
    pub assemblable_reads_per_bc: AsmReadsPerBcFormat,
    pub umi_summary: TsvFile<UmiSummaryRow>,
    pub barcode_support: CsvFile<BarcodeSupport>,
    pub barcode_brief: BarcodeDataBriefFile,
    pub barcode_full: BarcodeDataFile,
}

#[derive(Debug, Serialize, Deserialize, MartianStruct)]
pub struct SubsetAssemblyOutsStageOutputs {
    pub contig_annotations: JsonFile<Vec<ContigAnnotation>>,
    pub total_read_pairs: Option<i64>,
    pub corrected_barcode_counts: Option<JsonFile<SimpleHistogram<String>>>,
    pub assemblable_reads_per_bc: AsmReadsPerBcFormat,
    pub umi_summary: TsvFile<UmiSummaryRow>,
    pub barcode_support: Option<CsvFile<BarcodeSupport>>,
    pub barcode_brief: BarcodeDataBriefFile,
    pub barcode_full: BarcodeDataFile,
}

pub struct SubsetAssemblyOuts;

fn get_valid_barcodes_read_level_multiplexing(
    vdj_chemistry_def: &ChemistryDefs,
    fingerprints: &[Fingerprint],
    contig_annotations: &JsonFile<Vec<ContigAnnotation>>,
) -> TxHashSet<String> {
    let overhang_read_component = vdj_chemistry_def.overhang_read_barcode().unwrap();
    let overhang_range = overhang_read_component.offset()
        ..overhang_read_component.offset() + overhang_read_component.length();
    let overhang_seq_to_id = overhang_read_component.build_seq_to_id_map().unwrap();
    let valid_overhang_ids: HashSet<BarcodeId> = fingerprints
        .iter()
        .flat_map(Fingerprint::tag_names)
        .map(|tag_name| BarcodeId::pack(tag_name))
        .collect();
    let mut valid_bcs = TxHashSet::default();
    let contig_reader = contig_annotations.lazy_reader().unwrap();
    for ann in contig_reader {
        let ann: ContigAnnotation = ann.unwrap();
        let this_overhang =
            map_multiplexing_seq_to_id(&ann.barcode, &overhang_seq_to_id, &overhang_range);
        if valid_overhang_ids.contains(&this_overhang) {
            valid_bcs.insert(ann.barcode.clone());
        }
    }
    valid_bcs
}

fn get_valid_barcodes_cell_level_multiplexing(
    contig_annotations: &JsonFile<Vec<ContigAnnotation>>,
) -> TxHashSet<String> {
    let mut valid_bcs = TxHashSet::default();
    let contig_reader = contig_annotations.lazy_reader().unwrap();
    for ann in contig_reader {
        let ann: ContigAnnotation = ann.unwrap();
        // Was the barcode declared a cell by both the VDJ assembler and Gene Expression
        if ann.is_asm_cell.map_or(false, |a| a) && ann.is_gex_cell.map_or(false, |g| g) {
            valid_bcs.insert(ann.barcode.clone());
        }
    }
    valid_bcs
}

#[make_mro]
impl MartianMain for SubsetAssemblyOuts {
    type StageInputs = SubsetAssemblyOutsStageInputs;
    type StageOutputs = SubsetAssemblyOutsStageOutputs;
    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        let annot = match (args.contig_annotations, args.merged_annotations) {
            (Some(annot), None) => annot,
            (None, Some(annot)) => annot,
            (None, None) => unreachable!(),
            (Some(_), Some(_)) => unreachable!(),
        };

        if !args.per_sample || args.multiplexing_method.is_none() {
            return Ok(SubsetAssemblyOutsStageOutputs {
                contig_annotations: annot,
                total_read_pairs: Some(args.total_read_pairs),
                corrected_barcode_counts: Some(args.corrected_barcode_counts),
                assemblable_reads_per_bc: args.assemblable_reads_per_bc,
                umi_summary: args.umi_summary,
                barcode_support: Some(args.barcode_support),
                barcode_brief: args.barcode_brief,
                barcode_full: args.barcode_full,
            });
        }

        let (is_read_level_multiplexed, valid_bcs) = match args.multiplexing_method.unwrap() {
            BarcodeMultiplexingType::ReadLevel(ReadLevel::OH) => (
                true,
                get_valid_barcodes_read_level_multiplexing(
                    &args.vdj_chemistry_def,
                    &args.sample_fingerprint.unwrap().read()?,
                    &annot,
                ),
            ),
            BarcodeMultiplexingType::CellLevel(_) => {
                (false, get_valid_barcodes_cell_level_multiplexing(&annot))
            }
            BarcodeMultiplexingType::ReadLevel(ReadLevel::RTL) => {
                panic!("Unsupported multiplexing method!")
            }
        };

        // subset contig annotaions
        let contig_ann_json_file: JsonFormat<Json, Vec<ContigAnnotation>> =
            rover.make_path("contig_annotations.json");
        let mut contig_writer = contig_ann_json_file.lazy_writer()?;
        let contig_reader = annot.lazy_reader()?;
        for ann in contig_reader {
            let ann: ContigAnnotation = ann?;
            if valid_bcs.contains(&ann.barcode) {
                contig_writer.write_item(&ann)?;
            }
        }
        contig_writer.finish()?;

        // subset barcode data brief
        let mut total_read_pairs: i64 = 0;
        let barcode_data_brief_file: BarcodeDataBriefFile = rover.make_path("barcode_data_brief");
        let mut barcode_data_brief_writer = barcode_data_brief_file.lazy_writer()?;
        for brief in args.barcode_brief.lazy_reader()? {
            let brief = brief?;
            if valid_bcs.contains(&brief.barcode) {
                total_read_pairs += brief.read_pairs as i64;
                barcode_data_brief_writer.write_item(&brief)?;
            }
        }
        barcode_data_brief_writer.finish()?;

        // subset barcode data sum
        let barcode_full_file: BarcodeDataFile = rover.make_path("barcode_data");
        let mut barcode_full_writer = barcode_full_file.lazy_writer()?;
        for bc_data in args.barcode_full.lazy_reader()? {
            let bc_data = bc_data?;
            if valid_bcs.contains(&bc_data.barcode) {
                barcode_full_writer.write_item(&bc_data)?;
            }
        }

        // subset assemblable_reads_per_bc
        let mut bc_read_counts = args.assemblable_reads_per_bc.read()?;
        bc_read_counts.retain(|k, _| valid_bcs.contains(k));
        let assemblable_reads_per_bc_file: AsmReadsPerBcFormat =
            rover.make_path("assemblable_reads_per_bc");
        assemblable_reads_per_bc_file.write(&bc_read_counts)?;

        // subset corrected_barcode_counts
        let mut bc_counts_corrected = args.corrected_barcode_counts.read()?;
        bc_counts_corrected.retain(|k, _| valid_bcs.contains(k));
        let corrected_barcode_counts_file: JsonFile<_> =
            rover.make_path("corrected_barcode_counts");
        corrected_barcode_counts_file.write(&bc_counts_corrected)?;

        // subset umi_summary.tsv
        let umi_summary_file: TsvFile<UmiSummaryRow> = rover.make_path("umi_summary");
        let mut umi_summary_writer = umi_summary_file.lazy_writer()?;
        for row in args.umi_summary.lazy_reader()? {
            let row = row?;
            if valid_bcs.contains(&row.barcode) {
                umi_summary_writer.write_item(&row)?;
            }
        }
        umi_summary_writer.finish()?;

        // subset barcode_support.csv
        let barcode_support_file: CsvFile<BarcodeSupport> = rover.make_path("barcode_support");
        let mut barcode_support_writer = barcode_support_file.lazy_writer()?;
        for row in args.barcode_support.lazy_reader()? {
            let row = row?;
            if valid_bcs.contains(&row.barcode) {
                barcode_support_writer.write_item(&row)?;
            }
        }
        barcode_support_writer.finish()?;

        Ok(SubsetAssemblyOutsStageOutputs {
            contig_annotations: contig_ann_json_file,
            total_read_pairs: is_read_level_multiplexed.then_some(total_read_pairs),
            corrected_barcode_counts: is_read_level_multiplexed
                .then_some(corrected_barcode_counts_file),
            assemblable_reads_per_bc: assemblable_reads_per_bc_file,
            umi_summary: umi_summary_file,
            barcode_support: is_read_level_multiplexed.then_some(barcode_support_file),
            barcode_brief: barcode_data_brief_file,
            barcode_full: barcode_full_file,
        })
    }
}
