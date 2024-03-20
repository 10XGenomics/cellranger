//! SubsetAssemblyOuts stage code

use crate::assembly_types::{AsmReadsPerBcFormat, BarcodeSupport, ContigSummaryRow, UmiSummaryRow};
use crate::BarcodeDataBriefFile;
use anyhow::Result;
use barcode::whitelist::BarcodeId;
use cr_types::chemistry::{ChemistryDefs, ChemistryDefsExt};
use cr_types::{Fingerprint, FingerprintFile};
use itertools::Itertools;
use martian::prelude::*;
use martian_derive::{make_mro, MartianStruct};
use martian_filetypes::json_file::{Json, JsonFile, JsonFormat};
use martian_filetypes::tabular_file::{CsvFile, TsvFile};
use martian_filetypes::{FileTypeRead, FileTypeWrite, LazyFileTypeIO, LazyWrite};
use metric::{SimpleHistogram, TxHashSet};
use serde::{Deserialize, Serialize};
use std::collections::HashSet;
use vdj_ann::annotate::ContigAnnotation;
use vdj_asm_utils::filter_barcodes::map_multiplexing_seq_to_id;

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct SubsetAssemblyOutsStageInputs {
    pub per_sample: bool,
    pub count_chemistry_defs: Option<ChemistryDefs>,
    pub sample_fingerprint: Option<FingerprintFile>,
    pub contig_annotations: Option<JsonFile<Vec<ContigAnnotation>>>,
    pub merged_annotations: Option<JsonFile<Vec<ContigAnnotation>>>,
    pub total_read_pairs: i64,
    pub corrected_barcode_counts: JsonFile<SimpleHistogram<String>>,
    pub assemblable_reads_per_bc: AsmReadsPerBcFormat,
    pub contig_summary: TsvFile<ContigSummaryRow>,
    pub umi_summary: TsvFile<UmiSummaryRow>,
    pub barcode_support: CsvFile<BarcodeSupport>,
    pub barcode_brief: BarcodeDataBriefFile,
}

#[derive(Debug, Serialize, Deserialize, MartianStruct)]
pub struct SubsetAssemblyOutsStageOutputs {
    pub contig_annotations: JsonFile<Vec<ContigAnnotation>>,
    pub total_read_pairs: i64,
    pub corrected_barcode_counts: JsonFile<SimpleHistogram<String>>,
    pub assemblable_reads_per_bc: AsmReadsPerBcFormat,
    pub contig_summary: TsvFile<ContigSummaryRow>,
    pub umi_summary: TsvFile<UmiSummaryRow>,
    pub barcode_support: CsvFile<BarcodeSupport>,
    pub barcode_brief: BarcodeDataBriefFile,
}

pub struct SubsetAssemblyOuts;

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

        // Set-up for filtering barcodes based on overhang
        let maybe_overhang_read_component = args
            .count_chemistry_defs
            .iter()
            .filter_map(ChemistryDefs::overhang_read_barcode)
            .at_most_one()
            .unwrap();

        if !args.per_sample || maybe_overhang_read_component.is_none() {
            return Ok(SubsetAssemblyOutsStageOutputs {
                contig_annotations: annot,
                total_read_pairs: args.total_read_pairs,
                corrected_barcode_counts: args.corrected_barcode_counts,
                assemblable_reads_per_bc: args.assemblable_reads_per_bc,
                contig_summary: args.contig_summary,
                umi_summary: args.umi_summary,
                barcode_support: args.barcode_support,
                barcode_brief: args.barcode_brief,
            });
        }

        let overhang_read_component = maybe_overhang_read_component.unwrap();

        let overhang_range = overhang_read_component.offset()
            ..overhang_read_component.offset() + overhang_read_component.length();
        let overhang_seq_to_id = overhang_read_component.build_seq_to_id_map()?;
        let valid_overhang_ids: HashSet<BarcodeId> = args
            .sample_fingerprint
            .unwrap()
            .read()?
            .iter()
            .flat_map(Fingerprint::tag_names)
            .map(|tag_name| BarcodeId::pack(tag_name))
            .collect();

        // subset barcode data brief
        let mut total_read_pairs: i64 = 0;
        let mut valid_bcs = TxHashSet::default();
        let barcode_data_brief_file: BarcodeDataBriefFile = rover.make_path("barcode_data_brief");
        let mut barcode_data_brief_writer = barcode_data_brief_file.lazy_writer()?;
        for brief in args.barcode_brief.lazy_reader()? {
            let brief = brief?;
            let this_overhang =
                map_multiplexing_seq_to_id(&brief.barcode, &overhang_seq_to_id, &overhang_range);
            if valid_overhang_ids.contains(&this_overhang) {
                valid_bcs.insert(brief.barcode.clone());
                total_read_pairs += brief.read_pairs as i64;
                barcode_data_brief_writer.write_item(&brief)?;
            }
        }
        barcode_data_brief_writer.finish()?;

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

        // subset assemblable_reads_per_bc
        let mut bc_read_counts = args.assemblable_reads_per_bc.read()?;
        bc_read_counts.retain(|k, _| valid_bcs.contains(k));
        let assemblable_reads_per_bc_file: AsmReadsPerBcFormat =
            rover.make_path("assemblable_reads_per_bc");
        assemblable_reads_per_bc_file.write(&bc_read_counts)?;

        // subset corrected_barcode_counts
        let mut bc_counts_corrected = args.corrected_barcode_counts.read()?;
        bc_counts_corrected.retain(|k, _| valid_bcs.contains(k));
        let corrected_barcode_counts_file: JsonFile<SimpleHistogram<String>> =
            rover.make_path("corrected_barcode_counts");
        corrected_barcode_counts_file.write(&bc_counts_corrected)?;

        // subset contig_summary.tsv
        let contig_summary_file: TsvFile<ContigSummaryRow> = rover.make_path("contig_summary");
        let mut contig_summary_writer = contig_summary_file.lazy_writer()?;
        for row in args.contig_summary.lazy_reader()? {
            let row = row?;
            if valid_bcs.contains(&row.barcode) {
                contig_summary_writer.write_item(&row)?;
            }
        }
        contig_summary_writer.finish()?;

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
            total_read_pairs,
            corrected_barcode_counts: corrected_barcode_counts_file,
            assemblable_reads_per_bc: assemblable_reads_per_bc_file,
            contig_summary: contig_summary_file,
            umi_summary: umi_summary_file,
            barcode_support: barcode_support_file,
            barcode_brief: barcode_data_brief_file,
        })
    }
}
