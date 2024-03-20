//! AsmCallCells stage code

use crate::assembly::BarcodeDataBriefFile;
use anyhow::Result;
use barcode::whitelist::BarcodeId;
use cr_types::chemistry::{ChemistryDefs, ChemistryDefsExt};
use cr_types::{Fingerprint, FingerprintFile};
use martian::prelude::*;
use martian_derive::{make_mro, MartianStruct};
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::lz4_file::Lz4;
use martian_filetypes::{FileTypeRead, LazyFileTypeIO, LazyWrite};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeSet, HashMap, HashSet};
use std::path::PathBuf;
use vdj_ann::annotate::ContigAnnotation;
use vdj_asm_utils::barcode_data::analyze_barcode_data_brief;
use vdj_asm_utils::filter_barcodes::{
    cell_filter, confidence_filter, overhang_demux_filter, BarcodeCellInfo, BarcodeFilteringParams,
    Contigs,
};
use vdj_asm_utils::filter_log::{FilterLogEntry, FilterLogger, FilterSwitch};
use vdj_reference::VdjReceptor;
use vdj_types::VdjChain;
use vector_utils::bin_member;

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct AsmCallCellsStageInputs {
    pub receptor: Option<VdjReceptor>,
    pub denovo: bool,
    pub vdj_reference_path: Option<PathBuf>,
    pub count_chemistry_defs: Option<ChemistryDefs>,
    pub contig_annotations: JsonFile<Vec<ContigAnnotation>>,
    pub barcode_brief: BarcodeDataBriefFile,
    pub n50_n50_rpu: u32,
    pub filter_switch: FilterSwitch,
    pub sample_fingerprint: Option<FingerprintFile>,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct AsmCallCellsStageOutputs {
    pub contig_annotations: JsonFile<Vec<ContigAnnotation>>,
    #[mro_retain]
    pub filter_diagnostics: Lz4<JsonFile<Vec<FilterLogEntry>>>,
}

// This is our stage struct
pub struct AsmCallCells;

#[make_mro(mem_gb = 4)]
impl MartianMain for AsmCallCells {
    type StageInputs = AsmCallCellsStageInputs;
    type StageOutputs = AsmCallCellsStageOutputs;
    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        let n50_n50_rpu: i32 = args.n50_n50_rpu.try_into().unwrap();
        let mut contigs_by_bc = HashMap::<String, Contigs>::new();
        for ann in args.contig_annotations.lazy_reader()? {
            let ann: ContigAnnotation = ann?;

            let contig = contigs_by_bc.entry(ann.barcode.clone()).or_default();
            match ann.junction_support.is_some() {
                true => contig.good_contigs.push(ann),
                false => contig.reject_contigs.push(ann),
            };
        }
        let is_tcr =
            args.receptor == Some(VdjReceptor::TR) || args.receptor == Some(VdjReceptor::TRGD);
        let is_bcr = args.receptor == Some(VdjReceptor::IG);
        let gd_mode = args.receptor == Some(VdjReceptor::TRGD);
        let filter_diagnostics_file: Lz4<JsonFile<_>> = rover.make_path("filter_diagnostics");
        let mut filter_logger = FilterLogger::new(&filter_diagnostics_file)?;
        let mut barcode_cell_info = Vec::<BarcodeCellInfo>::new();

        let mut confident_bcs = HashSet::new();
        let reader = args.barcode_brief.lazy_reader()?;
        for brief in reader {
            let mut bc: BarcodeCellInfo = brief?.into();
            let contigs = contigs_by_bc.remove(&bc.barcode).unwrap_or_default();
            let filtering_params =
                BarcodeFilteringParams::build(&contigs, &bc, args.denovo, n50_n50_rpu, gd_mode)?;
            let mut low_confidence_reasons = Vec::new();
            bc.high_confidence =
                confidence_filter(&filtering_params, n50_n50_rpu, &mut low_confidence_reasons);

            if bc.high_confidence {
                confident_bcs.insert(bc.barcode.clone());
            };
            bc.now_a_cell = cell_filter(
                &filtering_params,
                &bc,
                args.denovo,
                is_tcr,
                is_bcr,
                n50_n50_rpu,
                Some(&mut filter_logger),
                low_confidence_reasons,
            );
            bc.paired = filtering_params.paired;
            bc.chimdata = contigs.build_chimdata(bc.now_a_cell, args.denovo)?;
            bc.jundata = contigs.build_jundata(bc.high_confidence)?;
            barcode_cell_info.push(bc);
        }

        let mut kills = Vec::<String>::new();
        let mut killsc = Vec::<String>::new();
        analyze_barcode_data_brief(
            &barcode_cell_info,
            args.filter_switch,
            &mut kills,
            &mut killsc,
            Some(&mut filter_logger),
        );

        // Update now_a_cell field in barcode_cell_info.
        for bc in &mut barcode_cell_info {
            if bc.now_a_cell && bin_member(&kills, &bc.barcode) {
                bc.now_a_cell = false;
            }
        }

        // OH demux filter
        if let (Some(count_chems), Some(fingerprint)) =
            (args.count_chemistry_defs, args.sample_fingerprint)
        {
            if let Some(overhang_read_barcode) = count_chems.overhang_read_barcode() {
                let valid_overhang_ids: HashSet<BarcodeId> = fingerprint
                    .read()?
                    .iter()
                    .flat_map(Fingerprint::tag_names)
                    .map(|tag_name| BarcodeId::pack(tag_name))
                    .collect();
                overhang_demux_filter(
                    overhang_read_barcode,
                    valid_overhang_ids,
                    &mut barcode_cell_info,
                    Some(&mut filter_logger),
                );
            }
        }

        let mut asm_cell_barcodes = BTreeSet::new();
        for x in &barcode_cell_info {
            if x.now_a_cell {
                asm_cell_barcodes.insert(x.barcode.clone());
            }
        }

        // In GD mode, take a pass to identify cells which have at least one productive G/D chain
        let mut gd_barcodes = HashSet::new();
        if gd_mode {
            for can in args.contig_annotations.lazy_reader()? {
                let can: ContigAnnotation = can?;
                if can.productive.unwrap_or(false)
                    && can
                        .annotations
                        .iter()
                        .any(|ann| matches!(ann.feature.chain, VdjChain::TRG | VdjChain::TRD))
                {
                    gd_barcodes.insert(can.barcode);
                }
            }
        }

        // Update contig annotation json file:
        // 1. Adjust high_confidence and is_cell values.
        // 2. Filter unrelated chains based on Gamma/Delta mode
        //      -> in GD mode, filter all AB (NOT VICE VERSA)
        //      -> in AB mode, no filter is applied to maintain consistency

        let contig_annotations_file: JsonFile<_> = rover.make_path("contig_annotations");
        let mut ann_writer = contig_annotations_file.lazy_writer()?;

        let reader = args.contig_annotations.lazy_reader()?;
        for can in reader {
            let mut can: ContigAnnotation = can?;
            can.high_confidence = confident_bcs.contains(&can.barcode);
            if bin_member(&killsc, &can.contig_name) || bin_member(&kills, &can.barcode) {
                can.high_confidence = false;
            }
            can.is_cell =
                asm_cell_barcodes.contains(&can.barcode) && !bin_member(&kills, &can.barcode);
            if gd_mode {
                can.is_cell = can.is_cell && gd_barcodes.contains(&can.barcode);
                can.productive = can.productive.map(|prod| {
                    prod && can
                        .annotations
                        .iter()
                        .any(|ann| matches!(ann.feature.chain, VdjChain::TRG | VdjChain::TRD))
                });
            }
            ann_writer.write_item(&can)?;
        }
        ann_writer.finish()?;

        Ok(AsmCallCellsStageOutputs {
            contig_annotations: contig_annotations_file,
            filter_diagnostics: filter_diagnostics_file,
        })
    }
}
