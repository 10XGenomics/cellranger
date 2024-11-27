//! AsmMetrics stage code

use crate::assembly::{BarcodeDataFile, VdjPrimers};
use anyhow::Result;
use cr_types::chemistry::ChemistryDef;
use cr_types::MetricsFile;
use io_utils::read_to_string_safe;
use martian::prelude::*;
use martian_derive::{make_mro, martian_filetype, MartianStruct};
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::{FileTypeRead, LazyFileTypeIO};
use metric::JsonReporter;
use serde::{Deserialize, Serialize};
use std::path::PathBuf;
use vdj_ann::refx::{make_vdj_ref_data_core, RefData};
use vdj_asm_utils::barcode_data::{metrics_json, BarcodeData, BarcodeDataSum};
use vdj_reference::VdjReceptor;

martian_filetype!(TxtFile, "txt");

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct AsmMetricsStageInputs {
    pub chemistry_def: ChemistryDef,
    pub vdj_reference_path: Option<PathBuf>,
    pub receptor: Option<VdjReceptor>,
    pub inner_enrichment_primers: Option<PathBuf>,
    pub total_read_pairs: Option<i64>,
    pub barcode_full: BarcodeDataFile,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct AsmMetricsStageOutputs {
    pub metrics_summary_json: MetricsFile,
    pub report: TxtFile,
}

// This is our stage struct
pub struct AsmMetrics;

#[make_mro(mem_gb = 8)]
impl MartianMain for AsmMetrics {
    type StageInputs = AsmMetricsStageInputs;
    type StageOutputs = AsmMetricsStageOutputs;
    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        let single_end = !args.chemistry_def.is_paired_end();
        let is_tcr =
            args.receptor == Some(VdjReceptor::TR) || args.receptor == Some(VdjReceptor::TRGD);
        let is_bcr = args.receptor == Some(VdjReceptor::IG);
        let mut refdata = RefData::new();
        if let Some(ref ref_path) = args.vdj_reference_path {
            let fasta = read_to_string_safe(format!("{}/fasta/regions.fa", ref_path.display()));
            make_vdj_ref_data_core(&mut refdata, &fasta, "", is_tcr, is_bcr, None);
        };
        let primers = VdjPrimers::new(
            args.vdj_reference_path.clone(),
            args.inner_enrichment_primers,
            args.receptor,
        )?;
        let bc_data: Vec<BarcodeData> = args
            .barcode_full
            .lazy_reader()?
            .map(std::result::Result::unwrap)
            .collect();
        let bc_data_sum = BarcodeDataSum::sum(&bc_data, &refdata);
        let report_file: TxtFile = rover.make_path("report");
        let mut report = report_file.buf_writer()?;
        let metrics_reporter = metrics_json(
            &bc_data_sum,
            single_end,
            &refdata,
            &primers.inner_primers,
            &primers.outer_primers,
            args.total_read_pairs,
            &mut report,
        );
        let reference_reporter = if let Some(ref_path) = args.vdj_reference_path {
            let report: JsonReporter = JsonFile::new(ref_path, "reference").read()?;
            report.add_prefix("vdj_reference")
        } else {
            JsonReporter::default()
        };
        let metrics_summary_json = MetricsFile::from_reporter(
            &rover,
            "metrics_summary_json",
            &(metrics_reporter + reference_reporter),
        )?;

        Ok(AsmMetricsStageOutputs {
            metrics_summary_json,
            report: report_file,
        })
    }
}
