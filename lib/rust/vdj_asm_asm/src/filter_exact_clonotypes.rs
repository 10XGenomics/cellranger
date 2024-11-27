//! FilterExactClonotypes stage code

use anyhow::Result;
use martian::prelude::*;
use martian_derive::{make_mro, MartianStruct};
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::{LazyFileTypeIO, LazyWrite};
use serde::{Deserialize, Serialize};
use std::collections::HashSet;
use vdj_ann::annotate::ContigAnnotation;
use vdj_asm_utils::exact_clonotyping::ExactClonotype;
use vdj_filter_barcodes::filter_clonotype_level::{
    build_wlcontaminfo_per_exact_clonotype, whitelist_contamination_filter,
};
use vdj_filter_barcodes::filter_log::{FilterLogEntry, FilterLogger, VdjFilterLogFormat};

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct FilterExactClonotypesStageInputs {
    pub exact_clonotypes: JsonFile<Vec<ExactClonotype>>,
    pub contig_annotations: JsonFile<Vec<ContigAnnotation>>,
    pub filter_diagnostics: VdjFilterLogFormat,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct FilterExactClonotypesStageOutputs {
    pub contig_annotations: JsonFile<Vec<ContigAnnotation>>,
    #[mro_retain]
    pub filter_diagnostics: VdjFilterLogFormat,
}

pub struct FilterExactClonotypes;

#[make_mro(mem_gb = 4)]
impl MartianMain for FilterExactClonotypes {
    type StageInputs = FilterExactClonotypesStageInputs;
    type StageOutputs = FilterExactClonotypesStageOutputs;
    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        // Set-up filter logging
        let filter_diagnostics_file: VdjFilterLogFormat = rover.make_path("filter_diagnostics");
        let mut filter_logger =
            extend_filter_log(args.filter_diagnostics, &filter_diagnostics_file);

        let exact_clonotypes: Vec<ExactClonotype> = args.exact_clonotypes.read_all()?;
        let contigs: Vec<ContigAnnotation> = args.contig_annotations.read_all()?;

        // Filter out gel bead whitelist contamination
        let whitelist_contam_info =
            build_wlcontaminfo_per_exact_clonotype(&contigs, exact_clonotypes);
        let mut filtered_bcs = HashSet::new();
        for info in whitelist_contam_info {
            let contam_bcs = whitelist_contamination_filter(&info, Some(&mut filter_logger));
            filtered_bcs.extend(contam_bcs);
        }

        // Update cell calls
        let contig_annotations: JsonFile<Vec<ContigAnnotation>> =
            rover.make_path("contig_annotations");
        let mut writer = contig_annotations.lazy_writer()?;
        for ann in contigs {
            writer.write_item(&ContigAnnotation {
                is_cell: if filtered_bcs.contains(&ann.barcode) {
                    false
                } else {
                    ann.is_cell
                },
                is_asm_cell: if filtered_bcs.contains(&ann.barcode) {
                    Some(false)
                } else {
                    ann.is_asm_cell
                },
                ..ann
            })?;
        }

        Ok(FilterExactClonotypesStageOutputs {
            contig_annotations,
            filter_diagnostics: filter_diagnostics_file,
        })
    }
}

fn extend_filter_log(in_fpath: VdjFilterLogFormat, out_fpath: &VdjFilterLogFormat) -> FilterLogger {
    let in_filter_log: Vec<FilterLogEntry> = in_fpath.read_all().unwrap();
    let mut filter_logger = FilterLogger::new(out_fpath).unwrap();
    for entry in in_filter_log {
        filter_logger.log(&entry);
    }
    filter_logger
}
