//! enclone_process
// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.
#![deny(missing_docs)]

mod core;
mod disintegrate;
mod doublets;
mod enclone;
mod filter_umi;
mod flag_defective;
mod merge_onesies;
mod populate_features;
mod process;
mod read_json;
mod some_filters;
mod split_candidate_clonotypes;
mod start;
mod weak_chains;

pub use crate::core::barcode_fate::BarcodeFate;
use crate::core::defs::EncloneControl;
pub use crate::core::defs::{ClonotypingConfig, ClonotypingFilterConfig, Dataset, InputSpec};
use anyhow::anyhow;
use enclone::innate::species;
use process::process_clonotypes::{process_clonotypes, write_fate};
use refx::{RefData, make_vdj_ref_data_core};
use start::{load_input_data, main_enclone_start};
use std::fs::File;
use std::io::{BufRead, BufReader};
use vdj_ann::refx;
pub use vdj_types::VdjReceptor;
use vector_utils::next_diff;

/// Run clonotyping.
pub fn run_enclone(cr_opt: ClonotypingConfig) -> anyhow::Result<()> {
    assert!(!cr_opt.dref_file.is_empty());
    assert!(cr_opt.max_cores.is_some());
    assert!(!cr_opt.proto.is_empty());
    assert!(!cr_opt.refname.is_empty());
    assert!(!cr_opt.input.origin_info.is_empty());

    let (ctl, refdata) = setup(cr_opt)?;
    let annotations = load_input_data(&ctl, &refdata)?;
    let (exacts, fate) =
        main_enclone_start(&ctl, annotations, &refdata, false).map_err(|e| anyhow!(e))?;
    write_fate(&ctl.cr_opt.fate_file, &fate).map_err(|e| anyhow!(e))?;
    process_clonotypes(&ctl, &refdata, &exacts)?;
    Ok(())
}

fn setup(cr_opt: ClonotypingConfig) -> anyhow::Result<(EncloneControl, RefData)> {
    let mut ctl = EncloneControl {
        origin_info: cr_opt.input.origin_info.clone(),
        cr_opt,
        trace_barcode: String::new(),
        species: String::new(),
    };

    validate_opts(&ctl).map_err(|e| anyhow!(e))?;

    if let Some(max_cores) = ctl.cr_opt.max_cores {
        let _ = rayon::ThreadPoolBuilder::new()
            .num_threads(max_cores)
            .build_global();
    }

    // Determine the reference sequence that is to be used.

    let mut refx = String::new();
    let fx = File::open(&ctl.cr_opt.refname);
    let f = BufReader::new(fx.unwrap());
    for line in f.lines() {
        let s = line.unwrap();
        refx += &s;
        refx += "\n";
    }

    // Build reference data.

    let refx2 = &refx;
    let mut refdata = RefData::new();
    let ext_refx = String::new();
    make_vdj_ref_data_core(
        &mut refdata,
        refx2,
        &ext_refx,
        ctl.is_tcr(),
        ctl.is_bcr(),
        None,
    );

    // Determine if the species is human or mouse or unknown.

    ctl.species = species(&refdata).to_string();

    // Return.

    Ok((ctl, refdata))
}

// TODO: this is assuredly a useless check that duplicates behavior in cellranger
fn validate_opts(ctl: &EncloneControl) -> Result<(), String> {
    // Check for duplicated directory paths.

    let mut dp = ctl
        .origin_info
        .iter()
        .map(|ds| ds.file.as_ref().to_str().unwrap().to_string())
        .collect::<Vec<String>>();
    dp.sort();
    let mut i = 0;
    while i < dp.len() {
        let j = next_diff(&dp, i);
        if j - i > 1 {
            return Err(format!("\nInput dataset path {} is duplicated.\n", dp[i]));
        }
        i = j;
    }

    Ok(())
}
