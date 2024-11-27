// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// This is a special entry point for cellranger, where we know that the arguments that could
// be passed are limited.  The code here is simplified and could be further simplified.

use self::refx::{make_vdj_ref_data_core, RefData};
use anyhow::anyhow;
use enclone::innate::species;
use enclone_args::proc_args_post::{proc_input_source, validate_opts};
use enclone_core::defs::EncloneControl;
pub use enclone_core::defs::{CellrangerFilterOpt, CellrangerOpt, InputSpec};
use enclone_process::process_clonotypes::{process_clonotypes, write_fate};
use enclone_stuff::start::{load_input_data, main_enclone_start};
use std::fs::File;
use std::io::{BufRead, BufReader};
use vdj_ann::refx;
pub use vdj_types::VdjReceptor;

pub fn main_enclone_ranger(cr_opt: CellrangerOpt) -> anyhow::Result<()> {
    assert!(!cr_opt.dref_file.is_empty());
    assert!(cr_opt.max_cores.is_some());
    assert!(cr_opt.pre.is_empty());
    assert!(!cr_opt.proto.is_empty());
    assert!(!cr_opt.refname.is_empty());

    let (ctl, refdata) = main_enclone_setup_ranger(cr_opt)?;
    let annotations = load_input_data(&ctl, &refdata)?;
    let (exacts, fate) =
        main_enclone_start(&ctl, annotations, &refdata, false).map_err(|e| anyhow!(e))?;
    write_fate(&ctl.cr_opt.fate_file, &fate).map_err(|e| anyhow!(e))?;
    process_clonotypes(&ctl, &refdata, &exacts)?;
    Ok(())
}

fn main_enclone_setup_ranger(cr_opt: CellrangerOpt) -> anyhow::Result<(EncloneControl, RefData)> {
    let mut ctl = EncloneControl {
        cr_opt,
        ..Default::default()
    };

    proc_input_source(&mut ctl).map_err(|e| anyhow!(e))?;
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
        ctl.gen_opt.is_tcr(),
        ctl.gen_opt.is_bcr(),
        None,
    );

    // Determine if the species is human or mouse or unknown.

    ctl.gen_opt.species = species(&refdata).to_string();

    // Return.

    Ok((ctl, refdata))
}
