// Copyright (c) 2022 10X Genomics, Inc. All rights reserved.

use crate::proc_args3::{get_path_fail, proc_meta, proc_xcr};
use enclone_core::defs::{EncloneControl, InputSpec};
use string_utils::TextUtils;
use vector_utils::{next_diff, unique_sort};

/// Process input source.
pub fn proc_input_source(ctl: &mut EncloneControl) -> Result<(), String> {
    match &ctl.cr_opt.input {
        None => {
            return Err(
                "No input data source provided; you must provide one of TCR, BCR, TCRGD, or META."
                    .to_string(),
            );
        }
        Some(InputSpec::MetaFile(path)) => {
            let f = get_path_fail(path, ctl, "META")?;
            if f.contains('/') {
                let d = f.rev_before("/").to_string();
                if !ctl.cr_opt.pre.contains(&d) {
                    ctl.cr_opt.pre.push(d);
                }
            }

            let (receptor, origin_info) = proc_meta(&f, ctl)?;
            ctl.origin_info = origin_info;
            ctl.gen_opt.receptor = receptor;
        }
        Some(InputSpec::Explicit(receptor, arg)) => {
            ctl.origin_info = proc_xcr(*receptor, arg, ctl)?;
            ctl.gen_opt.receptor = *receptor;
        }
    }

    let mut donors = Vec::<String>::new();
    for i in 0..ctl.origin_info.n() {
        donors.push(ctl.origin_info.donor_id[i].clone());
    }
    unique_sort(&mut donors);
    ctl.origin_info.donor_list = donors;
    Ok(())
}

pub fn validate_opts(ctl: &EncloneControl) -> Result<(), String> {
    // Check for duplicated directory paths.

    let mut dp = ctl.origin_info.dataset_path.clone();
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
