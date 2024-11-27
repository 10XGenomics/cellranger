// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// This file contains the two functions proc_xcr and proc_meta.

use enclone_core::defs::{EncloneControl, OriginInfo};
use io_utils::{dir_list, path_exists};
use itertools::Itertools;
use rayon::prelude::*;
use std::fmt::Write as _;
use std::fs::File;
use std::io::{BufRead, BufReader};
use string_utils::TextUtils;
use vdj_types::VdjReceptor;
use vector_utils::unique_sort;

// Functions to find the path to data.

pub fn get_path_fail(p: &str, ctl: &EncloneControl, source: &str) -> Result<String, String> {
    for x in &ctl.cr_opt.pre {
        let pp = format!("{x}/{p}");
        if path_exists(&pp) {
            return Ok(pp);
        }
    }
    if !path_exists(p) {
        if ctl.cr_opt.pre.is_empty() {
            let path = std::env::current_dir().unwrap();
            return Err(format!(
                "\nIn directory {}, unable to find the path {}.  This came from the {} argument.\n",
                path.display(),
                p,
                source
            ));
        }
        let path = std::env::current_dir().unwrap();
        let mut pre_msg = "Here are the number of entries in your PRE directories:\n".to_string();
        for x in &ctl.cr_opt.pre {
            let mut count = "(does not exist)".to_string();
            if path_exists(x) {
                count = dir_list(x).len().to_string();
            }
            writeln!(pre_msg, "{x}: {count}").unwrap();
        }
        return Err(format!(
            "\nIn directory {}, unable to find the\npath {},\n\
                even if prepended by any of the directories \
                in\nPRE={}.\nThis came from the {} argument.\n{}",
            path.display(),
            p,
            ctl.cr_opt.pre.iter().format(","),
            source,
            pre_msg
        ));
    }
    Ok(p.to_string())
}

fn get_path(p: &str, ctl: &EncloneControl, ok: &mut bool) -> String {
    *ok = false;
    for x in &ctl.cr_opt.pre {
        let pp = format!("{x}/{p}");
        if path_exists(&pp) {
            *ok = true;
            return pp;
        }
    }
    let pp = p.to_string();
    *ok = path_exists(&pp);
    pp
}

fn get_path_or_internal_id(p: &str, ctl: &EncloneControl, source: &str) -> Result<String, String> {
    let mut ok = false;
    let mut pp = get_path(p, ctl, &mut ok);
    if !ok {
        get_path_fail(&pp, ctl, source)?;
    }
    if !pp.ends_with("/outs") && path_exists(format!("{pp}/outs")) {
        pp = format!("{pp}/outs");
    }
    Ok(pp)
}

pub fn proc_xcr(
    receptor: VdjReceptor,
    val: &str,
    ctl: &EncloneControl,
) -> Result<OriginInfo, String> {
    let xcr = match receptor {
        VdjReceptor::IG => "BCR",
        VdjReceptor::TR => "TCR",
        VdjReceptor::TRGD => "TCRGD",
    };
    assert!(!val.is_empty());

    let donor_groups = vec![val];

    let mut origin_info = OriginInfo::default();

    for (id, d) in donor_groups.iter().enumerate() {
        let origin_groups = [&d[..]];

        for (is, s) in origin_groups.iter().enumerate() {
            let mut datasets = vec![&s[..]];
            for ds in &mut datasets {
                if ds.ends_with('/') {
                    *ds = ds.rev_before("/");
                }
            }
            for x in &datasets {
                let donor_name = format!("d{}", id + 1);
                let origin_name = format!("s{}", is + 1);
                origin_info.donor_id.push(donor_name);
                origin_info.origin_id.push(origin_name);
                let mut dataset_name = (*x).to_string();
                if dataset_name.contains('/') {
                    dataset_name = dataset_name.rev_after("/").to_string();
                }
                origin_info.dataset_id.push(dataset_name.clone());
            }
        }
    }

    // Get paths.  This will need to change when cellranger switches to multi.  This code is
    // parallelized because this code can indirectly make many calls to path_exists, and the wall
    // clock time for these can add up.  There should be a way to do this that does not involve
    // multithreading.

    let mut results = Vec::<(String, bool, String)>::new();
    for d in donor_groups {
        let origin_groups = (*d).split(':').collect::<Vec<&str>>();
        for s in origin_groups {
            let datasets = s.split(',').collect::<Vec<&str>>();
            for x in datasets {
                results.push((x.to_string(), false, String::new()));
            }
        }
    }

    results.par_iter_mut().for_each(|res| {
        let p = &mut res.0;
        let resx = get_path_or_internal_id(p, ctl, xcr);
        match resx {
            Err(resx) => res.2 = resx,
            Ok(resx) => {
                *p = resx;
                match receptor {
                    VdjReceptor::IG => {
                        if path_exists(format!("{p}/vdj_b")) {
                            *p = format!("{p}/vdj_b");
                        }
                        if path_exists(format!("{p}/multi/vdj_b")) {
                            *p = format!("{p}/multi/vdj_b");
                        }
                    }
                    VdjReceptor::TR => {
                        if path_exists(format!("{p}/vdj_t")) {
                            *p = format!("{p}/vdj_t");
                        }
                        if path_exists(format!("{p}/multi/vdj_t")) {
                            *p = format!("{p}/multi/vdj_t");
                        }
                    }
                    VdjReceptor::TRGD => (),
                }
            }
        }
    });
    for result in results {
        if !result.2.is_empty() {
            return Err(result.2);
        }
        origin_info.dataset_path.push(result.0);
    }

    Ok(origin_info)
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

fn proc_meta_core(
    lines: &[String],
    ctl: &EncloneControl,
) -> Result<(VdjReceptor, OriginInfo), String> {
    let fields: Vec<_> = lines
        .first()
        .expect("META file is empty")
        .split(',')
        .collect();
    let mut fields_sorted = fields.clone();
    unique_sort(&mut fields_sorted);
    if fields_sorted.len() < fields.len() {
        return Err("\nThe CSV file that you specified using the META argument \
                     has duplicate field names\nin its first line.\n"
            .to_string());
    }
    let allowed_fields = ["bcr", "donor", "origin", "tcr", "tcrgd", "color"];
    for x in &fields {
        if !allowed_fields.contains(x) {
            return Err(format!(
                "\nThe CSV file that you specified using the META argument \
                         has an illegal field name ({x}) in its first line.\n"
            ));
        }
    }
    let receptor = match (
        fields.contains(&"tcr"),
        fields.contains(&"tcrgd"),
        fields.contains(&"bcr"),
    ) {
        (true, false, false) => VdjReceptor::TR,
        (false, true, false) => VdjReceptor::TRGD,
        (false, false, true) => VdjReceptor::IG,
        (false, false, false) => {
            return Err("\nThe CSV file that you specified using the META argument \
                    has neither the field tcr, tcrgd, or bcr in its first line.\n"
                .to_string());
        }
        _ => {
            return Err("\nThe CSV file that you specified using the META argument \
                    has more than one of the fields tcr, bcr, and tcrgd in its first line.\n"
                .to_string());
        }
    };

    let mut origin_info = OriginInfo::default();

    for (count, s) in lines.iter().enumerate().skip(1) {
        if !s.starts_with('#') && !s.is_empty() {
            let val = s.split(',').collect::<Vec<&str>>();
            if val.len() != fields.len() {
                return Err(format!(
                    "\nMETA file line {} has a different number of fields than the \
                     first line of the file.\n",
                    count + 1
                ));
            }
            let mut path = String::new();
            let mut abbr = String::new();
            let mut origin = "s1".to_string();
            let mut donor = "d1".to_string();
            for i in 0..fields.len() {
                let x = &fields[i];
                let mut y = val[i].to_string();
                if y.starts_with('"') && y.ends_with('"') {
                    y = y.after("\"").rev_before("\"").to_string();
                }
                if *x == "tcr" || *x == "bcr" || *x == "tcrgd" {
                    if y.contains(':') {
                        path = y.after(":").to_string();
                        abbr = y.before(":").to_string();
                    } else {
                        path = y.to_string();
                        if path.contains('/') {
                            abbr = path.rev_after("/").to_string();
                        } else {
                            abbr.clone_from(&path);
                        }
                    }
                } else if *x == "origin" {
                    origin = y.to_string();
                } else if *x == "donor" {
                    donor = y.to_string();
                }
            }

            path = get_path_or_internal_id(&path, ctl, "META")?;
            match receptor {
                VdjReceptor::IG => {
                    if path_exists(format!("{path}/vdj_b")) {
                        path = format!("{path}/vdj_b");
                    }
                    if path_exists(format!("{path}/multi/vdj_b")) {
                        path = format!("{path}/multi/vdj_b");
                    }
                }
                VdjReceptor::TR => {
                    if path_exists(format!("{path}/vdj_t")) {
                        path = format!("{path}/vdj_t");
                    }
                    if path_exists(format!("{path}/multi/vdj_t")) {
                        path = format!("{path}/multi/vdj_t");
                    }
                }
                VdjReceptor::TRGD => {
                    if path_exists(format!("{path}/vdj_t_gd")) {
                        path = format!("{path}/vdj_t_gd");
                    }
                    if path_exists(format!("{path}/multi/vdj_t_gd")) {
                        path = format!("{path}/multi/vdj_t_gd");
                    }
                }
            }
            origin_info.dataset_path.push(path);
            origin_info.dataset_id.push(abbr);
            origin_info.donor_id.push(donor);
            origin_info.origin_id.push(origin);
        }
    }
    Ok((receptor, origin_info))
}

pub fn proc_meta(path: &str, ctl: &EncloneControl) -> Result<(VdjReceptor, OriginInfo), String> {
    let Ok(fx) = File::open(path) else {
        return Err(format!(
            "\nProblem with META: unable to read from the file\n\
                 \"{path}\".\nPlease check that that path makes sense and that you have read \
                 permission for it.\n"
        ));
    };
    let f = BufReader::new(fx);
    let mut lines = Vec::<String>::new();
    for line in f.lines() {
        let s = line.unwrap();
        lines.push(s);
    }
    proc_meta_core(&lines, ctl)
}
