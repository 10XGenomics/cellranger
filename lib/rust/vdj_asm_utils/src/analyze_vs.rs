// The function analyze_vs computes an error rate for contigs.  See the documentation
// that it prints for more information.  Also vdj_error_rate is a front end for it.
//
// Compare seg.rs, which does something similar.

#![allow(clippy::many_single_char_names)]
// TODO: fix these.
#![allow(clippy::needless_range_loop)]

use crate::barcode_data::{write_json_metric_f64, write_json_metric_ratio};
use align_tools::{affine_align, summary_less};
use debruijn::dna_string::DnaString;
use debruijn::Mer;
use io_utils::{fwriteln, open_for_read};
use itertools::Itertools;
use rayon::prelude::*;
use serde_json::Value;
use stats_utils::percent_ratio;
use std::collections::HashMap;
use std::fmt::Write as _;
use std::io::{BufRead, Write};
use std::path::Path;
use string_utils::{abbrev_list, strme, TextUtils};
use vdj_ann::annotate::annotate_seq;
use vdj_ann::refx::RefData;
use vdj_ann::transcript::is_valid;
use vector_utils::{bin_member, erase_if, next_diff1_6};

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
// CODE FOR STREAMING A JSON VECTOR
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Read an entry from a json file that represents a vector.  This is not completely
// general as it depends on assumptions about the formatting of the file.

fn read_vector_entry_from_json<R: BufRead>(json: &mut R) -> Option<Vec<u8>> {
    let mut line = String::new();
    if json.read_line(&mut line).is_err() || line == *"" || line == *"[]" {
        return None;
    }
    if line == *"[\n" {
        line.clear();
        assert!(json.read_line(&mut line).is_ok(), "json read failure 1");
    }
    let mut entry = Vec::<u8>::new();
    let (mut curlies, mut bracks, mut quotes) = (0_isize, 0_isize, 0_isize);
    let mut s = line.as_bytes();
    loop {
        if (s == b"]" || s == b"]\n") && curlies == 0 && bracks == 0 && quotes % 2 == 0 {
            return (!entry.is_empty()).then_some(entry);
        }
        let mut cpos = -1_isize;
        for i in (0..s.len() - 1).rev() {
            if s[i] == b',' {
                cpos = i as isize;
                break;
            }
            if s[i] != b' ' {
                break;
            }
        }
        let mut escaped = false;
        for i in 0..s.len() {
            if !escaped && s[i] == b'"' {
                quotes += 1;
            } else if !escaped && quotes % 2 == 0 {
                match s[i] {
                    b'{' => curlies += 1,
                    b'}' => curlies -= 1,
                    b'[' => bracks += 1,
                    b']' => bracks -= 1,
                    b',' => {
                        if i as isize == cpos && curlies == 0 && bracks == 0 && quotes % 2 == 0 {
                            return Some(entry);
                        }
                    }
                    _ => {}
                };
            }
            escaped = s[i] == b'\\' && !escaped;
            entry.push(s[i]);
        }
        line.clear();
        assert!(json.read_line(&mut line).is_ok(), "json read failure 2");
        s = line.as_bytes();
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
// ANALYZE V SEGMENTS
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
#[allow(clippy::too_many_arguments)]
pub fn analyze_vs(
    o: Vec<&str>,        // outs directories or lena ids
    species: &str,       // species (human or mouse or unknown)
    refdata: &RefData,   // reference data
    use_knowns: bool,    // use known alleles
    from_scratch: bool,  // compute annotations from scratch
    parallel: bool,      // parallelize
    do_aff: bool,        // show affine alignments of erroneous V segments (slow)
    context: bool,       // print context; also turns on json metrics
    novels: bool,        // print novel alleles
    json: &mut Vec<u8>,  // dump metrics here
    log: &mut Vec<u8>,   // dump logging here
    is_gd: Option<bool>, // is gamma/delta mode
) {
    let mut outs = Vec::<String>::new();
    let mut lenas = Vec::<String>::new();
    for i in 0..o.len() {
        lenas.push(o[i].to_owned());
        outs.push(o[i].to_owned());
    }
    for i in 0..lenas.len() {
        if lenas[i].contains('_') {
            lenas[i] = lenas[i].before("_").to_string();
        }
    }

    // Note CellRanger version.

    let version = "3.1";

    // Get known alleles.

    let mut knowns = Vec::<DnaString>::new();
    if use_knowns {
        let araw = include_str!("v_alleles.fasta");
        let fields: Vec<&str> = araw.split('\n').collect();
        for i in (1..fields.len()).step_by(2) {
            knowns.push(DnaString::from_dna_string(fields[i]));
        }
        knowns.sort();
    }

    // Find productive cells and the contigs associated with them.

    let mut tignames = Vec::<Vec<String>>::new();
    // contig seq, barcode, quals, tig_id, tig_start
    let mut tigs = Vec::<Vec<(DnaString, String, Vec<u8>, isize, isize)>>::new();
    for i in 0..outs.len() {
        //
        // From all_contig_annotations.json, find V segments that start at zero on the contig.
        // For these we record the feature id and the contig start position.

        let mut tigname_to_feature_id = HashMap::<String, isize>::new();
        let mut tigname_to_vstart = HashMap::<String, isize>::new();
        let f = {
            let base: &Path = outs[i].as_ref();
            let f = base.join("all_contig_annotations.json");
            if f.exists() {
                f
            } else {
                base.join("contig_annotations.json")
            }
        };
        let mut f = open_for_read![&f];
        let mut last_bc = String::default();
        let mut this_tignames = Vec::<String>::new();
        let mut tigname_to_id = HashMap::<String, usize>::new();
        let (mut n1, mut n2) = (0, 0);
        let mut qv = Vec::<Vec<u8>>::new();
        let mut dv = Vec::<DnaString>::new();
        let mut headers = Vec::<String>::new();
        loop {
            match read_vector_entry_from_json(&mut f) {
                None => break,
                Some(x) => {
                    let v: Value = serde_json::from_str(strme(&x)).unwrap();
                    let barcode = v["barcode"].to_string();
                    if barcode != last_bc {
                        if n1 >= 1 && n2 >= 1 {
                            for j in 0..this_tignames.len() {
                                tigname_to_id.insert(this_tignames[j].clone(), tignames.len());
                            }
                            tignames.push(this_tignames.clone());
                            tigs.push(Vec::<(DnaString, String, Vec<u8>, isize, isize)>::new());
                        }
                        this_tignames.clear();
                        n1 = 0;
                        n2 = 0;
                        last_bc = barcode.to_string();
                    }
                    let tigname = &v["contig_name"].to_string().between("\"", "\"").to_string();
                    let ann = v["annotations"].as_array().unwrap();
                    if v["productive"].as_bool().unwrap_or(false)
                        && v["is_cell"].as_bool().unwrap()
                        && v["high_confidence"].as_bool().unwrap()
                    {
                        this_tignames.push(tigname.clone());
                        let mut class = 0;
                        for i in 0..ann.len() {
                            let a = &ann[i];
                            let chain = a["feature"]["chain"]
                                .to_string()
                                .between("\"", "\"")
                                .to_string();
                            if chain == "TRA" || chain == "IGH" {
                                class = 1;
                            }
                            if chain == "TRB" || chain == "IGK" || chain == "IGL" {
                                class = 2;
                            }
                            let id = a["feature"]["feature_id"].as_i64().unwrap() as isize;
                            let region_type = &a["feature"]["region_type"];
                            let ref_start = a["annotation_match_start"].as_u64().unwrap() as usize;
                            let tig_start = a["contig_match_start"].as_i64().unwrap() as isize;
                            if region_type == "L-REGION+V-REGION" && ref_start == 0 {
                                tigname_to_vstart.insert(tigname.clone(), tig_start);
                                tigname_to_feature_id.insert(tigname.clone(), id);
                            }
                        }
                        if class == 1 {
                            n1 += 1;
                        } else if class == 2 {
                            n2 += 1;
                        }
                    }
                    let mut s = v["quals"].to_string().as_bytes().to_vec();
                    for i in 0..s.len() {
                        s[i] -= 33;
                    }
                    qv.push(s);
                    headers.push(v["contig_name"].to_string());
                    dv.push(DnaString::from_dna_string(&v["sequence"].to_string()));
                }
            }
        }
        if n1 >= 1 && n2 >= 1 {
            for j in 0..this_tignames.len() {
                tigname_to_id.insert(this_tignames[j].clone(), tignames.len());
            }
            tignames.push(this_tignames.clone());
            tigs.push(Vec::<(DnaString, String, Vec<u8>, isize, isize)>::new());
        }

        // Collect info.

        for j in 0..dv.len() {
            let tigname = &headers[j];
            if tigname_to_id.contains_key(tigname) {
                let id = tigname_to_id[tigname];
                let mut feature_id = -1_isize;
                let mut vstart = -1_isize;
                if tigname_to_feature_id.contains_key(tigname) {
                    feature_id = tigname_to_feature_id[tigname];
                }
                if tigname_to_vstart.contains_key(tigname) {
                    vstart = tigname_to_vstart[tigname];
                }
                tigs[id].push((
                    dv[j].clone(),
                    tigname.before("_").to_string(),
                    qv[j].clone(),
                    feature_id,
                    vstart,
                ));
            }
        }
    }

    // Flatten tigs to single-level vector and remove duplicate (contig,barcode) entries.

    let mut tigsx = Vec::<(DnaString, String, Vec<u8>, isize, isize)>::new();
    for i in 0..tigs.len() {
        tigsx.append(&mut tigs[i]);
    }
    tigsx.sort();
    let mut to_delete = vec![false; tigsx.len()];
    let mut i = 0;
    while i < tigsx.len() {
        // "let j = next_diff12_5(&tigsx, i)" but doesn't exist
        let mut j = i + 1;
        while j < tigsx.len() {
            if tigsx[j].0 != tigsx[i].0 || tigsx[j].1 != tigsx[i].1 {
                break;
            }
            j += 1;
        }
        for k in i + 1..j {
            to_delete[k] = true;
        }
        i = j;
    }
    erase_if(&mut tigsx, &to_delete);

    // Set up parallel "region".  This allows us to control the number of threads in a
    // parallel loop, without affecting parallism outside the region.

    let mut nt = 0;
    if !parallel {
        nt = 1;
    }
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(nt)
        .build()
        .unwrap();
    pool.install(|| {
        // Find the V regions.

        const VTRIM: usize = 30;
        // {(index, using, contig sequence, mismatches, error report, mismatch report, t, q60)}
        let mut rs = Vec::with_capacity(tigsx.len());
        for i in 0..tigsx.len() {
            rs.push((
                i,
                false,
                DnaString::new(),
                0,
                String::new(),
                String::new(),
                0,
                false,
            ));
        }
        let mut to_t = HashMap::<usize, usize>::new();
        for t in 0..refdata.refs.len() {
            to_t.insert(refdata.id[t] as usize, t);
        }
        #[allow(clippy::too_many_arguments)]
        fn make_report(
            refdata: &RefData,
            t: usize,
            l: usize,
            n: usize,
            i: usize,
            tig: &DnaString,
            do_aff: bool,
            tigsx: &[(DnaString, String, Vec<u8>, isize, isize)],
            r: &mut (usize, bool, DnaString, usize, String, String, usize, bool),
        ) {
            let mut mismatches = 0;
            let mut mis_report = refdata.rheaders[t]
                .after("|")
                .rev_before("|")
                .rev_before("|")
                .to_string();
            for s in 0..n - VTRIM {
                if tig.get(l + s) != refdata.refs[t].get(s) {
                    mismatches += 1;
                    write!(
                        mis_report,
                        " + {}{}",
                        s + 1,
                        tig.to_ascii_vec()[l + s] as char
                    )
                    .unwrap();
                }
            }
            let mut q60 = true;
            for m in l..l + n - VTRIM {
                if tigsx[i].2[m] < 60 {
                    q60 = false;
                }
            }
            r.1 = true;
            r.6 = t;
            r.7 = q60;
            if mismatches > 0 {
                let sub = tig.slice(l, l + n - VTRIM).to_owned();
                if do_aff {
                    let a = affine_align(&sub, &refdata.refs[t]);
                    let report = summary_less(&a);
                    r.4 = report;
                }
                r.2 = sub;
                r.3 = mismatches;
                r.5 = mis_report;
            }
        }
        rs.par_iter_mut().for_each(|r| {
            let i = r.0;
            let tig = &tigsx[i].0;
            let feature_id = tigsx[i].3;
            let vstart = tigsx[i].4;
            if !from_scratch {
                if feature_id >= 0 {
                    let t = to_t[&(feature_id as usize)];
                    let l = vstart as usize;
                    let n = refdata.refs[t].len();
                    assert!(
                        n >= 200,
                        "Reference V segment appears to have length {n}, which doesn't \
                             make sense.  Probably this function was called without the \
                             from_scratch option even though the reference has changed, or \
                             perhaps you supplied the wrong species."
                    );
                    if l + n <= tig.len() {
                        make_report(refdata, t, l, n, i, tig, do_aff, &tigsx, r);
                    }
                }
            } else {
                let mut ann = Vec::<(i32, i32, i32, i32, i32)>::new();
                annotate_seq(tig, refdata, &mut ann, true, false, true);
                let mut log = Vec::<u8>::new();
                if is_valid(tig, refdata, &ann, false, &mut log, is_gd) {
                    for l in 0..ann.len() {
                        let t = ann[l].2 as usize;
                        if refdata.is_v(t) && ann[l].3 == 0 {
                            let l = ann[l].0 as usize;
                            let n = refdata.refs[t].len();
                            if l + n <= tig.len() {
                                make_report(refdata, t, l, n, i, tig, do_aff, &tigsx, r);
                            }
                        }
                    }
                }
            }
        });

        let mut tiglets = Vec::<(DnaString, usize, String, String, usize, bool)>::new();
        let (mut total, mut total_q60) = (0, 0);
        let mut totals = vec![0; refdata.refs.len()];
        for i in 0..tigsx.len() {
            if rs[i].1 {
                total += 1;
                if rs[i].7 {
                    totals[rs[i].6] += 1;
                    total_q60 += 1;
                }
                if !rs[i].2.is_empty() {
                    tiglets.push((
                        rs[i].2.clone(),
                        rs[i].3,
                        rs[i].4.clone(),
                        rs[i].5.clone(),
                        rs[i].6,
                        rs[i].7, // q60
                    ));
                }
            }
        }
        tiglets.sort();

        // Look for common segments that are Q60.  We don't show known alleles, although
        // people might want to see them.

        let (mut common, mut common_q60) = (0, 0);
        const MIN_FREQ: usize = 5;
        const MIN_FRAC: f64 = 0.005;
        let mut i = 0;
        let mut common_clones = 0;
        let (mut errs0, mut errs0_q60) = (Vec::<usize>::new(), Vec::<usize>::new());
        let (mut errs, mut errs_q60) = (Vec::<String>::new(), Vec::<String>::new());
        let mut commons = Vec::<(String, String, String, DnaString)>::new();
        let mut ntiglets_q60 = 0;
        while i < tiglets.len() {
            let j = next_diff1_6(&tiglets, i as i32) as usize;
            let denom = totals[tiglets[i].4];
            let mut count = 0; // number that are q60
            for k in i..j {
                if tiglets[k].5 {
                    count += 1;
                    ntiglets_q60 += 1;
                }
            }
            if bin_member(&knowns, &tiglets[i].0)
                || (count >= MIN_FREQ && count as f64 / denom as f64 >= MIN_FRAC)
            {
                common += j - i;
                common_q60 += count;
                common_clones += 1;
                let descrip = format!(
                    "{}|mult={}/{}|lena={}",
                    tiglets[i].3,
                    count,
                    denom,
                    lenas.iter().format(",")
                );
                let gene_number = tiglets[i].3.before("|").to_string();
                let mut gene_name = tiglets[i].3.after("|").to_string();
                if gene_name.contains(' ') {
                    gene_name = gene_name.before(" ").to_string();
                }
                if !bin_member(&knowns, &tiglets[i].0) {
                    commons.push((gene_name, gene_number, descrip, tiglets[i].0.clone()));
                }
            } else {
                for k in i..j {
                    errs0.push(tiglets[k].1);
                    errs.push(tiglets[k].2.clone());
                    if tiglets[k].5 {
                        errs0_q60.push(tiglets[k].1);
                        errs_q60.push(tiglets[k].2.clone());
                    }
                }
            }
            i = j;
        }
        commons.sort();
        errs0.sort_unstable();
        errs.sort();
        errs0_q60.sort_unstable();
        errs_q60.sort();
        if context {
            fwriteln!(
                log,
                "\n▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓"
            );
            fwriteln!(log, "ERROR RATE ANALYSIS");
            fwriteln!(
                log,
                "▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓\n"
            );
            fwriteln!(
                log,
                "This analysis defines an error rate for contigs, which approximates the\n\
                 fraction of contigs that contain an error.  This is not a per-base error rate!\n\
                 THIS ANALYSIS IS ONLY ACCURATE FOR TCR DATA FROM HUMAN OR MOUSE B6 OR BALB/C.\n\
                 \n\
                 We only use contigs that appear in productive pairs.  Rather than use all the\n\
                 bases lying in contigs, we use only the bases lying in V segments, and trim off\n\
                 the last {} bases, which might be involved in recombination.  The remaining\n\
                 bases should match the germline sequence from which the segment originated.\n\
                 \n\
                 We look for exact matches either to the VDJ reference sequence or to a list of\n\
                 known alternate alleles.  All those that do not match are treated as errors.\n\
                 \n\
                 The error rate is computed as the number of mismatching segments, divided by\n\
                 the total number of segments.\n\
                 \n\
                 We also show the number of mismatches (as \"errs\"), in an abbreviated form,\n\
                 with exponents indicating multiplicity.  For example 1^10 refers to 10 V\n\
                 segments, each having one mismatch with the best reference.  Note that in some\n\
                 cases a large number of mismatches might be better represented as a small indel.\n\
                 \n\
                 Statistics are computed twice, once using all V segments, and a second time\n\
                 using only those whose bases all have quality 60 or higher.  For the most part\n\
                 these are the V segments that are fully covered by reads from at least two UMIs.\n\
                 \n\
                 When data from a new human donor are generated, new alleles may be found,\n\
                 provided that enough cells are assayed.  These alleles will be identified\n\
                 and printed here, and treated as known, but should also be added to the list of\n\
                 known alleles.",
                VTRIM
            );
        }
        fwriteln!(log, "\nALL V SEGMENTS");
        let mut err_rate = 0 as f64;
        if total > 0 {
            err_rate = percent_ratio(tiglets.len() - common, total);
        }
        fwriteln!(
            log,
            "{:.2}% = error rate = {} / {}",
            err_rate,
            tiglets.len() - common,
            total
        );
        fwriteln!(log, "errs = {}", abbrev_list(&errs0));
        if do_aff {
            fwriteln!(log, "     = {}", abbrev_list(&errs));
        }
        fwriteln!(log, "\nQ60 V SEGMENTS");
        let mut err_rate = 0 as f64;
        if total_q60 > 0 {
            err_rate = percent_ratio(ntiglets_q60 - common_q60, total_q60);
        }
        fwriteln!(
            log,
            "{:.2}% = error rate = {} / {}",
            err_rate,
            ntiglets_q60 - common_q60,
            total_q60
        );
        fwriteln!(log, "errs = {}", abbrev_list(&errs0_q60));
        if do_aff {
            fwriteln!(log, "     = {}", abbrev_list(&errs_q60));
        }
        fwriteln!(log, "\nNON-REFERENCE ALLELES");
        fwriteln!(log, "{} = known alleles", common_clones - commons.len());
        fwriteln!(log, "{} = novel alleles", commons.len());
        if !commons.is_empty() && novels {
            fwriteln!(log, "\nNOVEL ALLELES");
            for i in 0..commons.len() {
                fwriteln!(log, ">|CR{}|{}|{}|", version, species, commons[i].2);
                fwriteln!(log, "{}", commons[i].3.to_string());
            }
        }
        if context {
            fwriteln!(log, "\nMETRICS USED ABOVE");
            fwriteln!(log, "frac_of_v_segments_having_an_error");
            let metric = "frac_of_v_segments_having_an_error";
            write_json_metric_ratio(metric, tiglets.len() - common, total, json);
            fwriteln!(log, "frac_of_q60_v_segments_having_an_error");
            let metric = "frac_of_q60_v_segments_having_an_error";
            write_json_metric_ratio(metric, ntiglets_q60 - common_q60, total_q60, json);
            fwriteln!(log, "number_of_v_segments_having_an_error");
            let metric = "number_of_v_segments_having_an_error";
            write_json_metric_f64(metric, (tiglets.len() - common) as f64, json);
            fwriteln!(log, "number_of_q60_v_segments_having_an_error");
            let metric = "number_of_q60_v_segments_having_an_error";
            write_json_metric_f64(metric, (ntiglets_q60 - common_q60) as f64, json);
        }
        fwriteln!(log, "");
    });
}
