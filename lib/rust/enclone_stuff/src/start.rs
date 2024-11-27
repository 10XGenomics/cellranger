// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
//
// See README for documentation.

use crate::disintegrate::{disintegrate_onesies, DisintegrateOnesiesResult};
use crate::doublets::delete_doublets;
use crate::filter_umi::filter_umi;
use crate::flag_defective::flag_defective;
use crate::merge_onesies::merge_onesies;
use crate::populate_features::populate_features;
use crate::some_filters::{qual_filter, signature_filter};
use crate::split_candidate_clonotypes::split_candidate_clonotypes;
use crate::weak_chains::weak_chains;
use enclone::allele::{find_alleles, sub_alts};
use enclone::graph_filter::GraphFilter;
use enclone::info::build_info;
use enclone::join::join_exacts;
use enclone::misc1::{ArtificialFoursieFilter, CrossFilter};
use enclone::misc2::{check_for_barcode_reuse, find_exact_subclonotypes};
use enclone::misc3::sort_tig_bc;
use enclone::BarcodeFilter;
use enclone_args::read_json::{parse_json_annotations_files, Annotations};
use enclone_core::barcode_fate::BarcodeFate;
use enclone_core::defs::{AltRef, CloneInfo, EncloneControl, OriginInfo};
use enclone_core::enclone_structs::{BarcodeFates, CandidateClonotype, EncloneExacts};
use enclone_core::hcomp::heavy_complexity;
use enclone_process::define_mat::{define_mat, setup_define_mat};
use enclone_process::loupe::make_donor_refs;
use io_utils::fwriteln;
use itertools::Itertools;
use qd::dd;
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// This is a copy of stirling2_ratio_table from the stirling_numbers crate, that has been modified
// to use higher precision internal math.  This has also been speeded up, and in the process
// made less readable.
use qd::Double;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use vdj_ann::refx::RefData;
use vector_utils::{bin_member, erase_if};

pub fn stirling2_ratio_table_double(n_max: usize) -> Vec<Vec<Double>> {
    let mut s = Vec::<Vec<Double>>::new();
    let zero = dd![0.0];
    let one = dd![1.0];
    for n in 0..=n_max {
        s.push(vec![zero; n + 1]);
    }
    s[0][0] = one;
    let mut z = Vec::<Double>::new();
    let mut n2n1 = vec![dd![0.0]; n_max + 1];
    for (n, nn) in n2n1.iter_mut().enumerate().skip(2) {
        *nn = Double::from((n - 2) as u32) / Double::from((n - 1) as u32);
    }
    let mut k1k = vec![dd![0.0]; n_max];
    for (k, kk) in k1k.iter_mut().enumerate().skip(1) {
        *kk = Double::from((k - 1) as u32) / Double::from(k as u32);
    }
    let mut njn = Vec::<(usize, Double)>::new();
    for i in 0..n_max + 1 {
        njn.push((i, dd![0.0]));
    }
    njn.par_iter_mut().for_each(|res| {
        let n = res.0;
        if n >= 1 {
            let mut p = one;
            for j in 1..=n {
                p *= Double::from(j as u32) / Double::from(n as u32);
            }
            res.1 = p;
        }
    });

    // This is the slow part of the function.

    for n in 1..=n_max {
        s[n][0] = zero;
        for k in 1..n - 1 {
            z[k - 1] *= k1k[k];
        }
        if n >= 2 {
            z.push(n2n1[n].powi((n - 1) as i32));
        }
        for k in 1..n {
            let x = z[k - 1]; // = ((k-1)/k)^(n-1)
            s[n][k] = s[n - 1][k] + s[n - 1][k - 1] * x;
        }
        s[n][n] = njn[n].1;
    }
    s
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn load_input_data(ctl: &EncloneControl, refdata: &RefData) -> anyhow::Result<Annotations> {
    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Flag defective reference sequences.

    let mut log = Vec::<u8>::new();
    let mut broken = Vec::<bool>::new();
    flag_defective(ctl, refdata, &mut log, &mut broken);

    // Parse the json annotations file.

    let mut annotations = parse_json_annotations_files(ctl, refdata)?;

    // Populate features.

    let mut fr1_starts = Vec::<usize>::new();
    let mut fr2_starts = Vec::<Option<usize>>::new();
    let mut fr3_starts = Vec::<Option<usize>>::new();
    let mut cdr1_starts = Vec::<Option<usize>>::new();
    let mut cdr2_starts = Vec::<Option<usize>>::new();
    populate_features(
        refdata,
        &mut fr1_starts,
        &mut fr2_starts,
        &mut fr3_starts,
        &mut cdr1_starts,
        &mut cdr2_starts,
    );
    for tigi in &mut annotations.tig_bc {
        for x in tigi {
            x.fr1_start = fr1_starts[x.v_ref_id];
            x.fr2_start = fr2_starts[x.v_ref_id];
            x.fr3_start = fr3_starts[x.v_ref_id];
            x.cdr1_start = cdr1_starts[x.v_ref_id];
            x.cdr2_start = cdr2_starts[x.v_ref_id];
        }
    }

    // Test for no data.

    if ctl.origin_info.n() == 0 {
        anyhow::bail!("\nNo TCR or BCR data have been specified.\n");
    }

    Ok(annotations)
}

pub fn main_enclone_start(
    ctl: &EncloneControl,
    annotations: Annotations,
    refdata: &RefData,
    report_whitelist_contamination: bool,
) -> Result<(EncloneExacts, Vec<BarcodeFates>), String> {
    let Annotations {
        mut tig_bc,
        gex_cells,
        gex_cells_specified,
        mut fate,
    } = annotations;
    // Filter using light --> heavy graph.
    if ctl.cr_opt.filter.graph {
        GraphFilter::default().filter_items(&mut tig_bc, &mut fate);
    }

    // Sort tig_bc.
    sort_tig_bc(&mut tig_bc, refdata, ctl.cr_opt.mix_donors);

    // Cross filter.
    if ctl.cr_opt.filter.cross_dataset {
        CrossFilter {
            origin_info: &ctl.origin_info,
        }
        .filter_items(&mut tig_bc, &mut fate);
    }

    // Look for barcode reuse.
    check_for_barcode_reuse(&ctl.origin_info, &tig_bc)?;

    // Find exact subclonotypes.
    let mut exact_clonotypes = find_exact_subclonotypes(ctl, &tig_bc, refdata);

    if !ctl.gen_opt.trace_barcode.is_empty() {
        for ex in &exact_clonotypes {
            for clone in &ex.clones {
                if clone[0].barcode == ctl.gen_opt.trace_barcode {
                    println!(
                        "\nfound {} in an initial exact subclonotype having {} cells",
                        ctl.gen_opt.trace_barcode,
                        ex.ncells(),
                    );
                }
            }
        }
    }

    if ctl.cr_opt.filter.weak_foursies {
        ArtificialFoursieFilter.filter_items(&mut exact_clonotypes, &mut fate);
    }

    // Build info about clonotypes.  Note that this edits the V reference sequence to perform
    // an indel in some cases.
    // This may filter a barcode if it is found to be "improper".
    let mut info: Vec<CloneInfo> = build_info(refdata, &mut exact_clonotypes[..], &mut fate);

    // Derive consensus sequences for alternate alleles of V segments.  Then create donor
    // reference sequences for Loupe.
    let alt_refs = if ctl.gen_opt.no_alt_alleles {
        Vec::new()
    } else {
        find_alleles(refdata, ctl, &exact_clonotypes)
    };

    if !ctl.cr_opt.dref_file.is_empty() {
        write_donor_ref_file(&ctl.origin_info, &alt_refs, refdata, &ctl.cr_opt.dref_file);
    }
    let drefs = make_donor_refs(&alt_refs, refdata);

    // Update reference sequences for V segments by substituting in alt alleles if better.

    sub_alts(
        refdata,
        ctl,
        &alt_refs,
        &mut info[..],
        &mut exact_clonotypes[..],
    );

    // Compute to_bc, which maps (dataset_index, clonotype_id) to {barcodes}.
    // This is intended as a replacement for some old code below.

    let mut to_bc = HashMap::<(usize, usize), Vec<String>>::new();
    for (i, ex) in exact_clonotypes.iter().enumerate() {
        for clone in &ex.clones {
            let x = &clone[0];
            to_bc
                .entry((x.dataset_index, i))
                .or_default()
                .push(x.barcode.clone());
        }
    }

    // Make stirling ratio table.  Not sure that fixing the size of this is safe.

    let sr = stirling2_ratio_table_double(3000);

    // Compute complexity.

    if ctl.join_alg_opt.comp_filt < 1_000_000 {
        let jun = heavy_complexity(refdata, &exact_clonotypes, ctl, &drefs);
        for u in 0..exact_clonotypes.len() {
            let ex = &mut exact_clonotypes[u];
            for m in 0..ex.share.len() {
                if ex.share.len() == 2 && ex.share[m].left {
                    ex.share[m].jun = jun[u].clone();
                }
            }
        }
    }

    // Form equivalence relation on exact subclonotypes.  We also keep the raw joins, consisting
    // of pairs of info indices, that were originally joined.
    let (eq, raw_joins) = join_exacts(
        &to_bc,
        refdata,
        ctl,
        &exact_clonotypes,
        &info,
        &sr,
        report_whitelist_contamination,
    );

    // Disintegrate certain onesie clonotypes into single cell
    // clonotypes.  This requires editing of exact_clonotypes, info, eq, join_info and raw_joins.
    let DisintegrateOnesiesResult {
        eq,
        mut exact_clonotypes,
        info,
        mut raw_joins,
        disintegrated,
    } = disintegrate_onesies(eq, exact_clonotypes, info, raw_joins);

    // Lock info.
    let info = &info;

    // Update to_bc.

    let mut to_bc = HashMap::<(usize, usize), Vec<String>>::new();
    for (i, ex) in exact_clonotypes.iter().enumerate() {
        for clone in &ex.clones {
            let x = &clone[0];
            to_bc
                .entry((x.dataset_index, i))
                .or_default()
                .push(x.barcode.clone());
        }
    }

    // Restructure raw joins.

    raw_joins.sort_unstable();
    let raw_joins = {
        let mut raw_joins2 = vec![Vec::<usize>::new(); info.len()];
        for r in raw_joins {
            raw_joins2[r.0 as usize].push(r.1 as usize);
            raw_joins2[r.1 as usize].push(r.0 as usize);
        }
        raw_joins2
    };

    if !ctl.gen_opt.trace_barcode.is_empty() {
        for ex in &exact_clonotypes {
            for clone in &ex.clones {
                if clone[0].barcode == ctl.gen_opt.trace_barcode {
                    println!(
                        "\nfound {} in a pre-filter exact subclonotype having {} cells",
                        ctl.gen_opt.trace_barcode,
                        ex.ncells(),
                    );
                }
            }
        }
    }

    // Create candidate clonotypes.
    let mut candidate_clonotypes = eq.all_sets_iter_i32().collect();

    // Filter B cells based on UMI counts.
    if ctl.gen_opt.is_bcr() && (ctl.cr_opt.filter.umi_count || ctl.cr_opt.filter.umi_ratio) {
        // Conditionals on UMI filtering are processed internally.
        // We could lift them out if we delete the "clonotype marking" feature.
        candidate_clonotypes = filter_umi(
            ctl,
            &mut exact_clonotypes,
            info,
            candidate_clonotypes,
            &mut fate,
        );
    }

    if !ctl.gen_opt.trace_barcode.is_empty() {
        for ex in &exact_clonotypes {
            for clone in &ex.clones {
                if clone[0].barcode == ctl.gen_opt.trace_barcode {
                    println!(
                        "\nfound {} in an post-umi-filter exact subclonotype having {} cells",
                        ctl.gen_opt.trace_barcode,
                        ex.ncells(),
                    );
                }
            }
        }
    }

    // Remove cells that are not called cells by GEX.
    let candidate_clonotypes = {
        let mut candidate_clonotypes2 = Vec::<Vec<i32>>::new();
        for mut cc in candidate_clonotypes {
            let mut to_deletex = vec![false; cc.len()];
            for (&x, dx) in cc.iter().zip(to_deletex.iter_mut()) {
                let x: &CloneInfo = &info[x as usize];
                let ex = &mut exact_clonotypes[x.clonotype_index];
                let mut to_delete = vec![false; ex.ncells()];
                for (clone, d) in ex.clones.iter().take(ex.ncells()).zip(to_delete.iter_mut()) {
                    let li = clone[0].dataset_index;
                    let bc = &clone[0].barcode;
                    if gex_cells_specified[li] && !bin_member(&gex_cells[li], bc) {
                        *d = true;
                        fate[li].insert(bc.clone(), BarcodeFate::NotGexCell);
                    }
                }
                erase_if(&mut ex.clones, &to_delete);
                if ex.ncells() == 0 {
                    *dx = true;
                }
            }
            erase_if(&mut cc, &to_deletex);
            if !cc.is_empty() {
                candidate_clonotypes2.push(cc);
            }
        }
        candidate_clonotypes2
    };

    // Break up clonotypes containing a large number of chains. These are
    // very likely to be false merges
    let mut candidate_clonotypes: Vec<CandidateClonotype> = candidate_clonotypes
        .into_iter()
        .flat_map(|candidate_clonotype| {
            let (od, exacts) = setup_define_mat(&candidate_clonotype, info);
            let mat = define_mat(
                &to_bc,
                &sr,
                ctl,
                &exact_clonotypes,
                &exacts,
                &od,
                info,
                &raw_joins,
                refdata,
            );
            let num_chains = mat.len();
            if num_chains < ctl.cr_opt.split_max_chains {
                vec![candidate_clonotype]
            } else {
                let exacts_of_chains = mat
                    .iter()
                    .enumerate()
                    .flat_map(|(chain_num, chain_in_exact)| {
                        exacts
                            .iter()
                            .zip_eq(chain_in_exact.iter())
                            .filter_map(move |(e, chain)| chain.map(|_| (e, chain_num)))
                    })
                    .into_group_map()
                    .into_iter()
                    .map(|(k, v)| (v, *k))
                    .into_group_map();

                let mut group_of_exacts = HashMap::new();
                let mut group_num = 0;
                for (chains, chain_exacts) in exacts_of_chains.into_iter().sorted() {
                    if chains.len() == 1 {
                        for e in chain_exacts {
                            group_of_exacts.insert(e, group_num);
                            group_num += 1;
                        }
                    } else {
                        for e in chain_exacts {
                            group_of_exacts.insert(e, group_num);
                        }
                        group_num += 1;
                    }
                }

                let mut groups = vec![vec![]; group_num];

                for (_, exact_clonotype_id, val) in &od {
                    groups[group_of_exacts[exact_clonotype_id]].push(*val);
                }

                // To split every subclonotype
                // od
                //     .into_iter()
                //     .group_by(|o| o.1)
                //     .into_iter()
                //     .map(|(_, vals)| vals.map(|v| v.2).collect())
                //     .collect();
                groups
            }
        })
        .collect();

    // Delete exact subclonotypes that appear to represent doublets.
    if ctl.cr_opt.filter.doublet {
        delete_doublets(
            &mut candidate_clonotypes,
            &to_bc,
            &sr,
            ctl,
            &exact_clonotypes,
            info,
            &raw_joins,
            refdata,
            &mut fate,
        );
    }

    if ctl.cr_opt.filter.signature {
        signature_filter(
            &mut candidate_clonotypes,
            &to_bc,
            &sr,
            ctl,
            &exact_clonotypes,
            info,
            &raw_joins,
            &mut fate,
            refdata,
        );
    }

    // Merge onesies where totally unambiguous.
    merge_onesies(
        &mut candidate_clonotypes,
        ctl,
        &exact_clonotypes,
        info,
        &eq,
        &disintegrated,
    );

    // Check for disjoint candidate clonotypes.
    split_candidate_clonotypes(
        &mut candidate_clonotypes,
        &to_bc,
        &sr,
        ctl,
        &exact_clonotypes,
        info,
        &raw_joins,
        refdata,
    );

    // Test for weak chains.
    if ctl.cr_opt.filter.weak_chains {
        weak_chains(
            &mut candidate_clonotypes,
            &to_bc,
            &sr,
            ctl,
            &exact_clonotypes,
            info,
            &raw_joins,
            &mut fate,
            refdata,
        );
    }

    // Check for disjoint candidate clonotypes (again).
    split_candidate_clonotypes(
        &mut candidate_clonotypes,
        &to_bc,
        &sr,
        ctl,
        &exact_clonotypes,
        info,
        &raw_joins,
        refdata,
    );

    if ctl.cr_opt.filter.qual {
        qual_filter(
            &mut candidate_clonotypes,
            &to_bc,
            &sr,
            ctl,
            &exact_clonotypes,
            info,
            &raw_joins,
            &mut fate,
            refdata,
        );
    }

    // Check for disjoint candidate clonotypes (again again).
    split_candidate_clonotypes(
        &mut candidate_clonotypes,
        &to_bc,
        &sr,
        ctl,
        &exact_clonotypes,
        info,
        &raw_joins,
        refdata,
    );

    if !ctl.gen_opt.trace_barcode.is_empty() {
        for ex in &exact_clonotypes {
            for clone in &ex.clones {
                if clone[0].barcode == ctl.gen_opt.trace_barcode {
                    println!(
                        "\nfound {} in an intermediate exact subclonotype having {} cells",
                        ctl.gen_opt.trace_barcode,
                        ex.ncells(),
                    );
                }
            }
        }
    }
    Ok((
        EncloneExacts {
            to_bc,
            exact_clonotypes,
            raw_joins,
            info: info.clone(),
            candidate_clonotypes,
            drefs,
            sr,
        },
        fate,
    ))
}

fn write_donor_ref_file(
    origin_info: &OriginInfo,
    alt_refs: &[AltRef],
    refdata: &RefData,
    dref_file: &str,
) {
    let f = File::create(dref_file);
    if f.is_err() {
        eprintln!("\nError trying to write ctl.cr_opt.dref_file = {dref_file}.");
    }
    let mut f = BufWriter::new(f.unwrap());
    let mut count = 0;
    for i in 0..alt_refs.len() {
        let donor = alt_refs[i].donor;
        let ref_id = alt_refs[i].ref_id;
        if i > 0 && (donor != alt_refs[i - 1].donor || ref_id != alt_refs[i - 1].ref_id) {
            count = 0;
        }
        let alt_seq = &alt_refs[i].alt_seq;
        fwriteln!(
            f,
            ">{}:{}:{}:{} (reference record id : donor name : allele number : gene name)\n{}",
            refdata.id[ref_id],
            origin_info.donor_id[donor],
            count + 1,
            refdata.name[ref_id],
            alt_seq.to_string()
        );
        count += 1;
    }
}
