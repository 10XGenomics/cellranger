// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
//! This file provides the function process_clonotypes, which performs some final
//! filtering and writes out the barcode fate and the loupe clonotype files.
//! Additional behavior is injected by Enclone through the candidate clonotypeProcessor trait.
#![deny(missing_docs)]

use crate::core::defs::{CloneInfo, EncloneControl};
use crate::core::enclone_structs::{BarcodeFates, EncloneExacts};
use crate::process::define_column_info::define_column_info;
use crate::process::define_mat::{Od, define_mat};
use crate::process::delete_weaks::delete_weaks;
use crate::process::loupe::{loupe_out, make_loupe_clonotype};
use anyhow::Result;
use enclone_proto::types::Clonotype;
use itertools::Itertools;
use rayon::prelude::*;
use std::cmp::Reverse;
use std::fs::File;
use std::io::BufWriter;
use vdj_ann::refx::RefData;
use vector_utils::{erase_if, next_diff12_3};

/// Process clonotypes.
/// Filter out exact subclonotypes in candidate clonotypes that appear to be junk.
/// Write out barcode fates and loupe clonotype files.
pub(crate) fn process_clonotypes(
    ctl: &EncloneControl,
    refdata: &RefData,
    enclone_exacts: &EncloneExacts,
) -> Result<()> {
    let EncloneExacts {
        to_bc,
        exact_clonotypes,
        raw_joins,
        info,
        candidate_clonotypes,
        drefs: dref,
        sr,
    } = enclone_exacts;

    // Traverse the candidate clonotypes.
    let result_iter = candidate_clonotypes.par_iter().map(|cc| {
        let od: Vec<_> = cc
            .iter()
            .map(|id| {
                let x: &CloneInfo = &info[*id as usize];
                (x.origin.clone(), x.clonotype_index, *id)
            })
            .sorted()
            .collect();

        // Reconstruct the participating clones.  This is needed because most exact subclonotypes
        // having more than two chains have been split up.
        //
        // Capture these data into parallel data structures, one per exact subclonotype:
        // exacts: the exact subclonotype ids
        // mults:  number of cells [redundant, might remove]

        let mut exacts = Vec::<usize>::new();
        let mut mults = Vec::<usize>::new();
        let mut j = 0;
        while j < od.len() {
            let k = next_diff12_3(&od, j);
            let mut mult = 0_usize;
            for l in j..k {
                let x: &CloneInfo = &info[od[l].2 as usize];
                let m = x.clonotype_index;
                mult = exact_clonotypes[m].clones.len();
            }
            mults.push(mult);
            exacts.push(od[j].1);
            j = k;
        }

        // Identify the exact subclonotypes that are junk.

        sort_exact_clonotypes(ctl, refdata, enclone_exacts, &od, &mut exacts, &mut mults);

        // Define a matrix mat[col][ex] which is the column of the exact subclonotype
        // corresponding to the given column col of the clonotype, which may or may not be
        // defined.  Then define other information associated to each chain.  These are
        // reference sequence identifiers, CDR3 start positions, and the like.
        let mat = define_mat(
            to_bc,
            sr,
            ctl,
            exact_clonotypes,
            &exacts,
            &od,
            info,
            raw_joins,
            refdata,
        );
        let rsi = define_column_info(&exacts, exact_clonotypes, mat);
        let mat = &rsi.mat;

        // Filter.

        // Let n be the total number of cells in this pass.
        let n: usize = mults.iter().sum();

        let mut bads = vec![false; exacts.len()];

        if n > 0 {
            // Mark some weak exact subclonotypes for deletion.
            delete_weaks(&exacts, exact_clonotypes, mat, &mut bads);
        };

        // Delete weak exact subclonotypes.
        erase_if(&mut mults, &bads);
        erase_if(&mut exacts, &bads);

        sort_exact_clonotypes(ctl, refdata, enclone_exacts, &od, &mut exacts, &mut mults);

        // Define a matrix mat[col][ex] which is the column of the exact subclonotype
        // corresponding to the given column col of the clonotype, which may or may not be
        // defined.  Then define other information associated to each chain.  These are
        // reference sequence identifiers, CDR3 start positions, and the like.
        let mat = define_mat(
            to_bc,
            sr,
            ctl,
            exact_clonotypes,
            &exacts,
            &od,
            info,
            raw_joins,
            refdata,
        );
        let rsi = define_column_info(&exacts, exact_clonotypes, mat);

        // Generate Loupe data.

        let loupe_clonotype = (!ctl.cr_opt.proto.is_empty())
            .then(|| make_loupe_clonotype(exact_clonotypes, &exacts, &rsi, refdata, dref));

        // Let n be the total number of cells in this pass.

        let n: usize = mults.iter().sum();

        if n == 0 {
            return (0, loupe_clonotype);
        }

        let num_cells: usize = exacts
            .iter()
            .map(|exact| exact_clonotypes[*exact].ncells())
            .sum();

        (num_cells, loupe_clonotype)
    });
    let mut results: Vec<_> = result_iter.collect();

    // Sort results in descending order by number of cells.

    results.sort_by_key(|(num_cells, _)| Reverse(*num_cells));

    let mut all_loupe_clonotypes = Vec::<Clonotype>::new();

    for (_, loupe_clonotype) in results {
        all_loupe_clonotypes.extend(loupe_clonotype);
    }

    // Write loupe output.
    loupe_out(ctl, all_loupe_clonotypes, refdata, dref)?;

    Ok(())
}

/// Sort exact subclonotypes.
fn sort_exact_clonotypes(
    ctl: &EncloneControl,
    refdata: &RefData,
    enclone_exacts: &EncloneExacts,
    od: &[Od],
    exacts: &mut Vec<usize>,
    mults: &mut Vec<usize>,
) {
    let EncloneExacts {
        to_bc,
        exact_clonotypes,
        raw_joins,
        info,
        candidate_clonotypes: _,
        drefs: _,
        sr,
    } = enclone_exacts;
    let mat = define_mat(
        to_bc,
        sr,
        ctl,
        exact_clonotypes,
        exacts,
        od,
        info,
        raw_joins,
        refdata,
    );
    let priority = exacts
        .iter()
        .enumerate()
        .map(|(u, &exact)| {
            let typex = mat.iter().map(|col| col[u].is_some()).collect::<Vec<_>>();
            let clonotype_id = exact;
            let ex = &exact_clonotypes[clonotype_id];
            let mut utot0 = 0;
            if let Some(mid) = mat[0][u] {
                let ex = &exact_clonotypes[clonotype_id];
                for j in 0..ex.clones.len() {
                    utot0 += ex.clones[j][mid].umi_count;
                }
            }
            (typex, ex.ncells(), utot0)
        })
        .collect::<Vec<_>>();
    let permutation = permutation::sort(&priority[..]);
    *exacts = permutation.apply_slice(&exacts[..]);
    *mults = permutation.apply_slice(&mults[..]);
    exacts.reverse();
    mults.reverse();
}

pub(crate) fn write_fate(path: &str, fate: &[BarcodeFates]) -> Result<(), String> {
    // Write out the fate of each filtered barcode.
    if !path.is_empty() {
        let mut wtr =
            BufWriter::new(File::create(path).expect("Unable to open FATE_FILE for writing"));
        serde_json::to_writer_pretty(&mut wtr, &fate).map_err(|e| e.to_string())?;
    }
    Ok(())
}
