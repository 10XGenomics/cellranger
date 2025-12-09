// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
#![deny(missing_docs)]

// Miscellaneous functions.

use super::BarcodeFilter;
use crate::Dataset;
use crate::core::barcode_fate::BarcodeFate;
use crate::core::defs::{BarcodeContigs, ExactClonotype};
use itertools::Itertools;
use vector_utils::{bin_member, bin_position, erase_if, next_diff, next_diff1_3, unique_sort};

/// If a V..J segment appears in exactly one dataset, with frequency n, let x be the total
/// number of productive pairs for that dataset, and let y be the total number of productive
/// pairs for all datasets from the same origin.  If (x/y)^n <= 10^-6, i.e. the probability
/// that assuming even distribution, all instances of that V..J ended up in that one dataset,
/// delete all the productive pairs for that V..J segment that do not have at least 100
/// supporting UMIs.  (Note no attempt to do Bonferroni correction.)
///
/// For the case of two datasets for one origin, with equal numbers of productive pairs in
/// each, this corresponds roughly to the case n = 20.
///
/// Note that we could modify this to allow *some* occurrences in other datasets.
///
/// There are only certain ways that these misdistribution events could happen:
///
/// 1. A cell (and particularly a plasma cell or plasmablast) bursts after drawing cells to
///    make libraries, leaving behind cell fragments that seed separate GEMs
///    (probably most likely).
/// 2. Multiple gel beads end up in one GEM.
/// 3. Something involving like cells sticking together and subsequently separating.
/// 4. Physical contamination of libraries.
/// 5. Informatic mixup of libraries.
/// 6. Nothing other than a low probability event (unlikely).
///
/// Note that in case 1, we have evidence that a plasma cell or plasmablast existed in the
/// original cells that were drawn (perhaps breaking up in the process of drawing), and was
/// subsequently distintegrated.
pub(crate) struct CrossFilter<'a> {
    pub(crate) origin_info: &'a [Dataset],
}

impl BarcodeFilter<BarcodeContigs> for CrossFilter<'_> {
    fn fate_keys(&self, item: &BarcodeContigs) -> impl Iterator<Item = (usize, String)> {
        std::iter::once((item[0].dataset_index, item[0].barcode.clone()))
    }

    fn filter(&self, tig_bc: &[BarcodeContigs]) -> Vec<Option<BarcodeFate>> {
        // Get the list of dataset origins.  Here we allow the same origin name to have been used
        // for more than one donor, as we haven't explicitly prohibited that.

        let origins: Vec<_> = self
            .origin_info
            .iter()
            .map(|ds| (&ds.donor_id, &ds.origin_id))
            .sorted()
            .dedup()
            .collect();
        let to_origin = self
            .origin_info
            .iter()
            .map(|ds| bin_position(&origins, &(&ds.donor_id, &ds.origin_id)) as usize)
            .collect::<Vec<_>>();

        // For each dataset index, and each origin, compute the total number of productive pairs.

        let mut n_dataset_index = vec![0; self.origin_info.len()];
        let mut n_origin = vec![0; origins.len()];
        for tigi in tig_bc {
            for x in tigi {
                n_dataset_index[x.dataset_index] += 1;
                n_origin[to_origin[x.dataset_index]] += 1;
            }
        }

        // Find all the V..J segments, and for each, the number of times it appears in each
        // dataset ID.
        //
        // Note that there is no point running this unless we have at least two dataset IDs, and in
        // fact unless there is an origin with at least two dataset IDs.  Better: just gather data
        // for the origin for which there are at least two dataset IDs.  Also no point if NCROSS.

        let vjx = {
            let mut vjx = Vec::<(&[u8], usize, usize)>::new(); // (V..J, dataset index, count)
            for tigi in tig_bc {
                for x in tigi {
                    vjx.push((x.seq(), x.dataset_index, 1));
                }
            }
            vjx.sort();
            let mut to_delete = vec![false; vjx.len()];
            let mut i = 0;
            while i < vjx.len() {
                let j = next_diff(&vjx, i); // actually only need to check first two fields
                vjx[i].2 = j - i;
                for d in &mut to_delete[i + 1..j] {
                    *d = true;
                }
                i = j;
            }
            erase_if(&mut vjx, &to_delete);
            vjx
        };

        // Now do the cross filter.

        let mut blacklist = Vec::<&[u8]>::new();
        let mut i = 0;
        while i < vjx.len() {
            let j = next_diff1_3(&vjx, i);
            if j - i == 1 {
                let dataset_index = vjx[i].1;
                let n = vjx[i].2;
                let x = n_dataset_index[dataset_index];
                let y = n_origin[to_origin[dataset_index]];
                if y > 0 {
                    let p = (x as f64 / y as f64).powi(n as i32);
                    if p <= 1.0e-6 {
                        blacklist.push(vjx[i].0);
                    }
                }
            }
            i = j;
        }
        blacklist.sort();
        let mut to_delete = vec![None; tig_bc.len()];
        const UMIS_SAVE: usize = 100;
        for (i, tigi) in tig_bc.iter().enumerate() {
            for tig in tigi {
                if tig.umi_count < UMIS_SAVE && bin_member(&blacklist, &tig.seq()) {
                    to_delete[i] = Some(BarcodeFate::Cross);
                    break;
                }
            }
        }
        to_delete
    }
}

/// Filter out some foursie artifacts.
pub(crate) struct ArtificialFoursieFilter;

impl BarcodeFilter<ExactClonotype> for ArtificialFoursieFilter {
    fn fate_keys(&self, item: &ExactClonotype) -> impl Iterator<Item = (usize, String)> {
        item.clones
            .iter()
            .map(|clone| (clone[0].dataset_index, clone[0].barcode.clone()))
    }

    fn filter(&self, exact_clonotypes: &[ExactClonotype]) -> Vec<Option<BarcodeFate>> {
        let mut twosies = Vec::<(&[u8], &[u8])>::new();
        for ex in exact_clonotypes {
            if ex.share.len() == 2 && (ex.share[0].left ^ ex.share[1].left) && ex.ncells() >= 10 {
                twosies.push((ex.share[0].seq.as_ref(), ex.share[1].seq.as_ref()));
            }
        }
        unique_sort(&mut twosies);
        exact_clonotypes
            .iter()
            .map(|ex| {
                if ex.share.len() == 4 {
                    for (i1, s1) in ex.share.iter().enumerate() {
                        for s2 in &ex.share[i1 + 1..4] {
                            if s1.left ^ s2.left {
                                let p = (s1.seq.as_ref(), s2.seq.as_ref());
                                if bin_member(&twosies, &p) {
                                    return Some(BarcodeFate::FoursieKill);
                                }
                            }
                        }
                    }
                }
                None
            })
            .collect()
    }
}
