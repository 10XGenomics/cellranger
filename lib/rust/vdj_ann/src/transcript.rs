// Copyright (c) 2018 10X Genomics, Inc. All rights reserved.
#![expect(missing_docs)]

// Code that analyzes transcripts.

use crate::annotate::{Annotation, get_cdr3_using_ann};
use crate::refx::RefData;
use amino::{have_start, have_stop};
use debruijn::dna_string::DnaString;
use debruijn::kmer::Kmer20;
use debruijn::{Mer, Vmer};
use hyperbase::Hyper;
use itertools::iproduct;
use kmer_lookup::make_kmer_lookup_20_single;
use serde::{Deserialize, Serialize};
use std::cmp::max;
use std::str::FromStr;
use vdj_types::{VDJ_CHAINS, VdjChain};
use vector_utils::{lower_bound1_3, unique_sort};

const MIN_DELTA: i32 = -25;
const MIN_DELTA_IGH: i32 = -55;
const MAX_DELTA: i32 = 35;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
// TEST FOR VALID VDJ SEQUENCE
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

#[derive(Debug, Serialize, Deserialize, Clone, PartialEq, Default)]
pub struct ContigStatus {
    pub full_length: Option<bool>,
    pub has_v_start: Option<bool>,
    pub in_frame: Option<bool>,
    pub no_premature_stop: Option<bool>,
    pub has_cdr3: Option<bool>,
    pub has_expected_size: Option<bool>,
    pub correct_ann_order: Option<bool>,
}

impl ContigStatus {
    fn is_productive(&self) -> bool {
        match (
            self.full_length,
            self.has_v_start,
            self.in_frame,
            self.no_premature_stop,
            self.has_cdr3,
            self.has_expected_size,
            self.correct_ann_order,
        ) {
            (
                Some(true),
                Some(true),
                Some(true),
                Some(true),
                Some(true),
                Some(true),
                Some(true),
            ) => true,
            (_, _, _, _, _, _, _) => false,
        }
    }

    fn order_by(&self) -> u16 {
        match (
            self.full_length,
            self.has_v_start,
            self.in_frame,
            self.no_premature_stop,
            self.has_cdr3,
            self.has_expected_size,
            self.correct_ann_order,
        ) {
            (
                Some(true),
                Some(true),
                Some(true),
                Some(true),
                Some(true),
                Some(true),
                Some(true),
            ) => 0,
            (Some(_), Some(_), Some(_), Some(_), Some(_), Some(_), Some(_)) => 1,
            (Some(_), Some(_), Some(_), Some(_), Some(_), Some(_), _) => 2,
            (Some(_), Some(_), Some(_), Some(_), Some(_), _, _) => 3,
            (Some(_), Some(_), Some(_), Some(_), _, _, _) => 4,
            (Some(_), Some(_), Some(_), _, _, _, _) => 5,
            (Some(_), Some(_), _, _, _, _, _) => 6,
            (Some(_), _, _, _, _, _, _) => 7,
            _ => 8,
        }
    }
}

#[derive(Clone, Copy)]
pub struct Vstart {
    ref_id: usize,
    tig_start: usize,
}

#[derive(Clone, Copy)]
pub struct Jstop {
    ref_id: usize,
    tig_stop: usize,
}

fn find_inframe_vdj_pair(vstarts: Vec<Vstart>, jstops: Vec<Jstop>) -> Option<(Vstart, Jstop)> {
    let mut vj_combinations: Vec<(Vstart, Jstop, i32)> = iproduct!(vstarts, jstops)
        .map(|(v, j)| (v, j, j.tig_stop as i32 - v.tig_start as i32))
        .filter(|(_, _, n)| n > &0)
        .filter(|(_, _, n)| n % 3 == 1)
        .collect();
    vj_combinations.sort_by_key(|x| x.2);
    vj_combinations.last().map(|(v, j, _)| (*v, *j))
}

fn evaluate_contig_status(
    vdj_chain: VdjChain,
    ann: &[Annotation],
    reference: &RefData,
    contig: &DnaString,
) -> Option<ContigStatus> {
    let valid_v_type = format!("{vdj_chain}V");
    let valid_j_type = format!("{vdj_chain}J");
    let rheaders = &reference.rheaders;
    let refs = &reference.refs;
    let mut vstarts: Vec<Vstart> = ann
        .iter()
        .filter(|a| reference.is_v(a.ref_id as usize))
        .filter(|a| rheaders[a.ref_id as usize].contains(&valid_v_type))
        .filter(|a| a.ref_start == 0)
        .map(|a| Vstart {
            ref_id: a.ref_id as usize,
            tig_start: a.tig_start as usize,
        })
        .collect();
    let jstops: Vec<Jstop> = ann
        .iter()
        .filter(|a| reference.is_j(a.ref_id as usize))
        .filter(|a| rheaders[a.ref_id as usize].contains(&valid_j_type))
        .filter(|a| a.ref_start + a.match_len == refs[a.ref_id as usize].len() as i32)
        .map(|a| Jstop {
            ref_id: a.ref_id as usize,
            tig_stop: a.tig_start as usize + a.match_len as usize,
        })
        .collect();

    if vstarts.is_empty() && jstops.is_empty() {
        return None;
    }

    let mut contig_status = ContigStatus {
        full_length: Some(!vstarts.is_empty() && !jstops.is_empty()),
        ..Default::default()
    };

    // filter vstarts to require START codon
    vstarts.retain(|v| have_start(contig, v.tig_start));
    contig_status.has_v_start = match (contig_status.full_length, vstarts.is_empty()) {
        (Some(true), false) => Some(true),
        (_, _) => Some(false),
    };

    let inframe_pair = find_inframe_vdj_pair(vstarts, jstops);
    contig_status.in_frame = Some(inframe_pair.is_some());

    if let Some((vstart, jstop)) = inframe_pair {
        contig_status.no_premature_stop = Some(
            !(vstart.tig_start..jstop.tig_stop - 3)
                .step_by(3)
                .any(|j| have_stop(contig, j)),
        );
    };

    let found_cdr3s = get_cdr3_using_ann(contig, reference, ann);
    contig_status.has_cdr3 = Some(!found_cdr3s.is_empty());

    if let (Some((vstart, jstop)), Some(cdr3)) = (inframe_pair, found_cdr3s.first()) {
        let expected_len = (refs[vstart.ref_id].len() + refs[jstop.ref_id].len()) as i32
            + (3 * cdr3.aa_seq.len() as i32)
            - 20;
        let observed_len = jstop.tig_stop as i32 - vstart.tig_start as i32;
        let delta = expected_len - observed_len;
        let min_delta = if vdj_chain == VdjChain::IGH {
            MIN_DELTA_IGH
        } else {
            MIN_DELTA
        };
        contig_status.has_expected_size = if delta < min_delta || delta > MAX_DELTA {
            Some(false)
        } else {
            Some(true)
        }
    };

    let observed_order: Vec<i32> = ann
        .iter()
        .map(|a| reference.segtype[a.ref_id as usize])
        .map(|s| match s {
            b'U' => 0,
            b'V' => 1,
            b'D' => 2,
            b'J' => 3,
            b'C' => 4,
            _ => panic!("Invalid segtype"),
        })
        .collect();
    let mut expected_order = observed_order.clone();
    expected_order.sort();
    contig_status.correct_ann_order = Some(observed_order == expected_order);

    Some(contig_status)
}

pub fn is_productive_contig(
    b: &DnaString,
    refdata: &RefData,
    ann: &[Annotation],
) -> (bool, ContigStatus) {
    let contig_status = VDJ_CHAINS
        .iter()
        .map(|chain| VdjChain::from_str(chain).unwrap())
        .filter_map(|chain| evaluate_contig_status(chain, ann, refdata, b))
        .min_by_key(ContigStatus::order_by);
    if let Some(cs) = contig_status {
        return (cs.is_productive(), cs.clone());
    }
    (false, ContigStatus::default())
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
// JUNCTION REGION CODE
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Given a valid contig, find the junction sequence, which we define to be the 100
// bases ending where the right end of a J region aligns to the contig.

pub fn junction_seq(tig: &DnaString, refdata: &RefData, ann: &[Annotation], jseq: &mut DnaString) {
    let refs = &refdata.refs;
    let segtype = &refdata.segtype;
    const TAG: i32 = 100;
    let jstop = ann
        .iter()
        .filter_map(|a| {
            if segtype[a.ref_id as usize] == b'J'
                && (a.ref_start + a.match_len) as usize == refs[a.ref_id as usize].len()
                && a.tig_start + a.match_len >= TAG
            {
                Some(a.tig_start + a.match_len)
            } else {
                None
            }
        })
        .min()
        .unwrap_or_else(|| panic!("Could not find jstop for a valid contig!"));
    let jstart = jstop - TAG;
    *jseq = tig.slice(jstart as usize, jstop as usize).to_owned();
}

// Given a valid contig, find the read support for the junction sequence.  Return a
// pair (UMIs, nreads), consisting of the number of UMIs that cover the junction
// sequence, and the total number of reads on those UMIs that cover the junction
// sequence.
//
// This is a very restrictive definition of junction support!
//
// Note that it is possible to find zero UMIs, because each UMI individually may
// not cover the junction.

#[derive(Default, Clone)]
pub struct JunctionSupportCore {
    pub reads: i32,
    pub umis: i32,
    pub umi_ids: Vec<i32>,
}

pub fn junction_supp_core(
    reads: &[DnaString],
    x: &Hyper,
    umi_id: &[i32],
    jseq: &DnaString,
) -> JunctionSupportCore {
    let mut ids = Vec::<i32>::new();
    // ◼ What we're doing here is converting a Vec<u32> into a Vec<i32>.
    // ◼ There should be a function to do that.
    for e in 0..x.h.g.edge_count() {
        for id in &x.ids[e] {
            ids.push(*id as i32);
        }
    }
    unique_sort(&mut ids);
    let tigs = vec![jseq.clone()];
    let mut kmers_plus = Vec::<(Kmer20, i32, i32)>::new();
    make_kmer_lookup_20_single(&tigs, &mut kmers_plus);
    let mut idi = 0;
    let k = x.h.k as usize;
    let mut jsupp = JunctionSupportCore::default();
    while idi < ids.len() {
        let mut idj = idi + 1;
        while idj < ids.len() && umi_id[ids[idj] as usize] == umi_id[ids[idi] as usize] {
            idj += 1;
        }
        let mut mm = Vec::<(i32, i32)>::new();
        for ida in &ids[idi..idj] {
            let b = &reads[*ida as usize];
            if b.len() < k {
                continue;
            }
            for j in 0..b.len() - k + 1 {
                let z: Kmer20 = b.get_kmer(j);
                let low = lower_bound1_3(&kmers_plus, &z) as usize;
                for kmer in &kmers_plus[low..] {
                    if kmer.0 != z {
                        break;
                    }
                    let p = kmer.2 as usize;
                    if j > 0 && p > 0 && b.get(j - 1) == jseq.get(p - 1) {
                        continue;
                    }
                    let mut len = k;
                    loop {
                        if j + len == b.len() || p + len == jseq.len() {
                            break;
                        }
                        if b.get(j + len) != jseq.get(p + len) {
                            break;
                        }
                        len += 1;
                    }
                    mm.push((p as i32, (p + len) as i32));
                }
            }
        }
        mm.sort_unstable();
        let mut cov = true;
        if mm.is_empty() || mm[0].0 > 0 {
            cov = false;
        }
        let mut reach = 0;
        for m in &mm {
            if m.0 <= reach {
                reach = max(reach, m.1);
            }
        }
        if reach < jseq.len() as i32 {
            cov = false;
        }
        if cov {
            jsupp.umis += 1;
            jsupp.reads += mm.len() as i32;
            jsupp.umi_ids.push(umi_id[ids[idi] as usize]);
        }
        idi = idj;
    }
    jsupp
}
