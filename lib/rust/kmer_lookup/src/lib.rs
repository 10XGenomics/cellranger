//! kmer_lookup
// Copyright (c) 2018 10X Genomics, Inc. All rights reserved.
#![deny(missing_docs)]

// Kmer lookup.

use debruijn::dna_string::DnaString;
use debruijn::kmer::{Kmer12, Kmer20};
use debruijn::{Kmer, Mer, Vmer};
use vector_utils::{lower_bound1_3, upper_bound1_3};

/// Given a vector of DnaStrings dv, create a sorted vector whose entries are
/// (kmer, e, estart), where the kmer starts at position estart on dv[e].
fn make_kmer_lookup_single<K: Kmer>(dv: &[DnaString], x: &mut Vec<(K, i32, i32)>) {
    let sz = dv
        .iter()
        .filter(|b| b.len() >= K::k())
        .map(|b| b.len() - K::k() + 1)
        .sum();
    x.clear();
    x.reserve(sz);

    for (i, b) in dv.iter().enumerate() {
        for (j, kmer) in b.iter_kmers().enumerate() {
            x.push((kmer, i as i32, j as i32));
        }
    }

    x.sort();
}

/// Included for backward compatibility. Use make_kmer_lookup_single
pub fn make_kmer_lookup_20_single(dv: &[DnaString], x: &mut Vec<(Kmer20, i32, i32)>) {
    make_kmer_lookup_single(dv, x);
}

/// Included for backward compatibility. Use make_kmer_lookup_single
pub fn make_kmer_lookup_12_single(dv: &[DnaString], x: &mut Vec<(Kmer12, i32, i32)>) {
    make_kmer_lookup_single(dv, x);
}

/// Determine if a sequence perfectly matches in forward orientation.
pub fn match_12(b: &DnaString, dv: &[DnaString], x: &[(Kmer12, i32, i32)]) -> bool {
    let y: Kmer12 = b.get_kmer(0);
    let low = lower_bound1_3(x, &y);
    let high = upper_bound1_3(x, &y);
    for m in low..high {
        let mut l = 12;
        let t = x[m as usize].1 as usize;
        let mut p = x[m as usize].2 as usize + 12;
        while l < b.len() && p < dv[t].len() {
            if b.get(l) != dv[t].get(p) {
                break;
            }
            l += 1;
            p += 1;
        }
        if l == b.len() {
            return true;
        }
    }
    false
}
