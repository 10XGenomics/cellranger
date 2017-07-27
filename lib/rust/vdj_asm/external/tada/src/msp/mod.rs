//
// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
//

//! Compute the MSP of a string
//!
//! simple_scan method is based on:
//! Li, Yang. "MSPKmerCounter: a fast and memory efficient approach for k-mer counting." arXiv preprint arXiv:1505.06550 (2015).

use std::collections::HashSet;
use std::iter::FromIterator;
use std::cmp::min;
use rand::Rng;
use std::iter::Iterator;

#[inline(never)]
fn compute_pvals(p: usize, seq: &[u8]) -> Vec<u32> {
    if p > 16 {
        panic!("p can't be greater than 16");
    }

    let mut c: u32 = 0;
    let mut pvals: Vec<u32> = Vec::with_capacity(seq.len());

    if seq.len() < p {
        return pvals;
    }

    for i in 0..p {
        c = (c << 2) | ((seq[i as usize] as u32) & 0x3)
    }

    // v is the array of p-mers starting a each position of seq

    pvals.push(c);

    // Mask for the largest possible p-mer
    let mask: u32 = (1u32 << (2 * p)) - 1;

    for i in p..seq.len() {
        c = ((c << 2) | ((seq[i] as u32) & 0x3)) & mask;
        pvals.push(c);
    }

    pvals
}

#[inline(never)]
pub fn rc(seq: &[u8]) -> Vec<u8> {
    let rc: Vec<u8> = seq.iter().rev().map(|x| 3 - x).collect();
    rc
}

/// Determine MSP substrings of seq, for given k and p.
/// Returns a vector of tuples indicating the substrings, and the pmer values:
/// (p-mer value, min p-mer position, start position, end position)
/// permutation is a permutation of the lexicographically-sorted set of all pmers.
/// A permutation of pmers sorted by their inverse frequency in the dataset will give the
/// most even bucketing of MSPs over pmers.
pub fn simple_scan(k: usize,
                   p: usize,
                   seq: &[u8],
                   permutation: &Vec<usize>,
                   allow_rc: bool)
                   -> Vec<(u32, usize, usize, usize)> {
    let pfwd = compute_pvals(p, seq);
    let mut pvals: Vec<u32> = Vec::new();

    if allow_rc {
        let prev: Vec<u32> = compute_pvals(p, &rc(seq)[..]).into_iter().rev().collect();

        for (x, y) in pfwd.iter().zip(prev.iter()) {
            let mv = min(permutation[*x as usize] as u32,
                         permutation[*y as usize] as u32);
            pvals.push(mv);
        }
    } else {
        for x in pfwd.iter() {
            pvals.push(permutation[*x as usize] as u32);
        }
    }

    let pmin = |i: usize, j: usize| if pvals[i] <= pvals[j] {
        i
    } else {
        j
    };

    let m = seq.len();

    let find_min = |start, stop| {
        let mut pos = start;
        let mut min_pos = start;
        let mut min_val = pvals[start];

        while pos < stop {
            pos = pos + 1;
            if pvals[pos] < min_val {
                min_pos = pos;
                min_val = pvals[pos];
            }
        }

        min_pos
    };

    let mut min_positions = Vec::with_capacity(16);
    let mut min_pos = find_min(0, k - p);
    min_positions.push((0, min_pos));

    for i in 0..(m - k + 1) {
        if i > min_pos {
            min_pos = find_min(i, i + k - p);
            min_positions.push((i, min_pos));
        } else {
            let test_min = pmin(min_pos, i + k - p);

            if test_min != min_pos {
                min_pos = test_min;
                min_positions.push((i, min_pos));
            }
        }
    }

    let mut slices = Vec::with_capacity(min_positions.len());

    // Generate the slices of the final string
    for p in 0..min_positions.len() - 1 {
        let (start_pos, min_pos) = min_positions[p];
        let (next_pos, _) = min_positions[p + 1];
        slices.push((pvals[min_pos], min_pos, start_pos, next_pos + k - 1 - start_pos));
    }

    let (last_pos, min_pos) = min_positions[min_positions.len() - 1];
    slices.push((pvals[min_pos], min_pos, last_pos, m - last_pos));

    slices
}


fn all_kmers<T>(k: usize, seq: &[T]) -> Vec<&[T]> {
    (0..(seq.len() - k + 1)).map(|i| &seq[i..i + k]).collect()
}

fn test_all_kmers(k: usize, full_seq: &[u8], slices: Vec<(u32, usize, usize, usize)>) {
    let start_kmers = HashSet::from_iter(all_kmers(k, full_seq));

    let mut sliced_kmers = HashSet::new();

    for (_, _, slc_start, slc_len) in slices {
        let slc = &full_seq[slc_start..(slc_start + slc_len)];
        sliced_kmers.extend(all_kmers(k, slc));
    }

    if start_kmers != sliced_kmers {
        println!("start kmers: {:?}", start_kmers);
        println!("sliced kmers: {:?}", sliced_kmers);
        panic!("kmer sets not equal")
    }
}

fn random_dna<R: Rng>(sz: usize, r: &mut R) -> Vec<u8> {
    let mut vec = Vec::new();

    for _ in 0..sz {
        let v = (r.next_u32() % 4) as u8;
        vec.push(v);
    }
    vec
}


#[cfg(test)]
mod tests {
    use super::all_kmers;
    use rand;

    #[test]
    fn test1() {
        let v = vec![1u8, 2u8, 3u8, 4u8, 5u8, 6u8];
        let s = &v[..];
        let ak = all_kmers(2, s);

        if ak[0] < ak[1] {
            println!("sorts!")
        }

        println!("{:?}", ak);


        let v = vec![6u8, 5u8, 4u8, 3u8, 2u8, 1u8, 0u8];
        let s = &v[..];
        let mut ak = all_kmers(2, s);

        println!("{:?}", ak);
        ak.sort();


        if ak[0] < ak[1] {
            println!("sorts!")
        }
        println!("{:?}", ak);

    }

    #[test]
    fn test_slice() {
        let p = 8;
        let permutation: Vec<usize> = (0..(1 << 2 * p)).collect();

        let mut rng = rand::XorShiftRng::new_unseeded();
        for _ in 0..10 {
            let k = 50usize;
            let dna = super::random_dna(150, &mut rng);
            println!("{:?}", dna);
            let slices = super::simple_scan(k, p, &dna[..], &permutation, true);
            println!("Made {} slices from dna of length {:?}",
                     slices.len(),
                     dna.len());

            println!("slices: {:?}", slices);
            super::test_all_kmers(k, &dna[..], slices);
        }
    }

    #[test]
    fn test_sample() {
        // for testing MSP on specific error cases

        // let v1 : Vec<u8> = vec![3, 0, 3, 0, 3, 2, 3, 1, 0, 0, 2, 0, 0, 1, 1, 3, 0, 2, 0, 3, 2, 3, 3, 1, 2, 0, 2, 2, 3, 1, 0, 1, 2, 1, 2, 3, 2, 1, 0, 2, 1, 0, 2, 2, 1, 0, 2, 3, 0, 3, 2, 3, 0, 2, 0, 0, 1, 0, 2, 3, 1, 3, 2, 0, 2, 2, 2, 2, 1, 2, 0, 2, 1, 0, 0, 1, 3, 0, 0, 0, 2, 2, 3, 0, 0, 3, 2, 3, 1, 2, 0, 2, 0, 3, 2, 1, 3, 2, 2, 3, 3, 1, 2, 3, 0, 3, 1, 1, 2, 2, 2, 3, 1, 1, 2, 0, 2, 3, 3, 2, 0, 3, 1, 1, 1, 1, 3, 3, 3, 0, 1, 0, 0, 3, 0, 2, 3, 2, 2, 2, 3, 2, 3, 1, 0, 2, 3, 0, 2, 2, 2, 1, 0, 2, 3, 3, 3, 2, 2, 1, 3, 1, 0, 2, 0, 1, 0, 1, 2, 2, 2, 2, 2, 3, 0, 3, 3, 1, 2, 3, 2, 0, 2, 1, 0, 3, 1, 1, 0, 2, 1, 3, 0, 1, 2, 1, 1, 1, 3, 1, 0, 3, 0, 0, 2, 2, 0, 3, 2, 3, 2, 0, 2, 1, 2, 3, 0, 1, 3, 1, 1, 2, 2, 3, 0, 3, 3, 3, 3, 0, 3, 3, 0, 0, 0, 2, 3, 3, 1, 3, 3, 1, 2, 2, 0, 3, 2, 1, 0, 2];
        // let v2 : Vec<u8> = vec![3, 0, 3, 0, 3, 2, 3, 1, 0, 0, 2, 0, 0, 1, 1, 3, 0, 2, 0, 3, 2, 3, 3, 1, 2, 0, 2, 2, 3, 2, 0, 3, 2, 1, 2, 3, 2, 1, 0, 2, 1, 0, 2, 2, 1, 0, 2, 3, 0, 3, 2, 3, 0, 2, 3, 0, 1, 0, 2, 3, 1, 3, 2, 0, 2, 2, 0, 2, 1, 2, 0, 2, 1, 0, 0, 1, 3, 0, 0, 0, 2, 2, 3, 0, 0, 3, 2, 3, 1, 2, 0, 2, 0, 3, 2, 1, 3, 2, 2, 3, 3, 1, 2, 3, 0, 3, 1, 1, 2, 2, 2, 3, 1, 1, 2, 0, 0, 3, 3, 2, 0, 3, 1, 1, 1, 1, 3, 3, 3, 0, 1, 0, 0, 3, 0, 2, 3, 1, 2, 1, 3, 2, 3, 1, 0, 2, 3, 0, 2, 2, 2, 2, 1, 2, 3, 3, 3, 2, 2, 1, 3, 1, 0, 2, 0, 1, 0, 1, 2, 2, 2, 2, 2, 3, 0, 3, 3, 1, 2, 3, 2, 0, 2, 1, 0, 3, 1, 1, 0, 2, 1, 3, 0, 1, 2, 1, 2, 1, 3, 1, 0, 3, 0, 0, 2, 2, 0, 3, 2, 3, 2, 0, 2, 1, 2, 3, 0, 1, 3, 1, 1, 2, 2, 3, 0, 3, 3, 3, 3, 0, 3, 3, 0, 0, 0, 2, 3, 3, 1, 0, 3, 1, 2, 2, 0, 3, 2, 1, 0, 2];

        let v1: Vec<u8> = vec![3, 0, 3, 0, 1, 2, 3, 3, 0, 0, 0, 1, 1, 0, 3, 3, 1, 2, 0, 1, 1, 3,
                               2, 1, 1, 1, 2, 3, 3, 2, 1, 2, 2, 1, 2, 3, 2, 1, 0, 2, 1, 1, 2, 1,
                               0, 1, 2, 3, 0, 3, 2, 3, 0, 1, 3, 0, 1, 0, 2, 3, 3, 3, 2, 3, 2, 2,
                               0, 2, 1, 2, 2, 2, 0, 0, 1, 0, 1, 2, 1, 0, 3, 2, 3, 1, 0, 3, 2, 2,
                               2, 2, 1, 3, 0, 2, 1, 2, 1, 3, 3, 0, 1, 1, 2, 3, 0, 2, 3, 1, 3, 2,
                               3, 1, 1, 0, 2, 2, 1, 1, 1, 2, 0, 0, 1, 2, 1, 0, 3, 3, 3, 0, 1, 1,
                               0, 3, 0, 2, 2, 1, 2, 1, 3, 2, 3, 1, 0, 2, 3, 0, 2, 2, 0, 3, 3, 2,
                               3, 3, 0, 3, 0, 1, 1, 3, 0, 1, 0, 1, 0, 1, 2, 2, 0, 3, 3, 3, 1, 3,
                               1, 1, 0, 1, 2, 0, 2, 1, 0, 3, 1, 1, 1, 2, 0, 3, 0, 1, 2, 1, 1, 1,
                               3, 1, 2, 3, 0, 3, 2, 2, 0, 3, 2, 3, 3, 0, 2, 3, 2, 0, 0, 1, 2, 0,
                               3, 2, 2, 1, 2, 2, 3, 3, 3, 2, 3, 0, 3, 1, 0, 2, 3, 3, 0, 3, 1, 0,
                               3, 2, 1, 3, 2, 1, 1, 2];
        let v2: Vec<u8> = vec![3, 0, 3, 0, 1, 3, 3, 3, 0, 0, 0, 1, 1, 0, 3, 3, 1, 2, 0, 1, 1, 3,
                               2, 1, 1, 1, 2, 1, 3, 2, 1, 2, 2, 1, 2, 3, 2, 1, 0, 2, 1, 1, 2, 1,
                               1, 1, 2, 3, 0, 3, 2, 3, 0, 1, 3, 0, 1, 0, 2, 3, 3, 3, 2, 3, 2, 2,
                               0, 2, 1, 2, 2, 2, 0, 0, 1, 0, 1, 2, 1, 0, 3, 2, 3, 1, 0, 3, 2, 2,
                               2, 2, 1, 3, 0, 2, 1, 2, 1, 3, 3, 0, 1, 1, 2, 3, 0, 2, 3, 1, 3, 2,
                               3, 1, 1, 0, 2, 2, 1, 1, 0, 2, 0, 0, 1, 2, 1, 0, 3, 3, 3, 0, 1, 1,
                               0, 3, 0, 2, 2, 1, 2, 1, 3, 2, 3, 1, 0, 2, 3, 0, 2, 2, 0, 3, 3, 2,
                               3, 3, 0, 3, 0, 1, 1, 3, 0, 1, 0, 1, 0, 1, 2, 2, 0, 3, 3, 3, 1, 3,
                               1, 1, 0, 1, 2, 0, 2, 1, 0, 3, 1, 1, 1, 2, 0, 3, 0, 1, 2, 1, 3, 1,
                               3, 1, 2, 3, 0, 3, 2, 2, 0, 3, 2, 3, 3, 0, 2, 3, 2, 0, 0, 1, 2, 0,
                               3, 2, 2, 1, 2, 2, 3, 3, 3, 2, 3, 1, 3, 1, 0, 2, 2, 3, 0, 3, 1, 0,
                               3, 3, 1, 3, 2, 1, 1, 2];

        let p = 5;
        let permutation: Vec<usize> = (0..(1 << 2 * p)).collect();

        let s1 = super::simple_scan(35, 5, &v1[..], &permutation, true);
        let s2 = super::simple_scan(35, 5, &v2[..], &permutation, true);

        println!("{:?}", s1);
        println!("{:?}", s2);
        // assert!(false);
    }
}
