//
// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
//

/// Implementation of the FM-Index and BWT for BitEnc objects with 2bit encodings.

use std::iter::repeat;
use kmer;
use std::cmp::Ordering;
use bitenc;

use bincode;
use bincode::rustc_serialize::{encode_into, decode_from};

use std::fs::File;
use std::io::{BufWriter, BufReader};

use flate2::write::ZlibEncoder;
use flate2::read::ZlibDecoder;
use flate2::Compression;


// Inplace implementation of prescan over a slice.
pub fn prescan<T: Copy, F: Fn(T, T) -> T>(a: &mut [T], neutral: T, op: F) {
    let mut s = neutral;
    for v in a.iter_mut() {
        let t = *v;
        *v = s;
        s = op(s, t);
    }
}

/// An occurence array implementation.
#[derive(RustcEncodable, RustcDecodable)]
pub struct Occ {
    occ: Vec<Vec<usize>>,
    k: usize,
}


impl Occ {
    pub fn new(bwt: &bitenc::BitEnc, k: usize) -> Self {
        let n = bwt.len();
        let mut occ = Vec::with_capacity(n / k);
        let mut curr_occ: Vec<usize> = repeat(0).take(4).collect();
        for (i, c) in bwt.iter().enumerate() {
            curr_occ[c as usize] += 1;
            if i % k == 0 {
                occ.push(curr_occ.clone());
            }
        }

        Occ { occ: occ, k: k }
    }

    /// Get occurrence count of symbol a in BWT[..r+1].
    /// Complexity: O(k).
    pub fn get(&self, bwt: &bitenc::BitEnc, r: usize, a: u32) -> usize {
        let i = r / self.k;
        let mut count = self.occ[i][a as usize];
        for pos in i * self.k + 1..r + 1 {
            count += (bwt.get(pos).unwrap() as u32 == a) as usize;
        }
        count
    }
}


/// Calculate the less array for a given BWT. Complexity O(n).
pub fn less(bwt: &bitenc::BitEnc) -> Vec<usize> {
    let mut less: Vec<usize> = repeat(0).take(4).collect();
    for c in bwt.iter() {
        less[c as usize] += 1;
    }
    // calculate +-prescan
    prescan(&mut less[..], 0, |a, b| a + b);

    less
}

#[derive(RustcEncodable, RustcDecodable)]
pub struct FMIndex {
    bwt: bitenc::BitEnc,
    sa: Vec<usize>,
    sa_step: usize,
    less: Vec<usize>,
    occ: Occ,
}


impl FMIndex {
    pub fn empty() -> Self {
        FMIndex {
            bwt: bitenc::BitEnc::new(),
            sa: Vec::new(),
            sa_step: 1,
            less: Vec::new(),
            occ: Occ::new(&bitenc::BitEnc::new(), 1),
        }
    }

    pub fn new(bwt: bitenc::BitEnc, sa: Vec<usize>, sa_step: usize, k: usize) -> Self {
        // The next two lines will borrow the bwt, so they can't happen at
        // the same time as the creation of the FMIndex object (which will copy it).
        let less = less(&bwt);
        let occ = Occ::new(&bwt, k);
        FMIndex {
            bwt: bwt,
            sa: sa,
            sa_step: sa_step,
            less: less,
            occ: occ,
        }
    }

    pub fn len(&self) -> usize {
        self.bwt.len()
    }

    pub fn backward_search<'b, P: Iterator<Item = &'b u8> + DoubleEndedIterator>(&self,
                                                                                 pattern: P)
                                                                                 -> Vec<usize> {
        let (mut l, mut r) = (0, self.bwt.len() - 1);
        for &c in pattern.rev() {
            let a = kmer::base_to_bits(c);
            let less = self.less(a);
            l = less +
                if l > 0 {
                self.occ(l - 1, a)
            } else {
                0
            };
            if less + self.occ(r, a) < 1 {
                return Vec::new();
            }
            r = less + self.occ(r, a) - 1;
        }
        (l..r + 1).map(|x| self.get_sa_pos(x)).collect()
    }

    pub fn backward_search_bitenc(&self, bit_seq: bitenc::BitEnc) -> Vec<usize> {
        let (mut l, mut r) = (0, self.bwt.len() - 1);
        let vals: Vec<u8> = bit_seq.iter().collect();
        for &a in vals.iter().rev() {
            let less = self.less(a);
            l = less +
                if l > 0 {
                self.occ(l - 1, a)
            } else {
                0
            };
            if less + self.occ(r, a) < 1 {
                return Vec::new();
            }
            r = less + self.occ(r, a) - 1;
        }
        (l..r + 1).map(|x| self.get_sa_pos(x)).collect()
    }

    fn get_sa_pos(&self, i: usize) -> usize {
        let mut steps = 0;
        let mut pos = i;
        while pos > 0 && pos % self.sa_step != 0 && pos < self.bwt.len() {
            let a = self.bwt.get(pos).unwrap();
            pos = self.less(a) + self.occ(pos - 1, a);
            steps += 1;
        }
        (self.sa[pos / self.sa_step] + steps) % self.bwt.len()
    }

    fn occ(&self, r: usize, a: u8) -> usize {
        self.occ.get(&self.bwt, r, a as u32)
    }

    fn less(&self, a: u8) -> usize {
        self.less[a as usize]
    }
}

fn sort_bucket(bucket: &mut Vec<usize>,
               bitenc: &bitenc::BitEnc,
               bwt: &mut bitenc::BitEnc,
               sa: &mut Vec<usize>,
               sa_step: usize) {

    let len = bitenc.len();

    let suffix_cmp = |a: &usize, b: &usize| {
        let mut offset = 0;
        // will perform at most len comparisons
        for _ in 0..len {
            // Towards the end of the array we'll have to loop over.
            let cmp = bitenc.get((*a + offset) % len).cmp(&(bitenc.get((*b + offset) % len)));
            match cmp {
                Ordering::Greater => return cmp,
                Ordering::Less => return cmp,
                Ordering::Equal => {
                    offset += 1;
                }
            }
        }
        Ordering::Equal
    };

    bucket.sort_by(suffix_cmp);

    for bucket_pos in 0..bucket.len() {
        let sa_value = bucket[bucket_pos];
        let v = if sa_value > 0 {
            bitenc.get(sa_value - 1).expect("missing")
        } else {
            bitenc.get(bitenc.len() - 1).expect("missing")
        };
        if bwt.len() % sa_step == 0 {
            sa.push(bucket[bucket_pos]);
        }
        bwt.push(v);
    }
}

#[derive(Debug, PartialEq, Eq, RustcEncodable, RustcDecodable)]
pub struct BWT {
    pub bwt: bitenc::BitEnc,
    pub sa: Vec<usize>,
}

/// Computes the BWT of the given BitEnc.
/// npos : number of positions to use for bucketing
/// sa_step: will only store suffix array every sa_step positions
pub fn compute_bwt(bitenc: &bitenc::BitEnc,
                   npos: usize,
                   sa_step: usize)
                   -> BWT {
    let mut val = 0;
    let mask = (1 << (npos * bitenc.width())) - 1;
    let len = bitenc.len();

    let mut slots: Vec<Vec<usize>> = Vec::new();
    for _ in 0..(mask + 1) {
        slots.push(Vec::new());
    }
    // Bucket positions based on the starting n-mer
    for i in 0..len {
        if i == 0 {
            // Get the first npos values
            for i in 0..npos {
                val = val << bitenc.width();
                val = val | (bitenc.get(i).unwrap_or_default() as usize);
            }
        } else {
            val = (val << bitenc.width()) & mask;
            // Towards the end of the string we need to loop over.
            let next_pos = (i + npos - 1) % len;
            let next_val = bitenc.get(next_pos).unwrap_or_default() as usize;
            val = val | next_val;
        }
        slots[val as usize].push(i);
    }

    let mut bwt = bitenc::BitEnc::new();
    let mut sa: Vec<usize> = Vec::new();

    for bucket in slots.iter_mut() {
        sort_bucket(bucket, bitenc, &mut bwt, &mut sa, sa_step);
    }
    BWT { bwt: bwt, sa: sa }
}

pub fn compute_bwt_bucket(bitenc: &bitenc::BitEnc, npos: usize, bucket_val: usize) -> BWT {
    let mut val = 0;
    let mask = (1 << (npos * bitenc.width())) - 1;
    let len = bitenc.len();

    let mut bucket: Vec<usize> = Vec::new();

    // Bucket positions based on the starting n-mer
    for i in 0..len {
        if i == 0 {
            // Get the first npos values
            for i in 0..npos {
                val = val << bitenc.width();
                val = val | (bitenc.get(i).unwrap_or_default() as usize);
            }
        } else {
            val = (val << bitenc.width()) & mask;
            // Towards the end of the string we need to loop over.
            let next_pos = (i + npos - 1) % len;
            let next_val = bitenc.get(next_pos).unwrap_or_default() as usize;
            val = val | next_val;
        }
        if val == bucket_val {
            bucket.push(i);
        }
    }

    let mut bwt = bitenc::BitEnc::new();
    let mut sa: Vec<usize> = Vec::new();

    sort_bucket(&mut bucket, bitenc, &mut bwt, &mut sa, 1);
    BWT { bwt: bwt, sa: sa }
}

pub fn write_bwt(bwt: &BWT, filename: &str) {
    let f = File::create(filename).unwrap();
    let writer = BufWriter::new(f);
    let mut encoder = ZlibEncoder::new(writer, Compression::Best);
    encode_into(&bwt, &mut encoder, bincode::SizeLimit::Infinite).unwrap();
}

pub fn read_bwt(filename: &str) -> BWT {
    let f = File::open(filename).unwrap();
    let reader = BufReader::new(f);
    let mut decoder = ZlibDecoder::new(reader);
    let b: BWT = decode_from(&mut decoder, bincode::SizeLimit::Infinite).unwrap();
    b
}

pub fn merge_bwt_buckets(filenames: &Vec<String>, sa_step: usize) -> BWT {
    let mut bwt = bitenc::BitEnc::new();
    let mut sa: Vec<usize> = Vec::new();
    let mut count = 0;

    for filename in filenames {
        let f = File::open(filename).unwrap();
        let reader = BufReader::new(f);
        let mut decoder = ZlibDecoder::new(reader);
        let b: BWT = decode_from(&mut decoder, bincode::SizeLimit::Infinite).unwrap();

        for val in b.bwt.iter() {
            bwt.push(val);
        }

        for val in b.sa {
            if count % sa_step == 0 {
                sa.push(val);
            }
            count += 1;
        }
    }
    BWT { bwt: bwt, sa: sa }
}



#[cfg(test)]
mod tests {
    use super::*;
    use bitenc;
    use std::fs;
    use sim_tests;

    // Check for consistency between the BWT and SA when you have a full SA
    fn bwt_consistency(bwt: &BWT, dna: &Vec<u8>)
    {
        if bwt.bwt.len() == bwt.sa.len() {
            for i in 0..bwt.bwt.len() {
                let pos = (bwt.sa[i] + bwt.sa.len() - 1) % bwt.sa.len();
                assert_eq!(bwt.bwt.get(i).unwrap(), dna[pos])
            }
        }
    }

    fn check_bwt(dna: &Vec<u8>, npos: usize) {
        let bitenc = bitenc::BitEnc::from_bytes(dna);

        // Make sure that the BWT is identical irrespective of the number of buckets.
        let bwt1 = compute_bwt(&bitenc, 1, 1);
        let bwt2 = compute_bwt(&bitenc, npos, 1);
        assert_eq!(bwt1, bwt2);
        bwt_consistency(&bwt1, &dna);
    }

    #[test]
    fn test_bwt()
    {
        for npos in 2..9 {
            for iter in 0..2 {
                let dna = sim_tests::random_dna(5000);
                check_bwt(&dna, npos);
            }
        }
    }


    #[test]
    fn test_fm_index() {
        let dna = "ACGTACGT";
        // BWT is TTAACCGG
        let bitenc = bitenc::BitEnc::from_dna_string(dna);

        // Make sure that the BWT is identical irrespective of the number of buckets.
        let bwt1 = compute_bwt(&bitenc, 1, 1);
        let bwt2 = compute_bwt(&bitenc, 2, 1);
        assert_eq!(bwt1, bwt2);
        assert!(bwt1.sa == [4, 0, 5, 1, 6, 2, 7, 3] || bwt1.sa == [0, 4, 1, 5, 2, 6, 3, 7]);
        assert_eq!(&bwt1.bwt.to_dna_string(), "TTAACCGG");

        let fm = FMIndex::new(bwt1.bwt, bwt1.sa, 1, 1);
        // fm.less[i] is the number of characters that are smaller than i
        assert_eq!(fm.less, [0, 2, 4, 6]);
        assert_eq!(fm.occ.occ.len(), dna.len());
        // Number of occurrences of each letter up to position 3 of the BWT
        assert_eq!(fm.occ.occ[3], [2, 0, 0, 2]);
        let pattern = b"A";
        let pos = fm.backward_search(pattern.iter());
        assert!(pos == [0, 4] || pos == [4, 0]);

        let pattern = b"CGT";
        let pos = fm.backward_search(pattern.iter());
        assert!(pos == [1, 5] || pos == [5, 1]);

        // No sentinel used so occurrences should loop over the end of string.
        let pattern = b"GTA";
        let pos = fm.backward_search(pattern.iter());
        assert!(pos == [2, 6] || pos == [6, 2]);

        // Test different step sizes for the occ matrix
        // Make sure that the BWT is identical irrespective of the number of buckets.
        let bwt1 = compute_bwt(&bitenc, 1, 1);
        let fm2 = FMIndex::new(bwt1.bwt, bwt1.sa, 1, 3);
        let pattern = b"CGT";
        let pos = fm2.backward_search(pattern.iter());
        assert!(pos == [1, 5] || pos == [5, 1]);

        for sa_val in 1..10 {
            // Test runs of the same character
            let dna = "TGCATTAGAAAACTGCA";
            let bitenc = bitenc::BitEnc::from_dna_string(dna);
            let bwt1 = compute_bwt(&bitenc, 3, sa_val);
            let fm = FMIndex::new(bwt1.bwt, bwt1.sa, sa_val, 3);
            let pos = fm.backward_search(b"AAA".iter());
            assert!(pos == [8, 9] || pos == [9, 8]);
            let pos = fm.backward_search_bitenc(bitenc::BitEnc::from_dna_string("AAA"));
            assert!(pos == [8, 9] || pos == [9, 8]);

            // Search for something at the very beginning.
            let pos = fm.backward_search(b"TGCATT".iter());
            assert_eq!(pos, [0]);
            let pos = fm.backward_search_bitenc(bitenc::BitEnc::from_dna_string("TGCATT"));
            assert_eq!(pos, [0]);

            // Test for things that are not there.
            let pos = fm.backward_search(b"AAAAAAA".iter());
            assert_eq!(pos, []);

            // Test that the BWT properly loops over the end.
            let dna = "AATTAA";
            let bitenc = bitenc::BitEnc::from_dna_string(dna);
            let bwt1 = compute_bwt(&bitenc, 2, 1);
            let fm = FMIndex::new(bwt1.bwt, bwt1.sa, 1, 3);
            let pos = fm.backward_search(b"AAA".iter());
            assert!(pos == [4, 5] || pos == [5, 4]);
            let pos = fm.backward_search(b"AAAATT".iter());
            assert!(pos == [4]);
        }
    }

    #[test]
    fn test_merge_bwt_buckets()
    {
        for npos in 2..9 {
            for iter in 0..2 {
                let dna = sim_tests::random_dna(5000);
                check_merge_bwt_buckets(&dna, npos);
            }
        }
    }

    fn check_merge_bwt_buckets(dna: &Vec<u8>, npos: usize) {
        // Test runs of the same character
        let sa_val = 2;
        let bitenc = bitenc::BitEnc::from_bytes(dna);

        let nbuckets = (2 as usize).pow(npos as u32);

        let mut filenames = Vec::new();
        for i in 0..nbuckets {
            let f = format!("b{}.bin.gz", i);
            println!("fn: {}", f);
            filenames.push(f);
        }

        for i in 0..nbuckets {
            let bwt = compute_bwt_bucket(&bitenc, npos, i);
            write_bwt(&bwt, &filenames[i]);
        }
        println!("npos: {}", npos);

        let bwt = merge_bwt_buckets(&filenames, sa_val);

        let orig_bwt = compute_bwt(&bitenc, npos, sa_val);

        assert_eq!(orig_bwt.sa, bwt.sa);
        assert_eq!(orig_bwt.bwt, bwt.bwt);

        for filename in &filenames {
            fs::remove_file(filename).expect("remove");
        }
    }
}
