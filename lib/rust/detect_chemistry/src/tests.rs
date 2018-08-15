//
// Copyright (c) 2017 10x Genomics, Inc. All rights reserved.
//

use rand::Rng;
use rand;
use debruijn::bits_to_base;

#[allow(dead_code)]
pub fn random_seq(len: usize) -> Vec<u8> {
    // Random ascii string
    let mut r = rand::thread_rng();
    let mut seq = Vec::new();
    for _ in 0..len {
        seq.push(bits_to_base((r.gen::<u64>() % 4) as u8) as u8);
    }
    seq
}

#[allow(dead_code)]
pub fn random_seq_rng(len: usize, r: &mut impl Rng) -> Vec<u8> {
    // Random ascii string
    let mut seq = Vec::new();
    for _ in 0..len {
        seq.push(bits_to_base((r.gen::<u64>() % 4) as u8) as u8);
    }
    seq
}
