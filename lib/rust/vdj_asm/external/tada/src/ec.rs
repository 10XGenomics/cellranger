//
// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
//

//! Utilities for building DeBruijn graphs efficiently
use kmer::{K, Kmer};
//use std::collections::{VecDeque, HashSet};

use utils::{DMap};
use std::cmp::min;
/*
#[derive(Debug, Clone, PartialOrd, Ord, Eq, PartialEq)]
struct EcS {
    // Current kmer
    kmer: Kmer,

    // Previous state -- for traceback
    prev_state: usize,

    // Position of cursor
    pos: usize,
}


fn simple_options(state: EcS, prev_idx: usize, prev_score: isize, read: &Vec<u8>, buf: &mut Vec<(f32, EcS)>, kmers: &DMap<Kmer, u32>)
{
    let next_pos = state.pos + 1;
    if next_pos == read.len() {
        return;
    }

    let total_extend_counts = 0.0;
    for i in 0 .. 4
    {
        let next_kmer = state.kmer.extend_right(i);

        // Kmer score
        let (min_kmer, _) = next_kmer.min_rc();
        let kmer_counts = kmers.get(&min_kmer).map_or(0, |x| *x);
        total_extend_counts += kmer_counts as f32;
    }

    for i in 0 .. 4
    {
        let next_kmer = state.kmer.extend_right(i);
        let pos = state.pos + 1;

        // Penalties are calibrated so that a 'lone' kmer will be corrected to
        // a nearest-neighbor kmer if there are 4 or more observations of the neighbor
        let mut_penalty = if i == *read.get(pos).expect("oob") { -1.0 } else { 2.2 };

        // Kmer score
        let (min_kmer, _) = next_kmer.min_rc();
        let counts = kmers.get(&min_kmer).map_or(0, |x| *x) as f32;

        let count_score = (counts - 1.0 + 0.1).ln(); // - total_extend_counts.ln();

        let score = prev_score +  mut_penalty - count_score;

        let new_state = EcS { kmer: next_kmer, prev_state: prev_idx, pos: pos };
        buf.push((score, new_state));
    }
}

fn step(k: usize, states: &Vec<(isize, EcS)>, read: &Vec<u8>, kmers: &DMap<Kmer, Vec<u32>>) -> Vec<(isize, EcS)> {
    let mut new_states = Vec::new();

    for (idx, (prev_score, prev_state)) in states.iter().cloned().enumerate() {
        simple_options(prev_state, idx, prev_score, read, &mut new_states, kmers)
    }

    new_states.sort();
    new_states.truncate(k);
    new_states
}


/// Feed in seq in sequencer order -- starting with good quality reads is better!
pub fn correct_seq(read: &Vec<u8>, kmers: &DMap<Kmer, u32>) -> Vec<u8> {
    // find an edit of seq that makes it use only 'good' Kmers

    let mut first_kmer = Kmer::empty();
    for i in 0 .. K {
        first_kmer = first_kmer.set(i, read[i]);
    }

    println!("first kmer: {:?}", first_kmer);
    // Assume first kmer is correct
    if !kmers.contains_key(&first_kmer.min_rc().0) {
        println!("nope");
        return read.clone()
    }

    let init = EcS { kmer: first_kmer, pos: K-1, prev_state: 0 };

    // Initial state
    let mut states = vec![(0, init)];

    // trackback
    let mut history = Vec::new();

    loop {
        let new_states = step(64, &states, &read, kmers);
        history.push(states);

        if new_states.len() == 0 {
            break;
        }
        println!("{:?}", new_states[0]);
        states = new_states;
    }


    // reconstruct sequence
    let mut bp = Vec::new();
    let mut last_state = history[history.len() - 1].get(0).expect("history").clone();
    for i in 0 .. K {
        bp.push(last_state.1.kmer.get(K-1-i));
    }

    for i in (0..history.len()-1).rev() {
        let next_state = history[i].get(last_state.1.prev_state).expect("history");
        bp.push(next_state.1.kmer.get(0));
        last_state = next_state.clone();
    }

    let mut res = Vec::new();
    for b in bp.iter().rev() {
        res.push(*b);
    }
    res
}
*/
pub fn best_middle_base(qv: u8, kmer: Kmer, kmer_counts: &DMap<Kmer, u32>) -> u8 {
    let counts = kmer_counts.get(&kmer.min_rc().0).map_or(0, |x| *x);
    if counts >= 3 {
        return kmer.get(K/2)
    }

    let current = kmer.get(K/2);
    let mut ops : [f32; 4] = [0.0;4];
    for i in 0 .. 4 {
        let (cand_kmer, _) = kmer.set(K/2, i).min_rc();
        let cand_counts = kmer_counts.get(&cand_kmer).map_or(0, |x| *x);
        ops[i as usize] = cand_counts as f32;
    }

    let n = (ops[0] + ops[1] + ops[2] + ops[3]) as f32;


    let p_err = (10.0 as f32).powf(-(min(20, qv) as f32) / 10.0);
    let mut pmax = -1.0;
    let mut best_opt = 0;
    for i in 0 .. 4 {
        let posterior =
            if i == current {
                (1.0 -  p_err) * (ops[i as usize] - 1.0 + 0.05) / n
            }
            else
            {
                p_err * ops[i as usize] / n
            };

        if posterior > pmax {
            pmax = posterior;
            best_opt = i;
        }
    }

    best_opt
}


/// Feed in seq in sequencer order -- starting with good quality reads is better!
pub fn correct_seq(read: &Vec<u8>, qvs: &Vec<u8>, kmer_counts: &DMap<Kmer, u32>) -> Vec<u8> {
    // find an edit of seq that makes it use only 'good' Kmers

    let mut kmer = Kmer::empty();
    for i in 0 .. K {
        kmer = kmer.set(i, read[i]);
    }

    let mut new_read = read.clone();

    for i in (K/2) .. read.len() - K/2 {
        let best_option = best_middle_base(qvs[i], kmer, kmer_counts);
        new_read[i] = best_option;

        if i + K/2 + 1 < read.len() {
            kmer = kmer.extend_right(read[i + K/2 + 1]);
        }
    }

    new_read
}
