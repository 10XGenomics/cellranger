//
// Copyright (c) 2017 10x Genomics, Inc. All rights reserved.
//

use fxhash::FxHashMap;
use constants::UmiType;
use utils;
use std::collections::HashSet;

const VDJ_MAX_CDR3_LEN: usize = 80;
const VDJ_MIN_CDR3_LEN: usize = 26;
const NO_CDR3: &str = "None";
const NO_CDR3_ID: usize = 0;

const MIXING_THRESHOLD: f64 = 0.2;
const PURITY_FACTOR: usize = 5;
const WINNING_MARGIN: f64 = 1.0;
const MAX_BELIEF_CUTOFF: f64 = 0.3;

fn codon_to_aa(codon: &[u8]) -> u8 {
    assert!(codon.len() == 3);
    let aa = match codon {
        b"GGT" => b'G',
        b"GGC" => b'G',
        b"GGA" => b'G',
        b"GGG" => b'G',
        b"TGG" => b'W',
        b"TGT" => b'C',
        b"TGC" => b'C',
        b"TTT" => b'F',
        b"TTC" => b'F',
        b"TTA" => b'L',
        b"TTG" => b'L',
        b"CTT" => b'L',
        b"CTC" => b'L',
        b"CTA" => b'L',
        b"CTG" => b'L',
        b"ATT" => b'I',
        b"ATC" => b'I',
        b"ATA" => b'I',
        b"GTT" => b'V',
        b"GTC" => b'V',
        b"GTA" => b'V',
        b"GTG" => b'V',
        b"TCT" => b'S',
        b"TCC" => b'S',
        b"TCA" => b'S',
        b"TCG" => b'S',
        b"AGT" => b'S',
        b"AGC" => b'S',
        b"CCT" => b'P',
        b"CCC" => b'P',
        b"CCA" => b'P',
        b"CCG" => b'P',
        b"ACT" => b'T',
        b"ACC" => b'T',
        b"ACA" => b'T',
        b"ACG" => b'T',
        b"GCT" => b'A',
        b"GCC" => b'A',
        b"GCA" => b'A',
        b"GCG" => b'A',
        b"TAT" => b'Y',
        b"TAC" => b'Y',
        b"CAT" => b'H',
        b"CAC" => b'H',
        b"CAA" => b'Q',
        b"CAG" => b'Q',
        b"AAT" => b'N',
        b"AAC" => b'N',
        b"AAA" => b'K',
        b"AAG" => b'K',
        b"GAT" => b'D',
        b"GAC" => b'D',
        b"GAA" => b'E',
        b"GAG" => b'E',
        b"CGT" => b'R',
        b"CGC" => b'R',
        b"CGA" => b'R',
        b"CGG" => b'R',
        b"AGA" => b'R',
        b"AGG" => b'R',
        b"ATG" => b'M',
        b"TAG" => b'*',
        b"TAA" => b'*',
        b"TGA" => b'*',
        _ => panic!("Unexpected codon"),
    };
    return aa;
}

fn search_cdr3_in_seq(seq: &[u8]) -> Vec<String> {
    // Search the CDR3 signature in a sequence without guides from annotations.
    // This could lead to more false positive signature hits than the guided version.
    // Return value: A tuple (CDR3 amino-acid seq, start position in seq, end position)

    let min_cdr3_aas_len = VDJ_MIN_CDR3_LEN/3;
    let mut cdr3_aas = Vec::new();

    for frame in 0..3 {
        let mut aa: Vec<u8> = Vec::new();
        let mut i = frame;
        while i < seq.len().saturating_sub(3) {
            aa.push(codon_to_aa(&seq[i..(i+3)]));
            i += 3;
        }

        let mut fgxg_idx: usize;
        let mut fgxg_pos: usize;

        for idx in min_cdr3_aas_len..aa.len().saturating_sub(3) {
            // First try to find the end motif
            if (aa[idx] == b'F' || aa[idx] == b'W') && aa[idx + 1] == b'G' && aa[idx + 3] == b'G' {
                // The CDR3 includes the first F of the signature
                fgxg_idx = idx;
                fgxg_pos = frame + fgxg_idx * 3;

                let mut cys_idx = None;
                // Find the Cysteine closer to the end but past a minimum
                // number of amino acids
                let mut j: i32 = fgxg_idx as i32 - min_cdr3_aas_len as i32;
                while j >= 0 {
                    if aa[j as usize] == b'C' {
                        cys_idx = Some(j as usize);
                        break;
                    } else if aa[j as usize] == b'*' {
                        // break;
                    }
                    j -= 1;
                }
                match cys_idx {
                    Some(pos) => {
                        let cys_pos = frame + pos * 3;
                        if (fgxg_pos - cys_pos) < VDJ_MAX_CDR3_LEN {
                            // include the start of the FGXG
                            let mut tmp = Vec::new();
                            for k in pos..(fgxg_idx + 1) {
                                tmp.push(aa[k]);
                            }
                            cdr3_aas.push(tmp);
                        }
                    },
                    None => {},
                }
            }
        }
    }

    let mut result = Vec::new();
    for cdr3 in cdr3_aas {
        result.push(String::from_utf8(cdr3.clone()).unwrap());
    }
    result
}

fn map_seqs_to_cdr3(seqs: &Vec<String>) -> (Vec<String>, Vec<usize>, Vec<Vec<usize>>) {
    // Return list of CDR3s as a Vec<String>,
    // the index of cdr3 assigned to each seq as Vec<usize>
    // the list of paths assigned to each cdr3 as Vec<Vec<usize>>

    let mut all_cdr3s = Vec::new();
    all_cdr3s.push(NO_CDR3.to_string());
    let mut cdr3_count = Vec::new();
    cdr3_count.push(0usize);

    let mut all_cdr3s_of_path = vec![Vec::new(); seqs.len()];
    let mut unique_cdr3_of_path = vec![NO_CDR3_ID; seqs.len()];
    let mut index_of_cdr3 = FxHashMap::default();


    for (path_id, ref seq) in seqs.iter().enumerate() {
        let seq_as_bytes = seq.as_bytes();
        let mut cdr3s = search_cdr3_in_seq(&seq_as_bytes);
        cdr3s.sort();
        for &ref cdr3 in cdr3s.iter() {
            if index_of_cdr3.contains_key(cdr3) {
                let idx = *index_of_cdr3.get(cdr3).unwrap();
                all_cdr3s_of_path[path_id].push(idx);
                cdr3_count[idx] += 1usize;
            } else {
                all_cdr3s_of_path[path_id].push(all_cdr3s.len());
                index_of_cdr3.insert(cdr3.clone(), all_cdr3s.len());
                all_cdr3s.push(cdr3.clone());
                cdr3_count.push(1usize);
            }
        }
    }
    for (path_id, cdr3s) in all_cdr3s_of_path.iter().enumerate() {
        if cdr3s.is_empty() { continue; }
        let mut tmp = cdr3s.clone();
        tmp.sort_by_key(|x| -(cdr3_count[*x] as isize));
        unique_cdr3_of_path[path_id] = tmp[0];
    }


    let mut paths_of_cdr3 = vec![Vec::new(); all_cdr3s.len()];

    for (path_id, &cdr3_id) in unique_cdr3_of_path.iter().enumerate() {
        paths_of_cdr3[cdr3_id].push(path_id);
    }

    (all_cdr3s, unique_cdr3_of_path, paths_of_cdr3)
}

#[derive(Debug, Clone)]
struct Belief {
    total: f64,
    max: f64
}

impl Belief {
    pub fn new(_total: f64, _max: f64) -> Self {
        Belief {
            total: _total,
            max: _max
        }
    }
    pub fn max(&self) -> f64 {
        self.max
    }
    pub fn total(&self) -> f64 {
        self.total
    }
    pub fn update(&mut self, score: f64) {
        self.total += score;
        self.max = if score > self.max { score } else { self.max };
    }
}

fn compute_community_beliefs(umi_cdr3_score: &FxHashMap<UmiType, FxHashMap<usize, f64> >,
                                num_cdr3s: usize) -> Vec<Belief> {
    let mut comm_belief = vec![Belief::new(0.0f64, 0.0f64); num_cdr3s];
    for (_, cdr3_score) in umi_cdr3_score.iter() {
        for (&cdr3, &score) in cdr3_score {
            comm_belief[cdr3].update(score);
        }
    }
    comm_belief
}

fn compute_normalized_umi_cdr3_scores(all_umi_path_scores: &FxHashMap<UmiType, Vec<(usize, f64)>>,
                                      cdr3_of_path: &Vec<usize>
                                      ) -> FxHashMap<UmiType, FxHashMap<usize, f64>> {
    let mut umi_cdr3_scores: FxHashMap<UmiType, FxHashMap<usize, f64> > = FxHashMap::default();
    for (umi, ref path_scores) in all_umi_path_scores {
        let mut cdr3_score = FxHashMap::default();
        let mut sum_scores = 0.0f64;
        for &(path_id, score) in path_scores.iter() {
            let cdr3 = cdr3_of_path[path_id];
            if cdr3_score.contains_key(&cdr3) {
                if cdr3 == NO_CDR3_ID {
                    if *cdr3_score.get_mut(&cdr3).unwrap() < score {
                        *cdr3_score.get_mut(&cdr3).unwrap() = score;
                    }
                } else {
                    *cdr3_score.get_mut(&cdr3).unwrap() += score;
                }
            } else {
                cdr3_score.insert(cdr3, score);
            }
            sum_scores += score;
        }
        // Normalize the scores
        for (_, score) in &mut cdr3_score {
            *score /= sum_scores;
        }
        umi_cdr3_scores.insert(*umi, cdr3_score);
    }
    umi_cdr3_scores
}


// Returns a Vector of tuples (UMI, assigned path id)
pub fn assign_paths_to_umis(all_path_seqs: &Vec<String>,
            all_umi_path_scores: &FxHashMap<UmiType, Vec<(usize, f64)> >,
            score_factor: f64) -> Vec<(UmiType, usize)> {

    let (cdr3_list, cdr3_of_path, paths_of_cdr3) = map_seqs_to_cdr3(&all_path_seqs);
    assert!(cdr3_list[NO_CDR3_ID] == NO_CDR3);
    let umi_cdr3_scores = compute_normalized_umi_cdr3_scores(&all_umi_path_scores, &cdr3_of_path);
    // Community beliefs contain the total and max score of each CDR3 across all UMIs
    let community_beliefs = compute_community_beliefs(&umi_cdr3_scores, paths_of_cdr3.len());

    // We will consider only CDR3s which have a max community belief more that max_belief_cutoff
    let mut candidate_cdr3s = Vec::new();
    let mut blacklisted_cdr3s = HashSet::new();
    for (cdr3_id, belief) in community_beliefs.iter().enumerate() {
        if cdr3_id == NO_CDR3_ID { // Ignore the NO_CDR3 now
            continue;
        }
        if belief.max() > MAX_BELIEF_CUTOFF {
            candidate_cdr3s.push( (cdr3_id, belief.total()) );
        } else {
            blacklisted_cdr3s.insert(cdr3_id);
        }
    }

    // Sort candidates by decreasing total belief
    candidate_cdr3s.sort_by_key(|x| utils::NonNan::new(-x.1 as f32).unwrap());

    for (i, &(cdr3, _)) in candidate_cdr3s.iter().enumerate() {
        if blacklisted_cdr3s.contains(&cdr3) { continue; }
        for (j, &(other_cdr3, _)) in candidate_cdr3s.iter().enumerate() {
            if j <=i { continue; }
            let (is_mixed, winner) = check_mixing(cdr3, other_cdr3, &umi_cdr3_scores, &community_beliefs);
            if is_mixed {
                match winner {
                    Some(winner_cdr3) => {
                        let loser_cdr3 = if winner_cdr3==cdr3 { other_cdr3 } else { cdr3 };
                        blacklisted_cdr3s.insert(loser_cdr3);
                    },
                    None => {
                        blacklisted_cdr3s.insert(cdr3);
                        blacklisted_cdr3s.insert(other_cdr3);
                    }
                }
            }
        }
    }

    let mut assigned_cdr3_of_umi: FxHashMap<UmiType, usize> = FxHashMap::default();
    let mut blacklisted_umis: HashSet<UmiType> = HashSet::new();
    let mut max_path_score_of_umi: FxHashMap<UmiType, f64> = FxHashMap::default();

    for (&umi, path_scores) in all_umi_path_scores {
        let mut max_path_score = 0.0f64;

        for &(_, score) in path_scores.iter() {
            if score > max_path_score {
                max_path_score = score;
            }
        }
        max_path_score_of_umi.insert(umi, max_path_score);

        let mut na_assignable = false;
        let mut supports_blacklisted = false;
        let mut assignable_cdr3s = Vec::new();
        for &(path_id, score) in path_scores.iter() {
            if score > score_factor * max_path_score { // Path is assignable
                let cdr3 = cdr3_of_path[path_id];
                if cdr3 == NO_CDR3_ID {
                    na_assignable = true;
                } else {
                    if blacklisted_cdr3s.contains(&cdr3) {
                        supports_blacklisted = true;
                    } else if *umi_cdr3_scores.get(&umi).unwrap().get(&cdr3).unwrap() > 0.5 {
                        assignable_cdr3s.push(cdr3);
                    }
                }
            }
        }
        if assignable_cdr3s.is_empty() {
            if na_assignable && !supports_blacklisted {
                assigned_cdr3_of_umi.insert(umi, NO_CDR3_ID);
            } else {
                blacklisted_umis.insert(umi);
            }
        } else {
            assignable_cdr3s.sort_by_key(|x| utils::NonNan::new(-community_beliefs[*x].total() as f32).unwrap());
            assigned_cdr3_of_umi.insert(umi, assignable_cdr3s[0]);
        }
    }

    println!(" Mixture Filter blacklisted UMIs: {:?}", blacklisted_umis);

    let mut path_support = vec![0.0f64; all_path_seqs.len()];
    for (&umi, path_scores) in all_umi_path_scores {

        if blacklisted_umis.contains(&umi) { continue; }

        // All non-blacklisted UMIs have a cdr3 assigned
        let assigned_cdr3 = *assigned_cdr3_of_umi.get(&umi).unwrap();

        for &(path_id, score) in path_scores.iter() {
            if assigned_cdr3 == cdr3_of_path[path_id] {
                path_support[path_id] += score;
            }

        }
    }

    let mut assigned_paths_of_umis = Vec::new();
    'umiloop: for (&umi, path_scores) in all_umi_path_scores {

        if blacklisted_umis.contains(&umi) { continue; }

        let mut sorted_path_scores = path_scores.clone();
        sorted_path_scores.sort_by_key(|x| utils::NonNan::new(-path_support[x.0] as f32).unwrap());
        let assigned_cdr3 = *assigned_cdr3_of_umi.get(&umi).unwrap();
        let max_score = *max_path_score_of_umi.get(&umi).unwrap();

        for &(path_id, score) in sorted_path_scores.iter() {
            if (assigned_cdr3 == cdr3_of_path[path_id]) &&  (score > score_factor * max_score) {
                assigned_paths_of_umis.push((umi, path_id));
                continue 'umiloop;
            }
        }
    }

    assigned_paths_of_umis
}

fn check_mixing(cdr3_1: usize,
                cdr3_2: usize,
                umi_cdr3_scores: &FxHashMap<UmiType, FxHashMap<usize, f64> >,
                community_beliefs: &Vec<Belief>) -> (bool, Option<usize>) {

    let mut n1 = 0; // UMIs supporting only cdr3_1
    let mut n2 = 0; // UMIs supporting only cdr3_2
    let mut n12 = 0; // UMIs supporting both cdr3_1 and cdr3_2

    for (_, cdr3_score) in umi_cdr3_scores.iter() {
        let umi_supports_1: bool  = (cdr3_score.contains_key(&cdr3_1)) && (*cdr3_score.get(&cdr3_1).unwrap() > MIXING_THRESHOLD);
        let umi_supports_2: bool  = (cdr3_score.contains_key(&cdr3_2)) && (*cdr3_score.get(&cdr3_2).unwrap() > MIXING_THRESHOLD);
        if umi_supports_1 && umi_supports_2 {
            n12 += 1;
        } else if umi_supports_1 {
            n1 += 1;
        } else if umi_supports_2 {
            n2 += 1;
        }
    }
    let mixed = (n1 < PURITY_FACTOR * n12) || (n2 < PURITY_FACTOR * n12);
    let winner = {
        let margin = community_beliefs[cdr3_1].total() / community_beliefs[cdr3_2].total();
        if margin > WINNING_MARGIN {
            Some(cdr3_1)
        } else if margin < (1.0f64 / WINNING_MARGIN) {
            Some (cdr3_2)
        } else {
            None
        }
    };
    (mixed, winner)
}
