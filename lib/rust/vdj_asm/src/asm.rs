//
// Copyright (c) 2017 10x Genomics, Inc. All rights reserved.
//

use debruijn::Kmer;
use debruijn::fx::FxHashMap;
use graph::IndexedGraph;
use graph::Kmer1;
use graph::ReadDb;

use fastq;
use graph_read;
use graph;
use sw;
use utils;
use bam_utils;
use asm_helper;
use constants::{UmiType, ReadType, QUAL_OFFSET, MAX_NUM_KMERS, PROCESSED_UMI_TAG, KMER_LEN_BANDED_ALIGN, WINDOW_SIZE_BANDED_ALIGN};
use std::io::{Write};
use rust_htslib::bam;

use std::collections::{HashSet, HashMap};

fn count_support(rough_alignments: &FxHashMap<ReadType, (usize, f64)>,
                 read_db: &graph::ReadDb, min_sw_score: f64) ->
                 (HashMap<UmiType, (f64, usize, usize)>, HashSet<ReadType>) {

    let mut umi_path_scores = HashMap::new();

    let min_score = -1e6;
    let max_score = 1e6;

    let mut support_reads = HashSet::new();

    for (read_id, &(_, score)) in rough_alignments.iter() {
        if score < min_sw_score {
            continue;
        }
        let is_read1 = read_id % 2 == 0;

        let umi = read_db.reads.get(&read_id).unwrap().umi;
        if !umi_path_scores.contains_key(&umi) {
            umi_path_scores.insert(umi, (score, is_read1 as usize, !is_read1 as usize));
        } else {
            let val = umi_path_scores.get_mut(&umi).unwrap();
            let new_score;
            if (*val).0 > min_score && (*val).0 < max_score {
                new_score = (*val).0 + score;
            } else {
                new_score = (*val).0;
            }
            let new_read1 = (*val).1 + (is_read1 as usize);
            let new_read2 = (*val).2 + (!is_read1 as usize);

            *val = (new_score, new_read1, new_read2);
        }
        if score >= min_score {
            support_reads.insert(*read_id);
        }
    }

    (umi_path_scores, support_reads)
}


pub struct UmiCounter {
    pub map: HashMap<String, UmiType>,
    pub counts: HashMap<UmiType, usize>,
}

impl UmiCounter {
    pub fn new() -> UmiCounter {

        let mut umi_map = HashMap::new(); // UMI string -> UMI id
        umi_map.insert("".to_string(), 0 as UmiType);

        let mut umi_counts = HashMap::new();
        umi_counts.insert(0 as UmiType,0);

        UmiCounter { map: umi_map, counts: umi_counts }
    }

    pub fn len(&self) -> usize {
        self.counts.len()
    }

    pub fn get_umi(&mut self, header: &fastq::CellrangerFastqHeader) -> UmiType {

        let umi =
            match header.get_tag(PROCESSED_UMI_TAG) {
                Some(umi_str) => {
                    if !self.map.contains_key(umi_str) {
                        let new_id = self.map.len() as UmiType;
                        self.counts.insert(new_id, 0);
                        self.map.insert(umi_str.clone(), new_id);
                        new_id
                    } else {
                        self.map.get(umi_str).unwrap().clone()
                    }
                },
                None => 0 as UmiType
            };

        umi
    }

    pub fn count_umi(&mut self, umi: UmiType) {
        *(self.counts.get_mut(&umi).unwrap()) += 1;
    }

    pub fn get_good_umis(&self, min_umi_reads: usize) -> HashSet<UmiType> {
        self.counts.iter()
            .filter(|&(umi, count)| *count >= min_umi_reads && *umi != 0)
            .map(|(u, _)| *u)
            .collect()
    }

    pub fn get(&self, umi: &UmiType) -> Option<usize> {
        self.counts.get(umi).cloned()
    }
}


struct AssemblyParams {
    min_kmer_count: u32,
    qual_factor: f64,
    paths_per_component: usize,
    rt_error: f64,
    min_qual: u8,
    min_contig_len: usize,
    frac_path_reads: f64,
    min_sw_score: f64,
    score_factor: f64,
}


/// Assemble a set of reads and return assemblies and alignments of reads against them.
///
/// The algorithm works roughly as follows:
/// 1. Go through each connected component of the graph.
/// 2. For each component, sort edges by UMI support.
/// 3. Get the strongest edge, and get "high quality" paths from that edge (extending in both
/// directions). "High quality" means that branching points are followed only
/// if the branching base has at least some minimum quality.
/// 4. Invalidate all nodes in the returned paths.
/// 5. Get the next strongest (and not invalidated) edge and repeat 3.
/// 6. Repeat 3-5 npaths times.
/// 7. Sort paths by read support.
/// 8. Assign UMIs to paths: start from the strongest path and assign each UMI to that
/// path if the total read score of reads aligned against that path is within
/// score_factor from the total score of reads aligned againt the best path for that UMI.
/// 9. Remove paths that had no UMIs.
/// 10. Remove paths with <N99 of reads among all paths.
pub fn assemble_reads(_reads: &mut Vec<&graph_read::Read>,
                      umi_counts: &UmiCounter,
                      min_kmer_count: usize,
                      qual_factor: f64,
                      rt_error: f64,
                      min_qual: u8,
                      min_contig_len: usize,
                      frac_path_reads: f64,
                      min_align_score: f64,
                      score_factor: f64,
                      cons_only: bool,
                      scoring: sw::Scoring,
                      is_paired: bool,
                      frac_starts: f64,
                      ) -> (IndexedGraph<Kmer1, Vec<u32>>, Vec<(String, Vec<u8>, HashSet<UmiType>, HashMap<ReadType, sw::AlignmentPacket>)>) {

    let mut readsc: Vec<&graph_read::Read> = Vec::new();
    for r in _reads {
        if r.len() >= graph::Kmer1::k() {
            readsc.push(r);
        }
    }
    let reads = &mut readsc;
    reads.sort_by_key(|x| x.umi);

    let barcoded_graph = graph::build_graph(&reads, min_kmer_count, MAX_NUM_KMERS);
    let read_db = barcoded_graph.build_read_db(&reads);
    let all_read_ids = read_db.reads.keys().map(|x| *x).collect::<HashSet<ReadType>>();

    let mut good_contigs : Vec<(String, // path sequence after trimming
                                Vec<u8>, // base qualities
                                HashSet<UmiType>, // Set of UMIs on the path
                                HashMap<ReadType, sw::AlignmentPacket> // read id -> Alignment
                                )> = Vec::new();

    let components = barcoded_graph.connected_components();

    for ref mut component in components.iter() {
        if cons_only && good_contigs.len() == 1 {
            break;
        }

        // ========== PART 1 - get candidate paths ===========

        // UMI -> (score of path, # reads on path)
        let mut best_paths : HashMap<UmiType, (f64, usize)> = HashMap::new();

        let mut component_paths: Vec<(String, // path sequence
                                      HashSet<ReadType>, // Reads mapping on the path
                                      HashMap<UmiType, (f64, usize, usize)>, // From UMI to scores and read counts on this path
                                      f64 // total path score
                                      )> = Vec::new();

        let mut sorted_edges = component.clone();

        // Sort edges by increasing support -- initialize by read support instead of UMI support
        // this should avoid getting stuck on tips when assembling from 1 UMI
        sorted_edges.sort_by_key(|x| read_db.mappings[*x].len());

        let mut visited = HashSet::new(); // visited nodes
        let mut ext_qual_lookup: FxHashMap<(usize, i8), Vec<(usize, u8)> > = FxHashMap::default();

        let mut total_support: usize = 0;
        for edge in &sorted_edges {
            total_support += read_db.mappings[*edge].len();
        }
        let max_used_support = (frac_starts * (total_support as f64)) as usize;
        // Terminate path enumeration after starting from 
        // nodes containing frac_starts fraction of the reads
        let mut used_support = 0;

        while !sorted_edges.is_empty() && used_support < max_used_support {

            // Get the edge with the strongest support (will be at the end of sorted_edges).
            let strongest_edge = sorted_edges.pop().unwrap();
            used_support += read_db.mappings[strongest_edge].len();

            if visited.contains(&strongest_edge) { continue; }
            
            let new_paths = barcoded_graph.max_paths_by_qual(strongest_edge, &read_db,
                                                             qual_factor as f32, rt_error,
                                                             min_qual, &scoring, &mut ext_qual_lookup);

            // Mark the nodes along the path as "visited"
            for &(ref path, _) in new_paths.iter() {
                for &node in path.iter() {
                    visited.insert(node);
                }
            }

            for &(ref new_path, ref rough_alignments) in new_paths.iter() {
                let contig_seq = barcoded_graph.get_path_sequence(&new_path);
                if contig_seq.len() < min_contig_len {
                    continue;
                }

                let (umi_path_scores, rough_support_reads) = count_support(&rough_alignments, &read_db, min_align_score);

                let mut total_path_score = 0.0;

                // Keep track of the best score for each UMI. Also compute the total score of this
                // path across all good UMIs.
                for (&umi, &(score, _, _)) in umi_path_scores.iter() {

                    total_path_score = total_path_score + score;

                    if !best_paths.contains_key(&umi){
                        best_paths.insert(umi, (score, component_paths.len()));
                    } else if best_paths.get(&umi).unwrap().0 < score {
                        let val = best_paths.get_mut(&umi).unwrap();
                        *val = (score, component_paths.len());
                    }
                }

                component_paths.push((contig_seq, rough_support_reads, umi_path_scores, total_path_score));

            }
        }

        // ========== PART 2 - assign UMIs to paths ===========

        // Sort paths by decreasing total score
        component_paths.sort_by_key(|x| utils::NonNan::new(-x.3 as f32).unwrap());

        let mut assigned_umis : HashSet<UmiType> = HashSet::new();

        for &(ref contig_seq, ref rough_support_reads, ref umi_path_scores, _) in component_paths.iter() {
            let mut path_umis = HashSet::new();

            for (&umi, &(score, count1, count2)) in umi_path_scores.iter() {
                let total_umi_counts = umi_counts.get(&umi).unwrap();
                let count_frac1 = (count1 as f32) / (total_umi_counts as f32);
                let count_frac2 = (count2 as f32) / (total_umi_counts as f32);
                let total_count_frac = count_frac1 + count_frac2;

                // If this is not paired, require a fraction of all reads to map to the contig.
                // If this is paired, require a fraction of reads from each side to map to the contig.
                if (!is_paired && total_count_frac < frac_path_reads as f32) ||
                   (is_paired && (count_frac1 < 0.5 * (frac_path_reads as f32) || (count_frac2 < 0.5 * (frac_path_reads as f32)))) {
                    continue;
                }
                // If this is the optimal path for the UMI or close to the optimal, then assign it to this path.
                if score_factor * best_paths.get(&umi).unwrap().0 <= score && !assigned_umis.contains(&umi) {
                    path_umis.insert(umi);
                    // Can't use the same UMI again.
                    assigned_umis.insert(umi);
                }
            }

            if path_umis.is_empty() { continue; }

            let mut read_alignments = HashMap::new(); // read id -> Alignment object

            let refs = vec![contig_seq.clone()];
            let align_helper = sw::AlignHelper::new(&refs, scoring.clone(), KMER_LEN_BANDED_ALIGN, WINDOW_SIZE_BANDED_ALIGN);

            let mut good_alignments = Vec::new();
            let mut good_alignment_reads = Vec::new();

            let reads_to_align = match cons_only {
                true => &all_read_ids, // ALL reads
                false => rough_support_reads,
            };

            // Properly align the selected reads on the contig and compute quality scores.
            for read_id in reads_to_align.iter() {
                let read = read_db.reads.get(read_id).unwrap();
                if path_umis.contains(&read.umi) {
                    let wrapped_alignment = align_helper.find_read_matches(&read, min_align_score as i32);
                    match wrapped_alignment {
                        Some(al) => {
                            read_alignments.insert(*read_id, al.clone());
                            good_alignments.push((al, read));
                            good_alignment_reads.push(read_id);
                        },
                        None => {},
                    }
                }
            }

            if good_alignments.is_empty() { continue; }

            // Since we realign the reads and recompute the alignment scores, some reads that were
            // previously aligned might get dropped. So this can change the set of UMIs on this contig.
            let final_umis : HashSet<UmiType> = good_alignments.iter().map(|x| x.1.umi).collect();

            let pileup = align_helper.pileup(0, &good_alignments);

            let quals = align_helper.base_quals(0, &pileup, rt_error);
            assert!(quals.len() == contig_seq.len());

            let (start_pos, end_pos) = trim(&quals, 0);
            if end_pos <= 0 || start_pos >= contig_seq.len() {
                continue;
            }
            //let trimmed_quals = quals[start_pos..(end_pos + 1)].to_vec();
            let trimmed_seq = contig_seq[start_pos..(end_pos + 1)].to_string();
            if trimmed_seq.len() < min_contig_len {
                continue;
            }

            // Output the untrimmed sequence. We'd have to realign otherwise.
            good_contigs.push((contig_seq.to_string(), quals, final_umis.clone(), read_alignments));
            if cons_only && good_contigs.len() == 1 {
                break;
            }
        }
    }

    let nx_mapped_reads = utils::nx_count(&good_contigs.iter().map(|x| x.3.len() as usize).collect(), 99.0_f64);

    good_contigs = good_contigs.iter().filter(|&x| x.3.len() >= nx_mapped_reads).map(|x| x.clone()).collect();

    (barcoded_graph, good_contigs)
}


pub fn assemble_reads_with_mixture_filter(_reads: &mut Vec<&graph_read::Read>,
                      umi_counts: &UmiCounter,
                      min_kmer_count: usize,
                      qual_factor: f64,
                      rt_error: f64,
                      min_qual: u8,
                      min_contig_len: usize,
                      min_align_score: f64,
                      score_factor: f64,
                      scoring: sw::Scoring,
                      frac_starts: f64,
                      frac_path_reads: f64,
                      is_paired: bool
                      ) -> (IndexedGraph<Kmer1, Vec<u32>>, Vec<(String, Vec<u8>, HashSet<UmiType>, HashMap<ReadType, sw::AlignmentPacket>)>) {

    let mut readsc: Vec<&graph_read::Read> = Vec::new();
    for r in _reads {
        if r.len() >= graph::Kmer1::k() {
            readsc.push(r);
        }
    }
    let reads = &mut readsc;
    reads.sort_by_key(|x| x.umi);

    let barcoded_graph = graph::build_graph(&reads, min_kmer_count, MAX_NUM_KMERS);
    let read_db = barcoded_graph.build_read_db(&reads);

    // STEP 1 : Enumerate all paths and assign scores for valid (UMI, path) pairs
    let (all_path_seqs, all_umi_path_scores) = enumerate_paths(&read_db,
                                                    &umi_counts,
                                                    &barcoded_graph,
                                                    qual_factor,
                                                    rt_error,
                                                    min_qual,
                                                    min_contig_len,
                                                    min_align_score,
                                                    scoring.clone(),
                                                    frac_starts,
                                                    frac_path_reads,
                                                    is_paired);
    
    // STEP 2 : Group paths by CDR3 and assign scores for (UMI, CDR3) pairs
    // by summing over all the (UMI, path) scores for paths in the group
    // defined by the CDR3. Eliminate mixed CDR3 and assign paths to UMIs
    let assigned_paths_of_umis = asm_helper::assign_paths_to_umis(&all_path_seqs, &all_umi_path_scores, score_factor);

    let mut assigned_umis_of_paths: Vec<Vec<UmiType>> = vec![Vec::new(); all_path_seqs.len()];
    for &(umi, path_id) in assigned_paths_of_umis.iter() {
        assigned_umis_of_paths[path_id].push(umi);
    }

    // STEP 3 : Align reads from each UMI to the contig assigned to the UMI
    let good_contigs = align_reads(&read_db,
                            rt_error,
                            min_contig_len,
                            min_align_score,
                            scoring,
                            &all_path_seqs,
                            &assigned_umis_of_paths);
    
    (barcoded_graph, good_contigs)
}


fn enumerate_paths(read_db: &ReadDb,
                       umi_counts: &UmiCounter,
                       barcoded_graph: &IndexedGraph<Kmer1, Vec<u32>>,
                       qual_factor: f64,
                       rt_error: f64,
                       min_qual: u8,
                       min_contig_len: usize,
                       min_align_score: f64,
                       scoring: sw::Scoring,
                       frac_starts: f64,
                       frac_path_reads: f64,
                       is_paired: bool
                       ) -> (Vec<String>, FxHashMap<UmiType, Vec<(usize, f64)> >) {

    let components = barcoded_graph.connected_components();
    let mut path_id = 0;
    let mut all_path_seqs = Vec::new();
    let mut all_umi_path_scores: FxHashMap<UmiType, Vec<(usize, f64)> > = FxHashMap::default(); // UMI -> vector of (path_id, score)

    for ref mut component in components.iter() {

        // ========== PART 1 - Enumerate paths ===========

        let mut component_paths: Vec<(usize, // path sequence id
                                      HashSet<ReadType>, // Reads mapping on the path
                                      HashMap<UmiType, (f64, usize, usize)>, // From UMI to scores and read counts on this path
                                      )> = Vec::new();

        let mut sorted_edges = component.clone();

        // Sort edges by increasing support -- initialize by read support instead of UMI support
        // this should avoid getting stuck on tips when assembling from 1 UMI
        sorted_edges.sort_by_key(|x| read_db.mappings[*x].len());

        let mut visited = HashSet::new(); // visited nodes
        let mut ext_qual_lookup: FxHashMap<(usize, i8), Vec<(usize, u8)> > = FxHashMap::default();

        let mut total_support: usize = 0;
        for edge in &sorted_edges {
            total_support += read_db.mappings[*edge].len();
        }
        let max_used_support = (frac_starts * (total_support as f64)) as usize;
        // Terminate path enumeration after starting from 
        // nodes containing frac_starts fraction of the reads
        let mut used_support = 0;

        while !sorted_edges.is_empty() && used_support < max_used_support {

            // Get the edge with the strongest support (will be at the end of sorted_edges).
            let strongest_edge = sorted_edges.pop().unwrap();
            used_support += read_db.mappings[strongest_edge].len();

            if visited.contains(&strongest_edge) { continue; }

            let new_paths = barcoded_graph.max_paths_by_qual(strongest_edge, &read_db,
                                                             qual_factor as f32, rt_error,
                                                             min_qual, &scoring, &mut ext_qual_lookup);

            // Mark the nodes along the path as "visited"
            for &(ref path, _) in new_paths.iter() {
                for &node in path.iter() {
                    visited.insert(node);
                }
            }

            for &(ref new_path, ref rough_alignments) in new_paths.iter() {
                let contig_seq = barcoded_graph.get_path_sequence(&new_path);
                if contig_seq.len() < min_contig_len {
                    continue;
                }

                let (umi_path_scores, rough_support_reads) = count_support(&rough_alignments, &read_db, min_align_score);
                all_path_seqs.push(contig_seq);
                component_paths.push((path_id, rough_support_reads, umi_path_scores));
                path_id += 1;
            }
        }

        // ========== PART 2 - UMI vs Path Matrix ===========

        for &(pid, _, ref umi_path_scores) in component_paths.iter() {

            for (&umi, &(score, count1, count2)) in umi_path_scores.iter() {

                let total_umi_counts = umi_counts.get(&umi).unwrap();
                let count_frac1 = (count1 as f32) / (total_umi_counts as f32);
                let count_frac2 = (count2 as f32) / (total_umi_counts as f32);
                let total_count_frac = count_frac1 + count_frac2;

                // If this is not paired, require a fraction of all reads to map to the contig.
                // If this is paired, require a fraction of reads from each side to map to the contig.
                if (!is_paired && total_count_frac < frac_path_reads as f32) ||
                   (is_paired && (count_frac1 < 0.5 * (frac_path_reads as f32) || (count_frac2 < 0.5 * (frac_path_reads as f32)))) {
                    continue;
                }

                if all_umi_path_scores.contains_key(&umi) {
                    all_umi_path_scores.get_mut(&umi).unwrap().push((pid, score));
                } else {
                    all_umi_path_scores.insert(umi, vec![(pid, score)]);
                }
            }

        }
    }
    (all_path_seqs, all_umi_path_scores)
}

fn align_reads(read_db: &ReadDb,
                   rt_error: f64,
                   min_contig_len: usize,
                   min_align_score: f64,
                   scoring: sw::Scoring,
                   all_path_seqs: &Vec<String>,
                   assigned_umis_of_paths: &Vec<Vec<UmiType>>
                   ) -> Vec<(String, Vec<u8>, HashSet<UmiType>, HashMap<ReadType, sw::AlignmentPacket>)> {

    let all_read_ids = read_db.reads.keys().map(|x| *x).collect::<HashSet<ReadType>>();

    let mut good_contigs : Vec<(String, // path sequence after trimming
                                Vec<u8>, // base qualities
                                HashSet<UmiType>, // Set of UMIs on the path
                                HashMap<ReadType, sw::AlignmentPacket> // read id -> Alignment
                                )> = Vec::new();

    let reads_to_align = &all_read_ids;
    for (path_id, ref umis) in assigned_umis_of_paths.iter().enumerate() {

        if umis.is_empty() { continue; }

        let mut path_umis = HashSet::new();
        for &umi in umis.iter() {
            path_umis.insert(umi);
        }

        let mut read_alignments = HashMap::new(); // read id -> Alignment object
        let contig_seq = &all_path_seqs[path_id];
        let refs = vec![contig_seq.clone()];

        let align_helper = sw::AlignHelper::new(&refs, scoring.clone(), KMER_LEN_BANDED_ALIGN, WINDOW_SIZE_BANDED_ALIGN);

        let mut good_alignments = Vec::new();
        let mut good_alignment_reads = Vec::new();

        // Properly align the selected reads on the contig and compute quality scores.
        for read_id in reads_to_align.iter() {
            let read = read_db.reads.get(read_id).unwrap();
            if path_umis.contains(&read.umi) {
                let wrapped_alignment = align_helper.find_read_matches(&read, min_align_score as i32);
                match wrapped_alignment {
                    Some(al) => {
                        read_alignments.insert(*read_id, al.clone());
                        good_alignments.push((al, read));
                        good_alignment_reads.push(read_id);
                    },
                    None => {},
                }
            }
        }

        if good_alignments.is_empty() { continue; }

        // Since we realign the reads and recompute the alignment scores, some reads that were
        // previously aligned might get dropped. So this can change the set of UMIs on this contig.
        let final_umis : HashSet<UmiType> = good_alignments.iter().map(|x| x.1.umi).collect();

        let pileup = align_helper.pileup(0, &good_alignments);

        let quals = align_helper.base_quals(0, &pileup, rt_error);
        assert!(quals.len() == contig_seq.len());

        let (start_pos, end_pos) = trim(&quals, 0);
        if end_pos <= 0 || start_pos >= contig_seq.len() {
            continue;
        }
        //let trimmed_quals = quals[start_pos..(end_pos + 1)].to_vec();
        let trimmed_seq = contig_seq[start_pos..(end_pos + 1)].to_string();
        if trimmed_seq.len() < min_contig_len {
            continue;
        }

        // Output the untrimmed sequence. We'd have to realign otherwise.
        good_contigs.push((contig_seq.to_string(), quals, final_umis.clone(), read_alignments));
    }

    let nx_mapped_reads = utils::nx_count(&good_contigs.iter().map(|x| x.3.len() as usize).collect(), 99.0_f64);
    good_contigs = good_contigs.iter().filter(|&x| x.3.len() >= nx_mapped_reads).map(|x| x.clone()).collect();

    good_contigs
}

pub struct AssemblyOuts<T> {
    pub fasta_writer: T,
    pub fastq_writer: T,
    pub summary_writer: T,
    pub umi_summary_writer: T,
}

impl<T: Write> AssemblyOuts<T> {
    pub fn new(fasta_writer: T,
               fastq_writer: T,
               summary_writer: T,
               umi_summary_writer: T,
    ) -> AssemblyOuts<T> {
        AssemblyOuts {
            fasta_writer: fasta_writer,
            fastq_writer: fastq_writer,
            summary_writer: summary_writer,
            umi_summary_writer: umi_summary_writer,
        }
    }

    pub fn close(&mut self) {
        let _ = self.fasta_writer.flush();
        let _ = self.fastq_writer.flush();
        let _ = self.summary_writer.flush();
        let _ = self.umi_summary_writer.flush();
    }
}

pub fn write_assembly_results<T>(reads: &Vec<graph_read::Read>,
                                 umi_counts: &UmiCounter,
                                 good_umis: &HashSet<UmiType>,
                                 writers: &mut AssemblyOuts<T>,
                                 out_bam: &mut bam::Writer,
                                 contigs: &Vec<(String, Vec<u8>, HashSet<UmiType>, HashMap<ReadType, sw::AlignmentPacket>)>,
                                 contig_names: &Vec<String>,
                                 barcode: &str,
                                 read_cutoff: usize,
                                 single_end: bool) where T:Write {


    let mut all_alignments = HashMap::new();

    // From UMI id to UMI sequence. Just for reporting purposes.
    let mut rev_umi_map = HashMap::new();
    for (umi_str, umi_id) in umi_counts.map.iter() {
        rev_umi_map.insert(umi_id, umi_str);
    }

    // UMI -> contigs assigned to the UMI
    let mut umi_contigs = Vec::with_capacity(umi_counts.map.len());
    for _ in 0..umi_counts.map.len() {
        let tmp : Vec<String> = Vec::new();
        umi_contigs.push(tmp);
    }

    // Write assembled sequences to FASTQ
    for (contig_idx, &(ref contig_seq, ref quals, ref path_umis, ref alignments)) in contigs.iter().enumerate() {

        let mut good_pairs = HashSet::new();
        for (read_id, _) in alignments.iter() {
            if *read_id % 2 == 0 && alignments.contains_key(&(read_id + 1)) {
                good_pairs.insert(read_id);
            }
        }

        let ref contig_name = contig_names[contig_idx];

        for &umi in path_umis.iter() {
            umi_contigs[umi as usize].push(contig_name.clone());
        }

        writers.fasta_writer.write_fmt(format_args!(">{}\n{}\n", contig_name, contig_seq)).unwrap();
        writers.fastq_writer.write_fmt(format_args!("@{}\n{}\n+\n{}\n", contig_name, contig_seq,
                                       fastq::get_qual_string(&quals, QUAL_OFFSET))).unwrap();

        let mut umi_list : Vec<UmiType> = path_umis.iter().map(|x| *x).collect();
        umi_list.sort();
        writers.summary_writer.write_fmt(format_args!("{}\t{}\t{:?}\t{:?}\t{:?}\t{}\n",
                                         barcode, contig_name,
                                         alignments.len(), good_pairs.len(), path_umis.len(),
                                                      utils::vec_str(&umi_list))).unwrap();

        for (read_id, al) in alignments.iter() {
            all_alignments.insert(read_id, (contig_idx, al.clone()));
        }
    }


    for umi_id in 0..umi_counts.map.len() {
        let umi = umi_id as u32;
        if !umi_counts.get(&umi).is_some() {
            continue;
        }
        let umi_str = rev_umi_map.get(&umi).unwrap();
        // Note that the output counts will be after sampling.
        writers.umi_summary_writer.write_fmt(format_args!("{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                                                          barcode, umi, umi_str,
                                                          umi_counts.get(&umi).unwrap(),
                                                          read_cutoff,
                                                          good_umis.contains(&umi),
                                                          utils::vec_str(umi_contigs.get(umi as usize).unwrap()))).unwrap();
    }

    let ref mut read_iter = reads.iter();
    loop {
        match read_iter.next() {
            None => { break; }
            Some(read) => {
                let (contig_idx, mut alignment) = match all_alignments.get(&read.id) {
                    Some(&(idx, ref al)) => (idx, Some(al.clone())),
                    None => (0, None),
                };
                if !alignment.is_none() {
                    let mut unwrapped_alignment = alignment.unwrap();
                    unwrapped_alignment.ref_idx = contig_idx;
                    alignment = Some(unwrapped_alignment);
                }
                if single_end {
                    let mut rec = bam_utils::read_to_bam_record_opts(&read, &alignment, &None, true, true);
                    let _ = out_bam.write(&mut rec);
                } else {
                    let mate = read_iter.next().expect("Error while reading input reads - not paired?");
                    let (mate_idx, mut mate_alignment) = match all_alignments.get(&mate.id) {
                        Some(&(idx, ref al)) => (idx, Some(al.clone())),
                        None => (0, None),
                    };

                    if !mate_alignment.is_none() {
                        let mut unwrapped_alignment = mate_alignment.unwrap();
                        unwrapped_alignment.ref_idx = mate_idx;
                        mate_alignment = Some(unwrapped_alignment);
                    }
                    let mut rec = bam_utils::read_to_bam_record_opts(&read, &alignment, &mate_alignment, true, true);
                    let mut mate_rec = bam_utils::read_to_bam_record_opts(&mate, &mate_alignment, &alignment, true, true);

                    let _ = out_bam.write(&mut rec);
                    let _ = out_bam.write(&mut mate_rec);
                }
            }
        }
    }
}

pub fn trim<T>(quals: &Vec<T>, min_val: T) -> (usize, usize)
    where T: PartialOrd {
    let mut start_pos = 0;
    let mut end_pos = quals.len() - 1;
    while start_pos < quals.len() && quals[start_pos] <= min_val {
        start_pos += 1;
    }
    while quals[end_pos] <= min_val && end_pos > 0 {
        end_pos -= 1;
    }
    // Leave it to the caller to decide what to do when the whole sequence should be trimmed.
    (start_pos, end_pos)
}
