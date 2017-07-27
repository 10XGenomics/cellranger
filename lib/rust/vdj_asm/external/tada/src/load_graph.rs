//
// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
//

use std::io::prelude::*;
use std::cmp;
use std::fs::File;
use std::io::{BufWriter, BufReader};
use std::collections::{HashSet, HashMap};

use bincode;
use bincode::rustc_serialize::{EncodingResult, DecodingResult, encode_into, decode_from};
use rustc_serialize::json;

use flate2::write::ZlibEncoder;
use flate2::read::ZlibDecoder;
use flate2::Compression;

use std::collections::VecDeque;
use std::ascii::AsciiExt;

use bwt;
use bitenc;
use kmer;


pub const BWT_BUCKET_POS: usize = 3;
const OCC_STEP: usize = 64;
pub const SA_STEP: usize = 32;
// constants for computing stats
const MAX_EDGE_LEN: usize = 10000;
const LEN_BUCKET_SIZE: usize = 100;

// possible different layout to save memory.
// start_pos and length index into a shared BitEnd
#[derive(RustcEncodable, RustcDecodable, Clone)]
pub struct Edge {
    pub exts: kmer::Exts,

    // Represents the set of BCs that contributed kmers to this edge
    // Encoded as a index and length into a shared array of BC ids.
    // bcs_start_pos: usize,
    // bcs_length: u16,
    overloaded: bool,

    copy_number_estimate: u32,
}

pub struct Barcode {
    sequence: String,

    edge_start_pos: usize,
    edge_length: u32,
}


#[derive(RustcEncodable, RustcDecodable)]
pub struct Graph {
    // Concatenated edges
    pub all_seq: bitenc::BitEnc,
    edge_starts: Vec<usize>,
    // all_edge_bcs: Vec<u32>,
    // all_bc_edges: Vec<u32>,
    pub edges: Vec<Edge>,
    // bcs: Vec<Barcode>,
    pub first_kmer_idx: HashMap<kmer::Kmer, usize>,
    pub last_kmer_idx: HashMap<kmer::Kmer, usize>,

    fm: bwt::FMIndex,
}

impl Graph {
    pub fn new(all_seq: bitenc::BitEnc,
               edge_starts: Vec<usize>,
               edges: Vec<Edge>,
               first_kmer_idx: HashMap<kmer::Kmer, usize>,
               last_kmer_idx: HashMap<kmer::Kmer, usize>)
               -> Self {

        Graph {
            all_seq: all_seq,
            edge_starts: edge_starts,
            edges: edges,
            first_kmer_idx: first_kmer_idx,
            last_kmer_idx: last_kmer_idx,
            fm: bwt::FMIndex::empty(),
        }
    }

    pub fn from_graph_and_index(g: &Graph, bwt_file: &str) -> Self {
        let b = bwt::read_bwt(bwt_file);
        let fm = bwt::FMIndex::new(b.bwt, b.sa, SA_STEP, OCC_STEP);
        println!("{}", fm.len());
        Graph {
            all_seq: g.all_seq.clone(),
            edge_starts: g.edge_starts.clone(),
            edges: g.edges.clone(),
            first_kmer_idx: g.first_kmer_idx.clone(),
            last_kmer_idx: g.last_kmer_idx.clone(),
            fm: fm,
        }
    }

    pub fn len(&self) -> usize {
        self.all_seq.len()
    }

    pub fn num_edges(&self) -> usize {
        self.edges.len()
    }

    pub fn edge_length(&self, edge_idx: usize) -> usize {
        if edge_idx >= self.edges.len() {
            panic!("Invalid edge index");
        }
        if edge_idx < self.edges.len() - 1 {
            return self.edge_starts[edge_idx + 1] - self.edge_starts[edge_idx];
        }
        self.all_seq.len() - self.edge_starts[edge_idx]
    }

    /// Find the edge index and position in the edge from an index in the global sequence.
    pub fn edge_idx(&self, pos: usize) -> (usize, usize) {
        if pos >= self.all_seq.len() {
            panic!("Invalid position");
        }

        let edge_idx;
        if pos >= self.edge_starts[self.edges.len() - 1] {
            edge_idx = self.edges.len() - 1;
        } else {
            // binary_search returns Ok if the value is in the array
            // otherwise it returns Err(pos) where pos is the position
            // where the value should be inserted.
            edge_idx = match self.edge_starts.binary_search(&pos) {
                Ok(val) => val,
                Err(val) => val - 1,
            };
        }
        let edge_pos = pos - self.edge_starts[edge_idx];
        (edge_idx, edge_pos)
    }

    // Inverse of edge_idx()
    pub fn pos(&self, edge_idx: usize, edge_pos: usize) -> usize {
        self.edge_starts[edge_idx] + edge_pos
    }

    pub fn edge_seq(&self, edge_idx: usize) -> String {
        let start_pos = self.edge_starts[edge_idx];
        let edge_length = self.edge_length(edge_idx);
        let mut edge_seq = "".to_string();
        for i in start_pos..(start_pos + edge_length) {
            edge_seq = edge_seq + &kmer::bits_to_base(self.all_seq.get(i).unwrap()).to_string();
        }
        edge_seq
    }

    fn is_valid_match(&self, pos: usize, seq_len: usize) -> bool {
        let (edge_idx, edge_pos) = self.edge_idx(pos);
        edge_pos + seq_len <= self.edge_length(edge_idx)
    }

    /// Get the last k-mer of the corresponding edge
    fn last_kmer(&self, edge_idx: usize) -> kmer::Kmer {
        let mut kmer = kmer::Kmer::empty();
        let start_idx = self.edge_starts[edge_idx] + self.edge_length(edge_idx) - kmer::K;
        for i in 0..kmer::K {
            kmer = kmer.set(i, self.all_seq.get(i + start_idx).unwrap());
        }
        kmer
    }

    /// Construct a k-mer from all_seq[start_pos:end_pos].
    fn get_kmer(&self, start_pos: usize, end_pos: usize) -> kmer::Kmer {
        let mut kmer = kmer::Kmer::empty();
        for i in start_pos..end_pos {
            kmer = kmer.set(start_pos - i, self.all_seq.get(i).unwrap());
        }
        kmer
    }

    fn val_in_exts(&self, edge_idx: usize, dir: kmer::Dir, val: u8) -> bool {
        self.edges[edge_idx].exts.single_dir(dir).val & (1 << val) > 0
    }

    pub fn contains_seq(&self, bit_seq: &bitenc::BitEnc) -> bool {
        assert!(self.all_seq.len() == self.fm.len(), "Invalid index.");

        let min_len = cmp::min(bit_seq.len(), kmer::K);
        // find a match for the first min_len positions of bit_seq
        let pos: Vec<usize> = self.fm
                                  .backward_search_bitenc(bit_seq.pref(min_len))
                                  .iter()
                                  .map(|x| *x)
                                  .filter(|&x| self.is_valid_match(x, min_len))
                                  .collect();

        if pos.is_empty() {
            return false;
        }
        if min_len < kmer::K {
            return !pos.is_empty();
        }
        println!("Matching pos of prefix {}", pos[0]);
        // there shouldn't be more than one match
        let (mut edge_idx, mut edge_pos) = self.edge_idx(pos[0]);
        let mut first_kmer = self.last_kmer(edge_idx);
        // current position within bit_seq. All letters before that have been matched.
        let mut curr_pos = min_len;

        loop {
            let edge_length = self.edge_length(edge_idx);
            // how many positions on this edge to compare
            let cmp_pos = cmp::min(bit_seq.len() - curr_pos, edge_length - edge_pos - kmer::K);
            println!("pos {} edge {} edge pos {} cmp {}",
                     curr_pos,
                     edge_idx,
                     edge_pos,
                     cmp_pos);
            let cmp = compare_seqs(&self.all_seq,
                                   &bit_seq,
                                   self.pos(edge_idx, edge_pos) + kmer::K,
                                   curr_pos,
                                   cmp_pos);
            if !cmp {
                return false;
            }
            println!("matched some stuff");
            curr_pos += cmp_pos;
            if curr_pos >= bit_seq.len() {
                return true;
            }
            let val = bit_seq.get(curr_pos).unwrap();
            println!("Now need to extend by {}, exts are {}",
                     val,
                     self.edges[edge_idx].exts.val);
            if self.val_in_exts(edge_idx, kmer::Dir::Right, val) {
                first_kmer = first_kmer.extend_right(val);
                println!("Trying to find kmer {}", first_kmer.to_string());
                edge_idx = match self.first_kmer_idx.get(&first_kmer) {
                    Some(idx) => *idx,
                    None => panic!("Extension for edge {} not present in graph.", edge_idx),
                };
                first_kmer = self.last_kmer(edge_idx);
                edge_pos = 0;
                curr_pos += 1;
            } else {
                return false;
            }
        }
    }

    /// Returns a list of edges that would "spell out" the given seq.
    pub fn find_edges(&self, bit_seq: &bitenc::BitEnc) -> Vec<usize> {

        let mut edge_indices = Vec::new();
        let min_len = cmp::min(bit_seq.len(), kmer::K);

        // find a match for the first min_len positions of bit_seq
        let pos: Vec<usize> = self.fm
                                  .backward_search_bitenc(bit_seq.pref(min_len))
                                  .iter()
                                  .map(|x| *x)
                                  .filter(|&x| self.is_valid_match(x, min_len))
                                  .collect();

        if pos.is_empty() {
            return edge_indices;
        }

        if min_len < kmer::K {
            for p in pos.iter().map(|&x| self.edge_idx(x)).map(|(x, _)| x) {
                edge_indices.push(p);
            }
            return edge_indices;
        }

        // there shouldn't be more than one match
        assert_eq!(pos.len(), 1);
        let (mut edge_idx, mut edge_pos) = self.edge_idx(pos[0]);
        edge_indices.push(edge_idx);
        let mut curr_pos = min_len;
        edge_pos += min_len;

        while curr_pos < bit_seq.len() {
            let val = bit_seq.get(curr_pos).unwrap();
            if edge_pos < self.edge_length(edge_idx) {
                if val != self.all_seq.get(self.pos(edge_idx, edge_pos)).unwrap() {
                    return Vec::new();
                }
                edge_pos += 1;
            } else {
                if !self.val_in_exts(edge_idx, kmer::Dir::Right, val) {
                    return Vec::new();
                }
                let mut first_kmer = self.last_kmer(edge_idx);
                first_kmer = first_kmer.extend_right(val);
                edge_idx = match self.first_kmer_idx.get(&first_kmer) {
                    Some(idx) => *idx,
                    None => panic!("Extension for edge {} not present in graph.", edge_idx),
                };
                edge_pos = kmer::K;
                edge_indices.push(edge_idx);
            }
            curr_pos += 1;
        }
        edge_indices
    }

    /// Finds all the paths between two edges.
    pub fn find_paths(&self,
                      start_edge_idx: usize,
                      end_edge_idx: usize,
                      max_hops: usize,
                      good_edges: &HashSet<usize>)
                      -> Vec<GraphPath> {

        // List of all nodes visited during the search
        // Index in the list serves as node id.
        let mut nodes = Vec::<SearchNode>::new();
        // Queue of nodes to visit
        let mut node_queue = VecDeque::<usize>::new();
        // Set of visited graph edges to avoid cycles
        let mut visited = HashSet::<usize>::new();

        let first_node = SearchNode {
            edge_idx: start_edge_idx,
            node_idx: 0,
            parent: None,
            ext: None,
            depth: 0,
        };
        node_queue.push_back(first_node.node_idx);
        visited.insert(first_node.edge_idx);
        nodes.push(first_node);

        let mut last_nodes = Vec::<usize>::new();

        while !node_queue.is_empty() {
            let node_idx = node_queue.pop_front().unwrap();
            let node = nodes[node_idx].clone();
            // Since we're doing breadth first, at this point we've visited all nodes
            // at lower depth.
            if node.depth > max_hops {
                break;
            }
            if node.edge_idx == end_edge_idx {
                last_nodes.push(node_idx);
            }
            let last_kmer = self.last_kmer(node.edge_idx);
            for &val in self.edges[node.edge_idx].exts.get(kmer::Dir::Right).iter() {
                let next_kmer = last_kmer.extend_right(val);
                let edge_idx = match self.first_kmer_idx.get(&next_kmer) {
                    Some(idx) => *idx,
                    None => panic!("Extension for edge {} not present in graph.", node.edge_idx),
                };
                if good_edges.contains(&edge_idx) && !visited.contains(&edge_idx) {
                    let next_node = SearchNode {
                        edge_idx: edge_idx,
                        node_idx: nodes.len(),
                        parent: Some(node.node_idx),
                        ext: Some(val),
                        depth: node.depth + 1,
                    };
                    node_queue.push_back(next_node.node_idx);
                    visited.insert(edge_idx);
                    nodes.push(next_node);
                }
            }
        }

        let mut res = Vec::<GraphPath>::new();
        for node_idx in &last_nodes {
            res.push(traceback(&nodes, *node_idx));
        }
        res
    }

    pub fn compute_bwt_bucket(&self, bucket_val: usize) -> bwt::BWT {
        bwt::compute_bwt_bucket(&self.all_seq, BWT_BUCKET_POS, bucket_val)
    }

    pub fn graph_stats(&self,
                       summary_json: &str,
                       edge_len_tsv: &str,
                       seqs_file: &str,
                       gt_hits_tsv: &str) {
        let num_edges = self.num_edges();
        let seq_len = self.len();
        let mut uniq_seq_len = seq_len;
        let mut num_exts_l: usize = 0;
        let mut num_exts_r: usize = 0;

        let num_len_buckets = MAX_EDGE_LEN / LEN_BUCKET_SIZE;
        let mut len_counts: Vec<usize> = vec![0; num_len_buckets];

        for (eidx, edge) in self.edges.iter().enumerate() {
            num_exts_l += edge.exts.num_exts_l() as usize;
            num_exts_r += edge.exts.num_exts_r() as usize;
            uniq_seq_len -= (edge.exts.num_exts_r() as usize) * (kmer::K - 1);
            let edge_len = self.edge_length(eidx);
            let b = edge_len / LEN_BUCKET_SIZE;
            if b >= num_len_buckets {
                len_counts[num_len_buckets - 1] += 1
            } else {
                len_counts[b] += 1;
            }
        }
        let avg_exts_l: f32 = (num_exts_l as f32) / (num_edges as f32);
        let avg_exts_r: f32 = (num_exts_r as f32) / (num_edges as f32);

        // Write basic stats
        let stats = GraphStats {
            num_edges: num_edges,
            seq_len: seq_len,
            uniq_seq_len: uniq_seq_len,
            avg_exts_l: avg_exts_l,
            avg_exts_r: avg_exts_r,
        };
        stats.write_to_json(summary_json);

        // Write the edge length distribution
        let f = match File::create(edge_len_tsv) {
            Err(err) => panic!("couldn't create file {}: {}", edge_len_tsv, err),
            Ok(f) => f,
        };
        let mut writer = BufWriter::new(f);
        for (cidx, count) in len_counts.iter().enumerate() {
            write!(&mut writer, "{}\t{}\n", (cidx + 1) * LEN_BUCKET_SIZE, count).unwrap();
        }

        // Write the fraction of ground truth sequences found
        let f = match File::open(seqs_file) {
            Err(err) => panic!("couldn't open file {}: {}", gt_hits_tsv, err),
            Ok(f) => f,
        };

        let found_buckets = vec![0, 100, 200, 500, 1000, 5000];
        let mut nseqs = vec![0; found_buckets.len()];
        let mut nseqs_found = vec![0; found_buckets.len()];
        let mut nhops = vec![0; found_buckets.len()];

        let reader = BufReader::new(f);
        for line in reader.lines() {
            let line_str = line.unwrap();
            if line_str.len() < kmer::K {
                continue;
            }
            // binary_search returns Ok if the value is in the array
            // otherwise it returns Err(pos) where pos is the position
            // where the value should be inserted.
            let b = match found_buckets.binary_search(&line_str.len()) {
                Ok(val) => val,
                Err(val) => val - 1,
            };
            let s = bitenc::BitEnc::from_dna_string(&line_str.to_ascii_uppercase());
            let path = self.find_edges(&s);
            nseqs[b] += 1;
            nseqs_found[b] += (path.len() > 0) as usize;
            if path.len() > 0 {
                nhops[b] += path.len();
            }
        }

        let f = match File::create(gt_hits_tsv) {
            Err(err) => panic!("couldn't create file {}: {}", gt_hits_tsv, err),
            Ok(f) => f,
        };
        let mut writer2 = BufWriter::new(f);
        for (bidx, b) in found_buckets.iter().enumerate() {
            write!(&mut writer2,
                   "{}\t{}\t{}\t{}\n",
                   b,
                   nseqs[bidx],
                   nseqs_found[bidx],
                   nhops[bidx])
                .expect("write");
        }
    }
}

#[derive(RustcDecodable, RustcEncodable)]
struct GraphStats {
    num_edges: usize,
    seq_len: usize,
    uniq_seq_len: usize,
    avg_exts_l: f32,
    avg_exts_r: f32,
}

impl GraphStats {
    fn write_to_json(&self, json_file: &str) {
        let encoded = json::encode(&self).unwrap();
        let f = match File::create(json_file) {
            Err(err) => panic!("couldn't create file {}: {}", json_file, err),
            Ok(f) => f,
        };
        let mut writer = BufWriter::new(f);
        writer.write_all(encoded.as_bytes()).expect("write fail");
    }
}

fn traceback(nodes: &Vec<SearchNode>, start_idx: usize) -> GraphPath {
    let mut edges = Vec::new();
    edges.push(nodes[start_idx].edge_idx);
    let mut node_idx = start_idx;

    let mut exts = Vec::new();

    loop {
        match nodes[node_idx].parent {
            Some(idx) => {
                edges.push(nodes[idx].edge_idx);
                exts.push(nodes[node_idx].ext.unwrap());
                node_idx = idx;
            }
            None => break,
        }
    }
    assert_eq!(edges.len() - 1, exts.len());
    edges.reverse();
    exts.reverse();
    GraphPath {
        edges: edges,
        exts: exts,
    }
}

/// Auxiliary structure using for searching a graph.
#[derive(Clone)]
struct SearchNode {
    // Index of a graph edge
    edge_idx: usize,
    // Index of the node in a list of unique nodes
    node_idx: usize,
    // Index of parent node (i.e. how we got here during the search) or None if that's where
    // the search started.
    parent: Option<usize>,
    // Extension that takes us from the parent into this node.
    ext: Option<u8>,
    // # hops to get here since the beginning of the search
    depth: usize,
}

/// A path in the assembly graph.
pub struct GraphPath {
    pub edges: Vec<usize>,
    pub exts: Vec<u8>,
}

impl GraphPath {
    pub fn len(&self) -> usize {
        self.edges.len()
    }

    /// Returns the sequence spelled out by the path.
    pub fn seq(&self, g: &Graph) -> String {
        let mut seq = g.edge_seq(self.edges[0]).to_string();
        if self.edges.len() > 1 {
            for edge_idx in &self.edges[1..] {
                seq = seq + &g.edge_seq(*edge_idx)[(kmer::K - 1)..];
            }
        }
        seq
    }
}

fn compare_seqs(seq1: &bitenc::BitEnc,
                seq2: &bitenc::BitEnc,
                pos1: usize,
                pos2: usize,
                npos: usize)
                -> bool {
    for i in 0..npos {
        if seq1.get(pos1 + i) != seq2.get(pos2 + i) {
            return false;
        }
    }
    true
}


fn read_bytes<R>(reader: &mut R, num_bytes: u64) -> Vec<u8>
    where R: Read
{
    let mut seq = vec![];
    let mut chunk = reader.take(num_bytes as u64);
    let _ = chunk.read_to_end(&mut seq);
    seq
}
/*

/// Loads graph from a binary file.
/// The file is assumed to be in big endian with the following contents:
/// uint8 K
/// uint64 numEdges
/// List of edges
///
/// Each edge is:
/// uint32 length of seq byte array
/// []byte seq byte array
/// byte offset (offset is number of unused bases in seq byte array
/// i.e. if edge is AG, length of seq byte array = 1 and offset = 2)
/// byte context (bitmap of previous - first 4 bits - and next bases - last 4 bits - ordered A, C, G, T)
pub fn load_graph(filename: &str) -> Result<Graph, String> {
    let mut f = match File::open(filename) {
        Err(err) => panic!("couldn't open {}: {}", filename, err),
        Ok(f) => f,
    };

    // Read K
    let k: u8 = try!(read_u8(&mut f));
    assert_eq!(k as usize, kmer::K);

    // Read the number of edges
    let num_edges: u64 = try!(read_u64(&mut f));

    let mut edges = Vec::with_capacity(num_edges as usize);
    let mut seq_bitenc = bitenc::BitEnc::new(2);
    let mut edge_starts = Vec::with_capacity(num_edges as usize);
    let mut first_kmer_idx = HashMap::new();
    let mut last_kmer_idx = HashMap::new();

    for i in 0..num_edges {
        let num_bytes: usize = try!(read_u32(&mut f)) as usize;

        // let seq = read_bytes(&mut reader, num_bytes as u64);
        let mut seq: Vec<u8> = vec![0; num_bytes];
        for j in 0..num_bytes {
            seq[j] = try!(read_u8(&mut f));
        }

        let offset: usize = try!(read_u8(&mut f)) as usize;
        // Can't simply substract a u8 from a u32.
        let edge_length: usize = (num_bytes * 4 - offset) as usize;

        let context: u8 = try!(read_u8(&mut f));
        let exts = kmer::Exts { val: context };

        let mut tmp_bitenc = bitenc::BitEnc::new(2);
        tmp_bitenc.push_bytes(&seq, edge_length);

        first_kmer_idx.insert(tmp_bitenc.first_kmer(), i as usize);
        last_kmer_idx.insert(tmp_bitenc.last_kmer(), i as usize);

        edge_starts.push(seq_bitenc.len());

        let edge = Edge {
            exts: exts,
            overloaded: false,
            copy_number_estimate: 1,
        };

        seq_bitenc.push_bytes(&seq, edge_length);
        edges.push(edge);

        if i % 100000 == 0 && i > 0 {
            println!("Read {} edges", i);
        }
    }
    let g = Graph::new(seq_bitenc,
                       edge_starts,
                       edges,
                       first_kmer_idx,
                       last_kmer_idx);
    println!("Finished making graph: Edges {}, Total len {} bp",
             g.edges.len(),
             g.all_seq.len());
    Ok(g)
}*/

pub fn write_graph(g: &Graph, filename: &str) -> EncodingResult<()> {
    let f = match File::create(filename) {
        Err(err) => panic!("couldn't create file {}: {}", filename, err),
        Ok(f) => f,
    };
    let writer = BufWriter::new(f);
    let mut encoder = ZlibEncoder::new(writer, Compression::Fast);
    encode_into(&g, &mut encoder, bincode::SizeLimit::Infinite)
}

pub fn read_graph(filename: &str) -> DecodingResult<Graph> {
    let f = match File::open(filename) {
        Err(err) => panic!("couldn't open file {}: {}", filename, err),
        Ok(f) => f,
    };
    let reader = BufReader::new(f);
    let mut decoder = ZlibDecoder::new(reader);
    decode_from(&mut decoder, bincode::SizeLimit::Infinite)
}

#[cfg(test)]
mod tests {
    use super::*;
    use kmer;
    use bitenc;
    use std::collections::{HashSet, HashMap};
    use bwt;

    /*
    #[test]
    /// Text representation of graph:
    /// K = 48
    /// Number of edges = 4
    /// Edges:
    /// AAATTAAAGAAGTTTTCCCAACCTGGTTTTAAAAGCTGCATATTTGAAATTAAAAATATTTATTTATTTATTGTGTTAGTCCATTTTCACGCTGCTGATAAAAATATACCCAAGACTGGGAAGAAAA
    /// TTTTCTTCCCAGTCTTGGGTATATTTTTATCAGCAGCGTGAAAATGGACTAACACAATAAATAAATAAATATTTTTAATTTCAAATATGCAGCTTTTAAAACCAGGTTGGGAAAACTTCTTTAATTT
    /// GCTTCTTCCTTATTTTCTCTTGCTGCCACCCTGTAAGAAGTGCCTTTCACCTCCCACCATAATTCTGAGCCCTCCCCAGCCATATGGAACTGTAAGTCCAATTAAACCTCTTTTTCTTCCCAGTTTTGGGGATATGTTTATCAGCAGCGT
    /// ACGCTGCTGATAAACATATCCCCAAAACTGGGAAGAAAAAGAGGTTTAATTGGACTTACAGTTCCATATGGCTGGGGAGGGCTCAGAATTATGGTGGGAGGTGAAAGGCACTTCTTACAGGGTGGCAGCAAGAGAAAATAAGGAAGAAGC
    fn test_load_small_graph() {
        let filename = "test_data/test_small_graph";
        let mut g = match load_graph(filename) {
            Err(err) => panic!(err),
            Ok(g) => g,
        };

        // let encoded: Vec<u8> = encode(&g, SizeLimit::Infinite).unwrap();
        // let decoded: Graph2 = decode(&encoded[..]).unwrap();
        // assert!(g.first_kmer_idx == decoded.first_kmer_idx);

        assert_eq!(g.num_edges(), 4);
        assert_eq!(g.edge_length(0), 127);
        assert_eq!(g.edge_length(1), 127);
        assert_eq!(g.edge_length(2), 150);
        assert_eq!(g.edge_length(3), 150);

        assert_eq!(g.edge_idx(0), (0, 0));
        assert_eq!(g.edge_idx(120), (0, 120));
        assert_eq!(g.edge_idx(g.len() - 1), (3, 149));

        for edge in &g.edges {
            assert_eq!(edge.exts.val, 0);
        }

        assert_eq!(g.first_kmer_idx.len(), 4);
        assert_eq!(g.last_kmer_idx.len(), 4);

        let (sa, bwt) = bwt::compute_bwt(&g.all_seq, 1, 2);
        let fm = bwt::FMIndex::new(bwt, sa, 2, 2);
        g.fm = fm;

        // Search sequences smaller than K
        assert!(g.contains_seq(&bitenc::BitEnc::from_dna_string("AAATTAAAGAAGTTTTCC")));
        assert!(g.contains_seq(&bitenc::BitEnc::from_dna_string("CAGTTCCATATGGCT")));

        // Search for something with multiple occurrences
        assert!(g.contains_seq(&bitenc::BitEnc::from_dna_string("TTT")));
        let seq1 = "AAATTAAAGAAGTTTTCCCAACCTGGTTTTAAAAGCTGCATATTTGAAATTAAAAATATTTAT".to_string() +
                   "TTATTTATTGTGTTAGTCCATTTTCACGCTGCTGATAAAAATATACCCAAGACTGGGAAGAAAA";
        assert!(g.contains_seq(&bitenc::BitEnc::from_dna_string(&seq1)));
    }
    */

    //#[test]
    fn test_toy_graph() {
        // Edges:
        // Edge1: {A}_10 ACGT {A}_43
        // Edge2: ACGT {A}_43 AGGT
        // Edge3: {A}_43 AGGTTTAT
        // Edge4: ACGT {A}_43 TAAG
        // seq1 -> seq2 with ext A
        // seq1 -> seq4 with ext T
        // seq2 -> seq3 with ext T
        // Exts (remember right is higher order bits):
        // seq1 : 10010000 128+16
        // seq2 : 10000001
        // seq3 : 00001000
        // seq4 : 00000001
        let str0 = String::from_utf8(vec![b'A'; 43]).unwrap();
        let seq1 = String::from_utf8(vec![b'A'; 10]).unwrap() + "ACGT" + &str0;
        let seq2 = "ACGT".to_string() + &str0 + "AGGT";
        let seq3 = str0.clone() + "AGGTTTAT";
        let seq4 = "ACGT".to_string() + &str0 + "TAAG";
        let seqs = vec![seq1, seq2, seq3, seq4];

        let exts1 = kmer::Exts { val: 144u8 };
        let exts2 = kmer::Exts { val: 129u8 };
        let exts3 = kmer::Exts { val: 8u8 };
        let exts4 = kmer::Exts { val: 1u8 };
        let exts = vec![exts1, exts2, exts3, exts4];

        let num_edges = 4;
        let mut edges = Vec::with_capacity(num_edges);
        let mut seq_bitenc = bitenc::BitEnc::new();
        let mut edge_starts = Vec::with_capacity(num_edges);
        let mut first_kmer_idx = HashMap::new();
        let mut last_kmer_idx = HashMap::new();

        for i in 0..num_edges {
            let tmp_bitenc = bitenc::BitEnc::from_dna_string(&seqs[i]);

            first_kmer_idx.insert(tmp_bitenc.first_kmer(), i as usize);
            last_kmer_idx.insert(tmp_bitenc.last_kmer(), i as usize);

            edge_starts.push(seq_bitenc.len());

            let edge = Edge {
                exts: exts[i],
                overloaded: false,
                copy_number_estimate: 1,
            };

            for val in tmp_bitenc.iter() {
                seq_bitenc.push(val);
            }
            edges.push(edge);
        }

        let mut g = Graph::new(seq_bitenc,
                               edge_starts,
                               edges,
                               first_kmer_idx,
                               last_kmer_idx);
        let bwt = bwt::compute_bwt(&g.all_seq, 1, 1);
        let fm = bwt::FMIndex::new(bwt.bwt, bwt.sa, 1, 1);
        g.fm = fm;

        for i in 0..num_edges {
            assert_eq!(g.edge_seq(i), seqs[i]);
        }

        // Searches
        assert!(g.contains_seq(&bitenc::BitEnc::from_dna_string("AAA")));
        let e: HashSet<usize> = g.find_edges(&bitenc::BitEnc::from_dna_string("AAA"))
                                 .iter()
                                 .cloned()
                                 .collect();
        for i in 0..num_edges {
            assert!(e.contains(&i));
        }

        for i in 0..num_edges {
            assert!(g.contains_seq(&bitenc::BitEnc::from_dna_string(&seqs[i])));
            let e = g.find_edges(&bitenc::BitEnc::from_dna_string(&seqs[i]));
            assert!(e.contains(&i));
            assert_eq!(e.len(), 1);
        }
        // seq in the middle of an edge
        let s = "T".to_string() + &str0 + "A";
        assert!(g.contains_seq(&bitenc::BitEnc::from_dna_string(&s)));
        let e = g.find_edges(&bitenc::BitEnc::from_dna_string(&s));
        assert_eq!(e.into_iter().collect::<Vec<usize>>(), [1]);

        // seq across edges
        let s = "AAACGT".to_string() + &str0 + "AG";
        assert!(g.contains_seq(&bitenc::BitEnc::from_dna_string(&s)));
        let e = g.find_edges(&bitenc::BitEnc::from_dna_string(&s));
        assert!(e.contains(&0));
        assert!(e.contains(&1));
        assert_eq!(e.len(), 2);

        // seq that doesn't really exist in graph but exists in concatenation of edges
        let s = "AGGTA";
        assert!(!g.contains_seq(&bitenc::BitEnc::from_dna_string(&s)));
        let e = g.find_edges(&bitenc::BitEnc::from_dna_string(&s));
        assert_eq!(e.len(), 0);

        // Paths
        let good_edges: HashSet<usize> = vec![0, 1, 2, 3].iter().cloned().collect();
        let paths = g.find_paths(0, 2, 10, &good_edges);
        let path_seq = seqs[0].clone() + "AGGTTTAT";
        assert_eq!(paths[0].seq(&g), path_seq);
        assert_eq!(paths[0].edges, vec![0, 1, 2]);
        assert_eq!(paths[0].exts, vec![0, 3]);
        let paths = g.find_paths(0, 3, 10, &good_edges);
        assert_eq!(paths[0].edges, vec![0, 3]);
        assert_eq!(paths[0].exts, vec![3]);
        assert_eq!(g.find_paths(0, 2, 0, &good_edges).len(), 0);
        assert_eq!(g.find_paths(0, 2, 1, &good_edges).len(), 0);
        let paths = g.find_paths(0, 0, 0, &good_edges);
        assert_eq!(paths[0].edges, vec![0]);
        assert_eq!(paths[0].exts.len(), 0);
    }
}
