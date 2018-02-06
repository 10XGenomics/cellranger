//
// Copyright (c) 2017 10x Genomics, Inc. All rights reserved.
//

use std::io::{Read};
use std::hash::Hash;
use std::fmt::Debug;

use fnv::FnvHashSet;
use bio::io::{fasta};
use debruijn::{Kmer, Vmer};
use debruijn::dna_string::DnaString;
use itertools::Itertools;
use pdqsort;
use bincode;
use serde::de::DeserializeOwned;

use boomphf::Mphf;

#[derive(Serialize, Deserialize)]
pub enum KmerIndexType {
    HashSetIndex,
    BBHashIndex,
}

pub trait KmerPresenceQuery<K> {
    fn contains(&self, &K) -> bool;
}

#[derive(Serialize, Deserialize)]
pub struct HashSetKmerIndex<K: Hash + Eq> {
    index_type: KmerIndexType,
    pub kmers: FnvHashSet<K>,
}
impl<K: Hash + Eq> HashSetKmerIndex<K> {
    pub fn new(capacity: usize) -> HashSetKmerIndex<K> {
        HashSetKmerIndex {
            index_type: KmerIndexType::HashSetIndex,
            kmers: FnvHashSet::<K>::with_capacity_and_hasher(capacity, Default::default()),
        }
    }
}

impl<K: Hash + Eq> KmerPresenceQuery<K> for HashSetKmerIndex<K> {
    fn contains(&self, query: &K) -> bool{
        self.kmers.contains(query)
    }
}

#[derive(Serialize, Deserialize)]
pub struct BBHashKmerIndex<K: Hash + Eq + Clone + Debug> {
    index_type: KmerIndexType,
    pub kmers: Vec<K>,
    pub mphf: Mphf<K>,
}
impl<K: Hash + Eq + Clone + Debug> BBHashKmerIndex<K> {
    pub fn new(kmers: Vec<K>, mphf: Mphf<K>) -> BBHashKmerIndex<K> {
        BBHashKmerIndex {
            index_type: KmerIndexType::BBHashIndex,
            kmers: kmers,
            mphf: mphf,
        }
    }
}

impl<K: Hash + Eq + Debug + Clone> KmerPresenceQuery<K> for BBHashKmerIndex<K> {
    fn contains(&self, query: &K) -> bool{
        match self.mphf.try_hash(query) {
            Some(i) => self.kmers[i as usize] == *query,
            None => false,
        }
    }
}


pub fn index_transcripts_hashset<R: Read, K: Kmer + Hash>(fa_reader: fasta::Reader<R>,
                                                          skip_bases: usize) -> HashSetKmerIndex<K> {
    let mut idx = HashSetKmerIndex::<K>::new(1000000);

    info!("Building hashset");
    for rec in fa_reader.records()
        .map(|rec| rec.expect("Failed to parse FASTA record"))
        .filter(|rec| rec.seq().len() >= K::k()) {

            let dna_string = DnaString::from_acgt_bytes(&rec.seq());

            for kmer in dna_string.iter_kmers().step(1+skip_bases) {
                idx.kmers.insert(kmer);
            }
        }
    idx
}

pub fn index_transcripts_mphf<R: Read, K: Kmer + Hash>(fa_reader: fasta::Reader<R>,
                                                       est_n_kmers: usize,
                                                       skip_bases: usize,
                                                       gamma: f64) -> BBHashKmerIndex<K> {
    let mut kmers = Vec::<K>::with_capacity(est_n_kmers);

    for rec in fa_reader.records()
        .map(|rec| rec.expect("Failed to parse FASTA record"))
        .filter(|rec| rec.seq().len() >= K::k()) {

            let dna_string = DnaString::from_acgt_bytes(&rec.seq());

            for kmer in dna_string.iter_kmers().step(1+skip_bases) {
                kmers.push(kmer);
            }
        }

    info!("Sorting {} kmers", kmers.len());
    pdqsort::sort(&mut kmers);

    info!("Deduping {} kmers", kmers.len());
    kmers.dedup();

    info!("Hashing {} unique kmers", kmers.len());
    let mphf = Mphf::new(gamma, &kmers, None);

    info!("Arranging {} unique kmers", kmers.len());
    // Permuting in-place murders the cache, so copy
    let mut arranged_kmers = vec![K::empty(); kmers.len()];
    for kmer in kmers {
        arranged_kmers[mphf.hash(&kmer) as usize] = kmer;
    }

    BBHashKmerIndex::new(arranged_kmers, mphf)
}

// Get the type of kmer index stored in a file
pub fn get_index_type<R: Read>(mut reader: &mut R) -> KmerIndexType {
     bincode::deserialize_from(&mut reader, bincode::Infinite)
        .expect("Failed to read index type")
}

pub fn load_index<I: DeserializeOwned, R: Read>(mut reader: &mut R) -> I {
    let index: I = bincode::deserialize_from(&mut reader, bincode::Infinite)
        .expect("Failed to deserialize index");
    index
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_bbhash() {
        use index;
        use index::KmerPresenceQuery;
        use std::io::Cursor;
        use bio::io::fasta;
        use tests::random_seq_rng;
        use debruijn::Kmer;
        use debruijn::kmer::{IntKmer};
        use debruijn::{Vmer};
        use debruijn::dna_string::DnaString;
        use rand::{SeedableRng, StdRng};

        type MyKmer = IntKmer<u64>;

        let seed: &[_] = &[1, 2, 3, 4];
        let mut rng: StdRng = SeedableRng::from_seed(seed);

        let seq = random_seq_rng(100, &mut rng);
        let fa_str = format!(">1\n{}", String::from_utf8(seq.clone()).unwrap());

        let reader = fasta::Reader::new(Cursor::new(fa_str.as_bytes()));

        let index: index::BBHashKmerIndex<MyKmer> = index::index_transcripts_mphf(reader, 1000, 0, 2.0);
        //let index: index::HashSetKmerIndex<MyKmer> = index::index_transcripts_hashset(reader, 0);

        let dna_piece = DnaString::from_acgt_bytes(&seq[0..32].to_owned());
        let kmer: MyKmer = dna_piece.iter_kmers().nth(0).unwrap();

        assert!(index.contains(&kmer));
        assert!(!index.contains(&MyKmer::empty()));
    }
}
