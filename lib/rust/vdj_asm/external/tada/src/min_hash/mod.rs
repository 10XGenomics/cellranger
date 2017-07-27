//
// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
//

use std::hash::SipHasher;
use std::hash::Hash;
use std::marker::PhantomData;
use std::hash::Hasher;
use std::cmp::min;
use std::u64::MAX;
use std::collections::HashMap;

struct MinHasher<T:Hash>
{
    num_hashes: usize,
    hash_keys: Vec<(u64, u64)>,
    phantom: PhantomData<T>,
}

impl<T:Hash> MinHasher<T>
{
    pub fn new(n: usize) -> MinHasher<T>
    {
        let funcs = (0u64..(n as u64)).map(|i| (i,i+1)).collect();

        MinHasher
        {
            num_hashes: n,
            hash_keys: funcs,
            phantom: PhantomData,
        }
    }

    pub fn min_hash_signature(&self, values: &Vec<T>) -> Vec<u64>
    {
        let mut min_hash_vec = Vec::new();
        for &(k1, k2) in self.hash_keys.iter()
        {
            let ref mut hasher = SipHasher::new_with_keys(k1, k2);

            let mut min_hash : u64 = MAX;
            for v in values.iter()
            {
                v.hash(hasher);
                let hash_value = hasher.finish();
                min_hash = min(min_hash, hash_value);
            }

            min_hash_vec.push(min_hash);
        }

        min_hash_vec
    }
}

struct MinHashIndex<K, T:Hash>
{
    hasher: MinHasher<T>,
    items: Vec<K>,
    hash_dicts: Vec<HashMap<u64, Vec<usize>>>
}

impl<K, T:Hash> MinHashIndex<K, T>
{
    pub fn new(num_hashes: usize) -> MinHashIndex<K,T>
    {
        let hash_dicts = (0..num_hashes).map(|_| { HashMap::new() }).collect();

        MinHashIndex
        {
            hasher: MinHasher::new(num_hashes),
            items: Vec::new(),
            hash_dicts: hash_dicts
        }
    }

    pub fn insert(&mut self, key: K, set: &Vec<T>)
    {
        // Insert item
        let id = self.items.len();
        self.items.push(key);

        // Min-hash signature
        let sig = self.hasher.min_hash_signature(&set);

        for (mut h, sig_item) in self.hash_dicts.iter_mut().zip(sig)
        {
            let mut item_vec = h.entry(sig_item).or_insert_with(||{Vec::new()});
            item_vec.push(id);
        }
    }

    /// Query the min-hash index for items similar to the qset signature
    /// Returns items from the index with expected Jaccard similarity greater
    /// than min_similarity.
    pub fn query(&self, min_similarity: f64, query_set: &Vec<T>) -> Vec<(f64, &K)>
    {
        let mut hit_counts = HashMap::new();

        // Compute the signature of the query set
        let sig = self.hasher.min_hash_signature(query_set);

        for (h, sig_item) in self.hash_dicts.iter().zip(sig)
        {
            for i in h.get(&sig_item).unwrap_or(&Vec::new())
            {
                let counts = hit_counts.entry(*i).or_insert(0);
                *counts = *counts + 1;
            }
        }

        let mut res = Vec::new();
        for (item_id,hits) in hit_counts.iter()
        {
            let jacc = (*hits as f64) / (self.hasher.num_hashes as f64);
            if jacc > min_similarity
            {
                res.push((jacc, self.items.get(*item_id).unwrap()));
            }
        }

        res
    }
}
