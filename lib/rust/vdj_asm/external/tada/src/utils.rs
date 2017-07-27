//! Various useful methods (could be cleaned up)

use std::fs::File;
use std::path::{Path};
use std::collections::HashMap;
use std::io::{BufReader, BufWriter};
use std::io::BufRead;
use kmer::{Kmer, Bsp};

use bincode;
use rustc_serialize::{Encodable, Decodable};
use bincode::rustc_serialize::{EncodingResult, DecodingResult, encode_into, decode_from};

use std::collections::HashSet;
use linked_hash_map::LinkedHashMap;
use std::hash::BuildHasherDefault;
use std::default::Default;
use std::hash::{Hash, Hasher};
use itertools::{Itertools};

pub struct FnvHasher(u64);

impl Default for FnvHasher {
    #[inline]
    fn default() -> FnvHasher {
        FnvHasher(0xcbf29ce484222325)
    }
}

impl Hasher for FnvHasher {
    #[inline]
    fn finish(&self) -> u64 {
        self.0
    }

    #[inline]
    fn write(&mut self, bytes: &[u8]) {
        let FnvHasher(mut hash) = *self;

        for byte in bytes.iter() {
            hash = hash ^ (*byte as u64);
            hash = hash.wrapping_mul(0x100000001b3);
        }

        *self = FnvHasher(hash);
    }
}

pub type DSet<T> = HashSet<T, BuildHasherDefault<FnvHasher>>;

pub type DMap<K, T> = HashMap<K, T, BuildHasherDefault<FnvHasher>>;

pub type LMap<K, T> = LinkedHashMap<K, T, BuildHasherDefault<FnvHasher>>;


pub fn new_dset<T: Hash + Eq>() -> DSet<T> {
    let fnv = BuildHasherDefault::<FnvHasher>::default();
    HashSet::with_hasher(fnv)
}

pub fn new_dmap<K: Hash + Eq, T>() -> DMap<K, T> {
    let fnv = BuildHasherDefault::<FnvHasher>::default();
    HashMap::with_hasher(fnv)
}

pub fn new_lmap<K: Hash + Eq, T>() -> LMap<K, T> {
    let fnv = BuildHasherDefault::<FnvHasher>::default();
    LinkedHashMap::with_hasher(fnv)
}

pub fn write_obj<T: Encodable>(g: &T, filename: &Path) -> EncodingResult<()> {
    let f = match File::create(filename) {
        Err(err) => panic!("couldn't create file {:?}: {}", filename, err),
        Ok(f) => f,
    };
    let mut writer = BufWriter::new(f);
    encode_into(&g, &mut writer, bincode::SizeLimit::Infinite)
}

pub fn read_obj<T: Decodable>(filename: &Path) -> DecodingResult<T> {
    let f = match File::open(filename) {
        Err(err) => panic!("couldn't open file {:?}: {}", filename, err),
        Ok(f) => f,
    };
    let mut reader = BufReader::new(f);
    decode_from(&mut reader, bincode::SizeLimit::Infinite)
}




#[derive(Clone, RustcEncodable, RustcDecodable)]
pub struct MultiVec<T> {
    pub items: Vec<T>,
    pub start_pos: Vec<usize>,
    pub vec_len: Vec<u32>,
}

impl<T : Clone> MultiVec<T> {
    pub fn new() -> MultiVec<T> {
        MultiVec {
            items: Vec::new(),
            start_pos: Vec::new(),
            vec_len: Vec::new(),
        }
    }

    pub fn add<S: Iterator<Item = T>>(&mut self, items: S) {
        let start_pos = self.items.len();
        self.start_pos.push(start_pos);

        let mut n = 0;
        for i in items {
            self.items.push(i);
            n += 1;
        }

        self.vec_len.push(n);
    }

    pub fn add_slice(&mut self, items: &[T]) {
        let start_pos = self.items.len();
        self.start_pos.push(start_pos);

        let mut n = 0;
        for i in items {
            self.items.push(i.clone());
            n += 1;
        }

        self.vec_len.push(n);
    }

    pub fn len(&self) -> usize {
        return self.start_pos.len();
    }

    pub fn get_slice(&self, sub_vec: usize) -> &[T]
    {
        &self.items[(self.start_pos[sub_vec])..(self.start_pos[sub_vec] + self.vec_len[sub_vec] as usize)]
    }
}



pub fn load_barcode_whitelist(filename: &Path) -> HashMap<String, u32> {
    let f = File::open(filename).unwrap();
    let br = BufReader::new(f);

    let mut bc_map = HashMap::new();

    let mut i = 1u32;
    for l in br.lines() {
        bc_map.insert(l.unwrap() + "-1", i);
        i += 1;
    }

    bc_map
}


/// Read a shard and determine the valid kmers
#[inline(never)]
pub fn process_kmer_shard2(bsps: &Vec<Bsp>) -> DMap<Kmer, Vec<(u32, u16)>> {
    let mut kmer_dict = new_dmap();

    for b in bsps {

        //if (b.length as usize) < K || (b.length as usize) != b.sequence.length() || b.pos + b.length > 114
        //{
        //    println!("bad bsp: {:?}", b);
        //    panic!("bad bsp: {:?}", b);
        //}

        let pri = (b.partition, b.read);
        for k in b.kmers() {
            let (min_kmer, _) = k.min_rc();
            let e = kmer_dict.entry(min_kmer).or_insert_with(|| Vec::new());
            (*e).push(pri)
        }
    }

    // filter out kmers we don't like
    let mut final_kmers = new_dmap();
    let mut observed_kmers = 0;

    for (kmer, pris) in kmer_dict {
        observed_kmers += 1;
        let mut multiple_bcs = false;
        let first_bc = pris[0].0;

        for &(bc_id, _) in pris.iter() {
            if bc_id != first_bc {
                multiple_bcs = true;
                break;
            }
        }

        if multiple_bcs {
            final_kmers.insert(kmer, pris);
        }
    }

    println!("Kmers observed: {}. Kmers accepted: {}", observed_kmers, final_kmers.len());
    final_kmers
}

/// Read a shard and determine the valid kmers
#[inline(never)]
pub fn process_kmer_shard(bsps: &Vec<Bsp>) -> (DMap<Kmer, Vec<u32>>, DMap<Kmer, u32>) {

    let mut kmer_buckets : Vec<Vec<(Kmer, u32)>> = Vec::new();
    for _ in 0..256 {
        kmer_buckets.push(Vec::new());
    }
    //let mut kmers = Vec::new();


    info!("Enumerating kmers...");
    for b in bsps {
        for k in b.iter_kmers() {
            let (min_kmer, _) = k.min_rc();
            let bucket = min_kmer.bucket();
            kmer_buckets[bucket as usize].push((min_kmer, b.partition));
        }
    }

    let mut all_kmer_count = new_dmap();
    let mut final_kmers: DMap<Kmer, Vec<u32>> = new_dmap();
    let mut observed_kmers = 0;

    let mut s = 0;
    for bucket in kmer_buckets.iter() {
        s += bucket.len();
    }
    info!("Total size of buckets: {}", s);


    info!("Validating kmers...");
    for mut kmer_vec in kmer_buckets {

        kmer_vec.sort();

        for (kmer, group) in &kmer_vec.into_iter().group_by_lazy(|elt| elt.0) {
            let mut pris = Vec::new();
            observed_kmers = observed_kmers + 1;

            let mut num_bcs = 0;
            let mut last_bc = 0;
            let mut num_obs = 0;
            for (_, bc) in group {
                if bc != last_bc {
                    num_bcs += 1;
                    last_bc = bc;
                    pris.push(bc);
                }
                num_obs += 1
            }

            all_kmer_count.insert(kmer, num_obs);
            if num_bcs > 2 {
                final_kmers.insert(kmer, pris);
            }
        }
    }

    info!("Kmers observed: {}. Kmers accepted: {}", observed_kmers, final_kmers.len());
    (final_kmers, all_kmer_count)
}

pub fn read_fofn(fofn: &Path) -> Vec<String> {
    let f = File::open(fofn).unwrap();
    BufReader::new(f).lines().map(Result::unwrap).collect()
}
