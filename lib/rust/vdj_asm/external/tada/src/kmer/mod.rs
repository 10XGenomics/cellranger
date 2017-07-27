//
// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
//

//! Data structures for efficient handling of Kmers and related short DNA strings

use std::fmt;
use msp;

mod exts;
mod dir;
pub use self::dir::Dir;
pub use self::exts::Exts;


const KMER_BLOCKS: usize = 2;
const LMER_BLOCKS: usize = 3;
const BLOCK_SIZE: usize = 64;

/// Hard-coded K value for Kmers
pub const K: usize = 47;


/// Store a variable-length DNA sequence, up 92bp in length
#[derive(Copy, Clone, PartialEq, PartialOrd, Eq, Ord, Hash, RustcEncodable, RustcDecodable)]
pub struct Lmer {
    pub storage: [u64; LMER_BLOCKS],
}

pub fn complement(base: u8) -> u8 {
    (!base) & 0x3u8
}


fn block_set(kmer: u64, pos: usize, val: u8) -> u64 {
    let offset = (31 - pos) * 2;
    let mask = !(3 << offset);

    (kmer & mask) | ((val as u64) << offset)
}

fn block_get(kmer: u64, pos: usize) -> u8 {
    let offset = (31 - pos) * 2;
    ((kmer >> offset) & 3) as u8
}

impl Lmer {
    pub fn empty(len: usize) -> Lmer {
        let mut b = Lmer { storage: [0, 0, 0] };

        // Write the length into the last 8 bits
        b.storage[2] = b.storage[2] | ((len as u64) & 0xff);
        b
    }

    pub fn set(&self, pos: usize, val: u8) -> Lmer {
        let block = pos / 32;
        let offset = pos % 32;

        let block_val = block_set(self.storage[block], offset, val);
        let mut new_bases = self.storage;
        new_bases[block] = block_val;
        Lmer { storage: new_bases }
    }

    pub fn get(&self, pos: usize) -> u8 {
        let block = pos / 32;
        let offset = pos % 32;
        block_get(self.storage[block], offset)
    }

    pub fn len(&self) -> usize {
        (self.storage[2] & 0xff) as usize
    }

    // Slice out the kmer for this bsp
    pub fn get_kmer(&self, pos: usize) -> Kmer {
        let mut kmer = Kmer::empty();
        for i in 0..K {
            kmer = kmer.set(i, self.get(pos + i))
        }

        kmer
    }

    pub fn first_kmer(&self) -> Kmer {
        Kmer { storage: [self.storage[0], self.storage[1] & Kmer::lower_mask()] }
    }

    pub fn last_kmer(&self) -> Kmer {
        let s0 = self.storage[0];
        let s1 = self.storage[1];
        let s2 = self.storage[2];

        let l = self.len();

        if l < K {
            panic!("can't do less than K");
        } else if l == K {
            let e0 = s0;
            let e1 = s1;
            Kmer { storage: [e0, e1 & Kmer::lower_mask()] }
        } else if l <= 64 {
            let up_shift = 2 * (l - K);
            let down_shift = (64 - up_shift) % 64;
            let e0 = (s0 << up_shift) | (s1 >> down_shift);
            let e1 = s1 << up_shift;
            Kmer { storage: [e0, e1 & Kmer::lower_mask()] }
        } else if l - K < 32 {
            let up_shift = 2 * (l - K);
            let down_shift = (64 - up_shift) % 64;
            let e0 = (s0 << up_shift) | (s1 >> down_shift);
            let e1 = (s1 << up_shift) | (s2 >> down_shift);
            Kmer { storage: [e0, e1 & Kmer::lower_mask()] }
        } else if l - K - 32 == 0 {
            let e0 = s1;
            let e1 = s2;
            Kmer { storage: [e0, e1 & Kmer::lower_mask()]}
        } else {
            let up_shift = 2 * (l - K - 32);
            let down_shift = (64 - up_shift) % 64;
            let e0 = (s1 << up_shift) | (s2 >> down_shift);
            let e1 = s2 << up_shift;
            Kmer { storage: [e0, e1 & Kmer::lower_mask()]}
        }
    }

    pub fn term_kmer(&self, dir: Dir) -> Kmer {
        match dir {
            Dir::Left => self.first_kmer(),
            Dir::Right => self.last_kmer(),
        }
    }

    pub fn kmers(&self) -> Vec<Kmer> {
        let mut kmers = Vec::new();
        let mut k = self.get_kmer(0);
        kmers.push(k);

        for i in K .. self.len() {
            k = k.extend_right(self.get(i));
            kmers.push(k);
        }

        kmers
    }

    pub fn iter_kmers(&self) -> KmerIter
    {
        KmerIter {
            bases: &self,
            kmer: self.first_kmer(),
            pos: K,
        }
    }


    pub fn new<'a, I: Iterator<Item = &'a u8>>(s: I) -> Lmer {
        let mut b = Lmer { storage: [0, 0, 0] };
        let mut i = 0;
        for c in s {
            b = b.set(i, *c);
            i = i + 1;
        }

        if i > (32 * LMER_BLOCKS) - 4 {
            panic!("Sequence too big for Lmer");
        }

        // Write the length into the last 8 bits
        b.storage[2] = b.storage[2] | ((i as u64) & 0xff);
        b
    }

    pub fn new_from_storage(s: &[u64], k: usize) -> Lmer {
        assert!(s.len() == LMER_BLOCKS);
        let mut lmer = Lmer {storage: [s[0], s[1], s[2]]};
        lmer.storage[2] |= ((k as u64) & 0xff);
        lmer
    }

    fn new_from_slice_slow(s: &[u8]) -> Lmer {
        let mut b = Lmer { storage: [0, 0, 0] };
        let mut i = 0;
        for c in s {
            b = b.set(i, *c);
            i = i + 1;
        }

        // Write the length into the last 8 bits
        b.storage[2] = b.storage[2] | ((i as u64) & 0xff);
        b
    }

    fn new_from_slice(s: &[u8]) -> Lmer {
        let mut b = Lmer { storage: [0, 0, 0] };
        let mut i = 0;

        for blk in 0..3
        {
            let mut chnk : u64 = 0;

            for pos in 0..32
            {
                let offset = (31 - pos) * 2;
                let val = s[i];
                chnk = chnk | ((val as u64) << offset);
                i += 1;

                if i >= s.len(){
                    break;
                }
            }

            b.storage[blk] = chnk;

            if i >= s.len(){
                break;
            }
        }

        // Write length to last 8 bits
        b.storage[2] = b.storage[2] | ((i as u64) & 0xff);
        b
    }


    /// Generate the reverse complement Lmer
    pub fn rc(&self) -> Lmer {
        // not bits to get complement, then reverse order
        let s0 = reverse_by_twos(!self.storage[0]);
        let s1 = reverse_by_twos(!self.storage[1]);
        let s2 = reverse_by_twos(!self.storage[2]);

        let l = self.len();

        if l <= 32 {
            let e0 = s0 << (32 - l) * 2;
            let e1 = 0;
            let e2 = (l as u64) & 0xff;
            Lmer { storage: [e0, e1, e2] }
        } else if l < 64 {
            let up_shift = 2 * (64 - l);
            let down_shift = 64 - up_shift;
            let e0 = (s1 << up_shift) | (s0 >> down_shift);
            let e1 = s0 << up_shift;
            let e2 = (l as u64) & 0xff;
            Lmer { storage: [e0, e1, e2] }
        } else if l == 64 {
            let e0 = s1;
            let e1 = s0;
            let e2 = (l as u64) & 0xff;
            Lmer { storage: [e0, e1, e2] }
        } else {
            let up_shift = 2 * (96 - l);
            let down_shift = 64 - up_shift;
            let e0 = (s2 << up_shift) | (s1 >> down_shift);
            let e1 = (s1 << up_shift) | (s0 >> down_shift);
            let e2 = (s0 << up_shift) | ((l as u64) & 0xff);
            Lmer { storage: [e0, e1, e2] }
        }
    }
}

impl fmt::Debug for Lmer {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut s = String::new();
        for pos in 0..self.len() {
            s.push(bits_to_base(self.get(pos)))
        }
        write!(f, "{}", s)
    }
}

/// Iterate over the Kmers of an Lmer efficiently
pub struct KmerIter<'a>
{
    bases: &'a Lmer,
    kmer: Kmer,
    pos: usize
}

impl<'a> Iterator for KmerIter<'a> {
    type Item=Kmer;

    fn next(&mut self) -> Option<Kmer>
    {
        if self.pos <= self.bases.len()
        {
            let retval = self.kmer;
            self.kmer = self.kmer.extend_right(self.bases.get(self.pos));
            self.pos = self.pos + 1;
            Some(retval)
        }
        else
        {
            None
        }
    }
}




/// Convert an ASCII-encoded DNA base to a 2-bit representation
pub fn base_to_bits(c: u8) -> u8 {
    match c {
        b'A' => 0u8,
        b'C' => 1u8,
        b'G' => 2u8,
        b'T' => 3u8,
        _ => 0u8,
    }
}

/// Convert a 2-bit representation of a base to a char
pub fn bits_to_base(c: u8) -> char {
    match c {
        0u8 => 'A',
        1u8 => 'C',
        2u8 => 'G',
        3u8 => 'T',
        _ => 'X',
    }
}

/// Convert a 2-bit representation of a base to a char
pub fn bits_to_ascii(c: u8) -> u8 {
    match c {
        0u8 => 'A' as u8,
        1u8 => 'C' as u8,
        2u8 => 'G' as u8,
        3u8 => 'T' as u8,
        _ => 'X' as u8,
    }
}



/// Barcoded String Partition. Store an MSP-substring of read, tagged with the partition
/// number and read number within the partition.
#[derive(Debug, Copy, Clone)]
pub struct Bsp {
    /// MSP substring
    pub sequence: Lmer,
    /// Barcode partition
    pub partition: u32,
    /// read number from partition
    pub read: u16,
    /// position of substring in original read
    pub pos: u8,
    /// extensions
    pub extensions: Exts,
}


impl Bsp {
    #[inline(never)]
    pub fn new(read: &[u8], pos: usize, length: usize, partition_num: u32, read_num: u16) -> Bsp {
        let bases = Lmer::new_from_slice(&read[(pos as usize)..((pos + length) as usize)]);

        Bsp {
            sequence: bases,
            partition: partition_num,
            read: read_num,
            pos: pos as u8,
            extensions: Exts::from_slice_bounds(read, pos, length),
        }
    }

    pub fn kmers(&self) -> Vec<Kmer> {
        self.sequence.kmers()
    }


    pub fn iter_kmers(&self) -> KmerIter {
        self.sequence.iter_kmers()
    }


    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    pub fn msp_read(k: usize,
                    p: usize,
                    partition: u32,
                    read: u16,
                    seq: &[u8],
                    permutation: &Vec<usize>,
                    allow_rc: bool)
                    -> Vec<(u32, Bsp)> {
        if seq.len() < K {
            return Vec::new();
        }

        let msp_parts = msp::simple_scan(k, p, seq, permutation, allow_rc);

        let mut v = Vec::new();
        for (shard, _, start_pos, slice_length) in msp_parts {
            let b = Bsp::new(seq, start_pos, slice_length, partition, read);
            v.push((shard, b))
        }

        v
    }

    pub fn extend_kmer(&self, dir: Dir) -> Option<Kmer> {
        match self.extensions.get_unique_extension(dir) {
            Some(ext_base) => {
                let km = match dir {
                    Dir::Left => 0,
                    Dir::Right => self.len() - K,
                };
                Some(self.sequence.get_kmer(km).extend(ext_base, dir))
            }
            None => None,
        }
    }
}

/*
impl fmt::Debug for Bsp {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut s = String::new();
        for pos in 0..self.length {
            s.push(bits_to_base(self.sequence.get(pos as usize)))
        }

        write!(f, "{}", s)
    }
}
*/

#[inline]
pub fn reverse_by_twos(value: u64) -> u64 {
    // swap adjacent pairs
    let mut r = ((value & 0x3333333333333333u64) << 2) | ((value >> 2) & 0x3333333333333333u64);

    // swap nibbles
    r = ((r & 0x0F0F0F0F0F0F0F0Fu64) << 4) | ((r >> 4) & 0x0F0F0F0F0F0F0F0Fu64);

    // swap bytes
    r = ((r & 0x00FF00FF00FF00FFu64) << 8) | ((r >> 8) & 0x00FF00FF00FF00FFu64);

    // swap 2 bytes
    r = ((r & 0x0000FFFF0000FFFFu64) << 16) | ((r >> 16) & 0x0000FFFF0000FFFFu64);

    // swap 4 bytes
    r = ((r & 0x00000000FFFFFFFFu64) << 32) | ((r >> 32) & 0x00000000FFFFFFFFu64);

    r
}


/// A fixed-length Kmer sequence.
#[derive(Copy, Clone, PartialEq, PartialOrd, Eq, Ord, Hash, RustcEncodable, RustcDecodable)]
pub struct Kmer {
    pub storage: [u64; KMER_BLOCKS],
}

impl Kmer {
    pub fn empty() -> Kmer {
        Kmer { storage: [0, 0] }
    }

    pub fn from_bytes(bytes: &Vec<u8>) -> Vec<Kmer> {
        let mut r = Vec::new();

        let mut k0 = Kmer::empty();

        for i in 0..K {
            k0 = k0.set(i, bytes[i])
        }

        r.push(k0);

        for i in 1..(bytes.len() - K + 1) {
            k0 = k0.clone().extend_right(bytes[K + i - 1]);
            r.push(k0);
        }
        r
    }

    pub fn of_string(s: &str) -> Kmer {
        let sbytes: Vec<u8> = s.bytes().map(base_to_bits).collect();
        let mut k0 = Kmer::empty();

        for i in 0..K {
            k0 = k0.set(i, sbytes[i])
        }

        k0
    }

    fn addr(&self, pos: usize) -> (usize, usize) {
        let block = pos / 32;
        let bitpos = (31 - pos % 32) * 2;
        (block, bitpos)
    }

    /// Get the letter at the given position.
    pub fn get(&self, pos: usize) -> u8 {
        let (block, bit) = self.addr(pos);
        (self.storage[block] >> bit & 0x3) as u8
    }

    pub fn set(&self, pos: usize, v: u8) -> Kmer {
        let (block, bit) = self.addr(pos);
        let mask = !(0x3u64 << bit);

        let block_val = self.storage[block];
        let new_block_val = (block_val & mask) | ((v as u64) << bit);

        let mut s = self.storage;
        s[block] = new_block_val;

        Kmer { storage: s }
    }

    /// Return the reverse complement of this kmer
    pub fn rc(&self) -> Kmer {
        // not bits to get complement, then reverse order
        let lower = reverse_by_twos(!self.storage[1]);
        let upper = reverse_by_twos(!self.storage[0]);

        let up_shift = 2 * (64 - K);
        let down_shift = 64 - up_shift;
        let u2 = (lower << up_shift) | (upper >> down_shift);
        let l2 = upper << up_shift;

        Kmer { storage: [u2, l2] }
    }

    pub fn min_rc(&self) -> (Kmer, bool) {
        let kmer = self.clone();
        let rc = self.rc();
        if kmer < rc {
            (kmer, false)
        } else {
            (rc, true)
        }
    }

    pub fn lower_mask() -> u64 {
        if K < 64 {
            let l_bases = K - 32;
            let bits = (1 << 2 * l_bases) - 1;
            bits << (32 - l_bases) * 2
        } else {
            panic!("K >= 64 not implemented");
        }
    }
    /// Shift the base v into the left end of the kmer
    pub fn extend_left(&self, v: u8) -> Kmer {
        let upper = self.storage[0] >> 2 | ((v as u64) << 62);
        let lower = (self.storage[1] >> 2 | (self.storage[0] << 62)) & Kmer::lower_mask();
        Kmer { storage: [upper, lower] }
    }

    pub fn extend_right(&self, v: u8) -> Kmer {
        let upper = self.storage[0] << 2 | self.storage[1] >> 62;
        let lower = self.storage[1] << 2;
        let kmer = Kmer { storage: [upper, lower] };
        kmer.set(K - 1, v)
    }

    pub fn extend(&self, v: u8, dir: Dir) -> Kmer {
        match dir {
            Dir::Left => self.extend_left(v),
            Dir::Right => self.extend_right(v),
        }
    }

    pub fn get_extensions(&self, exts: Exts, dir: Dir) -> Vec<Kmer> {
        let ext_bases = exts.get(dir);
        ext_bases.iter().map(|b| self.extend(b.clone(), dir)).collect()
    }

    // Bucket kmer space into 256 buckets
    pub fn bucket(&self) -> u8
    {
        (self.storage[0] >> (64-8)) as u8
    }

    /// Generate all kmers from string
    pub fn kmers_from_string(str: &[u8]) -> Vec<Kmer> {
        // FIXME - Sofia: I don't think this is correct.
        // k0 is not copied at each iteration
        let mut r = Vec::new();

        if str.len() < K {
            return r;
        }

        let mut k0 = Kmer::empty();

        for i in 0..K {
            k0 = k0.set(i, str[i])
        }

        r.push(k0);

        for i in 1..(str.len() - K + 1) {
            k0 = k0.extend_right(str[K + i - 1]);
            r.push(k0);
        }
        r
    }

    pub fn to_string(&self) -> String {
        let mut s = String::new();
        for pos in 0..K {
            s.push(bits_to_base(self.get(pos)))
        }
        s
    }
}

impl fmt::Debug for Kmer {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut s = String::new();
        for pos in 0..K {
            s.push(bits_to_base(self.get(pos)))
        }

        write!(f, "{}", s)
    }
}


impl fmt::Debug for Pmer {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut s = String::new();
        for pos in 0..self.len() {
            s.push(bits_to_base(self.get(pos)))
        }

        write!(f, "{}", s)
    }
}

/// A variable-length DNA sequence, up to 28bp
#[derive(Copy, Clone, PartialEq, PartialOrd, Eq, Ord, Hash, RustcEncodable, RustcDecodable)]
pub struct Pmer {
    storage: u64,
}

impl Pmer {
    pub fn enumerate(p: usize) -> Vec<Pmer> {
        let mut pmers = Vec::new();

        for i in 0..(1 << 2 * p) {
            let pmer = Pmer { storage: i << (64 - 2 * p) | p as u64 };
            pmers.push(pmer);
        }

        pmers
    }

    pub fn empty(length: usize) -> Pmer {
        if length > 28 {
            panic!("pmer to big");
        }
        //
        Pmer { storage: 0 as u64 | length as u64 }
    }

    pub fn len(&self) -> usize {
        (self.storage & 0xFF) as usize
    }

    pub fn value(&self) -> usize {
        (self.storage >> (64 - 2 * self.len())) as usize
    }

    fn addr(&self, pos: usize) -> usize {
        let bitpos = (31 - pos) * 2;
        bitpos
    }

    /// Get the letter at the given position.
    pub fn get(&self, pos: usize) -> u8 {
        let bit = self.addr(pos);
        (self.storage >> bit & 0x3) as u8
    }

    pub fn set(&self, pos: usize, v: u8) -> Pmer {
        let bit = self.addr(pos);
        let mask = !(0x3u64 << bit);

        let val = self.storage;
        let new_val = (val & mask) | ((v as u64) << bit);

        Pmer { storage: new_val }
    }

    /// Return the reverse complement of this kmer
    pub fn rc(&self) -> Pmer {
        // not bits to get complement, then reverse order
        let len = self.len();
        let rev = reverse_by_twos(!self.storage);

        Pmer { storage: ((rev << (64-len*2)) & self.lmask()) | len as u64 }
    }

    pub fn min_rc(&self) -> (Pmer, bool) {
        let pmer = self.clone();
        let rc = pmer.rc();
        if pmer < rc {
            (pmer, false)
        } else {
            (rc, true)
        }
    }

    fn lmask(&self) -> u64
    {
        let all_ones = !0;
        let top_ones = all_ones << (64 - self.len()*2);
        top_ones
    }

    /// Shift the base v into the left end of the pmer
    pub fn extend_left(&self, v: u8) -> Pmer {
        let v = self.storage >> 2 | ((v as u64) << 62);
        Pmer { storage: (v & self.lmask()) | (self.len() as u64) }
    }

    pub fn extend_right(&self, v: u8) -> Pmer {
        let v = self.storage << 2 | (v as u64) << (64 - self.len() * 2);
        Pmer { storage: (v & self.lmask()) | (self.len() as u64) }
    }

    pub fn extend(&self, v: u8, dir: Dir) -> Pmer {
        match dir {
            Dir::Left => self.extend_left(v),
            Dir::Right => self.extend_right(v),
        }
    }

    pub fn get_extensions(&self, exts: Exts, dir: Dir) -> Vec<Pmer> {
        let ext_bases = exts.get(dir);
        ext_bases.iter().map(|b| self.extend(b.clone(), dir)).collect()
    }


    /// Generate all kmers from string
    pub fn pmers_from_string(str: &[u8], p: usize) -> Vec<Pmer> {

        let mut r = Vec::new();
        let mut k0 = Pmer::empty(p);

        for i in 0..p {
            k0 = k0.set(i, str[i])
        }

        r.push(k0);

        for i in 0..(str.len() - p) {
            k0 = k0.extend_right(str[p + i]);
            r.push(k0);
        }
        r
    }

    pub fn pmers_from_bytes(bytes: Vec<u8>, length: usize) -> Vec<Kmer> {
        let mut r = Vec::new();

        let mut k0 = Kmer::empty();

        for i in 0..length {
            k0 = k0.set(i, bytes[i])
        }

        r.push(k0);

        for i in 1..(bytes.len() - length + 1) {
            k0 = k0.clone().extend_right(bytes[length + i - 1]);
            r.push(k0);
        }
        r
    }


    pub fn of_string(s: &str) -> Pmer {
        let sbytes: Vec<u8> = s.bytes().map(base_to_bits).collect();
        let mut k0 = Pmer::empty(s.len());

        for i in 0..(s.len()) {
            k0 = k0.set(i, sbytes[i])
        }

        k0
    }

    pub fn to_string(&self) -> String {
        let mut s = String::new();
        for pos in 0..self.len() {
            s.push(bits_to_base(self.get(pos)))
        }
        s
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;
    use sim_tests::{random_dna, random_kmer, random_pmer, random_lmer, random_base};
    use rand::{Rng, XorShiftRng};

    fn check_kmer() {
        let km = random_kmer();
        let rc = km.rc();
        assert!(km != rc);

        let double_rc = rc.rc();
        assert!(km == double_rc);

        for i in 0..K {
            assert!(km.get(i) == (3 - rc.get(K - 1 - i)))
        }


        for i in 0..K {
            // Get and set
            let km2 = km.set(i, 0);
            assert!(km2.get(i) == 0);
        }

        let mut copy_kmer = Kmer::empty();
        for i in 0..K {
            copy_kmer = copy_kmer.set(i, km.get(i));
        }
        assert!(km == copy_kmer);


        // Extend right
        let nb = random_base();
        let ext_r = km.extend_right(nb);
        assert!(ext_r.get(K - 1) == nb);
        assert!(km.get(1) == ext_r.get(0));

        // Extend Left
        let nb = random_base();
        let ext_l = km.extend_left(nb);
        assert!(ext_l.get(0) == nb);
        assert!(ext_l.get(1) == km.get(0));

        // Shift twice
        let l_base = random_base();
        let r_base = random_base();
        let ts = km.set(0, l_base).set(K - 1, r_base);

        let double_shift = ts.extend_left(0).extend_right(r_base);
        assert!(ts == double_shift)
    }



    fn check_pmer() {

        let km = random_pmer();
        let p = km.len();

        let rc = km.rc();

        let double_rc = rc.rc();
        assert!(km == double_rc);

        for i in 0..p {
            //
            if km.get(i) != (3 - rc.get(p - 1 - i))
            {
                println!("km: {:?}, rc: {:?}", km, rc);
            }

            assert!(km.get(i) == (3 - rc.get(p - 1 - i)))
        }


        for i in 0..p {
            // Get and set
            let km2 = km.set(i, 0);
            assert!(km2.get(i) == 0);
        }

        let mut copy_kmer = Pmer::empty(p);
        for i in 0..p {
            copy_kmer = copy_kmer.set(i, km.get(i));
        }
        assert!(km == copy_kmer);


        // Extend right
        let nb = random_base();
        let ext_r = km.extend_right(nb);
        assert!(ext_r.get(p - 1) == nb);
        assert!(km.get(1) == ext_r.get(0));

        // Extend Left
        let nb = random_base();
        let ext_l = km.extend_left(nb);
        assert!(ext_l.get(0) == nb);
        assert!(ext_l.get(1) == km.get(0));

        // Shift twice
        let l_base = random_base();
        let r_base = random_base();
        let ts = km.set(0, l_base).set(p - 1, r_base);

        let double_shift = ts.extend_left(0).extend_right(r_base);
        assert!(ts == double_shift)
    }


    fn check_lmer() {

        let km = random_lmer();
        let p = km.len();

        let rc = km.rc();

        let double_rc = rc.rc();
        assert!(km == double_rc);

        for i in 0..p {
            //
            if km.get(i) != (3 - rc.get(p - 1 - i))
            {
                println!("km: {:?}, rc: {:?}", km, rc);
            }

            assert!(km.get(i) == (3 - rc.get(p - 1 - i)))
        }


        for i in 0..p {
            // Get and set
            let km2 = km.set(i, 0);
            assert!(km2.get(i) == 0);
        }

        let mut copy_kmer = Lmer::empty(p);
        for i in 0..p {
            copy_kmer = copy_kmer.set(i, km.get(i));
        }
        assert!(km == copy_kmer);


        let kmers = km.kmers();
        let iter_kmers : Vec<Kmer> = km.iter_kmers().collect();

        assert!(kmers == iter_kmers);

        if km.last_kmer() != kmers[kmers.len() - 1] {

            println!("lmer: {:?}, len: {}, laster: {:?},  iter_last: {:?}", km, km.len(), km.last_kmer(), kmers[kmers.len() - 1]);
        }

        assert!(km.last_kmer() == kmers[kmers.len() - 1]);
    }



    #[test]
    fn test_kmer() {
        for _ in 0..1000 {
            check_kmer();
        }
    }

    #[test]
    fn test_pmer() {
        for _ in 0..1000 {
            check_pmer();
        }
    }

    #[test]
    fn test_lmer() {
        for _ in 0..1000 {
            check_lmer();
        }
    }

    #[test]
    fn test_bsp() {
        let mut rng = XorShiftRng::new_unseeded();

        for _ in 0..100 {
            check_bsp_sequence(&mut rng);
        }
    }

    /// MSP a random string & verify that the BSPs match the original string
    fn check_bsp_sequence(r: &mut XorShiftRng) {

        let len = (r.next_u32() % 255) as usize;
        if len < K {
            return;
        }

        let _dna = random_dna(len as usize);
        let dna = &_dna[..];

        let p = 5;
        let permutation: Vec<usize> = (0..(1 << 2 * p)).collect();

        let bsps = Bsp::msp_read(K, p, 0, 0, dna, &permutation, true);

        println!("seq len: {}", _dna.len());
        for (_, b) in bsps {
            println!("bsp pos: {}, len: {}", b.pos, b.len());
            for i in 0..b.len() {
                assert!(b.sequence.get(i) == _dna[i + (b.pos as usize)]);
            }
        }
    }

    #[test]
    fn check_base_buf() {
        for i in K..((64 * 3 - 8) / 2) {
            for _ in 0..100 {
                let seq = random_dna(i);
                let base_buf = Lmer::new(seq.iter());

                // Check that double RC gives you back the same value
                let bb_rc = base_buf.rc();
                let bb2 = bb_rc.rc();
                assert!(base_buf == bb2);

                // Kmers from raw string match Lmer
                let mut raw_kmers = Kmer::kmers_from_string(&seq);
                raw_kmers.sort();
                let mut bb_kmers = base_buf.kmers();
                bb_kmers.sort();
                assert!(raw_kmers == bb_kmers);

                // RC'd kmers from Lmer match kmers from RC'd Lmer
                let mut rc_kmers: Vec<_> = base_buf.kmers().iter().map(|k| k.rc()).collect();
                rc_kmers.sort();
                let mut bb_rc_kmers = bb_rc.kmers();
                bb_rc_kmers.sort();
                assert!(rc_kmers == bb_rc_kmers);
            }
        }
    }

    /// MSP a bunch of similar strings and make sure each kmers always end up in the same shard
    #[test]
    fn check_consistent_shard() {
        let p = 5;
        let mut rng = XorShiftRng::new_unseeded();
        let len = 250;

        // make a random read and a bunch of edits of it
        let mut dna = random_dna(len);
        let mut dna_seqs = Vec::new();

        for _ in 0..250 {
            for _ in 0..16 {
                let pos = (rng.next_u32() as usize) % len;
                let bp = random_base();
                dna[pos] = bp;
            }

            dna_seqs.push(dna.clone());
        }

        let permutation: Vec<usize> = (0..(1 << 2 * p)).collect();


        let mut kmer_shard = HashMap::new();
        // shard each copy, tag each kmer by the
        // shard id.  fail if we ever get different shards for the same kmer

        for s in dna_seqs {
            for (shard, bsp) in Bsp::msp_read(K, p, 0, 0, &s[..], &permutation, true) {
                for km in bsp.kmers() {
                    let existing_shard = kmer_shard.insert(km, shard);
                    match existing_shard {
                        Some(orig_shard) => {
                            if shard != orig_shard {
                                println!("kmer: {:?}, orig_shard: {}, new_shard: {}",
                                         km,
                                         orig_shard,
                                         shard);
                            }
                            assert!(shard == orig_shard)
                        }

                        None => {}
                    }
                }
            }
        }
    }
}
