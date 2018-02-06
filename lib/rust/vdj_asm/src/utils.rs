//
// Copyright (c) 2017 10x Genomics, Inc. All rights reserved.
//

// Allow deprecated to suppress SipHasher warnings
#![allow(deprecated)]

use std::cmp::{Ord, Ordering};
use rand;
use std::hash::{Hash, SipHasher, Hasher};
use std::fmt::Display;
use std::io::Read;
use std::fs::File;
use std::path::{Path};
use constants::NUCS;
use lz4;

#[derive(PartialEq, PartialOrd, Clone, Debug)]
pub struct NonNan(f32);

impl NonNan {
    pub fn new(val: f32) -> Option<NonNan> {
        if val.is_nan() {
            None
        } else {
            Some(NonNan(val))
        }
    }
}

impl Eq for NonNan {}

impl Ord for NonNan {
    fn cmp(&self, other: &NonNan) -> Ordering {
        self.partial_cmp(other).unwrap()
    }
}

/// Computes log(exp(a) + exp(b)) in a way that tries to avoid over/underflows.
/// log(exp(a - b) + 1) + b = log((exp(a) + exp(b)) - exp(b)) + b = log(exp(a) + exp(b))
pub fn logaddexp(a: f64, b: f64, base: f64) -> f64 {
    let c: f64;
    let d: f64;

    // powf(a-b) doesn't do the right thing if a > b.
    if a > b {
        c = b;
        d = a;
    } else {
        c = a;
        d = b;
    }

    (base.powf(c - d) + 1.0).log(base) + d
}

// Computes
//
//     /    n        \
//     |  =====      |
//     |  \      x_i |
// log |   >    k    |
//    k|  /          |
//     |  =====      |
//     \  i = 1      /
//
// where k is the base and x_i are the elements in the vector vals
// using the function logaddexp() successively.
pub fn logaddexp_arr(vals: &Vec<f64>, base: f64) -> f64 {
    assert!(vals.len() > 1);
    let mut res = logaddexp(vals[0], vals[1], base);

    for &v in vals.iter().skip(2) {
        res = logaddexp(res, v, base);
    }
    res
}

// Computes ln(exp(a) + exp(b)) using a piecewise linear approximation
// See https://web.stanford.edu/~boyd/papers/rgp.html for details
//
// The coefficients are precomputed and listed in constants::LSE_COEFF
pub fn lnaddexp_approx(a: f64, b: f64) -> f64 {
    use constants::LSE_COEFF;
    let c: f64;
    let d: f64;
    if a > b {
        c = a;
        d = b;
    } else {
        c = b;
        d = a;
    }

    let order = LSE_COEFF.len();
    let mut result = c;
    for i in 0..order {
        let tmp = LSE_COEFF[i][0]*d + LSE_COEFF[i][1]*c + LSE_COEFF[i][2];
        if result < tmp {
            result = tmp;
        }
    }
    result

}

// Computes approximation of
//   /    n           \
//   |  =====         |
//   |  \             |
// ln|   >    exp(x_i)|
//   |  /             |
//   |  =====         |
//   \  i = 1         /
//
// x_i are the elements of the vector vals.
// The sum is computed by pairwise reduction using the fn lnaddexp_approx.
// As an exmaple the first reduction step in the case of 6 element vector would be
// ln(exp(x1) + exp(x2) + exp(x3) + exp(x4) + exp(x5) + exp(x6)) =
//       ln( exp(ln(exp(x1) + exp(x2))) + exp(ln(exp(x3) + exp(x4))) + exp(ln(exp(x5) + exp(x6))))
// and so on.
// Successive reduction as in the fn logaddexp_arr() would accumulate errors much faster
// compared to pairwise reduction,
pub fn lnaddexp_approx_arr(vals: &Vec<f64>) -> f64 {
    assert!(vals.len() > 1);
    let mut res = vals.clone();
    let mut step = 2;
    loop {
        let mut i = 0;
        while i < (vals.len() - step/2) {
            let tmp = lnaddexp_approx(res[i], res[i+step/2]);
            res[i] = tmp;
            i+=step;
        }
        if step >= vals.len() {
            break;
        }
        step *= 2;
    }
    return res[0];
}

// Computes approximation of
//
//     /    n        \
//     |  =====      |
//     |  \      x_i |
// log |   >    k    |
//    k|  /          |
//     |  =====      |
//     \  i = 1      /
//
// where k is the base and x_i are the elements in the vector vals
// using the lnaddexp_approx_arr() fn
pub fn logaddexp_approx_arr(vals: &Vec<f64>, base: f64) -> f64 {
    let lnbase = base.ln();
    let new_vals = vals.iter().map(|&x| x * lnbase).collect();
    lnaddexp_approx_arr(&new_vals)/lnbase
}

pub fn argmax<T: PartialOrd + Clone + Copy>(scores: &Vec<T>) -> Option<(usize, T)> {

    if scores.is_empty() { return None; }

    let mut max_index = 0;
    let mut max = scores.get(max_index).unwrap();

    for (i, v) in scores.iter().enumerate() {
        if max < v {
            max_index = i;
            max = v;
        }
    }

    Some((max_index, *max))
}

pub fn nx_count(nums: &Vec<usize>, x: f64) -> usize {

    if nums.is_empty() {
        return 0;
    }
    let frac = (x as f64) / 100.0;

    let mut nums_cp = nums.clone();
    nums_cp.sort();
    nums_cp.reverse();

    let mut sum = 0.0 as f64;
    for &n in nums_cp.iter() {
        sum = sum + (n as f64);
    }

    let mut curr_sum = 0.0 as f64;
    for &n in nums_cp.iter() {
        curr_sum = curr_sum + (n as f64);
        if (curr_sum / sum) >= frac {
            return n;
        }
    }
    nums_cp[nums_cp.len() - 1]
}

pub fn vec_str<T>(nums: &Vec<T>) -> String
    where T: Display {
    let mut out = "".to_string();
    for (idx, num) in nums.iter().enumerate() {
        if idx == 0 {
            out = format!("{}", num);
        } else {
            out = format!("{},{}", out, num);
        }
    }
    out.to_string()
}

pub fn rm_extension(path_name: &String) -> String {
    let path = Path::new(&path_name);
    let dir = path.parent().unwrap(); // dirname
    let name_no_ext = path.file_stem().unwrap(); // filename without the  extension
    dir.join(name_no_ext).into_os_string().into_string().unwrap()
}

pub fn replace_ns_rand(s: &String, rng: &mut rand::StdRng) -> String {
    let in_bytes = s.as_bytes();
    let mut out_bytes : Vec<u8> = Vec::with_capacity(in_bytes.len());
    for c in in_bytes.iter() {
        out_bytes.push(match *c {
            b'A' | b'C' | b'G' | b'T' => *c,
            _ => NUCS[rand::sample(rng, 0..4, 1)[0]],
        });
    }
    String::from_utf8(out_bytes).unwrap()
}

pub fn replace_ns(s: &String, name: &str) -> String {
    let in_bytes = s.as_bytes();
    let mut out_bytes : Vec<u8> = Vec::with_capacity(in_bytes.len());
    for (i, c) in in_bytes.iter().enumerate() {
        out_bytes.push(match *c {
            b'A' | b'C' | b'G' | b'T' => *c,
            _ => {
                let mut hasher = SipHasher::new();
                name.hash(&mut hasher);
                i.hash(&mut hasher);
                NUCS[(hasher.finish() % 4) as usize]
            },
        });
    }
    String::from_utf8(out_bytes).unwrap()
}

/// Convert an ASCII-encoded DNA base to a 2-bit representation
pub fn base_to_bits_rand(c: u8, rng: &mut rand::StdRng) -> u8 {
    match c {
        b'A' => 0u8,
        b'C' => 1u8,
        b'G' => 2u8,
        b'T' => 3u8,
        _ => rand::sample(rng, 0..4, 1)[0],
    }
}

pub fn base_to_bits_hash(c: u8, name: &str, pos: usize) -> u8 {

    match c {
        b'A' => 0u8,
        b'C' => 1u8,
        b'G' => 2u8,
        b'T' => 3u8,
        _ => {
                let mut hasher = SipHasher::new();
                name.hash(&mut hasher);
                pos.hash(&mut hasher);
                (hasher.finish() % 4) as u8
            },
    }
}

pub fn base_to_bits_safe(c: u8) -> u8 {
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
        _ => 'N',
    }
}

/// Convert a 2-bit representation of a base to a char
pub fn bits_to_ascii(c: u8) -> u8 {
    match c {
        0u8 => 'A' as u8,
        1u8 => 'C' as u8,
        2u8 => 'G' as u8,
        3u8 => 'T' as u8,
        _ => 'N' as u8,
    }
}

/// Pseudorandomly drop a read or read pair. Returns true to drop, false to keep.
pub fn drop_read(subsample_rate: f64, name: &str) -> bool {
    match subsample_rate {
        x if x == 0.0 => true,
        y if y == 1.0 => false,
        _ => {
            use rand::Rng;
            let mut hasher = SipHasher::new();

            // Strip the augmented data from this read name if there are any
            let key = match name.find("|||") {
                Some(i) => &name[0..i],
                None => name,
            };
            key.hash(&mut hasher);

            let hash = hasher.finish();
            let seed = [(hash >> 32) as u32, hash as u32, (hash >> 32) as u32, hash as u32];
            let mut rng: rand::XorShiftRng = rand::SeedableRng::from_seed(seed);
            rng.gen::<f64>() > subsample_rate
        }
    }
}

pub fn find_file_maybe_compressed(prefix: &str) -> Option<String> {
    if Path::new(prefix).exists() {
        return Some(prefix.to_owned())
    } else if (Path::new(&(prefix.clone().to_owned() + ".lz4"))).exists() {
        return Some(prefix.to_owned() + ".lz4")
    }
    None

}

pub fn open_lz4(filename: &str) -> lz4::Decoder<File> {
    let f = File::open(filename).expect("Failed to open file for reading");
    lz4::Decoder::new(f).expect("Failed to create lz4 decoder")
}

pub fn open_maybe_compressed(filename: &str) -> Box<Read> {
    match filename.ends_with(".lz4") {
        true => Box::new(open_lz4(filename)) as Box<Read>,
        false => Box::new(File::open(filename)
                          .expect("Failed to open file for reading")) as Box<Read>,
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use rand::{StdRng, SeedableRng};

    macro_rules! assert_delta {
        ($x:expr, $y:expr, $d:expr) => {
            if !($x - $y < $d || $y - $x < $d) { panic!(); }
        }
    }

    #[test]
    fn test_rm_extension() {
        assert_eq!(rm_extension(&"foo/bar/foobar.txt".to_string()), "foo/bar/foobar".to_string());
        assert_eq!(rm_extension(&"../foobar.txt".to_string()), "../foobar".to_string());
        assert_eq!(rm_extension(&"foo/bar/foo.bar.txt".to_string()), "foo/bar/foo.bar".to_string());
    }

    #[test]
    fn test_replace_ns() {
        let seed: &[_] = &[1, 2, 3, 4];
        let mut rng : StdRng = SeedableRng::from_seed(seed);
        assert_eq!(replace_ns_rand(&"ACGT".to_string(), &mut rng), "ACGT".to_string());
        assert!(replace_ns_rand(&"ACGTNNN".to_string(), &mut rng).starts_with("ACGT"));

        assert_eq!(replace_ns(&"ACGT".to_string(), "foo"), "ACGT".to_string());
        assert!(replace_ns(&"ACGTNNN".to_string(), "foo").starts_with("ACGT"));
    }

    #[test]
    fn test_drop_read() {
        assert_eq!(drop_read(0.0, "foo"), true);
        assert_eq!(drop_read(0.0000001, "foo"), true);
        assert_eq!(drop_read(0.9999999, "foo"), false);
        assert_eq!(drop_read(1.0, "foo"), false);

        let x = drop_read(0.5, "bar");
        for _ in 0..10 {
            assert_eq!(drop_read(0.5, "bar"), x);
        }

        let x = drop_read(0.5, "bar");
        for _ in 0..10 {
            assert_eq!(drop_read(0.5, "bar|||asdfasdf"), x);
        }
    }

    #[test]
    fn test_argmax() {
        let vec : Vec<usize> = vec![0, 1, 5, 3, 2];
        let (i, m) = argmax(&vec).unwrap();
        assert_eq!(i, 2);
        assert_eq!(m, 5);

        let vec : Vec<f32> = vec.iter().map(|x| *x as f32).collect();
        let (i, m) = argmax(&vec).unwrap();
        assert_eq!(i, 2);
        assert_eq!(m, 5.0);

        let vec : Vec<usize> = Vec::new();
        let m = argmax(&vec);
        assert!(m.is_none());
    }

    #[test]
    fn test_logaddexp() {
        let lse_err = 0.01;
        // Test against values computed using Python's numpy
        assert_delta!(logaddexp_arr(&vec![2.0, 3.0, 4.0, 5.0], 2.0), 5.90689059561, 1e-16);
        assert_delta!(logaddexp_arr(&vec![2.0, 3.0, 4.0, 5.0], 1.0_f64.exp()), 5.44018969856, 1e-16);

        // Test the approximation
        assert_delta!(logaddexp_arr(&vec![2.0, 3.0, 4.0, 5.0], 2.0), logaddexp_approx_arr(&vec![2.0, 3.0, 4.0, 5.0], 2.0), lse_err);
        assert_delta!(logaddexp_arr(&vec![2.0, 3.0, 4.0, 5.0], 1.0_f64.exp()), lnaddexp_approx_arr(&vec![2.0, 3.0, 4.0, 5.0]), lse_err);

        assert_eq!(logaddexp_arr(&vec![0.0, 0.0], 2.0), 1.0);
        assert_delta!(logaddexp_arr(&vec![1e-16, 1e-16], 2.0), 1.0, 1e-16);

        assert_delta!(logaddexp_arr(&vec![1.0, 1000000.0], 2.0), 1000000.0, 1e-16);
    }
}
