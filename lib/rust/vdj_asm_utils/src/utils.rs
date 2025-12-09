#![deny(missing_docs)]

#[cfg(test)]
use std::hash::{Hash, Hasher};
#[cfg(test)]
use std::path::{Path, PathBuf};

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
pub fn logaddexp_arr(vals: &[f64], base: f64) -> f64 {
    assert!(vals.len() > 1);
    let mut res = logaddexp(vals[0], vals[1], base);

    for &v in vals.iter().skip(2) {
        res = logaddexp(res, v, base);
    }
    res
}

#[cfg(test)]
pub fn rm_extension(path: &Path) -> PathBuf {
    path.with_extension("")
}

#[cfg(test)]
pub fn base_to_bits_hash(c: u8, name: &str, pos: usize) -> u8 {
    use metric::TxHasher;

    match c {
        b'A' => 0u8,
        b'C' => 1u8,
        b'G' => 2u8,
        b'T' => 3u8,
        _ => {
            let mut hasher = TxHasher::hasher();
            name.hash(&mut hasher);
            pos.hash(&mut hasher);
            (hasher.finish() % 4) as u8
        }
    }
}

#[cfg(test)]
pub mod tests {
    use super::*;

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
    pub fn lnaddexp_approx_arr(vals: &[f64]) -> f64 {
        assert!(vals.len() > 1);
        let mut res = vals.to_owned();
        let mut step = 2;
        loop {
            let mut i = 0;
            while i < (vals.len() - step / 2) {
                let tmp = lnaddexp_approx(res[i], res[i + step / 2]);
                res[i] = tmp;
                i += step;
            }
            if step >= vals.len() {
                break;
            }
            step *= 2;
        }
        res[0]
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
    fn logaddexp_approx_arr(vals: &[f64], base: f64) -> f64 {
        let lnbase = base.ln();
        let new_vals = vals.iter().map(|&x| x * lnbase).collect::<Vec<_>>();
        lnaddexp_approx_arr(&new_vals) / lnbase
    }

    // Computes ln(exp(a) + exp(b)) using a piecewise linear approximation
    // See https://web.stanford.edu/~boyd/papers/rgp.html for details
    //
    // The coefficients are precomputed and listed in constants::LSE_COEFF
    fn lnaddexp_approx(a: f64, b: f64) -> f64 {
        use crate::constants::LSE_COEFF;
        let c: f64;
        let d: f64;
        if a > b {
            c = a;
            d = b;
        } else {
            c = b;
            d = a;
        }

        LSE_COEFF.into_iter().fold(c, |result, lse| {
            result.max(lse[0] * d + lse[1] * c + lse[2])
        })
    }
    fn argmax<T: PartialOrd + Copy>(scores: &[T]) -> Option<(usize, T)> {
        if scores.is_empty() {
            return None;
        }
        let mut max_index = 0;
        let mut max = scores[max_index];
        for (i, &v) in scores.iter().enumerate() {
            if max < v {
                max_index = i;
                max = v;
            }
        }
        Some((max_index, max))
    }

    macro_rules! assert_delta {
        ($x:expr, $y:expr, $d:expr) => {
            if !($x - $y < $d || $y - $x < $d) {
                panic!();
            }
        };
    }

    #[test]
    fn test_rm_extension() {
        assert_eq!(
            rm_extension(Path::new("foo/bar/foobar.txt")).as_os_str(),
            "foo/bar/foobar"
        );
        assert_eq!(
            rm_extension(Path::new("../foobar.txt")).as_os_str(),
            "../foobar"
        );
        assert_eq!(
            rm_extension(Path::new("foo/bar/foo.bar.txt")).as_os_str(),
            "foo/bar/foo.bar"
        );
    }

    #[test]
    #[allow(clippy::float_cmp)]
    fn test_argmax() {
        let vec: Vec<usize> = vec![0, 1, 5, 3, 2];
        let (i, m) = argmax(&vec).unwrap();
        assert_eq!(i, 2);
        assert_eq!(m, 5);

        let vec: Vec<f32> = vec.iter().map(|x| *x as f32).collect();
        let (i, m) = argmax(&vec).unwrap();
        assert_eq!(i, 2);
        assert_eq!(m, 5.0);

        let vec: Vec<usize> = Vec::new();
        let m = argmax(&vec);
        assert!(m.is_none());
    }

    #[test]
    #[allow(clippy::float_cmp)]
    fn test_logaddexp() {
        let lse_err = 0.01;
        // Test against values computed using Python's numpy
        assert_delta!(
            logaddexp_arr(&[2.0, 3.0, 4.0, 5.0], 2.0),
            5.90689059561,
            1e-16
        );
        assert_delta!(
            logaddexp_arr(&[2.0, 3.0, 4.0, 5.0], 1.0_f64.exp()),
            5.44018969856,
            1e-16
        );

        // Test the approximation
        assert_delta!(
            logaddexp_arr(&[2.0, 3.0, 4.0, 5.0], 2.0),
            logaddexp_approx_arr(&[2.0, 3.0, 4.0, 5.0], 2.0),
            lse_err
        );
        assert_delta!(
            logaddexp_arr(&[2.0, 3.0, 4.0, 5.0], 1.0_f64.exp()),
            lnaddexp_approx_arr(&[2.0, 3.0, 4.0, 5.0]),
            lse_err
        );

        assert_eq!(logaddexp_arr(&[0.0, 0.0], 2.0), 1.0);
        assert_delta!(logaddexp_arr(&[1e-16, 1e-16], 2.0), 1.0, 1e-16);

        assert_delta!(logaddexp_arr(&[1.0, 1000000.0], 2.0), 1000000.0, 1e-16);
    }
}
