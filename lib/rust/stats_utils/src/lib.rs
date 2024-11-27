// Copyright (c) 2018 10X Genomics, Inc. All rights reserved.

// Compute some stats.
// ◼ Presumably much of this is available elsewhere.

use std::f64;

// Compute N50 of some numbers.

pub fn n50(v: &[i32]) -> i32 {
    if v.is_empty() {
        return 0;
    }
    for n in v {
        assert!(*n > 0);
    }
    let mut sum: i64 = 0;
    let mut half: i64 = 0;
    for n in v {
        sum += i64::from(*n);
    }
    let mut vs = v.to_owned();
    vs.sort_unstable();
    for i in 0..vs.len() {
        half += i64::from(vs[i]);
        if 2 * half == sum && i < vs.len() - 1 {
            return (vs[i] + vs[i + 1]) / 2;
        }
        if 2 * half >= sum {
            return vs[i];
        }
    }
    0 // never executed
}

pub fn n90(v: &[i32]) -> i32 {
    if v.is_empty() {
        return 0;
    }
    for n in v {
        assert!(*n > 0);
    }
    let mut sum: i64 = 0;
    let mut part: i64 = 0;
    for n in v {
        sum += i64::from(*n);
    }
    let mut vs = v.to_owned();
    vs.sort_unstable();
    for i in 0..vs.len() {
        part += i64::from(vs[i]);
        if 10 * part == 9 * sum && i < vs.len() - 1 {
            return (vs[i] + vs[i + 1]) / 2;
        }
        if part as f64 / sum as f64 >= 0.9_f64 {
            return vs[i];
        }
    }
    0 // never executed
}

// Compute mean of some numbers, returning zero on empty vector.

pub fn mean(v: &[i32]) -> f64 {
    let sum1 = v.len() as f64;
    let mut sum2 = 0_f64;
    for x in v {
        sum2 += f64::from(*x);
    }
    if sum1 == 0_f64 {
        return 0_f64;
    }
    sum2 / sum1
}

// Compute "length-weighted" mean of some numbers.  Normally this would be applied
// to positive numbers.  Returns 0 if the sum of the numbers is zero.
// ◼ Not sure what the right term for this is.  The word "length" is confusing.

pub fn len_weighted_mean(v: &[i32]) -> f64 {
    let mut sum1 = 0_f64;
    let mut sum2 = 0_f64;
    for x in v {
        sum1 += f64::from(*x);
        sum2 += f64::from(*x) * f64::from(*x);
    }
    if sum1 == 0_f64 {
        return 0_f64;
    }
    sum2 / sum1
}

// Compute absolute value of difference between two numbers.

pub fn abs_diff(a: usize, b: usize) -> usize {
    if a <= b {
        return b - a;
    }
    a - b
}

// Compute percent ratio.

pub fn percent_ratio(a: usize, b: usize) -> f64 {
    100_f64 * a as f64 / b as f64
}

// binomial_sum( n, k, p ): return sum_{i=0..k} choose(n,i) * p^i * (1-p)^(n-i)
//
// No attempt has been made to make this efficient or to pay attention to
// accuracy or overflow problems.

pub fn binomial_sum(n: usize, k: usize, p: f64) -> f64 {
    assert!(n >= 1);
    assert!(k <= n);
    let mut sum = 0.0;
    let mut choose = 1.0;
    for _ in 0..n {
        choose *= 1.0 - p;
    }
    let q = p / (1.0 - p);
    for i in 0..=k {
        sum += choose;
        choose *= (n - i) as f64;
        choose /= (i + 1) as f64;
        choose *= q;
    }
    sum
}
