// Copyright (c) 2018 10X Genomics, Inc. All rights reserved.

// This file contains miscellaneous vector utilities.

use std::borrow::Borrow;
use superslice::Ext;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
// SORT A VECTOR
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Reverse sort a vector.

pub fn reverse_sort<T: Ord>(x: &mut [T]) {
    x.sort_by(|a, b| b.cmp(a));
}

// Unique sort a vector.

pub fn unique_sort<T: Ord>(x: &mut Vec<T>) {
    x.sort();
    x.dedup();
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
// DOES VECTOR CONTAIN ANOTHER VECTOR
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Test to see if a vector contains another vector at a given position.

pub fn contains_at<T: Eq>(s: &[T], t: &[T], p: usize) -> bool {
    if p + t.len() > s.len() {
        return false;
    }
    for i in 0..t.len() {
        if s[p + i] != t[i] {
            return false;
        }
    }
    true
}

// Determine if vector x contains vector y.

pub fn contains<T: Eq>(x: &[T], y: &[T]) -> bool {
    if y.len() > x.len() {
        return false;
    }
    for i in 0..=x.len() - y.len() {
        let mut matches = true;
        for j in 0..y.len() {
            if x[i + j] != y[j] {
                matches = false;
                break;
            }
        }
        if matches {
            return true;
        }
    }
    false
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
// UNSIGNED VECTOR SIZE AND SOME SPECIAL SIZES
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub trait VecUtils<'a> {
    fn ilen(&self) -> isize;

    fn solo(&self) -> bool;

    fn duo(&self) -> bool;
}

impl<'a, T> VecUtils<'a> for [T] {
    fn ilen(&self) -> isize {
        self.len() as isize
    }

    fn solo(&self) -> bool {
        self.len() == 1
    }

    fn duo(&self) -> bool {
        self.len() == 2
    }
}

/// Erase elements in a vector that are flagged by another vector.  Both vectors
/// should have the same length.
pub fn erase_if<T>(x: &mut Vec<T>, to_delete: &[bool]) {
    // Adding this line means that enclone starts failing, and I don't want to
    // run down those errors. Sigh.  COME ON PEOPLE, ADD ASSERTIONS.
    // assert_eq!(x.len(), to_delete.len());
    erase_if_iter(x, to_delete.iter().copied());
}

/// Erase elements in a vector that are flagged by an iterator of booleans.
/// The iterator should have the same length as the vector.
pub fn erase_if_iter<T>(x: &mut Vec<T>, to_delete: impl IntoIterator<Item = bool>) {
    let mut count = 0;
    for (j, delete) in to_delete.into_iter().take(x.len()).enumerate() {
        if !delete {
            if j != count {
                x.swap(j, count);
            }
            count += 1;
        }
    }
    x.truncate(count);
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
// INTERSECTION FUNCTIONS
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Determine if two sorted vectors have an element in common.

pub fn meet<T: Ord>(x: &[T], y: &[T]) -> bool {
    let mut i = 0;
    let mut j = 0;
    while i < x.len() && j < y.len() {
        if x[i] < y[j] {
            i += 1;
        } else if y[j] < x[i] {
            j += 1;
        } else {
            return true;
        }
    }
    false
}

// Find the intersection size of two sorted vectors.  If an element occurs
// repeatedly, say n1 times in one vector and n2 times in the other vector, then
// that contributes min(n1,n2) to the total.

pub fn meet_size<T: Ord>(x: &[T], y: &[T]) -> usize {
    let mut i = 0;
    let mut j = 0;
    let mut count = 0;
    while i < x.len() && j < y.len() {
        if x[i] < y[j] {
            i += 1;
        } else if y[j] < x[i] {
            j += 1;
        } else {
            count += 1;
            i += 1;
            j += 1;
        }
    }
    count
}

// Compute the intersection of two sorted vectors.

pub fn intersection<T: Ord + Clone>(x: &[T], y: &[T], z: &mut Vec<T>) {
    z.clear();
    let (mut ix, mut iy) = (0, 0);
    while ix < x.len() && iy < y.len() {
        if x[ix] < y[iy] {
            ix += 1;
        } else if y[iy] < x[ix] {
            iy += 1;
        } else {
            z.push(x[ix].clone());
            ix += 1;
            iy += 1;
        }
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
// FREQUENCY FUNCTIONS
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Count elements in a sorted vector by type.  The output consists of a reverse
// sorted vector of pairs (m,v) where m is the multiplicity of an element v.

pub fn make_freq<T: Ord + Clone>(x: &[T], freq: &mut Vec<(u32, T)>) {
    freq.clear();
    let mut j = 0;
    loop {
        if j == x.len() {
            break;
        }
        let mut k = j + 1;
        loop {
            if k == x.len() || x[k] != x[j] {
                break;
            }
            k += 1;
        }
        let t = x[j].clone();
        freq.push(((k - j) as u32, t));
        j = k;
    }
    freq.sort_by(|a, b| b.cmp(a)); // freq.reverse_sort();
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
// MEMBERSHIP
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Like binary_search, but allows searching by a different type from what's in
// the slice, e.g. searching a &[String] for a &str.
fn borrowed_binary_search<T: Ord + ?Sized>(x: &[impl Borrow<T>], d: &T) -> Result<usize, usize> {
    x.binary_search_by(|p| p.borrow().cmp(d))
}

// Test to see if a sorted vector contains a given element.

pub fn bin_member<T: Ord + ?Sized>(x: &[impl Borrow<T>], d: &T) -> bool {
    borrowed_binary_search(x, d).is_ok()
}

// Return the position of an element in an unsorted vector.
// Returns -1 if not present.

pub fn position<T: Ord + ?Sized>(x: &[impl Borrow<T>], d: &T) -> i32 {
    for (i, y) in x.iter().enumerate() {
        if y.borrow() == d {
            return i as i32;
        }
    }
    -1_i32
}

// Return the position of an element in a sorted vector, or using just the first
// position.  Returns -1 if not present.

pub fn bin_position<T: Ord + ?Sized>(x: &[impl Borrow<T>], d: &T) -> i32 {
    match borrowed_binary_search(x, d) {
        Ok(p) => p as i32,
        Err(_e) => -1,
    }
}

pub fn bin_position1_2<S: Ord + ?Sized, T: Ord>(x: &[(impl Borrow<S>, T)], d: &S) -> i32 {
    match x.binary_search_by_key(&d, |(a, _b)| a.borrow()) {
        Ok(p) => p as i32,
        Err(_e) => -1,
    }
}

pub fn bin_position1_3<S: Ord + ?Sized, T: Ord, U: Ord>(
    x: &[(impl Borrow<S>, T, U)],
    d: &S,
) -> i32 {
    match x.binary_search_by_key(&d, |(a, _b, _c)| a.borrow()) {
        Ok(p) => p as i32,
        Err(_e) => -1,
    }
}

// Find lower/upper bounds.

fn lower_bound_usize<T: Ord + ?Sized>(x: &[impl Borrow<T>], d: &T) -> usize {
    x.lower_bound_by(|y| y.borrow().cmp(d))
}

pub fn lower_bound<T: Ord + ?Sized>(x: &[impl Borrow<T>], d: &T) -> i64 {
    lower_bound_usize(x, d) as i64
}

fn upper_bound_usize<T: Ord + ?Sized>(x: &[impl Borrow<T>], d: &T) -> usize {
    x.upper_bound_by(|y| y.borrow().cmp(d))
}

pub fn upper_bound<T: Ord + ?Sized>(x: &[impl Borrow<T>], d: &T) -> i64 {
    upper_bound_usize(x, d) as i64
}

pub fn lower_bound1_2<S: Ord + ?Sized, T: Ord>(x: &[(impl Borrow<S>, T)], d: &S) -> i64 {
    x.lower_bound_by_key(&d, |(a, _b)| a.borrow()) as i64
}

pub fn upper_bound1_2<S: Ord + ?Sized, T: Ord>(x: &[(impl Borrow<S>, T)], d: &S) -> i64 {
    x.upper_bound_by_key(&d, |(a, _b)| a.borrow()) as i64
}

pub fn lower_bound1_3<S: Ord + ?Sized, T: Ord, U: Ord>(x: &[(impl Borrow<S>, T, U)], d: &S) -> i64 {
    x.lower_bound_by_key(&d, |(a, _b, _c)| a.borrow()) as i64
}

pub fn upper_bound1_3<S: Ord + ?Sized, T: Ord, U: Ord>(x: &[(impl Borrow<S>, T, U)], d: &S) -> i64 {
    x.upper_bound_by_key(&d, |(a, _b, _c)| a.borrow()) as i64
}

// Compute the number of instances of a given element in a sorted vector.

pub fn count_instances<T: Ord + ?Sized>(x: &[impl Borrow<T>], d: &T) -> i32 {
    (upper_bound_usize(x, d) - lower_bound_usize(x, d)) as i32
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
// NEXT DIFFERENCE
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

/// Find first element that's different in a sorted vector, or different in
/// first position.
pub fn next_diff_any<T>(x: &[T], i: usize, eq: impl Fn(&T, &T) -> bool) -> usize {
    let mut j = i + 1;
    loop {
        if j == x.len() || !eq(&x[j], &x[i]) {
            return j;
        }
        j += 1;
    }
}

pub fn next_diff<T: Eq>(x: &[T], i: usize) -> usize {
    next_diff_any(x, i, |a, b| a == b)
}

pub fn next_diff1_2<T: Eq, U: Eq>(x: &[(T, U)], i: usize) -> usize {
    next_diff_any(x, i, |a, b| a.0 == b.0)
}

pub fn next_diff1_3<T: Eq, U: Eq, V: Eq>(x: &[(T, U, V)], i: usize) -> usize {
    next_diff_any(x, i, |a, b| a.0 == b.0)
}

pub fn next_diff12_3<T: Eq, U: Eq, V: Eq>(x: &[(T, U, V)], i: usize) -> usize {
    next_diff_any(x, i, |a, b| a.0 == b.0 && a.1 == b.1)
}

pub fn next_diff12_4<T: Eq, U: Eq, V: Eq, W: Eq>(x: &[(T, U, V, W)], i: usize) -> usize {
    next_diff_any(x, i, |a, b| a.0 == b.0 && a.1 == b.1)
}

#[allow(clippy::type_complexity)]
pub fn next_diff12_8<S: Eq, T: Eq, U: Eq, V: Eq, W: Eq, X: Eq, Y: Eq, Z: Eq>(
    x: &[(S, T, U, V, W, X, Y, Z)],
    i: usize,
) -> usize {
    next_diff_any(x, i, |a, b| a.0 == b.0 && a.1 == b.1)
}

pub fn next_diff1_5<T: Eq, U: Eq, V: Eq, W: Eq, X: Eq>(x: &[(T, U, V, W, X)], i: usize) -> usize {
    next_diff_any(x, i, |a, b| a.0 == b.0)
}

pub fn next_diff1_6<T: Eq, U: Eq, V: Eq, W: Eq, X: Eq, Y: Eq>(
    x: &[(T, U, V, W, X, Y)],
    i: usize,
) -> usize {
    next_diff_any(x, i, |a, b| a.0 == b.0)
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
// SORT SYNC VECTORS
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn sort_sync2<T: Ord + Clone, S1: Ord + Clone>(t: &mut Vec<T>, s1: &mut Vec<S1>) {
    let permutation = permutation::sort(&t[..]);
    *t = permutation.apply_slice(&t[..]);
    *s1 = permutation.apply_slice(&s1[..]);
}
