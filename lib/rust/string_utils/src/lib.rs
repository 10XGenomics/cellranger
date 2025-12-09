// Copyright (c) 2018 10X Genomics, Inc. All rights reserved.
//! This file contains some miscellaneous string utilities.
#![deny(missing_docs)]

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
// THINGS USED A LOT: SHORTHAND EXPRESSIONS FOR COMMON FUNCTIONALITY
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

/// Additional methods for &str.
pub trait TextUtils<'a> {
    /// Convert this string to usize or panic.
    fn force_usize(&self) -> usize;

    /// Convert this string to i32 or panic.
    fn force_i32(&self) -> i32;

    /// s.before(t): return the part of s before the first instance of t
    /// (or panic if t is not contained in s)
    fn before(&'a self, u: &str) -> &'a str;

    /// s.after(t): return the part of s after the first instance of t
    /// (or panic if t is not contained in s)
    fn after(&'a self, t: &str) -> &'a str;

    // s.between(t,u): return the part of s after the first instance of t and
    /// before the first instance of u after that
    fn between(&'a self, t: &str, u: &str) -> &'a str;

    /// s.rev_before(t): start from the end s, find the first instance of t, and
    /// return what's before that
    fn rev_before(&'a self, t: &str) -> &'a str;
}

impl<'a> TextUtils<'a> for str {
    fn force_usize(&self) -> usize {
        self.parse::<usize>()
            .unwrap_or_else(|_| panic!("could not convert \"{self}\" to usize"))
    }
    fn force_i32(&self) -> i32 {
        self.parse::<i32>()
            .unwrap_or_else(|_| panic!("could not convert \"{self}\" to i32"))
    }

    fn before(&'a self, u: &str) -> &'a str {
        let r = self
            .find(u)
            .unwrap_or_else(|| panic!("failed to find \"{u}\" in \"{self}\""));
        &self[0..r]
    }

    fn after(&'a self, t: &str) -> &'a str {
        let l = self
            .find(t)
            .unwrap_or_else(|| panic!("after failed to find \"{t}\" in \"{self}\""))
            + t.len();
        &self[l..self.len()]
    }

    fn between(&'a self, t: &str, u: &str) -> &'a str {
        let a = self.after(t);
        let r = a.find(u).unwrap_or_else(|| {
            panic!("between( \"{self}\", \"{t}\", \"{u}\" ) failed at second part")
        });
        &a[0..r]
    }

    fn rev_before(&'a self, t: &str) -> &'a str {
        let l = 0;
        let r = self.rfind(t).unwrap();
        &self[l..r]
    }
}
