// This file contains miscellaneous utilities for input and output.
// ◼  This is a total hodge-podge that should be pulled apart at some point.

use io_utils::{fwriteln, open_for_read};
use itertools::Itertools;
use std::collections::HashMap;
use std::io::prelude::*;
use std::iter::zip;
use std::path::Path;
use string_utils::parse_csv;
use vector_utils::next_diff;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
// PRINT COMPRESSED VEC OF INTEGERS
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Print vec of integers in a compressed form.
// E.g. 35, 40, 40, 40, 60, 60 is shown as 35, 40^3, 60^2.
// ◼ This is functionally the same as abbrev_list in string_utils.rs.

pub fn print_compressed(log: &mut Vec<u8>, q: &[u8]) {
    let mut i = 0;
    while i < q.len() {
        if i > 0 {
            write!(log, ", ").unwrap();
        }
        let j = next_diff(q, i);
        if j - i == 1 {
            write!(log, "{}", q[i]).unwrap();
        } else {
            write!(log, "{}^{}", q[i], j - i).unwrap();
        }
        i = j;
    }
    fwriteln!(log, "");
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
// GET METRIC VALUES
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Get all the metrics from a two-line customer CSV file "metrics_summary.csv".

pub fn metrics_from_csv(path: &Path) -> HashMap<String, String> {
    let Some((names, values)) = open_for_read!(path).lines().collect_tuple() else {
        panic!("Expected exactly two lines: {}", path.display());
    };
    zip(parse_csv(&names.unwrap()), parse_csv(&values.unwrap())).collect()
}
