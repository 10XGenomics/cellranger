//! Functions and traits useful for testing the martian stages

pub mod correctness;
pub mod tools;

pub use tools::{
    diff_metrics, ensure_identical_set_of_lines, ensure_no_diff, safe_copy, set_permissions,
};
