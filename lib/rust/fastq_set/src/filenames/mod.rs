//! Utilities for finding groups of FASTQ files on disk.
#![deny(missing_docs)]

pub mod bcl2fastq;
pub mod bcl_processor;
pub mod fastq_dir;
pub mod sbx;

use serde::{Deserialize, Serialize};
use std::collections::HashSet;

/// This enum stores the lane information associated with a fastq file:
#[derive(Deserialize, Serialize, Clone, Copy, PartialEq, Eq, Debug, PartialOrd, Ord)]
pub enum LaneMode {
    /// The fastq file contains data from all lanes.
    /// Such files are usually generated using the `--no-lane-splitting` option of `bcl2fastq`.
    NoLaneSplitting,
    /// The fastq file contains data from a single lane.
    SingleLane(usize),
}

impl From<usize> for LaneMode {
    fn from(lane: usize) -> Self {
        LaneMode::SingleLane(lane)
    }
}

/// Specification of which lanes to consider when finding FASTQ files.
#[derive(Deserialize, Serialize, Clone, PartialEq, Eq, Debug)]
pub enum LaneSpec {
    /// Consider all the lanes.
    Any,
    /// Only consider the given set of lanes.
    Lanes(HashSet<usize>),
}

impl LaneSpec {
    fn contains(&self, lane_mode: LaneMode) -> bool {
        match self {
            LaneSpec::Any => true,
            LaneSpec::Lanes(lanes) => match lane_mode {
                LaneMode::NoLaneSplitting => panic!(
                    "Fastq files are generated without splitting by lane. You cannot restrict the \
                    selection of fastq files by specifying lanes of interest in such a case."
                ),
                LaneMode::SingleLane(lane) => lanes.contains(&lane),
            },
        }
    }
}
