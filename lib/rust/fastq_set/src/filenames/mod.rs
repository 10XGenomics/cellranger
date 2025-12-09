//! Utilities for finding groups of FASTQ files on disk.
#![expect(missing_docs)]

pub mod bcl2fastq;
pub mod bcl_processor;
pub mod fastq_dir;

use crate::read_pair_iter::InputFastqs;
use anyhow::Error;
pub use bcl_processor::BclProcessorFastqDef;
use bcl_processor::SampleIndexSpec;
pub use bcl2fastq::Bcl2FastqDef;
use bcl2fastq::SampleNameSpec;
use serde::{Deserialize, Serialize};
use std::collections::HashSet;

/// A method to find a set of `InputFastqs` based on
/// some configuration information held by `self`,
/// and some conventions encoded in the implementing
/// type
pub trait FindFastqs {
    fn find_fastqs(&self) -> Result<Vec<InputFastqs>, Error>;
}

/// This enum stores the lane information associated with a fastq file:
/// - `NoLaneSplitting`: The fastq file contains data from all lanes. Such files are usually
///   generated using the --no-lane-splitting option to bcl2fastq
/// - `SingleLane(usize)`: The fastq file contains data from a single lane.
#[derive(Deserialize, Serialize, Clone, Copy, PartialEq, Eq, Debug, PartialOrd, Ord)]
pub enum LaneMode {
    NoLaneSplitting,
    SingleLane(usize),
}

impl From<usize> for LaneMode {
    fn from(lane: usize) -> Self {
        LaneMode::SingleLane(lane)
    }
}

#[derive(Deserialize, Serialize, Clone, PartialEq, Eq, Debug)]
pub enum LaneSpec {
    /// Consider all the lanes
    Any,
    /// Only consider the given set of lanes
    Lanes(HashSet<usize>),
}

impl LaneSpec {
    pub fn contains(&self, lane_mode: LaneMode) -> bool {
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

/// A pointer to FASTQ data on disk. Can be encoded in the standard Illumina
/// 'bcl2fastq' naming convention, or in the 10x-specific 'BclProcessor'
/// convention. Use the `find_fastqs()` method to find the concrete
/// `InputFastq` files corresponding to a `FastqDef`.
#[derive(Deserialize, Serialize, Clone, PartialEq, Eq, Debug)]
pub enum FastqDef {
    Bcl2Fastq(Bcl2FastqDef),
    BclProcessor(BclProcessorFastqDef),
}

impl FastqDef {
    pub fn bcl2fastq(
        fastq_path: String,
        sample_name_spec: SampleNameSpec,
        lane_spec: LaneSpec,
    ) -> FastqDef {
        FastqDef::Bcl2Fastq(Bcl2FastqDef {
            fastq_path,
            sample_name_spec,
            lane_spec,
        })
    }
    pub fn bcl_processor(
        fastq_path: String,
        sample_index_spec: SampleIndexSpec,
        lane_spec: LaneSpec,
    ) -> FastqDef {
        FastqDef::BclProcessor(BclProcessorFastqDef {
            fastq_path,
            sample_index_spec,
            lane_spec,
        })
    }
}

impl FindFastqs for FastqDef {
    fn find_fastqs(&self) -> Result<Vec<InputFastqs>, Error> {
        match self {
            FastqDef::Bcl2Fastq(d) => d.find_fastqs(),
            FastqDef::BclProcessor(d) => d.find_fastqs(),
        }
    }
}
