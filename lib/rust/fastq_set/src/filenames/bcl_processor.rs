//! Utilities for finding groups of FASTQ produced by the legacy `demux` pipeline from 10x Genomics.
#![expect(missing_docs)]

use crate::filenames::{LaneMode, LaneSpec};
use crate::read_pair_iter::FastqSet;
use crate::sseq::SSeq;
use anyhow::{Context, Result};
use itertools::Itertools;
use regex;
use serde::{Deserialize, Serialize};
use std::collections::HashSet;
use std::path::{Path, PathBuf};

/// Different ways to specify the sample index for `BclProcessorFastqDef`.
#[derive(Deserialize, Serialize, Clone, PartialEq, Eq, Debug)]
pub enum SampleIndexSpec {
    /// Allow all the sample indices **including X** irrespective of the number of N's in the SI
    Any,
    /// Allow sample indices from the given list of sequences. `indices` stores the list of allowed
    /// sequences and require that the SI match to within `max_n` Ns
    Sequences {
        indices: HashSet<SSeq>,
        max_n: usize,
    },
}

/// SampleIndexSpec from a single sequence
impl From<&str> for SampleIndexSpec {
    fn from(index: &str) -> Self {
        SampleIndexSpec::Sequences {
            indices: [index]
                .iter()
                .map(|ind| SSeq::from_bytes(ind.as_bytes()))
                .collect(),
            max_n: 1,
        }
    }
}

impl SampleIndexSpec {
    /// Does the given index match with the `SampleIndexSpec`
    pub fn matches(&self, index: &str) -> bool {
        match *self {
            SampleIndexSpec::Any => true,
            SampleIndexSpec::Sequences { ref indices, max_n } => indices
                .iter()
                .any(|target| match_seqs_with_n(target.as_bytes(), index.as_bytes(), max_n)),
        }
    }
}

/// A selector for a set of FASTQs emitted by the `BCL_PROCESSOR`
/// demultiplexing pipeline used internally at 10x, and in the
/// (now deprecated) `demux` customer command. `find_fastqs` will
/// select the set of FASTQs in `fastq_path` with the matching
/// sample index and lane values.
pub fn find_fastqs(
    read_path: &Path,
    sample_index_spec: &SampleIndexSpec,
    lane_spec: &LaneSpec,
) -> Result<Vec<FastqSet>> {
    // get all the files in the directory
    let all_fastq_sets = find_flowcell_fastqs(read_path)?;
    let mut res = Vec::new();

    for (bcl_proc, fastqs) in all_fastq_sets {
        if sample_index_spec.matches(&bcl_proc.si) && lane_spec.contains(bcl_proc.lane_mode()) {
            res.push(fastqs);
        }
    }

    res.sort();
    Ok(res)
}

fn match_seqs_with_n(target: &[u8], observed: &[u8], ns_allowed: usize) -> bool {
    if target.len() != observed.len() {
        return false;
    }

    let mut num_ns = 0;

    for (&tgt, &obs) in target.iter().zip_eq(observed.iter()) {
        if obs == b'N' {
            num_ns += 1;
        } else if obs != tgt {
            return false;
        }
    }

    num_ns <= ns_allowed
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct BclProcessorFileGroup {
    pub si: String,
    pub lane: usize,
    pub chunk: usize,
}

impl BclProcessorFileGroup {
    pub fn lane_mode(&self) -> LaneMode {
        LaneMode::SingleLane(self.lane)
    }
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct BclProcessorFile {
    pub group: BclProcessorFileGroup,
    pub read: String,
    pub full_path: PathBuf,
}

impl BclProcessorFile {
    /// Create a new `BclProcessorFile` from a full path to a FASTQ file.
    ///
    /// Returns `None` if the path is not a file or the filename does not
    /// match the expected pattern.
    fn new(full_path: &Path) -> Option<Self> {
        let re = "^read-([RI][A0-9])_si-([^_]+)_lane-([0-9]+)-chunk-([0-9]+).fastq(.gz|.lz4)?$";
        let re = regex::Regex::new(re).unwrap();
        full_path.file_name().and_then(|filename| {
            re.captures(filename.to_str().unwrap())
                .map(|caps| BclProcessorFile {
                    full_path: full_path.into(),
                    read: caps.get(1).unwrap().as_str().to_string(),
                    group: BclProcessorFileGroup {
                        si: caps.get(2).unwrap().as_str().to_string(),
                        lane: caps.get(3).unwrap().as_str().parse().unwrap(),
                        chunk: caps.get(4).unwrap().as_str().parse().unwrap(),
                    },
                })
        })
    }

    fn path_string(&self) -> String {
        self.full_path.to_str().unwrap().to_string()
    }
}

pub(crate) fn find_flowcell_fastqs(path: &Path) -> Result<Vec<(BclProcessorFileGroup, FastqSet)>> {
    let mut files = get_demux_files(path)?;
    files.sort();
    Ok(files
        .chunk_by(|x, y| x.group == y.group)
        .filter_map(|group_files| {
            group_files
                .iter()
                // ignore groups that don't have an RA file
                .find(|info| info.read == "RA")
                .map(|ra| {
                    let i1 = group_files.iter().find(|info| info.read == "I1");
                    let i2 = group_files.iter().find(|info| info.read == "I2");

                    let fastqs = FastqSet {
                        r1: ra.path_string(),
                        r2: None,
                        i1: i1.map(BclProcessorFile::path_string),
                        i2: i2.map(BclProcessorFile::path_string),
                        r1_interleaved: true,
                    };
                    (group_files[0].group.clone(), fastqs)
                })
        })
        .collect())
}

fn get_demux_files(path: &Path) -> Result<Vec<BclProcessorFile>> {
    std::fs::read_dir(path)
        .with_context(|| path.display().to_string())?
        .filter_map_ok(|x| BclProcessorFile::new(&x.path()))
        .try_collect()
        .with_context(|| path.display().to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bcl_fn() {
        let f1 = "read-RA_si-TTAGCGAT_lane-002-chunk-001.fastq.gz";

        let bcl = BclProcessorFile::new(Path::new(f1)).unwrap();

        let truth = BclProcessorFile {
            full_path: PathBuf::from(f1),
            read: "RA".to_string(),
            group: BclProcessorFileGroup {
                si: "TTAGCGAT".to_string(),
                lane: 2,
                chunk: 1,
            },
        };

        assert_eq!(bcl, truth);
    }

    #[test]
    fn load_dir() -> Result<()> {
        let path = Path::new("tests/filenames/bcl_processor");
        let fastqs = find_flowcell_fastqs(path)?;
        assert_eq!(fastqs.len(), 44);
        Ok(())
    }

    #[test]
    fn query_all_lanes() -> Result<()> {
        let path = Path::new("tests/filenames/bcl_processor");

        let fqs = find_fastqs(path, &"TCGAATGATC".into(), &LaneSpec::Any)?;

        assert_eq!(fqs.len(), 2);
        Ok(())
    }

    #[test]
    fn query_one_lane() -> Result<()> {
        let path = "tests/filenames/bcl_processor";

        let fqs = find_fastqs(
            Path::new(path),
            &"TCGAATGATC".into(),
            &LaneSpec::Lanes(vec![2].into_iter().collect()),
        )?;
        assert_eq!(fqs.len(), 1);
        assert_eq!(
            fqs[0].r1,
            "tests/filenames/bcl_processor/read-RA_si-TCGAATGATC_lane-002-chunk-001.fastq.gz"
        );
        Ok(())
    }

    #[test]
    fn test_si_any() {
        let fqs = find_fastqs(
            Path::new("tests/filenames/bcl_processor/"),
            &SampleIndexSpec::Any,
            &LaneSpec::Any,
        )
        .unwrap();
        assert_eq!(fqs.len(), 44);
    }
}

// The following tests are based on the tests in tenkit: lib/python/tenkit/test/test_fasta.py
#[cfg(test)]
mod tests_from_tenkit {
    use super::*;
    use pretty_assertions::assert_eq;

    #[test]
    fn test_find_input_fastq_files_10x_preprocess() -> Result<()> {
        let path = "tests/filenames/bcl_processor_2";

        let fqs = find_fastqs(
            Path::new(path),
            &"ACAGCAAC".into(),
            &LaneSpec::Lanes(vec![1, 2].into_iter().collect()),
        )?;

        let expected = vec![
            FastqSet {
                r1: format!("{path}/read-RA_si-ACAGCAAC_lane-001-chunk-001.fastq.gz"),
                r2: None,
                i1: Some(format!(
                    "{path}/read-I1_si-ACAGCAAC_lane-001-chunk-001.fastq.gz"
                )),
                i2: Some(format!(
                    "{path}/read-I2_si-ACAGCAAC_lane-001-chunk-001.fastq.gz"
                )),
                r1_interleaved: true,
            },
            FastqSet {
                r1: format!("{path}/read-RA_si-ACAGCAAC_lane-002-chunk-000.fastq.gz"),
                r2: None,
                i1: Some(format!(
                    "{path}/read-I1_si-ACAGCAAC_lane-002-chunk-000.fastq.gz"
                )),
                i2: Some(format!(
                    "{path}/read-I2_si-ACAGCAAC_lane-002-chunk-000.fastq.gz"
                )),
                r1_interleaved: true,
            },
        ];

        assert_eq!(fqs, expected);

        Ok(())
    }

    #[test]
    fn test_ra_missing() -> Result<()> {
        let path = "tests/filenames/bcl_processor_3";
        let fastqs = find_fastqs(
            Path::new(path),
            &"ACAGCAAC".into(),
            &LaneSpec::Lanes(vec![1, 2].into_iter().collect()),
        )?;
        let expected = vec![FastqSet {
            r1: format!("{path}/read-RA_si-ACAGCAAC_lane-001-chunk-001.fastq.gz"),
            r2: None,
            i1: Some(format!(
                "{path}/read-I1_si-ACAGCAAC_lane-001-chunk-001.fastq.gz"
            )),
            i2: Some(format!(
                "{path}/read-I2_si-ACAGCAAC_lane-001-chunk-001.fastq.gz"
            )),
            r1_interleaved: true,
        }];
        assert_eq!(fastqs, expected);

        Ok(())
    }
}
