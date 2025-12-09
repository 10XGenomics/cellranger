//! Utilities for finding groups of FASTQ produced by `bcl2fastq` from Illumina.
#![expect(missing_docs)]

use super::FindFastqs;
use crate::filenames::{LaneMode, LaneSpec};
use crate::read_pair_iter::InputFastqs;
use anyhow::{Context, Result};
use itertools::Itertools;
use lazy_static::lazy_static;
use regex::Regex;
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};

lazy_static! {
    static ref BCL2FASTQ_REGEX: Regex =
        Regex::new(r"^([\w_-]+)_S(\d+)_L(\d+)_([RI][A123])_(\d+).fastq(.gz|.lz4)?$").unwrap();
    static ref BCL2FASTQ_NO_LANE_SPLIT_REGEX: Regex =
        Regex::new(r"^([\w_-]+)_S(\d+)_([RI][A123])_(\d+).fastq(.gz|.lz4)?$").unwrap();
}

/// Different ways to specify sample names for the `Bcl2FastqDef`
#[derive(Deserialize, Serialize, Clone, PartialEq, Eq, Debug)]
pub enum SampleNameSpec {
    /// All the samples within the fastq directory
    Any,
    /// Only consider a known set of names
    Names(HashSet<String>),
}

/// SampleNameSpec for a single sample name
impl From<&str> for SampleNameSpec {
    fn from(sample_name: &str) -> Self {
        let mut names = HashSet::default();
        names.insert(sample_name.to_string());
        SampleNameSpec::Names(names)
    }
}

impl SampleNameSpec {
    pub fn contains(&self, sample_name: &str) -> bool {
        match self {
            SampleNameSpec::Any => true,
            SampleNameSpec::Names(names) => names.contains(sample_name),
        }
    }
}

/// A pointer to a set of FASTQ files on disk,
/// using the Illumina `bcl2fastq` naming conventions.
/// The `find_fastqs` method will find FASTQ files
/// of the form `heart_1k_v3_S1_L002_R2_001.fastq.gz`
/// with an optional `.gz` suffix.
#[derive(Deserialize, Serialize, Clone, PartialEq, Eq, Debug)]
pub struct Bcl2FastqDef {
    /// The path where to the demulitplexed FASTQ files
    pub fastq_path: String,

    /// Sample name(s) used for this sample
    pub sample_name_spec: SampleNameSpec,

    /// Lanes to include
    pub lane_spec: LaneSpec,
}

impl FindFastqs for Bcl2FastqDef {
    fn find_fastqs(&self) -> Result<Vec<InputFastqs>> {
        let all_fastqs = find_flowcell_fastqs(&self.fastq_path)?;

        let mut res = Vec::new();
        for (info, fastqs) in all_fastqs {
            if self.sample_name_spec.contains(&info.sample)
                && self.lane_spec.contains(info.lane_mode)
            {
                res.push(fastqs);
            }
        }

        res.sort();
        Ok(res)
    }
}

/// Metadata for a group of FASTQ files representing the different read components
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct IlmnFastqFileGroup {
    pub sample: String,
    pub s: usize,
    pub lane_mode: LaneMode,
    pub chunk: usize,
}

/// A parsed representation of an FASTQ file produced by
/// Illumina's bcl2fastq tool.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct IlmnFastqFile {
    pub group: IlmnFastqFileGroup,
    pub read: String,
    pub path: PathBuf,
}

/// Parse the ILMN fastq filename to get the read name, lane, read group,
/// and S field. We expect a filename of the form
/// <path>/<prefix>_S0_L001_R1_001.fastq
impl IlmnFastqFile {
    /// Attempt to parse `path` as an Illumina bcl2fastq-produced
    /// FASTQ file.
    pub fn new(path: impl AsRef<Path>) -> Option<IlmnFastqFile> {
        let filename = path.as_ref().file_name()?.to_str();

        if let Some(f) = filename {
            if let Some(cap) = BCL2FASTQ_REGEX.captures(f) {
                let sample = cap.get(1).unwrap().as_str().to_string();
                let s: usize = cap.get(2).unwrap().as_str().parse().unwrap();
                let lane: usize = cap.get(3).unwrap().as_str().parse().unwrap();
                let read = cap.get(4).unwrap().as_str().to_string();
                let chunk: usize = cap.get(5).unwrap().as_str().parse().unwrap();

                let r = Some(IlmnFastqFile {
                    group: IlmnFastqFileGroup {
                        sample,
                        s,
                        lane_mode: LaneMode::SingleLane(lane),
                        chunk,
                    },
                    read,

                    path: path.as_ref().into(),
                });

                return r;
            }

            // Try out the no lane split version next
            if let Some(cap) = BCL2FASTQ_NO_LANE_SPLIT_REGEX.captures(f) {
                let sample = cap.get(1).unwrap().as_str().to_string();
                let s: usize = cap.get(2).unwrap().as_str().parse().unwrap();
                let read = cap.get(3).unwrap().as_str().to_string();
                let chunk: usize = cap.get(4).unwrap().as_str().parse().unwrap();

                let r = Some(IlmnFastqFile {
                    group: IlmnFastqFileGroup {
                        sample,
                        s,
                        lane_mode: LaneMode::NoLaneSplitting,
                        chunk,
                    },
                    read,
                    path: path.as_ref().into(),
                });

                return r;
            }
        }

        None
    }
}

/// Find all the bcl2fastq FASTQ files present in `path`.
fn get_bcl2fastq_files(path: impl AsRef<Path>) -> Result<Vec<(IlmnFastqFile, PathBuf)>> {
    let path = path.as_ref();
    std::fs::read_dir(path)
        .with_context(|| path.display().to_string())?
        .map_ok(|f| f.path())
        .filter_map_ok(|f| IlmnFastqFile::new(&f).map(|x| (x, f)))
        .try_collect()
        .with_context(|| path.display().to_string())
}

/// Find all the sets of bcl2fastq FASTQ files present in `path` as well as directories directly
/// underneath `path`. Corresponding R1/R2/I1/I2 files are grouped together and reported in an
/// `InputFastqs`, along with a representative `IlmnFastqFile` struct.
pub fn find_flowcell_fastqs(
    path: impl AsRef<Path>,
) -> Result<Vec<(IlmnFastqFileGroup, InputFastqs)>> {
    let mut res = Vec::new();

    let mut files = get_bcl2fastq_files(&path)?;
    // Collect the files which are within the directories underneath `path`.
    // This typically means `path` corresponds to the project folder in the `mkfastq` outs
    for entry in std::fs::read_dir(&path)? {
        let entry = entry?.path();
        if entry.is_dir() {
            files.extend(get_bcl2fastq_files(entry)?);
        }
    }
    files.sort();

    for files in files.chunk_by(|x, y| x.0.group == y.0.group) {
        let mut my_files: HashMap<_, _> = files
            .iter()
            .map(|(info, path)| (info.read.as_str(), path.display().to_string()))
            .collect();

        let fastqs = match my_files.remove("R1") {
            Some(r1) => match my_files.remove("R3") {
                Some(r3) =>
                // Account for the confusing R1/R2/R3 convention that is used in SC3pv1, ATAC etc
                {
                    InputFastqs {
                        r1,
                        r2: Some(r3),
                        i1: my_files.remove("R2"),
                        i2: my_files.remove("I1"),
                        r1_interleaved: false,
                    }
                }
                None => InputFastqs {
                    r1,
                    r2: my_files.remove("R2"),
                    i1: my_files.remove("I1"),
                    i2: my_files.remove("I2"),
                    r1_interleaved: false,
                },
            },
            None => {
                // We will tolerate a missing R1 file here -- this
                // could be due to a copying error & unrelated to the
                // sample we're interested in, so we don't want to
                // throw an error just yet.
                continue;
            }
        };

        res.push((files[0].0.group.clone(), fastqs));
    }

    res.sort();
    Ok(res)
}

#[cfg(test)]
mod test {
    use super::*;
    use pretty_assertions::assert_eq;

    #[test]
    fn test_parse() {
        let filename = "heart_1k_v3_S1_L002_R2_001.fastq.gz";
        let r = IlmnFastqFile::new(filename);

        let expected = IlmnFastqFile {
            path: PathBuf::from(filename.to_string()),
            read: "R2".to_string(),
            group: IlmnFastqFileGroup {
                sample: "heart_1k_v3".to_string(),
                s: 1,
                lane_mode: 2.into(),
                chunk: 1,
            },
        };

        assert_eq!(r.unwrap(), expected);
    }

    #[test]
    fn test_parse_hyphen() {
        let filename = "heart-1k-v3_S1_L002_R2_001.fastq.gz";
        let r = IlmnFastqFile::new(filename);

        let expected = IlmnFastqFile {
            path: PathBuf::from(filename.to_string()),
            read: "R2".to_string(),
            group: IlmnFastqFileGroup {
                sample: "heart-1k-v3".to_string(),
                s: 1,
                lane_mode: 2.into(),
                chunk: 1,
            },
        };

        assert_eq!(r.unwrap(), expected);
    }

    #[test]
    fn test_parse_r3() {
        let filename = "heart_1k_v3_S1_L002_R3_001.fastq.gz";
        let r = IlmnFastqFile::new(filename);

        let expected = IlmnFastqFile {
            path: PathBuf::from(filename.to_string()),
            read: "R3".to_string(),
            group: IlmnFastqFileGroup {
                sample: "heart_1k_v3".to_string(),
                s: 1,
                lane_mode: 2.into(),
                chunk: 1,
            },
        };

        assert_eq!(r.unwrap(), expected);
    }

    #[test]
    fn test_parse_no_lane_split() {
        let filename = "tests/filenames/test_sample_S1_R1_001.fastq.gz";

        let r = IlmnFastqFile::new(filename);

        let expected = IlmnFastqFile {
            path: PathBuf::from(filename),
            read: "R1".to_string(),
            group: IlmnFastqFileGroup {
                sample: "test_sample".to_string(),
                s: 1,
                lane_mode: LaneMode::NoLaneSplitting,
                chunk: 1,
            },
        };

        assert_eq!(r.unwrap(), expected);
    }

    #[test]
    fn test_bad() {
        let filename = "heart_1k_v3_S1_LA_R2_001.fastq.gz";
        let r = IlmnFastqFile::new(filename);
        assert!(r.is_none());

        let filename = "heart_1k_v3_S1_L002_XX_001.fastq.gz";
        let r = IlmnFastqFile::new(filename);
        assert!(r.is_none());
    }

    #[test]
    fn query_bcl2fastq() -> Result<()> {
        let path = "tests/filenames/bcl2fastq";

        let query = Bcl2FastqDef {
            fastq_path: path.to_string(),
            sample_name_spec: "Infected".into(),
            lane_spec: LaneSpec::Any,
        };

        let fqs = query.find_fastqs()?;
        assert_eq!(fqs.len(), 2);
        assert_eq!(
            fqs[0].r1,
            "tests/filenames/bcl2fastq/Infected_S3_L001_R1_001.fastq"
        );
        assert_eq!(
            fqs[1].r1,
            "tests/filenames/bcl2fastq/Infected_S3_L002_R1_001.fastq"
        );
        Ok(())
    }

    #[test]
    fn query_bcl2fastq_lz4() -> Result<()> {
        let path = "tests/filenames/bcl2fastq_lz4";

        let query = Bcl2FastqDef {
            fastq_path: path.to_string(),
            sample_name_spec: "Infected".into(),
            lane_spec: LaneSpec::Any,
        };

        let fqs = query.find_fastqs()?;
        assert_eq!(fqs.len(), 2);
        assert_eq!(
            fqs[0].r1,
            "tests/filenames/bcl2fastq_lz4/Infected_S3_L001_R1_001.fastq.lz4"
        );
        assert_eq!(
            fqs[1].r1,
            "tests/filenames/bcl2fastq_lz4/Infected_S3_L002_R1_001.fastq.lz4"
        );
        Ok(())
    }

    #[test]
    fn query_bcl2fastq_lanes() -> Result<()> {
        let path = "tests/filenames/bcl2fastq";

        let query = Bcl2FastqDef {
            fastq_path: path.to_string(),
            sample_name_spec: "Infected".into(),
            lane_spec: LaneSpec::Lanes(vec![2].into_iter().collect()),
        };

        let fqs = query.find_fastqs()?;
        assert_eq!(fqs.len(), 1);
        assert_eq!(
            fqs[0].r1,
            "tests/filenames/bcl2fastq/Infected_S3_L002_R1_001.fastq"
        );
        Ok(())
    }

    #[test]
    fn test_bcl2fastq_no_lane_split() -> Result<()> {
        let path = "tests/filenames/bcl2fastq_no_lane_split";

        let query = Bcl2FastqDef {
            fastq_path: path.to_string(),
            sample_name_spec: "test_sample".into(),
            lane_spec: LaneSpec::Any,
        };

        let fqs = query.find_fastqs()?;
        // NOTE: R3 is not used
        let expected = vec![InputFastqs {
            r1: format!("{path}/test_sample_S1_R1_001.fastq.gz"),
            r2: Some(format!("{path}/test_sample_S1_R2_001.fastq.gz")),
            i1: Some(format!("{path}/test_sample_S1_I1_001.fastq.gz")),
            i2: None,
            r1_interleaved: false,
        }];
        assert_eq!(fqs, expected);
        Ok(())
    }
}

// The following tests are based on the tests in tenkit: lib/python/tenkit/test/test_fasta.py
#[cfg(test)]
mod tests_from_tenkit {
    use super::*;
    use pretty_assertions::assert_eq;

    #[test]
    fn test_find_input_fastq_files_bcl2fastq_demult() -> Result<()> {
        let path = "tests/filenames/bcl2fastq_2";

        let query = Bcl2FastqDef {
            fastq_path: path.to_string(),
            sample_name_spec: "test_sample".into(),
            lane_spec: LaneSpec::Any,
        };

        let fqs = query.find_fastqs()?;
        let expected = vec![InputFastqs {
            r1: format!("{path}/test_sample_S1_L001_R1_001.fastq.gz"),
            r2: Some(format!("{path}/test_sample_S1_L001_R3_001.fastq.gz")),
            i1: Some(format!("{path}/test_sample_S1_L001_R2_001.fastq.gz")),
            i2: Some(format!("{path}/test_sample_S1_L001_I1_001.fastq.gz")),
            r1_interleaved: false,
        }];
        assert_eq!(fqs, expected);
        Ok(())
    }

    #[test]
    fn test_find_input_fastq_files_bc2fastq_demult_project_1() -> Result<()> {
        let path = "tests/filenames/project_dir";

        let query = Bcl2FastqDef {
            fastq_path: path.to_string(),
            sample_name_spec: "test_sample".into(),
            lane_spec: LaneSpec::Any,
        };

        let fqs = query.find_fastqs()?;
        let expected = vec![
            InputFastqs {
                r1: format!("{path}/s1/test_sample_S1_L001_R1_001.fastq.gz"),
                r2: Some(format!("{path}/s1/test_sample_S1_L001_R3_001.fastq.gz")),
                i1: Some(format!("{path}/s1/test_sample_S1_L001_R2_001.fastq.gz")),
                i2: Some(format!("{path}/s1/test_sample_S1_L001_I1_001.fastq.gz")),
                r1_interleaved: false,
            },
            InputFastqs {
                r1: format!("{path}/s2/test_sample_S2_L001_R1_001.fastq.gz"),
                r2: Some(format!("{path}/s2/test_sample_S2_L001_R3_001.fastq.gz")),
                i1: Some(format!("{path}/s2/test_sample_S2_L001_R2_001.fastq.gz")),
                i2: Some(format!("{path}/s2/test_sample_S2_L001_I1_001.fastq.gz")),
                r1_interleaved: false,
            },
        ];
        assert_eq!(fqs, expected);
        Ok(())
    }

    #[test]
    fn test_find_input_fastq_files_bc2fastq_demult_project_2() -> Result<()> {
        let path = "tests/filenames/project_dir";

        let query = Bcl2FastqDef {
            fastq_path: path.to_string(),
            sample_name_spec: "test_sample2".into(),
            lane_spec: LaneSpec::Any,
        };

        let fqs = query.find_fastqs()?;
        let expected = vec![InputFastqs {
            r1: format!("{path}/s1/test_sample2_S1_L001_R1_001.fastq.gz"),
            r2: Some(format!("{path}/s1/test_sample2_S1_L001_R3_001.fastq.gz")),
            i1: Some(format!("{path}/s1/test_sample2_S1_L001_R2_001.fastq.gz")),
            i2: Some(format!("{path}/s1/test_sample2_S1_L001_I1_001.fastq.gz")),
            r1_interleaved: false,
        }];
        assert_eq!(fqs, expected);
        Ok(())
    }

    #[test]
    fn test_sample_name_verification() -> Result<()> {
        let path = "tests/filenames/tenkit91";
        for &s in &["test_sample", "test_sample_suffix"] {
            let query = Bcl2FastqDef {
                fastq_path: path.to_string(),
                sample_name_spec: s.into(),
                lane_spec: LaneSpec::Any,
            };
            let fqs = query.find_fastqs()?;
            let expected = vec![InputFastqs {
                r1: format!("{path}/{s}_S1_L001_R1_001.fastq.gz"),
                r2: Some(format!("{path}/{s}_S1_L001_R3_001.fastq.gz")),
                i1: Some(format!("{path}/{s}_S1_L001_R2_001.fastq.gz")),
                i2: Some(format!("{path}/{s}_S1_L001_I1_001.fastq.gz")),
                r1_interleaved: false,
            }];
            assert_eq!(fqs, expected);
        }
        Ok(())
    }

    #[test]
    fn test_sample_name_any() -> Result<()> {
        let path = "tests/filenames/tenkit91";

        let query = Bcl2FastqDef {
            fastq_path: path.to_string(),
            sample_name_spec: SampleNameSpec::Any,
            lane_spec: LaneSpec::Any,
        };
        let fqs = query.find_fastqs()?;
        let expected = vec![
            InputFastqs {
                r1: format!("{path}/test_sample_S1_L001_R1_001.fastq.gz"),
                r2: Some(format!("{path}/test_sample_S1_L001_R3_001.fastq.gz")),
                i1: Some(format!("{path}/test_sample_S1_L001_R2_001.fastq.gz")),
                i2: Some(format!("{path}/test_sample_S1_L001_I1_001.fastq.gz")),
                r1_interleaved: false,
            },
            InputFastqs {
                r1: format!("{path}/test_sample_suffix_S1_L001_R1_001.fastq.gz"),
                r2: Some(format!("{path}/test_sample_suffix_S1_L001_R3_001.fastq.gz")),
                i1: Some(format!("{path}/test_sample_suffix_S1_L001_R2_001.fastq.gz")),
                i2: Some(format!("{path}/test_sample_suffix_S1_L001_I1_001.fastq.gz")),
                r1_interleaved: false,
            },
        ];
        assert_eq!(fqs, expected);
        Ok(())
    }
}
