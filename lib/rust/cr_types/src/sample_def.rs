//!
//! There is code in fastq_set::sample_def, but having it in Cellranger makes
//! more sense

use crate::serde_helpers::NumberOrStr;
use crate::LibraryType;
use anyhow::{bail, ensure, Result};
use fastq_set::filenames::bcl2fastq::SampleNameSpec;
use fastq_set::filenames::bcl_processor::SampleIndexSpec;
use fastq_set::filenames::fastq_dir::{BclProcessorDir, FastqChecker};
use fastq_set::filenames::{FastqDef, FindFastqs, LaneSpec};
use fastq_set::read_pair_iter::InputFastqs;
use fastq_set::sample_index_map::SAMPLE_INDEX_MAP;
use fastq_set::SSeq;
use itertools::Itertools;
use martian_derive::MartianType;
use serde::{Deserialize, Deserializer, Serialize};
use std::collections::HashSet;
use std::convert::TryFrom;
use std::path::PathBuf;

pub const COUNT_HELP: &str = r#"No input FASTQs were found for the requested parameters.

If your files came from bcl2fastq or mkfastq:
 - Make sure you are specifying the correct --sample(s), i.e. matching the sample sheet
 - Make sure your files follow the correct naming convention, e.g. SampleName_S1_L001_R1_001.fastq.gz (and the R2 version)
 - Make sure your --fastqs points to the correct location.
 - Make sure your --lanes, if any, are correctly specified.

Refer to the "Specifying Input FASTQs" page at https://support.10xgenomics.com/ for more details.

"#;

#[allow(non_camel_case_types)]
#[derive(
    Serialize, Deserialize, Clone, Copy, PartialOrd, Ord, PartialEq, Eq, Hash, Debug, MartianType,
)]
pub enum FastqMode {
    BCL_PROCESSOR,
    ILMN_BCL2FASTQ,
}

pub fn parse_lanes<'de, D>(deserializer: D) -> Result<Option<Vec<usize>>, D::Error>
where
    D: Deserializer<'de>,
{
    let raw: Option<Vec<NumberOrStr>> = Option::deserialize(deserializer)?;
    match raw {
        Some(values) => Ok(Some(
            values
                .into_iter()
                .map(|v| usize::try_from(v).map_err(serde::de::Error::custom))
                .try_collect()?,
        )),
        None => Ok(None),
    }
}

/// Pointer to one logical set of FASTQ from a unique (library, gem_group) tuple
#[derive(Serialize, Deserialize, Clone, PartialOrd, PartialEq, Debug, MartianType)]
pub struct SampleDef {
    pub fastq_mode: FastqMode,
    pub gem_group: Option<u16>,
    #[serde(deserialize_with = "parse_lanes")]
    pub lanes: Option<Vec<usize>>,
    pub library_type: Option<LibraryType>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub r1_length: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub r2_length: Option<usize>,
    pub read_path: PathBuf,
    pub sample_indices: Option<Vec<String>>,
    pub sample_names: Option<Vec<String>>,
    pub subsample_rate: Option<f64>,
    pub fastq_id: Option<String>,
}

impl Default for SampleDef {
    fn default() -> Self {
        SampleDef {
            fastq_mode: FastqMode::BCL_PROCESSOR,
            gem_group: None,
            lanes: None,
            library_type: None,
            r1_length: None,
            r2_length: None,
            read_path: PathBuf::new(),
            sample_indices: None,
            sample_names: None,
            subsample_rate: None,
            fastq_id: None,
        }
    }
}

impl SampleDef {
    // Check for missing FASTQ files.
    #[allow(unused)]
    pub fn check_fastqs(&self, help_text: &str) -> Result<()> {
        let read_path = self.read_path.display();
        match self.fastq_mode {
            FastqMode::BCL_PROCESSOR => {
                let dir = BclProcessorDir::new(&self.read_path)?;
                if dir.is_empty() {
                    bail!(
                        "The directory {read_path} contains no 10X bcl processor compatible fastqs."
                    );
                }

                let lane_spec = self.lane_spec();

                // Check for missing sample indices.
                match self.sample_indices.as_deref() {
                    None | Some([]) => {
                        bail!(
                            "Expected non empty 'sample_indices' entry in sample_def when fastq_mode is BCL_PROCESSOR"
                        );
                    }
                    Some([sample_index]) if sample_index == "any" => {
                        // Any sample index is allowed.
                    }
                    Some(sample_indices) => {
                        let missing_sample_indices = sample_indices
                            .iter()
                            .filter(|si| !dir.contains_index(si))
                            .join(",");

                        ensure!(
                            missing_sample_indices.is_empty(),
                            "No input FASTQs were found for the sample indices \
                                {missing_sample_indices} in the path: {read_path}"
                        );

                        if let LaneSpec::Lanes(lanes) = &lane_spec {
                            for sample_index in sample_indices {
                                ensure!(
                                    lanes.iter().any(|lane| dir.contains_index_with_lane(sample_index, *lane)),
                                    "No input FASTQs were found for the sample index {sample_index} \
                                    and lanes {:?} in the path: {read_path}",
                                    lanes.iter().sorted().join(",")
                                );
                            }
                        }
                    }
                }

                // Check for missing lanes.
                if let LaneSpec::Lanes(lanes) = lane_spec {
                    let missing_lanes = lanes
                        .into_iter()
                        .filter(|&lane| !dir.contains_lane(lane))
                        .join(",");
                    ensure!(
                        missing_lanes.is_empty(),
                        "Lane(s) {missing_lanes} does not exist in the path: {read_path}"
                    );
                }
            }

            FastqMode::ILMN_BCL2FASTQ => {
                if matches!(self.sample_names.as_deref(), None | Some([])) {
                    bail!(
                        "Expected non empty 'sample_names' entry in sample_def when fastq_mode is ILMN_BCL2FASTQ"
                    );
                }
                let samples_found = FastqChecker::bcl2fastq_check_and_infer_sample_names(
                    &self.read_path,
                    dbg!(&self.parse_sample_names()),
                    &self.lanes,
                    help_text,
                )?;
                match self.sample_name_spec() {
                    SampleNameSpec::Any => {}
                    SampleNameSpec::Names(names) => {
                        let missing = names.difference(&samples_found).join("\n");
                        if !missing.is_empty() {
                            bail!(
                                "The following sample(s) are not found in path {read_path}:\n{missing}"
                            );
                        }
                    }
                }
            }
        }

        Ok(())
    }

    /// Sample Indices
    fn sample_index_spec(&self) -> SampleIndexSpec {
        let mut samp_indices = HashSet::new();
        let raw_indices = self.sample_indices.as_ref().unwrap();
        for index in raw_indices {
            if index == "any" {
                assert!(raw_indices.len() == 1);
                return SampleIndexSpec::Any;
            } else if let Some(seqs) = SAMPLE_INDEX_MAP.get(index.as_str()) {
                samp_indices.extend(seqs.iter().map(|s| SSeq::from_bytes(s.as_bytes())));
            } else if index.chars().all(|c| "ACGT".contains(c)) {
                samp_indices.insert(SSeq::from_bytes(index.as_bytes()));
            } else {
                panic!(
                    "Sample index '{index}' is not valid. Must be one of: any, SI-<number>, \
                     SI-<plate>-<well-coordinate>, 220<part-number>, or \
                     a nucleotide sequence."
                );
            }
        }
        SampleIndexSpec::Sequences {
            indices: samp_indices,
            max_n: 1,
        }
    }

    fn lane_spec(&self) -> LaneSpec {
        match &self.lanes {
            None => LaneSpec::Any,
            Some(lanes) => LaneSpec::Lanes(lanes.iter().copied().collect()),
        }
    }

    fn parse_sample_names(&self) -> Option<Vec<String>> {
        match self.sample_names {
            Some(ref names) => {
                if names.len() == 1 && names[0] == "any" {
                    None
                } else {
                    Some(names.clone())
                }
            }
            None => None,
        }
    }

    fn sample_name_spec(&self) -> SampleNameSpec {
        assert!(self.sample_names.is_some());
        match self.parse_sample_names() {
            Some(names) => SampleNameSpec::Names(names.into_iter().collect()),
            None => SampleNameSpec::Any,
        }
    }

    pub fn per_group_fastq_def(&self) -> Result<Vec<(Option<String>, FastqDef)>> {
        self.check_fastqs(COUNT_HELP)?;
        let mut result = Vec::new();
        match self.fastq_mode {
            FastqMode::BCL_PROCESSOR => {
                let raw_indices = self.sample_indices.as_ref().unwrap();
                for index in raw_indices {
                    if index == "any" {
                        assert!(raw_indices.len() == 1);
                        result.push((
                            None,
                            FastqDef::bcl_processor(
                                self.read_path.to_str().unwrap().to_owned(),
                                SampleIndexSpec::Any,
                                self.lane_spec(),
                            ),
                        ));
                    } else if let Some(seqs) = SAMPLE_INDEX_MAP.get(index.as_str()) {
                        result.push((
                            Some(index.to_string()),
                            FastqDef::bcl_processor(
                                self.read_path.to_str().unwrap().to_owned(),
                                SampleIndexSpec::Sequences {
                                    indices: seqs
                                        .iter()
                                        .map(|s| SSeq::from_bytes(s.as_bytes()))
                                        .collect(),
                                    max_n: 1,
                                },
                                self.lane_spec(),
                            ),
                        ));
                    } else if index.chars().all(|c| "ACGT".contains(c)) {
                        let mut indices = HashSet::new();
                        indices.insert(SSeq::from_bytes(index.as_bytes()));

                        result.push((
                            Some(index.to_string()),
                            FastqDef::bcl_processor(
                                self.read_path.to_str().unwrap().to_owned(),
                                SampleIndexSpec::Sequences { indices, max_n: 1 },
                                self.lane_spec(),
                            ),
                        ));
                    } else {
                        panic!(
                            "Sample index '{index}' is not valid. Must be one of: any, SI-<number>, \
                     SI-<plate>-<well-coordinate>, 220<part-number>, or \
                     a nucleotide sequence."
                        );
                    }
                }
            }
            FastqMode::ILMN_BCL2FASTQ => {
                let name_spec = self.sample_name_spec();
                match name_spec {
                    SampleNameSpec::Any => result.push((
                        None,
                        FastqDef::bcl2fastq(
                            self.read_path.to_str().unwrap().to_owned(),
                            SampleNameSpec::Any,
                            self.lane_spec(),
                        ),
                    )),
                    SampleNameSpec::Names(names) => {
                        for name in names {
                            result.push((
                                Some(name.clone()),
                                FastqDef::bcl2fastq(
                                    self.read_path.to_str().unwrap().to_owned(),
                                    name.as_str().into(),
                                    self.lane_spec(),
                                ),
                            ));
                        }
                    }
                }
            }
        }
        Ok(result)
    }

    /// Get a FastqDef corresponding to this SampleDef.
    pub fn get_fastq_def(&self) -> Result<FastqDef> {
        match self.fastq_mode {
            FastqMode::BCL_PROCESSOR => {
                if self.sample_indices.is_none() {
                    bail!(
                        "Expected 'sample_indices' entry in sample_def when fastq_mode is BCL_PROCESSOR"
                    );
                }

                Ok(FastqDef::bcl_processor(
                    self.read_path.to_str().unwrap().to_owned(),
                    self.sample_index_spec(),
                    self.lane_spec(),
                ))
            }

            FastqMode::ILMN_BCL2FASTQ => Ok(FastqDef::bcl2fastq(
                self.read_path.to_str().unwrap().to_owned(),
                self.sample_name_spec(),
                self.lane_spec(),
            )),
        }
    }

    /// Return a Vec of InputFastqs.
    /// Check if there's any data here and give a good error msg if not.
    /// This error should have already been caught if detect chemistry ran.
    pub fn find_fastqs(&self) -> Result<Vec<InputFastqs>> {
        self.check_fastqs(COUNT_HELP)?;
        let fastqs = self.get_fastq_def()?.find_fastqs()?;
        if fastqs.is_empty() {
            bail!("No input FASTQs were found for the sample def: {:#?}", self);
        }
        Ok(fastqs)
    }
}

#[cfg(test)]
mod sample_def_tests {
    use super::*;

    const GEX_SD: &str = r#"{
        "fastq_mode": "BCL_PROCESSOR",
        "gem_group": 1,
        "lanes": [1, 2],
        "library_type": "Antibody Capture",
        "read_path": "demultiplexed_fastq_path",
        "sample_indices": ["SI-P01-D10"],
        "target_set": null,
        "target_set_name": null
    }"#;

    const ANTIBODY_SD: &str = r#"{
        "fastq_mode": "BCL_PROCESSOR",
        "gem_group": 1,
        "lanes": [1],
        "library_type": "Gene Expression",
        "read_path": "demultiplexed_fastq_path",
        "sample_indices": ["SI-P2-D12"],
        "target_set": null,
        "target_set_name": null
    }"#;

    const TINY_MULTI_CS_SD: &str = r#"{
        "fastq_mode": "ILMN_BCL2FASTQ",
        "gem_group": 1,
        "lanes": [
            2
        ],
        "library_type": "Gene Expression",
        "read_path": "350_gex_fastqs",
        "sample_indices": [
            "any"
        ],
        "sample_names": [
            "350-tiny-test-gex"
        ],
        "target_set": null,
        "target_set_name": null
    }"#;

    const GEX_CS_SD: &str = r#"{
        "gem_group": null,
        "lanes": [
            "1"
        ],
        "sample_names": [
            "bamtofastq"
        ],
        "sample_indices": [
            "any"
        ],
        "target_set_name": "null",
        "fastq_mode": "ILMN_BCL2FASTQ",
        "read_path": "CR-300-01_chr21",
        "target_set": null
    }"#;

    #[test]
    fn read_gex_sample_def() {
        println!("{GEX_SD}");
        let gex_sd: SampleDef = serde_json::from_str(GEX_SD).unwrap();
        println!("{gex_sd:?}");
    }

    #[test]
    fn read_abs_sample_def() {
        let abs_sd: SampleDef = serde_json::from_str(ANTIBODY_SD).unwrap();
        println!("{abs_sd:?}");
    }

    #[test]
    fn read_gex_cs_sample_def() {
        let sd: SampleDef = serde_json::from_str(GEX_CS_SD).unwrap();
        println!("{sd:?}");
    }

    #[test]
    fn read_tiny_multi_sample_def() {
        let sd: SampleDef = serde_json::from_str(TINY_MULTI_CS_SD).unwrap();
        println!("{sd:?}");
    }
}

#[cfg(test)]
mod invalid_bcl2fastq_sample_def {
    use super::{FastqMode, SampleDef, COUNT_HELP};
    use std::path::Path;

    fn test_invalid_bcl2fastq_input(
        fastq_path: impl AsRef<Path>,
        lanes: Option<Vec<usize>>,
        sample_names: Vec<&str>,
    ) {
        let sdef = SampleDef {
            fastq_mode: FastqMode::ILMN_BCL2FASTQ,
            lanes,
            read_path: fastq_path.as_ref().to_path_buf(),
            sample_names: Some(sample_names.into_iter().map(String::from).collect()),
            ..SampleDef::default()
        };
        let check = sdef.check_fastqs(COUNT_HELP);
        println!("{check:#?}");
        assert!(check.is_err());
    }

    #[test]
    fn test_missing_only_sample() {
        test_invalid_bcl2fastq_input(
            "../dui_tests/test_resources/cellranger-count/too_few_reads",
            None,
            vec!["missing_sample"],
        );
    }

    #[test]
    fn test_missing_lane() {
        test_invalid_bcl2fastq_input(
            "../dui_tests/test_resources/cellranger-count/too_few_reads",
            Some(vec![2]),
            vec!["test_sample"],
        );
    }

    #[test]
    fn test_missing_one_sample() {
        test_invalid_bcl2fastq_input(
            "../dui_tests/test_resources/cellranger-count/too_few_reads",
            None,
            vec!["test_sample", "missing_sample"],
        );
    }

    #[test]
    fn test_invalid_folder() {
        test_invalid_bcl2fastq_input("/no/such/folder", None, vec!["test_sample"]);
    }

    #[test]
    fn test_empty_folder() {
        test_invalid_bcl2fastq_input(
            "../dui_tests/test_resources/cellranger-count/fastqs_empty",
            None,
            vec!["test_sample"],
        );
    }

    #[test]
    fn test_empty_sample_name() {
        test_invalid_bcl2fastq_input(
            "../dui_tests/test_resources/cellranger-count/too_few_reads",
            None,
            vec![],
        );
    }
}

#[cfg(test)]
mod invalid_bcl_proc_sample_def {
    use super::{FastqMode, SampleDef, COUNT_HELP};
    use std::path::Path;

    fn test_invalid_bcl_proc_input(
        fastq_path: impl AsRef<Path>,
        lanes: Option<Vec<usize>>,
        sample_indices: Option<Vec<&str>>,
    ) {
        let sdef = SampleDef {
            fastq_mode: FastqMode::BCL_PROCESSOR,
            lanes,
            read_path: fastq_path.as_ref().to_path_buf(),
            sample_indices: sample_indices.map(|xs| xs.into_iter().map(String::from).collect()),
            ..SampleDef::default()
        };
        let check = sdef.check_fastqs(COUNT_HELP);
        println!("{check:#?}");
        assert!(check.is_err());
    }

    #[test]
    fn test_none_sample_indices() {
        test_invalid_bcl_proc_input("", None, None);
    }

    #[test]
    fn test_empty_sample_indices() {
        test_invalid_bcl_proc_input("", None, Some(vec![]));
    }

    #[test]
    fn test_invalid_folder() {
        test_invalid_bcl_proc_input("/no/such/folder", None, Some(vec!["test_sample"]));
    }

    #[test]
    fn test_empty_folder() {
        test_invalid_bcl_proc_input(
            "../dui_tests/test_resources/cellranger-count/fastqs_empty",
            None,
            Some(vec!["SI-P2-D3"]),
        );
    }

    #[test]
    fn test_missing_si() {
        test_invalid_bcl_proc_input(
            "../dui_tests/test_resources/cellranger-count/cycle_failure_bc_v3",
            None,
            Some(vec!["SI-P2-D2"]),
        );
    }

    #[test]
    fn test_missing_one_si() {
        test_invalid_bcl_proc_input(
            "../dui_tests/test_resources/cellranger-count/cycle_failure_bc_v3",
            None,
            Some(vec!["SI-P2-D3", "SI-P2-D2"]),
        );
    }

    #[test]
    fn test_missing_lane() {
        test_invalid_bcl_proc_input(
            "../dui_tests/test_resources/cellranger-count/cycle_failure_bc_v3",
            Some(vec![2]),
            Some(vec!["SI-P2-D3"]),
        );
    }

    #[test]
    fn test_index_lane_pair_missing() {
        // Both the index and the lane are present, but the pair is missing
        test_invalid_bcl_proc_input(
            "../dui_tests/test_resources/cellranger-count/index_lane_pair_missing",
            Some(vec![1]),
            Some(vec!["TATCAGCCTA"]),
        );
    }
}
