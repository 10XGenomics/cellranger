#![expect(missing_docs)]
use crate::rna_read::AdapterTrimmingParameters;
use crate::{AlignerParam, LibraryType};
use anyhow::{Ok, Result, bail, ensure};
use fastq_set::SSeq;
use fastq_set::filenames::bcl_processor::SampleIndexSpec;
use fastq_set::filenames::bcl2fastq::{self, SampleNameSpec};
use fastq_set::filenames::fastq_dir::{BclProcessorDir, FastqChecker};
use fastq_set::filenames::{LaneSpec, bcl_processor, sbx};
use fastq_set::read_pair_iter::FastqSet;
use fastq_set::sample_index_map::SAMPLE_INDEX_MAP;
use itertools::Itertools;
use martian_derive::{MartianStruct, MartianType};
use serde::{Deserialize, Serialize};
use std::collections::HashSet;
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
    ONT,
    PACBIO,
    SBX,
}

impl FastqMode {
    /// Return the adapter trimming parameters required for this mode.
    pub fn adapter_trimming_parameters(
        self,
        aligner: Option<AlignerParam>,
    ) -> Option<AdapterTrimmingParameters> {
        match self {
            FastqMode::BCL_PROCESSOR | FastqMode::ILMN_BCL2FASTQ => None,
            FastqMode::SBX => Some(AdapterTrimmingParameters {
                is_directional: true,
                aligner: aligner.unwrap_or(AlignerParam::Minimap2),
            }),
            FastqMode::ONT | FastqMode::PACBIO => Some(AdapterTrimmingParameters {
                is_directional: false,
                aligner: aligner.unwrap_or(AlignerParam::Minimap2),
            }),
        }
    }
}

/// Pointer to one logical set of sequencing files from a unique (library, gem_group) tuple
#[derive(Serialize, Deserialize, Clone, PartialOrd, PartialEq, Debug, MartianStruct)]
pub struct SampleDef {
    pub fastq_mode: FastqMode,
    pub gem_group: Option<u16>,
    pub lanes: Option<Vec<usize>>,
    pub library_type: Option<LibraryType>,
    pub r1_length: Option<usize>,
    pub r2_length: Option<usize>,
    // Path to raw sequencing files in any allowed format
    pub read_path: PathBuf,
    pub sample_indices: Option<Vec<String>>,
    pub sample_names: Option<Vec<String>>,
    pub subsample_rate: Option<f64>,
    pub fastq_id: Option<String>,
}

impl SampleDef {
    /// Check for missing input sequencing data files.
    pub fn check_input_reads(&self, help_text: &str) -> Result<()> {
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
                                    lanes
                                        .iter()
                                        .any(|lane| dir
                                            .contains_index_with_lane(sample_index, *lane)),
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
            FastqMode::SBX => {
                if sbx::find_fastqs(&self.read_path)?.is_empty() {
                    bail!("The directory {read_path} contains no FASTQs.");
                }
            }
            FastqMode::ONT | FastqMode::PACBIO => {
                todo!()
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

    // Validates input reads and constructs FASTQ definitions if FASTQ paths are present.
    // Used in DETECT_CHEMISTRY
    pub fn per_group_fastq_def(&self) -> Result<Vec<(Option<String>, Vec<FastqSet>)>> {
        self.check_input_reads(COUNT_HELP)?;
        let mut result = Vec::new();
        match self.fastq_mode {
            FastqMode::BCL_PROCESSOR => {
                let raw_indices = self.sample_indices.as_ref().unwrap();
                for index in raw_indices {
                    if index == "any" {
                        assert!(raw_indices.len() == 1);
                        result.push((
                            None,
                            bcl_processor::find_fastqs(
                                &self.read_path,
                                &SampleIndexSpec::Any,
                                &self.lane_spec(),
                            )?,
                        ));
                    } else if let Some(seqs) = SAMPLE_INDEX_MAP.get(index.as_str()) {
                        result.push((
                            Some(index.clone()),
                            bcl_processor::find_fastqs(
                                &self.read_path,
                                &SampleIndexSpec::Sequences {
                                    indices: seqs
                                        .iter()
                                        .map(|s| SSeq::from_bytes(s.as_bytes()))
                                        .collect(),
                                    max_n: 1,
                                },
                                &self.lane_spec(),
                            )?,
                        ));
                    } else if index.chars().all(|c| "ACGT".contains(c)) {
                        let mut indices = HashSet::new();
                        indices.insert(SSeq::from_bytes(index.as_bytes()));

                        result.push((
                            Some(index.clone()),
                            bcl_processor::find_fastqs(
                                &self.read_path,
                                &SampleIndexSpec::Sequences { indices, max_n: 1 },
                                &self.lane_spec(),
                            )?,
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
                        bcl2fastq::find_fastqs(
                            &self.read_path,
                            &SampleNameSpec::Any,
                            &self.lane_spec(),
                        )?,
                    )),
                    SampleNameSpec::Names(names) => {
                        for name in names {
                            result.push((
                                Some(name.clone()),
                                bcl2fastq::find_fastqs(
                                    &self.read_path,
                                    &name.as_str().into(),
                                    &self.lane_spec(),
                                )?,
                            ));
                        }
                    }
                }
            }
            FastqMode::SBX => {
                // TODO: how should we break these up?
                result.push((None, sbx::find_fastqs(&self.read_path)?));
            }
            FastqMode::ONT | FastqMode::PACBIO => todo!(),
        }
        Ok(result)
    }

    /// Get sequencing units from this sample.
    pub fn find_fastqs(&self) -> Result<Vec<FastqSet>> {
        let fastqs = match self.fastq_mode {
            FastqMode::BCL_PROCESSOR => {
                if self.sample_indices.is_none() {
                    bail!(
                        "Expected 'sample_indices' entry in sample_def when fastq_mode is BCL_PROCESSOR"
                    );
                }

                bcl_processor::find_fastqs(
                    &self.read_path,
                    &self.sample_index_spec(),
                    &self.lane_spec(),
                )?
            }
            FastqMode::ILMN_BCL2FASTQ => bcl2fastq::find_fastqs(
                &self.read_path,
                &self.sample_name_spec(),
                &self.lane_spec(),
            )?,
            FastqMode::SBX => sbx::find_fastqs(&self.read_path)?,
            FastqMode::ONT | FastqMode::PACBIO => {
                todo!()
            }
        };
        if fastqs.is_empty() {
            bail!("No input FASTQs were found for the sample def: {self:#?}");
        }
        Ok(fastqs)
    }
}

pub struct SampleDefBuilder {
    pub fastq_mode: FastqMode,
    pub gem_group: Option<u16>,
    pub lanes: Option<Vec<usize>>,
    pub library_type: Option<LibraryType>,
    pub r1_length: Option<usize>,
    pub r2_length: Option<usize>,
    pub read_path: PathBuf,
    pub sample_indices: Option<Vec<String>>,
    pub sample_names: Option<Vec<String>>,
    pub subsample_rate: Option<f64>,
    pub fastq_id: Option<String>,
}

impl SampleDefBuilder {
    pub fn new(fastq_mode: FastqMode, read_path: PathBuf) -> Self {
        Self {
            fastq_mode,
            gem_group: None,
            lanes: None,
            library_type: None,
            r1_length: None,
            r2_length: None,
            read_path,
            sample_indices: None,
            sample_names: None,
            subsample_rate: None,
            fastq_id: None,
        }
    }

    pub fn lanes(mut self, lanes: Vec<usize>) -> Self {
        self.lanes = Some(lanes);
        self
    }

    pub fn library_type(mut self, library_type: LibraryType) -> Self {
        self.library_type = Some(library_type);
        self
    }

    pub fn r1_length(mut self, length: usize) -> Self {
        self.r1_length = Some(length);
        self
    }

    pub fn r2_length(mut self, length: usize) -> Self {
        self.r2_length = Some(length);
        self
    }

    pub fn sample_indices<T: AsRef<str>>(mut self, indices: Vec<T>) -> Self {
        self.sample_indices = Some(
            indices
                .into_iter()
                .map(|s| s.as_ref().to_string())
                .collect(),
        );
        self
    }

    pub fn sample_names<T: AsRef<str>>(mut self, names: Vec<T>) -> Self {
        self.sample_names = Some(names.into_iter().map(|s| s.as_ref().to_string()).collect());
        self
    }

    pub fn build(self) -> SampleDef {
        match self.fastq_mode {
            FastqMode::BCL_PROCESSOR | FastqMode::ILMN_BCL2FASTQ => SampleDef {
                fastq_mode: self.fastq_mode,
                gem_group: self.gem_group,
                lanes: self.lanes,
                library_type: self.library_type,
                r1_length: self.r1_length,
                r2_length: self.r2_length,
                read_path: self.read_path.clone(),
                sample_indices: self.sample_indices,
                sample_names: self.sample_names,
                subsample_rate: self.subsample_rate,
                fastq_id: self.fastq_id,
            },
            FastqMode::ONT | FastqMode::PACBIO | FastqMode::SBX => SampleDef {
                fastq_mode: self.fastq_mode,
                gem_group: self.gem_group,
                lanes: self.lanes,
                library_type: self.library_type,
                r1_length: self.r1_length,
                r2_length: self.r2_length,
                read_path: self.read_path,
                sample_indices: self.sample_indices,
                sample_names: self.sample_names,
                subsample_rate: self.subsample_rate,
                fastq_id: self.fastq_id,
            },
        }
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
        "fastq_path": "demultiplexed_fastq_path",
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
        "fastq_path": "demultiplexed_fastq_path",
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
        "fastq_path": "350_gex_fastqs",
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
            1
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
        "fastq_path": "CR-300-01_chr21",
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
    use super::{COUNT_HELP, FastqMode, SampleDefBuilder};
    use std::path::Path;

    fn test_invalid_bcl2fastq_input(
        fastq_path: impl AsRef<Path>,
        lanes: Option<Vec<usize>>,
        sample_names: Vec<&str>,
    ) {
        let mut builder =
            SampleDefBuilder::new(FastqMode::ILMN_BCL2FASTQ, fastq_path.as_ref().to_path_buf())
                .sample_names(sample_names);
        if let Some(lanes) = lanes {
            builder = builder.lanes(lanes);
        }
        let sdef = builder.build();
        let check = sdef.check_input_reads(COUNT_HELP);
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
    use super::{COUNT_HELP, FastqMode, SampleDefBuilder};
    use std::path::Path;

    fn test_invalid_bcl_proc_input(
        fastq_path: impl AsRef<Path>,
        lanes: Option<Vec<usize>>,
        sample_indices: Option<Vec<&str>>,
    ) {
        let mut builder =
            SampleDefBuilder::new(FastqMode::BCL_PROCESSOR, fastq_path.as_ref().to_path_buf());
        if let Some(lanes) = lanes {
            builder = builder.lanes(lanes);
        };
        if let Some(sample_indices) = sample_indices {
            builder = builder.sample_indices(sample_indices);
        }
        let sdef = builder.build();
        let check = sdef.check_input_reads(COUNT_HELP);
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
