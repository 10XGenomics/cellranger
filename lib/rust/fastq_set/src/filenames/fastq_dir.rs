//! Scan a directory of FASTQ file & provide access a subset of files based on lane or sample name.
#![expect(missing_docs)]

use crate::filenames::bcl_processor::{self, BclProcessorFileGroup};
use crate::filenames::bcl2fastq::{self, IlmnFastqFileGroup, SampleNameSpec};
use crate::filenames::{LaneMode, LaneSpec};
use crate::read_pair_iter::InputFastqs;
use crate::sample_index_map::SAMPLE_INDEX_MAP;
use anyhow::{Context, Result, bail};
use itertools::Itertools;
use std::collections::HashSet;
use std::path::{Path, PathBuf};

/// TODO: This should be moved to rust-utils
pub struct Bcl2FastqDir {
    fastq_path: PathBuf,
    fastq_data: Vec<(IlmnFastqFileGroup, InputFastqs)>,
    is_lane_split: bool,
    samples: HashSet<String>,
}

impl Bcl2FastqDir {
    pub fn new(fastq_path: impl AsRef<Path>) -> Result<Self> {
        let fastq_path = PathBuf::from(fastq_path.as_ref());
        if !fastq_path.exists() {
            bail!("{} does not exist", fastq_path.display());
        }

        let fastq_data = bcl2fastq::find_flowcell_fastqs(&fastq_path)?;

        if fastq_data.is_empty() {
            // The caller chooses what to do if the directory has no fastq files.
            return Ok(Bcl2FastqDir {
                fastq_path,
                fastq_data,
                is_lane_split: true,
                samples: HashSet::default(),
            });
        }

        let Ok(is_lane_split) = fastq_data
            .iter()
            .map(|(g, _)| match g.lane_mode {
                LaneMode::NoLaneSplitting => false,
                LaneMode::SingleLane(_) => true,
            })
            .dedup()
            .exactly_one()
        else {
            bail!(
                "Some files in the fastq path {} are split by lane, while some are not. \
                 This is not supported.",
                fastq_path.display()
            );
        };

        let samples = fastq_data.iter().map(|(g, _)| g.sample.clone()).collect();

        Ok(Bcl2FastqDir {
            fastq_path,
            fastq_data,
            is_lane_split,
            samples,
        })
    }

    pub fn fastq_count(&self) -> usize {
        self.fastq_data.len()
    }

    pub fn is_empty(&self) -> bool {
        self.fastq_data.is_empty()
    }

    pub fn contains_lane(&self, lane: usize) -> bool {
        assert!(
            self.is_lane_split,
            "Fastq files in {} are generated without splitting by lane.",
            self.fastq_path.display()
        );
        self.fastq_data
            .iter()
            .any(|(group, _)| group.lane_mode == lane.into())
    }

    pub fn samples(&self) -> &HashSet<String> {
        &self.samples
    }

    pub fn fastq_data(&self) -> &[(IlmnFastqFileGroup, InputFastqs)] {
        &self.fastq_data
    }

    pub fn filtered_fastq_data<'a>(
        &'a self,
        sample_name_spec: &'a SampleNameSpec,
        lane_spec: &'a LaneSpec,
    ) -> impl Iterator<Item = &'a (IlmnFastqFileGroup, InputFastqs)> + 'a {
        self.fastq_data.iter().filter(move |(g, _)| {
            sample_name_spec.contains(&g.sample) && lane_spec.contains(g.lane_mode)
        })
    }
}

pub struct BclProcessorDir {
    #[expect(dead_code)]
    fastq_path: PathBuf,
    fastq_data: Vec<(BclProcessorFileGroup, InputFastqs)>,
    sample_indices: HashSet<String>,
    lanes: HashSet<usize>,
}

impl BclProcessorDir {
    pub fn new(fastq_path: impl AsRef<Path>) -> Result<Self> {
        let fastq_path = PathBuf::from(fastq_path.as_ref());
        if !fastq_path.exists() {
            bail!("{} does not exist", fastq_path.display());
        }

        let fastq_data = bcl_processor::find_flowcell_fastqs(&fastq_path)?;
        let (sample_indices, lanes) = fastq_data
            .iter()
            .map(|(x, _)| (x.si.clone(), x.lane))
            .unzip();
        Ok(BclProcessorDir {
            fastq_path,
            fastq_data,
            sample_indices,
            lanes,
        })
    }

    pub fn is_empty(&self) -> bool {
        self.fastq_data.is_empty()
    }

    pub fn contains_index(&self, index: &str) -> bool {
        if let Some(seqs) = SAMPLE_INDEX_MAP.get(index) {
            seqs.iter().any(|&seq| self.sample_indices.contains(seq))
        } else {
            self.sample_indices.contains(index)
        }
    }

    pub fn contains_lane(&self, lane: usize) -> bool {
        self.lanes.contains(&lane)
    }

    pub fn contains_index_with_lane(&self, index: &str, lane: usize) -> bool {
        if let Some(seqs) = SAMPLE_INDEX_MAP.get(index) {
            seqs.iter().any(|&seq| {
                self.fastq_data
                    .iter()
                    .any(|(x, _)| x.si == seq && x.lane == lane)
            })
        } else {
            self.fastq_data
                .iter()
                .any(|(x, _)| x.si == index && x.lane == lane)
        }
    }
}

const FQ_HELP: &str = r#"No input FASTQs were found for the requested parameters.

If your files came from bcl2fastq or mkfastq:
 - Make sure you are specifying the correct --sample(s), i.e. matching the sample sheet
 - Make sure your files follow the correct naming convention, e.g. SampleName_S1_L001_R1_001.fastq.gz (and the R2 version)
 - Make sure your --fastqs points to the correct location.

Refer to the "Specifying Input FASTQs" page at https://support.10xgenomics.com/ for more details.

"#;

pub struct FastqChecker;

impl FastqChecker {
    /// Find the samples in `fastqs` for a single 'project resolved' fastq path that are
    /// compatible with `requested_samples`.
    ///
    /// Throw an error if none of the requested samples are found, or if no sample request
    /// was made, but multiple samples are present.
    pub fn bcl2fastq_check_and_infer_sample_names(
        fastq_path: impl AsRef<Path>,
        requested_samples: &Option<Vec<String>>,
        lanes: &Option<Vec<usize>>,
        help_text: &str,
    ) -> Result<HashSet<String>> {
        let fastq_path = fastq_path.as_ref();
        let bcl_dir = Bcl2FastqDir::new(fastq_path).context("Reading FASTQ directory")?;
        if bcl_dir.is_empty() {
            bail!(help_text.to_string());
        }

        // >=1 samples due to the check above
        let available_samples = bcl_dir.samples();
        let available_samples_str = available_samples.iter().sorted().join("\n");

        let sample_name_spec = match requested_samples {
            None => {
                // there should be exactly one sample present
                if available_samples.len() > 1 {
                    bail!(
                        "The --sample argument must be specified if multiple samples were \
                         demultiplexed in a run folder.  Available samples:\n{available_samples_str}",
                    );
                }
                SampleNameSpec::Names(available_samples.clone())
            }
            Some(names) => {
                let available_and_requested: HashSet<_> = names
                    .iter()
                    .filter(|&name| available_samples.contains(name))
                    .cloned()
                    .collect();
                if available_and_requested.is_empty() {
                    // Need to find a least one of the requested samples in this folder
                    bail!(
                        "Requested sample(s) not found in fastq directory {fastq_path:?}\n\
                         Available samples:\n{available_samples_str}",
                    );
                }
                SampleNameSpec::Names(available_and_requested)
            }
        };

        let lane_spec = match lanes {
            None => LaneSpec::Any,
            Some(l) => LaneSpec::Lanes(l.iter().copied().collect()),
        };

        if bcl_dir
            .filtered_fastq_data(&sample_name_spec, &lane_spec)
            .next()
            .is_none()
        {
            bail!(help_text.to_string());
        }

        Ok(match sample_name_spec {
            SampleNameSpec::Names(names) => names,
            SampleNameSpec::Any => unreachable!(),
        })
    }

    /// return the default help text for *_COUNTER_* pipelines
    pub fn count_help() -> &'static str {
        FQ_HELP
    }
}
