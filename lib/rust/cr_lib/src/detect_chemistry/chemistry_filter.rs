#![deny(missing_docs)]
use super::errors::DetectChemistryErrors;
use super::identity_check::check_read_identity;
use anyhow::Result;
use cr_types::LibraryType;
use cr_types::chemistry::{ChemistryDef, ChemistryName};
use cr_types::sample_def::{FastqMode, SampleDef};
use fastq_set::long_read_adapter_trimmer::{LongReadTrimmer, NoOpVisitor};
use fastq_set::read_pair::{ReadHeader, ReadPair};
use fastq_set::read_pair_iter::FastqSet;
use itertools::Itertools;
use metric::{TxHashMap, TxHashSet};
use stats::ReservoirSampler;
use std::collections::HashMap;
use std::fmt::{Display, Formatter};
use std::path::PathBuf;

const RANDOM_SEED: u64 = 124;

pub(crate) type DetectChemistryCandidates = TxHashMap<ChemistryName, ChemistryDef>;

pub(crate) type LoadedUnit = (DetectChemistryUnit, Vec<ReadPair>);

pub(crate) trait ChemistryFilter {
    const CONTEXT: &'static str;
    fn process_unit(
        &mut self,
        candidates: &DetectChemistryCandidates,
        unit: &DetectChemistryUnit,
        reads: &[ReadPair],
    ) -> Result<DetectChemistryCandidates, DetectChemistryErrors>;

    /// Filter a collection of explicit chemistry defs, indexed by their name.
    fn filter_chemistries(
        &mut self,
        candidates: &DetectChemistryCandidates,
        units: &[LoadedUnit],
    ) -> Result<DetectChemistryCandidates, DetectChemistryErrors> {
        let per_unit_chems = units
            .iter()
            .map(|(unit, read_pairs)| {
                println!("\n{} for {unit}", Self::CONTEXT);
                let compatible_chems = self.process_unit(candidates, unit, read_pairs)?;
                println!(
                    "Compatible chemistries: {}",
                    compatible_chems.keys().sorted().format(", ")
                );
                assert!(!compatible_chems.is_empty());
                Ok(compatible_chems)
            })
            .try_collect()?;

        // Chemistries that are compatible with all units.
        let mut result = intersect_many(&per_unit_chems);
        if result.is_empty() {
            println!("Result: none");
        } else {
            println!("Result: {}", result.keys().sorted().format(", "));
        }

        if result.is_empty() {
            return Err(DetectChemistryErrors::ConflictingChemistries {
                units: units.iter().map(|(unit, _)| unit.clone()).collect(),
                per_unit_chems,
                context: Self::CONTEXT,
            });
        }

        self.post_process(&mut result);

        Ok(result)
    }

    /// Perform any filter-specific post-processing on the candidates, in-place.
    fn post_process(&self, _candidates: &mut DetectChemistryCandidates) {}
}

/// Map a set of named chemistries into a mapping to their defs.
pub(crate) fn named_candidates(candidates: &TxHashSet<ChemistryName>) -> DetectChemistryCandidates {
    candidates
        .iter()
        .map(|&name| (name, ChemistryDef::named(name)))
        .collect()
}

/// Return a candidate collection which is the intersection of several.
///
/// Panic if two key intersections that appear in the final intersection do not
/// have identical chemistry defs.
fn intersect_many<'a>(
    items: impl IntoIterator<Item = &'a DetectChemistryCandidates>,
) -> DetectChemistryCandidates {
    let mut iter = items.into_iter();
    let Some(mut first) = iter.next().cloned() else {
        return HashMap::default();
    };
    let mut index_of_all: HashMap<_, _> =
        first.iter().map(|(&k, v)| (k, vec![v.clone()])).collect();
    for other in iter {
        // Collect all variations on each def for a given name,
        for (&name, def) in other {
            index_of_all.entry(name).or_default().push(def.clone());
        }
        first.retain(|item, _def| other.contains_key(item));
    }
    // For all the keys in the intersection, make sure all defs are identical.
    for name in first.keys() {
        for (c0, c1) in index_of_all[name].iter().tuples() {
            assert!(
                c0 == c1,
                "tried to take the intersection of two chemistries that are not identical matches:\n{c0:?}\n\n{c1:?}"
            );
        }
    }
    first
}

/// A "unit" of chemistry detection.
///
/// Defined for a specific library type, using a specific subset of the reads
/// sourced for that library (could be all reads, or a specific sample index).
///
/// For some cases, we load multiple units for a single library type to ensure
/// that we can detect cases where the user may have for example mixed up their
/// sample indexes and accidentally provided data from different libraries
/// sequenced on the same flow cell.
#[derive(Debug, Clone)]
pub struct DetectChemistryUnit {
    pub fastq_mode: FastqMode,
    pub fastqs: Vec<FastqSet>,
    pub group: Option<String>,
    pub read_path: PathBuf,
    pub r1_length: Option<usize>,
    pub r2_length: Option<usize>,
}

impl Display for DetectChemistryUnit {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}\"{}\"",
            match self.group {
                Some(ref g) => format!("Sample {g} in "),
                None => String::new(),
            },
            self.read_path.display()
        )
    }
}

impl DetectChemistryUnit {
    pub fn sampled_read_pairs(
        &self,
        sample_reads: usize,
        total_reads: usize,
        adapter_trimmer: Option<&LongReadTrimmer>,
    ) -> Result<Vec<ReadPair>> {
        let mut read_sampler = ReservoirSampler::new(sample_reads, RANDOM_SEED);
        for fastq in &self.fastqs {
            let read_pair_iter = fastq.read_iter(ReadHeader::KeepId, None)?;
            if read_pair_iter.get_is_single_ended() {
                println!("{self} is a single-end library!");
            } else {
                println!("{self} is a paired-end library!");
            }

            for read_pair in read_pair_iter {
                let mut read_pair = read_pair?;
                if let Some(trimmer) = adapter_trimmer {
                    let Some(trimmed_read) =
                        trimmer.trim_read(read_pair, ReadHeader::KeepId, &mut NoOpVisitor)?
                    else {
                        continue;
                    };
                    read_pair = trimmed_read;
                }
                read_sampler.add(read_pair);
                if read_sampler.num_items_seen() >= total_reads {
                    return Ok(read_sampler.done());
                }
            }
        }
        Ok(read_sampler.done())
    }

    pub fn check_read_identity(&self) -> Result<Vec<(&FastqSet, (u64, u64))>> {
        const MAX_READ_PAIRS: usize = 100;
        let mut result = vec![];
        for fastq in &self.fastqs {
            let read_pairs: Vec<_> = fastq
                .read_iter(ReadHeader::KeepId, None)?
                .take(MAX_READ_PAIRS)
                .try_collect()?;
            let (r1_hash, r2_hash) = check_read_identity(self, read_pairs.as_slice())?;
            result.push((fastq, (r1_hash, r2_hash)));
        }
        Ok(result)
    }
}

/// Generate per-library-type units from the input sample defs.
pub fn detect_chemistry_units(
    sample_defs: &[SampleDef],
    r1_length: Option<usize>,
    r2_length: Option<usize>,
) -> Result<HashMap<LibraryType, Vec<DetectChemistryUnit>>> {
    let mut units: HashMap<_, Vec<DetectChemistryUnit>> = HashMap::new();
    for sdef in sample_defs {
        for (group, fastqs) in sdef.per_group_fastq_def()? {
            // Trim lengths from global args take precedence over
            // length specifed in sample_def (see MAKE_SHARD)
            let r1_length = r1_length.or(sdef.r1_length);
            let r2_length = r2_length.or(sdef.r2_length);

            let library_type = sdef.library_type.unwrap_or_default();

            assert!(!fastqs.is_empty());
            units
                .entry(library_type)
                .or_default()
                .push(DetectChemistryUnit {
                    fastq_mode: sdef.fastq_mode,
                    fastqs,
                    group,
                    read_path: sdef.read_path.clone(),
                    r1_length,
                    r2_length,
                });
        }
    }
    Ok(units)
}

/// Fully-loaded units broken up by library type.
pub(crate) type UnitsPerLibrary = HashMap<LibraryType, Vec<LoadedUnit>>;
