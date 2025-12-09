#![deny(missing_docs)]
use super::errors::DetectChemistryErrors;
use super::identity_check::check_read_identity;
use anyhow::Result;
use cr_types::LibraryType;
use cr_types::chemistry::{ChemistryDef, ChemistryName};
use cr_types::sample_def::SampleDef;
use fastq_set::filenames::FindFastqs;
use fastq_set::read_pair::ReadPair;
use fastq_set::read_pair_iter::{InputFastqs, ReadPairIter};
use itertools::Itertools;
use metric::{TxHashMap, TxHashSet};
use stats::ReservoirSampler;
use std::fmt::{Display, Formatter};
use std::path::PathBuf;

const RANDOM_SEED: u64 = 124;

pub(crate) type DetectChemistryCandidates = TxHashMap<ChemistryName, ChemistryDef>;

pub(crate) trait ChemistryFilter {
    const CONTEXT: &'static str;
    fn process_unit(
        &mut self,
        candidates: &DetectChemistryCandidates,
        unit: &DetectChemistryUnit,
        reads: &[ReadPair],
    ) -> Result<DetectChemistryCandidates, DetectChemistryErrors>;

    /// Perform filtering on each individual input unit.
    /// Return a collection of compatible chemistries per input unit.
    ///
    /// The elements in the output vec correspond to the input units.
    fn filter_chemistries_per_unit(
        &mut self,
        candidates: &DetectChemistryCandidates,
        units: &[(DetectChemistryUnit, Vec<ReadPair>)],
    ) -> Result<Vec<DetectChemistryCandidates>, DetectChemistryErrors> {
        units
            .iter()
            .map(|unit| self.filter_chemistries_single_unit(candidates, unit))
            .try_collect()
    }

    /// Perform filtering on an individual input unit using explicit chemistries.
    fn filter_chemistries_single_unit(
        &mut self,
        candidates: &DetectChemistryCandidates,
        (unit, read_pairs): &(DetectChemistryUnit, Vec<ReadPair>),
    ) -> Result<DetectChemistryCandidates, DetectChemistryErrors> {
        println!("\n{} for {unit}", Self::CONTEXT);
        let compatible_chems = self.process_unit(candidates, unit, read_pairs)?;
        println!(
            "Compatible chemistries: {}",
            compatible_chems.keys().sorted().format(", ")
        );
        assert!(!compatible_chems.is_empty());
        Ok(compatible_chems)
    }

    /// Collapse a collection of per-unit chemistries into their intersection.
    ///
    /// Perform filter-specific post-processing.
    fn reduce_per_unit_chems(
        units: &[(DetectChemistryUnit, Vec<ReadPair>)],
        per_unit_chems: Vec<DetectChemistryCandidates>,
    ) -> Result<DetectChemistryCandidates, DetectChemistryErrors> {
        // Chemistries that are compatible with all libraries.
        let result = intersect_many(&per_unit_chems);
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

        Ok(result)
    }

    /// Filter a collection of explicit chemistry defs, indexed by their name.
    fn filter_chemistries(
        &mut self,
        candidates: &DetectChemistryCandidates,
        units: &[(DetectChemistryUnit, Vec<ReadPair>)],
    ) -> Result<DetectChemistryCandidates, DetectChemistryErrors> {
        let per_unit_chems = self.filter_chemistries_per_unit(candidates, units)?;

        Self::reduce_per_unit_chems(units, per_unit_chems)
    }

    /// Filter a collection of named default chemistries.
    fn filter_named_chemistries(
        &mut self,
        candidates: &TxHashSet<ChemistryName>,
        units: &[(DetectChemistryUnit, Vec<ReadPair>)],
    ) -> Result<TxHashSet<ChemistryName>, DetectChemistryErrors> {
        Ok(self
            .filter_chemistries(&named_candidates(candidates), units)?
            .into_keys()
            .collect())
    }
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
/// Assumes the chemistry values associated with each key are identical.
fn intersect_many<'a>(
    items: impl IntoIterator<Item = &'a DetectChemistryCandidates>,
) -> DetectChemistryCandidates {
    let mut iter = items.into_iter();
    let Some(mut first) = iter.next().cloned() else {
        return Default::default();
    };
    for other in iter {
        first.retain(|item, _def| other.contains_key(item));
    }
    first
}

#[derive(Debug, Clone)]
pub struct DetectChemistryUnit {
    pub fastqs: Vec<InputFastqs>,
    pub group: Option<String>,
    pub read_path: PathBuf,
    pub library_type: LibraryType,
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
    ) -> Result<Vec<ReadPair>> {
        let mut read_sampler = ReservoirSampler::new(sample_reads, RANDOM_SEED);
        for fastq in &self.fastqs {
            let read_pair_iter = ReadPairIter::from_fastq_files(fastq)?;
            if read_pair_iter.get_is_single_ended() {
                println!("{self} is a single-end library!");
            } else {
                println!("{self} is a paired-end library!");
            }
            for read_pair in read_pair_iter {
                read_sampler.add(read_pair?);
                if read_sampler.num_items_seen() >= total_reads {
                    return Ok(read_sampler.done());
                }
            }
        }
        Ok(read_sampler.done())
    }

    pub fn check_read_identity(&self) -> Result<Vec<(&InputFastqs, (u64, u64))>> {
        const MAX_READ_PAIRS: usize = 100;
        let mut result = vec![];
        for fastq in &self.fastqs {
            let read_pairs: Vec<_> = ReadPairIter::from_fastq_files(fastq)?
                .take(MAX_READ_PAIRS)
                .try_collect()?;
            let (r1_hash, r2_hash) = check_read_identity(self, read_pairs.as_slice())?;
            result.push((fastq, (r1_hash, r2_hash)));
        }
        Ok(result)
    }
}

pub fn detect_chemistry_units(
    sample_defs: &[SampleDef],
    r1_length: Option<usize>,
    r2_length: Option<usize>,
) -> Result<Vec<DetectChemistryUnit>> {
    let mut units = Vec::new();
    for sdef in sample_defs {
        for (group, fastq_def) in sdef.per_group_fastq_def()? {
            // Trim lengths from global args take precedence over
            // length specifed in sample_def (see MAKE_SHARD)
            let r1_length = r1_length.or(sdef.r1_length);
            let r2_length = r2_length.or(sdef.r2_length);

            let fastqs = fastq_def.find_fastqs()?;
            assert!(!fastqs.is_empty());
            units.push(DetectChemistryUnit {
                fastqs,
                group,
                read_path: sdef.read_path.clone(),
                library_type: sdef.library_type.unwrap_or_default(),
                r1_length,
                r2_length,
            });
        }
    }
    Ok(units)
}
