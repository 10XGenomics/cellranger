use super::errors::DetectChemistryErrors;
use super::identity_check::check_read_identity;
use anyhow::Result;
use cr_types::chemistry::ChemistryName;
use cr_types::sample_def::SampleDef;
use cr_types::LibraryType;
use fastq_set::filenames::FindFastqs;
use fastq_set::read_pair::ReadPair;
use fastq_set::read_pair_iter::{InputFastqs, ReadPairIter};
use itertools::Itertools;
use metric::TxHashSet;
use parameters_toml;
use stats::ReservoirSampler;
use std::fmt::{Display, Formatter};
use std::path::PathBuf;

const RANDOM_SEED: u64 = 124;

pub(crate) trait ChemistryFilter<'a> {
    fn context() -> &'static str;
    fn allowed_chemistries(&self) -> &'a TxHashSet<ChemistryName>;
    fn input_chemistries(&self) -> TxHashSet<ChemistryName>;
    fn process_unit(
        &mut self,
        unit: &DetectChemistryUnit,
        reads: &[ReadPair],
    ) -> Result<TxHashSet<ChemistryName>, DetectChemistryErrors>;

    fn filter_chemistries(
        &mut self,
        units: &[(DetectChemistryUnit, Vec<ReadPair>)],
    ) -> Result<TxHashSet<ChemistryName>, DetectChemistryErrors> {
        let per_unit_chems: Vec<_> = units
            .iter()
            .map(|(unit, read_pairs)| {
                println!("\n{} for {unit}", Self::context());
                let compatible_chems = self.process_unit(unit, read_pairs)?;
                println!(
                    "Compatible chemistries: {}",
                    compatible_chems.iter().sorted().format(", ")
                );
                assert!(!compatible_chems.is_empty());
                Ok(compatible_chems)
            })
            .try_collect()?;

        // Chemistries that are compatible with all libraries.
        let result = per_unit_chems
            .iter()
            .fold(self.input_chemistries(), |acc, x| {
                acc.intersection(x).copied().collect()
            });
        if result.is_empty() {
            println!("Result: none");
        } else {
            println!("Result: {}", result.iter().sorted().format(", "));
        }

        if result.is_empty() {
            return Err(DetectChemistryErrors::ConflictingChemistries {
                units: units.iter().map(|(unit, _)| unit.clone()).collect(),
                per_unit_chems,
                context: Self::context(),
            });
        }

        let result = self.post_process(result)?;

        // if we have more than 1 possible chemistry, remove the non-allowed ones
        Ok(if result.len() > 1 {
            result
                .intersection(self.allowed_chemistries())
                .copied()
                .collect()
        } else {
            result
        })
    }

    /// Perform post-filtering processing on the mixture of chemistries resulting
    /// from this filter.
    fn post_process(
        &self,
        result: TxHashSet<ChemistryName>,
    ) -> Result<TxHashSet<ChemistryName>, DetectChemistryErrors> {
        Ok(result)
    }
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
    pub fn sampled_read_pairs(&self) -> Result<Vec<ReadPair>> {
        let detect_chemistry_sample_reads = *parameters_toml::detect_chemistry_sample_reads()?;
        let detect_chemistry_total_reads = *parameters_toml::detect_chemistry_total_reads()?;

        let mut read_sampler = ReservoirSampler::new(detect_chemistry_sample_reads, RANDOM_SEED);
        for fastq in &self.fastqs {
            let read_pair_iter = ReadPairIter::from_fastq_files(fastq)?;
            if read_pair_iter.get_is_single_ended() {
                println!("{self} is a single-end library!");
            } else {
                println!("{self} is a paired-end library!");
            }
            for read_pair in read_pair_iter {
                read_sampler.add(read_pair?);
                if read_sampler.num_items_seen() >= detect_chemistry_total_reads {
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
