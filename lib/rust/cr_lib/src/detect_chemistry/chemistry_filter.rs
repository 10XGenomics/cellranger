use super::errors::DetectChemistryErrors;
use super::identity_check::check_read_identity;
use anyhow::Result;
use cr_types::chemistry::ChemistryName;
use cr_types::rna_read::LegacyLibraryType;
use cr_types::sample_def::SampleDef;
use fastq_set::filenames::FindFastqs;
use fastq_set::read_pair::ReadPair;
use fastq_set::read_pair_iter::{InputFastqs, ReadPairIter};
use itertools::Itertools;
use metric::{set, TxHashSet};
use parameters_toml;
use stats::ReservoirSampler;
use std::fmt::{Display, Formatter};
use std::iter::zip;
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
        units: &[DetectChemistryUnit],
        read_pairs_all_units: &[Vec<ReadPair>],
    ) -> Result<TxHashSet<ChemistryName>, DetectChemistryErrors> {
        use ChemistryName::{FivePrimePE, FivePrimeR2, MFRP, MFRP_47};

        let mut result: TxHashSet<_> = self.input_chemistries();
        let mut any_compatible: TxHashSet<_> = TxHashSet::default();
        let mut per_unit_chems = Vec::new(); // For error reporting
        for (unit, read_pairs) in zip(units, read_pairs_all_units) {
            println!("\n\n{} for {unit}", Self::context());
            let compatible_chems = self.process_unit(unit, read_pairs)?;
            println!("Compatible chemistries: {compatible_chems:?}");
            debug_assert!(!compatible_chems.is_empty());
            for chem in &compatible_chems {
                any_compatible.insert(*chem);
            }
            result = result.intersection(&compatible_chems).copied().collect();
            println!("Result: {result:?}");
            per_unit_chems.push(compatible_chems);
        }

        // we don't want to allow a situation where (SC5P-R2 + SC5P-PE -> SC5P-R2)
        if (any_compatible.contains(&FivePrimePE) && any_compatible.contains(&FivePrimeR2))
            && (!result.contains(&FivePrimePE) && result.contains(&FivePrimeR2))
        {
            return Err(DetectChemistryErrors::UnsupportedMixture {
                mixture: set![FivePrimeR2, FivePrimePE],
                specify: FivePrimeR2,
            });
        }
        if result.is_empty() {
            return Err(DetectChemistryErrors::ConflictingChemistries {
                units: units.to_vec(),
                per_unit_chems,
                context: Self::context(),
            });
        }

        Ok(self.finalize(if result == set![MFRP, MFRP_47] {
            // Use MFRP when both MFRP and MFRP-47 are compatible, because MFRP is a subset of MFRP-47.
            set![MFRP]
        } else {
            result
        }))
    }
    fn default_finalize(&self, chemistries: TxHashSet<ChemistryName>) -> TxHashSet<ChemistryName> {
        // if we have more than 1 possible chemistry, remove the non-allowed ones
        if chemistries.len() > 1 {
            chemistries
                .intersection(self.allowed_chemistries())
                .copied()
                .collect()
        } else {
            chemistries
        }
    }
    fn finalize(&self, chemistries: TxHashSet<ChemistryName>) -> TxHashSet<ChemistryName> {
        self.default_finalize(chemistries)
    }
}

#[derive(Debug, Clone)]
pub struct DetectChemistryUnit {
    pub fastqs: Vec<InputFastqs>,
    pub group: Option<String>,
    pub read_path: PathBuf,
    pub library_type: LegacyLibraryType,
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
                None => "".into(),
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
            for read_pair in ReadPairIter::from_fastq_files(fastq)? {
                read_sampler.add(read_pair?);
                if read_sampler.num_items_seen() >= detect_chemistry_total_reads {
                    return Ok(read_sampler.done());
                }
            }
        }
        Ok(read_sampler.done())
    }

    pub fn check_read_identity(&self) -> Result<Vec<(u64, u64)>> {
        const MAX_READ_PAIRS: usize = 100;
        let mut result = vec![];
        for fastq in &self.fastqs {
            let read_pairs: Vec<_> = ReadPairIter::from_fastq_files(fastq)?
                .take(MAX_READ_PAIRS)
                .try_collect()?;
            result.push(check_read_identity(self, read_pairs.as_slice())?);
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
            })
        }
    }
    Ok(units)
}
