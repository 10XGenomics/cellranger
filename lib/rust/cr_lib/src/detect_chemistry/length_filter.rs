use super::chemistry_filter::{ChemistryFilter, DetectChemistryUnit};
use super::errors::DetectChemistryErrors;
use anyhow::Result;
use cr_types::chemistry::{ChemistryDef, ChemistryName};
use cr_types::reference::feature_reference::{FeatureConfig, FeatureReferenceFile, FeatureType};
use fastq_set::read_pair::ReadPair;
use fastq_set::WhichRead;
use metric::{SimpleHistogram, TxHashMap, TxHashSet};
use std::fmt;

pub(crate) const WHICH_LEN_READS: [WhichRead; 3] = [WhichRead::R1, WhichRead::R2, WhichRead::I1];

pub(crate) struct LengthFilter<'a> {
    allowed_chems: &'a TxHashSet<ChemistryName>,
    chem_defs: Vec<ChemistryDef>,
    min_feature_read_lengths: TxHashMap<FeatureType, TxHashMap<WhichRead, usize>>,
}

#[derive(Debug)]
pub(crate) struct LengthStats {
    median: TxHashMap<WhichRead, usize>,
}

impl fmt::Display for LengthStats {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for which_read in &WHICH_LEN_READS {
            writeln!(
                f,
                " - {} median length = {}",
                which_read, self.median[which_read]
            )?;
        }
        Ok(())
    }
}

impl LengthStats {
    fn is_compatible(
        &self,
        chem_def: &ChemistryDef,
        skip: Option<&TxHashMap<WhichRead, usize>>,
    ) -> bool {
        WHICH_LEN_READS
            .iter()
            .filter(|r| !skip.is_some_and(|m| m.contains_key(r)))
            .all(|&r| self.median[&r] >= chem_def.min_read_length(r))
    }
}

impl<'a> ChemistryFilter<'a> for LengthFilter<'a> {
    fn context() -> &'static str {
        "Length based filtering"
    }
    fn allowed_chemistries(&self) -> &'a TxHashSet<ChemistryName> {
        self.allowed_chems
    }
    fn input_chemistries(&self) -> TxHashSet<ChemistryName> {
        self.chem_defs.iter().map(|def| def.name).collect()
    }
    fn process_unit(
        &mut self,
        unit: &DetectChemistryUnit,
        read_pairs: &[ReadPair],
    ) -> Result<TxHashSet<ChemistryName>, DetectChemistryErrors> {
        let (stats, compatible_chems, max_lengths) = self.process_reads(unit, read_pairs)?;
        println!("{stats:#?}");
        if compatible_chems.is_empty() {
            Err(DetectChemistryErrors::NotEnoughReadLength {
                stats,
                unit: Box::new(unit.clone()),
                chems: self.input_chemistries(),
                max_lengths,
            })
        } else {
            Ok(compatible_chems)
        }
    }
    fn post_process(
        &self,
        mut result: TxHashSet<ChemistryName>,
    ) -> Result<TxHashSet<ChemistryName>, DetectChemistryErrors> {
        use ChemistryName::{
            FivePrimePE, FivePrimePEV3, FivePrimeR2, FivePrimeR2V3, VdjPE, VdjPEV3, VdjR2, VdjR2V3,
        };
        if result.contains(&FivePrimePE) {
            result.remove(&FivePrimeR2);
        }
        if result.contains(&FivePrimePEV3) {
            result.remove(&FivePrimeR2V3);
        }
        if result.contains(&VdjPE) {
            result.remove(&VdjR2);
        }
        if result.contains(&VdjPEV3) {
            result.remove(&VdjR2V3);
        }
        Ok(result)
    }
}

impl<'a> LengthFilter<'a> {
    pub(crate) fn new<'b>(
        allowed_chems: &'a TxHashSet<ChemistryName>,
        chems: impl IntoIterator<Item = &'b ChemistryName>,
        feature_ref: Option<&FeatureReferenceFile>,
        feature_config: Option<&FeatureConfig>,
    ) -> Result<Self> {
        let min_feature_read_lengths = if let Some(ref_path) = feature_ref {
            ref_path.read(feature_config)?.min_feature_read_lengths()
        } else {
            TxHashMap::default()
        };
        Ok(LengthFilter {
            allowed_chems,
            chem_defs: chems.into_iter().map(|c| ChemistryDef::named(*c)).collect(),
            min_feature_read_lengths,
        })
    }

    #[allow(clippy::type_complexity)]
    pub(crate) fn process_reads(
        &self,
        unit: &DetectChemistryUnit,
        read_pairs: &[ReadPair],
    ) -> Result<
        (
            LengthStats,
            TxHashSet<ChemistryName>,
            TxHashMap<WhichRead, usize>,
        ),
        DetectChemistryErrors,
    > {
        let mut max_lengths = TxHashMap::default();
        if let Some(l) = unit.r1_length {
            max_lengths.insert(WhichRead::R1, l);
        }
        if let Some(l) = unit.r2_length {
            max_lengths.insert(WhichRead::R2, l);
        }

        let mut read_length_histograms = TxHashMap::default();
        for rp in read_pairs {
            for &which_read in &WHICH_LEN_READS {
                let max_len = max_lengths.get(&which_read).unwrap_or(&usize::MAX);
                read_length_histograms
                    .entry(which_read)
                    .or_insert_with(SimpleHistogram::default)
                    .observe(&rp.len(which_read).unwrap_or(0).min(*max_len));
            }
        }
        let stats = LengthStats {
            median: read_length_histograms
                .iter()
                .map(|(&k, v)| (k, v.percentiles([50.0]).map_or(0, |p| p[0] as usize)))
                .collect(),
        };

        if let Some(feature_type) = unit.library_type.feature_barcode_type() {
            // non-gene feature types get evaluated against the feature reference
            let Some(min_lengths) = self
                .min_feature_read_lengths
                .get(&FeatureType::Barcode(feature_type))
            else {
                return Err(DetectChemistryErrors::FeatureTypeNotInReference {
                    feature_type: FeatureType::Barcode(feature_type),
                });
            };
            for (&read, &obs_len) in &stats.median {
                let min_len = min_lengths.get(&read).copied().unwrap_or(0);
                if obs_len < min_len {
                    println!("read {read} with obs_len {obs_len} failed min_len {min_len}");
                    return Err(DetectChemistryErrors::FeatureTypeNotEnoughReadLength {
                        unit: Box::new(unit.clone()),
                        stats,
                        min_lengths: self.min_feature_read_lengths.clone(),
                        max_lengths,
                    });
                }
            }
            let compatible_chems = self
                .chem_defs
                .iter()
                .filter(|def| stats.is_compatible(def, Some(min_lengths)))
                .map(|def| def.name)
                .collect();
            return Ok((stats, compatible_chems, max_lengths));
        }

        let compatible_chems = self
            .chem_defs
            .iter()
            .filter(|def| stats.is_compatible(def, None))
            .map(|def| def.name)
            .collect();

        Ok((stats, compatible_chems, max_lengths))
    }
}
