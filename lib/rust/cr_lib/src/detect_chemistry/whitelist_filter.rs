#![deny(missing_docs)]
use super::chemistry_filter::{ChemistryFilter, DetectChemistryUnit};
use super::errors::DetectChemistryErrors;
use crate::detect_chemistry::chemistry_filter::DetectChemistryCandidates;
use anyhow::Result;
use barcode::{BarcodeConstruct, BcSegSeq, Whitelist, WhitelistSpec};
use cr_types::chemistry::{ChemistryDef, ChemistryName};
use cr_types::rna_read::extract_barcode;
use fastq_set::read_pair::ReadPair;
use itertools::Itertools;
use metric::{PercentMetric, TxHashMap};
use std::ops::Range;

#[derive(Default)]
pub(crate) struct WhitelistMatchFilter {
    /// On-demand cache of in-memory whitelists.
    cache: WhitelistCache,
}

impl ChemistryFilter for WhitelistMatchFilter {
    const CONTEXT: &'static str = "Whitelist based filtering";
    fn process_unit(
        &mut self,
        candidates: &DetectChemistryCandidates,
        unit: &DetectChemistryUnit,
        reads: &[ReadPair],
    ) -> Result<DetectChemistryCandidates, DetectChemistryErrors> {
        // Special case for HD CELLRANGER-7761
        if let Ok((name, _def)) = candidates.iter().exactly_one()
            && name.is_spatial()
        {
            return Ok(candidates.clone());
        }

        let min_frac_whitelist_match = 0.1;
        let mut all_wl_matches = TxHashMap::default();
        let mut best_wl_match = 0.0f64;
        for (&name, def) in candidates {
            let wl_matches = self.count_matches(def, reads)?.fraction();
            // TODO: Handle reads_with_bc = 0
            if wl_matches > best_wl_match {
                best_wl_match = wl_matches;
            }
            all_wl_matches.insert(name, (def, wl_matches));
            println!("{name}: {wl_matches}");
        }

        let sufficient_matches: TxHashMap<_, _> = all_wl_matches
            .iter()
            .filter(|&(_, &(_def, v))| v > min_frac_whitelist_match)
            .collect();
        if sufficient_matches.is_empty() {
            return Err(DetectChemistryErrors::NotEnoughWhitelistMatch {
                frac_matches: all_wl_matches
                    .iter()
                    .map(|(&name, &(_def, matches))| (name, matches))
                    .collect(),
                unit: Box::new(unit.clone()),
            });
        }

        let mut compatible_chems: DetectChemistryCandidates = sufficient_matches
            .iter()
            .filter_map(|(&&chem, &&(def, wl_matches))| {
                (wl_matches > 0.75f64 * best_wl_match).then_some((chem, def.clone()))
            })
            .collect();

        // Prefer the more specific chemistry when MFRP_47 and at least one other more specific
        // MFRP chemistry are compatible.
        if compatible_chems
            .keys()
            .any(|name| matches!(name, ChemistryName::MFRP_RNA | ChemistryName::MFRP_Ab))
        {
            compatible_chems.remove(&ChemistryName::MFRP_47);
        }

        // Use the chemistry with more valid barcodes when both R1 and R2 chemistries are compatible.
        let select_superior_matches_of_pair = [
            (ChemistryName::MFRP_RNA, ChemistryName::MFRP_RNA_R1),
            (ChemistryName::MFRP_Ab, ChemistryName::MFRP_Ab_R1),
            (ChemistryName::MFRP_Ab_R2pos50, ChemistryName::MFRP_Ab_R1),
            (ChemistryName::MFRP_96_RNA_R2, ChemistryName::MFRP_96_R1),
            (ChemistryName::Flex_v2_RNA_R2, ChemistryName::Flex_v2_R1),
            (ChemistryName::Flex_v2_Ab_R2_45, ChemistryName::Flex_v2_R1),
            (ChemistryName::Flex_v2_Ab_R2_64, ChemistryName::Flex_v2_R1),
            (
                ChemistryName::Flex_v2_CRISPR_R2_1,
                ChemistryName::Flex_v2_R1,
            ),
        ];
        for (r2_chem, r1_chem) in select_superior_matches_of_pair {
            if compatible_chems.contains_key(&r2_chem) && compatible_chems.contains_key(&r1_chem) {
                let inferior_chemistry =
                    if sufficient_matches[&r2_chem].1 < sufficient_matches[&r1_chem].1 {
                        r2_chem
                    } else {
                        r1_chem
                    };
                assert!(compatible_chems.remove(&inferior_chemistry).is_some());
            }
        }

        Ok(compatible_chems)
    }
}

impl WhitelistMatchFilter {
    fn count_matches(
        &mut self,
        chem_def: &ChemistryDef,
        read_pairs: &[ReadPair],
    ) -> Result<WhitelistMatchStats, DetectChemistryErrors> {
        let (whitelist, barcode_lengths) = self
            .cache
            .fetch(chem_def.barcode_whitelist_spec())
            .map_err(DetectChemistryErrors::WhitelistLoadError)?
            .unzip();

        let mut stats = WhitelistMatchStats::default();

        let barcode_components = chem_def.barcode_construct();
        let barcode_extraction = chem_def.barcode_extraction();

        fn contains(
            whitelist: BarcodeConstruct<&Whitelist>,
            sseq: BarcodeConstruct<&[u8]>,
        ) -> bool {
            // NOTE: This is robust to a single N cycle
            whitelist
                .zip(sseq)
                .map_option(|(wl, seq)| {
                    wl.match_to_whitelist_allow_one_n(BcSegSeq::from_bytes(seq), false)
                })
                .is_some()
        }

        for read_pair in read_pairs {
            stats.total_reads += 1;

            if let Ok((_bc_ranges, bc)) = extract_barcode(
                read_pair,
                barcode_components,
                barcode_extraction,
                whitelist,
                barcode_lengths,
                0,
            ) {
                stats.reads_with_bc += 1;
                if contains(whitelist, bc.sequence()) {
                    stats.reads_with_bc_in_wl += 1;
                }
            }
        }
        Ok(stats)
    }
}

#[derive(Default)]
struct WhitelistMatchStats {
    total_reads: i64,
    reads_with_bc: i64,
    reads_with_bc_in_wl: i64,
}

impl WhitelistMatchStats {
    fn fraction(self) -> f64 {
        PercentMetric::from((self.reads_with_bc_in_wl, self.reads_with_bc))
            .fraction()
            .unwrap_or(0.0f64) // Return 0 if the denominator is 0
    }
}

/// Cache of whitelists to avoid repeatedly loading the same one.
///
/// Also caches the whitelist barcode length range.
#[derive(Default)]
struct WhitelistCache {
    cache: TxHashMap<WhitelistSpec, (Whitelist, Range<usize>)>,
}

impl WhitelistCache {
    /// Fetch whitelists for specs, possibly loading if we don't have one already.
    pub(crate) fn fetch(
        &mut self,
        specs: BarcodeConstruct<&WhitelistSpec>,
    ) -> Result<BarcodeConstruct<(&Whitelist, &Range<usize>)>> {
        // Ensure we have all the specs loaded.
        for spec in specs.iter() {
            if !self.cache.contains_key(spec) {
                // Load this whitelist as a "plain" whitelist, regardless of type.
                //
                // This is sufficient for checking membership in the whitelist, but neglects
                // all translation information. Using this method instead of as_whitelist
                // reduces memory consumption since we are working in a context where
                // we only care about matching to the whitelist, but not using translation.
                let whitelist = Whitelist::Plain(spec.as_source()?.as_set()?);
                let lengths = whitelist.sequence_lengths();
                self.cache.insert(spec.clone(), (whitelist, lengths));
            }
        }
        Ok(specs.map(|spec| {
            let (wl, lens) = &self.cache[spec];
            (wl, lens)
        }))
    }
}
