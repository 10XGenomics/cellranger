#![deny(missing_docs)]
use super::chemistry_filter::DetectChemistryUnit;
use super::length_filter::{LengthStats, WHICH_LEN_READS};
use super::mapping_filter::MappingStats;
use crate::detect_chemistry::chemistry_filter::DetectChemistryCandidates;
use barcode::whitelist::BarcodeId;
use cr_types::ERROR_CODE_INFO;
use cr_types::chemistry::{AutoOrRefinedChemistry, ChemistryDef, ChemistryName};
use cr_types::reference::feature_reference::FeatureType;
use fastq_set::WhichRead;
use itertools::{Itertools, multizip};
use metric::{TxHashMap, TxHashSet};
use ordered_float::NotNan;
use std::cmp::Reverse;
use std::fmt;
use std::iter::Iterator;

#[derive(Debug, thiserror::Error)]
pub(crate) enum DetectChemistryErrors {
    NotEnoughReads {
        num_reads: usize,
        min_reads: usize,
        unit: Box<DetectChemistryUnit>,
    },
    ChemistryNotAllowed {
        input: AutoOrRefinedChemistry,
        allowed: Vec<AutoOrRefinedChemistry>,
    },
    NotEnoughWhitelistMatch {
        frac_matches: TxHashMap<ChemistryName, f64>,
        unit: Box<DetectChemistryUnit>,
    },
    ConflictingChemistries {
        units: Vec<DetectChemistryUnit>,
        per_unit_chems: Vec<DetectChemistryCandidates>,
        context: &'static str,
    },
    NotEnoughMapping {
        stats: MappingStats,
        unit: Box<DetectChemistryUnit>,
        chems: DetectChemistryCandidates,
    },
    NotEnoughReadLength {
        stats: LengthStats,
        unit: Box<DetectChemistryUnit>,
        chems: DetectChemistryCandidates,
        max_lengths: TxHashMap<WhichRead, usize>,
    },
    WhitelistLoadError(anyhow::Error),
    FeatureTypeNotEnoughReadLength {
        stats: LengthStats,
        unit: Box<DetectChemistryUnit>,
        min_lengths: TxHashMap<FeatureType, TxHashMap<WhichRead, usize>>,
        max_lengths: TxHashMap<WhichRead, usize>,
    },
    FeatureTypeNotInReference {
        feature_type: FeatureType,
    },
    ProbeBarcodeMixture {
        unit: Box<DetectChemistryUnit>,
        mixture: TxHashSet<BarcodeId>,
    },
}

impl fmt::Display for DetectChemistryErrors {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        #[allow(clippy::enum_glob_use)]
        use DetectChemistryErrors::*;
        let msg = match self {
            NotEnoughReads {
                num_reads,
                min_reads,
                unit,
            } => format!(
                "TXRNGR10001: There were not enough reads to auto detect the chemistry: {unit}\n\
                 Note that you can avoid auto detection by specifying the specific chemistry type \
                 and version.\n\
                 - Minimum number of required reads = {min_reads}\n\
                 - Number of reads available = {num_reads}\n\
                 {ERROR_CODE_INFO}"
            ),
            ChemistryNotAllowed { input, allowed } => format!(
                "The chemistry name '{input}' is not allowed in this pipeline. Allowed values:\n{}",
                allowed.iter().format("\n")
            ),
            NotEnoughWhitelistMatch { frac_matches, unit } => {
                if frac_matches.len() == 1 {
                    let chem = *frac_matches.keys().next().unwrap();
                    let chem_def = ChemistryDef::named(chem);
                    let matches_percent = 100.0 * frac_matches[&chem];
                    if chem.is_spatial() {
                        format!(
                            "The input data has a very low rate of barcodes ({matches_percent:.1}%) that match the selected \
                            slide type.\nPlease check that your input data and slide serial number or chosen \
                            slide version are correct.\nInput: {unit}"
                        )
                    } else if chem.is_vdj() {
                        format!(
                            "In the input data, an extremely low rate of correct barcodes was observed for this \
                            chemistry ({matches_percent:.1}%).\nPlease check your input data. Note: manual chemistry detection \
                            is not required in most cases.\nInput: {unit}"
                        )
                    } else if let Some(probe_barcode) = chem_def.barcode_read_type().option_probe()
                    {
                        format!(
                            "You selected chemistry {chem}, which expects the gel-bead barcode \
                             sequence in {} and the probe barcode sequence in {probe_barcode}.\n\
                             In the input data, an extremely low rate of correct barcodes was \
                             observed for this chemistry ({matches_percent:.1}%).\n\
                             Please check your input data and chemistry selection. \
                             Note: manual chemistry detection is not required in most cases.\n\
                             Input: {unit}",
                            chem_def.barcode_read_type().gel_bead(),
                        )
                    } else {
                        format!(
                            "You selected chemistry {}, which expects the cell barcode sequence in {}.\n\
                            In the input data, an extremely low rate of correct barcodes was observed for this \
                            chemistry ({:.1}%).\nPlease check your input data and chemistry selection. \
                            Note: manual chemistry detection is not required in most cases.\nInput: {}",
                            chem,
                            chem_def.barcode_read_type().gel_bead(),
                            matches_percent,
                            unit
                        )
                    }
                } else {
                    let frac_matches_str = frac_matches
                        .iter()
                        .sorted_by_key(|&(k, &v)| (Reverse(NotNan::new(v).unwrap()), k))
                        .map(|(k, v)| format!("- {:.1}% for chemistry {k}", 100.0 * v))
                        .format("\n");
                    format!(
                        "TXRNGR10002: An extremely low rate of correct barcodes was observed for \
                         all the candidate chemistry choices for the input: {unit}. \
                         Please check your input data.\n\
                         {frac_matches_str}\n\
                         {ERROR_CODE_INFO}"
                    )
                }
            }
            ConflictingChemistries {
                units,
                per_unit_chems,
                context,
            } => format!(
                "We detected conflicting chemistries for different set of input fastqs using {}.\n -{}",
                context,
                multizip((units, per_unit_chems))
                    .map(|(u, c)| format!(
                        "One of [{}] is compatible with {u}",
                        c.keys().sorted().format(", ")
                    ))
                    .format("\n -")
            ),
            NotEnoughMapping { stats, unit, chems } => format!(
                "TXRNGR10003: Unable to distinguish between [{chems}] chemistries based on the R2 \
                 read mapping for {unit}.\n\
                 {stats}\n\
                 {help_text}\n\n\
                 Please validate the inputs and/or specify the chemistry via \
                 the --chemistry argument. {ERROR_CODE_INFO}",
                chems = chems.keys().sorted().format(", "),
                help_text = MappingStats::help_text(),
            ),
            NotEnoughReadLength {
                stats,
                unit,
                chems,
                max_lengths,
            } => format!(
                "TXRNGR10004: The read lengths are incompatible with all \
                 the chemistries for {unit}.\n\
                 {stats}\n\
                 The minimum read length for different chemistries are:\n{min_read_lengths}\n\n\
                 We expect that at least 50% of the reads exceed the minimum length.\n\
                 {trim_lengths}\
                 {ERROR_CODE_INFO}",
                min_read_lengths = chems
                    .keys()
                    .sorted()
                    .map(|&c| {
                        let def = ChemistryDef::named(c);
                        let read_lengths = WHICH_LEN_READS
                            .iter()
                            .map(|&which| format!("{which}: {}", def.min_read_length(which)))
                            .format(", ");
                        format!("{c:8} - {read_lengths}")
                    })
                    .format("\n"),
                trim_lengths = if max_lengths.is_empty() {
                    String::new()
                } else {
                    format!(
                        "NOTE: You have specified the trim lengths explicitly:\n{}\n",
                        WHICH_LEN_READS
                            .iter()
                            .filter_map(|r| max_lengths.get(r).map(|l| format!("- {r}: {l} bases")))
                            .format("\n")
                    )
                }
            ),
            WhitelistLoadError(err) => format!("An error occurred when loading a whitelist: {err}"),
            FeatureTypeNotEnoughReadLength {
                stats,
                unit,
                min_lengths,
                max_lengths,
            } => format!(
                "TXRNGR10005: The read lengths are incompatible with all features described in the \
                 feature reference for {unit}.\n\
                 {stats}\n\
                 The minimum read length for different feature types are:\n{min_lengths}\n\n\
                 We expect that at least 50% of the reads exceed the minimum length.\n\
                 {trim_lengths}\
                 {ERROR_CODE_INFO}",
                min_lengths = min_lengths
                    .iter()
                    .sorted_by_key(|&(feature_type, _min_lengths)| feature_type)
                    .flat_map(|(feature_type, min_lengths)| {
                        min_lengths
                            .iter()
                            .sorted_by_key(|&(&which_read, _)| which_read as usize)
                            .map(move |(which_read, length)| {
                                format!(
                                    "{:26} - {length}",
                                    format!("{}/{which_read}", feature_type.as_str())
                                )
                            })
                    })
                    .format("\n"),
                trim_lengths = if max_lengths.is_empty() {
                    String::new()
                } else {
                    format!(
                        "NOTE: You have specified the trim lengths explicitly:\n{}\n",
                        WHICH_LEN_READS
                            .iter()
                            .filter_map(|r| max_lengths.get(r).map(|l| format!("- {r}: {l} bases")))
                            .format("\n")
                    )
                },
            ),
            FeatureTypeNotInReference { feature_type } => format!(
                "Reads with feature type {} were observed, but this feature type is not described \
                 anywhere in the provided feature reference.",
                feature_type.as_str()
            ),
            ProbeBarcodeMixture { mixture, unit } => format!(
                "We detected multiple probe barcodes in: {}\n\
                 Singleplex Flex chemistry is invalid with >1 probe barcode. If this is a \
                 multiplex Flex library please include a [samples] section defining the inputs.\n\
                 The following top probe barcodes were observed: {}.\n",
                unit,
                mixture.iter().sorted().format(", ")
            ),
        };
        write!(f, "{msg}")
    }
}
