//! Martian stage DETECT_CHEMISTRY
#![deny(missing_docs)]

use crate::barcode_overlap::FRPGemBarcodeOverlapRow;
use crate::detect_chemistry::adapter_filter::preprocess_units;
use crate::detect_chemistry::chemistry_filter::{
    ChemistryFilter, DetectChemistryCandidates, DetectChemistryUnit, LoadedUnit, UnitsPerLibrary,
    detect_chemistry_units, named_candidates,
};
use crate::detect_chemistry::errors::DetectChemistryErrors;
use crate::detect_chemistry::identity_check::check_fastq_identity;
use crate::detect_chemistry::length_filter::LengthFilter;
use crate::detect_chemistry::mapping_filter::ReadMappingFilter;
use crate::detect_chemistry::probe_bc_check::validate_no_probe_bc_mixture_in_sfrp;
use crate::detect_chemistry::probe_bc_pairing::{
    detect_probe_barcode_pairing, should_detect_probe_barcode_pairing,
};
use crate::detect_chemistry::whitelist_filter::WhitelistMatchFilter;
use anyhow::{Context, Result, bail, ensure};
use barcode::whitelist::BarcodeId;
use cr_types::chemistry::adapter::trimmer_for_chemistry;
use cr_types::chemistry::{
    AutoChemistryName, AutoOrRefinedChemistry, ChemistryDef, ChemistryDefs, ChemistryName,
    ChemistrySpecs, normalize_chemistry_def,
};
use cr_types::probe_set::ProbeSetReference;
use cr_types::reference::feature_reference::{FeatureConfig, FeatureReferenceFile};
use cr_types::reference::probe_set_reference::TargetSetFile;
use cr_types::sample_def::SampleDef;
use cr_types::{AlignerParam, BarcodeMultiplexingType, LibraryType, ReadLevel};
use fastq_set::long_read_adapter_trimmer::{AdapterLayout, TrimmerConfig};
use itertools::Itertools;
use martian::prelude::*;
use martian_derive::{MartianStruct, make_mro};
use martian_filetypes::FileTypeWrite;
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::tabular_file::CsvFile;
use metric::{TxHashMap, TxHashSet, join_metric_name};
use multi::config::multiconst::DEFAULT_INCLUDE_INTRONS;
use multi::config::{MultiConfigCsv, MultiConfigCsvFile, ProbeBarcodeIterationMode, SamplesCsv};
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::default::Default;
use std::path::{Path, PathBuf};

const MIN_READS_NEEDED: usize = 10_000;
const DEFAULT_CHEMISTRY_SAMPLE_READS: usize = 100_000;
const DEFAULT_CHEMISTRY_TOTAL_READS: usize = 2_000_000;
const DEFAULT_MIN_MAJOR_PROBE_BC_FRAC: f64 = 0.7;

#[derive(Clone, Serialize, Deserialize, MartianStruct)]
#[cfg_attr(test, derive(Default))]
pub struct DetectChemistryStageInputs {
    pub sample_def: Vec<SampleDef>,
    pub reference_path: Option<PathBuf>,
    pub feature_reference: Option<FeatureReferenceFile>,
    pub chemistry_specs: ChemistrySpecs,
    pub allowed_chems: Option<Vec<AutoOrRefinedChemistry>>,
    pub r1_length: Option<usize>,
    pub r2_length: Option<usize>,
    pub multi_config: Option<MultiConfigCsvFile>,
    pub is_pd: bool,
    pub custom_chemistry_defs: ChemistryDefs,
    pub feature_config: Option<FeatureConfig>,
    pub target_set: Option<TargetSetFile>,
    pub is_spatial: bool,
    pub aligner: Option<AlignerParam>,
}

impl DetectChemistryStageInputs {
    fn to_metadata(&self) -> Metadata<'_> {
        Metadata {
            allowed_chems: self.allowed_chems.as_deref(),
            sample_defs: &self.sample_def,
            reference_path: self.reference_path.as_deref(),
            feature_reference: self.feature_reference.as_ref(),
            feature_config: self.feature_config.as_ref(),
            is_pd: self.is_pd,
            is_spatial: self.is_spatial,
        }
    }
}

/// A subset of args data that is used across the stage.
struct Metadata<'a> {
    pub allowed_chems: Option<&'a [AutoOrRefinedChemistry]>,
    pub sample_defs: &'a [SampleDef],
    pub reference_path: Option<&'a Path>,
    pub feature_reference: Option<&'a FeatureReferenceFile>,
    pub feature_config: Option<&'a FeatureConfig>,
    pub is_pd: bool,
    pub is_spatial: bool,
}

pub type DetectedProbeBarcodePairingFile = JsonFile<TxHashMap<BarcodeId, BarcodeId>>;

#[derive(Clone, Serialize, Deserialize, MartianStruct)]
#[cfg_attr(test, derive(Debug))]
pub struct DetectChemistryStageOutputs {
    pub chemistry_defs: ChemistryDefs,
    pub library_trimmers: TxHashMap<LibraryType, TrimmerConfig>,
    pub is_antibody_only: bool,
    #[mro_retain]
    pub probe_barcode_overlap: Option<CsvFile<FRPGemBarcodeOverlapRow>>,
    #[mro_retain]
    pub detected_probe_barcode_pairing: Option<DetectedProbeBarcodePairingFile>,
}

/// Martian stage DETECT_CHEMISTRY
pub struct DetectChemistry;

#[make_mro(mem_gb = 20, volatile = strict)]
impl MartianMain for DetectChemistry {
    type StageInputs = DetectChemistryStageInputs;
    type StageOutputs = DetectChemistryStageOutputs;

    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        // Bail out if any of chemistry_specs are not in the optional list of allowed
        // chemistry names.
        for spec in args.chemistry_specs.values() {
            validate_chemistry_spec(*spec, args.allowed_chems.as_deref())?;
        }

        let multi_config_csv = args
            .multi_config
            .as_ref()
            .map(|csv| csv.read().with_context(|| csv.display().to_string()))
            .transpose()?;

        let (chemistry_defs, units, selected_layouts) =
            select_chemistries(&args, multi_config_csv.as_ref())?;

        let is_antibody_only = args
            .sample_def
            .iter()
            .all(|x| x.library_type == Some(LibraryType::Antibody));

        let library_trimmers: TxHashMap<LibraryType, TrimmerConfig> = chemistry_defs
            .iter()
            .filter_map(|(&lib_type, chem)| {
                let params = units[&lib_type][0]
                    .0
                    .fastq_mode
                    .adapter_trimming_parameters(args.aligner)?;
                let layout = selected_layouts.get(&lib_type).cloned().unwrap_or_else(|| {
                    panic!("No empirically determined adapter layout found for {lib_type}")
                });
                Some((lib_type, chem, params, layout))
            })
            .map(|(lib_type, chem, params, layout)| {
                let Some(trimmer) = trimmer_for_chemistry(chem, lib_type, layout, &params)? else {
                    bail!(
                        "unable to construct a read trimmer for {} {lib_type}",
                        chem.name
                    );
                };
                anyhow::Ok((lib_type, trimmer))
            })
            .try_collect()?;

        let mut outputs = DetectChemistryStageOutputs {
            library_trimmers,
            chemistry_defs,
            is_antibody_only,
            probe_barcode_overlap: None,
            detected_probe_barcode_pairing: None,
        };

        if let Some(multi_config_csv) = &multi_config_csv {
            handle_probe_barcode_translation(&rover, multi_config_csv, &units, &mut outputs)?;
        }

        // Remove RNA length if we're doing gapfill.
        if let Some(target_set) = &args.target_set
            && let Some(gex_chem) = outputs.chemistry_defs.get_mut(&LibraryType::Gex)
            && ProbeSetReference::from_path(target_set, None, 0)?.has_gap_probes
        {
            gex_chem.rna.length = None;
        }

        Ok(outputs)
    }
}

type AdaptersPerLibrary = TxHashMap<LibraryType, AdapterLayout>;

/// Top-level entry point to define the map of chemistries per library type.
fn select_chemistries(
    args: &DetectChemistryStageInputs,
    multi_config_csv: Option<&MultiConfigCsv>,
) -> Result<(ChemistryDefs, UnitsPerLibrary, AdaptersPerLibrary)> {
    let units = detect_chemistry_units(&args.sample_def, args.r1_length, args.r2_length)?;
    println!(
        "Number of fastq units = {}",
        units.values().map(Vec::len).sum::<usize>()
    );

    let (units, selected_layouts) = load_reads(units, multi_config_csv)?;
    ensure!(!units.is_empty(), "no reads loaded");

    let metadata = args.to_metadata();
    // Handle manually-specified chemistry.
    let final_defs = match unpack_chemistry_specs(&args.chemistry_specs)? {
        UnpackedChemistrySpecs::Manual(chems) => use_manual_chemistries(
            &chems,
            &metadata,
            multi_config_csv,
            &args.custom_chemistry_defs,
            &units,
        ),
        UnpackedChemistrySpecs::Auto(mode) => {
            // -----------------------------------------------------------------------------------------
            // If the sample def contains VDJ, none of the other library types should be present and
            // the chemistry def should be `SCVDJ_auto`.
            if metadata
                .sample_defs
                .iter()
                .any(|sd| sd.library_type.is_some_and(|lt| lt.is_vdj()))
            {
                ensure!(
                    metadata
                        .sample_defs
                        .iter()
                        .all(|sd| sd.library_type.is_some_and(|lt| lt.is_vdj())),
                    "VDJ cannot be mixed with other library types"
                );
                ensure!(
                    mode == AutoChemistryName::Vdj,
                    "VDJ chemistry_name_spec should be SCVDJ_auto"
                );
            }

            let candidates = {
                let mut selected_mode_chems = mode.allowed_chemistries(
                    metadata.is_pd,
                    multi_config_csv.is_some_and(|x| x.samples.is_some()),
                    multi_config_csv.is_some_and(is_rtl_uncollapsed),
                );
                if let Some(allowed_chems) = metadata.allowed_chems {
                    let allowed_chems: TxHashSet<_> = allowed_chems
                        .iter()
                        .filter_map(AutoOrRefinedChemistry::refined)
                        .collect();
                    // Filter down the possible chemistries using the allowed list.
                    selected_mode_chems.retain(|chem| allowed_chems.contains(chem));
                }
                named_candidates(&selected_mode_chems)
            };

            let is_overhang_multiplexed =
                multi_config_csv.is_some_and(MultiConfigCsv::is_overhang_multiplexed);

            let chems: ChemistryDefs =
                detect_chemistry_per_library_type(&candidates, &metadata, &units)?
                    .into_iter()
                    .map(|(library_type, chem)| {
                        let chem = validate_chemistry(
                            chem,
                            metadata.allowed_chems,
                            multi_config_csv,
                            is_overhang_multiplexed,
                        )?;
                        anyhow::Ok((library_type, chem))
                    })
                    .try_collect()?;
            assert!(!chems.is_empty());

            if let Some(config) = multi_config_csv
                && let Some(gex_chem) = chems.get(&LibraryType::Gex)
                && gex_chem.name.is_rtl() == Some(false)
            {
                let gex = config.gene_expression.as_ref();
                ensure!(
                    gex.is_some_and(|gex| gex.reference_path.is_some()),
                    "A reference transcriptome is required to analyze a Gene Expression library."
                );
            }

            let names: Vec<_> = chems.values().map(|chem| chem.name).unique().collect();

            // Ensure that we produced compatible chemistries. For earlier
            // products, this implies that we should have the same chemistry for
            // each library. For newer products, make sure that all chemistries
            // are from a single product family. For now we check identity based
            // on chemistry name, we may want to be more stringent in the future.

            if names.len() == 1 {
                // A single chemistry worked for all libraries. Do a bit of
                // validation.
                if names[0].is_sfrp() {
                    // Extract min_major_probe_bc_frac from multi config or use default
                    let min_major_probe_bc_frac = multi_config_csv
                        .and_then(|config| config.gene_expression.as_ref())
                        .and_then(|gex_params| gex_params.min_major_probe_bc_frac)
                        .unwrap_or(DEFAULT_MIN_MAJOR_PROBE_BC_FRAC);

                    // Bail out if a mixture of probe barcodes is observed for singleplex FRP chemistry
                    validate_no_probe_bc_mixture_in_sfrp(
                        chems.values().next().unwrap(),
                        &units,
                        metadata.feature_reference,
                        metadata.feature_config,
                        min_major_probe_bc_frac,
                    )?;
                }
                Ok(chems)
            }
            // Otherwise, check that we found chemistries from compatible
            // product families.
            else if names.iter().all(|chem| chem.is_rtl().unwrap_or(false))
                || names.iter().all(ChemistryName::is_sc_3p_v3)
                || names.iter().all(ChemistryName::is_sc_3p_v4)
            {
                Ok(chems)
            } else {
                // Chemistries seem to be mutually incompatible.
                bail!(DetectChemistryErrors::PerLibraryConflictingChemistries(
                    chems
                ))
            }
        }
    }?;
    Ok((final_defs, units, selected_layouts))
}

/// Represent the possible valid ways of unpacking chemistry specs.
enum UnpackedChemistrySpecs {
    Auto(AutoChemistryName),
    Manual(HashMap<LibraryType, ChemistryName>),
}

/// Unpack the provided chemistry specs into either auto or manual modes.
/// Expect a single auto chemistry mode, or a refined chemistry for every lib.
fn unpack_chemistry_specs(chemistry_specs: &ChemistrySpecs) -> Result<UnpackedChemistrySpecs> {
    let mut auto_modes: Vec<_> = chemistry_specs
        .values()
        .filter_map(AutoOrRefinedChemistry::auto)
        .collect();
    if !auto_modes.is_empty() {
        ensure!(
            auto_modes.len() == chemistry_specs.len(),
            "mixing auto and refined chemistries for different library types is not supported"
        );
        auto_modes.sort();
        auto_modes.dedup();
        ensure!(
            auto_modes.len() == 1,
            "multiple conflicting auto modes provided: {}",
            auto_modes.iter().format(", ")
        );
        return Ok(UnpackedChemistrySpecs::Auto(auto_modes[0]));
    }
    Ok(UnpackedChemistrySpecs::Manual(
        chemistry_specs
            .iter()
            .map(|(lib_type, spec)| (*lib_type, spec.refined().unwrap()))
            .collect(),
    ))
}

/// Validate manual chemistries for all units.
///
/// If any libraries are specified as custom chemistry, use the provided custom
/// chemistry for them.  No validation is performed for custom chemistries.
fn use_manual_chemistries(
    chems: &HashMap<LibraryType, ChemistryName>,
    metadata: &Metadata<'_>,
    multi_config_csv: Option<&MultiConfigCsv>,
    custom_chemistry_defs: &ChemistryDefs,
    all_units: &UnitsPerLibrary,
) -> Result<ChemistryDefs> {
    all_units
        .iter()
        .map(|(&library_type, units)| {
            let chem_name = chems[&library_type];
            let chem_def = if chem_name == ChemistryName::Custom {
                use_custom_chemistry(custom_chemistry_defs.get(&library_type))?
            } else {
                use_manual_chemistry(library_type, chem_name, metadata, multi_config_csv, units)?
            };
            Ok((library_type, chem_def))
        })
        .try_collect()
}

/// If custom chemistry is specified, extract it from the args.
fn use_custom_chemistry(custom_chemistry_def: Option<&ChemistryDef>) -> Result<ChemistryDef> {
    let Some(custom_def) = custom_chemistry_def else {
        bail!(
            "custom chemistry def should be present if the input chemistry is '{}'",
            ChemistryName::Custom
        );
    };
    ensure!(
        custom_def.name == ChemistryName::Custom,
        "expected a custom chemistry but found {}",
        custom_def.name
    );
    Ok(normalize_chemistry_def(custom_def.clone()).unwrap_or_else(|| custom_def.clone()))
}

/// If manual chemistry is specified, extract and validate it.
///
/// Check that it passes the whitelist filter. Manual chemistry is used as an
/// escape hatch by customers for data which is usually of lower quality or by
/// QA to run non-standard fuzzed FASTQs. No minimum number of reads is enforced
/// here. Emit a warning if there are few valid barcodes.
fn use_manual_chemistry(
    library_type: LibraryType,
    chem: ChemistryName,
    metadata: &Metadata<'_>,
    multi_config_csv: Option<&MultiConfigCsv>,
    units: &[LoadedUnit],
) -> Result<ChemistryDef> {
    let mut chem = ChemistryDef::named(chem);
    if library_type == LibraryType::Crispr
        && let Some(crispr_probe_barcode_location) = multi_config_csv
            .and_then(|config| config.feature.as_ref())
            .and_then(|feature| feature.crispr_probe_barcode_location)
    {
        // Config validation should have made sure of this.
        assert!(chem.name.is_mfrp());
        chem.override_probe_barcode_extraction(crispr_probe_barcode_location);
    }
    let chems: TxHashMap<_, _> = [(chem.name, chem.clone())].into_iter().collect();

    if chem.name.is_mfrp() {
        // Check the read length to ensure that the probe barcode is sequenced.
        let _matching_chemistries = LengthFilter::new(
            library_type,
            metadata.feature_reference,
            metadata.feature_config,
        )?
        .filter_chemistries(&chems, units)?;
    }

    if let Err(err) =
        WhitelistMatchFilter::new(metadata.is_spatial).filter_chemistries(&chems, units)
    {
        if metadata.is_spatial || matches!(err, DetectChemistryErrors::WhitelistLoadError(_)) {
            bail!(err);
        }
        println!("WARNING: {err:#}");
    }

    validate_rtl(multi_config_csv, &chem)?;
    Ok(chem)
}

type CandidatesPerLibrary = HashMap<LibraryType, DetectChemistryCandidates>;

/// Load a subsample of reads from each input unit.
///
/// Potentially perform post-processing of the reads if necessary for the
/// sequencing platform in use. Returns the loaded units and the empirically
/// selected adapter layout per library type (for libraries that require adapter
/// trimming).
fn load_reads(
    units: HashMap<LibraryType, Vec<DetectChemistryUnit>>,
    multi_config_csv: Option<&MultiConfigCsv>,
) -> Result<(UnitsPerLibrary, AdaptersPerLibrary)> {
    // Check for duplicate R1 and R2 files amongst units.
    check_fastq_identity(units.values().flatten())?;

    // Get chemistry detection parameters from multi config or use defaults
    let sample_reads = multi_config_csv
        .as_ref()
        .and_then(|config| config.gene_expression.as_ref())
        .and_then(|gex_params| gex_params.detect_chemistry_sample_reads)
        .unwrap_or(DEFAULT_CHEMISTRY_SAMPLE_READS);

    let total_reads = multi_config_csv
        .as_ref()
        .and_then(|config| config.gene_expression.as_ref())
        .and_then(|gex_params| gex_params.detect_chemistry_total_reads)
        .unwrap_or(DEFAULT_CHEMISTRY_TOTAL_READS);

    let combined: Vec<(LibraryType, Vec<LoadedUnit>, Option<AdapterLayout>)> = units
        .into_iter()
        .sorted_by_key(|(lib_type, _)| *lib_type)
        .map(|(lib_type, units)| {
            let loaded: Vec<LoadedUnit> = units
                .iter()
                .map(|unit| {
                    println!("Sampling reads from: {unit}");
                    let reads = unit.sampled_read_pairs(sample_reads, total_reads, None)?;
                    anyhow::Ok((unit.clone(), reads))
                })
                .try_collect()?;

            // STAR aligner is used by the downstream ReadMappingFilter.
            let Some(trimmer_params) = units
                .iter()
                .map(|u| {
                    u.fastq_mode
                        .adapter_trimming_parameters(Some(AlignerParam::Star))
                })
                .unique()
                .exactly_one()
                .unwrap()
            else {
                return anyhow::Ok((lib_type, loaded, None));
            };

            // We need to do adapter trimming on these units.
            let (reloaded, selected_layout) =
                preprocess_units(loaded, &trimmer_params, lib_type, sample_reads, total_reads)?;

            anyhow::Ok((lib_type, reloaded, Some(selected_layout)))
        })
        .try_collect()?;

    let mut units_per_library = UnitsPerLibrary::new();
    let mut selected_layouts = AdaptersPerLibrary::default();
    for (lib_type, units, layout) in combined {
        units_per_library.insert(lib_type, units);
        if let Some(layout) = layout {
            selected_layouts.insert(lib_type, layout);
        }
    }
    Ok((units_per_library, selected_layouts))
}

/// Determine which chemistry to use for all units of a single library type.
///
/// For most cases this will resolve to a single compatible chemistry.
/// For certain cases such as Antigen Capture, we cannot select a unique
/// chemistry, and we require a follow-on consensus step to resolve to a single
/// chemistry.
fn detect_chemistry_for_library(
    library_type: LibraryType,
    candidates: &DetectChemistryCandidates,
    metadata: &Metadata<'_>,
    units: &[LoadedUnit],
) -> Result<DetectChemistryCandidates> {
    #[allow(clippy::enum_glob_use)]
    use ChemistryName::*;

    for (unit, read_pairs) in units {
        ensure!(
            read_pairs.len() >= MIN_READS_NEEDED,
            DetectChemistryErrors::NotEnoughReads {
                num_reads: read_pairs.len(),
                min_reads: MIN_READS_NEEDED,
                unit: Box::new(unit.clone()),
            }
        );
    }

    println!(
        "Potential chemistries: {}",
        candidates.keys().sorted().format(", ")
    );

    // Read length based filtering.
    let candidates = LengthFilter::new(
        library_type,
        metadata.feature_reference,
        metadata.feature_config,
    )?
    .filter_chemistries(candidates, units)?;

    println!(
        "Candidates after length filter: {}",
        candidates.keys().sorted().format(", ")
    );

    // -----------------------------------------------------------------------------------------
    // For each unit of fastqs, ensure that a sufficient fraction of reads contain barcodes
    // which match the whitelist for at least one of the possible chemistries.
    // TODO: we should share the whitelist cache between different invocations of
    // this function to avoid repeatedly loading them for each library type.
    let mut candidates =
        WhitelistMatchFilter::new(metadata.is_spatial).filter_chemistries(&candidates, units)?;

    // If we're down to one candidate, we're done.
    if candidates.len() == 1 {
        return Ok(candidates);
    }

    // Downselect to regular throughput kits if HT and regular are both compatible.
    if candidates_are(&candidates, &[ThreePrimeV3PolyA, ThreePrimeV3HTPolyA]) {
        candidates.remove(&ThreePrimeV3HTPolyA).unwrap();
        return Ok(candidates);
    } else if candidates_are(&candidates, &[ThreePrimeV3CS1, ThreePrimeV3HTCS1]) {
        candidates.remove(&ThreePrimeV3HTCS1).unwrap();
        return Ok(candidates);
    } else if candidates_are(&candidates, &[FivePrimeR2, FivePrimeHT]) {
        candidates.remove(&FivePrimeHT).unwrap();
        return Ok(candidates);
    }

    println!(
        "After whitelist filter: {}",
        candidates.keys().sorted().format(", ")
    );

    let is_antibody_only = metadata
        .sample_defs
        .iter()
        .all(|sd| sd.library_type == Some(LibraryType::Antibody));
    let candidates = if is_antibody_only {
        let mut result = candidates;
        // we define a new chemistry named "SC-FB" because this could be 3' v2 polyA capture
        // antibody library or a 5' antibody library
        let redundant_ab_chems = [ThreePrimeV2, FivePrimeR2, FivePrimePE];

        if redundant_ab_chems
            .iter()
            .any(|chem| result.contains_key(chem))
        {
            // FIXME: this is a hole in ensuring that chemistry defs flow through
            // Any state updated on the input defs is lost here... may not be important
            // if this is only for old products.
            result.insert(FeatureBarcodingOnly, FeatureBarcodingOnly.def());
            for chem in &redundant_ab_chems {
                result.remove(chem);
            }
        }
        result
    } else if library_type.is_gex() {
        ReadMappingFilter::new(metadata.reference_path.unwrap())?
            .filter_chemistries(&candidates, units)?
    } else {
        candidates
    };

    println!(
        "After mapping filter: {}",
        candidates.keys().sorted().format(", ")
    );

    if candidates.len() != 1 {
        println!(
            "Could not distinguish between {} for {library_type}",
            candidates.keys().sorted().format(", ")
        );
    }

    Ok(candidates)
}

/// Return true if the candidates are an exact match to the provided set of names.
fn candidates_are(candidates: &DetectChemistryCandidates, are: &[ChemistryName]) -> bool {
    if candidates.len() != are.len() {
        return false;
    }
    are.iter().all(|name| candidates.contains_key(name))
}

/// Run chemistry detection on every library independently.
fn detect_chemistry_per_library_type(
    candidates: &DetectChemistryCandidates,
    metadata: &Metadata<'_>,
    all_units: &UnitsPerLibrary,
) -> Result<ChemistryDefs> {
    let candidates_per_library: CandidatesPerLibrary = all_units
        .iter()
        .sorted_by_key(|&(&lib_type, _)| lib_type)
        .map(|(&library_type, units)| {
            println!("\nDetecting chemistry for {library_type}");
            let candidates =
                detect_chemistry_for_library(library_type, candidates, metadata, units)?;

            println!(
                "\nDetected {} for {library_type}",
                candidates.keys().format(", ")
            );
            anyhow::Ok((library_type, candidates))
        })
        .try_collect()?;

    // If we have one chemistry per library type, we're done.
    if candidates_per_library
        .values()
        .all(|chems| chems.len() == 1)
    {
        return Ok(candidates_per_library
            .into_iter()
            .map(|(library_type, chems)| (library_type, chems.into_values().exactly_one().unwrap()))
            .collect());
    }

    // If we can find a consensus chemistry that works for all libraries, use it.
    // Note that we take the consensus based on chemistry name but then retain
    // the individual defs for each library, in case a library-specific customization
    // was applied.
    let consensus: Vec<_> = candidates_per_library
        .values()
        .map(|candidates| {
            candidates
                .keys()
                .copied()
                .collect::<HashSet<ChemistryName>>()
        })
        .reduce(|consensus, names| &consensus & &names)
        .unwrap()
        .into_iter()
        .sorted()
        .collect();
    match consensus.as_slice() {
        [] => {
            // Arbitrarily report one library with multiple chems to fail with.
            let (lib_type, candidates) = candidates_per_library
                .into_iter()
                .sorted_by_key(|(lib_type, _)| *lib_type)
                .find(|(_, candidates)| candidates.len() > 1)
                .unwrap();
            bail!(
                "Could not distinguish between {} for {lib_type}",
                candidates.keys().sorted().format(", ")
            );
        }
        [consensus] => Ok(candidates_per_library
            .into_iter()
            .map(|(library_type, mut chems)| (library_type, chems.remove(consensus).unwrap()))
            .collect()),
        too_many => {
            bail!(
                "Could not distinguish between {}",
                too_many.iter().format(", ")
            );
        }
    }
}

/// Validate that the specified chemistry is allowed.
fn validate_chemistry_spec(
    chemistry_spec: AutoOrRefinedChemistry,
    allowed_chems: Option<&[AutoOrRefinedChemistry]>,
) -> Result<()> {
    let Some(allowed_chems) = allowed_chems else {
        return Ok(());
    };

    ensure!(
        allowed_chems.contains(&chemistry_spec),
        DetectChemistryErrors::ChemistryNotAllowed {
            input: chemistry_spec,
            allowed: allowed_chems.to_vec(),
        }
    );

    Ok(())
}

/// Validate a chemistry, and do some post-processing.
fn validate_chemistry(
    chemistry: ChemistryDef,
    allowed_chems: Option<&[AutoOrRefinedChemistry]>,
    multi_config_csv: Option<&MultiConfigCsv>,
    overhang_multiplexing: bool,
) -> Result<ChemistryDef> {
    if chemistry.name == ChemistryName::ArcV1 {
        bail!(
            "Cell Ranger detected the chemistry {}, which may indicate a workflow \
            error during sample preparation. Please check the reagents used to prepare this \
            sample and contact 10x Genomics support for further assistance. If this workflow is \
            intentional, you can force Cell Ranger to process this data by manually specifying \
            the {} chemistry in your analysis configuration.",
            chemistry.name,
            chemistry.name
        );
    }

    validate_rtl(multi_config_csv, &chemistry)?;

    let chemistry_with_overhang = if overhang_multiplexing {
        // FIXME: this will discard any dynamic changes made to the chemistry
        ChemistryDef::named(chemistry.name.get_overhang_version()?)
    } else {
        chemistry
    };

    validate_chemistry_spec(
        AutoOrRefinedChemistry::Refined(chemistry_with_overhang.name),
        allowed_chems,
    )?;

    Ok(chemistry_with_overhang)
}

/// Validate the chemistry with RTL-related parameters.
fn validate_rtl(multi_config_csv: Option<&MultiConfigCsv>, chemistry: &ChemistryDef) -> Result<()> {
    let Some(config) = multi_config_csv else {
        // cellranger count does not support RTL chemistries.
        return Ok(());
    };

    if !chemistry.name.is_rtl().expect("is_rtl is None only for custom and spatial chemistries, which do not use chemistry detection") {
        let Some(samples) = &config.samples else {
            return Ok(());
        };
        ensure!(
            samples.multiplexing_type!=BarcodeMultiplexingType::ReadLevel(ReadLevel::RTL),
            "A non-Flex chemistry {} was detected, and the [samples] \
             section has a probe_barcode_ids column. The probe_barcode_ids column may only be \
             specified with Flex chemistries.",
             chemistry.name
        );
        return Ok(());
    }

    if let Some(gex) = &config.gene_expression {
        if config.libraries.has_gene_expression() {
            ensure!(gex.has_probe_set(), "Flex chemistries require a probe-set.");
        }
        ensure!(
            gex.include_introns == DEFAULT_INCLUDE_INTRONS,
            "The [gene-expression] section specifies the parameter include-introns, \
             which is not valid for Flex chemistries."
        );
    }

    if chemistry.name.is_sfrp() {
        ensure!(
            config.samples.is_none(),
            "We detected singleplex Flex chemistry from the data. \
             Sample definitions are unsupported for singleplex inputs. \
             To process this data as multiplex Flex you will need to specify `MFRP` \
             as the chemistry in the config.csv."
        );
    }

    Ok(())
}

/// Return true if probe barcode IDs are uncollapsed, like BC001A. For PD use.
fn is_rtl_uncollapsed(multi_config_csv: &MultiConfigCsv) -> bool {
    multi_config_csv
        .sample_barcode_ids_used_in_experiment(ProbeBarcodeIterationMode::All)
        .into_iter()
        .all(|x| x.ends_with(['A', 'B', 'C', 'D']))
}

/// Identify probe barcode pairings if necessary.
/// If pairings are present, modify the translation whitelist to include those pairings.
fn handle_probe_barcode_translation(
    rover: &MartianRover,
    multi_config_csv: &MultiConfigCsv,
    units: &UnitsPerLibrary,
    outputs: &mut DetectChemistryStageOutputs,
) -> Result<()> {
    let detected_pairing_id_translation = if !should_detect_probe_barcode_pairing(multi_config_csv)
    {
        None
    } else {
        let (overlaps, pairings) = detect_probe_barcode_pairing(&outputs.chemistry_defs, units)?;

        let probe_barcode_overlap: CsvFile<_> = rover.make_path("probe_barcode_overlap");
        probe_barcode_overlap.write(&overlaps)?;
        outputs.probe_barcode_overlap = Some(probe_barcode_overlap);

        let pairings = pairings.into_iter().collect();
        let detected_probe_barcode_pairing: JsonFile<_> =
            rover.make_path("detected_probe_barcode_pairing");
        detected_probe_barcode_pairing.write(&pairings)?;
        outputs.detected_probe_barcode_pairing = Some(detected_probe_barcode_pairing);

        if pairings.is_empty() {
            None
        } else {
            // The pairings are generated as RTL: AB; invert the mapping so they
            // serve to translate an AB barcode into the paired RTL barcode.
            Some(pairings.into_iter().map(|(rtl, ab)| (ab, rtl)).collect())
        }
    };

    let explicit_pairing = multi_config_csv
        .samples
        .as_ref()
        .map(SamplesCsv::get_translated_probe_barcodes)
        .and_then(|pairing| {
            if pairing.is_empty() {
                None
            } else {
                Some(pairing)
            }
        });

    // Remap the probe barcode whitelist(s) if we have a probe barcode pairing.
    // Use the explicit pairing if one is provided.
    // Otherwise, use the detected pairing if there are any.
    if let Some(pairing) = explicit_pairing.or(detected_pairing_id_translation) {
        // Use the GEX probe barcode whitelist as the target.
        let target_probe_bc_whitelist = outputs.chemistry_defs[&LibraryType::Gex]
            .barcode_whitelist_source()?
            .probe();
        for (&library_type, chemistry_def) in &mut outputs.chemistry_defs {
            chemistry_def.translate_probe_barcode_whitelist_with_id_map(
                &pairing,
                &target_probe_bc_whitelist,
                rover.make_path(join_metric_name(
                    library_type,
                    "probe_barcode_translation_whitelist.tsv",
                )),
            )?;
        }
    };
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use cr_types::chemistry::adapter::{layouts, trimmer_for_chemistry};
    use cr_types::chemistry::{ChemistryDefsExt, known_chemistry_defs};
    use cr_types::reference::probe_set_reference::TargetSetFile;
    use cr_types::sample_def::{FastqMode, SampleDefBuilder};
    use cr_types::{BarcodeMultiplexingType, LibraryType, ReadLevel, SequencingPlatform};
    use dui_tests::stage_test::StageFailTest;
    use dui_tests::{DuiTest, stage_fail_dui_test};
    use metric::set;
    use multi::config::{
        DataType, GeneExpressionParams, Lanes, LibrariesCsv, Library, SampleBarcodeColumn,
        SampleRow,
    };
    use std::io::Write;
    use std::path::Path;
    use std::string::ToString;

    #[test]
    fn test_validate_rtl() {
        use ChemistryName::{MFRP_Ab, MFRP_RNA, SFRP, ThreePrimeV3CS1, ThreePrimeV3PolyA};

        assert!(validate_rtl(None, &ThreePrimeV3PolyA.def()).is_ok());

        let gex_library = Library {
            read_path: PathBuf::default(),
            physical_library_id: String::default(),
            feature_type: LibraryType::Gex,
            subsample_rate: None,
            chemistry: None,
            data_type: DataType::BclProcessor {
                sample_indices: Vec::default(),
                lanes: Lanes::Any,
            },
            platform: SequencingPlatform::Illumina,
        };

        let antibody_library = Library {
            read_path: PathBuf::default(),
            physical_library_id: String::default(),
            feature_type: LibraryType::Antibody,
            subsample_rate: None,
            chemistry: None,
            data_type: DataType::BclProcessor {
                sample_indices: Vec::default(),
                lanes: Lanes::Any,
            },
            platform: SequencingPlatform::Illumina,
        };

        let gex_without_probe_set = MultiConfigCsv {
            libraries: LibrariesCsv::new(vec![gex_library]),
            gene_expression: Some(GeneExpressionParams {
                probe_set: vec![],
                include_introns: DEFAULT_INCLUDE_INTRONS,
                ..GeneExpressionParams::default()
            }),
            feature: None,
            vdj: None,
            samples: None,
            antigen_specificity: None,
            functional_map: None,
        };
        assert!(validate_rtl(Some(&gex_without_probe_set), &MFRP_RNA.def()).is_err());
        assert!(validate_rtl(Some(&gex_without_probe_set), &SFRP.def()).is_err());
        assert!(validate_rtl(Some(&gex_without_probe_set), &ThreePrimeV3PolyA.def()).is_ok());

        let gex_with_probe_set = MultiConfigCsv {
            gene_expression: Some(GeneExpressionParams {
                probe_set: vec![TargetSetFile::from("probe_set")],
                include_introns: DEFAULT_INCLUDE_INTRONS,
                ..GeneExpressionParams::default()
            }),
            ..gex_without_probe_set
        };
        assert!(validate_rtl(Some(&gex_with_probe_set), &MFRP_RNA.def()).is_ok());
        assert!(validate_rtl(Some(&gex_with_probe_set), &SFRP.def()).is_ok());
        assert!(validate_rtl(Some(&gex_with_probe_set), &ThreePrimeV3PolyA.def()).is_ok()); // to aggr with RTL

        let antibody_only = MultiConfigCsv {
            libraries: LibrariesCsv::new(vec![antibody_library]),
            gene_expression: Some(GeneExpressionParams {
                probe_set: vec![],
                include_introns: DEFAULT_INCLUDE_INTRONS,
                ..GeneExpressionParams::default()
            }),
            feature: None,
            vdj: None,
            samples: None,
            antigen_specificity: None,
            functional_map: None,
        };
        // This assertion might seem non-sensical, but the validation function
        // under test doesn't care about the library type, so it should not fail
        // for Gex even if the config is for antibody-only.
        assert!(validate_rtl(Some(&antibody_only), &MFRP_RNA.def()).is_ok());
        assert!(validate_rtl(Some(&antibody_only), &MFRP_Ab.def()).is_ok());
        assert!(validate_rtl(Some(&antibody_only), &SFRP.def()).is_ok());
        assert!(validate_rtl(Some(&antibody_only), &ThreePrimeV3CS1.def()).is_ok());

        let include_introns_false = MultiConfigCsv {
            gene_expression: Some(GeneExpressionParams {
                probe_set: vec![TargetSetFile::from("probe_set")],
                include_introns: false,
                ..GeneExpressionParams::default()
            }),
            libraries: LibrariesCsv::new(vec![]),
            feature: None,
            vdj: None,
            samples: None,
            antigen_specificity: None,
            functional_map: None,
        };
        assert!(validate_rtl(Some(&include_introns_false), &MFRP_RNA.def()).is_err());
        assert!(validate_rtl(Some(&include_introns_false), &SFRP.def()).is_err());
        assert!(validate_rtl(Some(&include_introns_false), &ThreePrimeV3CS1.def()).is_ok()); // to aggr with RTL

        let samples = MultiConfigCsv {
            samples: Some(SamplesCsv {
                rows: vec![SampleRow::from_probe_barcode_id("BC001")],
                multiplexing_type: BarcodeMultiplexingType::ReadLevel(ReadLevel::RTL),
                sample_barcode_column_name: SampleBarcodeColumn::ProbeBarcodeIds,
            }),
            ..gex_with_probe_set
        };
        assert!(validate_rtl(Some(&samples), &MFRP_RNA.def()).is_ok());
        assert!(validate_rtl(Some(&samples), &SFRP.def()).is_err());
        assert!(validate_rtl(Some(&samples), &ThreePrimeV3CS1.def()).is_err());
    }

    macro_rules! err_snapshot {
        ($val: expr) => {{
            let outs = DetectChemistry.test_run_tmpdir($val);
            if let Err(err) = &outs {
                let e = format!("{:#}", err);
                insta::assert_snapshot!(e);
            };
            assert!(outs.is_err());
        }};
    }

    fn specs_input_gex(chem: impl Into<AutoOrRefinedChemistry>) -> ChemistrySpecs {
        specs_input(LibraryType::Gex, chem.into())
    }

    fn specs_input(lib_type: LibraryType, chem: AutoOrRefinedChemistry) -> ChemistrySpecs {
        [(lib_type, chem)].into_iter().collect()
    }

    fn test_invalid_bcl2fastq_input(
        fastq_path: impl AsRef<Path>,
        lanes: Option<Vec<usize>>,
        sample_names: Vec<&str>,
    ) -> DetectChemistryStageInputs {
        let mut builder =
            SampleDefBuilder::new(FastqMode::ILMN_BCL2FASTQ, fastq_path.as_ref().to_path_buf())
                .sample_names(sample_names);
        if let Some(lanes) = lanes {
            builder = builder.lanes(lanes);
        }
        DetectChemistryStageInputs {
            sample_def: vec![builder.build()],
            chemistry_specs: specs_input_gex(AutoChemistryName::Count),
            ..DetectChemistryStageInputs::default()
        }
    }

    // ---------------------------------------------------------------------------------------------
    // Some of these tests are also done in the context of a sample def.

    #[test]
    fn test_missing_only_sample() {
        err_snapshot!(test_invalid_bcl2fastq_input(
            "../dui_tests/test_resources/cellranger-count/too_few_reads",
            None,
            vec!["missing_sample"],
        ));
    }

    #[test]
    fn test_missing_lane() {
        err_snapshot!(test_invalid_bcl2fastq_input(
            "../dui_tests/test_resources/cellranger-count/too_few_reads",
            Some(vec![2]),
            vec!["test_sample"],
        ));
    }

    #[test]
    fn test_missing_one_sample() {
        err_snapshot!(test_invalid_bcl2fastq_input(
            "../dui_tests/test_resources/cellranger-count/too_few_reads",
            None,
            vec!["test_sample", "missing_sample"],
        ));
    }

    #[test]
    fn test_invalid_folder() {
        err_snapshot!(test_invalid_bcl2fastq_input(
            "/no/such/folder",
            None,
            vec!["test_sample"]
        ));
    }

    #[test]
    fn test_empty_folder() {
        err_snapshot!(test_invalid_bcl2fastq_input(
            "../dui_tests/test_resources/cellranger-count/fastqs_empty",
            None,
            vec!["test_sample"],
        ));
    }

    #[test]
    fn test_imbalanced_fastq() {
        err_snapshot!(test_invalid_bcl2fastq_input(
            "../dui_tests/test_resources/cellranger-count/fastqs_imbalanced",
            None,
            vec!["test_sample"],
        ));
    }

    #[test]
    fn test_fastqs_invalid_chars() {
        err_snapshot!(test_invalid_bcl2fastq_input(
            "../dui_tests/test_resources/cellranger-count/fastqs_invalid_chars",
            None,
            vec!["test_sample"],
        ));
    }

    #[test]
    fn test_fastqs_incomplete() {
        err_snapshot!(test_invalid_bcl2fastq_input(
            "../dui_tests/test_resources/cellranger-count/fastqs_incomplete",
            None,
            vec!["test_sample"],
        ));
    }

    #[test]
    fn test_mismatched_header() {
        err_snapshot!(test_invalid_bcl2fastq_input(
            "../dui_tests/test_resources/cellranger-count/fastqs_mismatched_header",
            None,
            vec!["test_sample"],
        ));
    }

    #[test]
    fn test_duplicate_fastqs() {
        err_snapshot!(test_invalid_bcl2fastq_input(
            "../dui_tests/test_resources/cellranger-count/fastqs_duplicate_inputs",
            None,
            vec!["test_sample"],
        ));
    }

    #[test]
    fn test_duplicate_fastqs_recursive() {
        err_snapshot!(test_invalid_bcl2fastq_input(
            "../cr_lib/test/fastq/duplicate_fastqs_recursive",
            None,
            vec!["test_sample"],
        ));
    }

    #[test]
    fn test_duplicate_r1_r2() {
        err_snapshot!(test_invalid_bcl2fastq_input(
            "../dui_tests/test_resources/cellranger-count/fastqs_duplicate_r1_r2",
            None,
            vec!["test_sample"],
        ));
    }

    #[test]
    fn test_corrupted_gzip() {
        err_snapshot!(test_invalid_bcl2fastq_input(
            "../dui_tests/test_resources/cellranger-count/corrupted_gzip",
            None,
            vec!["corrupted"],
        ));
    }

    #[test]
    fn test_too_few_reads() {
        err_snapshot!(test_invalid_bcl2fastq_input(
            "../dui_tests/test_resources/cellranger-count/too_few_reads",
            None,
            vec!["test_sample"],
        ));
    }

    #[test]
    fn test_name_spec_not_in_allowed() {
        let args = DetectChemistryStageInputs {
            chemistry_specs: specs_input_gex(ChemistryName::VdjPE),
            allowed_chems: Some(vec![ChemistryName::ThreePrimeV3PolyA.into()]),
            ..Default::default()
        };
        err_snapshot!(args);
    }

    #[test]
    fn test_custom_chemistry_spec() {
        let mut def: ChemistryDef = ChemistryDef::named(ChemistryName::ThreePrimeV3PolyA);
        def.name = ChemistryName::Custom;
        def.rna.offset = 13;
        let args = DetectChemistryStageInputs {
            chemistry_specs: specs_input_gex(ChemistryName::Custom),
            custom_chemistry_defs: HashMap::from_iter([(LibraryType::Gex, def.clone())]),
            sample_def: vec![
                SampleDefBuilder::new(
                    FastqMode::ILMN_BCL2FASTQ,
                    "../dui_tests/test_resources/cellranger-count/CR-300-01_SC3pv3_15k".into(),
                )
                .sample_names(vec!["CR-300-01"])
                .build(),
            ],
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(outs.chemistry_defs.values().exactly_one().unwrap(), &def);
    }

    #[test]
    fn test_mix_of_custom_and_manual_chemistry_spec() {
        let mut def: ChemistryDef = ChemistryDef::named(ChemistryName::ThreePrimeV3CS1);
        def.name = ChemistryName::Custom;
        def.rna.offset = 13;
        let args = DetectChemistryStageInputs {
            chemistry_specs: [
                (
                    LibraryType::Gex,
                    AutoOrRefinedChemistry::Refined(ChemistryName::MFRP_RNA),
                ),
                (
                    LibraryType::Antibody,
                    AutoOrRefinedChemistry::Refined(ChemistryName::Custom),
                ),
            ]
            .into_iter()
            .collect(),
            custom_chemistry_defs: HashMap::from_iter([(LibraryType::Antibody, def.clone())]),
            sample_def: vec![
                SampleDefBuilder::new(
                    FastqMode::ILMN_BCL2FASTQ,
                    "../dui_tests/test_resources/cellranger-count/CR-300-01_SC3pv3_15k".into(),
                )
                .sample_names(vec!["CR-300-01"])
                .build(),
                SampleDefBuilder::new(
                    FastqMode::ILMN_BCL2FASTQ,
                    "../dui_tests/test_resources/cellranger-count/pbmc_1k_protein_v3_antibody"
                        .into(),
                )
                .sample_names(vec!["pbmc_1k_protein_v3_antibody"])
                .library_type(LibraryType::Antibody)
                .build(),
            ],
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(
            outs.chemistry_defs[&LibraryType::Gex].name,
            ChemistryName::MFRP_RNA
        );
        assert_eq!(outs.chemistry_defs[&LibraryType::Antibody], def);
    }

    #[test]
    fn test_custom_chemistry_changed_to_known_chemistry() {
        // First determine all of the unique known chemistries.
        let unique_known_defs: Vec<_> = known_chemistry_defs()
            .values()
            .cloned()
            .filter_map(normalize_chemistry_def)
            .collect();
        assert!(!unique_known_defs.is_empty());
        for known_def in unique_known_defs {
            let mut custom_def = known_def.clone();
            custom_def.name = ChemistryName::Custom;
            custom_def.description = "I'm so custom I'm not even custom at all".to_string();

            let args = DetectChemistryStageInputs {
                chemistry_specs: specs_input_gex(ChemistryName::Custom),
                custom_chemistry_defs: HashMap::from_iter([(LibraryType::Gex, custom_def)]),
                sample_def: vec![
                    SampleDefBuilder::new(
                        FastqMode::ILMN_BCL2FASTQ,
                        "../dui_tests/test_resources/cellranger-count/CR-300-01_SC3pv3_15k".into(),
                    )
                    .sample_names(vec!["CR-300-01"])
                    .build(),
                ],
                ..Default::default()
            };
            let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
            // Ensure that the output is mapped to a known chemistry
            assert_eq!(
                outs.chemistry_defs.values().exactly_one().unwrap(),
                &known_def
            );
        }
    }

    /// User supplies an incorrect option to `chemistry`, and it runs to completion.
    /// We supply `chemistry=SC3Pv2` for a datasets which is `SC3Pv3`.
    #[test]
    fn test_incorrect_chemistry_spec() {
        let args = DetectChemistryStageInputs {
            sample_def: vec![
                SampleDefBuilder::new(
                    FastqMode::ILMN_BCL2FASTQ,
                    "../dui_tests/test_resources/cellranger-count/CR-300-01_SC3pv3_15k".into(),
                )
                .sample_names(vec!["CR-300-01"])
                .build(),
            ],
            chemistry_specs: specs_input_gex(ChemistryName::ThreePrimeV2),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(
            outs.chemistry_defs.primary().name,
            ChemistryName::ThreePrimeV2
        );
    }

    #[test]
    fn test_correct_chemistry_spec() {
        let args = DetectChemistryStageInputs {
            sample_def: vec![
                SampleDefBuilder::new(
                    FastqMode::ILMN_BCL2FASTQ,
                    "../dui_tests/test_resources/cellranger-count/CR-300-01_SC3pv3_15k".into(),
                )
                .sample_names(vec!["CR-300-01"])
                .build(),
            ],
            chemistry_specs: specs_input_gex(ChemistryName::ThreePrimeV3PolyA),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(
            outs.chemistry_defs.primary().name,
            ChemistryName::ThreePrimeV3PolyA
        );
        assert!(!outs.is_antibody_only);
    }

    #[test]
    fn test_auto_chemistry_spec_3pv3() {
        let args = DetectChemistryStageInputs {
            sample_def: vec![
                SampleDefBuilder::new(
                    FastqMode::ILMN_BCL2FASTQ,
                    "../dui_tests/test_resources/cellranger-count/CR-300-01_SC3pv3_15k".into(),
                )
                .sample_names(vec!["CR-300-01"])
                .build(),
            ],
            chemistry_specs: specs_input_gex(AutoChemistryName::Count),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(
            outs.chemistry_defs.primary().name,
            ChemistryName::ThreePrimeV3PolyA
        );
        assert!(!outs.is_antibody_only);
    }

    #[test]
    fn test_correct_chemistry_spec_3pv3ht() {
        // should detect ThreePrimeV3HT when SCPV3HT is specified
        let args = DetectChemistryStageInputs {
            sample_def: vec![
                SampleDefBuilder::new(
                    FastqMode::ILMN_BCL2FASTQ,
                    "../dui_tests/test_resources/cellranger-count/HT_chemistry".into(),
                )
                .sample_names(vec!["1107522"])
                .build(),
            ],
            chemistry_specs: specs_input_gex(ChemistryName::ThreePrimeV3HTPolyA),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(
            outs.chemistry_defs.primary().name,
            ChemistryName::ThreePrimeV3HTPolyA
        );
        assert!(!outs.is_antibody_only);
    }

    #[test]
    fn test_auto_chemistry_spec_ht() {
        // should detect ThreePrimeV3PolyA when HT fastqs are run
        let args = DetectChemistryStageInputs {
            sample_def: vec![
                SampleDefBuilder::new(
                    FastqMode::ILMN_BCL2FASTQ,
                    "../dui_tests/test_resources/cellranger-count/HT_chemistry".into(),
                )
                .sample_names(vec!["1107522"])
                .build(),
            ],
            chemistry_specs: specs_input_gex(AutoChemistryName::Count),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(
            outs.chemistry_defs.primary().name,
            ChemistryName::ThreePrimeV3PolyA
        );
        assert!(!outs.is_antibody_only);
    }

    #[test]
    fn test_mixed_chemistries_3pv3_5p() -> Result<()> {
        stage_fail_dui_test!(
            MixedChemistries3p5p,
            description: "User inputs two sets of fastqs, one is 3' v3, the other is 5'",
            stage: DetectChemistry,
            args: DetectChemistryStageInputs {
                sample_def: vec![
                    SampleDefBuilder::new(
                        FastqMode::ILMN_BCL2FASTQ,
                        "../dui_tests/test_resources/cellranger-count/CR-300-01_SC3pv3_15k".into()
                    ).sample_names(vec!["CR-300-01"]).build(),
                    SampleDefBuilder::new(
                        FastqMode::ILMN_BCL2FASTQ,
                        "../dui_tests/test_resources/cellranger-count/CR-CUSTOM-22_SC5pr2_15k".into()
                    ).sample_names(vec!["CR-CUSTOM-22"]).build(),
                ],
                chemistry_specs: specs_input_gex(AutoChemistryName::Count),
                ..Default::default()
            },
        );
        Ok(())
    }

    #[test]
    fn test_cycle_failure() {
        let args = DetectChemistryStageInputs {
            sample_def: vec![
                SampleDefBuilder::new(
                    FastqMode::BCL_PROCESSOR,
                    "../dui_tests/test_resources/cellranger-count/cycle_failure_bc_v3".into(),
                )
                .sample_indices(vec!["SI-P2-D3"])
                .build(),
            ],
            chemistry_specs: specs_input_gex(AutoChemistryName::Count),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(
            outs.chemistry_defs.primary().name,
            ChemistryName::ThreePrimeV3PolyA
        );
        assert!(!outs.is_antibody_only);
    }

    #[test]
    fn test_vdj_pe() {
        let args = DetectChemistryStageInputs {
            sample_def: vec![
                SampleDefBuilder::new(
                    FastqMode::ILMN_BCL2FASTQ,
                    "../dui_tests/test_resources/cellranger-vdj/vdj_v1_hs_pbmc3_t".into(),
                )
                .sample_names(vec!["vdj_v1_hs_pbmc3_t"])
                .build(),
            ],
            chemistry_specs: specs_input_gex(AutoChemistryName::Vdj),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(outs.chemistry_defs.primary().name, ChemistryName::VdjPE);
        assert!(!outs.is_antibody_only);
    }

    #[test]
    fn test_vdj_r2() {
        let args = DetectChemistryStageInputs {
            sample_def: vec![
                SampleDefBuilder::new(
                    FastqMode::ILMN_BCL2FASTQ,
                    "../dui_tests/test_resources/cellranger-vdj/vdj_v1_mm_pbmc4_b_26x91".into(),
                )
                .sample_names(vec!["vdj_v1_mm_pbmc4_b_26x91"])
                .build(),
            ],
            chemistry_specs: specs_input_gex(AutoChemistryName::Vdj),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(outs.chemistry_defs.primary().name, ChemistryName::VdjR2);
        assert!(!outs.is_antibody_only);
    }

    #[test]
    fn test_ab_only_3pv3() {
        let args = DetectChemistryStageInputs {
            sample_def: vec![
                SampleDefBuilder::new(
                    FastqMode::ILMN_BCL2FASTQ,
                    "../dui_tests/test_resources/cellranger-count/pbmc_1k_protein_v3_antibody".into()
                )
                .sample_names(vec!["pbmc_1k_protein_v3_antibody"])
                .lanes(vec![4])
                .library_type(LibraryType::Antibody)
                .build(),
            ],
            chemistry_specs: specs_input_gex(AutoChemistryName::Count),
            feature_reference: Some("../dui_tests/test_resources/reference/feature/Biolegend_x22_TotalSeqB_14Ab_Isotype_Controls_20180601.csv".into()),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(
            outs.chemistry_defs.primary().name,
            ChemistryName::ThreePrimeV3CS1
        );
        assert!(outs.is_antibody_only);
    }

    #[test]
    fn test_ab_short_only_3pv3() {
        let args = DetectChemistryStageInputs {
            sample_def: vec![
                SampleDefBuilder::new(
                    FastqMode::ILMN_BCL2FASTQ,
                    "../dui_tests/test_resources/cellranger-count/pbmc_1k_protein_v3_antibody_len24".into()
                )
                .sample_names(vec!["pbmc_1k_protein_v3_antibody"])
                .lanes(vec![4])
                .library_type(LibraryType::Antibody)
                .build(),
                ],
            chemistry_specs: specs_input_gex(AutoChemistryName::Count),
            feature_reference: Some("../dui_tests/test_resources/reference/feature/Biolegend_x22_TotalSeqB_14Ab_Isotype_Controls_len24_20210721.csv".into()),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(
            outs.chemistry_defs.primary().name,
            ChemistryName::ThreePrimeV3CS1
        );
        assert!(outs.is_antibody_only);
    }

    #[test]
    fn test_ab_only_5p() {
        let args = DetectChemistryStageInputs {
            sample_def: vec![
                SampleDefBuilder::new(
                    FastqMode::ILMN_BCL2FASTQ,
                   "../dui_tests/test_resources/cellranger-count/vdj_v1_hs_pbmc2_antibody".into()
                )
                .sample_names(vec!["vdj_v1_hs_pbmc2_antibody"])
                .lanes(vec![1])
                .library_type(LibraryType::Antibody)
                .build(),
            ],
            chemistry_specs: specs_input_gex(AutoChemistryName::Count),
            feature_reference: Some("../dui_tests/test_resources/reference/feature/Biolegend_x22_TotalSeqB_14Ab_Isotype_Controls_20180601.csv".into()),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(
            outs.chemistry_defs.primary().name,
            ChemistryName::FeatureBarcodingOnly
        );
        assert!(outs.is_antibody_only);
    }

    #[test]
    fn test_ab_and_gex_3pv3() {
        let args = DetectChemistryStageInputs {
            sample_def: vec![SampleDefBuilder::new(
                FastqMode::ILMN_BCL2FASTQ,
                "../dui_tests/test_resources/cellranger-count/pbmc_1k_protein_v3_antibody".into(),
            ).sample_names(vec!["pbmc_1k_protein_v3_antibody"]).lanes(vec![4]).library_type(LibraryType::Antibody).build(),
            SampleDefBuilder::new(
                FastqMode::ILMN_BCL2FASTQ,
                "../dui_tests/test_resources/cellranger-count/pbmc_1k_protein_v3_gex".into(),
            ).sample_names(vec!["pbmc_1k_protein_v3_gex"]).lanes(vec![4]).library_type(LibraryType::Gex).build(),
            ],
            chemistry_specs: specs_input_gex(AutoChemistryName::Count),
            feature_reference: Some("../dui_tests/test_resources/reference/feature/Biolegend_x22_TotalSeqB_14Ab_Isotype_Controls_20180601.csv".into()),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(
            outs.chemistry_defs
                .values()
                .map(|def| def.name)
                .collect::<TxHashSet<_>>(),
            set![
                ChemistryName::ThreePrimeV3PolyA,
                ChemistryName::ThreePrimeV3CS1
            ]
        );
        assert!(!outs.is_antibody_only);
    }

    #[test]
    fn test_ab_short_and_gex_3pv3() {
        let args = DetectChemistryStageInputs {
            sample_def: vec![
                SampleDefBuilder::new(
                    FastqMode::ILMN_BCL2FASTQ,
                    "../dui_tests/test_resources/cellranger-count/pbmc_1k_protein_v3_antibody_len24".into(),
                )
                .sample_names(vec!["pbmc_1k_protein_v3_antibody"])
                .lanes(vec![4])
                .library_type(LibraryType::Antibody)
                .build(),
                SampleDefBuilder::new(
                    FastqMode::ILMN_BCL2FASTQ,
                    "../dui_tests/test_resources/cellranger-count/pbmc_1k_protein_v3_gex".into(),
                )
                .sample_names(vec!["pbmc_1k_protein_v3_gex"])
                .lanes(vec![4])
                .library_type(LibraryType::Gex)
                .build(),
            ],
            chemistry_specs: specs_input_gex(AutoChemistryName::Count),
            feature_reference: Some("../dui_tests/test_resources/reference/feature/Biolegend_x22_TotalSeqB_14Ab_Isotype_Controls_len24_20210721.csv".into()),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(
            outs.chemistry_defs
                .values()
                .map(|def| def.name)
                .collect::<TxHashSet<_>>(),
            set![
                ChemistryName::ThreePrimeV3PolyA,
                ChemistryName::ThreePrimeV3CS1
            ]
        );
        assert!(!outs.is_antibody_only);
    }

    #[test]
    fn test_ab_too_short_and_gex_3pv3() {
        let args = DetectChemistryStageInputs {
            sample_def: vec![
                SampleDefBuilder::new(
                    FastqMode::ILMN_BCL2FASTQ,
                    "../dui_tests/test_resources/cellranger-count/pbmc_1k_protein_v3_antibody"
                        .into(),
                )
                .sample_names(vec!["pbmc_1k_protein_v3_antibody"])
                .lanes(vec![4])
                .library_type(LibraryType::Antibody)
                .build(),
                SampleDefBuilder::new(
                    FastqMode::ILMN_BCL2FASTQ,
                    "../dui_tests/test_resources/cellranger-count/pbmc_1k_protein_v3_gex".into(),
                )
                .sample_names(vec!["pbmc_1k_protein_v3_gex"])
                .lanes(vec![4])
                .library_type(LibraryType::Gex)
                .build(),
            ],
            chemistry_specs: specs_input_gex(AutoChemistryName::Count),
            feature_reference: Some(
                "../dui_tests/test_resources/reference/feature/too_long_features.csv".into(),
            ),
            ..Default::default()
        };
        err_snapshot!(args);
    }

    #[test]
    fn test_crispr_and_gex_3pv3() {
        let args = DetectChemistryStageInputs {
            sample_def: vec![
                SampleDefBuilder::new(
                    FastqMode::ILMN_BCL2FASTQ,
                    "../dui_tests/test_resources/cellranger-count/K562_5k_crispr_v3_crispr".into(),
                )
                .sample_names(vec!["K562_5k_crispr_v3_crispr"])
                .library_type(LibraryType::Crispr)
                .build(),
                SampleDefBuilder::new(
                    FastqMode::ILMN_BCL2FASTQ,
                    "../dui_tests/test_resources/cellranger-count/K562_5k_crispr_v3_gex".into(),
                )
                .sample_names(vec!["K562_5k_crispr_v3_gex"])
                .library_type(LibraryType::Gex)
                .build(),
            ],
            chemistry_specs: specs_input_gex(AutoChemistryName::Count),
            feature_reference: Some(
                "../dui_tests/test_resources/reference/feature/RAB1A_NT_Aug2020.csv".into(),
            ),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(
            outs.chemistry_defs
                .values()
                .map(|def| def.name)
                .collect::<TxHashSet<_>>(),
            set![
                ChemistryName::ThreePrimeV3PolyA,
                ChemistryName::ThreePrimeV3CS1
            ]
        );
        assert!(!outs.is_antibody_only);
    }

    #[test]
    fn test_crispr_too_short_and_gex_3pv3() {
        let args = DetectChemistryStageInputs {
            sample_def: vec![
                SampleDefBuilder::new(
                    FastqMode::ILMN_BCL2FASTQ,
                    "../dui_tests/test_resources/cellranger-count/K562_5k_crispr_v3_crispr".into(),
                )
                .sample_names(vec!["K562_5k_crispr_v3_crispr"])
                .library_type(LibraryType::Crispr)
                .build(),
                SampleDefBuilder::new(
                    FastqMode::ILMN_BCL2FASTQ,
                    "../dui_tests/test_resources/cellranger-count/K562_5k_crispr_v3_gex".into(),
                )
                .sample_names(vec!["K562_5k_crispr_v3_gex"])
                .library_type(LibraryType::Gex)
                .build(),
            ],
            chemistry_specs: specs_input_gex(AutoChemistryName::Count),
            feature_reference: Some(
                "../dui_tests/test_resources/reference/feature/too_long_features.csv".into(),
            ),
            ..Default::default()
        };
        err_snapshot!(args);
    }

    #[test]
    fn test_sc3pv1() {
        let args = DetectChemistryStageInputs {
            sample_def: vec![
                SampleDefBuilder::new(
                    FastqMode::ILMN_BCL2FASTQ,
                    "../dui_tests/test_resources/cellranger-count/CR-CUSTOM-02_SC3pv1_15k".into(),
                )
                .sample_names(vec!["test_sample"])
                .build(),
            ],
            chemistry_specs: specs_input_gex(AutoChemistryName::Count),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(
            outs.chemistry_defs.primary().name,
            ChemistryName::ThreePrimeV1
        );
        assert!(!outs.is_antibody_only);
    }

    #[test]
    fn test_mixed_chemistries_3pv1_3pv2() {
        let args = DetectChemistryStageInputs {
            sample_def: vec![
                SampleDefBuilder::new(
                    FastqMode::ILMN_BCL2FASTQ,
                    "../dui_tests/test_resources/cellranger-count/CR-CUSTOM-02_SC3pv1_15k".into(),
                )
                .sample_names(vec!["test_sample"])
                .build(),
                SampleDefBuilder::new(
                    FastqMode::ILMN_BCL2FASTQ,
                    "../dui_tests/test_resources/cellranger-count/CR-120-01_SC3pv2_15k".into(),
                )
                .sample_names(vec!["CR-120-01"])
                .build(),
            ],
            chemistry_specs: specs_input_gex(AutoChemistryName::Count),
            ..Default::default()
        };
        err_snapshot!(args);
    }

    #[test]
    fn test_no_whitelist_match_auto() -> Result<()> {
        stage_fail_dui_test!(
            NoWhitelistMatch,
            description: "The input data does not match any of the whitelists.",
            stage: DetectChemistry,
            args: DetectChemistryStageInputs {
                sample_def: vec![SampleDefBuilder::new(
                    FastqMode::ILMN_BCL2FASTQ,
                    "../dui_tests/test_resources/cellranger-count/no_barcode_div_v3".into(),
                ).sample_names(vec!["test_sample"]).build()],
                chemistry_specs: specs_input_gex(AutoChemistryName::Count),
                ..Default::default()
            },
        );
        Ok(())
    }

    #[test]
    fn test_reads_too_short() -> Result<()> {
        stage_fail_dui_test!(
            ReadsTooShort,
            description: "In the input data, Read 1 is only 24 bases & Read 2 is only 30 bases.",
            stage: DetectChemistry,
            args: DetectChemistryStageInputs {
                sample_def: vec![SampleDefBuilder::new(
                    FastqMode::ILMN_BCL2FASTQ,
                    "../dui_tests/test_resources/cellranger-count/reads_too_short_v2".into()
                ).sample_names(vec!["CR-120-01"]).build()],
                chemistry_specs: specs_input_gex(AutoChemistryName::Count),
                ..Default::default()
            },
        );
        Ok(())
    }

    #[test]
    fn test_vdj_pe_with_r1_r2_length() {
        let args = DetectChemistryStageInputs {
            sample_def: vec![
                SampleDefBuilder::new(
                    FastqMode::ILMN_BCL2FASTQ,
                    "../dui_tests/test_resources/cellranger-vdj/vdj_v1_hs_pbmc3_t".into(),
                )
                .sample_names(vec!["vdj_v1_hs_pbmc3_t"])
                .build(),
            ],
            chemistry_specs: specs_input_gex(AutoChemistryName::Vdj),
            r1_length: Some(26),
            r2_length: Some(91),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(outs.chemistry_defs.primary().name, ChemistryName::VdjR2);
        assert!(!outs.is_antibody_only);
    }

    #[test]
    fn test_r1_too_short_3pv3() {
        let args = DetectChemistryStageInputs {
            sample_def: vec![
                SampleDefBuilder::new(
                    FastqMode::ILMN_BCL2FASTQ,
                    "../dui_tests/test_resources/cellranger-count/CR-300-01_SC3pv3_15k".into(),
                )
                .sample_names(vec!["CR-300-01"])
                .build(),
            ],
            chemistry_specs: specs_input_gex(AutoChemistryName::Count),
            r1_length: Some(25),
            ..Default::default()
        };
        err_snapshot!(args);
    }

    #[test]
    fn test_r2_too_short_3pv3() -> Result<()> {
        stage_fail_dui_test!(
            R2LengthTooShort,
            description: "User inputs an r2-length that is too small.",
            stage: DetectChemistry,
            args: DetectChemistryStageInputs {
                sample_def: vec![SampleDefBuilder::new(
                    FastqMode::ILMN_BCL2FASTQ,
                    "../dui_tests/test_resources/cellranger-count/CR-300-01_SC3pv3_15k".into(),
                )
                .sample_names(vec!["CR-300-01"])
                .build()],
                chemistry_specs: specs_input_gex(AutoChemistryName::Count),
                r2_length: Some(14),
                ..Default::default()
            },
        );
        Ok(())
    }

    #[test]
    fn test_sample_def_too_much_trimming() -> Result<()> {
        stage_fail_dui_test!(
            NotEnoughReadLength,
            description: "Sample def contains an r1-length that is too small.",
            stage: DetectChemistry,
            args: DetectChemistryStageInputs {
                sample_def: vec![SampleDefBuilder::new(
                    FastqMode::ILMN_BCL2FASTQ,
                    "../dui_tests/test_resources/cellranger-count/CR-300-01_SC3pv3_15k".into(),
                ).sample_names(vec!["CR-300-01"]).r1_length(25).build()],
                chemistry_specs: specs_input_gex(AutoChemistryName::Count),
                ..Default::default()
            },
        );
        Ok(())
    }

    #[test]
    fn test_r1_r2_okay_3pv3() {
        let args = DetectChemistryStageInputs {
            sample_def: vec![
                SampleDefBuilder::new(
                    FastqMode::ILMN_BCL2FASTQ,
                    "../dui_tests/test_resources/cellranger-count/CR-300-01_SC3pv3_15k".into(),
                )
                .sample_names(vec!["CR-300-01"])
                .build(),
            ],
            chemistry_specs: specs_input_gex(AutoChemistryName::Count),
            r1_length: Some(26),
            r2_length: Some(75),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(
            outs.chemistry_defs.primary().name,
            ChemistryName::ThreePrimeV3PolyA
        );
        assert!(!outs.is_antibody_only);
    }

    #[test]
    fn test_few_reads_explicit_chem() {
        let args = DetectChemistryStageInputs {
            sample_def: vec![
                SampleDefBuilder::new(
                    FastqMode::ILMN_BCL2FASTQ,
                    "../dui_tests/test_resources/cellranger-count/sc3pv2_5k".into(),
                )
                .sample_names(vec!["tinygex"])
                .build(),
            ],
            chemistry_specs: specs_input_gex(ChemistryName::ThreePrimeV2),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(
            outs.chemistry_defs.primary().name,
            ChemistryName::ThreePrimeV2
        );
        assert!(!outs.is_antibody_only);
    }

    #[test]
    // When chemistry is not specified, check for probe bc mixtures is on
    // if the input is from count or config is missing [samples] section.
    fn test_rtl_singleplex_probe_bc_mixture_auto_chem() {
        let args = DetectChemistryStageInputs {
            sample_def: vec![SampleDefBuilder::new(
                FastqMode::ILMN_BCL2FASTQ,
                "../dui_tests/test_resources/cellranger-count/RTL_single_sample_probe_bcs_mixture"
                    .into(),
            ).sample_names(vec!["test_sample"]).build()],
            chemistry_specs: specs_input_gex(AutoChemistryName::Count),
            ..Default::default()
        };
        err_snapshot!(args);
    }

    #[test]
    // When refined SFRP chemistry is specified, check for probe bc mixtures is skipped
    fn test_rtl_singleplex_probe_bc_mixture_explicit_chem() {
        let args = DetectChemistryStageInputs {
            sample_def: vec![SampleDefBuilder::new(
                FastqMode::ILMN_BCL2FASTQ,
                "../dui_tests/test_resources/cellranger-count/RTL_single_sample_probe_bcs_mixture"
                    .into(),
            )
            .sample_names(vec!["test_sample"])
            .library_type(LibraryType::Gex)
            .build()],
            chemistry_specs: specs_input_gex(ChemistryName::SFRP),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(outs.chemistry_defs.primary().name, ChemistryName::SFRP);
    }

    #[test]
    fn test_rtl_no_probe_set() {
        let multi_config = b"\
            [gene-expression]\n\
            ref,/ref\n\
            create-bam,false\n\
            [libraries]\n\
            fastq_id,fastqs,feature_types\n\
            rtl,/fastqs,Gene Expression\n";
        let multi_config_csv = tempfile::Builder::new().suffix(".csv").tempfile().unwrap();
        multi_config_csv.as_file().write_all(multi_config).unwrap();

        let args = DetectChemistryStageInputs {
            sample_def: vec![
                SampleDefBuilder::new(
                    FastqMode::ILMN_BCL2FASTQ,
                    "../dui_tests/test_resources/cellranger-count/CR-300-01_SC3pv3_15k".into(),
                )
                .sample_names(vec!["CR-300-01"])
                .build(),
            ],
            chemistry_specs: specs_input_gex(ChemistryName::SFRP),
            multi_config: Some(MultiConfigCsvFile::from(multi_config_csv.path())),
            ..Default::default()
        };
        err_snapshot!(args);
    }

    #[test]
    fn test_overhang_chemistry() {
        let args = DetectChemistryStageInputs {
            sample_def: vec![
                SampleDefBuilder::new(
                    FastqMode::ILMN_BCL2FASTQ,
                    "../dui_tests/test_resources/cellranger-count/CR-300-01_SC3pv3_15k".into(),
                )
                .sample_names(vec!["CR-300-01"])
                .build(),
            ],
            chemistry_specs: specs_input_gex(AutoChemistryName::Count),
            multi_config: Some(
                "../dui_tests/test_resources/cellranger-multi/overhang_3pv3.csv".into(),
            ),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(
            outs.chemistry_defs.primary().name,
            ChemistryName::ThreePrimeV3PolyAOCM
        );
    }

    #[test]
    fn test_auto_chemistry_spec_sc5pv3_sbx() {
        let args = DetectChemistryStageInputs {
            sample_def: vec![
                SampleDefBuilder::new(FastqMode::SBX, "test/fastq/sbx/sbx_fiveprime_frag".into())
                    .build(),
            ],
            chemistry_specs: specs_input_gex(AutoChemistryName::Count),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        let chem = outs.chemistry_defs.primary();
        assert_eq!(chem.name, ChemistryName::FivePrimeR2V3);
        assert_eq!(
            outs.library_trimmers[&LibraryType::Gex],
            trimmer_for_chemistry(
                chem,
                LibraryType::Gex,
                layouts::FIVE_PRIME_GEX,
                &FastqMode::SBX.adapter_trimming_parameters(None).unwrap()
            )
            .unwrap()
            .unwrap()
        );
        assert!(!outs.is_antibody_only);
    }

    #[test]
    fn test_auto_chemistry_spec_sc3pv4_sbx() {
        let args = DetectChemistryStageInputs {
            sample_def: vec![
                SampleDefBuilder::new(FastqMode::SBX, "test/fastq/sbx/sbx_threeprime_frag".into())
                    .build(),
            ],
            chemistry_specs: specs_input_gex(AutoChemistryName::Count),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        let chem = outs.chemistry_defs.primary();
        assert_eq!(chem.name, ChemistryName::ThreePrimeV4PolyA);
        assert_eq!(
            outs.library_trimmers[&LibraryType::Gex],
            trimmer_for_chemistry(
                chem,
                LibraryType::Gex,
                layouts::FIVE_PRIME_GEX,
                &FastqMode::SBX.adapter_trimming_parameters(None).unwrap()
            )
            .unwrap()
            .unwrap()
        );
        assert!(!outs.is_antibody_only);
    }

    #[test]
    fn test_auto_chemistry_spec_flex_sbx() {
        let multi_config = b"\
            [gene-expression]\n\
            ref,/ref\n\
            probe-set,foobar.csv\n\
            create-bam,false\n\
            [libraries]\n\
            fastq_id,fastqs,feature_types\n\
            rtl,/fastqs,Gene Expression\n\
            [samples]\n\
            sample_id,probe_barcode_ids
            test,A-A01
            ";
        let multi_config_csv = tempfile::Builder::new().suffix(".csv").tempfile().unwrap();
        multi_config_csv.as_file().write_all(multi_config).unwrap();
        let args = DetectChemistryStageInputs {
            sample_def: vec![
                SampleDefBuilder::new(FastqMode::SBX, "test/fastq/sbx/sbx_flex_gex".into()).build(),
            ],
            chemistry_specs: specs_input_gex(AutoChemistryName::Count),
            multi_config: Some(MultiConfigCsvFile::from(multi_config_csv.path())),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        let chem = outs.chemistry_defs.primary();
        assert_eq!(chem.name, ChemistryName::Flex_v2_R1);
        assert_eq!(
            outs.library_trimmers[&LibraryType::Gex],
            trimmer_for_chemistry(
                chem,
                LibraryType::Gex,
                layouts::FLEX_GEX,
                &FastqMode::SBX.adapter_trimming_parameters(None).unwrap()
            )
            .unwrap()
            .unwrap()
        );
        assert!(!outs.is_antibody_only);
    }

    #[cfg(feature = "slow_tests")]
    mod tests_requiring_alignment {
        use super::*;
        use cr_types::sample_def::FastqMode;
        use dui_tests::stage_test::StageFailTest;
        use dui_tests::{DuiTest, stage_fail_dui_test};
        use test_refdata::{refdata_path, sere_path, testdata_path};

        #[test]
        fn test_auto_chemistry_spec_3pv2() {
            let args = DetectChemistryStageInputs {
                sample_def: vec![
                    SampleDefBuilder::new(
                        FastqMode::ILMN_BCL2FASTQ,
                        "../dui_tests/test_resources/cellranger-count/CR-120-01_SC3pv2_15k".into(),
                    )
                    .sample_names(vec!["CR-120-01"])
                    .build(),
                ],
                reference_path: Some(refdata_path("GRCh38")),
                chemistry_specs: specs_input_gex(AutoChemistryName::Count),
                ..Default::default()
            };
            let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
            assert_eq!(
                outs.chemistry_defs.primary().name,
                ChemistryName::ThreePrimeV2
            );
            assert!(!outs.is_antibody_only);
        }

        #[test]
        fn test_auto_chemistry_spec_5p_r2() {
            let args = DetectChemistryStageInputs {
                sample_def: vec![
                    SampleDefBuilder::new(
                        FastqMode::ILMN_BCL2FASTQ,
                        "../dui_tests/test_resources/cellranger-count/CR-CUSTOM-22_SC5pr2_15k"
                            .into(),
                    )
                    .sample_names(vec!["CR-CUSTOM-22"])
                    .build(),
                ],
                reference_path: Some(refdata_path("GRCh38")),
                chemistry_specs: specs_input_gex(AutoChemistryName::Count),
                ..Default::default()
            };
            let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
            assert_eq!(
                outs.chemistry_defs.primary().name,
                ChemistryName::FivePrimeR2
            );
            assert!(!outs.is_antibody_only);
        }

        #[test]
        fn test_5p_auto_chemistry_spec_5p_r2() {
            let args = DetectChemistryStageInputs {
                sample_def: vec![
                    SampleDefBuilder::new(
                        FastqMode::ILMN_BCL2FASTQ,
                        "../dui_tests/test_resources/cellranger-count/CR-CUSTOM-22_SC5pr2_15k"
                            .into(),
                    )
                    .sample_names(vec!["CR-CUSTOM-22"])
                    .build(),
                ],
                reference_path: Some(refdata_path("GRCh38")),
                chemistry_specs: specs_input_gex(AutoChemistryName::FivePrime),
                ..Default::default()
            };
            let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
            assert_eq!(
                outs.chemistry_defs.primary().name,
                ChemistryName::FivePrimeR2
            );
            assert!(!outs.is_antibody_only);
        }

        #[test]
        fn test_auto_chemistry_spec_5p_pe() {
            let args = DetectChemistryStageInputs {
                sample_def: vec![
                    SampleDefBuilder::new(
                        FastqMode::ILMN_BCL2FASTQ,
                        "../dui_tests/test_resources/cellranger-count/CR-CUSTOM-17_SC5p-PE_15k"
                            .into(),
                    )
                    .sample_names(vec!["CR-CUSTOM-17"])
                    .build(),
                ],
                reference_path: Some(refdata_path("hg19_and_mm10-3.0.0")),
                chemistry_specs: specs_input_gex(AutoChemistryName::Count),
                ..Default::default()
            };
            let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
            assert_eq!(
                outs.chemistry_defs.primary().name,
                ChemistryName::FivePrimePE
            );
            assert!(!outs.is_antibody_only);
        }

        #[test]
        fn test_unsupported_auto_chemistry_spec_5p_r2_and_pe() {
            let args = DetectChemistryStageInputs {
                sample_def: vec![
                    SampleDefBuilder::new(
                        FastqMode::ILMN_BCL2FASTQ,
                        "../dui_tests/test_resources/cellranger-count/CR-CUSTOM-22_SC5pr2_15k"
                            .into(),
                    )
                    .sample_names(vec!["CR-CUSTOM-22"])
                    .build(),
                    SampleDefBuilder::new(
                        FastqMode::ILMN_BCL2FASTQ,
                        "../dui_tests/test_resources/cellranger-count/CR-CUSTOM-17_SC5p-PE_15k"
                            .into(),
                    )
                    .sample_names(vec!["CR-CUSTOM-17"])
                    .build(),
                ],
                reference_path: Some(refdata_path("hg19_and_mm10-3.0.0")),
                chemistry_specs: specs_input_gex(AutoChemistryName::Count),
                ..Default::default()
            };
            let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
            assert_eq!(
                outs.chemistry_defs.primary().name,
                ChemistryName::FivePrimeR2
            );
        }

        #[test]
        fn test_unsupported_auto_chemistry_spec_5p_pe_and_r2() {
            let args = DetectChemistryStageInputs {
                sample_def: vec![
                    SampleDefBuilder::new(
                        FastqMode::ILMN_BCL2FASTQ,
                        "../dui_tests/test_resources/cellranger-count/CR-CUSTOM-17_SC5p-PE_15k"
                            .into(),
                    )
                    .sample_names(vec!["CR-CUSTOM-17"])
                    .build(),
                    SampleDefBuilder::new(
                        FastqMode::ILMN_BCL2FASTQ,
                        "../dui_tests/test_resources/cellranger-count/CR-CUSTOM-22_SC5pr2_15k"
                            .into(),
                    )
                    .sample_names(vec!["CR-CUSTOM-22"])
                    .build(),
                ],
                reference_path: Some(refdata_path("hg19_and_mm10-3.0.0")),
                chemistry_specs: specs_input_gex(AutoChemistryName::Count),
                ..Default::default()
            };
            let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
            assert_eq!(
                outs.chemistry_defs.primary().name,
                ChemistryName::FivePrimeR2
            );
        }

        #[test]
        fn test_ab_and_gex_5p() {
            let args = DetectChemistryStageInputs {
                sample_def: vec![
                    SampleDefBuilder::new(
                        FastqMode::ILMN_BCL2FASTQ,
                        sere_path("vdj_nextgem_hs_pbmc3_5gex_protein/vdj_nextgem_hs_pbmc3_5gex_protein_antibody"),
                    )
                    .sample_names(vec!["vdj_nextgem_hs_pbmc3_5gex_protein_antibody"])
                    .lanes(vec![1])
                    .library_type(LibraryType::Antibody)
                    .build(),
                    SampleDefBuilder::new(
                        FastqMode::ILMN_BCL2FASTQ,
                        sere_path("vdj_nextgem_hs_pbmc3_5gex_protein/vdj_nextgem_hs_pbmc3_5gex_protein_gex"),
                    )
                    .sample_names(vec!["vdj_nextgem_hs_pbmc3_5gex_protein_gex"])
                    .lanes(vec![1])
                    .library_type(LibraryType::Gex)
                    .build(),
                ],
                chemistry_specs: specs_input_gex(AutoChemistryName::Count),
                reference_path: Some(refdata_path("GRCh38")),
                feature_reference: Some("../dui_tests/test_resources/reference/feature/Biolegend_x22_TotalSeqB_14Ab_Isotype_Controls_20180601.csv".into()),
                ..Default::default()
            };
            let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
            assert_eq!(
                outs.chemistry_defs.primary().name,
                ChemistryName::FivePrimeR2
            );
            assert!(!outs.is_antibody_only);
        }

        #[test]
        fn test_gex_and_ab_sfrp() {
            let args = DetectChemistryStageInputs {
                sample_def: vec![
                    SampleDefBuilder::new(
                        FastqMode::ILMN_BCL2FASTQ,
                        testdata_path("fastqs/cellranger/multi/1247954_1270332_frp_gex_and_ab_fastqs"),
                    )
                    .sample_names(vec!["1247954_GEX"])
                    .lanes(vec![4])
                    .library_type(LibraryType::Gex)
                    .build(),
                    SampleDefBuilder::new(
                        FastqMode::ILMN_BCL2FASTQ,
                        testdata_path("fastqs/cellranger/multi/1247954_1270332_frp_gex_and_ab_fastqs"),
                    )
                    .sample_names(vec!["1247954_AB"])
                    .lanes(vec![4])
                    .library_type(LibraryType::Antibody)
                    .build()
                ],
                multi_config: Some("../dui_tests/test_resources/cellranger-multi/1247954_1270332_rtl_sfrp_gex_and_ab.csv".into()),
                chemistry_specs: specs_input_gex(AutoChemistryName::Count),
                feature_reference: Some(testdata_path("cellranger/multi/feature_refs/BioLegend_x22_TotalSeqB_Human_PBMC_20210114.csv").into()),
                ..Default::default()
            };
            let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
            assert_eq!(outs.chemistry_defs.primary().name, ChemistryName::SFRP);
            assert!(!outs.is_antibody_only);
        }

        #[test]
        fn test_gex_only_sfrp() {
            let args = DetectChemistryStageInputs {
                sample_def: vec![
                    SampleDefBuilder::new(
                        FastqMode::ILMN_BCL2FASTQ,
                        testdata_path("fastqs/cellranger/multi/1245033_sfrp_fastqs"),
                    )
                    .sample_names(vec!["CTACTGAATT"])
                    .lanes(vec![4])
                    .library_type(LibraryType::Gex)
                    .build(),
                ],
                multi_config: Some(
                    "../dui_tests/test_resources/cellranger-multi/1245033_sfrp.csv".into(),
                ),
                chemistry_specs: specs_input_gex(AutoChemistryName::Count),
                ..Default::default()
            };
            let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
            assert_eq!(outs.chemistry_defs.primary().name, ChemistryName::SFRP);
        }

        #[test]
        fn test_gex_only_mfrp() {
            let args = DetectChemistryStageInputs {
                sample_def: vec![
                    SampleDefBuilder::new(
                        FastqMode::ILMN_BCL2FASTQ,
                        testdata_path("fastqs/cellranger/multi/1216790_1245489_mfrp_fastqs"),
                    )
                    .sample_names(vec!["AATTTCGGGT"])
                    .lanes(vec![4])
                    .library_type(LibraryType::Gex)
                    .build(),
                ],
                multi_config: Some(
                    "../dui_tests/test_resources/cellranger-multi/1216790_1245489.csv".into(),
                ),
                chemistry_specs: specs_input_gex(AutoChemistryName::Count),
                ..Default::default()
            };
            let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
            assert_eq!(outs.chemistry_defs.primary().name, ChemistryName::MFRP_RNA);
        }

        #[test]
        fn test_mfrp_r1() {
            let args = DetectChemistryStageInputs {
                sample_def: vec![
                    SampleDefBuilder::new(
                        FastqMode::ILMN_BCL2FASTQ,
                        testdata_path("fastqs/cellranger/multi/mfrp_r1_50_90"),
                    )
                    .sample_names(vec!["tiny_gex"])
                    .lanes(vec![1])
                    .library_type(LibraryType::Gex)
                    .build(),
                ],
                multi_config: Some(
                    "../dui_tests/test_resources/cellranger-multi/1320608_mfrp_r1_50_90.csv".into(),
                ),
                chemistry_specs: specs_input_gex(AutoChemistryName::Count),
                ..Default::default()
            };
            let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
            assert_eq!(
                outs.chemistry_defs.primary().name,
                ChemistryName::MFRP_RNA_R1
            );
        }

        #[test]
        fn test_ab_only_mfrp() {
            let args = DetectChemistryStageInputs {
                sample_def: vec![
                    SampleDefBuilder::new(
                        FastqMode::ILMN_BCL2FASTQ,
                        testdata_path("fastqs/cellranger/multi/1461767_mfrp_gex_ab/ab"),
                    )
                    .sample_names(vec!["tiny_ab"])
                    .lanes(vec![1])
                    .library_type(LibraryType::Gex)
                    .build(),
                ],
                multi_config: Some(
                    "../dui_tests/test_resources/cellranger-multi/1461767_mfrp_ab_only.csv".into(),
                ),
                chemistry_specs: specs_input_gex(AutoChemistryName::Count),
                ..Default::default()
            };
            let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
            assert_eq!(outs.chemistry_defs.primary().name, ChemistryName::MFRP_Ab);
        }

        #[test]
        fn test_ab_only_mfrp_r1() {
            let args = DetectChemistryStageInputs {
                sample_def: vec![
                    SampleDefBuilder::new(
                        FastqMode::ILMN_BCL2FASTQ,
                        testdata_path("fastqs/cellranger/multi/1464930_mfrp_r1_gex_ab/ab"),
                    )
                    .sample_names(vec!["tiny_ab"])
                    .lanes(vec![1])
                    .library_type(LibraryType::Gex)
                    .build(),
                ],
                multi_config: Some(
                    "../dui_tests/test_resources/cellranger-multi/1464930_mfrp_r1_ab_only.csv"
                        .into(),
                ),
                chemistry_specs: specs_input_gex(AutoChemistryName::Count),
                ..Default::default()
            };
            let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
            assert_eq!(
                outs.chemistry_defs.primary().name,
                ChemistryName::MFRP_Ab_R1
            );
        }

        #[test]
        fn test_per_lib_chem_mfrp() {
            let args = DetectChemistryStageInputs{
                sample_def: vec![
                    SampleDefBuilder::new(
                        FastqMode::ILMN_BCL2FASTQ,
                        testdata_path("fastqs/cellranger/multi/1548572/tiny_ab"),
                    )
                    .r1_length(28)
                    .sample_names(vec!["tiny_ab"])
                    .lanes(vec![1])
                    .library_type(LibraryType::Antibody)
                    .build(),
                    SampleDefBuilder::new(
                        FastqMode::ILMN_BCL2FASTQ,
                        testdata_path("fastqs/cellranger/multi/1548572/tiny_gex"),
                    )
                    .sample_names(vec!["tiny_gex"])
                    .lanes(vec![1])
                    .library_type(LibraryType::Gex)
                    .build(),
                ],
                multi_config: Some("../dui_tests/test_resources/cellranger-multi/1548572_mfrp_r1_ab_r2pos50.csv".into()),
                feature_reference: Some("../dui_tests/test_resources/reference/feature/PTG_Multipro_Human_Immune_Intracellular_Antibodies_PROT002.csv".into()),
                chemistry_specs: specs_input_gex(AutoChemistryName::Count),
                ..Default::default()
            };
            let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
            assert_eq!(
                outs.chemistry_defs[&LibraryType::Gex].name,
                ChemistryName::MFRP_RNA_R1
            );
            assert_eq!(
                outs.chemistry_defs[&LibraryType::Antibody].name,
                ChemistryName::MFRP_Ab_R2pos50
            );
        }

        #[test]
        fn test_ab_only_sfrp() {
            let args = DetectChemistryStageInputs {
                sample_def: vec![
                    SampleDefBuilder::new(
                        FastqMode::ILMN_BCL2FASTQ,
                        testdata_path("fastqs/cellranger/multi/1247954_1270332_frp_ab_fastqs/ab"),
                    )
                    .sample_names(vec!["ab"])
                    .lanes(vec![1])
                    .library_type(LibraryType::Gex)
                    .build(),
                ],
                multi_config: Some(
                    "../dui_tests/test_resources/cellranger-multi/1247954_sfrp_ab_only.csv".into(),
                ),
                chemistry_specs: specs_input_gex(AutoChemistryName::Count),
                ..Default::default()
            };
            let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
            assert_eq!(outs.chemistry_defs.primary().name, ChemistryName::SFRP);
        }

        #[test]
        fn test_few_aligned_reads_3pv2() -> Result<()> {
            stage_fail_dui_test!(
                FewAlignedReads3pv2,
                description: "The input data is 3pv2, but few reads map confidently so that it \
                cannot be distinguished from 5'.",
                stage: DetectChemistry,
                args: DetectChemistryStageInputs {
                    sample_def: vec![
                        SampleDefBuilder::new(
                            FastqMode::ILMN_BCL2FASTQ,
                            "../dui_tests/test_resources/cellranger-count/few_aligned_reads_v2/".into(),
                        )
                        .sample_names(vec!["test_sample"])
                        .build(),
                    ],
                    reference_path: Some(refdata_path("GRCh38")),
                    chemistry_specs: specs_input_gex(AutoChemistryName::Count),
                    ..Default::default()
                },
            );
            Ok(())
        }

        #[test]
        fn test_mixed_chemistries_3pv2_5p_r2() {
            let args = DetectChemistryStageInputs {
                sample_def: vec![
                    SampleDefBuilder::new(
                        FastqMode::ILMN_BCL2FASTQ,
                        "../dui_tests/test_resources/cellranger-count/CR-120-01_SC3pv2_15k".into(),
                    )
                    .sample_names(vec!["CR-120-01"])
                    .build(),
                    SampleDefBuilder::new(
                        FastqMode::ILMN_BCL2FASTQ,
                        "../dui_tests/test_resources/cellranger-count/CR-CUSTOM-22_SC5pr2_15k"
                            .into(),
                    )
                    .sample_names(vec!["CR-CUSTOM-22"])
                    .build(),
                ],
                chemistry_specs: specs_input_gex(AutoChemistryName::Count),
                reference_path: Some(refdata_path("GRCh38")),
                ..Default::default()
            };
            err_snapshot!(args);
        }

        #[test]
        fn test_mixed_chemistries_3pv2_5p_pe() {
            let args = DetectChemistryStageInputs {
                sample_def: vec![
                    SampleDefBuilder::new(
                        FastqMode::ILMN_BCL2FASTQ,
                        "../dui_tests/test_resources/cellranger-count/CR-120-01_SC3pv2_15k".into(),
                    )
                    .sample_names(vec!["CR-120-01"])
                    .build(),
                    SampleDefBuilder::new(
                        FastqMode::ILMN_BCL2FASTQ,
                        "../dui_tests/test_resources/cellranger-count/CR-CUSTOM-17_SC5p-PE_15k"
                            .into(),
                    )
                    .sample_names(vec!["CR-CUSTOM-17"])
                    .build(),
                ],
                chemistry_specs: specs_input_gex(AutoChemistryName::Count),
                reference_path: Some(refdata_path("GRCh38")),
                ..Default::default()
            };
            err_snapshot!(args);
        }

        #[test]
        fn test_5pr2v3() {
            let args = DetectChemistryStageInputs {
                sample_def: vec![
                    SampleDefBuilder::new(
                        FastqMode::ILMN_BCL2FASTQ,
                        testdata_path("fastqs/cellranger/multi/1523899_gex_5pr2v3"),
                    )
                    .sample_names(vec!["tiny_gex"])
                    .lanes(vec![4])
                    .library_type(LibraryType::Gex)
                    .build(),
                ],
                multi_config: Some(
                    "../dui_tests/test_resources/cellranger-multi/1523899_gex_5pr2v3.csv".into(),
                ),
                chemistry_specs: specs_input_gex(AutoChemistryName::Count),
                ..Default::default()
            };
            let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
            assert_eq!(
                outs.chemistry_defs.primary().name,
                ChemistryName::FivePrimeR2V3
            );
        }

        #[test]
        fn test_5pr2ohv3() {
            let args = DetectChemistryStageInputs {
                sample_def: vec![
                    SampleDefBuilder::new(
                        FastqMode::ILMN_BCL2FASTQ,
                        testdata_path("fastqs/cellranger/multi/1523899_gex_5pr2v3"),
                    )
                    .sample_names(vec!["tiny_gex"])
                    .lanes(vec![4])
                    .library_type(LibraryType::Gex)
                    .build(),
                ],
                multi_config: Some(
                    "../dui_tests/test_resources/cellranger-multi/1523899_gex_5pr2ohv3.csv".into(),
                ),
                chemistry_specs: specs_input_gex(AutoChemistryName::Count),
                ..Default::default()
            };
            let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
            assert_eq!(
                outs.chemistry_defs.primary().name,
                ChemistryName::FivePrimeR2OCMV3
            );
        }

        #[test]
        fn test_3pv4() {
            let args = DetectChemistryStageInputs {
                sample_def: vec![
                    SampleDefBuilder::new(
                        FastqMode::ILMN_BCL2FASTQ,
                        testdata_path("fastqs/cellranger/multi/1561255_gex/tiny_gex"),
                    )
                    .sample_names(vec!["gex"])
                    .lanes(vec![1])
                    .library_type(LibraryType::Gex)
                    .build(),
                ],
                chemistry_specs: specs_input_gex(AutoChemistryName::Count),
                ..Default::default()
            };
            let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
            assert_eq!(
                outs.chemistry_defs.primary().name,
                ChemistryName::ThreePrimeV4PolyA
            );
        }
    }
}
