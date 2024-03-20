//! Martian stage DETECT_CHEMISTRY

use crate::barcode_overlap::FRPGemBarcodeOverlapRow;
use crate::detect_chemistry::chemistry_filter::{
    detect_chemistry_units, ChemistryFilter, DetectChemistryUnit,
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
use anyhow::{bail, ensure, Context, Result};
use barcode::whitelist::BarcodeId;
use cr_types::chemistry::{
    normalize_chemistry_def, AutoChemistryName, AutoOrRefinedChemistry, ChemistryDef,
    ChemistryDefs, ChemistryName, ChemistrySpecs,
};
use cr_types::reference::feature_reference::{FeatureConfig, FeatureReferenceFile};
use cr_types::sample_def::SampleDef;
use cr_types::LibraryType;
use fastq_set::read_pair::ReadPair;
use itertools::Itertools;
use martian::prelude::*;
use martian_derive::{make_mro, MartianStruct};
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::tabular_file::CsvFile;
use martian_filetypes::FileTypeWrite;
use metric::{join_metric_name, set, TxHashMap, TxHashSet};
use multi::config::{
    multiconst, MultiConfigCsv, MultiConfigCsvFile, ProbeBarcodeIterationMode, SamplesCsv,
};
use parameters_toml::{fiveprime_multiplexing, threeprime_lt_multiplexing};
use serde::{Deserialize, Serialize};
use slice_group_by::GroupBy;
use std::collections::HashMap;
use std::path::{Path, PathBuf};
#[allow(clippy::enum_glob_use)]
use ChemistryName::*;

const MIN_READS_NEEDED: usize = 10_000;

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
    pub custom_chemistry_def: Option<ChemistryDef>,
    pub feature_config: Option<FeatureConfig>,
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
}

pub type DetectedProbeBarcodePairingFile = JsonFile<TxHashMap<BarcodeId, BarcodeId>>;

#[derive(Clone, Serialize, Deserialize, MartianStruct)]
#[cfg_attr(test, derive(Debug))]
pub struct DetectChemistryStageOutputs {
    pub chemistry_defs: ChemistryDefs,
    pub is_antibody_only: bool,
    #[mro_retain]
    pub probe_barcode_overlap: Option<CsvFile<FRPGemBarcodeOverlapRow>>,
    #[mro_retain]
    pub detected_probe_barcode_pairing: Option<DetectedProbeBarcodePairingFile>,
}

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

        let units = &load_reads(&args)?;
        ensure!(!units.is_empty(), "no reads loaded");

        let chemistry_defs = select_chemistries(&args, multi_config_csv.as_ref(), units)?;

        let is_antibody_only = args
            .sample_def
            .iter()
            .all(|x| x.library_type == Some(LibraryType::Antibody));

        let mut outputs = DetectChemistryStageOutputs {
            chemistry_defs,
            is_antibody_only,
            probe_barcode_overlap: None,
            detected_probe_barcode_pairing: None,
        };

        if let Some(multi_config_csv) = &multi_config_csv {
            handle_probe_barcode_translation(&rover, multi_config_csv, units, &mut outputs)?;
        }

        Ok(outputs)
    }
}

/// Top-level entry point to define the map of chemistries per library type.
fn select_chemistries(
    args: &DetectChemistryStageInputs,
    multi_config_csv: Option<&MultiConfigCsv>,
    units: &[(DetectChemistryUnit, Vec<ReadPair>)],
) -> Result<ChemistryDefs> {
    let metadata = args.to_metadata();
    // Handle manually-specified chemistry.
    Ok(match unpack_chemistry_specs(&args.chemistry_specs)? {
        UnpackedChemistrySpecs::Manual(chems) => use_manual_chemistries(
            &chems,
            &metadata,
            multi_config_csv,
            args.custom_chemistry_def.as_ref(),
            units,
        )?,

        UnpackedChemistrySpecs::Auto(mode) => {
            match detect_chemistry(mode, &metadata, multi_config_csv, units) {
                Ok(chemistry_def) => {
                    // One chemistry was compatible with all units.
                    clone_chemistry_for_libraries(&chemistry_def, units)
                }

                Err(err) => {
                    // No single chemistry was compatible with all units.
                    // Attempt to detect chemistry per library type for RTL experiments.
                    let Some(DetectChemistryErrors::ConflictingChemistries {
                        per_unit_chems, ..
                    }) = err.downcast_ref::<DetectChemistryErrors>()
                    else {
                        bail!(err);
                    };

                    let is_rtl = per_unit_chems
                        .iter()
                        .flatten()
                        .all(|chem| chem.is_rtl().unwrap_or(false));
                    ensure!(is_rtl, err);

                    detect_chemistry_per_library_type(mode, &metadata, multi_config_csv, units)?
                }
            }
        }
    })
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
    custom_chemistry_def: Option<&ChemistryDef>,
    all_units: &[(DetectChemistryUnit, Vec<ReadPair>)],
) -> Result<ChemistryDefs> {
    divide_units_by_library_type(all_units)
        .map(|(library_type, units)| {
            let chem_name = chems[&library_type];
            if chem_name == ChemistryName::Custom {
                use_custom_chemistry(custom_chemistry_def)
            } else {
                use_manual_chemistry(chem_name, metadata, multi_config_csv, units)
            }
            .map(|chem_def| (library_type, chem_def))
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
    chem: ChemistryName,
    metadata: &Metadata<'_>,
    multi_config_csv: Option<&MultiConfigCsv>,
    units: &[(DetectChemistryUnit, Vec<ReadPair>)],
) -> Result<ChemistryDef> {
    let chems = &set![chem];

    if chem.is_rtl_multiplexed() {
        // Check the read length to ensure that the probe barcode is sequenced.
        let _matching_chemistries = LengthFilter::new(
            chems,
            chems,
            metadata.feature_reference,
            metadata.feature_config,
        )?
        .filter_chemistries(units)?;
    }

    if let Err(err) = WhitelistMatchFilter::new(chems, chems)?.filter_chemistries(units) {
        if chem.is_spatial() {
            bail!(err);
        }
        println!("WARNING: {err:#}");
    }

    validate_multiplexing(chem, metadata.sample_defs)?;
    validate_rtl(multi_config_csv, chem)?;
    Ok(ChemistryDef::named(chem))
}

/// Populate a ChemistryDefs map with a copy of the chemistry for each library type.
fn clone_chemistry_for_libraries(
    chemistry_def: &ChemistryDef,
    units: &[(DetectChemistryUnit, Vec<ReadPair>)],
) -> ChemistryDefs {
    units
        .iter()
        .map(|(unit, _reads)| unit.library_type)
        .unique()
        .zip(std::iter::repeat(chemistry_def.clone()))
        .collect()
}

/// Load a subsample of reads from each input unit.
fn load_reads(
    args: &DetectChemistryStageInputs,
) -> Result<Vec<(DetectChemistryUnit, Vec<ReadPair>)>> {
    let units = detect_chemistry_units(&args.sample_def, args.r1_length, args.r2_length)?;
    println!("Number of fastq units = {}", units.len());

    // Check for duplicate R1 and R2 files amongst units.
    check_fastq_identity(&units)?;

    units
        .into_iter()
        .sorted_by_key(|unit| unit.library_type)
        .map(|unit| {
            println!("Sampling reads from: {unit}");
            let reads = unit.sampled_read_pairs()?;
            Ok((unit, reads))
        })
        .try_collect()
}

/// Determine which chemistry to use for the specified auto detection mode.
fn detect_chemistry(
    mode: AutoChemistryName,
    metadata: &Metadata<'_>,
    multi_config_csv: Option<&MultiConfigCsv>,
    units: &[(DetectChemistryUnit, Vec<ReadPair>)],
) -> Result<ChemistryDef> {
    // Check for overhang multiplexing
    let is_overhang_multiplexed =
        multi_config_csv.is_some_and(MultiConfigCsv::is_overhang_multiplexed);

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

    let possible_chemistries = mode.allowed_chemistries(
        metadata.is_pd,
        multi_config_csv.is_some_and(|x| x.samples.is_some()),
        multi_config_csv.is_some_and(is_rtl_uncollapsed),
    );

    let allowed_chemistries = metadata.allowed_chems.map_or_else(
        || possible_chemistries.clone(),
        |chems| {
            chems
                .iter()
                .filter_map(AutoOrRefinedChemistry::refined)
                .collect()
        },
    );

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
        possible_chemistries.iter().sorted().format(", ")
    );

    // Read length based filtering
    let length_matching_chemistries = LengthFilter::new(
        &allowed_chemistries,
        &possible_chemistries,
        metadata.feature_reference,
        metadata.feature_config,
    )?
    .filter_chemistries(units)?;

    println!(
        "After length filter: {}",
        length_matching_chemistries.iter().sorted().format(", ")
    );

    // -----------------------------------------------------------------------------------------
    // For each unit of fastqs, ensure that a sufficient fraction of reads contain barcodes
    // which match the whitelist for at least one of the possible chemistries
    let wl_matching_chemistries =
        WhitelistMatchFilter::new(&allowed_chemistries, &length_matching_chemistries)?
            .filter_chemistries(units)?;

    let chosen_chemistry_def = if wl_matching_chemistries.len() == 1 {
        let &chem = wl_matching_chemistries.iter().exactly_one().unwrap();
        if chem == SFRP {
            // Bail out if a mixture of probe barcodes is observed for singleplex FRP chemistry
            validate_no_probe_bc_mixture_in_sfrp(
                units,
                metadata.feature_reference,
                metadata.feature_config,
            )?;
        }
        Some(validate_chemistry(
            chem,
            metadata.allowed_chems,
            metadata.sample_defs,
            multi_config_csv,
            is_overhang_multiplexed,
        )?)
    } else if wl_matching_chemistries == set![ThreePrimeV3, ThreePrimeV3HT, ThreePrimeV3LT] {
        Some(validate_chemistry(
            ThreePrimeV3LT,
            metadata.allowed_chems,
            metadata.sample_defs,
            multi_config_csv,
            is_overhang_multiplexed,
        )?)
    } else if wl_matching_chemistries == set![ThreePrimeV3, ThreePrimeV3HT] {
        Some(validate_chemistry(
            ThreePrimeV3,
            metadata.allowed_chems,
            metadata.sample_defs,
            multi_config_csv,
            is_overhang_multiplexed,
        )?)
    } else if wl_matching_chemistries == set![ThreePrimeV4, ThreePrimeV4HT] {
        Some(validate_chemistry(
            ThreePrimeV4,
            metadata.allowed_chems,
            metadata.sample_defs,
            multi_config_csv,
            is_overhang_multiplexed,
        )?)
    } else if wl_matching_chemistries == set![FivePrimeR2, FivePrimeHT] {
        Some(validate_chemistry(
            FivePrimeR2,
            metadata.allowed_chems,
            metadata.sample_defs,
            multi_config_csv,
            is_overhang_multiplexed,
        )?)
    } else if wl_matching_chemistries == set![FivePrimeR2V3, FivePrimeHTV3] {
        Some(validate_chemistry(
            FivePrimeR2V3,
            metadata.allowed_chems,
            metadata.sample_defs,
            multi_config_csv,
            is_overhang_multiplexed,
        )?)
    } else {
        None
    };
    if let Some(chosen_chemistry_def) = chosen_chemistry_def {
        return Ok(chosen_chemistry_def);
    }

    println!(
        "After whitelist filter: {}",
        wl_matching_chemistries.iter().sorted().format(", ")
    );

    let is_antibody_only = metadata
        .sample_defs
        .iter()
        .all(|sd| sd.library_type == Some(LibraryType::Antibody));
    let expected_mapping_chemistries = if is_antibody_only {
        let mut result = wl_matching_chemistries;
        // we define a new chemistry named "SC-FB" because this could be 3' v2 polyA capture
        // antibody library or a 5' antibody library
        let redundant_ab_chems = [ThreePrimeV2, FivePrimeR2, FivePrimePE];
        if redundant_ab_chems.iter().any(|chem| result.contains(chem)) {
            result.insert(FeatureBarcodingOnly);
            for chem in &redundant_ab_chems {
                result.remove(chem);
            }
        }
        result
    } else {
        let mut mapper = ReadMappingFilter::new(
            metadata.reference_path.unwrap(),
            &allowed_chemistries,
            wl_matching_chemistries,
        )?;
        mapper.filter_chemistries(units)?
    };

    println!(
        "After mapping filter: {}",
        expected_mapping_chemistries.iter().sorted().format(", ")
    );

    ensure!(
        expected_mapping_chemistries.len() == 1,
        "Could not distinguish between {}",
        expected_mapping_chemistries.iter().sorted().format(", ")
    );

    validate_chemistry(
        expected_mapping_chemistries
            .into_iter()
            .exactly_one()
            .unwrap(),
        metadata.allowed_chems,
        metadata.sample_defs,
        multi_config_csv,
        is_overhang_multiplexed,
    )
}

/// Run chemistry detection on every library independently.
/// This variant produces an independent result for each library without any
/// constraint on the individual chemistries being mutually compatible.
fn detect_chemistry_per_library_type(
    mode: AutoChemistryName,
    metadata: &Metadata<'_>,
    multi_config_csv: Option<&MultiConfigCsv>,
    all_units: &[(DetectChemistryUnit, Vec<ReadPair>)],
) -> Result<ChemistryDefs> {
    divide_units_by_library_type(all_units)
        .map(|(library_type, units)| {
            println!("\nDetecting chemistry for {library_type}");
            let chemistry = detect_chemistry(mode, metadata, multi_config_csv, units)?;
            ensure!(
                chemistry.name.compatible_with_library_type(library_type),
                "The chemistry {} was detected for {library_type} but they are not compatible; \
                 please check that your library configurations are associated with the correct \
                 library type.",
                chemistry.name,
            );
            println!("\nDetected chemistry {} for {library_type}", chemistry.name);
            anyhow::Ok((library_type, chemistry))
        })
        .try_collect()
}

/// Divide up the read units into slices based on library type.
fn divide_units_by_library_type(
    units: &[(DetectChemistryUnit, Vec<ReadPair>)],
) -> impl Iterator<Item = (LibraryType, &[(DetectChemistryUnit, Vec<ReadPair>)])> {
    units
        .linear_group_by_key(|(unit, _reads)| unit.library_type)
        .map(|units| (units[0].0.library_type, units))
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

/// Validate a chemistry.
fn validate_chemistry(
    chemistry: ChemistryName,
    allowed_chems: Option<&[AutoOrRefinedChemistry]>,
    sample_defs: &[SampleDef],
    multi_config_csv: Option<&MultiConfigCsv>,
    overhang_multiplexing: bool,
) -> Result<ChemistryDef> {
    if chemistry == ChemistryName::ThreePrimeV3LT {
        bail!("The chemistry SC3Pv3LT (Single Cell 3'v3 LT) is no longer supported. To analyze this data, use Cell Ranger 7.2 or earlier.");
    }
    if chemistry == ChemistryName::ArcV1 {
        bail!(
            "Cell Ranger detected the chemistry {chemistry}, which may indicate a workflow \
            error during sample preparation. Please check the reagents used to prepare this \
            sample and contact 10x Genomics support for further assistance. If this workflow is \
            intentional, you can force Cell Ranger to process this data by manually specifying \
            the {chemistry} chemistry in your analysis configuration."
        );
    }
    if chemistry == ChemistryName::ThreePrimeV4
        && sample_defs
            .iter()
            .any(|sdef| sdef.library_type == Some(LibraryType::Crispr))
    {
        bail!("The chemistry SC3Pv4 (Single Cell 3'v4) is not supported with CRISPR Guide Capture libraries.")
    }

    validate_multiplexing(chemistry, sample_defs)?;
    validate_rtl(multi_config_csv, chemistry)?;

    let chemistry_with_overhang = if overhang_multiplexing {
        chemistry.get_overhang_version()?
    } else {
        chemistry
    };

    validate_chemistry_spec(
        AutoOrRefinedChemistry::Refined(chemistry_with_overhang),
        allowed_chems,
    )?;

    Ok(ChemistryDef::named(chemistry_with_overhang))
}

/// Validate combinations of chemistry and types of multiplexing.
fn validate_multiplexing(chemistry_type: ChemistryName, sample_defs: &[SampleDef]) -> Result<()> {
    if !sample_defs
        .iter()
        .any(|sdef| sdef.library_type == Some(LibraryType::Cellplex))
    {
        return Ok(());
    }

    match chemistry_type {
        ThreePrimeV3LT => ensure!(
            *threeprime_lt_multiplexing()?,
            "Multiplexing Capture libraries are not supported with Single Cell 3' v3 LT chemistry"
        ),
        FivePrimeR1 | FivePrimeR2 | FivePrimePE => ensure!(
            *fiveprime_multiplexing()?,
            "Multiplexing Capture libraries are not supported with Single Cell 5' chemistries"
        ),
        _ => (),
    }
    Ok(())
}

/// Validate the chemistry with RTL-related parameters.
fn validate_rtl(multi_config_csv: Option<&MultiConfigCsv>, chemistry: ChemistryName) -> Result<()> {
    let Some(config) = multi_config_csv else {
        // cellranger count does not support RTL chemistries.
        return Ok(());
    };

    if !chemistry.is_rtl().expect("is_rtl is None only for custom and spatial chemistries, which do not use chemistry detection") {
        let Some(samples) = &config.samples else {
            return Ok(());
        };
        ensure!(
            !samples.has_probe_barcode_ids(),
            "A non-Fixed RNA Profiling chemistry {chemistry} was detected, and the [samples] \
             section has a probe_barcode_ids column. The probe_barcode_ids column may only be \
             specified with Fixed RNA Profiling chemistries.",
        );
        return Ok(());
    }

    if let Some(gex) = &config.gene_expression {
        if config.libraries.has_gene_expression() {
            ensure!(
                gex.probe_set.is_some(),
                "Fixed RNA Profiling chemistries require a probe-set."
            );
        }
        ensure!(
            gex.include_introns == multiconst::DEFAULT_INCLUDE_INTRONS,
            "The [gene-expression] section specifies the parameter include-introns, \
             which is not valid for Fixed RNA Profiling chemistries."
        );
    }

    if chemistry == SFRP {
        ensure!(
            config.samples.is_none(),
            "We detected singleplex Fixed RNA Profiling chemistry from the data. \
             Sample definitions are unsupported for singleplex inputs. \
             To process this data as multiplex Fixed RNA Profiling you will need to specify `MFRP` \
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
        .all(|x| x.ends_with(&['A', 'B', 'C', 'D']))
}

/// Identify probe barcode pairings if necessary.
/// If pairings are present, modify the translation whitelist to include those pairings.
fn handle_probe_barcode_translation(
    rover: &MartianRover,
    multi_config_csv: &MultiConfigCsv,
    units: &[(DetectChemistryUnit, Vec<ReadPair>)],
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
            .barcode_whitelist()
            .probe()
            .as_source(true)?;
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
    use cr_types::chemistry::{known_chemistry_defs, ChemistryDefsExt};
    use cr_types::sample_def::FastqMode;
    use cr_types::LibraryType;
    use dui_tests::stage_test::StageFailTest;
    use dui_tests::{stage_fail_dui_test, DuiTest};
    use multi::config::{GemWell, GeneExpressionParams, Lanes, LibrariesCsv, Library, SampleRow};
    use std::io::Write;
    use std::path::Path;
    use std::string::ToString;

    #[test]
    fn test_validate_rtl() {
        use ChemistryName::{MFRP_Ab, ThreePrimeV3, MFRP_RNA, SFRP};

        assert!(validate_rtl(None, ThreePrimeV3).is_ok());

        let gex_library = Library::BclProcessor {
            fastq_path: PathBuf::default(),
            sample_indices: Vec::default(),
            lanes: Lanes::Any,
            physical_library_id: String::default(),
            feature_type: LibraryType::Gex,
            gem_well: GemWell(0),
            subsample_rate: None,
            chemistry: None,
        };

        let antibody_library = Library::BclProcessor {
            fastq_path: PathBuf::default(),
            sample_indices: Vec::default(),
            lanes: Lanes::Any,
            physical_library_id: String::default(),
            feature_type: LibraryType::Antibody,
            gem_well: GemWell(0),
            subsample_rate: None,
            chemistry: None,
        };

        let gex_without_probe_set = MultiConfigCsv {
            libraries: LibrariesCsv(vec![gex_library]),
            gene_expression: Some(GeneExpressionParams {
                probe_set: None,
                include_introns: multiconst::DEFAULT_INCLUDE_INTRONS,
                ..GeneExpressionParams::default()
            }),
            ..MultiConfigCsv::default()
        };
        assert!(validate_rtl(Some(&gex_without_probe_set), MFRP_RNA).is_err());
        assert!(validate_rtl(Some(&gex_without_probe_set), SFRP).is_err());
        assert!(validate_rtl(Some(&gex_without_probe_set), ThreePrimeV3).is_ok());

        let gex_with_probe_set = MultiConfigCsv {
            gene_expression: Some(GeneExpressionParams {
                probe_set: Some(PathBuf::new()),
                include_introns: multiconst::DEFAULT_INCLUDE_INTRONS,
                ..GeneExpressionParams::default()
            }),
            ..gex_without_probe_set
        };
        assert!(validate_rtl(Some(&gex_with_probe_set), MFRP_RNA).is_ok());
        assert!(validate_rtl(Some(&gex_with_probe_set), SFRP).is_ok());
        assert!(validate_rtl(Some(&gex_with_probe_set), ThreePrimeV3).is_ok()); // to aggr with RTL

        let antibody_only = MultiConfigCsv {
            libraries: LibrariesCsv(vec![antibody_library]),
            gene_expression: Some(GeneExpressionParams {
                probe_set: None,
                include_introns: multiconst::DEFAULT_INCLUDE_INTRONS,
                ..GeneExpressionParams::default()
            }),
            ..MultiConfigCsv::default()
        };
        // This assertion might seem non-sensical, but the validation function
        // under test doesn't care about the library type, so it should not fail
        // for Gex even if the config is for antibody-only.
        assert!(validate_rtl(Some(&antibody_only), MFRP_RNA).is_ok());
        assert!(validate_rtl(Some(&antibody_only), MFRP_Ab).is_ok());
        assert!(validate_rtl(Some(&antibody_only), SFRP).is_ok());
        assert!(validate_rtl(Some(&antibody_only), ThreePrimeV3).is_ok());

        let include_introns_false = MultiConfigCsv {
            gene_expression: Some(GeneExpressionParams {
                probe_set: Some(PathBuf::new()),
                include_introns: false,
                ..GeneExpressionParams::default()
            }),
            ..MultiConfigCsv::default()
        };
        assert!(validate_rtl(Some(&include_introns_false), MFRP_RNA).is_err());
        assert!(validate_rtl(Some(&include_introns_false), SFRP).is_err());
        assert!(validate_rtl(Some(&include_introns_false), ThreePrimeV3).is_ok()); // to aggr with RTL

        let samples = MultiConfigCsv {
            samples: Some(SamplesCsv(vec![SampleRow::from_probe_barcode_id("BC001")])),
            ..gex_with_probe_set
        };
        assert!(validate_rtl(Some(&samples), MFRP_RNA).is_ok());
        assert!(validate_rtl(Some(&samples), SFRP).is_err());
        assert!(validate_rtl(Some(&samples), ThreePrimeV3).is_err());
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
        DetectChemistryStageInputs {
            sample_def: vec![SampleDef {
                fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                lanes,
                read_path: fastq_path.as_ref().to_path_buf(),
                sample_names: Some(sample_names.into_iter().map(String::from).collect()),
                ..SampleDef::default()
            }],
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
            allowed_chems: Some(vec![ChemistryName::ThreePrimeV3.into()]),
            ..Default::default()
        };
        err_snapshot!(args);
    }

    #[test]
    fn test_custom_chemistry_spec() {
        let mut def: ChemistryDef = ChemistryDef::named(ChemistryName::ThreePrimeV3);
        def.name = ChemistryName::Custom;
        def.rna.offset = 13;
        let args = DetectChemistryStageInputs {
            chemistry_specs: specs_input_gex(ChemistryName::Custom),
            custom_chemistry_def: Some(def.clone()),
            sample_def: vec![SampleDef {
                fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                read_path: "../dui_tests/test_resources/cellranger-count/CR-300-01_SC3pv3_15k"
                    .into(),
                sample_names: Some(vec!["CR-300-01".into()]),
                ..Default::default()
            }],
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(outs.chemistry_defs.values().exactly_one().unwrap(), &def);
    }

    #[test]
    fn test_mix_of_custom_and_manual_chemistry_spec() {
        let mut def: ChemistryDef = ChemistryDef::named(ChemistryName::ThreePrimeV3);
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
            custom_chemistry_def: Some(def.clone()),
            sample_def: vec![
                SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path: "../dui_tests/test_resources/cellranger-count/CR-300-01_SC3pv3_15k"
                        .into(),
                    sample_names: Some(vec!["CR-300-01".into()]),
                    ..Default::default()
                },
                SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path:
                        "../dui_tests/test_resources/cellranger-count/pbmc_1k_protein_v3_antibody"
                            .into(),
                    sample_names: Some(vec!["pbmc_1k_protein_v3_antibody".into()]),
                    library_type: Some(LibraryType::Antibody),
                    ..Default::default()
                },
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
                custom_chemistry_def: Some(custom_def),
                sample_def: vec![SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path: "../dui_tests/test_resources/cellranger-count/CR-300-01_SC3pv3_15k"
                        .into(),
                    sample_names: Some(vec!["CR-300-01".into()]),
                    ..Default::default()
                }],
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
            sample_def: vec![SampleDef {
                fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                read_path: "../dui_tests/test_resources/cellranger-count/CR-300-01_SC3pv3_15k"
                    .into(),
                sample_names: Some(vec!["CR-300-01".into()]),
                ..Default::default()
            }],
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
            sample_def: vec![SampleDef {
                fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                read_path: "../dui_tests/test_resources/cellranger-count/CR-300-01_SC3pv3_15k"
                    .into(),
                sample_names: Some(vec!["CR-300-01".into()]),
                ..Default::default()
            }],
            chemistry_specs: specs_input_gex(ChemistryName::ThreePrimeV3),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(
            outs.chemistry_defs.primary().name,
            ChemistryName::ThreePrimeV3
        );
        assert!(!outs.is_antibody_only);
    }

    #[test]
    fn test_auto_chemistry_spec_3pv3() {
        let args = DetectChemistryStageInputs {
            sample_def: vec![SampleDef {
                fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                read_path: "../dui_tests/test_resources/cellranger-count/CR-300-01_SC3pv3_15k"
                    .into(),
                sample_names: Some(vec!["CR-300-01".into()]),
                ..Default::default()
            }],
            chemistry_specs: specs_input_gex(AutoChemistryName::Count),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(
            outs.chemistry_defs.primary().name,
            ChemistryName::ThreePrimeV3
        );
        assert!(!outs.is_antibody_only);
    }

    #[test]
    fn test_correct_chemistry_spec_3pv3ht() {
        // should detect ThreePrimeV3HT when SCPV3HT is specified
        let args = DetectChemistryStageInputs {
            sample_def: vec![SampleDef {
                fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                read_path: "../dui_tests/test_resources/cellranger-count/HT_chemistry".into(),
                sample_names: Some(vec!["1107522".into()]),
                ..Default::default()
            }],
            chemistry_specs: specs_input_gex(ChemistryName::ThreePrimeV3HT),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(
            outs.chemistry_defs.primary().name,
            ChemistryName::ThreePrimeV3HT
        );
        assert!(!outs.is_antibody_only);
    }

    #[test]
    fn test_auto_chemistry_spec_ht() {
        // should detect ThreePrimeV3 when HT fastqs are run
        let args = DetectChemistryStageInputs {
            sample_def: vec![SampleDef {
                fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                read_path: "../dui_tests/test_resources/cellranger-count/HT_chemistry".into(),
                sample_names: Some(vec!["1107522".into()]),
                ..Default::default()
            }],
            chemistry_specs: specs_input_gex(AutoChemistryName::Count),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(
            outs.chemistry_defs.primary().name,
            ChemistryName::ThreePrimeV3
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
                    SampleDef {
                        fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                        read_path: "../dui_tests/test_resources/cellranger-count/CR-300-01_SC3pv3_15k"
                            .into(),
                        sample_names: Some(vec!["CR-300-01".into()]),
                        ..Default::default()
                    },
                    SampleDef {
                        fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                        read_path:
                            "../dui_tests/test_resources/cellranger-count/CR-CUSTOM-22_SC5pr2_15k"
                                .into(),
                        sample_names: Some(vec!["CR-CUSTOM-22".into()]),
                        ..Default::default()
                    },
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
            sample_def: vec![SampleDef {
                fastq_mode: FastqMode::BCL_PROCESSOR,
                read_path: "../dui_tests/test_resources/cellranger-count/cycle_failure_bc_v3"
                    .into(),
                sample_indices: Some(vec!["SI-P2-D3".into()]),
                ..Default::default()
            }],
            chemistry_specs: specs_input_gex(AutoChemistryName::Count),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(
            outs.chemistry_defs.primary().name,
            ChemistryName::ThreePrimeV3
        );
        assert!(!outs.is_antibody_only);
    }

    #[test]
    fn test_vdj_pe() {
        let args = DetectChemistryStageInputs {
            sample_def: vec![SampleDef {
                fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                read_path: "../dui_tests/test_resources/cellranger-vdj/vdj_v1_hs_pbmc3_t".into(),
                sample_names: Some(vec!["vdj_v1_hs_pbmc3_t".into()]),
                ..Default::default()
            }],
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
            sample_def: vec![SampleDef {
                fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                read_path: "../dui_tests/test_resources/cellranger-vdj/vdj_v1_mm_pbmc4_b_26x91"
                    .into(),
                sample_names: Some(vec!["vdj_v1_mm_pbmc4_b_26x91".into()]),
                ..Default::default()
            }],
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
            sample_def: vec![SampleDef {
                fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                read_path:
                    "../dui_tests/test_resources/cellranger-count/pbmc_1k_protein_v3_antibody"
                        .into(),
                sample_names: Some(vec!["pbmc_1k_protein_v3_antibody".into()]),
                lanes: Some(vec![4]),
                library_type: Some(LibraryType::Antibody),
                ..Default::default()
            }],
            chemistry_specs: specs_input_gex(AutoChemistryName::Count),
            feature_reference: Some("../dui_tests/test_resources/reference/feature/Biolegend_x22_TotalSeqB_14Ab_Isotype_Controls_20180601.csv".into()),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(
            outs.chemistry_defs.primary().name,
            ChemistryName::ThreePrimeV3
        );
        assert!(outs.is_antibody_only);
    }

    #[test]
    fn test_ab_short_only_3pv3() {
        let args = DetectChemistryStageInputs {
            sample_def: vec![SampleDef {
                fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                read_path:
                    "../dui_tests/test_resources/cellranger-count/pbmc_1k_protein_v3_antibody_len24"
                        .into(),
                sample_names: Some(vec!["pbmc_1k_protein_v3_antibody".into()]),
                lanes: Some(vec![4]),
                library_type: Some(LibraryType::Antibody),
                ..Default::default()
            }],
            chemistry_specs: specs_input_gex(AutoChemistryName::Count),
            feature_reference: Some("../dui_tests/test_resources/reference/feature/Biolegend_x22_TotalSeqB_14Ab_Isotype_Controls_len24_20210721.csv".into()),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(
            outs.chemistry_defs.primary().name,
            ChemistryName::ThreePrimeV3
        );
        assert!(outs.is_antibody_only);
    }

    #[test]
    fn test_ab_only_5p() {
        let args = DetectChemistryStageInputs {
            sample_def: vec![SampleDef {
                fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                read_path: "../dui_tests/test_resources/cellranger-count/vdj_v1_hs_pbmc2_antibody"
                    .into(),
                sample_names: Some(vec!["vdj_v1_hs_pbmc2_antibody".into()]),
                lanes: Some(vec![1]),
                library_type: Some(LibraryType::Antibody),
                ..Default::default()
            }],
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
            sample_def: vec![
                SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path:
                        "../dui_tests/test_resources/cellranger-count/pbmc_1k_protein_v3_antibody"
                            .into(),
                    sample_names: Some(vec!["pbmc_1k_protein_v3_antibody".into()]),
                    lanes: Some(vec![4]),
                    library_type: Some(LibraryType::Antibody),
                    ..Default::default()
                },
                SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path:
                        "../dui_tests/test_resources/cellranger-count/pbmc_1k_protein_v3_gex".into(),
                    sample_names: Some(vec!["pbmc_1k_protein_v3_gex".into()]),
                    lanes: Some(vec![4]),
                    library_type: Some(LibraryType::Gex),
                    ..Default::default()
                },
            ],
            chemistry_specs: specs_input_gex(AutoChemistryName::Count),
            feature_reference: Some("../dui_tests/test_resources/reference/feature/Biolegend_x22_TotalSeqB_14Ab_Isotype_Controls_20180601.csv".into()),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(
            outs.chemistry_defs.primary().name,
            ChemistryName::ThreePrimeV3
        );
        assert!(!outs.is_antibody_only);
    }

    #[test]
    fn test_ab_short_and_gex_3pv3() {
        let args = DetectChemistryStageInputs {
            sample_def: vec![
                SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path:
                        "../dui_tests/test_resources/cellranger-count/pbmc_1k_protein_v3_antibody_len24"
                            .into(),
                    sample_names: Some(vec!["pbmc_1k_protein_v3_antibody".into()]),
                    lanes: Some(vec![4]),
                    library_type: Some(LibraryType::Antibody),
                    ..Default::default()
                },
                SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path:
                        "../dui_tests/test_resources/cellranger-count/pbmc_1k_protein_v3_gex".into(),
                    sample_names: Some(vec!["pbmc_1k_protein_v3_gex".into()]),
                    lanes: Some(vec![4]),
                    library_type: Some(LibraryType::Gex),
                    ..Default::default()
                },
            ],
            chemistry_specs: specs_input_gex(AutoChemistryName::Count),
            feature_reference: Some("../dui_tests/test_resources/reference/feature/Biolegend_x22_TotalSeqB_14Ab_Isotype_Controls_len24_20210721.csv".into()),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(
            outs.chemistry_defs.primary().name,
            ChemistryName::ThreePrimeV3
        );
        assert!(!outs.is_antibody_only);
    }

    #[test]
    fn test_ab_too_short_and_gex_3pv3() {
        let args = DetectChemistryStageInputs {
            sample_def: vec![
                SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path:
                        "../dui_tests/test_resources/cellranger-count/pbmc_1k_protein_v3_antibody"
                            .into(),
                    sample_names: Some(vec!["pbmc_1k_protein_v3_antibody".into()]),
                    lanes: Some(vec![4]),
                    library_type: Some(LibraryType::Antibody),
                    ..Default::default()
                },
                SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path:
                        "../dui_tests/test_resources/cellranger-count/pbmc_1k_protein_v3_gex".into(),
                    sample_names: Some(vec!["pbmc_1k_protein_v3_gex".into()]),
                    lanes: Some(vec![4]),
                    library_type: Some(LibraryType::Gex),
                    ..Default::default()
                },
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
                SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path:
                        "../dui_tests/test_resources/cellranger-count/K562_5k_crispr_v3_crispr"
                            .into(),
                    sample_names: Some(vec!["K562_5k_crispr_v3_crispr".into()]),
                    library_type: Some(LibraryType::Crispr),
                    ..Default::default()
                },
                SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path: "../dui_tests/test_resources/cellranger-count/K562_5k_crispr_v3_gex"
                        .into(),
                    sample_names: Some(vec!["K562_5k_crispr_v3_gex".into()]),
                    library_type: Some(LibraryType::Gex),
                    ..Default::default()
                },
            ],
            chemistry_specs: specs_input_gex(AutoChemistryName::Count),
            feature_reference: Some(
                "../dui_tests/test_resources/reference/feature/RAB1A_NT_Aug2020.csv".into(),
            ),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(
            outs.chemistry_defs.primary().name,
            ChemistryName::ThreePrimeV3
        );
        assert!(!outs.is_antibody_only);
    }

    #[test]
    fn test_crispr_too_short_and_gex_3pv3() {
        let args = DetectChemistryStageInputs {
            sample_def: vec![
                SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path:
                        "../dui_tests/test_resources/cellranger-count/K562_5k_crispr_v3_crispr"
                            .into(),
                    sample_names: Some(vec!["K562_5k_crispr_v3_crispr".into()]),
                    library_type: Some(LibraryType::Crispr),
                    ..Default::default()
                },
                SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path: "../dui_tests/test_resources/cellranger-count/K562_5k_crispr_v3_gex"
                        .into(),
                    sample_names: Some(vec!["K562_5k_crispr_v3_gex".into()]),
                    library_type: Some(LibraryType::Gex),
                    ..Default::default()
                },
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
            sample_def: vec![SampleDef {
                fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                read_path: "../dui_tests/test_resources/cellranger-count/CR-CUSTOM-02_SC3pv1_15k"
                    .into(),
                sample_names: Some(vec!["test_sample".into()]),
                ..Default::default()
            }],
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
                SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path:
                        "../dui_tests/test_resources/cellranger-count/CR-CUSTOM-02_SC3pv1_15k"
                            .into(),
                    sample_names: Some(vec!["test_sample".into()]),
                    ..Default::default()
                },
                SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path: "../dui_tests/test_resources/cellranger-count/CR-120-01_SC3pv2_15k"
                        .into(),
                    sample_names: Some(vec!["CR-120-01".into()]),
                    ..Default::default()
                },
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
                sample_def: vec![SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path: "../dui_tests/test_resources/cellranger-count/no_barcode_div_v3".into(),
                    sample_names: Some(vec!["test_sample".into()]),
                    ..SampleDef::default()
                }],
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
                sample_def: vec![SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path: "../dui_tests/test_resources/cellranger-count/reads_too_short_v2".into(),
                    sample_names: Some(vec!["CR-120-01".into()]),
                    ..SampleDef::default()
                }],
                chemistry_specs: specs_input_gex(AutoChemistryName::Count),
                ..Default::default()
            },
        );
        Ok(())
    }

    #[test]
    fn test_vdj_pe_with_r1_r2_length() {
        let args = DetectChemistryStageInputs {
            sample_def: vec![SampleDef {
                fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                read_path: "../dui_tests/test_resources/cellranger-vdj/vdj_v1_hs_pbmc3_t".into(),
                sample_names: Some(vec!["vdj_v1_hs_pbmc3_t".into()]),
                ..Default::default()
            }],
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
            sample_def: vec![SampleDef {
                fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                read_path: "../dui_tests/test_resources/cellranger-count/CR-300-01_SC3pv3_15k"
                    .into(),
                sample_names: Some(vec!["CR-300-01".into()]),
                ..Default::default()
            }],
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
                sample_def: vec![SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path: "../dui_tests/test_resources/cellranger-count/CR-300-01_SC3pv3_15k"
                        .into(),
                    sample_names: Some(vec!["CR-300-01".into()]),
                    ..Default::default()
                }],
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
                sample_def: vec![SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path: "../dui_tests/test_resources/cellranger-count/CR-300-01_SC3pv3_15k"
                        .into(),
                    sample_names: Some(vec!["CR-300-01".into()]),
                    r1_length: Some(25), // Trimmed too much in sample_def
                    ..Default::default()
                }],
                chemistry_specs: specs_input_gex(AutoChemistryName::Count),
                ..Default::default()
            },
        );
        Ok(())
    }

    #[test]
    fn test_r1_r2_okay_3pv3() {
        let args = DetectChemistryStageInputs {
            sample_def: vec![SampleDef {
                fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                read_path: "../dui_tests/test_resources/cellranger-count/CR-300-01_SC3pv3_15k"
                    .into(),
                sample_names: Some(vec!["CR-300-01".into()]),
                ..Default::default()
            }],
            chemistry_specs: specs_input_gex(AutoChemistryName::Count),
            r1_length: Some(26),
            r2_length: Some(75),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(
            outs.chemistry_defs.primary().name,
            ChemistryName::ThreePrimeV3
        );
        assert!(!outs.is_antibody_only);
    }

    #[test]
    fn test_few_reads_explicit_chem() {
        let args = DetectChemistryStageInputs {
            sample_def: vec![SampleDef {
                fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                read_path: "../dui_tests/test_resources/cellranger-count/sc3pv2_5k".into(),
                sample_names: Some(vec!["tinygex".into()]),
                ..Default::default()
            }],
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
            sample_def: vec![SampleDef {
                fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                read_path:
                    "../dui_tests/test_resources/cellranger-count/RTL_single_sample_probe_bcs_mixture"
                        .into(),
                sample_names: Some(vec!["test_sample".into()]),
                ..Default::default()
            }],
            chemistry_specs: specs_input_gex(AutoChemistryName::Count),
            ..Default::default()
        };
        err_snapshot!(args);
    }

    #[test]
    // When refined SFRP chemistry is specified, check for probe bc mixtures is skipped
    fn test_rtl_singleplex_probe_bc_mixture_explicit_chem() {
        let args = DetectChemistryStageInputs {
            sample_def: vec![SampleDef {
                fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                read_path:
                    "../dui_tests/test_resources/cellranger-count/RTL_single_sample_probe_bcs_mixture"
                        .into(),
                sample_names: Some(vec!["test_sample".into()]),
                library_type: Some(LibraryType::Gex),
                target_set: Some("../dui_tests/test_resources/cellranger-count/RTL_single_sample_probe_bcs_mixture/probe_set.csv".into()),
                ..Default::default()
            }],
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
            sample_def: vec![SampleDef {
                fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                read_path: "../dui_tests/test_resources/cellranger-count/CR-300-01_SC3pv3_15k"
                    .into(),
                sample_names: Some(vec!["CR-300-01".into()]),
                ..Default::default()
            }],
            chemistry_specs: specs_input_gex(ChemistryName::SFRP),
            multi_config: Some(MultiConfigCsvFile::from(multi_config_csv.path())),
            ..Default::default()
        };
        err_snapshot!(args);
    }

    #[test]
    fn test_overhang_chemistry() {
        let args = DetectChemistryStageInputs {
            sample_def: vec![SampleDef {
                fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                read_path: "../dui_tests/test_resources/cellranger-count/CR-300-01_SC3pv3_15k"
                    .into(),
                sample_names: Some(vec!["CR-300-01".into()]),
                ..Default::default()
            }],
            chemistry_specs: specs_input_gex(AutoChemistryName::Count),
            multi_config: Some(
                "../dui_tests/test_resources/cellranger-multi/overhang_3pv3.csv".into(),
            ),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(
            outs.chemistry_defs.primary().name,
            ChemistryName::ThreePrimeV3OH
        );
    }

    #[cfg(feature = "slow_tests")]
    mod tests_requiring_alignment {
        use super::*;
        use cr_types::sample_def::FastqMode;
        use dui_tests::stage_test::StageFailTest;
        use dui_tests::{stage_fail_dui_test, DuiTest};
        use test_refdata::{refdata_path, sere_path, testdata_path};

        #[test]
        fn test_auto_chemistry_spec_3pv2() {
            let args = DetectChemistryStageInputs {
                sample_def: vec![SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path: "../dui_tests/test_resources/cellranger-count/CR-120-01_SC3pv2_15k"
                        .into(),
                    sample_names: Some(vec!["CR-120-01".into()]),
                    ..Default::default()
                }],
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
                sample_def: vec![SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path:
                        "../dui_tests/test_resources/cellranger-count/CR-CUSTOM-22_SC5pr2_15k"
                            .into(),
                    sample_names: Some(vec!["CR-CUSTOM-22".into()]),
                    ..Default::default()
                }],
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
                sample_def: vec![SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path:
                        "../dui_tests/test_resources/cellranger-count/CR-CUSTOM-22_SC5pr2_15k"
                            .into(),
                    sample_names: Some(vec!["CR-CUSTOM-22".into()]),
                    ..Default::default()
                }],
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
                sample_def: vec![SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path:
                        "../dui_tests/test_resources/cellranger-count/CR-CUSTOM-17_SC5p-PE_15k"
                            .into(),
                    sample_names: Some(vec!["CR-CUSTOM-17".into()]),
                    ..Default::default()
                }],
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
                    SampleDef {
                        fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                        read_path:
                            "../dui_tests/test_resources/cellranger-count/CR-CUSTOM-22_SC5pr2_15k"
                                .into(),
                        sample_names: Some(vec!["CR-CUSTOM-22".into()]),
                        ..Default::default()
                    },
                    SampleDef {
                        fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                        read_path:
                            "../dui_tests/test_resources/cellranger-count/CR-CUSTOM-17_SC5p-PE_15k"
                                .into(),
                        sample_names: Some(vec!["CR-CUSTOM-17".into()]),
                        ..Default::default()
                    },
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
                    SampleDef {
                        fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                        read_path:
                            "../dui_tests/test_resources/cellranger-count/CR-CUSTOM-17_SC5p-PE_15k"
                                .into(),
                        sample_names: Some(vec!["CR-CUSTOM-17".into()]),
                        ..Default::default()
                    },
                    SampleDef {
                        fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                        read_path:
                            "../dui_tests/test_resources/cellranger-count/CR-CUSTOM-22_SC5pr2_15k"
                                .into(),
                        sample_names: Some(vec!["CR-CUSTOM-22".into()]),
                        ..Default::default()
                    },
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
            sample_def: vec![SampleDef {
                fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                read_path: sere_path("vdj_nextgem_hs_pbmc3_5gex_protein/vdj_nextgem_hs_pbmc3_5gex_protein_antibody"),
                sample_names: Some(vec!["vdj_nextgem_hs_pbmc3_5gex_protein_antibody".into()]),
                lanes: Some(vec![1]),
                library_type: Some(LibraryType::Antibody),
                ..Default::default()
            }, SampleDef {
                fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                read_path: sere_path("vdj_nextgem_hs_pbmc3_5gex_protein/vdj_nextgem_hs_pbmc3_5gex_protein_gex"),
                sample_names: Some(vec!["vdj_nextgem_hs_pbmc3_5gex_protein_gex".into()]),
                lanes: Some(vec![1]),
                library_type: Some(LibraryType::Gex),
                ..Default::default()
            }],
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
                SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path: testdata_path("fastqs/cellranger/multi/1247954_1270332_frp_gex_and_ab_fastqs"),
                    sample_names: Some(vec!["1247954_GEX".into()]),
                    target_set: Some(refdata_path("targeted_panels/human_wta_RTL.Fletcher_v7.0.all_manufacturing.mismatch_base_pairing_revision.csv")),
                    lanes: Some(vec![4]),
                    library_type: Some(LibraryType::Gex),
                    ..Default::default()
                },
                SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path: testdata_path("fastqs/cellranger/multi/1247954_1270332_frp_gex_and_ab_fastqs"),
                    sample_names: Some(vec!["1247954_AB".into()]),
                    lanes: Some(vec![4]),
                    library_type: Some(LibraryType::Antibody),
                    ..Default::default()
                }
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
            sample_def: vec![SampleDef {
                fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                read_path: testdata_path("fastqs/cellranger/multi/1245033_sfrp_fastqs"),
                sample_names: Some(vec!["CTACTGAATT".into()]),
                target_set: Some(refdata_path("targeted_panels/human_wta_RTL.Fletcher_v7.0.all_manufacturing.mismatch_base_pairing_revision.csv")),
                lanes: Some(vec![4]),
                library_type: Some(LibraryType::Gex),
                ..Default::default()
            }],
            multi_config: Some("../dui_tests/test_resources/cellranger-multi/1245033_sfrp.csv".into()),
            chemistry_specs: specs_input_gex(AutoChemistryName::Count),
            ..Default::default()
        };
            let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
            assert_eq!(outs.chemistry_defs.primary().name, ChemistryName::SFRP);
        }

        #[test]
        fn test_gex_only_mfrp() {
            let args = DetectChemistryStageInputs {
            sample_def: vec![SampleDef {
                fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                read_path: testdata_path("fastqs/cellranger/multi/1216790_1245489_mfrp_fastqs"),
                sample_names: Some(vec!["AATTTCGGGT".into()]),
                target_set: Some(refdata_path("targeted_panels/human_wta_RTL.Fletcher_v7.0.all_manufacturing.mismatch_base_pairing_revision.csv")),
                lanes: Some(vec![4]),
                library_type: Some(LibraryType::Gex),
                ..Default::default()
            }],
            multi_config: Some("../dui_tests/test_resources/cellranger-multi/1216790_1245489.csv".into()),
            chemistry_specs: specs_input_gex(AutoChemistryName::Count),
            ..Default::default()
        };
            let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
            assert_eq!(outs.chemistry_defs.primary().name, ChemistryName::MFRP_RNA);
        }

        #[test]
        fn test_mfrp_r1() {
            let args = DetectChemistryStageInputs{
            sample_def: vec![SampleDef {
                fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                read_path: testdata_path("fastqs/cellranger/multi/mfrp_r1_50_90"),
                sample_names: Some(vec!["tiny_gex".into()]),
                target_set: Some(refdata_path("targeted_panels/human_wta_RTL.Fletcher_v7.0.all_manufacturing.mismatch_base_pairing_revision.csv")),
                lanes: Some(vec![1]),
                library_type: Some(LibraryType::Gex),
                ..Default::default()
            }],
            multi_config: Some("../dui_tests/test_resources/cellranger-multi/1320608_mfrp_r1_50_90.csv".into()),
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
            let args = DetectChemistryStageInputs{
            sample_def: vec![SampleDef {
                fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                read_path: testdata_path("fastqs/cellranger/multi/1461767_mfrp_gex_ab/ab"),
                sample_names: Some(vec!["tiny_ab".into()]),
                target_set: Some(refdata_path("targeted_panels/human_wta_RTL.Fletcher_v7.0.all_manufacturing.mismatch_base_pairing_revision.csv")),
                lanes: Some(vec![1]),
                library_type: Some(LibraryType::Gex),
                ..Default::default()
            }],
            multi_config: Some("../dui_tests/test_resources/cellranger-multi/1461767_mfrp_ab_only.csv".into()),
            chemistry_specs: specs_input_gex(AutoChemistryName::Count),
            ..Default::default()
            };
            let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
            assert_eq!(outs.chemistry_defs.primary().name, ChemistryName::MFRP_Ab);
        }

        #[test]
        fn test_ab_only_mfrp_r1() {
            let args = DetectChemistryStageInputs{
            sample_def: vec![SampleDef {
                fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                read_path: testdata_path("fastqs/cellranger/multi/1464930_mfrp_r1_gex_ab/ab"),
                sample_names: Some(vec!["tiny_ab".into()]),
                target_set: Some(refdata_path("targeted_panels/human_wta_RTL.Fletcher_v7.0.all_manufacturing.mismatch_base_pairing_revision.csv")),
                lanes: Some(vec![1]),
                library_type: Some(LibraryType::Gex),
                ..Default::default()
            }],
            multi_config: Some("../dui_tests/test_resources/cellranger-multi/1464930_mfrp_r1_ab_only.csv".into()),
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
                SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path: testdata_path("fastqs/cellranger/multi/1548572/tiny_ab"),
                    r1_length: Some(28),
                    sample_names: Some(vec!["tiny_ab".into()]),
                    target_set: Some(refdata_path("targeted_panels/human_wta_RTL.Fletcher_v7.0.all_manufacturing.mismatch_base_pairing_revision.csv")),
                    lanes: Some(vec![1]),
                    library_type: Some(LibraryType::Antibody),
                    ..Default::default()
                },
                SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path: testdata_path("fastqs/cellranger/multi/1548572/tiny_gex"),
                    sample_names: Some(vec!["tiny_gex".into()]),
                    target_set: Some(refdata_path("targeted_panels/human_wta_RTL.Fletcher_v7.0.all_manufacturing.mismatch_base_pairing_revision.csv")),
                    lanes: Some(vec![1]),
                    library_type: Some(LibraryType::Gex),
                    ..Default::default()
                }
            ],
            multi_config: Some("../dui_tests/test_resources/cellranger-multi/1548572_mfrp_r1_ab_r2pos50.csv".into()),
            feature_reference: Some("../dui_tests/test_resources/reference/feature/PTG_Multipro_Human_Immune_Intracellular_Antibodies_PROT002.csv".into()),
            chemistry_specs: specs_input_gex(AutoChemistryName::Count),
            ..Default::default()
        };
            let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
            assert_eq!(outs.chemistry_defs[&LibraryType::Gex].name, MFRP_RNA_R1);
            assert_eq!(
                outs.chemistry_defs[&LibraryType::Antibody].name,
                MFRP_Ab_R2pos50
            );
        }

        #[test]
        fn test_ab_only_sfrp() {
            let args = DetectChemistryStageInputs{
            sample_def: vec![SampleDef {
                fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                read_path: testdata_path("fastqs/cellranger/multi/1247954_1270332_frp_ab_fastqs/ab"),
                sample_names: Some(vec!["ab".into()]),
                target_set: Some(refdata_path("targeted_panels/human_wta_RTL.Fletcher_v7.0.all_manufacturing.mismatch_base_pairing_revision.csv")),
                lanes: Some(vec![1]),
                library_type: Some(LibraryType::Gex),
                ..Default::default()
            }],
            multi_config: Some("../dui_tests/test_resources/cellranger-multi/1247954_sfrp_ab_only.csv".into()),
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
                    sample_def: vec![SampleDef {
                        fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                        read_path: "../dui_tests/test_resources/cellranger-count/few_aligned_reads_v2/"
                            .into(),
                        sample_names: Some(vec!["test_sample".into()]),
                        ..Default::default()
                    }],
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
                    SampleDef {
                        fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                        read_path:
                            "../dui_tests/test_resources/cellranger-count/CR-120-01_SC3pv2_15k"
                                .into(),
                        sample_names: Some(vec!["CR-120-01".into()]),
                        ..Default::default()
                    },
                    SampleDef {
                        fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                        read_path:
                            "../dui_tests/test_resources/cellranger-count/CR-CUSTOM-22_SC5pr2_15k"
                                .into(),
                        sample_names: Some(vec!["CR-CUSTOM-22".into()]),
                        ..Default::default()
                    },
                ],
                chemistry_specs: specs_input_gex(AutoChemistryName::Count),
                reference_path: Some(refdata_path("GRCh38")),
                ..Default::default()
            };
            err_snapshot!(args)
        }

        #[test]
        fn test_mixed_chemistries_3pv2_5p_pe() {
            let args = DetectChemistryStageInputs {
                sample_def: vec![
                    SampleDef {
                        fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                        read_path:
                            "../dui_tests/test_resources/cellranger-count/CR-120-01_SC3pv2_15k"
                                .into(),
                        sample_names: Some(vec!["CR-120-01".into()]),
                        ..Default::default()
                    },
                    SampleDef {
                        fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                        read_path:
                            "../dui_tests/test_resources/cellranger-count/CR-CUSTOM-17_SC5p-PE_15k"
                                .into(),
                        sample_names: Some(vec!["CR-CUSTOM-17".into()]),
                        ..Default::default()
                    },
                ],
                chemistry_specs: specs_input_gex(AutoChemistryName::Count),
                reference_path: Some(refdata_path("GRCh38")),
                ..Default::default()
            };
            err_snapshot!(args)
        }

        #[test]
        fn test_5pr2v3() {
            let args = DetectChemistryStageInputs {
                sample_def: vec![SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path: testdata_path("fastqs/cellranger/multi/1523899_gex_5pr2v3"),
                    sample_names: Some(vec!["tiny_gex".into()]),
                    lanes: Some(vec![4]),
                    library_type: Some(LibraryType::Gex),
                    ..Default::default()
                }],
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
                sample_def: vec![SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path: testdata_path("fastqs/cellranger/multi/1523899_gex_5pr2v3"),
                    sample_names: Some(vec!["tiny_gex".into()]),
                    lanes: Some(vec![4]),
                    library_type: Some(LibraryType::Gex),
                    ..Default::default()
                }],
                multi_config: Some(
                    "../dui_tests/test_resources/cellranger-multi/1523899_gex_5pr2ohv3.csv".into(),
                ),
                chemistry_specs: specs_input_gex(AutoChemistryName::Count),
                ..Default::default()
            };
            let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
            assert_eq!(
                outs.chemistry_defs.primary().name,
                ChemistryName::FivePrimeR2OHV3
            );
        }

        #[test]
        fn test_3pv4() {
            let args = DetectChemistryStageInputs {
                sample_def: vec![SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path: testdata_path("fastqs/cellranger/multi/1561255_gex/tiny_gex"),
                    sample_names: Some(vec!["gex".into()]),
                    lanes: Some(vec![1]),
                    library_type: Some(LibraryType::Gex),
                    ..Default::default()
                }],
                chemistry_specs: specs_input_gex(AutoChemistryName::Count),
                ..Default::default()
            };
            let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
            assert_eq!(
                outs.chemistry_defs.primary().name,
                ChemistryName::ThreePrimeV4
            );
        }

        #[test]
        fn test_crispr_and_gex_3pv4_fails() {
            let args = DetectChemistryStageInputs {
                sample_def: vec![
                    SampleDef {
                        fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                        read_path: testdata_path("fastqs/cellranger/multi/detectchem_unittest_tiny_gex_crispr_3pv4/crispr"),
                        sample_names: Some(vec!["crispr".into()]),
                        library_type: Some(LibraryType::Crispr),
                        ..Default::default()
                    },
                    SampleDef {
                        fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                        read_path: testdata_path("fastqs/cellranger/multi/detectchem_unittest_tiny_gex_crispr_3pv4/gex"),
                        sample_names: Some(vec!["gex".into()]),
                        library_type: Some(LibraryType::Gex),
                        ..Default::default()
                    },
                ],
                chemistry_specs: specs_input_gex(AutoChemistryName::Count),
                feature_reference: Some(
                    "../dui_tests/test_resources/reference/feature/RAB1A_NT_Aug2020.csv".into(),
                ),
                ..Default::default()
            };
            err_snapshot!(args)
        }
    }
}
