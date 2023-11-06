//! Martian stage DETECT_CHEMISTRY

use crate::barcode_overlap::FRPGemBarcodeOverlapRow;
use crate::detect_chemistry::chemistry_filter::{
    detect_chemistry_units, ChemistryFilter, DetectChemistryUnit,
};
use crate::detect_chemistry::errors::DetectChemistryErrors;
use crate::detect_chemistry::identity_check::check_fastq_identity;
use crate::detect_chemistry::length_filter::LengthFilter;
use crate::detect_chemistry::mapping_filter::ReadMappingFilter;
use crate::detect_chemistry::probe_bc_check::check_probe_bc;
use crate::detect_chemistry::probe_bc_pairing::{
    detect_probe_barcode_pairing, should_detect_probe_barcode_pairing,
};
use crate::detect_chemistry::whitelist_filter::WhitelistMatchFilter;
use anyhow::{ensure, Context, Result};
use barcode::whitelist::BarcodeId;
use cr_types::chemistry::{
    normalize_chemistry_def, AutoChemistryName, AutoOrRefinedChemistry, ChemistryDef, ChemistryName,
};
use cr_types::reference::feature_reference::{FeatureConfig, FeatureReferenceFile};
use cr_types::rna_read::LegacyLibraryType;
use cr_types::sample_def::SampleDef;
use fastq_set::read_pair::ReadPair;
use itertools::Itertools;
use lazy_static::lazy_static;
use martian::prelude::*;
use martian_derive::{make_mro, MartianStruct};
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::tabular_file::CsvFile;
use martian_filetypes::FileTypeWrite;
use metric::{set, TxHashMap, TxHashSet};
use multi::config::{
    multiconst, MultiConfigCsv, MultiConfigCsvFile, ProbeBarcodeIterationMode, SamplesCsv,
};
use parameters_toml::{fiveprime_multiplexing, threeprime_lt_multiplexing};
use serde::{Deserialize, Serialize};
use std::borrow::Cow;
use std::iter::zip;
use std::path::PathBuf;
use AutoChemistryName::{Count, FivePrime, ThreePrime, Vdj};
use AutoOrRefinedChemistry::{Auto, Refined};
#[allow(clippy::enum_glob_use)]
use ChemistryName::*;

const MIN_READS_NEEDED: usize = 10_000;

lazy_static! {
    static ref MIX_SC5P_R2_PE: TxHashSet<ChemistryName> = {
        let mut set = TxHashSet::default();
        set.insert(FivePrimeR2);
        set.insert(FivePrimePE);
        set
    };
    static ref MIX_SC5P_R2: TxHashSet<ChemistryName> = {
        let mut set = TxHashSet::default();
        set.insert(FivePrimeR2);
        set
    };
}

#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct DetectChemistryStageInputs {
    pub sample_def: Vec<SampleDef>,
    pub reference_path: Option<PathBuf>,
    pub feature_reference: Option<FeatureReferenceFile>,
    pub chemistry_name_spec: AutoOrRefinedChemistry,
    pub allowed_chems: Option<Vec<AutoOrRefinedChemistry>>,
    pub r1_length: Option<usize>,
    pub r2_length: Option<usize>,
    pub multi_config: Option<MultiConfigCsvFile>,
    pub is_pd: bool,
    pub custom_chemistry_def: Option<ChemistryDef>,
    pub feature_config: Option<FeatureConfig>,
}

pub type DetectedProbeBarcodePairingFile = JsonFile<TxHashMap<BarcodeId, BarcodeId>>;

#[derive(Clone, Serialize, Deserialize, MartianStruct)]
#[cfg_attr(test, derive(Debug))]
pub struct DetectChemistryStageOutputs {
    pub chemistry_def: ChemistryDef,
    pub is_antibody_only: bool,
    #[mro_retain]
    pub probe_barcode_overlap: Option<CsvFile<FRPGemBarcodeOverlapRow>>,
    #[mro_retain]
    pub detected_probe_barcode_pairing: Option<DetectedProbeBarcodePairingFile>,
}

/// Return true if probe barcode IDs are uncollapsed, like BC001A. For PD use.
fn is_rtl_uncollapsed(multi_config_csv: &MultiConfigCsv) -> bool {
    multi_config_csv
        .sample_barcode_ids_used_in_experiment(ProbeBarcodeIterationMode::All)
        .into_iter()
        .all(|x| x.ends_with(&['A', 'B', 'C', 'D']))
}

/// Validate combinations of chemistry and types of multiplexing.
fn validate_multiplexing(chemistry_type: ChemistryName, sample_defs: &[SampleDef]) -> Result<()> {
    if !sample_defs
        .iter()
        .any(|sdef| sdef.library_type == Some(LegacyLibraryType::Multiplexing))
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

/// Validate a manually-specified chemistry.
///
/// If the chemistry is explicitly specified as a refined chemistry, we will check that it passes
/// the whitelist filter. Manual chemistry is used as an escape hatch by customers for data which is
/// usually of lower quality or by QA to run non-standard fuzzed FASTQs. No minimum number of reads
/// is enforced here. Emit a warning if there are few valid barcodes.
fn validate_manual_chemistry(
    chem: ChemistryName,
    sample_defs: &[SampleDef],
    multi_config_csv: Option<&MultiConfigCsv>,
    feature_reference: Option<&FeatureReferenceFile>,
    feature_config: Option<&FeatureConfig>,
    units: &[DetectChemistryUnit],
    read_pairs_all_units: &[Vec<ReadPair>],
) -> Result<()> {
    let chems = &set![chem];

    if chem.is_rtl_multiplexed() {
        // Check the read length to ensure that the probe barcode is sequenced.
        let _matching_chemistries =
            LengthFilter::new(chems, chems, feature_reference, feature_config)?
                .filter_chemistries(units, read_pairs_all_units)?;
    }

    if let Err(err) =
        WhitelistMatchFilter::new(chems, chems)?.filter_chemistries(units, read_pairs_all_units)
    {
        if chem.is_spatial() {
            return Err(anyhow::Error::from(err));
        } else {
            println!("WARNING: {err:#}");
        }
    }

    validate_multiplexing(chem, sample_defs)?;
    validate_rtl(multi_config_csv, chem)?;
    Ok(())
}

// This is our stage struct
pub struct DetectChemistry;

#[make_mro(mem_gb = 20, volatile = strict)]
impl MartianMain for DetectChemistry {
    type StageInputs = DetectChemistryStageInputs;
    type StageOutputs = DetectChemistryStageOutputs;

    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        let multi_config_csv = if let Some(csv) = &args.multi_config {
            Some(csv.read().with_context(|| csv.display().to_string())?)
        } else {
            None
        };
        let multi_config_csv = multi_config_csv.as_ref();

        let (chemistry_def, reads) = Self::main_inner(&args, multi_config_csv)?;

        let is_antibody_only = args
            .sample_def
            .iter()
            .all(|sd| sd.library_type == Some(LegacyLibraryType::AntibodyCapture));

        let mut outputs = DetectChemistryStageOutputs {
            chemistry_def,
            is_antibody_only,
            probe_barcode_overlap: None,
            detected_probe_barcode_pairing: None,
        };

        if let Some(multi_config_csv) = multi_config_csv {
            handle_probe_barcode_translation(&args, &rover, multi_config_csv, &mut outputs, reads)?;
        }
        Ok(outputs)
    }
}

impl DetectChemistry {
    /// All logic that determines which chemistry to use.
    fn main_inner(
        args: &DetectChemistryStageInputs,
        multi_config_csv: Option<&MultiConfigCsv>,
    ) -> Result<(ChemistryDef, Option<ReadUnits>)> {
        let allowed_chems = args.allowed_chems.as_deref();
        let sample_defs = &args.sample_def;

        // -----------------------------------------------------------------------------------------
        // Bail out if the input chemistry_name_spec is not in the optional list of allowed
        // chemistry names
        validate_chemistry_spec(args.chemistry_name_spec, allowed_chems)?;

        // -----------------------------------------------------------------------------------------
        // The simplest case is when the input chemistry is custom. **Nothing** is checked if the
        // chemistry is "custom". This is used only for PD purposes.
        if let Refined(ChemistryName::Custom) = args.chemistry_name_spec {
            let custom_def = args.custom_chemistry_def.clone().expect(
                "Custom chemistry def should be present if the input chemistry is 'Custom'",
            );
            assert_eq!(custom_def.name, ChemistryName::Custom);
            return Ok((
                normalize_chemistry_def(custom_def.clone()).unwrap_or(custom_def),
                None,
            ));
        }

        // Check for overhang multiplexing
        let is_overhang_multiplexed =
            multi_config_csv.map_or(false, MultiConfigCsv::is_overhang_multiplexed);

        // -----------------------------------------------------------------------------------------
        // If the sample def contains VDJ, none of the other library types should be present and
        // the chemistry def should be `SCVDJ_auto`.
        if args
            .sample_def
            .iter()
            .any(|sd| sd.library_type == Some(LegacyLibraryType::Vdj))
        {
            ensure!(
                args.sample_def
                    .iter()
                    .all(|sd| sd.library_type == Some(LegacyLibraryType::Vdj)),
                "VDJ cannot be mixed with other library types"
            );
            ensure!(
                args.chemistry_name_spec == Auto(Vdj),
                "VDJ chemistry_name_spec should be SCVDJ_auto"
            );
        }

        let (ref units, read_pairs_all_units) = load_read_units(args)?;

        // Define groups of refined chemistries
        let mut default_chems = TxHashSet::default();
        let three_prime_chems: TxHashSet<ChemistryName> = set![
            ThreePrimeV1,
            ThreePrimeV2,
            ThreePrimeV3,
            ThreePrimeV3LT,
            ThreePrimeV3HT
        ];
        let five_prime_chems: TxHashSet<ChemistryName> = set![FivePrimeR2, FivePrimePE];
        default_chems.extend(three_prime_chems.clone());
        default_chems.extend(five_prime_chems.clone());

        // Singleplex and multiplex FRP chemistries are indistinguishable if the probe BC is sequenced
        // Singleplex FRP is valid with either count or a multi config without a [samples] section,
        // while multiplex FRP is only valid with a multi config that has a [samples] section.
        let samples = multi_config_csv.and_then(|csv| csv.samples.as_ref());
        default_chems.extend(if samples.is_none() {
            [SFRP].iter()
        } else if args.is_pd {
            if multi_config_csv.map_or(false, is_rtl_uncollapsed) {
                [MFRP_uncollapsed, MFRP_R1_48_uncollapsed].iter()
            } else {
                [MFRP, MFRP_47, MFRP_R1].iter()
            }
        } else {
            [MFRP, MFRP_R1].iter()
        });

        let possible_chemistries = match args.chemistry_name_spec {
            Refined(chem) => {
                validate_manual_chemistry(
                    chem,
                    sample_defs,
                    multi_config_csv,
                    args.feature_reference.as_ref(),
                    args.feature_config.as_ref(),
                    units,
                    &read_pairs_all_units,
                )?;
                return Ok((ChemistryDef::named(chem), Some(read_pairs_all_units)));
            }
            Auto(Count) => default_chems,
            Auto(ThreePrime) => three_prime_chems,
            Auto(FivePrime) => five_prime_chems,
            Auto(Vdj) => set![VdjPE, VdjR2],
        };

        let allowed_chemistries = match allowed_chems {
            Some(chems) => Cow::Owned(
                chems
                    .iter()
                    .filter_map(|x| {
                        if let Refined(chem) = x {
                            Some(*chem)
                        } else {
                            None
                        }
                    })
                    .collect(),
            ),
            None => Cow::Borrowed(&possible_chemistries),
        };

        for (unit, read_pairs) in zip(units, &read_pairs_all_units) {
            ensure!(
                read_pairs.len() >= MIN_READS_NEEDED,
                DetectChemistryErrors::NotEnoughReads {
                    num_reads: read_pairs.len(),
                    min_reads: MIN_READS_NEEDED,
                    unit: Box::new(unit.clone()),
                }
            );
        }

        println!("Potential chemistries: {possible_chemistries:?}");

        // Read length based filtering
        let length_matching_chemistries = LengthFilter::new(
            &allowed_chemistries,
            &possible_chemistries,
            args.feature_reference.as_ref(),
            args.feature_config.as_ref(),
        )?
        .filter_chemistries(units, &read_pairs_all_units)?;

        println!("After length filter: {length_matching_chemistries:?}");

        // -----------------------------------------------------------------------------------------
        // For each unit of fastqs, ensure that a sufficient fraction of reads contain barcodes
        // which match the whitelist for at least one of the possible chemistries
        let wl_matching_chemistries =
            WhitelistMatchFilter::new(&allowed_chemistries, &length_matching_chemistries)?
                .filter_chemistries(units, &read_pairs_all_units)?;

        let chosen_chemistry_def = if wl_matching_chemistries.len() == 1 {
            let chem = *wl_matching_chemistries.iter().next().unwrap();
            if chem == SFRP {
                // Bail out if a mixture of probe barcodes is observed for singleplex FRP chemistry
                check_probe_bc(
                    units,
                    &read_pairs_all_units,
                    args.feature_reference.as_ref(),
                    args.feature_config.as_ref(),
                )?;
            }
            Some(validate_chemistry(
                chem,
                allowed_chems,
                sample_defs,
                multi_config_csv,
                is_overhang_multiplexed,
            )?)
        } else if wl_matching_chemistries == set![ThreePrimeV3, ThreePrimeV3HT, ThreePrimeV3LT] {
            Some(validate_chemistry(
                ThreePrimeV3LT,
                allowed_chems,
                sample_defs,
                multi_config_csv,
                is_overhang_multiplexed,
            )?)
        } else if wl_matching_chemistries == set![ThreePrimeV3, ThreePrimeV3HT] {
            Some(validate_chemistry(
                ThreePrimeV3,
                allowed_chems,
                sample_defs,
                multi_config_csv,
                is_overhang_multiplexed,
            )?)
        } else if wl_matching_chemistries == set![FivePrimeR2, FivePrimeHT] {
            Some(validate_chemistry(
                FivePrimeR2,
                allowed_chems,
                sample_defs,
                multi_config_csv,
                is_overhang_multiplexed,
            )?)
        } else {
            None
        };
        if let Some(chosen_chemistry_def) = chosen_chemistry_def {
            return Ok((chosen_chemistry_def, Some(read_pairs_all_units)));
        }

        println!("After whitelist filter: {wl_matching_chemistries:?}");

        let is_antibody_only = sample_defs
            .iter()
            .all(|sd| sd.library_type == Some(LegacyLibraryType::AntibodyCapture));
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
                args.reference_path.as_deref().unwrap(),
                &allowed_chemistries,
                wl_matching_chemistries,
            )?;
            mapper.filter_chemistries(units, &read_pairs_all_units)?
        };

        println!("After mapping filter: {expected_mapping_chemistries:?}");

        ensure!(
            expected_mapping_chemistries.len() == 1,
            "Could not distinguish between {expected_mapping_chemistries:?}"
        );

        Ok((
            validate_chemistry(
                expected_mapping_chemistries.into_iter().next().unwrap(),
                allowed_chems,
                sample_defs,
                multi_config_csv,
                is_overhang_multiplexed,
            )?,
            Some(read_pairs_all_units),
        ))
    }
}

/// A nested collection of reads from a detect chemistry unit, per unit.
type ReadUnits = Vec<Vec<ReadPair>>;

/// Load a subsample of reads from each input unit.
fn load_read_units(
    args: &DetectChemistryStageInputs,
) -> Result<(Vec<DetectChemistryUnit>, ReadUnits)> {
    let units = detect_chemistry_units(&args.sample_def, args.r1_length, args.r2_length)?;
    println!("Number of fastq units = {}", units.len());

    // check for duplicate R1 and R2 files amongst units
    check_fastq_identity(&units)?;

    let read_pairs_all_units = units
        .iter()
        .map(|unit| {
            println!("Sampling reads from: {unit}");
            unit.sampled_read_pairs()
        })
        .try_collect()?;
    Ok((units, read_pairs_all_units))
}

/// Identify probe barcode pairings if necessary.
/// If pairings are present, modify the translation whitelist to include those
/// pairings.
fn handle_probe_barcode_translation(
    args: &DetectChemistryStageInputs,
    rover: &MartianRover,
    multi_config_csv: &MultiConfigCsv,
    outputs: &mut DetectChemistryStageOutputs,
    reads: Option<ReadUnits>,
) -> Result<()> {
    let detected_pairing_id_translation = if !should_detect_probe_barcode_pairing(multi_config_csv)
    {
        None
    } else {
        let reads = match reads {
            Some(reads) => reads,
            None => {
                // Load reads if we haven't already loaded them.
                let (_, reads) = load_read_units(args)?;
                reads
            }
        };
        let (overlaps, pairings) = detect_probe_barcode_pairing(&outputs.chemistry_def, &reads)?;

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

    // Remap the probe barcode whitelist if we have a probe barcode pairing.
    // Use the explicit pairing if one is provided.
    // Otherwise, use the detected pairing if there are any.
    if let Some(pairing) = explicit_pairing.or(detected_pairing_id_translation) {
        let path: PathBuf = rover.make_path("probe_barcode_translation_whitelist.tsv");
        outputs
            .chemistry_def
            .translate_probe_barcode_whitelist_with_id_map(&pairing, path)?;
    };
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use cr_types::chemistry::known_chemistry_defs;
    use cr_types::sample_def::FastqMode;
    use dui_tests::stage_test::StageFailTest;
    use dui_tests::{stage_fail_dui_test, DuiTest};
    use multi::config::{
        FeatureType, GemWell, GeneExpressionParams, Lanes, LibrariesCsv, Library, SampleRow,
    };
    use std::io::Write;
    use std::path::Path;
    use std::string::ToString;

    #[test]
    fn test_validate_rtl() {
        use ChemistryName::{ThreePrimeV3, MFRP, SFRP};

        assert!(validate_rtl(None, ThreePrimeV3).is_ok());

        let gex_library = Library::BclProcessor {
            fastq_path: PathBuf::default(),
            sample_indices: Vec::default(),
            lanes: Lanes::Any,
            physical_library_id: String::default(),
            feature_types: vec![FeatureType::GeneExpression],
            gem_well: GemWell(0),
            subsample_rate: None,
        };

        let antibody_library = Library::BclProcessor {
            fastq_path: PathBuf::default(),
            sample_indices: Vec::default(),
            lanes: Lanes::Any,
            physical_library_id: String::default(),
            feature_types: vec![FeatureType::AntibodyCapture],
            gem_well: GemWell(0),
            subsample_rate: None,
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
        assert!(validate_rtl(Some(&gex_without_probe_set), MFRP).is_err());
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
        assert!(validate_rtl(Some(&gex_with_probe_set), MFRP).is_ok());
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
        assert!(validate_rtl(Some(&antibody_only), MFRP).is_ok());
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
        assert!(validate_rtl(Some(&include_introns_false), MFRP).is_err());
        assert!(validate_rtl(Some(&include_introns_false), SFRP).is_err());
        assert!(validate_rtl(Some(&include_introns_false), ThreePrimeV3).is_ok()); // to aggr with RTL

        let samples = MultiConfigCsv {
            samples: Some(SamplesCsv(vec![SampleRow::from_probe_barcode_id("BC001")])),
            ..gex_with_probe_set
        };
        assert!(validate_rtl(Some(&samples), MFRP).is_ok());
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

    impl Default for DetectChemistryStageInputs {
        fn default() -> Self {
            DetectChemistryStageInputs {
                sample_def: vec![],
                reference_path: None,
                chemistry_name_spec: AutoChemistryName::Count.into(),
                allowed_chems: None,
                r1_length: None,
                r2_length: None,
                feature_reference: None,
                multi_config: None,
                is_pd: false,
                custom_chemistry_def: None,
                feature_config: None,
            }
        }
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
            chemistry_name_spec: AutoChemistryName::Count.into(),
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
        ))
    }

    #[test]
    fn test_duplicate_r1_r2() {
        err_snapshot!(test_invalid_bcl2fastq_input(
            "../dui_tests/test_resources/cellranger-count/fastqs_duplicate_r1_r2",
            None,
            vec!["test_sample"],
        ))
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
            chemistry_name_spec: ChemistryName::VdjPE.into(),
            allowed_chems: Some(vec![ChemistryName::ThreePrimeV3.into()]),
            ..Default::default()
        };
        err_snapshot!(args)
    }

    #[test]
    fn test_custom_chemistry_spec() {
        let mut def: ChemistryDef = ChemistryDef::named(ChemistryName::ThreePrimeV3);
        def.name = ChemistryName::Custom;
        def.rna.offset = 13;
        let args = DetectChemistryStageInputs {
            chemistry_name_spec: ChemistryName::Custom.into(),
            custom_chemistry_def: Some(def.clone()),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(outs.chemistry_def, def);
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
                chemistry_name_spec: ChemistryName::Custom.into(),
                custom_chemistry_def: Some(custom_def),
                ..Default::default()
            };
            let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
            // Ensure that the output is mapped to a known chemistry
            assert_eq!(outs.chemistry_def, known_def);
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
                ..SampleDef::default()
            }],
            chemistry_name_spec: ChemistryName::ThreePrimeV2.into(),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(outs.chemistry_def.name, ChemistryName::ThreePrimeV2);
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
            chemistry_name_spec: ChemistryName::ThreePrimeV3.into(),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(outs.chemistry_def.name, ChemistryName::ThreePrimeV3);
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
            chemistry_name_spec: AutoChemistryName::Count.into(),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(outs.chemistry_def.name, ChemistryName::ThreePrimeV3);
        assert!(!outs.is_antibody_only);
    }

    #[test]
    fn test_correct_chemistry_spec_lt() {
        let args = DetectChemistryStageInputs {
            sample_def: vec![SampleDef {
                fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                read_path: "../dui_tests/test_resources/cellranger-count/LT_chemistry".into(),
                sample_names: Some(vec!["1077080".into()]),
                ..Default::default()
            }],
            chemistry_name_spec: ChemistryName::ThreePrimeV3LT.into(),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(outs.chemistry_def.name, ChemistryName::ThreePrimeV3LT);
        assert!(!outs.is_antibody_only);
    }

    #[test]
    fn test_auto_chemistry_spec_3pv3lt() {
        let args = DetectChemistryStageInputs {
            sample_def: vec![SampleDef {
                fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                read_path: "../dui_tests/test_resources/cellranger-count/LT_chemistry".into(),
                sample_names: Some(vec!["1077080".into()]),
                ..Default::default()
            }],
            chemistry_name_spec: AutoChemistryName::Count.into(),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(outs.chemistry_def.name, ChemistryName::ThreePrimeV3LT);
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
            chemistry_name_spec: ChemistryName::ThreePrimeV3HT.into(),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(outs.chemistry_def.name, ChemistryName::ThreePrimeV3HT);
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
            chemistry_name_spec: AutoChemistryName::Count.into(),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(outs.chemistry_def.name, ChemistryName::ThreePrimeV3);
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
                chemistry_name_spec: AutoChemistryName::Count.into(),
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
            chemistry_name_spec: AutoChemistryName::Count.into(),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(outs.chemistry_def.name, ChemistryName::ThreePrimeV3);
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
            chemistry_name_spec: AutoChemistryName::Vdj.into(),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(outs.chemistry_def.name, ChemistryName::VdjPE);
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
            chemistry_name_spec: AutoChemistryName::Vdj.into(),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(outs.chemistry_def.name, ChemistryName::VdjR2);
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
                library_type: Some(LegacyLibraryType::AntibodyCapture),
                ..Default::default()
            }],
            chemistry_name_spec: AutoChemistryName::Count.into(),
            feature_reference: Some("../dui_tests/test_resources/reference/feature/Biolegend_x22_TotalSeqB_14Ab_Isotype_Controls_20180601.csv".into()),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(outs.chemistry_def.name, ChemistryName::ThreePrimeV3);
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
                library_type: Some(LegacyLibraryType::AntibodyCapture),
                ..Default::default()
            }],
            chemistry_name_spec: AutoChemistryName::Count.into(),
            feature_reference: Some("../dui_tests/test_resources/reference/feature/Biolegend_x22_TotalSeqB_14Ab_Isotype_Controls_len24_20210721.csv".into()),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(outs.chemistry_def.name, ChemistryName::ThreePrimeV3);
        assert!(outs.is_antibody_only);
    }

    #[test]
    fn test_ab_only_3pv3lt() {
        let args = DetectChemistryStageInputs {
            sample_def: vec![SampleDef {
                fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                read_path: "../dui_tests/test_resources/cellranger-count/LT_chemistry_ab".into(),
                sample_names: Some(vec!["1077130".into()]),
                lanes: Some(vec![1]),
                library_type: Some(LegacyLibraryType::AntibodyCapture),
                ..Default::default()
            }],
            chemistry_name_spec: AutoChemistryName::Count.into(),
            feature_reference: Some("../dui_tests/test_resources/reference/feature/Biolegend_x22_TotalSeqB_14Ab_Isotype_Controls_20180601.csv".into()),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(outs.chemistry_def.name, ChemistryName::ThreePrimeV3LT);
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
                library_type: Some(LegacyLibraryType::AntibodyCapture),
                ..Default::default()
            }],
            chemistry_name_spec: AutoChemistryName::Count.into(),
            feature_reference: Some("../dui_tests/test_resources/reference/feature/Biolegend_x22_TotalSeqB_14Ab_Isotype_Controls_20180601.csv".into()),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(outs.chemistry_def.name, ChemistryName::FeatureBarcodingOnly);
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
                    library_type: Some(LegacyLibraryType::AntibodyCapture),
                    ..Default::default()
                },
                SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path:
                        "../dui_tests/test_resources/cellranger-count/pbmc_1k_protein_v3_gex".into(),
                    sample_names: Some(vec!["pbmc_1k_protein_v3_gex".into()]),
                    lanes: Some(vec![4]),
                    library_type: Some(LegacyLibraryType::GeneExpression),
                    ..Default::default()
                },
            ],
            chemistry_name_spec: AutoChemistryName::Count.into(),
            feature_reference: Some("../dui_tests/test_resources/reference/feature/Biolegend_x22_TotalSeqB_14Ab_Isotype_Controls_20180601.csv".into()),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(outs.chemistry_def.name, ChemistryName::ThreePrimeV3);
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
                    library_type: Some(LegacyLibraryType::AntibodyCapture),
                    ..Default::default()
                },
                SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path:
                        "../dui_tests/test_resources/cellranger-count/pbmc_1k_protein_v3_gex".into(),
                    sample_names: Some(vec!["pbmc_1k_protein_v3_gex".into()]),
                    lanes: Some(vec![4]),
                    library_type: Some(LegacyLibraryType::GeneExpression),
                    ..Default::default()
                },
            ],
            chemistry_name_spec: AutoChemistryName::Count.into(),
            feature_reference: Some("../dui_tests/test_resources/reference/feature/Biolegend_x22_TotalSeqB_14Ab_Isotype_Controls_len24_20210721.csv".into()),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(outs.chemistry_def.name, ChemistryName::ThreePrimeV3);
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
                    library_type: Some(LegacyLibraryType::AntibodyCapture),
                    ..Default::default()
                },
                SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path:
                        "../dui_tests/test_resources/cellranger-count/pbmc_1k_protein_v3_gex".into(),
                    sample_names: Some(vec!["pbmc_1k_protein_v3_gex".into()]),
                    lanes: Some(vec![4]),
                    library_type: Some(LegacyLibraryType::GeneExpression),
                    ..Default::default()
                },
            ],
            chemistry_name_spec: AutoChemistryName::Count.into(),
            feature_reference: Some(
                "../dui_tests/test_resources/reference/feature/too_long_features.csv".into(),
            ),
            ..Default::default()
        };
        err_snapshot!(args)
    }

    #[test]
    fn test_ab_and_gex_3pv3lt() {
        let args = DetectChemistryStageInputs {
            sample_def: vec![
                SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path: "../dui_tests/test_resources/cellranger-count/LT_chemistry_gex_ab"
                        .into(),
                    sample_names: Some(vec!["1077130_ab".into()]),
                    lanes: Some(vec![1]),
                    library_type: Some(LegacyLibraryType::AntibodyCapture),
                    ..Default::default()
                },
                SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path: "../dui_tests/test_resources/cellranger-count/LT_chemistry_gex_ab"
                        .into(),
                    sample_names: Some(vec!["1077130_gex".into()]),
                    lanes: Some(vec![1]),
                    library_type: Some(LegacyLibraryType::GeneExpression),
                    ..Default::default()
                },
            ],
            chemistry_name_spec: AutoChemistryName::Count.into(),
            feature_reference: Some("../dui_tests/test_resources/reference/feature/Biolegend_x22_TotalSeqB_14Ab_Isotype_Controls_20180601.csv".into()),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(outs.chemistry_def.name, ChemistryName::ThreePrimeV3LT);
        assert!(!outs.is_antibody_only);
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
                    library_type: Some(LegacyLibraryType::CrisprGuideCapture),
                    ..Default::default()
                },
                SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path: "../dui_tests/test_resources/cellranger-count/K562_5k_crispr_v3_gex"
                        .into(),
                    sample_names: Some(vec!["K562_5k_crispr_v3_gex".into()]),
                    library_type: Some(LegacyLibraryType::GeneExpression),
                    ..Default::default()
                },
            ],
            chemistry_name_spec: AutoChemistryName::Count.into(),
            feature_reference: Some(
                "../dui_tests/test_resources/reference/feature/RAB1A_NT_Aug2020.csv".into(),
            ),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(outs.chemistry_def.name, ChemistryName::ThreePrimeV3);
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
                    library_type: Some(LegacyLibraryType::CrisprGuideCapture),
                    ..Default::default()
                },
                SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path: "../dui_tests/test_resources/cellranger-count/K562_5k_crispr_v3_gex"
                        .into(),
                    sample_names: Some(vec!["K562_5k_crispr_v3_gex".into()]),
                    library_type: Some(LegacyLibraryType::GeneExpression),
                    ..Default::default()
                },
            ],
            chemistry_name_spec: AutoChemistryName::Count.into(),
            feature_reference: Some(
                "../dui_tests/test_resources/reference/feature/too_long_features.csv".into(),
            ),
            ..Default::default()
        };
        err_snapshot!(args)
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
            chemistry_name_spec: AutoChemistryName::Count.into(),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(outs.chemistry_def.name, ChemistryName::ThreePrimeV1);
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
            chemistry_name_spec: AutoChemistryName::Count.into(),
            ..Default::default()
        };
        err_snapshot!(args)
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
                chemistry_name_spec: AutoChemistryName::Count.into(),
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
                chemistry_name_spec: AutoChemistryName::Count.into(),
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
            chemistry_name_spec: AutoChemistryName::Vdj.into(),
            r1_length: Some(26),
            r2_length: Some(91),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(outs.chemistry_def.name, ChemistryName::VdjR2);
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
            chemistry_name_spec: AutoChemistryName::Count.into(),
            r1_length: Some(25),
            ..Default::default()
        };
        err_snapshot!(args)
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
                chemistry_name_spec: AutoChemistryName::Count.into(),
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
                chemistry_name_spec: AutoChemistryName::Count.into(),
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
            chemistry_name_spec: AutoChemistryName::Count.into(),
            r1_length: Some(26),
            r2_length: Some(75),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(outs.chemistry_def.name, ChemistryName::ThreePrimeV3);
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
            chemistry_name_spec: ChemistryName::ThreePrimeV2.into(),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(outs.chemistry_def.name, ChemistryName::ThreePrimeV2);
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
            chemistry_name_spec: AutoChemistryName::Count.into(),
            ..Default::default()
        };
        err_snapshot!(args)
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
                library_type: Some(LegacyLibraryType::GeneExpression),
                target_set: Some("../dui_tests/test_resources/cellranger-count/RTL_single_sample_probe_bcs_mixture/probe_set.csv".into()),
                ..Default::default()
            }],
            chemistry_name_spec: ChemistryName::SFRP.into(),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(outs.chemistry_def.name, ChemistryName::SFRP);
    }

    #[test]
    fn test_rtl_no_probe_set() {
        let multi_config = b"\
            [gene-expression]\n\
            ref,/ref\n\
            [libraries]\n\
            fastq_id,fastqs,feature_types\n\
            rtl,/fastqs,Gene Expression\n";
        let multi_config_csv = tempfile::Builder::new().suffix(".csv").tempfile().unwrap();
        multi_config_csv.as_file().write_all(multi_config).unwrap();

        let args = DetectChemistryStageInputs {
            chemistry_name_spec: ChemistryName::SFRP.into(),
            multi_config: Some(MultiConfigCsvFile::from(multi_config_csv.path())),
            ..Default::default()
        };
        err_snapshot!(args)
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
            chemistry_name_spec: AutoChemistryName::Count.into(),
            multi_config: Some(
                "../dui_tests/test_resources/cellranger-multi/overhang_3pv3.csv".into(),
            ),
            ..Default::default()
        };
        let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
        assert_eq!(outs.chemistry_def.name, ChemistryName::ThreePrimeV3OH);
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
                chemistry_name_spec: AutoChemistryName::Count.into(),
                ..Default::default()
            };
            let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
            assert_eq!(outs.chemistry_def.name, ChemistryName::ThreePrimeV2);
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
                chemistry_name_spec: AutoChemistryName::Count.into(),
                ..Default::default()
            };
            let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
            assert_eq!(outs.chemistry_def.name, ChemistryName::FivePrimeR2);
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
                chemistry_name_spec: AutoChemistryName::FivePrime.into(),
                ..Default::default()
            };
            let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
            assert_eq!(outs.chemistry_def.name, ChemistryName::FivePrimeR2);
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
                chemistry_name_spec: AutoChemistryName::Count.into(),
                ..Default::default()
            };
            let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
            assert_eq!(outs.chemistry_def.name, ChemistryName::FivePrimePE);
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
                chemistry_name_spec: AutoChemistryName::Count.into(),
                ..Default::default()
            };
            err_snapshot!(args)
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
                chemistry_name_spec: AutoChemistryName::Count.into(),
                ..Default::default()
            };
            err_snapshot!(args)
        }

        #[test]
        fn test_ab_and_gex_5p() {
            let args = DetectChemistryStageInputs {
            sample_def: vec![SampleDef {
                fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                read_path: sere_path("vdj_nextgem_hs_pbmc3_5gex_protein/vdj_nextgem_hs_pbmc3_5gex_protein_antibody"),
                sample_names: Some(vec!["vdj_nextgem_hs_pbmc3_5gex_protein_antibody".into()]),
                lanes: Some(vec![1]),
                library_type: Some(LegacyLibraryType::AntibodyCapture),
                ..Default::default()
            }, SampleDef {
                fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                read_path: sere_path("vdj_nextgem_hs_pbmc3_5gex_protein/vdj_nextgem_hs_pbmc3_5gex_protein_gex"),
                sample_names: Some(vec!["vdj_nextgem_hs_pbmc3_5gex_protein_gex".into()]),
                lanes: Some(vec![1]),
                library_type: Some(LegacyLibraryType::GeneExpression),
                ..Default::default()
            }],
            chemistry_name_spec: AutoChemistryName::Count.into(),
            reference_path: Some(refdata_path("GRCh38")),
            feature_reference: Some("../dui_tests/test_resources/reference/feature/Biolegend_x22_TotalSeqB_14Ab_Isotype_Controls_20180601.csv".into()),
            ..Default::default()
        };
            let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
            assert_eq!(outs.chemistry_def.name, ChemistryName::FivePrimeR2);
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
                    library_type: Some(LegacyLibraryType::GeneExpression),
                    ..Default::default()
                },
                SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    read_path: testdata_path("fastqs/cellranger/multi/1247954_1270332_frp_gex_and_ab_fastqs"),
                    sample_names: Some(vec!["1247954_AB".into()]),
                    lanes: Some(vec![4]),
                    library_type: Some(LegacyLibraryType::AntibodyCapture),
                    ..Default::default()
                }
            ],
            multi_config: Some("../dui_tests/test_resources/cellranger-multi/1247954_1270332_rtl_sfrp_gex_and_ab.csv".into()),
            chemistry_name_spec: AutoChemistryName::Count.into(),
            feature_reference: Some(testdata_path("cellranger/multi/feature_refs/BioLegend_x22_TotalSeqB_Human_PBMC_20210114.csv").into()),
            ..Default::default()
        };
            let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
            assert_eq!(outs.chemistry_def.name, ChemistryName::SFRP);
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
                library_type: Some(LegacyLibraryType::GeneExpression),
                ..Default::default()
            }],
            multi_config: Some("../dui_tests/test_resources/cellranger-multi/1245033_sfrp.csv".into()),
            chemistry_name_spec: AutoChemistryName::Count.into(),
            ..Default::default()
        };
            let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
            assert_eq!(outs.chemistry_def.name, ChemistryName::SFRP);
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
                library_type: Some(LegacyLibraryType::GeneExpression),
                ..Default::default()
            }],
            multi_config: Some("../dui_tests/test_resources/cellranger-multi/1216790_1245489.csv".into()),
            chemistry_name_spec: AutoChemistryName::Count.into(),
            ..Default::default()
        };
            let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
            assert_eq!(outs.chemistry_def.name, ChemistryName::MFRP);
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
                library_type: Some(LegacyLibraryType::GeneExpression),
                ..Default::default()
            }],
            multi_config: Some("../dui_tests/test_resources/cellranger-multi/1320608_mfrp_r1_50_90.csv".into()),
            chemistry_name_spec: AutoChemistryName::Count.into(),
            ..Default::default()
            };
            let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
            assert_eq!(outs.chemistry_def.name, ChemistryName::MFRP_R1);
        }

        #[test]
        fn test_ab_only_mfrp() {
            let args = DetectChemistryStageInputs{
            sample_def: vec![SampleDef {
                fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                read_path: testdata_path("fastqs/cellranger/multi/1402853_mfrp_ab/ab"),
                sample_names: Some(vec!["ab".into()]),
                target_set: Some(refdata_path("targeted_panels/human_wta_RTL.Fletcher_v7.0.all_manufacturing.mismatch_base_pairing_revision.csv")),
                lanes: Some(vec![1]),
                library_type: Some(LegacyLibraryType::GeneExpression),
                ..Default::default()
            }],
            multi_config: Some("../dui_tests/test_resources/cellranger-multi/1402853_mfrp_ab_only.csv".into()),
            chemistry_name_spec: AutoChemistryName::Count.into(),
            ..Default::default()
            };
            let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
            assert_eq!(outs.chemistry_def.name, ChemistryName::MFRP);
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
                library_type: Some(LegacyLibraryType::GeneExpression),
                ..Default::default()
            }],
            multi_config: Some("../dui_tests/test_resources/cellranger-multi/1247954_sfrp_ab_only.csv".into()),
            chemistry_name_spec: AutoChemistryName::Count.into(),
            ..Default::default()
            };
            let outs = DetectChemistry.test_run_tmpdir(args).unwrap();
            assert_eq!(outs.chemistry_def.name, ChemistryName::SFRP);
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
                    chemistry_name_spec: AutoChemistryName::Count.into(),
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
                chemistry_name_spec: AutoChemistryName::Count.into(),
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
                chemistry_name_spec: AutoChemistryName::Count.into(),
                reference_path: Some(refdata_path("GRCh38")),
                ..Default::default()
            };
            err_snapshot!(args)
        }
    }
}
