#![expect(missing_docs)]
use super::{
    AntigenSpecificityRow, FunctionalMapRow, Library, MultiConfigCsv, ProbeBarcodeIterationMode,
    SampleRow, create_feature_config, multiconst,
};
use crate::cmo_set::load_default_cmo_set;
use crate::config::{ChemistrySet, PROBE_BARCODE_ID_GROUPING, get_default_overhang_set, libsconst};
use anyhow::{Context, Result, bail, ensure};
use barcode::whitelist::{
    BarcodeId, RTLMultiplexingBarcodeType, categorize_rtl_multiplexing_barcode_id,
};
use cr_types::chemistry::{AutoOrRefinedChemistry, ChemistryDef, ChemistryName};
use cr_types::reference::feature_extraction::FeatureExtractor;
use cr_types::reference::feature_reference::{BeamMode, FeatureConfig, FeatureReference};
use cr_types::{ERROR_CODE_INFO, FeatureBarcodeType, LibraryType, VdjChainType};
use fastq_set::WhichRead;
use itertools::Itertools;
use metric::{TxHashMap, TxHashSet};
use std::collections::HashSet;
use std::fmt::{Display, Formatter};
use std::fs::File;
use std::io::BufReader;
use std::path::Path;
use std::sync::Arc;
use strsim::levenshtein;
use transcriptome::Transcriptome;
use vdj_reference::VdjReference;

#[derive(Debug, Clone, Copy)]
pub struct SectionCtx {
    pub section: &'static str,
    pub field: &'static str,
}

impl Display for SectionCtx {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "[{}] {}", self.section, self.field)
    }
}

pub fn check_duplicate_libraries(libs: &[Library]) -> Result<()> {
    if libs.is_empty() {
        bail!(
            "[{}] is empty, no libraries provided! Please specify 1 or more libraries.",
            multiconst::LIBRARIES
        );
    }
    for i in 0..(libs.len() - 1) {
        for j in (i + 1)..libs.len() {
            if libs[i].overlaps(&libs[j]) {
                bail!(
                    "[{}] has overlapping library definitions between rows {} and {}",
                    multiconst::LIBRARIES,
                    i + 1,
                    j + 1
                );
            }
        }
    }
    Ok(())
}

fn find_approximate_matches<'a>(
    input: &'a str,
    set: &'a TxHashSet<String>,
    max_distance: usize,
) -> impl Iterator<Item = &'a str> {
    set.iter()
        .map(String::as_str)
        .filter(move |item| levenshtein(input, item) <= max_distance)
}

struct InvalidIdsResult<'a> {
    invalid_ids: Vec<&'a str>,
    possible_typos: Vec<(&'a str, Vec<&'a str>)>,
}

fn get_invalid_id_entries<'a>(
    valid_options: &'a TxHashSet<String>,
    entries: impl IntoIterator<Item = &'a str>,
    search_for_typos: bool,
) -> InvalidIdsResult<'a> {
    let invalid_ids: Vec<&'a str> = entries
        .into_iter()
        .filter(|x| !valid_options.contains(*x))
        .collect();
    let mut possible_typos: Vec<_> = vec![];
    if search_for_typos {
        possible_typos = invalid_ids
            .iter()
            .filter_map(|word| {
                let approx: Vec<_> = find_approximate_matches(word, valid_options, 4).collect();
                (!approx.is_empty()).then_some((*word, approx))
            })
            .collect();
    }

    InvalidIdsResult {
        invalid_ids,
        possible_typos,
    }
}

fn test_reserved_words(words: &TxHashSet<String>) -> Vec<String> {
    words
        .iter()
        .filter(|word| {
            matches!(
                word.to_ascii_lowercase().as_str(),
                "blank" | "unassigned" | "multiplet"
            )
        })
        .cloned()
        .collect()
}

pub const MULTI_HELP: &str = r#"No input FASTQs were found for the requested parameters.

Please ensure the following for each entry in the [libraries] section of your multi CSV:
 - Make sure you are specifying the correct `fastq_id`, i.e. matching the sample sheet
 - Make sure your files follow the correct naming convention, e.g. SampleName_S1_L001_R1_001.fastq.gz (and the R2 version)
 - Make sure your `fastqs` points to the correct location.
 - Make sure your `lanes`, if any, are correctly specified.

Refer to the "Specifying Input FASTQs" page at https://support.10xgenomics.com/ for more details.

"#;

/// Ensure that the library types are compatible with RTL.
pub fn check_libraries_rtl(cfg: &MultiConfigCsv) -> Result<()> {
    let (supported, unsupported): (Vec<_>, Vec<_>) = cfg
        .libraries
        .library_types()
        .sorted()
        .dedup()
        .partition(|&library_type| {
            matches!(
                library_type,
                LibraryType::Gex | LibraryType::Antibody | LibraryType::Crispr
            )
        });
    ensure!(
        unsupported.is_empty(),
        "These library types are not compatible with Flex: {}",
        unsupported.iter().format(", ")
    );
    if supported.contains(&LibraryType::Crispr) {
        ensure!(
            supported.contains(&LibraryType::Gex),
            "{} libraries must be paired with a {} library",
            LibraryType::Crispr,
            LibraryType::Gex
        );
    }
    Ok(())
}

pub fn check_libraries(
    cfg: &MultiConfigCsv,
    fref: Option<Arc<FeatureReference>>,
    is_pd: bool,
    hostname: &str,
) -> Result<()> {
    use Library::{Bcl2Fastq, BclProcessor};

    for lib in &cfg.libraries.0 {
        let fastq_path = match &lib {
            Bcl2Fastq { fastqs, .. } => fastqs,
            BclProcessor { fastq_path, .. } => fastq_path,
        };
        if !fastq_path.is_absolute() {
            bail!(
                "Specified FASTQ folder must be an absolute path: {}",
                fastq_path.display(),
            );
        }
        if !fastq_path.exists() {
            bail!(
                "On machine: {}, specified FASTQ folder does not exist: {}",
                hostname,
                fastq_path.display(),
            );
        }
        if !fastq_path.is_dir() {
            bail!(
                "Specified FASTQ directory is not a directory: {}",
                fastq_path.display(),
            );
        }

        ensure!(
            std::fs::read_dir(fastq_path)
                .with_context(|| {
                    format!(
                        "On machine: {}, permissions forbid opening specified FASTQ directory: {}",
                        hostname,
                        fastq_path.display()
                    )
                })?
                .next()
                .is_some(),
            "Specified FASTQ directory is empty: {}",
            fastq_path.display(),
        );

        // Make sure out reference actually contains this feature barcode type.
        if let Some(feature_type) = lib.library_type().feature_barcode_type()
            && let Some(fref) = &fref
        {
            ensure!(
                fref.feature_maps.contains_key(&feature_type),
                "TXRNGR10014: You declared a library with feature_type = '{feature_type}', but \
                     there are no features with that feature_type in the feature reference. \
                     {ERROR_CODE_INFO}",
                feature_type = feature_type.as_str(),
            );
        }

        // traditionally, no error is thrown if feature reference is not provided

        match lib.to_sample_def() {
            Err(_) => bail!("{MULTI_HELP}"),
            Ok(sdef) => sdef.check_fastqs(MULTI_HELP)?,
        }
        // check fastq_path exists, is a folder, is not empty, and sample indices are valid

        ensure!(
            lib.chemistry() != Some(AutoOrRefinedChemistry::Custom) || is_pd,
            "Unknown chemistry {} in library {}.",
            ChemistryName::Custom,
            lib.physical_library_id(),
        );
    }

    // Detect if we are in Beam-T or Beam-Ab mode
    let mut beam_mode: Option<BeamMode> = None;
    if cfg.libraries.has_antigen_capture() && cfg.libraries.has_vdj_t_or_gd() {
        beam_mode = Some(BeamMode::BeamT);
    } else if cfg.libraries.has_antigen_capture() && cfg.libraries.has_vdj_b() {
        beam_mode = Some(BeamMode::BeamAB);
    }

    if let Some(beam_mode) = beam_mode {
        fref.as_ref()
            .unwrap()
            .validate_beam_feature_ref(beam_mode)?;
    }

    // Ensure that we don't only have multiplexing libraries.
    let library_types: Vec<_> = cfg.libraries.library_types().unique().collect();
    ensure!(
        library_types.as_slice()
            != [LibraryType::FeatureBarcodes(
                FeatureBarcodeType::Multiplexing
            )],
        "Only Multiplexing Capture feature libraries were provided. \
         Please provide additional non-Multiplexing Capture libraries."
    );

    if cfg.is_rtl() && !is_pd {
        check_libraries_rtl(cfg)?;
    }
    if cfg.libraries.has_multiplexing() && cfg.libraries.has_vdj() {
        ensure!(
            is_pd,
            "The combination of VDJ libraries and Multiplexing Capture libraries is not supported."
        );
    }

    Ok(())
}

pub fn check_samples(
    cfg: &MultiConfigCsv,
    fref: Option<Arc<FeatureReference>>,
    is_pd: bool,
) -> Result<()> {
    let Some(samples) = &cfg.samples else {
        return Ok(());
    };
    // Validate sample names against reserved words
    let sample_names: TxHashSet<_> = samples
        .0
        .iter()
        .map(|sample| sample.sample_id.clone())
        .collect();
    let invalid_sample_names = test_reserved_words(&sample_names);
    ensure!(
        invalid_sample_names.is_empty(),
        "Invalid sample_ids ('{}') provided, please ensure you are not using the reserved words \
        'blank', 'multiplet', or 'unassigned'.",
        invalid_sample_names.join("', '")
    );

    // Validate CMO IDs
    let multiplexing_ids = fref
        .as_ref()
        .map_or_else(TxHashSet::default, |x| x.multiplexing_ids());

    let invalid_multiplexing_ids = test_reserved_words(&multiplexing_ids);
    ensure!(
        invalid_multiplexing_ids.is_empty(),
        "Invalid cmo_ids ('{}') provided, please ensure you are not using the reserved words \
        'blank', 'multiplet', or 'unassigned'.",
        invalid_multiplexing_ids.join("', '")
    );

    for sample in &samples.0 {
        let invalid_result = get_invalid_id_entries(
            &multiplexing_ids,
            sample.cmo_ids.iter().flatten().map(String::as_str),
            true,
        );
        if !invalid_result.invalid_ids.is_empty() {
            let mut err_msg = format!(
                "Unknown cmo_ids ('{}') provided for sample '{}', please ensure you are either \
                 using valid 10x CMO IDs or are providing the correct [gene-expression] cmo-set.",
                invalid_result.invalid_ids.join("', '"),
                sample.sample_id,
            );

            err_msg = format!(
                "{err_msg}\n\nValid IDs are currently:\n{}",
                multiplexing_ids.iter().sorted().join("\n")
            );

            if !invalid_result.possible_typos.is_empty() {
                err_msg = format!(
                    "{err_msg}\n\n\
                     Did you perhaps mean any of the following? (Input -> Intended Value)?\n{}",
                    invalid_result
                        .possible_typos
                        .iter()
                        .map(|z| format!("{} -> {}\n", z.0, z.1.join(" or ")))
                        .join("\n")
                );
            }
            bail!(err_msg);
        }
    }

    // Validate HASHTAG IDs
    if samples.has_hashtag_ids() {
        let hashtag_ids: TxHashSet<String> = samples
            .0
            .iter()
            .filter_map(|s| s.sample_barcode_ids(ProbeBarcodeIterationMode::All))
            .flatten()
            .map(ToString::to_string)
            .collect();
        let invalid_hashtag_ids = test_reserved_words(&hashtag_ids);
        ensure!(
            invalid_hashtag_ids.is_empty(),
            "Invalid hashtag_ids ('{}') provided, please ensure you are not using the reserved words \
            'blank', 'multiplet', or 'unassigned'.",
            invalid_hashtag_ids.join("', '")
        );
        let antibody_ids = fref.map_or_else(TxHashSet::default, |x| x.antibody_ids());
        for sample in &samples.0 {
            let invalid_result = get_invalid_id_entries(
                &antibody_ids,
                sample.hashtag_ids.iter().flatten().map(String::as_str),
                false,
            );
            ensure!(
                invalid_result.invalid_ids.is_empty(),
                "Unknown hashtag_ids ('{}') provided for sample '{}', please ensure that you are using \
                IDs that correspond to a valid Antibody Capture feature.",
                invalid_result.invalid_ids.join("', '"),
                sample.sample_id
            );
        }
    }

    // Validate probe barcode IDs
    if samples.has_probe_barcode_ids() && !is_pd {
        let probe_barcode_whitelist_ids: Vec<BarcodeId> = cfg
            .chemistry_specs()?
            .into_iter()
            .flat_map(|(library_type, spec)| {
                match spec {
                    AutoOrRefinedChemistry::Refined(chem) => vec![chem],
                    AutoOrRefinedChemistry::Auto(_) => {
                        vec![
                            ChemistrySet::Mfrp
                                .chemistry_for_library_type(library_type)
                                // Unclear if this should never fail at this point in validation,
                                // so just use a default in case we have an invalid library type combination.
                                .unwrap_or(ChemistryName::MFRP_RNA),
                            ChemistryName::Flex_v2_R1,
                        ]
                    }
                }
            })
            .map(|chem| {
                ChemistryDef::named(chem)
                    .barcode_whitelist_source()?
                    .probe()
                    .get_ids()
            })
            .flatten_ok()
            .try_collect()?;

        for sample in &samples.0 {
            {
                let valid_options: TxHashSet<_> = probe_barcode_whitelist_ids
                    .iter()
                    .map(ToString::to_string)
                    .collect();
                let invalid_result = get_invalid_id_entries(
                    &valid_options,
                    sample
                        .sample_barcode_ids(ProbeBarcodeIterationMode::All)
                        .unwrap(),
                    false,
                );
                ensure!(
                    invalid_result.invalid_ids.is_empty(),
                    "Unknown probe_barcode_ids ('{}') provided for sample '{}', please ensure you \
                     are using IDs from the following list of valid 10x probe barcode IDs: {}",
                    invalid_result.invalid_ids.join("', '"),
                    sample.sample_id,
                    probe_barcode_whitelist_ids
                        .into_iter()
                        .sorted()
                        .dedup()
                        .format(", "),
                );
            }

            for grouping in sample
                .sample_barcode_id_groupings(ProbeBarcodeIterationMode::All)
                .unwrap()
            {
                check_probe_barcode_id_grouping(
                    &grouping,
                    cfg.libraries.has_gene_expression(),
                    cfg.libraries.has_antibody_capture(),
                    cfg.libraries.has_crispr_guide_capture(),
                )
                .with_context(|| {
                    format!(
                        "Invalid probe_barcode_ids grouping ('{}') provided for sample '{}'",
                        grouping.join(PROBE_BARCODE_ID_GROUPING),
                        sample.sample_id
                    )
                })?;
            }
        }
    }

    // Validate OH IDs
    if samples.has_overhang_ids() {
        let default_overhang_ids: TxHashSet<String> = get_default_overhang_set()
            .iter()
            .map(ToString::to_string)
            .collect();

        for sample in &samples.0 {
            let invalid_result = get_invalid_id_entries(
                &default_overhang_ids,
                sample
                    .sample_barcode_ids(ProbeBarcodeIterationMode::All)
                    .unwrap(),
                true,
            );
            ensure!(
                // allow invalid overhang IDs for PD-samples only
                invalid_result.invalid_ids.is_empty() || is_pd,
                "Unknown overhang_ids ('{}') provided for sample '{}'.",
                invalid_result.invalid_ids.iter().join("', '"),
                sample.sample_id
            );
        }
    }
    Ok(())
}

/// Validate that a probe barcode ID grouping is valid.
fn check_probe_barcode_id_grouping(
    barcode_ids: &[&str],
    has_gex_lib: bool,
    has_ab_lib: bool,
    has_cr_lib: bool,
) -> Result<()> {
    use RTLMultiplexingBarcodeType::{Antibody, Crispr, Gene};

    let barcode_types: Vec<_> = barcode_ids
        .iter()
        .map(|x| categorize_rtl_multiplexing_barcode_id(x))
        .collect();
    if barcode_types == [None] {
        // This barcode ID is compatible with all library types.
        return Ok(());
    }

    let barcode_types: Vec<_> = barcode_types.into_iter().flatten().collect();
    match (has_gex_lib, has_ab_lib, has_cr_lib) {
        (false, false, false) => Ok(()),
        (true, false, false) => match_one(&barcode_types, Gene),
        (false, true, false) => match_one(&barcode_types, Antibody),
        (true, true, false) => match_two(&barcode_types, Antibody),
        (true, false, true) => match_two(&barcode_types, Crispr),
        (true, true, true) => match barcode_types.as_slice() {
            // Can match one, two, or three.
            // Order should always be RTL+AB+CR
            [] => bail!("no barcode ID provided"),
            [Gene] | [Gene, Antibody] | [Gene, Crispr] | [Gene, Antibody, Crispr] => Ok(()),
            [other_type] => bail!("expected a {Gene} barcode ID, not {other_type}"),
            anything_else => bail!(
                "expected either a single {Gene} barcode ID, \
                 or a pair of {Gene} + {Antibody} or {Gene} + {Crispr} barcode IDs, \
                 or a triple of {Gene} + {Antibody} + {Crispr} barcode IDs, not {}",
                anything_else.iter().join(" + "),
            ),
        },
        // We should be validating this case somewhere upstack.
        (false, _, true) => panic!("CRISPR library present without GEX"),
    }
}

/// Match a single barcode of the provided type.
fn match_one(
    barcode_types: &[RTLMultiplexingBarcodeType],
    expected: RTLMultiplexingBarcodeType,
) -> Result<()> {
    match barcode_types {
        [] => bail!("no barcode ID provided"),
        [one] if *one == expected => Ok(()),
        [other_type] => bail!("expected {expected} barcode ID, not {other_type}"),
        too_many => bail!(
            "expected a single {expected} barcode ID, not {}",
            too_many.iter().join(" + "),
        ),
    }
}

/// Match one or two barcodes; one must be RTL.
fn match_two(
    barcode_types: &[RTLMultiplexingBarcodeType],
    expected: RTLMultiplexingBarcodeType,
) -> Result<()> {
    use RTLMultiplexingBarcodeType::Gene;
    match barcode_types {
        [] => bail!("no barcode ID provided"),
        [Gene] => Ok(()),
        [other_type] => bail!("expected a {Gene} barcode ID, not {other_type}"),
        [Gene, other] if *other == expected => Ok(()),

        [other, Gene] if *other == expected => bail!(
            "when pairing {Gene} and {expected} barcode IDs, \
             provide the {Gene} barcode ID first"
        ),
        anything_else => bail!(
            "expected either a {Gene} barcode ID or a pair of \
             {Gene} + {expected} barcode IDs, not {}",
            anything_else.iter().join(" + "),
        ),
    }
}

/// Validate the [samples] section, check for re-used cmo_ids or probe_barcode_ids.
pub fn check_duplicate_sample_barcode_ids(samples: &[SampleRow]) -> Result<()> {
    // Check for duplicate sample barcode IDs.
    // doing this here instead of when building the MultiConfigCsv,
    // during which we try to construct a CrMultiGraph and enforce
    // any of its rules
    let mut used: TxHashMap<&str, &SampleRow> = TxHashMap::default();
    for sample in samples {
        if let Some(sample_barcode_ids) = sample.sample_barcode_ids(ProbeBarcodeIterationMode::All)
        {
            for sample_barcode_id in sample_barcode_ids {
                if let Some(orig) = used.insert(sample_barcode_id, sample) {
                    bail!(
                        "Re-used {} ('{}') provided for sample '{}' (already provided for sample '{}').",
                        sample.sample_barcode_ids_column_name(),
                        sample_barcode_id,
                        sample.sample_id,
                        orig.sample_id
                    );
                }
            }
        }
    }
    Ok(())
}

/// Validate the [samples] section. Check for duplicate sample_ids
pub fn check_duplicate_samples(samples: &[SampleRow]) -> Result<()> {
    check_unique(samples.iter().map(|s| &s.sample_id), "samples", "sample_id")?;
    Ok(())
}

fn check_unique<'a, T: 'static + Eq + Display>(
    vals: impl Iterator<Item = &'a T>,
    section: &str,
    field: &str,
) -> Result<()> {
    let mut seen_values = Vec::new();
    for v in vals {
        ensure!(
            !seen_values.contains(&v),
            "In the [{section}] section, there are two different entries with {field} = '{v}'. \
             Please give each entry a unique value in the {field} column. If a sample was tagged \
             with multiple CMOs, these values should be entered as pipe-separated entries in the \
             config file (e.g. CMO1|CMO2).",
        );
        seen_values.push(v);
    }
    Ok(())
}

pub fn check_gem_wells(libs: &[Library]) -> Result<()> {
    assert!(!libs.is_empty());
    let gem_wells: Vec<_> = libs
        .iter()
        .map(|x| x.gem_well().0)
        .sorted()
        .dedup()
        .collect();
    ensure!(
        gem_wells[0] == 1,
        "[{}] gem_well numbering must start at 1",
        multiconst::LIBRARIES
    );

    for (i, gw) in gem_wells.into_iter().enumerate() {
        let expected = i as u16 + 1;
        ensure!(
            gw == expected,
            "[{}] gem_wells must be numbered contiguously starting from 1, missing gem_well: {}",
            multiconst::LIBRARIES,
            expected
        );
    }
    Ok(())
}

pub fn check_physical_library_ids(libs: &[Library]) -> Result<()> {
    use multiconst::LIBRARIES;

    let mut mapping = TxHashMap::default();
    let mut num_vdj_auto_libs = 0;
    for (i, lib) in libs.iter().enumerate() {
        if lib.library_type() == LibraryType::VdjAuto {
            num_vdj_auto_libs += 1;
            continue;
        }
        if let Some((j, pli)) = mapping.insert(
            (lib.library_type(), lib.gem_well()),
            (i, lib.physical_library_id()),
        ) {
            ensure!(
                pli == lib.physical_library_id(),
                "[{LIBRARIES}] invalid physical_library_id on row {row0}: '{library0}'. \
                 This feature_type '{library_type}' has already been given a physical_library_id: \
                 '{library1}' on row {row1}",
                library_type = lib.library_type(),
                library0 = lib.physical_library_id(),
                library1 = pli,
                row0 = i + 2,
                row1 = j + 2
            );
        }
    }
    ensure!(
        num_vdj_auto_libs <= 3,
        "[{LIBRARIES}] invalid physical_library_id(s): found {num_vdj_auto_libs} VDJ libraries, \
         but we expect at most 3 VDJ libraries",
    );
    Ok(())
}

const GEX_REF_FILES: &[&str; 2] = &["reference.json", "fasta/genome.fa"];
const GEX_STAR_FILES: &[&str; 8] = &[
    "star/chrLength.txt",
    "star/chrNameLength.txt",
    "star/chrName.txt",
    "star/chrStart.txt",
    "star/Genome",
    "star/genomeParameters.txt",
    "star/SA",
    "star/SAindex",
];

pub fn check_gex_reference<D: Display>(
    ref_ctx: &D,
    ref_path: &Path,
    hostname: &str,
) -> Result<Transcriptome> {
    for file in GEX_REF_FILES {
        let path = &ref_path.join(file);
        if check_file(ref_ctx, path).is_err() {
            bail!(
                "Your {ref_ctx} does not contain the expected files, or they are not readable. Please check your reference folder on {hostname}."
            );
        }
    }

    if check_file(ref_ctx, &ref_path.join("genes/genes.gtf")).is_err()
        && check_file(ref_ctx, &ref_path.join("genes/genes.gtf.gz")).is_err()
    {
        bail!(
            "Your {} is missing gene annotations that should be present at {}/genes/genes.gtf[.gz], or they are not readable. Please check your reference folder on {}.",
            ref_ctx,
            ref_path.display(),
            hostname,
        );
    }

    for file in GEX_STAR_FILES {
        let path = &ref_path.join(file);
        if check_file(ref_ctx, path).is_err() {
            bail!("Your {ref_ctx} does not appear to be indexed. Please run `mkref`");
        }
    }

    Transcriptome::from_reference_path(ref_path)
}

pub fn check_feature_reference<D: Display>(
    ref_ctx: &D,
    ref_path: &Path,
    hostname: &str,
    feature_config: Option<&FeatureConfig>,
) -> Result<Arc<FeatureReference>> {
    if check_file(ref_ctx, ref_path).is_err() {
        bail!(
            "Your {} is either missing or not readable from {}: {}",
            ref_ctx,
            hostname,
            ref_path.display(),
        );
    }
    let rdr = BufReader::new(File::open(ref_path)?);
    let fref = FeatureReference::new(&[], None, Some(rdr), None, None, None, feature_config)?;
    let fref = Arc::new(fref);
    let _ = FeatureExtractor::new(fref.clone(), None)?;
    Ok(fref)
}

/// Return a feature reference complete with CMOs, if needed for the analysis.
/// The boolean returned in second position is true if the feature reference
/// is using 10X internal CMOs.  If the analysis isn't using CMO multiplexing,
/// returns None in second position.
pub fn build_feature_reference_with_cmos(
    cfg: &MultiConfigCsv,
    is_pd: bool,
    hostname: &str,
) -> Result<(Option<Arc<FeatureReference>>, Option<bool>)> {
    let feature_reference = cfg
        .feature
        .as_ref()
        .and_then(|f| f.reference_path.as_ref())
        .map(|feature_reference_path| {
            check_feature_reference(
                &SectionCtx {
                    section: "feature",
                    field: "reference",
                },
                feature_reference_path,
                hostname,
                create_feature_config(
                    cfg.antigen_specificity.as_ref(),
                    cfg.functional_map.as_ref(),
                    cfg.libraries.beam_mode(),
                    cfg.samples.as_ref(),
                )
                .as_ref(),
            )
        })
        .transpose()?;

    let cmo_set = cfg
        .gene_expression
        .as_ref()
        .and_then(|gex| gex.cmo_set.as_ref())
        .map(|cmo_set| {
            check_cmo_set(
                &SectionCtx {
                    section: "gene-expression",
                    field: "cmo-set",
                },
                cmo_set,
                hostname,
            )
        })
        .transpose()?;

    if !cfg.libraries.has_multiplexing() && cmo_set.is_none() {
        // Not CMO multiplexed, return the provided feature reference.
        return Ok((feature_reference, None));
    }

    let (feature_reference, tenx_cmos) =
        validate_cmo_set(feature_reference, cmo_set)?.build_feature_reference(is_pd)?;
    Ok((Some(feature_reference), Some(tenx_cmos)))
}

fn check_cmo_set<D: Display>(
    ref_ctx: &D,
    ref_path: &Path,
    hostname: &str,
) -> Result<Arc<FeatureReference>> {
    let fref = check_feature_reference(ref_ctx, ref_path, hostname, None)?;
    ensure!(
        fref.feature_defs
            .iter()
            .all(|x| x.feature_type.is_multiplexing()),
        "All CMO set definitions must be of feature_type 'Multiplexing Capture'"
    );
    Ok(fref)
}

pub fn check_vdj_reference<D: Display>(
    ref_ctx: &D,
    ref_path: &Path,
    hostname: &str,
) -> Result<VdjReference> {
    let path = ref_path.join("fasta/regions.fa");
    if check_file(ref_ctx, &path).is_err() {
        bail!(
            "Your {} does not contain the expected file ({}), or they are not readable. Please check your reference folder on {}.",
            ref_ctx,
            path.display(),
            hostname,
        );
    }
    VdjReference::check(ref_path).with_context(|| {
        format!(
            "Failed while verifying the contents of the V(D)J reference ({}) at \"{}\"",
            ref_ctx,
            ref_path.display(),
        )
    })?;
    let vdj_ref = VdjReference::from_reference_fasta(&path)?;
    Ok(vdj_ref)
}

pub fn check_library_combinations(libraries: &[Library]) -> Result<()> {
    let mut has_vdj = false;
    let mut has_antigen = false;
    let mut has_gex = false;
    let mut has_crispr = false;
    let mut uniq_vdj_features: HashSet<VdjChainType> = HashSet::new();
    for lib in libraries {
        has_vdj |= lib.is_vdj();
        has_antigen |= lib.is_antigen();
        has_gex |= lib.is_gex();
        has_crispr |= lib.is_crispr();
        if let Some(vdj_chain_type) = lib.library_type().vdj_chain_type() {
            uniq_vdj_features.insert(vdj_chain_type);
        }
    }
    if has_antigen && !(has_gex && has_vdj) {
        bail!(
            "[{}] Antigen Capture library requires paired VDJ and Gene Expression libraries.",
            multiconst::LIBRARIES
        );
    }
    if has_antigen && uniq_vdj_features.len() > 1 {
        bail!(
            "[{}] Antigen Capture library cannot be accompanied by more than one VDJ feature type.",
            multiconst::LIBRARIES
        );
    }
    if has_antigen && has_crispr {
        bail!(
            "[{}] The combination of Antigen Capture and CRISPR libraries is not supported.",
            multiconst::LIBRARIES
        );
    }
    Ok(())
}

/// Enforce valid combinations of library chemistries.
/// Currently only Flex libraries can use this feature, so make sure that anything
/// specified here is a Flex chemistry.
pub fn check_library_chemistries(libraries: &[Library]) -> Result<()> {
    use libsconst::CHEMISTRY;
    use multiconst::LIBRARIES;

    let chems: Vec<_> = libraries.iter().filter_map(Library::chemistry).collect();

    // All None is valid.
    if chems.is_empty() {
        return Ok(());
    }

    // Auto chemistry is not supported in the libraries section.
    let auto_chems: Vec<_> = chems
        .iter()
        .filter_map(AutoOrRefinedChemistry::auto)
        .collect();
    ensure!(
        auto_chems.is_empty(),
        "[{LIBRARIES}] Specifying auto chemistries at the library level is not supported: ({}).",
        auto_chems.iter().format(", ")
    );

    // Ensure that if they provided one, they provided all.
    ensure!(
        chems.len() == libraries.len(),
        "[{LIBRARIES}] A chemistry name must be provided for all libraries or none."
    );

    // Ensure that everything is a flex chemistry.
    let non_flex_chems: Vec<_> = chems
        .iter()
        .filter(|&&chem| chem.is_rtl() != Some(true) && chem != AutoOrRefinedChemistry::Custom)
        .collect();
    ensure!(
        non_flex_chems.is_empty(),
        "[{LIBRARIES}] Only Flex assays may specify chemistry at the per-library level; \
         invalid chemistries: {}",
        non_flex_chems.iter().unique().format(", "),
    );

    // Ensure that we have the same chemistry if multiple libraries are of the same type.
    for (lib_type, libs) in libraries
        .iter()
        .into_group_map_by(|lib| lib.library_type())
        .into_iter()
        .sorted_by_key(|(lib_type, _libs)| *lib_type)
    {
        let unique_chems: Vec<_> = libs
            .into_iter()
            .filter_map(Library::chemistry)
            .unique()
            .collect();
        ensure!(
            unique_chems.len() <= 1,
            "[{LIBRARIES}] Conflicting {CHEMISTRY} for {lib_type} libraries ({chems}); \
             manual chemistry must be the same for all libraries of the same type.",
            chems = unique_chems.iter().format(", ")
        );

        if let Some(chem) = unique_chems.into_iter().at_most_one().unwrap() {
            // Already validated that this is a explicit chem.
            let chem = chem.refined().unwrap();
            ensure!(
                chem.compatible_with_library_type(lib_type),
                "[{LIBRARIES}] The {CHEMISTRY} {chem} is not valid for the feature type {lib_type}",
            );
        }
    }
    Ok(())
}

pub fn check_antigen_specificity(controls: &[AntigenSpecificityRow]) -> Result<()> {
    let mut control_ids = TxHashSet::default();
    for control in controls {
        // If the set did have this value present, insert returns false
        if !control_ids.insert(&control.control_id) {
            bail!("Duplicate {} 'control_id'!", control.control_id);
        }
    }

    let alleles: TxHashSet<_> = controls
        .iter()
        .filter_map(|x| x.mhc_allele.as_deref())
        .collect();
    if alleles.is_empty() {
        ensure!(
            controls.len() < 2,
            "More than one value for 'control_id' specified without 'mhc_allele' information. \
             For BCR Antigen Capture experiments only a single control is valid, \
             for TCR Antigen Capture both 'control_id' and 'mhc_allele' are required."
        );
    } else {
        ensure!(
            alleles.len() == controls.len(),
            "For TCR Antigen Capture experiments values in both 'control_id' and 'mhc_allele' \
             are required and only a single 'control_id' can be specified for an 'mhc_allele'."
        );
    }
    Ok(())
}

pub fn check_feature_functional_map(func_map_rows: &[FunctionalMapRow]) -> Result<()> {
    let mut seen_functional_names = TxHashSet::default();
    let mut seen_feature_ids = TxHashSet::default();
    for row in func_map_rows {
        if !seen_functional_names.insert(&row.functional_name) {
            bail!(
                "Duplicate {} 'functional_name' in [feature-functional-map] section!",
                row.functional_name
            );
        }
        for feature_id in &row.feature_ids {
            if !seen_feature_ids.insert(feature_id) {
                bail!("Duplicate {feature_id} 'feature_id' in [feature-functional-map] section!");
            }
        }
    }
    Ok(())
}

pub fn check_file<D: Display>(ctx: D, path: &Path) -> Result<()> {
    if path.is_file() && File::open(path).is_ok() {
        return Ok(());
    }
    bail!(
        "{} is not readable, please check file permissions: {}",
        ctx,
        path.display()
    )
}

fn check_tenx_cmos(feature_ref: &FeatureReference) -> bool {
    let tenx_cmos = load_default_cmo_set().unwrap();
    let tenx_cmos: TxHashSet<_> = tenx_cmos
        .iter()
        .map(|x| (x.read, &x.pattern, &x.sequence))
        .collect();
    for fdef in &feature_ref.feature_defs {
        if !fdef.feature_type.is_multiplexing() {
            continue;
        }
        if !tenx_cmos.contains(&(fdef.read, &fdef.pattern, &fdef.sequence)) {
            return false;
        }
    }
    true
}

#[derive(Debug, Clone)]
pub enum CmoSetAction {
    Builtin,
    MergeBuiltin(Arc<FeatureReference>),
    InCmoSet(Option<Arc<FeatureReference>>, Arc<FeatureReference>),
    InFeatureRef(Arc<FeatureReference>),
}

const MAX_MULTIPLEXING_TAGS: usize = 12;

impl CmoSetAction {
    /// Construct the feature reference given the policy and data in self.
    pub fn build_feature_reference(self, is_pd: bool) -> Result<(Arc<FeatureReference>, bool)> {
        use CmoSetAction::{Builtin, InCmoSet, InFeatureRef, MergeBuiltin};
        let tenx_cmos: bool;
        let (feature_ref, tenx_cmos) = match self {
            Builtin => {
                let mut fref = FeatureReference::default();
                let multi_defs: Vec<_> = load_default_cmo_set()?
                    .into_iter()
                    .map(|x| x.into_feature_def(fref.feature_defs.len()))
                    .collect();
                fref.append_feature_defs(&multi_defs);
                (Arc::new(fref), true)
            }
            MergeBuiltin(mut feature_ref) => {
                let fref = Arc::make_mut(&mut feature_ref);
                let multi_defs: Vec<_> = load_default_cmo_set()?
                    .into_iter()
                    .map(|x| x.into_feature_def(fref.feature_defs.len()))
                    .collect();
                fref.append_feature_defs(&multi_defs);
                (feature_ref, true)
            }
            InCmoSet(feature_ref, cmo_set) => {
                if let Some(mut feature_ref) = feature_ref {
                    let fref = Arc::make_mut(&mut feature_ref);
                    fref.append_feature_defs(&cmo_set.feature_defs);
                    tenx_cmos = check_tenx_cmos(fref);
                    (feature_ref, tenx_cmos)
                } else {
                    tenx_cmos = check_tenx_cmos(&cmo_set);
                    (cmo_set, tenx_cmos)
                }
            }
            InFeatureRef(feature_ref) => {
                tenx_cmos = check_tenx_cmos(&feature_ref);
                (feature_ref, tenx_cmos)
            }
        };
        let num_multiplexing_tags = feature_ref.feature_defs.iter().fold(0, |acc, fdef| {
            if fdef.feature_type.is_multiplexing() {
                acc + 1
            } else {
                acc
            }
        });

        if !is_pd {
            ensure!(
                num_multiplexing_tags <= MAX_MULTIPLEXING_TAGS,
                "More than the maximum number of supported multiplexing tags \
                 ({num_multiplexing_tags} > {MAX_MULTIPLEXING_TAGS}) have been provided.",
            );
        }

        for fdef in &feature_ref.feature_defs {
            if !fdef.feature_type.is_multiplexing() {
                continue;
            }
            if fdef.id.is_empty() {
                bail!("Multiplexing Capture feature id field cannot be empty.");
            } else if let Some(pos) = fdef.id.find('*') {
                bail!(
                    r"Multiplexing Capture feature id field contains illegal character at position {}: '{}'
Multiplexing Capture feature ids may only ASCII characters, and must not use whitespace, asterisk, slash, quote or comma characters.",
                    pos + 1,
                    fdef.id,
                );
            } else if fdef.name.is_empty() {
                bail!("Multiplexing Capture feature name field cannot be empty.");
            }
        }
        Ok((feature_ref, tenx_cmos))
    }
}

/// This function returns the behavior required to prepare the feature_ref for possible
/// multiplexing libraries. If the configuration is somehow invalid, it returns an Error,
/// otherwise if None there is no need to do anything, otherwise it needs to call
/// CmoSetAction::build_feature_reference in order to construct a unified FeatureReference for
/// processing by the pipeline, which requires a unified representation of all feature barcode
/// constructs. You may only call this function after feature reference validation has already
/// taken place.
pub fn validate_cmo_set(
    feature_ref: Option<Arc<FeatureReference>>,
    cmo_set: Option<Arc<FeatureReference>>,
) -> Result<CmoSetAction> {
    // there _are_ multiplexing libraries, figure out how to proceed
    let fref_has_cmos = feature_ref.as_ref().is_some_and(|fref| {
        fref.feature_defs
            .iter()
            .any(|fdef| fdef.feature_type.is_multiplexing())
    });
    if fref_has_cmos && let Some(cmo_set) = cmo_set {
        // determine if there are any that are CMOs represented in one and not the other
        fn get_cmos(x: &FeatureReference) -> TxHashSet<(WhichRead, &str, &str)> {
            x.feature_defs
                .iter()
                .filter(|x| x.feature_type.is_multiplexing())
                .map(|x| (x.read, x.pattern.as_str(), x.sequence.as_str()))
                .collect()
        }
        let fref = get_cmos(feature_ref.as_ref().unwrap());
        let cmos = get_cmos(cmo_set.as_ref());
        if fref.symmetric_difference(&cmos).count() > 0 {
            bail!(
                "Multiplexing Capture feature definitions provided in both the [feature] \
                 reference-path and [gene-expression] cmo-set. You must put your complete list of \
                 Multiplexing Capture feature definitions in only one of these feature references."
            );
        }
        Ok(CmoSetAction::InFeatureRef(feature_ref.unwrap()))
    } else if let Some(cmo_set) = cmo_set {
        Ok(CmoSetAction::InCmoSet(feature_ref, cmo_set))
    } else if fref_has_cmos {
        Ok(CmoSetAction::InFeatureRef(feature_ref.unwrap()))
    } else if let Some(feature_ref) = feature_ref {
        Ok(CmoSetAction::MergeBuiltin(feature_ref))
    } else {
        Ok(CmoSetAction::Builtin)
    }
}

#[cfg(test)]
mod test {
    use crate::config::preflight::check_probe_barcode_id_grouping;

    #[test]
    fn test_check_probe_barcode_id_grouping() {
        let should_be_ok = vec![
            (vec!["BC001"], (true, false, false)),
            (vec!["AB001"], (false, true, false)),
            (vec!["BC025"], (false, true, false)), // pre-AB way of specifying
            (vec!["BC001"], (true, true, false)),
            (vec!["BC001", "AB001"], (true, true, false)),
            (vec!["BC001", "CR001"], (true, false, true)),
            (vec!["BC001", "AB001", "CR001"], (true, true, true)),
        ];
        for (bcs, (has_gex, has_ab, has_cr)) in should_be_ok {
            assert!(check_probe_barcode_id_grouping(&bcs, has_gex, has_ab, has_cr).is_ok());
        }

        let should_fail = vec![
            (vec![], (false, true, false)),
            (vec![], (true, false, false)),
            (vec![], (true, true, false)),
            (vec!["BC001"], (false, true, false)),
            (vec!["AB001"], (true, false, false)),
            (vec!["BC025"], (true, false, false)), // pre-AB way of specifying
            (vec!["AB001"], (true, true, false)),  // only provided a AB barcode
            (vec!["CR001"], (true, false, true)),  // only provided a CR barcode
            (vec!["AB001"], (true, false, true)),  // wrong type
            (vec!["AB001", "BC001"], (true, true, false)), // wrong order
            (vec!["CR001", "BC001"], (true, false, true)), // wrong order
            (vec!["BC001", "CR001"], (true, true, false)), // wrong type
            (vec!["BC001", "AB001", "BC002"], (true, true, false)), // too many
            (vec!["BC001", "AB001", "BC002"], (true, false, true)), // too many
            (vec!["BC001", "CR001", "AB002"], (true, true, true)), // wrong order
            (vec!["BC001", "AB001", "BC002"], (true, true, true)), // wrong type
        ];
        for (bcs, (has_gex, has_ab, has_cr)) in should_fail {
            assert!(
                check_probe_barcode_id_grouping(&bcs, has_gex, has_ab, has_cr).is_err(),
                "grouping didn't fail: {bcs:?}"
            );
        }
    }
}
