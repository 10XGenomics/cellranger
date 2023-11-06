#![allow(dead_code, unused_variables)]

use super::{
    create_feature_config, multiconst, AntigenSpecificityRow, FunctionalMapRow, Library,
    MultiConfigCsv, ProbeBarcodeIterationMode, SampleRow,
};
use crate::cmo_set::load_default_cmo_set;
use crate::config::{get_default_overhang_set, FeatureType, PROBE_BARCODE_ID_GROUPING};
use anyhow::{bail, ensure, Context, Result};
use barcode::whitelist::{categorize_multiplexing_barcode_id, MultiplexingBarcodeType};
use barcode::WhitelistSource;
use cr_types::chemistry::{ChemistryDef, ChemistryName};
use cr_types::reference::feature_extraction::FeatureExtractor;
use cr_types::reference::feature_reference::{BeamMode, FeatureConfig, FeatureReference};
use cr_types::reference::reference_info::ReferenceInfo;
use fastq_set::WhichRead;
use itertools::Itertools;
use metric::{TxHashMap, TxHashSet};
use std::collections::HashSet;
use std::fmt::{Display, Formatter};
use std::fs::File;
use std::io::BufReader;
use std::path::Path;
use std::str::FromStr;
use std::sync::Arc;
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

fn get_invalid_id_entries<'a>(
    valid_options: &TxHashSet<String>,
    entries: impl IntoIterator<Item = &'a str>,
) -> Vec<&'a str> {
    entries
        .into_iter()
        .filter(|x| !valid_options.contains(*x))
        .collect()
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
    use FeatureType::{AntibodyCapture, GeneExpression};
    let unsupported_feature_types: Vec<_> = cfg
        .libraries
        .feature_types()
        .filter(|x| !matches!(x, GeneExpression | AntibodyCapture))
        .sorted()
        .dedup()
        .collect();
    ensure!(
        unsupported_feature_types.is_empty(),
        "These library types are not compatible with Flex: {}",
        unsupported_feature_types.iter().join(", ")
    );
    Ok(())
}

pub fn check_libraries(
    cfg: &MultiConfigCsv,
    fref: Option<Arc<FeatureReference>>,
    is_pd: bool,
    hostname: &str,
    max_multiplexing_tags: usize,
) -> Result<()> {
    use Library::{Bcl2Fastq, BclProcessor};

    for (i, lib) in cfg.libraries.0.iter().enumerate() {
        let row = i + 1;
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
        let dirents =
            std::fs::read_dir(fastq_path).and_then(Iterator::collect::<std::io::Result<Vec<_>>>);
        if let Ok(dirents) = dirents {
            if dirents.is_empty() {
                bail!(
                    "Specified FASTQ directory is empty: {}",
                    fastq_path.display(),
                );
            }
        } else {
            bail!(
                "On machine: {}, permissions forbid opening specified FASTQ directory: {}",
                hostname,
                fastq_path.display(),
            );
        }
        for &feature_type in lib.feature_types() {
            use super::FeatureType::{
                AntibodyCapture, AntigenCapture, CrisprGuideCapture, Custom, GeneExpression,
                MultiplexingCapture, VDJ, VDJ_B, VDJ_T, VDJ_T_GD,
            };
            // convert from our FeatureType to the internal type used by FeatureReference
            let ftype = match feature_type {
                GeneExpression | VDJ | VDJ_T | VDJ_T_GD | VDJ_B => continue,
                AntibodyCapture => cr_types::types::FeatureType::Antibody,
                MultiplexingCapture => cr_types::types::FeatureType::Multiplexing,
                CrisprGuideCapture => cr_types::types::FeatureType::CRISPR,
                AntigenCapture => cr_types::types::FeatureType::Antigen,
                Custom => cr_types::types::FeatureType::Custom,
            };
            if let Some(fref) = &fref {
                if !fref.feature_maps.contains_key(&ftype) {
                    bail!(
                        "You declared a library with feature_type = '{}', but there are no features with that feature_type in the feature reference.",
                        feature_type,
                    );
                }
            }
            // traditionally, no error is thrown if feature reference is not provided
        }
        match lib.to_sample_def(cfg) {
            Err(_) => bail!("{}", MULTI_HELP),
            Ok(sdef) => sdef.check_fastqs(MULTI_HELP)?,
        }
        // check fastq_path exists, is a folder, is not empty, and sample indices are valid
    }

    // Check that Antigen Capture feature definitions match
    // the 10x allowed whitelist in case of _CS runs
    if cfg.libraries.has_antigen_capture() && !is_pd {
        fref.as_ref().unwrap().check_tenx_beam()?
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
            .validate_beam_feature_ref(beam_mode)?
    }

    // check that we don't _only_ have multiplexing libraries
    let library_types = cfg
        .libraries
        .0
        .iter()
        .flat_map(super::Library::feature_types)
        .collect::<TxHashSet<_>>()
        .into_iter()
        .collect::<Vec<_>>();
    if let [super::FeatureType::MultiplexingCapture] = library_types.as_slice() {
        bail!(
            "Only Multiplexing Capture feature libraries were provided, please provide additional non-Multiplexing Capture libraries."
        );
    }

    if cfg.is_rtl() && !is_pd {
        check_libraries_rtl(cfg)?;
    }

    Ok(())
}

pub fn check_samples(
    cfg: &MultiConfigCsv,
    fref: Option<Arc<FeatureReference>>,
    is_pd: bool,
    max_multiplexing_tags: usize,
) -> Result<()> {
    let Some(samples) = &cfg.samples else {
        return Ok(());
    };
    // Validate sample names against reserved words
    let sample_names = samples
        .0
        .iter()
        .map(|sample| sample.sample_id.clone())
        .collect::<TxHashSet<_>>();
    let invalid_sample_names = test_reserved_words(&sample_names);
    ensure!(
        invalid_sample_names.is_empty(),
        "Invalid sample_ids ('{}') provided, please ensure you are not using the reserved words \
        'blank', 'multiplet', or 'unassigned'.",
        invalid_sample_names.join("', '")
    );

    // Validate CMO IDs
    let multiplexing_ids = fref.map_or_else(TxHashSet::default, |x| x.multiplexing_ids());

    let invalid_multiplexing_ids = test_reserved_words(&multiplexing_ids);
    ensure!(
        invalid_multiplexing_ids.is_empty(),
        "Invalid cmo_ids ('{}') provided, please ensure you are not using the reserved words \
        'blank', 'multiplet', or 'unassigned'.",
        invalid_multiplexing_ids.join("', '")
    );

    for sample in &samples.0 {
        let invalid_entries = get_invalid_id_entries(
            &multiplexing_ids,
            sample.cmo_ids.iter().flatten().map(String::as_str),
        );
        ensure!(
            invalid_entries.is_empty(),
            "Unknown cmo_ids ('{}') provided for sample '{}', please ensure you are either using \
            valid 10x CMO IDs or are providing the correct [gene-expression] cmo-set.",
            invalid_entries.join("', '"),
            sample.sample_id
        );
    }

    // Validate probe barcode IDs
    if samples.has_probe_barcode_ids() && !is_pd {
        let gex = cfg.gene_expression.as_ref().unwrap();
        let chemistry_name = gex
            .chemistry
            .name()
            .map(ChemistryName::from_str)
            .transpose()
            .ok()
            .flatten()
            .unwrap_or(ChemistryName::MFRP);
        let chemistry_def = ChemistryDef::named(chemistry_name);
        let probe_barcode_whitelist = chemistry_def.barcode_whitelist().probe();
        let probe_barcode_whitelist_ids: Vec<_> =
            WhitelistSource::from_spec(probe_barcode_whitelist, true, None)?.get_ids()?;
        for sample in &samples.0 {
            let invalid_entries = get_invalid_id_entries(
                &probe_barcode_whitelist_ids
                    .iter()
                    .map(ToString::to_string)
                    .collect(),
                sample
                    .sample_barcode_ids(ProbeBarcodeIterationMode::All)
                    .unwrap(),
            );
            ensure!(
                invalid_entries.is_empty(),
                "Unknown probe_barcode_ids ('{}') provided for sample '{}', please ensure you are \
                using IDs from the following list of valid 10x probe barcode IDs: {}",
                invalid_entries.join("', '"),
                sample.sample_id,
                probe_barcode_whitelist_ids.join(", "),
            );
            for grouping in sample
                .sample_barcode_id_groupings(ProbeBarcodeIterationMode::All)
                .unwrap()
            {
                check_probe_barcode_id_grouping(
                    &grouping,
                    cfg.libraries.has_gene_expression(),
                    cfg.libraries.has_antibody_capture(),
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

    // CELLRANGER-7549
    if samples.has_overhang_ids() && !is_pd {
        bail!("[samples] section contains an invalid column: overhang_ids");
    }

    // Validate OH IDs
    if samples.has_overhang_ids() {
        let default_overhang_ids = get_default_overhang_set()
            .iter()
            .map(ToString::to_string)
            .collect();
        for sample in &samples.0 {
            let invalid_entries = get_invalid_id_entries(
                &default_overhang_ids,
                sample
                    .sample_barcode_ids(ProbeBarcodeIterationMode::All)
                    .unwrap(),
            );
            ensure!(
                invalid_entries.is_empty(),
                "Unknown overhang_ids ('{}') provided for sample '{}'.",
                invalid_entries.join("', '"),
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
) -> Result<()> {
    use MultiplexingBarcodeType::{Antibody, RTL};

    let barcode_types = barcode_ids
        .iter()
        .copied()
        .map(categorize_multiplexing_barcode_id)
        .collect::<Vec<_>>();
    match (has_gex_lib, has_ab_lib) {
        (false, false) => Ok(()),
        (true, false) => match_one(&barcode_types, RTL),
        (false, true) => match_one(&barcode_types, Antibody),
        (true, true) => {
            // can match either a single RTL barcode or a pair of two
            match &barcode_types[..] {
                [] => bail!("no barcode ID provided"),
                [RTL] | [RTL, Antibody] => Ok(()),
                [other_type] => bail!("expected a {RTL} barcode ID, not {other_type}"),
                [Antibody, RTL] => bail!(
                    "when pairing {RTL} and {Antibody} barcode IDs, \
                     provide the {RTL} barcode ID first"
                ),
                anything_else => bail!(
                    "expected either a {RTL} barcode ID or a pair of \
                     {RTL} + {Antibody} barcode IDs, not {}",
                    anything_else.iter().join(" + "),
                ),
            }
        }
    }
}

/// Match a single barcode of the provided type.
fn match_one(
    barcode_types: &[MultiplexingBarcodeType],
    expected: MultiplexingBarcodeType,
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
        if seen_values.contains(&v) {
            bail!(
                "In the [{}] section, there are two different entries with {} = '{}'. Please give each entry a unique value in the {} column. If a sample was tagged with multiple CMOs, these values should be entered as pipe-separated entries in the config file (e.g. CMO1|CMO2).",
                section,
                field,
                v,
                field,
            );
        } else {
            seen_values.push(v);
        }
    }

    Ok(())
}

pub fn check_gem_wells(libs: &[Library]) -> Result<()> {
    assert!(!libs.is_empty());
    let gem_wells = libs
        .iter()
        .map(|x| x.gem_well().0)
        .sorted()
        .dedup()
        .collect::<Vec<_>>();
    if gem_wells[0] != 1 {
        bail!(
            "[{}] gem_well numbering must start at 1",
            multiconst::LIBRARIES
        );
    }
    for (i, &gw) in gem_wells.iter().enumerate() {
        let expected = i as u16 + 1;
        if gw != expected {
            bail!(
                "[{}] gem_wells must be numbered contiguously starting from 1, missing gem_well: {}",
                multiconst::LIBRARIES,
                expected
            );
        }
    }
    Ok(())
}

pub fn check_physical_library_ids(libs: &[Library]) -> Result<()> {
    use super::FeatureType::VDJ;
    let mut mapping = TxHashMap::default();
    let mut num_vdj_libs = 0;
    for (i, lib) in libs.iter().enumerate() {
        if lib.feature_types() == [VDJ] {
            num_vdj_libs += 1;
            continue;
        }
        if let Some((j, pli)) = mapping.insert(
            (lib.feature_types(), lib.gem_well()),
            (i, lib.physical_library_id()),
        ) {
            if pli != lib.physical_library_id() {
                // TODO: update this message when gem_well comes along
                bail!(
                    r#"[{}] invalid physical_library_id on row {}: '{}'
This feature_type ({}) has already been given a physical_library_id: '{}' on row {}"#,
                    multiconst::LIBRARIES,
                    i + 2,
                    lib.physical_library_id(),
                    lib.feature_types()
                        .iter()
                        .map(|ft| format!("'{ft}'"))
                        .join(", "),
                    pli,
                    j + 2
                );
            }
        }
    }
    if num_vdj_libs > 3 {
        bail!(
            r#"[{}] invalid physical_library_id(s): found {} VDJ libraries, but we expect at most 3 VDJ libraries"#,
            multiconst::LIBRARIES,
            num_vdj_libs,
        );
    }
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
                "Your {} does not contain the expected files, or they are not readable. Please check your reference folder on {}.",
                ref_ctx,
                hostname,
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
            bail!(
                "Your {} does not appear to be indexed. Please run `mkref`",
                ref_ctx
            );
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
    let fref = FeatureReference::new(
        &ReferenceInfo::default(),
        &Transcriptome::dummy(),
        Some(rdr),
        None,
        None,
        None,
        feature_config,
    )?;
    let fref = Arc::new(fref);
    let _ = FeatureExtractor::new(fref.clone(), None, None)?;
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
    max_multiplexing_tags: usize,
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

    let (feature_reference, tenx_cmos) = validate_cmo_set(feature_reference, cmo_set)?
        .build_feature_reference(is_pd, max_multiplexing_tags)?;
    Ok((Some(feature_reference), Some(tenx_cmos)))
}

fn check_cmo_set<D: Display>(
    ref_ctx: &D,
    ref_path: &Path,
    hostname: &str,
) -> Result<Arc<FeatureReference>> {
    let fref = check_feature_reference(ref_ctx, ref_path, hostname, None)?;
    for (i, fdef) in fref.feature_defs.iter().enumerate() {
        if fdef.feature_type != cr_types::types::FeatureType::Multiplexing {
            bail!("All CMO set definitions must be of feature_type 'Multiplexing Capture'");
        }
    }
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
    use FeatureType::{VDJ, VDJ_B, VDJ_T, VDJ_T_GD};
    let mut has_multiplexing = false;
    let mut has_vdj = false;
    let mut has_antigen = false;
    let mut has_gex = false;
    let mut has_crispr = false;
    let mut uniq_vdj_features: HashSet<&FeatureType> = HashSet::new();
    for lib in libraries {
        has_multiplexing |= lib.is_multiplexing();
        has_vdj |= lib.is_vdj();
        has_antigen |= lib.is_antigen();
        has_gex |= lib.is_gex();
        has_crispr |= lib.is_crispr();
        for feature_type in lib.feature_types() {
            match feature_type {
                VDJ | VDJ_T | VDJ_B | VDJ_T_GD => _ = uniq_vdj_features.insert(feature_type),
                _ => continue,
            }
        }
    }
    if has_multiplexing && has_vdj {
        bail!(
            "[{}] The combination of VDJ libraries and Multiplexing Capture libraries is not supported.",
            multiconst::LIBRARIES
        );
    }
    if has_antigen && !(has_gex & has_vdj) {
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

pub fn check_antigen_specificity(controls: &[AntigenSpecificityRow]) -> Result<()> {
    let mut control_ids = TxHashSet::default();
    for control in controls {
        // If the set did have this value present, insert returns false
        if !control_ids.insert(&control.control_id) {
            bail!("Duplicate {} 'control_id'!", control.control_id);
        }
    }
    let alleles = controls
        .iter()
        .filter_map(|x| x.mhc_allele.as_ref())
        .collect::<TxHashSet<_>>();
    if alleles.is_empty() && controls.len() > 1 {
        bail!(
            "More than one value for 'control_id' specified without 'mhc_allele' information. \
            For BCR Antigen Capture experiments only a single control is valid, for TCR Antigen Capture both 'control_id' and 'mhc_allele' are required."
        );
    }
    if !alleles.is_empty() && (alleles.len() != controls.len()) {
        bail!(
            "For TCR Antigen Capture experiments values in both 'control_id' and 'mhc_allele' are required \
            and only a single 'control_id' can be specified for an 'mhc_allele'. "
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
                bail!(
                    "Duplicate {} 'feature_id' in [feature-functional-map] section!",
                    feature_id
                );
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
    let tenx_cmos = tenx_cmos
        .iter()
        .map(|x| (x.read, &x.pattern, &x.sequence))
        .collect::<TxHashSet<_>>();
    for fdef in &feature_ref.feature_defs {
        if fdef.feature_type != cr_types::types::FeatureType::Multiplexing {
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

impl CmoSetAction {
    /// Construct the feature reference given the policy and data in self.
    pub fn build_feature_reference(
        self,
        is_pd: bool,
        max_multiplexing_tags: usize,
    ) -> Result<(Arc<FeatureReference>, bool)> {
        use CmoSetAction::{Builtin, InCmoSet, InFeatureRef, MergeBuiltin};
        let tenx_cmos: bool;
        let (feature_ref, tenx_cmos) = match self {
            Builtin => {
                let mut fref = FeatureReference::default();
                let multi_defs = load_default_cmo_set()?
                    .into_iter()
                    .map(|x| x.into_feature_def(fref.feature_defs.len()))
                    .collect::<Vec<_>>();
                fref.append_feature_defs(&multi_defs);
                (Arc::new(fref), true)
            }
            MergeBuiltin(mut feature_ref) => {
                let fref = Arc::make_mut(&mut feature_ref);
                let multi_defs = load_default_cmo_set()?
                    .into_iter()
                    .map(|x| x.into_feature_def(fref.feature_defs.len()))
                    .collect::<Vec<_>>();
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
            if fdef.feature_type == cr_types::types::FeatureType::Multiplexing {
                acc + 1
            } else {
                acc
            }
        });
        if num_multiplexing_tags > max_multiplexing_tags && !is_pd {
            bail!(
                "More than the maximum number of supported multiplexing tags ({} > {}) have been provided.",
                num_multiplexing_tags,
                max_multiplexing_tags
            );
        }
        for fdef in &feature_ref.feature_defs {
            if fdef.feature_type != cr_types::types::FeatureType::Multiplexing {
                continue;
            }
            if fdef.id.is_empty() {
                bail!("Multiplexing Capture feature id field cannot be empty.");
            } else if let Some(pos) = fdef.id.find('*') {
                bail!(
                    r#"Multiplexing Capture feature id field contains illegal character at position {}: '{}'
Multiplexing Capture feature ids may only ASCII characters, and must not use whitespace, asterisk, slash, quote or comma characters."#,
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
    use cr_types::types::FeatureType::Multiplexing;
    // there _are_ multiplexing libraries, figure out how to proceed
    let fref_has_cmos = feature_ref
        .as_ref()
        .map(|fref| {
            fref.feature_defs
                .iter()
                .any(|fdef| fdef.feature_type == Multiplexing)
        })
        .unwrap_or(false);
    if fref_has_cmos && cmo_set.is_some() {
        // determine if there are any that are CMOs represented in one and not the other
        fn get_cmos(x: &FeatureReference) -> TxHashSet<(WhichRead, &str, &str)> {
            x.feature_defs
                .iter()
                .filter(|x| x.feature_type == Multiplexing)
                .map(|x| (x.read, x.pattern.as_str(), x.sequence.as_str()))
                .collect::<TxHashSet<_>>()
        }
        let fref = get_cmos(feature_ref.as_ref().unwrap());
        let cmos = get_cmos(cmo_set.as_ref().unwrap());
        if fref.symmetric_difference(&cmos).count() > 0 {
            bail!(
                "Multiplexing Capture feature definitions provided in both the [feature] \
                 reference-path and [gene-expression] cmo-set. You must put your complete list of \
                 Multiplexing Capture feature definitions in only one of these feature references."
            );
        }
        Ok(CmoSetAction::InFeatureRef(feature_ref.unwrap()))
    } else if cmo_set.is_some() {
        Ok(CmoSetAction::InCmoSet(feature_ref, cmo_set.unwrap()))
    } else if fref_has_cmos {
        Ok(CmoSetAction::InFeatureRef(feature_ref.unwrap()))
    } else if feature_ref.is_some() {
        Ok(CmoSetAction::MergeBuiltin(feature_ref.unwrap()))
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
            (vec!["BC001"], (true, false)),
            (vec!["AB001"], (false, true)),
            (vec!["BC025"], (false, true)), // pre-AB way of specifying
            (vec!["BC001"], (true, true)),
            (vec!["BC001", "AB001"], (true, true)),
        ];
        for (bcs, (has_gex, has_ab)) in should_be_ok {
            assert!(check_probe_barcode_id_grouping(&bcs, has_gex, has_ab).is_ok());
        }

        let should_fail = vec![
            (vec![], (false, true)),
            (vec![], (true, false)),
            (vec![], (true, true)),
            (vec!["BC001"], (false, true)),
            (vec!["AB001"], (true, false)),
            (vec!["BC025"], (true, false)), // pre-AB way of specifying
            (vec!["AB001"], (true, true)),  // only provided a AB barcode
            (vec!["OH001"], (true, false)), // wrong type
            (vec!["OH001"], (false, true)), // wrong type
            (vec!["OH001"], (true, true)),  // wrong type
            (vec!["AB001", "BC001"], (true, true)), // wrong order
            (vec!["BC001", "OH001"], (true, true)), // wrong type
            (vec!["BC001", "AB001", "BC002"], (true, true)), // too many
        ];
        for (bcs, (has_gex, has_ab)) in should_fail {
            assert!(
                check_probe_barcode_id_grouping(&bcs, has_gex, has_ab).is_err(),
                "grouping didn't fail: {bcs:?}"
            );
        }
    }
}
