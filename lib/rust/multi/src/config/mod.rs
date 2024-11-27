/// everything in a MultiConfigCsv
pub(crate) mod csv;
pub(crate) mod parse;
pub mod preflight;
pub(crate) mod scsv;

use self::csv::CsvParser;
use self::parse::{parse_prefixed_range, parse_range, parse_vec, Parse, ParseCtx};
use self::preflight::{
    check_antigen_specificity, check_duplicate_libraries, check_duplicate_sample_barcode_ids,
    check_duplicate_samples, check_feature_functional_map, check_gem_wells,
    check_library_combinations, check_physical_library_ids,
};
use self::scsv::{section_csv, Section, SectionHdr, Span, XtraData};
use crate::config::multiconst::FUNCTIONAL_MAP;
use crate::config::preflight::check_library_chemistries;
use crate::config::samplesconst::GLOBAL_MINIMUM_UMIS;
use anyhow::{anyhow, bail, ensure, Context, Result};
use barcode::whitelist::BarcodeId;
use barcode::WhitelistSource;
use cloud_utils;
use cr_types::cell_annotation::CellAnnotationModel;
use cr_types::chemistry::{
    AutoChemistryName, AutoOrRefinedChemistry, ChemistryName, ChemistrySpecs,
};
use cr_types::constants::DEFAULT_MIN_CRISPR_UMI_THRESHOLD;
use cr_types::reference::feature_reference::{
    BeamMode, FeatureConfig, SpecificityControls, MHC_ALLELE, NO_ALLELE,
};
use cr_types::reference::probe_set_reference::TargetSetFile;
use cr_types::sample_def::{FastqMode, SampleDef};
use cr_types::types::{BarcodeMultiplexingType, CellLevel, CrMultiGraph, ReadLevel};
use cr_types::{AlignerParam, FeatureBarcodeType, LibraryType, TargetingMethod, VdjChainType};
use fastq_set::filenames::FastqDef;
use itertools::{process_results, Itertools};
use martian::{AsMartianPrimaryType, MartianPrimaryType};
use martian_derive::martian_filetype;
use metric::{TxHashMap, TxHashSet};
use regex::Regex;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::convert::{AsRef, TryFrom};
use std::fmt::{Display, Formatter};
use std::fs::File;
use std::hash::Hash;
use std::io::{BufReader, Read};
use std::iter::FromIterator;
use std::path::{Path, PathBuf};
use std::str::FromStr;
use strum_macros::{Display as EnumDisplay, EnumString};

const MIN_FORCE_CELLS: usize = 10;

const ERROR_INCLUDE_INTRONS_WITH_RTL: &str = "The [gene-expression] section specifies the parameter include-introns, which is not valid for Flex chemistries.";

const DEFAULT_OVERHANG_WL: &str = "overhang";

type NomErr<'a> = nom::Err<nom::error::Error<Span<'a>>>;

pub const CONTROL_ID: &str = "control_id";
pub const AG_SPEC_REQ_HDRS: &[&str] = &[CONTROL_ID];
pub const AG_SPEC_OPT_HDRS: &[&str] = &[MHC_ALLELE];
pub const FUNCTIONAL_NAME: &str = "functional_name";
pub const FEATURE_IDS: &str = "feature_ids";
pub const FUNC_MAP_REQ_HDRS: &[&str] = &[FUNCTIONAL_NAME, FEATURE_IDS];
pub const FUNC_MAP_OPT_HDRS: &[&str] = &[];
pub const SEPARATOR: &str = "|";

/// Return an iterator of any characters in the provided string that do not match the pattern.
/// The invalid characters are deduplicated.
macro_rules! invalid_chars {
    ($input:expr, $valid:pat) => {
        $input
            .chars()
            .filter(|c: &char| !matches!(c, $valid))
            .unique()
            .collect::<Vec<_>>()
    };
}

pub struct Bool(bool);

impl FromStr for Bool {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self> {
        let s = s.to_lowercase();
        match s.as_str() {
            "0" | "f" | "false" => Ok(Bool(false)),
            "1" | "t" | "true" => Ok(Bool(true)),
            _ => bail!("invalid boolean, must be one of [true, false]"),
        }
    }
}

impl From<Bool> for bool {
    fn from(b: Bool) -> bool {
        b.0
    }
}

pub struct Ident(String);

impl AsRef<str> for Ident {
    fn as_ref(&self) -> &str {
        &self.0
    }
}

impl Display for Ident {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        f.pad(self.0.as_str())
    }
}

impl From<Ident> for String {
    fn from(v: Ident) -> Self {
        v.0
    }
}

impl FromStr for Ident {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self> {
        ensure!(!s.is_empty(), "must be non-empty");
        let invalid = invalid_chars!(s, 'A'..='Z' | 'a'..='z' | '0'..='9' | '_' | '-');
        ensure!(
            invalid.is_empty(),
            "invalid character(s): '{}', must contain only \
             letters (A-Z and a-z), digits (0-9), underscore (_), and hyphen (-)",
            invalid.iter().join("', '")
        );
        Ok(Ident((*s).to_string()))
    }
}

pub struct SampleId(String);

impl Display for SampleId {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        f.pad(self.0.as_str())
    }
}

impl From<SampleId> for String {
    fn from(v: SampleId) -> Self {
        v.0
    }
}

impl FromStr for SampleId {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self> {
        const MAX_ID_LEN: usize = 64;
        let Ident(t) = Ident::from_str(s)
            .map_err(|err| anyhow!("{err} and be no more than {MAX_ID_LEN} characters long"))?;
        ensure!(
            t.len() <= MAX_ID_LEN,
            "must be no more than {MAX_ID_LEN} characters long"
        );
        Ok(SampleId(t))
    }
}

pub struct MhcAllele(String);

impl From<MhcAllele> for String {
    fn from(v: MhcAllele) -> Self {
        v.0
    }
}

impl FromStr for MhcAllele {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self> {
        let invalid = invalid_chars!(s, 'A'..='Z' | 'a'..='z' | '0'..='9' | '_' | '-' | '*' | ':');
        ensure!(
            invalid.is_empty(),
            "invalid character(s): '{}', must match [A-Za-z0-9_-*:]+",
            invalid.iter().join("', '")
        );
        Ok(MhcAllele((*s).to_string()))
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct AtLeastOne(pub usize);

impl FromStr for AtLeastOne {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self> {
        let v: usize = s.parse()?;
        ensure!(v >= 1, "must be >= 1");
        Ok(AtLeastOne(v))
    }
}

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Probability(pub f64);

impl FromStr for Probability {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self> {
        let p: f64 = s.parse()?;
        ensure!((0.0..=1.0).contains(&p), "must be in [0, 1]");
        Ok(Probability(p))
    }
}

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Percent(pub f64);

impl FromStr for Percent {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self> {
        let p: f64 = s.parse()?;
        ensure!((0.0..=100.0).contains(&p), "must be in [0, 100]");
        Ok(Percent(p))
    }
}

/// A named collection of chemistries, zero or one per library type.
/// This concept doesn't escape the configuration parsing stage, thus why it is
/// declared here instead of with the rest of the chemistries.
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize, EnumDisplay, EnumString,
)]
#[strum(ascii_case_insensitive)]
pub enum ChemistrySet {
    #[serde(rename = "MFRP")]
    #[strum(to_string = "MFRP")]
    Mfrp,
    #[serde(rename = "MFRP-R1")]
    #[strum(to_string = "MFRP-R1")]
    MfrpR1,
    #[serde(rename = "SC3Pv3")]
    #[strum(to_string = "SC3Pv3")]
    ThreePrimeV3,
    #[serde(rename = "SC3Pv3HT")]
    #[strum(to_string = "SC3Pv3HT")]
    ThreePrimeV3HT,
    #[serde(rename = "SC3Pv4")]
    #[strum(to_string = "SC3Pv4")]
    ThreePrimeV4,
    #[serde(rename = "SC3Pv4OCM")]
    #[strum(to_string = "SC3Pv4OCM")]
    ThreePrimeV4OCM,
}

impl ChemistrySet {
    /// Return true if all of the chemistries in this set are RTL.
    /// False if not. None if ambiguous (though it is unlikely that any
    /// meaningful chemistry set is ambiguous in this regard).
    fn is_rtl(&self) -> Option<bool> {
        match self {
            Self::Mfrp | Self::MfrpR1 => Some(true),
            Self::ThreePrimeV3
            | Self::ThreePrimeV3HT
            | Self::ThreePrimeV4
            | Self::ThreePrimeV4OCM => Some(false),
        }
    }

    fn is_mfrp(&self) -> bool {
        match self {
            Self::Mfrp | Self::MfrpR1 => true,
            Self::ThreePrimeV3
            | Self::ThreePrimeV3HT
            | Self::ThreePrimeV4
            | Self::ThreePrimeV4OCM => false,
        }
    }

    /// Return the specific chemistry for the provided library type.
    pub fn chemistry_for_library_type(&self, library_type: LibraryType) -> Result<ChemistryName> {
        match self {
            Self::Mfrp => match library_type {
                LibraryType::Gex => Some(ChemistryName::MFRP_RNA),
                LibraryType::Antibody => Some(ChemistryName::MFRP_Ab),
                LibraryType::Crispr => Some(ChemistryName::MFRP_CRISPR),
                _ => None,
            },
            Self::MfrpR1 => match library_type {
                LibraryType::Gex => Some(ChemistryName::MFRP_RNA_R1),
                LibraryType::Antibody => Some(ChemistryName::MFRP_Ab_R1),
                _ => None,
            },
            Self::ThreePrimeV3 => match library_type {
                LibraryType::Gex => Some(ChemistryName::ThreePrimeV3PolyA),
                LibraryType::Antibody
                | LibraryType::Crispr
                | LibraryType::Cellplex
                | LibraryType::Custom => Some(ChemistryName::ThreePrimeV3CS1),
                _ => None,
            },
            Self::ThreePrimeV3HT => match library_type {
                LibraryType::Gex => Some(ChemistryName::ThreePrimeV3HTPolyA),
                LibraryType::Antibody
                | LibraryType::Crispr
                | LibraryType::Custom
                | LibraryType::Cellplex => Some(ChemistryName::ThreePrimeV3HTCS1),
                _ => None,
            },
            Self::ThreePrimeV4 => match library_type {
                LibraryType::Gex => Some(ChemistryName::ThreePrimeV4PolyA),
                LibraryType::Antibody | LibraryType::Custom => Some(ChemistryName::ThreePrimeV4CS1),
                _ => None,
            },
            Self::ThreePrimeV4OCM => match library_type {
                LibraryType::Gex => Some(ChemistryName::ThreePrimeV4PolyAOCM),
                LibraryType::Antibody | LibraryType::Custom => {
                    Some(ChemistryName::ThreePrimeV4CS1OCM)
                }
                _ => None,
            },
        }
        .ok_or_else(|| {
            anyhow!("The chemistry set {self} does not support the library type {library_type}.")
        })
    }
}

/// Specify a chemistry.
/// Can be an auto detection mode, a specific manual chemistry, or a
/// chemistry set to be unpacked later.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(untagged)]
pub enum ChemistryParam {
    AutoOrRefined(AutoOrRefinedChemistry),
    Set(ChemistrySet),
}

impl AsMartianPrimaryType for ChemistryParam {
    fn as_martian_primary_type() -> MartianPrimaryType {
        MartianPrimaryType::Str
    }
}

impl ChemistryParam {
    #[allow(non_upper_case_globals)]
    /// Alias for custom chemistry.
    pub const Custom: Self = Self::AutoOrRefined(AutoOrRefinedChemistry::Custom);

    fn is_rtl(&self) -> Option<bool> {
        match self {
            Self::AutoOrRefined(spec) => spec.is_rtl(),
            Self::Set(set) => set.is_rtl(),
        }
    }

    /// Return true if this chemistry is in the MFRP family.
    fn is_mfrp(&self) -> bool {
        match self {
            Self::AutoOrRefined(spec) => spec.is_mfrp(),
            Self::Set(set) => set.is_mfrp(),
        }
    }
}

impl std::fmt::Display for ChemistryParam {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::AutoOrRefined(chem) => chem.fmt(f),
            Self::Set(set) => set.fmt(f),
        }
    }
}

impl FromStr for ChemistryParam {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if let Ok(set) = ChemistrySet::from_str(s) {
            return Ok(Self::Set(set));
        }
        Ok(Self::AutoOrRefined(
            ParseAutoOrRefinedChemistry::from_str(s)?.0,
        ))
    }
}

/// Wrapper type needed to interface with our parsing machinery.
struct ParseAutoOrRefinedChemistry(AutoOrRefinedChemistry);

impl FromStr for ParseAutoOrRefinedChemistry {
    type Err = anyhow::Error;

    #[allow(clippy::enum_glob_use)]
    fn from_str(chemistry: &str) -> Result<Self> {
        use AutoChemistryName::*;
        use AutoOrRefinedChemistry::*;
        use ChemistryName::*;
        Ok(Self(match chemistry.to_ascii_lowercase().as_str() {
            // TODO: it would be ideal to avoid repeating the string definitions
            // of the chemistries here, but the aliases and case-insensitive
            // behavior are specific to this parsing and it should remain
            // confined to this module.

            // NOTE: if you add or modify options here, you must update the
            // equivalent matcher in the telemetry configuration.
            "auto" => Auto(Count),
            "custom" => Refined(Custom),
            "threeprime" => Auto(ThreePrime),
            "fiveprime" => Auto(FivePrime),
            "sc3pv1" => Refined(ThreePrimeV1),
            "sc3pv2" => Refined(ThreePrimeV2),
            "sc3pv3-polya" => Refined(ThreePrimeV3PolyA),
            "sc3pv3-cs1" => Refined(ThreePrimeV3CS1),
            "sc3pv4-polya" => Refined(ThreePrimeV4PolyA),
            "sc3pv4-cs1" => Refined(ThreePrimeV4CS1),
            "sc3pv4-polya-ocm" => Refined(ThreePrimeV4PolyAOCM),
            "sc3pv4-cs1-ocm" => Refined(ThreePrimeV4CS1OCM),
            // Removed from CS - once disabled for PD, remove this as a valid
            // parse option and move the preflight check here.
            "sc3pv3lt" => Refined(ThreePrimeV3LT),
            "sc3pv3ht-polya" => Refined(ThreePrimeV3HTPolyA),
            "sc3pv3ht-cs1" => Refined(ThreePrimeV3HTCS1),
            "sc5p-pe" | "sc5ppe" => Refined(FivePrimePE),
            "sc5p-pe-v3" | "sc5ppev3" => Refined(FivePrimePEV3),
            "sc5p-pe-ocm-v3" | "sc5ppeocmv3" => Refined(FivePrimePEOCMV3),
            "sc5p-r2" | "sc5pr2" => Refined(FivePrimeR2),
            "sc5p-r2-v3" | "sc5pr2v3" => Refined(FivePrimeR2V3),
            "sc5p-r2-ocm-v3" | "sc5pr2ocmv3" => Refined(FivePrimeR2OCMV3),
            "sc5pht" => Refined(FivePrimeHT),
            "sc-fb" | "scfb" => Refined(FeatureBarcodingOnly),
            "sfrp" => Refined(SFRP),
            "sfrp-no-trim-r2" => Refined(SfrpNoTrimR2),
            "mfrp-rna" => Refined(MFRP_RNA),
            "mfrp-ab" => Refined(MFRP_Ab),
            "mfrp-crispr" => Refined(MFRP_CRISPR),
            "mfrp-47" => Refined(MFRP_47),
            "mfrp-uncollapsed" => Refined(MFRP_uncollapsed),
            "mfrp-rna-r1" => Refined(MFRP_RNA_R1),
            "mfrp-ab-r1" => Refined(MFRP_Ab_R1),
            "mfrp-ab-r2pos50" => Refined(MFRP_Ab_R2pos50),
            "mfrp-r1-48-uncollapsed" => Refined(MFRP_R1_48_uncollapsed),
            "mfrp-r1-no-trim-r2" => Refined(MfrpR1NoTrimR2),
            "arc-v1" => Refined(ArcV1),
            _ => bail!("unknown chemistry: {chemistry}"),
        }))
    }
}

/// The gene-expression parameters in the experiment CSV
#[derive(Debug, Default, Serialize)]
pub struct GeneExpressionParams {
    pub reference_path: Option<PathBuf>,
    pub probe_set: Vec<TargetSetFile>,
    pub emptydrops_minimum_umis: Option<usize>,
    pub global_minimum_umis: Option<usize>,
    pub max_mito_percent: Option<f64>,
    pub r1_length: Option<usize>,
    pub r2_length: Option<usize>,
    pub chemistry: Option<ChemistryParam>,
    pub expect_cells: Option<usize>,
    pub force_cells: Option<usize>,
    pub no_secondary_analysis: bool,
    pub include_introns: bool,
    pub check_library_compatibility: bool,
    pub aligner: Option<AlignerParam>,
    pub create_bam: bool,
    pub filter_probes: Option<bool>,
    pub filter_high_occupancy_gems: bool,
    pub cmo_set: Option<PathBuf>,
    pub min_assignment_confidence: Option<f64>,
    pub barcode_sample_assignment: Option<PathBuf>,
    /// Select which model is used by the cell annotation service. The default string is "default" but will be validated
    pub cell_annotation_model: Option<CellAnnotationModel>,
    pub skip_cell_annotation: bool,
    pub tenx_cloud_token_path: Option<String>,
}

impl GeneExpressionParams {
    /// Return the probe set.
    pub fn probe_set(&self) -> &[TargetSetFile] {
        &self.probe_set
    }

    pub fn has_probe_set(&self) -> bool {
        !self.probe_set.is_empty()
    }

    /// Return the targeting_method.
    pub fn targeting_method(&self) -> Option<TargetingMethod> {
        self.has_probe_set()
            .then_some(TargetingMethod::TemplatedLigation)
    }

    fn has_force_cells(&self) -> bool {
        self.force_cells.is_some()
    }

    fn has_expect_cells(&self) -> bool {
        self.expect_cells.is_some()
    }

    pub fn has_expect_or_force_cells(&self) -> bool {
        self.has_expect_cells() || self.has_force_cells()
    }

    pub fn has_global_minimum_umis(&self) -> bool {
        self.global_minimum_umis.is_some()
    }

    pub fn invalid_parameter_with_antigen_capture(&self) -> Option<String> {
        if self.has_probe_set() {
            Some("probe-set".to_owned())
        } else if self.cmo_set.is_some() {
            Some("cmo-set".to_owned())
        } else if self.min_assignment_confidence.is_some() {
            Some("min-assignment-confidence".to_owned())
        } else if self.barcode_sample_assignment.is_some() {
            Some("barcode-sample-assignment".to_owned())
        } else {
            None
        }
    }

    /// Return whether the config has a probe-set or RTL chemistry.
    pub fn is_rtl(&self) -> bool {
        self.targeting_method() == Some(TargetingMethod::TemplatedLigation)
            || self.chemistry.as_ref().and_then(ChemistryParam::is_rtl) == Some(true)
    }

    /// Return either the specified token path or the default one
    pub fn get_tenx_cloud_token_path(&self) -> Option<String> {
        if self.tenx_cloud_token_path.is_none() {
            if let Ok(default_token_path) = cloud_utils::default_token_path() {
                return Some(default_token_path);
            }
            eprintln!("{}", cloud_utils::CELL_ANNOTATION_HOMEDIR_MSG);
            return None;
        }
        self.tenx_cloud_token_path.clone()
    }
}

/// A combinator to turn Some("") into None, useful for parsing logic
fn empty_is_none<'a, 'b>(s: &'a Span<'b>) -> Option<&'a Span<'b>> {
    if s.fragment().is_empty() {
        None
    } else {
        Some(s)
    }
}

impl<'a> TryFrom<&Section<'a>> for GeneExpressionParams {
    type Error = anyhow::Error;

    fn try_from(sec: &Section<'a>) -> Result<Self> {
        let ctx = ParseCtx::Hdr(sec.name);
        let mut reference_path: Option<PathBuf> = None;
        let mut probe_set: Vec<TargetSetFile> = Vec::new();
        let mut filter_probes: Option<bool> = None;
        let mut emptydrops_minimum_umis = None;
        let mut global_minimum_umis = None;
        let mut max_mito_percent = None;
        let mut r1_length: Option<usize> = None;
        let mut r2_length: Option<usize> = None;
        let mut chemistry = None;
        let mut expect_cells: Option<usize> = None;
        let mut force_cells: Option<usize> = None;
        let mut no_secondary_analysis = false;
        let mut include_introns: Option<bool> = None;
        let mut check_library_compatibility = true;
        let mut aligner: Option<AlignerParam> = None;
        let mut create_bam = None;
        let mut cmo_set: Option<PathBuf> = None;
        let mut min_assignment_confidence: Option<f64> = None;
        let mut barcode_sample_assignment: Option<PathBuf> = None;
        let mut cell_annotation_model: Option<CellAnnotationModel> = None;
        let mut skip_cell_annotation: bool = false;
        let mut tenx_cloud_token_path: Option<String> = None;
        let mut filter_high_occupancy_gems = true;
        for row in &sec.rows {
            if row.is_empty() {
                continue;
            }
            let param = row[0].fragment().to_ascii_lowercase();
            let ctx = ctx.with_col(param.as_str());
            match param.as_str() {
                "ref" | "reference" | "reference-path" => {
                    if let Some(path) = row.get(1).and_then(empty_is_none) {
                        reference_path = Some(path.parse::<PathBuf>(ctx)?);
                    }
                }
                "probe-set" => {
                    if let Some(path) = row.get(1).and_then(empty_is_none) {
                        probe_set.push(TargetSetFile::from(path.parse::<PathBuf>(ctx)?));
                    }
                }
                "filter-probes" => {
                    if let Some(val) = row.get(1).and_then(empty_is_none) {
                        filter_probes = Some(val.parse::<Bool>(ctx)?.into());
                    }
                }
                "emptydrops-minimum-umis" => {
                    if let Some(val) = row.get(1).and_then(empty_is_none) {
                        emptydrops_minimum_umis = val.parse::<usize>(ctx)?.into();
                    }
                }
                "global-minimum-umis" => {
                    if let Some(val) = row.get(1).and_then(empty_is_none) {
                        global_minimum_umis = val.parse::<usize>(ctx)?.into();
                    }
                }
                "max-mito-percent" => {
                    if let Some(val) = row.get(1).and_then(empty_is_none) {
                        max_mito_percent = Some(val.parse::<Percent>(ctx)?.0);
                    }
                }
                "r1-length" => {
                    if let Some(val) = row.get(1).and_then(empty_is_none) {
                        r1_length = Some(val.parse::<AtLeastOne>(ctx)?.0);
                    }
                }
                "r2-length" => {
                    if let Some(val) = row.get(1).and_then(empty_is_none) {
                        r2_length = Some(val.parse::<AtLeastOne>(ctx)?.0);
                    }
                }
                "chemistry" => {
                    if let Some(val) = row.get(1).and_then(empty_is_none) {
                        chemistry = Some(val.parse::<ChemistryParam>(ctx)?);
                    }
                }
                "expect-cells" => {
                    if let Some(val) = row.get(1).and_then(empty_is_none) {
                        expect_cells = Some(val.parse::<usize>(ctx)?);
                    }
                }
                "force-cells" => {
                    if let Some(val) = row.get(1).and_then(empty_is_none) {
                        let parsed_val = val.parse::<usize>(ctx)?;
                        ensure!(
                            parsed_val >= MIN_FORCE_CELLS,
                            "The 'force-cells' parameter specified under the '[gene-expression]' \
                             section needs to be at least {}. The value you have specified is {} \
                             which is too low.",
                            MIN_FORCE_CELLS,
                            parsed_val
                        );
                        force_cells = Some(parsed_val);
                    }
                }
                "introns" | "includeintrons" | "include-introns" => {
                    if let Some(val) = row.get(1).and_then(empty_is_none) {
                        include_introns = Some(val.parse::<Bool>(ctx)?.into());
                    }
                }
                "check_library_compatibility" | "check-library-compatibility" => {
                    if let Some(val) = row.get(1).and_then(empty_is_none) {
                        check_library_compatibility = val.parse::<Bool>(ctx)?.into();
                    }
                }
                "aligner" => {
                    if let Some(val) = row.get(1).and_then(empty_is_none) {
                        aligner = Some(val.parse::<AlignerParam>(ctx)?);
                    }
                }
                "nosecondary" | "no-secondary" | "no-secondary-analysis" => {
                    if let Some(val) = row.get(1).and_then(empty_is_none) {
                        no_secondary_analysis = val.parse::<Bool>(ctx)?.into();
                    }
                }
                "createbam" | "create-bam" => {
                    if let Some(val) = row.get(1).and_then(empty_is_none) {
                        create_bam = Some(val.parse::<Bool>(ctx)?.0);
                    }
                }
                "cmoset" | "cmo-set" => {
                    if let Some(val) = row.get(1).and_then(empty_is_none) {
                        cmo_set = Some(val.parse::<PathBuf>(ctx)?);
                    }
                }
                "min-assignment-confidence" | "min-assignment-conf" | "min-asgmt-conf" => {
                    if let Some(val) = row.get(1).and_then(empty_is_none) {
                        min_assignment_confidence = Some(val.parse::<Probability>(ctx)?.0);
                    }
                }
                "barcode-sample-assignment" => {
                    if let Some(val) = row.get(1).and_then(empty_is_none) {
                        barcode_sample_assignment = Some(val.parse::<PathBuf>(ctx)?);
                    }
                }
                "cell-annotation-model" => {
                    if let Some(val) = row.get(1).and_then(empty_is_none) {
                        cell_annotation_model = Some(val.parse::<CellAnnotationModel>(ctx)?);
                    }
                }
                "skip-cell-annotation" => {
                    if let Some(val) = row.get(1).and_then(empty_is_none) {
                        skip_cell_annotation = val.parse::<Bool>(ctx)?.into();
                    }
                }
                "tenx-cloud-token-path" => {
                    if let Some(val) = row.get(1).and_then(empty_is_none) {
                        tenx_cloud_token_path = Some(val.parse::<String>(ctx)?);
                    }
                }
                "filter-high-occupancy-gems" => {
                    if let Some(val) = row.get(1).and_then(empty_is_none) {
                        filter_high_occupancy_gems = val.parse::<Bool>(ctx)?.into();
                    }
                }
                "nobam" | "no-bam" => bail!("{ctx} no-bam has been replaced with create-bam"),
                _ => {
                    bail!(
                        "{ctx} unknown parameter '{}' provided at line: {}, col: {}",
                        row[0].fragment(),
                        row[0].location_line(),
                        row[0].get_utf8_column(),
                    );
                }
            }
        }

        ensure!(
            !(expect_cells.is_some() && force_cells.is_some()),
            "{ctx} only one of force-cells or expect-cells is allowed.",
        );

        if filter_probes.is_some() {
            ensure!(
                !probe_set.is_empty(),
                "{ctx} filter-probes requires a probe-set.",
            );
        }

        if !probe_set.is_empty() {
            ensure!(
                cmo_set.is_none(),
                "{ctx} When probe-set is specified, cmo-set is an invalid parameter.",
            );
            ensure!(
                min_assignment_confidence.is_none(),
                "{ctx} When probe-set is specified, min-assignment-confidence \
                 is an invalid parameter.",
            );
            ensure!(
                barcode_sample_assignment.is_none(),
                "{ctx} When probe-set is specified, barcode-sample-assignment \
                 is an invalid parameter.",
            );
        }

        // Disallow setting include-introns for RTL chemistries.
        if chemistry.as_ref().and_then(ChemistryParam::is_rtl) == Some(true) {
            ensure!(include_introns.is_none(), ERROR_INCLUDE_INTRONS_WITH_RTL);
        }

        // Filter probes by default.
        if !probe_set.is_empty() && filter_probes.is_none() {
            filter_probes = Some(true);
        }

        let Some(create_bam) = create_bam else {
            bail!("{ctx} create-bam is a required parameter")
        };

        Ok(GeneExpressionParams {
            reference_path,
            probe_set,
            filter_probes,
            emptydrops_minimum_umis,
            global_minimum_umis,
            max_mito_percent,
            r1_length,
            r2_length,
            chemistry,
            expect_cells,
            force_cells,
            no_secondary_analysis,
            include_introns: include_introns.unwrap_or(true),
            check_library_compatibility,
            aligner,
            create_bam,
            cmo_set,
            min_assignment_confidence,
            barcode_sample_assignment,
            cell_annotation_model,
            skip_cell_annotation,
            tenx_cloud_token_path,
            filter_high_occupancy_gems,
        })
    }
}

#[derive(Debug, Serialize, Deserialize)]
pub struct FeatureParams {
    pub reference_path: Option<PathBuf>,
    pub r1_length: Option<usize>,
    pub r2_length: Option<usize>,
    pub filter_aggregates: bool,
    pub min_crispr_umi: usize,
}

impl<'a> TryFrom<&Section<'a>> for FeatureParams {
    type Error = anyhow::Error;

    fn try_from(sec: &Section<'a>) -> Result<Self> {
        let ctx = ParseCtx::Hdr(sec.name);
        let mut reference_path: Option<PathBuf> = None;
        let mut r1_length: Option<usize> = None;
        let mut r2_length: Option<usize> = None;
        let mut filter_aggregates = true;
        let mut min_crispr_umi = DEFAULT_MIN_CRISPR_UMI_THRESHOLD;

        for row in &sec.rows {
            if row.is_empty() {
                continue;
            }
            let param = row[0].fragment().to_ascii_lowercase();
            let ctx = ctx.with_col(param.as_str());
            match param.as_str() {
                "ref" | "reference" | "reference-path" => {
                    if let Some(val) = row.get(1).and_then(empty_is_none) {
                        reference_path = Some(val.parse::<PathBuf>(ctx)?);
                    }
                }
                "r1-length" => {
                    if let Some(val) = row.get(1).and_then(empty_is_none) {
                        r1_length = Some(val.parse::<AtLeastOne>(ctx)?.0);
                    }
                }
                "r2-length" => {
                    if let Some(val) = row.get(1).and_then(empty_is_none) {
                        r2_length = Some(val.parse::<AtLeastOne>(ctx)?.0);
                    }
                }
                "filter-aggregates" => {
                    if let Some(val) = row.get(1).and_then(empty_is_none) {
                        filter_aggregates = val.parse::<Bool>(ctx)?.into();
                    }
                }
                "min-crispr-umi" => {
                    if let Some(val) = row.get(1).and_then(empty_is_none) {
                        min_crispr_umi = val.parse::<AtLeastOne>(ctx)?.0;
                    }
                }
                _ => {
                    bail!(
                        "{ctx} unknown parameter '{}' provided at line: {}, col: {}",
                        row[0].fragment(),
                        row[0].location_line(),
                        row[0].get_utf8_column(),
                    );
                }
            }
        }
        Ok(FeatureParams {
            reference_path,
            r1_length,
            r2_length,
            filter_aggregates,
            min_crispr_umi,
        })
    }
}

#[derive(Debug, Serialize, Deserialize)]
pub struct VdjParams {
    pub reference_path: PathBuf,
    pub inner_enrichment_primers: Option<PathBuf>,
    pub r1_length: Option<usize>,
    pub r2_length: Option<usize>,
    pub multiplet_filter: Option<bool>,
    pub shared_contig_filter: Option<bool>,
    pub umi_baseline_filter: Option<bool>,
    pub min_contig_length: Option<usize>,
    pub skip_clonotyping: Option<bool>,
}

impl<'a> TryFrom<&Section<'a>> for VdjParams {
    type Error = anyhow::Error;

    fn try_from(sec: &Section<'a>) -> Result<Self> {
        let ctx = ParseCtx::Hdr(sec.name);
        let mut reference_path = PathBuf::new();
        let mut inner_enrichment_primers: Option<PathBuf> = None;
        let mut r1_length: Option<usize> = None;
        let mut r2_length: Option<usize> = None;
        let mut multiplet_filter = None;
        let mut shared_contig_filter = None;
        let mut umi_baseline_filter = None;
        let mut min_contig_length = None;
        let mut skip_clonotyping = None;
        for row in &sec.rows {
            if row.is_empty() {
                continue;
            }
            let param = row[0].fragment().to_ascii_lowercase();
            let ctx = ctx.with_col(row[0].fragment());
            match param.as_str() {
                "ref" | "reference" | "reference-path" => {
                    reference_path = row
                        .get(1)
                        .ok_or_else(|| {
                            anyhow!("{ctx} no value provided for '{}'", row[0].fragment())
                        })?
                        .parse::<PathBuf>(ctx)?;
                }
                "inner-enrichment-primers" => {
                    if let Some(path) = row.get(1).and_then(empty_is_none) {
                        inner_enrichment_primers = Some(path.parse::<PathBuf>(ctx)?);
                    }
                }
                "r1-length" => {
                    if let Some(val) = row.get(1).and_then(empty_is_none) {
                        r1_length = Some(val.parse::<AtLeastOne>(ctx)?.0);
                    }
                }
                "r2-length" => {
                    if let Some(val) = row.get(1).and_then(empty_is_none) {
                        r2_length = Some(val.parse::<AtLeastOne>(ctx)?.0);
                    }
                }
                "multiplet-filter" => {
                    if let Some(val) = row.get(1).and_then(empty_is_none) {
                        multiplet_filter = Some(val.parse::<Bool>(ctx)?.into());
                    }
                }
                "shared-contig-filter" => {
                    if let Some(val) = row.get(1).and_then(empty_is_none) {
                        shared_contig_filter = Some(val.parse::<Bool>(ctx)?.into());
                    }
                }
                "umi-baseline-filter" => {
                    if let Some(val) = row.get(1).and_then(empty_is_none) {
                        umi_baseline_filter = Some(val.parse::<Bool>(ctx)?.into());
                    }
                }
                "min-contig-length" => {
                    if let Some(val) = row.get(1).and_then(empty_is_none) {
                        min_contig_length = Some(val.parse::<usize>(ctx)?);
                    }
                }
                "skip-clonotyping" => {
                    if let Some(val) = row.get(1).and_then(empty_is_none) {
                        skip_clonotyping = Some(val.parse::<Bool>(ctx)?.into());
                    }
                }
                _ => {
                    bail!(
                        "{} unknown parameter '{}' provided at line: {}, col: {}",
                        sec.name,
                        row[0].fragment(),
                        row[0].location_line(),
                        row[0].get_utf8_column(),
                    );
                }
            }
        }
        ensure!(
            !reference_path.as_os_str().is_empty(),
            "{ctx} reference is missing"
        );
        Ok(VdjParams {
            reference_path,
            inner_enrichment_primers,
            r1_length,
            r2_length,
            multiplet_filter,
            shared_contig_filter,
            umi_baseline_filter,
            min_contig_length,
            skip_clonotyping,
        })
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(transparent)]
pub struct GemWell(pub u16);

impl FromStr for GemWell {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self> {
        Ok(GemWell(s.parse::<u16>()?))
    }
}

/// Wrapper type for case-insensitive parsing.
struct ParseFeatureType(LibraryType);

impl FromStr for ParseFeatureType {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(ParseFeatureType(match s.to_ascii_lowercase().as_str() {
            "gene expression" => LibraryType::Gex,
            "vdj" => LibraryType::VdjAuto,
            "vdj-t" => LibraryType::Vdj(VdjChainType::VdjT),
            "vdj-t-gd" => LibraryType::Vdj(VdjChainType::VdjTGD),
            "vdj-b" => LibraryType::Vdj(VdjChainType::VdjB),
            "antibody capture" => LibraryType::Antibody,
            "multiplexing capture" => LibraryType::Cellplex,
            "crispr guide capture" => LibraryType::Crispr,
            "antigen capture" => LibraryType::Antigen,
            "custom" => LibraryType::Custom,
            _ => bail!("unknown feature_type '{s}'"),
        }))
    }
}

// TODO: we have LaneSpec, use it here
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(from = "Option<Vec<usize>>", into = "Option<Vec<usize>>")]
pub enum Lanes {
    Any,
    Lanes(Vec<usize>),
}

impl Lanes {
    pub fn overlaps(&self, other: &Lanes) -> bool {
        use self::Lanes::{Any, Lanes};
        match (self, other) {
            (Any, Any) | (Any, _) | (_, Any) => true,
            (Lanes(l1), Lanes(l2)) => {
                let l1: TxHashSet<_> = l1.iter().copied().collect();
                let l2: TxHashSet<_> = l2.iter().copied().collect();
                // if the intersection isn't empty, we overlap
                l1.intersection(&l2).count() > 0
            }
        }
    }
}

impl From<Option<Vec<usize>>> for Lanes {
    fn from(o: Option<Vec<usize>>) -> Self {
        match o {
            None => Lanes::Any,
            Some(ls) => Lanes::Lanes(ls),
        }
    }
}

impl From<Lanes> for Option<Vec<usize>> {
    fn from(val: Lanes) -> Option<Vec<usize>> {
        match val {
            Lanes::Any => None,
            Lanes::Lanes(ls) => Some(ls),
        }
    }
}

// this is why we have Lanes implemented here instead of LaneSpec:
//   I cannot define TryFrom (from std) for an external type
impl<'a, 'b, 'c> TryFrom<(ParseCtx<'a, SectionHdr<'b>>, Span<'c>)> for Lanes {
    type Error = anyhow::Error;

    fn try_from((ctx, span): (ParseCtx<'a, SectionHdr<'b>>, Span<'c>)) -> Result<Self> {
        if &span.fragment().to_ascii_lowercase() == "any" {
            return Ok(Lanes::Any);
        }

        let Ok(lanes) = parse_vec(span.clone()) else {
            bail!(
                "{ctx} has invalid {} '{}' at line: {}, col: {}",
                libsconst::LANES,
                span.fragment(),
                span.location_line(),
                span.get_utf8_column(),
            );
        };

        Ok(Lanes::Lanes(process_results(
            lanes.into_iter().map(parse_range).flatten_ok(),
            |iter| iter.sorted().dedup().collect(),
        )?))
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub enum Library {
    Bcl2Fastq {
        fastq_id: String,
        fastqs: PathBuf,
        lanes: Lanes,
        physical_library_id: String,
        feature_type: LibraryType,
        gem_well: GemWell,
        subsample_rate: Option<f64>,
        ilmn_fastq_id: Option<String>, // Fastq ID might be edited to handle case of duplicate fastq IDs referring to different data
        chemistry: Option<AutoOrRefinedChemistry>,
    },
    BclProcessor {
        fastq_path: PathBuf,
        sample_indices: Vec<String>,
        lanes: Lanes,
        physical_library_id: String,
        feature_type: LibraryType,
        gem_well: GemWell,
        subsample_rate: Option<f64>,
        chemistry: Option<AutoOrRefinedChemistry>,
    },
}

impl Library {
    pub fn lanes(&self) -> &Lanes {
        match self {
            Self::Bcl2Fastq { lanes, .. } => lanes,
            Self::BclProcessor { lanes, .. } => lanes,
        }
    }

    pub fn library_type(&self) -> LibraryType {
        match self {
            Self::Bcl2Fastq { feature_type, .. } => *feature_type,
            Self::BclProcessor { feature_type, .. } => *feature_type,
        }
    }

    pub fn gem_well(&self) -> GemWell {
        match self {
            Self::Bcl2Fastq { gem_well, .. } => *gem_well,
            Self::BclProcessor { gem_well, .. } => *gem_well,
        }
    }

    pub fn physical_library_id(&self) -> &str {
        match self {
            Self::Bcl2Fastq {
                physical_library_id,
                ..
            } => physical_library_id.as_str(),
            Self::BclProcessor {
                physical_library_id,
                ..
            } => physical_library_id.as_str(),
        }
    }

    pub fn chemistry(&self) -> Option<AutoOrRefinedChemistry> {
        match self {
            Self::Bcl2Fastq { chemistry, .. } => *chemistry,
            Self::BclProcessor { chemistry, .. } => *chemistry,
        }
    }

    pub fn is_multiplexing(&self) -> bool {
        self.library_type() == LibraryType::Cellplex
    }

    pub fn is_vdj(&self) -> bool {
        self.library_type().is_vdj()
    }

    pub fn is_antigen(&self) -> bool {
        self.library_type().is_fb_type(FeatureBarcodeType::Antigen)
    }

    pub fn is_gex(&self) -> bool {
        self.library_type().is_gex()
    }

    pub fn is_crispr(&self) -> bool {
        self.library_type().is_fb_type(FeatureBarcodeType::Crispr)
    }

    pub fn overlaps(&self, other: &Library) -> bool {
        match self {
            Self::Bcl2Fastq {
                fastq_id: fastq_id1,
                fastqs: fastq_path1,
                lanes: lanes1,
                ilmn_fastq_id: ilmn_fastq_id1,
                ..
            } => {
                if let Self::Bcl2Fastq {
                    fastq_id: fastq_id2,
                    fastqs: fastq_path2,
                    lanes: lanes2,
                    ilmn_fastq_id: ilmn_fastq_id2,
                    ..
                } = other
                {
                    let fastq_id1 = ilmn_fastq_id1.as_ref().unwrap_or(fastq_id1);
                    let fastq_id2 = ilmn_fastq_id2.as_ref().unwrap_or(fastq_id2);
                    if fastq_id1 != fastq_id2 || fastq_path1 != fastq_path2 {
                        return false;
                    }
                    return lanes1.overlaps(lanes2);
                }
                false
            }
            Self::BclProcessor {
                fastq_path: fastq_path1,
                sample_indices,
                lanes: lanes1,
                ..
            } => {
                // something bogus is happening here requiring me to provide State type parameter
                let sample_indices1 = TxHashSet::from_iter(sample_indices);
                if let Self::BclProcessor {
                    fastq_path: fastq_path2,
                    sample_indices,
                    lanes: lanes2,
                    ..
                } = other
                {
                    if fastq_path1 != fastq_path2 {
                        return false;
                    }
                    let sample_indices2 = TxHashSet::from_iter(sample_indices);
                    if sample_indices1.intersection(&sample_indices2).count() == 0 {
                        return false;
                    }
                    return lanes1.overlaps(lanes2);
                }
                false
            }
        }
    }

    pub fn to_fastq_def(&self) -> Result<FastqDef> {
        self.to_sample_def()?.get_fastq_def()
    }

    pub fn to_sample_def(&self) -> Result<SampleDef> {
        let lanes = self.lanes().clone().into();
        let library_type = Some({
            // CELLRANGER-7889: this was previously an domain boundary where all
            // VDJ chain types would be collapsed down into a single "VDJ" type.
            // Retain this invariant for the moment by collapsing all VDJ chain types
            // down to Auto.
            let library_type = self.library_type();
            if library_type.is_vdj() {
                LibraryType::VdjAuto
            } else {
                library_type
            }
        });
        match self {
            Self::Bcl2Fastq {
                fastqs,
                fastq_id,
                gem_well,
                subsample_rate,
                ilmn_fastq_id,
                ..
            } => {
                let sample_name = match ilmn_fastq_id {
                    Some(orig) => orig.clone(),
                    None => fastq_id.clone(),
                };
                Ok(SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    gem_group: Some(gem_well.0),
                    lanes,
                    library_type,
                    r1_length: None,
                    r2_length: None,
                    read_path: fastqs.clone(),
                    sample_indices: None,
                    sample_names: Some(vec![sample_name]),
                    subsample_rate: *subsample_rate,
                    fastq_id: Some(fastq_id.clone()),
                })
            }
            Self::BclProcessor {
                fastq_path,
                sample_indices,
                gem_well,
                subsample_rate,
                ..
            } => {
                // Heuristic to use the flowcell id as the fastq_id based on how internal fastqs
                // are along paths like:
                // /somepath/pipestances/HKWWVDSXC/BCL_PROCESSOR_PD/...
                let pattern = r"([^/]+)/BCL_PROCESSOR_PD";
                let re = Regex::new(pattern).unwrap();

                let flow_cell_id = re
                    .captures(fastq_path.to_str().unwrap())
                    .and_then(|caps| caps.get(1).map(|m| String::from(m.as_str())));
                Ok(SampleDef {
                    fastq_mode: FastqMode::BCL_PROCESSOR,
                    gem_group: Some(gem_well.0),
                    lanes,
                    library_type,
                    r1_length: None,
                    r2_length: None,
                    read_path: fastq_path.clone(),
                    sample_indices: Some(sample_indices.clone()),
                    sample_names: None,
                    subsample_rate: *subsample_rate,
                    fastq_id: flow_cell_id,
                })
            }
        }
    }
}

mod libsconst {
    pub const FASTQ_ID: &str = "fastq_id";
    pub const FASTQS: &str = "fastqs";
    pub const PHYSICAL_LIBRARY_ID: &str = "physical_library_id";
    pub const FEATURE_TYPES: &str = "feature_types";
    pub const LIBS_REQ_HDRS: &[&str] = &[FASTQ_ID, FASTQS, FEATURE_TYPES];

    pub const LANES: &str = "lanes";
    pub const GEM_WELL: &str = "gem_well";
    pub const SUBSAMPLE_RATE: &str = "subsample_rate";
    pub const CHEMISTRY: &str = "chemistry";
    #[cfg(not(test))]
    pub const LIBS_OPT_HDRS: &[&str] = &[PHYSICAL_LIBRARY_ID, LANES, SUBSAMPLE_RATE, CHEMISTRY];
    #[cfg(test)]
    pub const LIBS_OPT_HDRS: &[&str] = &[
        PHYSICAL_LIBRARY_ID,
        LANES,
        GEM_WELL,
        SUBSAMPLE_RATE,
        CHEMISTRY,
    ];

    pub const FASTQ_PATH: &str = "fastq_path";
    pub const SAMPLE_INDICES: &str = "sample_indices";
    pub const LIBS_INT_REQ_HDRS: &[&str] = &[FASTQ_PATH, SAMPLE_INDICES, FEATURE_TYPES];
    #[cfg(not(test))]
    pub const LIBS_INT_OPT_HDRS: &[&str] = &[PHYSICAL_LIBRARY_ID, LANES, SUBSAMPLE_RATE, CHEMISTRY];
    #[cfg(test)]
    pub const LIBS_INT_OPT_HDRS: &[&str] = &[
        PHYSICAL_LIBRARY_ID,
        LANES,
        GEM_WELL,
        SUBSAMPLE_RATE,
        CHEMISTRY,
    ];
}

// TODO: should these just be SampleDef?
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
#[serde(transparent)]
pub struct LibrariesCsv(pub Vec<Library>);

impl LibrariesCsv {
    /// Return an iterator over the library types.
    pub fn library_types(&self) -> impl Iterator<Item = LibraryType> + '_ {
        self.0.iter().map(Library::library_type)
    }

    /// Return true if any libraries are gene expression.
    pub fn has_gene_expression(&self) -> bool {
        self.has_library_type(LibraryType::Gex)
    }

    /// Return true if any libraries are Antibody/Crispr/Antigen/Custom.
    /// Return false for CMO multiplexing, which comes with a built-in list.
    pub fn has_feature_barcode(&self) -> bool {
        self.library_types().any(|ft| {
            let Some(fb_type) = ft.feature_barcode_type() else {
                return false;
            };
            fb_type != FeatureBarcodeType::Multiplexing
        })
    }

    /// Return true if any libraries are of the provided type.
    pub fn has_library_type(&self, library_type: LibraryType) -> bool {
        self.library_types().any(|ft| ft == library_type)
    }

    pub fn has_feature_barcode_type(&self, feature_barcode_type: FeatureBarcodeType) -> bool {
        self.has_library_type(LibraryType::FeatureBarcodes(feature_barcode_type))
    }

    /// Return true if any libraries are antibody capture.
    pub fn has_antibody_capture(&self) -> bool {
        self.has_feature_barcode_type(FeatureBarcodeType::Antibody)
    }

    /// Return true if any libraries are CRISPR guide capture.
    pub fn has_crispr_guide_capture(&self) -> bool {
        self.has_feature_barcode_type(FeatureBarcodeType::Crispr)
    }

    /// Return true if any libraries are antigen capture.
    pub fn has_antigen_capture(&self) -> bool {
        self.has_feature_barcode_type(FeatureBarcodeType::Antigen)
    }

    /// Return true if any libraries are CMO multiplexing.
    pub fn has_multiplexing(&self) -> bool {
        self.has_feature_barcode_type(FeatureBarcodeType::Multiplexing)
    }

    pub fn has_vdj(&self) -> bool {
        self.library_types().any(|x| x.is_vdj())
    }

    pub fn has_vdj_t_or_gd(&self) -> bool {
        self.has_library_type(LibraryType::Vdj(VdjChainType::VdjT))
            || self.has_library_type(LibraryType::Vdj(VdjChainType::VdjTGD))
    }

    pub fn has_vdj_b(&self) -> bool {
        self.has_library_type(LibraryType::Vdj(VdjChainType::VdjB))
    }

    pub fn number_of_vdj(&self) -> usize {
        self.library_types()
            .filter_map(|x| x.vdj_chain_type())
            .unique()
            .count()
    }

    pub fn beam_mode(&self) -> Option<BeamMode> {
        if self.has_vdj_t_or_gd() {
            return Some(BeamMode::BeamT);
        } else if self.has_vdj_b() {
            return Some(BeamMode::BeamAB);
        }
        None
    }

    /// Return the map of chemistry specs, for any libraries that provided one.
    /// Return None if chemistry was not provided at the library level.
    ///
    /// Validation must ensure that this is either None, or a map containing
    /// an entry for every library. Validation also must ensure that different
    /// libraries of the same feature type have matching chemistries.
    pub fn chemistry_specs(&self) -> Option<ChemistrySpecs> {
        let specs: ChemistrySpecs = self
            .0
            .iter()
            .filter_map(|lib| lib.chemistry().map(|chem| (lib.library_type(), chem)))
            .collect();
        if specs.is_empty() {
            return None;
        }
        Some(specs)
    }
}

fn library_type_short(ft: LibraryType) -> &'static str {
    match ft {
        LibraryType::Gex => "GEX",
        LibraryType::FeatureBarcodes(fb) => {
            #[allow(clippy::enum_glob_use)]
            use FeatureBarcodeType::*;
            match fb {
                Antibody => "ABC",
                Antigen => "AGC",
                Crispr => "CGC",
                Custom => "CUST",
                Multiplexing => "CMO",
            }
        }
        LibraryType::Vdj(ct) => {
            #[allow(clippy::enum_glob_use)]
            use VdjChainType::*;
            match ct {
                Auto => "VDJ",
                VdjT => "VDJT",
                VdjB => "VDJB",
                VdjTGD => "VDJTGD",
            }
        }
        LibraryType::Atac => unreachable!(),
    }
}

impl<'a> TryFrom<&Section<'a>> for LibrariesCsv {
    type Error = anyhow::Error;

    fn try_from(sec: &Section<'a>) -> Result<Self> {
        use libsconst::{
            CHEMISTRY, FASTQS, FASTQ_ID, FASTQ_PATH, FEATURE_TYPES, GEM_WELL, LANES,
            LIBS_INT_OPT_HDRS, LIBS_INT_REQ_HDRS, LIBS_OPT_HDRS, LIBS_REQ_HDRS,
            PHYSICAL_LIBRARY_ID, SAMPLE_INDICES, SUBSAMPLE_RATE,
        };
        let hdr = sec.name;
        let (parser, is_internal) = {
            if let Ok(parser) = CsvParser::new(sec.clone(), LIBS_INT_REQ_HDRS, LIBS_INT_OPT_HDRS) {
                (parser, true)
            } else {
                (
                    CsvParser::new(sec.clone(), LIBS_REQ_HDRS, LIBS_OPT_HDRS)?,
                    false,
                )
            }
        };
        let mut data = vec![];
        let mut num_vdj_libs = 0;

        // Keep track of how many times we've seen a fastq-id.
        let mut fastq_id_counts = TxHashMap::default();

        for row in parser.rows() {
            let ctx = ParseCtx::HdrRow(hdr, row + 1);
            let feature_type = parser
                .find_req(row, FEATURE_TYPES)?
                .parse::<ParseFeatureType>(ctx.with_col(FEATURE_TYPES))?
                .0;
            let lanes = parser
                .find_opt(row, LANES)?
                .and_then(empty_is_none)
                .map(|span| Lanes::try_from((ctx, span.clone())))
                .transpose()?
                .unwrap_or(Lanes::Any);
            let gem_well = if cfg!(test) {
                parser
                    .find_opt(row, GEM_WELL)?
                    .and_then(empty_is_none)
                    .map(|gw| gw.parse::<GemWell>(ctx.with_col(GEM_WELL)))
                    .transpose()?
                    .unwrap_or(GemWell(1))
            } else {
                GemWell(1)
            };
            let physical_library_id = if let Some(pli) = parser
                .find_opt(row, PHYSICAL_LIBRARY_ID)?
                .and_then(empty_is_none)
            {
                pli.parse::<Ident>(ctx.with_col(PHYSICAL_LIBRARY_ID))?
                    .to_string()
            } else {
                if feature_type == LibraryType::VdjAuto {
                    num_vdj_libs += 1;
                    ensure!(
                        num_vdj_libs == 1,
                        "{} is invalid: no more than 1 VDJ library may rely upon an \
                         auto-generated physical_library_id. Please either provide \
                         physical_libary_ids for all VDJ libraries or use more \
                         specific feature_types like VDJ-T and VDJ-B.",
                        ctx.with_col(PHYSICAL_LIBRARY_ID)
                    );
                }
                format!("{}_{}", library_type_short(feature_type), gem_well.0)
            };
            let subsample_rate = parser
                .find_opt(row, SUBSAMPLE_RATE)?
                .and_then(empty_is_none)
                .map(|sr| sr.parse::<Probability>(ctx.with_col(SUBSAMPLE_RATE)))
                .transpose()?
                .map(|sr| sr.0);
            let chemistry = parser
                .find_opt(row, CHEMISTRY)?
                .and_then(empty_is_none)
                .map(|sr| sr.parse::<ParseAutoOrRefinedChemistry>(ctx.with_col(CHEMISTRY)))
                .transpose()?
                .map(|param| param.0);

            if is_internal {
                let fastq_path = PathBuf::from(parser.find_req(row, FASTQ_PATH)?.fragment());
                let sample_indices = parser.find_req(row, SAMPLE_INDICES)?;
                let sample_indices: Vec<_> = parse_vec(sample_indices.clone())
                    .map_err(|_: NomErr<'_>| {
                        anyhow!(
                            "{ctx} has invalid {SAMPLE_INDICES} '{}' at line: {}, col: {}",
                            sample_indices.fragment(),
                            sample_indices.location_line(),
                            sample_indices.get_utf8_column(),
                        )
                    })?
                    .into_iter()
                    .map(|si| si.to_string())
                    .collect();

                data.push(Library::BclProcessor {
                    fastq_path,
                    sample_indices,
                    lanes,
                    feature_type,
                    physical_library_id,
                    gem_well,
                    subsample_rate,
                    chemistry,
                });
            } else {
                let mut fastq_id = parser
                    .find_req(row, FASTQ_ID)?
                    .parse::<Ident>(ctx.with_col(FASTQ_ID))?
                    .to_string();

                // if the same fastq ID has been used already, add a count afterward
                // hold on to the original fastq ID to use as "sample name" and pass bcl2fastq checks
                let mut ilmn_fastq_id = None;
                let fastq_count = fastq_id_counts.entry(fastq_id.clone()).or_insert(0);
                *fastq_count += 1;

                if *fastq_count > 1 {
                    ilmn_fastq_id = Some(fastq_id.clone());
                    fastq_id = format!("{fastq_id} ({fastq_count})");
                }

                let fastqs = PathBuf::from(parser.find_req(row, FASTQS)?.fragment());

                data.push(Library::Bcl2Fastq {
                    fastq_id,
                    fastqs,
                    lanes,
                    feature_type,
                    physical_library_id,
                    gem_well,
                    subsample_rate,
                    ilmn_fastq_id,
                    chemistry,
                });
            }
        }
        check_duplicate_libraries(&data)?;
        check_gem_wells(&data)?;
        check_physical_library_ids(&data)?;
        check_library_combinations(&data)?;
        check_library_chemistries(&data)?;
        Ok(LibrariesCsv(data))
    }
}

#[derive(Debug)]
pub struct GemWells(pub Vec<GemWell>);

impl<'a> TryFrom<Span<'a>> for GemWells {
    type Error = anyhow::Error;

    fn try_from(s: Span<'a>) -> Result<Self> {
        Ok(GemWells(parse_range(s)?.map(GemWell).collect::<Vec<_>>()))
    }
}

/// Multiple probe barcodes may be grouped explicitly using this character.
/// This enables cross-library declarations of different probe barcodes for
/// different library types.
pub const PROBE_BARCODE_ID_GROUPING: &str = "+";

#[derive(Debug, Default, Serialize, Deserialize)]
pub struct SampleRow {
    pub sample_id: String,
    cmo_ids: Option<Vec<String>>,
    hashtag_ids: Option<Vec<String>>,
    /// NOTE: these may be +-concatenated groupings. Access via the
    /// sample_barcode_ids method to control how they are unpacked.
    probe_barcode_ids: Option<Vec<String>>,
    overhang_ids: Option<Vec<String>>,
    pub description: String,
    pub expect_cells: Option<usize>,
    pub force_cells: Option<usize>,
    pub emptydrops_minimum_umis: Option<usize>,
    pub global_minimum_umis: Option<usize>,
    pub max_mito_percent: Option<f64>,
}

impl SampleRow {
    /// Create a new SampleRow with the given probe_barcode_id, for testing.
    pub fn from_probe_barcode_id(probe_barcode_id: &str) -> Self {
        Self {
            probe_barcode_ids: Some(vec![probe_barcode_id.to_string()]),
            ..Self::default()
        }
    }

    /// Return an iterator of the sample barcode IDs, either CMO or probe barcode IDs.
    /// Control probe barcode iteration using the provided option.
    pub(crate) fn sample_barcode_ids(
        &self,
        probe_barcode_iteration_mode: ProbeBarcodeIterationMode,
    ) -> Option<impl Iterator<Item = &str> + '_> {
        self.sample_barcode_id_groupings(probe_barcode_iteration_mode)
            .map(|bcs| bcs.into_iter().flatten())
    }

    /// Return all of the sample barcode IDs, either CMO or probe barcode IDs.
    /// Retain any existing |-separated groupings in the result.
    /// Control probe barcode iteration using the provided option.
    pub(crate) fn sample_barcode_id_groupings(
        &self,
        probe_barcode_iteration_mode: ProbeBarcodeIterationMode,
    ) -> Option<Vec<Vec<&str>>> {
        match (
            self.cmo_ids.as_ref(),
            self.hashtag_ids.as_ref(),
            self.probe_barcode_ids.as_ref(),
            self.overhang_ids.as_ref(),
        ) {
            (Some(x), None, None, None)
            | (None, Some(x), None, None)
            | (None, None, None, Some(x)) => {
                // cmo_ids or hashtag_ids or overhang_ids
                Some(x.iter().map(|v| vec![v.as_str()]).collect())
            }
            (None, None, Some(x), None) => {
                // probe_barcode_ids
                Some(
                    x.iter()
                        .map(move |bcid| {
                            // Conditionally take only the first piece of the barcode if we're excluding
                            // mapped barcodes.
                            bcid.split(PROBE_BARCODE_ID_GROUPING)
                                .take(match probe_barcode_iteration_mode {
                                    ProbeBarcodeIterationMode::All => usize::MAX,
                                    ProbeBarcodeIterationMode::Mapped => 1,
                                })
                                .collect()
                        })
                        .collect(),
                )
            }
            (None, None, None, None) => None,
            _ => panic!("found more than one source of sample barcode IDs"),
        }
    }

    /// Return the type of cell multiplexing.
    pub fn barcode_multiplexing_type(&self) -> BarcodeMultiplexingType {
        match (
            self.cmo_ids.is_some(),
            self.hashtag_ids.is_some(),
            self.probe_barcode_ids.is_some(),
            self.overhang_ids.is_some(),
        ) {
            (true, false, false, false) => BarcodeMultiplexingType::CellLevel(CellLevel::CMO),
            (false, true, false, false) => BarcodeMultiplexingType::CellLevel(CellLevel::Hashtag),
            (false, false, true, false) => BarcodeMultiplexingType::ReadLevel(ReadLevel::RTL),
            (false, false, false, true) => BarcodeMultiplexingType::ReadLevel(ReadLevel::OH),
            _ => unreachable!(),
        }
    }

    /// Return the appropriate column name.
    pub fn sample_barcode_ids_column_name(&self) -> &str {
        match self.barcode_multiplexing_type() {
            BarcodeMultiplexingType::CellLevel(CellLevel::CMO) => samplesconst::CMO_IDS,
            BarcodeMultiplexingType::CellLevel(CellLevel::Hashtag) => samplesconst::HASHTAG_IDS,
            BarcodeMultiplexingType::ReadLevel(ReadLevel::RTL) => samplesconst::PROBE_BARCODE_IDS,
            BarcodeMultiplexingType::ReadLevel(ReadLevel::OH) => samplesconst::OH_IDS,
        }
    }
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(transparent)]
pub struct SamplesCsv(pub Vec<SampleRow>);

impl SamplesCsv {
    pub fn has_probe_barcode_ids(&self) -> bool {
        self.0
            .iter()
            .any(|sample| sample.probe_barcode_ids.is_some())
    }

    pub fn has_overhang_ids(&self) -> bool {
        self.0.iter().any(|sample| sample.overhang_ids.is_some())
    }

    pub fn has_hashtag_ids(&self) -> bool {
        self.0.iter().any(|sample| sample.hashtag_ids.is_some())
    }

    pub fn has_cmo_ids(&self) -> bool {
        self.0.iter().any(|sample| sample.cmo_ids.is_some())
    }

    fn has_force_cells(&self) -> bool {
        self.0.iter().any(|sample| sample.force_cells.is_some())
    }

    fn has_expect_cells(&self) -> bool {
        self.0.iter().any(|sample| sample.expect_cells.is_some())
    }

    pub fn has_expect_or_force_cells(&self) -> bool {
        self.has_expect_cells() || self.has_force_cells()
    }

    pub fn has_emptydrops_minimum_umis(&self) -> bool {
        self.0
            .iter()
            .any(|sample| sample.emptydrops_minimum_umis.is_some())
    }

    pub fn has_global_minimum_umis(&self) -> bool {
        self.0
            .iter()
            .any(|sample| sample.global_minimum_umis.is_some())
    }

    pub fn get_expect_cells(&self) -> TxHashMap<String, Option<f64>> {
        self.0
            .iter()
            .map(|sample| {
                (
                    sample.sample_id.clone(),
                    sample.expect_cells.map(|x| x as f64),
                )
            })
            .collect()
    }

    pub fn get_force_cells(&self) -> TxHashMap<String, Option<f64>> {
        self.0
            .iter()
            .map(|sample| {
                (
                    sample.sample_id.clone(),
                    sample.force_cells.map(|x| x as f64),
                )
            })
            .collect()
    }

    pub fn get_emptydrops_minimum_umis(&self) -> TxHashMap<String, Option<f64>> {
        self.0
            .iter()
            .map(|sample| {
                (
                    sample.sample_id.clone(),
                    sample.emptydrops_minimum_umis.map(|x| x as f64),
                )
            })
            .collect()
    }

    pub fn get_global_minimum_umis(&self) -> TxHashMap<String, Option<f64>> {
        self.0
            .iter()
            .map(|sample| {
                (
                    sample.sample_id.clone(),
                    (sample.global_minimum_umis.map(|x| x as f64)),
                )
            })
            .collect()
    }

    pub fn get_max_mito_percent(&self) -> TxHashMap<String, Option<f64>> {
        self.0
            .iter()
            .map(|sample| (sample.sample_id.clone(), sample.max_mito_percent))
            .collect()
    }

    pub fn is_rtl_multiplexed(&self) -> bool {
        self.0.iter().all(|sample| {
            sample.barcode_multiplexing_type() == BarcodeMultiplexingType::ReadLevel(ReadLevel::RTL)
        })
    }

    pub fn get_cmo_sample_map(&self) -> TxHashMap<String, String> {
        let mut res = TxHashMap::default();
        for row in &self.0 {
            if let Some(ref cmo_ids) = row.cmo_ids {
                for cmo_id in cmo_ids {
                    assert!(
                        res.insert(cmo_id.clone(), row.sample_id.clone()).is_none(),
                        "cmo_id {cmo_id} used for multiple samples"
                    );
                }
            }
        }
        res
    }

    pub fn get_hashtag_sample_map(&self) -> TxHashMap<String, String> {
        let mut res = TxHashMap::default();
        for row in &self.0 {
            if let Some(ref hashtag_ids) = row.hashtag_ids {
                for hashtag_id in hashtag_ids {
                    assert!(
                        res.insert(hashtag_id.clone(), row.sample_id.clone())
                            .is_none(),
                        "hashtag_id {hashtag_id} used for multiple samples"
                    );
                }
            }
        }
        res
    }

    /// Return a mapping from source to translated probe barcode.
    /// This method handles interpreting probe barcodes concatenated with a +
    /// which indicate to the pipeline that all of the IDs following the first
    /// ID should be dynamically mapped into the barcode associated with the
    /// first ID.
    ///
    /// This is used to merge samples from different library types that
    /// have been created from the same GEM well.
    pub fn get_translated_probe_barcodes(&self) -> TxHashMap<BarcodeId, BarcodeId> {
        // Here we assume that we have already validated that there are no
        // duplicates in the probe barcodes in the entire config.
        self.iter_translated_probe_barcodes()
            .flat_map(|(target, sources)| sources.map(move |s| (s, target)))
            .collect()
    }

    /// Return an iterator over tuples of target barcode and barcodes to be translated into it.
    fn iter_translated_probe_barcodes(
        &self,
    ) -> impl Iterator<Item = (BarcodeId, impl Iterator<Item = BarcodeId> + '_)> + '_ {
        self.0
            .iter()
            .filter_map(|row| row.probe_barcode_ids.as_ref())
            .flatten()
            .filter_map(|bc| {
                let mut pieces_iter = bc.split(PROBE_BARCODE_ID_GROUPING);
                pieces_iter
                    .next()
                    .map(|target| (BarcodeId::pack(target), pieces_iter.map(BarcodeId::pack)))
            })
    }
}

/// Directive to determine how we iterate over probe barcodes.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum ProbeBarcodeIterationMode {
    /// Include all probe barcodes, including those dynamically mapped into others.
    All,
    /// Ignore probe barcodes that have been mapped into another.
    Mapped,
}

pub mod samplesconst {
    pub const SAMPLE_ID: &str = "sample_id";
    pub const CMO_IDS: &str = "cmo_ids";
    pub const HASHTAG_IDS: &str = "hashtag_ids";
    pub const PROBE_BARCODE_IDS: &str = "probe_barcode_ids";
    pub const OH_IDS: &str = "ocm_barcode_ids";
    pub const _GEM_WELLS: &str = "gem_wells";
    pub const DESCRIPTION: &str = "description";
    pub const EXPECT_CELLS: &str = "expect_cells";
    pub const FORCE_CELLS: &str = "force_cells";
    pub const EMPTYDROPS_MINIMUM_UMIS: &str = "emptydrops_minimum_umis";
    pub const GLOBAL_MINIMUM_UMIS: &str = "global_minimum_umis";
    pub const MAX_MITO_FRAC: &str = "max_mito_percent";
    pub const SAMP_REQ_HDRS: &[&str] = &[SAMPLE_ID];
    #[cfg(not(test))]
    pub const SAMP_OPT_HDRS: &[&str] = &[
        CMO_IDS,
        HASHTAG_IDS,
        OH_IDS,
        PROBE_BARCODE_IDS,
        DESCRIPTION,
        EXPECT_CELLS,
        FORCE_CELLS,
        EMPTYDROPS_MINIMUM_UMIS,
        GLOBAL_MINIMUM_UMIS,
        MAX_MITO_FRAC,
    ];
    #[cfg(test)]
    pub const SAMP_OPT_HDRS: &[&str] = &[
        CMO_IDS,
        HASHTAG_IDS,
        PROBE_BARCODE_IDS,
        OH_IDS,
        _GEM_WELLS,
        DESCRIPTION,
        EXPECT_CELLS,
        FORCE_CELLS,
        EMPTYDROPS_MINIMUM_UMIS,
        GLOBAL_MINIMUM_UMIS,
        MAX_MITO_FRAC,
    ];
}

impl<'a> TryFrom<(&TxHashSet<GemWell>, &Section<'a>)> for SamplesCsv {
    type Error = anyhow::Error;

    fn try_from((valid_gws, sec): (&TxHashSet<GemWell>, &Section<'a>)) -> Result<Self> {
        use samplesconst::{
            CMO_IDS, DESCRIPTION, EMPTYDROPS_MINIMUM_UMIS, EXPECT_CELLS, FORCE_CELLS, HASHTAG_IDS,
            MAX_MITO_FRAC, OH_IDS, PROBE_BARCODE_IDS, SAMPLE_ID, SAMP_OPT_HDRS, SAMP_REQ_HDRS,
            _GEM_WELLS,
        };
        let hdr = sec.name;
        let parser = CsvParser::new(sec.clone(), SAMP_REQ_HDRS, SAMP_OPT_HDRS)?;
        let mut data = vec![];
        for row in parser.rows() {
            let ctx = ParseCtx::HdrRow(hdr, row + 1);
            let sample_id: String = parser
                .find_req(row, SAMPLE_ID)?
                .parse::<SampleId>(ctx.with_col(SAMPLE_ID))?
                .into();

            let cmo_ids = parser.find_opt(row, CMO_IDS)?;
            let hashtag_ids = parser.find_opt(row, HASHTAG_IDS)?;
            let probe_barcode_ids = parser.find_opt(row, PROBE_BARCODE_IDS)?;
            let overhang_ids = parser.find_opt(row, OH_IDS)?;
            let sample_barcode_ids = cmo_ids
                .or(hashtag_ids)
                .or(probe_barcode_ids)
                .or(overhang_ids);
            let barcode_multiplexing_type = match (
                cmo_ids,
                hashtag_ids,
                probe_barcode_ids,
                overhang_ids,
            ) {
                (Some(_), None, None, None) => BarcodeMultiplexingType::CellLevel(CellLevel::CMO),
                (None, Some(_), None, None) => {
                    BarcodeMultiplexingType::CellLevel(CellLevel::Hashtag)
                }
                (None, None, Some(_), None) => BarcodeMultiplexingType::ReadLevel(ReadLevel::RTL),
                (None, None, None, Some(_)) => BarcodeMultiplexingType::ReadLevel(ReadLevel::OH),
                (None, None, None, None) => {
                    bail!(
                        "{ctx} requires either {CMO_IDS} or {HASHTAG_IDS} or {PROBE_BARCODE_IDS} or {OH_IDS} column \
                         to be specified"
                    )
                }
                (_, _, _, _) => {
                    bail!(
                        "{ctx} has more than one of the following mutually exclusive columns: \
                    {CMO_IDS}, {HASHTAG_IDS}, {PROBE_BARCODE_IDS} and/or {OH_IDS}"
                    )
                }
            };
            let expect_cells = parser.find_opt(row, EXPECT_CELLS)?;
            let force_cells = parser.find_opt(row, FORCE_CELLS)?;
            let emptydrops_minimum_umis = parser.find_opt(row, EMPTYDROPS_MINIMUM_UMIS)?;
            let global_minimum_umis = parser.find_opt(row, GLOBAL_MINIMUM_UMIS)?;
            let max_mito_percent = parser.find_opt(row, MAX_MITO_FRAC)?;
            if probe_barcode_ids.is_none() && overhang_ids.is_none() {
                ensure!(
                    expect_cells.is_none(),
                    "{ctx} has {EXPECT_CELLS} column specified \
                     without the required {PROBE_BARCODE_IDS} column",
                );
                ensure!(
                    force_cells.is_none(),
                    "{ctx} has {FORCE_CELLS} column specified \
                     without the required {PROBE_BARCODE_IDS} column",
                );
            };
            let sample_barcode_column_name = match barcode_multiplexing_type {
                BarcodeMultiplexingType::CellLevel(CellLevel::CMO) => CMO_IDS,
                BarcodeMultiplexingType::CellLevel(CellLevel::Hashtag) => HASHTAG_IDS,
                BarcodeMultiplexingType::ReadLevel(ReadLevel::RTL) => PROBE_BARCODE_IDS,
                BarcodeMultiplexingType::ReadLevel(ReadLevel::OH) => OH_IDS,
            };

            let sample_barcode_ids: Option<Vec<_>> =
                if let Some(sample_barcode_ids) = sample_barcode_ids.and_then(empty_is_none) {
                    let Ok(sample_barcode_ids_in) = parse_vec(sample_barcode_ids.clone()) else {
                        bail!(
                            "{ctx} has invalid {} '{}' at line: {}, col: {}",
                            sample_barcode_column_name,
                            sample_barcode_ids.fragment(),
                            sample_barcode_ids.location_line(),
                            sample_barcode_ids.get_utf8_column(),
                        )
                    };

                    let sample_barcode_ids_iter = sample_barcode_ids_in
                        .into_iter()
                        .map(parse_prefixed_range)
                        .flatten_ok()
                        .map(|multiplexing_id| {
                            let multiplexing_id = multiplexing_id?;
                            match barcode_multiplexing_type {
                                BarcodeMultiplexingType::ReadLevel(ReadLevel::RTL) => {
                                    validate_probe_barcode_id(multiplexing_id)
                                }
                                BarcodeMultiplexingType::CellLevel(CellLevel::CMO) => {
                                    multiplexing_id.parse::<Ident>().map(|_| multiplexing_id)
                                }
                                BarcodeMultiplexingType::CellLevel(CellLevel::Hashtag) => {
                                    multiplexing_id.parse::<Ident>().map(|_| multiplexing_id)
                                }
                                BarcodeMultiplexingType::ReadLevel(ReadLevel::OH) => {
                                    Ok(multiplexing_id)
                                }
                            }
                            .with_context(|| {
                                format!(
                                    "{ctx} has invalid {} '{}' at line: {}, col: {}",
                                    sample_barcode_column_name,
                                    sample_barcode_ids.fragment(),
                                    sample_barcode_ids.location_line(),
                                    sample_barcode_ids.get_utf8_column(),
                                )
                            })
                        });

                    Some(process_results(sample_barcode_ids_iter, |iter| {
                        iter.sorted().dedup().collect()
                    })?)
                } else {
                    None
                };

            let gem_wells_frag = None;
            let gem_wells = gem_wells_frag
                .and_then(empty_is_none)
                .map(|span| {
                    parse_vec(span.clone())
                        .map_err(|_: NomErr<'_>| {
                            anyhow!(
                                "{ctx} has invalid {_GEM_WELLS} '{}' at line: {}, col: {}",
                                span.fragment(),
                                span.location_line(),
                                span.get_utf8_column(),
                            )
                        })?
                        .into_iter()
                        .map(|gws| anyhow::Ok(parse_range::<u16>(gws)?.map(GemWell)))
                        .flatten_ok()
                        .collect()
                })
                .transpose()?
                .unwrap_or_else(|| vec![GemWell(1)]);
            for gem_well in &gem_wells {
                if let Some(gem_wells_frag) = gem_wells_frag {
                    ensure!(
                        valid_gws.contains(gem_well),
                        "{ctx} has invalid gem_well '{}' at line: {}, col: {}: \
                         please check for consistency with [{}]",
                        gem_well.0,
                        gem_wells_frag.location_line(),
                        gem_wells_frag.get_utf8_column(),
                        multiconst::LIBRARIES,
                    );
                } else {
                    ensure!(
                        valid_gws.contains(gem_well),
                        "{ctx} has invalid gem_well '{}': please check for consistency with [{}]",
                        gem_well.0,
                        multiconst::LIBRARIES
                    );
                }
            }
            let description = parser
                .find_opt(row, DESCRIPTION)?
                .map_or_else(String::default, |x| String::from(x as &str));
            let force_cells = force_cells
                .and_then(empty_is_none)
                .map(|fc| fc.parse::<usize>(ctx.with_col(FORCE_CELLS)))
                .transpose()?;
            let expect_cells = expect_cells
                .and_then(empty_is_none)
                .map(|ec| ec.parse::<usize>(ctx.with_col(EXPECT_CELLS)))
                .transpose()?;
            ensure!(
                !(force_cells.is_some() && expect_cells.is_some()),
                "{ctx} has both force-cells and expect-cells specified, only one is allowed.",
            );
            let emptydrops_minimum_umis = emptydrops_minimum_umis
                .and_then(empty_is_none)
                .map(|ec| ec.parse::<usize>(ctx.with_col(EMPTYDROPS_MINIMUM_UMIS)))
                .transpose()?;
            let global_minimum_umis = global_minimum_umis
                .and_then(empty_is_none)
                .map(|ec| ec.parse::<usize>(ctx.with_col(GLOBAL_MINIMUM_UMIS)))
                .transpose()?;
            let max_mito_percent = max_mito_percent
                .and_then(empty_is_none)
                .map(|ec| ec.parse::<Percent>(ctx.with_col(MAX_MITO_FRAC)))
                .transpose()?
                .map(|mm| mm.0);

            // Check to make sure that sample barcode IDs were provided and that all entries are non-empty
            ensure!(
                sample_barcode_ids
                    .as_ref()
                    .is_some_and(|ids| ids.iter().all(|x| !x.is_empty())),
                "{ctx} has an empty entry in the {sample_barcode_column_name} column. \
                 All entries must be non-empty.",
            );

            data.push(match barcode_multiplexing_type {
                BarcodeMultiplexingType::CellLevel(CellLevel::CMO) => SampleRow {
                    sample_id,
                    cmo_ids: sample_barcode_ids,
                    hashtag_ids: None,
                    probe_barcode_ids: None,
                    overhang_ids: None,
                    description,
                    force_cells: None,
                    expect_cells: None,
                    emptydrops_minimum_umis,
                    global_minimum_umis,
                    max_mito_percent,
                },
                BarcodeMultiplexingType::CellLevel(CellLevel::Hashtag) => SampleRow {
                    sample_id,
                    cmo_ids: None,
                    hashtag_ids: sample_barcode_ids,
                    probe_barcode_ids: None,
                    overhang_ids: None,
                    description,
                    force_cells: None,
                    expect_cells: None,
                    emptydrops_minimum_umis,
                    global_minimum_umis,
                    max_mito_percent,
                },
                BarcodeMultiplexingType::ReadLevel(ReadLevel::RTL) => SampleRow {
                    sample_id,
                    cmo_ids: None,
                    hashtag_ids: None,
                    probe_barcode_ids: sample_barcode_ids,
                    overhang_ids: None,
                    description,
                    force_cells,
                    expect_cells,
                    emptydrops_minimum_umis,
                    global_minimum_umis,
                    max_mito_percent,
                },
                BarcodeMultiplexingType::ReadLevel(ReadLevel::OH) => SampleRow {
                    sample_id,
                    cmo_ids: None,
                    hashtag_ids: None,
                    probe_barcode_ids: None,
                    overhang_ids: sample_barcode_ids,
                    description,
                    force_cells,
                    expect_cells,
                    emptydrops_minimum_umis,
                    global_minimum_umis,
                    max_mito_percent,
                },
            });
        }
        check_duplicate_sample_barcode_ids(&data)?;
        check_duplicate_samples(&data)?;
        Ok(SamplesCsv(data))
    }
}

/// Check that a probe barcode ID is valid.
/// Handles the possibility that the barcode ID is a concatenated pairing.
fn validate_probe_barcode_id(bcid: String) -> Result<String> {
    ensure!(!bcid.is_empty(), "must be non-empty");
    ensure!(
        !bcid.split(PROBE_BARCODE_ID_GROUPING).any(str::is_empty),
        "must be non-empty"
    );

    let invalid = invalid_chars!(bcid, 'A'..='Z' | '0'..='9' | '+');
    ensure!(
        invalid.is_empty(),
        "invalid character(s): '{}', must contain only \
         uppercase letters (A-Z), digits (0-9), and plus (+)",
        invalid.iter().join("', '")
    );
    Ok(bcid)
}

/// Returns the set of default overhang IDs
fn get_default_overhang_set() -> TxHashSet<BarcodeId> {
    WhitelistSource::named(DEFAULT_OVERHANG_WL, true)
        .unwrap()
        .get_ids()
        .unwrap()
        .into_iter()
        .collect()
}

#[derive(Debug, Serialize, Deserialize)]
pub struct GemWellParams {
    pub gem_well: GemWell,
    pub force_cells: Option<usize>,
    pub expect_cells: Option<usize>,
    pub vdj_force_cells: Option<usize>,
}

mod gemwellconst {
    pub const GEM_WELL: &str = "gem_well";
    pub const FORCE_CELLS: &str = "force_cells";
    pub const EXPECT_CELLS: &str = "expect_cells";
    pub const VDJ_FORCE_CELLS: &str = "vdj_force_cells";
    pub const GEM_WELL_REQ_HDRS: &[&str] = &[GEM_WELL];
    pub const GEM_WELL_OPT_HDRS: &[&str] = &[FORCE_CELLS, EXPECT_CELLS, VDJ_FORCE_CELLS];
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(transparent)]
pub struct GemWellsCsv(pub TxHashMap<GemWell, GemWellParams>);

impl<'a> TryFrom<(&TxHashSet<GemWell>, &Section<'a>)> for GemWellsCsv {
    type Error = anyhow::Error;

    fn try_from((valid_gws, sec): (&TxHashSet<GemWell>, &Section<'a>)) -> Result<Self> {
        use gemwellconst::{
            EXPECT_CELLS, FORCE_CELLS, GEM_WELL, GEM_WELL_OPT_HDRS, GEM_WELL_REQ_HDRS,
            VDJ_FORCE_CELLS,
        };
        let hdr = sec.name;
        let parser = CsvParser::new(sec.clone(), GEM_WELL_REQ_HDRS, GEM_WELL_OPT_HDRS)?;
        let mut data = TxHashMap::default();
        for row in parser.rows() {
            let ctx = ParseCtx::HdrRow(hdr, row + 1);
            let gem_well = parser
                .find_req(row, GEM_WELL)?
                .parse::<GemWell>(ctx.with_col(GEM_WELL))?;
            ensure!(
                valid_gws.contains(&gem_well),
                "{ctx} has invalid gem_well: {}, please check for consistency with [{}]",
                gem_well.0,
                multiconst::LIBRARIES,
            );
            let force_cells = parser
                .find_opt(row, FORCE_CELLS)?
                .and_then(empty_is_none)
                .map(|fc| fc.parse::<usize>(ctx.with_col(FORCE_CELLS)))
                .transpose()?;
            let expect_cells = parser
                .find_opt(row, EXPECT_CELLS)?
                .and_then(empty_is_none)
                .map(|ec| ec.parse::<usize>(ctx.with_col(EXPECT_CELLS)))
                .transpose()?;
            let vdj_force_cells = parser
                .find_opt(row, VDJ_FORCE_CELLS)?
                .and_then(empty_is_none)
                .map(|vfc| vfc.parse::<usize>(ctx.with_col(VDJ_FORCE_CELLS)))
                .transpose()?;
            let duplicate_entry = data.insert(
                gem_well,
                GemWellParams {
                    gem_well,
                    force_cells,
                    expect_cells,
                    vdj_force_cells,
                },
            );
            ensure!(
                duplicate_entry.is_none(),
                "{ctx} duplicate entry detected: {}",
                gem_well.0
            );
        }
        Ok(GemWellsCsv(data))
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct AntigenSpecificityRow {
    pub control_id: String,
    pub mhc_allele: Option<String>,
}

#[derive(Debug, Serialize, Deserialize, Clone)]
#[serde(transparent)]
pub struct AntigenSpecificityCsv(pub Vec<AntigenSpecificityRow>);

impl<'a> TryFrom<&Section<'a>> for AntigenSpecificityCsv {
    type Error = anyhow::Error;

    fn try_from(sec: &Section<'a>) -> Result<Self> {
        let hdr = sec.name;
        let parser = CsvParser::new(sec.clone(), AG_SPEC_REQ_HDRS, AG_SPEC_OPT_HDRS)?;
        let mut data = vec![];
        for row in parser.rows() {
            let ctx = ParseCtx::HdrRow(hdr, row + 1);
            let control_id = parser.find_req(row, CONTROL_ID)?.to_string();
            let mhc_allele = parser
                .find_opt(row, MHC_ALLELE)?
                .and_then(empty_is_none)
                .map(|a| a.parse::<MhcAllele>(ctx.with_col(MHC_ALLELE)))
                .transpose()?
                .map(String::from);
            data.push(AntigenSpecificityRow {
                control_id,
                mhc_allele,
            });
        }
        check_antigen_specificity(&data)?;
        Ok(AntigenSpecificityCsv(data))
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct FunctionalMapRow {
    pub functional_name: String,
    pub feature_ids: Vec<String>,
}

#[derive(Debug, Serialize, Deserialize, Clone)]
#[serde(transparent)]
pub struct FunctionalMapCsv(pub Vec<FunctionalMapRow>);

impl<'a> TryFrom<&Section<'a>> for FunctionalMapCsv {
    type Error = anyhow::Error;

    fn try_from(sec: &Section<'a>) -> Result<Self> {
        let parser = CsvParser::new(sec.clone(), FUNC_MAP_REQ_HDRS, FUNC_MAP_OPT_HDRS)?;
        let mut data = vec![];
        for row in parser.rows() {
            let functional_name = parser.find_req(row, FUNCTIONAL_NAME)?.to_string();
            let feature_ids: Vec<_> = parser
                .find_req(row, FEATURE_IDS)?
                .split(SEPARATOR)
                .map(String::from)
                .collect();
            data.push(FunctionalMapRow {
                functional_name,
                feature_ids,
            });
        }
        check_feature_functional_map(&data)?;
        Ok(FunctionalMapCsv(data))
    }
}
pub fn create_feature_config(
    antigen_specificity_csv: Option<&AntigenSpecificityCsv>,
    functional_map_csv: Option<&FunctionalMapCsv>,
    beam_mode: Option<BeamMode>,
    samples_csv: Option<&SamplesCsv>,
) -> Option<FeatureConfig> {
    let specificity_controls = create_specificity_controls(antigen_specificity_csv);

    let mut functional_map = HashMap::new();

    if let Some(functional_map_csv) = functional_map_csv {
        for functional_map_row in &functional_map_csv.0 {
            let functional_name = functional_map_row.functional_name.clone();
            for feature_id in functional_map_row.feature_ids.clone() {
                if let Some(existing_value) =
                    functional_map.insert(feature_id, functional_name.clone())
                {
                    assert_ne!(existing_value, functional_name);
                }
            }
        }
    }

    let functional_map = if !functional_map.is_empty() {
        Some(functional_map)
    } else {
        None
    };

    let hashtag_ids: Option<Vec<String>> = samples_csv.and_then(|samples| {
        if samples.has_hashtag_ids() {
            Some(
                samples
                    .0
                    .iter()
                    .filter_map(|s| s.sample_barcode_ids(ProbeBarcodeIterationMode::All))
                    .flatten()
                    .map(std::string::ToString::to_string)
                    .collect(),
            )
        } else {
            None
        }
    });

    if specificity_controls.is_some() || functional_map_csv.is_some() || hashtag_ids.is_some() {
        Some(FeatureConfig {
            specificity_controls,
            beam_mode,
            functional_map,
            hashtag_ids,
        })
    } else {
        None
    }
}

pub fn create_specificity_controls(
    antigen_specificity_csv: Option<&AntigenSpecificityCsv>,
) -> Option<SpecificityControls> {
    let mut ag_spec_controls = SpecificityControls {
        control_for_allele: HashMap::new(),
        has_mhc_allele_column: false,
    };
    if let Some(antigen_specificity_csv) = antigen_specificity_csv {
        for antigen_row in &antigen_specificity_csv.0 {
            let allele = if let Some(ref mhc_allele) = antigen_row.mhc_allele {
                ag_spec_controls.has_mhc_allele_column = true;
                mhc_allele
            } else {
                NO_ALLELE
            };
            ag_spec_controls
                .control_for_allele
                .insert(allele.to_string(), antigen_row.control_id.clone());
        }
        // Some(ag_spec_controls)
    }

    if antigen_specificity_csv.is_some() {
        Some(ag_spec_controls)
    } else {
        None
    }
}

pub mod multiconst {
    use lazy_static::lazy_static;
    use metric::TxHashSet;

    pub const GENE_EXPRESSION: &str = "gene-expression";
    pub const GEX: &str = "gex";
    pub const FEATURE: &str = "feature";
    pub const VDJ: &str = "vdj";
    pub const LIBRARIES: &str = "libraries";
    pub const LIBS: &str = "libs";
    pub const SAMPLES: &str = "samples";
    pub const GEM_WELLS: &str = "gem-wells";
    pub const GWS: &str = "gws";
    pub const ANTIGEN_SPECIFICITY: &str = "antigen-specificity";
    pub const FUNCTIONAL_MAP: &str = "feature-functional-map";

    lazy_static! {
        pub static ref VALID_SECTIONS: TxHashSet<&'static str> = {
            let mut s = TxHashSet::default();
            s.insert(GENE_EXPRESSION);
            s.insert(GEX);
            s.insert(FEATURE);
            s.insert(VDJ);
            s.insert(LIBRARIES);
            s.insert(LIBS);
            s.insert(ANTIGEN_SPECIFICITY);
            // TODO: bring this back one day
            #[cfg(test)]
            {
                s.insert(SAMPLES);
                s.insert(GEM_WELLS);
                s.insert(GWS);
            }
            s
        };
    }

    // Default value of include-introns for non-hybcap GEX analyses.
    pub const DEFAULT_INCLUDE_INTRONS: bool = true;
}

martian_filetype! { MultiConfigCsvFile, "csv" }

impl MultiConfigCsvFile {
    /// Read this multi config CSV.
    pub fn read(&self) -> Result<MultiConfigCsv> {
        MultiConfigCsv::from_csv(self)
    }
}

/// A container for the contents of an MultiConfigCsv
#[derive(Debug, Default, Serialize)]
pub struct MultiConfigCsv {
    pub gene_expression: Option<GeneExpressionParams>,
    pub feature: Option<FeatureParams>,
    pub vdj: Option<VdjParams>,
    pub libraries: LibrariesCsv,
    pub samples: Option<SamplesCsv>,
    pub gem_wells: Option<GemWellsCsv>,
    pub antigen_specificity: Option<AntigenSpecificityCsv>,
    pub functional_map: Option<FunctionalMapCsv>,
}

/// A representation of a `cellranger multi` configuration
impl MultiConfigCsv {
    /// Load a MultiConfigCsv from a path
    fn from_csv(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        MultiConfigCsv::from_reader(
            BufReader::new(File::open(path).with_context(|| path.display().to_string())?),
            XtraData::new(path),
        )
    }

    /// Load an MultiConfigCsv from any `impl Read` along with parsing `XtraData`
    pub fn from_reader<R: Read, X: Into<XtraData>>(mut reader: R, xtra: X) -> Result<Self> {
        use multiconst::{
            ANTIGEN_SPECIFICITY, FEATURE, GEM_WELLS, GENE_EXPRESSION, GEX, GWS, LIBRARIES, LIBS,
            SAMPLES, VDJ,
        };
        let mut buf = String::new();
        let _ = reader.read_to_string(&mut buf)?;
        let s = if buf.starts_with('\u{feff}') {
            &buf.as_str()[3..]
        } else {
            buf.as_str()
        };
        let input = Span::new_extra(s, xtra.into());
        let (_, mut sections) = section_csv(input).map_err(|e| match e {
            nom::Err::Error(e) | nom::Err::Failure(e) => anyhow!(
                "failed to parse CSV at line: {}, col: {}",
                e.input.location_line(),
                e.input.get_utf8_column()
            ),
            nom::Err::Incomplete(_) => {
                anyhow!("failed to parse CSV, incomplete information available to pinpoint error")
            }
        })?;
        sections.sort_by_key(|s| {
            let name = s.name.fragment().to_ascii_lowercase();
            match name.as_str() {
                GENE_EXPRESSION | GEX | VDJ | FEATURE | ANTIGEN_SPECIFICITY => 0,
                LIBRARIES | LIBS => 1,
                SAMPLES | GEM_WELLS | GWS => 2,
                _ => 3,
            }
        });
        // I tried to implement TryFrom, but the compiler complained of overlapping impls
        //   and I gave up before figuring out why
        let mut builder = MultiConfigCsvBuilder::new();
        for section in &sections {
            builder.push(section)?;
        }
        let cfg = builder.build()?;

        // this forces additional consistency checking of [libraries] section
        let _ = cfg.to_multi_graph("dummy_sample", "dummy_desc", None)?;

        Ok(cfg)
    }

    /// Produce the multi graph for this config.
    /// If this config does not have an explicit probe barcode pairing and one
    /// is provided, it will be used to generate the tag IDs.
    /// The explicit pairings should have the form target_bc:mapped_bc, such that
    /// the value is the barcode ID that will be mapped into the key.
    pub fn to_multi_graph(
        &self,
        sample_id: &str,
        sample_desc: &str,
        detected_probe_barcode_pairing: Option<&TxHashMap<BarcodeId, BarcodeId>>,
    ) -> Result<CrMultiGraph> {
        use cr_types::types::{CrMultiGraphBuilder, Fingerprint, GemWell};
        let mut builder = CrMultiGraphBuilder::new();
        for lib in &self.libraries.0 {
            let physical_library_id = lib.physical_library_id().to_string();
            let fastq = lib.to_fastq_def()?;
            let gem_well = lib.gem_well().0.into();
            builder.push_library(physical_library_id, fastq, lib.library_type(), gem_well)?;
        }
        if let Some(samples) = &self.samples {
            // Generate a reverse lookup for translated barcodes.
            let translated_source_map = {
                let explicit_mapping: TxHashMap<_, Vec<_>> = samples
                    .iter_translated_probe_barcodes()
                    .map(|(target, sources)| (target, sources.collect()))
                    .collect();
                // If we have a detected pairing and there is no explicit mapping,
                // use the detected pairing.
                let have_explicit_mapping = explicit_mapping
                    .values()
                    .any(|translated_bcs| !translated_bcs.is_empty());
                match (have_explicit_mapping, detected_probe_barcode_pairing) {
                    (true, _) | (false, None) => explicit_mapping,
                    (false, Some(detected)) => detected
                        .iter()
                        .map(|(target_bc, mapped_bc)| (*target_bc, vec![*mapped_bc]))
                        .collect(),
                }
            };
            for sample in &samples.0 {
                if let Some(sample_barcode_ids) =
                    sample.sample_barcode_ids(ProbeBarcodeIterationMode::Mapped)
                {
                    for sample_barcode_id in sample_barcode_ids {
                        builder.push_fingerprint(
                            sample.sample_id.clone(),
                            sample.description.clone(),
                            Fingerprint::tagged(
                                GemWell::default(),
                                sample_barcode_id.to_string(),
                                translated_source_map
                                    .get(sample_barcode_id)
                                    .map(|bcs| bcs.iter().map(ToString::to_string).collect())
                                    .unwrap_or_default(),
                                sample.barcode_multiplexing_type(),
                            ),
                        )?;
                    }
                } else {
                    builder.push_fingerprint(
                        sample.sample_id.clone(),
                        sample.description.clone(),
                        Fingerprint::untagged(GemWell::default()),
                    )?;
                };
            }
        } else {
            builder.push_fingerprint(
                (*sample_id).to_string(),
                (*sample_desc).to_string(),
                Fingerprint::untagged(GemWell::default()),
            )?;
        };

        Ok(builder.build())
    }

    /// Return all sample barcode IDs that we expect to appear in the analysis.
    ///
    /// Control probe barcode iteration using the provided option.
    pub fn sample_barcode_ids_used_in_experiment(
        &self,
        probe_barcode_iteration_mode: ProbeBarcodeIterationMode,
    ) -> TxHashSet<String> {
        if let Some(ref samps) = self.samples {
            samps
                .0
                .iter()
                .filter_map(|r| r.sample_barcode_ids(probe_barcode_iteration_mode))
                .flatten()
                .map(String::from)
                .collect()
        } else {
            TxHashSet::default()
        }
    }

    /// Return true if multiplexed using overhangs.
    pub fn is_overhang_multiplexed(&self) -> bool {
        self.samples
            .as_ref()
            .is_some_and(SamplesCsv::has_overhang_ids)
    }

    /// Return true if multiplexed using RTL probe barcodes.
    pub fn is_rtl_multiplexed(&self) -> bool {
        self.samples
            .as_ref()
            .is_some_and(SamplesCsv::is_rtl_multiplexed)
    }

    /// Return whether either [gene-expression] or [samples] is RTL.
    pub fn is_rtl(&self) -> bool {
        self.gene_expression
            .as_ref()
            .is_some_and(GeneExpressionParams::is_rtl)
            || self.is_rtl_multiplexed()
    }

    /// Return the mapping of chemistry specs, derived from a combination of
    /// the gex section and the libraries.
    pub fn chemistry_specs(&self) -> Result<ChemistrySpecs> {
        let gex_chemistry_spec = self.gene_expression.as_ref().and_then(|gex| gex.chemistry);
        let Some(specs) = self.libraries.chemistry_specs() else {
            // No per-library specs, use gex or default
            let spec = gex_chemistry_spec.unwrap_or(ChemistryParam::AutoOrRefined(
                AutoOrRefinedChemistry::Auto(AutoChemistryName::Count),
            ));
            return self
                .libraries
                .0
                .iter()
                .map(|lib| {
                    let library_type = lib.library_type();
                    match spec {
                        ChemistryParam::AutoOrRefined(chem) => Ok((library_type, chem)),
                        ChemistryParam::Set(chem_set) => Ok((
                            library_type,
                            AutoOrRefinedChemistry::Refined(
                                chem_set.chemistry_for_library_type(library_type)?,
                            ),
                        )),
                    }
                })
                .collect();
        };
        assert!(gex_chemistry_spec.is_none());
        Ok(specs)
    }
}

#[derive(Debug, Default)]
struct MultiConfigCsvBuilder {
    pub gene_expression: Option<GeneExpressionParams>,
    pub feature: Option<FeatureParams>,
    pub vdj: Option<VdjParams>,
    pub libraries: Option<LibrariesCsv>,
    pub samples: Option<SamplesCsv>,
    pub gem_wells: Option<GemWellsCsv>,
    pub antigen_specificity: Option<AntigenSpecificityCsv>,
    pub functional_map: Option<FunctionalMapCsv>,
}

macro_rules! setter {
    ($name:ident, $typ:ident) => {
        fn $name(&mut self, section: &Section<'_>) -> Result<()> {
            ensure!(
                self.$name.is_none(),
                "failed to parse CSV, duplicate [{}] at line: {}, col: {}",
                section.name.fragment(),
                section.name.location_line(),
                section.name.get_utf8_column(),
            );
            let $name = $typ::try_from(section)?;
            self.$name = Some($name);
            Ok(())
        }
    };
}

macro_rules! setter_validate_gws {
    ($name:ident, $typ:ident, $constant:expr) => {
        fn $name(&mut self, section: &Section<'_>) -> Result<()> {
            ensure!(
                self.$name.is_none(),
                "failed to parse CSV, duplicate [{}] at line: {}, col: {}",
                section.name.fragment(),
                section.name.location_line(),
                section.name.get_utf8_column(),
            );

            let libraries = self.libraries.as_ref().ok_or_else(|| {
                anyhow!(
                    "failed to parse [{}] before [{}]",
                    $constant,
                    multiconst::LIBRARIES,
                )
            })?;

            let valid_gws: TxHashSet<_> = libraries.0.iter().map(Library::gem_well).collect();
            let $name = $typ::try_from((&valid_gws, section))?;
            self.$name = Some($name);
            return Ok(());
        }
    };
}

impl MultiConfigCsvBuilder {
    fn new() -> Self {
        Self::default()
    }

    fn push(&mut self, section: &Section<'_>) -> Result<()> {
        use multiconst::{
            ANTIGEN_SPECIFICITY, FEATURE, GEM_WELLS, GENE_EXPRESSION, GEX, GWS, LIBRARIES, LIBS,
            SAMPLES, VDJ,
        };
        let name = section.name.fragment().to_ascii_lowercase();
        match name.as_str() {
            GENE_EXPRESSION | GEX => self.gene_expression(section),
            FEATURE => self.feature(section),
            VDJ => self.vdj(section),
            LIBRARIES | LIBS => self.libraries(section),
            SAMPLES => self.samples(section),
            GEM_WELLS | GWS => self.gem_wells(section),
            ANTIGEN_SPECIFICITY => self.antigen_specificity(section),
            FUNCTIONAL_MAP => self.functional_map(section),
            _ => bail!(
                "failed to parse CSV, unknown section [{}] at line: {}, col: {}",
                section.name.fragment(),
                section.name.location_line(),
                section.name.get_utf8_column(),
            ),
        }
    }

    setter!(gene_expression, GeneExpressionParams);
    setter!(feature, FeatureParams);
    setter!(vdj, VdjParams);
    setter!(libraries, LibrariesCsv);
    setter!(antigen_specificity, AntigenSpecificityCsv);
    setter!(functional_map, FunctionalMapCsv);

    setter_validate_gws!(samples, SamplesCsv, multiconst::SAMPLES);
    setter_validate_gws!(gem_wells, GemWellsCsv, multiconst::GEM_WELLS);

    fn build(self) -> Result<MultiConfigCsv> {
        let MultiConfigCsvBuilder {
            gene_expression: gene_expression_owned,
            feature,
            vdj,
            libraries,
            samples: samples_owned,
            gem_wells,
            antigen_specificity,
            functional_map,
        } = self;

        let Some(libraries) = libraries else {
            bail!(
                "failed to parse CSV, [{}] section not provided",
                multiconst::LIBRARIES
            );
        };

        let gene_expression = gene_expression_owned.as_ref();
        let samples = samples_owned.as_ref();
        let gex_section_chemistry = gene_expression.and_then(|x| x.chemistry);
        let library_chemistries = libraries.chemistry_specs();

        // Cannot mix specifying chemistry in gex section and libraries.
        ensure!(
            gex_section_chemistry.is_none() || library_chemistries.is_none(),
            "failed to parse CSV: chemistry specified in both the [{}] and [{}] sections",
            multiconst::GENE_EXPRESSION,
            multiconst::LIBRARIES
        );

        let any_source_chemistries: TxHashSet<_> = library_chemistries
            .into_iter()
            .flat_map(HashMap::into_values)
            .map(ChemistryParam::AutoOrRefined)
            .chain(gex_section_chemistry)
            .collect();

        let has_probe_barcode_ids = samples.is_some_and(SamplesCsv::has_probe_barcode_ids);
        let is_rtl = has_probe_barcode_ids
            || any_source_chemistries
                .iter()
                .any(|x| x.is_rtl() == Some(true));

        if libraries.has_gene_expression() {
            ensure!(
                gene_expression.is_some(),
                "failed to parse CSV: [{}] section omitted but gene expression libraries provided",
                multiconst::GENE_EXPRESSION,
            );
        }

        if libraries.has_feature_barcode() {
            ensure!(
                feature.as_ref().is_some_and(|x| x.reference_path.is_some()),
                "failed to parse CSV: [{}] reference is missing \
                 but feature barcode libraries provided",
                multiconst::FEATURE,
            );
        }

        if libraries.has_vdj() {
            ensure!(
                vdj.is_some(),
                "failed to parse CSV: [{}] section omitted but VDJ libraries provided",
                multiconst::VDJ,
            );
        }

        if libraries.has_multiplexing() {
            ensure!(
                samples.is_some(),
                "failed to parse CSV: [{}] section omitted \
                 but Multiplexing Capture libraries provided",
                multiconst::SAMPLES,
            );
        }

        let has_hashtag_ids = samples.is_some_and(SamplesCsv::has_hashtag_ids);
        if has_hashtag_ids {
            ensure!(
                libraries.has_gene_expression() && libraries.has_antibody_capture(),
                "failed to parse CSV: [{}] section specified {} parameter \
                 which require both Gene Expression and Antibody Capture libraries.",
                multiconst::SAMPLES,
                samplesconst::HASHTAG_IDS,
            );
        }

        if antigen_specificity.is_some() {
            ensure!(
                libraries.has_antigen_capture(),
                "failed to parse CSV: [{}] section is provided but no Antigen Capture libraries",
                multiconst::ANTIGEN_SPECIFICITY,
            );
        }

        if functional_map.is_some() {
            ensure!(
                libraries.has_feature_barcode(),
                "failed to parse CSV: [{}] section is provided \
                 but no feature barcode libraries provided",
                multiconst::FUNCTIONAL_MAP,
            );
        }

        if libraries.has_antigen_capture() {
            let invalid_parameter = gene_expression
                .unwrap()
                .invalid_parameter_with_antigen_capture();
            ensure!(
                invalid_parameter.is_none(),
                "failed to parse CSV: [{}] section specified {} parameter but [{}] section \
                 includes Antigen Capture libraries. Invalid parameter in this context",
                multiconst::GENE_EXPRESSION,
                invalid_parameter.unwrap(),
                multiconst::LIBRARIES,
            );
        }

        if any_source_chemistries.iter().any(|chem| {
            *chem
                == ChemistryParam::AutoOrRefined(AutoOrRefinedChemistry::Refined(
                    ChemistryName::SFRP,
                ))
        }) {
            ensure!(
                samples.is_none(),
                "failed to parse CSV: specified `{}` as the `chemistry` but a \
                 [{}] section is also provided",
                ChemistryName::SFRP,
                multiconst::SAMPLES,
            );
        }

        if let Some(mfrp_chem) = any_source_chemistries.iter().find(|spec| spec.is_mfrp()) {
            ensure!(
                has_probe_barcode_ids,
                "failed to parse CSV: specified `{}` as the `chemistry` but a \
                 [{}] section with {} column is not provided",
                mfrp_chem,
                multiconst::SAMPLES,
                samplesconst::PROBE_BARCODE_IDS,
            );

            ensure!(
                !gene_expression.is_some_and(GeneExpressionParams::has_expect_or_force_cells),
                "failed to parse CSV: specified `{}` as the `chemistry` and \
                 cell calling parameters in {} section. For multiplex Flex libraries \
                 `expect-cells` or `force-cells' parameter is valid only in [{}] section.",
                mfrp_chem,
                multiconst::GENE_EXPRESSION,
                multiconst::SAMPLES,
            );
        }

        if is_rtl {
            if libraries.has_gene_expression() {
                ensure!(
                    gene_expression.is_some_and(GeneExpressionParams::has_probe_set),
                    "failed to parse CSV: [{}] section is missing `probe-set` a required \
                     parameter for Flex",
                    multiconst::GENE_EXPRESSION,
                );
            }

            // Disallow setting include-introns for RTL chemistries.
            ensure!(
                gene_expression.map_or(true, |x| x.include_introns
                    == multiconst::DEFAULT_INCLUDE_INTRONS),
                ERROR_INCLUDE_INTRONS_WITH_RTL,
            );
        }

        if has_probe_barcode_ids {
            // NOTE: non-multiplexed chemistry + having probe_barcode_ids column(invalid) passes here but would fail in checks above
            ensure!(
                any_source_chemistries
                    .iter()
                    .all(|chem| chem.is_rtl() != Some(false)),
                "failed to parse CSV: [{}] section manually specifies a \
                 non-Flex chemistry but [{}] section has a {} column. \
                 The {} column may only be specified with Flex chemistries",
                multiconst::GENE_EXPRESSION,
                multiconst::SAMPLES,
                samplesconst::PROBE_BARCODE_IDS,
                samplesconst::PROBE_BARCODE_IDS,
            );
        }

        ensure!(
            !(gene_expression.is_some_and(GeneExpressionParams::has_expect_or_force_cells)
                && has_probe_barcode_ids),
            "failed to parse CSV: [{}] section has a {} column indicating multiplex \
             Flex chemistry, however a cell calling parameter is specified \
             in [{}] section. The parameters `expect-cells` or `force-cells' are valid \
             in [{}] section for singleplex Flex and in [{}] section for multiplex Flex.",
            multiconst::SAMPLES,
            samplesconst::PROBE_BARCODE_IDS,
            multiconst::GENE_EXPRESSION,
            multiconst::GENE_EXPRESSION,
            multiconst::SAMPLES,
        );

        ensure!(
            !(gene_expression.is_some_and(GeneExpressionParams::has_expect_or_force_cells)
                && samples.is_some_and(SamplesCsv::has_expect_or_force_cells)),
            "failed to parse CSV: 'expect-cells' or 'force-cells' parameters are specified \
             in both [{}] and [{}] sections",
            multiconst::GENE_EXPRESSION,
            multiconst::SAMPLES,
        );

        ensure!(
            !(gene_expression.is_some_and(GeneExpressionParams::has_global_minimum_umis)
                && samples.is_some_and(SamplesCsv::has_global_minimum_umis)),
            "failed to parse CSV: 'global-minimum-umis' parameter is specified \
             in both [{}] and [{}] sections",
            multiconst::GENE_EXPRESSION,
            multiconst::SAMPLES,
        );

        Ok(MultiConfigCsv {
            gene_expression: gene_expression_owned,
            feature,
            vdj,
            libraries,
            samples: samples_owned,
            gem_wells,
            antigen_specificity,
            functional_map,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::scsv::XtraData;
    use super::{MultiConfigCsv, ProbeBarcodeIterationMode, SampleRow};
    use crate::config::{ChemistryParam, ChemistrySet};
    use anyhow::Result;
    use barcode::whitelist::BarcodeId;
    use cr_types::chemistry::{AutoOrRefinedChemistry, ChemistryName};
    use cr_types::reference::probe_set_reference::TargetSetFile;
    use cr_types::{BarcodeMultiplexingType, CellLevel, Fingerprint, GemWell, LibraryType};
    use insta::{assert_debug_snapshot, assert_json_snapshot, assert_snapshot};
    use itertools::Itertools;
    use metric::TxHashMap;
    use std::collections::HashMap;
    use std::path::{Path, PathBuf};

    // initialize insta test harness
    #[ctor::ctor]
    fn init() {
        // this ensures insta knows where to find its snap tests
        let cwd = std::env::current_dir().unwrap();
        let workspace_root = cwd.parent().unwrap();
        std::env::set_var("INSTA_WORKSPACE_ROOT", workspace_root);
    }

    /// Assert that parsing the provided CSV contents fails with a matching error.
    fn assert_parse_fails(contents: &str, expected_error: &str) {
        assert_eq!(
            expected_error,
            MultiConfigCsv::from_reader(contents.as_bytes(), XtraData::new("tests"))
                .unwrap_err()
                .to_string(),
        );
    }

    #[test]
    fn load_simple_external_cmos() -> Result<()> {
        let csv = r#"
[gene-expression]
ref,/path/to/gex/ref
chemistry,
create-bam,true

[feature]
reference,/path/to/feature/ref

[vdj]
ref,/path/to/vdj/ref
inner-enrichment-primers,

[libraries]
fastq_id,fastqs,physical_library_id,feature_types,gem_well,subsample_rate,
mygex,/path/to/fastqs,a,Gene Expression,,
mycmo,/path/to/fastqs,1,Multiplexing Capture,,,

[samples]
sample_id,cmo_ids,gem_wells,description
whoami,cmo_2,1-1,hidad!,

[gem-wells]
gem_well,force_cells,expect_cells,vdj_force_cells,
1,,,
"#;
        let xtra = XtraData::new("tests::load_simple_external_cmos");
        let exp = MultiConfigCsv::from_reader(csv.as_bytes(), xtra)?;
        println!("{exp:?}");
        Ok(())
    }

    #[test]
    fn load_simple_external_probe_bcs() -> Result<()> {
        let csv = r#"
[gene-expression]
ref,/path/to/gex/ref
probe-set,/path/to/probe_set
create-bam,true

[libraries]
fastq_id,fastqs,physical_library_id,feature_types,gem_well,subsample_rate,
mygex,/path/to/fastqs,a,Gene Expression,,
mycmo,/path/to/fastqs,1,Multiplexing Capture,,,

[samples]
sample_id,probe_barcode_ids,gem_wells,description,force_cells,expect_cells,
whoami,BC1,,hidad!,,100
"#;
        let xtra = XtraData::new("tests::load_simple_external_probe_bcs");
        let exp = MultiConfigCsv::from_reader(csv.as_bytes(), xtra)?;
        println!("{exp:?}");
        Ok(())
    }

    #[test]
    fn load_include_introns_default() -> Result<()> {
        let csv = r#"
[gene-expression]
ref,/path/to/gex/ref
probe-set,/path/to/probe_set
create-bam,true

[libraries]
fastq_id,fastqs,physical_library_id,feature_types,gem_well,subsample_rate,
mygex,/path/to/fastqs,a,Gene Expression,,
mycmo,/path/to/fastqs,1,Multiplexing Capture,,,

[samples]
sample_id,probe_barcode_ids,gem_wells,description,force_cells,expect_cells,
whoami,BC1,,hidad!,,100
"#;
        let xtra = XtraData::new("tests::load_introns_for_rtl");
        let exp = MultiConfigCsv::from_reader(csv.as_bytes(), xtra)?;
        assert!(exp.gene_expression.unwrap().include_introns);

        let csv_with_introns = r#"
[gene-expression]
ref,/path/to/gex/ref
chemistry,
create-bam,true

[feature]
reference,/path/to/feature/ref

[vdj]
ref,/path/to/vdj/ref
inner-enrichment-primers,

[libraries]
fastq_id,fastqs,physical_library_id,feature_types,gem_well,subsample_rate,
mygex,/path/to/fastqs,a,Gene Expression,,
mycmo,/path/to/fastqs,1,Multiplexing Capture,,,

[samples]
sample_id,cmo_ids,gem_wells,description
whoami,cmo_2,1-1,hidad!,

[gem-wells]
gem_well,force_cells,expect_cells,vdj_force_cells,
1,,,
"#;
        let xtra2 = XtraData::new("tests::include_introns");
        let exp2 = MultiConfigCsv::from_reader(csv_with_introns.as_bytes(), xtra2)?;
        assert!(exp2.gene_expression.unwrap().include_introns);
        let csv_with_no_introns = r#"
[gene-expression]
ref,/path/to/gex/ref
include-introns,false
create-bam,true

[feature]
reference,/path/to/feature/ref

[vdj]
ref,/path/to/vdj/ref
inner-enrichment-primers,

[libraries]
fastq_id,fastqs,physical_library_id,feature_types,gem_well,subsample_rate,
mygex,/path/to/fastqs,a,Gene Expression,,
mycmo,/path/to/fastqs,1,Multiplexing Capture,,,

[samples]
sample_id,cmo_ids,gem_wells,description
whoami,cmo_2,1-1,hidad!,

[gem-wells]
gem_well,force_cells,expect_cells,vdj_force_cells,
1,,,
"#;
        let xtra3 = XtraData::new("tests::no_include_introns");
        let exp3 = MultiConfigCsv::from_reader(csv_with_no_introns.as_bytes(), xtra3)?;
        assert!(!exp3.gene_expression.unwrap().include_introns);
        Ok(())
    }

    #[test]
    fn load_simple_internal() -> Result<()> {
        let csv = r#"
[gene-expression]
reference-path,/path/to/gex/ref
r1-length,28
r2-length,51
chemistry,SC3Pv3
create-bam,true

[feature]
ref,test/feature/cmo_features.csv
r1-length,28
r2-length,51

[vdj]
ref,/path/to/vdj/ref
inner-enrichment-primers,/path/to/primers
r1-length,28
r2-length,51

[libraries]
fastq_path,sample_indices,lanes,physical_library_id,feature_types,gem_well,subsample_rate,
/path/to/fastqs,SI-001,1-3|5,whatsit_gex,Gene Expression,1,0.5,
/path/to/fastqs,SI-002,4,whatsit_cmo,Multiplexing Capture,1,1,

[samples]
sample_id,cmo_ids,gem_wells,description
abc,cmo_1-3|cmo_5,1,"hi mom",
def,cmo_4,1,hi dad!

[gem-wells]
gem_well,force_cells,expect_cells,vdj_force_cells,
1,9,99,999
"#;
        let xtra = XtraData::new("tests::load_simple_internal");
        let exp = MultiConfigCsv::from_reader(csv.as_bytes(), xtra)?;
        println!("{exp:?}");
        Ok(())
    }

    #[test]
    fn invalid_sample_id_chars() {
        let csv = r#"
[gene-expression]
ref,GRCh38-2020-A-chr21
create-bam,true

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types
bamtofastq,fastqs/cellranger/multi/VDJ_GEX_small_multi/small_gex_fastqs_chr21_new,any,gex_1,gene expression
pbmc_1k_protein_v3_antibody,../dui_tests/test_resources/cellranger-count/pbmc_1k_protein_v3_antibody,any,cmo1,Multiplexing Capture

[samples]
sample_id,cmo_ids,gem_wells,description
IAmABad*SampleId,CMO1,1,some cells
"#;
        let xtra = XtraData::new("tests::invalid_sample_id_chars");
        let res = MultiConfigCsv::from_reader(csv.as_bytes(), xtra);
        assert_snapshot!(res.unwrap_err());
    }

    #[test]
    fn invalid_sample_id_chars_2() {
        let csv = r#"
[gene-expression]
ref,GRCh38-2020-A-chr21
create-bam,true

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types
bamtofastq,fastqs/cellranger/multi/VDJ_GEX_small_multi/small_gex_fastqs_chr21_new,any,gex_1,gene expression
pbmc_1k_protein_v3_antibody,../dui_tests/test_resources/cellranger-count/pbmc_1k_protein_v3_antibody,any,cmo1,Multiplexing Capture

[samples]
sample_id,cmo_ids,gem_wells,description
IAmABad#SampleId,CMO1,1,some cells
"#;
        let xtra = XtraData::new("tests::invalid_sample_id_chars_2");
        let res = MultiConfigCsv::from_reader(csv.as_bytes(), xtra);
        assert_snapshot!(res.unwrap_err());
    }

    #[test]
    fn too_long_sample_id() {
        let csv = r#"
[gene-expression]
ref,GRCh38-2020-A-chr21
create-bam,true

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types
bamtofastq,fastqs/cellranger/multi/VDJ_GEX_small_multi/small_gex_fastqs_chr21_new,any,gex_1,gene expression
pbmc_1k_protein_v3_antibody,../dui_tests/test_resources/cellranger-count/pbmc_1k_protein_v3_antibody,any,cmo1,Multiplexing Capture

[samples]
sample_id,cmo_ids,gem_wells,description
IAmToooooooooooooooooooooooooooooooooooooooooooooooooooooooooLong,CMO1,1,some cells
"#;
        let xtra = XtraData::new("tests::too_long_sample_id");
        let res = MultiConfigCsv::from_reader(csv.as_bytes(), xtra);
        assert_snapshot!(res.unwrap_err());
    }

    #[test]
    fn load_trailing_whitespace() -> Result<()> {
        let csv = r#"
[gene-expression]
ref,GRCh38-2020-A-chr21
chemistry,SC5P-R2
create-bam,true

[feature]
ref,test/feature/cmo_features.csv

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types
bamtofastq,fastqs/cellranger/multi/VDJ_GEX_small_multi/small_gex_fastqs_chr21_new,any,gex_1,gene expression
pbmc_1k_protein_v3_antibody,../dui_tests/test_resources/cellranger-count/pbmc_1k_protein_v3_antibody,any,cmo1,Multiplexing Capture

[samples]
sample_id,cmo_ids,gem_wells,description
HUMAN_T,CMO1,1,some cells
"#;
        let xtra = XtraData::new("tests::load_trailing_whitespace");
        let exp = MultiConfigCsv::from_reader(csv.as_bytes(), xtra)?;
        println!("{exp:?}");
        Ok(())
    }

    #[test]
    fn test_check_library_compatibility() -> Result<()> {
        let csv = r#"
[gene-expression]
check_library_compatibility,false
reference,GRCh38-2020-A-chr21
create-bam,true

[feature]
ref,test/feature/cmo_features.csv

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types
bamtofastq,fastqs/cellranger/multi/VDJ_GEX_small_multi/small_gex_fastqs_chr21_new,any,gex_1,gene expression
"#;
        let xtra = XtraData::new("tests::test_check_library_compatibility");
        let exp = MultiConfigCsv::from_reader(csv.as_bytes(), xtra)?;
        println!("{exp:?}");
        assert!(!exp.gene_expression.unwrap().check_library_compatibility);
        let csv2 = r#"
[gene-expression]
reference,GRCh38-2020-A-chr21
create-bam,true

[feature]
ref,test/feature/cmo_features.csv

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types
bamtofastq,fastqs/cellranger/multi/VDJ_GEX_small_multi/small_gex_fastqs_chr21_new,any,gex_1,gene expression

"#;
        let xtra2 = XtraData::new("tests::test_check_library_compatibility");
        let exp2 = MultiConfigCsv::from_reader(csv2.as_bytes(), xtra2)?;
        assert!(exp2.gene_expression.unwrap().check_library_compatibility);
        let csv3 = r#"
[gene-expression]
ref,GRCh38-2020-A-chr21
check-library-compatibility,false
create-bam,true

[feature]
ref,test/feature/cmo_features.csv

[vdj]
ref,vdj/vdj_GRCh38_alts_ensembl-4.0.0

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
vdj_t,fastqs/cellranger/multi/vdj_t_tiny,any,lib2,vdj-t,
vdj_ig,fastqs/cellranger/multi/vdj_ig_tiny,any,lib4,vdj-b,
FB_fastqs_id,path/to/FB_fastqs,1|2,any,Antibody Capture,
"#;
        let xtra3 = XtraData::new("tests::test_check_library_compatibility");
        let exp3 = MultiConfigCsv::from_reader(csv3.as_bytes(), xtra3)?;
        assert!(!exp3.gene_expression.unwrap().check_library_compatibility);
        Ok(())
    }

    #[test]
    fn no_physical_library_id() -> Result<()> {
        let csv = r#"
[gene-expression]
reference-path,/path/to/gex/ref
r1-length,28
r2-length,51
chemistry,SC3Pv3
create-bam,true

[feature]
ref,/path/to/feature/ref
r1-length,28
r2-length,51

[vdj]
ref,/path/to/vdj/ref
inner-enrichment-primers,/path/to/primers
r1-length,28
r2-length,51

[libraries]
fastq_path,sample_indices,lanes,feature_types,gem_well,subsample_rate,
/path/to/fastqs,SI-001,1-3|5,Gene Expression,1,0.5,
/path/to/fastqs,SI-002,4,Multiplexing Capture,1,1,

[samples]
sample_id,cmo_ids
sample,CMO1
"#;
        let xtra = XtraData::new("tests::no_physical_library_id");
        let exp = MultiConfigCsv::from_reader(csv.as_bytes(), xtra)?;
        println!("{exp:?}");
        Ok(())
    }

    #[test]
    fn test_vdj_typo() -> Result<()> {
        let csv = r#"
[gene-expression]
ref,GRCh38-2020-A-chr21
chemistry,SC5P-R2
create-bam,true

[  vdj       ]
ref,vdj/vdj_GRCh38_alts_ensembl-4.0.0

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
bamtofastq,fastqs/cellranger/multi/targeted_tiny,any,lib1,gene expression,
vdj_t,fastqs/cellranger/multi/vdj_t_tiny,any,lib2,vdj-t,
vdj_ig,fastqs/cellranger/multi/vdj_ig_tiny,any,lib4,vdj-b,
"#;
        let xtra = XtraData::new("test::test_vdj_typo");
        let exp = MultiConfigCsv::from_reader(csv.as_bytes(), xtra)?;
        println!("{exp:?}");
        Ok(())
    }

    /// Ensure that extra commas in the multi CSV are handled nicely
    #[test]
    fn test_extra_commas() -> Result<()> {
        let csv = r#"
# this is a comment,,,
[gene-expression],,,
reference,/home/labs/bioservices/services/expression_references/refdata-gex-mm10-2020-A,,
chemistry,auto,,
include-introns,FALSE,,
no-secondary,FALSE,,
create-bam,TRUE,,
,,,
[feature],,,
reference,/home/labs/abramson/Collaboration/Shir_cite-seq/Shir_feature_ref.csv,,
,,,
,,,
[libraries],,,
fastq_id,fastqs,lanes,feature_types
A_36_cDNA,/home/labs/abramson/Collaboration/Shir_cite-seq/210106_A00929_0240_AHY77NDRXX/HY77NDRXX/outs/fastq_path/,any,gene expression
A_36_Ab,/home/labs/abramson/Collaboration/Shir_cite-seq/210106_A00929_0240_AHY77NDRXX/HY77NDRXX/outs/fastq_path/,any,antibody capture
"#;

        let xtra = XtraData::new("test::extra_commas");
        let exp = MultiConfigCsv::from_reader(csv.as_bytes(), xtra)?;
        println!("{exp:?}");
        Ok(())
    }

    #[test]
    fn test_autogen_phys_lib_id_two_vdj_libs() {
        let csv = r#"
[vdj]
ref,vdj/vdj_GRCh38_alts_ensembl-4.0.0

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
vdj_t,fastqs/cellranger/multi/vdj_t_tiny,any,,vdj,
vdj_ig,fastqs/cellranger/multi/vdj_ig_tiny,any,,vdj,
"#;
        let xtra = XtraData::new("test::test_autogen_phys_lib_id_two_vdj_libs");
        let res = MultiConfigCsv::from_reader(csv.as_bytes(), xtra);
        assert_snapshot!(res.unwrap_err());
    }

    /// Ensure that antigen capture is accompanied with a VDJ library
    #[test]
    fn test_antigen_requires_vdj() -> Result<()> {
        let csv = r#"
    [gene-expression]
    ref,mm10-2020-A-chr19
    create-bam,true

    [feature]
    ref,cellranger/multi/feature_refs/20211122_v1.1.csv

    [libraries]
    fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
    tiny_gex,fastqs/cellranger/multi/1245140_gex_vdj_beam_ab/gex,any,gex,gene expression,0.5
    tiny_an_b,fastqs/cellranger/multi/1245140_gex_vdj_beam_ab/an_b,any,ab,antigen capture,0.5
    "#;

        let xtra = XtraData::new("test::antigen_required_vdj");
        let res = MultiConfigCsv::from_reader(csv.as_bytes(), xtra);
        assert_snapshot!(res.unwrap_err());
        Ok(())
    }

    #[test]
    fn test_antigen_specificity_requires_ag() -> Result<()> {
        let csv = r#"
    [gene-expression]
    ref,mm10-2020-A-chr19
    create-bam,true

    [feature]
    ref,cellranger/multi/feature_refs/20211122_v1.1.csv

    [antigen-specificity]
    control_id
    feature_01

    [libraries]
    fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
    tiny_gex,fastqs/cellranger/multi/1245140_gex_vdj_beam_ab/gex,any,gex,gene expression,0.5
    "#;

        let xtra = XtraData::new("test::antigen_specificity_requires_ag");
        let res = MultiConfigCsv::from_reader(csv.as_bytes(), xtra);
        assert_snapshot!(res.unwrap_err());
        Ok(())
    }

    #[test]
    fn test_antigen_specificity_multiple_control_ids() -> Result<()> {
        let csv = r#"
    [gene-expression]
    ref,mm10-2020-A-chr19
    create-bam,true

    [feature]
    ref,cellranger/multi/feature_refs/20211122_v1.1.csv

    [antigen-specificity]
    control_id, mhc_allele
    feature_01,
    feature_02,

    [libraries]
    fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
    tiny_gex,fastqs/cellranger/multi/1245140_gex_vdj_beam_ab/gex,any,gex,gene expression,0.5
    "#;

        let xtra = XtraData::new("test::antigen_specificity_multiple_control_ids");
        let res = MultiConfigCsv::from_reader(csv.as_bytes(), xtra);
        assert_snapshot!(res.unwrap_err());
        Ok(())
    }

    #[test]
    fn test_antigen_specificity_missing_allele() -> Result<()> {
        let csv = r#"
    [gene-expression]
    ref,mm10-2020-A-chr19
    create-bam,true

    [feature]
    ref,cellranger/multi/feature_refs/20211122_v1.1.csv

    [antigen-specificity]
    control_id, mhc_allele
    feature_01, a1
    feature_02,

    [libraries]
    fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
    tiny_gex,fastqs/cellranger/multi/1245140_gex_vdj_beam_ab/gex,any,gex,gene expression,0.5
    "#;

        let xtra = XtraData::new("test::antigen_specificity_missing_allele");
        let res = MultiConfigCsv::from_reader(csv.as_bytes(), xtra);
        assert_snapshot!(res.unwrap_err());
        Ok(())
    }

    #[test]
    fn test_invalid_mhc_allele_character() -> Result<()> {
        let csv = r#"
    [gene-expression]
    ref,mm10-2020-A-chr19
    create-bam,true

    [feature]
    ref,cellranger/multi/feature_refs/20211122_v1.1.csv

    [vdj]
    ref,vdj/vdj_GRCh38_alts_ensembl-4.0.0

    [antigen-specificity]
    control_id, mhc_allele
    feature_01, :??@

    [libraries]
    fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
    tiny_gex,fastqs/cellranger/multi/1245140_gex_vdj_beam_ab/gex,any,gex,gene expression,0.5
    tiny_vdj,fastqs/cellranger/multi/1245140_gex_vdj_beam_ab/vdj,any,vdj,vdj,0.5
    tiny_an_b,fastqs/cellranger/multi/1245140_gex_vdj_beam_ab/an_b,any,ab,antigen capture,0.5
    "#;

        let xtra = XtraData::new("test::invalid_mhc_allele_character");
        let res = MultiConfigCsv::from_reader(csv.as_bytes(), xtra);
        assert_snapshot!(res.unwrap_err());
        Ok(())
    }

    #[test]
    fn test_functional_map_requires_feature() -> Result<()> {
        let csv = r#"
    [gene-expression]
    ref,mm10-2020-A-chr19
    create-bam,true

    [feature-functional-map]
    functional_name,feature_ids
    flu,feature01|feature02
    cmv,feature03

    [libraries]
    fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
    tiny_gex,fastqs/cellranger/multi/1245140_gex_vdj_beam_ab/gex,any,gex,gene expression,0.5
    "#;

        let xtra = XtraData::new("test::test_functional_map_requires_feature");
        let res = MultiConfigCsv::from_reader(csv.as_bytes(), xtra);
        assert_snapshot!(res.unwrap_err());
        Ok(())
    }

    #[test]
    fn test_invalid_functional_map_value() -> Result<()> {
        // Duplicate feature ids
        let csv = r#"
    [gene-expression]
    ref,mm10-2020-A-chr19
    create-bam,true

    [feature]
    ref,cellranger/multi/feature_refs/20211122_v1.1.csv

    [vdj]
    ref,vdj/vdj_GRCh38_alts_ensembl-4.0.0

    [feature-functional-map]
    functional_name,feature_ids
    flu,feature01|feature02
    cmv,feature01

    [libraries]
    fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
    tiny_gex,fastqs/cellranger/multi/1245140_gex_vdj_beam_ab/gex,any,gex,gene expression,0.5
    tiny_vdj,fastqs/cellranger/multi/1245140_gex_vdj_beam_ab/vdj,any,vdj,vdj,0.5
    tiny_an_b,fastqs/cellranger/multi/1245140_gex_vdj_beam_ab/an_b,any,ab,antigen capture,0.5
    "#;

        let xtra = XtraData::new("test::test_invalid_functional_map_value");
        let res = MultiConfigCsv::from_reader(csv.as_bytes(), xtra);
        assert_snapshot!(res.unwrap_err());

        // Duplicate functional name
        let csv = r#"
    [gene-expression]
    ref,mm10-2020-A-chr19
    create-bam,true

    [feature]
    ref,cellranger/multi/feature_refs/20211122_v1.1.csv

    [vdj]
    ref,vdj/vdj_GRCh38_alts_ensembl-4.0.0

    [feature-functional-map]
    functional_name,feature_ids
    flu,feature01|feature02
    flu,feature03

    [libraries]
    fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
    tiny_gex,fastqs/cellranger/multi/1245140_gex_vdj_beam_ab/gex,any,gex,gene expression,0.5
    tiny_vdj,fastqs/cellranger/multi/1245140_gex_vdj_beam_ab/vdj,any,vdj,vdj,0.5
    tiny_an_b,fastqs/cellranger/multi/1245140_gex_vdj_beam_ab/an_b,any,ab,antigen capture,0.5
    "#;

        let xtra = XtraData::new("test::test_invalid_functional_map_value");
        let res = MultiConfigCsv::from_reader(csv.as_bytes(), xtra);
        assert_snapshot!(res.unwrap_err());

        Ok(())
    }

    #[test]
    fn test_blank_lines() -> Result<()> {
        let csv = r#"
    [gene-expression]
    ref,mm10-2020-A-chr19
    create-bam,true

    probe-set,/path/to/probe_set

    [libraries]
    fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
    mygex,/path/to/fastqs,any,gex,gene expression,0.5

    mycmo,/path/to/fastqs,any,cmo,Multiplexing Capture,

    [samples]
    sample_id,cmo_ids
    sample1,1

    sample2,2
    "#;

        let xtra = XtraData::new("test::blank_lines");
        let res = MultiConfigCsv::from_reader(csv.as_bytes(), xtra)?;

        let gex = res
            .gene_expression
            .expect("gene expression section not present");
        assert_eq!(
            gex.probe_set,
            vec![TargetSetFile::from("/path/to/probe_set")]
        );
        assert_eq!(gex.reference_path, Some(PathBuf::from("mm10-2020-A-chr19")));

        let lib_count = res.libraries.0.len();
        assert_eq!(lib_count, 2, "expected 2 libraries, found {lib_count}");

        let sample_count = res.samples.expect("samples section not present").0.len();
        assert_eq!(sample_count, 2);

        Ok(())
    }

    #[test]
    fn test_multiple_probe_sets() -> Result<()> {
        let csv = r#"
    [gene-expression]
    ref,mm10-2020-A-chr19
    create-bam,true
    probe-set,/path/to/probe_set1
    probe-set,/path/to/probe_set2

    [libraries]
    fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
    mygex,/path/to/fastqs,any,gex,gene expression,0.5

    mycmo,/path/to/fastqs,any,cmo,Multiplexing Capture,

    [samples]
    sample_id,cmo_ids
    sample1,1

    sample2,2
    "#;

        let xtra = XtraData::new("test::blank_lines");
        let res = MultiConfigCsv::from_reader(csv.as_bytes(), xtra)?;

        let gex = res
            .gene_expression
            .expect("gene expression section not present");
        assert_eq!(
            gex.probe_set,
            vec![
                TargetSetFile::from("/path/to/probe_set1"),
                TargetSetFile::from("/path/to/probe_set2"),
            ]
        );
        assert_eq!(
            gex.reference_path.as_deref(),
            Some(Path::new("mm10-2020-A-chr19"))
        );

        let lib_count = res.libraries.0.len();
        assert_eq!(lib_count, 2, "expected 2 libraries, found {lib_count}");

        let sample_count = res.samples.expect("samples section not present").0.len();
        assert_eq!(sample_count, 2);

        Ok(())
    }

    #[test]
    fn test_emptydrops_per_sample() -> Result<()> {
        let csv = r#"
    [gene-expression]
    ref,mm10-2020-A-chr19
    create-bam,true

    [libraries]
    fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
    mygex,/path/to/fastqs,any,gex,gene expression,0.5
    mycmo,/path/to/fastqs,any,cmo,Multiplexing Capture,

    [samples]
    sample_id,cmo_ids,emptydrops_minimum_umis
    sample1,1,100
    sample2,2,200
    "#;

        let xtra = XtraData::new("test::");
        let res = MultiConfigCsv::from_reader(csv.as_bytes(), xtra)?;
        res.samples
            .expect("samples section not present")
            .0
            .iter()
            .for_each(|sample| match sample.sample_id.as_str() {
                "sample1" => assert_eq!(sample.emptydrops_minimum_umis, Some(100)),
                "sample2" => assert_eq!(sample.emptydrops_minimum_umis, Some(200)),
                _ => panic!("unexpected sample found: {}", sample.sample_id),
            });
        Ok(())
    }

    #[test]
    fn test_emptydrops_present() -> Result<()> {
        let csv = r#"
    [gene-expression]
    ref,mm10-2020-A-chr19
    emptydrops-minimum-umis,75
    create-bam,true

    [libraries]
    fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
    tiny_gex,fastqs/cellranger/multi/1245140_gex_vdj_beam_ab/gex,any,gex,gene expression,0.5
    "#;

        let xtra = XtraData::new("test::");
        let res = MultiConfigCsv::from_reader(csv.as_bytes(), xtra)?;
        assert_eq!(
            res.gene_expression
                .expect("gene expression missing")
                .emptydrops_minimum_umis,
            Some(75),
        );
        Ok(())
    }

    fn create_bc_pairing(
        pairs: &[(&'static str, &'static str)],
    ) -> TxHashMap<BarcodeId, BarcodeId> {
        pairs
            .iter()
            .map(|&(x, y)| (BarcodeId::pack(x), BarcodeId::pack(y)))
            .collect()
    }

    #[test]
    fn test_unpack_probe_barcode_ids() -> Result<()> {
        let csv = r#"
    [gene-expression]
    ref,mm10-2020-A-chr19
    probe-set,/path/to/probe_set
    create-bam,true

    [libraries]
    fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
    mygex,/path/to/fastqs,any,gex,gene expression,0.5

    [samples]
    sample_id,probe_barcode_ids
    sample1,BC1|BC2+BC3
    sample2,BC4|BC5+BC6+BC7
    "#;

        let xtra = XtraData::new("test::");
        let res = MultiConfigCsv::from_reader(csv.as_bytes(), xtra)?;
        let samples = res.samples.as_ref().expect("samples section not present");

        let mapping = samples.get_translated_probe_barcodes();
        let expected = create_bc_pairing(&[("BC3", "BC2"), ("BC6", "BC5"), ("BC7", "BC5")]);
        assert_eq!(mapping, expected);

        let expect_unpack_bcs = |row: &SampleRow, include_mapped, expected| {
            assert_eq!(
                row.sample_barcode_ids(include_mapped)
                    .unwrap()
                    .sorted()
                    .collect::<Vec<_>>(),
                expected
            );
        };
        expect_unpack_bcs(
            &samples.0[0],
            ProbeBarcodeIterationMode::All,
            vec!["BC1", "BC2", "BC3"],
        );
        expect_unpack_bcs(
            &samples.0[0],
            ProbeBarcodeIterationMode::Mapped,
            vec!["BC1", "BC2"],
        );
        expect_unpack_bcs(
            &samples.0[1],
            ProbeBarcodeIterationMode::All,
            vec!["BC4", "BC5", "BC6", "BC7"],
        );
        expect_unpack_bcs(
            &samples.0[1],
            ProbeBarcodeIterationMode::Mapped,
            vec!["BC4", "BC5"],
        );

        // Test that we generate a multi graph with expected content.
        assert_json_snapshot!(
            "expected_probe_barcode_pairing_graph",
            res.to_multi_graph("test", "test description", None)?
        );

        // Test that we ignore a detected probe barcode pairing if we have an
        // explicit one.
        let different_pairing = create_bc_pairing(&[("BC1", "BC10"), ("BC2", "BC7")]);
        assert_json_snapshot!(
            "expected_probe_barcode_pairing_graph",
            res.to_multi_graph("test", "test description", Some(&different_pairing))?
        );
        Ok(())
    }

    #[test]
    fn test_no_translated_probe_barcodes() -> Result<()> {
        let csv = r#"
    [gene-expression]
    ref,mm10-2020-A-chr19
    probe-set,/path/to/probe_set
    create-bam,true

    [libraries]
    fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
    mygex,/path/to/fastqs,any,gex,gene expression,0.5

    [samples]
    sample_id,probe_barcode_ids
    sample1,BC1|BC2
    sample2,BC3
    "#;

        let xtra = XtraData::new("test::");
        let res = MultiConfigCsv::from_reader(csv.as_bytes(), xtra)?;
        assert!(res
            .samples
            .expect("samples section not present")
            .get_translated_probe_barcodes()
            .is_empty());
        Ok(())
    }

    #[test]
    fn test_multi_graph_detected_probe_barcode_pairing() -> Result<()> {
        let csv = r#"
    [gene-expression]
    ref,mm10-2020-A-chr19
    probe-set,/path/to/probe_set
    create-bam,true

    [libraries]
    fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
    mygex,/path/to/fastqs,any,gex,gene expression,0.5

    [samples]
    sample_id,probe_barcode_ids
    sample1,BC1|BC2
    sample2,BC3
    "#;

        let xtra = XtraData::new("test::");
        let res = MultiConfigCsv::from_reader(csv.as_bytes(), xtra)?;
        let detected_pairing = create_bc_pairing(&[("BC1", "AB1"), ("BC3", "AB3")]);
        assert_json_snapshot!(res.to_multi_graph(
            "test",
            "test description",
            Some(&detected_pairing)
        )?);
        Ok(())
    }

    #[test]
    fn test_cannot_mix_gex_and_lib_chemistry() {
        assert_parse_fails(
            r#"
            [gene-expression]
            ref,mm10-2020-A-chr19
            chemistry,mfrp
            create-bam,true

            [libraries]
            fastq_id,fastqs,lanes,physical_library_id,feature_types,chemistry
            mygex,/path/to/fastqs,any,gex,gene expression,mfrp-rna
            "#,
            "failed to parse CSV: chemistry specified in both the [gene-expression] and [libraries] sections",
        );
    }

    #[test]
    fn test_cannot_use_auto_lib_chem() {
        assert_parse_fails(
            r#"
            [gene-expression]
            ref,mm10-2020-A-chr19
            create-bam,true

            [libraries]
            fastq_id,fastqs,lanes,physical_library_id,feature_types,chemistry
            mygex,/path/to/fastqs,any,gex,gene expression,auto
            "#,
            "[libraries] Specifying auto chemistries at the library level is not supported: (auto).",
        );
    }

    #[test]
    fn test_cannot_provide_partial_lib_chems() {
        assert_parse_fails(
            r#"
            [gene-expression]
            ref,mm10-2020-A-chr19
            create-bam,true

            [libraries]
            fastq_id,fastqs,lanes,physical_library_id,feature_types,chemistry
            mygex,/path/to/fastqs1,any,gex,gene expression,mfrp-rna
            myab,/path/to/fastqs2,any,ab,antibody capture,
            "#,
            "[libraries] A chemistry name must be provided for all libraries or none.",
        );
    }

    #[test]
    fn test_only_flex_chems_in_libs() {
        assert_parse_fails(
            r#"
            [gene-expression]
            ref,mm10-2020-A-chr19
            create-bam,true

            [libraries]
            fastq_id,fastqs,lanes,physical_library_id,feature_types,chemistry
            mygex,/path/to/fastqs1,any,gex,gene expression,SC3Pv1
            myab,/path/to/fastqs2,any,ab,antibody capture,SC3Pv1
            "#,
            "[libraries] Only Flex assays may specify chemistry at the per-library level; invalid chemistries: SC3Pv1",
        );
    }

    #[test]
    fn test_no_conflicting_chems_for_lib_type() {
        assert_parse_fails(
            r#"
            [gene-expression]
            ref,mm10-2020-A-chr19
            create-bam,true

            [libraries]
            fastq_id,fastqs,lanes,physical_library_id,feature_types,chemistry
            mygex,/path/to/fastqs1,any,gex,gene expression,MFRP-RNA
            mygex1,/path/to/fastqs2,any,gex,gene expression,SFRP
            "#,
            "[libraries] Conflicting chemistry for Gene Expression libraries (MFRP-RNA, SFRP); manual chemistry must be the same for all libraries of the same type.",
        );
    }

    #[test]
    fn test_duplicate_sample_ids() {
        let res = MultiConfigCsv::from_csv("test/invalid_csvs/duplicate_sample_ids.csv");
        assert_snapshot!(res.unwrap_err());
    }

    #[test]
    fn test_mismatched_lib_feature_types() {
        let res = MultiConfigCsv::from_csv("test/invalid_csvs/mismatched_lib_feature_types.csv");
        assert_snapshot!(res.unwrap_err());
    }

    #[test]
    fn test_redefined_physical_library_id() {
        let res = MultiConfigCsv::from_csv("test/invalid_csvs/redefined_physical_library_id.csv");
        assert_snapshot!(res.unwrap_err());
    }

    #[test]
    fn test_force_cells_too_low() {
        assert_snapshot!(
            MultiConfigCsv::from_csv("test/invalid_csvs/force_cells_too_low.csv").unwrap_err()
        );
    }

    #[test]
    fn test_duplicate_libraries() {
        assert_snapshot!(
            MultiConfigCsv::from_csv("test/invalid_csvs/duplicate_libraries.csv").unwrap_err()
        );
    }

    #[test]
    fn test_multiplexing_no_samples() {
        assert_snapshot!(
            MultiConfigCsv::from_csv("test/invalid_csvs/multiplexing_no_samples.csv").unwrap_err()
        );
    }

    #[test]
    fn test_multiple_vdjt_libraries() {
        assert_debug_snapshot!(MultiConfigCsv::from_csv("test/multiple_vdj_t.csv"));
    }

    #[test]
    fn test_crispr_with_umi_thresh() {
        insta::assert_debug_snapshot!(MultiConfigCsv::from_csv("test/crispr_with_umi_thresh.csv"));
    }

    #[test]
    fn test_expect_cells_and_force_cells_samples() {
        assert_snapshot!(MultiConfigCsv::from_csv(
            "test/invalid_csvs/expect_cells_and_force_cells_samples.csv"
        )
        .unwrap_err());
    }

    #[test]
    fn test_frp_chemistry_with_samples() {
        assert_snapshot!(
            MultiConfigCsv::from_csv("test/invalid_csvs/frp_chem_with_samples.csv").unwrap_err()
        );
    }

    #[test]
    fn test_mfrp_chemistry_without_samples() {
        assert_snapshot!(
            MultiConfigCsv::from_csv("test/invalid_csvs/mfrp_chem_no_samples.csv").unwrap_err()
        );
    }

    #[test]
    fn test_mfrp_chemistry_with_missing_probe_barcode_entry() {
        assert_snapshot!(MultiConfigCsv::from_csv(
            "test/invalid_csvs/mfrp_chemistry_with_missing_probe_barcode_entry.csv"
        )
        .unwrap_err());
    }

    #[test]
    fn test_use_mfrp_chemistry_set() -> Result<()> {
        let csv = r#"
    [gene-expression]
    ref,mm10-2020-A-chr19
    probe-set,/path/to/probe_set
    chemistry,mfrp
    create-bam,false

    [feature]
    ref,cellranger/multi/feature_refs/20211122_v1.1.csv

    [libraries]
    fastq_id,fastqs,lanes,physical_library_id,feature_types
    mygex,/path/to/fastqs,any,gex,gene expression
    myab,/path/to/fastqs,any,ab,antibody capture

    [samples]
    sample_id,probe_barcode_ids
    sample1,BC1|BC2
    sample2,BC3
    "#;

        let res = MultiConfigCsv::from_reader(csv.as_bytes(), XtraData::new("test::"))?;
        assert_eq!(
            res.gene_expression.as_ref().unwrap().chemistry,
            Some(ChemistryParam::Set(ChemistrySet::Mfrp))
        );
        assert_eq!(
            res.chemistry_specs().unwrap(),
            [
                (
                    LibraryType::Gex,
                    AutoOrRefinedChemistry::Refined(ChemistryName::MFRP_RNA)
                ),
                (
                    LibraryType::Antibody,
                    AutoOrRefinedChemistry::Refined(ChemistryName::MFRP_Ab)
                )
            ]
            .into_iter()
            .collect::<HashMap<_, _>>(),
        );

        Ok(())
    }

    #[test]
    fn test_hashtag_id_missing_gex_lib() {
        let csv = r#"
        [gene-expression]
        ref,GRCh38-2020
        create-bam,true

        [feature]
        ref,path/to/feature_ref.csv

        [libraries]
        fastq_id,fastqs,physical_library_id,feature_types
        myab,/path/to/fastqs,ab,antibody capture

        [samples]
        sample_id, hashtag_ids, description
        1,CD4|CD8a,t_cells
        2,CD19, b_cells
        "#;

        assert_snapshot!(
            MultiConfigCsv::from_reader(csv.as_bytes(), XtraData::new("test::")).unwrap_err()
        );
    }

    #[test]
    fn test_multi_graph_hashtag_fingerprint() -> Result<()> {
        let csv = r#"
    [gene-expression]
    ref,GRCh38-2020
    create-bam,true

    [feature]
    ref,path/to/feature_ref.csv

    [libraries]
    fastq_id,fastqs,physical_library_id,feature_types
    myab,/path/to/fastqs,ab,antibody capture
    myab,/path/to/other/fastqs,gex,gene expression

    [samples]
    sample_id, hashtag_ids, description
    1,CD4|CD8a,t_cells
    2,CD19, b_cells
    "#;

        let xtra = XtraData::new("test::");
        let res = MultiConfigCsv::from_reader(csv.as_bytes(), xtra)?;
        let sample1_fprint = vec![
            Fingerprint::Tagged {
                gem_well: GemWell(1),
                tag_name: "CD4".to_owned(),
                barcode_multiplexing_type: BarcodeMultiplexingType::CellLevel(CellLevel::Hashtag),
                translated_tag_names: vec![],
            },
            Fingerprint::Tagged {
                gem_well: GemWell(1),
                tag_name: "CD8a".to_owned(),
                barcode_multiplexing_type: BarcodeMultiplexingType::CellLevel(CellLevel::Hashtag),
                translated_tag_names: vec![],
            },
        ];
        let mgraph_fprint = res
            .to_multi_graph("sample_id", "sample_desc", None)?
            .samples
            .first()
            .unwrap()
            .fingerprints
            .clone();
        assert_eq!(sample1_fprint, mgraph_fprint);
        Ok(())
    }
}
