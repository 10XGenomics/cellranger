//! Multi config CSV
#![deny(missing_docs)]

pub(super) mod csv;
pub(super) mod parse;
pub mod preflight;
pub(super) mod scsv;

use self::csv::CsvParser;
use self::multiconst::FUNCTIONAL_MAP;
use self::parse::{
    Parse, ParseCtx, parse_prefixed_plate_range, parse_prefixed_range, parse_range, parse_vec,
};
use self::preflight::{
    check_antigen_specificity, check_duplicate_libraries, check_duplicate_sample_barcode_ids,
    check_duplicate_samples, check_feature_functional_map, check_library_chemistries,
    check_library_combinations, check_physical_library_ids, check_platform_restrictions,
};
use self::scsv::{Section, SectionHdr, Span, XtraData, section_csv};
use anyhow::{Context, Result, anyhow, bail, ensure};
use barcode::whitelist::{
    BarcodeId, RTLMultiplexingBarcodeType, ReqStrand, categorize_rtl_multiplexing_barcode_id,
};
use barcode::{BarcodeConstruct, WhitelistSource};
use cloud_utils;
use cr_types::cell_annotation::CellAnnotationModel;
use cr_types::chemistry::{
    AutoChemistryName, AutoOrRefinedChemistry, ChemistryDef, ChemistryName, ChemistrySpecs,
};
use cr_types::constants::DEFAULT_MIN_CRISPR_UMI_THRESHOLD;
use cr_types::reference::feature_reference::{
    BeamMode, FeatureConfig, MHC_ALLELE, NO_ALLELE, SpecificityControls,
};
use cr_types::reference::probe_set_reference::TargetSetFile;
use cr_types::sample_def::{FastqMode, SampleDef};
use cr_types::types::{
    BarcodeMultiplexingType, CellLevel, CrMultiGraph, CrMultiGraphBuilder, Fingerprint, ReadLevel,
};
use cr_types::{
    AlignerParam, FeatureBarcodeType, LibraryType, MinimapPreset, SequencingPlatform,
    TargetingMethod, VdjChainType, enum_maker,
};
use itertools::Itertools;
use martian::{AsMartianPrimaryType, MartianPrimaryType};
use metric::{TxHashMap, TxHashSet};
use regex::Regex;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::convert::{AsRef, TryFrom};
use std::fmt::{Display, Formatter};
use std::fs::File;
use std::hash::Hash;
use std::io::{BufReader, Read};
use std::path::{Path, PathBuf};
use std::str::FromStr;
use strum::{Display as EnumDisplay, EnumString, IntoEnumIterator, IntoStaticStr};

const MIN_FORCE_CELLS: usize = 10;

const DEFAULT_OVERHANG_WL: &str = "overhang";

type NomErr<'a> = nom::Err<nom::error::Error<Span<'a>>>;

const CONTROL_ID: &str = "control_id";
const AG_SPEC_REQ_HDRS: &[&str] = &[CONTROL_ID];
const AG_SPEC_OPT_HDRS: &[&str] = &[MHC_ALLELE];
const FUNCTIONAL_NAME: &str = "functional_name";
const FEATURE_IDS: &str = "feature_ids";
const FUNC_MAP_REQ_HDRS: &[&str] = &[FUNCTIONAL_NAME, FEATURE_IDS];
const FUNC_MAP_OPT_HDRS: &[&str] = &[];
const SEPARATOR: &str = "|";

// Allowed whitelists for multiplexing types.
const OCM_WHITELISTS: &[&str] = &[DEFAULT_OVERHANG_WL];
const CMX_WHITELISTS: &[&str] = &["cmx_rt"];

/// Barcode IDs for multiplexing type detection
pub struct MultiplexingBarcodeIds {
    rtl_barcode_ids: TxHashSet<BarcodeId>,
    ocm_barcode_ids: TxHashSet<BarcodeId>,
    cmx_barcode_ids: TxHashSet<BarcodeId>,
}

impl MultiplexingBarcodeIds {
    /// Load all barcode IDs from all chemistries' whitelists.
    pub fn load() -> Result<Self> {
        let mut rtl_whitelist_specs = TxHashSet::default();

        let ocm_barcode_ids = Self::load_whitelist(OCM_WHITELISTS)?;
        let cmx_barcode_ids = Self::load_whitelist(CMX_WHITELISTS)?;

        // Iterate over all chemistries to collect RTL whitelist specs
        for chem_name in ChemistryName::iter() {
            if chem_name.is_rtl() == Some(true) {
                let chem_def = ChemistryDef::named(chem_name);
                if let BarcodeConstruct::GelBeadAndProbe(construct) = chem_def.barcode_construct() {
                    rtl_whitelist_specs.insert(construct.probe.whitelist().clone());
                }
            }
        }

        let rtl_barcode_ids = Self::load_whitelist_specs(&rtl_whitelist_specs)?;
        Ok(Self {
            rtl_barcode_ids,
            ocm_barcode_ids,
            cmx_barcode_ids,
        })
    }

    /// Load barcode IDs from a list of named whitelists into a single set.
    fn load_whitelist(whitelists: &[&str]) -> Result<TxHashSet<BarcodeId>> {
        let mut ids = TxHashSet::default();
        for &wl in whitelists {
            ids.extend(WhitelistSource::named(wl, true, ReqStrand::Forward)?.get_ids()?);
        }
        Ok(ids)
    }

    /// Load barcode IDs from a set of whitelist specs.
    fn load_whitelist_specs(
        specs: &TxHashSet<barcode::WhitelistSpec>,
    ) -> Result<TxHashSet<BarcodeId>> {
        let mut ids = TxHashSet::default();
        for spec in specs {
            ids.extend(spec.as_source()?.get_ids()?);
        }
        Ok(ids)
    }

    /// Infer multiplexing type from barcode ID
    pub fn infer_type(&self, barcode_id: &BarcodeId) -> Result<ReadLevel> {
        if self.ocm_barcode_ids.contains(barcode_id) {
            Ok(ReadLevel::OH)
        } else if self.cmx_barcode_ids.contains(barcode_id) {
            Ok(ReadLevel::CMX)
        } else if self.rtl_barcode_ids.contains(barcode_id) {
            Ok(ReadLevel::RTL)
        } else {
            bail!("Unknown barcode ID '{barcode_id}'. Must be a valid OCM or RTL barcode ID.")
        }
    }

    /// Validate all IDs in sample are same type
    pub fn validate_consistent_type(&self, barcode_ids: &[BarcodeId]) -> Result<ReadLevel> {
        ensure!(!barcode_ids.is_empty(), "No barcode IDs provided");

        let first_type = self.infer_type(&barcode_ids[0])?;

        for id in &barcode_ids[1..] {
            let id_type = self.infer_type(id)?;
            ensure!(
                id_type == first_type,
                "Mixed types: '{}' is {:?} but '{}' is {:?}",
                barcode_ids[0],
                first_type,
                id,
                id_type
            );
        }

        Ok(first_type)
    }
}
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

enum_maker! {
    /// Sample barcode column names used in the [samples] CSV section.
    #[derive(strum::EnumIter, EnumDisplay)]
    #[expect(missing_docs)]
    SampleBarcodeColumn,
    (CmoIds, "cmo_ids"),
    (HashtagIds, "hashtag_ids"),
    (SampleBarcodeIds, "sample_barcode_ids"),
    (ProbeBarcodeIds, "probe_barcode_ids"),
    (OcmBarcodeIds, "ocm_barcode_ids"),
    (CmxBarcodeIds, "cmx_barcode_ids")
}

struct Bool(bool);

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

pub(crate) struct Ident(String);

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

struct SampleId(String);

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

struct MhcAllele(String);

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
struct AtLeastOne(usize);

impl FromStr for AtLeastOne {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self> {
        let v: usize = s.parse()?;
        ensure!(v >= 1, "must be >= 1");
        Ok(AtLeastOne(v))
    }
}

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
struct Probability(f64);

impl FromStr for Probability {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self> {
        let p: f64 = s.parse()?;
        ensure!((0.0..=1.0).contains(&p), "must be in [0, 1]");
        Ok(Probability(p))
    }
}

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
struct Percent(f64);

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
#[expect(missing_docs, non_camel_case_types)]
pub enum ChemistrySet {
    #[serde(rename = "MFRP")]
    #[strum(to_string = "MFRP")]
    Mfrp,
    #[serde(rename = "MFRP-R1")]
    #[strum(to_string = "MFRP-R1")]
    MfrpR1,
    #[serde(rename = "Flex-v2-R2")]
    #[strum(to_string = "Flex-v2-R2")]
    Flex_v2_R2,
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
    fn is_rtl(&self) -> bool {
        match self {
            Self::Mfrp | Self::MfrpR1 | Self::Flex_v2_R2 => true,
            Self::ThreePrimeV3
            | Self::ThreePrimeV3HT
            | Self::ThreePrimeV4
            | Self::ThreePrimeV4OCM => false,
        }
    }

    fn is_mfrp(&self) -> bool {
        match self {
            Self::Mfrp | Self::MfrpR1 | Self::Flex_v2_R2 => true,
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
            Self::Flex_v2_R2 => match library_type {
                LibraryType::Gex | LibraryType::Crispr | LibraryType::Custom => {
                    Some(ChemistryName::Flex_v2_RNA_R2)
                }
                LibraryType::Antibody => Some(ChemistryName::Flex_v2_Ab_R2_64),
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
                LibraryType::Antibody | LibraryType::Crispr | LibraryType::Custom => {
                    Some(ChemistryName::ThreePrimeV4CS1)
                }
                _ => None,
            },
            Self::ThreePrimeV4OCM => match library_type {
                LibraryType::Gex => Some(ChemistryName::ThreePrimeV4PolyAOCM),
                LibraryType::Antibody | LibraryType::Crispr | LibraryType::Custom => {
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
#[expect(missing_docs)]
pub enum ChemistryParam {
    AutoOrRefined(AutoOrRefinedChemistry),
    Set(ChemistrySet),
}

/// The cell-calling mode (filter method) selected via the `cell-calling-mode` parameter
/// in the `[gene-expression]` section of the multi config CSV.
///
/// Values mirror `FilterMethod` in `lib/python/cellranger/cell_calling_helpers.py`.
/// These values are intentionally excluded:
/// - `manual` is selected by configuring `barcode-sample-assignment`
/// - `topn` is selected by configuring `force-cells`
/// - `targeted` is dead code and equivalent to `gradient`
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize, EnumDisplay, EnumString,
)]
#[serde(rename_all = "snake_case")]
#[strum(serialize_all = "snake_case")]
pub enum CellCallingMode {
    /// Take barcodes within an order of magnitude of the max by count.
    Ordmag,
    /// `Ordmag`, then find barcodes that differ from the ambient profile.
    OrdmagNonambient,
    /// `Ordmag`, then find barcodes above the steepest gradient in the log-log rank plot.
    Gradient,
}

impl AsMartianPrimaryType for CellCallingMode {
    fn as_martian_primary_type() -> MartianPrimaryType {
        MartianPrimaryType::Str
    }
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
            Self::Set(set) => Some(set.is_rtl()),
        }
    }

    /// Return true if this chemistry is in the MFRP family.
    fn is_mfrp(&self) -> bool {
        match self {
            Self::AutoOrRefined(spec) => spec.is_mfrp(),
            Self::Set(set) => set.is_mfrp(),
        }
    }

    /// Return true if this chemistry is in the SFRP family.
    fn is_sfrp(&self) -> bool {
        match self {
            Self::AutoOrRefined(AutoOrRefinedChemistry::Refined(chem)) => chem.is_sfrp(),
            Self::AutoOrRefined(AutoOrRefinedChemistry::Auto(_)) | Self::Set(_) => false,
        }
    }
}

impl Display for ChemistryParam {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
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
            "flex-v2-r1" => Refined(Flex_v2_R1),
            "flex-v2-rna-r2" => Refined(Flex_v2_RNA_R2),
            "flex-v2-ab-r2:45" => Refined(Flex_v2_Ab_R2_45),
            "flex-v2-ab-r2:64" => Refined(Flex_v2_Ab_R2_64),
            "flex-v2-crispr-r2:1" => Refined(Flex_v2_CRISPR_R2_1),
            "flex-v2-singleplex" => Refined(Flex_v2_singleplex),
            "arc-v1" => Refined(ArcV1),
            "combinatorialmultiplexingv1" => Refined(CombinatorialMultiplexingV1),
            _ => bail!("unknown chemistry: {chemistry}"),
        }))
    }
}

/// The gene-expression parameters in the experiment CSV
#[derive(Debug, Default, Serialize)]
#[expect(missing_docs)]
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
    pub cell_calling_mode: Option<CellCallingMode>,
    pub enable_tsne: bool,
    pub no_secondary_analysis: bool,
    pub include_introns: bool,
    pub check_library_compatibility: bool,
    pub aligner: Option<AlignerParam>,
    pub minimap_preset: Option<MinimapPreset>,
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
    /// Hidden parameter: Number of reads to sample for chemistry detection
    pub detect_chemistry_sample_reads: Option<usize>,
    /// Hidden parameter: Total reads to consider for chemistry detection
    pub detect_chemistry_total_reads: Option<usize>,
    /// Hidden parameter: Command line parameters of STAR
    pub star_parameters: Option<String>,
    /// Hidden parameter: Minimum fraction of the single major probe barcode for singleplex FRP libraries
    pub min_major_probe_bc_frac: Option<f64>,
    /// Hidden parameter: minimum UMI length for customers with too-short seq configs
    pub umi_min_read_length: Option<usize>,
}

impl GeneExpressionParams {
    /// Return the probe set.
    pub fn probe_set(&self) -> &[TargetSetFile] {
        &self.probe_set
    }

    /// Return whether a probe set is specified.
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

    /// Return whether expect-cells or force-cells is specified.
    pub fn has_expect_or_force_cells(&self) -> bool {
        self.has_expect_cells() || self.has_force_cells()
    }

    /// Return whether global-minimum-umis is specified.
    pub fn has_global_minimum_umis(&self) -> bool {
        self.global_minimum_umis.is_some()
    }

    /// Return the name of an invalid parameter when antigen capture is used, if any.
    pub fn invalid_parameter_with_antigen_capture(&self) -> Option<&'static str> {
        if self.has_probe_set() {
            Some("probe-set")
        } else if self.cmo_set.is_some() {
            Some("cmo-set")
        } else if self.min_assignment_confidence.is_some() {
            Some("min-assignment-confidence")
        } else if self.barcode_sample_assignment.is_some() {
            Some("barcode-sample-assignment")
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
            eprintln!("{}", cloud_utils::cell_annotation_homedir_msg());
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
        let mut cell_calling_mode: Option<CellCallingMode> = None;
        let mut enable_tsne = false;
        let mut no_secondary_analysis = false;
        let mut include_introns: Option<bool> = None;
        let mut check_library_compatibility = true;
        let mut aligner: Option<AlignerParam> = None;
        let mut minimap_preset: Option<MinimapPreset> = None;
        let mut create_bam = None;
        let mut cmo_set: Option<PathBuf> = None;
        let mut min_assignment_confidence: Option<f64> = None;
        let mut barcode_sample_assignment: Option<PathBuf> = None;
        let mut cell_annotation_model: Option<CellAnnotationModel> = None;
        let mut skip_cell_annotation: bool = false;
        let mut tenx_cloud_token_path: Option<String> = None;
        let mut filter_high_occupancy_gems = true;
        let mut detect_chemistry_sample_reads: Option<usize> = None;
        let mut detect_chemistry_total_reads: Option<usize> = None;
        let mut star_parameters: Option<String> = None;
        let mut min_major_probe_bc_frac: Option<f64> = None;
        let mut umi_min_read_length: Option<usize> = None;
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
                        let probe_set_path = path.parse::<PathBuf>(ctx)?;
                        ensure!(
                            probe_set_path.extension().is_some_and(|ext| ext == "csv"),
                            MultiConfigCsvError::InvalidProbeSetPathExtension { probe_set_path }
                        );
                        probe_set.push(TargetSetFile::from(probe_set_path));
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
                             section needs to be at least {MIN_FORCE_CELLS}. The value you have specified is {parsed_val} \
                             which is too low."
                        );
                        force_cells = Some(parsed_val);
                    }
                }
                "cell-calling-mode" => {
                    if let Some(val) = row.get(1).and_then(empty_is_none) {
                        cell_calling_mode = Some(val.parse::<CellCallingMode>(ctx)?);
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
                "minimap-preset" => {
                    if let Some(val) = row.get(1).and_then(empty_is_none) {
                        minimap_preset = Some(val.parse::<MinimapPreset>(ctx)?);
                    }
                }
                "enable-tsne" => {
                    if let Some(val) = row.get(1).and_then(empty_is_none) {
                        enable_tsne = val.parse::<Bool>(ctx)?.into();
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
                "detect-chemistry-sample-reads" => {
                    if let Some(val) = row.get(1).and_then(empty_is_none) {
                        detect_chemistry_sample_reads = Some(val.parse::<usize>(ctx)?);
                    }
                }
                "detect-chemistry-total-reads" => {
                    if let Some(val) = row.get(1).and_then(empty_is_none) {
                        detect_chemistry_total_reads = Some(val.parse::<usize>(ctx)?);
                    }
                }
                "star-parameters" => {
                    if let Some(val) = row.get(1).and_then(empty_is_none) {
                        star_parameters = Some((*val.fragment()).to_string());
                    }
                }
                "min-major-probe-bc-frac" => {
                    if let Some(val) = row.get(1).and_then(empty_is_none) {
                        min_major_probe_bc_frac = Some(val.parse::<ParseFraction>(ctx)?.0);
                    }
                }
                "umi-min-read-length" => {
                    if let Some(val) = row.get(1).and_then(empty_is_none) {
                        umi_min_read_length = Some(val.parse::<usize>(ctx)?);
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
            ensure!(
                include_introns.is_none(),
                MultiConfigCsvError::FlexIncludeIntrons
            );
        }

        // Filter probes by default.
        if !probe_set.is_empty() && filter_probes.is_none() {
            filter_probes = Some(true);
        }

        // Require minimap-preset if aligner is set to minimap2
        if aligner == Some(AlignerParam::Minimap2) {
            ensure!(
                minimap_preset.is_some(),
                "{ctx} with minimap2 aligner minimap-preset is a required parameter"
            );
        }

        // Validate that minimap-preset requires aligner to be set to minimap2
        if minimap_preset.is_some() {
            ensure!(
                aligner == Some(AlignerParam::Minimap2),
                "{ctx} minimap-preset requires aligner to be set to minimap2"
            );
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
            cell_calling_mode,
            enable_tsne,
            no_secondary_analysis,
            include_introns: include_introns.unwrap_or(true),
            check_library_compatibility,
            aligner,
            minimap_preset,
            create_bam,
            cmo_set,
            min_assignment_confidence,
            barcode_sample_assignment,
            cell_annotation_model,
            skip_cell_annotation,
            tenx_cloud_token_path,
            filter_high_occupancy_gems,
            detect_chemistry_sample_reads,
            detect_chemistry_total_reads,
            star_parameters,
            min_major_probe_bc_frac,
            umi_min_read_length,
        })
    }
}

struct ParseFraction(f64);

impl FromStr for ParseFraction {
    type Err = anyhow::Error;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let val: f64 = s.parse()?;
        ensure!((0.0..=1.0).contains(&val), "must be in the range [0, 1]");
        Ok(Self(val))
    }
}

/// Section [feature] parameters.
#[derive(Debug, Serialize, Deserialize)]
pub struct FeatureParams {
    /// ref or reference or reference-path
    pub reference_path: Option<PathBuf>,
    /// r1-length
    pub r1_length: Option<usize>,
    /// r2-length
    pub r2_length: Option<usize>,
    /// filter-aggregates
    pub filter_aggregates: bool,
    /// min-crispr-umi
    pub min_crispr_umi: usize,
    /// crispr-probe-barcode-location
    /// Optional location for CRISPR probe barcode on R2.
    /// See: CELLRANGER-9564
    pub crispr_probe_barcode_location: Option<usize>,
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
        let mut crispr_probe_barcode_location: Option<usize> = None;

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
                "crispr-probe-barcode-location" => {
                    if let Some(val) = row.get(1).and_then(empty_is_none) {
                        crispr_probe_barcode_location = Some(val.parse::<usize>(ctx)?);
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
            crispr_probe_barcode_location,
        })
    }
}

/// Section [vdj] parameters.
#[derive(Debug, Serialize, Deserialize)]
pub struct VdjParams {
    /// ref or reference or reference-path
    pub reference_path: Option<PathBuf>,
    /// denovo
    pub denovo: Option<bool>,
    /// inner-enrichment-primers
    pub inner_enrichment_primers: Option<PathBuf>,
    /// r1-length
    pub r1_length: Option<usize>,
    /// r2-length
    pub r2_length: Option<usize>,
    /// multiplet-filter
    pub multiplet_filter: Option<bool>,
    /// shared-contig-filter
    pub shared_contig_filter: Option<bool>,
    /// umi-baseline-filter
    pub umi_baseline_filter: Option<bool>,
    /// min-contig-length
    pub min_contig_length: Option<usize>,
    /// skip-clonotyping
    pub skip_clonotyping: Option<bool>,
    /// max-reads-per-barcode
    pub max_reads_per_barcode: Option<usize>,
}

impl<'a> TryFrom<&Section<'a>> for VdjParams {
    type Error = anyhow::Error;

    fn try_from(sec: &Section<'a>) -> Result<Self> {
        let ctx = ParseCtx::Hdr(sec.name);
        let mut reference_path: Option<PathBuf> = None;
        let mut denovo = None;
        let mut inner_enrichment_primers: Option<PathBuf> = None;
        let mut r1_length: Option<usize> = None;
        let mut r2_length: Option<usize> = None;
        let mut multiplet_filter = None;
        let mut shared_contig_filter = None;
        let mut umi_baseline_filter = None;
        let mut min_contig_length = None;
        let mut skip_clonotyping = None;
        let mut max_reads_per_barcode = None;
        for row in &sec.rows {
            if row.is_empty() {
                continue;
            }
            let param = row[0].fragment().to_ascii_lowercase();
            let ctx = ctx.with_col(row[0].fragment());
            match param.as_str() {
                "ref" | "reference" | "reference-path" => {
                    if let Some(path) = row.get(1).and_then(empty_is_none) {
                        reference_path = Some(path.parse::<PathBuf>(ctx)?);
                    }
                }
                "denovo" => {
                    if let Some(val) = row.get(1).and_then(empty_is_none) {
                        denovo = Some(val.parse::<Bool>(ctx)?.into());
                    }
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
                "max-reads-per-barcode" => {
                    if let Some(val) = row.get(1).and_then(empty_is_none) {
                        max_reads_per_barcode = Some(val.parse::<usize>(ctx)?);
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
        if denovo != Some(true) {
            ensure!(reference_path.is_some(), "{ctx} reference is missing");
        }
        if reference_path.is_none() {
            ensure!(
                inner_enrichment_primers.is_some(),
                "{ctx} requires inner_enrichment_primers to be specified if reference is missing"
            );
        }

        Ok(VdjParams {
            reference_path,
            denovo,
            inner_enrichment_primers,
            r1_length,
            r2_length,
            multiplet_filter,
            shared_contig_filter,
            umi_baseline_filter,
            min_contig_length,
            skip_clonotyping,
            max_reads_per_barcode,
        })
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

/// Specifes which lanes to use.
/// TODO: Use `LaneSpec` here.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(from = "Option<Vec<usize>>", into = "Option<Vec<usize>>")]
pub enum Lanes {
    /// Use all lanes.
    Any,
    /// Use the specified lanes.
    Lanes(Vec<usize>),
}

impl Lanes {
    fn overlaps(&self, other: &Lanes) -> bool {
        use self::Lanes::{Any, Lanes};
        match (self, other) {
            (Any, _) | (_, Any) => true,
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

        Ok(Lanes::Lanes(
            lanes
                .into_iter()
                .map(parse_range)
                .flatten_ok()
                .process_results(|iter| iter.sorted())?
                .dedup()
                .collect(),
        ))
    }
}

/// subsample_rate
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
#[serde(transparent)]
pub struct SubsampleRate(f64);

/// How should we interpret the data in the source directory?
#[derive(Debug, Serialize, Deserialize, Clone, PartialEq)]
#[expect(missing_docs)]
pub enum DataType {
    /// FASTQ files from `bcl2fastq`
    Bcl2Fastq {
        fastq_id: String,
        lanes: Lanes,
        /// Fastq ID might be edited to handle case of duplicate fastq IDs referring to different data.
        ilmn_fastq_id: Option<String>,
    },
    /// FASTQ files from `BclProcessor`
    BclProcessor {
        sample_indices: Vec<String>,
        lanes: Lanes,
    },
    /// Long-read files with additional specifications.
    LongRead { format: SequencingDataFormat },
}

/// Known sequencing data formats.
#[derive(Debug, Serialize, Deserialize, Clone, PartialEq)]
pub enum SequencingDataFormat {
    /// Expect to read FASTQ files.
    Fastq,
    /// Expect to read BAM files.
    Bam,
}

/// One row of the section [libraries].
#[derive(Debug, Serialize, Deserialize, Clone, PartialEq)]
#[expect(missing_docs)]
pub struct Library {
    pub physical_library_id: String,
    pub feature_type: LibraryType,
    pub subsample_rate: Option<SubsampleRate>,
    pub chemistry: Option<AutoOrRefinedChemistry>,
    pub read_path: PathBuf,
    pub platform: SequencingPlatform,
    pub data_type: DataType,
}

impl Library {
    /// Return true if the library is CMO.
    pub fn is_multiplexing(&self) -> bool {
        self.feature_type == LibraryType::Cellplex
    }

    /// Return true if the library is VDJ.
    pub fn is_vdj(&self) -> bool {
        self.feature_type.is_vdj()
    }

    fn is_antigen(&self) -> bool {
        self.feature_type.is_fb_type(FeatureBarcodeType::Antigen)
    }

    /// Return true if the library is GEX.
    pub fn is_gex(&self) -> bool {
        self.feature_type.is_gex()
    }

    fn is_crispr(&self) -> bool {
        self.feature_type.is_fb_type(FeatureBarcodeType::Crispr)
    }

    /// Return true if the sequencing data from these two libraries overlap.
    fn overlaps(&self, other: &Library) -> bool {
        if self.read_path != other.read_path {
            return false;
        }
        match &self.data_type {
            DataType::Bcl2Fastq {
                fastq_id,
                lanes,
                ilmn_fastq_id,
            } => {
                let DataType::Bcl2Fastq {
                    fastq_id: fastq_id1,
                    lanes: lanes1,
                    ilmn_fastq_id: ilmn_fastq_id1,
                } = &other.data_type
                else {
                    return false;
                };
                let fastq_id1 = ilmn_fastq_id1.as_ref().unwrap_or(fastq_id1);
                let fastq_id = ilmn_fastq_id.as_ref().unwrap_or(fastq_id);
                if fastq_id1 != fastq_id {
                    return false;
                }
                lanes1.overlaps(lanes)
            }
            DataType::BclProcessor {
                sample_indices,
                lanes,
            } => {
                let DataType::BclProcessor {
                    sample_indices: sample_indices1,
                    lanes: lanes1,
                } = &other.data_type
                else {
                    return false;
                };
                // Check for overlap between sample indices.
                if sample_indices.iter().all(|i| !sample_indices1.contains(i)) {
                    return false;
                }
                // Some index overlap - only a problem if there is lanes overlap.
                lanes1.overlaps(lanes)
            }
            DataType::LongRead { format } => {
                let DataType::LongRead { format: format1 } = &other.data_type else {
                    return false;
                };
                format == format1
            }
        }
    }

    /// Convert to `SampleDef`.
    pub fn to_sample_def(&self) -> SampleDef {
        let library_type = Some({
            // CELLRANGER-7889: this was previously an domain boundary where all
            // VDJ chain types would be collapsed down into a single "VDJ" type.
            // Retain this invariant for the moment by collapsing all VDJ chain types
            // down to Auto.
            let library_type = self.feature_type;
            if library_type.is_vdj() {
                LibraryType::VdjAuto
            } else {
                library_type
            }
        });
        let gem_group = Some(1);
        match &self.data_type {
            DataType::Bcl2Fastq {
                fastq_id,
                lanes,
                ilmn_fastq_id,
            } => {
                let sample_name = match ilmn_fastq_id {
                    Some(orig) => orig.clone(),
                    None => fastq_id.clone(),
                };
                SampleDef {
                    fastq_mode: FastqMode::ILMN_BCL2FASTQ,
                    gem_group,
                    lanes: lanes.clone().into(),
                    library_type,
                    r1_length: None,
                    r2_length: None,
                    read_path: self.read_path.clone(),
                    sample_indices: None,
                    sample_names: Some(vec![sample_name]),
                    subsample_rate: self.subsample_rate.map(|sr| sr.0),
                    fastq_id: Some(fastq_id.clone()),
                }
            }
            DataType::BclProcessor {
                sample_indices,
                lanes,
            } => {
                // Heuristic to use the flowcell id as the fastq_id based on how internal fastqs
                // are along paths like:
                // /somepath/pipestances/HKWWVDSXC/BCL_PROCESSOR_PD/...
                let pattern = r"([^/]+)/BCL_PROCESSOR_PD";
                let re = Regex::new(pattern).unwrap();

                let flow_cell_id = re
                    .captures(self.read_path.to_str().unwrap())
                    .and_then(|caps| caps.get(1).map(|m| String::from(m.as_str())));
                SampleDef {
                    fastq_mode: FastqMode::BCL_PROCESSOR,
                    gem_group,
                    lanes: lanes.clone().into(),
                    library_type,
                    r1_length: None,
                    r2_length: None,
                    read_path: self.read_path.clone(),
                    sample_indices: Some(sample_indices.clone()),
                    sample_names: None,
                    subsample_rate: self.subsample_rate.map(|sr| sr.0),
                    fastq_id: flow_cell_id,
                }
            }
            DataType::LongRead { format: _format } => SampleDef {
                // TODO: use the format spec here instead of platform, SampleDef
                // shouldn't know about platforms
                fastq_mode: match self.platform {
                    SequencingPlatform::PacBio => FastqMode::PACBIO,
                    SequencingPlatform::OxfordNanopore => FastqMode::ONT,
                    SequencingPlatform::Sbx => FastqMode::SBX,
                    SequencingPlatform::Illumina | SequencingPlatform::Ultima => unreachable!(),
                },
                gem_group,
                lanes: None,
                library_type,
                r1_length: None,
                r2_length: None,
                read_path: self.read_path.clone(),
                sample_indices: None,
                sample_names: None,
                subsample_rate: self.subsample_rate.map(|sr| sr.0),
                fastq_id: None,
            },
        }
    }
}

mod libsconst {
    pub(super) const FASTQ_ID: &str = "fastq_id";
    pub(super) const FASTQS: &str = "fastqs";
    pub(super) const PHYSICAL_LIBRARY_ID: &str = "physical_library_id";
    pub(super) const FEATURE_TYPES: &str = "feature_types";
    pub(super) const LIBS_REQ_HDRS: &[&str] = &[FASTQ_ID, FASTQS, FEATURE_TYPES];

    pub(super) const LANES: &str = "lanes";
    pub(super) const SUBSAMPLE_RATE: &str = "subsample_rate";
    pub(super) const CHEMISTRY: &str = "chemistry";
    pub(super) const LIBS_OPT_HDRS: &[&str] =
        &[PHYSICAL_LIBRARY_ID, LANES, SUBSAMPLE_RATE, CHEMISTRY];

    pub(super) const FASTQ_PATH: &str = "fastq_path";
    pub(super) const SAMPLE_INDICES: &str = "sample_indices";
    pub(super) const LIBS_INT_REQ_HDRS: &[&str] = &[FASTQ_PATH, SAMPLE_INDICES, FEATURE_TYPES];
    pub(super) const LIBS_INT_OPT_HDRS: &[&str] =
        &[PHYSICAL_LIBRARY_ID, LANES, SUBSAMPLE_RATE, CHEMISTRY];

    pub(super) const PLATFORM: &str = "platform";
    pub(super) const READ_PATH: &str = "read_path";
    pub(super) const LIBS_ANY_REQ_HDRS: &[&str] = &[PLATFORM, READ_PATH, FEATURE_TYPES];
    pub(super) const LIBS_ANY_OPT_HDRS: &[&str] = &[
        FASTQ_ID,
        PHYSICAL_LIBRARY_ID,
        LANES,
        SUBSAMPLE_RATE,
        CHEMISTRY,
    ];
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq, EnumDisplay, EnumString)]
enum LibrariesCsvFormat {
    /// format: fastq_path,sample_indices,lanes,physical_library_id,feature_types,subsample_rate,chemistry
    Internal,
    /// format: fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate,chemistry
    DefaultCs,
    /// format: platform,read_path,feature_types,fastq_id,physical_library_id,lanes,subsample_rate,chemistry
    Universal,
}

/// Section [libraries].
/// TODO: should these just be SampleDef?
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct LibrariesCsv {
    /// The rows of the section [libraries].
    pub libs: Vec<Library>,
}

impl LibrariesCsv {
    /// Create a new LibrariesCsv with the specified libraries and format
    pub fn new(libraries: Vec<Library>) -> Self {
        Self { libs: libraries }
    }

    /// Return an iterator over the library types.
    fn library_types(&self) -> impl Iterator<Item = LibraryType> + '_ {
        self.libs.iter().map(|l| l.feature_type)
    }

    /// Return true if any libraries are gene expression.
    pub fn has_gene_expression(&self) -> bool {
        self.has_library_type(LibraryType::Gex)
    }

    /// Return true if any libraries are Antibody/Crispr/Antigen/Custom.
    /// Return false for CMO multiplexing, which comes with a built-in list.
    fn has_feature_barcode(&self) -> bool {
        self.library_types().any(|ft| {
            let Some(fb_type) = ft.feature_barcode_type() else {
                return false;
            };
            fb_type != FeatureBarcodeType::Multiplexing
        })
    }

    /// Return true if any libraries are of the provided type.
    fn has_library_type(&self, library_type: LibraryType) -> bool {
        self.library_types().any(|ft| ft == library_type)
    }

    fn has_feature_barcode_type(&self, feature_barcode_type: FeatureBarcodeType) -> bool {
        self.has_library_type(LibraryType::FeatureBarcodes(feature_barcode_type))
    }

    /// Return true if any libraries are antibody capture.
    pub fn has_antibody_capture(&self) -> bool {
        self.has_feature_barcode_type(FeatureBarcodeType::Antibody)
    }

    /// Return true if any libraries are CRISPR guide capture.
    fn has_crispr_guide_capture(&self) -> bool {
        self.has_feature_barcode_type(FeatureBarcodeType::Crispr)
    }

    /// Return true if any libraries are custom capture.
    fn has_custom_capture(&self) -> bool {
        self.has_feature_barcode_type(FeatureBarcodeType::Custom)
    }

    /// Return true if any libraries are antigen capture.
    fn has_antigen_capture(&self) -> bool {
        self.has_feature_barcode_type(FeatureBarcodeType::Antigen)
    }

    /// Return true if any libraries are CMO multiplexing.
    fn has_multiplexing(&self) -> bool {
        self.has_feature_barcode_type(FeatureBarcodeType::Multiplexing)
    }

    fn has_vdj(&self) -> bool {
        self.library_types().any(|x| x.is_vdj())
    }

    fn has_vdj_t_or_gd(&self) -> bool {
        self.has_library_type(LibraryType::Vdj(VdjChainType::VdjT))
            || self.has_library_type(LibraryType::Vdj(VdjChainType::VdjTGD))
    }

    fn has_vdj_b(&self) -> bool {
        self.has_library_type(LibraryType::Vdj(VdjChainType::VdjB))
    }

    /// Return the BEAM mode, if any VDJ libraries are present.
    pub fn beam_mode(&self) -> Option<BeamMode> {
        if self.has_vdj_t_or_gd() {
            Some(BeamMode::BeamT)
        } else if self.has_vdj_b() {
            Some(BeamMode::BeamAB)
        } else {
            None
        }
    }

    /// Return the map of chemistry specs, for any libraries that provided one.
    /// Return None if chemistry was not provided at the library level.
    ///
    /// Validation must ensure that this is either None, or a map containing
    /// an entry for every library. Validation also must ensure that different
    /// libraries of the same feature type have matching chemistries.
    fn chemistry_specs(&self) -> Option<ChemistrySpecs> {
        let specs: ChemistrySpecs = self
            .libs
            .iter()
            .filter_map(|lib| lib.chemistry.map(|chem| (lib.feature_type, chem)))
            .collect();
        if specs.is_empty() {
            return None;
        }
        Some(specs)
    }

    /// Return the sequencing platform for GEX libraries, if any are present.
    pub fn gex_platform(&self) -> Option<SequencingPlatform> {
        self.libs
            .iter()
            .find(|lib| lib.is_gex())
            .map(|l| l.platform)
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
            CHEMISTRY, FASTQ_ID, FASTQ_PATH, FASTQS, FEATURE_TYPES, LANES, LIBS_ANY_OPT_HDRS,
            LIBS_ANY_REQ_HDRS, LIBS_INT_OPT_HDRS, LIBS_INT_REQ_HDRS, LIBS_OPT_HDRS, LIBS_REQ_HDRS,
            PHYSICAL_LIBRARY_ID, SAMPLE_INDICES, SUBSAMPLE_RATE,
        };
        let hdr = sec.name;
        let (parser, format) = {
            if let Ok(parser) = CsvParser::new(sec.clone(), LIBS_ANY_REQ_HDRS, LIBS_ANY_OPT_HDRS) {
                (parser, LibrariesCsvFormat::Universal)
            } else if let Ok(parser) =
                CsvParser::new(sec.clone(), LIBS_INT_REQ_HDRS, LIBS_INT_OPT_HDRS)
            {
                (parser, LibrariesCsvFormat::Internal)
            } else {
                (
                    CsvParser::new(sec.clone(), LIBS_REQ_HDRS, LIBS_OPT_HDRS)?,
                    LibrariesCsvFormat::DefaultCs,
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
                format!(
                    "{}_{gem_group}",
                    library_type_short(feature_type),
                    gem_group = 1,
                )
            };
            let subsample_rate = parser
                .find_opt(row, SUBSAMPLE_RATE)?
                .and_then(empty_is_none)
                .map(|sr| sr.parse::<Probability>(ctx.with_col(SUBSAMPLE_RATE)))
                .transpose()?
                .map(|sr| sr.0)
                .map(SubsampleRate);
            let chemistry = parser
                .find_opt(row, CHEMISTRY)?
                .and_then(empty_is_none)
                .map(|sr| sr.parse::<ParseAutoOrRefinedChemistry>(ctx.with_col(CHEMISTRY)))
                .transpose()?
                .map(|param| param.0);

            if format == LibrariesCsvFormat::Internal {
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

                data.push(Library {
                    physical_library_id,
                    feature_type,
                    subsample_rate,
                    chemistry,
                    read_path: PathBuf::from(parser.find_req(row, FASTQ_PATH)?.fragment()),
                    platform: SequencingPlatform::Illumina,
                    data_type: DataType::BclProcessor {
                        sample_indices,
                        lanes,
                    },
                });
            } else if format == LibrariesCsvFormat::DefaultCs {
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

                data.push(Library {
                    physical_library_id,
                    feature_type,
                    subsample_rate,
                    chemistry,
                    read_path: PathBuf::from(parser.find_req(row, FASTQS)?.fragment()),
                    platform: SequencingPlatform::Illumina,
                    data_type: DataType::Bcl2Fastq {
                        fastq_id,
                        ilmn_fastq_id,
                        lanes,
                    },
                });
            } else {
                let platform = parser
                    .find_req(row, libsconst::PLATFORM)?
                    .parse::<SequencingPlatform>(ctx.with_col(libsconst::PLATFORM))?;
                let read_path =
                    PathBuf::from(parser.find_req(row, libsconst::READ_PATH)?.fragment());
                let data_type = match platform {
                    SequencingPlatform::Illumina => {
                        let mut fastq_id = parser
                            .find_opt(row, FASTQ_ID)?
                            .ok_or_else(|| anyhow!("{ctx} is missing required field '{FASTQ_ID}'"))?
                            .parse::<Ident>(ctx.with_col(FASTQ_ID))?
                            .to_string();
                        let mut ilmn_fastq_id = None;
                        let fastq_count = fastq_id_counts.entry(fastq_id.clone()).or_insert(0);
                        *fastq_count += 1;

                        if *fastq_count > 1 {
                            ilmn_fastq_id = Some(fastq_id.clone());
                            fastq_id = format!("{fastq_id} ({fastq_count})");
                        }
                        DataType::Bcl2Fastq {
                            fastq_id,
                            ilmn_fastq_id,
                            lanes,
                        }
                    }
                    SequencingPlatform::OxfordNanopore | SequencingPlatform::PacBio => {
                        DataType::LongRead {
                            format: SequencingDataFormat::Bam,
                        }
                    }
                    SequencingPlatform::Sbx => DataType::LongRead {
                        format: SequencingDataFormat::Fastq,
                    },
                    SequencingPlatform::Ultima => {
                        todo!("library specifications for additional sequencers");
                    }
                };
                data.push(Library {
                    physical_library_id,
                    feature_type,
                    subsample_rate,
                    chemistry,
                    read_path,
                    platform,
                    data_type,
                });
            }
        }
        check_duplicate_libraries(&data)?;
        check_platform_restrictions(&data)?;
        check_physical_library_ids(&data)?;
        check_library_combinations(&data)?;
        check_library_chemistries(&data)?;
        Ok(LibrariesCsv::new(data))
    }
}

/// Multiple probe barcodes may be grouped explicitly using this character.
/// This enables cross-library declarations of different probe barcodes for
/// different library types.
pub const PROBE_BARCODE_ID_GROUPING: &str = "+";

/// A row of the section [samples].
#[derive(Debug, Default, Serialize, Deserialize)]
pub struct SampleRow {
    sample_id: String,
    cmo_ids: Option<Vec<String>>,
    hashtag_ids: Option<Vec<String>>,
    /// NOTE: these may be +-concatenated groupings. Access via the
    /// sample_barcode_ids method to control how they are unpacked.
    probe_barcode_ids: Option<Vec<String>>,
    overhang_ids: Option<Vec<String>>,
    cmx_barcode_ids: Option<Vec<String>>,
    sample_barcode_ids: Option<Vec<String>>, // Unified column for sample multiplexing barcode IDs (RTL/OH/CMX).
    description: String,
    expect_cells: Option<usize>,
    force_cells: Option<usize>,
    emptydrops_minimum_umis: Option<usize>,
    global_minimum_umis: Option<usize>,
    max_mito_percent: Option<f64>,
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
    pub(super) fn sample_barcode_ids(
        &self,
        probe_barcode_iteration_mode: ProbeBarcodeIterationMode,
    ) -> Option<impl Iterator<Item = &str> + '_> {
        self.sample_barcode_id_groupings(probe_barcode_iteration_mode)
            .map(|bcs| bcs.into_iter().flatten())
    }

    /// Return all of the sample barcode IDs, either CMO or probe barcode IDs.
    /// Retain any existing |-separated groupings in the result.
    /// Control probe barcode iteration using the provided option.
    pub(super) fn sample_barcode_id_groupings(
        &self,
        probe_barcode_iteration_mode: ProbeBarcodeIterationMode,
    ) -> Option<Vec<Vec<&str>>> {
        match (
            self.cmo_ids.as_ref(),
            self.hashtag_ids.as_ref(),
            self.probe_barcode_ids.as_ref(),
            self.overhang_ids.as_ref(),
            self.cmx_barcode_ids.as_ref(),
            self.sample_barcode_ids.as_ref(),
        ) {
            (Some(x), None, None, None, None, None)
            | (None, Some(x), None, None, None, None)
            | (None, None, None, Some(x), None, None)
            | (None, None, None, None, Some(x), None) => {
                // cmo_ids or hashtag_ids or overhang_ids or cmx_barcode_ids
                Some(x.iter().map(|v| vec![v.as_str()]).collect())
            }
            (None, None, Some(x), None, None, None) | (None, None, None, None, None, Some(x)) => {
                // probe_barcode_ids or sample_barcode_ids.
                // We already validate that overhang_ids and cmx_barcode_ids won't contain the PROBE_BARCODE_ID_GROUPING character.
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

            (None, None, None, None, None, None) => None,
            _ => panic!("found more than one source of sample barcode IDs"),
        }
    }
}

/// Section [samples].
#[derive(Debug, Serialize)]
pub struct SamplesCsv {
    /// The rows of the [samples] section.
    pub rows: Vec<SampleRow>,
    /// The multiplexing type for this section.
    pub multiplexing_type: BarcodeMultiplexingType,
    /// The actual column name used for barcode IDs in the CSV.
    pub sample_barcode_column_name: SampleBarcodeColumn,
}

impl SamplesCsv {
    /// Infer the ReadLevel subtype (RTL/OH/CMX) from barcode IDs in a cell value.
    /// Uses only the first barcode ID for inference; full validation happens in row parsing.
    fn infer_read_level_multiplexing_type<'a>(
        barcode_ids_span: Span<'a>,
        multiplexing_barcode_ids: &MultiplexingBarcodeIds,
    ) -> Result<ReadLevel> {
        let ids_vec = parse_vec(barcode_ids_span.clone())
            .map_err(|e| anyhow!("invalid barcode IDs format: {e}"))?;
        let first_entry = ids_vec
            .first()
            .ok_or_else(|| anyhow!("barcode IDs list is empty"))?
            .clone();
        let ids = parse_prefixed_range(first_entry)?; // Get first entry from a range '-'
        // Split on '+' to strip any cross-library groupings (e.g. "BC001+AB001" → "BC001").
        let first_bc_id = ids[0]
            .split(PROBE_BARCODE_ID_GROUPING)
            .next()
            .unwrap_or(&ids[0]);

        let barcode_id = BarcodeId::pack(first_bc_id);
        multiplexing_barcode_ids.infer_type(&barcode_id)
    }

    /// Determine the barcode multiplexing type and column name from a specific row of a samples section.
    fn infer_multiplexing_type_from_row<'a>(
        parser: &CsvParser<'a, SectionHdr<'_>>,
        row: usize,
        hdr: SectionHdr<'_>,
        multiplexing_barcode_ids: Option<&MultiplexingBarcodeIds>,
    ) -> Result<(BarcodeMultiplexingType, SampleBarcodeColumn)> {
        let ctx = ParseCtx::HdrRow(hdr, row + 1);

        // whitelist based inference is needed based on the column contents.
        let sample_barcode =
            parser.find_opt(row, SampleBarcodeColumn::SampleBarcodeIds.as_str())?;
        // Columns for which type is explicit from the column name, no whitelist inference needed.
        let cmo = parser
            .find_opt(row, SampleBarcodeColumn::CmoIds.as_str())?
            .is_some();
        let hashtag = parser
            .find_opt(row, SampleBarcodeColumn::HashtagIds.as_str())?
            .is_some();
        let probe = parser
            .find_opt(row, SampleBarcodeColumn::ProbeBarcodeIds.as_str())?
            .is_some();
        let oh = parser
            .find_opt(row, SampleBarcodeColumn::OcmBarcodeIds.as_str())?
            .is_some();
        let cmx = parser
            .find_opt(row, SampleBarcodeColumn::CmxBarcodeIds.as_str())?
            .is_some();

        match (cmo, hashtag, sample_barcode.is_some(), probe, oh, cmx) {
            (true, false, false, false, false, false) => Ok((
                BarcodeMultiplexingType::CellLevel(CellLevel::CMO),
                SampleBarcodeColumn::CmoIds,
            )),
            (false, true, false, false, false, false) => Ok((
                BarcodeMultiplexingType::CellLevel(CellLevel::Hashtag),
                SampleBarcodeColumn::HashtagIds,
            )),
            (false, false, false, true, false, false) => Ok((
                BarcodeMultiplexingType::ReadLevel(ReadLevel::RTL),
                SampleBarcodeColumn::ProbeBarcodeIds,
            )),
            (false, false, false, false, true, false) => Ok((
                BarcodeMultiplexingType::ReadLevel(ReadLevel::OH),
                SampleBarcodeColumn::OcmBarcodeIds,
            )),
            (false, false, false, false, false, true) => Ok((
                BarcodeMultiplexingType::ReadLevel(ReadLevel::CMX),
                SampleBarcodeColumn::CmxBarcodeIds,
            )),
            // New unified column: infer type from barcode whitelist content.
            (false, false, true, false, false, false) => {
                let ids = multiplexing_barcode_ids.ok_or_else(|| {
                    anyhow!("multiplexing_barcode_ids required for read-level multiplexing")
                })?;
                let read_level_type =
                    Self::infer_read_level_multiplexing_type(sample_barcode.unwrap().clone(), ids)
                        .with_context(|| {
                            format!(
                                "{ctx} invalid {}",
                                SampleBarcodeColumn::SampleBarcodeIds.as_str()
                            )
                        })?;
                Ok((
                    BarcodeMultiplexingType::ReadLevel(read_level_type),
                    SampleBarcodeColumn::SampleBarcodeIds,
                ))
            }
            (false, false, false, false, false, false) => {
                // TODO(CR11): add sample_barcode_ids back to error message
                bail!(
                    "{ctx} requires one of {}, {}, {} or {} columns to be specified",
                    SampleBarcodeColumn::CmoIds.as_str(),
                    SampleBarcodeColumn::HashtagIds.as_str(),
                    SampleBarcodeColumn::ProbeBarcodeIds.as_str(),
                    SampleBarcodeColumn::OcmBarcodeIds.as_str()
                )
            }
            _ => {
                // TODO(CR11): add sample_barcode_ids back to error message
                bail!(
                    "{ctx} has more than one of the following mutually exclusive columns: \
                     {}, {}, {}, and/or {}",
                    SampleBarcodeColumn::CmoIds.as_str(),
                    SampleBarcodeColumn::HashtagIds.as_str(),
                    SampleBarcodeColumn::ProbeBarcodeIds.as_str(),
                    SampleBarcodeColumn::OcmBarcodeIds.as_str()
                )
            }
        }
    }

    fn has_force_cells(&self) -> bool {
        self.rows.iter().any(|sample| sample.force_cells.is_some())
    }

    fn has_expect_cells(&self) -> bool {
        self.rows.iter().any(|sample| sample.expect_cells.is_some())
    }

    fn has_expect_or_force_cells(&self) -> bool {
        self.has_expect_cells() || self.has_force_cells()
    }

    fn has_global_minimum_umis(&self) -> bool {
        self.rows
            .iter()
            .any(|sample| sample.global_minimum_umis.is_some())
    }

    /// Return the column `expect_cells`.
    pub fn get_expect_cells(&self) -> TxHashMap<String, Option<f64>> {
        self.rows
            .iter()
            .map(|sample| {
                (
                    sample.sample_id.clone(),
                    sample.expect_cells.map(|x| x as f64),
                )
            })
            .collect()
    }

    /// Return the column `force_cells`.
    pub fn get_force_cells(&self) -> TxHashMap<String, Option<f64>> {
        self.rows
            .iter()
            .map(|sample| {
                (
                    sample.sample_id.clone(),
                    sample.force_cells.map(|x| x as f64),
                )
            })
            .collect()
    }

    /// Return the column `emptydrops_minimum_umis`.
    pub fn get_emptydrops_minimum_umis(&self) -> TxHashMap<String, Option<f64>> {
        self.rows
            .iter()
            .map(|sample| {
                (
                    sample.sample_id.clone(),
                    sample.emptydrops_minimum_umis.map(|x| x as f64),
                )
            })
            .collect()
    }

    /// Return the column `global_minimum_umis`.
    pub fn get_global_minimum_umis(&self) -> TxHashMap<String, Option<f64>> {
        self.rows
            .iter()
            .map(|sample| {
                (
                    sample.sample_id.clone(),
                    (sample.global_minimum_umis.map(|x| x as f64)),
                )
            })
            .collect()
    }

    /// Return the column `max_mito_percent`.
    pub fn get_max_mito_percent(&self) -> TxHashMap<String, Option<f64>> {
        self.rows
            .iter()
            .map(|sample| (sample.sample_id.clone(), sample.max_mito_percent))
            .collect()
    }

    fn is_rtl_multiplexed(&self) -> bool {
        self.multiplexing_type == BarcodeMultiplexingType::ReadLevel(ReadLevel::RTL)
    }

    pub(super) fn get_cmo_sample_map(&self) -> TxHashMap<String, String> {
        let mut res = TxHashMap::default();
        for row in &self.rows {
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

    pub(super) fn get_hashtag_sample_map(&self) -> TxHashMap<String, String> {
        let mut res = TxHashMap::default();
        for row in &self.rows {
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

    /// Iterate over tuples of sample ID and multiplexing barcode IDs.
    /// sample_id,BC001+AB001|BC002|AB002
    /// [(sample_id, [[BC001, AB001], [BC002, AB002]])]
    fn iter_sample_and_probe_barcodes(
        &self,
    ) -> impl Iterator<
        Item = (
            SampleId,
            impl Iterator<Item = impl Iterator<Item = BarcodeId> + '_> + '_,
        ),
    > + '_ {
        self.rows.iter().map(|row| {
            (
                SampleId::from_str(&row.sample_id).unwrap(),
                row.probe_barcode_ids
                    .iter()
                    .flatten()
                    .map(|x| x.split(PROBE_BARCODE_ID_GROUPING).map(BarcodeId::pack)),
            )
        })
    }

    /// Return all distinct RTLMultiplexingBarcodeType in the probe_id column.
    fn rtl_multiplexing_barcode_types(
        &self,
    ) -> impl Iterator<Item = Option<RTLMultiplexingBarcodeType>> + '_ {
        self.iter_sample_and_probe_barcodes()
            .flat_map(|(_, ids)| ids)
            .flatten()
            .map(|x| categorize_rtl_multiplexing_barcode_id(x.as_str()))
            .unique()
    }

    /// Return whether these multiplexing barcode types need translation.
    fn has_translated_multiplexing_barcode_types(&self) -> bool {
        self.rtl_multiplexing_barcode_types().any(|x| x.is_some())
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
        self.iter_sample_and_probe_barcodes()
            .flat_map(|(_sample_id, barcode_ids)| {
                barcode_ids.filter_map(|mut pieces| pieces.next().map(|target| (target, pieces)))
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

pub(super) mod samplesconst {
    use super::SampleBarcodeColumn;
    use strum::IntoEnumIterator;

    pub(super) const SAMPLE_ID: &str = "sample_id";
    pub(super) const DESCRIPTION: &str = "description";
    pub(super) const EXPECT_CELLS: &str = "expect_cells";
    pub(super) const FORCE_CELLS: &str = "force_cells";
    pub(super) const EMPTYDROPS_MINIMUM_UMIS: &str = "emptydrops_minimum_umis";
    pub(super) const GLOBAL_MINIMUM_UMIS: &str = "global_minimum_umis";
    pub(super) const MAX_MITO_FRAC: &str = "max_mito_percent";
    pub(super) const SAMP_REQ_HDRS: &[&str] = &[SAMPLE_ID];

    pub(super) fn samp_opt_hdrs() -> Vec<&'static str> {
        let mut headers: Vec<&'static str> = SampleBarcodeColumn::iter()
            .map(|variant| variant.as_str())
            .collect();

        headers.extend_from_slice(&[
            DESCRIPTION,
            EXPECT_CELLS,
            FORCE_CELLS,
            EMPTYDROPS_MINIMUM_UMIS,
            GLOBAL_MINIMUM_UMIS,
            MAX_MITO_FRAC,
        ]);

        headers
    }
}

impl<'a> TryFrom<&Section<'a>> for SamplesCsv {
    type Error = anyhow::Error;

    fn try_from(sec: &Section<'a>) -> Result<Self> {
        use samplesconst::{
            DESCRIPTION, EMPTYDROPS_MINIMUM_UMIS, EXPECT_CELLS, FORCE_CELLS, GLOBAL_MINIMUM_UMIS,
            MAX_MITO_FRAC, SAMP_REQ_HDRS, SAMPLE_ID, samp_opt_hdrs,
        };

        let hdr = sec.name;
        let opt_hdrs = samp_opt_hdrs();
        let parser = CsvParser::new(sec.clone(), SAMP_REQ_HDRS, &opt_hdrs)?;
        let mut data = vec![];

        // Check if we need to load whitelists (only for sample_barcode_ids column)
        let first_row = parser
            .rows()
            .next()
            .ok_or_else(|| anyhow!("no rows in samples section"))?;
        let needs_whitelist_loading = parser
            .find_opt(first_row, SampleBarcodeColumn::SampleBarcodeIds.as_str())?
            .is_some();

        // Load barcode whitelists once if needed
        let multiplexing_barcode_ids = if needs_whitelist_loading {
            Some(MultiplexingBarcodeIds::load()?)
        } else {
            None
        };

        // Determine multiplexing type and column name based on the first row
        let (barcode_multiplexing_type, sample_barcode_column_name) =
            SamplesCsv::infer_multiplexing_type_from_row(
                &parser,
                first_row,
                hdr,
                multiplexing_barcode_ids.as_ref(),
            )?;

        for row in parser.rows() {
            let ctx = ParseCtx::HdrRow(hdr, row + 1);

            // Validate that this row has the same multiplexing type as the first row
            let (row_multiplexing_type, _) = SamplesCsv::infer_multiplexing_type_from_row(
                &parser,
                row,
                hdr,
                multiplexing_barcode_ids.as_ref(),
            )?;
            ensure!(
                row_multiplexing_type == barcode_multiplexing_type,
                "{ctx} has multiplexing type {row_multiplexing_type:?} but expected {barcode_multiplexing_type:?} based on the first row"
            );

            let sample_id: String = parser
                .find_req(row, SAMPLE_ID)?
                .parse::<SampleId>(ctx.with_col(SAMPLE_ID))?
                .into();

            let cmo_ids = parser.find_opt(row, SampleBarcodeColumn::CmoIds.as_str())?;
            let hashtag_ids = parser.find_opt(row, SampleBarcodeColumn::HashtagIds.as_str())?;
            let sample_barcode_ids_col =
                parser.find_opt(row, SampleBarcodeColumn::SampleBarcodeIds.as_str())?;
            let probe_barcode_ids =
                parser.find_opt(row, SampleBarcodeColumn::ProbeBarcodeIds.as_str())?;
            let overhang_ids = parser.find_opt(row, SampleBarcodeColumn::OcmBarcodeIds.as_str())?;
            let cmx_barcode_ids =
                parser.find_opt(row, SampleBarcodeColumn::CmxBarcodeIds.as_str())?;

            // Ensure only one sample barcode column is present
            let sample_barcode_ids = match (
                cmo_ids,
                hashtag_ids,
                sample_barcode_ids_col,
                probe_barcode_ids,
                overhang_ids,
                cmx_barcode_ids,
            ) {
                (Some(ids), None, None, None, None, None)
                | (None, Some(ids), None, None, None, None)
                | (None, None, Some(ids), None, None, None)
                | (None, None, None, Some(ids), None, None)
                | (None, None, None, None, Some(ids), None)
                | (None, None, None, None, None, Some(ids)) => Some(ids),
                (None, None, None, None, None, None) => None,
                // TODO(CR11): add sample_barcode_ids back to error message
                _ => bail!(
                    "{ctx} has more than one of the following mutually exclusive columns: \
                     {}, {}, {}, and/or {}",
                    SampleBarcodeColumn::CmoIds.as_str(),
                    SampleBarcodeColumn::HashtagIds.as_str(),
                    SampleBarcodeColumn::ProbeBarcodeIds.as_str(),
                    SampleBarcodeColumn::OcmBarcodeIds.as_str()
                ),
            };
            let expect_cells = parser.find_opt(row, EXPECT_CELLS)?;
            let force_cells = parser.find_opt(row, FORCE_CELLS)?;
            let emptydrops_minimum_umis = parser.find_opt(row, EMPTYDROPS_MINIMUM_UMIS)?;
            let global_minimum_umis = parser.find_opt(row, GLOBAL_MINIMUM_UMIS)?;
            let max_mito_percent = parser.find_opt(row, MAX_MITO_FRAC)?;
            if probe_barcode_ids.is_none()
                && overhang_ids.is_none()
                && cmx_barcode_ids.is_none()
                && sample_barcode_ids_col.is_none()
            {
                // TODO(CR11): add sample_barcode_ids back to error message
                ensure!(
                    expect_cells.is_none(),
                    "{ctx} has {EXPECT_CELLS} column specified \
                     without the required {} or {} column",
                    SampleBarcodeColumn::ProbeBarcodeIds.as_str(),
                    SampleBarcodeColumn::OcmBarcodeIds.as_str()
                );
                // TODO(CR11): add sample_barcode_ids back to error message
                ensure!(
                    force_cells.is_none(),
                    "{ctx} has {FORCE_CELLS} column specified \
                     without the required {} or {} column",
                    SampleBarcodeColumn::ProbeBarcodeIds.as_str(),
                    SampleBarcodeColumn::OcmBarcodeIds.as_str()
                );
            };
            let sample_barcode_column_name = match barcode_multiplexing_type {
                BarcodeMultiplexingType::CellLevel(CellLevel::CMO) => SampleBarcodeColumn::CmoIds,
                BarcodeMultiplexingType::CellLevel(CellLevel::Hashtag) => {
                    SampleBarcodeColumn::HashtagIds
                }
                BarcodeMultiplexingType::ReadLevel(_) => {
                    if sample_barcode_ids_col.is_some() {
                        SampleBarcodeColumn::SampleBarcodeIds
                    } else {
                        match barcode_multiplexing_type {
                            BarcodeMultiplexingType::ReadLevel(ReadLevel::RTL) => {
                                SampleBarcodeColumn::ProbeBarcodeIds
                            }
                            BarcodeMultiplexingType::ReadLevel(ReadLevel::OH) => {
                                SampleBarcodeColumn::OcmBarcodeIds
                            }
                            BarcodeMultiplexingType::ReadLevel(ReadLevel::CMX) => {
                                SampleBarcodeColumn::CmxBarcodeIds
                            }
                            _ => unreachable!(),
                        }
                    }
                }
            };

            let sample_barcode_ids: Option<Vec<_>> = if let Some(sample_barcode_ids) =
                sample_barcode_ids.and_then(empty_is_none)
            {
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
                        .map(|s| match barcode_multiplexing_type {
                            BarcodeMultiplexingType::ReadLevel(ReadLevel::CMX) => {
                                parse_prefixed_plate_range(s)
                            }
                            _ => parse_prefixed_range(s),
                        })
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
                                    ensure!(
                                        !multiplexing_id.contains(PROBE_BARCODE_ID_GROUPING),
                                        "must not contain '{PROBE_BARCODE_ID_GROUPING}' (barcode grouping is only supported for probe barcodes)"
                                    );
                                    Ok(multiplexing_id)
                                }
                                BarcodeMultiplexingType::ReadLevel(ReadLevel::CMX) => {
                                    ensure!(
                                        !multiplexing_id.contains(PROBE_BARCODE_ID_GROUPING),
                                        "must not contain '{PROBE_BARCODE_ID_GROUPING}' (barcode grouping is only supported for probe barcodes)"
                                    );
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

                Some(
                    sample_barcode_ids_iter
                        .process_results(|iter| iter.sorted())?
                        .dedup()
                        .collect(),
                )
            } else {
                None
            };

            // When using the unified sample_barcode_ids column, validate that all barcodes
            // in this row (including all parts of '+' groupings) are the same ReadLevel type.
            if sample_barcode_ids_col.is_some()
                && let (Some(parsed_ids), Some(mux_ids)) =
                    (&sample_barcode_ids, &multiplexing_barcode_ids)
            {
                let bc_ids: Vec<BarcodeId> = parsed_ids
                    .iter()
                    .flat_map(|id| id.split(PROBE_BARCODE_ID_GROUPING).map(BarcodeId::pack))
                    .collect();
                mux_ids.validate_consistent_type(&bc_ids).with_context(|| {
                    format!(
                        "{ctx} has mixed barcode types in {}",
                        SampleBarcodeColumn::SampleBarcodeIds.as_str()
                    )
                })?;
            }

            let description = parser
                .find_opt(row, DESCRIPTION)?
                .map_or_else(String::default, ToString::to_string);
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

            // Populate fields based on the sample barcode column that we are using.
            let mut final_cmo_ids = None;
            let mut final_hashtag_ids = None;
            let mut final_read_level_ids = None;
            let mut final_probe_ids = None;
            let mut final_oh_ids = None;
            let mut final_cmx_ids = None;

            match sample_barcode_column_name {
                SampleBarcodeColumn::CmoIds => final_cmo_ids = sample_barcode_ids,
                SampleBarcodeColumn::HashtagIds => final_hashtag_ids = sample_barcode_ids,
                SampleBarcodeColumn::SampleBarcodeIds => final_read_level_ids = sample_barcode_ids,
                SampleBarcodeColumn::ProbeBarcodeIds => final_probe_ids = sample_barcode_ids,
                SampleBarcodeColumn::OcmBarcodeIds => final_oh_ids = sample_barcode_ids,
                SampleBarcodeColumn::CmxBarcodeIds => final_cmx_ids = sample_barcode_ids,
            }

            data.push(SampleRow {
                sample_id,
                cmo_ids: final_cmo_ids,
                hashtag_ids: final_hashtag_ids,
                sample_barcode_ids: final_read_level_ids,
                probe_barcode_ids: final_probe_ids,
                overhang_ids: final_oh_ids,
                cmx_barcode_ids: final_cmx_ids,
                description,
                force_cells: if matches!(
                    barcode_multiplexing_type,
                    BarcodeMultiplexingType::ReadLevel(_)
                ) {
                    force_cells
                } else {
                    None
                },
                expect_cells: if matches!(
                    barcode_multiplexing_type,
                    BarcodeMultiplexingType::ReadLevel(_)
                ) {
                    expect_cells
                } else {
                    None
                },
                emptydrops_minimum_umis,
                global_minimum_umis,
                max_mito_percent,
            });
        }
        check_duplicate_sample_barcode_ids(&data, sample_barcode_column_name)?;
        check_duplicate_samples(&data)?;
        Ok(SamplesCsv {
            rows: data,
            multiplexing_type: barcode_multiplexing_type,
            sample_barcode_column_name,
        })
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

    let invalid = invalid_chars!(bcid, 'A'..='Z' | '0'..='9' | '-' | '+');
    ensure!(
        invalid.is_empty(),
        "invalid character(s): '{}', must contain only \
         uppercase letters (A-Z), digits (0-9), hyphen (-), and plus (+)",
        invalid.iter().format("', '")
    );
    Ok(bcid)
}

/// Returns the set of default overhang IDs
fn get_default_overhang_set() -> TxHashSet<BarcodeId> {
    WhitelistSource::named(DEFAULT_OVERHANG_WL, true, ReqStrand::Forward)
        .unwrap()
        .get_ids()
        .unwrap()
        .into_iter()
        .collect()
}

/// One row of the section [antigen-specificity].
#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct AntigenSpecificityRow {
    control_id: String,
    mhc_allele: Option<String>,
}

/// The section [antigen-specificity].
#[derive(Debug, Serialize, Deserialize, Clone)]
#[serde(transparent)]
pub struct AntigenSpecificityCsv(Vec<AntigenSpecificityRow>);

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

/// One row of the section [feature-functional-map].
#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct FunctionalMapRow {
    functional_name: String,
    feature_ids: Vec<String>,
}

/// The section [feature-functional-map].
#[derive(Debug, Serialize, Deserialize, Clone)]
#[serde(transparent)]
pub struct FunctionalMapCsv(Vec<FunctionalMapRow>);

impl<'a> TryFrom<&Section<'a>> for FunctionalMapCsv {
    type Error = anyhow::Error;

    fn try_from(sec: &Section<'a>) -> Result<Self> {
        let parser = CsvParser::new(sec.clone(), FUNC_MAP_REQ_HDRS, FUNC_MAP_OPT_HDRS)?;
        let mut data = vec![];
        for row in parser.rows() {
            let functional_name = parser.find_req(row, FUNCTIONAL_NAME)?.to_string();
            let feature_ids = parser
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

/// Create a `FeatureConfig` from the various feature-barcoding related CSV sections.
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
        if samples.multiplexing_type == BarcodeMultiplexingType::CellLevel(CellLevel::Hashtag) {
            Some(
                samples
                    .rows
                    .iter()
                    .filter_map(|s| s.sample_barcode_ids(ProbeBarcodeIterationMode::All))
                    .flatten()
                    .map(ToString::to_string)
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

fn create_specificity_controls(
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

#[expect(missing_docs)]
pub mod multiconst {
    pub const GENE_EXPRESSION: &str = "gene-expression";
    pub(super) const GEX: &str = "gex";
    pub const FEATURE: &str = "feature";
    pub const VDJ: &str = "vdj";
    pub(super) const LIBRARIES: &str = "libraries";
    pub(super) const LIBS: &str = "libs";
    pub(crate) const SAMPLES: &str = "samples";
    pub(super) const ANTIGEN_SPECIFICITY: &str = "antigen-specificity";
    pub(super) const FUNCTIONAL_MAP: &str = "feature-functional-map";

    /// Default value of include-introns for non-hybcap GEX analyses.
    pub const DEFAULT_INCLUDE_INTRONS: bool = true;
}

mod types {
    #![expect(missing_docs)]
    use martian_derive::martian_filetype;
    use serde::{Deserialize, Serialize};
    martian_filetype!(MultiConfigCsvFile, "csv");
}

pub use types::MultiConfigCsvFile;

impl MultiConfigCsvFile {
    /// Read this multi config CSV.
    pub fn read(&self) -> Result<MultiConfigCsv> {
        MultiConfigCsv::from_csv(self)
    }
}

/// A multi config CSV file.
#[derive(Debug, Serialize)]
pub struct MultiConfigCsv {
    /// Section [gene-expression]
    pub gene_expression: Option<GeneExpressionParams>,
    /// Section [feature]
    pub feature: Option<FeatureParams>,
    /// Section [vdj]
    pub vdj: Option<VdjParams>,
    /// Section [libraries]
    pub libraries: LibrariesCsv,
    /// Section [samples]
    pub samples: Option<SamplesCsv>,
    /// Section [antigen-specificity]
    pub antigen_specificity: Option<AntigenSpecificityCsv>,
    /// Section [feature-functional-map]
    pub functional_map: Option<FunctionalMapCsv>,
}

/// A representation of a `cellranger multi` configuration
impl MultiConfigCsv {
    /// Load a MultiConfigCsv from a path
    pub fn from_csv(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        MultiConfigCsv::from_reader(
            BufReader::new(File::open(path).with_context(|| path.display().to_string())?),
            XtraData::new(path),
        )
    }

    /// Load an MultiConfigCsv from any `impl Read` along with parsing `XtraData`
    pub fn from_reader<R: Read, X: Into<XtraData>>(mut reader: R, xtra: X) -> Result<Self> {
        use multiconst::{
            ANTIGEN_SPECIFICITY, FEATURE, GENE_EXPRESSION, GEX, LIBRARIES, LIBS, SAMPLES, VDJ,
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
                SAMPLES => 2,
                _ => 3,
            }
        });
        // I tried to implement TryFrom, but the compiler complained of overlapping impls
        //   and I gave up before figuring out why
        let mut builder = MultiConfigCsvBuilder::default();
        for section in &sections {
            builder.push(section)?;
        }
        let cfg = builder.build()?;

        // TODO: Define multi_graph_representation for Ont/Pacbio, curently fails to find FastqDef
        if !cfg.libraries.libs.iter().any(|lib| {
            matches!(
                lib.platform,
                SequencingPlatform::OxfordNanopore
                    | SequencingPlatform::PacBio
                    | SequencingPlatform::Sbx
            )
        }) {
            // this forces additional consistency checking of [libraries] section
            let _ = cfg.to_multi_graph("dummy_sample", "dummy_desc", None)?;
        }
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
        let mut builder = CrMultiGraphBuilder::new();
        for lib in &self.libraries.libs {
            builder.push_library(lib.physical_library_id.clone(), lib.feature_type)?;
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
            for sample in &samples.rows {
                if let Some(sample_barcode_ids) =
                    sample.sample_barcode_ids(ProbeBarcodeIterationMode::Mapped)
                {
                    for sample_barcode_id in sample_barcode_ids {
                        builder.push_fingerprint(
                            sample.sample_id.clone(),
                            sample.description.clone(),
                            Fingerprint::tagged(
                                sample_barcode_id.to_string(),
                                translated_source_map
                                    .get(sample_barcode_id)
                                    .map(|bcs| bcs.iter().map(ToString::to_string).collect())
                                    .unwrap_or_default(),
                                samples.multiplexing_type,
                            ),
                        )?;
                    }
                } else {
                    builder.push_fingerprint(
                        sample.sample_id.clone(),
                        sample.description.clone(),
                        Fingerprint::Untagged,
                    )?;
                };
            }
        } else {
            builder.push_fingerprint(
                (*sample_id).to_string(),
                (*sample_desc).to_string(),
                Fingerprint::Untagged,
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
    ) -> TxHashSet<&str> {
        if let Some(samples) = &self.samples {
            samples
                .rows
                .iter()
                .filter_map(|x| x.sample_barcode_ids(probe_barcode_iteration_mode))
                .flatten()
                .collect()
        } else {
            TxHashSet::default()
        }
    }

    /// Return true if multiplexed using overhangs.
    pub fn is_overhang_multiplexed(&self) -> bool {
        self.samples.as_ref().is_some_and(|s| {
            s.multiplexing_type == BarcodeMultiplexingType::ReadLevel(ReadLevel::OH)
        })
    }

    /// Return true if multiplexed using CMX barcodes.
    pub fn is_cmx_multiplexed(&self) -> bool {
        self.samples.as_ref().is_some_and(|s| {
            s.multiplexing_type == BarcodeMultiplexingType::ReadLevel(ReadLevel::CMX)
        })
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

    /// Return whether these multiplexing barcode types need translation.
    pub fn has_translated_multiplexing_barcode_types(&self) -> bool {
        self.samples
            .as_ref()
            .is_some_and(SamplesCsv::has_translated_multiplexing_barcode_types)
    }

    /// Return the [samples] multiplexing type, if any.
    pub fn barcode_multiplexing_type(&self) -> Option<BarcodeMultiplexingType> {
        self.samples.as_ref().map(|s| s.multiplexing_type)
    }

    /// Return true if any barcode multiplexing is configured.
    pub fn is_multiplexed(&self) -> bool {
        self.barcode_multiplexing_type().is_some()
    }

    /// Return true if cell-level barcode multiplexing is configured.
    pub fn is_cell_multiplexed(&self) -> bool {
        self.barcode_multiplexing_type()
            .is_some_and(|t| matches!(t, BarcodeMultiplexingType::CellLevel(_)))
    }

    /// Return true if read-level barcode multiplexing is configured.
    pub fn is_read_multiplexed(&self) -> bool {
        self.barcode_multiplexing_type()
            .is_some_and(|t| matches!(t, BarcodeMultiplexingType::ReadLevel(_)))
    }

    /// Return an iterator over the library types.
    pub fn library_types(&self) -> impl Iterator<Item = LibraryType> + '_ {
        self.libraries.library_types()
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
                .libs
                .iter()
                .map(|lib| match spec {
                    ChemistryParam::AutoOrRefined(chem) => Ok((lib.feature_type, chem)),
                    ChemistryParam::Set(chem_set) => Ok((
                        lib.feature_type,
                        AutoOrRefinedChemistry::Refined(
                            chem_set.chemistry_for_library_type(lib.feature_type)?,
                        ),
                    )),
                })
                .collect();
        };
        assert!(gex_chemistry_spec.is_none());
        Ok(specs)
    }
}

#[derive(Debug, Default)]
struct MultiConfigCsvBuilder {
    gene_expression: Option<GeneExpressionParams>,
    feature: Option<FeatureParams>,
    vdj: Option<VdjParams>,
    libraries: Option<LibrariesCsv>,
    samples: Option<SamplesCsv>,
    antigen_specificity: Option<AntigenSpecificityCsv>,
    functional_map: Option<FunctionalMapCsv>,
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

#[derive(thiserror::Error, Debug, PartialEq)]
enum MultiConfigCsvError {
    #[error(
        "Conflicting chemistry specifications: presence of Flex libraries is implied by the \
         presence of probe-barcode-ids but all libraries have specified non-Flex chemistries."
    )]
    ConflictingFlexChemistry,
    #[error("A reference transcriptome is required to analyze a non-Flex Gene Expression library.")]
    MissingRefPathForGex,
    #[error(
        "[{}] section is missing `probe-set` a required parameter for Flex",
        multiconst::GENE_EXPRESSION
    )]
    FlexMissingProbeSet,
    #[error(
        "The [{}] section is missing both `probe-set` and `ref` parameters; at least one is required \
         to analyze a `Gene Expression` library.",
        multiconst::GENE_EXPRESSION
    )]
    GexMissingProbeSetAndRefPath,
    #[error(
        "The [{}] section specifies the parameter include-introns, which is not valid for Flex \
         chemistries.",
        multiconst::GENE_EXPRESSION
    )]
    FlexIncludeIntrons,
    #[error(
        "The 'probe-set' parameter specified under the '[gene-expression]' section must have a 'csv' filename extension: {probe_set_path:?}"
    )]
    InvalidProbeSetPathExtension { probe_set_path: PathBuf },
}

impl MultiConfigCsvBuilder {
    fn push(&mut self, section: &Section<'_>) -> Result<()> {
        use multiconst::{
            ANTIGEN_SPECIFICITY, FEATURE, GENE_EXPRESSION, GEX, LIBRARIES, LIBS, SAMPLES, VDJ,
        };
        let name = section.name.fragment().to_ascii_lowercase();
        match name.as_str() {
            GENE_EXPRESSION | GEX => self.gene_expression(section),
            FEATURE => self.feature(section),
            VDJ => self.vdj(section),
            LIBRARIES | LIBS => self.libraries(section),
            SAMPLES => self.samples(section),
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

    fn samples(&mut self, section: &Section<'_>) -> Result<()> {
        ensure!(
            self.samples.is_none(),
            "failed to parse CSV, duplicate [{}] at line: {}, col: {}",
            section.name.fragment(),
            section.name.location_line(),
            section.name.get_utf8_column(),
        );
        self.samples = Some(SamplesCsv::try_from(section)?);
        Ok(())
    }

    fn build(self) -> Result<MultiConfigCsv> {
        let MultiConfigCsvBuilder {
            gene_expression: gene_expression_owned,
            feature,
            vdj: vdj_owned,
            libraries,
            samples: samples_owned,
            antigen_specificity,
            functional_map,
            ..
        } = self;

        let Some(libraries) = libraries else {
            bail!(
                "failed to parse CSV, [{}] section not provided",
                multiconst::LIBRARIES
            );
        };

        let gene_expression = gene_expression_owned.as_ref();
        let vdj = vdj_owned.as_ref();
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

        let has_probe_barcode_ids = samples.is_some_and(|s| {
            s.multiplexing_type == BarcodeMultiplexingType::ReadLevel(ReadLevel::RTL)
        });
        let chem_rtls: TxHashSet<_> = any_source_chemistries
            .iter()
            .map(ChemistryParam::is_rtl)
            .collect();
        // If any chemistry is RTL or if probe barcode IDs are provided then rtl_is_implied is true.
        let rtl_is_implied = has_probe_barcode_ids || chem_rtls.contains(&Some(true));
        // If all chemistries are specified as not RTL, then no_rtl_is_implied is true.
        let no_rtl_is_implied = chem_rtls
            .iter()
            .exactly_one()
            .is_ok_and(|&x| x == Some(false));
        let is_rtl = match (rtl_is_implied, no_rtl_is_implied) {
            (false, false) => None, // Could be RTL, because at least one chemistry is unspecified
            (false, true) => Some(false), // Not RTL because all chemistries are not RTL
            (true, false) => Some(true), // RTL because has_probe_barcode_ids
            (true, true) => bail!(MultiConfigCsvError::ConflictingFlexChemistry),
        };
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

            if libraries.has_library_type(LibraryType::Vdj(VdjChainType::Auto)) {
                ensure!(
                    vdj.unwrap().reference_path.is_some(),
                    "failed to parse CSV: if [{}] reference is omitted, then the feature-types in the [{}] section are required to be specified as: {}, {} or {}",
                    multiconst::VDJ,
                    multiconst::LIBRARIES,
                    LibraryType::Vdj(VdjChainType::VdjT),
                    LibraryType::Vdj(VdjChainType::VdjB),
                    LibraryType::Vdj(VdjChainType::VdjTGD),
                );
            }
        }

        if libraries.has_multiplexing() {
            ensure!(
                samples.is_some(),
                "failed to parse CSV: [{}] section omitted \
                 but Multiplexing Capture libraries provided",
                multiconst::SAMPLES,
            );
        }

        let has_hashtag_ids = samples.is_some_and(|s| {
            s.multiplexing_type == BarcodeMultiplexingType::CellLevel(CellLevel::Hashtag)
        });
        if has_hashtag_ids {
            ensure!(
                libraries.has_gene_expression() && libraries.has_antibody_capture(),
                "failed to parse CSV: [{}] section specified {} parameter \
                 which require both Gene Expression and Antibody Capture libraries.",
                multiconst::SAMPLES,
                SampleBarcodeColumn::HashtagIds.as_str(),
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
                "failed to parse CSV: [{FUNCTIONAL_MAP}] section is provided \
                 but no feature barcode libraries provided"
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

        if any_source_chemistries.iter().any(ChemistryParam::is_sfrp) {
            ensure!(
                samples.is_none(),
                "failed to parse CSV: specified `{chems}` as the `chemistry` but a \
                 [{samples}] section is also provided",
                chems = any_source_chemistries
                    .iter()
                    .filter(|chem| chem.is_sfrp())
                    .format(", "),
                samples = multiconst::SAMPLES,
            );
        }

        if let Some(mfrp_chem) = any_source_chemistries.iter().find(|spec| spec.is_mfrp()) {
            ensure!(
                has_probe_barcode_ids,
                "failed to parse CSV: specified `{}` as the `chemistry` but a \
                 [{}] section with {} column is not provided",
                mfrp_chem,
                multiconst::SAMPLES,
                SampleBarcodeColumn::ProbeBarcodeIds.as_str(),
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

        if libraries.has_gene_expression() && is_rtl == Some(false) {
            ensure!(
                gene_expression.is_some_and(|x| x.reference_path.is_some()),
                MultiConfigCsvError::MissingRefPathForGex,
            );
        }

        if is_rtl == Some(true) {
            if libraries.has_gene_expression() {
                ensure!(
                    gene_expression.is_some_and(GeneExpressionParams::has_probe_set),
                    MultiConfigCsvError::FlexMissingProbeSet,
                );
            }

            // Disallow setting include-introns for RTL chemistries.
            ensure!(
                gene_expression
                    .is_none_or(|x| x.include_introns == multiconst::DEFAULT_INCLUDE_INTRONS),
                MultiConfigCsvError::FlexIncludeIntrons,
            );
        }

        if libraries.has_gene_expression() {
            ensure!(
                gene_expression.is_some_and(GeneExpressionParams::has_probe_set)
                    || gene_expression.is_some_and(|x| x.reference_path.is_some()),
                MultiConfigCsvError::GexMissingProbeSetAndRefPath,
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
                SampleBarcodeColumn::ProbeBarcodeIds.as_str(),
                SampleBarcodeColumn::ProbeBarcodeIds.as_str(),
            );
        }

        // Ensure we only use the crispr-probe-barcode-location feature with
        // appropriate settings (multiplexed Flex with a manual CRISPR chemistry)
        if feature
            .as_ref()
            .is_some_and(|feat| feat.crispr_probe_barcode_location.is_some())
        {
            ensure!(
                has_probe_barcode_ids,
                "`crispr-probe-barcode-location` may only be used with a \
                 multiplexed Flex analysis but no [{}] section is present",
                SampleBarcodeColumn::ProbeBarcodeIds.as_str()
            );
            ensure!(
                libraries.has_crispr_guide_capture(),
                "`crispr-probe-barcode-location` is set but no CRISPR Guide \
                 Capture library is present."
            );
            ensure!(
                !any_source_chemistries.is_empty()
                    && any_source_chemistries
                        .iter()
                        .all(|chem| chem.is_rtl() == Some(true)),
                "`crispr-probe-barcode-location` may only be used with \
                 manually specified RTL chemistries"
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
            SampleBarcodeColumn::ProbeBarcodeIds.as_str(),
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
            vdj: vdj_owned,
            libraries,
            samples: samples_owned,
            antigen_specificity,
            functional_map,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::scsv::XtraData;
    use super::*;
    use anyhow::Result;
    use barcode::whitelist::BarcodeId;
    use cr_types::chemistry::{AutoOrRefinedChemistry, ChemistryName};
    use cr_types::reference::probe_set_reference::TargetSetFile;
    use cr_types::{BarcodeMultiplexingType, CellLevel, Fingerprint, LibraryType};
    use insta::{assert_debug_snapshot, assert_json_snapshot, assert_snapshot};
    use itertools::Itertools;
    use metric::TxHashMap;
    use std::collections::HashMap;
    use std::path::{Path, PathBuf};

    // initialize insta test harness
    #[ctor::ctor(unsafe)]
    fn init() {
        // this ensures insta knows where to find its snap tests
        let cwd = std::env::current_dir().unwrap();
        let workspace_root = cwd.parent().unwrap();
        unsafe { std::env::set_var("INSTA_WORKSPACE_ROOT", workspace_root) }
    }

    /// Assert that `observed` contains `expected`, and panic with a message showing both strings if not.
    pub(super) fn assert_contains(observed: &str, expected: &str) {
        assert!(
            observed.contains(expected),
            "assertion `observed.contains(expected)` failed\n\
             observed: {observed}\n\
             expected: {expected}",
        );
    }

    /// Assert that parsing the provided CSV contents fails with a matching error.
    fn assert_parse_fails(contents: &str, expected_error: &str) {
        assert_contains(
            &MultiConfigCsv::from_reader(contents.as_bytes(), XtraData::new("tests"))
                .unwrap_err()
                .to_string(),
            expected_error,
        );
    }

    /// Assert that parsing the provided CSV succeeds.
    fn assert_parse_ok(contents: &str) -> MultiConfigCsv {
        MultiConfigCsv::from_reader(contents.as_bytes(), XtraData::new("tests")).unwrap()
    }

    #[test]
    fn test_gex_without_transcriptome_no_chem() {
        let csv = "
            [gene-expression]
            create-bam,false
            [libraries]
            fastq_id,fastqs,feature_types
            gex,/path/to/fastqs,Gene Expression"
            .as_bytes();
        let exp_err = MultiConfigCsvError::GexMissingProbeSetAndRefPath;
        let act_err = MultiConfigCsv::from_reader(csv, XtraData::new(""))
            .unwrap_err()
            .downcast::<MultiConfigCsvError>()
            .unwrap();
        assert_eq!(act_err, exp_err);
    }

    #[test]
    fn test_gex_without_transcriptome_sfrp_chem() {
        let csv = "
            [gene-expression]
            create-bam,false
            [libraries]
            fastq_id,fastqs,feature_types,chemistry
            gex,/path/to/fastqs,Gene Expression,SFRP"
            .as_bytes();
        let exp_err = MultiConfigCsvError::FlexMissingProbeSet;
        let act_err = MultiConfigCsv::from_reader(csv, XtraData::new(""))
            .unwrap_err()
            .downcast::<MultiConfigCsvError>()
            .unwrap();
        assert_eq!(act_err, exp_err);
    }

    #[test]
    fn test_gex_without_transcriptome_3pv3_chem() {
        let csv = "
            [gene-expression]
            create-bam,false
            chemistry,SC3Pv3
            [libraries]
            fastq_id,fastqs,feature_types
            gex,/path/to/fastqs,Gene Expression"
            .as_bytes();
        let exp_err = MultiConfigCsvError::MissingRefPathForGex;
        let act_err = MultiConfigCsv::from_reader(csv, XtraData::new(""))
            .unwrap_err()
            .downcast::<MultiConfigCsvError>()
            .unwrap();
        assert_eq!(act_err, exp_err);
    }

    #[test]
    fn test_conflicting_rtl_chem() {
        let csv = "
            [gene-expression]
            create-bam,false
            chemistry,SC3Pv3

            [libraries]
            fastq_id,fastqs,feature_types
            gex,/path/to/fastqs,Gene Expression

            [samples]
            sample_id,probe_barcode_ids
            whoami,BC1
            "
        .as_bytes();
        let exp_err = MultiConfigCsvError::ConflictingFlexChemistry;
        let act_err = MultiConfigCsv::from_reader(csv, XtraData::new(""))
            .unwrap_err()
            .downcast::<MultiConfigCsvError>()
            .unwrap();
        assert_eq!(act_err, exp_err);
    }

    #[test]
    fn load_simple_external_cmos() -> Result<()> {
        let csv = r"
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
fastq_id,fastqs,physical_library_id,feature_types
mygex,/path/to/fastqs,a,Gene Expression
mycmo,/path/to/fastqs,1,Multiplexing Capture

[samples]
sample_id,cmo_ids,description
whoami,cmo_2,hidad!
";
        let xtra = XtraData::new("tests::load_simple_external_cmos");
        let exp = MultiConfigCsv::from_reader(csv.as_bytes(), xtra)?;
        println!("{exp:?}");
        Ok(())
    }

    #[test]
    fn load_simple_external_probe_bcs() -> Result<()> {
        let csv = r"
[gene-expression]
ref,/path/to/gex/ref
probe-set,/path/to/probe_set.csv
create-bam,true

[libraries]
fastq_id,fastqs,physical_library_id,feature_types
mygex,/path/to/fastqs,a,Gene Expression
mycmo,/path/to/fastqs,1,Multiplexing Capture

[samples]
sample_id,probe_barcode_ids,description,force_cells,expect_cells,
whoami,BC1,hidad!,,100
";
        let xtra = XtraData::new("tests::load_simple_external_probe_bcs");
        let exp = MultiConfigCsv::from_reader(csv.as_bytes(), xtra)?;
        println!("{exp:?}");
        Ok(())
    }

    #[test]
    fn test_probe_set_tsv_ext() {
        let csv: &'static str = r"
[gene-expression]
probe-set,/path/to/probe_set.tsv
create-bam,true

[libraries]
fastq_id,fastqs,physical_library_id,feature_types
mygex,/path/to/fastqs,a,Gene Expression

[samples]
sample_id,probe_barcode_ids,description,force_cells,expect_cells,
whoami,BC1,hidad!,,100
";
        let exp_err = MultiConfigCsvError::InvalidProbeSetPathExtension {
            probe_set_path: PathBuf::from("/path/to/probe_set.tsv"),
        };
        let act_err = MultiConfigCsv::from_reader(csv.as_bytes(), XtraData::new(""))
            .unwrap_err()
            .downcast::<MultiConfigCsvError>()
            .unwrap();
        assert_eq!(act_err, exp_err);
    }

    #[test]
    fn test_probe_set_no_ext() {
        let csv: &'static str = r"
[gene-expression]
probe-set,/path/to/probe_set
create-bam,true

[libraries]
fastq_id,fastqs,physical_library_id,feature_types
mygex,/path/to/fastqs,a,Gene Expression

[samples]
sample_id,probe_barcode_ids,description,force_cells,expect_cells,
whoami,BC1,hidad!,,100
";
        let exp_err = MultiConfigCsvError::InvalidProbeSetPathExtension {
            probe_set_path: PathBuf::from("/path/to/probe_set"),
        };
        let act_err = MultiConfigCsv::from_reader(csv.as_bytes(), XtraData::new(""))
            .unwrap_err()
            .downcast::<MultiConfigCsvError>()
            .unwrap();
        assert_eq!(act_err, exp_err);
    }

    #[test]
    fn load_include_introns_default() -> Result<()> {
        let csv = r"
[gene-expression]
ref,/path/to/gex/ref
probe-set,/path/to/probe_set.csv
create-bam,true

[libraries]
fastq_id,fastqs,physical_library_id,feature_types
mygex,/path/to/fastqs,a,Gene Expression
mycmo,/path/to/fastqs,1,Multiplexing Capture

[samples]
sample_id,probe_barcode_ids,description,force_cells,expect_cells,
whoami,BC1,hidad!,,100
";
        let xtra = XtraData::new("tests::load_introns_for_rtl");
        let exp = MultiConfigCsv::from_reader(csv.as_bytes(), xtra)?;
        assert!(exp.gene_expression.unwrap().include_introns);

        let csv_with_introns = r"
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
fastq_id,fastqs,physical_library_id,feature_types
mygex,/path/to/fastqs,a,Gene Expression
mycmo,/path/to/fastqs,1,Multiplexing Capture

[samples]
sample_id,cmo_ids,description
whoami,cmo_2,hidad!
";
        let xtra2 = XtraData::new("tests::include_introns");
        let exp2 = MultiConfigCsv::from_reader(csv_with_introns.as_bytes(), xtra2)?;
        assert!(exp2.gene_expression.unwrap().include_introns);
        let csv_with_no_introns = r"
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
fastq_id,fastqs,physical_library_id,feature_types
mygex,/path/to/fastqs,a,Gene Expression
mycmo,/path/to/fastqs,1,Multiplexing Capture

[samples]
sample_id,cmo_ids,description
whoami,cmo_2,hidad!
";
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
fastq_path,sample_indices,lanes,physical_library_id,feature_types,subsample_rate,
/path/to/fastqs,SI-001,1-3|5,whatsit_gex,Gene Expression,0.5
/path/to/fastqs,SI-002,4,whatsit_cmo,Multiplexing Capture,1

[samples]
sample_id,cmo_ids,description
abc,cmo_1-3|cmo_5,"hi mom"
def,cmo_4,hi dad!
"#;
        let xtra = XtraData::new("tests::load_simple_internal");
        let exp = MultiConfigCsv::from_reader(csv.as_bytes(), xtra)?;
        println!("{exp:?}");
        Ok(())
    }

    #[test]
    fn invalid_sample_id_chars() {
        let csv = r"
[gene-expression]
ref,GRCh38-2020-A-chr21
create-bam,true

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types
bamtofastq,fastqs/cellranger/multi/VDJ_GEX_small_multi/small_gex_fastqs_chr21_new,any,gex_1,gene expression
pbmc_1k_protein_v3_antibody,../dui_tests/test_resources/cellranger-count/pbmc_1k_protein_v3_antibody,any,cmo1,Multiplexing Capture

[samples]
sample_id,cmo_ids,description
IAmABad*SampleId,CMO1,some cells
";
        let xtra = XtraData::new("tests::invalid_sample_id_chars");
        let res = MultiConfigCsv::from_reader(csv.as_bytes(), xtra);
        assert_snapshot!(res.unwrap_err());
    }

    #[test]
    fn invalid_sample_id_chars_2() {
        let csv = r"
[gene-expression]
ref,GRCh38-2020-A-chr21
create-bam,true

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types
bamtofastq,fastqs/cellranger/multi/VDJ_GEX_small_multi/small_gex_fastqs_chr21_new,any,gex_1,gene expression
pbmc_1k_protein_v3_antibody,../dui_tests/test_resources/cellranger-count/pbmc_1k_protein_v3_antibody,any,cmo1,Multiplexing Capture

[samples]
sample_id,cmo_ids,description
IAmABad#SampleId,CMO1,some cells
";
        let xtra = XtraData::new("tests::invalid_sample_id_chars_2");
        let res = MultiConfigCsv::from_reader(csv.as_bytes(), xtra);
        assert_snapshot!(res.unwrap_err());
    }

    #[test]
    fn too_long_sample_id() {
        let csv = r"
[gene-expression]
ref,GRCh38-2020-A-chr21
create-bam,true

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types
bamtofastq,fastqs/cellranger/multi/VDJ_GEX_small_multi/small_gex_fastqs_chr21_new,any,gex_1,gene expression
pbmc_1k_protein_v3_antibody,../dui_tests/test_resources/cellranger-count/pbmc_1k_protein_v3_antibody,any,cmo1,Multiplexing Capture

[samples]
sample_id,cmo_ids,description
IAmToooooooooooooooooooooooooooooooooooooooooooooooooooooooooLong,CMO1,some cells
";
        let xtra = XtraData::new("tests::too_long_sample_id");
        let res = MultiConfigCsv::from_reader(csv.as_bytes(), xtra);
        assert_snapshot!(res.unwrap_err());
    }

    #[test]
    fn load_trailing_whitespace() -> Result<()> {
        let csv = r"
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
sample_id,cmo_ids,description
HUMAN_T,CMO1,some cells
";
        let xtra = XtraData::new("tests::load_trailing_whitespace");
        let exp = MultiConfigCsv::from_reader(csv.as_bytes(), xtra)?;
        println!("{exp:?}");
        Ok(())
    }

    #[test]
    fn test_check_library_compatibility() -> Result<()> {
        let csv = r"
[gene-expression]
check_library_compatibility,false
reference,GRCh38-2020-A-chr21
create-bam,true

[feature]
ref,test/feature/cmo_features.csv

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types
bamtofastq,fastqs/cellranger/multi/VDJ_GEX_small_multi/small_gex_fastqs_chr21_new,any,gex_1,gene expression
";
        let xtra = XtraData::new("tests::test_check_library_compatibility");
        let exp = MultiConfigCsv::from_reader(csv.as_bytes(), xtra)?;
        println!("{exp:?}");
        assert!(!exp.gene_expression.unwrap().check_library_compatibility);
        let csv2 = r"
[gene-expression]
reference,GRCh38-2020-A-chr21
create-bam,true

[feature]
ref,test/feature/cmo_features.csv

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types
bamtofastq,fastqs/cellranger/multi/VDJ_GEX_small_multi/small_gex_fastqs_chr21_new,any,gex_1,gene expression

";
        let xtra2 = XtraData::new("tests::test_check_library_compatibility");
        let exp2 = MultiConfigCsv::from_reader(csv2.as_bytes(), xtra2)?;
        assert!(exp2.gene_expression.unwrap().check_library_compatibility);
        let csv3 = r"
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
";
        let xtra3 = XtraData::new("tests::test_check_library_compatibility");
        let exp3 = MultiConfigCsv::from_reader(csv3.as_bytes(), xtra3)?;
        assert!(!exp3.gene_expression.unwrap().check_library_compatibility);
        Ok(())
    }

    #[test]
    fn no_physical_library_id() -> Result<()> {
        let csv = r"
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
fastq_path,sample_indices,lanes,feature_types,subsample_rate
/path/to/fastqs,SI-001,1-3|5,Gene Expression,0.5
/path/to/fastqs,SI-002,4,Multiplexing Capture,1

[samples]
sample_id,cmo_ids
sample,CMO1
";
        let xtra = XtraData::new("tests::no_physical_library_id");
        let exp = MultiConfigCsv::from_reader(csv.as_bytes(), xtra)?;
        println!("{exp:?}");
        Ok(())
    }

    #[test]
    fn test_vdj_typo() -> Result<()> {
        let csv = r"
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
";
        let xtra = XtraData::new("test::test_vdj_typo");
        let exp = MultiConfigCsv::from_reader(csv.as_bytes(), xtra)?;
        println!("{exp:?}");
        Ok(())
    }

    /// Ensure that extra commas in the multi CSV are handled nicely
    #[test]
    fn test_extra_commas() -> Result<()> {
        let csv = r"
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
";

        let xtra = XtraData::new("test::extra_commas");
        let exp = MultiConfigCsv::from_reader(csv.as_bytes(), xtra)?;
        println!("{exp:?}");
        Ok(())
    }

    #[test]
    fn test_autogen_phys_lib_id_two_vdj_libs() {
        let csv = r"
[vdj]
ref,vdj/vdj_GRCh38_alts_ensembl-4.0.0

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
vdj_t,fastqs/cellranger/multi/vdj_t_tiny,any,,vdj,
vdj_ig,fastqs/cellranger/multi/vdj_ig_tiny,any,,vdj,
";
        let xtra = XtraData::new("test::test_autogen_phys_lib_id_two_vdj_libs");
        let res = MultiConfigCsv::from_reader(csv.as_bytes(), xtra);
        assert_snapshot!(res.unwrap_err());
    }

    /// Ensure that antigen capture is accompanied with a VDJ library
    #[test]
    fn test_antigen_requires_vdj() {
        let csv = r"
    [gene-expression]
    ref,mm10-2020-A-chr19
    create-bam,true

    [feature]
    ref,cellranger/multi/feature_refs/20211122_v1.1.csv

    [libraries]
    fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
    tiny_gex,fastqs/cellranger/multi/1245140_gex_vdj_beam_ab/gex,any,gex,gene expression,0.5
    tiny_an_b,fastqs/cellranger/multi/1245140_gex_vdj_beam_ab/an_b,any,ab,antigen capture,0.5
    ";

        let xtra = XtraData::new("test::antigen_required_vdj");
        let res = MultiConfigCsv::from_reader(csv.as_bytes(), xtra);
        assert_snapshot!(res.unwrap_err());
    }

    #[test]
    fn test_antigen_specificity_requires_ag() {
        let csv = r"
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
    ";

        let xtra = XtraData::new("test::antigen_specificity_requires_ag");
        let res = MultiConfigCsv::from_reader(csv.as_bytes(), xtra);
        assert_snapshot!(res.unwrap_err());
    }

    #[test]
    fn test_antigen_specificity_multiple_control_ids() {
        let csv = r"
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
    ";

        let xtra = XtraData::new("test::antigen_specificity_multiple_control_ids");
        let res = MultiConfigCsv::from_reader(csv.as_bytes(), xtra);
        assert_snapshot!(res.unwrap_err());
    }

    #[test]
    fn test_antigen_specificity_missing_allele() {
        let csv = r"
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
    ";

        let xtra = XtraData::new("test::antigen_specificity_missing_allele");
        let res = MultiConfigCsv::from_reader(csv.as_bytes(), xtra);
        assert_snapshot!(res.unwrap_err());
    }

    #[test]
    fn test_invalid_mhc_allele_character() {
        let csv = r"
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
    ";

        let xtra = XtraData::new("test::invalid_mhc_allele_character");
        let res = MultiConfigCsv::from_reader(csv.as_bytes(), xtra);
        assert_snapshot!(res.unwrap_err());
    }

    #[test]
    fn test_functional_map_requires_feature() {
        let csv = r"
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
    ";

        let xtra = XtraData::new("test::test_functional_map_requires_feature");
        let res = MultiConfigCsv::from_reader(csv.as_bytes(), xtra);
        assert_snapshot!(res.unwrap_err());
    }

    #[test]
    fn test_invalid_functional_map_value() {
        // Duplicate feature ids
        let csv = r"
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
    ";

        let xtra = XtraData::new("test::test_invalid_functional_map_value");
        let res = MultiConfigCsv::from_reader(csv.as_bytes(), xtra);
        assert_snapshot!(res.unwrap_err());

        // Duplicate functional name
        let csv = r"
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
    ";

        let xtra = XtraData::new("test::test_invalid_functional_map_value");
        let res = MultiConfigCsv::from_reader(csv.as_bytes(), xtra);
        assert_snapshot!(res.unwrap_err());
    }

    #[test]
    fn test_blank_lines() -> Result<()> {
        let csv = r"
    [gene-expression]
    ref,mm10-2020-A-chr19
    create-bam,true

    probe-set,/path/to/probe_set.csv

    [libraries]
    fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
    mygex,/path/to/fastqs,any,gex,gene expression,0.5

    mycmo,/path/to/fastqs,any,cmo,Multiplexing Capture,

    [samples]
    sample_id,cmo_ids
    sample1,1

    sample2,2
    ";

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

        let lib_count = res.libraries.libs.len();
        assert_eq!(lib_count, 2, "expected 2 libraries, found {lib_count}");

        let sample_count = res.samples.expect("samples section not present").rows.len();
        assert_eq!(sample_count, 2);

        Ok(())
    }

    #[test]
    fn test_multiple_probe_sets() -> Result<()> {
        let csv = r"
    [gene-expression]
    ref,mm10-2020-A-chr19
    create-bam,true
    probe-set,/path/to/probe_set1.csv
    probe-set,/path/to/probe_set2.csv

    [libraries]
    fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
    mygex,/path/to/fastqs,any,gex,gene expression,0.5

    mycmo,/path/to/fastqs,any,cmo,Multiplexing Capture,

    [samples]
    sample_id,cmo_ids
    sample1,1

    sample2,2
    ";

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

        let lib_count = res.libraries.libs.len();
        assert_eq!(lib_count, 2, "expected 2 libraries, found {lib_count}");

        let sample_count = res.samples.expect("samples section not present").rows.len();
        assert_eq!(sample_count, 2);

        Ok(())
    }

    #[test]
    fn test_emptydrops_per_sample() -> Result<()> {
        let csv = r"
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
    ";

        let xtra = XtraData::new("test::");
        let res = MultiConfigCsv::from_reader(csv.as_bytes(), xtra)?;
        res.samples
            .expect("samples section not present")
            .rows
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
        let csv = r"
    [gene-expression]
    ref,mm10-2020-A-chr19
    emptydrops-minimum-umis,75
    create-bam,true

    [libraries]
    fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
    tiny_gex,fastqs/cellranger/multi/1245140_gex_vdj_beam_ab/gex,any,gex,gene expression,0.5
    ";

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

    #[test]
    fn test_cell_calling_mode() {
        let config = assert_parse_ok(
            r"
            [gene-expression]
            ref,GRCh38-2024-A
            create-bam,true
            cell-calling-mode,gradient

            [libraries]
            fastq_id,fastqs,feature_types
            tiny_gex,/path/to/fastqs,Gene Expression",
        );
        assert_eq!(
            config.gene_expression.unwrap().cell_calling_mode,
            Some(CellCallingMode::Gradient),
        );
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
    fn test_plate_range_expansion_with_pipe() -> Result<()> {
        let csv = r"
    [gene-expression]
    ref,ref_path
    create-bam, false

    [libraries]
    fastq_id,fastqs,feature_types
    dummy_id,dummy_path,gene expression

    [samples]
    sample_id,cmx_barcode_ids
    sample1,QRT01A|QRT02B-02D
        ";

        let xtra = XtraData::new("test");
        let res = MultiConfigCsv::from_reader(csv.as_bytes(), xtra)?;
        let samples = res.samples.as_ref().unwrap();
        assert_eq!(
            samples.rows[0].cmx_barcode_ids.clone().unwrap(),
            vec!["QRT01A", "QRT02B", "QRT02C", "QRT02D"]
        );

        Ok(())
    }
    #[test]
    fn test_unpack_probe_barcode_ids() -> Result<()> {
        let csv = r"
    [gene-expression]
    ref,mm10-2020-A-chr19
    probe-set,/path/to/probe_set.csv
    create-bam,true

    [libraries]
    fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
    mygex,/path/to/fastqs,any,gex,gene expression,0.5

    [samples]
    sample_id,probe_barcode_ids
    sample1,BC1|BC2+BC3
    sample2,BC4|BC5+BC6+BC7
    ";

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
            &samples.rows[0],
            ProbeBarcodeIterationMode::All,
            vec!["BC1", "BC2", "BC3"],
        );
        expect_unpack_bcs(
            &samples.rows[0],
            ProbeBarcodeIterationMode::Mapped,
            vec!["BC1", "BC2"],
        );
        expect_unpack_bcs(
            &samples.rows[1],
            ProbeBarcodeIterationMode::All,
            vec!["BC4", "BC5", "BC6", "BC7"],
        );
        expect_unpack_bcs(
            &samples.rows[1],
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
        let csv = r"
    [gene-expression]
    ref,mm10-2020-A-chr19
    probe-set,/path/to/probe_set.csv
    create-bam,true

    [libraries]
    fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
    mygex,/path/to/fastqs,any,gex,gene expression,0.5

    [samples]
    sample_id,probe_barcode_ids
    sample1,BC1|BC2
    sample2,BC3
    ";

        let xtra = XtraData::new("test::");
        let res = MultiConfigCsv::from_reader(csv.as_bytes(), xtra)?;
        assert!(
            res.samples
                .expect("samples section not present")
                .get_translated_probe_barcodes()
                .is_empty()
        );
        Ok(())
    }

    #[test]
    fn test_multi_graph_detected_probe_barcode_pairing() -> Result<()> {
        let csv = r"
    [gene-expression]
    ref,mm10-2020-A-chr19
    probe-set,/path/to/probe_set.csv
    create-bam,true

    [libraries]
    fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
    mygex,/path/to/fastqs,any,gex,gene expression,0.5

    [samples]
    sample_id,probe_barcode_ids
    sample1,BC1|BC2
    sample2,BC3
    ";

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
            r"
            [gene-expression]
            ref,mm10-2020-A-chr19
            chemistry,mfrp
            create-bam,true

            [libraries]
            fastq_id,fastqs,lanes,physical_library_id,feature_types,chemistry
            mygex,/path/to/fastqs,any,gex,gene expression,mfrp-rna
            ",
            "failed to parse CSV: chemistry specified in both the [gene-expression] and [libraries] sections",
        );
    }

    #[test]
    fn test_cannot_use_auto_lib_chem() {
        assert_parse_fails(
            r"
            [gene-expression]
            ref,mm10-2020-A-chr19
            create-bam,true

            [libraries]
            fastq_id,fastqs,lanes,physical_library_id,feature_types,chemistry
            mygex,/path/to/fastqs,any,gex,gene expression,auto
            ",
            "[libraries] Specifying auto chemistries at the library level is not supported: (auto).",
        );
    }

    #[test]
    fn test_cannot_provide_partial_lib_chems() {
        assert_parse_fails(
            r"
            [gene-expression]
            ref,mm10-2020-A-chr19
            create-bam,true

            [libraries]
            fastq_id,fastqs,lanes,physical_library_id,feature_types,chemistry
            mygex,/path/to/fastqs1,any,gex,gene expression,mfrp-rna
            myab,/path/to/fastqs2,any,ab,antibody capture,
            ",
            "[libraries] A chemistry name must be provided for all libraries or none.",
        );
    }

    #[test]
    fn test_only_flex_chems_in_libs() {
        assert_parse_fails(
            r"
            [gene-expression]
            ref,mm10-2020-A-chr19
            create-bam,true

            [libraries]
            fastq_id,fastqs,lanes,physical_library_id,feature_types,chemistry
            mygex,/path/to/fastqs1,any,gex,gene expression,SC3Pv1
            myab,/path/to/fastqs2,any,ab,antibody capture,SC3Pv1
            ",
            "[libraries] Only Flex assays may specify chemistry at the per-library level; invalid chemistries: SC3Pv1.",
        );
    }

    #[test]
    fn test_crispr_probe_barcode_position() {
        assert_parse_fails(
            r"
            [gene-expression]
            ref,mm10-2020-A-chr19
            create-bam,false
            probe-set,/path/to/probe_set.csv

            [feature]
            crispr-probe-barcode-location,1
            reference,/path/to/reference

            [libraries]
            fastq_id,fastqs,lanes,physical_library_id,feature_types
            mygex,/path/to/fastqs1,any,gex,gene expression
            mycrispr,/path/to/fastqs2,any,ab,crispr guide capture

            [samples]
            sample_id,probe_barcode_ids
            sample1,BC1|BC2
            sample2,BC3
            ",
            "`crispr-probe-barcode-location` may only be used with manually specified RTL chemistries",
        );
    }

    #[test]
    fn test_no_conflicting_chems_for_lib_type() {
        assert_parse_fails(
            r"
            [gene-expression]
            ref,mm10-2020-A-chr19
            create-bam,true

            [libraries]
            fastq_id,fastqs,lanes,physical_library_id,feature_types,chemistry
            mygex,/path/to/fastqs1,any,gex,gene expression,MFRP-RNA
            mygex1,/path/to/fastqs2,any,gex,gene expression,SFRP
            ",
            "[libraries] Conflicting chemistry for Gene Expression libraries (MFRP-RNA, SFRP); manual chemistry must be the same for all libraries of the same type.",
        );
    }

    #[test]
    fn test_frac_parameters() {
        assert_parse_fails(
            r"
            [gene-expression]
            ref,mm10-2020-A-chr19
            create-bam,true
            min-major-probe-bc-frac,10

            [libraries]
            fastq_id,fastqs,lanes,physical_library_id,feature_types
            mygex1,/path/to/fastqs2,any,gex,gene expression
            ",
            "[gene-expression] has invalid min-major-probe-bc-frac '10' at line: 5, col: 37: must be in the range [0, 1]",
        );

        let csv = assert_parse_ok(
            r"
            [gene-expression]
            ref,mm10-2020-A-chr19
            create-bam,true
            min-major-probe-bc-frac,0.3

            [libraries]
            fastq_id,fastqs,lanes,physical_library_id,feature_types
            mygex1,/path/to/fastqs2,any,gex,gene expression
            ",
        );
        assert_eq!(
            Some(0.3),
            csv.gene_expression.unwrap().min_major_probe_bc_frac
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
        assert_snapshot!(
            MultiConfigCsv::from_csv("test/invalid_csvs/expect_cells_and_force_cells_samples.csv")
                .unwrap_err()
        );
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
        assert_snapshot!(
            MultiConfigCsv::from_csv(
                "test/invalid_csvs/mfrp_chemistry_with_missing_probe_barcode_entry.csv"
            )
            .unwrap_err()
        );
    }

    #[test]
    fn test_use_mfrp_chemistry_set() -> Result<()> {
        let csv = r"
    [gene-expression]
    ref,mm10-2020-A-chr19
    probe-set,/path/to/probe_set.csv
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
    ";

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
        let csv = r"
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
        ";

        assert_snapshot!(
            MultiConfigCsv::from_reader(csv.as_bytes(), XtraData::new("test::")).unwrap_err()
        );
    }

    #[test]
    fn test_multi_graph_hashtag_fingerprint() -> Result<()> {
        let csv = r"
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
    ";

        let xtra = XtraData::new("test::");
        let res = MultiConfigCsv::from_reader(csv.as_bytes(), xtra)?;
        let sample1_fprint = vec![
            Fingerprint::Tagged {
                tag_name: "CD4".to_string(),
                barcode_multiplexing_type: BarcodeMultiplexingType::CellLevel(CellLevel::Hashtag),
                translated_tag_names: vec![],
            },
            Fingerprint::Tagged {
                tag_name: "CD8a".to_string(),
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

    #[test]
    fn test_vdj_denovo_missing_chaintype() {
        let csv = r"
        [vdj]
        denovo, true
        inner-enrichment-primers,path/to/primers.txt

        [libraries]
        fastq_id,fastqs,physical_library_id,feature_types
        myvdj,/path/to/fastqs,vdj,VDJ
        ";

        assert_snapshot!(
            MultiConfigCsv::from_reader(csv.as_bytes(), XtraData::new("test::")).unwrap_err()
        );
    }

    #[test]
    fn test_any_platform_libs_req_headers() {
        assert_parse_ok(
            r"
[gene-expression]
ref,/path/to/gex/ref
create-bam,true

[libraries]
platform,read_path,feature_types
ont,/path/to/bams,Gene Expression",
        );
        assert_parse_ok(
            r"
[gene-expression]
ref,/path/to/gex/ref
create-bam,true

[libraries]
platform,read_path,feature_types
pb,/path/to/bams,Gene Expression",
        );
        assert_parse_ok(
            r"
[gene-expression]
ref,/path/to/gex/ref
create-bam,true

[libraries]
platform,read_path,feature_types, fastq_id
ilmn,/path/to/fastqs,Gene Expression,mygex",
        );
    }

    #[test]
    fn test_any_platform_libs_opt_headers() {
        let csv = assert_parse_ok(
            r"
[gene-expression]
ref,/path/to/gex/ref
create-bam,true

[libraries]
platform,read_path,feature_types,fastq_id,physical_library_id,lanes,subsample_rate,chemistry
ont,/path/to/bams,Gene Expression,mygex,a,any,0.5,",
        );
        assert_eq!(
            csv.libraries.libs.first().unwrap(),
            &Library {
                read_path: PathBuf::from("/path/to/bams"),
                physical_library_id: "a".to_string(),
                feature_type: LibraryType::Gex,
                subsample_rate: Some(SubsampleRate(0.5)),
                chemistry: None,
                platform: SequencingPlatform::OxfordNanopore,
                data_type: DataType::LongRead {
                    format: SequencingDataFormat::Bam
                },
            }
        );
        let csv = assert_parse_ok(
            r"
[gene-expression]
ref,/path/to/gex/ref
create-bam,true

[libraries]
platform,read_path,feature_types,fastq_id,physical_library_id,lanes,subsample_rate,chemistry
pb,/path/to/bams,Gene Expression,mygex,a,any,0.5,",
        );
        assert_eq!(
            csv.libraries.libs.first().unwrap(),
            &Library {
                read_path: PathBuf::from("/path/to/bams"),
                physical_library_id: "a".to_string(),
                feature_type: LibraryType::Gex,
                subsample_rate: Some(SubsampleRate(0.5)),
                chemistry: None,
                platform: SequencingPlatform::PacBio,
                data_type: DataType::LongRead {
                    format: SequencingDataFormat::Bam
                },
            }
        );
        let csv = assert_parse_ok(
            r"
[gene-expression]
ref,/path/to/gex/ref
create-bam,true

[libraries]
platform,read_path,feature_types,fastq_id,physical_library_id,lanes,subsample_rate,chemistry
ilmn,/path/to/fastqs,Gene Expression,mygex,a,any,0.5,",
        );
        assert_eq!(
            csv.libraries.libs.first().unwrap(),
            &Library {
                read_path: PathBuf::from("/path/to/fastqs"),
                physical_library_id: "a".to_string(),
                feature_type: LibraryType::Gex,
                subsample_rate: Some(SubsampleRate(0.5)),
                chemistry: None,
                platform: SequencingPlatform::Illumina,
                data_type: DataType::Bcl2Fastq {
                    fastq_id: "mygex".to_string(),
                    lanes: Lanes::Any,
                    ilmn_fastq_id: None
                },
            }
        );
    }

    #[test]
    fn test_gex_multi_platforms() {
        assert_parse_fails(
            r"
[gene-expression]
ref,/path/to/gex/ref
create-bam,true

[feature]
ref,path/to/feature_ref.csv

[libraries]
platform,read_path,feature_types,fastq_id
ont,/path/to/bams,Gene Expression,
ilmn,/path/to/fastqs,Gene Expression,mygex",
            "[libraries] All Gene Expression libraries must use the same sequencing platform.",
        );
    }

    #[test]
    fn test_non_gex_lr_platforms() {
        assert_parse_fails(
            r"
[gene-expression]
ref,/path/to/gex/ref
create-bam,true

[feature]
ref,path/to/feature_ref.csv

[libraries]
platform,read_path,feature_types, fastq_id
ont,/path/to/bams,Gene Expression,
ont,/path/to/fastqs,Antibody Capture,",
            "[libraries] Only Gene Expression libraries are supported with ONT or PacBio platforms.",
        );
    }

    #[test]
    fn test_platform_libs_mix() {
        assert_parse_ok(
            r"
[gene-expression]
ref,/path/to/gex/ref
create-bam,true

[feature]
ref,path/to/feature_ref.csv

[libraries]
platform,read_path,feature_types, fastq_id
ont,/path/to/bams,Gene Expression,
ilmn,/path/to/fastqs,Antibody Capture,myab",
        );
    }
    #[test]
    fn test_minimap_preset_aligner_conflict() {
        let csv = r"
        [gene-expression]
        aligner,star
        minimap-preset,splice
        create-bam,True
        reference-path,/path/to/reference

        [libraries]
        fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate,chemistry
        ONT,/path/to/fastqs,1,lib1,Gene Expression,,
        ";
        assert_snapshot!(
            MultiConfigCsv::from_reader(csv.as_bytes(), XtraData::new("test::")).unwrap_err()
        );
    }
    #[test]
    fn test_invalid_minimap_preset() {
        let csv = r"
        [gene-expression]
        aligner,minimap2
        minimap-preset,map-ont
        create-bam,True
        reference-path,/path/to/reference

        [libraries]
        fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate,chemistry
        ONT,/path/to/fastqs,1,lib1,Gene Expression,,
        ";
        assert_snapshot!(
            MultiConfigCsv::from_reader(csv.as_bytes(), XtraData::new("test::")).unwrap_err()
        );
    }

    #[test]
    fn test_whitelists_are_disjoint() -> Result<()> {
        let barcode_ids = MultiplexingBarcodeIds::load()?;
        let sets = [
            ("RTL", &barcode_ids.rtl_barcode_ids),
            ("OCM", &barcode_ids.ocm_barcode_ids),
            ("CMX", &barcode_ids.cmx_barcode_ids),
        ];
        for i in 0..sets.len() {
            for j in (i + 1)..sets.len() {
                let (name_i, set_i) = sets[i];
                let (name_j, set_j) = sets[j];
                let intersection: Vec<_> = set_i.intersection(set_j).collect();
                assert!(
                    intersection.is_empty(),
                    "{name_i} and {name_j} whitelists share {} barcode IDs: {intersection:?}",
                    intersection.len()
                );
            }
        }

        Ok(())
    }
}
