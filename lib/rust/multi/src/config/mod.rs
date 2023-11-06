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
use anyhow::{anyhow, bail, ensure, Context, Result};
use barcode::whitelist::{find_whitelist, BarcodeId};
use barcode::{WhitelistSource, WhitelistSpec};
use cr_types::reference::feature_reference::{
    BeamMode, FeatureConfig, SpecificityControls, MHC_ALLELE, NO_ALLELE,
};
use cr_types::rna_read::LegacyLibraryType;
use cr_types::sample_def::{FastqMode, SampleDef};
use cr_types::types::{CellMultiplexingType, CrMultiGraph};
use cr_types::{AlignerParam, TargetingMethod};
use fastq_set::filenames::FastqDef;
use itertools::Itertools;
use martian_derive::martian_filetype;
use metric::{TxHashMap, TxHashSet};
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::convert::{AsRef, TryFrom};
use std::fmt::{Display, Formatter};
use std::fs::File;
use std::hash::Hash;
use std::io::{BufReader, Read};
use std::iter::FromIterator;
use std::path::{Path, PathBuf};
use std::str::FromStr;

const MIN_FORCE_CELLS: usize = 10;

const ERROR_INCLUDE_INTRONS_WITH_RTL: &str = "The [gene-expression] section specifies the parameter include-introns, which is not valid for Fixed RNA Profiling chemistries.";

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
        let t = Ident::from_str(s)
            .map_err(|err| anyhow!("{err} and be no more than {MAX_ID_LEN} characters long"))?;
        if t.0.len() > MAX_ID_LEN {
            bail!("must be no more {MAX_ID_LEN} characters long");
        }
        Ok(SampleId(t.0))
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
        let v = s.parse::<usize>()?;
        if v >= 1 {
            Ok(AtLeastOne(v))
        } else {
            bail!("must be >= 1")
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Probability(pub f64);

impl FromStr for Probability {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self> {
        let p = s.parse::<f64>()?;
        if (0.0..=1.0).contains(&p) {
            Ok(Probability(p))
        } else {
            bail!("must be in [0, 1]")
        }
    }
}

// TODO: unify this with cr_types/src/chemistry/mod.rs and
//   cr_wrap/src/chemistry_arg.rs
/// A enum for the different possible ways to specify gene-expression chemistry
#[allow(non_camel_case_types)]
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, Serialize, Deserialize)]
pub enum ChemistryParam {
    #[serde(rename = "auto")]
    #[default]
    Auto,
    #[serde(rename = "custom")]
    Custom,
    #[serde(rename = "threeprime")]
    ThreePrime,
    #[serde(rename = "fiveprime")]
    FivePrime,
    SC3Pv1,
    SC3Pv2,
    SC3Pv3,
    SC3Pv3LT,
    SC3Pv3HT,
    #[serde(rename = "SC5P-PE")]
    SC5P_PE,
    #[serde(rename = "SC5P-R2")]
    SC5P_R2,
    SC5PHT,
    #[serde(rename = "SC-FB")]
    SC_FB,
    SFRP, // singleplex fixed RNA profiling
    MFRP, // multiplex fixed RNA profiling
    #[serde(rename = "MFRP-47")]
    MFRP_47, // multiplex fixed RNA profiling with 47 probe barcodes
    #[serde(rename = "MFRP-uncollapsed")]
    MFRP_uncollapsed, // multiplex fixed RNA profiling with 64 non-based-balanced probe barcodes
    #[serde(rename = "MFRP-R1")]
    MFRP_R1, // multiplex fixed RNA profiling (probeBC on R1)
    #[serde(rename = "MFRP-R1-48-uncollapsed")]
    MFRP_R1_48_uncollapsed, // multiplex fixed RNA profiling (probeBC on R1) with 192 non-based-balanced probe barcodes
    ARCv1, // Multiome GEX chemistry
}

impl ChemistryParam {
    /// Return whether the chemistry is RTL, and None when the chemistry is auto or custom.
    fn is_rtl(&self) -> Option<bool> {
        #[allow(clippy::enum_glob_use)]
        use ChemistryParam::*;
        match self {
            Auto | Custom => None,
            SFRP | MFRP | MFRP_47 | MFRP_uncollapsed | MFRP_R1 | MFRP_R1_48_uncollapsed => {
                Some(true)
            }
            ARCv1 | FivePrime | SC3Pv1 | SC3Pv2 | SC3Pv3 | SC3Pv3HT | SC3Pv3LT | SC5PHT
            | SC5P_PE | SC5P_R2 | SC_FB | ThreePrime => Some(false),
        }
    }

    /// Return Some if this chemistry is explicit, or None if it is auto.
    pub fn name(&self) -> Option<&str> {
        #[allow(clippy::enum_glob_use)]
        use ChemistryParam::*;
        match self {
            Auto => None,
            Custom => Some("custom"),
            ThreePrime => Some("threeprime"),
            FivePrime => Some("fiveprime"),
            SC3Pv1 => Some("SC3Pv1"),
            SC3Pv2 => Some("SC3Pv2"),
            SC3Pv3 => Some("SC3Pv3"),
            SC3Pv3LT => Some("SC3Pv3LT"),
            SC3Pv3HT => Some("SC3Pv3HT"),
            SC5P_PE => Some("SC5P-PE"),
            SC5P_R2 => Some("SC5P-R2"),
            SC5PHT => Some("SC5PHT"),
            SC_FB => Some("SC-FB"),
            SFRP => Some("SFRP"),
            MFRP => Some("MFRP"),
            MFRP_47 => Some("MFRP-47"),
            MFRP_uncollapsed => Some("MFRP-uncollapsed"),
            MFRP_R1 => Some("MFRP-R1"),
            MFRP_R1_48_uncollapsed => Some("MFRP-R1-48-uncollapsed"),
            ARCv1 => Some("ARC-v1"),
        }
    }
}

impl FromStr for ChemistryParam {
    type Err = anyhow::Error;

    fn from_str(chemistry: &str) -> Result<Self> {
        #[allow(clippy::enum_glob_use)]
        use ChemistryParam::*;
        match chemistry.to_ascii_lowercase().as_str() {
            "auto" => Ok(Auto),
            "custom" => Ok(Custom),
            "threeprime" => Ok(ThreePrime),
            "fiveprime" => Ok(FivePrime),
            "sc3pv1" => Ok(SC3Pv1),
            "sc3pv2" => Ok(SC3Pv2),
            "sc3pv3" => Ok(SC3Pv3),
            "sc3pv3lt" => Ok(SC3Pv3LT),
            "sc3pv3ht" => Ok(SC3Pv3HT),
            "sc5p-pe" | "sc5ppe" => Ok(SC5P_PE),
            "sc5p-r2" | "sc5pr2" => Ok(SC5P_R2),
            "sc5pht" => Ok(SC5PHT),
            "sc-fb" | "scfb" => Ok(SC_FB),
            "sfrp" => Ok(SFRP),
            "mfrp" => Ok(MFRP),
            "mfrp-47" => Ok(MFRP_47),
            "mfrp-uncollapsed" => Ok(MFRP_uncollapsed),
            "mfrp-r1" => Ok(MFRP_R1),
            "mfrp-r1-48-uncollapsed" => Ok(MFRP_R1_48_uncollapsed),
            "arc-v1" => Ok(ARCv1),
            _ => bail!("unknown chemistry: {chemistry}"),
        }
    }
}

impl Display for ChemistryParam {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result<(), std::fmt::Error> {
        let mut json_str = serde_json::to_string(&self).unwrap();
        json_str.remove(0);
        json_str.pop();
        write!(f, "{json_str}")
    }
}

/// The gene-expression parameters in the experiment CSV
#[derive(Debug, Default, Serialize, Deserialize)]
pub struct GeneExpressionParams {
    pub reference_path: PathBuf,
    pub probe_set: Option<PathBuf>,
    pub emptydrops_minimum_umis: Option<usize>,
    pub r1_length: Option<usize>,
    pub r2_length: Option<usize>,
    pub chemistry: ChemistryParam,
    pub expect_cells: Option<usize>,
    pub force_cells: Option<usize>,
    pub no_secondary_analysis: bool,
    pub include_introns: bool,
    pub check_library_compatibility: bool,
    pub aligner: Option<AlignerParam>,
    pub no_bam: bool,
    pub filter_probes: Option<bool>,
    pub cmo_set: Option<PathBuf>,
    pub min_assignment_confidence: Option<f64>,
    pub barcode_sample_assignment: Option<PathBuf>,
}

impl GeneExpressionParams {
    /// Return the probe set.
    pub fn probe_set(&self) -> Option<&Path> {
        self.probe_set.as_deref()
    }

    /// Return the targeting_method.
    pub fn targeting_method(&self) -> Option<TargetingMethod> {
        self.probe_set
            .is_some()
            .then_some(TargetingMethod::TemplatedLigation)
    }

    pub fn has_force_cells(&self) -> bool {
        self.force_cells.is_some()
    }

    pub fn has_expect_cells(&self) -> bool {
        self.expect_cells.is_some()
    }

    pub fn invalid_parameter_with_antigen_capture(&self) -> Option<String> {
        if self.probe_set.is_some() {
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
            || self.chemistry.is_rtl().unwrap_or(false)
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
        let mut reference_path = PathBuf::new();
        let mut probe_set: Option<PathBuf> = None;
        let mut filter_probes: Option<bool> = None;
        let mut emptydrops_minimum_umis = None;
        let mut r1_length: Option<usize> = None;
        let mut r2_length: Option<usize> = None;
        let mut chemistry = ChemistryParam::default();
        let mut expect_cells: Option<usize> = None;
        let mut force_cells: Option<usize> = None;
        let mut no_secondary_analysis = false;
        let mut include_introns: Option<bool> = None;
        let mut check_library_compatibility = true;
        let mut aligner: Option<AlignerParam> = None;
        let mut no_bam = false;
        let mut cmo_set: Option<PathBuf> = None;
        let mut min_assignment_confidence: Option<f64> = None;
        let mut barcode_sample_assignment: Option<PathBuf> = None;
        for row in &sec.rows {
            if row.is_empty() {
                continue;
            }
            let param = row[0].fragment().to_ascii_lowercase();
            let ctx = ctx.with_col(param.as_str());
            match param.as_str() {
                "ref" | "reference" | "reference-path" => {
                    reference_path = row
                        .get(1)
                        .ok_or_else(|| {
                            anyhow!("{ctx} no value provided for '{}'", row[0].fragment())
                        })?
                        .parse::<PathBuf>(ctx)?;
                }
                "probe-set" => {
                    if let Some(path) = row.get(1).and_then(empty_is_none) {
                        probe_set = Some(path.parse::<PathBuf>(ctx)?);
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
                        chemistry = val.parse::<ChemistryParam>(ctx)?;
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
                        if parsed_val < MIN_FORCE_CELLS {
                            bail!(
                                "The 'force-cells' parameter specified under the '[gene-expression]' \
                                section needs to be at least {}. The value you have specified is {} \
                                which is too low.",
                                MIN_FORCE_CELLS,
                                parsed_val
                            )
                        }
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
                "nobam" | "no-bam" => {
                    if let Some(val) = row.get(1).and_then(empty_is_none) {
                        no_bam = val.parse::<Bool>(ctx)?.into();
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
                _ => {
                    bail!(
                        "{} unknown parameter '{}' provided at line: {}, col: {}",
                        ctx,
                        row[0].fragment(),
                        row[0].location_line(),
                        row[0].get_utf8_column(),
                    );
                }
            }
        }

        if reference_path.as_os_str().is_empty() {
            bail!("{} reference is missing", ctx);
        }
        if expect_cells.is_some() && force_cells.is_some() {
            bail!(
                "{} only one of force-cells or expect-cells is allowed.",
                ctx
            );
        }
        if filter_probes.is_some() && probe_set.is_none() {
            bail!("{} filter-probes requires a probe-set.", ctx);
        }
        if probe_set.is_some()
            && (cmo_set.is_some()
                | min_assignment_confidence.is_some()
                | barcode_sample_assignment.is_some())
        {
            bail!(
                "{} When probe-set is specified, cmo-set, min-assignment-confidence & barcode-sample-assignment are invalid parameters.",
                ctx
            );
        }

        // Disallow setting include-introns for RTL chemistries.
        if chemistry.is_rtl() == Some(true) {
            ensure!(include_introns.is_none(), ERROR_INCLUDE_INTRONS_WITH_RTL);
        }

        // Filter probes by default.
        if probe_set.is_some() && filter_probes.is_none() {
            filter_probes = Some(true);
        }

        Ok(GeneExpressionParams {
            reference_path,
            probe_set,
            filter_probes,
            emptydrops_minimum_umis,
            r1_length,
            r2_length,
            chemistry,
            expect_cells,
            force_cells,
            no_secondary_analysis,
            include_introns: include_introns.unwrap_or(true),
            check_library_compatibility,
            aligner,
            no_bam,
            cmo_set,
            min_assignment_confidence,
            barcode_sample_assignment,
        })
    }
}

#[derive(Debug, Serialize, Deserialize)]
pub struct FeatureParams {
    pub reference_path: Option<PathBuf>,
    pub r1_length: Option<usize>,
    pub r2_length: Option<usize>,
}

impl<'a> TryFrom<&Section<'a>> for FeatureParams {
    type Error = anyhow::Error;

    fn try_from(sec: &Section<'a>) -> Result<Self> {
        let ctx = ParseCtx::Hdr(sec.name);
        let mut reference_path: Option<PathBuf> = None;
        let mut r1_length: Option<usize> = None;
        let mut r2_length: Option<usize> = None;
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
                _ => {
                    bail!(
                        "{} unknown parameter '{}' provided at line: {}, col: {}",
                        ctx,
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
    pub r2_revcomp: Option<bool>,
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
        let mut r2_revcomp = None;
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
                "r2-revcomp" => {
                    if let Some(val) = row.get(1).and_then(empty_is_none) {
                        r2_revcomp = Some(val.parse::<Bool>(ctx)?.into());
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
        if reference_path.as_os_str().is_empty() {
            bail!("{} reference is missing", ctx);
        }
        Ok(VdjParams {
            reference_path,
            inner_enrichment_primers,
            r1_length,
            r2_length,
            multiplet_filter,
            shared_contig_filter,
            umi_baseline_filter,
            r2_revcomp,
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

#[allow(non_camel_case_types)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub enum FeatureType {
    #[serde(rename = "Gene Expression")]
    GeneExpression,
    VDJ,
    #[serde(rename = "VDJ-T")]
    VDJ_T,
    #[serde(rename = "VDJ-T-GD")]
    VDJ_T_GD,
    #[serde(rename = "VDJ-B")]
    VDJ_B,
    #[serde(rename = "Antibody Capture")]
    AntibodyCapture,
    #[serde(rename = "Multiplexing Capture")]
    MultiplexingCapture,
    #[serde(rename = "CRISPR Guide Capture")]
    CrisprGuideCapture,
    #[serde(rename = "Antigen Capture")]
    AntigenCapture,
    Custom,
}

impl Display for FeatureType {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        use FeatureType::{
            AntibodyCapture, AntigenCapture, CrisprGuideCapture, Custom, GeneExpression,
            MultiplexingCapture, VDJ, VDJ_B, VDJ_T, VDJ_T_GD,
        };
        write!(
            f,
            "{}",
            match self {
                GeneExpression => "Gene Expression",
                VDJ => "VDJ",
                VDJ_T => "VDJ-T",
                VDJ_T_GD => "VDJ-T-GD",
                VDJ_B => "VDJ-B",
                AntibodyCapture => "Antibody Capture",
                AntigenCapture => "Antigen Capture",
                MultiplexingCapture => "Multiplexing Capture",
                CrisprGuideCapture => "CRISPR Guide Capture",
                Custom => "Custom",
            }
        )
    }
}

fn legacy_library_type(fts: &[FeatureType]) -> Result<LegacyLibraryType> {
    use FeatureType::{
        AntibodyCapture, AntigenCapture, CrisprGuideCapture, Custom, GeneExpression,
        MultiplexingCapture, VDJ, VDJ_B, VDJ_T, VDJ_T_GD,
    };
    match fts {
        [GeneExpression] => Ok(LegacyLibraryType::GeneExpression),
        [VDJ] | [VDJ_T] | [VDJ_B] | [VDJ_T_GD] => Ok(LegacyLibraryType::Vdj),
        [Custom] => Ok(LegacyLibraryType::Custom),
        [CrisprGuideCapture] => Ok(LegacyLibraryType::CrisprGuideCapture),
        [AntibodyCapture] => Ok(LegacyLibraryType::AntibodyCapture),
        [MultiplexingCapture] => Ok(LegacyLibraryType::Multiplexing),
        [AntigenCapture] => Ok(LegacyLibraryType::AntigenCapture),
        _ => bail!("unimplemented feature_type: '{:?}'", fts),
    }
}

impl<'a> TryFrom<Span<'a>> for FeatureType {
    type Error = anyhow::Error;

    fn try_from(s: Span<'a>) -> Result<Self> {
        use FeatureType::{
            AntibodyCapture, AntigenCapture, CrisprGuideCapture, Custom, GeneExpression,
            MultiplexingCapture, VDJ, VDJ_B, VDJ_T, VDJ_T_GD,
        };
        let lcs = s.fragment().to_ascii_lowercase();
        match lcs.as_str() {
            "gene expression" => Ok(GeneExpression),
            "vdj" => Ok(VDJ),
            "vdj-t" => Ok(VDJ_T),
            "vdj-t-gd" => Ok(VDJ_T_GD),
            "vdj-b" => Ok(VDJ_B),
            "antibody capture" => Ok(AntibodyCapture),
            "multiplexing capture" => Ok(MultiplexingCapture),
            "crispr guide capture" => Ok(CrisprGuideCapture),
            "antigen capture" => Ok(AntigenCapture),
            "custom" => Ok(Custom),
            _ => bail!(
                "unknown feature_type '{}' at line: {}, col: {}",
                s.fragment(),
                s.location_line(),
                s.get_utf8_column()
            ),
        }
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
        if span.fragment().to_ascii_lowercase().as_str() == "any" {
            Ok(Lanes::Any)
        } else {
            let mut lanes: Vec<_> = parse_vec(span.clone())
                .map_err(|_: NomErr<'_>| {
                    anyhow!(
                        "{ctx} has invalid {} '{}' at line: {}, col: {}",
                        libsconst::LANES,
                        span.fragment(),
                        span.location_line(),
                        span.get_utf8_column(),
                    )
                })?
                .into_iter()
                .map(|l| Ok(parse_range::<usize>(l.clone())?.collect::<Vec<_>>()))
                .collect::<Result<Vec<_>>>()?
                .into_iter()
                .flatten()
                .collect();
            lanes.sort_unstable();
            lanes.dedup();
            Ok(Lanes::Lanes(lanes))
        }
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub enum Library {
    Bcl2Fastq {
        fastq_id: String,
        fastqs: PathBuf,
        lanes: Lanes,
        physical_library_id: String,
        feature_types: Vec<FeatureType>,
        gem_well: GemWell,
        subsample_rate: Option<f64>,
        ilmn_fastq_id: Option<String>, // Fastq ID might be edited to handle case of duplicate fastq IDs referring to different data
    },
    BclProcessor {
        fastq_path: PathBuf,
        sample_indices: Vec<String>,
        lanes: Lanes,
        physical_library_id: String,
        feature_types: Vec<FeatureType>,
        gem_well: GemWell,
        subsample_rate: Option<f64>,
    },
}

impl Library {
    pub fn lanes(&self) -> &Lanes {
        use Library::{Bcl2Fastq, BclProcessor};
        match self {
            Bcl2Fastq { lanes, .. } => lanes,
            BclProcessor { lanes, .. } => lanes,
        }
    }

    pub fn feature_types(&self) -> &[FeatureType] {
        use Library::{Bcl2Fastq, BclProcessor};
        match self {
            Bcl2Fastq { feature_types, .. } => feature_types.as_slice(),
            BclProcessor { feature_types, .. } => feature_types.as_slice(),
        }
    }

    pub fn gem_well(&self) -> GemWell {
        use Library::{Bcl2Fastq, BclProcessor};
        match self {
            Bcl2Fastq { gem_well, .. } => *gem_well,
            BclProcessor { gem_well, .. } => *gem_well,
        }
    }

    pub fn physical_library_id(&self) -> &str {
        use Library::{Bcl2Fastq, BclProcessor};
        match self {
            Bcl2Fastq {
                physical_library_id,
                ..
            } => physical_library_id.as_str(),
            BclProcessor {
                physical_library_id,
                ..
            } => physical_library_id.as_str(),
        }
    }

    pub fn is_count(&self) -> Result<bool> {
        use FeatureType::{
            AntibodyCapture, AntigenCapture, CrisprGuideCapture, Custom, GeneExpression,
            MultiplexingCapture, VDJ, VDJ_B, VDJ_T, VDJ_T_GD,
        };
        match self.feature_types() {
            [GeneExpression]
            | [AntibodyCapture]
            | [AntigenCapture]
            | [CrisprGuideCapture]
            | [Custom]
            | [MultiplexingCapture] => Ok(true),
            [VDJ] | [VDJ_T] | [VDJ_B] | [VDJ_T_GD] => Ok(false),
            _ => bail!("unimplemented feature_type: '{:?}'", self.feature_types()),
        }
    }

    pub fn is_multiplexing(&self) -> bool {
        use FeatureType::MultiplexingCapture;
        matches!(self.feature_types(), [MultiplexingCapture])
    }

    pub fn is_vdj(&self) -> bool {
        use FeatureType::{VDJ, VDJ_B, VDJ_T, VDJ_T_GD};
        matches!(self.feature_types(), [VDJ] | [VDJ_T] | [VDJ_B] | [VDJ_T_GD])
    }

    pub fn is_antigen(&self) -> bool {
        use FeatureType::AntigenCapture;
        matches!(self.feature_types(), [AntigenCapture])
    }

    pub fn is_gex(&self) -> bool {
        use FeatureType::GeneExpression;
        matches!(self.feature_types(), [GeneExpression])
    }

    pub fn is_crispr(&self) -> bool {
        use FeatureType::CrisprGuideCapture;
        matches!(self.feature_types(), [CrisprGuideCapture])
    }

    pub fn overlaps(&self, other: &Library) -> bool {
        use Library::{Bcl2Fastq, BclProcessor};
        match self {
            Bcl2Fastq {
                fastq_id: fastq_id1,
                fastqs: fastq_path1,
                lanes: lanes1,
                ilmn_fastq_id: ilmn_fastq_id1,
                ..
            } => {
                if let Bcl2Fastq {
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
            BclProcessor {
                fastq_path: fastq_path1,
                sample_indices,
                lanes: lanes1,
                ..
            } => {
                // something bogus is happening here requiring me to provide State type parameter
                let sample_indices1 = TxHashSet::from_iter(sample_indices);
                if let BclProcessor {
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

    pub fn to_fastq_def(&self, config: &MultiConfigCsv) -> Result<FastqDef> {
        self.to_sample_def(config)?.get_fastq_def()
    }

    pub fn to_sample_def(&self, config: &MultiConfigCsv) -> Result<SampleDef> {
        use LegacyLibraryType::GeneExpression;
        use Library::{Bcl2Fastq, BclProcessor};
        let lanes = self.lanes().clone().into();
        let library_type = Some(legacy_library_type(self.feature_types())?);
        let target_set = match library_type {
            Some(GeneExpression) => config
                .gene_expression
                .as_ref()
                .and_then(GeneExpressionParams::probe_set)
                .map(PathBuf::from),
            _ => None,
        };
        let target_set_name = target_set
            .as_ref()
            .map(|x| x.file_stem().unwrap().to_string_lossy().into_owned());
        match self {
            Bcl2Fastq {
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
                    target_set,
                    target_set_name,
                    fastq_id: Some(fastq_id.clone()),
                })
            }
            BclProcessor {
                fastq_path,
                sample_indices,
                gem_well,
                subsample_rate,
                ..
            } => Ok(SampleDef {
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
                target_set,
                target_set_name,
                fastq_id: None,
            }),
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
    #[cfg(not(test))]
    pub const LIBS_OPT_HDRS: &[&str] = &[PHYSICAL_LIBRARY_ID, LANES, SUBSAMPLE_RATE];
    #[cfg(test)]
    pub const LIBS_OPT_HDRS: &[&str] = &[PHYSICAL_LIBRARY_ID, LANES, GEM_WELL, SUBSAMPLE_RATE];

    pub const FASTQ_PATH: &str = "fastq_path";
    pub const SAMPLE_INDICES: &str = "sample_indices";
    pub const LIBS_INT_REQ_HDRS: &[&str] = &[FASTQ_PATH, SAMPLE_INDICES, FEATURE_TYPES];
    #[cfg(not(test))]
    pub const LIBS_INT_OPT_HDRS: &[&str] = &[PHYSICAL_LIBRARY_ID, LANES, SUBSAMPLE_RATE];
    #[cfg(test)]
    pub const LIBS_INT_OPT_HDRS: &[&str] = &[PHYSICAL_LIBRARY_ID, LANES, GEM_WELL, SUBSAMPLE_RATE];
}

// TODO: should these just be SampleDef?
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
#[serde(transparent)]
pub struct LibrariesCsv(pub Vec<Library>);

impl LibrariesCsv {
    /// Return an iterator over the feature types.
    pub fn feature_types(&self) -> impl Iterator<Item = FeatureType> + '_ {
        self.0.iter().flat_map(Library::feature_types).copied()
    }

    /// Return true if any libraries are gene expression.
    pub fn has_gene_expression(&self) -> bool {
        self.has_feature_type(FeatureType::GeneExpression)
    }

    /// Return true if any libraries are antibody capture.
    pub fn has_antibody_capture(&self) -> bool {
        self.has_feature_type(FeatureType::AntibodyCapture)
    }

    /// Return true if any libraries are Antibody/Crispr/Antigen/Custom.
    /// Return false for CMO multiplexing, which comes with a built-in list.
    pub fn has_feature_barcode(&self) -> bool {
        use FeatureType::{AntibodyCapture, AntigenCapture, CrisprGuideCapture, Custom};
        self.0.iter().flat_map(Library::feature_types).any(|ft| {
            matches!(
                ft,
                AntibodyCapture | CrisprGuideCapture | AntigenCapture | Custom
            )
        })
    }

    /// Return true if any libraries are of the provided feature type.
    pub fn has_feature_type(&self, feature_type: FeatureType) -> bool {
        self.0
            .iter()
            .flat_map(Library::feature_types)
            .any(|ft| *ft == feature_type)
    }

    /// Return true if any libraries are antigen capture.
    pub fn has_antigen_capture(&self) -> bool {
        self.has_feature_type(FeatureType::AntigenCapture)
    }

    /// Return true if any libraries are CMO multiplexing.
    pub fn has_multiplexing(&self) -> bool {
        self.0
            .iter()
            .flat_map(Library::feature_types)
            .any(|&ft| ft == FeatureType::MultiplexingCapture)
    }

    pub fn has_vdj(&self) -> bool {
        use FeatureType::{VDJ, VDJ_B, VDJ_T, VDJ_T_GD};
        self.0
            .iter()
            .flat_map(Library::feature_types)
            .any(|ft| matches!(ft, VDJ | VDJ_T | VDJ_B | VDJ_T_GD))
    }

    pub fn has_vdj_t_or_gd(&self) -> bool {
        use FeatureType::{VDJ_T, VDJ_T_GD};
        for lib in &self.0 {
            for feature_type in lib.feature_types() {
                match feature_type {
                    VDJ_T | VDJ_T_GD => return true,
                    _ => continue,
                }
            }
        }
        false
    }

    pub fn has_vdj_b(&self) -> bool {
        use FeatureType::VDJ_B;
        for lib in &self.0 {
            for feature_type in lib.feature_types() {
                match feature_type {
                    VDJ_B => return true,
                    _ => continue,
                }
            }
        }
        false
    }

    pub fn number_of_vdj(&self) -> usize {
        use FeatureType::{VDJ, VDJ_B, VDJ_T, VDJ_T_GD};
        let mut uniq_feature_types: HashSet<&FeatureType> = HashSet::new();
        for lib in &self.0 {
            for feature_type in lib.feature_types() {
                match feature_type {
                    VDJ | VDJ_T | VDJ_B | VDJ_T_GD => _ = uniq_feature_types.insert(feature_type),
                    _ => continue,
                }
            }
        }
        uniq_feature_types.len()
    }

    pub fn beam_mode(&self) -> Option<BeamMode> {
        if self.has_vdj_t_or_gd() {
            return Some(BeamMode::BeamT);
        } else if self.has_vdj_b() {
            return Some(BeamMode::BeamAB);
        }
        None
    }
}

fn feature_types_short(ft: &[FeatureType]) -> Result<&'static str> {
    use FeatureType::{
        AntibodyCapture, AntigenCapture, CrisprGuideCapture, Custom, GeneExpression,
        MultiplexingCapture, VDJ, VDJ_B, VDJ_T, VDJ_T_GD,
    };
    match ft {
        [GeneExpression] => Ok("GEX"),
        [AntibodyCapture] => Ok("ABC"),
        [AntigenCapture] => Ok("AGC"),
        [CrisprGuideCapture] => Ok("CGC"),
        [Custom] => Ok("CUST"),
        [MultiplexingCapture] => Ok("CMO"),
        [VDJ] => Ok("VDJ"),
        [VDJ_T] => Ok("VDJT"),
        [VDJ_B] => Ok("VDJB"),
        [VDJ_T_GD] => Ok("VDJTGD"),
        _ => bail!("unimplemented feature_type: '{:?}'", ft),
    }
}

impl<'a> TryFrom<&Section<'a>> for LibrariesCsv {
    type Error = anyhow::Error;

    fn try_from(sec: &Section<'a>) -> Result<Self> {
        use libsconst::{
            FASTQS, FASTQ_ID, FASTQ_PATH, FEATURE_TYPES, GEM_WELL, LANES, LIBS_INT_OPT_HDRS,
            LIBS_INT_REQ_HDRS, LIBS_OPT_HDRS, LIBS_REQ_HDRS, PHYSICAL_LIBRARY_ID, SAMPLE_INDICES,
            SUBSAMPLE_RATE,
        };
        let hdr = sec.name;
        match CsvParser::new(sec.clone(), LIBS_INT_REQ_HDRS, LIBS_INT_OPT_HDRS) {
            Ok(parser) => {
                let mut data = vec![];
                let mut num_vdj_libs = 0;
                for row in parser.rows() {
                    let ctx = ParseCtx::HdrRow(hdr, row + 1);
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
                    let feature_types = parser.find_req(row, FEATURE_TYPES)?;
                    let feature_types: Vec<_> = parse_vec(feature_types.clone())
                        .map_err(|_: NomErr<'_>| {
                            anyhow!(
                                "{ctx} has invalid {FEATURE_TYPES} '{}' at line: {}, col: {}",
                                feature_types.fragment(),
                                feature_types.location_line(),
                                feature_types.get_utf8_column()
                            )
                        })?
                        .into_iter()
                        .map(FeatureType::try_from)
                        .try_collect()?;
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
                    let physical_library_id = parser
                        .find_opt(row, PHYSICAL_LIBRARY_ID)?
                        .and_then(empty_is_none)
                        .map(|pli| {
                            pli.parse::<Ident>(ctx.with_col(PHYSICAL_LIBRARY_ID))
                                .map(String::from)
                        })
                        .unwrap_or_else(|| {
                            if feature_types == [FeatureType::VDJ] {
                                num_vdj_libs += 1;
                                if num_vdj_libs > 1 {
                                    bail!(
                                        r#"{} is invalid: no more than 1 VDJ library may rely upon an auto-generated physical_library_id
Please either provide physical_libary_ids for all VDJ libraries or use more specific feature_types like VDJ-T and VDJ-B."#,
                                        ctx.with_col(PHYSICAL_LIBRARY_ID)
                                    );
                                }
                            }
                            feature_types_short(&feature_types)
                                .map(|x| format!("{x}_{}", gem_well.0))
                        })?;
                    let subsample_rate = parser
                        .find_opt(row, SUBSAMPLE_RATE)?
                        .and_then(empty_is_none)
                        .map(|sr| sr.parse::<Probability>(ctx.with_col(SUBSAMPLE_RATE)))
                        .transpose()?
                        .map(|sr| sr.0);
                    data.push(Library::BclProcessor {
                        fastq_path,
                        sample_indices,
                        lanes,
                        feature_types,
                        physical_library_id,
                        gem_well,
                        subsample_rate,
                    });
                }
                check_duplicate_libraries(&data)?;
                check_gem_wells(&data)?;
                check_physical_library_ids(&data)?;
                check_library_combinations(&data)?;
                Ok(LibrariesCsv(data))
            }
            Err(_) => {
                let parser = CsvParser::new(sec.clone(), LIBS_REQ_HDRS, LIBS_OPT_HDRS)?;
                let mut data = vec![];
                let mut num_vdj_libs = 0;
                let mut fastq_id_counts = TxHashMap::default(); // keep track of how many times we've seen a fastq-id.

                for row in parser.rows() {
                    let ctx = ParseCtx::HdrRow(hdr, row + 1);
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
                    let feature_types = parser.find_req(row, FEATURE_TYPES)?;
                    let feature_types: Vec<_> = parse_vec(feature_types.clone())
                        .map_err(|_: NomErr<'_>| {
                            anyhow!(
                                "{ctx} has invalid {FEATURE_TYPES} '{}' at line: {}, col: {}",
                                feature_types.fragment(),
                                feature_types.location_line(),
                                feature_types.get_utf8_column()
                            )
                        })?
                        .into_iter()
                        .map(FeatureType::try_from)
                        .try_collect()?;
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
                    let physical_library_id = parser
                        .find_opt(row, PHYSICAL_LIBRARY_ID)?
                        .and_then(empty_is_none)
                        .map(|pli| {
                            pli.parse::<Ident>(ctx.with_col(PHYSICAL_LIBRARY_ID))
                                .map(String::from)
                        })
                        .unwrap_or_else(|| {
                            if feature_types == [FeatureType::VDJ] {
                                num_vdj_libs += 1;
                                if num_vdj_libs > 1 {
                                    bail!(
                                        r#"{} is invalid: no more than 1 VDJ library may rely upon auto-generated physical_library_id
Please either provide physical_libary_ids for all VDJ libraries or use more specific feature_types like VDJ-T and VDJ-B."#,
                                        ctx.with_col(PHYSICAL_LIBRARY_ID)
                                    );
                                }
                            }
                            feature_types_short(&feature_types)
                                .map(|x| format!("{x}_{}", gem_well.0))
                        })?;
                    let subsample_rate = parser
                        .find_opt(row, SUBSAMPLE_RATE)?
                        .and_then(empty_is_none)
                        .map(|sr| sr.parse::<Probability>(ctx.with_col(SUBSAMPLE_RATE)))
                        .transpose()?
                        .map(|sr| sr.0);
                    data.push(Library::Bcl2Fastq {
                        fastq_id,
                        fastqs,
                        lanes,
                        feature_types,
                        physical_library_id,
                        gem_well,
                        subsample_rate,
                        ilmn_fastq_id,
                    });
                }
                check_duplicate_libraries(&data)?;
                check_gem_wells(&data)?;
                check_physical_library_ids(&data)?;
                check_library_combinations(&data)?;
                Ok(LibrariesCsv(data))
            }
        }
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
    /// NOTE: these may be +-concatenated groupings. Access via the
    /// sample_barcode_ids method to control how they are unpacked.
    probe_barcode_ids: Option<Vec<String>>,
    overhang_ids: Option<Vec<String>>,
    pub description: String,
    pub expect_cells: Option<usize>,
    pub force_cells: Option<usize>,
    pub emptydrops_minimum_umis: Option<usize>,
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
            self.probe_barcode_ids.as_ref(),
            self.overhang_ids.as_ref(),
        ) {
            (Some(x), None, None) | (None, None, Some(x)) => {
                // cmo_ids or overhang_ids
                Some(x.iter().map(|v| vec![v.as_str()]).collect())
            }
            (None, Some(x), None) => {
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
            (None, None, None) => None,
            _ => panic!("found more than one source of sample barcode IDs"),
        }
    }

    /// Return the type of cell multiplexing.
    pub fn cell_multiplexing_type(&self) -> CellMultiplexingType {
        match (
            self.cmo_ids.is_some(),
            self.probe_barcode_ids.is_some(),
            self.overhang_ids.is_some(),
        ) {
            (true, false, false) => CellMultiplexingType::CMO,
            (false, true, false) => CellMultiplexingType::RTL,
            (false, false, true) => CellMultiplexingType::OH,
            _ => unreachable!(),
        }
    }

    /// Return the appropriate column name, either cmo_ids or probe_barcode_ids.
    pub fn sample_barcode_ids_column_name(&self) -> &str {
        match self.cell_multiplexing_type() {
            CellMultiplexingType::CMO => samplesconst::CMO_IDS,
            CellMultiplexingType::RTL => samplesconst::PROBE_BARCODE_IDS,
            CellMultiplexingType::OH => samplesconst::OH_IDS,
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

    pub fn has_force_cells(&self) -> bool {
        self.0.iter().any(|sample| sample.force_cells.is_some())
    }

    pub fn has_expect_cells(&self) -> bool {
        self.0.iter().any(|sample| sample.expect_cells.is_some())
    }

    pub fn has_emptydrops_minimum_umis(&self) -> bool {
        self.0
            .iter()
            .any(|sample| sample.emptydrops_minimum_umis.is_some())
    }

    pub fn get_expect_cells(&self) -> TxHashMap<String, Option<usize>> {
        self.0
            .iter()
            .map(|sample| (sample.sample_id.clone(), sample.expect_cells))
            .collect()
    }

    pub fn get_force_cells(&self) -> TxHashMap<String, Option<usize>> {
        self.0
            .iter()
            .map(|sample| (sample.sample_id.clone(), sample.force_cells))
            .collect()
    }

    pub fn get_emptydrops_minimum_umis(&self) -> TxHashMap<String, Option<usize>> {
        self.0
            .iter()
            .map(|sample| (sample.sample_id.clone(), sample.emptydrops_minimum_umis))
            .collect()
    }

    pub fn is_rtl_multiplexed(&self) -> bool {
        self.0
            .iter()
            .all(|sample| sample.cell_multiplexing_type() == CellMultiplexingType::RTL)
    }

    pub fn get_cmo_sample_map(&self) -> TxHashMap<String, String> {
        let mut res = TxHashMap::default();
        for row in &self.0 {
            if let Some(ref cmo_ids) = row.cmo_ids {
                for cmo_id in cmo_ids {
                    if res.insert(cmo_id.clone(), row.sample_id.clone()).is_some() {
                        panic!("cmo_id {cmo_id} used for multiple samples");
                    }
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

mod samplesconst {
    pub const SAMPLE_ID: &str = "sample_id";
    pub const CMO_IDS: &str = "cmo_ids";
    pub const PROBE_BARCODE_IDS: &str = "probe_barcode_ids";
    pub const OH_IDS: &str = "overhang_ids";
    pub const _GEM_WELLS: &str = "gem_wells";
    pub const DESCRIPTION: &str = "description";
    pub const EXPECT_CELLS: &str = "expect_cells";
    pub const FORCE_CELLS: &str = "force_cells";
    pub const EMPTYDROPS_MINIMUM_UMIS: &str = "emptydrops_minimum_umis";
    pub const SAMP_REQ_HDRS: &[&str] = &[SAMPLE_ID];
    #[cfg(not(test))]
    pub const SAMP_OPT_HDRS: &[&str] = &[
        CMO_IDS,
        OH_IDS,
        PROBE_BARCODE_IDS,
        DESCRIPTION,
        EXPECT_CELLS,
        FORCE_CELLS,
        EMPTYDROPS_MINIMUM_UMIS,
    ];
    #[cfg(test)]
    pub const SAMP_OPT_HDRS: &[&str] = &[
        CMO_IDS,
        PROBE_BARCODE_IDS,
        OH_IDS,
        _GEM_WELLS,
        DESCRIPTION,
        EXPECT_CELLS,
        FORCE_CELLS,
        EMPTYDROPS_MINIMUM_UMIS,
    ];
}

impl<'a> TryFrom<(&TxHashSet<GemWell>, &Section<'a>)> for SamplesCsv {
    type Error = anyhow::Error;

    fn try_from((valid_gws, sec): (&TxHashSet<GemWell>, &Section<'a>)) -> Result<Self> {
        use samplesconst::{
            CMO_IDS, DESCRIPTION, EMPTYDROPS_MINIMUM_UMIS, EXPECT_CELLS, FORCE_CELLS, OH_IDS,
            PROBE_BARCODE_IDS, SAMPLE_ID, SAMP_OPT_HDRS, SAMP_REQ_HDRS, _GEM_WELLS,
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
            let probe_barcode_ids = parser.find_opt(row, PROBE_BARCODE_IDS)?;
            let overhang_ids = parser.find_opt(row, OH_IDS)?;
            let sample_barcode_ids = cmo_ids.or(probe_barcode_ids).or(overhang_ids);
            let cell_multiplexing_type = match (cmo_ids, probe_barcode_ids, overhang_ids) {
                (Some(_), None, None) => CellMultiplexingType::CMO,
                (None, Some(_), None) => CellMultiplexingType::RTL,
                (None, None, Some(_)) => CellMultiplexingType::OH,
                (None, None, None) => {
                    bail!(
                        // CELLRANGER-7549
                        "{} requires either {} or {} column to be specified",
                        ctx,
                        CMO_IDS,
                        PROBE_BARCODE_IDS
                    )
                }
                (Some(_), Some(_), _) => {
                    bail!(
                        "{} has mutually exclusive columns {} and {}",
                        ctx,
                        CMO_IDS,
                        PROBE_BARCODE_IDS
                    )
                }
                (None, Some(_), Some(_)) => {
                    bail!(
                        "{} has mutually exclusive columns {} and {}",
                        ctx,
                        PROBE_BARCODE_IDS,
                        OH_IDS
                    )
                }
                (Some(_), None, Some(_)) => {
                    bail!(
                        "{} has mutually exclusive columns {} and {}",
                        ctx,
                        CMO_IDS,
                        OH_IDS
                    )
                }
            };
            let expect_cells = parser.find_opt(row, EXPECT_CELLS)?;
            let force_cells = parser.find_opt(row, FORCE_CELLS)?;
            let emptydrops_minimum_umis = parser.find_opt(row, EMPTYDROPS_MINIMUM_UMIS)?;
            if probe_barcode_ids.is_none() && overhang_ids.is_none() {
                if expect_cells.is_some() {
                    bail!(
                        "{} has {} column specified without the required {} column",
                        ctx,
                        EXPECT_CELLS,
                        PROBE_BARCODE_IDS
                    );
                };
                if force_cells.is_some() {
                    bail!(
                        "{} has {} column specified without the required {} column",
                        ctx,
                        FORCE_CELLS,
                        PROBE_BARCODE_IDS
                    );
                };
            };
            let sample_barcode_column_name = match cell_multiplexing_type {
                CellMultiplexingType::CMO => CMO_IDS,
                CellMultiplexingType::RTL => PROBE_BARCODE_IDS,
                CellMultiplexingType::OH => OH_IDS,
            };

            let sample_barcode_ids = sample_barcode_ids
                .and_then(empty_is_none)
                .map(|sample_barcode_ids| {
                    parse_vec(sample_barcode_ids.clone())
                        .map_err(|_: NomErr<'_>| {
                            anyhow!(
                                "{ctx} has invalid {} '{}' at line: {}, col: {}",
                                sample_barcode_column_name,
                                sample_barcode_ids.fragment(),
                                sample_barcode_ids.location_line(),
                                sample_barcode_ids.get_utf8_column(),
                            )
                        })?
                        .into_iter()
                        .map(parse_prefixed_range)
                        .collect::<Result<Vec<_>>>()?
                        .into_iter()
                        .flatten()
                        .map(|multiplexing_id| {
                            match cell_multiplexing_type {
                                CellMultiplexingType::RTL => {
                                    validate_probe_barcode_id(multiplexing_id)
                                }
                                CellMultiplexingType::CMO => {
                                    multiplexing_id.parse::<Ident>().map(|_| multiplexing_id)
                                }
                                // CELLRANGER-7549
                                CellMultiplexingType::OH => Ok(multiplexing_id),
                            }
                            .with_context(|| {
                                format!(
                                    "{ctx} has invalid {} '{}' at line: {}, col: {}",
                                    sample_barcode_column_name,
                                    sample_barcode_ids.fragment(),
                                    sample_barcode_ids.location_line(),
                                    sample_barcode_ids.get_utf8_column()
                                )
                            })
                        })
                        .collect::<Result<TxHashSet<_>>>()
                })
                .transpose()?
                .map(|sample_barcode_ids| sample_barcode_ids.into_iter().collect::<Vec<_>>());

            let gem_wells_frag = None;
            let _gem_wells = gem_wells_frag
                .and_then(empty_is_none)
                .map(|gem_wells| {
                    parse_vec(gem_wells.clone())
                        .map_err(|_: NomErr<'_>| {
                            anyhow!(
                                "{ctx} has invalid {_GEM_WELLS} '{}' at line: {}, col: {}",
                                gem_wells.fragment(),
                                gem_wells.location_line(),
                                gem_wells.get_utf8_column(),
                            )
                        })?
                        .into_iter()
                        .map(|gws| {
                            Ok(parse_range::<u16>(gws.clone())?
                                .map(GemWell)
                                .collect::<Vec<_>>())
                        })
                        .collect::<Result<Vec<_>>>()
                        .map(|gws| gws.into_iter().flatten().collect::<Vec<_>>())
                })
                .unwrap_or_else(|| Ok(vec![GemWell(1)]))?;
            for gem_well in &_gem_wells {
                if !valid_gws.contains(gem_well) {
                    if let Some(gem_wells_frag) = gem_wells_frag {
                        bail!(
                            "{} has invalid gem_well '{}' at line: {}, col: {}: please check for consistency with [{}]",
                            ctx,
                            gem_well.0,
                            gem_wells_frag.location_line(),
                            gem_wells_frag.get_utf8_column(),
                            multiconst::LIBRARIES,
                        );
                    } else {
                        bail!(
                            "{} has invalid gem_well '{}': please check for consistency with [{}]",
                            ctx,
                            gem_well.0,
                            multiconst::LIBRARIES
                        );
                    }
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
            if force_cells.is_some() && expect_cells.is_some() {
                bail!(
                    "{} has both force-cells and expect-cells specified, only one is allowed.",
                    ctx,
                );
            };
            let emptydrops_minimum_umis = emptydrops_minimum_umis
                .and_then(empty_is_none)
                .map(|ec| ec.parse::<usize>(ctx.with_col(EMPTYDROPS_MINIMUM_UMIS)))
                .transpose()?;

            // Check to make sure that sample barcode IDs were provided and that all entries are non-empty
            let empty_sample_barcodes_ids = match sample_barcode_ids {
                Some(ref ids) => ids.iter().any(String::is_empty),
                None => true,
            };

            if empty_sample_barcodes_ids {
                bail!(
                    "{} has an empty entry in the {} column. All entries must be non-empty.",
                    ctx,
                    sample_barcode_column_name,
                );
            }

            data.push(match cell_multiplexing_type {
                CellMultiplexingType::CMO => SampleRow {
                    sample_id,
                    cmo_ids: sample_barcode_ids,
                    probe_barcode_ids: None,
                    overhang_ids: None,
                    description,
                    force_cells: None,
                    expect_cells: None,
                    emptydrops_minimum_umis,
                },
                CellMultiplexingType::RTL => SampleRow {
                    sample_id,
                    cmo_ids: None,
                    probe_barcode_ids: sample_barcode_ids,
                    overhang_ids: None,
                    description,
                    force_cells,
                    expect_cells,
                    emptydrops_minimum_umis,
                },
                CellMultiplexingType::OH => SampleRow {
                    sample_id,
                    cmo_ids: None,
                    probe_barcode_ids: None,
                    overhang_ids: sample_barcode_ids,
                    description,
                    force_cells,
                    expect_cells,
                    emptydrops_minimum_umis,
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
    let overhang_wl_spec = WhitelistSpec::DynamicTranslation {
        translation_whitelist_path: find_whitelist(DEFAULT_OVERHANG_WL, true, None).unwrap(),
    };
    WhitelistSource::from_spec(&overhang_wl_spec, true, None)
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
        use std::collections::hash_map::Entry::{Occupied, Vacant};
        let hdr = sec.name;
        let parser = CsvParser::new(sec.clone(), GEM_WELL_REQ_HDRS, GEM_WELL_OPT_HDRS)?;
        let mut data = TxHashMap::default();
        for row in parser.rows() {
            let ctx = ParseCtx::HdrRow(hdr, row + 1);
            let gem_well = parser
                .find_req(row, GEM_WELL)?
                .parse::<GemWell>(ctx.with_col(GEM_WELL))?;
            if !valid_gws.contains(&gem_well) {
                bail!(
                    "{} has invalid gem_well: {}, please check for consistency with [{}]",
                    ctx,
                    gem_well.0,
                    multiconst::LIBRARIES,
                );
            }
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
            match data.entry(gem_well) {
                Vacant(entry) => {
                    entry.insert(GemWellParams {
                        gem_well,
                        force_cells,
                        expect_cells,
                        vdj_force_cells,
                    });
                }
                Occupied(_) => {
                    bail!("{} duplicate entry detected: {}", ctx, gem_well.0);
                }
            }
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
                    assert_ne!(existing_value, functional_name)
                }
            }
        }
    }

    let functional_map = if !functional_map.is_empty() {
        Some(functional_map)
    } else {
        None
    };

    if specificity_controls.is_some() || functional_map_csv.is_some() {
        Some(FeatureConfig {
            specificity_controls,
            beam_mode,
            functional_map,
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
#[derive(Debug, Default, Serialize, Deserialize)]
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
    fn from_csv<P: AsRef<Path>>(path: P) -> Result<Self> {
        let f = File::open(path.as_ref())?;
        let reader = BufReader::new(f);
        let xtra = XtraData::new(path);
        let config = MultiConfigCsv::from_reader(reader, xtra)?;
        Ok(config)
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
            _ => anyhow!("failed to parse CSV, incomplete information available to pinpoint error"),
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
            let fastq = lib.to_fastq_def(self)?;
            let library_features = legacy_library_type(lib.feature_types())?.into();
            let gem_well = lib.gem_well().0.into();
            builder.push_library(physical_library_id, fastq, library_features, gem_well)?;
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
                                sample.cell_multiplexing_type(),
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
            .map_or(false, SamplesCsv::has_overhang_ids)
    }

    /// Return true if multiplexed using RTL probe barcodes.
    pub fn is_rtl_multiplexed(&self) -> bool {
        self.samples
            .as_ref()
            .map_or(false, SamplesCsv::is_rtl_multiplexed)
    }

    /// Return whether either [gene-expression] or [samples] is RTL.
    pub fn is_rtl(&self) -> bool {
        self.gene_expression
            .as_ref()
            .is_some_and(GeneExpressionParams::is_rtl)
            || self.is_rtl_multiplexed()
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
            if self.$name.is_some() {
                bail!(
                    "failed to parse CSV, duplicate [{}] at line: {}, col: {}",
                    section.name.fragment(),
                    section.name.location_line(),
                    section.name.get_utf8_column(),
                );
            }
            let $name = $typ::try_from(section)?;
            self.$name = Some($name);
            Ok(())
        }
    };
}

macro_rules! setter_validate_gws {
    ($name:ident, $typ:ident, $constant:expr) => {
        fn $name(&mut self, section: &Section<'_>) -> Result<()> {
            if self.$name.is_some() {
                bail!(
                    "failed to parse CSV, duplicate [{}] at line: {}, col: {}",
                    section.name.fragment(),
                    section.name.location_line(),
                    section.name.get_utf8_column(),
                );
            }
            if let Some(libraries) = self.libraries.as_ref() {
                let valid_gws: TxHashSet<_> = libraries.0.iter().map(Library::gem_well).collect();
                let $name = $typ::try_from((&valid_gws, section))?;
                self.$name = Some($name);
                return Ok(());
            }
            bail!(
                "failed to parse [{}] before [{}]",
                $constant,
                multiconst::LIBRARIES,
            )
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
            feature: feature_owned,
            vdj,
            libraries,
            samples: samples_owned,
            gem_wells,
            antigen_specificity,
            functional_map,
        } = self;
        let gene_expression = gene_expression_owned.as_ref();
        let chemistry = gene_expression.map(|x| x.chemistry);
        let feature = feature_owned.as_ref();
        let samples = samples_owned.as_ref();
        let has_probe_barcode_ids = samples.map_or(false, SamplesCsv::has_probe_barcode_ids);
        let is_rtl = has_probe_barcode_ids || chemistry.and_then(|x| x.is_rtl()) == Some(true);
        let libraries = if let Some(libraries) = libraries {
            libraries
        } else {
            bail!(
                "failed to parse CSV, [{}] section not provided",
                multiconst::LIBRARIES,
            );
        };
        if libraries.has_gene_expression() && gene_expression.is_none() {
            bail!(
                "failed to parse CSV: [{}] section omitted but gene expression libraries provided",
                multiconst::GENE_EXPRESSION
            );
        }
        let has_feature_ref = feature.map(|x| x.reference_path.is_some()).unwrap_or(false);
        if libraries.has_feature_barcode() && !has_feature_ref {
            bail!(
                "failed to parse CSV: [{}] reference is missing but feature barcode libraries provided",
                multiconst::FEATURE
            );
        }
        if libraries.has_vdj() && vdj.is_none() {
            bail!(
                "failed to parse CSV: [{}] section omitted but VDJ libraries provided",
                multiconst::VDJ
            );
        }
        if libraries.has_multiplexing() && samples.is_none() {
            bail!(
                "failed to parse CSV: [{}] section omitted but Multiplexing Capture libraries provided",
                multiconst::SAMPLES
            );
        }
        if antigen_specificity.is_some() && !libraries.has_antigen_capture() {
            bail!(
                "failed to parse CSV: [{}] section is provided but no Antigen Capture libraries",
                multiconst::ANTIGEN_SPECIFICITY,
            );
        }
        if functional_map.is_some() && !libraries.has_feature_barcode() {
            bail!(
                "failed to parse CSV: [{}] section is provided but no feature barcode libraries provided",
                multiconst::FUNCTIONAL_MAP,
            );
        }
        if libraries.has_antigen_capture()
            && gene_expression
                .unwrap()
                .invalid_parameter_with_antigen_capture()
                .is_some()
        {
            bail!(
                "failed to parse CSV: [{}] section specified {} parameter but [{}] section \
                 includes Antigen Capture libraries. Invalid parameter in this context",
                multiconst::GENE_EXPRESSION,
                gene_expression
                    .unwrap()
                    .invalid_parameter_with_antigen_capture()
                    .unwrap(),
                multiconst::LIBRARIES
            );
        }
        if chemistry == Some(ChemistryParam::SFRP) && samples.is_some() {
            bail!(
                "failed to parse CSV: [{}] section specified `SFRP` as the `chemistry` but a \
                 [{}] section is also provided",
                multiconst::GENE_EXPRESSION,
                multiconst::SAMPLES
            );
        }
        if chemistry == Some(ChemistryParam::MFRP) && (samples.is_none() || !has_probe_barcode_ids)
        {
            bail!(
                "failed to parse CSV: [{}] section specified `MFRP` as the `chemistry` but a \
                 [{}] section with {} column is not provided",
                multiconst::GENE_EXPRESSION,
                multiconst::SAMPLES,
                samplesconst::PROBE_BARCODE_IDS,
            );
        }

        if chemistry == Some(ChemistryParam::MFRP)
            && gene_expression.map_or(false, |gex| gex.has_expect_cells() || gex.has_force_cells())
        {
            bail!(
                "failed to parse CSV: [{}] section specified `MFRP` as the `chemistry` and cell calling parameters. \
                 For multiplex Fixed RNA Profiling libraries `expect-cells` or `force-cells' parameter is valid only in [{}] section.",
                multiconst::GENE_EXPRESSION,
                multiconst::SAMPLES
            );
        }

        let probe_set = gene_expression.and_then(|gex| gex.probe_set.as_ref());
        if is_rtl && libraries.has_gene_expression() && probe_set.is_none() {
            bail!(
                "failed to parse CSV: [{}] section is missing `probe-set` a required parameter for \
                 Fixed RNA Profiling",
                multiconst::GENE_EXPRESSION,
            );
        }

        // Disallow setting include-introns for RTL chemistries.
        if is_rtl
            && gene_expression.map_or(false, |x| {
                x.include_introns != multiconst::DEFAULT_INCLUDE_INTRONS
            })
        {
            bail!(ERROR_INCLUDE_INTRONS_WITH_RTL);
        }

        if has_probe_barcode_ids {
            // NOTE: non-multiplexed chemistry + having probe_barcode_ids column(invalid) passes here but would fail in checks above
            if chemistry.and_then(|x| x.is_rtl()) == Some(false) {
                bail!(
                    "failed to parse CSV: [{}] section manually specifies a non-Fixed RNA Profiling chemistry but [{}] section has a {} column. The {} column may only be specified with Fixed RNA Profiling chemistries",
                    multiconst::GENE_EXPRESSION,
                    multiconst::SAMPLES,
                    samplesconst::PROBE_BARCODE_IDS,
                    samplesconst::PROBE_BARCODE_IDS,
                );
            }
        }

        if gene_expression.map_or(false, |gex| gex.has_expect_cells() || gex.has_force_cells())
            && samples.map_or(false, SamplesCsv::has_probe_barcode_ids)
        {
            bail!(
                "failed to parse CSV: [{}] section has a {} column indicating multiplex Fixed RNA Profiling chemistry, \
                however a cell calling parameter is specified in [{}] section. The parameters `expect-cells` or `force-cells' \
                are valid in [{}] section for singleplex Fixed RNA Profiling and in [{}] section for multiplex Fixed RNA Profiling.",
                multiconst::SAMPLES,
                samplesconst::PROBE_BARCODE_IDS,
                multiconst::GENE_EXPRESSION,
                multiconst::GENE_EXPRESSION,
                multiconst::SAMPLES,
            );
        }

        if gene_expression.map_or(false, |gex| gex.has_expect_cells() || gex.has_force_cells())
            && samples.map_or(false, |s| s.has_expect_cells() || s.has_force_cells())
        {
            bail!(
                "failed to parse CSV: 'expect-cells' or 'force-cells' parameters are specified in both [{}] and [{}] sections",
                multiconst::GENE_EXPRESSION,
                multiconst::SAMPLES
            );
        }

        Ok(MultiConfigCsv {
            gene_expression: gene_expression_owned,
            feature: feature_owned,
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
    use anyhow::Result;
    use barcode::whitelist::BarcodeId;
    use itertools::Itertools;
    use metric::TxHashMap;
    use std::path::Path;

    // initialize insta test harness
    #[ctor::ctor]
    fn init() {
        // this ensures insta knows where to find its snap tests
        let cwd = std::env::current_dir().unwrap();
        let workspace_root = cwd.parent().unwrap();
        std::env::set_var("INSTA_WORKSPACE_ROOT", workspace_root);
    }

    #[test]
    fn load_simple_external_cmos() -> Result<()> {
        let csv = r#"
[gene-expression]
ref,/path/to/gex/ref
chemistry,

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
        insta::assert_display_snapshot!(res.unwrap_err());
    }

    #[test]
    fn invalid_sample_id_chars_2() {
        let csv = r#"
[gene-expression]
ref,GRCh38-2020-A-chr21

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
        insta::assert_display_snapshot!(res.unwrap_err());
    }

    #[test]
    fn too_long_sample_id() {
        let csv = r#"
[gene-expression]
ref,GRCh38-2020-A-chr21

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
        insta::assert_display_snapshot!(res.unwrap_err());
    }

    #[test]
    fn load_trailing_whitespace() -> Result<()> {
        let csv = r#"
[gene-expression]
ref,GRCh38-2020-A-chr21
chemistry,SC5P-R2

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
no-bam,FALSE,,
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
        insta::assert_display_snapshot!(res.unwrap_err());
    }

    /// Ensure that antigen capture is accompanied with a VDJ library
    #[test]
    fn test_antigen_requires_vdj() -> Result<()> {
        let csv = r#"
    [gene-expression]
    ref,mm10-2020-A-chr19

    [feature]
    ref,cellranger/multi/feature_refs/20211122_v1.1.csv

    [libraries]
    fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
    tiny_gex,fastqs/cellranger/multi/1245140_gex_vdj_beam_ab/gex,any,gex,gene expression,0.5
    tiny_an_b,fastqs/cellranger/multi/1245140_gex_vdj_beam_ab/an_b,any,ab,antigen capture,0.5
    "#;

        let xtra = XtraData::new("test::antigen_required_vdj");
        let res = MultiConfigCsv::from_reader(csv.as_bytes(), xtra);
        insta::assert_display_snapshot!(res.unwrap_err());
        Ok(())
    }

    #[test]
    fn test_antigen_specificity_requires_ag() -> Result<()> {
        let csv = r#"
    [gene-expression]
    ref,mm10-2020-A-chr19

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
        insta::assert_display_snapshot!(res.unwrap_err());
        Ok(())
    }

    #[test]
    fn test_antigen_specificity_multiple_control_ids() -> Result<()> {
        let csv = r#"
    [gene-expression]
    ref,mm10-2020-A-chr19

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
        insta::assert_display_snapshot!(res.unwrap_err());
        Ok(())
    }

    #[test]
    fn test_antigen_specificity_missing_allele() -> Result<()> {
        let csv = r#"
    [gene-expression]
    ref,mm10-2020-A-chr19

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
        insta::assert_display_snapshot!(res.unwrap_err());
        Ok(())
    }

    #[test]
    fn test_invalid_mhc_allele_character() -> Result<()> {
        let csv = r#"
    [gene-expression]
    ref,mm10-2020-A-chr19

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
        insta::assert_display_snapshot!(res.unwrap_err());
        Ok(())
    }

    #[test]
    fn test_functional_map_requires_feature() -> Result<()> {
        let csv = r#"
    [gene-expression]
    ref,mm10-2020-A-chr19

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
        insta::assert_display_snapshot!(res.unwrap_err());
        Ok(())
    }

    #[test]
    fn test_invalid_functional_map_value() -> Result<()> {
        // Duplicate feature ids
        let csv = r#"
    [gene-expression]
    ref,mm10-2020-A-chr19

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
        insta::assert_display_snapshot!(res.unwrap_err());

        // Duplicate functional name
        let csv = r#"
    [gene-expression]
    ref,mm10-2020-A-chr19

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
        insta::assert_display_snapshot!(res.unwrap_err());

        Ok(())
    }

    #[test]
    fn test_blank_lines() -> Result<()> {
        let csv = r#"
    [gene-expression]
    ref,mm10-2020-A-chr19

    probe-set,/path/to/probe_set

    [libraries]
    fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
    mygex,/path/to/fastqs,any,gex,gene expression,0.5,100

    mycmo,/path/to/fastqs,any,cmo,Multiplexing Capture,,

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
            gex.probe_set.as_deref(),
            Some(Path::new("/path/to/probe_set"))
        );
        assert_eq!(&gex.reference_path, Path::new("mm10-2020-A-chr19"));

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

    [libraries]
    fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
    mygex,/path/to/fastqs,any,gex,gene expression,0.5,100
    mycmo,/path/to/fastqs,any,cmo,Multiplexing Capture,,

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

    [libraries]
    fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
    mygex,/path/to/fastqs,any,gex,gene expression,0.5,100

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
        insta::assert_json_snapshot!(
            "expected_probe_barcode_pairing_graph",
            res.to_multi_graph("test", "test description", None)?
        );

        // Test that we ignore a detected probe barcode pairing if we have an
        // explicit one.
        let different_pairing = create_bc_pairing(&[("BC1", "BC10"), ("BC2", "BC7")]);
        insta::assert_json_snapshot!(
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

    [libraries]
    fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
    mygex,/path/to/fastqs,any,gex,gene expression,0.5,100

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

    [libraries]
    fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
    mygex,/path/to/fastqs,any,gex,gene expression,0.5,100

    [samples]
    sample_id,probe_barcode_ids
    sample1,BC1|BC2
    sample2,BC3
    "#;

        let xtra = XtraData::new("test::");
        let res = MultiConfigCsv::from_reader(csv.as_bytes(), xtra)?;
        let detected_pairing = create_bc_pairing(&[("BC1", "AB1"), ("BC3", "AB3")]);
        insta::assert_json_snapshot!(res.to_multi_graph(
            "test",
            "test description",
            Some(&detected_pairing)
        )?);
        Ok(())
    }

    #[test]
    fn test_duplicate_sample_ids() {
        let res = MultiConfigCsv::from_csv("test/invalid_csvs/duplicate_sample_ids.csv");
        insta::assert_display_snapshot!(res.unwrap_err());
    }

    #[test]
    fn test_mismatched_lib_feature_types() {
        let res = MultiConfigCsv::from_csv("test/invalid_csvs/mismatched_lib_feature_types.csv");
        insta::assert_display_snapshot!(res.unwrap_err());
    }

    #[test]
    fn test_redefined_physical_library_id() {
        let res = MultiConfigCsv::from_csv("test/invalid_csvs/redefined_physical_library_id.csv");
        insta::assert_display_snapshot!(res.unwrap_err());
    }

    #[test]
    fn test_force_cells_too_low() {
        insta::assert_display_snapshot!(MultiConfigCsv::from_csv(
            "test/invalid_csvs/force_cells_too_low.csv"
        )
        .unwrap_err());
    }

    #[test]
    fn test_duplicate_libraries() {
        insta::assert_display_snapshot!(MultiConfigCsv::from_csv(
            "test/invalid_csvs/duplicate_libraries.csv"
        )
        .unwrap_err());
    }

    #[test]
    fn test_multiplexing_no_samples() {
        insta::assert_display_snapshot!(MultiConfigCsv::from_csv(
            "test/invalid_csvs/multiplexing_no_samples.csv"
        )
        .unwrap_err());
    }

    #[test]
    fn test_multiple_vdjt_libraries() {
        insta::assert_debug_snapshot!(MultiConfigCsv::from_csv("test/multiple_vdj_t.csv"));
    }

    #[test]
    fn test_expect_cells_and_force_cells_samples() {
        insta::assert_display_snapshot!(MultiConfigCsv::from_csv(
            "test/invalid_csvs/expect_cells_and_force_cells_samples.csv"
        )
        .unwrap_err());
    }

    #[test]
    fn test_frp_chemistry_with_samples() {
        insta::assert_display_snapshot!(MultiConfigCsv::from_csv(
            "test/invalid_csvs/frp_chem_with_samples.csv"
        )
        .unwrap_err());
    }

    #[test]
    fn test_mfrp_chemistry_without_samples() {
        insta::assert_display_snapshot!(MultiConfigCsv::from_csv(
            "test/invalid_csvs/mfrp_chem_no_samples.csv"
        )
        .unwrap_err());
    }

    #[test]
    fn test_mfrp_chemistry_with_missing_probe_barcode_entry() {
        insta::assert_display_snapshot!(MultiConfigCsv::from_csv(
            "test/invalid_csvs/mfrp_chemistry_with_missing_probe_barcode_entry.csv"
        )
        .unwrap_err());
    }
}
