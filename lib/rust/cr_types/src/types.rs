#![expect(missing_docs)]
use anyhow::{Context, Result, bail, ensure};
use barcode::{Barcode, BarcodeConstruct, BarcodeFromString, BcSegSeq};
use fastq_set::filenames::FastqDef;
use fastq_set::filenames::bcl2fastq::SampleNameSpec;
use itertools::Itertools;
use martian::{AsMartianPrimaryType, MartianFileType, MartianPrimaryType};
use martian_derive::{MartianStruct, MartianType, martian_filetype};
use martian_filetypes::FileTypeRead;
use martian_filetypes::bin_file::BinaryFormat;
use martian_filetypes::json_file::{JsonFile, JsonFormat};
use metric::{AsMetricPrefix, Histogram, SimpleHistogram, TxHashMap, TxHashSet};
use rust_htslib::bam::{Read as _, Reader as BamReader};
use serde::de::{IgnoredAny, Visitor};
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use serde_with::DeserializeFromStr;
use shardio::SortKey;
use std::borrow::Cow;
use std::convert::TryFrom;
use std::fmt::{self, Display, Formatter};
use std::fs::read_to_string;
use std::hash::Hash;
use std::path::PathBuf;
use std::str::FromStr;
use strum::{EnumString, IntoStaticStr};
use umi::UmiType;
// Used to indicate a UMI count does not come from a probe
// in datasets that are mixed probe/no-probe (e.g. RTL + FB)
pub const PROBE_IDX_SENTINEL_VALUE: i32 = -1;

// File Types

martian_filetype!(BaiIndexFile, "bam.bai");

martian_filetype!(CsiIndexFile, "bam.csi");

martian_filetype!(H5File, "h5");

// Shardio file containing FeatureBarcodeCount records
martian_filetype!(CountShardFile, "csf");

martian_filetype! { BcSegmentCountFile, "bsc" }
pub type BcSegmentCountFormat = BinaryFormat<
    BcSegmentCountFile,
    TxHashMap<LibraryType, BarcodeConstruct<SimpleHistogram<BcSegSeq>>>,
>;

martian_filetype!(FeatureCountFile, "fbc");
pub type FeatureCountFormat = BinaryFormat<FeatureCountFile, Vec<i64>>;

martian_filetype!(BarcodeSetFile, "blf");
pub type BarcodeSetFormat = JsonFormat<BarcodeSetFile, TxHashSet<Barcode>>;

martian_filetype!(_FingerprintFile, "fprint");
pub type FingerprintFile = JsonFormat<_FingerprintFile, Vec<Fingerprint>>;

// End File Types

/// A feature ID.
pub type FeatureID = String;

/// A feature name.
pub type FeatureName = String;

/// A genome name.
#[derive(
    Clone,
    Debug,
    Default,
    Hash,
    Eq,
    PartialEq,
    Ord,
    PartialOrd,
    Serialize,
    Deserialize,
    derive_more::Deref,
    derive_more::Display,
)]
pub struct GenomeName(String);

impl GenomeName {
    /// Return this genome name as a string slice.
    pub fn as_str(&self) -> &str {
        &self.0
    }
}

impl AsRef<[u8]> for GenomeName {
    fn as_ref(&self) -> &[u8] {
        self.0.as_ref()
    }
}

impl From<&str> for GenomeName {
    fn from(value: &str) -> Self {
        GenomeName(value.to_string())
    }
}

impl From<GenomeName> for String {
    fn from(value: GenomeName) -> Self {
        value.0
    }
}

impl AsMetricPrefix for GenomeName {
    fn as_metric_prefix(&self) -> Option<&str> {
        Some(self.0.as_str())
    }
}

impl AsMartianPrimaryType for GenomeName {
    fn as_martian_primary_type() -> MartianPrimaryType {
        MartianPrimaryType::Str
    }
}

/// Count of the number of UMIs observed for one feature in one barcode.  Corresponds
/// to a single entry in the feature x barcode matrix.
#[derive(Serialize, Deserialize, Debug, Clone, Copy, Ord, PartialOrd, Eq, PartialEq)]
pub struct FeatureBarcodeCount {
    pub barcode: Barcode,
    pub feature_idx: u32,
    pub umi_count: u32,
}

/// Sort by barcode and then by feature.
pub struct BarcodeThenFeatureOrder;

impl SortKey<FeatureBarcodeCount> for BarcodeThenFeatureOrder {
    type Key = (Barcode, u32);

    fn sort_key(fbc: &FeatureBarcodeCount) -> Cow<'_, Self::Key> {
        Cow::Owned((fbc.barcode, fbc.feature_idx))
    }
}

/// Count of the number of UMIs observed for one probe in one barcode.  Corresponds
/// to a single entry in the probe x barcode matrix.
#[derive(Serialize, Deserialize, Clone, Ord, PartialOrd, Eq, PartialEq)]
pub struct ProbeBarcodeCount {
    pub barcode: Barcode,
    pub probe_idx: u32,
    pub umi_count: u32,
}

/// A single UMI count - corresponds to a single entry
/// in the `molecule_info.h5` file. Must be be nested
/// inside a `BcUmiInfo` to know the barcode associated
/// with this count.
#[derive(Clone, Copy, PartialOrd, Ord, PartialEq, Eq, Serialize, Deserialize, Debug)]
pub struct UmiCount {
    pub library_idx: u16,
    pub feature_idx: u32,
    pub umi: u32,
    pub read_count: u32,
    pub utype: UmiType,
    pub probe_idx: Option<i32>,
}

/// UMI count information for one barcode
#[derive(Serialize, Deserialize, Debug)]
pub struct BcUmiInfo {
    pub barcode: Barcode,
    pub umi_counts: Vec<UmiCount>,
}

impl SortKey<BcUmiInfo> for BcUmiInfo {
    type Key = Barcode;

    fn sort_key(bc_info: &BcUmiInfo) -> Cow<'_, Self::Key> {
        Cow::Owned(bc_info.barcode)
    }
}

impl BcUmiInfo {
    /// Convert a `BcUmiInfo`, which contains one entry per UMI, into a FeatureBarcodeCount
    /// which contains one entry per non-zero feature.
    pub fn feature_counts(&self) -> impl Iterator<Item = FeatureBarcodeCount> + '_ {
        SimpleHistogram::from_iter_owned(self.umi_counts.iter().map(|c| c.feature_idx))
            .into_iter()
            .map(|(feature_idx, umi_count)| FeatureBarcodeCount {
                barcode: self.barcode,
                umi_count: u32::try_from(umi_count.count()).unwrap(),
                feature_idx,
            })
    }

    /// Convert a `BcUmiInfo`, which contains one entry per UMI, into a `ProbeBarcodeCount`
    /// which contains one entry per non-zero probe.
    pub fn probe_counts(&self) -> impl Iterator<Item = ProbeBarcodeCount> + '_ {
        SimpleHistogram::from_iter_owned(self.umi_counts.iter().filter_map(|x| match x.probe_idx {
            None | Some(PROBE_IDX_SENTINEL_VALUE) => None,
            Some(probe_index) => Some(u32::try_from(probe_index).unwrap()),
        }))
        .into_iter()
        .map(|(probe_idx, umi_count)| ProbeBarcodeCount {
            barcode: self.barcode,
            umi_count: u32::try_from(umi_count.count()).unwrap(),
            probe_idx,
        })
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct FileOrBytes {
    #[mro_type = "file"]
    pub file: Option<PathBuf>,
    pub bytes: Option<String>,
}

macro_rules! enum_maker {
    ($(#[$meta:meta])* $name:ident, $( ($field:ident, $lit: literal $(, $alias: literal)*) ),*) => {
        #[derive(
            Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize, Deserialize,
            EnumString, IntoStaticStr,
        )]
        $(#[$meta])*
        pub enum $name {
            $(
                #[serde(rename = $lit)]
                #[strum(to_string = $lit)]
                $(
                    #[serde(alias = $alias)]
                    #[strum(serialize = $alias)]
                )*
                $field,
            )*
        }

        impl $name {
            /// Return the string representation of this variant.
            pub fn as_str(&self) -> &'static str {
                self.into()
            }
        }
    };
}

/// A sample barcode ID string
pub type SampleBarcodeID = String;

/// A sample ID string
pub type SampleId = String;

/// A SampleAssignment refers to the sample assignment status of a barcode or group of barcodes:
/// - `Unassigned`: not assigned to any sample
/// - `Assigned(String)`: assigned to the specified sample_id (a string)
#[derive(Debug, Clone, PartialEq, Eq, Hash, Ord, PartialOrd)]
pub enum SampleAssignment {
    NonMultiplexed,
    Unassigned,
    Assigned(SampleId),
}

impl SampleAssignment {
    pub fn as_str(&self) -> &str {
        match self {
            SampleAssignment::NonMultiplexed => "non_multiplexed",
            SampleAssignment::Unassigned => "unassigned",
            SampleAssignment::Assigned(sample_id) => sample_id.as_str(),
        }
    }
}

impl From<String> for SampleAssignment {
    fn from(s: String) -> Self {
        match s.as_str() {
            "unassigned" | "Unassigned" => SampleAssignment::Unassigned,
            "non_multiplexed" | "Non_multiplexed" => SampleAssignment::NonMultiplexed,
            _ => SampleAssignment::Assigned(s),
        }
    }
}

impl AsMartianPrimaryType for SampleAssignment {
    fn as_martian_primary_type() -> MartianPrimaryType {
        MartianPrimaryType::Str
    }
}

impl Display for SampleAssignment {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

impl Serialize for SampleAssignment {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        match self {
            SampleAssignment::Assigned(sample) => serializer.serialize_str(sample),
            SampleAssignment::Unassigned => serializer.serialize_str("unassigned"),
            SampleAssignment::NonMultiplexed => serializer.serialize_str("non_multiplexed"),
        }
    }
}

impl<'de> Deserialize<'de> for SampleAssignment {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        struct SampleAssignmentVisitor;
        impl Visitor<'_> for SampleAssignmentVisitor {
            type Value = SampleAssignment;
            fn expecting(&self, f: &mut Formatter<'_>) -> fmt::Result {
                f.write_str(
                    "A string of either the assigned sample ID, \"unassigned\", or \"non_multiplexed\"",
                )
            }

            fn visit_str<E: serde::de::Error>(self, s: &str) -> Result<Self::Value, E> {
                Ok(SampleAssignment::from(s.to_string()))
            }
        }
        deserializer.deserialize_any(SampleAssignmentVisitor)
    }
}

pub type SampleBarcodesFile = JsonFile<SampleBarcodes>;

/// A map of sample IDs to their barcodes.
#[derive(Debug, Eq, PartialEq)]
pub struct SampleBarcodes(TxHashMap<SampleAssignment, Vec<Barcode>>);

impl SampleBarcodes {
    /// Read a sample barcodes JSON file and return the samples.
    pub fn read_samples(
        sample_barcodes_json: Option<&SampleBarcodesFile>,
    ) -> Result<Vec<SampleAssignment>> {
        let Some(sample_barcodes_json) = sample_barcodes_json else {
            return Ok(vec![SampleAssignment::NonMultiplexed]);
        };

        let f: JsonFile<TxHashMap<SampleAssignment, IgnoredAny>> =
            JsonFile::from_path(sample_barcodes_json);
        Ok(f.read()?.into_keys().collect())
    }

    /// Read a sample barcodes JSON file and return the samples and NonMultiplexed.
    pub fn read_samples_with_unassigned(
        sample_barcodes_json: Option<&SampleBarcodesFile>,
    ) -> Result<Vec<SampleAssignment>> {
        if sample_barcodes_json.is_none() {
            return Ok(vec![SampleAssignment::NonMultiplexed]);
        };

        let mut samples = SampleBarcodes::read_samples(sample_barcodes_json)?;
        samples.push(SampleAssignment::Unassigned);
        Ok(samples)
    }

    /// Read a sample barcodes JSON file.
    pub fn read_from_json(sample_barcodes_json: Option<&SampleBarcodesFile>) -> Result<Self> {
        let Some(sample_barcodes_json) = sample_barcodes_json else {
            return Ok(Self(TxHashMap::default()));
        };

        let loaded: TxHashMap<SampleAssignment, Vec<BarcodeFromString>> = serde_json::from_str(
            &read_to_string(sample_barcodes_json)
                .with_context(|| sample_barcodes_json.display().to_string())?,
        )?;

        let sample_to_barcodes = loaded
            .into_iter()
            .map(|(sample, barcodes)| (sample, barcodes.into_iter().map(Barcode::from).collect()))
            .collect();
        Ok(SampleBarcodes(sample_to_barcodes))
    }

    /// Return whether this analysis is multiplexed.
    /// "multiplexed" in this context describes any multi run where we are writing a per-sample
    /// molecule info/bams/etc rather than the count or raw (uber) versions of files even if it's
    /// "non-multiplexed" multi run, that single sample is treated like a single muxed sample for
    /// simplicity internally.
    pub fn is_multiplexed(&self) -> bool {
        !self.0.is_empty()
    }

    /// Return the barcodes for one sample.
    pub fn get_barcodes(&self, sample: &SampleAssignment) -> Option<&[Barcode]> {
        if self.is_multiplexed() {
            Some(self.0[sample].as_slice())
        } else {
            None
        }
    }

    /// Return the barcodes for one sample.
    pub fn into_barcodes(mut self, sample: &SampleAssignment) -> Option<Vec<Barcode>> {
        if self.is_multiplexed() {
            let barcodes = self.0.remove(sample);
            assert!(barcodes.is_some());
            barcodes
        } else {
            None
        }
    }

    /// Return the samples.
    pub fn get_samples(&self) -> Vec<&SampleAssignment> {
        if self.is_multiplexed() {
            self.0.keys().sorted().collect()
        } else {
            vec![&SampleAssignment::NonMultiplexed]
        }
    }

    /// Return the samples and Unassigned.
    pub fn get_samples_with_unassigned(&self) -> Vec<&SampleAssignment> {
        if self.is_multiplexed() {
            let mut samples = self.get_samples();
            samples.push(&SampleAssignment::Unassigned);
            samples
        } else {
            vec![&SampleAssignment::NonMultiplexed]
        }
    }
}

/// A map of barcodes to sample assignments.
pub struct BarcodeToSample<'a> {
    barcode_to_sample: Option<TxHashMap<&'a Barcode, &'a SampleAssignment>>,
}

impl BarcodeToSample<'_> {
    /// Construct a map of barcodes to samples.
    pub fn construct(sample_barcodes: &SampleBarcodes) -> BarcodeToSample<'_> {
        if !sample_barcodes.is_multiplexed() {
            return BarcodeToSample {
                barcode_to_sample: None,
            };
        }

        let barcode_to_sample = sample_barcodes
            .0
            .iter()
            .flat_map(|(sample, barcodes)| barcodes.iter().map(move |barcode| (barcode, sample)))
            .collect();
        BarcodeToSample {
            barcode_to_sample: Some(barcode_to_sample),
        }
    }

    /// Return the sample assignment of a barcode.
    pub fn get_sample(&self, barcode: &Barcode) -> &SampleAssignment {
        let Some(barcode_to_sample) = &self.barcode_to_sample else {
            return &SampleAssignment::NonMultiplexed;
        };
        barcode_to_sample
            .get(barcode)
            .unwrap_or(&&SampleAssignment::Unassigned)
    }
}

enum_maker! {
    /// Type of targeting
    #[derive(MartianType)]
    TargetingMethod,
    (HybridCapture, "hybrid_capture"),
    (TemplatedLigation, "templated_ligation")
}

enum_maker! {
    /// Multiplexing information encoded at the level of reads.
    ReadLevel,
    (RTL, "RTL"),
    (OH, "OH")
}

enum_maker! {
    /// Multiplexing information encoded at the level of cells.
    CellLevel,
    (CMO, "CMO"),
    (Hashtag, "Hashtag")
}

#[derive(Debug, PartialEq, Eq, Clone, Copy, Ord, PartialOrd, Hash, Serialize, Deserialize)]
#[serde(untagged)]
pub enum BarcodeMultiplexingType {
    CellLevel(CellLevel),
    ReadLevel(ReadLevel),
}

impl BarcodeMultiplexingType {
    pub const CMO: Self = Self::CellLevel(CellLevel::CMO);
    #[allow(non_upper_case_globals)]
    pub const Hashtag: Self = Self::CellLevel(CellLevel::Hashtag);
    pub const RTL: Self = Self::ReadLevel(ReadLevel::RTL);
    pub const OH: Self = Self::ReadLevel(ReadLevel::OH);

    pub fn as_str(&self) -> &'static str {
        match self {
            BarcodeMultiplexingType::CellLevel(cl) => cl.as_str(),
            BarcodeMultiplexingType::ReadLevel(rl) => rl.as_str(),
        }
    }

    pub fn is_rtl(&self) -> bool {
        *self == Self::RTL
    }

    pub fn is_oh(&self) -> bool {
        *self == Self::OH
    }

    pub fn is_cmo(&self) -> bool {
        *self == Self::CMO
    }

    pub fn is_hashtag(&self) -> bool {
        *self == Self::Hashtag
    }
}

impl AsMartianPrimaryType for BarcodeMultiplexingType {
    fn as_martian_primary_type() -> MartianPrimaryType {
        MartianPrimaryType::Str
    }
}

enum_maker! {
    VdjChainType,
    (VdjT, "VDJ-T"),
    (VdjB, "VDJ-B"),
    (VdjTGD, "VDJ-T-GD"),
    (Auto, "VDJ")
}

#[allow(clippy::derivable_impls)]
impl Default for VdjChainType {
    fn default() -> Self {
        VdjChainType::Auto
    }
}

enum_maker! {
    #[derive(strum::Display, strum::EnumIter)]
    FeatureBarcodeType,
    (Antibody, "Antibody Capture"),
    (Antigen, "Antigen Capture"),
    (Crispr, "CRISPR Guide Capture"),
    (Multiplexing, "Multiplexing Capture", "FEATURETEST", "Multiplexing Tag Capture"),
    (Custom, "Custom")
}

impl FeatureBarcodeType {
    const ANTIBODY_SNAKE_CASE: &'static str = "antibody_capture";
    const ANTIGEN_SNAKE_CASE: &'static str = "antigen_capture";
    const CRISPR_SNAKE_CASE: &'static str = "crispr_guide_capture";
    const MUTLIPLEXING_SNAKE_CASE: &'static str = "multiplexing_capture";
    const CUSTOM_SNAKE_CASE: &'static str = "custom";

    /// Return a snake_case representation of this feature barcode type.
    pub fn as_snake_case(&self) -> &'static str {
        match self {
            Self::Antibody => Self::ANTIBODY_SNAKE_CASE,
            Self::Antigen => Self::ANTIGEN_SNAKE_CASE,
            Self::Crispr => Self::CRISPR_SNAKE_CASE,
            Self::Multiplexing => Self::MUTLIPLEXING_SNAKE_CASE,
            Self::Custom => Self::CUSTOM_SNAKE_CASE,
        }
    }

    /// Parse a snake-case representation of this feature barcode type.
    pub fn from_snake_case(s: &str) -> Result<Self> {
        Ok(match s {
            Self::ANTIBODY_SNAKE_CASE => Self::Antibody,
            Self::ANTIGEN_SNAKE_CASE => Self::Antigen,
            Self::CRISPR_SNAKE_CASE => Self::Crispr,
            Self::MUTLIPLEXING_SNAKE_CASE => Self::Multiplexing,
            Self::CUSTOM_SNAKE_CASE => Self::Custom,
            _ => {
                bail!("unable to parse '{s}' as a feature barcode type");
            }
        })
    }

    // Return the metric prefix for this feature barcode type as a static str.
    pub fn as_metric_prefix_static(&self) -> Option<&'static str> {
        match self {
            Self::Antibody => Some("ANTIBODY"),
            Self::Antigen => Some("ANTIGEN"),
            Self::Crispr => Some("CRISPR"),
            Self::Custom => Some("Custom"),
            Self::Multiplexing => Some("MULTIPLEXING"),
        }
    }
}

impl AsMetricPrefix for FeatureBarcodeType {
    /// Return the metric prefix for this feature barcode type.
    fn as_metric_prefix(&self) -> Option<&str> {
        self.as_metric_prefix_static()
    }
}

/// An enum which encapsulates both a library type and all the features
/// associated with it. All the gem wells in a single multi run should have
/// an identical set of `LibraryType`
#[derive(
    Debug, Default, Copy, Clone, PartialOrd, Ord, PartialEq, Eq, Hash, Serialize, DeserializeFromStr,
)]
#[serde(into = "&str")]
pub enum LibraryType {
    #[default]
    GeneExpression,
    Vdj(VdjChainType),
    FeatureBarcodes(FeatureBarcodeType),
    Atac,
}

impl Display for LibraryType {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

impl From<LibraryType> for &'static str {
    fn from(library_type: LibraryType) -> &'static str {
        library_type.as_str()
    }
}

impl FromStr for LibraryType {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(match s {
            Self::GENE_EXPRESSION_STR => Self::GeneExpression,
            Self::ATAC_STR => Self::Atac,
            _ => {
                if let Ok(fb) = FeatureBarcodeType::from_str(s) {
                    Self::FeatureBarcodes(fb)
                } else if let Ok(ct) = VdjChainType::from_str(s) {
                    Self::Vdj(ct)
                } else {
                    bail!("unable to parse {s} as LibraryType");
                }
            }
        })
    }
}

impl AsMartianPrimaryType for LibraryType {
    fn as_martian_primary_type() -> MartianPrimaryType {
        MartianPrimaryType::Str
    }
}

impl From<FeatureBarcodeType> for LibraryType {
    fn from(feature_type: FeatureBarcodeType) -> LibraryType {
        LibraryType::FeatureBarcodes(feature_type)
    }
}

#[allow(non_upper_case_globals)]
impl LibraryType {
    /// Alias for GeneExpression.
    pub const Gex: Self = Self::GeneExpression;
    /// Alias for antibody capture.
    pub const Antibody: Self = Self::FeatureBarcodes(FeatureBarcodeType::Antibody);
    /// Alias for antigen capture.
    pub const Antigen: Self = Self::FeatureBarcodes(FeatureBarcodeType::Antigen);
    /// Alias for CRISPR guide capture.
    pub const Crispr: Self = Self::FeatureBarcodes(FeatureBarcodeType::Crispr);
    /// Alias for cell multiplexing.
    pub const Cellplex: Self = Self::FeatureBarcodes(FeatureBarcodeType::Multiplexing);
    /// Alias for custom feature barcoding.
    pub const Custom: Self = Self::FeatureBarcodes(FeatureBarcodeType::Custom);
    /// Alias for VDJ Auto chain.
    pub const VdjAuto: Self = Self::Vdj(VdjChainType::Auto);

    const GENE_EXPRESSION_STR: &'static str = "Gene Expression";
    const ATAC_STR: &'static str = "Chromatin Accessibility";

    /// Return true if this is a gene expression library.
    pub fn is_gex(&self) -> bool {
        *self == Self::GeneExpression
    }

    /// Return true if this is a VDJ library.
    pub fn is_vdj(&self) -> bool {
        matches!(self, Self::Vdj(_))
    }

    /// Return true if this is a feature barcoding library.
    pub fn is_fb(&self) -> bool {
        matches!(self, Self::FeatureBarcodes(_))
    }

    /// Return true if this library type is feature barcoding of the specified type.
    pub fn is_fb_type(self, feature_barcode_type: FeatureBarcodeType) -> bool {
        self.feature_barcode_type() == Some(feature_barcode_type)
    }

    /// Return the feature barcoding type, if this is a feature barcoding library.
    pub fn feature_barcode_type(&self) -> Option<FeatureBarcodeType> {
        match self {
            Self::FeatureBarcodes(fb) => Some(*fb),
            Self::GeneExpression | Self::Vdj(_) | Self::Atac => None,
        }
    }

    /// Return the VDJ chain type, if this is a VDJ library.
    pub fn vdj_chain_type(&self) -> Option<VdjChainType> {
        match self {
            Self::Vdj(ct) => Some(*ct),
            Self::GeneExpression | Self::FeatureBarcodes(_) | Self::Atac => None,
        }
    }

    /// Return the string representation of this library type.
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::GeneExpression => Self::GENE_EXPRESSION_STR,
            Self::Atac => Self::ATAC_STR,
            Self::FeatureBarcodes(x) => x.as_str(),
            Self::Vdj(x) => x.as_str(),
        }
    }

    /// Return the metrix prefix for this library type as a static str.
    pub fn as_metric_prefix_static(&self) -> Option<&'static str> {
        match self {
            Self::GeneExpression => None,
            Self::FeatureBarcodes(fbt) => fbt.as_metric_prefix_static(),
            Self::Vdj(_) => Some("VDJ"),
            Self::Atac => Some("atac"),
        }
    }
}

impl AsMetricPrefix for LibraryType {
    fn as_metric_prefix(&self) -> Option<&str> {
        self.as_metric_prefix_static()
    }
}

/// A group of GEMs that were processed together, usually derived from a single
/// microfluidic channel
#[derive(Debug, Clone, Copy, PartialOrd, Ord, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct GemWell(pub u16);
impl AsMartianPrimaryType for GemWell {
    fn as_martian_primary_type() -> MartianPrimaryType {
        MartianPrimaryType::Int
    }
}
impl Default for GemWell {
    fn default() -> Self {
        GemWell(1)
    }
}
impl From<u16> for GemWell {
    fn from(gem_well: u16) -> Self {
        GemWell(gem_well)
    }
}
impl GemWell {
    pub fn inner(self) -> u16 {
        self.0
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Fastq {
    id: String,
    subsample_rate: Option<f64>,
    def: FastqDef,
}

impl Fastq {
    pub fn id(&self) -> &str {
        &self.id
    }
}

/// Pointer to all the FASTQs from a unique (library, gem_group) tuple.
/// Entries with the same `library_name` must be from the same physical sequencing
/// library, and must have been configured with the same `gem_well`, `library_type`,
/// and `feature_types`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct LibraryDef {
    pub physical_library_id: String,
    pub gem_well: GemWell,
    pub library_type: LibraryType,
    /// In CS mode, this will be `FastqDef::Bcl2Fastq(..)`. The `fastq_id` column corresponds to the
    /// `sample_name` in Bcl2Fastq
    pub fastqs: Vec<Fastq>, // Should be non empty
}

impl LibraryDef {
    fn push_fastq(&mut self, fastq_def: FastqDef, subsample_rate: Option<f64>) {
        let fastq_id = match fastq_def {
            FastqDef::Bcl2Fastq(ref def) => match &def.sample_name_spec {
                SampleNameSpec::Any => unreachable!(),
                SampleNameSpec::Names(names) => names.iter().exactly_one().unwrap().clone(),
            },
            FastqDef::BclProcessor(_) => {
                // TODO: Not sure if there is any other obvious name
                format!("{}_fastq_{}", self.physical_library_id, self.fastqs.len())
            }
        };
        self.fastqs.push(Fastq {
            id: fastq_id,
            subsample_rate,
            def: fastq_def,
        });
    }
}

/// A fingerprint refers to either:
/// - `Untagged(GemWell)`: All the cells in a gem well where the cells
///   are not tagged
/// - `Tagged(GemWell, SampleBarcodeID)`: A subset of cells in a gem well
///   which are tagged with a CMO with the given feature name
///   We can assume that no two samples within a multi run can share the same
///   `Fingerprint`.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub enum Fingerprint {
    Untagged {
        gem_well: GemWell,
    },
    Tagged {
        gem_well: GemWell,
        tag_name: SampleBarcodeID,
        /// Sample barcodes that were dynamically translated into this tag.
        translated_tag_names: Vec<SampleBarcodeID>,
        barcode_multiplexing_type: BarcodeMultiplexingType,
    },
}

impl Fingerprint {
    pub fn untagged(gem_well: GemWell) -> Self {
        Self::Untagged { gem_well }
    }

    pub fn tagged(
        gem_well: GemWell,
        tag_name: SampleBarcodeID,
        translated_tag_names: Vec<SampleBarcodeID>,
        barcode_multiplexing_type: BarcodeMultiplexingType,
    ) -> Self {
        Self::Tagged {
            gem_well,
            tag_name,
            translated_tag_names,
            barcode_multiplexing_type,
        }
    }

    /// Return true if this fingerprint is tagged.
    fn is_tagged(&self) -> bool {
        match *self {
            Fingerprint::Untagged { .. } => false,
            Fingerprint::Tagged { .. } => true,
        }
    }

    /// Return an iterator over all tag names.
    pub fn tag_names(&self) -> Box<dyn Iterator<Item = &SampleBarcodeID> + '_> {
        match self {
            Fingerprint::Untagged { .. } => Box::new(std::iter::empty()),
            Fingerprint::Tagged {
                tag_name,
                translated_tag_names,
                ..
            } => Box::new(std::iter::once(tag_name).chain(translated_tag_names)),
        }
    }

    pub fn gem_well(&self) -> GemWell {
        match *self {
            Fingerprint::Untagged { gem_well } | Fingerprint::Tagged { gem_well, .. } => gem_well,
        }
    }

    pub fn barcode_multiplexing_type(&self) -> Option<BarcodeMultiplexingType> {
        match self {
            Fingerprint::Untagged { .. } => None,
            Fingerprint::Tagged {
                barcode_multiplexing_type: cell_multiplexing_type,
                ..
            } => Some(*cell_multiplexing_type),
        }
    }

    pub fn is_cmo_multiplexed(&self) -> bool {
        self.barcode_multiplexing_type() == Some(BarcodeMultiplexingType::CellLevel(CellLevel::CMO))
    }

    pub fn is_hashtag_multiplexed(&self) -> bool {
        self.barcode_multiplexing_type()
            == Some(BarcodeMultiplexingType::CellLevel(CellLevel::Hashtag))
    }

    pub fn is_rtl_multiplexed(&self) -> bool {
        self.barcode_multiplexing_type() == Some(BarcodeMultiplexingType::ReadLevel(ReadLevel::RTL))
    }

    pub fn is_overhang_multiplexed(&self) -> bool {
        self.barcode_multiplexing_type() == Some(BarcodeMultiplexingType::ReadLevel(ReadLevel::OH))
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Sample {
    pub sample_id: String,
    pub description: String,
    /// set of places where I can find cells for this sample.
    /// it is invalid to have an empty fingerprint
    /// In 4.0, either all the fingerprints are untagged or all are cmo tagged
    /// No two samples can share the same fingerprint
    pub fingerprints: Vec<Fingerprint>, // Should not be empty
}

impl Sample {
    /// Returns the cell multiplexing type in use.
    /// Panics if more than one is found, which should be impossible.
    pub fn barcode_multiplexing_type(&self) -> Option<BarcodeMultiplexingType> {
        self.fingerprints
            .iter()
            .map(Fingerprint::barcode_multiplexing_type)
            .dedup()
            .exactly_one()
            .unwrap()
    }

    /// Return an iterator over all tag names in this sample.
    pub fn tag_names(&self) -> impl Iterator<Item = &SampleBarcodeID> {
        self.fingerprints.iter().flat_map(Fingerprint::tag_names)
    }
}

/// Incrementally construct the `CrMultiGraph` while tracking metadata to
/// perform validations as new entries as added.
pub struct CrMultiGraphBuilder {
    libraries: TxHashMap<String, LibraryDef>,
    // Set of library features of each gem well
    // In 4.0, all gem wells need to have identical set of library features
    gem_well_feat: TxHashMap<GemWell, TxHashSet<LibraryType>>,
    samples: TxHashMap<String, Sample>,
    fingerprints: TxHashSet<Fingerprint>,
    multiplexed_gem_wells: TxHashSet<GemWell>,
}

impl CrMultiGraphBuilder {
    pub fn new() -> Self {
        CrMultiGraphBuilder {
            libraries: TxHashMap::default(),
            gem_well_feat: TxHashMap::default(),
            samples: TxHashMap::default(),
            fingerprints: TxHashSet::default(),
            multiplexed_gem_wells: TxHashSet::default(),
        }
    }
    pub fn push_library(
        &mut self,
        physical_library_id: String,
        fastq: FastqDef,
        library_type: LibraryType,
        gem_well: GemWell,
    ) -> Result<()> {
        if library_type.is_fb_type(FeatureBarcodeType::Multiplexing) {
            self.multiplexed_gem_wells.insert(gem_well);
        }

        match self.libraries.entry(physical_library_id.clone()) {
            std::collections::hash_map::Entry::Vacant(e) => {
                let mut lib_def = LibraryDef {
                    physical_library_id,
                    gem_well,
                    library_type,
                    fastqs: vec![],
                };
                lib_def.push_fastq(fastq, None);
                e.insert(lib_def);
            }
            std::collections::hash_map::Entry::Occupied(mut e) => {
                let entry = e.get_mut();

                // Consistent checks: make sure there aren't conflicting entries for the same library
                if entry.gem_well != gem_well {
                    bail!(
                        "You supplied two different 'gem_well' values for physical_library_id = '{}'.\nEach physical_library_id must be associated with a single 'gem_well'. Supplied 'gem_well' values: '{}' and '{}'",
                        physical_library_id,
                        entry.gem_well.0,
                        gem_well.0
                    );
                }

                if entry.library_type != library_type {
                    bail!(
                        "You supplied two different 'feature_type' values for physical_library_id = {physical_library_id}.\nEach physical_library_id must be associated with a single 'feature_type'.",
                    );
                }

                // duplicate FASTQ def
                if entry.fastqs.iter().any(|fq| fq.def == fastq) {
                    bail!(
                        "You supplied duplicate FASTQ inforation different for physical_library_id = {physical_library_id}.\n You may have a redundant entry in your library specificaiton. Duplicate FASTQ definition:\n'{fastq:#?}",
                    );
                }

                entry.push_fastq(fastq, None);
            }
        }
        self.gem_well_feat
            .entry(gem_well)
            .or_default()
            .insert(library_type);

        Ok(())
    }

    // This needs to be called after populating all the libraries. This could be enforced
    // by a simple state machine, but it is currently not.
    pub fn push_fingerprint(
        &mut self,
        sample_id: String,
        description: String,
        fingerprint: Fingerprint,
    ) -> Result<()> {
        ensure!(
            !self.fingerprints.contains(&fingerprint),
            "Sample with sample_id = '{sample_id}' contains duplicated multiplexing information in the [samples] section: {fingerprint:?}"
        );
        if fingerprint.is_tagged() {
            // TODO: Good error message
            // Check that all the fingerprints are cmo tagged
            assert!(self.fingerprints.iter().all(Fingerprint::is_tagged));

            // make sure the GemWell containing this sample is multiplexed
            if fingerprint.is_cmo_multiplexed()
                && !self.multiplexed_gem_wells.contains(&fingerprint.gem_well())
            {
                bail!(
                    "Sample with sample_id = '{}' in gem_well = '{}' was declared as a multiplexed sample, but no 'Multiplexing Capture' library was provided for that gem_well in the [libraries] section",
                    sample_id,
                    &fingerprint.gem_well().0
                );
            }
        } else {
            // Check that either all the fingerprints are untagged
            assert!(self.fingerprints.iter().all(|fp| !fp.is_tagged()));
        }
        self.fingerprints.insert(fingerprint.clone());

        // Make sure that we have see this gem well
        assert!(self.gem_well_feat.contains_key(&fingerprint.gem_well())); // TODO: Good error message

        // TODO: Check that the description is consistent between two entries of a sample
        self.samples
            .entry(sample_id.clone())
            .or_insert(Sample {
                sample_id,
                description,
                fingerprints: Vec::new(),
            })
            .fingerprints
            .push(fingerprint);

        Ok(())
    }

    pub fn build(self) -> CrMultiGraph {
        assert!(
            self.multiplexed_gem_wells.is_empty()
                || self
                    .multiplexed_gem_wells
                    .iter()
                    .all(|gw| self.gem_well_feat.contains_key(gw))
        ); // TODO: Good error message

        // Make sure that there is at least one sample per gem well
        let sample_gem_wells: TxHashSet<_> = self
            .fingerprints
            .iter()
            .map(Fingerprint::gem_well)
            .collect();
        assert!(sample_gem_wells.len() == self.gem_well_feat.len());

        let mut libraries: Vec<_> = self.libraries.into_values().collect();
        libraries.sort_by(|a, b| {
            (a.gem_well, &a.physical_library_id).cmp(&(b.gem_well, &b.physical_library_id))
        });

        let mut samples: Vec<_> = self.samples.into_values().collect();
        samples.sort_by(|a, b| a.sample_id.cmp(&b.sample_id));

        CrMultiGraph { libraries, samples }
    }
}

impl Default for CrMultiGraphBuilder {
    fn default() -> Self {
        Self::new()
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CrMultiGraph {
    /// Incoming sequencing data
    pub libraries: Vec<LibraryDef>,

    /// Sample metadata table
    pub samples: Vec<Sample>,
}

impl CrMultiGraph {
    /// Return true if this analysis is using any form of multiplexing.
    pub fn is_multiplexed(&self) -> bool {
        self.barcode_multiplexing_type().is_some()
    }

    /// Return true if one or more libraries is of the specified type.
    pub fn has_library_type(&self, t: LibraryType) -> bool {
        self.libraries.iter().any(|lib| lib.library_type == t)
    }

    /// Return an iterator over all the library types in the experiment
    pub fn library_types(&self) -> impl Iterator<Item = LibraryType> + '_ {
        self.libraries.iter().map(|lib| lib.library_type)
    }

    /// Return the cell multiplexing type in use, if there is one.
    pub fn barcode_multiplexing_type(&self) -> Option<BarcodeMultiplexingType> {
        self.samples
            .iter()
            .map(Sample::barcode_multiplexing_type)
            .dedup()
            .exactly_one()
            .unwrap()
    }

    /// Returns if the multiplexing type is read-level or not.
    pub fn is_read_level_multiplexed(&self) -> bool {
        matches!(
            self.barcode_multiplexing_type(),
            Some(BarcodeMultiplexingType::ReadLevel(_))
        )
    }

    /// Returns if the multiplexing type is cell-level or not.
    pub fn is_cell_level_multiplexed(&self) -> bool {
        matches!(
            self.barcode_multiplexing_type(),
            Some(BarcodeMultiplexingType::CellLevel(_))
        )
    }

    /// Returns a map from demux tag name (e.g. OCM barcode ID or probe barcode ID) to sample ID information.
    pub fn get_tag_name_to_sample_id_map(&self) -> Result<TxHashMap<&str, (&str, &[String])>> {
        self.samples
            .iter()
            .flat_map(|sample| {
                sample.fingerprints.iter().map(|fingerprint| {
                    let Fingerprint::Tagged {
                        tag_name,
                        translated_tag_names,
                        ..
                    } = fingerprint
                    else {
                        bail!("found an untagged fingerprint when creating tag name map: {fingerprint:?}");
                    };
                    Ok((
                        tag_name.as_str(),
                        (sample.sample_id.as_str(), translated_tag_names.as_slice()),
                    ))
                })
            })
            .collect()
    }
}

#[derive(Deserialize, Serialize, Debug, PartialOrd, Ord, PartialEq, Eq)]
#[serde(untagged)]
pub enum LibraryInfo {
    Count(crate::rna_read::LibraryInfo),
    Aggr(crate::aggr::LibraryInfo),
}

impl From<crate::aggr::LibraryInfo> for LibraryInfo {
    fn from(v: crate::aggr::LibraryInfo) -> Self {
        LibraryInfo::Aggr(v)
    }
}

impl From<crate::rna_read::LibraryInfo> for LibraryInfo {
    fn from(v: crate::rna_read::LibraryInfo) -> Self {
        LibraryInfo::Count(v)
    }
}

impl TryFrom<LibraryInfo> for crate::rna_read::LibraryInfo {
    type Error = anyhow::Error;

    fn try_from(value: LibraryInfo) -> Result<Self> {
        match value {
            LibraryInfo::Count(v) => Ok(v),
            LibraryInfo::Aggr(_) => todo!(),
        }
    }
}

impl TryFrom<LibraryInfo> for crate::aggr::LibraryInfo {
    type Error = anyhow::Error;

    fn try_from(value: LibraryInfo) -> Result<Self> {
        match value {
            LibraryInfo::Aggr(v) => Ok(v),
            LibraryInfo::Count(_) => todo!(),
        }
    }
}

/// The aligner used to align the reads to the reference.
#[derive(Clone, Copy, Debug, PartialEq, Eq, MartianType, EnumString, Deserialize, Serialize)]
#[serde(rename_all = "snake_case")]
#[strum(serialize_all = "snake_case")]
pub enum AlignerParam {
    /// 10x probe aligner used for RTL chemistries
    Hurtle,
    /// STAR (Spliced Transcripts Alignment to a Reference)
    Star,
    /// Minimap2
    Minimap2,
}

impl AlignerParam {
    /// The minimum MAPQ threshold for a confidently-mapped read.
    pub fn high_conf_mapq(&self) -> u8 {
        match self {
            AlignerParam::Hurtle => 255,
            AlignerParam::Star => 255,
            AlignerParam::Minimap2 => 60,
        }
    }
}

/// GeneticDemuxParams struct
#[derive(Clone, Serialize, Deserialize, MartianStruct, Debug)]
pub struct GeneticDemuxParams {
    /// The number of samples to demux into
    pub num_clusters: usize,
    /// Candidate SNP VCF file
    pub candidate_snps: Vcf,
    /// known genotypes
    pub known_genotypes: Option<KnownGenotypes>,
    /// The minimum number of UMIs supporting the reference allele
    pub min_ref_umi: Option<usize>,
    /// The minimum number of UMIs supporting the alternate allele
    pub min_alt_umi: Option<usize>,
    /// The minimum number of cells containing ref allele for the variant to be used for clustering
    pub min_ref_cells: Option<usize>,
    /// The minimum number of cells containing alt allele for the variant to be used for clustering
    pub min_alt_cells: Option<usize>,
    /// The number of random seedings
    pub restarts: Option<usize>,
}

/// Known genotype information
#[derive(Clone, Serialize, Deserialize, MartianStruct, Debug)]
pub struct KnownGenotypes {
    /// Known genotypes VCF file
    pub vcf: Vcf,
    /// Sample names
    pub sample_names: Vec<String>,
}

// VCF file
martian_filetype!(Vcf, "vcf");

// Matrix Market file
martian_filetype!(Mtx, "mtx");

// FASTA file
martian_filetype!(FastaFile, "fa");

// Shardio file containing BAM Record objects from rust_htslib,
// in position-sorted order
martian_filetype!(AlignShardFile, "asf");

// BAM file's SAM header
martian_filetype!(SamHeaderFile, "sam");

/// BAM contig information
pub struct BamContigInfo {
    /// Contig name
    pub name: String,
    /// Contig target ID (0-based index in the BAM header)
    pub tid: u32,
    /// Contig length
    pub length: u64,
}

impl SamHeaderFile {
    /// Get a vector of contig information
    pub fn contig_info(&self) -> impl Iterator<Item = BamContigInfo> {
        let bam_reader = BamReader::from_path(self).unwrap();
        let header = bam_reader.header();
        #[expect(clippy::needless_collect)]
        let bam_contig_infos: Vec<_> = header
            .target_names()
            .into_iter()
            .map(move |s| {
                let tid = header.tid(s).unwrap();
                let name = String::from_utf8(s.to_vec()).unwrap();
                let length = header.target_len(tid).unwrap();
                BamContigInfo { name, tid, length }
            })
            .collect();
        bam_contig_infos.into_iter()
    }

    /// Get a map of each contig name to its target ID and length
    pub fn contig_info_by_name(&self) -> TxHashMap<String, BamContigInfo> {
        self.contig_info()
            .map(|info| (info.name.clone(), info))
            .collect()
    }

    /// Get a map of each target ID to contig name and length
    pub fn contig_info_by_tid(&self) -> TxHashMap<u32, BamContigInfo> {
        self.contig_info().map(|info| (info.tid, info)).collect()
    }
}

/// Wrap around API calls to read/write/parse
impl Vcf {
    // TODO: Implement parse/read/write functions
}

/// Wrap around API calls to read/write/parse
impl Mtx {
    // TODO: Implement parse/read/write functions
}

#[cfg(test)]
mod tests {
    use super::*;
    use barcode::BcSeq;
    use fastq_set::filenames::LaneSpec;
    use insta::assert_json_snapshot;
    use std::fs::write;
    use tempfile::NamedTempFile;

    fn barcode(seq: &[u8]) -> Barcode {
        Barcode::with_seq(1, BcSeq::from_bytes(seq), true)
    }

    #[test]
    /// Test SampleBarcodes::read_from_json
    fn test_sample_barcodes_read_from_json() {
        let tempfile = NamedTempFile::with_suffix(".json").unwrap();
        let file_path = tempfile.path();
        let json_string = r#"{
            "Sample1": [
                "AAACTGGGTAAAGGCCATCCCAAC-1",
                "AAAGCGAAGACCATGGATCCCAAC-1"
            ],
            "Sample2": [
                "TTGAGCGAGTCAACTTAACGCCGA-1",
                "TTGAGGTAGGGATGAAAACGCCGA-1",
                "TTGAGGAGTTAGTGAGAACGCCGA-1"
            ]
        }"#;
        write(file_path, json_string).unwrap();

        let sample_barcode_expected = SampleBarcodes(TxHashMap::from_iter([
            (
                SampleAssignment::from("Sample1".to_string()),
                vec![
                    barcode(b"AAACTGGGTAAAGGCCATCCCAAC"),
                    barcode(b"AAAGCGAAGACCATGGATCCCAAC"),
                ],
            ),
            (
                SampleAssignment::from("Sample2".to_string()),
                vec![
                    barcode(b"TTGAGCGAGTCAACTTAACGCCGA"),
                    barcode(b"TTGAGGTAGGGATGAAAACGCCGA"),
                    barcode(b"TTGAGGAGTTAGTGAGAACGCCGA"),
                ],
            ),
        ]));

        assert_eq!(
            SampleBarcodes::read_from_json(Some(&SampleBarcodesFile::from(file_path))).unwrap(),
            sample_barcode_expected,
        );
    }

    fn untagged_single_gw_graph() -> Result<CrMultiGraph> {
        // 1 sample, 1 gem well, GEX + VDJ + Ab, 1 sequencing per library
        let mut builder = CrMultiGraphBuilder::new();
        builder.push_library(
            "GEX_LIB".into(),
            FastqDef::bcl2fastq(
                "/path/to/gex/fastq/".into(),
                "my_gex_data".into(),
                LaneSpec::Any,
            ),
            LibraryType::Gex,
            GemWell(1),
        )?;

        builder.push_library(
            "VDJ_LIB".into(),
            FastqDef::bcl2fastq(
                "/path/to/vdj/fastq/".into(),
                "my_vdj_data".into(),
                LaneSpec::Any,
            ),
            LibraryType::VdjAuto,
            GemWell(1),
        )?;

        builder.push_library(
            "AB_LIB".into(),
            FastqDef::bcl2fastq(
                "/path/to/ab/fastq/".into(),
                "my_ab_data".into(),
                LaneSpec::Any,
            ),
            LibraryType::Antibody,
            GemWell(1),
        )?;

        builder.push_fingerprint(
            "PBMC".into(),
            "10k Human PBMC".into(),
            Fingerprint::untagged(GemWell(1)),
        )?;

        Ok(builder.build())
    }

    fn tagged_multi_gw_graph() -> Result<CrMultiGraph> {
        // 2 samples, 1 tag per sample, 2 gem wells, GEX + CMO, 1 sequencing per library
        let mut builder = CrMultiGraphBuilder::new();
        builder.push_library(
            "GEX_LIB_1".into(),
            FastqDef::bcl2fastq(
                "/path/to/gex_1/fastq/".into(),
                "my_gex_1_data".into(),
                LaneSpec::Any,
            ),
            LibraryType::Gex,
            GemWell(1),
        )?;
        builder.push_library(
            "GEX_LIB_2".into(),
            FastqDef::bcl2fastq(
                "/path/to/gex_2/fastq/".into(),
                "my_gex_2_data".into(),
                LaneSpec::Any,
            ),
            LibraryType::Gex,
            GemWell(2),
        )?;
        builder.push_library(
            "CMO_LIB_1".into(),
            FastqDef::bcl2fastq(
                "/path/to/cmo_1/fastq/".into(),
                "my_cmo_1_data".into(),
                LaneSpec::Any,
            ),
            LibraryType::Cellplex,
            GemWell(1),
        )?;
        builder.push_library(
            "CMO_LIB_2".into(),
            FastqDef::bcl2fastq(
                "/path/to/cmo_2/fastq/".into(),
                "my_cmo_2_data".into(),
                LaneSpec::Any,
            ),
            LibraryType::Cellplex,
            GemWell(2),
        )?;

        builder.push_fingerprint(
            "PBMC_donor_A".into(),
            "Human PBMC from donor A".into(),
            Fingerprint::tagged(
                GemWell(1),
                "CMO500".into(),
                Vec::default(),
                BarcodeMultiplexingType::CellLevel(CellLevel::CMO),
            ),
        )?;

        builder.push_fingerprint(
            "PBMC_donor_A".into(),
            "Human PBMC from donor A".into(),
            Fingerprint::tagged(
                GemWell(2),
                "CMO500".into(),
                Vec::default(),
                BarcodeMultiplexingType::CellLevel(CellLevel::CMO),
            ),
        )?;

        builder.push_fingerprint(
            "PBMC_donor_B".into(),
            "Human PBMC from donor B".into(),
            Fingerprint::tagged(
                GemWell(1),
                "CMO501".into(),
                Vec::default(),
                BarcodeMultiplexingType::CellLevel(CellLevel::CMO),
            ),
        )?;

        builder.push_fingerprint(
            "PBMC_donor_B".into(),
            "Human PBMC from donor B".into(),
            Fingerprint::tagged(
                GemWell(2),
                "CMO501".into(),
                Vec::default(),
                BarcodeMultiplexingType::CellLevel(CellLevel::CMO),
            ),
        )?;

        Ok(builder.build())
    }

    #[test]
    fn test_py_api_graph_snapshots() {
        insta::with_settings!({snapshot_path => "snapshots/py_api"}, {
            assert_json_snapshot!("multi_graph_untagged", untagged_single_gw_graph().unwrap());
            assert_json_snapshot!("multi_graph_tagged", tagged_multi_gw_graph().unwrap());
        });
    }
}
