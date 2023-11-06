use crate::barcode_index::BarcodeIndex;
use crate::bit_encode::BitEncoded;
use crate::rna_read::LegacyLibraryType;
use anyhow::{bail, ensure, Context, Result};
use barcode::{Barcode, BarcodeConstruct, BarcodeFromString, BcSegSeq};
use fastq_set::filenames::bcl2fastq::SampleNameSpec;
use fastq_set::filenames::FastqDef;
use itertools::Itertools;
use martian::{AsMartianPrimaryType, MartianFileType, MartianPrimaryType};
use martian_derive::{martian_filetype, MartianStruct, MartianType};
use martian_filetypes::bin_file::BinaryFormat;
use martian_filetypes::json_file::{JsonFile, JsonFormat};
use martian_filetypes::FileTypeRead;
use metric::{JsonReport, JsonReporter, SimpleHistogram, TxHashMap, TxHashSet};
use serde::de::{IgnoredAny, Visitor};
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use shardio::SortKey;
use std::borrow::Cow;
use std::convert::TryFrom;
use std::fmt::{self, Formatter};
use std::fs::read_to_string;
use std::hash::Hash;
use std::path::PathBuf;
use std::str::FromStr;
use strum::IntoEnumIterator;
use strum_macros::{Display, EnumCount, EnumIter, EnumString};
use umi::UmiType;

// Used to indicate a UMI count does not come from a probe
// in datasets that are mixed probe/no-probe (e.g. RTL + FB)
pub const PROBE_IDX_SENTINEL_VALUE: i32 = -1;

// File Types

martian_filetype! { BcSegmentCountFile, "bsc" }
pub type BcSegmentCountFormat = BinaryFormat<
    BcSegmentCountFile,
    TxHashMap<LibraryFeatures, BarcodeConstruct<SimpleHistogram<BcSegSeq>>>,
>;

martian_filetype! { BcCountFile, "bcc" }
pub type BcCountDataType = TxHashMap<LibraryFeatures, SimpleHistogram<Barcode>>;
pub type BcCountFormat = BinaryFormat<BcCountFile, BcCountDataType>;

martian_filetype! { TotalBcCountFile, "tbcc" }
pub type TotalBcCountDataType = SimpleHistogram<Barcode>;
pub type TotalBcCountFormat = BinaryFormat<TotalBcCountFile, TotalBcCountDataType>;

martian_filetype!(FeatureCountFile, "fbc");
pub type FeatureCountFormat = BinaryFormat<FeatureCountFile, Vec<i64>>;

martian_filetype!(BarcodeSetFile, "blf");
pub type BarcodeSetFormat = JsonFormat<BarcodeSetFile, TxHashSet<Barcode>>;

martian_filetype!(BarcodeIndexFile, "bi");
pub type BarcodeIndexFormat<B> = BinaryFormat<BarcodeIndexFile, BarcodeIndex<B>>;

// End File Tyes
/// A trait for reads that may have a sample index.
pub trait HasSampleIndex {
    fn si_seq(&self) -> Option<&[u8]>;
    fn si_qual(&self) -> Option<&[u8]>;
}

/// Count of the number of UMIs observed for one feature in one barcode.  Corresponds
/// to a single entry in the feature x barcode matrix.
#[derive(Serialize, Deserialize, Clone, Copy, Ord, PartialOrd, Eq, PartialEq)]
pub struct FeatureBarcodeCount<B> {
    pub barcode: B,
    pub feature_idx: u32,
    pub umi_count: u32,
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
    pub fn feature_counts(&self) -> impl Iterator<Item = FeatureBarcodeCount<Barcode>> + '_ {
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

#[derive(Debug, Clone, Hash, PartialEq, Eq, Ord, PartialOrd, Copy, Serialize, Deserialize)]
pub enum ReqStrand {
    #[serde(rename = "+")]
    Forward,
    #[serde(rename = "-")]
    Reverse,
}

pub trait HasMultiplexing {
    fn has_multiplexing(&self) -> bool;
}

macro_rules! enum_maker {
    ($(#[$meta:meta])* $name:ident, $( ($field:ident, $lit: literal $(, $alias: literal)*) ),*) => {
        #[derive(EnumCount, EnumString, EnumIter, Display, Debug, PartialEq, Eq, Hash, Clone, Copy, PartialOrd, Ord, Serialize, Deserialize, MartianType)]
        pub enum $name {
            $(
                #[strum(to_string = $lit)]
                #[serde(rename = $lit)]
                $(#[serde(alias = $alias)]
                #[strum(serialize = $alias)])*
                $field,
            )*
        }
        impl From<u8> for $name {
            fn from(val: u8) -> Self {
                match $name::iter().skip(val as usize).next() {
                    Some(v) => v,
                    None => panic!("Cannot convert {} to {}", val, stringify!($name)),
                }
            }
        }
        impl From<$name> for u8 {
            fn from(val: $name) -> u8 {
                val as u8
            }
        }
    };
}

pub type SampleBarcodeID = String;
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

impl AsMartianPrimaryType for SampleAssignment {
    fn as_martian_primary_type() -> MartianPrimaryType {
        MartianPrimaryType::Str
    }
}

impl fmt::Display for SampleAssignment {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = match self {
            SampleAssignment::NonMultiplexed => "non_multiplexed".to_string(),
            SampleAssignment::Unassigned => "unassigned".to_string(),
            SampleAssignment::Assigned(sample_id) => sample_id.to_string(),
        };

        write!(f, "{s}")
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
        impl<'de> Visitor<'de> for SampleAssignmentVisitor {
            type Value = SampleAssignment;
            fn expecting(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
                f.write_str(
                    "A string of either the assigned sample ID, \"unassigned\", or \"non_multiplexed\"",
                )
            }

            fn visit_str<E>(self, s: &str) -> Result<Self::Value, E>
            where
                E: serde::de::Error,
            {
                match s {
                    "unassigned" | "Unassigned" => Ok(SampleAssignment::Unassigned),
                    "non_multiplexed" | "Non_multiplexed" => Ok(SampleAssignment::NonMultiplexed),
                    _ => Ok(SampleAssignment::Assigned(s.to_string())),
                }
            }
        }
        deserializer.deserialize_any(SampleAssignmentVisitor)
    }
}

/// describes whether a molecule info is raw (uber), sample, or count
/// stored in molecule info metrics and used to set aggr preflights and etc
#[derive(Debug, Clone, Hash, PartialEq, Eq, Ord, PartialOrd, Copy, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum MoleculeInfoType {
    Raw,
    PerSample,
    Count,
}

pub type SampleBarcodesFile = JsonFile<SampleBarcodes>;

/// A map of sample IDs to their barcodes.
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
    TargetingMethod,
    (HybridCapture, "hybrid_capture"),
    (TemplatedLigation, "templated_ligation")
}

enum_maker! {
    /// Type of multiplexing, CMO or RTL.
    CellMultiplexingType,
    (CMO, "CMO"),
    (RTL, "RTL"),
    (OH, "OH")
}

enum_maker! {
    /// We want the users to have the option to explicitly specify whether
    /// a VDJ library was TCR or IG enriched. We would auto detect the enrichment by
    /// default. However, in some rare cases when the quality of the library is
    /// poor, the auto detection fails and the pipeline would be unable to proceed.
    /// So we are allowing the users to optionally specify `T Cell` or `B Cell` as a
    /// `feature_type` when the `library_type` is VDJ. If unspecified, the default value
    /// of `Auto` will be used.
    VdjChainType,
    (Tcr, "T Cell Receptor"),
    (Ig, "B Cell Receptor"),
    (Auto, "Auto")
}

#[allow(clippy::derivable_impls)]
impl Default for VdjChainType {
    fn default() -> Self {
        VdjChainType::Auto
    }
}

enum_maker! {
    /// The feature types supported in our matrix
    FeatureType,
    (Antibody, "Antibody Capture"),
    (Antigen, "Antigen Capture"),
    (CRISPR, "CRISPR Guide Capture"),
    (Multiplexing, "Multiplexing Capture", "FEATURETEST", "Multiplexing Tag Capture"),
    (Custom, "Custom"),
    (Gene, "Gene Expression")
}

impl FeatureType {
    /// Return an iterator over all feature barcode library types.
    pub fn feature_barcodes() -> impl Iterator<Item = FeatureType> {
        FeatureType::iter().filter(|&f| f != FeatureType::Gene)
    }

    /// Return the metric prefix for this feature type, and None for GEX.
    pub fn metric_prefix(self) -> Option<&'static str> {
        #[allow(clippy::enum_glob_use)]
        use FeatureType::*;
        match self {
            Gene => None,
            Antibody => Some("ANTIBODY"),
            Antigen => Some("ANTIGEN"),
            CRISPR => Some("CRISPR"),
            Custom => Some("Custom"),
            Multiplexing => Some("MULTIPLEXING"),
        }
    }

    /// Join the metric prefix and the metric name, separated by an underscore.
    pub fn join<'a>(&self, metric_name: &'a str) -> Cow<'a, str> {
        if let Some(prefix) = self.metric_prefix() {
            Cow::Owned(format!("{prefix}_{metric_name}"))
        } else {
            Cow::Borrowed(metric_name)
        }
    }
}

impl TryFrom<LegacyLibraryType> for FeatureType {
    type Error = LegacyLibraryType;

    fn try_from(legacy_library_type: LegacyLibraryType) -> Result<FeatureType, LegacyLibraryType> {
        match legacy_library_type {
            LegacyLibraryType::GeneExpression => Ok(FeatureType::Gene),
            LegacyLibraryType::AntibodyCapture => Ok(FeatureType::Antibody),
            LegacyLibraryType::AntigenCapture => Ok(FeatureType::Antigen),
            LegacyLibraryType::CrisprGuideCapture => Ok(FeatureType::CRISPR),
            LegacyLibraryType::Custom => Ok(FeatureType::Custom),
            LegacyLibraryType::Multiplexing => Ok(FeatureType::Multiplexing),
            LegacyLibraryType::ATAC | LegacyLibraryType::Vdj => Err(legacy_library_type),
        }
    }
}

impl HasMultiplexing for FeatureType {
    fn has_multiplexing(&self) -> bool {
        *self == FeatureType::Multiplexing
    }
}

/// An enum which encapsulates both a library type and all the features
/// associated with it. All the gem wells in a single multi run should have
/// an identical set of `LibraryFeatures`
#[derive(
    Debug,
    Display,
    Copy,
    Clone,
    PartialOrd,
    Ord,
    PartialEq,
    Eq,
    Hash,
    Serialize,
    Deserialize,
    MartianType,
)]
pub enum LibraryFeatures {
    /// The inner feature type should be `FeatureType::Gene`. It is natural to ask
    /// why we even need to store a `FeatureType` within this enum variant. The reason is
    /// doing so would ensure that all the variants in the enum would serialize into
    /// a map. This makes `LibraryFeatures` a compatible type for martian stage args/outs
    GeneExpression(FeatureType),
    Vdj(VdjChainType),
    FeatureBarcodes(BitEncoded<FeatureType>),
}

impl From<FeatureType> for LibraryFeatures {
    fn from(feature_type: FeatureType) -> LibraryFeatures {
        use FeatureType::{Antibody, Antigen, Custom, Gene, Multiplexing, CRISPR};
        match feature_type {
            Gene => LibraryFeatures::GeneExpression(feature_type),
            Antibody | Antigen | CRISPR | Multiplexing | Custom => {
                LibraryFeatures::FeatureBarcodes(feature_type.into())
            }
        }
    }
}

impl From<LegacyLibraryType> for LibraryFeatures {
    fn from(legacy_library_type: LegacyLibraryType) -> LibraryFeatures {
        match legacy_library_type {
            LegacyLibraryType::GeneExpression => LibraryFeatures::GeneExpression(FeatureType::Gene),
            LegacyLibraryType::AntibodyCapture => {
                LibraryFeatures::FeatureBarcodes(FeatureType::Antibody.into())
            }
            LegacyLibraryType::AntigenCapture => {
                LibraryFeatures::FeatureBarcodes(FeatureType::Antigen.into())
            }
            LegacyLibraryType::CrisprGuideCapture => {
                LibraryFeatures::FeatureBarcodes(FeatureType::CRISPR.into())
            }
            LegacyLibraryType::Vdj => LibraryFeatures::Vdj(VdjChainType::Auto),
            LegacyLibraryType::Custom => {
                LibraryFeatures::FeatureBarcodes(FeatureType::Custom.into())
            }
            LegacyLibraryType::Multiplexing => {
                LibraryFeatures::FeatureBarcodes(FeatureType::Multiplexing.into())
            }
            LegacyLibraryType::ATAC => panic!("Unexpected library type 'Chromatin Accessibility'"),
        }
    }
}

impl HasMultiplexing for LibraryFeatures {
    fn has_multiplexing(&self) -> bool {
        match *self {
            LibraryFeatures::FeatureBarcodes(f) => f.iter().any(|f| f.has_multiplexing()),
            _ => false,
        }
    }
}

impl LibraryFeatures {
    pub fn gex() -> Self {
        LibraryFeatures::GeneExpression(FeatureType::Gene)
    }
    pub fn has_feature(self, feature_type: FeatureType) -> bool {
        match self {
            LibraryFeatures::GeneExpression(f) => f == feature_type,
            LibraryFeatures::FeatureBarcodes(feats) => feats.iter().any(|f| f == feature_type),
            _ => false,
        }
    }
    pub fn library_type(self) -> LibraryType {
        match self {
            LibraryFeatures::GeneExpression(_) => LibraryType::GeneExpression,
            LibraryFeatures::Vdj(_) => LibraryType::ImmuneProfiling,
            LibraryFeatures::FeatureBarcodes(_) => LibraryType::FeatureBarcoding,
        }
    }
    pub fn legacy_library_type(self) -> LegacyLibraryType {
        match self {
            LibraryFeatures::GeneExpression(_) => LegacyLibraryType::GeneExpression,
            LibraryFeatures::Vdj(_) => LegacyLibraryType::Vdj,
            LibraryFeatures::FeatureBarcodes(features) => {
                match features.iter().exactly_one().unwrap() {
                    FeatureType::Antibody => LegacyLibraryType::AntibodyCapture,
                    FeatureType::Antigen => LegacyLibraryType::AntigenCapture,
                    FeatureType::CRISPR => LegacyLibraryType::CrisprGuideCapture,
                    FeatureType::Custom => LegacyLibraryType::Custom,
                    FeatureType::Multiplexing => LegacyLibraryType::Multiplexing,
                    _ => panic!(),
                }
            }
        }
    }
}

enum_maker! {
    /// TODO: Should targeting be a different library type?
    LibraryType,
    (GeneExpression, "Gene Expression"),
    (ImmuneProfiling, "Immune Profiling"),
    (FeatureBarcoding, "Feature Barcoding")
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

/// Pointer to all the FASTQs from a unique (library, gem_group) tuple.
/// Entries with the same `library_name` must be from the same physical sequencing
/// library, and must have been configured with the same `gem_well`, `library_type`,
/// and `feature_types`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct LibraryDef {
    pub physical_library_id: String,
    pub gem_well: GemWell,
    pub library_features: LibraryFeatures,
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
    pub fn legacy_library_type(&self) -> LegacyLibraryType {
        self.library_features.legacy_library_type()
    }
}

impl HasMultiplexing for LibraryDef {
    /// Return true if this is a FeatureBarcoding library
    /// with `Multiplexing` feature type
    fn has_multiplexing(&self) -> bool {
        self.library_features.has_multiplexing()
    }
}

impl LibraryDef {
    pub fn library_type(&self) -> LibraryType {
        self.library_features.library_type()
    }
    // pub fn feature_types(&self) -> Option<impl Iterator<Item = FeatureType>> {
    //     self.library_features.feature_types()
    // }
}

/// A fingerprint refers to either:
/// - `Untagged(GemWell)`: All the cells in a gem well where the cells
///   are not tagged
/// - `Tagged(GemWell, SampleBarcodeID)`: A subset of cells in a gem well
///   which are tagged with a CMO with the given feature name
/// We can assume that no two samples within a multi run can share the same
/// `Fingerprint`.
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
        cell_multiplexing_type: CellMultiplexingType,
    },
}

impl HasMultiplexing for Fingerprint {
    /// Return true if this fingerprint is CMO or RTL tagged
    fn has_multiplexing(&self) -> bool {
        match *self {
            Fingerprint::Untagged { .. } => false,
            Fingerprint::Tagged { .. } => true,
        }
    }
}

impl Fingerprint {
    pub fn untagged(gem_well: GemWell) -> Self {
        Self::Untagged { gem_well }
    }

    pub fn tagged(
        gem_well: GemWell,
        tag_name: SampleBarcodeID,
        translated_tag_names: Vec<SampleBarcodeID>,
        cell_multiplexing_type: CellMultiplexingType,
    ) -> Self {
        Self::Tagged {
            gem_well,
            tag_name,
            translated_tag_names,
            cell_multiplexing_type,
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

    pub fn cell_multiplexing_type(&self) -> Option<CellMultiplexingType> {
        match *self {
            Fingerprint::Untagged { .. } => None,
            Fingerprint::Tagged {
                cell_multiplexing_type,
                ..
            } => Some(cell_multiplexing_type),
        }
    }

    pub fn is_cmo_multiplexed(&self) -> bool {
        self.cell_multiplexing_type() == Some(CellMultiplexingType::CMO)
    }

    pub fn is_rtl_multiplexed(&self) -> bool {
        self.cell_multiplexing_type() == Some(CellMultiplexingType::RTL)
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

impl HasMultiplexing for Sample {
    /// Returns `true` if all the `Fingerprint`s are CMO tagged
    /// Returns `false` if all the `Fingerprint`s are untagged
    /// Panics otherwise
    fn has_multiplexing(&self) -> bool {
        self.fingerprints
            .iter()
            .map(HasMultiplexing::has_multiplexing)
            .dedup()
            .exactly_one()
            .unwrap()
    }
}

impl Sample {
    /// Returns the cell multiplexing type in use.
    /// Panics if more than one is found, which should be impossible.
    pub fn cell_multiplexing_type(&self) -> Option<CellMultiplexingType> {
        self.fingerprints
            .iter()
            .map(Fingerprint::cell_multiplexing_type)
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
    gem_well_feat: TxHashMap<GemWell, TxHashSet<LibraryFeatures>>,
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
        library_features: LibraryFeatures,
        gem_well: GemWell,
    ) -> Result<()> {
        if library_features.has_multiplexing() {
            self.multiplexed_gem_wells.insert(gem_well);
        }

        match self.libraries.entry(physical_library_id.clone()) {
            std::collections::hash_map::Entry::Vacant(e) => {
                let mut lib_def = LibraryDef {
                    physical_library_id,
                    gem_well,
                    library_features,
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

                if entry.library_features != library_features {
                    bail!(
                        "You supplied two different 'feature_type' values for physical_library_id = {}.\nEach physical_library_id must be associated with a single 'feature_type'.",
                        physical_library_id
                    );
                }

                // duplicate FASTQ def
                if entry.fastqs.iter().any(|fq| fq.def == fastq) {
                    bail!(
                        "You supplied duplicate FASTQ inforation different for physical_library_id = {}.\n You may have a redundant entry in your library specificaiton. Duplicate FASTQ definition:\n'{:#?}",
                        physical_library_id,
                        fastq
                    );
                }

                entry.push_fastq(fastq, None);
            }
        }
        self.gem_well_feat
            .entry(gem_well)
            .or_insert_with(TxHashSet::default)
            .insert(library_features);

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
            "Sample with sample_id = '{}' contains duplicated multiplexing information in the [samples] section: {:?}",
            sample_id,
            fingerprint,
        );
        if fingerprint.has_multiplexing() {
            // TODO: Good error message
            // Check that all the fingerprints are cmo tagged
            assert!(self
                .fingerprints
                .iter()
                .all(HasMultiplexing::has_multiplexing));

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
            assert!(self.fingerprints.iter().all(|fp| !fp.has_multiplexing()));
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
    pub fn get_multiplexing_method(&self) -> Option<CellMultiplexingType> {
        match self.has_multiplexing() {
            false => None,
            true => Some(
                self.samples
                    .iter()
                    .flat_map(|sample| {
                        sample
                            .fingerprints
                            .iter()
                            .map(|fingerprint| fingerprint.cell_multiplexing_type().unwrap())
                            .collect::<Vec<CellMultiplexingType>>()
                    })
                    .dedup()
                    .exactly_one()
                    .unwrap(),
            ),
        }
    }

    /// Return true if one or more libraries is of the specified type.
    pub fn has_legacy_library_type(&self, t: LegacyLibraryType) -> bool {
        self.libraries
            .iter()
            .any(|lib| lib.legacy_library_type() == t)
    }
}

impl HasMultiplexing for CrMultiGraph {
    fn has_multiplexing(&self) -> bool {
        let lib_mult = self.libraries.iter().any(HasMultiplexing::has_multiplexing);
        let samp_mult = self
            .samples
            .iter()
            .map(HasMultiplexing::has_multiplexing)
            .dedup()
            .exactly_one()
            .unwrap();
        if lib_mult {
            assert!(samp_mult);
        }
        samp_mult
    }
}

impl CrMultiGraph {
    /// Return the cell multiplexing type in use, if there is one.
    pub fn cell_multiplexing_type(&self) -> Option<CellMultiplexingType> {
        self.samples
            .iter()
            .map(Sample::cell_multiplexing_type)
            .dedup()
            .exactly_one()
            .unwrap()
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct CrMultiParams {
    // TODO: Reference, Target sets, force cells etc
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CrMultiConfig {
    graph: CrMultiGraph,
    params: CrMultiParams,
}

#[derive(
    EnumString,
    Display,
    Debug,
    PartialEq,
    Eq,
    Hash,
    Clone,
    Copy,
    PartialOrd,
    Ord,
    Serialize,
    Deserialize,
)]
pub enum NormalizationMode {
    #[strum(to_string = "mapped")]
    MappedReadsPerCell,
    #[strum(to_string = "none")]
    None,
    #[strum(to_string = "reads_per_umi")]
    ReadsPerUmi,
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
            _ => todo!(),
        }
    }
}

impl TryFrom<LibraryInfo> for crate::aggr::LibraryInfo {
    type Error = anyhow::Error;

    fn try_from(value: LibraryInfo) -> Result<Self> {
        match value {
            LibraryInfo::Aggr(v) => Ok(v),
            _ => todo!(),
        }
    }
}

// TODO: Enable this when we extend to aggr groups
// We need to keep track of per sample metadata
/*
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct AggregationReference {
    pub aggr_id: String,

    /// Samples to include in this aggregation. Each string
    /// must correspond to a declared sample
    pub sample_ids: Vec<String>,

    /// algorithm to use for count normalization,
    pub normalization: NormalizationMode,

    /// Sample metadata column to for creating batch
    /// labels for chemistry batch correction.
    /// None disables CBC.
    pub batch_variable: Option<String>,
}

// pointer to one group of cells in a sample. If fingerprint is null, take all cells from the GemWell.
// It is invalid to point to a GemWell with Multiplexing with a fingerprint of None.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct SampleReference {
    fingerprint: Option<Fingerprint>,
    gem_well: GemWell,
}

impl SampleReference {
    pub fn new(fingerprint: Option<Fingerprint>, gem_well: GemWell) -> SampleReference {
        SampleReference {
            fingerprint,
            gem_well,
        }
    }
}

*/

/// The aligner used to align the reads to the reference.
#[derive(Clone, Copy, Debug, PartialEq, Eq, MartianType, Deserialize, Serialize)]
pub enum AlignerParam {
    /// The 10x probe aligner used for RTL chemistries.
    #[serde(rename = "hurtle")]
    Hurtle,
    /// STAR (Spliced Transcripts Alignment to a Reference)
    #[serde(rename = "star")]
    Star,
}

impl FromStr for AlignerParam {
    type Err = anyhow::Error;

    fn from_str(aligner: &str) -> Result<Self> {
        match aligner.to_ascii_lowercase().as_str() {
            "hurtle" => Ok(AlignerParam::Hurtle),
            "star" => Ok(AlignerParam::Star),
            _ => bail!("Invalid aligner: {aligner}"),
        }
    }
}

impl fmt::Display for AlignerParam {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result<(), std::fmt::Error> {
        let aligner = match self {
            Self::Hurtle => "hurtle",
            Self::Star => "star",
        };
        write!(f, "{aligner}",)
    }
}

impl JsonReport for AlignerParam {
    fn to_json_reporter(&self) -> JsonReporter {
        std::iter::once(("", self)).collect()
    }
}

#[cfg(test)]
mod py_api_tests {
    use super::*;
    use fastq_set::filenames::bcl_processor::SampleIndexSpec;
    use fastq_set::filenames::LaneSpec;
    use insta::assert_json_snapshot;
    use std::collections::BTreeMap;
    use strum::EnumCount;

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
            LibraryFeatures::gex(),
            GemWell(1),
        )?;

        builder.push_library(
            "VDJ_LIB".into(),
            FastqDef::bcl2fastq(
                "/path/to/vdj/fastq/".into(),
                "my_vdj_data".into(),
                LaneSpec::Any,
            ),
            LibraryFeatures::Vdj(VdjChainType::Auto),
            GemWell(1),
        )?;

        builder.push_library(
            "AB_LIB".into(),
            FastqDef::bcl2fastq(
                "/path/to/ab/fastq/".into(),
                "my_ab_data".into(),
                LaneSpec::Any,
            ),
            LibraryFeatures::FeatureBarcodes(FeatureType::Antibody.into()),
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
            LibraryFeatures::gex(),
            GemWell(1),
        )?;
        builder.push_library(
            "GEX_LIB_2".into(),
            FastqDef::bcl2fastq(
                "/path/to/gex_2/fastq/".into(),
                "my_gex_2_data".into(),
                LaneSpec::Any,
            ),
            LibraryFeatures::gex(),
            GemWell(2),
        )?;
        builder.push_library(
            "CMO_LIB_1".into(),
            FastqDef::bcl2fastq(
                "/path/to/cmo_1/fastq/".into(),
                "my_cmo_1_data".into(),
                LaneSpec::Any,
            ),
            LibraryFeatures::FeatureBarcodes(FeatureType::Multiplexing.into()),
            GemWell(1),
        )?;
        builder.push_library(
            "CMO_LIB_2".into(),
            FastqDef::bcl2fastq(
                "/path/to/cmo_2/fastq/".into(),
                "my_cmo_2_data".into(),
                LaneSpec::Any,
            ),
            LibraryFeatures::FeatureBarcodes(FeatureType::Multiplexing.into()),
            GemWell(2),
        )?;

        builder.push_fingerprint(
            "PBMC_donor_A".into(),
            "Human PBMC from donor A".into(),
            Fingerprint::tagged(
                GemWell(1),
                "CMO500".into(),
                Vec::default(),
                CellMultiplexingType::CMO,
            ),
        )?;

        builder.push_fingerprint(
            "PBMC_donor_A".into(),
            "Human PBMC from donor A".into(),
            Fingerprint::tagged(
                GemWell(2),
                "CMO500".into(),
                Vec::default(),
                CellMultiplexingType::CMO,
            ),
        )?;

        builder.push_fingerprint(
            "PBMC_donor_B".into(),
            "Human PBMC from donor B".into(),
            Fingerprint::tagged(
                GemWell(1),
                "CMO501".into(),
                Vec::default(),
                CellMultiplexingType::CMO,
            ),
        )?;

        builder.push_fingerprint(
            "PBMC_donor_B".into(),
            "Human PBMC from donor B".into(),
            Fingerprint::tagged(
                GemWell(2),
                "CMO501".into(),
                Vec::default(),
                CellMultiplexingType::CMO,
            ),
        )?;

        Ok(builder.build())
    }

    #[test]
    fn test_py_api_graph_snapshots() {
        let lib_def_base = LibraryDef {
            physical_library_id: "GEX_1".to_string(),
            gem_well: GemWell(2),
            library_features: LibraryFeatures::gex(),
            fastqs: vec![],
        };

        let mut lib_def_1 = lib_def_base.clone();
        lib_def_1.push_fastq(
            FastqDef::bcl2fastq(
                "/path/to/fastq".into(),
                "my_gex_sample".into(),
                LaneSpec::Any,
            ),
            None,
        );

        let mut lib_def_2 = lib_def_base;
        lib_def_2.push_fastq(
            FastqDef::bcl_processor("/path/to/fastq".into(), SampleIndexSpec::Any, LaneSpec::Any),
            None,
        );

        let feat_type_map: BTreeMap<_, _> = (0..FeatureType::COUNT)
            .map(|i| (FeatureType::from(i as u8).to_string(), i as u8))
            .collect();

        let cell_multiplexing_type_map: BTreeMap<_, _> = (0..CellMultiplexingType::COUNT)
            .map(|i| (CellMultiplexingType::from(i as u8).to_string(), i as u8))
            .collect();

        insta::with_settings!({snapshot_path => "snapshots/py_api"}, {
            assert_json_snapshot!("feat_type_map", feat_type_map);
            assert_json_snapshot!("cell_multiplexing_type_map", cell_multiplexing_type_map);
            assert_json_snapshot!("lib_feat_genes", LibraryFeatures::gex());
            assert_json_snapshot!(
                "lib_feat_vdj_auto",
                LibraryFeatures::Vdj(VdjChainType::Auto)
            );
            assert_json_snapshot!("lib_feat_vdj_tcr", LibraryFeatures::Vdj(VdjChainType::Tcr));
            assert_json_snapshot!("lib_feat_vdj_bcr", LibraryFeatures::Vdj(VdjChainType::Ig));
            assert_json_snapshot!(
                "lib_feat_fb_ab",
                LibraryFeatures::FeatureBarcodes(FeatureType::Antibody.into())
            );
            assert_json_snapshot!(
                "lib_feat_fb_all",
                LibraryFeatures::FeatureBarcodes(FeatureType::feature_barcodes().collect_vec().into())
            );

            assert_json_snapshot!("lib_def_1", lib_def_1);
            assert_json_snapshot!("lib_def_2", lib_def_2);

            assert_json_snapshot!("fingerprint_untagged", Fingerprint::untagged(GemWell(1)));
            assert_json_snapshot!("fingerprint_tagged", Fingerprint::tagged(GemWell(2), "CMO500".into(), Vec::default(), CellMultiplexingType::CMO));
            assert_json_snapshot!("fingerprint_rtl_tagged", Fingerprint::tagged(GemWell(2), "CMO500".into(), Vec::default(), CellMultiplexingType::RTL));

            assert_json_snapshot!("sample_untagged", Sample {
                sample_id: "100".into(),
                description: "My sample".into(),
                fingerprints: vec![Fingerprint::untagged(GemWell(1))]
            });
            assert_json_snapshot!("sample_tagged", Sample {
                sample_id: "101".into(),
                description: "My tagged sample".into(),
                fingerprints: vec![
                    Fingerprint::tagged(GemWell(1), "CMO500".into(), Vec::default(), CellMultiplexingType::CMO),
                    Fingerprint::tagged(GemWell(1), "CMO501".into(), Vec::default(), CellMultiplexingType::CMO),
                    Fingerprint::tagged(GemWell(2), "CMO500".into(), Vec::default(), CellMultiplexingType::CMO)
                ]
            });
            assert_json_snapshot!("multi_graph_untagged", untagged_single_gw_graph().unwrap());
            assert_json_snapshot!("multi_graph_tagged", tagged_multi_gw_graph().unwrap());
        });
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use strum::EnumCount;

    #[test]
    fn test_legacy_conversion() {
        for legacy_type in LegacyLibraryType::iter() {
            if legacy_type == LegacyLibraryType::ATAC {
                continue;
            }
            let lib_feat = LibraryFeatures::from(legacy_type);
            assert_eq!(legacy_type, lib_feat.legacy_library_type());
        }
    }

    #[test]
    #[should_panic]
    fn test_invalid_legacy_conversion() {
        let lib_feat = LibraryFeatures::FeatureBarcodes(
            vec![FeatureType::Antibody, FeatureType::Antigen].into(),
        );
        let _ = lib_feat.legacy_library_type();
    }

    #[test]
    fn test_iter_order() {
        for (i, feat) in FeatureType::iter().enumerate() {
            assert_eq!(i as u8, feat as u8);
        }
    }

    #[test]
    fn test_u8_roundtrip() {
        for i in 0..FeatureType::COUNT {
            let f = FeatureType::from(i as u8);
            let j: u8 = f.into();
            assert_eq!(i as u8, j);
        }
    }

    #[test]
    fn test_bit_encoded() {
        let mut encoded = BitEncoded::new();
        encoded.push(FeatureType::Antibody);
        assert_eq!(encoded.inspect_bits(), 1u8);
        encoded.push(FeatureType::Multiplexing);
        assert_eq!(encoded.inspect_bits(), 9u8);
        assert_eq!(
            encoded.iter().collect_vec(),
            vec![FeatureType::Antibody, FeatureType::Multiplexing]
        );
    }

    #[test]
    fn test_multiplexing_aliases() {
        let feat: FeatureType = serde_json::from_str(r#""FEATURETEST""#).unwrap();
        assert_eq!(feat, FeatureType::Multiplexing);
        let feat: FeatureType = serde_json::from_str(r#""Multiplexing Capture""#).unwrap();
        assert_eq!(feat, FeatureType::Multiplexing);
        let feat: FeatureType = serde_json::from_str(r#""Multiplexing Tag Capture""#).unwrap();
        assert_eq!(feat, FeatureType::Multiplexing);
    }
}
