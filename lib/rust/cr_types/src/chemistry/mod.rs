//! Specification of analysis chemistries.

mod chemistry_defs;

use crate::types::ReqStrand;
use crate::LibraryType;
use anyhow::{bail, Result};
use barcode::whitelist::{find_slide_design, BarcodeId};
use barcode::{
    BarcodeConstruct, BcSegSeq, GelBeadAndProbeConstruct, Segments, WhitelistSource, WhitelistSpec,
};
use chemistry_defs::get_chemistry_def;
pub use chemistry_defs::{known_chemistry_defs, normalize_chemistry_def};
use fastq_set::read_pair::{RpRange, WhichRead};
use fastq_set::WhichEnd;
use itertools::{chain, Itertools};
use martian::{AsMartianPrimaryType, MartianPrimaryType};
use martian_derive::{MartianStruct, MartianType};
use metric::{JsonReport, JsonReporter, TxHashMap, TxHashSet};
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::collections::HashMap;
use std::fmt::{Display, Formatter};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::iter::FromIterator;
use std::path::PathBuf;
use std::str::FromStr;
use strum_macros::{Display, EnumIter, EnumString};

// "Auto" chemistries which are not fully specified and need to be further
// refined into one of
#[derive(
    Copy,
    Clone,
    Debug,
    Eq,
    Hash,
    PartialEq,
    Ord,
    PartialOrd,
    Display,
    EnumIter,
    EnumString,
    Serialize,
    Deserialize,
)]
pub enum AutoChemistryName {
    #[serde(rename = "auto")]
    #[strum(serialize = "auto")]
    Count,
    #[serde(rename = "threeprime", alias = "SC3P_auto")]
    #[strum(to_string = "threeprime", serialize = "SC3P_auto")]
    ThreePrime,
    #[serde(rename = "fiveprime", alias = "SC5P_auto")]
    #[strum(to_string = "fiveprime", serialize = "SC5P_auto")]
    FivePrime,
    #[serde(rename = "SCVDJ_auto")]
    #[strum(serialize = "SCVDJ_auto")]
    Vdj,
}

const THREE_PRIME_AUTO_CHEMS: [ChemistryName; 7] = [
    ChemistryName::ThreePrimeV1,
    ChemistryName::ThreePrimeV2,
    ChemistryName::ThreePrimeV3,
    // LT is no longer supported, but we retain it here so that we can issue
    // an accurate error message in DETECT_CHEMISTRY if we detect an LT kit.
    ChemistryName::ThreePrimeV3LT,
    ChemistryName::ThreePrimeV3HT,
    ChemistryName::ThreePrimeV4,
    // ARC is not supported by count, but is included to improve debugging
    // of reagent mix-ups.
    ChemistryName::ArcV1,
];

const FIVE_PRIME_AUTO_CHEMS: [ChemistryName; 4] = [
    ChemistryName::FivePrimeR2,
    ChemistryName::FivePrimePE,
    ChemistryName::FivePrimeR2V3,
    ChemistryName::FivePrimePEV3,
];

const VDJ_AUTO_CHEMS: [ChemistryName; 6] = [
    ChemistryName::VdjPE,
    ChemistryName::VdjR2,
    ChemistryName::VdjPEV3,
    ChemistryName::VdjR2V3,
    ChemistryName::VdjR2FRP,
    ChemistryName::ThreePrimeV3,
];

impl AutoChemistryName {
    /// Return the allowed chemistries for this auto chemistry mode.
    /// Which RTL chemistries are allowed is determined somewhat dynamically,
    /// thus we require passing some configuration parameters.
    pub fn allowed_chemistries(
        &self,
        is_pd: bool,
        is_rtl_multiplexed: bool,
        using_rtl_uncollapsed_probe_barcodes: bool,
    ) -> TxHashSet<ChemistryName> {
        match self {
            Self::ThreePrime => THREE_PRIME_AUTO_CHEMS.iter().copied().collect(),
            Self::FivePrime => FIVE_PRIME_AUTO_CHEMS.iter().copied().collect(),
            Self::Vdj => VDJ_AUTO_CHEMS.iter().copied().collect(),
            Self::Count => {
                let frp_chems = if is_rtl_multiplexed {
                    if is_pd {
                        // PD
                        if using_rtl_uncollapsed_probe_barcodes {
                            // Uncollapsed probe barcode chemistries
                            [
                                ChemistryName::MFRP_uncollapsed,
                                ChemistryName::MFRP_R1_48_uncollapsed,
                            ]
                            .as_slice()
                        } else {
                            // Collapsed probe barcode chemistries
                            [
                                ChemistryName::MFRP_RNA,
                                ChemistryName::MFRP_Ab,
                                ChemistryName::MFRP_47,
                                ChemistryName::MFRP_RNA_R1,
                                ChemistryName::MFRP_Ab_R1,
                                ChemistryName::MFRP_Ab_R2pos50,
                                ChemistryName::MFRP_CRISPR,
                            ]
                            .as_slice()
                        }
                    } else {
                        // CS
                        [
                            ChemistryName::MFRP_RNA,
                            ChemistryName::MFRP_Ab,
                            ChemistryName::MFRP_RNA_R1,
                            ChemistryName::MFRP_Ab_R1,
                            ChemistryName::MFRP_Ab_R2pos50,
                            ChemistryName::MFRP_CRISPR,
                        ]
                        .as_slice()
                    }
                } else {
                    // Without a [samples] section
                    [ChemistryName::SFRP].as_slice()
                };
                chain!(&THREE_PRIME_AUTO_CHEMS, &FIVE_PRIME_AUTO_CHEMS, frp_chems)
                    .copied()
                    .collect()
            }
        }
    }
}

/// Chemistry names supported by Cellranger excluding auto chemistries
/// - PE stands for "Paired End"
#[derive(
    Copy,
    Clone,
    Debug,
    Eq,
    PartialEq,
    Hash,
    Ord,
    PartialOrd,
    Display,
    EnumIter,
    EnumString,
    Serialize,
    Deserialize,
    MartianType,
)]
#[allow(non_camel_case_types)]
pub enum ChemistryName {
    #[serde(rename = "custom")]
    #[strum(serialize = "custom")]
    Custom,

    /// Singleplex fixed RNA profiling
    SFRP,

    /// Multiplex fixed RNA profiling
    #[serde(rename = "MFRP-RNA")]
    #[strum(serialize = "MFRP-RNA")]
    MFRP_RNA,

    /// Multiplex fixed RNA profiling, antibody capture
    #[serde(rename = "MFRP-Ab")]
    #[strum(serialize = "MFRP-Ab")]
    MFRP_Ab,

    /// Multiplex fixed RNA profiling with 47 probe barcodes
    #[serde(rename = "MFRP-47")]
    #[strum(serialize = "MFRP-47")]
    MFRP_47,

    /// Multiplex fixed RNA profiling without collapsing base-balanced barcodes
    #[serde(rename = "MFRP-uncollapsed")]
    #[strum(serialize = "MFRP-uncollapsed")]
    MFRP_uncollapsed,

    /// Multiplex fixed RNA profiling (probeBC on R1)
    #[serde(rename = "MFRP-RNA-R1")]
    #[strum(serialize = "MFRP-RNA-R1")]
    MFRP_RNA_R1,

    /// Multiplex fixed RNA profiling, antibody capture (probeBC on R1)
    #[serde(rename = "MFRP-Ab-R1")]
    #[strum(serialize = "MFRP-Ab-R1")]
    MFRP_Ab_R1,

    /// Multiplex fixed RNA profiling (probeBC on R1) with 192 non-base-balanced probe barcodes
    #[serde(rename = "MFRP-R1-48-uncollapsed")]
    #[strum(serialize = "MFRP-R1-48-uncollapsed")]
    MFRP_R1_48_uncollapsed,

    /// Multiplex fixed RNA Profiling (probe barcode at R2:50)
    #[serde(rename = "MFRP-Ab-R2pos50")]
    #[strum(serialize = "MFRP-Ab-R2pos50")]
    MFRP_Ab_R2pos50,

    /// Multiplex fixed RNA Profiling (CRISPR)
    #[serde(rename = "MFRP-CRISPR")]
    #[strum(serialize = "MFRP-CRISPR")]
    MFRP_CRISPR,

    #[serde(rename = "SC-FB")]
    #[strum(serialize = "SC-FB")]
    FeatureBarcodingOnly,

    #[serde(rename = "SCVDJ")]
    #[strum(serialize = "SCVDJ")]
    VdjPE,
    #[serde(rename = "SCVDJ-R2")]
    #[strum(serialize = "SCVDJ-R2")]
    VdjR2,
    #[serde(rename = "SCVDJ-v3")]
    #[strum(serialize = "SCVDJ-v3")]
    VdjPEV3,
    #[serde(rename = "SCVDJ-Splint-R2-FRP")]
    #[strum(serialize = "SCVDJ-Splint-R2-FRP")]
    VdjR2FRP,
    #[serde(rename = "SCVDJ-R2-v3")]
    #[strum(serialize = "SCVDJ-R2-v3")]
    VdjR2V3,
    #[serde(rename = "SCVDJ-R1")]
    #[strum(serialize = "SCVDJ-R1")]
    VdjR1, // Deprecated now, does not have a V2 version

    #[serde(rename = "SC5P-R1")]
    #[strum(serialize = "SC5P-R1")]
    FivePrimeR1,
    #[serde(rename = "SC5P-R2")]
    #[strum(serialize = "SC5P-R2")]
    FivePrimeR2,
    #[serde(rename = "SC5P-R2-OH")]
    #[strum(serialize = "SC5P-R2-OH")]
    FivePrimeR2OH,
    #[serde(rename = "SC5PHT")]
    #[strum(serialize = "SC5PHT")]
    FivePrimeHT,
    #[serde(rename = "SC5P-PE")]
    #[strum(serialize = "SC5P-PE")]
    FivePrimePE,

    #[serde(rename = "SC5P-R1-v3")]
    #[strum(serialize = "SC5P-R1-v3")]
    FivePrimeR1V3,
    #[serde(rename = "SC5P-R2-v3")]
    #[strum(serialize = "SC5P-R2-v3")]
    FivePrimeR2V3,
    #[serde(rename = "SC5P-R2-OH-v3")]
    #[strum(serialize = "SC5P-R2-OH-v3")]
    FivePrimeR2OHV3,
    #[serde(rename = "SC5P-PE-v3")]
    #[strum(serialize = "SC5P-PE-v3")]
    FivePrimePEV3,

    #[serde(rename = "SC3Pv1")]
    #[strum(serialize = "SC3Pv1")]
    ThreePrimeV1,
    #[serde(rename = "SC3Pv2")]
    #[strum(serialize = "SC3Pv2")]
    ThreePrimeV2,
    #[serde(rename = "SC3Pv3")]
    #[strum(serialize = "SC3Pv3")]
    ThreePrimeV3,
    #[serde(rename = "SC3Pv3-OH")]
    #[strum(serialize = "SC3Pv3-OH")]
    ThreePrimeV3OH,
    #[serde(rename = "SC3Pv3LT")]
    #[strum(serialize = "SC3Pv3LT")]
    ThreePrimeV3LT, // deprecated, to be removed
    #[serde(rename = "SC3Pv3HT")]
    #[strum(serialize = "SC3Pv3HT")]
    ThreePrimeV3HT,
    #[serde(rename = "SC3Pv4")]
    #[strum(serialize = "SC3Pv4")]
    ThreePrimeV4,
    #[serde(rename = "SC3Pv4-OH")]
    #[strum(serialize = "SC3Pv4-OH")]
    ThreePrimeV4OH,

    #[serde(rename = "SPATIAL3Pv1")]
    #[strum(serialize = "SPATIAL3Pv1")]
    SpatialThreePrimeV1,
    #[serde(rename = "SPATIAL3Pv2")]
    #[strum(serialize = "SPATIAL3Pv2")]
    SpatialThreePrimeV2,
    #[serde(rename = "SPATIAL3Pv3")]
    #[strum(serialize = "SPATIAL3Pv3")]
    SpatialThreePrimeV3,
    #[serde(rename = "SPATIAL3Pv4")]
    #[strum(serialize = "SPATIAL3Pv4")]
    SpatialThreePrimeV4,
    #[serde(rename = "SPATIAL3Pv5")]
    #[strum(serialize = "SPATIAL3Pv5")]
    SpatialThreePrimeV5,

    #[serde(rename = "SPATIAL-HD-v1")]
    #[strum(serialize = "SPATIAL-HD-v1")]
    SpatialHdV1,

    #[serde(rename = "ARC-v1")]
    #[strum(serialize = "ARC-v1")]
    ArcV1,

    #[serde(rename = "ATAC-v1")]
    #[strum(serialize = "ATAC-v1")]
    AtacV1,
}

impl ChemistryName {
    /// Return Some(true) if an RTL chemistry, Some(false) if not an RTL chemistry, and None if unknown.
    /// Custom and spatial chemistries may be either 3'GEX or RTL.
    pub fn is_rtl(&self) -> Option<bool> {
        #[allow(clippy::enum_glob_use)]
        use ChemistryName::*;
        match self {
            // Custom and spatial chemistries may be either 3'GEX or RTL.
            Custom | SpatialThreePrimeV1 | SpatialThreePrimeV2 | SpatialThreePrimeV3
            | SpatialThreePrimeV4 | SpatialThreePrimeV5 => None,

            // Visium HD v1 is RTL
            SpatialHdV1 => Some(true),

            // FRP chemistries are RTL.
            SFRP
            | MFRP_RNA
            | MFRP_Ab
            | MFRP_47
            | MFRP_uncollapsed
            | MFRP_RNA_R1
            | MFRP_Ab_R1
            | MFRP_R1_48_uncollapsed
            | MFRP_Ab_R2pos50
            | MFRP_CRISPR => Some(true),

            // All other chemistries are not RTL.
            // These are explicitly listed for future exhaustiveness checks.
            FeatureBarcodingOnly | VdjPE | VdjR2 | VdjPEV3 | VdjR2FRP | VdjR2V3 | VdjR1
            | FivePrimeR1 | FivePrimeR2 | FivePrimeR2OH | FivePrimeHT | FivePrimePE
            | FivePrimeR1V3 | FivePrimeR2V3 | FivePrimeR2OHV3 | FivePrimePEV3 | ThreePrimeV1
            | ThreePrimeV2 | ThreePrimeV3 | ThreePrimeV3OH | ThreePrimeV3LT | ThreePrimeV3HT
            | ThreePrimeV4 | ThreePrimeV4OH | ArcV1 | AtacV1 => Some(false),
        }
    }

    /// Return true if a multiplexed RTL chemistry.
    pub fn is_rtl_multiplexed(&self) -> bool {
        self.is_rtl().unwrap_or(false) && self != &ChemistryName::SFRP
    }

    /// Return true if this is an MFRP chemistry.
    pub fn is_mfrp(&self) -> bool {
        #[allow(clippy::enum_glob_use)]
        use ChemistryName::*;
        match self {
            MFRP_RNA
            | MFRP_Ab
            | MFRP_47
            | MFRP_uncollapsed
            | MFRP_RNA_R1
            | MFRP_Ab_R1
            | MFRP_R1_48_uncollapsed
            | MFRP_Ab_R2pos50
            | MFRP_CRISPR => true,
            Custom | SFRP | FeatureBarcodingOnly | VdjPE | VdjR2 | VdjPEV3 | VdjR2FRP | VdjR2V3
            | VdjR1 | FivePrimeR1 | FivePrimeR2 | FivePrimeR2OH | FivePrimeHT | FivePrimePE
            | FivePrimeR1V3 | FivePrimeR2V3 | FivePrimeR2OHV3 | FivePrimePEV3 | ThreePrimeV1
            | ThreePrimeV2 | ThreePrimeV3 | ThreePrimeV3OH | ThreePrimeV3LT | ThreePrimeV3HT
            | ThreePrimeV4 | ThreePrimeV4OH | SpatialThreePrimeV1 | SpatialThreePrimeV2
            | SpatialThreePrimeV3 | SpatialThreePrimeV4 | SpatialThreePrimeV5 | SpatialHdV1
            | ArcV1 | AtacV1 => false,
        }
    }

    /// Return true if this chemistry is compatible with the provided library type.
    /// NOTE: this method is currently only used for Flex, and the full implementation
    /// would be verbose. This can be expanded in the future if other products
    /// need to answer this question.
    pub fn compatible_with_library_type(&self, library_type: LibraryType) -> bool {
        assert!(self.is_mfrp() || *self == SFRP);
        #[allow(clippy::enum_glob_use)]
        use ChemistryName::*;
        match library_type {
            LibraryType::Gex => matches!(
                self,
                SFRP | MFRP_RNA | MFRP_47 | MFRP_uncollapsed | MFRP_RNA_R1 | MFRP_R1_48_uncollapsed
            ),
            LibraryType::Antibody => matches!(
                self,
                SFRP | MFRP_Ab
                    | MFRP_47
                    | MFRP_uncollapsed
                    | MFRP_Ab_R1
                    | MFRP_R1_48_uncollapsed
                    | MFRP_Ab_R2pos50
            ),
            LibraryType::Crispr => matches!(
                self,
                SFRP | MFRP_47 | MFRP_uncollapsed | MFRP_R1_48_uncollapsed | MFRP_CRISPR
            ),
            _ => panic!("the library type {library_type} is not supported by Flex"),
        }
    }

    /// Return true if a spatial chemistry.
    pub fn is_spatial(&self) -> bool {
        #[allow(clippy::enum_glob_use)]
        use ChemistryName::*;
        matches!(
            self,
            SpatialThreePrimeV1
                | SpatialThreePrimeV2
                | SpatialThreePrimeV3
                | SpatialThreePrimeV4
                | SpatialThreePrimeV5
                | SpatialHdV1
        )
    }

    /// Return true if a spatial chemistry that supports CytAssist
    pub fn is_cytassist_compatible(&self) -> bool {
        #[allow(clippy::enum_glob_use)]
        use ChemistryName::*;
        matches!(
            self,
            SpatialThreePrimeV4 | SpatialThreePrimeV5 | SpatialHdV1
        )
    }

    // Returns true if a spatial HD chemistry is an RTL chemistry
    pub fn is_spatial_hd_rtl(&self) -> bool {
        #[allow(clippy::enum_glob_use)]
        use ChemistryName::*;
        matches!(self, SpatialHdV1)
    }

    /// Return true if a spatial chemistry that supports feature barcoding
    pub fn is_spatial_fb(&self) -> bool {
        #[allow(clippy::enum_glob_use)]
        use ChemistryName::*;
        matches!(
            self,
            SpatialThreePrimeV3 | SpatialThreePrimeV4 | SpatialThreePrimeV5
        )
    }

    /// Return overhang version of chemistry.
    pub fn get_overhang_version(&self) -> Result<ChemistryName> {
        #[allow(clippy::enum_glob_use)]
        use ChemistryName::*;
        match self {
            ThreePrimeV3 => Ok(ThreePrimeV3OH),
            ThreePrimeV4 => Ok(ThreePrimeV4OH),
            FivePrimeR2 => Ok(FivePrimeR2OH),
            FivePrimeR2V3 => Ok(FivePrimeR2OHV3),

            _ => bail!("Overhang chemistry undefined for {self}"),
        }
    }

    /// Return true if a VDJ chemistry.
    pub fn is_vdj(&self) -> bool {
        #[allow(clippy::enum_glob_use)]
        use ChemistryName::*;
        matches!(
            self,
            VdjPE | VdjR1 | VdjR2 | VdjPEV3 | VdjR2V3 | VdjR2FRP | ThreePrimeV3
        )
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(untagged)] // Ensures that this (de)serializes from/to a string
pub enum AutoOrRefinedChemistry {
    Auto(AutoChemistryName),
    Refined(ChemistryName),
}

impl AsMartianPrimaryType for AutoOrRefinedChemistry {
    fn as_martian_primary_type() -> MartianPrimaryType {
        MartianPrimaryType::Str
    }
}

impl FromStr for AutoOrRefinedChemistry {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<AutoOrRefinedChemistry> {
        Ok(if let Ok(auto_chem) = AutoChemistryName::try_from(s) {
            AutoOrRefinedChemistry::Auto(auto_chem)
        } else if let Ok(chem) = ChemistryName::try_from(s) {
            AutoOrRefinedChemistry::Refined(chem)
        } else {
            bail!("Could not parse {s} as AutoOrRefinedChemistry")
        })
    }
}

impl From<AutoChemistryName> for AutoOrRefinedChemistry {
    fn from(auto_chem: AutoChemistryName) -> AutoOrRefinedChemistry {
        AutoOrRefinedChemistry::Auto(auto_chem)
    }
}

impl From<ChemistryName> for AutoOrRefinedChemistry {
    fn from(chem: ChemistryName) -> AutoOrRefinedChemistry {
        AutoOrRefinedChemistry::Refined(chem)
    }
}

impl Display for AutoOrRefinedChemistry {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result<(), std::fmt::Error> {
        match self {
            AutoOrRefinedChemistry::Auto(x) => write!(f, "{x}"),
            AutoOrRefinedChemistry::Refined(x) => write!(f, "{x}"),
        }
    }
}

impl AutoOrRefinedChemistry {
    #[allow(non_upper_case_globals)]
    /// Alias for custom chemistry.
    pub const Custom: Self = Self::Refined(ChemistryName::Custom);

    /// Return whether the chemistry is RTL.
    /// Return None if the chemistry is auto.
    pub fn is_rtl(&self) -> Option<bool> {
        match self {
            Self::Auto(_) => None,
            Self::Refined(name) => name.is_rtl(),
        }
    }

    /// Return true if this chemistry is specifically an MFRP variant.
    pub fn is_mfrp(&self) -> bool {
        match self {
            Self::Auto(_) => false,
            Self::Refined(name) => name.is_mfrp(),
        }
    }

    /// Return true if this is custom chemistry.
    pub fn is_custom(&self) -> bool {
        *self == Self::Custom
    }

    /// Return Some if the chemistry is refined, None if auto.
    pub fn refined(&self) -> Option<ChemistryName> {
        if let Self::Refined(name) = *self {
            Some(name)
        } else {
            None
        }
    }

    /// Return true if this is an auto chemistry.
    pub fn is_auto(&self) -> bool {
        self.auto().is_some()
    }

    /// Return Some if the chemistry is auto, None if refined.
    pub fn auto(&self) -> Option<AutoChemistryName> {
        if let Self::Auto(mode) = *self {
            Some(mode)
        } else {
            None
        }
    }
}

impl ChemistryName {
    // figure out the headers for bamtofastq
    pub fn bamtofastq_headers(&self) -> &'static [&'static str] {
        #[allow(clippy::enum_glob_use)]
        use ChemistryName::*;

        match self {
            // standard 'modern' chemistries
            ThreePrimeV2 | ThreePrimeV3 | ThreePrimeV3OH | ThreePrimeV3LT | ThreePrimeV3HT
            | ThreePrimeV4 | ThreePrimeV4OH | SpatialThreePrimeV1 | SpatialThreePrimeV2
            | SpatialThreePrimeV3 | SpatialThreePrimeV4 | SpatialThreePrimeV5 | FivePrimeR2
            | FivePrimeR2OH | FivePrimeHT | FivePrimeR2V3 | FivePrimeR2OHV3
            | FeatureBarcodingOnly | SFRP | ArcV1 => &[
                "10x_bam_to_fastq:R1(CR:CY,UR:UY)",
                "10x_bam_to_fastq:R2(SEQ:QUAL)",
            ],

            // Multiplexed fixed RNA profiling
            MFRP_RNA
            | MFRP_Ab
            | MFRP_47
            | MFRP_uncollapsed
            | MFRP_RNA_R1
            | MFRP_Ab_R1
            | MFRP_R1_48_uncollapsed
            | MFRP_Ab_R2pos50
            | MFRP_CRISPR => &[
                "10x_bam_to_fastq:R1(GR:GY,UR:UY,1R:1Y)",
                "10x_bam_to_fastq:R2(SEQ:QUAL,2R:2Y)",
            ],

            // old school 3' v1 chemistry
            ThreePrimeV1 => &[
                "10x_bam_to_fastq:I1(CR:CY)",
                "10x_bam_to_fastq:R1(SEQ:QUAL)",
                "10x_bam_to_fastq:R2(UR:UY)",
            ],

            // VDJ chemistry should not be coming through here - may need changing if
            // we decide to make bamtofastq compatible BAMs from VDJ in the future.
            VdjR1 | VdjR2 | VdjPE | VdjR2V3 | VdjPEV3 | VdjR2FRP => {
                panic!("VDJ don't yet make a normal BAM file")
            }

            // this chemistry supported for customers & bamtofastq doesn't work, so don't emit any
            // bamtofastq headers
            FivePrimeR1 | FivePrimeR1V3 => &[],

            // 5' PE chem
            FivePrimePE | FivePrimePEV3 => &[
                "10x_bam_to_fastq:R1(CR:CY,UR:UY,SEQ:QUAL)",
                "10x_bam_to_fastq:R2(SEQ:QUAL)",
            ],

            // for custom chemistry, don't bother with bamtofastq headers
            Custom => &[],

            // for ATAC BAMs generated by ATAC-v1 or ARC-v1
            AtacV1 => &[
                "10x_bam_to_fastq:R1(SEQ:QUAL,TR:TQ)",
                "10x_bam_to_fastq:R2(SEQ:QUAL,TR:TQ)",
                "10x_bam_to_fastq:I1(BC:QT)",
                "10x_bam_to_fastq:I2(CR:CY)",
                "10x_bam_to_fastq_seqnames:R1,R3,I1,R2",
            ],
            SpatialHdV1 => &[
                "10x_bam_to_fastq:R1(1R:1Y)",
                "10x_bam_to_fastq:R2(SEQ:QUAL,2R:2Y)",
            ],
        }
    }
}

#[derive(Copy, Serialize, Deserialize, Clone, PartialEq, Eq, Debug, MartianStruct)]
pub struct RnaReadComponent {
    #[mro_type = "string"]
    pub read_type: WhichRead,
    pub offset: usize,
    pub length: Option<usize>,
    pub min_length: Option<usize>,
}

impl From<RnaReadComponent> for RpRange {
    fn from(v: RnaReadComponent) -> RpRange {
        RpRange::new(v.read_type, v.offset, v.length)
    }
}

// TODO: Should this be in the `umi` crate?
#[derive(Copy, Serialize, Deserialize, Clone, PartialEq, Eq, Debug, MartianType)]
#[serde(rename_all = "snake_case")]
pub enum UmiTranslation {
    /// Translate the UMI part to a single base. Used in multi barcode chemistries
    /// to map the 4bp splints between barcode parts which have a diversity of 4 to a
    /// single nucleotide
    SingleBase,
}

// TODO: Should this be in the `umi` crate?
#[derive(Serialize, Deserialize, Clone, PartialEq, Eq, Debug, MartianStruct)]
pub struct UmiWhitelistSpec {
    pub slide: String,
    #[mro_type = "string"]
    pub part: slide_design::OligoPart,
    pub translation: UmiTranslation,
}

impl UmiWhitelistSpec {
    pub fn sequences(&self) -> Vec<String> {
        let slide_path = find_slide_design(&self.slide).expect("Failed to find slide design file");
        slide_design::load_oligos(&slide_path, self.part).unwrap()
    }
}

#[derive(Serialize, Deserialize, Clone, PartialEq, Eq, Debug, MartianStruct)]
pub struct UmiReadComponent {
    #[mro_type = "string"]
    pub read_type: WhichRead,
    pub offset: usize,
    /// The length of the UMI. At most this number of bases will be extracted for use as a UMI.
    pub length: usize,
    /// If a shorter UMI can be used, add it here. None indicates that the full
    /// `length` is required.
    pub min_length: Option<usize>,
    pub whitelist: Option<UmiWhitelistSpec>,
}

#[derive(
    Serialize, Deserialize, Clone, Copy, PartialOrd, Ord, PartialEq, Eq, Hash, Debug, MartianType,
)]
#[serde(rename_all = "snake_case")]
pub(crate) enum BarcodeKind {
    GelBead,
    LeftProbe,
    RightProbe,
    SpotSegment,
    // To enable OH multiplexing
    Overhang,
}

/// One component of a segmented barcode.
#[derive(Serialize, Deserialize, Clone, PartialEq, Eq, Hash, Debug, MartianStruct)]
pub struct BarcodeReadComponent {
    #[mro_type = "string"]
    pub(crate) read_type: WhichRead,
    pub(crate) kind: BarcodeKind,
    pub(crate) offset: usize,
    pub(crate) length: usize,
    pub(crate) whitelist: WhitelistSpec,
}

impl BarcodeReadComponent {
    /// Return the whitelist.
    pub fn whitelist(&self) -> &WhitelistSpec {
        &self.whitelist
    }

    /// Return the length of this barcode segment.
    pub fn length(&self) -> usize {
        self.length
    }

    /// Return the offset of this barcode segment.
    pub fn offset(&self) -> usize {
        self.offset
    }

    pub fn rprange(&self) -> RpRange {
        RpRange::new(self.read_type, self.offset, Some(self.length))
    }

    /// Return whether this barcode componenet is a gel bead barcode.
    pub fn is_gel_bead(&self) -> bool {
        self.kind == BarcodeKind::GelBead
    }

    /// Return whether this barcode componenet is a probe barcode.
    pub fn is_probe(&self) -> bool {
        matches!(self.kind, BarcodeKind::LeftProbe | BarcodeKind::RightProbe)
    }

    /// Replace the whitelist with a dynamically generated translation.
    /// Use the provided other whitelist for the translated sequenced.
    /// The whitelist will be written to the provided path.
    pub fn translate_whitelist_with_id_map(
        &mut self,
        id_map: &TxHashMap<BarcodeId, BarcodeId>,
        translate_to: &WhitelistSource,
        path: PathBuf,
    ) -> Result<()> {
        let source = self.whitelist.as_source(true)?;
        let mut file = BufWriter::new(File::create(&path)?);
        for row in source.create_translation_from_id_map(id_map, translate_to)? {
            writeln!(file, "{}\t{}\t{}", row.0, row.1, row.2)?;
        }
        file.flush()?;
        self.whitelist = WhitelistSpec::DynamicTranslation {
            translation_whitelist_path: path,
        };
        Ok(())
    }

    /// For a translated whitelist build a sequence to id map
    pub fn build_seq_to_id_map(&self) -> Result<TxHashMap<BcSegSeq, BarcodeId>> {
        self.whitelist().as_source(true)?.as_translation_seq_to_id()
    }
}

fn bc_vec_to_construct(comps: &[BarcodeReadComponent]) -> BarcodeConstruct<&BarcodeReadComponent> {
    use BarcodeKind::{GelBead, LeftProbe, Overhang, RightProbe, SpotSegment};
    match comps {
        [] => panic!("No barcode segments"),
        [gel_bead] => {
            assert_eq!(gel_bead.kind, GelBead);
            BarcodeConstruct::GelBeadOnly(gel_bead)
        }
        [first, other] => match (first.kind, other.kind) {
            (GelBead, LeftProbe) | (GelBead, RightProbe) => {
                BarcodeConstruct::GelBeadAndProbe(GelBeadAndProbeConstruct {
                    gel_bead: first,
                    probe: other,
                })
            }
            (GelBead, Overhang) => BarcodeConstruct::GelBeadOnly(first),
            (SpotSegment, SpotSegment) => BarcodeConstruct::Segmented(Segments::from_iter(comps)),
            _ => unimplemented!("Combination of two other kinds segments are not supported"),
        },
        [_, _, _] => unimplemented!("Barcodes with three segments are not supported"),
        [_, _, _, _] => {
            assert!(comps.iter().all(|x| x.kind == SpotSegment));
            BarcodeConstruct::Segmented(Segments::from_iter(comps))
        }
        [_, _, _, _, _, ..] => {
            unimplemented!("Barcodes with five or more segments are not supported")
        }
    }
}

fn bc_vec_to_overhang_read_component(
    comps: &[BarcodeReadComponent],
) -> Option<&BarcodeReadComponent> {
    use BarcodeKind::{GelBead, LeftProbe, Overhang, RightProbe, SpotSegment};
    match comps {
        [] => panic!("No barcode segments"),
        [_] => None,
        [first, other] => match (first.kind, other.kind) {
            (GelBead, LeftProbe) | (GelBead, RightProbe) => None,
            (GelBead, Overhang) => Some(other),
            (SpotSegment, SpotSegment) => None,
            _ => unimplemented!("Combination of two other kinds segments are not supported"),
        },
        [_, _, _] => unimplemented!("Barcodes with three segments are not supported"),
        [_, _, _, _] => None,
        [_, _, _, _, _, ..] => {
            unimplemented!("Barcodes with five or more segments are not supported")
        }
    }
}

/// Methods by which barcodes are extracted from the read
#[derive(Serialize, Deserialize, Clone, Default, PartialEq, Eq, Debug)]
#[serde(rename_all = "snake_case")]
#[serde(tag = "method", content = "params")]
pub enum BarcodeExtraction {
    /// Extract each part of the barcode independently using the offset
    /// and length specified in the `BarcodeReadComponent`. This works
    /// for all chemistries where barcodes are present at a fixed position in read1
    #[default]
    Independent,
    /// A two part barcode chemistry, where the two parts are adjacent to
    /// each other. The position of bc1 is flexible due to variable length
    /// regions upstream of bc1 and the possibility of indels. Instead of an
    /// independent extraction, extract both bc1 and bc2 simultaneously.
    /// Barcode correction would allow both indels and mismatches.
    JointBc1Bc2 {
        /// The minimum possible starting position of barcode part 1
        min_offset: usize,
        /// The maximum possible starting position of barcode part 1 (inclusive)
        max_offset: usize,
    },
}

impl AsMartianPrimaryType for BarcodeExtraction {
    fn as_martian_primary_type() -> MartianPrimaryType {
        MartianPrimaryType::Map
    }
}

/// Define a chemistry supported by our RNA products.
///
/// A chemistry tells you where & how to look for various read components
/// (cell barcode, cDNA sequence, UMI, sample index) from the FASTQ
/// cluster data. It is convenient to specify the chemistry as a JSON
/// file and let `ChemistryDef` deserialize it.
///
/// As an example, consider the Single Cell V(D)J read layout below:
/// ![Plot](../../../../doc-media/fastq/scvdj_chem.png)
/// The corresponding chemistry definition would be:
///
/// - Barcode is present in the first 16 bases of Read 1.
/// This translates to a `read_type` "R1", `read_offset` 0
/// and `read_length` 16. Valid options for `read_type` are
/// "R1", "R2", "I1", "I2"
/// ``` text
/// {
///     "barcode_read_length": 16,
///     "barcode_read_offset": 0,
///     "barcode_read_type": "R1",
///     "barcode_whitelist": "737K-august-2016",
/// ```
/// - Description and name for the chemistry
/// ``` text
///     "description": "Single Cell V(D)J",
///     "name": "SCVDJ",
/// ```
/// - Specify the `endedness` of the product. This would be `three_prime` for Single Cell
/// 3' Gene expression and `five_prime` for VDJ and 5' Gene expression
/// ``` text
///     "endedness": "five_prime",
/// ```
/// - Every RNA product will have cDNA sequences in
/// one or both of read1 and read2. Hence an `rna_read`
/// definition is necessary for each chemistry. `rna_read2`
/// is optional. For V(D)J, `rna_read` would be the bases in read1
/// beyond the spacer sequence and `rna_read2` would be all the
/// bases in read2. It is not necessary that `rna_read` correspond to
/// read1. For example, in Single Cell 5' R2-only mode, `rna_read`
/// would be all of read2 and `rna_read2` would be empty.
/// ``` text
///     "rna": {
///        "length": null,
///        "offset": 41,
///        "read_type": "R1",
///     }
/// ```
/// - Optionally specify `rna_read2`
/// ``` text
///     "rna2": {
///        "length": null,
///        "offset": 0,
///        "read_type": "R2",
///     }
/// ```
/// - Strandedness is `+` when the rna_read and the transcript are
/// expected to be in the same orientation and `-` otherwise.
/// ``` text
///     "strandedness": "+",
/// ```
/// - UMI is present in the bases 16-25 of read1
/// ``` text
///     "umi": {
///         "length": 10,
///         "min_length": null,
///         "offset": 16,
///         "read_type": "R1"
///     }
/// }
/// ```
#[derive(Serialize, Deserialize, Clone, PartialEq, Eq, Debug, MartianStruct)]
#[serde(deny_unknown_fields)]
pub struct ChemistryDef {
    pub name: ChemistryName,
    pub description: String,
    #[mro_type = "string"]
    pub endedness: Option<WhichEnd>,
    #[mro_type = "string"]
    pub strandedness: ReqStrand,
    barcode: Vec<BarcodeReadComponent>,
    pub umi: Vec<UmiReadComponent>,
    pub rna: RnaReadComponent,
    pub rna2: Option<RnaReadComponent>,
    barcode_extraction: Option<BarcodeExtraction>,
}

impl JsonReport for ChemistryDef {
    fn to_json_reporter(&self) -> JsonReporter {
        let Value::Object(map) = serde_json::to_value(self).unwrap() else {
            unreachable!();
        };
        map.into_iter()
            // Only report the newly introduced key if it's not null
            .filter(|(k, v)| !(k == "barcode_extraction" && v == &Value::Null))
            .map(|(k, v)| (format!("chemistry_{k}"), v))
            .collect()
    }
}

impl ChemistryDef {
    /// Return whether this chemistry is RTL.
    pub fn is_rtl(&self) -> Option<bool> {
        self.name.is_rtl()
    }

    /// Select a statically-defined chemistry definition by name.
    /// Panics if there is no chemistry def with this name.
    pub fn named(name: ChemistryName) -> Self {
        get_chemistry_def(name)
            .cloned()
            .unwrap_or_else(|| panic!("no chemistry definition found for {name}"))
    }

    pub fn is_paired_end(&self) -> bool {
        self.rna2.is_some()
    }

    pub fn barcode_construct(&self) -> BarcodeConstruct<&BarcodeReadComponent> {
        bc_vec_to_construct(&self.barcode)
    }

    /// Return the overhang multiplexing barcode read component,
    /// if one exists for this chemistry
    pub fn overhang_read_barcode(&self) -> Option<&BarcodeReadComponent> {
        bc_vec_to_overhang_read_component(&self.barcode)
    }

    pub fn barcode_range(&self) -> BarcodeConstruct<RpRange> {
        if matches!(
            self.barcode_extraction,
            Some(BarcodeExtraction::JointBc1Bc2 { .. })
        ) {
            panic!("Barcode range cannot be extracted for joint extraction using this function")
        }
        self.barcode_construct().map(BarcodeReadComponent::rprange)
    }

    pub fn barcode_read_type(&self) -> BarcodeConstruct<WhichRead> {
        self.barcode_construct().map(|bc| bc.read_type)
    }

    pub fn barcode_whitelist(&self) -> BarcodeConstruct<&WhitelistSpec> {
        self.barcode_construct()
            .map(BarcodeReadComponent::whitelist)
    }

    pub fn min_read_length(&self, which: WhichRead) -> usize {
        // Minimum number of bases required for alignment/assembly
        const MIN_RNA_BASES: usize = 25;

        // this is a special case because there's 15bp of TSO at the beginning
        //   of R1 which is not really alignable
        if self.name == ChemistryName::FivePrimePE && which == WhichRead::R1 {
            return 81;
        }

        if self.name == ChemistryName::FivePrimePEV3 && which == WhichRead::R1 {
            return 83;
        }

        let mut result = 0;
        for bc in &self.barcode {
            if bc.read_type == which {
                result = result.max(bc.offset + bc.length);
            }
        }

        for umi in &self.umi {
            if umi.read_type == which {
                result = result.max(umi.offset + umi.min_length.unwrap_or(umi.length));
            }
        }

        if self.rna.read_type == which {
            result = result.max(self.rna.offset + self.rna.min_length.unwrap_or(MIN_RNA_BASES));
        }
        if let Some(rna2) = self.rna2 {
            if rna2.read_type == which {
                result = result.max(rna2.offset + rna2.length.unwrap_or(MIN_RNA_BASES));
            }
        }

        result
    }

    pub fn barcode_extraction(&self) -> Option<&BarcodeExtraction> {
        self.barcode_extraction.as_ref()
    }

    /// Translate the probe barcode whitelist using the provided ID mapping.
    /// Exactly one barcode read component should be of the appropriate type.
    /// The translated whitelist will be written to the provided path.
    pub fn translate_probe_barcode_whitelist_with_id_map(
        &mut self,
        id_map: &TxHashMap<BarcodeId, BarcodeId>,
        translate_to: &WhitelistSource,
        path: PathBuf,
    ) -> Result<()> {
        self.barcode
            .iter_mut()
            .filter(|x| x.is_probe())
            .exactly_one()
            .unwrap()
            .translate_whitelist_with_id_map(id_map, translate_to, path)
    }
}

/// One chemistry def per library type.
pub type ChemistryDefs = HashMap<LibraryType, ChemistryDef>;
/// One chemistry spec per library type.
pub type ChemistrySpecs = HashMap<LibraryType, AutoOrRefinedChemistry>;

pub trait ChemistryDefsExt {
    /// Return the primary chemistry definition.
    fn primary(&self) -> &ChemistryDef;

    /// Return the chemistry name.
    fn name(&self) -> String;

    /// Return the chemistry description.
    fn description(&self) -> String;

    /// Return whether this chemistry is RTL.
    fn is_rtl(&self) -> Option<bool>;

    /// Return a single overhang read barcode component, if one is defined.
    /// Panic if more than one is found.
    fn overhang_read_barcode(&self) -> Option<&BarcodeReadComponent>;

    /// Return the endedness of this chemistry.
    fn endedness(&self) -> Option<WhichEnd>;
}

impl ChemistryDefsExt for ChemistryDefs {
    /// Return the only chemistry definition when there is only one.
    /// Return the GEX chemistry definition when there are multiple.
    /// Fail when there are multiple chemistry definitions and no GEX.
    fn primary(&self) -> &ChemistryDef {
        if let Ok((_library_type, chemistry_def)) = self.iter().exactly_one() {
            chemistry_def
        } else {
            self.get(&LibraryType::Gex)
                .expect("no gene expression chemistry found")
        }
    }

    /// Return a comma-separated list of the chemistry names.
    fn name(&self) -> String {
        self.iter()
            .sorted_by_key(|(&library_type, _)| library_type)
            .map(|(_, chem)| &chem.name)
            .unique()
            .join(", ")
    }

    /// Return a comma-separated list of the chemistry descriptions.
    fn description(&self) -> String {
        self.iter()
            .sorted_by_key(|(&library_type, _)| library_type)
            .map(|(_, chem)| &chem.description)
            .unique()
            .join(", ")
    }

    /// Return whether this chemistry is RTL.
    fn is_rtl(&self) -> Option<bool> {
        self.values()
            .filter_map(ChemistryDef::is_rtl)
            .dedup()
            .at_most_one()
            .unwrap()
    }

    /// Return a single overhang read barcode component, if one is defined.
    /// Panic if more than one is found.
    fn overhang_read_barcode(&self) -> Option<&BarcodeReadComponent> {
        self.values()
            .filter_map(ChemistryDef::overhang_read_barcode)
            .unique()
            .at_most_one()
            .unwrap()
    }

    /// Return the endedness of the chemistry.
    fn endedness(&self) -> Option<WhichEnd> {
        self.values()
            .filter_map(|x| x.endedness)
            .dedup()
            .at_most_one()
            .unwrap()
    }
}

#[derive(Serialize, Deserialize, Clone, PartialEq, Eq, Debug, MartianType)]
struct ReadDefs {
    strand: ReqStrand,
    endedness: Option<WhichEnd>,
    rna_read_type: WhichRead,
    rna_read_offset: usize,
    rna_read_length: Option<usize>,
    rna_read_min_length: Option<usize>,
    rna_read2_type: Option<WhichRead>,
    rna_read2_offset: Option<usize>,
    rna_read2_length: Option<usize>,
    rna_read2_min_length: Option<usize>,
    umi_read_type: WhichRead,
    umi_read_offset: usize,
    umi_read_length: usize,
    umi_read_min_length: Option<usize>,
    barcode_read_type: WhichRead,
    barcode_read_offset: usize,
    barcode_read_length: usize,
    right_probe_barcode_read_type: Option<WhichRead>,
    right_probe_barcode_read_offset: Option<usize>,
    right_probe_barcode_read_length: Option<usize>,
    right_probe_barcode_whitelist: Option<String>,
}

#[derive(Serialize, Deserialize, Clone, PartialEq, Eq, Debug, MartianStruct)]
pub struct IndexScheme {
    name: String,
    read_defs: ReadDefs,
}

impl IndexScheme {
    /// Match logic from lib/python/puppy/argshim/lena.py, function `make_custom_chemistry_def`
    pub fn to_chemistry_def(self, barcode_whitelist: &str) -> ChemistryDef {
        let name = ChemistryName::Custom;
        let description = format!("custom: {}", self.name);
        let strandedness = self.read_defs.strand;
        let endedness = self.read_defs.endedness;

        let gel_bead_barcode = BarcodeReadComponent {
            read_type: self.read_defs.barcode_read_type,
            kind: BarcodeKind::GelBead,
            offset: self.read_defs.barcode_read_offset,
            length: self.read_defs.barcode_read_length,
            whitelist: WhitelistSpec::TxtFile {
                name: barcode_whitelist.into(),
            },
        };

        let probe_barcode = if let Some(read_type) = self.read_defs.right_probe_barcode_read_type {
            let right_probe_barcode_whitelist =
                self.read_defs.right_probe_barcode_whitelist.unwrap();
            assert!(!right_probe_barcode_whitelist.is_empty());
            Some(BarcodeReadComponent {
                read_type,
                kind: BarcodeKind::RightProbe,
                offset: self.read_defs.right_probe_barcode_read_offset.unwrap(),
                length: self.read_defs.right_probe_barcode_read_length.unwrap(),
                whitelist: WhitelistSpec::TxtFile {
                    name: right_probe_barcode_whitelist,
                },
            })
        } else {
            None
        };

        let barcode = std::iter::once(gel_bead_barcode)
            .chain(probe_barcode)
            .collect();

        let umi = vec![UmiReadComponent {
            read_type: self.read_defs.umi_read_type,
            offset: self.read_defs.umi_read_offset,
            length: self.read_defs.umi_read_length,
            min_length: self.read_defs.umi_read_min_length,
            whitelist: None, // TODO: build it
        }];

        let rna = RnaReadComponent {
            read_type: self.read_defs.rna_read_type,
            offset: self.read_defs.rna_read_offset,
            length: self.read_defs.rna_read_length,
            min_length: self.read_defs.rna_read_min_length,
        };

        let rna2 = self
            .read_defs
            .rna_read2_type
            .map(|read_type| RnaReadComponent {
                read_type,
                offset: self.read_defs.rna_read2_offset.unwrap(),
                length: self.read_defs.rna_read2_length,
                min_length: self.read_defs.rna_read2_min_length,
            });

        ChemistryDef {
            name,
            description,
            endedness,
            strandedness,
            barcode,
            umi,
            rna,
            rna2,
            barcode_extraction: None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use strum::IntoEnumIterator;

    #[test]
    fn test_parse_auto() {
        assert_eq!(
            "threeprime".parse::<AutoChemistryName>().unwrap(),
            AutoChemistryName::ThreePrime
        );
        assert_eq!(
            "SC3P_auto".parse::<AutoChemistryName>().unwrap(),
            AutoChemistryName::ThreePrime
        );
    }

    #[test]
    fn test_standard_chemistry_defs() {
        for chem_name in ChemistryName::iter().filter(|name| {
            *name != ChemistryName::Custom
                && *name != ChemistryName::VdjR1
                && *name != ChemistryName::AtacV1
        }) {
            println!("Testing {chem_name}");
            let _ = ChemistryDef::named(chem_name);
        }
    }

    #[test]
    fn test_index_scheme_to_chemistry_def() {
        let index_scheme = IndexScheme {
            name: "RTL R1 test".into(),
            read_defs: ReadDefs {
                strand: ReqStrand::Reverse,
                endedness: Some(WhichEnd::ThreePrime),
                barcode_read_type: WhichRead::R1,
                barcode_read_length: 16,
                barcode_read_offset: 0,
                right_probe_barcode_read_type: Some(WhichRead::R1),
                right_probe_barcode_whitelist: Some("probe-barcodes-whitelist".into()),
                right_probe_barcode_read_offset: Some(40),
                right_probe_barcode_read_length: Some(8),
                rna_read_type: WhichRead::R2,
                rna_read_offset: 0,
                rna_read_length: Some(50),
                rna_read_min_length: None,
                rna_read2_type: None,
                rna_read2_offset: None,
                rna_read2_length: None,
                rna_read2_min_length: None,
                umi_read_type: WhichRead::R1,
                umi_read_offset: 16,
                umi_read_length: 12,
                umi_read_min_length: None,
            },
        };

        let chemistry_def = index_scheme.to_chemistry_def("737k-whitelist");

        assert_eq!(chemistry_def.description, "custom: RTL R1 test");
        assert_eq!(chemistry_def.name, ChemistryName::Custom);
        assert_eq!(chemistry_def.strandedness, ReqStrand::Reverse);
        assert_eq!(chemistry_def.endedness, Some(WhichEnd::ThreePrime));

        assert_eq!(chemistry_def.rna2, None);

        assert_eq!(chemistry_def.umi.len(), 1);
        assert_eq!(chemistry_def.umi[0].read_type, WhichRead::R1);
        assert_eq!(chemistry_def.umi[0].offset, 16);
        assert_eq!(chemistry_def.umi[0].length, 12);
        assert_eq!(chemistry_def.umi[0].min_length, None);

        assert_eq!(chemistry_def.rna.read_type, WhichRead::R2);
        assert_eq!(chemistry_def.rna.offset, 0);
        assert_eq!(chemistry_def.rna.length, Some(50));
        assert_eq!(chemistry_def.rna.min_length, None);

        assert_eq!(chemistry_def.barcode.len(), 2);
        assert_eq!(
            chemistry_def
                .barcode
                .iter()
                .filter(|b| b.kind == BarcodeKind::GelBead)
                .count(),
            1
        );
        assert_eq!(
            chemistry_def
                .barcode
                .iter()
                .filter(|b| b.kind == BarcodeKind::RightProbe)
                .count(),
            1
        );
        for barcode in chemistry_def.barcode {
            match barcode.kind {
                BarcodeKind::GelBead => {
                    assert_eq!(barcode.read_type, WhichRead::R1);
                    assert_eq!(barcode.offset, 0);
                    assert_eq!(barcode.length, 16);
                    assert_eq!(
                        barcode.whitelist,
                        WhitelistSpec::TxtFile {
                            name: "737k-whitelist".into()
                        }
                    );
                }
                BarcodeKind::RightProbe => {
                    assert_eq!(barcode.read_type, WhichRead::R1);
                    assert_eq!(barcode.offset, 40);
                    assert_eq!(barcode.length, 8);
                    assert_eq!(
                        barcode.whitelist,
                        WhitelistSpec::TxtFile {
                            name: "probe-barcodes-whitelist".into()
                        }
                    );
                }
                _ => panic!("Invalid barcode detected"),
            }
        }
    }

    #[test]
    fn test_display_impl() {
        assert_eq!(AutoChemistryName::Count.to_string(), "auto");
        assert_eq!(ChemistryName::VdjPE.to_string(), "SCVDJ");
    }

    #[test]
    fn test_auto_or_refined_chemistry() {
        assert_eq!(
            "auto".parse::<AutoOrRefinedChemistry>().unwrap(),
            AutoOrRefinedChemistry::Auto(AutoChemistryName::Count)
        );
        assert_eq!(
            "SCVDJ".parse::<AutoOrRefinedChemistry>().unwrap(),
            AutoOrRefinedChemistry::Refined(ChemistryName::VdjPE)
        );
        assert!("boo".parse::<AutoOrRefinedChemistry>().is_err());
    }

    #[test]
    fn test_serde_auto_or_refined_chem() {
        use AutoOrRefinedChemistry::{Auto, Refined};
        assert_eq!(
            serde_json::to_string(&Auto(AutoChemistryName::Count)).unwrap(),
            serde_json::to_string("auto").unwrap()
        );
        assert_eq!(
            serde_json::to_string(&Refined(ChemistryName::VdjPE)).unwrap(),
            serde_json::to_string("SCVDJ").unwrap()
        );
        assert_eq!(
            serde_json::from_str::<AutoOrRefinedChemistry>(r#""threeprime""#).unwrap(),
            Auto(AutoChemistryName::ThreePrime)
        );
        assert_eq!(
            serde_json::from_str::<AutoOrRefinedChemistry>(r#""SC3Pv3""#).unwrap(),
            Refined(ChemistryName::ThreePrimeV3)
        );
    }

    #[test]
    fn test_martian_type_auto_or_refined_chem() {
        assert_eq!(
            AutoOrRefinedChemistry::as_martian_primary_type(),
            MartianPrimaryType::Str
        );
    }

    #[test]
    fn test_from_auto() {
        for auto_chem in AutoChemistryName::iter() {
            let c: AutoOrRefinedChemistry = auto_chem.into();
            assert_eq!(c, AutoOrRefinedChemistry::Auto(auto_chem));
        }
    }

    #[test]
    fn test_from_refined() {
        for chem in ChemistryName::iter() {
            let c: AutoOrRefinedChemistry = chem.into();
            assert_eq!(c, AutoOrRefinedChemistry::Refined(chem));
        }
    }

    #[test]
    fn test_min_read_length() {
        use ChemistryName::{
            FivePrimeHT, FivePrimePE, FivePrimeR1, FivePrimeR2, ThreePrimeV1, ThreePrimeV2,
            ThreePrimeV3, ThreePrimeV3HT, ThreePrimeV3LT, VdjPE, VdjR2,
        };
        use WhichRead::{I1, R1, R2};
        assert_eq!(ChemistryDef::named(ThreePrimeV1).min_read_length(I1), 14);
        assert_eq!(ChemistryDef::named(ThreePrimeV1).min_read_length(R1), 25);
        assert_eq!(ChemistryDef::named(ThreePrimeV1).min_read_length(R2), 10);

        assert_eq!(ChemistryDef::named(ThreePrimeV2).min_read_length(R1), 26);
        assert_eq!(ChemistryDef::named(ThreePrimeV2).min_read_length(R2), 25);

        assert_eq!(ChemistryDef::named(ThreePrimeV3).min_read_length(R1), 26);
        assert_eq!(ChemistryDef::named(ThreePrimeV3).min_read_length(R2), 25);

        assert_eq!(ChemistryDef::named(ThreePrimeV3LT).min_read_length(R1), 26);
        assert_eq!(ChemistryDef::named(ThreePrimeV3LT).min_read_length(R2), 25);

        assert_eq!(ChemistryDef::named(ThreePrimeV3HT).min_read_length(R1), 26);
        assert_eq!(ChemistryDef::named(ThreePrimeV3HT).min_read_length(R2), 25);

        assert_eq!(ChemistryDef::named(FivePrimeR1).min_read_length(R1), 66);
        assert_eq!(ChemistryDef::named(FivePrimeR1).min_read_length(R2), 0);

        assert_eq!(ChemistryDef::named(FivePrimeR2).min_read_length(R1), 26);
        assert_eq!(ChemistryDef::named(FivePrimeR2).min_read_length(R2), 25);

        assert_eq!(ChemistryDef::named(FivePrimeHT).min_read_length(R1), 26);
        assert_eq!(ChemistryDef::named(FivePrimeHT).min_read_length(R2), 25);

        assert_eq!(ChemistryDef::named(FivePrimePE).min_read_length(R1), 81);
        assert_eq!(ChemistryDef::named(FivePrimePE).min_read_length(R2), 25);

        assert_eq!(ChemistryDef::named(VdjPE).min_read_length(R1), 66);
        assert_eq!(ChemistryDef::named(VdjPE).min_read_length(R2), 25);

        assert_eq!(ChemistryDef::named(VdjR2).min_read_length(R1), 26);
        assert_eq!(ChemistryDef::named(VdjR2).min_read_length(R2), 25);
    }

    #[test]
    fn test_chem_json_report() {
        insta::assert_json_snapshot!(
            ChemistryDef::named(ChemistryName::ThreePrimeV3).to_json_reporter()
        );
    }

    #[test]
    fn test_vec_to_construct_1() {
        let barcode = vec![BarcodeReadComponent {
            read_type: WhichRead::R1,
            kind: BarcodeKind::GelBead,
            offset: 0,
            length: 16,
            whitelist: WhitelistSpec::TxtFile {
                name: "737K-august-2016".to_string(),
            },
        }];
        assert_eq!(
            bc_vec_to_construct(&barcode),
            BarcodeConstruct::GelBeadOnly(&barcode[0])
        );
    }

    #[test]
    fn test_vec_to_construct_2() {
        let barcode = vec![
            BarcodeReadComponent {
                read_type: WhichRead::R1,
                kind: BarcodeKind::GelBead,
                offset: 0,
                length: 16,
                whitelist: WhitelistSpec::TxtFile {
                    name: "737K-august-2016".to_string(),
                },
            },
            BarcodeReadComponent {
                read_type: WhichRead::R2,
                kind: BarcodeKind::RightProbe,
                offset: 40,
                length: 8,
                whitelist: WhitelistSpec::TxtFile {
                    name: "turtle-scheme1a-right".to_string(),
                },
            },
        ];
        assert_eq!(
            bc_vec_to_construct(&barcode),
            BarcodeConstruct::new_gel_bead_and_probe(&barcode[0], &barcode[1])
        );
    }

    #[test]
    #[should_panic]
    fn test_vec_to_construct_3() {
        let barcode = vec![
            BarcodeReadComponent {
                read_type: WhichRead::R1,
                kind: BarcodeKind::GelBead,
                offset: 0,
                length: 16,
                whitelist: WhitelistSpec::TxtFile {
                    name: "737K-august-2016".to_string(),
                },
            },
            BarcodeReadComponent {
                read_type: WhichRead::R2,
                kind: BarcodeKind::LeftProbe,
                offset: 0,
                length: 8,
                whitelist: WhitelistSpec::TxtFile {
                    name: "turtle-scheme1a-left".to_string(),
                },
            },
            BarcodeReadComponent {
                read_type: WhichRead::R2,
                kind: BarcodeKind::RightProbe,
                offset: 40,
                length: 8,
                whitelist: WhitelistSpec::TxtFile {
                    name: "turtle-scheme1a-right".to_string(),
                },
            },
        ];
        bc_vec_to_construct(&barcode);
    }

    #[test]
    fn test_vec_to_construct_4() {
        let barcode = vec![
            BarcodeReadComponent {
                read_type: WhichRead::R1,
                kind: BarcodeKind::SpotSegment,
                offset: 35,
                length: 7,
                whitelist: WhitelistSpec::TxtFile {
                    name: "omni-part1.txt".to_string(),
                },
            },
            BarcodeReadComponent {
                read_type: WhichRead::R1,
                kind: BarcodeKind::SpotSegment,
                offset: 49,
                length: 7,
                whitelist: WhitelistSpec::TxtFile {
                    name: "omni-part2.txt".to_string(),
                },
            },
            BarcodeReadComponent {
                read_type: WhichRead::R1,
                kind: BarcodeKind::SpotSegment,
                offset: 63,
                length: 7,
                whitelist: WhitelistSpec::TxtFile {
                    name: "omni-part3.txt".to_string(),
                },
            },
            BarcodeReadComponent {
                read_type: WhichRead::R1,
                kind: BarcodeKind::SpotSegment,
                offset: 77,
                length: 7,
                whitelist: WhitelistSpec::TxtFile {
                    name: "omni-part4.txt".to_string(),
                },
            },
        ];

        assert_eq!(
            bc_vec_to_construct(&barcode),
            BarcodeConstruct::Segmented(Segments {
                segment1: &barcode[0],
                segment2: &barcode[1],
                segment3: Some(&barcode[2]),
                segment4: Some(&barcode[3]),
            })
        );

        assert_eq!(
            bc_vec_to_construct(&barcode[0..2]),
            BarcodeConstruct::Segmented(Segments {
                segment1: &barcode[0],
                segment2: &barcode[1],
                segment3: None,
                segment4: None,
            })
        );
    }

    #[test]
    #[should_panic]
    fn test_vec_to_construct_illegal_1() {
        let barcode = vec![];
        bc_vec_to_construct(&barcode);
    }

    #[test]
    #[should_panic]
    fn test_vec_to_construct_illegal_2() {
        let barcode = vec![
            BarcodeReadComponent {
                read_type: WhichRead::R1,
                kind: BarcodeKind::GelBead,
                offset: 0,
                length: 16,
                whitelist: WhitelistSpec::TxtFile {
                    name: "737K-august-2016".to_string(),
                },
            },
            BarcodeReadComponent {
                read_type: WhichRead::R1,
                kind: BarcodeKind::GelBead,
                offset: 16,
                length: 16,
                whitelist: WhitelistSpec::TxtFile {
                    name: "737K-august-2016".to_string(),
                },
            },
        ];
        bc_vec_to_construct(&barcode);
    }

    #[test]
    #[should_panic]
    fn test_vec_to_construct_illegal_3() {
        let barcode = vec![
            BarcodeReadComponent {
                read_type: WhichRead::R2,
                kind: BarcodeKind::RightProbe,
                offset: 0,
                length: 16,
                whitelist: WhitelistSpec::TxtFile {
                    name: "turtle-scheme1-right".to_string(),
                },
            },
            BarcodeReadComponent {
                read_type: WhichRead::R2,
                kind: BarcodeKind::RightProbe,
                offset: 16,
                length: 16,
                whitelist: WhitelistSpec::TxtFile {
                    name: "turtle-scheme1-right".to_string(),
                },
            },
        ];
        bc_vec_to_construct(&barcode);
    }

    #[test]
    fn test_whitelist_spec() {
        assert_eq!(
            serde_json::from_str::<WhitelistSpec>(
                r#"{
                    "name": "3M-february-2018"
                }"#,
            )
            .unwrap(),
            WhitelistSpec::TxtFile {
                name: "3M-february-2018".into(),
            }
        );

        assert_eq!(
            serde_json::from_str::<WhitelistSpec>(
                r#"{
                    "slide": "heildelberg_build2",
                    "part": "bc1"
                }"#,
            )
            .unwrap(),
            WhitelistSpec::SlideFile {
                slide: "heildelberg_build2".into(),
                part: slide_design::OligoPart::Bc1,
            }
        );
    }

    #[test]
    fn test_umi_spec() {
        assert_eq!(
            serde_json::from_str::<UmiReadComponent>(
                r#"{
                    "length": 12,
                    "min_length": 10,
                    "offset": 16,
                    "read_type": "R1"
                }"#,
            )
            .unwrap(),
            UmiReadComponent {
                read_type: WhichRead::R1,
                length: 12,
                min_length: Some(10),
                offset: 16,
                whitelist: None,
            }
        );
        assert_eq!(
            serde_json::from_str::<UmiReadComponent>(
                r#"{
                    "length": 4,
                    "offset": 7,
                    "read_type": "R1",
                    "whitelist": {
                        "slide": "test",
                        "part": "bc1",
                        "translation": "single_base"
                    }
                }"#,
            )
            .unwrap(),
            UmiReadComponent {
                read_type: WhichRead::R1,
                length: 4,
                min_length: None,
                offset: 7,
                whitelist: Some(UmiWhitelistSpec {
                    slide: "test".into(),
                    part: slide_design::OligoPart::Bc1,
                    translation: UmiTranslation::SingleBase,
                }),
            }
        );
    }
}
