use anyhow::Result;
use martian_derive::MartianType;
use serde::{Deserialize, Serialize};
use strum::EnumString;

#[derive(
    Clone,
    Copy,
    Debug,
    Deserialize,
    EnumString,
    Eq,
    Hash,
    MartianType,
    PartialEq,
    PartialOrd,
    Ord,
    Serialize,
)]
pub(crate) enum FeatureType {
    #[strum(serialize = "Gene Expression")]
    #[strum(serialize = "gene_expression")]
    Gene,
    #[strum(serialize = "Antibody Capture")]
    #[strum(serialize = "antibody_capture")]
    Antibody,
    #[strum(serialize = "Antigen Capture")]
    #[strum(serialize = "antigen_capture")]
    Antigen,
    #[strum(serialize = "Multiplexing Capture")]
    #[strum(serialize = "FEATURETEST")]
    #[strum(serialize = "Multiplexing Tag Capture")]
    #[strum(serialize = "multiplexing_capture")]
    Multiplexing,
    #[strum(serialize = "CRISPR Guide Capture")]
    #[strum(serialize = "crispr_guide_capture")]
    Crispr,
    #[strum(serialize = "Custom")]
    #[strum(serialize = "custom")]
    Custom,
}

impl FeatureType {
    /// Return a space-separated string representation of this feature type.
    pub(crate) fn as_str(&self) -> &'static str {
        #[allow(clippy::enum_glob_use)]
        use FeatureType::*;
        match self {
            Antibody => "Antibody Capture",
            Antigen => "Antigen Capture",
            Crispr => "CRISPR Guide Capture",
            Custom => "Custom",
            Gene => "Gene Expression",
            Multiplexing => "Multiplexing Capture",
        }
    }

    /// Return an underscore-separated lowercase string representation of this feature type.
    pub(crate) fn lc(&self) -> &'static str {
        #[allow(clippy::enum_glob_use)]
        use FeatureType::*;
        match self {
            Gene => "gene_expression",
            Antibody => "antibody_capture",
            Antigen => "antigen_capture",
            Multiplexing => "multiplexing_capture",
            Crispr => "crispr_guide_capture",
            Custom => "custom",
        }
    }
}

impl std::fmt::Display for FeatureType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
        write!(f, "{}", self.as_str())
    }
}
