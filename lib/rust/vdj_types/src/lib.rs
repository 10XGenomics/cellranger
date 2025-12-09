//! vdj_types
// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.
#![expect(missing_docs)]

use serde::{Deserialize, Serialize};
use std::convert::TryFrom;
use std::fmt;
use std::str::FromStr;

// From https://danielkeep.github.io/tlborm/book/blk-counting.html
macro_rules! replace_expr {
    ($_t:tt $sub:expr) => {
        $sub
    };
}

macro_rules! count_tts {
    ($($tts:tt)*) => {0usize $(+ replace_expr!($tts 1usize))*};
}

macro_rules! make_enum {
    (
        name: $name:ident,
        variants:[$( ($field:ident, $lit: literal) ,)*],
        const_var_name: $const_var_name:ident,
    ) => {
        pub const $const_var_name: [&str; count_tts!($($field)*)] = [
            $($lit,)*
        ];

        #[derive(
            Debug,
            Copy,
            Clone,
            PartialEq,
            Eq,
            PartialOrd,
            Ord,
            Serialize,
            Deserialize,
            Hash,
        )]
        pub enum $name {
            $(
                #[serde(rename = $lit)]
                $field,
            )*
        }

        impl $name {
            pub fn all() -> [Self; count_tts!($($field)*)] {
                [
                    $($name::$field,)*
                ]
            }
        }

        impl std::fmt::Display for $name {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> fmt::Result {
                write!(f, "{}", match self {
                    $(
                        $name::$field => $lit,
                    )*
                })
            }
        }

        impl From<$name> for &'static str {
            fn from(src: $name) -> &'static str {
                match src {
                    $(
                        $name::$field => $lit,
                    )*
                }
            }
        }

        impl std::str::FromStr for $name {
            type Err = String;

            fn from_str(s: &str) -> Result<Self, Self::Err> {
                match s {
                    $(
                        $lit => Ok($name::$field),
                    )*
                    unknown => Err(
                        format!("Unknown variant '{}' for {}. Supported variants are: [{}]", unknown, stringify!($name), $const_var_name.join(", "))
                    )
                }
            }
        }
    };
}

// FIXME: collapse with similar types in Cellranger once codebases are merged.
#[derive(Debug, PartialEq, Default, Clone, Copy)]
pub enum VdjReceptor {
    #[default]
    TR,
    TRGD,
    IG,
}

make_enum! {
    name: VdjChain,
    variants: [
        (IGH, "IGH"),
        (IGK, "IGK"),
        (IGL, "IGL"),
        (TRA, "TRA"),
        (TRB, "TRB"),
        (TRD, "TRD"),
        (TRG, "TRG"),
    ],
    const_var_name: VDJ_CHAINS,
}

/// A contig could contain different regions that belong to different
/// chains. For e.g TRD V-region with TRB J-region
#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize, Hash)]
#[serde(into = "String", try_from = "&str")]
pub enum VdjContigChain {
    Single(VdjChain),
    Multi, // Could store a bitvec here listing the chains
}

impl fmt::Display for VdjContigChain {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}",
            match self {
                VdjContigChain::Single(chain) => chain.to_string(),
                VdjContigChain::Multi => "Multi".to_string(),
            }
        )
    }
}

impl From<VdjContigChain> for String {
    fn from(contig_chain: VdjContigChain) -> String {
        contig_chain.to_string()
    }
}

impl FromStr for VdjContigChain {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "Multi" => Ok(VdjContigChain::Multi),
            single => Ok(VdjContigChain::Single(single.parse()?)),
        }
    }
}
impl TryFrom<&str> for VdjContigChain {
    type Error = String;
    fn try_from(s: &str) -> Result<Self, Self::Error> {
        s.parse()
    }
}

make_enum! {
    name: VdjRegion,
    variants: [
        (UTR, "5'UTR"), // 5′ untranslated region (5′ UTR)
        (V, "L-REGION+V-REGION"), // Variable region
        (D, "D-REGION"), // Diversity region
        (J, "J-REGION"), // Joining region
        (C, "C-REGION"), // Constant region
    ],
    const_var_name: VDJ_REGIONS,
}

// FIXME: this isn't really a great place for this, but the current dep tree
// just doesn't have one.
/// Return the number of read pairs per barcode we should process.
///
/// In paired-end mode, one read pair counts as two reads.
/// Use the default `MAX_READS_PER_BARCODE` if no override is provided.
pub fn get_max_read_pairs_per_barcode(
    is_paired_end: bool,
    max_reads_per_barcode: Option<usize>,
) -> usize {
    let limit = max_reads_per_barcode.unwrap_or(80_000);
    if is_paired_end { limit / 2 } else { limit }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::str::FromStr;

    #[test]
    fn test_vdj_region_invalid_from_str() {
        assert_eq!(
            VdjRegion::from_str("V-REGION").unwrap_err(),
            "Unknown variant 'V-REGION' for VdjRegion. Supported variants are: [5'UTR, L-REGION+V-REGION, D-REGION, J-REGION, C-REGION]"
        );
    }

    #[test]
    fn test_vdj_chain_invalid_from_str() {
        assert_eq!(
            VdjChain::from_str("").unwrap_err(),
            "Unknown variant '' for VdjChain. Supported variants are: [IGH, IGK, IGL, TRA, TRB, TRD, TRG]"
        );
    }

    #[test]
    fn test_vdj_region_from_str() {
        assert_eq!(VdjRegion::from_str("5'UTR"), Ok(VdjRegion::UTR));
        assert_eq!(VdjRegion::from_str("L-REGION+V-REGION"), Ok(VdjRegion::V));
        assert_eq!(VdjRegion::from_str("D-REGION"), Ok(VdjRegion::D));
        assert_eq!(VdjRegion::from_str("J-REGION"), Ok(VdjRegion::J));
        assert_eq!(VdjRegion::from_str("C-REGION"), Ok(VdjRegion::C));

        assert_eq!(
            serde_json::from_str::<VdjRegion>("\"5'UTR\"").unwrap(),
            VdjRegion::UTR
        );
        assert_eq!(
            serde_json::from_str::<VdjRegion>("\"L-REGION+V-REGION\"").unwrap(),
            VdjRegion::V
        );
        assert_eq!(
            serde_json::from_str::<VdjRegion>("\"D-REGION\"").unwrap(),
            VdjRegion::D
        );
        assert_eq!(
            serde_json::from_str::<VdjRegion>("\"J-REGION\"").unwrap(),
            VdjRegion::J
        );
        assert_eq!(
            serde_json::from_str::<VdjRegion>("\"C-REGION\"").unwrap(),
            VdjRegion::C
        );

        assert_eq!(serde_json::to_string(&VdjRegion::UTR).unwrap(), "\"5'UTR\"");
        assert_eq!(
            serde_json::to_string(&VdjRegion::V).unwrap(),
            "\"L-REGION+V-REGION\""
        );
        assert_eq!(
            serde_json::to_string(&VdjRegion::D).unwrap(),
            "\"D-REGION\""
        );
        assert_eq!(
            serde_json::to_string(&VdjRegion::J).unwrap(),
            "\"J-REGION\""
        );
        assert_eq!(
            serde_json::to_string(&VdjRegion::C).unwrap(),
            "\"C-REGION\""
        );
    }

    #[test]
    fn test_vdj_chain_from_str() {
        assert_eq!(VdjChain::from_str("IGH"), Ok(VdjChain::IGH));
        assert_eq!(VdjChain::from_str("IGK"), Ok(VdjChain::IGK));
        assert_eq!(VdjChain::from_str("IGL"), Ok(VdjChain::IGL));
        assert_eq!(VdjChain::from_str("TRA"), Ok(VdjChain::TRA));
        assert_eq!(VdjChain::from_str("TRB"), Ok(VdjChain::TRB));
        assert_eq!(VdjChain::from_str("TRD"), Ok(VdjChain::TRD));
        assert_eq!(VdjChain::from_str("TRG"), Ok(VdjChain::TRG));
    }

    #[test]
    fn test_vdj_contig_chain() {
        for chain in VdjChain::all() {
            let contig_chain = VdjContigChain::Single(chain);
            assert_eq!(contig_chain.to_string(), chain.to_string());
            let chain_str = serde_json::to_string(&chain).unwrap();
            assert_eq!(serde_json::to_string(&contig_chain).unwrap(), chain_str);
            assert_eq!(
                serde_json::from_str::<VdjContigChain>(&chain_str).unwrap(),
                contig_chain,
            );
            assert_eq!(
                chain.to_string().parse::<VdjContigChain>().unwrap(),
                contig_chain
            );
        }
        assert_eq!(VdjContigChain::Multi.to_string(), "Multi");
        assert_eq!(
            "Multi".parse::<VdjContigChain>().unwrap(),
            VdjContigChain::Multi
        );
        assert_eq!(
            serde_json::to_string(&VdjContigChain::Multi).unwrap(),
            "\"Multi\""
        );
    }
}
