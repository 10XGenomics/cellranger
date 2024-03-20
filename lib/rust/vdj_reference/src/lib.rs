//!
//! Data structure for parsing the VDJ reference and other helper function
//!
// Warning groups (as of rust 1.55)
#![deny(
    future_incompatible,
    nonstandard_style,
    rust_2018_compatibility,
    rust_2021_compatibility,
    rust_2018_idioms,
    unused
)]
// Other warnings (as of rust 1.55)
#![deny(
    asm_sub_register,
    bad_asm_style,
    bindings_with_variant_name,
    clashing_extern_declarations,
    confusable_idents,
    const_item_mutation,
    deprecated,
    deref_nullptr,
    drop_bounds,
    dyn_drop,
    elided_lifetimes_in_paths,
    exported_private_dependencies,
    function_item_references,
    improper_ctypes,
    improper_ctypes_definitions,
    incomplete_features,
    inline_no_sanitize,
    invalid_value,
    irrefutable_let_patterns,
    large_assignments,
    mixed_script_confusables,
    non_shorthand_field_patterns,
    no_mangle_generic_items,
    overlapping_range_endpoints,
    renamed_and_removed_lints,
    stable_features,
    temporary_cstring_as_ptr,
    trivial_bounds,
    type_alias_bounds,
    uncommon_codepoints,
    unconditional_recursion,
    unknown_lints,
    unnameable_test_items,
    unused_comparisons,
    while_true
)]

use anyhow::{bail, Result};
use bio::io::fasta::{self, Record};
use errors::{ErrorContext, HeaderErrors, VdjReferenceErrors};
use martian_derive::MartianType;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::convert::TryFrom;
use std::fmt::Debug;
use std::fs::File;
use std::io::BufRead;
use std::ops::Deref;
use std::path::{Path, PathBuf};
use std::str::{self, FromStr};
use strum::IntoEnumIterator;
use strum_macros::{Display, EnumIter};
use vdj_types::VDJ_CHAINS;

pub mod errors;
pub mod lookup;
pub use lookup::{KmerClassify, KmerClassifyStrategy};
pub use vdj_types::{VdjChain, VdjRegion};

#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize, Hash)]
pub struct FeatureId(u32);

/// Enum to denote T-cell receptors or Immunoglobulin receptors
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
    Display,
    EnumIter,
    Hash,
    MartianType,
)]

// TODO(CELLRANGER-7889) collapse this with VdjChainType in cr_types.
pub enum VdjReceptor {
    #[serde(rename = "TR")]
    TR,
    #[serde(rename = "TR_GD")]
    #[strum(to_string = "TR_GD")]
    TRGD,
    IG,
}

impl FromStr for VdjReceptor {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<VdjReceptor> {
        Ok(match s {
            "TR" | "TCR" => VdjReceptor::TR,
            "IG" => VdjReceptor::IG,
            "TR_GD" => VdjReceptor::TRGD,
            _ => bail!("Unknown variant '{s}' for chain type. Supported variants are: [TR, IG]"),
        })
    }
}

/// Create a `VdjReceptor` from a `VdjChain`
///
/// # Example
/// ```rust
/// use vdj_reference::{VdjChain, VdjReceptor};
/// let chain = VdjChain::IGK;
/// let receptor = VdjReceptor::from(chain);
/// assert_eq!(receptor, VdjReceptor::IG);
/// ```
impl From<VdjChain> for VdjReceptor {
    fn from(chain: VdjChain) -> Self {
        use VdjChain::{IGH, IGK, IGL, TRA, TRB, TRD, TRG};
        match chain {
            IGH | IGK | IGL => VdjReceptor::IG,
            TRD | TRG => VdjReceptor::TRGD,
            TRA | TRB => VdjReceptor::TR,
        }
    }
}

/// Any VDJ chain can be categorized into a `Light` chain or a `Heavy` chain.
/// The term `Heavy` refers to any chain which has a D-Region.
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
    Display,
    EnumIter,
    Hash,
)]
pub enum VdjChainCategory {
    Light,
    Heavy,
}

/// Create a `VdjChainCategory` from a `VdjChain`
///
/// # Example
/// ```rust
/// use vdj_reference::{VdjChain, VdjChainCategory};
/// let chain = VdjChain::IGH;
/// let category = VdjChainCategory::from(chain);
/// assert_eq!(category, VdjChainCategory::Heavy);
/// let category: VdjChainCategory = VdjChain::TRA.into();
/// assert_eq!(category, VdjChainCategory::Light);
/// ```
impl From<VdjChain> for VdjChainCategory {
    fn from(chain: VdjChain) -> Self {
        use VdjChain::{IGH, TRB, TRD};
        match chain {
            IGH | TRB | TRD => VdjChainCategory::Heavy, // These have a D-Region
            _ => VdjChainCategory::Light,
        }
    }
}

/// A parsed single header from the reference fasta.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct VdjReferenceHeader {
    feature_id: FeatureId,       // Globally unique, 1-based integer; e.g., 3
    record_id: String,           // Originating transcript or cDNA id; e.g., ENST00000390547
    display_name: String, // Fully qualified name suitable for display (gene*allele), e.g., TRAV1-1*01. If the allele is absent this is just gene_name.
    gene_name: String,    // E.g., TRAV1-1
    region: VdjRegion,    // U, V, D, J or C
    receptor: VdjReceptor, // E.g TR
    chain: VdjChain,      // Eg IGH
    subclass: Option<String>, // Not the same across species, so just use a string instead of Enum
    allele_name: Option<String>, // Allele (according to IMGT, for example); e.g., 01.
}

impl FromStr for VdjReferenceHeader {
    type Err = HeaderErrors;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let words: Vec<&str> = s.splitn(2, char::is_whitespace).collect();
        if words.len() != 2 {
            return Err(HeaderErrors::NoDescription);
        }

        VdjReferenceHeader::from_words(words[0], words[1])
    }
}

pub(crate) const GENE_NAME_DISALLOWED: [&str; 9] = [
    "segment", "before", "after", "UTR", "V-REGION", "D-REGION", "J-REGION", "C-REGION", "L+V",
];

impl VdjReferenceHeader {
    pub fn from_record(rec: &Record) -> Result<Self, HeaderErrors> {
        let word0 = rec.id();
        let Some(word1) = rec.desc() else {
            return Err(HeaderErrors::NoDescription);
        };
        Self::from_words(word0, word1)
    }

    fn from_words(word0: &str, word1: &str) -> Result<Self, HeaderErrors> {
        let parts_0: Vec<&str> = word0.split('|').collect();
        if parts_0.len() != 2 {
            return Err(HeaderErrors::UnexpectedIdFormat {
                num_parts: parts_0.len(),
                id: word0.to_string(),
            });
        }

        let parts_1: Vec<&str> = word1.split('|').collect();
        if parts_1.len() != 7 {
            return Err(HeaderErrors::UnexpectedDescFormat {
                num_parts: parts_1.len(),
                desc: word1.to_string(),
            });
        }

        let (display_name, allele_name) = if parts_1[6] == "00" {
            (parts_0[1].into(), None)
        } else {
            (
                format!("{}*{}", parts_0[1], parts_1[6]),
                Some(parts_1[6].into()),
            )
        };

        let find_range = |part1_idx| {
            let start =
                word0.len() + 2 + (0..part1_idx).map(|i| parts_1[i].len() + 1).sum::<usize>();
            start..start + parts_1[part1_idx].len()
        };

        const GENE_NAME_IDX: usize = 1;
        const REGION_IDX: usize = 2;
        const RECEPTOR_IDX: usize = 3;
        const CHAIN_IDX: usize = 4;

        let mut header = VdjReferenceHeader {
            feature_id: FeatureId(if parts_0[0].starts_with('>') {
                parts_0[0][1..]
                    .parse()
                    .map_err(|_| HeaderErrors::FeatureIdNotAnInteger {
                        feature_id: parts_0[0][1..].to_string(),
                    })?
            } else {
                parts_0[0]
                    .parse()
                    .map_err(|_| HeaderErrors::FeatureIdNotAnInteger {
                        feature_id: parts_0[0].to_string(),
                    })?
            }),
            record_id: parts_1[0].into(),
            display_name,
            gene_name: parts_1[GENE_NAME_IDX].into(),
            region: parts_1[REGION_IDX].parse().map_err(|e: String| {
                HeaderErrors::UnknownVariant {
                    parse_error: e.replace("VdjRegion", "region"),
                    range: find_range(REGION_IDX),
                }
            })?,
            receptor: parts_1[RECEPTOR_IDX]
                .parse()
                .map_err(|parse_error: anyhow::Error| HeaderErrors::UnknownVariant {
                    parse_error: parse_error.to_string(),
                    range: find_range(RECEPTOR_IDX),
                })?,
            chain: parts_1[CHAIN_IDX].parse().map_err(|e: String| {
                HeaderErrors::UnknownVariant {
                    parse_error: e.replace("VdjChain", "chain name"),
                    range: find_range(CHAIN_IDX),
                }
            })?,
            subclass: if parts_1[5] == "None" {
                None
            } else {
                Some(parts_1[5].into())
            },
            allele_name,
        };

        // Update receptor based on chain
        header.receptor = header.chain.into();

        // enclone/vdj_ann has additional assumptions which we will check here
        if !header.gene_name.starts_with(parts_1[CHAIN_IDX]) {
            return Err(HeaderErrors::GeneNameDoesNotStartWithChainName {
                range: find_range(GENE_NAME_IDX),
                gene_name: header.gene_name,
                chain_name: parts_1[CHAIN_IDX].to_owned(),
            });
        }
        if GENE_NAME_DISALLOWED
            .iter()
            .any(|disallowed| header.gene_name.contains(disallowed))
        {
            return Err(HeaderErrors::GeneNameNotAllowed {
                gene_name: header.gene_name,
                range: find_range(GENE_NAME_IDX),
            });
        }

        if ["TR", "IG"].iter().all(|s| header.gene_name.contains(s)) {
            return Err(HeaderErrors::GeneNameAmbiguous {
                gene_name: header.gene_name,
                range: find_range(GENE_NAME_IDX),
                exclusive_options: vec!["TR", "IG"],
            });
        }

        let chain_matches: Vec<_> = VDJ_CHAINS
            .into_iter()
            .filter(|&chain| header.gene_name.contains(chain))
            .collect();
        if chain_matches.len() != 1 {
            return Err(HeaderErrors::GeneNameAmbiguous {
                gene_name: header.gene_name,
                range: find_range(GENE_NAME_IDX),
                exclusive_options: chain_matches,
            });
        }

        Ok(header)
    }
}

/// A single entry (header + sequence) in the reference fasta
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct VdjReferenceEntry {
    header: VdjReferenceHeader,
    pub sequence: Vec<u8>,
}

impl Deref for VdjReferenceEntry {
    type Target = VdjReferenceHeader;
    fn deref(&self) -> &Self::Target {
        &self.header
    }
}

impl TryFrom<Record> for VdjReferenceEntry {
    type Error = HeaderErrors;

    fn try_from(record: Record) -> Result<Self, Self::Error> {
        Ok(VdjReferenceEntry {
            header: VdjReferenceHeader::from_record(&record)?,
            sequence: record.seq().to_vec(), // Calling seq() does not work here for multiline fasta sequence
        })
    }
}

impl VdjReferenceEntry {
    pub fn seq(&self) -> &[u8] {
        &self.sequence
    }
}

fn _regions_fa_path(ref_folder: &Path) -> PathBuf {
    ref_folder.join("fasta/regions.fa")
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct VdjReference {
    data: Vec<VdjReferenceEntry>,
}

pub(crate) const ALLOWED_NUCLEOTIDES: &str = "ACGTURYKMSWBDHVN";

impl VdjReference {
    pub fn check(ref_folder: &Path) -> Result<()> {
        let fa_file = _regions_fa_path(ref_folder);
        let fa_reader = fasta::Reader::from_file(&fa_file)?;
        let mut context = ErrorContext::new(&fa_file)?;

        let mut header_errors = Vec::new();
        let mut line_of_id = HashMap::new();

        for record in fa_reader.records() {
            let record = record.map_err(|e| VdjReferenceErrors::CannotReadFastaRecord {
                fa_file: fa_file.clone(),
                state: context.state.clone(),
                error: e,
            })?;

            context.advance(&record)?;
            for (i, base) in record.seq().iter().enumerate() {
                if !ALLOWED_NUCLEOTIDES
                    .as_bytes()
                    .iter()
                    .any(|allowed| allowed == base)
                {
                    return Err(VdjReferenceErrors::InvalidBaseInSequence {
                        fa_file,
                        state: context.state.clone(),
                        base: *base as char,
                        position: i + 1,
                    }
                    .into());
                }
            }
            match VdjReferenceEntry::try_from(record) {
                Ok(entry) => match line_of_id.get(&entry.feature_id) {
                    Some(&last_line) => {
                        header_errors.push((
                            context.state.clone(),
                            HeaderErrors::DuplicateId {
                                id: entry.feature_id.0,
                                last_line,
                            },
                        ));
                    }
                    None => {
                        line_of_id.insert(entry.feature_id, context.state.line_num());
                    }
                },
                Err(e) => header_errors.push((context.state.clone(), e)),
            }
        }
        if !header_errors.is_empty() {
            return Err(VdjReferenceErrors::InvalidHeader {
                fa_file,
                errors: header_errors,
            }
            .into());
        }
        Ok(())
    }

    pub fn from_reference_folder(ref_folder: &Path) -> Result<Self> {
        VdjReference::from_reference_fasta(&_regions_fa_path(ref_folder))
    }

    pub fn from_reference_fasta(ref_path: &Path) -> Result<Self> {
        VdjReference::from_fasta_reader(fasta::Reader::from_file(ref_path)?)
    }

    pub fn from_fasta_reader<R: BufRead>(reader: fasta::Reader<R>) -> Result<Self> {
        let mut data = Vec::new();
        let records = reader.records();
        for record in records {
            data.push(VdjReferenceEntry::try_from(record?)?);
        }
        Ok(VdjReference { data })
    }
    pub fn iter(&self) -> impl Iterator<Item = &VdjReferenceEntry> {
        self.data.iter()
    }
    pub fn iter_seq(&self) -> impl Iterator<Item = &[u8]> {
        self.data.iter().map(VdjReferenceEntry::seq)
    }
    pub fn iter_receptor_filtered(
        &self,
        receptor: VdjReceptor,
    ) -> impl Iterator<Item = &'_ VdjReferenceEntry> {
        self.data
            .iter()
            .filter(move |entry| entry.receptor == receptor)
    }
    pub fn iter_region_filtered(
        &self,
        region: VdjRegion,
    ) -> impl Iterator<Item = &'_ VdjReferenceEntry> {
        self.data.iter().filter(move |entry| entry.region == region)
    }
    pub fn iter_chain_filtered(
        &self,
        chain: VdjChain,
    ) -> impl Iterator<Item = &'_ VdjReferenceEntry> {
        self.data.iter().filter(move |entry| entry.chain == chain)
    }
}

impl<'a> IntoIterator for &'a VdjReference {
    type IntoIter = std::slice::Iter<'a, VdjReferenceEntry>;
    type Item = &'a VdjReferenceEntry;

    fn into_iter(self) -> Self::IntoIter {
        self.data.iter()
    }
}

pub trait HierarchyChoices: Sized {
    fn choices() -> Vec<Self>;
}

impl HierarchyChoices for VdjReceptor {
    fn choices() -> Vec<Self> {
        VdjReceptor::iter().collect()
    }
}

impl HierarchyChoices for VdjChainCategory {
    fn choices() -> Vec<Self> {
        VdjChainCategory::iter().collect()
    }
}

impl HierarchyChoices for VdjChain {
    fn choices() -> Vec<Self> {
        VdjChain::all().to_vec()
    }
}

impl HierarchyChoices for VdjRegion {
    fn choices() -> Vec<Self> {
        VdjRegion::all().to_vec()
    }
}

pub trait VdjHierarchy
where
    Self: Copy + PartialEq + HierarchyChoices,
{
    fn from_header(header: &VdjReferenceHeader) -> Self;
    fn from_entry(entry: &VdjReferenceEntry) -> Self {
        Self::from_header(&entry.header)
    }
}

impl VdjHierarchy for VdjReceptor {
    fn from_header(header: &VdjReferenceHeader) -> Self {
        header.receptor
    }
}

impl VdjHierarchy for VdjChain {
    fn from_header(header: &VdjReferenceHeader) -> Self {
        header.chain
    }
}

impl VdjHierarchy for VdjRegion {
    fn from_header(header: &VdjReferenceHeader) -> Self {
        header.region
    }
}

impl VdjHierarchy for VdjChainCategory {
    fn from_header(header: &VdjReferenceHeader) -> Self {
        header.chain.into()
    }
}

/// Metadata associated with the reference from the `reference.json` file in the reference directory
///
/// Only the `genomes` and `version` fields are populated. There are other fields in the json as
/// which can be added to this struct later if needed.
#[derive(Deserialize, Debug, Default)]
pub struct VdjReferenceInfo {
    pub genomes: String,
    pub version: Option<String>,
}

impl VdjReferenceInfo {
    pub fn from_reference_folder(ref_folder: &Path) -> Result<Self> {
        VdjReferenceInfo::from_reference_json(&ref_folder.join("reference.json"))
    }

    pub fn from_reference_json(json_path: &Path) -> Result<Self> {
        Ok(serde_json::from_reader(File::open(json_path)?)?)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use test_refdata::{refdata_available, refdata_path};

    #[test]
    fn test_vdj_receptor() {
        assert_eq!(VdjReceptor::from_str("TR").unwrap(), VdjReceptor::TR);
        assert_eq!(VdjReceptor::from_str("TCR").unwrap(), VdjReceptor::TR);
        assert_eq!(VdjReceptor::from_str("IG").unwrap(), VdjReceptor::IG);
        assert_eq!(serde_json::to_string(&VdjReceptor::TR).unwrap(), r#""TR""#);
    }

    #[test]
    fn test_chain_category() {
        assert_eq!(
            VdjChainCategory::from(VdjChain::IGH),
            VdjChainCategory::Heavy
        );
        assert_eq!(
            VdjChainCategory::from(VdjChain::IGK),
            VdjChainCategory::Light
        );
        assert_eq!(
            VdjChainCategory::from(VdjChain::IGL),
            VdjChainCategory::Light
        );
        assert_eq!(
            VdjChainCategory::from(VdjChain::TRA),
            VdjChainCategory::Light
        );
        assert_eq!(
            VdjChainCategory::from(VdjChain::TRB),
            VdjChainCategory::Heavy
        );
        assert_eq!(
            VdjChainCategory::from(VdjChain::TRG),
            VdjChainCategory::Light
        );
        assert_eq!(
            VdjChainCategory::from(VdjChain::TRD),
            VdjChainCategory::Heavy
        );
    }

    #[test]
    fn test_chain_receptor() {
        assert_eq!(VdjReceptor::from(VdjChain::IGH), VdjReceptor::IG);
        assert_eq!(VdjReceptor::from(VdjChain::IGK), VdjReceptor::IG);
        assert_eq!(VdjReceptor::from(VdjChain::IGL), VdjReceptor::IG);
        assert_eq!(VdjReceptor::from(VdjChain::TRA), VdjReceptor::TR);
        assert_eq!(VdjReceptor::from(VdjChain::TRB), VdjReceptor::TR);
        assert_eq!(VdjReceptor::from(VdjChain::TRG), VdjReceptor::TRGD);
        assert_eq!(VdjReceptor::from(VdjChain::TRD), VdjReceptor::TRGD);
    }

    #[test]
    fn test_parse_header_1() {
        let header_str = ">10|IGHD1-26 ENST00000390567|IGHD1-26|D-REGION|IG|IGH|None|00";
        let expected_header = VdjReferenceHeader {
            feature_id: FeatureId(10),
            record_id: "ENST00000390567".into(),
            display_name: "IGHD1-26".into(),
            gene_name: "IGHD1-26".into(),
            region: VdjRegion::D,
            receptor: VdjReceptor::IG,
            chain: VdjChain::IGH,
            subclass: None,
            allele_name: None,
        };
        let actual_header: VdjReferenceHeader = header_str.parse().unwrap();
        assert_eq!(actual_header, expected_header);
    }

    #[test]
    fn test_parse_header_2() {
        let header_str = ">40|IGHG1 ENST00000390542|IGHG1|C-REGION|IG|IGH|G1|00";
        let expected_header = VdjReferenceHeader {
            feature_id: FeatureId(40),
            record_id: "ENST00000390542".into(),
            display_name: "IGHG1".into(),
            gene_name: "IGHG1".into(),
            region: VdjRegion::C,
            receptor: VdjReceptor::IG,
            chain: VdjChain::IGH,
            subclass: Some("G1".into()),
            allele_name: None,
        };
        let actual_header: VdjReferenceHeader = header_str.parse().unwrap();
        assert_eq!(actual_header, expected_header);
    }

    #[test]
    fn test_parse_header_3() {
        let header_str = ">430|TRAJ4 ENST00000390533|TRAJ4|J-REGION|TR|TRA|None|01";
        let expected_header = VdjReferenceHeader {
            feature_id: FeatureId(430),
            record_id: "ENST00000390533".into(),
            display_name: "TRAJ4*01".into(),
            gene_name: "TRAJ4".into(),
            region: VdjRegion::J,
            receptor: VdjReceptor::TR,
            chain: VdjChain::TRA,
            subclass: None,
            allele_name: Some("01".into()),
        };
        let actual_header: VdjReferenceHeader = header_str.parse().unwrap();
        assert_eq!(actual_header, expected_header);
    }

    #[test]
    fn test_parse_header_4() {
        let header_str = ">382|IGLV5-45 ENST00000390296|IGLV5-45|L-REGION+V-REGION|IG|IGL|None|00";
        let expected_header = VdjReferenceHeader {
            feature_id: FeatureId(382),
            record_id: "ENST00000390296".into(),
            display_name: "IGLV5-45".into(),
            gene_name: "IGLV5-45".into(),
            region: VdjRegion::V,
            receptor: VdjReceptor::IG,
            chain: VdjChain::IGL,
            subclass: None,
            allele_name: None,
        };
        let actual_header: VdjReferenceHeader = header_str.parse().unwrap();
        assert_eq!(actual_header, expected_header);
    }

    #[test]
    fn test_parse_header_5() {
        let header_str = "382|IGLV5-45 ENST00000390296|IGLV5-45|L-REGION+V-REGION|IG|IGL|None|00";
        let expected_header = VdjReferenceHeader {
            feature_id: FeatureId(382),
            record_id: "ENST00000390296".into(),
            display_name: "IGLV5-45".into(),
            gene_name: "IGLV5-45".into(),
            region: VdjRegion::V,
            receptor: VdjReceptor::IG,
            chain: VdjChain::IGL,
            subclass: None,
            allele_name: None,
        };
        let actual_header: VdjReferenceHeader = header_str.parse().unwrap();
        assert_eq!(actual_header, expected_header);
    }

    #[test]
    fn test_vdj_reference_entry() {
        let fasta =
            b">9|IGHD1-20 ENST00000450276|IGHD1-20|D-REGION|IG|IGH|None|00\nGGTATAACTGGAACGAC";
        let mut reader = fasta::Reader::new(&fasta[..]).records();
        let record = reader.next().unwrap().unwrap();
        let entry = VdjReferenceEntry::try_from(record).unwrap();
        let expected_entry = VdjReferenceEntry {
            header: VdjReferenceHeader {
                feature_id: FeatureId(9),
                record_id: "ENST00000450276".into(),
                display_name: "IGHD1-20".into(),
                gene_name: "IGHD1-20".into(),
                region: VdjRegion::D,
                receptor: VdjReceptor::IG,
                chain: VdjChain::IGH,
                subclass: None,
                allele_name: None,
            },
            sequence: b"GGTATAACTGGAACGAC".to_vec(),
        };
        assert_eq!(entry, expected_entry);
    }

    #[test]
    fn test_vdj_reference_entry_multiline() {
        let fasta = b">457|TRAV1-2 ENST00000390423|TRAV1-2|5'UTR|TR|TRA|None|00\nCCCACATGAAGTGTCTACCTTCTG\nCAGACTCCAATGGCTCAGGAACTGGGAATGCAGTG\nCCAGGCTCGTGGTATCCTGCAGCAG";
        let mut reader = fasta::Reader::new(&fasta[..]).records();
        let record = reader.next().unwrap().unwrap();
        let entry = VdjReferenceEntry::try_from(record).unwrap();
        let expected_entry = VdjReferenceEntry {
            header: VdjReferenceHeader {
                feature_id: FeatureId(457),
                record_id: "ENST00000390423".into(),
                display_name: "TRAV1-2".into(),
                gene_name: "TRAV1-2".into(),
                region: VdjRegion::UTR,
                receptor: VdjReceptor::TR,
                chain: VdjChain::TRA,
                subclass: None,
                allele_name: None,
            },
            sequence: b"CCCACATGAAGTGTCTACCTTCTGCAGACTCCAATGGCTCAGGAACTGGGAATGCAGTGCCAGGCTCGTGGTATCCTGCAGCAG".to_vec(),
        };
        println!("{}", str::from_utf8(&entry.sequence).unwrap());
        assert_eq!(entry, expected_entry);
    }

    #[test]
    fn test_vdj_reference_entry_deref() {
        let fasta = b">457|TRAV1-2 ENST00000390423|TRAV1-2|5'UTR|TR|TRA|None|00\nCCCACATGAAGTGTCTACCTTCTG\nCAGACTCCAATGGCTCAGGAACTGGGAATGCAGTG\nCCAGGCTCGTGGTATCCTGCAGCAG";
        let mut reader = fasta::Reader::new(&fasta[..]).records();
        let record = reader.next().unwrap().unwrap();
        let entry = VdjReferenceEntry::try_from(record).unwrap();
        assert_eq!(entry.feature_id, FeatureId(457));
        assert_eq!(entry.record_id, "ENST00000390423".to_string());
        assert_eq!(entry.display_name, "TRAV1-2".to_string());
        assert_eq!(entry.gene_name, "TRAV1-2".to_string());
        assert_eq!(entry.region, VdjRegion::UTR);
        assert_eq!(entry.receptor, VdjReceptor::TR);
        assert_eq!(entry.chain, VdjChain::TRA);
        assert_eq!(entry.subclass, None);
        assert_eq!(entry.allele_name, None);
        assert_eq!(
            entry.sequence,
            b"CCCACATGAAGTGTCTACCTTCTGCAGACTCCAATGGCTCAGGAACTGGGAATGCAGTGCCAGGCTCGTGGTATCCTGCAGCAG"
                .to_vec()
        );
    }

    #[test]
    fn test_vdj_reference_entry_tcr_gd() {
        let fasta = b">629|TRGJ1 ENSMUST00000200495|TRGJ1|J-REGION|TR|TRG|None|00\nATAGCTCAGGTTTTCACAAGGTATTTGCAGAAGGAACTAAGCTCATAGTAATTCCCTCTG";
        let mut reader = fasta::Reader::new(&fasta[..]).records();
        let record = reader.next().unwrap().unwrap();
        let entry = VdjReferenceEntry::try_from(record).unwrap();
        let expected_entry = VdjReferenceEntry {
            header: VdjReferenceHeader {
                feature_id: FeatureId(629),
                record_id: "ENSMUST00000200495".into(),
                display_name: "TRGJ1".into(),
                gene_name: "TRGJ1".into(),
                region: VdjRegion::J,
                receptor: VdjReceptor::TRGD,
                chain: VdjChain::TRG,
                subclass: None,
                allele_name: None,
            },
            sequence: b"ATAGCTCAGGTTTTCACAAGGTATTTGCAGAAGGAACTAAGCTCATAGTAATTCCCTCTG".to_vec(),
        };
        println!("{}", str::from_utf8(&entry.sequence).unwrap());
        assert_eq!(entry, expected_entry);
    }

    #[test]
    fn test_vdj_reference_folder() {
        let reference = VdjReference::from_reference_folder(Path::new(
            "../dui_tests/test_resources/reference/vdj_GRCh38_alts_ensembl-4.0.0",
        ))
        .unwrap();
        assert_eq!(reference.data.len(), 743);
        assert_eq!(reference.iter().count(), 743);
        assert_eq!(reference.iter_seq().count(), 743);
        assert_eq!(
            reference.iter_receptor_filtered(VdjReceptor::TR).count(),
            308
        );
        assert_eq!(
            reference.iter_receptor_filtered(VdjReceptor::TRGD).count(),
            37
        );
        assert_eq!(reference.iter_region_filtered(VdjRegion::V).count(), 312);
        assert_eq!(reference.iter_chain_filtered(VdjChain::IGH).count(), 213);
    }

    #[test]
    fn test_ref_json() {
        let info: VdjReferenceInfo = serde_json::from_str(r#"{
    "fasta_hash": "10aa54c58952f920c072270fc459d1e304c2ba3c2da20ac80528064f638dc1f1",
    "genomes": "vdj_GRCh38_alts_ensembl",
    "gtf_hash": "b2b50cc12da4d2bda69207aa7fd51bf648826d0d2f39199e87922bf107d81ed0",
    "input_fasta_files": "release-94/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa",
    "input_gtf_files": "release-94/gtf/homo_sapiens/Homo_sapiens.GRCh38.94.chr_patch_hapl_scaff.gtf",
    "mkref_version": "",
    "type": "V(D)J Reference",
    "version": "4.0.0"
}"#).unwrap();
        assert_eq!(info.genomes, "vdj_GRCh38_alts_ensembl");
        assert_eq!(info.version, Some("4.0.0".into()));

        let info = VdjReferenceInfo::from_reference_folder(Path::new(
            "../dui_tests/test_resources/reference/vdj_GRCh38_alts_ensembl-4.0.0",
        ))
        .unwrap();
        assert_eq!(info.genomes, "vdj_GRCh38_alts_ensembl-4.0.0");
        assert_eq!(info.version, Some("4.0.0".into()));
    }

    #[test]
    fn test_receptor_name() {
        assert_eq!(VdjReceptor::from_str("TR_GD").unwrap(), VdjReceptor::TRGD);
    }

    #[test]
    fn test_invalid_ref() {
        for folder in [
            "../dui_tests/test_resources/reference/GRCh38-2020-A-chrM/",
            "../dui_tests/test_resources/reference/vdj_invalid_1/",
            "../dui_tests/test_resources/reference/vdj_invalid_2/",
            "../dui_tests/test_resources/reference/vdj_invalid_3/",
        ] {
            println!("----------------------------------------------------");
            _error_log(VdjReference::check(Path::new(folder)).unwrap_err());
        }
    }

    #[test]
    fn test_valid_ref() {
        if !refdata_available() {
            return;
        }
        for folder in [
            refdata_path("vdj/vdj_GRCh38_alts_ensembl-1.0.0"),
            refdata_path("vdj/vdj_GRCh38_alts_ensembl-1.0.1"),
            refdata_path("vdj/vdj_GRCh38_alts_ensembl-1.0.2"),
            refdata_path("vdj/vdj_GRCh38_alts_ensembl-2.0.0"),
            refdata_path("vdj/vdj_GRCh38_alts_ensembl-3.1.0"),
            refdata_path("vdj/vdj_GRCh38_alts_ensembl-3.1.0_OLD"),
            refdata_path("vdj/vdj_GRCh38_alts_ensembl-4.0.0"),
            refdata_path("vdj/vdj_GRCh38_alts_ensembl-5.0.0"),
            refdata_path("vdj/vdj_GRCm38_alts_ensembl-1.0.0"),
            refdata_path("vdj/vdj-GRCm38-alts-ensembl-2.2.0"),
            refdata_path("vdj/vdj_GRCm38_alts_ensembl-3.1.0"),
            refdata_path("vdj/vdj_GRCm38_alts_ensembl-4.0.0"),
            refdata_path("vdj/vdj_GRCm38_alts_ensembl-5.0.0"),
            refdata_path("vdj/vdj_IMGT_20170916-2.1.0"),
            refdata_path("vdj/vdj_IMGT_human_20200415-0.0.0"),
            refdata_path("vdj/vdj_IMGT_mouse_20171012-2.1.0"),
            refdata_path("vdj/vdj_IMGT_mouse_20180723-2.2.0"),
        ] {
            VdjReference::check(&folder).unwrap();
        }
    }

    fn _error_log(e: anyhow::Error) {
        println!("ERROR: {e}");
        for c in e.chain().skip(1) {
            println!("\tCaused by: {c}");
        }
    }
}
