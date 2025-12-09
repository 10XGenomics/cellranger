// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
#![expect(missing_docs)]

use debruijn::dna_string::DnaString;
use itertools::Itertools;
use martian_derive::MartianStruct;
use martian_filetypes::json_file::JsonFile;
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use vdj_ann::annotate::{Annotation, ContigAnnotation};
use vdj_types::VdjReceptor;

/// Metadata for a single clonotyping dataset.
#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct Dataset {
    /// Path to contig annotations for this dataset.
    pub file: JsonFile<Vec<ContigAnnotation>>,
    /// Short name for this donor.
    pub donor_id: String,
    /// Sample short name.
    pub origin_id: String,
}

#[derive(Debug, Default)]
pub struct InputSpec {
    pub receptor: VdjReceptor,
    pub origin_info: Vec<Dataset>,
}

/// Configurable clonotyping parameters.
#[derive(Debug)]
pub struct ClonotypingConfig {
    /// Where we should fetch input data.
    /// This also specifies which receptor we're working with.
    pub input: InputSpec,

    /// Path to reference.
    pub refname: String,
    /// Max number of threads to use for parallel processing.
    pub max_cores: Option<usize>,

    /// Path to donor reference output file.
    pub dref_file: String,
    /// Path to protobuf output file.
    pub proto: String,
    /// Optional path to a json file containing metadata.
    pub proto_metadata: String,
    /// Optional path to write out barcode fate.
    pub fate_file: String,

    /// Cell/clonotype filtering configuration.
    pub filter: ClonotypingFilterConfig,

    /// allow cells from different donors to be placed in the same clonotype
    pub mix_donors: bool,

    /// Break up clonotypes than have `split_max_chains` chains or more
    pub split_max_chains: usize,
}

#[derive(Debug, PartialEq)]
pub struct ClonotypingFilterConfig {
    /// heavy -> weak chain graph filtering
    pub graph: bool,
    /// cross-dataset filtering
    pub cross_dataset: bool,
    /// umi count filter
    pub umi_count: bool,
    /// umi ratio filter
    pub umi_ratio: bool,
    /// filter weak chains from clonotypes
    pub weak_chains: bool,
    /// filter weak foursies
    pub weak_foursies: bool,
    /// filter putative doublets
    pub doublet: bool,
    /// signature filtering
    pub signature: bool,
    /// filter out exact subclonotypes having a weak base
    pub qual: bool,
}

impl Default for ClonotypingFilterConfig {
    fn default() -> Self {
        Self {
            graph: true,
            cross_dataset: true,
            umi_count: true,
            umi_ratio: true,
            weak_chains: true,
            weak_foursies: true,
            doublet: true,
            signature: true,
            qual: true,
        }
    }
}

impl EncloneControl {
    pub fn is_bcr(&self) -> bool {
        self.cr_opt.input.receptor == VdjReceptor::IG
    }

    pub fn is_tcr(&self) -> bool {
        matches!(
            self.cr_opt.input.receptor,
            VdjReceptor::TR | VdjReceptor::TRGD
        )
        // if self.bcr {
        //     assert!(!self.tcr);
        //     return false;
        // }
        // // The original logic for computing this value was based on the confusing
        // // assumption that this is always true if self.bcr is false.
        // // I've preserved this for now, but this line might make sense to add:
        // // assert!(self.tcr);
        // true
    }

    /// Return sorted/unique donor names.
    pub fn donor_list(&self) -> Vec<String> {
        self.origin_info
            .iter()
            .map(|d| d.donor_id.clone())
            .sorted()
            .dedup()
            .collect()
    }
}

/// Specification of an alternative reference sequence.
#[derive(Clone, PartialOrd, Ord, PartialEq, Eq)]
pub struct AltRef {
    pub donor: usize,
    pub ref_id: usize,
    pub alt_seq: DnaString,
    pub support: usize,
    pub is_ref: bool,
}

/// Set up control datastructure (EncloneControl).  This is stuff that is constant for a given
/// run of enclone.  If you add something to this, be sure to update the "changed" section in
/// enclone_server.rs, if needed.
pub(crate) struct EncloneControl {
    /// If non-empty, trace the provided barcode through clonotyping.
    pub trace_barcode: String,
    /// Human or mouse or unknown, determined from the reference sequence
    pub species: String,
    /// Config options used by cellranger.
    pub cr_opt: ClonotypingConfig,
    /// origin (sample) info
    pub origin_info: Vec<Dataset>,
}

pub type BarcodeContigs = Vec<Contig>;

/// Set up data structure to track clonotype data.  A Contig is for one contig;
/// a Vec<Contig> is for one barcode, and an ExactClonotype is for an exact subclonotype.
#[derive(PartialEq, Eq, Clone)] // not sure these are all needed
pub struct Contig {
    /// CDR3 DNA sequence
    pub cdr3_dna: String,
    /// length of V..J sequence
    pub len: usize,
    /// start of V on full contig sequence
    pub v_start: usize,
    /// stop of aligned V on full contig sequence
    pub v_stop: usize,
    /// stop of aligned V on reference V
    pub v_stop_ref: usize,
    /// start of aligned D on full contig sequence
    pub d_start: Option<usize>,
    /// start of aligned J on full contig sequence
    pub j_start: usize,
    /// start of aligned J on reference J
    pub j_start_ref: usize,
    /// stop of J on full contig sequence
    pub j_stop: usize,
    /// start of C on full contig sequence
    pub c_start: Option<usize>,
    /// full contig sequence
    pub full_seq: Vec<u8>,
    /// index of 5'-UTR in ref file if found
    pub u_ref_id: Option<usize>,
    /// index of V segment reference sequence in ref file
    pub v_ref_id: usize,
    /// index of D segment reference sequence in ref file
    pub d_ref_id: Option<usize>,
    /// index of J segment reference sequence in ref file
    pub j_ref_id: usize,
    /// index of C segment reference sequence in ref file
    pub c_ref_id: Option<usize>,
    /// name of C segment reference sequence in ref file
    pub c_ref_name: Option<String>,
    /// start position in bases of FWR1 on V..J
    pub fr1_start: usize,
    /// start position in bases of CDR1 on V..J
    pub cdr1_start: Option<usize>,
    /// start position in bases of FWR2 on V..J
    pub fr2_start: Option<usize>,
    /// start position in bases of CDR2 on V..J
    pub cdr2_start: Option<usize>,
    /// start position in bases of FWR3 on V..J
    pub fr3_start: Option<usize>,
    /// CDR3 amino acid sequence
    pub cdr3_aa: String,
    /// start position in bases of CDR3 on V..J
    pub cdr3_start: usize,
    /// quality scores, truncated to V..J
    pub quals: Vec<u8>,
    /// quality scores
    pub full_quals: Vec<u8>,
    /// barcode
    pub barcode: String,
    /// name of contig
    pub tigname: String,
    /// true if this is IGH or TRB (or TRD in gamma/delta mode)
    pub left: bool,
    /// index of dataset
    pub dataset_index: usize,
    /// index of donor
    pub donor_index: Option<usize>,
    /// number of UMIs supporting contig
    pub umi_count: usize,
    /// number of reads supporting contig
    pub read_count: usize,
    /// e.g. IGH
    pub chain_type: String,
    /// V annotation (one or two entries), for V..J
    pub annv: Vec<Annotation>,
    // ctl.mix_donors for sorting
    pub mix_donors: bool,
}

impl Contig {
    pub fn seq(&self) -> &[u8] {
        &self.full_seq[self.v_start..self.j_stop]
    }
}

impl PartialOrd for Contig {
    fn partial_cmp(&self, other: &Contig) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Contig {
    fn cmp(&self, other: &Contig) -> Ordering {
        assert_eq!(self.mix_donors, other.mix_donors);

        // Order by cdr3_dna.

        if self.cdr3_dna < other.cdr3_dna {
            return Ordering::Less;
        } else if self.cdr3_dna > other.cdr3_dna {
            return Ordering::Greater;

        // Order by chain length.
        } else if self.len < other.len {
            return Ordering::Less;
        } else if self.len > other.len {
            return Ordering::Greater;

        // Order by chain sequence.
        } else if self.seq() < other.seq() {
            return Ordering::Less;
        } else if self.seq() > other.seq() {
            return Ordering::Greater;
        }

        // order by constant ref region presence
        if self.c_ref_id.is_none() && other.c_ref_id.is_some() {
            return Ordering::Less;
        } else if other.c_ref_id.is_none() && self.c_ref_id.is_some() {
            return Ordering::Greater;

        // Order by constant region name.
        } else if self.c_ref_id.is_some()
            && other.c_ref_id.is_some()
            && self.c_ref_name.clone().unwrap() < other.c_ref_name.clone().unwrap()
        {
            return Ordering::Less;
        } else if self.c_ref_id.is_some()
            && other.c_ref_id.is_some()
            && self.c_ref_name.clone().unwrap() > other.c_ref_name.clone().unwrap()
        {
            return Ordering::Greater;

        // Order by JC delta.
        } else if self.c_start.is_some()
            && other.c_start.is_some()
            && self.c_start.unwrap() + other.j_stop < other.c_start.unwrap() + self.j_stop
        {
            return Ordering::Less;
        } else if self.c_start.is_some()
            && other.c_start.is_some()
            && self.c_start.unwrap() + other.j_stop > other.c_start.unwrap() + self.j_stop
        {
            return Ordering::Greater;

        // Order by donor if MIX_DONORS option used.
        } else if self.mix_donors && self.donor_index < other.donor_index {
            return Ordering::Less;
        } else if self.mix_donors && self.donor_index > other.donor_index {
            return Ordering::Greater;
        }

        Ordering::Equal
    }
}

/// The ExactClonotype data structure stores information that could be exhibited as a
/// Vec<Vec<Vec<Contig>>>, but it avoids repetition of identical data.
///
/// CellContig: data for each cell
/// SharedContig: shared data
#[derive(Clone)]
pub struct CellContig {
    /// quality scores, truncated to V..J
    pub quals: Vec<u8>,
    /// barcode
    pub barcode: String,
    /// name of contig
    pub tigname: String,
    /// index of dataset
    pub dataset_index: usize,
    /// index of donor
    pub donor_index: Option<usize>,
    /// number of UMIs supporting contig
    pub umi_count: usize,
    /// number of reads supporting contig
    pub read_count: usize,
}

#[derive(Clone, Default)]
pub struct Junction {
    /// junction alignment complexity
    pub hcomp: usize,
}

/// The ExactClonotype data structure stores information that could be exhibited as a
/// Vec<Vec<Vec<Contig>>>, but it avoids repetition of identical data.
///
/// CellContig: data for each cell
/// SharedContig: shared data
#[derive(Clone)]
pub struct SharedContig {
    /// CDR3 DNA sequence
    pub cdr3_dna: String,
    /// V..J contig subsequence
    pub seq: Vec<u8>,
    /// V..J, possibly with mod 3 del
    pub seq_del: Vec<u8>,
    /// V..J, possibly with mod 3 del at mod 3 start
    pub seq_del_amino: Vec<u8>,
    /// amino acid sequence, after removing indel if present
    pub aa_mod_indel: Vec<u8>,
    /// insertions in V..J (currently at most one) = {(pos, seq)}
    /// **before** the given position
    pub ins: Vec<(usize, Vec<u8>)>,

    /// full contig sequence (consensus)
    pub full_seq: Vec<u8>,
    /// start of V on full contig sequence
    pub v_start: usize,
    /// stop of aligned V on full contig sequence
    pub v_stop: usize,
    /// stop of aligned V on reference V
    pub v_stop_ref: usize,
    /// start of aligned J on full contig sequence
    pub j_start: usize,
    /// start of aligned J on reference J
    pub j_start_ref: usize,
    /// stop of J on full contig sequence
    pub j_stop: usize,
    /// index of 5'-UTR in ref file if found
    pub u_ref_id: Option<usize>,
    /// index of V segment reference sequence in ref file
    pub v_ref_id: usize,
    /// optional index into alt_refs
    pub v_ref_id_donor: Option<usize>,
    /// donor id for v_ref_id_donor
    pub v_ref_id_donor_donor: Option<usize>,
    /// alt ref id for donor id for v_ref_id_donor
    pub v_ref_id_donor_alt_id: Option<usize>,
    /// index of D segment reference sequence in ref file
    pub d_ref_id: Option<usize>,
    /// index of J segment reference sequence in ref file
    pub j_ref_id: usize,
    /// index of C segment reference sequence in ref file
    pub c_ref_id: Option<usize>,
    /// start position in bases of FWR1 on V..J
    pub fr1_start: usize,
    /// start position in bases of CDR1 on V..J
    pub cdr1_start: Option<usize>,
    /// start position in bases of FWR2 on V..J
    pub fr2_start: Option<usize>,
    /// start position in bases of CDR2 on V..J
    pub cdr2_start: Option<usize>,
    /// start position in bases of FWR3 on V..J
    pub fr3_start: Option<usize>,
    /// CDR3 amino acid sequence
    pub cdr3_aa: String,
    /// start position in bases of CDR3 on V..J
    pub cdr3_start: usize,
    /// true if this is IGH or TRB (or TRD in gamma/delta mode)
    pub left: bool,
    /// e.g. IGH
    pub chain_type: String,
    /// V annotation (one or two entries), for V..J
    pub annv: Vec<Annotation>,
    /// reference V segment (possibly donor allele)
    pub vs: DnaString,
    /// reference J segment
    pub js: DnaString,
    pub inkt_alpha_chain_gene_match: bool,
    pub inkt_alpha_chain_junction_match: bool,
    pub inkt_beta_chain_gene_match: bool,
    pub inkt_beta_chain_junction_match: bool,
    pub mait_alpha_chain_gene_match: bool,
    pub mait_alpha_chain_junction_match: bool,
    pub mait_beta_chain_gene_match: bool,
    pub mait_beta_chain_junction_match: bool,
    pub jun: Junction,
}

#[derive(Clone)]
pub struct ExactClonotype {
    /// clone info that is shared
    pub share: Vec<SharedContig>,
    /// clone info, excluding shared stuff
    pub clones: Vec<Vec<CellContig>>,
}

impl ExactClonotype {
    pub fn ncells(&self) -> usize {
        self.clones.len()
    }
}

/// Define clonotype info data structure.  The fact that we have multiple data structures
/// encapsulating the same info is legacy and should be cleaned up.
///
/// The vectors in a CloneInfo object mostly have length two.  The exceptions are in
/// improper clones (not having chains of both types).
#[derive(Eq, Ord, PartialEq, PartialOrd, Default, Clone)] // not sure we need all these
pub struct CloneInfo {
    // CAUTION: this struct is sorted by a derived impl, and this sort behavior
    // should probably be maintained if we refactor this struct to be a vec of
    // structs instead of a struct of vecs.
    /// V..J contig lengths (will sort by this)
    pub lens: Vec<usize>,
    /// contigs, truncated to V..J (with possible - chars inserted)
    /// note only need tigs in has_del case, could improve
    pub tigs: Vec<Vec<u8>>,
    /// same as tigs, but deletion shifted to mod 3 position
    /// if there is one (rare, so wasteful, should be Option)
    pub tigs_amino: Vec<Vec<u8>>,
    /// contigs, truncated to V..J, packed (doesn't show - chars)
    pub tigsp: Vec<DnaString>,
    /// if - chars inserted to represent deletion
    pub has_del: Vec<bool>,
    /// the columns of the exact_clonotype that were extracted (used?)
    pub exact_cols: Vec<usize>,
    /// index into vector of all exact subclonotypes (across origins)
    pub clonotype_index: usize,
    /// origin indices
    pub origin: Vec<usize>,
    /// reference V segments (possibly donor allele)
    pub vs: Vec<DnaString>,
    /// reference J segments
    pub js: Vec<DnaString>,
    /// ids of V segments
    pub vsids: Vec<usize>,
    /// ids of J segments
    pub jsids: Vec<usize>,
    /// cdr3 nucleotide seqs
    pub cdr3s: Vec<String>,
}

/// Every entry in a ColInfo is a vector whose number of entries is the number of chains
/// in a clonotype.
#[derive(Clone, Default)]
pub struct ColInfo {
    pub uids: Vec<Option<usize>>,
    pub vids: Vec<usize>,
    pub vpids: Vec<Option<usize>>,
    pub dids: Vec<Option<usize>>,
    pub jids: Vec<usize>,
    pub cids: Vec<Option<usize>>,
    pub seqss: Vec<Vec<Vec<u8>>>,
    pub mat: Vec<Vec<Option<usize>>>,
}

// Potential join structure.

#[derive(Default)]
pub struct PotentialJoin {
    pub k1: usize,
    pub k2: usize,
    pub cd: isize,
    pub shares: Vec<isize>,
    pub err: bool,
}
