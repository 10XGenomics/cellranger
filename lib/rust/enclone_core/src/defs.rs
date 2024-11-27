// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use debruijn::dna_string::DnaString;
use vdj_ann::annotate::Annotation;
use vdj_types::VdjReceptor;

/// Clonotyping algorithm heuristics.
#[derive(PartialEq, Eq)]
pub struct ClonotypeHeuristics {
    pub max_diffs: usize,
    pub max_degradation: usize,
    pub ref_v_trim: usize,
    pub ref_j_trim: usize,
}

impl Default for ClonotypeHeuristics {
    fn default() -> Self {
        Self {
            max_diffs: 1_000_000,
            max_degradation: 2,
            ref_v_trim: 15,
            ref_j_trim: 15,
        }
    }
}

/// Origin info data structure.
#[derive(Default, PartialEq, Eq)]
pub struct OriginInfo {
    // parallel vectors
    /// map dataset index to vdj path
    pub dataset_path: Vec<String>,
    /// map dataset index to dataset short name
    pub dataset_id: Vec<String>,
    /// map dataset index to donor short name
    pub donor_id: Vec<String>,
    /// map dataset id to origin (sample) short name
    pub origin_id: Vec<String>,
    /// unique-sorted list of donor short names
    pub donor_list: Vec<String>,
}

impl OriginInfo {
    // number of datasets
    pub fn n(&self) -> usize {
        self.dataset_path.len()
    }
}

#[derive(Debug, PartialEq)]
pub enum InputSpec {
    /// Explicitly-provided pair of receptor and datasets.
    Explicit(VdjReceptor, String),
    /// Path to a META file containing input spec.
    MetaFile(String),
}

/// The subset of configuration options used by Cellranger.
#[derive(Debug, PartialEq)]
pub struct CellrangerOpt {
    /// Where we should fetch input data.
    /// This also specifies which receptor we're working with.
    // FIXME: we should be able to eliminate optionality.
    pub input: Option<InputSpec>,

    // FIXME: cellranger actually always sets this to an empty vec
    // should eliminate this from this config.
    pub pre: Vec<String>,

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

    /// Cell/clonotype filtering options that may be configured from Cellranger.
    pub filter: CellrangerFilterOpt,

    /// allow cells from different donors to be placed in the same clonotype
    pub mix_donors: bool,

    // Clonotype joining options.
    /// Break up clonotypes than have `split_max_chains` chains or more
    pub split_max_chains: usize,
}

impl Default for CellrangerOpt {
    fn default() -> Self {
        Self {
            input: Default::default(),
            pre: Default::default(),
            refname: Default::default(),
            max_cores: Default::default(),
            dref_file: Default::default(),
            proto: Default::default(),
            proto_metadata: Default::default(),
            fate_file: Default::default(),
            filter: Default::default(),
            mix_donors: Default::default(),
            split_max_chains: usize::MAX,
        }
    }
}

#[derive(Debug, PartialEq)]
pub struct CellrangerFilterOpt {
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

impl Default for CellrangerFilterOpt {
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

impl CellrangerFilterOpt {
    /// Return a filter configuration that should disable all filters.
    pub fn disabled() -> Self {
        Self {
            graph: false,
            cross_dataset: false,
            umi_count: false,
            umi_ratio: false,
            weak_chains: false,
            weak_foursies: false,
            doublet: false,
            signature: false,
            qual: false,
        }
    }
}

/// Miscellaneous general options.
#[derive(PartialEq)]
pub struct GeneralOpt {
    pub receptor: VdjReceptor,
    pub trace_barcode: String,
    /// human or mouse or unknown, determined from the reference sequence
    pub species: String,
    pub jscore_match: i32,
    pub jscore_mismatch: i32,
    pub jscore_bits_multiplier: f64,
    pub jscore_gap_open: i32,
    pub jscore_gap_extend: i32,
    pub no_alt_alleles: bool,
}

impl Default for GeneralOpt {
    fn default() -> Self {
        Self {
            receptor: Default::default(),
            trace_barcode: Default::default(),
            species: Default::default(),
            jscore_match: 20,
            jscore_mismatch: -20,
            jscore_bits_multiplier: 2.2,
            jscore_gap_open: -120,
            jscore_gap_extend: -20,
            no_alt_alleles: Default::default(),
        }
    }
}

impl GeneralOpt {
    pub fn is_bcr(&self) -> bool {
        self.receptor == VdjReceptor::IG
    }

    pub fn is_tcr(&self) -> bool {
        matches!(self.receptor, VdjReceptor::TR | VdjReceptor::TRGD)
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
}

/// Allele-finding algorithmic options.
#[derive(PartialEq, Eq)]
pub struct AlleleAlgOpt {
    pub min_mult: usize,
    pub min_alt: usize,
}

impl Default for AlleleAlgOpt {
    fn default() -> Self {
        Self {
            min_mult: 4,
            min_alt: 4,
        }
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

/// Data about alleles
#[derive(Clone, Default)]
pub struct AlleleData {
    pub alt_refs: Vec<AltRef>,
}

/// Join algorithmic options.
#[derive(PartialEq)]
pub struct JoinAlgOpt {
    /// max score for join
    pub max_score: f64,
    pub max_cdr3_diffs: usize,
    pub join_cdr3_ident: f64,
    pub fwr1_cdr12_delta: f64,
    pub cdr3_normal_len: usize,
    pub auto_share: usize,
    pub comp_filt: usize,
    pub comp_filt_bound: usize,
}

impl Default for JoinAlgOpt {
    fn default() -> Self {
        Self {
            max_score: 100_000.0,
            max_cdr3_diffs: 1000,
            join_cdr3_ident: 85.0,
            fwr1_cdr12_delta: 20.0,
            cdr3_normal_len: 42,
            auto_share: 15,
            comp_filt: 8,
            comp_filt_bound: 80,
        }
    }
}

/// Set up control datastructure (EncloneControl).  This is stuff that is constant for a given
/// run of enclone.  If you add something to this, be sure to update the "changed" section in
/// enclone_server.rs, if needed.
#[derive(Default)]
pub struct EncloneControl {
    /// miscellaneous general options
    pub gen_opt: GeneralOpt,
    /// Config options used by cellranger.
    pub cr_opt: CellrangerOpt,
    /// algorithmic heuristics
    pub heur: ClonotypeHeuristics,
    /// origin (sample) info
    pub origin_info: OriginInfo,
    /// algorithmic options for allele finding
    pub allele_alg_opt: AlleleAlgOpt,
    /// algorithmic options for join
    pub join_alg_opt: JoinAlgOpt,
}

/// Set up data structure to track clonotype data.  A TigData is for one contig;
/// a Vec<TigData> is for one barcode, and an ExactClonotype is for an exact subclonotype.
#[derive(Eq, Ord, PartialEq, PartialOrd, Clone)] // not sure these are all needed
pub struct TigData {
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
}

impl TigData {
    pub fn seq(&self) -> &[u8] {
        &self.full_seq[self.v_start..self.j_stop]
    }
}

/// The ExactClonotype data structure stores information that could be exhibited as a
/// Vec<Vec<Vec<TigData>>>, but it avoids repetition of identical data.
///
/// TigData0: data for each cell
/// TigData1: shared data
#[derive(Clone)]
pub struct TigData0 {
    /// quality scores, truncated to V..J
    pub quals: Vec<u8>,
    /// start of V on full contig sequence
    pub v_start: usize,
    /// stop of J on full contig sequence
    pub j_stop: usize,
    /// start of C on full contig sequence
    pub c_start: Option<usize>,
    /// full contig sequence
    pub full_seq: Vec<u8>,
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
    /// index of V segment reference sequence in ref file
    pub v_ref_id: usize,
}

#[derive(Clone, Default)]
pub struct Junction {
    /// junction alignment complexity
    pub hcomp: usize,
}

/// The ExactClonotype data structure stores information that could be exhibited as a
/// Vec<Vec<Vec<TigData>>>, but it avoids repetition of identical data.
///
/// TigData0: data for each cell
/// TigData1: shared data
#[derive(Clone)]
pub struct TigData1 {
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
    /// start of aligned D on full contig sequence
    pub d_start: Option<usize>,
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
    pub share: Vec<TigData1>,
    /// clone info, excluding shared stuff
    pub clones: Vec<Vec<TigData0>>,
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
