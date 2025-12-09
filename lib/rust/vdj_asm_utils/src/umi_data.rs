//! This captures and analyzes data about vdj umis.
#![deny(missing_docs)]

use barcode::Barcode;
use debruijn::dna_string::DnaString;
use serde::{Deserialize, Serialize};
use umi::Umi;
use vdj_types::VdjChain;

/// Strongly-typed foreign key index into a vector of contigs
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ContigID(pub usize);

impl ContigID {
    /// Return a string name for this contig in the context of the provided barcode
    pub fn tigname(self, barcode: &String) -> String {
        format!("{}_contig_{}", barcode, self.0)
    }
}

#[expect(missing_docs)]
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct UmiInfo {
    pub barcode: Barcode,
    pub umi: Umi,
    pub umi_id: i32,
    pub contig: Option<ContigID>,
    pub chain_type: Option<VdjChain>,
    /// Has more than 1 non-empty read
    pub nonsolo: bool,
    /// Survives assembly graph
    pub surviving: bool,
    /// Supports contig junction sequence
    pub junction_support: bool,
    pub reads: usize,
    pub read_pairs: usize,
    pub nonempty_reads: usize,
}

impl UmiInfo {
    /// Create a new UmiInfo object.
    pub fn new(
        umi: &String,
        umi_id: i32,
        barcode: &String,
        reads: Vec<&DnaString>,
        single_end: bool,
    ) -> UmiInfo {
        let umi = Umi::new(umi.as_bytes());
        let barcode = Barcode::from_bytes(barcode.as_bytes()).unwrap();
        let num_reads = reads.len();
        let read_pairs = if single_end { num_reads } else { num_reads / 2 };
        let num_nonempty_reads = reads
            .into_iter()
            .filter(|read| !DnaString::is_empty(read))
            .count();

        UmiInfo {
            barcode,
            umi,
            umi_id,
            reads: num_reads,
            nonempty_reads: num_nonempty_reads,
            read_pairs,
            contig: None,
            chain_type: None,
            nonsolo: false,
            surviving: false,
            junction_support: false,
        }
    }
}
