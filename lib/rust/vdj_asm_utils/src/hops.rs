// Data structure for storing information about cross-library contamination on a
// flowcell, created by find_hops.rs.  Such contamination might arise from Illumina
// index hopping or direct contamination between libraries e.g. from adjacent wells.
//
// On a flowcell, there may be multiple lena ids ("a group") that refer to the same
// sample indices, e.g. at different levels of downsampling.  For purposes here,
// these are treated as synonymous, and tracked by a data structure lena_index.
//
// This data structure has entries (barcode,umi,lid,mult) where lid is a lena
// group index and mult is the number of read pairs for the given barcode/umi/lid.
//
// Singleton entries for given (barcode,umi,lid) are not shown, nor are entries
// for invalid barcodes, although these conditions are not enforced by the data
// structure.  Barcodes and umis containing Ns are intended to be ignored.
//
// Packing of barcode and umi: see pack_dna.rs.
//
// ISSUES
//
// ◼ Has hardcoded barcode and umi lengths.
//
// ◼ Depending on the lengths of the barcode and umi, it could be slightly more
//   efficient to store the barcode + umi as a single packed entry.
//
// ◼ If there is ever a case where the sample indices for two lena ids overlap but
//   are not equal, this will create a mess.

use io_utils::{fwriteln, open_for_write_new};
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use std::io::Write;
use string_utils::strme;
use tenkit2::pack_dna::{unpack_bases_10, unpack_bases_16};

#[derive(Serialize, Deserialize)]
pub struct FlowcellCrossEntry {
    pub bc: [u8; 4],  // packed barcode, assumed to be 16 bases
    pub umi: [u8; 3], // packed umi, assumed to be 10 bases
    pub lid: u16,     // lena group index (unsafe to store as u8)
    pub mult: u32,    // multiplicity
}

#[derive(Serialize, Deserialize)]
pub struct FlowcellContam {
    pub lena_index: Vec<Vec<usize>>, // maps index to sorted group of lena ids
    pub entry: Vec<FlowcellCrossEntry>, // entries
}

impl FlowcellContam {
    pub fn new() -> FlowcellContam {
        FlowcellContam {
            lena_index: Vec::new(),
            entry: Vec::new(),
        }
    }

    #[allow(clippy::many_single_char_names)]
    pub fn print(&self, out: &str) {
        let mut f = open_for_write_new![&out];
        let (mut i, mut xcount) = (0, 0);
        while i < self.entry.len() {
            let mut j = i + 1;
            while j < self.entry.len() {
                if self.entry[j].bc != self.entry[i].bc || self.entry[j].umi != self.entry[i].umi {
                    break;
                }
                j += 1;
            }
            let (mut b, mut u) = ([0_u8; 16], [0_u8; 10]);
            unpack_bases_16(&self.entry[i].bc, &mut b);
            unpack_bases_10(&self.entry[i].umi, &mut u);
            xcount += 1;
            fwriteln!(f, "\n[{}] b = {}, u = {}", xcount, strme(&b), strme(&u));
            for k in i..j {
                let lenas = &self.lena_index[self.entry[k].lid as usize];
                fwriteln!(
                    f,
                    "mult = {}, lenas = {}",
                    self.entry[k].mult,
                    lenas.iter().format(",")
                );
            }
            i = j;
        }
    }
}

impl Default for FlowcellContam {
    fn default() -> Self {
        Self::new()
    }
}
