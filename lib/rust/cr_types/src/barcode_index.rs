///
/// Types for dealing with 10x barcodes
///
use crate::BcCountFormat;
use anyhow::Result;
use barcode::Barcode;
use itertools::Itertools;
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::FileTypeRead;
use metric::{TxHashMap, TxHashSet};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BarcodeIndex {
    sorted_barcodes: Vec<Barcode>,
    index: TxHashMap<Barcode, usize>, // Perhaps binary search is enough?
}

impl BarcodeIndex {
    pub fn new(
        bc_count_file: &BcCountFormat,
        barcodes_under_tissue: Option<&JsonFile<Vec<String>>>,
    ) -> Result<Self> {
        let barcodes = bc_count_file
            .read()?
            .into_values()
            .flat_map(|h| h.into_iter().map(|(barcode, _)| barcode));
        Ok(if let Some(json_file) = barcodes_under_tissue {
            let additional_barcodes = json_file.read()?;
            barcodes
                .chain(additional_barcodes.into_iter().map(|x| x.parse().unwrap()))
                .collect()
        } else {
            barcodes.collect()
        })
    }
}

impl FromIterator<Barcode> for BarcodeIndex {
    /// Construct a barcode index from an iterator of barcodes.
    /// The iterator may contain duplicate barcodes, which will be ignored.
    fn from_iter<I: IntoIterator<Item = Barcode>>(barcodes: I) -> Self {
        let sorted_barcodes: Vec<Barcode> = barcodes.into_iter().sorted().dedup().collect();
        let index = sorted_barcodes
            .iter()
            .enumerate()
            .map(|(i, barcode)| (*barcode, i))
            .collect();
        BarcodeIndex {
            index,
            sorted_barcodes,
        }
    }
}

impl BarcodeIndex {
    /// Return the numerical index of the specified barcode.
    pub fn get_index(&self, barcode: &Barcode) -> usize {
        self.index[barcode]
    }

    /// Return the slice of sorted barcodes.
    pub fn sorted_barcodes(&self) -> &[Barcode] {
        &self.sorted_barcodes
    }

    /// Return an indicator vector corresponding to a set of barcodes
    /// A vector of bools returned where index i is true if
    /// sorted_barcodes[i] is in the set. Else false.
    /// If the set contains barcodes not in the BarcodeIndex, this is
    /// not caught
    pub fn into_indicator_vec(&self, filter_set: &TxHashSet<Barcode>) -> Vec<bool> {
        self.sorted_barcodes
            .iter()
            .map(|x| filter_set.contains(x))
            .collect()
    }

    /// Return the sorted barcodes.
    pub fn into_sorted_barcodes(self) -> Vec<Barcode> {
        self.sorted_barcodes
    }

    /// Return whether this index contains no barcodes.
    pub fn is_empty(&self) -> bool {
        self.sorted_barcodes.is_empty()
    }

    /// Return the number of barcodes.
    pub fn len(&self) -> usize {
        self.sorted_barcodes.len()
    }

    /// Checks if a barcode is in the barcode index
    pub fn contains_barcode(&self, barcode: &Barcode) -> bool {
        self.index.contains_key(barcode)
    }
}

#[cfg(test)]

mod tests {
    use super::*;
    use barcode::BcSeq;

    fn barcode(seq: &[u8]) -> Barcode {
        Barcode::with_seq(1, BcSeq::from_bytes(seq), true)
    }

    #[test]
    fn test_into_indicator_vec() {
        // Not in sorted order. BC0 and 2 are the same
        // In order it is BC3, BC0=BC2, BC1
        let vec_of_barcodes = vec![
            barcode(b"TTGAGGAGTTAGTGAGAACGCCGA"),
            barcode(b"TTGAGGTAGGGATGAAAACGCCGA"),
            barcode(b"TTGAGGAGTTAGTGAGAACGCCGA"),
            barcode(b"TTGAGCGAGTCAACTTAACGCCGA"),
        ];

        // These are BC3, BC1, BCnotInPreviousVec
        // In BCIndex generated this would be BC0, BC2, BCNotInIndex
        let query_vec_of_barcodes = vec![
            barcode(b"TTGAGCGAGTCAACTTAACGCCGA"),
            barcode(b"TTGAGGTAGGGATGAAAACGCCGA"),
            barcode(b"TTGAGGTAGGGATGATAACGCCGA"),
        ];
        let b_index = BarcodeIndex::from_iter(vec_of_barcodes);
        let set_of_barcodes = TxHashSet::from_iter(query_vec_of_barcodes);
        assert_eq!(
            b_index.into_indicator_vec(&set_of_barcodes),
            vec![true, false, true]
        );
    }
}
