// Stage to quickly get the counts of UMIs and Reads per Barcode for valid
// cell barcodes for a given feature type (e.g. AB, CMO, GEX), etc.
//use strum_macros::EnumString;
use cr_h5::molecule_info::{
    BarcodeIdxType, FullUmiCount, LibraryIdxType, MoleculeInfoIterator, MoleculeInfoReader,
};
use cr_types::reference::feature_reference::{FeatureDef, FeatureType};
use ndarray::Axis;
use pyo3::prelude::*;
use std::collections::{HashMap, HashSet};
use std::iter::zip;
use std::path::Path;
use std::str::FromStr;

//  A struct that can be used to filter out barcodes from a Molecule Info that are not associated
// with a spot and/or cell
pub(crate) struct FilteredBarcodeFilter {
    filtered_barcodes: HashSet<(BarcodeIdxType, LibraryIdxType)>,
    // last valid bc_lib_group
    last_valid_bc_lib: (BarcodeIdxType, LibraryIdxType),
    // last invalid bc_lib_group
    last_invalid_bc_lib: (BarcodeIdxType, LibraryIdxType),
}

impl FilteredBarcodeFilter {
    pub(crate) fn new(path: &Path) -> Self {
        // Get Cell-Associated Barcodes
        let bcs = MoleculeInfoReader::read_barcode_info(path)
            .expect("Could not read barcode info")
            .0;
        let mut filtered_barcodes = HashSet::with_capacity(bcs.nrows());
        for cell_bc in bcs.axis_iter(Axis(0)) {
            filtered_barcodes.insert((cell_bc[0], cell_bc[1] as u16));
        }
        FilteredBarcodeFilter {
            filtered_barcodes,
            last_valid_bc_lib: (BarcodeIdxType::MAX, LibraryIdxType::MAX),
            last_invalid_bc_lib: (BarcodeIdxType::MAX, LibraryIdxType::MAX),
        }
    }
    // Returns true if in the filtered barcodes
    #[allow(clippy::wrong_self_convention)]
    pub(crate) fn is_cell_or_spot(&mut self, bc_lib: (BarcodeIdxType, LibraryIdxType)) -> bool {
        let is_cell: bool;
        if bc_lib == self.last_valid_bc_lib {
            is_cell = true;
        } else if bc_lib == self.last_invalid_bc_lib {
            is_cell = false;
        } else {
            is_cell = self.filtered_barcodes.contains(&bc_lib);
            if is_cell {
                self.last_valid_bc_lib = bc_lib;
            } else {
                self.last_invalid_bc_lib = bc_lib;
            }
        }
        is_cell
    }
}

/// Data on counts per barcode for UMI and Reads
#[derive(Clone)]
struct BarcodeCounts {
    umi_cnts: Vec<u32>,
    read_cnts: Vec<u32>,
}

impl BarcodeCounts {
    fn new(num_barcodes: BarcodeIdxType) -> Self {
        BarcodeCounts {
            umi_cnts: vec![0; num_barcodes as usize],
            read_cnts: vec![0; num_barcodes as usize],
        }
    }
}

/// A struct to accumulate data for the given feature type, for each
/// feature of interest, we keep two arrays of size num_barcodes
/// to count how many reads and umis each barcode had
struct FeatureCountsPerBarcode {
    // first cmo feature index
    idx_start: usize,
    // last cmo feature index
    idx_end: usize,
    // name of all features in order
    features: Vec<FeatureDef>,
    // counts of umis and reads per barcode
    cnts: Vec<BarcodeCounts>,

    // Cell associated barcode/gem groups
    cell_associated_filter: FilteredBarcodeFilter,
}

impl FeatureCountsPerBarcode {
    fn new(cur_path: &Path, features: Vec<FeatureDef>) -> Self {
        assert!(
            !features.is_empty(),
            "No features of this type were found in the FeatureRef."
        );
        // Assuming all tags are in
        let mut index = features[0].index;
        for feature in &features[1..] {
            index += 1;
            assert!(
                index == feature.index,
                "The feature entries in the FeatureRef were not in serial order"
            );
        }
        // Get valid cell-barcodes
        let num_barcodes = MoleculeInfoReader::read_barcodes(cur_path).unwrap().len();
        let cell_associated_filter = FilteredBarcodeFilter::new(cur_path);
        // Create the vectors to hold the barcode read/umi counts
        let size = index - features[0].index + 1;
        let cnts = vec![BarcodeCounts::new(num_barcodes as BarcodeIdxType); size];

        FeatureCountsPerBarcode {
            idx_start: features[0].index,
            idx_end: features.last().unwrap().index,
            features,
            cnts,
            cell_associated_filter,
        }
    }
    fn consume_umicount(&mut self, cnt: FullUmiCount) {
        let idx = cnt.umi_data.feature_idx as usize;
        if idx >= self.idx_start && idx <= self.idx_end {
            let bc_lib = (cnt.barcode_idx, cnt.umi_data.library_idx);
            let is_cell = self.cell_associated_filter.is_cell_or_spot(bc_lib);
            if is_cell {
                let n_reads = cnt.umi_data.read_count;
                let barcode_idx = cnt.barcode_idx as usize;
                let index = idx - self.idx_start;
                let cnts = &mut self.cnts[index];
                cnts.umi_cnts[barcode_idx] += 1;
                cnts.read_cnts[barcode_idx] += n_reads;
            }
        }
    }
}

/// Method to count the number of UMIs and Reads associated with each cell-associated barcode
/// for a given library type.  It returns a dictionary with a string key of feature_ids and a
/// value of a tuple with vectors of reads per barcode, and umis per barcode.  Prior to returning, any
/// cell-associated barcodes with 0 counts are removed, and it is in general not guaranteed what
/// index in the count arrays corresponds to a given barcode.
#[allow(clippy::type_complexity)]
#[pyfunction]
pub(crate) fn counts_per_barcode(
    _py: Python<'_>,
    mol_info_path: String,
    feature_type_str: String,
) -> PyResult<HashMap<String, (Vec<u32>, Vec<u32>)>> {
    // Sum up all the read counts for a list of features, returning a tuple of
    // cnts, feature_id, feature_type
    // method is inefficient because it parses all data columns in mol_info but only uses 2,
    // also is efficient because it doesn't use a lot of memory of once and streams the mol_info
    let cur_path = Path::new(&mol_info_path);
    let iterator = MoleculeInfoIterator::new(cur_path).unwrap();
    let feature_type = FeatureType::from_str(&feature_type_str).unwrap();
    let features: Vec<FeatureDef> = iterator
        .feature_ref
        .feature_defs
        .iter()
        .filter(|&x| x.feature_type == feature_type)
        .cloned()
        .collect();

    let mut counter = FeatureCountsPerBarcode::new(cur_path, features);

    // Count up Per-barcode types
    iterator.for_each(|x| {
        counter.consume_umicount(x);
    });

    // Now translate back for Python
    Ok(zip(
        counter.features.into_iter().map(|feat| feat.id),
        counter.cnts,
    )
    .map(|(feature, mut cnt)| {
        cnt.read_cnts.retain(|&x| x != 0);
        cnt.umi_cnts.retain(|&x| x != 0);
        (feature, (cnt.umi_cnts, cnt.read_cnts))
    })
    .collect())
}
