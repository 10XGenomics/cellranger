// Various functions to iterate over molecule info files
#![deny(missing_docs)]

// Code for downsampling and collecting filters, to move a downsampler out of pandas and into Rust

use cr_h5::molecule_info::{
    BarcodeIdxType, FeatureIdxType, FullUmiCount, LibraryIdxType, MoleculeInfoIterator,
    MoleculeInfoReader, PerLibrarySubSampler,
};
use cr_types::{LibraryInfo, LibraryType};
use itertools::Itertools;
use pyo3::prelude::*;
use std::collections::{HashMap, HashSet};
use std::path::Path;

/// Read and downsample a molecule_info.h5 file, and produce a feature-barcode matrix.
/// Exclude barcodes that are not in the filtered barcodes dataset /barcode_info/pass_filter.
#[pyfunction]
pub(super) fn downsample_molinfo(
    _py: Python<'_>,
    mol_info_path: String,
    subsample_rates: Vec<f64>,
) -> (Vec<FeatureIdxType>, Vec<BarcodeIdxType>, Vec<u32>) {
    let mut downsampler = PerLibrarySubSampler::new(subsample_rates, 0);
    let (features, barcodes, counts) = MoleculeInfoIterator::new(Path::new(&mol_info_path))
        .unwrap()
        .cell_barcodes_only(true)
        .unwrap()
        .map(|x| downsampler.down_sample(x))
        .filter_map(|x| {
            (x.umi_data.read_count >= 1).then_some((x.umi_data.feature_idx, x.barcode_idx))
        })
        .dedup_with_count()
        .map(|(count, (feature, barcode))| {
            let count: u32 = count.try_into().unwrap();
            (feature, barcode, count)
        })
        .multiunzip();
    (features, barcodes, counts)
}

struct TotalReadsCounter {
    total_reads: usize,
    total_reads_in_cells: usize,
    cell_barcodes: HashSet<BarcodeIdxType>,
}

impl TotalReadsCounter {
    fn observe_umi(&mut self, umi: &FullUmiCount) {
        self.total_reads += umi.umi_data.read_count as usize;
        if self.cell_barcodes.contains(&umi.barcode_idx) {
            self.total_reads_in_cells += umi.umi_data.read_count as usize;
        }
    }
}

/// Count the reads and reads associated with cell barcodes.
#[pyfunction]
pub(super) fn count_reads_and_reads_in_cells(
    _py: Python<'_>,
    mol_info_path: String,
    library_idxs: Vec<LibraryIdxType>,
    cell_barcodes: Vec<BarcodeIdxType>,
) -> PyResult<(usize, usize)> {
    let mut reporter = TotalReadsCounter {
        total_reads: 0,
        total_reads_in_cells: 0,
        cell_barcodes: cell_barcodes.into_iter().collect(),
    };

    for umi in MoleculeInfoIterator::new(Path::new(&mol_info_path))
        .unwrap()
        .filter(|umi| library_idxs.contains(&umi.umi_data.library_idx))
    {
        reporter.observe_umi(&umi);
    }

    Ok((reporter.total_reads, reporter.total_reads_in_cells))
}

/// Count on-target usable reads for each library. GEX reads must be on-target.
#[pyfunction]
pub(super) fn count_usable_reads(_py: Python<'_>, mol_info_path: String) -> HashMap<u16, usize> {
    let mol_info_path = Path::new(&mol_info_path);
    let library_types: HashMap<_, _> = MoleculeInfoReader::read_library_info(mol_info_path)
        .unwrap()
        .into_iter()
        .map(|x| match x {
            LibraryInfo::Count(x) => (x.library_id, x.library_type),
            LibraryInfo::Aggr(_) => unreachable!(),
        })
        .collect();

    let molinfo = MoleculeInfoIterator::new(mol_info_path).unwrap();
    let gex_target_set = &molinfo.feature_ref.target_set.clone();
    molinfo
        .cell_barcodes_only(true)
        .unwrap()
        .map(|x| x.umi_data)
        .filter(|x| {
            if library_types[&x.library_idx] != LibraryType::Gex {
                return true;
            };
            let Some(gex_target_set) = gex_target_set else {
                return true;
            };
            gex_target_set.is_on_target(x.feature_idx)
        })
        .map(|x| (x.library_idx, x.read_count as usize))
        .into_grouping_map()
        .sum()
}

pub(super) struct LibraryAndFeatureIdxFilter {
    pub(super) library_idxs: Vec<LibraryIdxType>,
    pub(super) feature_idxs: HashSet<FeatureIdxType>,
}

impl LibraryAndFeatureIdxFilter {
    fn passes_filter(&self, umi: &FullUmiCount) -> bool {
        self.library_idxs.contains(&umi.umi_data.library_idx)
            && self.feature_idxs.contains(&umi.umi_data.feature_idx)
    }
}

#[pyfunction]
pub(super) fn get_num_umis_per_barcode(
    _py: Python<'_>,
    mol_info_path: String,
    library_idxs: Vec<LibraryIdxType>,
    target_feature_indices: Vec<FeatureIdxType>,
) -> PyResult<Vec<usize>> {
    let cur_path = Path::new(&mol_info_path);
    let filter = LibraryAndFeatureIdxFilter {
        library_idxs,
        feature_idxs: target_feature_indices.into_iter().collect(),
    };
    let counts = MoleculeInfoIterator::new(cur_path)
        .expect("Could not open molecule info")
        .filter_map(|x| filter.passes_filter(&x).then_some(x.barcode_idx))
        .dedup_with_count()
        .map(|(count, _)| count)
        .sorted()
        .collect();
    Ok(counts)
}
