// Various functions to iterate over molecule info files

// Code for downsampling and collecting filters, to move a downsampler out of pandas and into Rust

use cr_h5::molecule_info::{
    BarcodeIdxType, FeatureIdxType, FullUmiCount, GemGroupType, LibraryIdxType,
    MoleculeInfoIterator, MoleculeInfoReader, MoleculeInfoWriter, PerLibrarySubSampler,
};
use cr_types::{LibraryInfo, LibraryType};
use itertools::Itertools;
use numpy::PyArray1;
use pyo3::exceptions::PyIOError;
use pyo3::prelude::*;
use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};

/// Read and downsample a molecule_info.h5 file, and produce a feature-barcode matrix.
/// Exclude barcodes that are not in the filtered barcodes dataset /barcode_info/pass_filter.
#[pyfunction]
pub(crate) fn downsample_molinfo(
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
pub(crate) fn count_reads_and_reads_in_cells(
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
pub(crate) fn count_usable_reads(_py: Python<'_>, mol_info_path: String) -> HashMap<u16, usize> {
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

/// Method to concatenate molecule. Takes in a molecule info file as a path where all the top
/// level data (feature_ref, barcodes, metrics, library_info, etc). has been set (by python code in
/// MERGE_MOLECULES) and proceeds to append the mol_info_columns to the input file by concatentating
/// and modifying the inputs.
/// bc_idx_offsets: a vector of values to offset each barcode index by
/// gg_map: a vector of arrays where the index is the old gem group, and the value at that index the
/// new
/// lib_idx_maps: same as gg_map but for library indices
#[pyfunction]
pub(crate) fn concatenate_molecule_infos(
    _py: Python<'_>,
    out_path: String,
    sources: Vec<PathBuf>,
    bc_idx_offsets: Vec<BarcodeIdxType>,
    gg_map: Vec<&PyArray1<GemGroupType>>,
    lib_idx_maps: Vec<&PyArray1<LibraryIdxType>>,
) -> PyResult<()> {
    let out_path = Path::new(&out_path);
    let dest_result = MoleculeInfoWriter::from_file(out_path);
    match dest_result {
        Ok(mut dest) => {
            let gg_maps: Vec<_> = gg_map.into_iter().map(|gg| gg.to_vec().unwrap()).collect();
            let lib_idx_maps: Vec<_> = lib_idx_maps
                .into_iter()
                .map(|lib| lib.to_vec().unwrap())
                .collect();
            dest.concatenate_many(sources, bc_idx_offsets, gg_maps, lib_idx_maps)
                .expect("Error concatenating molecule infos");

            dest.flush().expect("Could not flush buffers");
            Ok(())
        }
        Err(_e) => Err(PyIOError::new_err("Could not open output file")),
    }
}

pub(crate) struct LibraryAndFeatureIdxFilter {
    pub(crate) library_idxs: Vec<LibraryIdxType>,
    pub(crate) feature_idxs: HashSet<FeatureIdxType>,
}

impl LibraryAndFeatureIdxFilter {
    fn passes_filter(&self, umi: &FullUmiCount) -> bool {
        self.library_idxs.contains(&umi.umi_data.library_idx)
            && self.feature_idxs.contains(&umi.umi_data.feature_idx)
    }
}

#[pyfunction]
pub(crate) fn get_num_umis_per_barcode(
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
