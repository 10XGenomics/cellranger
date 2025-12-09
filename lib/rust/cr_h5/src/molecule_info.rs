#![expect(missing_docs)]
use crate::iter::H5Iterator;
use crate::{
    ChunkedWriter, ColumnAction, feature_reference_io, probe_reference_io, scalar_attribute,
};
use anyhow::{Context, Result, bail};
use barcode::{Barcode, BarcodeContent, MAX_BARCODE_LENGTH};
use cr_types::chemistry::BarcodeReadComponent;
use cr_types::probe_set::Probe;
use cr_types::reference::feature_reference::{FeatureDef, FeatureReference};
use cr_types::{
    BarcodeIndex, BcUmiInfo, GenomeName, LibraryInfo, LibraryType, PROBE_IDX_SENTINEL_VALUE,
    UmiCount,
};
use hdf5::types::{FixedAscii, TypeDescriptor, VarLenAscii, VarLenUnicode};
use hdf5::{Dataset, Extents, File, Group};
use itertools::{Itertools, izip, process_results};
use metric::{TxHashMap, TxHashSet};
use ndarray::{Array1, Array2, Axis, s};
use rand::SeedableRng;
use rand::rngs::SmallRng;
use rand_distr::{Binomial, Distribution};
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::cmp::min;
use std::collections::BTreeSet;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::str::FromStr;
use umi::UmiType;
pub use visium_hd::bin_barcodes;

const H5_FILETYPE_KEY: &str = "filetype";
const METRICS_JSON: &str = "metrics_json";
const MOLECULE_H5_FILETYPE: &str = "molecule";
const FILE_VERSION_KEY: &str = "file_version";
const BARCODE_INFO_GROUP_NAME: &str = "barcode_info";
const PASS_FILTER_DATASET_NAME: &str = "pass_filter";
const GENOMES_DATASET_NAME: &str = "genomes";
const BARCODE_DATASET_NAME: &str = "barcodes";
const PROBE_GROUP_NAME: &str = "probes";
const LIBRARY_INFO: &str = "library_info";

const BARCODE_IDX_COL_NAME: &str = "barcode_idx";
const LIBRARY_IDX_COL_NAME: &str = "library_idx";
const FEATURE_IDX_COL_NAME: &str = "feature_idx";
const COUNT_COL_NAME: &str = "count";
const GEM_GROUP_COL_NAME: &str = "gem_group";
const UMI_COL_NAME: &str = "umi";
const UMI_TYPE_COL_NAME: &str = "umi_type";
const PROBE_IDX_COL_NAME: &str = "probe_idx";

const LIBRARY_METRICS_JSON: &str = "libraries";
const RAW_READS_IN_LIBRARY_METRICS_JSON: &str = "raw_read_pairs";

// const MOL_INFO_COLUMNS: [&'static str; 7] = [BARCODE_IDX_COL_NAME, LIBRARY_IDX_COL_NAME,
// FEATURE_IDX_COL_NAME, COUNT_COL_NAME, GEM_GROUP_COL_NAME, UMI_COL_NAME, UMI_TYPE_COL_NAME];

// Version history of molecule_info.h5 files.  We will only support parsing >=3
//
// Version 6:
//    adds `probe_idx` int32 dataset and `probes` group to annotate RTL reads with the probe
//    they came from.
//
// Version 5: PR #3815
//    adds `umi_type` uint32 dataset to distinguish between transcriptomic and non-transcriptomic
//    umi's in intron mode. We only use 1 bit of information and the remaining 31 bits will be used
//    in the future to further annotate umis.
//
// Version 4: PR #2466
//   revises how metrics required for aggr are stored in the molecule info.
//   They used to be stored as python pickle strings, which is difficult to interoperate with
//   outside python. Change to storing all the extra metrics in a json string.
//
// Version 3: PR #774
//     Columns about read mapping read mapping removed
//        * nonconf_mapped_reads
//        * unmapped_reads
//        * umi_corrected_reads
//        * barcode_corrected_reads
//        * conf_mapped_unique_read_pos
//     Barcodes changed from 2 bit encoded `barcode` into `barcode_idx`
//     `feaure_idx` and feature ref added, 'gene', 'genome' removed
//     `library_idx` added
//     `reads` renamed `count`
//     `barcode_info` added
//     tables was switched with h5py for writing
//
//
//  Version 2: commit 1feb9c0ccf0efeb19729ea954ac8f58876c1ee9d
//     First time a version key is used, everything before considered V1 by absence of key
//
const CURRENT_VERSION: i64 = 6;

type FA256 = FixedAscii<256>;
type FALibraryInfo = FixedAscii<{ 256 * 1024 }>;
type FABc = FixedAscii<MAX_BARCODE_LENGTH>;
pub type BarcodeIdxType = u64;
pub type FeatureIdxType = u32;
pub type GemGroupType = u16;
pub type LibraryIdxType = u16;
pub type ProbeIdxType = i32;
type UmiTypeType = u32;
// Capital T to distinguish names from other variables
type UmiTypeT = u32;
type CountType = u32;

/// Read umi count data from a molecule_info.h5 file.
pub struct MoleculeInfoReader;

pub struct MolInfoDatasets {
    barcode_idx: Dataset,
    feature_idx: Dataset,
    probe_idx: Option<Dataset>,
    gem_group: Dataset,
    library_idx: Dataset,
    umi: Dataset,
    count: Dataset,
    /// Added in molecule_info.h5 version 5.
    umi_type: Option<Dataset>,
}

/// A single UMI count with minimal BC information (bc index and gem group)
/// used when we don't want to pass around a bunch of barcode sequences
#[derive(Clone, Copy, PartialOrd, Ord, PartialEq, Eq, Debug, Serialize, Deserialize)]
pub struct FullUmiCount {
    pub barcode_idx: BarcodeIdxType,
    pub gem_group: GemGroupType,
    pub umi_data: UmiCount,
}

fn binomial_sample(umi: &FullUmiCount, sample_rate: f64, rng: &mut SmallRng) -> u32 {
    let count = umi.umi_data.read_count as u64;
    if sample_rate == 1.0 {
        count as u32
    } else {
        let bin = Binomial::new(count, sample_rate).unwrap();
        let new_count = bin.sample(rng);
        new_count as u32
    }
}

pub struct PerLibrarySubSampler {
    rate_per_lib: Vec<f64>,
    rng: SmallRng,
}

impl PerLibrarySubSampler {
    pub fn new(rate_per_lib: Vec<f64>, seed: u64) -> Self {
        Self {
            rate_per_lib,
            rng: SmallRng::seed_from_u64(seed),
        }
    }

    pub fn down_sample(&mut self, mut umi: FullUmiCount) -> FullUmiCount {
        let lib_index = umi.umi_data.library_idx as usize;
        let sample_rate = self.rate_per_lib[lib_index];
        // Get the number of molecules and make a Binomial draw using the subsample rate.
        // If we keep at least one molecule, record the barcode and feature indices
        let new_count = binomial_sample(&umi, sample_rate, &mut self.rng);
        umi.umi_data.read_count = new_count;
        umi
    }
}

pub struct UniformDownSampler {
    downsample_rate: f64,
    rng: SmallRng,
}

impl UniformDownSampler {
    pub fn new(downsample_rate: f64, seed: u64) -> Self {
        Self {
            downsample_rate,
            rng: SmallRng::seed_from_u64(seed),
        }
    }

    pub fn down_sample(&mut self, mut umi: FullUmiCount) -> FullUmiCount {
        let new_count = binomial_sample(&umi, self.downsample_rate, &mut self.rng);
        umi.umi_data.read_count = new_count;
        umi
    }
}

// Workaround for the rust hdf5 bindings not seeming to have a
// an easy way to iterate through the data.
pub struct MolInfoCache {
    // Start of slice in molecule info this section is from
    pub start: usize,
    // End of slice in molecule info this section is from
    pub end: usize,
    // Data from between start and end
    pub data: Vec<FullUmiCount>,
}

impl MolInfoCache {
    // Load a cache of umi counts
    fn new(mol_info_dsets: &MolInfoDatasets, start: usize, end: usize) -> Result<MolInfoCache> {
        let barcode_idx = mol_info_dsets.barcode_idx.read_slice_1d(start..end)?;
        let feature_idx = mol_info_dsets.feature_idx.read_slice_1d(start..end)?;
        let gem_group = mol_info_dsets.gem_group.read_slice_1d(start..end)?;
        let library_idx = mol_info_dsets.library_idx.read_slice_1d(start..end)?;
        let umi = mol_info_dsets.umi.read_slice_1d(start..end)?;
        let probe_idx: Array1<Option<ProbeIdxType>> =
            if let Some(probe_idx) = &mol_info_dsets.probe_idx {
                Array1::from_iter(probe_idx.read_slice_1d(start..end)?.map(|x| Some(*x)))
            } else {
                Array1::from_elem(umi.len(), None)
            };
        let count = mol_info_dsets.count.read_slice_1d(start..end)?;
        let umi_type = if let Some(umi_type) = &mol_info_dsets.umi_type {
            umi_type.read_slice_1d(start..end)?
        } else {
            Array1::from_elem(umi.len(), UmiType::default().to_u32())
        };

        let vec: Vec<FullUmiCount> = izip!(
            &barcode_idx,
            &gem_group,
            &library_idx,
            &feature_idx,
            &probe_idx,
            &umi,
            &count,
            &umi_type
        )
        .map(|x| FullUmiCount {
            barcode_idx: *x.0,
            gem_group: *x.1,
            umi_data: UmiCount {
                library_idx: *x.2,
                feature_idx: *x.3,
                probe_idx: *x.4,
                umi: *x.5,
                read_count: *x.6,
                utype: UmiType::from_u32(x.7),
            },
        })
        .collect();
        Ok(MolInfoCache {
            start,
            end,
            data: vec,
        })
    }
}

impl MolInfoDatasets {
    fn new(file: &File) -> Result<MolInfoDatasets> {
        let probe_idx = file.dataset(PROBE_IDX_COL_NAME).ok();
        let d = MolInfoDatasets {
            barcode_idx: file.dataset("barcode_idx")?,
            feature_idx: file.dataset("feature_idx")?,
            gem_group: file.dataset("gem_group")?,
            library_idx: file.dataset("library_idx")?,
            probe_idx,
            umi: file.dataset("umi")?,
            count: file.dataset("count")?,
            umi_type: file.dataset("umi_type").ok(),
        };
        Ok(d)
    }
}

const ITERATOR_CHUNK_SIZE: usize = 262144;
// Iterate through the file producing BcUmiInfo
pub struct MoleculeInfoIterator {
    pub feature_ref: FeatureReference,
    pub probes: Option<Vec<Probe>>,
    //Barcodes should probably be next, but not currently used
    datasets: MolInfoDatasets,
    cache: Option<MolInfoCache>,
    index: usize,
    // Total elements in the dset
    size: usize,
    path: PathBuf,
    // If not None, exclude molecules whose barcodes are not in the set
    barcode_ids: Option<TxHashSet<BarcodeIdxType>>,
    // If not None, exclude molecules whose features are not in the set
    feature_ids: Option<TxHashSet<FeatureIdxType>>,
}

fn check_version(f: &File) -> Result<i64> {
    let version = f.attr(FILE_VERSION_KEY)?;
    Ok(version.as_reader().read_scalar()?)
}

impl MoleculeInfoIterator {
    pub fn new(path: &Path) -> Result<MoleculeInfoIterator> {
        let file = File::open(path).with_context(|| path.display().to_string())?;
        let mol_info_dsets = MolInfoDatasets::new(&file)?;
        let size = mol_info_dsets.barcode_idx.size();
        let features = feature_reference_io::from_h5(&file.group("features")?)?;
        let probes = match &file.group(PROBE_GROUP_NAME) {
            Ok(value) => Some(probe_reference_io::from_h5(
                value,
                features
                    .feature_defs
                    .iter()
                    .map(|x| &x.genome)
                    .unique()
                    .cloned()
                    .collect(),
            )?),
            Err(_) => None,
        };
        Ok(MoleculeInfoIterator {
            feature_ref: features,
            probes,
            datasets: mol_info_dsets,
            cache: None,
            index: 0,
            size,
            path: path.to_path_buf(),
            barcode_ids: None,
            feature_ids: None,
        })
    }

    pub fn cell_barcodes_only(mut self, enable: bool) -> Result<Self> {
        if !enable {
            return Ok(self);
        }
        let barcode_ids = MoleculeInfoReader::read_filtered_barcode_ids(&self.path)?;
        self.barcode_ids = Some(barcode_ids);
        Ok(self)
    }

    pub fn filter_features<F: Fn(&FeatureDef) -> bool>(mut self, f: F) -> Result<Self> {
        let feature_ids: TxHashSet<_> = self
            .feature_ref
            .feature_defs
            .iter()
            .filter_map(|fdef| {
                if f(fdef) {
                    Some(fdef.index as FeatureIdxType)
                } else {
                    None
                }
            })
            .collect();
        self.feature_ids = Some(feature_ids);
        Ok(self)
    }

    fn allow_entry(&self, entry: &FullUmiCount) -> bool {
        self.barcode_ids
            .as_ref()
            .is_none_or(|ids| ids.contains(&entry.barcode_idx))
            && self
                .feature_ids
                .as_ref()
                .is_none_or(|ids| ids.contains(&entry.umi_data.feature_idx))
    }
}

impl Iterator for MoleculeInfoIterator {
    type Item = FullUmiCount;
    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if self.index >= self.size {
                break None;
            }

            // initialize
            if self.cache.is_none() {
                self.cache = Some(
                    MolInfoCache::new(
                        &self.datasets,
                        self.index,
                        min(self.index + ITERATOR_CHUNK_SIZE, self.size),
                    )
                    .unwrap(),
                );
            }
            match &self.cache {
                Some(x) => {
                    let cur_index = self.index % ITERATOR_CHUNK_SIZE;
                    let cur_data = x.data[cur_index];
                    self.index += 1;
                    // Clear the cache if we've read past it
                    if self.index >= x.end {
                        self.cache = None;
                    }
                    if self.allow_entry(&cur_data) {
                        break Some(cur_data);
                    }
                }
                None => panic!("Logical flaw in iteration"),
            }
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = self.size.saturating_sub(self.index);
        (remaining, Some(remaining))
    }
}

impl MoleculeInfoReader {
    /// Read a molecule_info.h5 file, yield the per-barcode umi info and the feature reference
    pub fn read(path: &Path) -> Result<(Vec<BcUmiInfo>, FeatureReference)> {
        let file = File::open(path)?;

        let mol_info_dsets = MolInfoDatasets::new(&file)?;

        let barcode_idx = mol_info_dsets.barcode_idx.read_1d::<BarcodeIdxType>()?;
        let feature_idx = mol_info_dsets.feature_idx.read_1d::<FeatureIdxType>()?;
        let gem_group = mol_info_dsets.gem_group.read_1d::<GemGroupType>()?;
        let library_idx = mol_info_dsets.library_idx.read_1d::<LibraryIdxType>()?;

        let probe_idx = match &mol_info_dsets.probe_idx {
            Some(probe_idx) => Some(probe_idx.read_1d::<ProbeIdxType>()?),
            None => None,
        };
        let umi = mol_info_dsets.umi.read_1d::<UmiTypeT>()?;
        let count = mol_info_dsets.count.read_1d::<CountType>()?;
        let umi_type = if let Some(umi_type) = &mol_info_dsets.umi_type {
            umi_type.read_1d()?
        } else {
            Array1::from_elem(umi.len(), UmiType::default().to_u32())
        };

        let barcode_seqs = file.dataset("barcodes")?.read_1d::<FABc>()?;
        let features = feature_reference_io::from_h5(&file.group("features")?)?;

        let mut cur_bc = barcode_idx[0];
        let mut cur_gg = gem_group[0];
        let mut start = 0;
        let mut i = 0;

        let mut all_bc_info = Vec::new();

        loop {
            if start >= barcode_idx.len() {
                break;
            }

            loop {
                if i >= barcode_idx.len() || (barcode_idx[i], gem_group[i]) != (cur_bc, cur_gg) {
                    break;
                }
                i += 1;
            }

            let end = i;
            let mut umi_counts = Vec::with_capacity(end - start);

            // construct the fastq_10x barcode
            let barcode = Barcode::with_content(
                gem_group[start],
                BarcodeContent::from_bytes(barcode_seqs[barcode_idx[start] as usize].as_bytes())?,
                true,
            );

            for j in start..end {
                let cur_probe_idx = probe_idx.as_ref().map(|indices| indices[j]);
                let c = UmiCount {
                    library_idx: library_idx[j],
                    feature_idx: feature_idx[j],
                    probe_idx: cur_probe_idx,
                    umi: umi[j],
                    read_count: count[j],
                    utype: UmiType::from_u32(&umi_type[j]),
                };

                assert_eq!(gem_group[j], barcode.gem_group());
                umi_counts.push(c);
            }

            let bc_umi_info = BcUmiInfo {
                barcode,
                umi_counts,
            };

            start = end;
            if start < barcode_idx.len() {
                cur_bc = barcode_idx[start];
                cur_gg = gem_group[start];
            }

            all_bc_info.push(bc_umi_info);
        }

        assert_eq!(all_bc_info.len(), barcode_idx.iter().unique().count());

        let total_umis: usize = all_bc_info.iter().map(|x| x.umi_counts.len()).sum();
        assert_eq!(barcode_idx.len(), total_umis);

        Ok((all_bc_info, features))
    }

    fn open(path: &Path) -> Result<File> {
        File::open(path).with_context(|| path.display().to_string())
    }

    /// Read the barcodes.
    pub fn read_barcodes(path: &Path) -> Result<impl Iterator<Item = Result<BarcodeContent>>> {
        Ok(H5Iterator::<FABc>::new(
            Self::open(path)?.dataset(BARCODE_DATASET_NAME)?,
            ITERATOR_CHUNK_SIZE,
        )
        .map(|barcode| BarcodeContent::from_bytes(barcode?.as_bytes())))
    }

    /// Return the number of barcodes.
    pub fn read_barcodes_size(path: &Path) -> Result<usize> {
        Ok(Self::open(path)?.dataset(BARCODE_DATASET_NAME)?.size())
    }

    /// Return the number of barcodes with non-zero counts.
    ///
    /// This requires streaming the barcode IDs.
    /// For a file that has been pre-trimmed to only include non-zero barcodes,
    /// this should be the same as the number of barcodes.
    pub fn count_non_zero_bcs(path: &Path) -> Result<usize> {
        process_results(
            H5Iterator::<BarcodeIdxType>::new(
                Self::open(path)?.dataset(BARCODE_IDX_COL_NAME)?,
                ITERATOR_CHUNK_SIZE,
            ),
            |bc_idx_iter| bc_idx_iter.unique().count(),
        )
    }

    pub fn read_gem_groups_size(path: &Path) -> Result<usize> {
        Ok(Self::open(path)?.dataset(GEM_GROUP_COL_NAME)?.size())
    }

    pub fn read_filtered_barcode_ids(path: &Path) -> Result<TxHashSet<u64>> {
        let barcode_info_pass_filter = MoleculeInfoReader::read_barcode_info_pass_filter(path)?;
        Ok(barcode_info_pass_filter
            .index_axis(Axis(1), 0)
            .into_iter()
            .copied()
            .collect())
    }

    pub fn read_feature_ref(path: &Path) -> Result<FeatureReference> {
        feature_reference_io::from_h5(&Self::open(path)?.group("features")?)
    }

    /// Read a molecule_info.h5 file and return the datasets barcode_info/pass_filter.
    pub fn read_barcode_info_pass_filter(path: &Path) -> Result<Array2<u64>> {
        Ok(Self::open(path)?
            .group(BARCODE_INFO_GROUP_NAME)?
            .dataset(PASS_FILTER_DATASET_NAME)?
            .read_2d::<u64>()?)
    }

    /// Read a molecule_info.h5 file and return the datasets barcode_info/genomes.
    pub fn read_barcode_info_genomes(path: &Path) -> Result<Vec<GenomeName>> {
        Ok(Self::open(path)?
            .group(BARCODE_INFO_GROUP_NAME)?
            .dataset(GENOMES_DATASET_NAME)?
            .read_1d::<FA256>()?
            .into_iter()
            .map(|x| GenomeName::from(x.as_str()))
            .collect())
    }

    /// Read a molecule_info.h5 file and return a vector of LibraryInfo structs
    pub fn read_library_info(path: &Path) -> Result<Vec<LibraryInfo>> {
        let dataset = Self::open(path)?.dataset(LIBRARY_INFO)?;
        assert_eq!(dataset.size(), 1);
        let dtype = dataset.dtype()?;
        Ok(match dtype.to_descriptor()? {
            TypeDescriptor::FixedAscii(_) => {
                serde_json::from_str(dataset.read_1d::<FALibraryInfo>()?[0].as_str())?
            }
            TypeDescriptor::VarLenAscii => {
                serde_json::from_str(dataset.read_1d::<VarLenAscii>()?[0].as_str())?
            }
            _ => bail!("unexpected library_info dtype: {dtype:?}"),
        })
    }

    /// Read a molecule_info.h5 file and return the molecule_info JSON string.
    pub fn read_metrics(path: &Path) -> Result<String> {
        let ds = Self::open(path)?.dataset(METRICS_JSON)?;
        match ds.dtype()?.to_descriptor()? {
            TypeDescriptor::VarLenAscii => Ok(ds.read_scalar::<VarLenAscii>()?.to_string()),
            TypeDescriptor::VarLenUnicode => Ok(ds.read_scalar::<VarLenUnicode>()?.to_string()),
            _ => bail!("Can't convert metrics to string"),
        }
    }

    /// Whether this molecule info corresponds to a visium HD sample
    pub fn is_visium_hd(path: &Path) -> Result<bool> {
        let Some(chemistry_barcode) = Value::from_str(&Self::read_metrics(path)?)?
            .as_object_mut()
            .unwrap()
            .remove("chemistry_barcode")
        else {
            return Ok(false);
        };
        let Ok(chemistry_barcode) =
            serde_json::from_value::<Vec<BarcodeReadComponent>>(chemistry_barcode)
        else {
            return Ok(false);
        };
        Ok(chemistry_barcode
            .iter()
            .all(|x| x.whitelist().slide_name().is_some()))
    }

    /// Read a molecule_info.h5 file and return the gem groups
    pub fn read_gem_groups(path: &Path) -> Result<Array1<GemGroupType>> {
        Ok(Self::open(path)?
            .dataset(GEM_GROUP_COL_NAME)?
            .read_1d::<GemGroupType>()?)
    }

    /// Return an iterator over gem groups.
    pub fn iter_gem_groups(path: &Path) -> Result<H5Iterator<GemGroupType>> {
        Ok(H5Iterator::new(
            Self::open(path)?.dataset(GEM_GROUP_COL_NAME)?,
            ITERATOR_CHUNK_SIZE,
        ))
    }

    pub fn nrows(path: &Path) -> Result<usize> {
        Ok(Self::open(path)?.dataset(GEM_GROUP_COL_NAME)?.size())
    }

    // Return the filtered barcodes of the specified library type.
    pub fn read_filtered_barcodes(path: &Path, library_type: LibraryType) -> Result<Vec<Barcode>> {
        let lib_to_gem_group: TxHashMap<LibraryIdxType, GemGroupType> =
            MoleculeInfoReader::read_library_info(path)?
                .into_iter()
                .filter_map(|x| match x {
                    LibraryInfo::Count(x) if x.library_type == library_type => {
                        Some((x.library_id, x.gem_group))
                    }
                    LibraryInfo::Count(_) => None,
                    LibraryInfo::Aggr(_) => unimplemented!(),
                })
                .collect();

        let filtered_barcode_indices = MoleculeInfoReader::read_filtered_barcode_ids(path)?;
        let index_to_barcode: TxHashMap<BarcodeIdxType, BarcodeContent> = H5Iterator::new(
            Self::open(path)?.dataset(BARCODE_DATASET_NAME)?,
            ITERATOR_CHUNK_SIZE,
        )
        .enumerate()
        .map(|(index, barcode): (usize, Result<FABc>)| {
            BarcodeContent::from_str(&barcode?).map(|barcode| (index as BarcodeIdxType, barcode))
        })
        .filter_ok(|(index, _barcode)| filtered_barcode_indices.contains(index))
        .try_collect()?;

        Ok(MoleculeInfoReader::read_barcode_info_pass_filter(path)?
            .outer_iter()
            .filter_map(|barcode_library| {
                let &[barcode_index, library_index, _genome_index] =
                    barcode_library.as_slice().unwrap()
                else {
                    unreachable!();
                };
                lib_to_gem_group
                    .get(&(library_index as LibraryIdxType))
                    .map(|&gem_group| {
                        Barcode::with_content(gem_group, index_to_barcode[&barcode_index], true)
                    })
            })
            .collect())
    }

    /// Returns the library IDs in count GEX libraries in the molecule info
    pub fn get_count_gex_library_ids(path: &Path) -> Result<TxHashSet<LibraryIdxType>> {
        Ok(Self::read_library_info(path)?
            .iter()
            .filter_map(|lib_info| {
                if let LibraryInfo::Count(rna_lib) = lib_info {
                    rna_lib.library_type.is_gex().then_some(rna_lib.library_id)
                } else {
                    None
                }
            })
            .collect())
    }

    /// Returns the total raw reads in count GEX libraries in the molecule info
    pub fn get_raw_reads_in_count_gex_libraries(path: &Path) -> Result<i64> {
        let gex_lib_ids = Self::get_count_gex_library_ids(path)?;
        let total_gex_reads = Value::from_str(&Self::read_metrics(path)?)?[LIBRARY_METRICS_JSON]
            .as_object()
            .unwrap()
            .iter()
            .filter(|(lib_id, _lib_data)| gex_lib_ids.contains(&lib_id.parse().unwrap()))
            .map(|(_lib_id, lib_metrics)| {
                lib_metrics[RAW_READS_IN_LIBRARY_METRICS_JSON]
                    .as_i64()
                    .unwrap()
            })
            .sum();
        Ok(total_gex_reads)
    }
}

/// Write per-bc UMI count data progressively to a molecule_info.h5 file
pub struct MoleculeInfoWriter {
    file: File,
    writers: ColumnWriters,
}

const MOL_INFO_BUFFER_SZ: usize = 1 << 20;

impl MoleculeInfoWriter {
    /// Create a molecule_info.h5 writer.
    pub fn new(
        path: &Path,
        feature_ref: &FeatureReference,
        probes: Option<&[Probe]>,
        filtered_probes: Option<&[bool]>,
        barcodes: impl IntoIterator<Item = BarcodeContent>,
        library_info: &[LibraryInfo],
    ) -> Result<MoleculeInfoWriter> {
        let f = File::create(path)?;

        let library_info_json = serde_json::to_string(library_info)?;
        if library_info_json.len() < FALibraryInfo::capacity() {
            // Use FixedAscii for compatibility with CR prior to version 10.
            f.new_dataset_builder()
                .with_data(&[FALibraryInfo::from_ascii(&library_info_json)?])
                .create(LIBRARY_INFO)?;
        } else {
            f.new_dataset_builder()
                .with_data(&[VarLenAscii::from_ascii(&library_info_json)?])
                .create(LIBRARY_INFO)?;
        }

        // h5 headers
        scalar_attribute(
            &f,
            H5_FILETYPE_KEY,
            VarLenUnicode::from_str(MOLECULE_H5_FILETYPE)?,
        )?;

        // Version 5:
        // - adds `umi_type` uint32 dataset to distinguish between transcriptomic and non-transcriptomic
        //   umi's in intron mode. We only use 1 bit of information and the remaining 31 bits will be used
        //   in the future to further annotate umis.
        scalar_attribute(&f, FILE_VERSION_KEY, CURRENT_VERSION)?;

        Self::write_barcodes(&f, barcodes)?;

        feature_reference_io::to_h5(feature_ref, &mut f.create_group("features")?)?;
        if let Some(probes) = probes {
            let mut probe_group = f.create_group(PROBE_GROUP_NAME)?;
            probe_reference_io::to_h5(probes, filtered_probes.unwrap(), &mut probe_group)?;
        }

        Ok(MoleculeInfoWriter {
            writers: ColumnWriters::new(&f, probes.is_some(), ColumnAction::CreateNew)?,
            file: f,
        })
    }

    /// Write a metrics JSON string to /metrics_json dataset.
    pub fn write_metrics(&mut self, metrics: &str) -> Result<()> {
        self.file
            .new_dataset::<VarLenAscii>()
            .shape(())
            .create(METRICS_JSON)?
            .write_scalar(&VarLenAscii::from_ascii(metrics)?)?;
        Ok(())
    }

    /// write the /barcode_info/pass_filter and /barcode_info/genomes datsets
    /// to h5 file.
    pub fn write_barcode_info(
        &mut self,
        pass_filter: &Array2<u64>,
        genomes: &[GenomeName],
    ) -> Result<()> {
        let barcode_info_group = self.file.create_group(BARCODE_INFO_GROUP_NAME)?;
        write_barcode_info(&barcode_info_group, pass_filter, genomes)
    }

    /// write /barcodes to h5 file.
    fn write_barcodes(f: &File, barcodes: impl IntoIterator<Item = BarcodeContent>) -> Result<()> {
        // TODO: deduplicate with write_barcodes_column helper function in count_matrix

        // Buffer for writing out the fixed-width string representation of each barcode.
        let mut buf = Vec::with_capacity(FABc::capacity());
        let formatted_bc_iter = barcodes.into_iter().map(|bc| {
            write!(&mut buf, "{bc}")?;
            let formatted = FABc::from_ascii(&buf)?;
            buf.clear();
            anyhow::Ok(formatted)
        });
        formatted_bc_iter.process_results(|bc_iter| {
            ChunkedWriter::write_all(
                f,
                BARCODE_DATASET_NAME,
                MOL_INFO_BUFFER_SZ,
                ColumnAction::CreateNew,
                bc_iter,
            )
        })?
    }

    /// Write umi count data for one barcode to the molecule_info.h5 file.
    pub fn fill(&mut self, barcode_index: &BarcodeIndex, count_data: &BcUmiInfo) -> Result<()> {
        let barcode_idx = barcode_index.must_get(&count_data.barcode) as u64;
        for c in &count_data.umi_counts {
            self.write(count_data.barcode.gem_group(), barcode_idx, c)?;
        }
        Ok(())
    }

    /// Write data for a single molecule.
    pub fn write(&mut self, gem_group: u16, barcode_idx: u64, count: &UmiCount) -> Result<()> {
        self.writers.write(gem_group, barcode_idx, count)
    }

    /// Flush buffered data to disk. This method must be called when all data has been written.
    pub fn flush(&mut self) -> Result<()> {
        self.writers.flush()
    }

    /// Concatendate multiple molecule info files together.
    ///
    /// Remaps barcode indices, gem groups, and library IDs using the provided
    /// mappings
    pub fn concatenate_many(
        &mut self,
        sources: Vec<PathBuf>,
        bc_idx_offsets: Vec<BarcodeIdxType>,
        gg_maps: Vec<Vec<GemGroupType>>,
        lib_idx_maps: Vec<Vec<LibraryIdxType>>,
    ) -> Result<()> {
        for (idx, f) in sources.iter().enumerate() {
            let bc_offset = bc_idx_offsets[idx];
            let lib_map = &lib_idx_maps[idx];
            let gg_map = &gg_maps[idx];

            let src = MoleculeInfoIterator::new(f)?;

            for mut count in src {
                count.umi_data.library_idx = lib_map[usize::from(count.umi_data.library_idx)];
                self.write(
                    gg_map[usize::from(count.gem_group)],
                    count.barcode_idx + bc_offset,
                    &count.umi_data,
                )?;
            }
        }
        self.flush()
    }

    pub fn concatenate_metrics(
        metrics_list: Vec<serde_json::Map<String, Value>>,
    ) -> Result<String> {
        let mut combined_metrics: Option<serde_json::Map<String, Value>> = None;
        let mut gg_metrics: serde_json::Map<String, Value> = Default::default();
        let mut lib_metrics: serde_json::Map<String, Value> = Default::default();
        let mut targeted_metrics: Vec<serde_json::Map<String, Value>> = Default::default();

        for mut single_metrics in metrics_list {
            if combined_metrics.is_none() {
                combined_metrics = Some(
                    single_metrics
                        .iter()
                        .filter(|(k, _)| !k.starts_with("target"))
                        .map(|(k, v)| (k.clone(), v.clone()))
                        .collect(),
                );
            }

            // Collect the targeted metrics for each input.
            targeted_metrics.push(
                single_metrics
                    .iter()
                    .filter(|(k, _)| k.starts_with("target"))
                    .map(|(k, v)| (k.clone(), v.clone()))
                    .collect(),
            );

            // concatenate new gem groups to the metrics. if it collides with an existing
            // gem group, the old one will be overwritten.
            if let Some(single_gg_metrics) = single_metrics.remove("gem_groups") {
                let Value::Object(mut single_gg_metrics) = single_gg_metrics else {
                    unreachable!();
                };

                if let Some(analysis_parameters) = single_metrics.remove("analysis_parameters") {
                    let analysis_parameters = analysis_parameters.as_object().unwrap();
                    for metrics in single_gg_metrics.values_mut() {
                        metrics
                            .as_object_mut()
                            .unwrap()
                            .extend(analysis_parameters.clone());
                    }
                }
                gg_metrics.extend(single_gg_metrics);
            }

            lib_metrics.extend(single_metrics["libraries"].as_object().unwrap().clone());
        }

        let mut combined_metrics = combined_metrics.unwrap();
        combined_metrics.extend([
            ("is_aggregated".to_string(), Value::Bool(true)),
            ("gem_groups".to_string(), serde_json::to_value(gg_metrics)?),
            ("libraries".to_string(), serde_json::to_value(lib_metrics)?),
        ]);
        combined_metrics.remove("analysis_parameters");

        // Pass through the targeting metrics if all inputs are identical.
        if let Ok(targeted_metrics) = targeted_metrics.into_iter().all_equal_value() {
            combined_metrics.extend(targeted_metrics);
        }

        Ok(serde_json::to_string(&combined_metrics)?)
    }

    pub fn merge_barcode_infos(
        mut bc_infos: Vec<(Array2<u64>, Vec<GenomeName>)>,
    ) -> (Array2<u64>, Vec<GenomeName>) {
        assert!(!bc_infos.is_empty());
        let (last_pf, genomes) = bc_infos.pop().unwrap();

        let mut pfs = Vec::with_capacity(bc_infos.len() + 1);
        for (pf, genome) in &bc_infos {
            assert_eq!(pf.shape()[1], 3);
            assert_eq!(genome, &genomes);
            pfs.push(pf.view());
        }
        pfs.push(last_pf.view());

        let new_pf = ndarray::concatenate(Axis(0), &pfs).unwrap();

        (new_pf, genomes)
    }
}

/// Manage writers for individual molecule info columns.
struct ColumnWriters {
    gem_group: ChunkedWriter<GemGroupType>,
    barcode_idx: ChunkedWriter<BarcodeIdxType>,
    feature_idx: ChunkedWriter<FeatureIdxType>,
    library_idx: ChunkedWriter<LibraryIdxType>,
    probe_idx: Option<ChunkedWriter<ProbeIdxType>>,
    umi: ChunkedWriter<UmiTypeT>,
    count: ChunkedWriter<CountType>,
    umi_type: ChunkedWriter<UmiTypeType>,
}

impl ColumnWriters {
    /// Initialize buffer writers for new columns.
    pub fn new(group: &Group, with_probes: bool, action: ColumnAction) -> Result<Self> {
        let probe_idx = if with_probes {
            Some(ChunkedWriter::in_group(
                group,
                PROBE_IDX_COL_NAME,
                MOL_INFO_BUFFER_SZ,
                action,
            )?)
        } else {
            None
        };
        Ok(Self {
            gem_group: ChunkedWriter::in_group(
                group,
                GEM_GROUP_COL_NAME,
                MOL_INFO_BUFFER_SZ,
                action,
            )?,
            barcode_idx: ChunkedWriter::in_group(
                group,
                BARCODE_IDX_COL_NAME,
                MOL_INFO_BUFFER_SZ,
                action,
            )?,
            feature_idx: ChunkedWriter::in_group(
                group,
                FEATURE_IDX_COL_NAME,
                MOL_INFO_BUFFER_SZ,
                action,
            )?,
            library_idx: ChunkedWriter::in_group(
                group,
                LIBRARY_IDX_COL_NAME,
                MOL_INFO_BUFFER_SZ,
                action,
            )?,
            umi: ChunkedWriter::in_group(group, UMI_COL_NAME, MOL_INFO_BUFFER_SZ, action)?,
            count: ChunkedWriter::in_group(group, COUNT_COL_NAME, MOL_INFO_BUFFER_SZ, action)?,
            umi_type: ChunkedWriter::in_group(
                group,
                UMI_TYPE_COL_NAME,
                MOL_INFO_BUFFER_SZ,
                action,
            )?,
            probe_idx,
        })
    }

    /// Write info for a single molecule.
    pub fn write(&mut self, gem_group: u16, barcode_idx: u64, count: &UmiCount) -> Result<()> {
        self.gem_group.write(gem_group)?;
        self.barcode_idx.write(barcode_idx)?;
        self.feature_idx.write(count.feature_idx)?;
        self.library_idx.write(count.library_idx)?;
        self.umi.write(count.umi)?;
        self.count.write(count.read_count)?;
        self.umi_type.write(count.utype.to_u32())?;
        if let Some(probe_idx_writer) = &mut self.probe_idx {
            // If FB data has no probe_idx, we want to keep the column lengths consistent,
            // so use a placeholder value.
            probe_idx_writer.write(count.probe_idx.unwrap_or(PROBE_IDX_SENTINEL_VALUE))?;
        }
        Ok(())
    }

    /// Flush all columns to disk.
    pub fn flush(&mut self) -> Result<()> {
        self.gem_group.flush()?;
        self.barcode_idx.flush()?;
        self.feature_idx.flush()?;
        self.library_idx.flush()?;
        if let Some(probe_idx) = &mut self.probe_idx {
            probe_idx.flush()?;
        }
        self.umi.flush()?;
        self.count.flush()?;
        self.umi_type.flush()?;
        Ok(())
    }
}

/// Open a molecule info file in read/write mode for modification.
///
/// Validates that the file is of the correct version.
pub fn open_for_modification(path: &Path) -> Result<File> {
    let f = File::open_rw(path)?;
    let version = check_version(&f)?;
    if version < CURRENT_VERSION {
        bail!(
            "Molecule info file {} was produced by an older software version.",
            path.display()
        );
    } else if version > CURRENT_VERSION {
        bail!(
            "Molecule info file {} was produced by a newer software version.",
            path.display()
        );
    }
    Ok(f)
}

// Specific to Visium HD samples
mod visium_hd {
    use super::{
        BARCODE_DATASET_NAME, BARCODE_INFO_GROUP_NAME, BarcodeIdxType, ColumnWriters, FA256, FABc,
        FullUmiCount, GENOMES_DATASET_NAME, GemGroupType, ITERATOR_CHUNK_SIZE,
        MoleculeInfoIterator, MoleculeInfoWriter, PASS_FILTER_DATASET_NAME, PROBE_IDX_COL_NAME,
        open_for_modification, write_barcode_info,
    };
    use crate::ColumnAction;
    use crate::iter::H5Iterator;
    use anyhow::Result;
    use barcode::BarcodeContent;
    use cr_types::GenomeName;
    use hdf5::File;
    use itertools::Itertools;
    use ndarray::Array2;
    use shardio::{ShardReader, ShardWriter, SortKey};
    use std::borrow::Cow;
    use std::collections::HashMap;
    use std::path::Path;

    struct GemGroupBarcodeOrder;
    impl SortKey<FullUmiCount> for GemGroupBarcodeOrder {
        type Key = (GemGroupType, BarcodeIdxType);
        fn sort_key(t: &FullUmiCount) -> Cow<'_, Self::Key> {
            Cow::Owned((t.gem_group, t.barcode_idx))
        }
    }

    /// For Visim HD samples, bin the barcodes at the given bin scale and update the
    /// molecule info datasets.
    ///
    /// The molecule info file is modified in place. We also sort the columns in the molecule info file
    /// by barcode as the barcode order is changed after binning.
    pub fn bin_barcodes(molecule_h5: &Path, bin_scale: u32, tmp_shard_path: &Path) -> Result<()> {
        let (binned_barcodes_unique_sorted, binned_barcode_idx_per_molecule) =
            bin_spatial_barcodes(molecule_h5, bin_scale)?;

        // use shardio for out-of-memory sorting
        write_shard_file(
            tmp_shard_path,
            molecule_h5,
            &binned_barcode_idx_per_molecule,
        )?;

        let file = open_for_modification(molecule_h5)?;
        update_barcode_info_with_binned_barcodes(&file, &binned_barcode_idx_per_molecule)?;

        // Write binned barcodes
        file.unlink(BARCODE_DATASET_NAME)?;
        MoleculeInfoWriter::write_barcodes(&file, binned_barcodes_unique_sorted.into_iter())?;
        // Write the columns in sorted order
        let mut column_writers = ColumnWriters::new(
            &file,
            file.link_exists(PROBE_IDX_COL_NAME),
            ColumnAction::ReplaceExisting,
        )?;
        for umi_count in
            ShardReader::<FullUmiCount, GemGroupBarcodeOrder>::open(tmp_shard_path)?.iter()?
        {
            let umi_count = umi_count?;
            column_writers.write(
                umi_count.gem_group,
                umi_count.barcode_idx,
                &umi_count.umi_data,
            )?;
        }
        column_writers.flush()?;
        file.close()?;

        Ok(())
    }

    /// Write a shard file with the molecules in sorted order by the binned barcode.
    ///
    /// We are using the shard file to avoid loading all the data into memory. We will
    /// instead buffer the data in chunks and write them to the shard file.
    ///
    /// The shard file is then read back in sorted order and the molecule info datasets
    /// are updated in place.
    fn write_shard_file(
        tmp_shard_path: &Path,
        molecule_h5: &Path,
        binned_barcode_idx_per_molecule: &[u64],
    ) -> Result<()> {
        let mut shard_writer: ShardWriter<FullUmiCount, GemGroupBarcodeOrder> =
            ShardWriter::new(tmp_shard_path, 256, 8192, 1_048_576)?;
        let mut sender = shard_writer.get_sender();
        for umi_count in MoleculeInfoIterator::new(molecule_h5)? {
            let new_barcode_idx = binned_barcode_idx_per_molecule[umi_count.barcode_idx as usize];
            sender.send(FullUmiCount {
                barcode_idx: new_barcode_idx,
                ..umi_count
            })?;
        }
        sender.finished()?;
        shard_writer.finish()?;
        Ok(())
    }

    /// Bin the barcodes at the given bin scale.
    ///
    /// Returns the binned barcodes and the binned barcode index of each molecule
    /// in the molecule info file.
    fn bin_spatial_barcodes(
        molecule_h5: &Path,
        bin_scale: u32,
    ) -> Result<(Vec<BarcodeContent>, Vec<u64>)> {
        let file = File::open(molecule_h5)?;
        let binned_barcodes =
            H5Iterator::<FABc>::new(file.dataset(BARCODE_DATASET_NAME)?, ITERATOR_CHUNK_SIZE)
                .map(|x| x.and_then(|y| BarcodeContent::from_bytes(y.as_bytes())))
                .process_results(|iter| {
                    iter.map(|x| match x {
                        BarcodeContent::SpatialIndex(index) => {
                            BarcodeContent::SpatialIndex(index.binned(bin_scale))
                        }
                        _ => unreachable!("bin_barcodes called on non-spatial barcode"),
                    })
                    .collect::<Vec<_>>()
                })?;

        let binned_barcodes_unique_sorted: Vec<BarcodeContent> =
            binned_barcodes.iter().unique().copied().sorted().collect();

        let binned_barcodes_index: HashMap<_, _> = binned_barcodes_unique_sorted
            .iter()
            .enumerate()
            .map(|(a, &b)| (b, a as u64))
            .collect();

        let binned_barcode_idx_per_molecule: Vec<u64> = binned_barcodes
            .into_iter()
            .map(|x| binned_barcodes_index[&x])
            .collect();
        file.flush()?;
        file.close()?;
        Ok((
            binned_barcodes_unique_sorted,
            binned_barcode_idx_per_molecule,
        ))
    }

    /// Update the barcode info datasets with the binned barcodes.
    fn update_barcode_info_with_binned_barcodes(
        file: &File,
        binned_barcode_idx_per_molecule: &[u64],
    ) -> Result<()> {
        let old_barcode_info = file.group(BARCODE_INFO_GROUP_NAME)?;
        let genomes: Vec<GenomeName> = old_barcode_info
            .dataset(GENOMES_DATASET_NAME)?
            .read_1d::<FA256>()?
            .into_iter()
            .map(|x| GenomeName::from(x.as_str()))
            .collect();
        let pass_filter = Array2::from(
            old_barcode_info
                .dataset(PASS_FILTER_DATASET_NAME)?
                .read_2d::<u64>()?
                .outer_iter()
                .map(|x| {
                    let &[barcode_index, library_index, genome_index] = x.as_slice().unwrap()
                    else {
                        unreachable!();
                    };
                    [
                        binned_barcode_idx_per_molecule[barcode_index as usize],
                        library_index,
                        genome_index,
                    ]
                })
                .unique()
                .sorted()
                .collect_vec(),
        );
        file.unlink(BARCODE_INFO_GROUP_NAME)?;
        write_barcode_info(
            &file.create_group(BARCODE_INFO_GROUP_NAME)?,
            &pass_filter,
            &genomes,
        )?;
        Ok(())
    }
}

/// Trim barcodes
///
/// Changes the data in the dataset so that only barcodes that were
/// pass filtered (and optionally those with >0 counts) are stored
/// in the file.
pub fn trim_barcodes(path: &Path, pass_only: bool) -> Result<()> {
    let file = File::open_rw(path)?;
    // TODO: check barcodes here
    let old_barcodes = file.dataset(BARCODE_DATASET_NAME)?;

    let old_barcode_info = file.group(BARCODE_INFO_GROUP_NAME)?;
    let mut new_pass_filter = old_barcode_info
        .dataset(PASS_FILTER_DATASET_NAME)?
        .read_2d::<u64>()?;

    // Collect all barcodes in pass_filter
    let mut bc_idx_to_retain: BTreeSet<usize> = new_pass_filter
        .slice(s![.., 0])
        .iter()
        .map(|&x| x as usize)
        .collect();

    let barcode_idx_ds = file.dataset(BARCODE_IDX_COL_NAME)?;

    if !pass_only {
        // Also keep any barcode with count > 0
        for bc in H5Iterator::<BarcodeIdxType>::new(
            file.dataset(BARCODE_IDX_COL_NAME)?,
            ITERATOR_CHUNK_SIZE,
        ) {
            bc_idx_to_retain.insert(bc? as usize);
        }
    }

    // Implementation note: this is an optimization to avoid checking a hashmap for
    // new indices for the old barcodes, but 0 is a valid position and so can result
    // in wrong positions.
    let mut new_positions = vec![None; old_barcodes.size()];
    for (pos, &v) in bc_idx_to_retain.iter().enumerate() {
        new_positions[v] = Some(pos);
    }

    // update chunks
    let size = barcode_idx_ds.size();
    let mut index = 0;

    while index < size {
        let end = min(size, index + ITERATOR_CHUNK_SIZE);
        let mut future: Array1<usize> = barcode_idx_ds.read_slice_1d(index..end)?;
        for x in &mut future {
            *x = new_positions[*x].expect("Error accessing invalid barcode index");
        }
        barcode_idx_ds.write_slice(future.view(), index..end)?;
        index += ITERATOR_CHUNK_SIZE;
    }

    // Update barcode_info
    for x in new_pass_filter.slice_mut(s![.., 0]) {
        *x = new_positions[*x as usize].expect("Error accessing invalid barcode index") as u64;
    }

    // TODO: Jira: CELLRANGER-9147: Consider using reference_genomes list instead
    let genomes: Vec<GenomeName> = old_barcode_info
        .dataset(GENOMES_DATASET_NAME)?
        .read_1d::<FA256>()?
        .into_iter()
        .map(|x| GenomeName::from(x.as_str()))
        .collect();

    // Write new barcode info.
    file.unlink(BARCODE_INFO_GROUP_NAME)?;
    write_barcode_info(
        &file.create_group(BARCODE_INFO_GROUP_NAME)?,
        &new_pass_filter,
        &genomes,
    )?;

    // Write new barcodes by streaming the old dataset directly into a new dataset.
    let temp_new_bc_ds_name = "temp_new_bcs";
    process_results(
        H5Iterator::<FABc>::new(old_barcodes, ITERATOR_CHUNK_SIZE)
            .enumerate()
            .filter_map(|(i, bc_read_result)| {
                if !bc_idx_to_retain.contains(&i) {
                    return None;
                }
                Some(bc_read_result)
            }),
        |new_barcodes_iter| {
            ChunkedWriter::write_all(
                &file,
                temp_new_bc_ds_name,
                MOL_INFO_BUFFER_SZ,
                ColumnAction::CreateNew,
                new_barcodes_iter,
            )
        },
    )??;
    file.unlink(BARCODE_DATASET_NAME)?;
    file.relink(temp_new_bc_ds_name, BARCODE_DATASET_NAME)?;
    Ok(())
}

/// write the /barcode_info/pass_filter and /barcode_info/genomes datsets
/// to h5 file.
fn write_barcode_info(
    barcode_info_group: &Group,
    pass_filter: &Array2<u64>,
    genomes: &[GenomeName],
) -> Result<()> {
    let ds = barcode_info_group
        .new_dataset::<u64>()
        .shuffle()
        .deflate(1)
        .chunk((1 << 16, 16))
        .shape(Extents::resizable(pass_filter.shape().into()))
        .create(PASS_FILTER_DATASET_NAME)?;

    ds.as_writer().write(pass_filter.view())?;

    let genome_array =
        Array1::from_shape_fn(genomes.len(), |i| FA256::from_ascii(&genomes[i]).unwrap());

    let ds = barcode_info_group
        .new_dataset::<FA256>()
        .shuffle()
        .deflate(1)
        .chunk(1 << 16)
        .shape(Extents::resizable(genome_array.shape().into()))
        .create(GENOMES_DATASET_NAME)?;

    ds.as_writer().write(genome_array.view())?;

    Ok(())
}

#[cfg(test)]
mod molecule_info_tests {

    use super::*;
    use cr_types::reference::feature_reference::FeatureType;
    use cr_types::types::FeatureBarcodeType;
    use std::fs::{Permissions, copy, set_permissions};
    use std::os::unix::fs::PermissionsExt;
    use tempfile::NamedTempFile;

    #[test]
    fn concatenate_metrics_1() {
        let m1 = serde_json::json!({
        "cellranger_version":"2021.0428.1",
        "chemistry_barcode":[{"kind":"gel_bead",
        "length":16,
        "offset":0,
        "read_type":"R1",
        "whitelist":"3M-february-2018"}],
        "chemistry_description":"Single Cell 3' v3",
        "chemistry_endedness":"three_prime",
        "chemistry_name":"SC3Pv3",
        "chemistry_rna":{"length":Value::Null,
        "min_length":15,
        "offset":0,
        "read_type":"R2"},
        "chemistry_rna2":Value::Null,
        "chemistry_strandedness":"+",
        "chemistry_umi":{"length":12,
        "min_length":10,
        "offset":16,
        "read_type":"R1"},
        "gem_groups":{"1":{"force_cells":Value::Null,
        "include_introns":false,
        "recovered_cells":2000}},
        "is_aggregated":true,
        "libraries":{"0":{"feature_read_pairs":980404509,
        "raw_read_pairs":1004464914,
        "usable_read_pairs":371030925}},
        "molecule_info_type":"count",
        "reference_fasta_hash":"b6f131840f9f337e7b858c3d1e89d7ce0321b243",
        "reference_gtf_hash":"d9aa710e1eab4e2bdb1a87c25d4cc8a9397db121",
        "reference_gtf_hash.gz":"",
        "reference_mkref_version":""
           });
        let m1 = m1.as_object().expect("Error parsing input data");

        let mut m2 = m1.clone();
        m2["gem_groups"] = serde_json::json!({
          "2": {
            "force_cells":Value::Null,
            "include_introns":false,
            "recovered_cells":2000
        }});

        m2["libraries"] = serde_json::json!({
          "1":{
            "feature_read_pairs":980404509,
            "raw_read_pairs":1004464914,
            "usable_read_pairs":371030925
        }});

        let metrics = Value::from_str(
            &MoleculeInfoWriter::concatenate_metrics(vec![m1.clone(), m2]).unwrap(),
        )
        .unwrap();
        let metrics = metrics.as_object().unwrap();

        assert_eq!(
            metrics["gem_groups"],
            serde_json::json!(
            {
              "1": {
                "force_cells":Value::Null,
                "include_introns":false,
                "recovered_cells":2000},
              "2": {
                "force_cells":Value::Null,
                "include_introns":false,
                "recovered_cells":2000},
            })
        );

        assert_eq!(
            metrics["libraries"],
            serde_json::json!(
            {
                "0":{
                    "feature_read_pairs":980404509,
                    "raw_read_pairs":1004464914,
                    "usable_read_pairs":371030925},
                "1":{
                    "feature_read_pairs":980404509,
                    "raw_read_pairs":1004464914,
                    "usable_read_pairs":371030925}
            })
        );
    }

    #[test]
    fn test_mol_info_reader() {
        let path = Path::new("test/h5/pbmc_1k_v2_molecule_info.h5");
        let filtered_barcodes =
            MoleculeInfoReader::read_filtered_barcodes(path, LibraryType::Gex).unwrap();
        let filtered_barcode_ids = MoleculeInfoReader::read_filtered_barcode_ids(path).unwrap();
        dbg!(filtered_barcode_ids.len());
        assert_eq!(filtered_barcodes.len(), 996);
        assert_eq!(filtered_barcodes.len(), filtered_barcode_ids.len());

        let fref = MoleculeInfoReader::read_feature_ref(path).unwrap();
        assert_eq!(fref.num_features(), 33538);

        for (i, fdef) in fref.feature_defs.iter().enumerate() {
            assert_eq!(i, fdef.index);
        }

        let num_bcs = MoleculeInfoReader::read_barcodes_size(path).unwrap();
        assert_eq!(737280, num_bcs);
        let num_nz_bcs = MoleculeInfoReader::count_non_zero_bcs(path).unwrap();
        assert_eq!(142324, num_nz_bcs);
    }

    #[test]
    fn test_mol_info_iter() -> Result<()> {
        let mol_info_path = Path::new("test/h5/pbmc_1k_v2_molecule_info.h5");

        assert_eq!(MoleculeInfoIterator::new(mol_info_path)?.count(), 4303635);

        assert_eq!(
            MoleculeInfoIterator::new(mol_info_path)?
                .cell_barcodes_only(true)?
                .count(),
            3524098
        );

        assert_eq!(
            MoleculeInfoIterator::new(mol_info_path)?
                .filter_features(|fdef| fdef.index == 500)?
                .count(),
            79
        );

        assert_eq!(
            MoleculeInfoIterator::new(mol_info_path)?
                .cell_barcodes_only(true)?
                .filter_features(|fdef| fdef.index == 500)?
                .count(),
            69
        );
        Ok(())
    }

    #[test]
    fn test_mol_info_iter_antigen() -> Result<()> {
        let mol_info_path = Path::new("test/h5/antigen_tiny.h5");
        assert_eq!(
            MoleculeInfoIterator::new(mol_info_path)?
                .cell_barcodes_only(true)?
                .filter_features(
                    |fdef| fdef.feature_type == FeatureType::Barcode(FeatureBarcodeType::Antigen)
                )?
                .count(),
            55193
        );
        assert_eq!(
            MoleculeInfoReader::read_filtered_barcode_ids(mol_info_path)
                .unwrap()
                .len(),
            115
        );
        Ok(())
    }

    #[test]
    fn test_trim_barcodes() -> Result<()> {
        let template = Path::new("test/h5/pbmc_1k_v2_molecule_info.h5");
        let tempfile = NamedTempFile::new()?;
        let path = tempfile.path();
        copy(template, path)?;
        set_permissions(path, Permissions::from_mode(0o644))?;

        let num_bcs = MoleculeInfoReader::read_barcodes_size(path)?;
        assert_eq!(737280, num_bcs);
        let num_nz_bcs = MoleculeInfoReader::count_non_zero_bcs(path)?;
        assert_eq!(142324, num_nz_bcs);

        trim_barcodes(path, false)?;

        assert_eq!(num_nz_bcs, MoleculeInfoReader::read_barcodes_size(path)?);
        Ok(())
    }
}
