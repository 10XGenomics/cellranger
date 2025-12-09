#![deny(missing_docs)]
use barcode::BcSeq;
use cr_h5::molecule_info::{
    BarcodeIdxType, FullUmiCount, GemGroupType, LibraryIdxType, MoleculeInfoIterator,
    MoleculeInfoReader,
};
use cr_types::reference::feature_reference::{FeatureDef, FeatureType, HASHTAG};
use cr_types::rna_read::LibraryInfo;
use cr_types::types::FeatureBarcodeType;
use cr_types::{BarcodeMultiplexingType, CellLevel, LibraryType};
use itertools::Itertools;
use pyo3::prelude::*;
use std::collections::{HashMap, HashSet};
use std::convert::TryInto;
use std::fs::File;
use std::path::Path;

/// The file non_singlet_barcodes.json has three top-level keys: Blank, Multiplet, Unassigned.
const MULTIPLET_JSON_KEY: &str = "Multiplet";

// A struct to accumulate data needed by the COMPUTE_EXTRA_MULTIPLEXING_METRICS stage
struct TagReadCounts {
    // Maps the feature index to the local index within the vectors defined below
    feature_index_to_local_index: HashMap<usize, usize>,
    // name of all cmos
    tag_names: Vec<String>,
    // Vectors to count reads for various types of barcodes
    cell_associated_counts: Vec<usize>,
    multiplet_counts: Vec<usize>,
    singlet_counts: Vec<usize>,
    total_counts: Vec<usize>,
    // HashSets to identiyf different types of barcodes
    singlets: HashSet<(BarcodeIdxType, GemGroupType)>,
    multiplets: HashSet<(BarcodeIdxType, GemGroupType)>,
    cell_associated: HashSet<(BarcodeIdxType, GemGroupType)>,
}

impl TagReadCounts {
    fn new(
        features: Vec<FeatureDef>,
        singlets: HashSet<(BarcodeIdxType, GemGroupType)>,
        multiplets: HashSet<(BarcodeIdxType, GemGroupType)>,
        cell_associated: HashSet<(BarcodeIdxType, GemGroupType)>,
    ) -> Self {
        assert!(
            !features.is_empty(),
            "No multiplexing tags were found in the FeatureRef."
        );
        // Assuming all tags are in
        let feature_index_to_local_index = features
            .iter()
            .enumerate()
            .map(|(i, feature_def)| (feature_def.index, i))
            .collect();
        let size = features.len();
        TagReadCounts {
            feature_index_to_local_index,
            singlet_counts: vec![0; size],
            multiplet_counts: vec![0; size],
            cell_associated_counts: vec![0; size],
            total_counts: vec![0; size],
            singlets,
            multiplets,
            cell_associated,
            tag_names: features.into_iter().map(|f| f.id).collect(),
        }
    }

    fn consume_umicount(&mut self, cnt: FullUmiCount) {
        let idx = cnt.umi_data.feature_idx as usize;
        if let Some(&index) = self.feature_index_to_local_index.get(&idx) {
            let n_reads = cnt.umi_data.read_count as usize;
            self.total_counts[index] += n_reads;
            let bc_gg = (cnt.barcode_idx, cnt.gem_group);
            if self.singlets.contains(&bc_gg) {
                self.singlet_counts[index] += n_reads;
                self.cell_associated_counts[index] += n_reads;
            } else if self.multiplets.contains(&bc_gg) {
                self.multiplet_counts[index] += n_reads;
                self.cell_associated_counts[index] += n_reads;
            } else if self.cell_associated.contains(&bc_gg) {
                self.cell_associated_counts[index] += n_reads;
            }
        }
    }
}

fn load_tag_json(
    cur_path: &Path,
    bc_to_ind: &HashMap<BcSeq, usize>,
    select_tag: Option<&str>,
) -> HashSet<(BarcodeIdxType, GemGroupType)> {
    // Given a JSON file that is a serialized CellsPerFeature python object
    // decode it to a set of barcode_idx and gem group values
    let json_data = File::open(cur_path).unwrap();
    let data: HashMap<String, Vec<String>> = serde_json::from_reader(json_data).unwrap();
    let mut bc_gg: HashSet<(BarcodeIdxType, GemGroupType)> = HashSet::new();
    for (tag, barcodes) in data {
        if matches!(select_tag, Some(t) if t != tag) {
            continue;
        }
        for bc in barcodes {
            let sp: Vec<&str> = bc.split('-').collect();
            let gg: GemGroupType = sp[1].parse().unwrap();
            let bc = BcSeq::from_bytes(sp[0].as_bytes());
            let bc_ind = bc_to_ind[&bc] as u64;
            bc_gg.insert((bc_ind, gg));
        }
    }
    bc_gg
}

fn make_barcode_to_index_map(mol_h5_fname: &Path) -> anyhow::Result<HashMap<BcSeq, usize>> {
    // Gets a map translating a barcode sequence to a index value
    MoleculeInfoReader::read_barcodes(mol_h5_fname)?
        .enumerate()
        .map(|(index, barcode)| Ok((*barcode?.sequence(), index)))
        .try_collect()
}

// Method counts the total number of reads for each feature
#[allow(clippy::type_complexity)]
#[pyfunction]
pub(super) fn tag_read_counts(
    _py: Python<'_>,
    mol_info_path: String,
    singlets_path: String,
    non_singlets_path: String,
    multiplexing_method: String,
) -> PyResult<(Vec<String>, Vec<usize>, Vec<usize>, Vec<usize>, Vec<usize>)> {
    // Calculates metrics for the COMPUTE_EXTRA_MULTIPLEXING_METRICS stage
    // Given a molecule info file, it will go through and for each tag CMO or Hashtag count how many reads are
    // 1 - In "singlet" barcodes
    // 2 - In "multiplet" barcodes
    // 3-  In cell-associated barcodes
    // 4 - Total Reads

    let cur_path = Path::new(&mol_info_path);

    // Load up singlet and multiplet barcodes
    let bc_to_ind = make_barcode_to_index_map(cur_path).unwrap();
    let singlets = load_tag_json(singlets_path.as_ref(), &bc_to_ind, None);
    let multiplets = load_tag_json(
        non_singlets_path.as_ref(),
        &bc_to_ind,
        Some(MULTIPLET_JSON_KEY),
    );
    let (library_type, feature_type, feature_tag) =
        match serde_json::from_str(&multiplexing_method).unwrap() {
            BarcodeMultiplexingType::CellLevel(CellLevel::CMO) => (
                LibraryType::Cellplex,
                FeatureType::Barcode(FeatureBarcodeType::Multiplexing),
                None,
            ),
            BarcodeMultiplexingType::CellLevel(CellLevel::Hashtag) => (
                LibraryType::Antibody,
                FeatureType::Barcode(FeatureBarcodeType::Antibody),
                Some(HASHTAG),
            ),
            _ => panic!("Unsupported multiplexing method!"),
        };

    // Load up cell associated barcodes for multiplexing libraries
    let library_info: Vec<LibraryInfo> = MoleculeInfoReader::read_library_info(cur_path)
        .unwrap()
        .into_iter()
        .map(TryInto::try_into)
        .try_collect()
        .unwrap();

    let mut lib_to_gg: HashMap<LibraryIdxType, GemGroupType> = HashMap::new();
    library_info
        .into_iter()
        .filter(|x| x.library_type == library_type)
        .for_each(|z| {
            lib_to_gg.insert(z.library_id, z.gem_group);
        });
    let barcodes = MoleculeInfoReader::read_barcode_info_pass_filter(cur_path).unwrap();
    let mut cell_associated_bcs: HashSet<(BarcodeIdxType, GemGroupType)> =
        HashSet::with_capacity(barcodes.dim().0 / 2); // divide by 2 as we typically have GEX + CMO at a minimum, and we only add in CMO
    barcodes.outer_iter().for_each(|z| {
        let gg = lib_to_gg.get(&(z[1] as u16));
        if let Some(x) = gg {
            cell_associated_bcs.insert((z[0], *x)); // (BC, GG)
        }
    });

    // Now to tally up different types

    let iterator = MoleculeInfoIterator::new(cur_path).unwrap();
    // What are the multiplexing features?
    let multiplexing_features = iterator
        .feature_ref
        .feature_defs
        .iter()
        .filter(|&x| x.feature_type == feature_type)
        .filter(|&x| match feature_tag {
            Some(tag) => x.tags.contains_key(tag),
            None => true,
        })
        .cloned()
        .collect();

    // Get ready to tally
    let mut cnt_consumer = TagReadCounts::new(
        multiplexing_features,
        singlets,
        multiplets,
        cell_associated_bcs,
    );

    // Counts up reads per feature
    //let float_cnts: Vec<f64> = Vec::FromIterator(cnts.map(|x| *x as f64));
    iterator.for_each(|x| {
        cnt_consumer.consume_umicount(x);
    });
    Ok((
        cnt_consumer.tag_names,
        cnt_consumer.singlet_counts,
        cnt_consumer.multiplet_counts,
        cnt_consumer.cell_associated_counts,
        cnt_consumer.total_counts,
    ))
}
