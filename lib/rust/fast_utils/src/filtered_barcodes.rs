use barcode::Barcode;
use cr_types::filtered_barcodes::{FilteredBarcodesCsv, FilteredBarcodesCsvRow};
use cr_types::GenomeName;
use itertools::Itertools;
use martian_filetypes::bin_file::BincodeFile;
use martian_filetypes::{FileTypeRead, FileTypeWrite};
use numpy::PyArray1;
use pyanyhow::Result;
use pyo3::prelude::*;
use pyo3::types::PyBytes;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeSet, HashMap};
use std::path::PathBuf;

type GenomeNameStr = str;
type GenomeNameString = String;

#[pyclass]
pub struct FilteredBarcodes {
    /// Sorted list of unique barcodes.
    sorted_barcodes: Vec<Barcode>,

    /// List of unique genomes.
    genomes: Vec<GenomeName>,

    /// The (barcode, genome) pairs are stored as a CSR sparse matrix.
    /// In case of barnyard samples, the same barcode can appear multiple times
    /// in the filtered barcodes csv with different genomes.
    barcodes_indptr: Vec<usize>,
    genome_indices: Vec<usize>,
}

enum GenomeChoice {
    All,
    One(usize),
    Unobserved,
}

impl GenomeChoice {
    fn contains(&self, index: usize) -> bool {
        match self {
            GenomeChoice::All => true,
            GenomeChoice::One(genome_index) => *genome_index == index,
            GenomeChoice::Unobserved => false,
        }
    }
}

#[pymethods]
impl FilteredBarcodes {
    #[new]
    pub fn new(csv_file: PathBuf) -> Result<Self> {
        let mut contents = FilteredBarcodesCsv::from(csv_file).read()?;
        contents.sort_unstable_by_key(|f| f.barcode);

        let mut genomes = vec![];
        let mut sorted_barcodes = Vec::with_capacity(contents.len());
        let mut genome_indices = Vec::with_capacity(contents.len());
        let mut barcodes_indptr = Vec::with_capacity(contents.len());

        for FilteredBarcodesCsvRow { genome, barcode } in contents {
            // Cheap linear lookup because the number of genomes is small.
            let genome_index = match genomes.iter().position(|g| g == &genome) {
                Some(index) => index,
                None => {
                    genomes.push(genome);
                    genomes.len() - 1
                }
            };
            if sorted_barcodes.last() != Some(&barcode) {
                barcodes_indptr.push(genome_indices.len());
                sorted_barcodes.push(barcode);
            }
            genome_indices.push(genome_index);
        }
        barcodes_indptr.push(genome_indices.len());
        Ok(FilteredBarcodes {
            sorted_barcodes,
            genomes,
            genome_indices,
            barcodes_indptr,
        })
    }

    pub fn num_cells(&self) -> usize {
        self.sorted_barcodes.len()
    }

    pub fn sorted_barcode_indices<'py>(
        &self,
        py: Python<'py>,
        genome: Option<&GenomeNameStr>,
    ) -> &'py PyArray1<usize> {
        let genome_choice = self.genome_choice(genome);
        PyArray1::from_iter(
            py,
            self.iter_barcode_genome_indices()
                .filter_map(|(bc_idx, genome_idx)| {
                    genome_choice.contains(genome_idx).then_some(bc_idx)
                })
                .dedup(),
        )
    }

    pub fn index_of_barcode(&self, barcode: &PyAny) -> usize {
        self._index_of_barcode(barcode)
            .unwrap_or_else(|| panic!("Could not find index for barcode: {barcode}"))
    }

    pub fn per_genome_barcodes(
        &self,
        py: Python<'_>,
    ) -> HashMap<GenomeNameString, Vec<Py<PyBytes>>> {
        let genome_to_barcodes = self
            .iter_barcode_genome_indices()
            .map(|(bc_idx, genome_idx)| {
                (
                    genome_idx,
                    PyBytes::new(py, self.sorted_barcodes[bc_idx].to_string().as_bytes()).into(),
                )
            })
            .fold(
                vec![Vec::new(); self.genomes.len()],
                |mut acc, (genome_idx, barcode)| {
                    acc[genome_idx].push(barcode);
                    acc
                },
            );

        genome_to_barcodes
            .into_iter()
            .enumerate()
            .map(|(genome_idx, barcodes)| (self.genomes[genome_idx].to_string(), barcodes))
            .collect()
    }

    pub fn sorted_barcodes(&self, py: Python<'_>) -> Vec<Py<PyBytes>> {
        self.sorted_barcodes
            .iter()
            .map(|bc| PyBytes::new(py, bc.to_string().as_bytes()).into())
            .collect()
    }

    pub fn sorted_barcodes_string(&self) -> Vec<String> {
        self.sorted_barcodes
            .iter()
            .map(ToString::to_string)
            .collect()
    }

    /// Returns true if the barcode is in the list of filtered barcodes for the given genome.
    ///
    /// NOTE: In the python code genome = "" is used to refer to all genomes.
    pub fn contains(&self, barcode: &PyAny, genome: Option<&GenomeNameStr>) -> bool {
        self._index_of_barcode(barcode).is_some_and(|bc_idx| {
            let genome_choice = self.genome_choice(genome);
            self.iter_genome_indices(bc_idx)
                .any(|genome_idx| genome_choice.contains(genome_idx))
        })
    }

    pub fn cells_per_gem_group(&self) -> HashMap<u16, usize> {
        self.sorted_barcodes
            .iter()
            .map(|bc| bc.gem_group())
            .counts()
    }
}

impl FilteredBarcodes {
    // Returns an iterator over all the (barcode index, genome index) pairs.
    fn iter_barcode_genome_indices(&self) -> impl Iterator<Item = (usize, usize)> + '_ {
        (0..self.sorted_barcodes.len()).flat_map(move |bc_idx| {
            self.iter_genome_indices(bc_idx)
                .map(move |genome_idx| (bc_idx, genome_idx))
        })
    }

    // Returns an iterator over the genome indices for the given barcode index.
    fn iter_genome_indices(&self, bc_idx: usize) -> impl Iterator<Item = usize> + '_ {
        (self.barcodes_indptr[bc_idx]..self.barcodes_indptr[bc_idx + 1])
            .map(move |idx| self.genome_indices[idx])
    }

    fn genome_choice(&self, genome: Option<&GenomeNameStr>) -> GenomeChoice {
        match genome {
            Some(genome) if !genome.is_empty() => {
                // In OCM barnyard, multiplexed samples can have cells called exclusively from one species
                self.index_of_genome(genome)
                    .map_or(GenomeChoice::Unobserved, GenomeChoice::One)
            }
            Some(_) | None => GenomeChoice::All,
        }
    }

    fn _index_of_barcode(&self, barcode: &PyAny) -> Option<usize> {
        let barcode: Barcode = Self::parse_barcode(barcode).unwrap();
        self.sorted_barcodes.binary_search(&barcode).ok()
    }

    fn parse_barcode(barcode_str_or_bytes: &PyAny) -> Result<Barcode> {
        // Check if the input is a bytes object
        if let Ok(bytes) = barcode_str_or_bytes.extract::<&[u8]>() {
            return Ok(Barcode::from_bytes(bytes)?);
        }

        // Check if the input is a string
        if let Ok(string) = barcode_str_or_bytes.extract::<String>() {
            return Ok(string.parse().unwrap());
        }

        // Otherwise return an error
        Err(anyhow::anyhow!("barcode must be bytes or str, got {barcode_str_or_bytes:?}").into())
    }

    fn index_of_genome(&self, genome: &GenomeNameStr) -> Option<usize> {
        self.genomes.iter().position(|g| g.as_str() == genome)
    }
}

#[derive(Serialize, Deserialize, PartialEq, Eq, Hash)]
struct Key {
    gem_group: u16,
    genome: String,
}

type FilteredBcsGroupsFile = BincodeFile<Vec<(Key, Vec<Barcode>)>>;

#[pyfunction]
pub fn save_filtered_bcs_groups(
    py: Python<'_>,
    filtered_bcs_groups: HashMap<(u16, String), Vec<Py<PyBytes>>>,
    file_name: PathBuf,
) -> Result<()> {
    FilteredBcsGroupsFile::from(file_name).write(
        &filtered_bcs_groups
            .into_iter()
            .map(|(k, barcodes)| {
                (
                    Key {
                        gem_group: k.0,
                        genome: k.1,
                    },
                    barcodes
                        .into_iter()
                        .map(|bc| Barcode::from_bytes(bc.as_bytes(py)).unwrap())
                        .collect(),
                )
            })
            .collect(),
    )?;
    Ok(())
}

#[pyfunction]
#[allow(clippy::type_complexity)]
pub fn load_filtered_bcs_groups(
    py: Python<'_>,
    files: Vec<PathBuf>,
) -> Result<(Vec<Py<PyBytes>>, HashMap<String, Vec<Py<PyBytes>>>)> {
    let mut all_barcodes = BTreeSet::new();
    let mut per_genome_barcodes: HashMap<String, BTreeSet<Barcode>> = HashMap::new();
    for file in files {
        let contents = FilteredBcsGroupsFile::from(file).read()?;
        for (Key { genome, .. }, barcodes) in contents {
            for barcode in barcodes {
                all_barcodes.insert(barcode);
                per_genome_barcodes
                    .entry(genome.clone())
                    .or_default()
                    .insert(barcode);
            }
        }
    }

    fn convert(py: Python<'_>, barcodes: BTreeSet<Barcode>) -> Vec<Py<PyBytes>> {
        barcodes
            .into_iter()
            .map(|bc| PyBytes::new(py, bc.to_string().as_bytes()).into())
            .collect()
    }
    Ok((
        convert(py, all_barcodes),
        per_genome_barcodes
            .into_iter()
            .map(|(genome, barcodes)| (genome, convert(py, barcodes)))
            .collect(),
    ))
}
