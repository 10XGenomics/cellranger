//! Wrapping the Rust handling of the multi config.
use cr_types::{CellMultiplexingType, CrMultiGraph, Fingerprint};
use itertools::Itertools;
use martian::MartianFileType;
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::FileTypeRead;
use pyo3::exceptions::PyException;
use pyo3::prelude::*;
use pyo3::types::PyType;
use std::collections::{HashMap, HashSet};
use std::path::Path;

/// Wrapper type for interacting with the multi graph from Python.
/// Add methods below as needed for extracting data.
/// TODO: it might make sense to lift the pyo3 annotations directly up to the
/// CrMultiGraph itself, rather than wrapping for Python consumption.
#[pyclass]
pub struct MultiGraph(CrMultiGraph);

#[pymethods]
impl MultiGraph {
    /// Load the multi graph from the provided path.
    #[classmethod]
    pub fn from_path(_cls: &PyType, path: &str) -> PyResult<Self> {
        Ok(Self(JsonFile::from_path(Path::new(path)).read().map_err(
            |err| PyErr::new::<PyException, _>(format!("{err:#}")),
        )?))
    }

    /// Load the multi graph from the provided JSON-encoded string.
    #[classmethod]
    pub fn from_str(_cls: &PyType, contents: &str) -> PyResult<Self> {
        Ok(Self(serde_json::from_str(contents).map_err(|err| {
            PyErr::new::<PyException, _>(format!("{err:#}"))
        })?))
    }

    pub fn sample_ids(&self) -> Vec<String> {
        self.0
            .samples
            .iter()
            .map(|sample| sample.sample_id.clone())
            .collect()
    }

    pub fn sample_tag_ids(&self) -> HashMap<String, Vec<String>> {
        self.0
            .samples
            .iter()
            .map(|sample| {
                (
                    sample.sample_id.clone(),
                    sample
                        .fingerprints
                        .iter()
                        .filter_map(|f| match f {
                            Fingerprint::Tagged { tag_name, .. } => Some(tag_name.clone()),
                            Fingerprint::Untagged { .. } => None,
                        })
                        .collect(),
                )
            })
            .collect()
    }

    /// Extract the names of all feature barcoding types present.
    pub fn feature_types(&self) -> HashSet<String> {
        self.0
            .libraries
            .iter()
            .filter_map(|lib| lib.library_type.feature_barcode_type())
            .unique()
            .map(|fbt| fbt.to_string())
            .collect()
    }

    /// Return true if data is multiplexed using any method.
    pub fn is_multiplexed(&self) -> bool {
        self.0.is_multiplexed()
    }

    /// Return true if this data is using CMO multiplexing.
    pub fn is_cmo_multiplexed(&self) -> bool {
        self.0.cell_multiplexing_type() == Some(CellMultiplexingType::CMO)
    }

    /// Return true if this data is using RTL multiplexing.
    pub fn is_rtl_multiplexed(&self) -> bool {
        self.0.cell_multiplexing_type() == Some(CellMultiplexingType::RTL)
    }

    /// Return true if this data is using OH multiplexing.
    pub fn is_oh_multiplexed(&self) -> bool {
        self.0.cell_multiplexing_type() == Some(CellMultiplexingType::OH)
    }

    /// Return tuples of sample ID and finerprints.
    /// THIS IS ONLY USED FOR CREATION OF GRAPHVIZ.
    pub fn sample_info_for_graphviz(&self) -> Vec<(String, Vec<PyFingerprint>)> {
        self.0
            .samples
            .iter()
            .map(|sample| {
                (
                    sample.sample_id.clone(),
                    sample
                        .fingerprints
                        .iter()
                        .map(PyFingerprint::from_fingerprint)
                        .collect::<Vec<_>>(),
                )
            })
            .collect()
    }

    /// Return tuples of physical library ID, gem well, associated FASTQ IDs.
    /// THIS IS ONLY USED FOR CREATION OF GRAPHVIZ.
    pub fn library_info_for_graphviz(&self) -> Vec<(String, u16, Vec<String>)> {
        self.0
            .libraries
            .iter()
            .map(|lib| {
                (
                    lib.physical_library_id.clone(),
                    lib.gem_well.0,
                    lib.fastqs
                        .iter()
                        .map(|fastq| fastq.id().to_string())
                        .collect::<Vec<_>>(),
                )
            })
            .collect()
    }
}

#[pyclass]
pub struct PyFingerprint {
    #[pyo3(get)]
    gem_well: u16,
    #[pyo3(get)]
    tag_names: Vec<String>,
}

impl PyFingerprint {
    fn from_fingerprint(fingerprint: &Fingerprint) -> Self {
        match fingerprint {
            Fingerprint::Tagged {
                gem_well,
                tag_name,
                translated_tag_names,
                ..
            } => Self {
                gem_well: gem_well.0,
                tag_names: std::iter::once(tag_name)
                    .chain(translated_tag_names)
                    .cloned()
                    .collect(),
            },
            Fingerprint::Untagged { gem_well } => Self {
                gem_well: gem_well.0,
                tag_names: Default::default(),
            },
        }
    }
}
