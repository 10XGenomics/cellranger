//! Wrapping the Rust handling of the multi config.
use cr_types::{CellMultiplexingType, CrMultiGraph, Fingerprint, HasMultiplexing, LibraryFeatures};
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
            .filter_map(|lib| match lib.library_features {
                LibraryFeatures::GeneExpression(_) | LibraryFeatures::Vdj(_) => None,
                LibraryFeatures::FeatureBarcodes(features) => {
                    Some(features.iter().map(|f| f.to_string()))
                }
            })
            .flatten()
            .collect()
    }

    /// Return true if data is multiplexed.
    pub fn is_multiplexed(&self) -> bool {
        self.0.has_multiplexing()
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
}
