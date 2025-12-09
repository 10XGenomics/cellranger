//! Wrapping the Rust handling of the multi config.
#![deny(missing_docs)]
use cr_types::{BarcodeMultiplexingType, CellLevel, CrMultiGraph, Fingerprint, ReadLevel};
use martian::MartianFileType;
use martian_filetypes::FileTypeRead;
use martian_filetypes::json_file::JsonFile;
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
    pub fn from_path(_cls: &Bound<'_, PyType>, path: &str) -> PyResult<Self> {
        Ok(Self(JsonFile::from_path(Path::new(path)).read().map_err(
            |err| PyErr::new::<PyException, _>(format!("{err:#}")),
        )?))
    }

    /// Load the multi graph from the provided JSON-encoded string.
    #[classmethod]
    pub fn from_str(_cls: &Bound<'_, PyType>, contents: &str) -> PyResult<Self> {
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

    /// Get the names of all library types.
    pub fn library_types(&self) -> HashSet<&'static str> {
        self.0
            .libraries
            .iter()
            .map(|x| x.library_type.as_str())
            .collect()
    }

    /// Return true if data is multiplexed using any method.
    pub fn is_multiplexed(&self) -> bool {
        self.0.is_multiplexed()
    }

    /// Return true if this data is using CMO multiplexing.
    pub fn is_cmo_multiplexed(&self) -> bool {
        self.0.barcode_multiplexing_type()
            == Some(BarcodeMultiplexingType::CellLevel(CellLevel::CMO))
    }

    /// Return true if this data is using HASHTAG multiplexing.
    pub fn is_hashtag_multiplexed(&self) -> bool {
        self.0.barcode_multiplexing_type()
            == Some(BarcodeMultiplexingType::CellLevel(CellLevel::Hashtag))
    }

    /// Return true if this data is using RTL multiplexing.
    pub fn is_rtl_multiplexed(&self) -> bool {
        self.0.barcode_multiplexing_type()
            == Some(BarcodeMultiplexingType::ReadLevel(ReadLevel::RTL))
    }

    /// Return true if this data is using OH multiplexing.
    pub fn is_oh_multiplexed(&self) -> bool {
        self.0.barcode_multiplexing_type()
            == Some(BarcodeMultiplexingType::ReadLevel(ReadLevel::OH))
    }

    /// Return the multiplexing type as String. Returns None if not multiplexed.
    pub fn get_barcode_multiplexing_type(&self) -> PyResult<Option<String>> {
        Ok(self
            .0
            .barcode_multiplexing_type()
            .map(|x| x.as_str().to_string()))
    }
}
