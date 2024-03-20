use itertools::Itertools;
use metric::TxHashMap;
use pyanyhow::Result;
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyDict};
use serde::{Deserialize, Serialize};

#[pyclass(module = "cellranger.fast_utils")]
#[derive(Clone, Serialize, Deserialize)]
pub struct MatrixBarcodeIndex(TxHashMap<String, usize>);

fn strip_null_suffix(bytes: &[u8]) -> &[u8] {
    let end = bytes
        .iter()
        .position(|&b| b == b'\x00')
        .unwrap_or(bytes.len());
    &bytes[0..end]
}

#[pymethods]
impl MatrixBarcodeIndex {
    #[new]
    fn new(inner: TxHashMap<String, usize>) -> Self {
        MatrixBarcodeIndex(inner)
    }

    #[staticmethod]
    fn from_raw_bytes(py: Python<'_>, barcode_bytes: Py<PyBytes>, chunk_size: usize) -> Self {
        MatrixBarcodeIndex(
            barcode_bytes
                .as_bytes(py)
                .chunks(chunk_size)
                .enumerate()
                .map(|(i, bc)| {
                    (
                        String::from_utf8(strip_null_suffix(bc).to_vec()).unwrap(),
                        i,
                    )
                })
                .collect(),
        )
    }

    fn __deepcopy__(&self, _memo: &PyDict) -> Self {
        self.clone()
    }
    fn __setstate__(&mut self, state: &PyBytes) -> Result<()> {
        *self = bincode::deserialize(state.as_bytes())?;
        Ok(())
    }
    fn __getstate__<'py>(&self, py: Python<'py>) -> Result<&'py PyBytes> {
        Ok(PyBytes::new(py, &bincode::serialize(&self)?))
    }
    fn __getnewargs__(&self) -> (TxHashMap<String, usize>,) {
        (self.0.clone(),)
    }
    fn is_equal(&self, other: &Self) -> bool {
        self.0 == other.0
    }

    fn bc_to_int(&self, _py: Python<'_>, barcode: &PyAny) -> PyResult<usize> {
        let barcode = if let Ok(bytes) = barcode.extract::<&[u8]>() {
            std::str::from_utf8(strip_null_suffix(bytes)).unwrap()
        } else if let Ok(string) = barcode.extract::<&str>() {
            string
        } else {
            panic!("barcode must be bytes or str, got {barcode:?}");
        };
        self.0.get(barcode).copied().ok_or_else(|| {
            PyErr::new::<pyo3::exceptions::PyKeyError, _>(format!(
                "Barcode {barcode} not found in index",
            ))
        })
    }

    fn bcs_to_ints(
        &self,
        py: Python<'_>,
        barcodes: Vec<&PyAny>,
        sort: bool,
    ) -> PyResult<Vec<usize>> {
        let mut indices: Vec<_> = barcodes
            .into_iter()
            .map(|barcode| self.bc_to_int(py, barcode))
            .try_collect()?;
        if sort {
            indices.sort_unstable();
        }
        Ok(indices)
    }
}
