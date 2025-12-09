//! pyanyhow
#![deny(missing_docs)]

use pyo3::PyErr;
use pyo3::exceptions::PyException;

/// A Result of PyAnyhowError
pub type Result<T, E = PyAnyhowError> = core::result::Result<T, E>;

/// Wrap a anyhow::Error
#[derive(Debug)]
pub struct PyAnyhowError(pub anyhow::Error);

impl From<PyAnyhowError> for PyErr {
    fn from(error: PyAnyhowError) -> Self {
        PyException::new_err(format!("{:#}", error.0))
    }
}

impl<T> From<T> for PyAnyhowError
where
    anyhow::Error: From<T>,
{
    fn from(error: T) -> Self {
        Self(error.into())
    }
}
