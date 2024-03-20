use pyo3::exceptions::PyException;
use pyo3::PyErr;

pub type Result<T, E = PyAnyhowError> = core::result::Result<T, E>;

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
