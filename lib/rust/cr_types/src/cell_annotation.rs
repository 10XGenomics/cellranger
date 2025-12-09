#![expect(missing_docs)]
use anyhow::{Error, Result};
use serde::Serialize;
use std::str::FromStr;

#[derive(Clone, Debug, Serialize)]
pub enum CellAnnotationModel {
    Auto,
    Custom(String),
}

impl CellAnnotationModel {
    /// Converts cell annotation model to ppln inputs
    /// CellAnnotationModel::Custom("my model".to_string())
    /// gets mapped to Some("my model".to_string())
    /// CellAnnotationModel::Auto gets mapped to None
    pub fn to_pipeline_inputs(&self) -> Option<String> {
        if let Self::Custom(model_name) = self {
            Some(model_name.clone())
        } else {
            None
        }
    }
}

impl FromStr for CellAnnotationModel {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self> {
        if s.to_lowercase().trim() == "auto" {
            Ok(CellAnnotationModel::Auto)
        } else {
            Ok(CellAnnotationModel::Custom(s.to_string()))
        }
    }
}
