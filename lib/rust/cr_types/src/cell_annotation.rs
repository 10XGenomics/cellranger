#![expect(missing_docs)]
use anyhow::{Error, Result, bail};
use serde::Serialize;
use std::str::FromStr;

const MODELS_NOT_ALLOWED: [&str; 2] = ["azimuth", "panhuman"];

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
        } else if s.matches(',').count() > 0 {
            bail!("Please provide only one model name in --cell-annotation-model.")
        } else if MODELS_NOT_ALLOWED
            .iter()
            .any(|&model| s.to_lowercase().contains(model))
        {
            bail!(
                "Local models Azimuth-panhuman is run by default. Only enter cloud models in --cell-annotation-model."
            )
        } else {
            Ok(CellAnnotationModel::Custom(s.to_string()))
        }
    }
}
