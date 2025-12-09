//! Module for supporting pipeline - 10x Cloud Analysis interaction.
#![deny(missing_docs)]

use anyhow::Result;
use std::env;
use std::path::Path;

/// An error message
pub const CELL_ANNOTATION_HOMEDIR_MSG: &str = "Could not find a 10x cloud token in the default user directory.  In order to enable cell annotation, run cellranger cloud auth setup, or please supply a --tenx-cloud-token-path argument.";

/// Return the path to the default cloud token
pub fn default_token_path() -> Result<String> {
    let homedir = env::var("HOME")?;
    let credentials_path = format!("{homedir}/.config/txg/credentials");
    let credentials_path_obj = Path::new(credentials_path.as_str());
    if !credentials_path_obj.exists() {
        return Err(anyhow::anyhow!("Default credentials file not found."));
    }
    Ok(credentials_path)
}
