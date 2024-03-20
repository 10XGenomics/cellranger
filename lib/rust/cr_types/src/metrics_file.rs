//! Martian filetype for representing metrics JSON files.
//! Use this type for any metrics JSON file whose contents is not represented
//! by a static type.
use anyhow::{Context, Result};
use martian::MartianRover;
use martian_derive::martian_filetype;
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::FileTypeRead;
use metric::{JsonReport, JsonReporter, TxHashMap};
use serde::{Deserialize, Serialize};
use serde_json::{Map, Value};
use std::collections::HashMap;
use std::fs::File;
use std::io::BufWriter;

martian_filetype! { MetricsFile, "json" }

impl MetricsFile {
    /// Write the provided serializable object into this JSON file.
    pub fn write(&self, value: &impl Serialize) -> Result<()> {
        Ok(serde_json::to_writer_pretty(
            BufWriter::new(File::create(self).with_context(|| self.display().to_string())?),
            value,
        )?)
    }

    /// Create a metrics file and write the provided reporter directly into it.
    /// Convenience constructor for a very common case.
    pub fn from_reporter<R: JsonReport>(
        rover: &MartianRover,
        filename: &str,
        reporter: &R,
    ) -> Result<Self> {
        let file: Self = rover.make_path(filename);
        reporter.report(&file)?;
        Ok(file)
    }
}

/// Read a metrics file as a dynamic hash map.
impl FileTypeRead<TxHashMap<String, Value>> for MetricsFile {
    fn read_from<R: std::io::Read>(reader: R) -> Result<TxHashMap<String, Value>> {
        JsonFile::read_from(reader)
    }
}

/// Read a metrics file as a dynamic hash map.
impl FileTypeRead<HashMap<String, Value>> for MetricsFile {
    fn read_from<R: std::io::Read>(reader: R) -> Result<HashMap<String, Value>> {
        JsonFile::read_from(reader)
    }
}

/// Read a metrics file as a dynamic hash map.
impl FileTypeRead<Map<String, Value>> for MetricsFile {
    fn read_from<R: std::io::Read>(reader: R) -> Result<Map<String, Value>> {
        JsonFile::read_from(reader)
    }
}

/// Read a metrics file directly into a JsonReporter.
impl FileTypeRead<JsonReporter> for MetricsFile {
    fn read_from<R: std::io::Read>(reader: R) -> Result<JsonReporter> {
        JsonFile::read_from(reader)
    }
}

/// Read a metrics file into a dynamic JSON Value.
impl FileTypeRead<Value> for MetricsFile {
    fn read_from<R: std::io::Read>(reader: R) -> Result<Value> {
        JsonFile::read_from(reader)
    }
}
