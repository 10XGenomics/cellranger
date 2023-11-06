//! Martian filetype for representing metrics JSON files.
//! Use this type for any metrics JSON file whose contents is not represented
//! by a static type.
use anyhow::Result;
use martian::MartianRover;
use martian_derive::martian_filetype;
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::FileTypeRead;
use metric::{JsonReport, JsonReporter, TxHashMap};
use serde::{Deserialize, Serialize};
use serde_json::{Map, Value};
use std::collections::HashMap;

martian_filetype! { MetricsFile, "json" }

impl MetricsFile {
    /// Write the provided JSON report into this file.
    pub fn write_report<R: JsonReport>(&self, reporter: &R) -> Result<()> {
        reporter.report(self)
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
