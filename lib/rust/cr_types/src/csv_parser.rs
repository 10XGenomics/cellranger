use anyhow::{anyhow, bail, Context, Result};
use csv::StringRecord;
use itertools::Itertools;
use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};
use std::str::FromStr;

/// Helper class for parsing configuration CSV files, validating them & providing good error messages
pub struct CsvParser {
    filetype: String,
    filename: PathBuf,
    headers: Vec<String>,
    rows: Vec<StringRecord>,
    col_map: HashMap<String, usize>,
    line: usize,
}

impl CsvParser {
    /// Create a CSV parser / validator for file `filename`. `required_headers` are checked and an error will be
    /// returned if they're not present. `filetype` is a readable description of the kind of CSV file being parsed and
    /// will be used in error messages.
    pub fn new<T: AsRef<str>>(
        filename: &Path,
        required_headers: impl IntoIterator<Item = T>,
        filetype: &str,
    ) -> Result<CsvParser> {
        let file = File::open(filename).with_context(|| filename.display().to_string())?;
        let buf_rdr = BufReader::new(file);
        let mut rdr = csv::Reader::from_reader(buf_rdr);

        let mut headers = rdr.headers()?.clone();
        headers.trim();
        let headers: Vec<_> = headers.iter().map(String::from).collect();

        let mut rows = Vec::new();
        for result in rdr.records() {
            let mut record = result?;
            record.trim();
            rows.push(record);
        }

        let col_map = CsvParser::check_headers(filename, required_headers, &headers)?;

        Ok(CsvParser {
            filetype: filetype.to_string(),
            filename: filename.to_path_buf(),
            headers,
            rows,
            col_map,
            line: 0,
        })
    }

    /// Number of data rows
    pub fn len(&self) -> usize {
        self.rows.len()
    }
    pub fn is_empty(&self) -> bool {
        self.rows.is_empty()
    }

    /// File name
    pub fn filename(&self) -> &Path {
        &self.filename
    }

    /// Headers found in CSV
    pub fn headers(&self) -> &[String] {
        &self.headers
    }

    /// Set the line number (not including the header) to pull data from
    pub fn set_line(&mut self, line: usize) {
        self.line = line;
    }

    /// Get map of extra key-value pairs in CSV in columns
    /// that aren't listed in `ignore_cols`.
    pub fn get_extra_data(&self) -> BTreeMap<String, String> {
        let mut res = BTreeMap::new();

        for h in &self.headers {
            res.insert(h.clone(), self.get_string(h));
        }

        res
    }

    /// Get a value of type `T` from column `col` on the current line.
    /// Returns None is the column doesn't exist. Returns an error if the
    /// contents of the csv cannot be parsed as `T`. Panics if the field doesn't exist.
    /// Returns Ok(None) for an empty field
    pub fn try_parse_field<T>(&self, col: &str, expected: &str) -> Result<Option<T>>
    where
        T: FromStr,
        Result<T, <T as FromStr>::Err>: anyhow::Context<T, <T as FromStr>::Err>,
    {
        let i = self.col_map[col];
        let v = self.rows[self.line][i].trim();
        if v.is_empty() {
            return Ok(None);
        }

        Ok(Some(v.parse::<T>().with_context(|| {
            format!(
                "Error in {} file '{}'. On line {} in '{col}' column: \
                 Expected a {expected} but received '{v}'",
                self.filetype,
                self.filename.display(),
                self.line
            )
        })?))
    }

    /// Get a value of type `T` from column `col` on the current line.
    /// Returns None is the column doesn't exist. Returns an error if the
    /// contents of the csv cannot be parsed as `T`. Panics if the field doesn't exist.
    /// Returns Ok(None) for an empty field
    pub fn parse_field<T>(&self, col: &str, expected: &str) -> Result<T>
    where
        T: FromStr,
        Result<T, <T as FromStr>::Err>: anyhow::Context<T, <T as FromStr>::Err>,
    {
        match self.try_parse_field(col, expected) {
            Err(e) => Err(e),
            Ok(Some(v)) => Ok(v),
            Ok(None) => {
                bail!(
                    "Error in {} file '{}'. On line {} in '{col}' column: \
                     Expected a {expected}, but got empty value",
                    self.filetype,
                    self.filename.display(),
                    self.line,
                );
            }
        }
    }

    /// Get a usize from column `col` on the current line.
    pub fn try_parse_many<T>(&self, col: &str, expected: &str) -> Result<Option<Vec<T>>>
    where
        T: FromStr,
        Result<T, <T as FromStr>::Err>: anyhow::Context<T, <T as FromStr>::Err>,
    {
        let col = self.col_map.get(col);
        if col.is_none() {
            return Ok(None);
        }

        let col = col.unwrap();
        let val = self.rows[self.line][*col].trim();

        let xs: Vec<_> = val
            .split('|')
            .filter_map(|v| {
                let v = v.trim();
                if v.is_empty() {
                    return None;
                }
                Some(v.parse::<T>().with_context(|| {
                    format!(
                        "Error in {} file: {}. On line {} in {col} column: \
                         Expected a {expected} but received '{v}'",
                        self.filetype,
                        self.filename.display(),
                        self.line
                    )
                }))
            })
            .try_collect()?;
        Ok(Some(xs))
    }

    /// Attempt to parse a set of integers. Accepted forms are "8", "1|2", "1-4".
    pub fn try_parse_int_range(&self, col: &str, expected: &str) -> Result<Option<Vec<usize>>> {
        let Some(&icol) = self.col_map.get(col) else {
            return Ok(None);
        };
        let vals: Vec<&str> = self.rows[self.line][icol].split('-').collect();

        // attempt entries of the form '1-4' and convert into a vec of integers
        if vals.len() == 2 {
            let v1 = vals[0].trim();
            let v2 = vals[1].trim();

            if let Ok(start) = v1.parse::<usize>() {
                if let Ok(end) = v2.parse::<usize>() {
                    let range: Vec<usize> = (start..(end + 1)).collect();
                    return Ok(Some(range));
                }
            }
        }

        // fallback to
        self.try_parse_many(col, expected)
    }

    /// Get a string from column `col` on the current line.
    /// Returns an error on an empty string, Panics if column doesn't exist
    pub fn require_string(&self, col: &str) -> Result<String> {
        self.try_get_string(col).ok_or_else(|| {
            anyhow!(
                "Error in {} file '{}'. On line {} in '{col}' column: \
                 Value required but cell is empty.",
                self.filetype,
                self.filename.display(),
                self.line,
            )
        })
    }

    /// Get a string from column `col` on the current line.
    /// Returns "" for an empty string, Panics if column doesn't exist
    pub fn get_string(&self, col: &str) -> String {
        self.try_get_string(col).unwrap_or_default()
    }

    /// Get a string from column `col` on the current line.
    /// Returns None for an empty string, Panics if column doesn't exist
    pub fn try_get_string(&self, col: &str) -> Option<String> {
        let col = self.col_map.get(col);
        let col = col.unwrap();

        let val = self.rows[self.line][*col].trim().to_string();

        if val.as_str() == "" {
            None
        } else {
            Some(val)
        }
    }

    fn check_headers<T: AsRef<str>>(
        file_arg: &Path,
        required: impl IntoIterator<Item = T>,
        headers: &[String],
    ) -> Result<HashMap<String, usize>> {
        let mut result = HashMap::new();

        // check that we have required headers
        for r in required {
            if !headers.contains(&r.as_ref().into()) {
                bail!(
                    "The input file '{}' must contain a column named '{}', but it was not found. \
                    Please check the headers in the CSV file.",
                    file_arg.display(),
                    r.as_ref()
                );
            }
        }

        // column name to column index map
        for (i, h) in headers.iter().enumerate() {
            result.insert(h.to_string(), i);
        }

        Ok(result)
    }
}
