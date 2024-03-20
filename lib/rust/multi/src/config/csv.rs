use super::parse::ParseCtx;
use super::scsv::{Section, Span, XtraData};
use anyhow::{anyhow, bail, Context, Result};
use itertools::Itertools;
use nom_locate::LocatedSpan;
use std::collections::HashMap;
use std::convert::TryFrom;
use std::fmt::Display;
use std::ops::Range;
use std::str::FromStr;

#[derive(Debug)]
pub struct CsvParser<'a, T: Display> {
    section: Section<'a, T>,
    required_headers: HashMap<String, usize>,
    optional_headers: HashMap<String, Option<usize>>,
}

#[derive(Clone)]
pub enum HeaderReq<'a> {
    Single(&'a str),
    Either(&'a [&'a str]),
}

pub struct HeaderReqs<'a>(Vec<HeaderReq<'a>>);

impl<'a> From<&'a [&'a str]> for HeaderReqs<'a> {
    fn from(reqs: &'a [&'a str]) -> Self {
        HeaderReqs(reqs.iter().map(|s| HeaderReq::Single(s)).collect())
    }
}

impl<'a> From<&'a [HeaderReq<'a>]> for HeaderReqs<'a> {
    fn from(reqs: &'a [HeaderReq<'a>]) -> Self {
        HeaderReqs(reqs.to_vec())
    }
}

impl<'a, T: Display> CsvParser<'a, T> {
    pub fn new(
        section: Section<'a, T>,
        required_headers: impl Into<HeaderReqs<'a>>,
        optional_headers: &[&str],
    ) -> Result<Self> {
        let Section { ref name, ref rows } = section;
        let ctx = ParseCtx::Hdr(name);
        let hdr = rows.iter().next().ok_or_else(
            #[cold]
            || anyhow!("{ctx} is empty"),
        )?;
        // Make sure that we don't have any empty cells in the header, besides
        // possible trailing empty cells.
        validate_no_intermediate_empty_cells(hdr).with_context(|| ctx.to_string())?;

        // Ensure that the header and the rest of the rows are the same length.
        // Take into account that the header and the rows might have some trailing empty cells.
        let trimmed_header_len = trimmed_len(hdr.as_slice());
        for (i, row) in rows.iter().enumerate().skip(1) {
            let trimmed_row_len = trimmed_len(row.as_slice());
            let too_short = row.len() < trimmed_header_len;
            let too_long = trimmed_row_len > trimmed_header_len;
            if too_short || too_long {
                let disp_row_len = if too_short {
                    row.len()
                } else {
                    trimmed_row_len
                };
                bail!(
                    "{} has {disp_row_len} column(s) but the header has {trimmed_header_len}",
                    ParseCtx::HdrRow(name, i + 1),
                );
            }
        }
        // TODO: use a reduce here to detect duplicate header entries
        let mut hdr: HashMap<String, usize> = hdr
            .iter()
            .enumerate()
            .map(|(i, v)| (v.fragment().to_ascii_lowercase(), i))
            .collect();
        let required_headers = required_headers
            .into()
            .0
            .iter()
            .map(|options| {
                // First detect tabs in header, likely a TSV file instead of CSV?
                if hdr.keys().any(|k| k.contains('\t')) {
                    bail!("Tab delimiter detected in header, provide a CSV file instead");
                }
                let k = match options {
                    HeaderReq::Single(k) => k,
                    HeaderReq::Either(reqs) => {
                        // Filter required keys for validation
                        let keys: Vec<_> = reqs.iter().filter(|k| hdr.contains_key(**k)).collect();
                        if keys.len() == 1 {
                            // Only one key? All good, use it later
                            *keys[0]
                        } else if keys.is_empty() {
                            bail!("Missing either one of {reqs:?} required headers")
                        } else {
                            bail!(
                                "Should have only one of {reqs:?} required header, found {keys:?}"
                            )
                        }
                    }
                };

                hdr.remove_entry(k).ok_or_else(
                    #[cold]
                    || anyhow!("{ctx} is missing required header '{k}'"),
                )
            })
            .try_collect()?;
        let optional_headers = optional_headers
            .iter()
            .map(|&k| (k.to_string(), hdr.remove(k)))
            .collect();
        for (k, _) in hdr {
            if !k.is_empty() {
                bail!("{ctx} has extra column '{k}'");
            }
        }
        Ok(CsvParser {
            section,
            required_headers,
            optional_headers,
        })
    }

    pub fn find_req(&self, row: usize, hdr: &str) -> Result<&Span<'a>> {
        let Section { ref name, ref rows } = &self.section;
        let ctx = ParseCtx::HdrRow(name, row + 1);
        if let Some(vals) = rows.get(row) {
            return self
                .required_headers
                .get(hdr)
                .ok_or_else(
                    #[cold]
                    || anyhow!("{ctx} has no required field '{hdr}'"),
                )
                .and_then(|&col| {
                    vals.get(col).ok_or_else(
                        #[cold]
                        || anyhow!("{ctx} is missing required field '{hdr}'"),
                    )
                });
        }
        bail!("[{name}] has no row {}", row + 1)
    }

    pub fn find_opt(&self, row: usize, hdr: &str) -> Result<Option<&Span<'a>>> {
        let Section { ref name, ref rows } = &self.section;
        let ctx = ParseCtx::HdrRow(name, row + 1);
        if let Some(vals) = rows.get(row) {
            return self
                .optional_headers
                .get(hdr)
                .ok_or_else(
                    #[cold]
                    || anyhow!("{ctx} has no optional field '{hdr}'"),
                )
                .map(|col| col.map(|col| vals.get(col)).flatten());
        }
        bail!("{ctx} does not exist")
    }

    #[allow(dead_code)]
    pub fn parse_req<R>(&self, row: usize, hdr: &str) -> Result<R>
    where
        R: FromStr,
        <R as FromStr>::Err: std::error::Error + Sync + Send + 'static,
    {
        let Section { ref name, .. } = &self.section;
        let ctx = ParseCtx::HdrRow(name, row + 1);
        self.find_req(row, hdr).and_then(|val| {
            val.fragment().parse::<R>().with_context(|| {
                format!(
                    "{ctx} has invalid {hdr} '{}' at line: {}, col: {}",
                    val.fragment(),
                    val.location_line(),
                    val.naive_get_utf8_column(),
                )
            })
        })
    }

    #[allow(dead_code)]
    pub fn try_into_req<R>(&self, row: usize, hdr: &str) -> Result<R>
    where
        R: TryFrom<Span<'a>>,
        <R as TryFrom<Span<'a>>>::Error: std::error::Error + Sync + Send + 'static,
    {
        self.find_req(row, hdr)
            .and_then(|val| Ok(R::try_from(val.clone())?))
    }

    #[allow(dead_code)]
    pub fn parse_opt<R>(&self, row: usize, hdr: &str) -> Result<Option<R>>
    where
        R: FromStr,
        <R as FromStr>::Err: std::error::Error + Send + Sync + 'static,
    {
        let Section { ref name, .. } = &self.section;
        let ctx = ParseCtx::HdrRow(name, row + 1);
        self.find_opt(row, hdr).and_then(|val| {
            val.map(|val| {
                val.fragment().parse::<R>().with_context(|| {
                    format!(
                        "{ctx} has invalid {hdr} '{}' at line: {}, col: {}",
                        val.fragment(),
                        val.location_line(),
                        val.naive_get_utf8_column(),
                    )
                })
            })
            .transpose()
        })
    }

    #[allow(dead_code)]
    pub fn try_into_opt<R>(&self, row: usize, hdr: &str) -> Result<Option<R>>
    where
        R: TryFrom<Span<'a>>,
        <R as TryFrom<Span<'a>>>::Error: std::error::Error + Sync + Send + 'static,
    {
        self.find_opt(row, hdr)
            .and_then(|val| val.map(|val| Ok(R::try_from(val.clone())?)).transpose())
    }

    pub fn rows(&self) -> Range<usize> {
        // skip over the hdr row
        1..self.section.rows.len()
    }
}

/// Return the length of this row after omitting all trailing empty cells.
fn trimmed_len(row: &[LocatedSpan<&str, XtraData>]) -> usize {
    row.iter().rev().skip_while(|cell| cell.is_empty()).count()
}

/// Check that the provided row has no intermediate empty cells.
/// Trailing empty cells are OK.
fn validate_no_intermediate_empty_cells(row: &[LocatedSpan<&str, XtraData>]) -> Result<()> {
    let mut first_empty_cell = None;
    for (i, cell) in row.iter().enumerate() {
        if cell.is_empty() && first_empty_cell.is_none() {
            first_empty_cell = Some(i);
            continue;
        }
        if let Some(empty_pos) = first_empty_cell {
            if !cell.is_empty() {
                bail!("intermediate empty cell at col {}", empty_pos + 1);
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::config::scsv::section_csv;

    #[test]
    fn test_intermediate_empty_cell_check() {
        let check = |row_str: &str| {
            let row: Vec<_> = row_str
                .split(',')
                .map(|cell| LocatedSpan::new_extra(cell, XtraData::new("::test")))
                .collect();
            validate_no_intermediate_empty_cells(&row)
        };
        assert!(check("foo,bar").is_ok());
        assert!(check("foo,bar,").is_ok());
        assert!(check("foo,bar,,").is_ok());
        assert!(check("foo,bar,,nope").is_err());
        assert!(check("foo,bar,,nope,,").is_err());
    }

    fn check_header_row_len_match(
        hdr: &str,
        row: &str,
        expected_header_len: usize,
        expected_row_len: usize,
    ) {
        let contents = format!("[test]\n{hdr}\n{row}\n");
        let span = Span::new_extra(&contents, XtraData::new("::test"));
        let section = section_csv(span)
            .unwrap()
            .1
            .into_iter()
            .exactly_one()
            .unwrap();
        let res = CsvParser::new(section, &["foo"][..], &["bar"]);

        if expected_row_len == expected_header_len {
            assert!(res.is_ok());
        } else {
            assert_eq!(
                format!("[test] row 2 has {expected_row_len} column(s) but the header has {expected_header_len}"),
                res.unwrap_err().to_string()
            );
        }
    }

    #[test]
    fn test_header_row_len_match() {
        check_header_row_len_match("foo", "baz", 1, 1);
        // check that we handle trailing extra cells correctly
        check_header_row_len_match("foo,,", "baz", 1, 1);
        check_header_row_len_match("foo", "baz,,", 1, 1);
        check_header_row_len_match("foo,bar", "baz", 2, 1);
        check_header_row_len_match("foo,bar,,", "baz", 2, 1);
        check_header_row_len_match("foo,bar,,", "baz,,bad", 2, 3);
        // trailing empty cells in the data row should be fine
        check_header_row_len_match("foo,bar,,", "baz,,,", 3, 3);
    }
}
