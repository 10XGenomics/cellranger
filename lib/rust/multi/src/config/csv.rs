use super::parse::ParseCtx;
use super::scsv::{Section, Span};
use anyhow::{anyhow, bail, Context, Result};
use itertools::Itertools;
use std::collections::HashMap;
use std::convert::TryFrom;
use std::fmt::Display;
use std::ops::Range;
use std::str::FromStr;

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
