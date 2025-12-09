#![expect(missing_docs)]
use crate::config::csv::CsvParser;
use crate::config::parse::{Parse, ParseCtx};
use crate::config::scsv::{Section, Span, XtraData, plain_csv};
use crate::config::{Ident, MultiConfigCsv, multiconst, samplesconst};
use anyhow::{Context, Result, anyhow, bail};
use itertools::FoldWhile::{Continue, Done};
use itertools::Itertools;
use metric::{TxHashMap, TxHashSet};
use nom_locate::LocatedSpan;
use std::convert::TryFrom;
use std::fmt::Display;
use std::fs::File;
use std::io::BufWriter;
use std::path::Path;
use std::str::FromStr;
use std::string::ToString;

/// A Barcode with sequence and gem group
#[derive(Debug)]
struct Barcode((String, u16));

impl FromStr for Barcode {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self> {
        let mut parts = s.split('-');
        let barcode = parts.next().unwrap().to_string();
        if !barcode.trim_matches(&['A', 'C', 'G', 'T'][..]).is_empty() {
            bail!("invalid barcode: {s}");
        }
        let gem_group = parts.next().unwrap_or("1").parse::<u16>()?;
        if gem_group != 1 {
            bail!("invalid gem group {gem_group}, must be 1");
        }
        if parts.next().is_some() {
            bail!("invalid barcode: {s}");
        }
        Ok(Barcode((barcode, gem_group)))
    }
}

impl Display for Barcode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}-{}", self.0.0, self.0.1)
    }
}

#[derive(Debug)]
enum SampleAssignment {
    Sample(String),
    SampleTag(String, String),
    Unassigned,
    Multiplet,
    Blank,
}

impl SampleAssignment {
    pub fn is_non_singlet(&self) -> bool {
        use SampleAssignment::{Sample, SampleTag};
        !matches!(self, Sample(_) | SampleTag(_, _))
    }
}

impl Display for SampleAssignment {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        use SampleAssignment::{Blank, Multiplet, Sample, SampleTag, Unassigned};
        match self {
            Sample(sample_id) => f.write_str(sample_id),
            SampleTag(sample_id, _) => f.write_str(sample_id),
            Blank => f.write_str("Blank"),
            Multiplet => f.write_str("Multiplet"),
            Unassigned => f.write_str("Unassigned"),
        }
    }
}

#[derive(Debug)]
struct SampleAssignmentRow {
    barcode: Barcode,
    assignment: SampleAssignment,
}

#[derive(Debug)]
pub struct SampleAssignmentCsv {
    rows: Vec<SampleAssignmentRow>,
    cmo_sample_map: TxHashMap<String, String>,
    samples: TxHashSet<String>,
}

impl SampleAssignmentCsv {
    pub fn from_file(path: &Path, cfg: &MultiConfigCsv) -> Result<Self> {
        let buf = std::fs::read_to_string(path)?;
        let s = if buf.starts_with('\u{feff}') {
            &buf.as_str()[3..]
        } else {
            buf.as_str()
        };
        let input = Span::new_extra(s, path.into());
        let name = format!("{}", path.display());

        Self::from_span(input, name, cfg).with_context(|| {
            anyhow!(
                "Error parsing barcode sample assignment file '{}'",
                path.display()
            )
        })
    }

    pub fn from_span<'a, X: Into<LocatedSpan<&'a str, XtraData>>>(
        span: X,
        name: String,
        cfg: &MultiConfigCsv,
    ) -> Result<Self> {
        let (_, sec) = plain_csv(name.clone(), span.into()).map_err(|e| match e {
            nom::Err::Error(e) | nom::Err::Failure(e) => anyhow!(
                "failed to parse barcode sample assignment CSV {} at line: {}, col: {}",
                e.input.extra,
                e.input.location_line(),
                e.input.naive_get_utf8_column()
            ),
            nom::Err::Incomplete(_) => anyhow!(
                "failed to parse barcode sample assignment CSV {name}, \
                 incomplete information available to pinpoint error",
            ),
        })?;
        SampleAssignmentCsv::try_from((&sec, cfg))
    }

    pub fn to_sample_barcodes_json(&self, path: &Path) -> Result<()> {
        use SampleAssignment::{Sample, SampleTag};
        let sample_barcodes = self.samples.iter().fold(
            TxHashMap::<String, Vec<String>>::default(),
            |mut acc, sample_id| {
                acc.insert(sample_id.clone(), vec![]);
                acc
            },
        );
        let sample_barcodes =
            self.rows
                .iter()
                .fold(sample_barcodes, |mut acc, row| match row.assignment {
                    Sample(ref sample_id) | SampleTag(ref sample_id, _) => {
                        acc.entry(sample_id.clone())
                            .or_default()
                            .push(row.barcode.to_string());
                        acc
                    }
                    _ => acc,
                });
        let mut writer = BufWriter::new(File::create(path)?);
        serde_json::to_writer_pretty(&mut writer, &sample_barcodes)?;
        Ok(())
    }

    pub fn to_cell_barcodes_json(&self, path: &Path) -> Result<()> {
        let barcodes: Vec<_> = self
            .rows
            .iter()
            .map(|row| row.barcode.to_string())
            .collect();
        let mut writer = BufWriter::new(File::create(path)?);
        serde_json::to_writer_pretty(&mut writer, &barcodes)?;
        Ok(())
    }

    pub fn to_non_singlet_barcodes_json(&self, path: &Path) -> Result<()> {
        let assignments = {
            let mut asgmts = TxHashMap::<String, Vec<String>>::default();
            asgmts.insert(SampleAssignment::Blank.to_string(), vec![]);
            asgmts.insert(SampleAssignment::Multiplet.to_string(), vec![]);
            asgmts.insert(SampleAssignment::Unassigned.to_string(), vec![]);
            asgmts
        };
        let assignments = self.rows.iter().fold(assignments, |mut acc, row| {
            if row.assignment.is_non_singlet() {
                acc.entry(row.assignment.to_string())
                    .or_default()
                    .push(row.barcode.to_string());
            }
            acc
        });
        let mut writer = BufWriter::new(File::create(path)?);
        serde_json::to_writer_pretty(&mut writer, &assignments)?;
        Ok(())
    }
    pub fn to_cells_per_tag_json(&self, path: &Path) -> Result<bool> {
        use SampleAssignment::{Multiplet, Sample, SampleTag};
        let cells_per_tag = self.cmo_sample_map.iter().fold(
            TxHashMap::<String, Vec<String>>::default(),
            |mut acc, (tag, _)| {
                acc.insert(tag.clone(), vec![]);
                acc
            },
        );
        let cells_per_tag = self
            .rows
            .iter()
            .fold_while(Some(cells_per_tag), |mut acc, row| match row.assignment {
                SampleTag(_, ref tag) => {
                    if let Some(acc) = acc.as_mut() {
                        acc.entry(tag.clone())
                            .or_insert_with(Vec::new)
                            .push(row.barcode.to_string());
                    }
                    Continue(acc)
                }
                Sample(_) => Done(None),
                Multiplet => {
                    // we don't know which tags a multiplet was assigned to,
                    // so assign each multiplet to every tag
                    if let Some(acc) = acc.as_mut() {
                        for x in acc.values_mut() {
                            x.push(row.barcode.to_string());
                        }
                    }
                    Continue(acc)
                }
                _ => Continue(acc),
            })
            .into_inner();
        if let Some(cells_per_tag) = cells_per_tag {
            let mut writer = BufWriter::new(File::create(path)?);
            serde_json::to_writer_pretty(&mut writer, &cells_per_tag)?;
            Ok(true)
        } else {
            Ok(false)
        }
    }
}

mod hdrs {
    use crate::config::csv::HeaderReq;

    pub(super) const BARCODE: &str = "barcode";
    pub(super) const SAMPLE_ID: &str = "sample_id";
    pub(super) const ASSIGNMENT: &str = "assignment";
    pub(super) const REQ_HDRS: &[HeaderReq<'_>] = &[
        HeaderReq::Single(BARCODE),
        HeaderReq::Either(&[SAMPLE_ID, ASSIGNMENT]),
    ];
    pub(super) const OPT_HDRS: &[&str] = &[];
}

impl<'a> TryFrom<(&Section<'a, String>, &MultiConfigCsv)> for SampleAssignmentCsv {
    type Error = anyhow::Error;

    fn try_from((sec, cfg): (&Section<'a, String>, &MultiConfigCsv)) -> Result<Self> {
        use hdrs::{ASSIGNMENT, BARCODE, OPT_HDRS, REQ_HDRS, SAMPLE_ID};

        let samples = cfg.samples.as_ref().ok_or_else(
            #[cold]
            || {
                anyhow!(
                    "[{}] barcode-sample-assignment provided but no samples defined in [{}]",
                    multiconst::GENE_EXPRESSION,
                    multiconst::SAMPLES,
                )
            },
        )?;
        let (sample_map, sample_barcode_id) = match (
            samples.has_cmo_ids(),
            samples.has_hashtag_ids(),
        ) {
            (true, false) => (samples.get_cmo_sample_map(), samplesconst::CMO_IDS),
            (false, true) => (samples.get_hashtag_sample_map(), samplesconst::HASHTAG_IDS),
            (_, _) => bail!(
                "[{}] barcode-sample-assignment requires samples definition in [{}] to use the mutually exclusive {} or {}",
                multiconst::GENE_EXPRESSION,
                multiconst::SAMPLES,
                samplesconst::CMO_IDS,
                samplesconst::HASHTAG_IDS
            ),
        };
        let samples: TxHashSet<_> = sample_map.values().cloned().collect();
        let parser = CsvParser::new(sec.clone(), REQ_HDRS, OPT_HDRS)?;
        let mut rows = vec![];
        let hdr = &sec.name;
        for row in parser.rows() {
            let ctx = ParseCtx::HdrRow(hdr, row + 1);
            let barcode = parser
                .find_req(row, BARCODE)?
                .parse::<Barcode>(ctx.with_col(BARCODE))?;

            let sample_id = parser.find_req(row, SAMPLE_ID).map_or_else(
                |_| Ok(None),
                |si| {
                    si.parse::<Ident>(ctx.with_col(SAMPLE_ID))
                        .map(String::from)
                        .map(Some)
                },
            )?;

            let assignment = parser.find_req(row, ASSIGNMENT).map_or_else(
                |_| Ok(None),
                |si| {
                    si.parse::<Ident>(ctx.with_col(ASSIGNMENT))
                        .map(String::from)
                        .map(Some)
                },
            )?;

            let assignment = match (sample_id, assignment) {
                (Some(sample_id), None) => {
                    let lc_sample_id = sample_id.to_ascii_lowercase();
                    match lc_sample_id.as_str() {
                        "blank" => SampleAssignment::Blank,
                        "multiplet" => SampleAssignment::Multiplet,
                        "unassigned" => SampleAssignment::Unassigned,
                        _ => {
                            if !samples.contains(&sample_id) {
                                bail!(
                                    "invalid {} on row {}, sample {} not found in [samples]",
                                    SAMPLE_ID,
                                    row + 1,
                                    sample_id
                                );
                            }
                            SampleAssignment::Sample(sample_id)
                        }
                    }
                }
                (None, Some(assignment)) => {
                    let lc_assignment = assignment.to_ascii_lowercase();
                    match lc_assignment.as_str() {
                        "blank" => SampleAssignment::Blank,
                        "multiplet" => SampleAssignment::Multiplet,
                        "unassigned" => SampleAssignment::Unassigned,
                        _ => {
                            let sample_id = sample_map
                                .get(&assignment)
                                .ok_or_else(
                                    #[cold]
                                    || {
                                        anyhow!(
                                            "invalid {} on row {}, {} {} not found in [samples]",
                                            ASSIGNMENT,
                                            row + 1,
                                            sample_barcode_id,
                                            assignment
                                        )
                                    },
                                )?
                                .clone();
                            SampleAssignment::SampleTag(sample_id, assignment)
                        }
                    }
                }
                _ => {
                    bail!("exactly one column of {SAMPLE_ID} or {ASSIGNMENT} must be provided!");
                }
            };
            rows.push(SampleAssignmentRow {
                barcode,
                assignment,
            });
        }

        // assert no duplicated barcodes are present
        let duplicated_barcodes: Vec<_> = rows
            .iter()
            .map(|r| r.barcode.to_string())
            .duplicates()
            .collect();
        if !duplicated_barcodes.is_empty() {
            bail!("Barcodes with multiple assignments found: {duplicated_barcodes:?}");
        }

        Ok(SampleAssignmentCsv {
            rows,
            samples,
            cmo_sample_map: sample_map,
        })
    }
}

#[cfg(test)]
mod tests {
    use crate::barcode_sample_assignment::SampleAssignmentCsv;
    use crate::config::MultiConfigCsv;
    use crate::config::scsv::{Span, XtraData};
    use anyhow::Result;

    const BASE_CONFIG: &str = r"
[gene-expression]
ref,/path/to/gex/ref
barcode-sample-assignment,sample_bc_assignment_disallowed_sample_id.csv
create-bam,true

[libraries]
fastq_id,fastqs,lanes,physical_library_id,feature_types,subsample_rate
tiny_gex,/path_to_fastqs,any,gex,Gene Expression,
tiny_cmo,/path_to_fastqs,any,cmo,Multiplexing Capture,

[samples]
sample_id,cmo_ids,description
raji,CMO301,raji
jurkat,CMO302,jurkat
";

    #[test]
    fn incorrect_sample_id() -> Result<()> {
        let xtra = XtraData::new("tests::incorrect_sample_id");
        let exp = MultiConfigCsv::from_reader(BASE_CONFIG.as_bytes(), xtra)?;

        let sample_assign = r"Barcode,Sample_ID
AAAGGTAGTGTTAGCT-1,jurkat
AACCATGAGCACACCC-1,jurkat
AAGATAGCAGAGATGC-1,raji raji
";
        let input = Span::new_extra(sample_assign, "tests::incorrect_sample_id".into());
        let name = "tests::incorrect_sample_id".into();
        let exp = SampleAssignmentCsv::from_span(input, name, &exp);
        assert!(exp.is_err());
        let err_msg = exp.unwrap_err().to_string();
        assert!(err_msg.contains("tests::incorrect_sample_id row 4 has invalid sample_id 'raji raji' at line: 4, col: 20: invalid character(s): ' '"),
        "{}", err_msg);

        Ok(())
    }

    #[test]
    fn duplicate_barcodes() -> Result<()> {
        let xtra = XtraData::new("tests::duplicate_barcodes");
        let exp = MultiConfigCsv::from_reader(BASE_CONFIG.as_bytes(), xtra)?;

        let sample_assign = r"Barcode,Sample_ID
AAAGGTAGTGTTAGCT-1,jurkat
AAAGGTAGTGTTAGCT-1,jurkat
AAGATAGCAGAGATGC-1,raji
";
        let input = Span::new_extra(sample_assign, "tests::duplicate_barcodes".into());
        let name = "tests::duplicate_barcodes".into();
        let exp = SampleAssignmentCsv::from_span(input, name, &exp);
        assert!(exp.is_err());
        let err_msg = exp.unwrap_err().to_string();
        assert!(
            err_msg.contains(r#"Barcodes with multiple assignments found: ["AAAGGTAGTGTTAGCT-1"]"#),
            "{}",
            err_msg
        );

        Ok(())
    }

    #[test]
    fn invalid_gem_group() -> Result<()> {
        let xtra = XtraData::new("tests::invalid_gem_group");
        let exp = MultiConfigCsv::from_reader(BASE_CONFIG.as_bytes(), xtra)?;

        let sample_assign = r"Barcode,Sample_ID
AAAGGTAGTGTTAGCT-1,jurkat
AAGATAGCAGAGATGC-2,raji
";
        let input = Span::new_extra(sample_assign, "tests::invalid_gem_group".into());
        let name = "tests::invalid_gem_group".into();
        let exp = SampleAssignmentCsv::from_span(input, name, &exp);
        assert!(exp.is_err());
        let err_msg = exp.unwrap_err().to_string();
        assert!(
            err_msg.contains(r"invalid gem group 2, must be 1"),
            "{}",
            err_msg
        );

        Ok(())
    }
}
