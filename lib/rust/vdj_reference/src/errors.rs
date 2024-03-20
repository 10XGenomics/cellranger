#![allow(clippy::enum_glob_use)]
use anyhow::Result;
use bio::io::fasta::Record;
use itertools::Itertools;
use std::fmt::{self, Write};
use std::fs::File;
use std::io::{BufRead, BufReader, Lines};
use std::iter::Enumerate;
use std::ops::Range;
use std::path::{Path, PathBuf};
pub struct ErrorContext {
    reader: Enumerate<Lines<BufReader<File>>>,
    pub state: FastaState,
}

impl ErrorContext {
    pub fn new(fa_file: &Path) -> Result<Self> {
        Ok(ErrorContext {
            reader: BufReader::new(File::open(fa_file)?).lines().enumerate(),
            state: FastaState::NoRecordProcessed,
        })
    }
    pub fn advance(&mut self, record: &Record) -> Result<()> {
        let current_header = match record.desc() {
            Some(desc) => format!(">{} {desc}", record.id()),
            None => format!(">{}", record.id()),
        };

        self.state = loop {
            match self.reader.next() {
                Some((line_num, line)) => {
                    let header = line?;
                    if header.trim_end() == current_header.trim_end() {
                        break FastaState::ProcessedRecord { line_num, header };
                    }
                }
                None => break FastaState::Unknown, // Should be unreachable!
            }
        };
        Ok(())
    }
}

#[derive(Clone, Debug)]
pub enum FastaState {
    NoRecordProcessed,
    ProcessedRecord { line_num: usize, header: String },
    Unknown,
}

impl FastaState {
    fn is_processed(&self) -> bool {
        matches!(self, &FastaState::ProcessedRecord { .. })
    }
    pub fn line_num(&self) -> Option<usize> {
        if let FastaState::ProcessedRecord { line_num, .. } = self {
            Some(*line_num)
        } else {
            None
        }
    }
}

impl fmt::Display for FastaState {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            FastaState::NoRecordProcessed => {
                write!(f, "No fasta records were processed successfully!")
            }
            FastaState::ProcessedRecord { line_num, header } => {
                write!(f, "LINE {:5}: {header}", line_num + 1)
            }
            FastaState::Unknown => Ok(()),
        }
    }
}

fn _last_record_print(last: Option<&(usize, String)>) -> String {
    match last.as_ref() {
        Some((n, line)) => {
            format!(
                "The last successfully processed record was:\nLine {:5}:{}",
                n + 1,
                line
            )
        }
        None => "No fasta records were processed successfully!".into(),
    }
}

#[derive(Debug)]
pub enum VdjReferenceErrors {
    CannotReadFastaRecord {
        fa_file: PathBuf,
        state: FastaState,
        error: std::io::Error,
    },
    InvalidBaseInSequence {
        fa_file: PathBuf,
        state: FastaState,
        base: char,
        position: usize,
    },
    InvalidHeader {
        fa_file: PathBuf,
        errors: Vec<(FastaState, HeaderErrors)>,
    },
}

impl fmt::Display for VdjReferenceErrors {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        use VdjReferenceErrors::*;

        fn _common_error(fa_file: &Path) -> String {
            format!(
                "Error reading the fasta file \"{}\" due to:",
                fa_file.display()
            )
        }

        let msg = match self {
            CannotReadFastaRecord {
                fa_file,
                state,
                error,
            } => {
                let state_msg = if state.is_processed() {
                    format!("The last successfully processed record is:\n{state}\n")
                } else {
                    state.to_string()
                };
                format!(
                    "{}\n\n- {}\n\n{}\n",
                    _common_error(fa_file),
                    error,
                    state_msg
                )
            }
            InvalidBaseInSequence {
                fa_file,
                state,
                base,
                position,
            } => {
                format!(
                    "{}\n\n\
                    - Invalid character '{}' in the sequence at position {}. Only {} characters are allowed. \
                    Please correct the sequence associated with the following header:\n\n{}\n",
                    _common_error(fa_file),
                    base,
                    position,
                    crate::ALLOWED_NUCLEOTIDES,
                    state
                )
            }
            InvalidHeader { fa_file, errors } => {
                let num_errors = errors.len();
                let sep = "....................";

                let error_msg = errors
                    .iter()
                    .enumerate()
                    .map(|(i, (state, err))| {
                        format!(
                            "{}\nError {}/{}\n{}\n\n{}\n{}\n",
                            sep,
                            i + 1,
                            num_errors,
                            sep,
                            state,
                            match state {
                                FastaState::ProcessedRecord { header, .. } =>
                                    err.with_highlight(header),
                                _ => err.to_string(),
                            }
                        )
                    })
                    .join("\n");

                format!("{}\n\n{error_msg}", _common_error(fa_file))
            }
        };

        write!(f, "{msg}")
    }
}

impl std::error::Error for VdjReferenceErrors {}

#[derive(Debug)]
pub enum HeaderErrors {
    NoDescription,
    UnexpectedIdFormat {
        num_parts: usize,
        id: String,
    },
    UnexpectedDescFormat {
        num_parts: usize,
        desc: String,
    },
    FeatureIdNotAnInteger {
        feature_id: String,
    },
    UnknownVariant {
        range: Range<usize>,
        parse_error: String,
    },
    GeneNameDoesNotStartWithChainName {
        range: Range<usize>,
        gene_name: String,
        chain_name: String,
    },
    GeneNameNotAllowed {
        range: Range<usize>,
        gene_name: String,
    },
    GeneNameAmbiguous {
        range: Range<usize>,
        gene_name: String,
        exclusive_options: Vec<&'static str>,
    },
    DuplicateId {
        id: u32,
        last_line: Option<usize>,
    },
}

impl HeaderErrors {
    fn with_highlight(&self, header: &str) -> String {
        let mut result = String::new();
        let hl_range = self.highlight_range(header);
        const LINE_NUM_PREFIX_LEN: usize = 12;
        for _ in 0..LINE_NUM_PREFIX_LEN {
            result.push(' ');
        }
        for i in 0..header.len() {
            result.push(if hl_range.contains(&i) { '^' } else { ' ' });
        }
        result.push('\n');
        for _ in 0..LINE_NUM_PREFIX_LEN + hl_range.start {
            result.push(' ');
        }
        write!(&mut result, "{self}").unwrap();
        result
    }

    fn highlight_range(&self, header: &str) -> Range<usize> {
        use HeaderErrors::*;
        let n = header.len();
        match self {
            NoDescription => 0..n,
            UnexpectedIdFormat { id, .. } => 0..id.len() + 1,
            UnexpectedDescFormat { desc, .. } => n - desc.len()..n,
            FeatureIdNotAnInteger { feature_id } => 1..feature_id.len() + 1,
            UnknownVariant { range, .. } => range.clone(),
            GeneNameDoesNotStartWithChainName { range, .. } => range.clone(),
            GeneNameNotAllowed { range, .. } => range.clone(),
            GeneNameAmbiguous { range, .. } => range.clone(),
            DuplicateId { id, .. } => 1..id.to_string().len() + 1,
        }
    }
}

impl std::error::Error for HeaderErrors {}

impl fmt::Display for HeaderErrors {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        use HeaderErrors::*;
        let msg = match self {
            NoDescription => {
                "Expected two strings (id & description) separated by a space in the FASTA header. Found no description".to_string()
            }
            UnexpectedIdFormat { num_parts, id } => {
                format!(
                    "Record ID in the FASTA header must consist of 2 fields \
                    ('feature_id', 'display_name') separated by a '|' (for e.g '10|IGLV5-52'). \
                    Found {num_parts} values in {id}"
                )
            }
            UnexpectedDescFormat { num_parts, desc } => {
                format!(
                    "Record description in the FASTA header must consist of 7 fields \
                    ('record_id', 'gene_name', 'region_type', 'chain_type', 'chain', 'isotype', \
                    'allele_name') separated by '|' (for e.g 'ENST00000390539|IGHA2|C-REGION|IG|IGH|A2|00'). \
                    Found {num_parts} values in '{desc}'"
                )
            }
            FeatureIdNotAnInteger { feature_id } => {
                format!("Expected an integer for feature_id. Found '{feature_id}'")
            }
            UnknownVariant { parse_error, .. } => parse_error.clone(),
            GeneNameDoesNotStartWithChainName {
                gene_name,
                chain_name,
                ..
            } => {
                format!(
                    "The gene name '{gene_name}' does not begin with the chain name '{chain_name}'."
                )
            }
            GeneNameNotAllowed { gene_name, .. } => {
                format!(
                    "The gene name '{}' contains one of the disallowed words: [{}].",
                    gene_name,
                    crate::GENE_NAME_DISALLOWED.join(", ")
                )
            }
            GeneNameAmbiguous {
                gene_name,
                exclusive_options,
                ..
            } => {
                format!(
                    "We found all of the phrases [{}] in the gene name '{}'. \
                    Only one such phrase is allowed as a substring in the gene name.",
                    exclusive_options.iter().join(", "),
                    gene_name
                )
            }
            DuplicateId { id, last_line } => format!(
                "The feature ID '{}' is duplicated. {}",
                id,
                last_line.map_or(String::new(), |line| format!(
                    "The same ID was encountered in line {}.",
                    line + 1
                ))
            ),
        };
        write!(f, "{msg}")
    }
}
