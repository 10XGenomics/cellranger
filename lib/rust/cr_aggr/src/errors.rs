use crate::parse_aggr_csv::MULTI_OUTS;
use itertools::Itertools;
use std::path::PathBuf;

#[derive(Debug, PartialEq, Eq, thiserror::Error)]
pub enum FieldResolutionErrors {
    #[error(
        "Could not resolve path {path:?}. The paths specified in the CSV needs to be one of:\n\
        - An absolute path\n\
        - A path relative to the pipestance root {pipestance_root:?}\n\
        - A path relative to the folder containing the aggregation csv {agg_csv_path:?}"
    )]
    RelPathDoesntExist {
        path: PathBuf,
        agg_csv_path: PathBuf,
        pipestance_root: PathBuf,
    },

    #[error("Expected a folder, but got a file: {path:?}")]
    ExpectedFolderGotFile { path: PathBuf },

    #[error("Expected a file, but got a folder: {path:?}")]
    ExpectedFileGotFolder { path: PathBuf },

    #[error("The path {path:?} does not exist.")]
    AbsPathDoesntExist { path: PathBuf },

    #[error(
        "Found an invalid character {invalid} in the text {text} at position {position}. The allowed \
         characters are upper case letters A-Z, lower case letters a-z, numbers 0-9 \
         and the following special characters {allowed}"
    )]
    InvalidTextCharacter {
        invalid: char,
        position: usize,
        text: String,
        allowed: &'static str,
    },

    #[error(
        "Could not resolve multi '{MULTI_OUTS}' path. A correct multi '{MULTI_OUTS}' path \
         needs to have at least one of the following files present: \
         path/count/molecule_info.h5, path/vdj_b/vdj_contig_info.pb, or \
         path/vdj_b/vdj_contig_info.pb"
    )]
    CouldNotResolveMultiPath {},

    // Specified top level multi out instead of sample outs in CR 6.0 or greater
    #[error(
        "You have specified the top level outs folder generated using the multi pipeline. \
         Please specify one of the sample output folders inside the 'outs/per_sample_outs' \
         folder."
    )]
    TopLevelMultiOut,
}

#[derive(Debug, thiserror::Error)]
pub enum ParseAggrCsvErrors {
    #[error("The provided CSV file is empty: {path:?}")]
    EmptyCsv { path: PathBuf },

    #[error(
        "Duplicate entry detected in column '{col}' of the CSV file {path:?}. \
            The column '{col}' needs to contain unique values, however the value '{duplicate_value}' in line {line} \
            was encountered already."
    )]
    DuplicateInUniqueColumn {
        path: PathBuf,
        col: String,
        line: usize,
        duplicate_value: String,
    },

    #[error("Parsing column '{col}' in line {line} of the CSV file {path:?}")]
    ErrorParsingField {
        path: PathBuf,
        col: String,
        line: usize,
        source: anyhow::Error,
    },

    #[error(
        "The column '{col}' in line {line} of the CSV file {path:?} is empty, which is not allowed."
    )]
    EmptyField {
        path: PathBuf,
        col: String,
        line: usize,
    },

    #[error(
        "The column at position {col} of the CSV file {path:?} is empty, which is not allowed."
    )]
    EmptyColumnHeader { path: PathBuf, col: usize },

    #[error(
        "Unable to determine the input mode in the aggr CSV file based on the headers: \
        '{}'.\n\nWe support one of the following modes:\n{}",
        headers.join(","),
        supported_modes.iter()
            .map(|(name, specify, headers)| {
                format!(
                    "- If you are looking to aggr outputs from {} by specifying the {}, \
                    the following headers are required: '{}'",
                    name,
                    specify,
                    headers.join(",")
                )
            })
            .join("\n")
    )]
    ErrorDetectingPipeCase {
        supported_modes: Vec<(&'static str, &'static str, Vec<String>)>,
        headers: Vec<String>,
    },

    #[error(
        "Unable to parse sample '{sample_id}' in aggr CSV file. The sample has a VDJ library \
         but the csv is missing required donor or origin headers."
    )]
    ErrorParsingMultiVdjSample { sample_id: String },

    #[error(
        "Origin {origin} is shared between donors {donor1} and {donor2}. Origin needs to be unique identifier \
        for a set of cells coming from a donor and it cannot be shared across donors."
    )]
    DuplicateOriginAcrossDonors {
        origin: String,
        donor1: String,
        donor2: String,
    },
}
