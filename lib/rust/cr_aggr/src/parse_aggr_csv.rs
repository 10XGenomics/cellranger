//! Martian stage PARSE_AGGR_CSV

use crate::errors::{FieldResolutionErrors, ParseAggrCsvErrors};
use anyhow::{bail, Result};
use cr_types::csv_parser::CsvParser;
use itertools::Itertools;
use martian::prelude::*;
use martian_derive::{make_mro, MartianStruct, MartianType};
use martian_filetypes::tabular_file::CsvFile;
use path_clean::PathClean;
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};

/// Different flavours of fields allowed in the CSV file
#[derive(PartialEq, Clone, Copy)]
enum FieldKind {
    FilePath,
    MultiFolderPath,
    Text,
}

#[derive(PartialEq, Clone)]
enum FieldConstraint {
    RequiredUnique,
    Required,
}

struct CsvField {
    name: String,
    kind: FieldKind,
    constraint: FieldConstraint,
}

impl CsvField {
    fn new(name: impl ToString, kind: FieldKind, constraint: FieldConstraint) -> Self {
        CsvField {
            name: name.to_string(),
            kind,
            constraint,
        }
    }
}

struct MultiFilePaths {
    count_v5: PathBuf,
    count_v6: PathBuf,
    vdjt: PathBuf,
    vdjtgd: PathBuf,
    vdjb: PathBuf,
}

impl MultiFilePaths {
    fn new(multi_path: &Path) -> Self {
        MultiFilePaths {
            count_v5: multi_path.join("count/molecule_info.h5"),
            count_v6: multi_path.join("count/sample_molecule_info.h5"),
            vdjt: multi_path.join("vdj_t/vdj_contig_info.pb"),
            vdjtgd: multi_path.join("vdj_t_gd/vdj_contig_info.pb"),
            vdjb: multi_path.join("vdj_b/vdj_contig_info.pb"),
        }
    }
    fn sub_libraries(&self) -> Vec<MultiSubLibraries> {
        let mut ret_vec: Vec<MultiSubLibraries> = vec![];
        if self.count_v5.exists() || self.count_v6.exists() {
            ret_vec.push(MultiSubLibraries::Count);
        }
        if self.vdjt.exists() {
            ret_vec.push(MultiSubLibraries::VdjT);
        }
        if self.vdjtgd.exists() {
            ret_vec.push(MultiSubLibraries::VdjTGD);
        }
        if self.vdjb.exists() {
            ret_vec.push(MultiSubLibraries::VdjB);
        }
        ret_vec
    }
    fn mol_info(&self) -> Option<&Path> {
        if self.count_v5.exists() {
            Some(self.count_v5.as_ref())
        } else if self.count_v6.exists() {
            Some(self.count_v6.as_ref())
        } else {
            None
        }
    }
}

fn process_multi_libraries(
    multi_libs: Vec<MultiAggrLibrary>,
    out_csv: CsvFile<()>,
) -> Result<ParseAggrCsvStageOutputs> {
    let mut count_libraries: Vec<CountLibrary> = vec![];
    let mut vdj_aggr_inputs: Vec<VdjAggrInput> = vec![];
    let mut disable_count_aggr = true;
    let mut disable_vdj_aggr = true;

    // Ensure that all libraries contain the same MultiSubLibraries
    let any_sub_lib: Vec<MultiSubLibraries> = multi_libs
        .iter()
        .flat_map(|mlib| &mlib.libraries)
        .unique()
        .cloned()
        .collect();

    if multi_libs
        .iter()
        .any(|mlib| !any_sub_lib.iter().all(|lib| mlib.libraries.contains(lib)))
    {
        let desc = multi_libs
            .iter()
            .map(|mlib| {
                format!(
                    "- '{}' contains [{}]",
                    mlib.sample_id,
                    mlib.libraries
                        .iter()
                        .map(|lib| match lib {
                            MultiSubLibraries::Count => "count",
                            MultiSubLibraries::VdjT => "vdj_t",
                            MultiSubLibraries::VdjTGD => "vdj_t_gd",
                            MultiSubLibraries::VdjB => "vdj_b",
                        })
                        .join(", "),
                )
            })
            .join("\n");
        bail!(
            "The multi outs folders supplied as inputs to aggr contain an inconsistent set of \
            libraries:\n{}\n\nIn order to aggr multi outs, we require that all the individual \
            inputs contains the same set of library types. If you would like to aggr just one \
            library type, you can supply the `molecule_info` (for count) or the `vdj_contig_info` \
            (for vdj)",
            desc
        );
    }

    // Gather TCR libraries.

    let mut vdj_libs: Vec<VdjAggrCsvLibrary> = vec![];
    for mlib in &multi_libs {
        // For each library:
        // Detect what MultiSubLibraries are present (vdj_t)
        let sublibs = &mlib.libraries;
        if sublibs.contains(&MultiSubLibraries::VdjT) {
            disable_vdj_aggr = false;
            // Check if meta contains origin and donor
            let mut vdj_meta = mlib.meta.clone();
            let donor = vdj_meta.remove(DONOR_HEADER);
            let origin = vdj_meta.remove(ORIGIN_HEADER);
            if donor.is_none() || origin.is_none() {
                return Err(ParseAggrCsvErrors::ErrorParsingMultiVdjSample {
                    sample_id: mlib.sample_id.clone(),
                }
                .into());
            }
            if sublibs.contains(&MultiSubLibraries::VdjT) {
                vdj_libs.push(VdjAggrCsvLibrary {
                    library_id: mlib.sample_id.clone(),
                    vdj_contig_info: MultiFilePaths::new(&mlib.multi_outs).vdjt,
                    donor: donor.clone().unwrap(),
                    origin: origin.clone().unwrap(),
                    meta: vdj_meta.clone(),
                });
            }
        }
    }
    if !vdj_libs.is_empty() {
        vdj_aggr_inputs.push(VdjAggrInput {
            libraries: vdj_libs,
        });
    }

    // Gather BCR libraries.

    let mut vdj_libs: Vec<VdjAggrCsvLibrary> = vec![];
    for mlib in &multi_libs {
        // For each library:
        // Detect what MultiSubLibraries are present (vdj_b)
        let sublibs = &mlib.libraries;
        if sublibs.contains(&MultiSubLibraries::VdjB) {
            disable_vdj_aggr = false;
            // Check if meta contains origin and donor
            let mut vdj_meta = mlib.meta.clone();
            let donor = vdj_meta.remove(DONOR_HEADER);
            let origin = vdj_meta.remove(ORIGIN_HEADER);
            if donor.is_none() || origin.is_none() {
                return Err(ParseAggrCsvErrors::ErrorParsingMultiVdjSample {
                    sample_id: mlib.sample_id.clone(),
                }
                .into());
            }
            if sublibs.contains(&MultiSubLibraries::VdjB) {
                vdj_libs.push(VdjAggrCsvLibrary {
                    library_id: mlib.sample_id.clone(),
                    vdj_contig_info: MultiFilePaths::new(&mlib.multi_outs).vdjb,
                    donor: donor.clone().unwrap(),
                    origin: origin.clone().unwrap(),
                    meta: vdj_meta.clone(),
                });
            }
        }
    }
    if !vdj_libs.is_empty() {
        vdj_aggr_inputs.push(VdjAggrInput {
            libraries: vdj_libs,
        });
    }

    // Gather count libraries.

    for mut mlib in multi_libs {
        // For each library:
        // Detect what MultiSubLibraries are present (count)
        let sublibs = &mlib.libraries;
        if sublibs.contains(&MultiSubLibraries::Count) {
            disable_count_aggr = false;
            count_libraries.push(CountLibrary {
                library_id: mlib.sample_id.clone(),
                molecule_h5: MultiFilePaths::new(&mlib.multi_outs)
                    .mol_info()
                    .unwrap()
                    .to_owned(),
                batch: mlib.meta.remove(BATCH_HEADER),
                meta: mlib.meta.clone(),
            });
        }
    }

    // Return.

    Ok(ParseAggrCsvStageOutputs {
        aggregation_csv: Some(out_csv),
        count_libraries,
        vdj_aggr_inputs,
        disable_count_aggr,
        disable_vdj_aggr,
    })
}

fn abs_path(p: PathBuf) -> std::io::Result<PathBuf> {
    Ok(if p.is_absolute() {
        p
    } else {
        std::env::current_dir()?.join(p)
    }
    .clean())
}

/// Resolve `to_be_resolved` as an absolute path or a path relative to the pipestance_root or
/// a path relative to the folder which contains `agg_csv_path`
fn resolve_path(
    to_be_resolved: &Path,
    pipestance_root: &Path,
    agg_csv_path: &Path,
) -> Result<PathBuf> {
    if to_be_resolved.is_absolute() {
        // Check if path exists
        if to_be_resolved.exists() {
            Ok(to_be_resolved.into())
        } else {
            Err(FieldResolutionErrors::AbsPathDoesntExist {
                path: to_be_resolved.into(),
            }
            .into())
        }
    } else {
        // Check if is relative to pipestance root
        let rel_root_path = pipestance_root.join(to_be_resolved);
        if rel_root_path.exists() {
            return Ok(abs_path(rel_root_path)?);
        }

        // Check if is relative to csv file
        let rel_csv_path = agg_csv_path.parent().unwrap().join(to_be_resolved);
        if rel_csv_path.exists() {
            return Ok(abs_path(rel_csv_path)?);
        }
        Err(FieldResolutionErrors::RelPathDoesntExist {
            path: to_be_resolved.into(),
            agg_csv_path: agg_csv_path.into(),
            pipestance_root: pipestance_root.into(),
        }
        .into())
    }
}

impl FieldKind {
    fn validate_and_resolve(
        self,
        to_be_resolved: &str,
        pipestance_root: &Path,
        agg_csv_path: &Path,
    ) -> Result<String> {
        match self {
            FieldKind::Text => validate_text_value(to_be_resolved),
            FieldKind::FilePath => {
                let to_be_resolved = Path::new(to_be_resolved);
                let resolved = resolve_path(to_be_resolved, pipestance_root, agg_csv_path)?;
                if resolved.is_file() {
                    Ok(resolved.to_str().unwrap().into())
                } else {
                    Err(FieldResolutionErrors::ExpectedFileGotFolder { path: resolved }.into())
                }
            }
            FieldKind::MultiFolderPath => {
                let to_be_resolved = Path::new(to_be_resolved);
                let resolved = resolve_path(to_be_resolved, pipestance_root, agg_csv_path)?;
                let multi_paths = MultiFilePaths::new(&resolved);

                if resolved.is_dir() {
                    if multi_paths.count_v5.is_file()
                        || multi_paths.count_v6.is_file()
                        || multi_paths.vdjt.is_file()
                        || multi_paths.vdjtgd.is_file()
                        || multi_paths.vdjb.is_file()
                    {
                        Ok(resolved.to_str().unwrap().into())
                    } else if resolved.join("per_sample_outs").exists() {
                        Err(FieldResolutionErrors::TopLevelMultiOut.into())
                    } else {
                        Err(FieldResolutionErrors::CouldNotResolveMultiPath {}.into())
                    }
                } else {
                    Err(FieldResolutionErrors::ExpectedFolderGotFile { path: resolved }.into())
                }
            }
        }
    }
}

const ALLOWED_SPECIAL_CHARS: &str = r#"@%^&* ()-_+=[]|:;'."#;
fn validate_text_value(value: &str) -> Result<String> {
    let result = value.to_string();
    for (i, c) in result.chars().enumerate() {
        let valid = c.is_ascii_alphanumeric() || ALLOWED_SPECIAL_CHARS.contains(c);
        if !valid {
            return Err(FieldResolutionErrors::InvalidTextCharacter {
                invalid: c,
                position: i + 1, // 1-based index
                text: result,
                allowed: ALLOWED_SPECIAL_CHARS,
            }
            .into());
        }
    }
    Ok(result)
}

trait AggrCsvParser {
    type SampleDef;
    fn known_csv_fields() -> Vec<CsvField>;
    fn build_sample_def(
        &self,
        valid_row: HashMap<String, String>,
    ) -> Result<Self::SampleDef, ParseAggrCsvErrors>;

    fn check<T: AsRef<str>>(headers: impl IntoIterator<Item = T>) -> bool {
        let headers: Vec<_> = headers.into_iter().map(|x| x.as_ref().into()).collect();
        Self::known_csv_fields()
            .iter()
            .all(|f| headers.contains(&f.name))
    }
    fn known_csv_field_names() -> Vec<String> {
        Self::known_csv_fields()
            .into_iter()
            .map(|f| f.name)
            .collect()
    }
    fn parse_csv(&self, pipestance_root: &Path, csv_file: &Path) -> Result<Vec<Self::SampleDef>> {
        // Open the CSV
        // Read the headers and check them
        // Build all CSV fields, any custom fields become required here
        // Iterate through lines, check for uniqueness, validate field contents
        // Pass a validated row to build_sample_def
        // Collect and return the SampleDef

        let mut csv = {
            let dummy_csv = CsvParser::new(csv_file, std::iter::empty::<String>(), "")?;
            let headers = dummy_csv.headers();
            CsvParser::new(csv_file, headers, "")?
        };
        let num_lines = csv.len();
        let mut ret_vec: Vec<Self::SampleDef> = vec![];
        if num_lines == 0 {
            return Err(ParseAggrCsvErrors::EmptyCsv {
                path: csv.filename().into(),
            }
            .into());
        }
        // TODO Create a new vec for all headers
        let known_headers = Self::known_csv_field_names();
        let all_headers: Vec<_> = csv
            .headers()
            .iter()
            .filter(|&h| !known_headers.contains(h)) // Filter the custom headers
            .map(|h| CsvField::new(h, FieldKind::Text, FieldConstraint::Required)) // Those are required now
            .chain(Self::known_csv_fields()) // Add the known headers
            .collect();

        let mut dup_map = HashMap::<String, HashSet<String>>::new();
        for i in 0..num_lines {
            csv.set_line(i);
            let mut row = HashMap::new();
            for h in &all_headers {
                let Some(val) = csv.try_parse_field::<String>(&h.name, "string")? else {
                    return Err(ParseAggrCsvErrors::EmptyField {
                        path: csv.filename().into(),
                        col: h.name.clone(),
                        line: i + 2, // 1-based index, add one for header as well
                    }
                    .into());
                };
                let resolved_val = h
                    .kind
                    .validate_and_resolve(&val, pipestance_root, csv_file)
                    .map_err(|e| ParseAggrCsvErrors::ErrorParsingField {
                        path: csv.filename().into(),
                        col: h.name.clone(),
                        line: i + 2, // 1-based index, add one for header as well
                        source: e,
                    })?;

                if h.constraint == FieldConstraint::RequiredUnique {
                    if !dup_map.contains_key(&h.name) {
                        dup_map.insert(h.name.clone(), HashSet::<String>::new());
                    }
                    if !dup_map[&h.name].contains(&resolved_val) {
                        dup_map
                            .get_mut(&h.name)
                            .unwrap()
                            .insert(resolved_val.clone());
                    } else {
                        return Err(ParseAggrCsvErrors::DuplicateInUniqueColumn {
                            path: csv.filename().into(),
                            col: h.name.clone(),
                            line: i + 2,
                            duplicate_value: resolved_val,
                        }
                        .into());
                    }
                }
                row.insert(h.name.clone(), resolved_val);
            }
            // row is a valid row now
            ret_vec.push(self.build_sample_def(row)?);
        }
        Ok(ret_vec)
    }
}

const SAMPLE_ID_HEADER: &str = "sample_id";
const MOLECULE_H5_HEADER: &str = "molecule_h5";
const BATCH_HEADER: &str = "batch";
const VDJ_CONTIG_INFO_HEADER: &str = "vdj_contig_info";
pub const DONOR_HEADER: &str = "donor";
pub const ORIGIN_HEADER: &str = "origin";
pub const MULTI_OUTS: &str = "sample_outs";

#[derive(Copy, Clone)]
struct CountParser;
impl AggrCsvParser for CountParser {
    type SampleDef = CountLibrary;
    fn known_csv_fields() -> Vec<CsvField> {
        use FieldConstraint::RequiredUnique;
        use FieldKind::{FilePath, Text};
        vec![
            CsvField::new(SAMPLE_ID_HEADER, Text, RequiredUnique),
            CsvField::new(MOLECULE_H5_HEADER, FilePath, RequiredUnique),
        ]
    }
    fn build_sample_def(
        &self,
        mut valid_row: HashMap<String, String>,
    ) -> Result<Self::SampleDef, ParseAggrCsvErrors> {
        Ok(CountLibrary {
            library_id: valid_row.remove(SAMPLE_ID_HEADER).unwrap(),
            molecule_h5: PathBuf::from(valid_row.remove(MOLECULE_H5_HEADER).unwrap()),
            batch: valid_row.remove(BATCH_HEADER),
            meta: valid_row,
        })
    }
}

struct VdjParser;
impl AggrCsvParser for VdjParser {
    type SampleDef = VdjAggrCsvLibrary;
    fn known_csv_fields() -> Vec<CsvField> {
        use FieldConstraint::{Required, RequiredUnique};
        use FieldKind::{FilePath, Text};
        vec![
            CsvField::new(SAMPLE_ID_HEADER, Text, RequiredUnique),
            CsvField::new(VDJ_CONTIG_INFO_HEADER, FilePath, RequiredUnique),
            CsvField::new(DONOR_HEADER, Text, Required),
            CsvField::new(ORIGIN_HEADER, Text, Required),
        ]
    }
    fn build_sample_def(
        &self,
        mut valid_row: HashMap<String, String>,
    ) -> Result<Self::SampleDef, ParseAggrCsvErrors> {
        Ok(VdjAggrCsvLibrary {
            library_id: valid_row.remove(SAMPLE_ID_HEADER).unwrap(),
            vdj_contig_info: PathBuf::from(valid_row.remove(VDJ_CONTIG_INFO_HEADER).unwrap()),
            donor: valid_row.remove(DONOR_HEADER).unwrap(),
            origin: valid_row.remove(ORIGIN_HEADER).unwrap(),
            meta: valid_row,
        })
    }
}

struct MultiParser;
impl AggrCsvParser for MultiParser {
    type SampleDef = MultiAggrLibrary;
    fn known_csv_fields() -> Vec<CsvField> {
        use FieldConstraint::RequiredUnique;
        use FieldKind::{MultiFolderPath, Text};
        vec![
            CsvField::new(SAMPLE_ID_HEADER, Text, RequiredUnique),
            CsvField::new(MULTI_OUTS, MultiFolderPath, RequiredUnique),
        ]
    }
    fn build_sample_def(
        &self,
        mut valid_row: HashMap<String, String>,
    ) -> Result<Self::SampleDef, ParseAggrCsvErrors> {
        let multi_outs = PathBuf::from(valid_row.remove(MULTI_OUTS).unwrap());
        let libraries = MultiFilePaths::new(&multi_outs).sub_libraries();
        Ok(MultiAggrLibrary {
            sample_id: valid_row.remove(SAMPLE_ID_HEADER).unwrap(),
            multi_outs,
            libraries,
            meta: valid_row,
        })
    }
}

#[derive(Debug, Serialize, Clone, PartialEq, Eq, Hash)]
enum MultiSubLibraries {
    VdjT,
    VdjTGD,
    VdjB,
    Count,
}

#[derive(Debug, Serialize, Clone, PartialEq)]
struct MultiAggrLibrary {
    sample_id: String,
    multi_outs: PathBuf,
    meta: HashMap<String, String>,
    libraries: Vec<MultiSubLibraries>,
}

/// Types of aggr config CSVs
#[derive(Debug, PartialEq)]
enum PipeCase {
    Case1Count,
    Case2Vdj,
    Case3Multi,
}

//TODO address conflicts (when two or more cases match)
fn detect_pipe_case<T: AsRef<str>>(
    headers: impl IntoIterator<Item = T> + Copy,
) -> Result<PipeCase, ParseAggrCsvErrors> {
    let count_match = CountParser::check(headers);
    let vdj_match = VdjParser::check(headers);
    let multi_match = MultiParser::check(headers);
    match (count_match, vdj_match, multi_match) {
        (true, false, false) => Ok(PipeCase::Case1Count),
        (false, true, false) => Ok(PipeCase::Case2Vdj),
        (false, false, true) => Ok(PipeCase::Case3Multi),
        _ => Err(ParseAggrCsvErrors::ErrorDetectingPipeCase {
            supported_modes: vec![
                (
                    "count",
                    "molecule_info.h5",
                    CountParser::known_csv_field_names(),
                ),
                (
                    "vdj",
                    "vdj_contig_info.pb",
                    VdjParser::known_csv_field_names(),
                ),
                (
                    "multi",
                    "output folder",
                    MultiParser::known_csv_field_names(),
                ),
            ],
            headers: headers.into_iter().map(|h| h.as_ref().into()).collect(),
        }),
    }
}

//TODO: count, vdj, multi -> only one of them should be satisfied. otherwise error.CsvParser
// if it's multi, check if it's multi vdj.
// batch is optional and not forbidden for vdj.

#[derive(Debug, Clone, Deserialize, MartianStruct)]
pub struct ParseAggrCsvStageInputs {
    pipestance_root: PathBuf,
    aggregation_csv: Option<CsvFile<()>>,
}

/// Maintains the existing count sample def used in aggr
#[derive(Debug, Serialize, Deserialize, MartianType, Clone, PartialEq)]
struct CountLibrary {
    library_id: String,
    molecule_h5: PathBuf,
    #[serde(skip_serializing_if = "Option::is_none")]
    batch: Option<String>,
    #[serde(flatten)]
    meta: HashMap<String, String>,
}

/// All the libraries defined in the aggr input Csv
#[derive(Debug, Serialize, Deserialize, MartianStruct, Clone, PartialEq, Eq)]
pub struct VdjAggrInput {
    pub libraries: Vec<VdjAggrCsvLibrary>,
}

/// A single library in the aggr csv input. This corresponds to data in a single row
#[derive(Debug, Serialize, Deserialize, MartianStruct, Clone, PartialEq, Eq)]
pub struct VdjAggrCsvLibrary {
    pub library_id: String,
    pub vdj_contig_info: PathBuf,
    pub donor: String,
    pub origin: String,
    pub meta: HashMap<String, String>,
}

fn check_outputs(outs: ParseAggrCsvStageOutputs) -> Result<ParseAggrCsvStageOutputs> {
    // check for duplicate origin between separate donors
    let mut dup_map = HashMap::<String, String>::new();
    for vdj_inp in &outs.vdj_aggr_inputs {
        for library in &vdj_inp.libraries {
            if !dup_map.contains_key(&library.origin) {
                dup_map.insert(library.origin.clone(), library.donor.clone());
            } else if dup_map[&library.origin] != library.donor {
                return Err(ParseAggrCsvErrors::DuplicateOriginAcrossDonors {
                    origin: library.origin.clone(),
                    donor1: dup_map[&library.origin].clone(),
                    donor2: library.donor.clone(),
                }
                .into());
            }
        }
    }
    Ok(outs)
}
#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct ParseAggrCsvStageOutputs {
    /// Input file copied, because it is a pipeline output as well
    aggregation_csv: Option<CsvFile<()>>,
    count_libraries: Vec<CountLibrary>,
    /// In the general case of multi, we could have upto two sets of inputs
    /// one for TCR and the other for IG
    vdj_aggr_inputs: Vec<VdjAggrInput>,
    disable_count_aggr: bool,
    disable_vdj_aggr: bool,
}

// This is our stage struct
pub struct ParseAggrCsv;

fn handle_deprecated_headers(csv: &Path) -> Result<()> {
    let req_headers: &[&str] = &[];
    let parser = CsvParser::new(csv, req_headers, "aggr csv")?;
    let headers = parser.headers();
    let has_lib_id = headers.contains(&"library_id".into());
    let has_lib_outs = headers.contains(&"library_outs".into());
    if has_lib_id && has_lib_outs {
        bail!(
            "The columns 'library_id' and 'library_outs' are deprecated. Please use the column \
            name 'sample_id' instead of 'library_id' and 'sample_outs' instead of 'library_outs'"
        )
    } else if has_lib_id {
        bail!(
            "The column 'library_id' is deprecated. Please use the column name 'sample_id' \
            instead of 'library_id'."
        )
    }
    Ok(())
}

#[make_mro]
impl MartianMain for ParseAggrCsv {
    type StageInputs = ParseAggrCsvStageInputs;
    type StageOutputs = ParseAggrCsvStageOutputs;
    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        // Need to handle parsing in three cases:
        // - Count only aggr
        // - Vdj only aggr
        // - Count + Vdj aggr (multi)

        let Some(aggregation_csv) = args.aggregation_csv.as_ref() else {
            return Ok(ParseAggrCsvStageOutputs {
                aggregation_csv: None,
                count_libraries: Vec::new(),
                vdj_aggr_inputs: Vec::new(),
                disable_count_aggr: true,
                disable_vdj_aggr: true,
            });
        };

        handle_deprecated_headers(aggregation_csv)?;

        // Check if csv file exists and find the headers
        let csv = CsvParser::new(aggregation_csv, [SAMPLE_ID_HEADER], "aggr csv")?;

        // Error out if empty header.
        for (i, h) in csv.headers().iter().enumerate() {
            if h.is_empty() {
                // +1 for 1 based indexing
                return Err(ParseAggrCsvErrors::EmptyColumnHeader {
                    path: csv.filename().into(),
                    col: i + 1,
                }
                .into());
            }
        }

        // Make a copy of the csv
        let out_csv: CsvFile<()> = rover.make_path("aggregation.csv");
        std::fs::copy(aggregation_csv, &out_csv)?;

        // step 1: detect which parser to use
        let pipe = detect_pipe_case(csv.headers())?;
        // step 2: run parser and build the libraries
        let outs = match pipe {
            PipeCase::Case1Count => ParseAggrCsvStageOutputs {
                aggregation_csv: Some(out_csv),
                count_libraries: CountParser.parse_csv(&args.pipestance_root, aggregation_csv)?,
                vdj_aggr_inputs: Vec::new(),
                disable_count_aggr: false,
                disable_vdj_aggr: true,
            },
            PipeCase::Case2Vdj => ParseAggrCsvStageOutputs {
                aggregation_csv: Some(out_csv),
                count_libraries: Vec::new(),
                vdj_aggr_inputs: vec![VdjAggrInput {
                    libraries: VdjParser.parse_csv(&args.pipestance_root, aggregation_csv)?,
                }],
                disable_count_aggr: true,
                disable_vdj_aggr: false,
            },
            PipeCase::Case3Multi => process_multi_libraries(
                MultiParser.parse_csv(&args.pipestance_root, aggregation_csv)?,
                out_csv,
            )?,
        };
        check_outputs(outs)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use insta::{assert_debug_snapshot, assert_json_snapshot, assert_snapshot};
    use martian_filetypes::json_file::JsonFile;
    use martian_filetypes::tabular_file::CsvFile;
    use martian_filetypes::FileTypeRead;
    use pretty_assertions::assert_eq;

    impl ParseAggrCsvStageInputs {
        #[cfg(test)]
        fn test(csv: &str) -> Self {
            ParseAggrCsvStageInputs {
                pipestance_root: PathBuf::new(),
                aggregation_csv: Some(CsvFile::from(csv)),
            }
        }
    }

    fn test_roundtrip_and_representation(input: &str) -> Result<()> {
        let input_file = JsonFile::from(input);
        let sample_defs: Vec<CountLibrary> = input_file.read()?;
        assert_debug_snapshot!(sample_defs);
        let out = serde_json::to_string(&sample_defs)?;
        assert_eq!(
            serde_json::from_str::<serde_json::Value>(&out)?,
            serde_json::from_str::<serde_json::Value>(&std::fs::read_to_string(&input_file)?)?
        );
        Ok(())
    }

    fn err_string(csv: &str) -> String {
        let output = ParseAggrCsv.test_run_tmpdir(ParseAggrCsvStageInputs::test(csv));
        format!("{:#}", output.unwrap_err())
    }

    fn error_to_string(err: anyhow::Error) -> String {
        format!("{err:#}")
    }

    macro_rules! err_snapshot {
        ($fname: expr) => {
            insta::assert_snapshot!(err_string($fname));
        };
    }

    #[test]
    fn test_sample_def_roundtrip() -> Result<()> {
        test_roundtrip_and_representation("test_resources/parse_aggr_csv/sample_defs_1.json")?; // library_id,molecule_info
        test_roundtrip_and_representation("test_resources/parse_aggr_csv/sample_defs_2.json")?; // library_id,molecule_info,batch
        test_roundtrip_and_representation("test_resources/parse_aggr_csv/sample_defs_3.json")?; // library_id,molecule_info,key1
        test_roundtrip_and_representation("test_resources/parse_aggr_csv/sample_defs_4.json")?; // library_id,molecule_info,key1,batch
        Ok(())
    }

    #[test]
    fn test_csv_doesnt_exist() {
        err_snapshot!("non/existant/file.csv");
    }

    #[test]
    fn test_missing_lib_id() {
        // only extra columns
        err_snapshot!(
            "test_resources/parse_aggr_csv/pipes/csvs/missing_h5_id_column_three_user_def_fields.csv"
        );
    }

    #[test]
    fn test_missing_lib_id_vdj() {
        err_snapshot!("test_resources/parse_aggr_csv/pipes/csvs/vdj_proto_no_lib.csv");
    }

    #[test]
    fn test_missing_donor() {
        let csv = "test_resources/parse_aggr_csv/pipes/csvs/vdj_proto_no_donor.csv";
        err_snapshot!(csv);
    }
    #[test]
    fn test_missing_origin() {
        let csv = "test_resources/parse_aggr_csv/pipes/csvs/vdj_proto_no_origin.csv";
        err_snapshot!(csv);
    }

    #[test]
    fn test_duplicate_h5() {
        let csv = "test_resources/parse_aggr_csv/pipes/csvs/duplicate_h5.csv";
        assert_eq!(
            err_string(csv),
            ParseAggrCsvErrors::DuplicateInUniqueColumn {
                path: csv.into(),
                col: MOLECULE_H5_HEADER.to_string(),
                line: 3,
                duplicate_value: std::env::current_dir()
                    .unwrap()
                    .join("test_resources/parse_aggr_csv/pipes/run1_outs/molecule_info.h5")
                    .to_str()
                    .unwrap()
                    .into(),
            }
            .to_string()
        );
    }
    #[test]
    fn test_duplicate_vdj_proto() {
        let csv = "test_resources/parse_aggr_csv/pipes/csvs/vdj_duplicate_proto.csv";
        assert_eq!(
            err_string(csv),
            ParseAggrCsvErrors::DuplicateInUniqueColumn {
                path: csv.into(),
                col: VDJ_CONTIG_INFO_HEADER.to_string(),
                line: 3,
                duplicate_value: std::env::current_dir()
                    .unwrap()
                    .join("test_resources/parse_aggr_csv/pipes/run3_outs/vdj_contig_info.pb")
                    .to_str()
                    .unwrap()
                    .into(),
            }
            .to_string()
        );
    }
    #[test]
    fn test_duplicate_sample_id() {
        let csv = "test_resources/parse_aggr_csv/pipes/csvs/duplicate_sample_id.csv";
        assert_eq!(
            err_string(csv),
            ParseAggrCsvErrors::DuplicateInUniqueColumn {
                path: csv.into(),
                col: SAMPLE_ID_HEADER.to_string(),
                line: 3,
                duplicate_value: "sample1".to_string()
            }
            .to_string()
        );
    }

    #[test]
    fn test_duplicate_sample_id_vdj() {
        let csv = "test_resources/parse_aggr_csv/pipes/csvs/vdj_duplicate_sample_id.csv";
        assert_eq!(
            err_string(csv),
            ParseAggrCsvErrors::DuplicateInUniqueColumn {
                path: csv.into(),
                col: SAMPLE_ID_HEADER.to_string(),
                line: 3,
                duplicate_value: "sample1".to_string()
            }
            .to_string()
        );
    }

    #[test]
    fn test_duplicate_line() {
        let csv = "test_resources/parse_aggr_csv/pipes/csvs/duplicate_line.csv";
        assert_eq!(
            err_string(csv),
            ParseAggrCsvErrors::DuplicateInUniqueColumn {
                path: csv.into(),
                col: SAMPLE_ID_HEADER.to_string(),
                line: 3,
                duplicate_value: "sample1".to_string()
            }
            .to_string()
        );
    }

    #[test]
    fn test_empty_lib_id() {
        let csv = "test_resources/parse_aggr_csv/pipes/csvs/empty_lib_id.csv";
        assert_eq!(
            err_string(csv),
            ParseAggrCsvErrors::EmptyField {
                path: csv.into(),
                col: SAMPLE_ID_HEADER.to_string(),
                line: 2
            }
            .to_string()
        );
    }

    #[test]
    fn test_empty_vdj_proto() {
        let csv = "test_resources/parse_aggr_csv/pipes/csvs/vdj_empty_proto.csv";
        assert_eq!(
            err_string(csv),
            ParseAggrCsvErrors::EmptyField {
                path: csv.into(),
                col: VDJ_CONTIG_INFO_HEADER.to_string(),
                line: 2
            }
            .to_string()
        );
    }

    #[test]
    fn test_empty_donor() {
        let csv = "test_resources/parse_aggr_csv/pipes/csvs/vdj_empty_donor.csv";
        assert_eq!(
            err_string(csv),
            ParseAggrCsvErrors::EmptyField {
                path: csv.into(),
                col: DONOR_HEADER.to_string(),
                line: 3
            }
            .to_string()
        );
    }

    #[test]
    fn test_missing_column_title() {
        let csv = "test_resources/parse_aggr_csv/bad_csvs/missing_header.csv";
        assert_eq!(
            err_string(csv),
            ParseAggrCsvErrors::EmptyColumnHeader {
                path: csv.into(),
                col: 3,
            }
            .to_string()
        );
    }

    #[test]
    fn test_empty_extra() {
        let csv = "test_resources/parse_aggr_csv/pipes/csvs/vdj_empty_extra_field.csv";
        assert_eq!(
            err_string(csv),
            ParseAggrCsvErrors::EmptyField {
                path: csv.into(),
                col: "key".to_string(),
                line: 3
            }
            .to_string()
        );
    }

    #[test]
    fn test_invalid_h5() {
        let csv = "test_resources/parse_aggr_csv/pipes/csvs/missing_h5_file.csv";
        assert_eq!(
            err_string(csv),
            error_to_string(anyhow::Error::from(ParseAggrCsvErrors::ErrorParsingField {
                col: MOLECULE_H5_HEADER.to_string(),
                line: 2,
                path: PathBuf::from(csv),
                source: anyhow::Error::from(FieldResolutionErrors::AbsPathDoesntExist {
                    path: PathBuf::from("/no/such/file/molecule_info.h5"),
                })
            }))
        );
    }

    #[test]
    fn test_invalid_vdj_proto() {
        let csv = "test_resources/parse_aggr_csv/pipes/csvs/vdj_invalid_proto.csv";
        assert_eq!(
            err_string(csv),
            error_to_string(anyhow::Error::from(ParseAggrCsvErrors::ErrorParsingField {
                col: VDJ_CONTIG_INFO_HEADER.to_string(),
                line: 2,
                path: PathBuf::from(csv),
                source: anyhow::Error::from(FieldResolutionErrors::RelPathDoesntExist {
                    path: PathBuf::from("../invalid/dir/vdj_contig_info.pb"),
                    agg_csv_path: PathBuf::from(csv),
                    pipestance_root: PathBuf::new()
                })
            }))
        );
    }

    #[test]
    fn test_invalid_vdj_dup_origin() {
        let csv = "test_resources/parse_aggr_csv/pipes/csvs/invalid_vdj_1.csv";
        assert_eq!(
            err_string(csv),
            "Origin o1 is shared between donors d1 and d2. Origin needs to be unique \
            identifier for a set of cells coming from a donor and it cannot be shared across donors."
                .to_string()
        );
    }

    #[test]
    fn test_only_lib_id() {
        let csv = "test_resources/parse_aggr_csv/pipes/csvs/only_lib_id_and_extra.csv";
        err_snapshot!(csv);
    }
    #[test]
    fn test_missing_header() {
        err_snapshot!("test_resources/parse_aggr_csv/pipes/csvs/no_header.csv");
    }

    #[test]
    fn test_valid_count_csv_1() -> Result<()> {
        let inp = ParseAggrCsvStageInputs::test(
            "test_resources/parse_aggr_csv/pipes/csvs/valid_count_1.csv",
        );
        let out = ParseAggrCsv.test_run_tmpdir(inp)?;
        let expected_count_sample_def = &[
            CountLibrary {
                library_id: "sample1".to_string(),
                molecule_h5: abs_path(PathBuf::from(
                    "test_resources/parse_aggr_csv/pipes/run1_outs/molecule_info.h5",
                ))?,
                batch: None,
                meta: HashMap::new(),
            },
            CountLibrary {
                library_id: "sample2".to_string(),
                molecule_h5: abs_path(PathBuf::from(
                    "test_resources/parse_aggr_csv/pipes/run2_outs/molecule_info.h5",
                ))?,
                batch: None,
                meta: HashMap::new(),
            },
        ];

        assert_eq!(&out.count_libraries, expected_count_sample_def);
        assert_eq!(&out.vdj_aggr_inputs, &[]);
        assert!(!out.disable_count_aggr);
        assert!(out.disable_vdj_aggr);
        Ok(())
    }

    macro_rules! hashmap {
        ($( $key: expr => $val: expr ),*) => {{
             let mut map = ::std::collections::HashMap::new();
             $( map.insert($key.into(), $val.into()); )*
             map
        }}
    }

    #[test]
    fn test_valid_count_csv_2() -> Result<()> {
        let inp = ParseAggrCsvStageInputs::test(
            "test_resources/parse_aggr_csv/pipes/csvs/valid_count_2.csv",
        );
        let out = ParseAggrCsv.test_run_tmpdir(inp)?;
        let expected_count_sample_def = vec![
            CountLibrary {
                library_id: "sample1".to_string(),
                molecule_h5: abs_path(PathBuf::from(
                    "test_resources/parse_aggr_csv/pipes/run1_outs/molecule_info.h5",
                ))?,
                batch: None,
                meta: hashmap!["key2" => "val", "key1" => "foo"],
            },
            CountLibrary {
                library_id: "sample2".to_string(),
                molecule_h5: abs_path(PathBuf::from(
                    "test_resources/parse_aggr_csv/pipes/run2_outs/molecule_info.h5",
                ))?,
                batch: None,
                meta: hashmap!["key2" => "val", "key1" => "bar"],
            },
        ];
        assert_eq!(out.count_libraries, expected_count_sample_def);
        assert!(out.vdj_aggr_inputs.is_empty());
        assert!(!out.disable_count_aggr);
        assert!(out.disable_vdj_aggr);
        Ok(())
    }

    #[test]
    fn test_valid_count_csv_3() -> Result<()> {
        let inp = ParseAggrCsvStageInputs::test(
            "test_resources/parse_aggr_csv/pipes/csvs/valid_count_3.csv",
        );
        let out = ParseAggrCsv.test_run_tmpdir(inp)?;
        let expected_count_sample_def = vec![
            CountLibrary {
                library_id: "sample1".to_string(),
                molecule_h5: abs_path(PathBuf::from(
                    "test_resources/parse_aggr_csv/pipes/run1_outs/molecule_info.h5",
                ))?,
                batch: Some("b1".to_string()),
                meta: hashmap!["key1" => "val"],
            },
            CountLibrary {
                library_id: "sample2".to_string(),
                molecule_h5: abs_path(PathBuf::from(
                    "test_resources/parse_aggr_csv/pipes/run2_outs/molecule_info.h5",
                ))?,
                batch: Some("b2".to_string()),
                meta: hashmap!["key1" => "val"],
            },
        ];

        assert_eq!(out.count_libraries, expected_count_sample_def);
        assert!(out.vdj_aggr_inputs.is_empty());
        assert!(!out.disable_count_aggr);
        assert!(out.disable_vdj_aggr);
        Ok(())
    }

    #[test]
    fn test_valid_vdj_csv_1() -> Result<()> {
        let inp = ParseAggrCsvStageInputs::test(
            "test_resources/parse_aggr_csv/pipes/csvs/valid_vdj_1.csv",
        );
        let out = ParseAggrCsv.test_run_tmpdir(inp)?;
        let expected_vdj_sample_def = vec![
            VdjAggrCsvLibrary {
                library_id: "sample1".to_string(),
                vdj_contig_info: abs_path(PathBuf::from(
                    "test_resources/parse_aggr_csv/pipes/run3_outs/vdj_contig_info.pb",
                ))?,
                donor: "d1".to_string(),
                origin: "o1".to_string(),
                meta: HashMap::new(),
            },
            VdjAggrCsvLibrary {
                library_id: "sample2".to_string(),
                vdj_contig_info: abs_path(PathBuf::from(
                    "test_resources/parse_aggr_csv/pipes/run4_outs/vdj_contig_info.pb",
                ))?,
                donor: "d1".to_string(),
                origin: "o2".to_string(),
                meta: HashMap::new(),
            },
        ];

        assert!(out.count_libraries.is_empty());
        assert_eq!(
            out.vdj_aggr_inputs,
            vec![VdjAggrInput {
                libraries: expected_vdj_sample_def
            }]
        );
        assert!(out.disable_count_aggr);
        assert!(!out.disable_vdj_aggr);

        Ok(())
    }
    #[test]
    fn test_valid_vdj_csv_2() -> Result<()> {
        let inp = ParseAggrCsvStageInputs::test(
            "test_resources/parse_aggr_csv/pipes/csvs/valid_vdj_2.csv",
        );
        let out = ParseAggrCsv.test_run_tmpdir(inp)?;
        let expected_vdj_sample_def = vec![
            VdjAggrCsvLibrary {
                library_id: "sample1".to_string(),
                vdj_contig_info: abs_path(PathBuf::from(
                    "test_resources/parse_aggr_csv/pipes/run3_outs/vdj_contig_info.pb",
                ))?,
                donor: "d1".to_string(),
                origin: "o1".to_string(),
                meta: hashmap!["key2" => "val", "key1" => "bar1"],
            },
            VdjAggrCsvLibrary {
                library_id: "sample2".to_string(),
                vdj_contig_info: abs_path(PathBuf::from(
                    "test_resources/parse_aggr_csv/pipes/run4_outs/vdj_contig_info.pb",
                ))?,
                donor: "d1".to_string(),
                origin: "o1".to_string(),
                meta: hashmap!["key2" => "val", "key1" => "bar2"],
            },
        ];

        assert!(out.count_libraries.is_empty());
        assert_eq!(
            out.vdj_aggr_inputs,
            vec![VdjAggrInput {
                libraries: expected_vdj_sample_def
            }]
        );
        assert!(out.disable_count_aggr);
        assert!(!out.disable_vdj_aggr);
        Ok(())
    }

    #[test]
    fn test_valid_multi_1() -> Result<()> {
        let inp = ParseAggrCsvStageInputs::test(
            "test_resources/parse_aggr_csv/pipes/csvs/valid_multi_1.csv",
        );
        let out = ParseAggrCsv.test_run_tmpdir(inp)?;
        let expected_count_sample_def = vec![
            CountLibrary {
                library_id: "sample1".to_string(),
                molecule_h5: abs_path(PathBuf::from(
                    "test_resources/parse_aggr_csv/pipes/run7_outs/count/molecule_info.h5",
                ))?,
                batch: Some("b1".into()),
                meta: hashmap![
                    "key2" => "val1",
                    "key1" => "bar1",
                    "origin" => "o1",
                    "donor" => "d1"
                ],
            },
            CountLibrary {
                library_id: "sample2".to_string(),
                molecule_h5: abs_path(PathBuf::from(
                    "test_resources/parse_aggr_csv/pipes/run9_outs/count/molecule_info.h5",
                ))?,
                batch: Some("b1".into()),
                meta: hashmap![
                    "key2" => "val2",
                    "key1" => "bar2",
                    "origin" => "o2",
                    "donor" => "d2"
                ],
            },
        ];
        let expected_vdj_t_sample_def = vec![
            VdjAggrCsvLibrary {
                library_id: "sample1".to_string(),
                vdj_contig_info: abs_path(PathBuf::from(
                    "test_resources/parse_aggr_csv/pipes/run7_outs/vdj_t/vdj_contig_info.pb",
                ))?,
                donor: "d1".to_string(),
                origin: "o1".to_string(),
                meta: hashmap!["key2" => "val1", "key1" => "bar1", "batch" => "b1"],
            },
            VdjAggrCsvLibrary {
                library_id: "sample2".to_string(),
                vdj_contig_info: abs_path(PathBuf::from(
                    "test_resources/parse_aggr_csv/pipes/run9_outs/vdj_t/vdj_contig_info.pb",
                ))?,
                donor: "d2".to_string(),
                origin: "o2".to_string(),
                meta: hashmap!["key2" => "val2", "key1" => "bar2", "batch" => "b1"],
            },
        ];

        let expected_vdj_b_sample_def = vec![
            VdjAggrCsvLibrary {
                library_id: "sample1".to_string(),
                vdj_contig_info: abs_path(PathBuf::from(
                    "test_resources/parse_aggr_csv/pipes/run7_outs/vdj_b/vdj_contig_info.pb",
                ))?,
                donor: "d1".to_string(),
                origin: "o1".to_string(),
                meta: hashmap!["key2" => "val1", "key1" => "bar1", "batch" => "b1"],
            },
            VdjAggrCsvLibrary {
                library_id: "sample2".to_string(),
                vdj_contig_info: abs_path(PathBuf::from(
                    "test_resources/parse_aggr_csv/pipes/run9_outs/vdj_b/vdj_contig_info.pb",
                ))?,
                donor: "d2".to_string(),
                origin: "o2".to_string(),
                meta: hashmap!["key2" => "val2", "key1" => "bar2", "batch" => "b1"],
            },
        ];

        assert_eq!(out.count_libraries, expected_count_sample_def);
        assert_eq!(
            out.vdj_aggr_inputs,
            vec![
                VdjAggrInput {
                    libraries: expected_vdj_t_sample_def
                },
                VdjAggrInput {
                    libraries: expected_vdj_b_sample_def
                },
            ]
        );
        assert!(!out.disable_count_aggr);
        assert!(!out.disable_vdj_aggr);
        Ok(())
    }

    #[test]
    fn test_valid_multi_3() -> Result<()> {
        let inp = ParseAggrCsvStageInputs::test(
            "test_resources/parse_aggr_csv/pipes/csvs/valid_multi_3.csv",
        );
        let out = ParseAggrCsv.test_run_tmpdir(inp)?;
        let expected_count_sample_def = vec![CountLibrary {
            library_id: "sample1".to_string(),
            molecule_h5: abs_path(PathBuf::from(
                "test_resources/parse_aggr_csv/pipes/run5_outs/count/molecule_info.h5",
            ))?,
            batch: Some("b1".to_string()),
            meta: hashmap![
            "key1" => "bar1",
            "key2" => "val1"],
        }];

        assert_eq!(out.count_libraries, expected_count_sample_def);
        assert!(out.vdj_aggr_inputs.is_empty());
        assert!(!out.disable_count_aggr);
        assert!(out.disable_vdj_aggr);
        Ok(())
    }

    #[test]
    fn test_invalid_multi() -> Result<()> {
        // Path is not a multi outs
        err_snapshot!("test_resources/parse_aggr_csv/pipes/csvs/invalid_multi_1.csv");

        // Path is a multi outs, but non of expected files (count/molecule_info.h5,
        // vdj_b/vdj_contig_info.pb, vdj_t/vdj_contig_info.pb) are available
        err_snapshot!("test_resources/parse_aggr_csv/pipes/csvs/invalid_multi_2.csv");

        // Path is a multi outs which includes a vdj library, but donor and origin
        // are not included in the csv
        err_snapshot!("test_resources/parse_aggr_csv/pipes/csvs/invalid_multi_3.csv");

        // Origin is shared across two donors
        err_snapshot!("test_resources/parse_aggr_csv/pipes/csvs/invalid_multi_4.csv");
        Ok(())
    }

    #[test]
    fn test_invalid_multi_5() -> Result<()> {
        // Libraries are not consistent between different output folders
        err_snapshot!("test_resources/parse_aggr_csv/pipes/csvs/invalid_multi_5.csv");
        Ok(())
    }

    #[test]
    fn test_check_headers_count() -> Result<()> {
        assert!(CountParser::check([SAMPLE_ID_HEADER, MOLECULE_H5_HEADER]),);
        assert!(CountParser::check([
            SAMPLE_ID_HEADER,
            MOLECULE_H5_HEADER,
            BATCH_HEADER
        ]));
        assert!(!CountParser::check([SAMPLE_ID_HEADER]));
        assert!(!CountParser::check([
            SAMPLE_ID_HEADER,
            VDJ_CONTIG_INFO_HEADER
        ]));
        // Matches two cases, must be resolved in detect_pipe
        assert!(CountParser::check([
            SAMPLE_ID_HEADER,
            MOLECULE_H5_HEADER,
            VDJ_CONTIG_INFO_HEADER
        ]));

        Ok(())
    }

    #[test]
    fn test_detect_pipe() -> Result<()> {
        assert_eq!(
            detect_pipe_case([SAMPLE_ID_HEADER, MOLECULE_H5_HEADER])?,
            PipeCase::Case1Count
        );
        assert_eq!(
            detect_pipe_case([
                SAMPLE_ID_HEADER,
                VDJ_CONTIG_INFO_HEADER,
                DONOR_HEADER,
                ORIGIN_HEADER
            ])?,
            PipeCase::Case2Vdj
        );
        assert_eq!(
            detect_pipe_case([SAMPLE_ID_HEADER, MULTI_OUTS])?,
            PipeCase::Case3Multi
        );
        assert_eq!(
            detect_pipe_case([SAMPLE_ID_HEADER, MULTI_OUTS, DONOR_HEADER, ORIGIN_HEADER])?,
            PipeCase::Case3Multi
        );

        Ok(())
    }

    #[test]
    fn test_resolve_file_path() -> Result<()> {
        let cargo_dir = std::env::current_dir()?;
        let abs_path =
            cargo_dir.join("test_resources/parse_aggr_csv/path_resolution/pipe_root/molecule.h5");

        let pipestance_root =
            cargo_dir.join("test_resources/parse_aggr_csv/path_resolution/pipe_root/");
        // test 1:
        // test_resources/parse_aggr_csv/path_resolution/
        //      agg.csv
        //      pipe_root/
        //          molecule.h5
        let agg_csv_path = cargo_dir.join("test_resources/parse_aggr_csv/path_resolution/agg.csv");

        // absolute path
        assert_eq!(
            abs_path.to_str().unwrap(),
            FieldKind::FilePath.validate_and_resolve(
                abs_path.to_str().unwrap(),
                &pipestance_root,
                &agg_csv_path,
            )?
        );

        // absolute path
        assert_eq!(
            abs_path,
            resolve_path(&abs_path, &pipestance_root, &agg_csv_path)?
        );

        // File instead of folder
        assert!(FieldKind::MultiFolderPath
            .validate_and_resolve(abs_path.to_str().unwrap(), &pipestance_root, &agg_csv_path,)
            .is_err());

        assert_eq!(
            abs_path.to_str().unwrap(),
            FieldKind::FilePath.validate_and_resolve(
                "molecule.h5", // relative to pipestance root
                &pipestance_root,
                &agg_csv_path,
            )?
        );

        assert_eq!(
            abs_path,
            // relative to pipestance root
            resolve_path(Path::new("molecule.h5"), &pipestance_root, &agg_csv_path,)?
        );

        assert_eq!(
            abs_path.to_str().unwrap(),
            FieldKind::FilePath
                .validate_and_resolve(
                    "pipe_root/molecule.h5", // relative to agg.csv dir
                    &pipestance_root,
                    &agg_csv_path,
                )
                .unwrap()
        );

        assert_eq!(
            abs_path,
            // relative to agg.csv dir
            resolve_path(
                Path::new("pipe_root/molecule.h5"),
                &pipestance_root,
                &agg_csv_path,
            )?
        );

        // test 2:
        // test_resources/parse_aggr_csv/path_resolution/
        //      csv_dir/
        //          agg.csv
        //      pipe_root/
        //          molecule.h5
        assert_eq!(
            abs_path.to_str().unwrap(),
            FieldKind::FilePath.validate_and_resolve(
                "../pipe_root/molecule.h5", // relative to agg.csv dir
                &pipestance_root,
                &cargo_dir.join("test_resources/parse_aggr_csv/path_resolution/csv_dir/agg.csv"),
            )?
        );
        assert_eq!(
            abs_path,
            // relative to agg.csv dir
            resolve_path(
                Path::new("../pipe_root/molecule.h5"),
                &pipestance_root,
                &cargo_dir.join("test_resources/parse_aggr_csv/path_resolution/csv_dir/agg.csv"),
            )?
        );

        // test 3: Incorrect path
        // absolute path
        assert!(resolve_path(Path::new("/no/such/path"), &pipestance_root, &agg_csv_path).is_err());

        let wrong_abs_path =
            cargo_dir.join("test_resources/parse_aggr_csv/path_resolution/pipe_root/WRONG.h5");
        assert!(FieldKind::FilePath
            .validate_and_resolve(
                wrong_abs_path.to_str().unwrap(),
                &pipestance_root,
                &agg_csv_path,
            )
            .is_err());

        // relative path
        assert!(resolve_path(Path::new("WRONG.h5"), &pipestance_root, &agg_csv_path).is_err());
        assert!(FieldKind::FilePath
            .validate_and_resolve("WRONG.h5", &pipestance_root, &agg_csv_path,)
            .is_err());

        assert!(FieldKind::FilePath
            .validate_and_resolve("../WRONG/molecule.h5", &pipestance_root, &agg_csv_path,)
            .is_err());

        Ok(())
    }

    #[test]
    fn test_resolve_multi_folder_path() -> Result<()> {
        let cargo_dir = std::env::current_dir()?;
        // resolve a directory
        // test_resources/parse_aggr_csv/path_resolution/
        //      csv_dir/
        //          agg.csv
        //      pipe_root/
        //          another_dir/
        //              (not a multi dir!)
        let abs_path =
            cargo_dir.join("test_resources/parse_aggr_csv/path_resolution/pipe_root/another_dir");
        let pipestance_root =
            cargo_dir.join("test_resources/parse_aggr_csv/path_resolution/pipe_root");
        let agg_csv_path = cargo_dir.join("test_resources/parse_aggr_csv/path_resolution/agg.csv");

        // absolute path (can be resolved)
        assert_eq!(
            abs_path,
            resolve_path(&abs_path, &pipestance_root, &agg_csv_path)?
        );
        // absolute path is not a correct multi out
        let abs_path_str = abs_path.to_str().unwrap();
        assert!(FieldKind::MultiFolderPath
            .validate_and_resolve(abs_path_str, &pipestance_root, &agg_csv_path,)
            .is_err());

        // Expected a file, for a folder
        assert!(FieldKind::FilePath
            .validate_and_resolve(abs_path_str, &pipestance_root, &agg_csv_path,)
            .is_err());

        // path can be resolved
        assert_eq!(
            abs_path,
            resolve_path(
                Path::new("another_dir"), // relative to pipestance root
                &pipestance_root,
                &agg_csv_path
            )?
        );
        // path can be resolved but it's not a correct multi out
        assert!(FieldKind::MultiFolderPath
            .validate_and_resolve(
                "another_dir", // relative to pipestance root
                &pipestance_root,
                &agg_csv_path,
            )
            .is_err());

        // path can be resolved
        assert_eq!(
            abs_path,
            resolve_path(
                Path::new("pipe_root/another_dir"), // relative to agg.csv dir
                &pipestance_root,
                &agg_csv_path,
            )?
        );
        // path can be resolved but it's not a correct multi out
        assert!(FieldKind::MultiFolderPath
            .validate_and_resolve(
                "pipe_root/another_dir", // relative to agg.csv dir
                &pipestance_root,
                &agg_csv_path,
            )
            .is_err());

        // test 2: resolve a multi directory
        // test/parse_aggr_csv/path_resolution/
        //      csv_dir/
        //          agg.csv
        //      pipe_root/
        //          multi_dir_w_count/
        //              count/
        //                  molecule_info.h5
        let abs_path = cargo_dir
            .join("test_resources/parse_aggr_csv/path_resolution/pipe_root/multi_dir_w_count");
        // absolute path
        assert_eq!(
            abs_path.to_str().unwrap(),
            FieldKind::MultiFolderPath
                .validate_and_resolve(abs_path.to_str().unwrap(), &pipestance_root, &agg_csv_path)
                .unwrap()
        );
        // relative to pipestance root
        assert_eq!(
            abs_path.to_str().unwrap(),
            FieldKind::MultiFolderPath
                .validate_and_resolve("multi_dir_w_count", &pipestance_root, &agg_csv_path)
                .unwrap()
        );
        // relative to agg.csv dir
        assert_eq!(
            abs_path.to_str().unwrap(),
            FieldKind::MultiFolderPath
                .validate_and_resolve(
                    "pipe_root/multi_dir_w_count",
                    &pipestance_root,
                    &agg_csv_path
                )
                .unwrap()
        );

        // test 3: resolve a multi directory
        // test/parse_aggr_csv/path_resolution/
        //      csv_dir/
        //          agg.csv
        //      pipe_root/
        //          multi_dir_w_vdj/
        //              vdj_t/
        //                  vdj_contig_info.pb
        //              vdj_b/
        //                  vdj_contig_info.pb
        let abs_path = cargo_dir
            .join("test_resources/parse_aggr_csv/path_resolution/pipe_root/multi_dir_w_vdj");
        // absolute path
        assert_eq!(
            abs_path.to_str().unwrap(),
            FieldKind::MultiFolderPath
                .validate_and_resolve(abs_path.to_str().unwrap(), &pipestance_root, &agg_csv_path)
                .unwrap()
        );
        // relative to pipestance root
        assert_eq!(
            abs_path.to_str().unwrap(),
            FieldKind::MultiFolderPath
                .validate_and_resolve("multi_dir_w_vdj/", &pipestance_root, &agg_csv_path)
                .unwrap()
        );
        // relative to agg.csv dir
        assert_eq!(
            abs_path.to_str().unwrap(),
            FieldKind::MultiFolderPath
                .validate_and_resolve(
                    "pipe_root/multi_dir_w_vdj/",
                    &pipestance_root,
                    &agg_csv_path
                )
                .unwrap()
        );
        // test 4: resolve a multi directory
        // test/parse_aggr_csv/path_resolution/
        //      csv_dir/
        //          agg.csv
        //      pipe_root/
        //          multi_dir_w_empty_subdirs/
        //              count/
        //              vdj_b/
        //              vdj_t/
        let abs_path = cargo_dir.join(
            "test_resources/parse_aggr_csv/path_resolution/pipe_root/multi_dir_w_empty_subdirs",
        );
        // absolute path
        assert!(FieldKind::MultiFolderPath
            .validate_and_resolve(abs_path.to_str().unwrap(), &pipestance_root, &agg_csv_path)
            .is_err());
        // relative to pipestance root
        assert!(FieldKind::MultiFolderPath
            .validate_and_resolve(
                "multi_dir_w_empty_subdirs/",
                &pipestance_root,
                &agg_csv_path
            )
            .is_err());
        // relative to agg.csv dir
        assert!(FieldKind::MultiFolderPath
            .validate_and_resolve(
                "pipe_root/multi_dir_w_empty_subdirs/",
                &pipestance_root,
                &agg_csv_path
            )
            .is_err());
        Ok(())
    }

    #[inline]
    fn build_sample_def_snapshot<T>(parser: T, row: HashMap<String, String>)
    where
        T: AggrCsvParser,
        T::SampleDef: Serialize,
    {
        // Settings for insta to force sorting of maps
        let mut settings = insta::Settings::clone_current();
        settings.set_sort_maps(true);
        settings.bind(|| assert_json_snapshot!(parser.build_sample_def(row).unwrap()));
    }

    #[test]
    fn test_build_sample_def() {
        // Testing only with valid rows

        // Count

        // Hashmap with batch
        build_sample_def_snapshot(
            CountParser,
            hashmap![
                SAMPLE_ID_HEADER => "ABC1",
                MOLECULE_H5_HEADER => "AAA.h5",
                BATCH_HEADER => "b1"
            ],
        );

        // Hashmap without batch
        build_sample_def_snapshot(
            CountParser,
            hashmap![SAMPLE_ID_HEADER => "ABC1", MOLECULE_H5_HEADER => "AAA.h5"],
        );

        // Hashmap with batch and extra keys
        build_sample_def_snapshot(
            CountParser,
            hashmap![
                SAMPLE_ID_HEADER => "ABC1",
                MOLECULE_H5_HEADER => "AAA.h5",
                BATCH_HEADER => "b1",
                "key1" => "val1",
                "key2" => "val2"
            ],
        );

        // Hashmap with only extra keys and no batch
        build_sample_def_snapshot(
            CountParser,
            hashmap![
                SAMPLE_ID_HEADER => "ABC1",
                MOLECULE_H5_HEADER => "AAA.h5",
                "key1" => "val1",
                "key2" => "val2"
            ],
        );

        // VDJ

        // Hashmap without extra keys
        build_sample_def_snapshot(
            VdjParser,
            hashmap![
                SAMPLE_ID_HEADER => "ABC1",
                VDJ_CONTIG_INFO_HEADER => "AAA.proto",
                DONOR_HEADER => "d1",
                ORIGIN_HEADER => "o1"
            ],
        );

        // Hashmap with extra keys
        build_sample_def_snapshot(
            VdjParser,
            hashmap![
                SAMPLE_ID_HEADER => "ABC1",
                VDJ_CONTIG_INFO_HEADER => "AAA.proto",
                DONOR_HEADER => "d1",
                ORIGIN_HEADER => "o1",
                "key1" => "v1",
                "key2" => "v2"
            ],
        );

        // Multi

        // Hashmap for multi with only count
        build_sample_def_snapshot(
            MultiParser,
            hashmap![
                SAMPLE_ID_HEADER => "ABC1",
                MULTI_OUTS => "multi_out/",
                DONOR_HEADER => "d1",
                ORIGIN_HEADER => "o1"
            ],
        );

        // Hashmap for multi with vdj, count, and extra keys
        build_sample_def_snapshot(
            MultiParser,
            hashmap![
                SAMPLE_ID_HEADER => "ABC1",
                MULTI_OUTS => "multi_out/",
                DONOR_HEADER => "d1",
                ORIGIN_HEADER => "o1",
                BATCH_HEADER => "b1",
                "key1" => "v1",
                "key2" => "v2"
            ],
        );
    }

    #[test]
    fn test_text_field_validation() {
        assert_snapshot!(validate_text_value("Hello>world").unwrap_err());
        assert_eq!(validate_text_value("hEllO_123").unwrap(), "hEllO_123");
        assert!(validate_text_value("hE`llO{}123").is_err());
    }

    #[test]
    fn test_legacy_library_id_error() {
        // The error user would see if they use library_id
        err_snapshot!("test_resources/parse_aggr_csv/legacy_library_id.csv");
    }

    #[test]
    fn test_legacy_library_outs_error() {
        // The error user would see if they use library_id and library_outs
        err_snapshot!("test_resources/parse_aggr_csv/legacy_library_outs.csv");
    }

    #[test]
    fn test_invalid_multi_v6_outs() {
        // User specifies the outs folder of the 6.0 multi pipeline
        err_snapshot!("test_resources/parse_aggr_csv/pipes/csvs/invalid_multi_6.csv");
    }

    #[test]
    fn test_valid_multi_v6() -> Result<()> {
        let inp = ParseAggrCsvStageInputs::test(
            "test_resources/parse_aggr_csv/pipes/csvs/valid_multi_2.csv",
        );
        let out = ParseAggrCsv.test_run_tmpdir(inp)?;
        let expected_count_sample_def = vec![
            CountLibrary {
                library_id: "sample1".to_string(),
                molecule_h5: std::env::current_dir()?
                .join(
                    "test_resources/parse_aggr_csv/pipes/run10_outs/per_sample_outs/sample1/count/sample_molecule_info.h5",
                ),
                batch: None,
                meta: HashMap::default(),
            },
            CountLibrary {
                library_id: "sample2".to_string(),
                molecule_h5: std::env::current_dir()?
                .join(
                    "test_resources/parse_aggr_csv/pipes/run10_outs/per_sample_outs/sample2/count/sample_molecule_info.h5",
                ),
                batch: None,
                meta: HashMap::default(),
            },
        ];

        assert_eq!(out.count_libraries, expected_count_sample_def);
        assert!(out.vdj_aggr_inputs.is_empty());
        assert!(!out.disable_count_aggr);
        assert!(out.disable_vdj_aggr);
        Ok(())
    }
}
