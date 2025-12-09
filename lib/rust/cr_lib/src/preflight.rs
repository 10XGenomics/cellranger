#![deny(missing_docs)]

use anyhow::{Context, Result, bail, ensure};
use bio::alphabets::dna::revcomp;
use cr_h5::probe_reference_io::PROBE_DATA_LEN;
use cr_types::TargetingMethod;
use cr_types::probe_set::ProbeSetReferenceMetadata;
use cr_types::reference::feature_reference::{FeatureReference, FeatureType};
use cr_types::reference::probe_set_reference::TargetSetFile;
use cr_types::reference::reference_info::{MULTI_GENOME_SEPARATOR, ReferenceInfo};
use cr_types::types::FeatureBarcodeType;
use itertools::Itertools;
use lazy_static::lazy_static;
use martian::MartianRover;
use metric::{TxHashMap, TxHashSet};
use multi::config::preflight::check_file;
use regex::bytes::Regex;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use transcriptome::Transcriptome;
use vdj_reference::VdjReference;
use vdj_types::VdjRegion;

const MIN_NOFILE: u64 = 1024;

pub fn hostname() -> String {
    if let Ok(n) = ::hostname::get() {
        return n.to_string_lossy().into_owned();
    }
    "<unknown hostname>".to_string()
}

/// Check the global open file handle limit on Linux. Do nothing on macOS.
#[cfg(not(target_os = "linux"))]
fn check_global_file_handle_limit() -> Result<()> {
    Ok(())
}
#[cfg(target_os = "linux")]
fn check_global_file_handle_limit() -> Result<()> {
    const GLOBAL_NOFILE_PATH: &str = "/proc/sys/fs/file-max";
    const MIN_GLOBAL_NOFILE: i64 = 32_768;

    ensure!(
        Path::new(GLOBAL_NOFILE_PATH).is_file(),
        "On machine: {}, {GLOBAL_NOFILE_PATH} does not exist.",
        hostname()
    );

    let files_str = std::fs::read_to_string(GLOBAL_NOFILE_PATH)?;
    let files: i64 = files_str
        .split_ascii_whitespace()
        .next()
        .unwrap_or_default()
        .parse()
        .map_err(|_| {
            anyhow::anyhow!(
                "On machine: {}, {GLOBAL_NOFILE_PATH} contains a non-integer global \
                 open file handle limit: {files_str}",
                hostname(),
            )
        })?;
    ensure!(
        files >= MIN_GLOBAL_NOFILE,
        "On machine: {}, global open file handle limit ({files}) is less than {MIN_GLOBAL_NOFILE}. \
         Please set the global file handle limit to {MIN_GLOBAL_NOFILE} before restarting the pipeline.",
        hostname(),
    );
    Ok(())
}

pub fn check_resource_limits() -> Result<()> {
    if let Some(nofile) = rustix::process::getrlimit(rustix::process::Resource::Nofile).maximum {
        ensure!(
            nofile >= MIN_NOFILE,
            "On machine: {}, process open file handle hard limit ({nofile}) is less than \
             {MIN_NOFILE}. Please run `ulimit -n {MIN_NOFILE}` before restarting the pipeline.",
            hostname(),
        );
    }

    check_global_file_handle_limit()
}

/// Required probe set CSV metadata.
const REQUIRED_METADATA_ARRAY: [&str; 4] = [
    "panel_name",
    "panel_type",
    "reference_genome",
    "reference_version",
];

/// Probe set CSV headers.
const HEADERS: [&str; 8] = [
    "gene_id",
    "bait_seq",
    "bait_id",
    "probe_seq",
    "probe_id",
    "included",
    "region",
    "gene_name",
];

/// Minimum number of records in a target set.
const MIN_TARGET_PANEL_LENGTH: usize = 10;

lazy_static! {
    static ref REQUIRED_METADATA: TxHashSet<&'static str> =
        REQUIRED_METADATA_ARRAY.iter().copied().collect();
    static ref ALLOWED_HEADERS: TxHashSet<&'static str> = HEADERS.iter().copied().collect();
    static ref BAIT_SEQ: Regex = Regex::new("^[ACGTN-]+$").unwrap(); // "-" used as separator between probe halves (and if present gap)
}

fn make_required_headers(method: TargetingMethod) -> TxHashSet<&'static str> {
    if method == TargetingMethod::TemplatedLigation {
        ["gene_id", "probe_seq", "probe_id"]
            .iter()
            .copied()
            .collect()
    } else {
        ["gene_id"].iter().copied().collect()
    }
}

/// Validate the length and content of the gene_id and probe_id columns of the probe set reference CSV.
fn validate_probe_field(field_name: &str, field_value: &str) -> Result<()> {
    ensure!(
        field_value.len() <= PROBE_DATA_LEN,
        "{field_name} must not be longer than {PROBE_DATA_LEN} characters: \"{field_value}\""
    );
    ensure!(
        field_value.chars().all(|c| c.is_ascii_graphic()),
        "{field_name} must contain only ASCII characters: \"{field_value}\""
    );
    Ok(())
}

/// Validate that a target set is valid.
///
/// If enforce_min_panel_size is true, check that the user has supplied enough
/// probes for a valid analysis.
///
/// Return a set of target gene IDs that are present in the probe set CSV file.
pub fn validate_target_panel(
    rover: &MartianRover,
    transcriptome: Option<&Transcriptome>,
    ref_info: Option<&ReferenceInfo>,
    probe_set: &TargetSetFile,
    targeting_method: TargetingMethod,
    enforce_min_panel_size: bool,
) -> Result<TxHashSet<String>> {
    use TargetingMethod::{HybridCapture, TemplatedLigation};

    ensure!(
        check_file("The probe set CSV file", probe_set).is_ok(),
        "The probe set CSV file is either missing or not readable from {}: {}",
        hostname(),
        probe_set.display(),
    );

    // Read the transcriptome GTF to map gene names to gene IDs, used for hybcap.
    let gene_ids: TxHashSet<_> = transcriptome
        .map(|transcriptome| transcriptome.genes.iter().map(|x| x.id.as_str()).collect())
        .unwrap_or_default();
    let gene_name_to_id: TxHashMap<_, _> = if let Some(transcriptome) = transcriptome {
        transcriptome
            .genes
            .iter()
            .map(|x| (x.name.as_str(), x.id.as_str()))
            .collect()
    } else {
        TxHashMap::default()
    };

    let metadata = ProbeSetReferenceMetadata::load_from(probe_set)?;
    let metadata_fields = metadata.keys().map(String::as_str).collect();
    let missing_metadata_fields: Vec<_> = REQUIRED_METADATA
        .difference(&metadata_fields)
        .sorted()
        .collect();
    ensure!(
        missing_metadata_fields.is_empty(),
        "The following metadata fields are required in the probe set CSV header: \"{}\", but were \
         not found. Please include these fields in #field=value format at the top of the file.",
        missing_metadata_fields.into_iter().join("\", \""),
    );

    let panel_name = metadata.panel_name();
    ensure!(
        !panel_name.contains('/'),
        "The character \"/\" cannot appear in the probe set CSV panel name: {panel_name}"
    );

    let (file_format, file_format_version) = match (
        metadata.get("target_panel_file_format"),
        metadata.get("probe_set_file_format"),
    ) {
        (None, None) => bail!(
            "The probe_set_file_format metadata field is required in the probe set CSV file header."
        ),
        (Some(version), None) => ("target_panel_file_format", version),
        (None, Some(version)) => ("probe_set_file_format", version),
        (Some(_), Some(_)) => bail!(
            "The probe set CSV file specifes both target_panel_file_format and probe_set_file_format."
        ),
    };

    match (file_format, targeting_method) {
        ("target_panel_file_format", TemplatedLigation) => bail!(
            "The probe set CSV requires the metadata field `probe_set_file_format` in the header. \
             Please include this field in #field=value format at the top of the probe set file."
        ),
        ("probe_set_file_format", HybridCapture) => bail!(
            "The target panel CSV requires the metadata field `target_panel_file_format` in the header. \
             Please include this field in #field=value format at the top of the target set file."
        ),
        _ => (),
    }
    match (
        file_format_version.parse::<f64>(),
        file_format_version.contains('.'),
    ) {
        (Ok(v), true) if v > 3.0 => bail!(
            "The probe set CSV file contains unknown #{file_format}={file_format_version}. Must be 3.0 or less."
        ),
        (Err(_), true) | (Ok(_), false) => bail!(
            "The probe set CSV file contains an invalid value for the #{file_format}={file_format_version}. \
             It must conform to the format X.Y, such as 1.0."
        ),
        _ => (),
    }

    if let (Some(ref_info), TemplatedLigation) = (ref_info, targeting_method) {
        for probe_set_genome in metadata.reference_genome().split(MULTI_GENOME_SEPARATOR) {
            ensure!(
                ref_info
                    .genomes
                    .iter()
                    .any(|g| g.as_str() == probe_set_genome),
                "Probe set CSV reference genome \"{probe_set_genome}\" (split from \
                 \"{}\") does not match reference genome(s) in the reference: {:?}",
                metadata.reference_genome(),
                ref_info.genomes,
            );
        }

        if let Some(reference_version) = &ref_info.version {
            let probe_set_version = metadata.reference_version();
            if probe_set_version != reference_version {
                let msg = format!(
                    "Reference version \"{reference_version}\" does not match \
                     probe set CSV reference version \"{probe_set_version}\""
                );
                if ref_info.genomes.len() >= 2 {
                    rover.alarm(&msg).unwrap();
                } else {
                    bail!("{msg}");
                }
            }
        }
    }

    // Read the probe set CSV file.
    let mut reader = csv::ReaderBuilder::new()
        .comment(Some(b'#'))
        .from_path(probe_set)
        .with_context(|| probe_set.display().to_string())?;

    let headers: Vec<_> = {
        let headers = reader.headers()?;
        let duplicates: Vec<_> = headers.iter().duplicates().collect();
        ensure!(
            duplicates.is_empty(),
            "The probe set CSV file header contains a duplicate field: \"{}\". \
             Please ensure that the file has no duplicated fields.",
            duplicates.join(", ")
        );
        headers.iter().collect()
    };

    let required_headers = make_required_headers(targeting_method);
    let headers_set = headers.iter().copied().collect();
    let missing_headers: Vec<_> = required_headers.difference(&headers_set).sorted().collect();
    ensure!(
        missing_headers.is_empty(),
        "The probe set CSV file header does not contain one or more required comma-separated \
         fields: \"{}\". The following required fields were found: \"{}\". \
         Please check that the file is in CSV format and has the required field names.",
        missing_headers.into_iter().join("\", \""),
        required_headers
            .intersection(&headers_set)
            .sorted()
            .join("\", \""),
    );

    let extra_headers: Vec<_> = headers_set.difference(&ALLOWED_HEADERS).sorted().collect();
    ensure!(
        extra_headers.is_empty(),
        "The probe set CSV file header contains invalid columns: \"{}\". \
         Only the following comma-delimited fields may be used in the file header: \"{}\"",
        extra_headers.into_iter().join("\", \""),
        HEADERS.join("\", \""),
    );

    if file_format == "probe_set_file_format" {
        let header_row = &headers[0..3];
        match header_row {
            ["gene_id", "probe_seq", "probe_id"] => (),
            ["gene_id", "bait_seq", "bait_id"] => (),
            _ => {
                bail!(
                    "The probe set CSV file headers are incorrectly specified: \
                    \"{}\". Expected order of required first three headers is: \
                    \"gene_id, probe_seq, probe_id\".",
                    header_row.join(", ")
                )
            }
        }
    }

    let hdr_to_idx: TxHashMap<_, _> = headers
        .into_iter()
        .enumerate()
        .map(|(i, h)| (h.to_string(), i))
        .collect();

    let gene_id_idx = hdr_to_idx["gene_id"];
    let mut total_records = 0;
    let mut excluded_records = 0;
    let mut probe_ids = TxHashSet::default();
    let mut target_gene_ids = TxHashSet::default();
    for (i, record) in reader.records().enumerate() {
        let row = 2 + metadata.len() + i; // 1 for 0-indexing, 1 for header
        if let Err(ref error) = record {
            match error.kind() {
                &csv::ErrorKind::UnequalLengths {
                    pos: _,
                    expected_len: _,
                    len,
                } => {
                    if len as usize > hdr_to_idx.len() {
                        bail!(
                            "The probe set CSV file contains more columns than the header on row \
                             {row}. Please use a CSV file with a header for each column and a \
                             value for each column in each row. You may have a comma character in \
                             a field. Commas are permitted in some fields, but fields containing \
                             commas must be enclosed in quotes.",
                        );
                    }
                    bail!(
                        "The probe set CSV file contains an empty column or fewer columns than \
                         the header on row {row}. You might have a missing comma. \
                         Please use a CSV file with a header for each column and a value for \
                         each column in each row.",
                    );
                }
                _ => {
                    bail!("The probe set CSV file failed to parse on row {row}.")
                }
            }
        }

        let record = record?;
        if let Some(gene_id) = record.get(gene_id_idx) {
            validate_probe_field("gene_id", gene_id)?;
            let targeted_gene_id = match targeting_method {
                // Custom probes may target genes not found in the reference transcriptome.
                TemplatedLigation => gene_id,
                HybridCapture => {
                    let gene = gene_id.split('.').next().unwrap();
                    if transcriptome.is_none() {
                        gene_id
                    } else if gene_ids.contains(gene) {
                        gene
                    } else if let Some(&gene_id) = gene_name_to_id.get(gene) {
                        gene_id
                    } else {
                        bail!(
                            "The target panel CSV file contains a record to a gene not seen \
                             in the reference: \"{gene}\""
                        );
                    }
                }
            };
            target_gene_ids.insert(targeted_gene_id.to_string());
        }

        match (hdr_to_idx.get("bait_seq"), hdr_to_idx.get("probe_seq")) {
            (None, None) => (),
            (Some(&idx), None) | (None, Some(&idx)) => {
                if let Some(probe_seq) = record.get(idx) {
                    ensure!(
                        BAIT_SEQ.is_match(probe_seq.as_bytes()),
                        "The probe set CSV file contains an invalid probe sequence on row {row}: \
                         \"{probe_seq}\". May only contain \"ACGTN-\"",
                    );
                }
            }
            (Some(_), Some(_)) => {
                bail!("The probe set CSV file contains both bait_seq and probe_seq.")
            }
        }

        match (hdr_to_idx.get("bait_id"), hdr_to_idx.get("probe_id")) {
            (None, None) => (),
            (Some(&idx), None) | (None, Some(&idx)) => {
                if let Some(probe_id) = record.get(idx) {
                    validate_probe_field("probe_id", probe_id)?;
                    ensure!(
                        probe_ids.insert(probe_id.to_string()),
                        "The probe set CSV file contains a duplicate probe ID on row {row}: \
                         \"{probe_id}\". All entries in the probe_id column must be unique.",
                    );
                }
            }
            (Some(_), Some(_)) => {
                bail!("The probe set CSV file contains both bait_id and probe_id.")
            }
        }

        match hdr_to_idx.get("included") {
            None => (),
            Some(&idx) => {
                if let Some(included) = record.get(idx) {
                    if !["true", "false"].contains(&included.to_lowercase().as_str()) {
                        bail!(
                            r#"The column "included" must be "true" or "false" but saw "{included}""#
                        );
                    }
                    if included.to_lowercase().as_str() == "false" {
                        excluded_records += 1;
                    };
                }
            }
        }

        if let Some(&gene_name_index) = hdr_to_idx.get("gene_name") {
            let gene_id = &record[gene_id_idx];
            let Some(gene_name) = record.get(gene_name_index) else {
                bail!("gene_name must not be empty for {gene_id} on row {row}");
            };
            ensure!(
                !gene_name.is_empty(),
                "gene_name must not be empty for {gene_id} on row {row}",
            );
            validate_probe_field("gene_name", gene_name)?;
            if let Some(transcriptome) = transcriptome
                && let Some(gene_index) = transcriptome.gene_id_to_idx.get(gene_id)
            {
                let transcriptome_gene_name =
                    transcriptome.genes[gene_index.0 as usize].name.as_str();
                ensure!(
                    gene_name == transcriptome_gene_name,
                    "The gene_name {gene_name} of gene_id {gene_id} on row {row} does not \
                         match the gene name {transcriptome_gene_name} in the transcriptome",
                );
            }
        }

        total_records += 1;
    }

    if enforce_min_panel_size {
        ensure!(
            total_records >= MIN_TARGET_PANEL_LENGTH,
            "{MIN_TARGET_PANEL_LENGTH} or more records must be specified in the probe set CSV file \
            for compatibility with downstream analysis. Number of records found: {total_records}."
        );

        ensure!(
            (total_records - excluded_records) >= MIN_TARGET_PANEL_LENGTH,
            "{MIN_TARGET_PANEL_LENGTH} or more records must be specified in the probe set CSV file \
            for compatibility with downstream analysis. Number of records found: {total_records}. \
            Number of records set to false for the `included` field: {excluded_records}",
        );
    }

    Ok(target_gene_ids)
}

lazy_static! {
    static ref SPECIAL_IDS: TxHashSet<&'static str> = {
        let mut s = TxHashSet::default();
        s.insert("none");
        s.insert("non-targeting");
        s.insert("ignore");
        s
    };
}

pub fn check_crispr_target_genes(
    transcriptome: &Transcriptome,
    feature_ref: &FeatureReference,
    target_gene_ids: Option<&TxHashSet<String>>,
) -> Result<()> {
    let gene_id_to_name: TxHashMap<_, _> = transcriptome
        .genes
        .iter()
        .map(|x| (x.id.as_str(), x.name.as_str()))
        .collect();
    for feat in &feature_ref.feature_defs {
        if feat.feature_type != FeatureType::Barcode(FeatureBarcodeType::Crispr) {
            continue;
        }
        if let Some(target_id) = feat.tags.get("target_gene_id").map(String::as_str) {
            if SPECIAL_IDS.contains(target_id.to_ascii_lowercase().as_str()) {
                continue;
            }
            if let Some(gene_name) = gene_id_to_name.get(target_id) {
                target_gene_ids
                    .map(|ids| {
                        ensure!(
                            ids.contains(target_id),
                            "CRISPR: You specified target_gene_id = \"{target_id}\" as a guide RNA \
                             in the feature reference, but this gene is not specified in the \
                             gene_id column of the probe set CSV file.",
                        );
                        Ok(())
                    })
                    .transpose()?;
                match feat.tags.get("target_gene_name") {
                    Some(target_name) if !target_name.is_empty() => {
                        if target_name != gene_name {
                            bail!(
                                "CRISPR: You specified target_gene_id = \"{target_id}\" and \
                                 target_gene_name = \"{target_name}\" in the feature reference, but the \
                                 transcriptome has gene_id = \"{target_id}\" with name = \"{gene_name}\". Please \
                                 ensure the target_gene_name field has the correct gene name that \
                                 matches the transcriptome."
                            );
                        }
                    }
                    _ => bail!(
                        "CRISPR: No target_gene_name specified for target_gene_id = \"{target_id}\" in the feature reference."
                    ),
                }
            } else {
                bail!(
                    "CRISPR: A target_gene_id (\"{target_id}\") declared for one or more guide RNAs in the \
                     feature reference does not exist in the transcriptome. Please specify a \
                     target_gene_id that exists in the reference, or use the string \
                     \"Non-Targeting\" to indicate a control guide."
                );
            }
        }
    }
    Ok(())
}

lazy_static! {
    static ref PRIMER_SEQ: Regex = Regex::new(r"^[ACGT]+$").unwrap();
}

pub fn check_vdj_inner_enrichment_primers(
    vdj_ref_path: Option<&Path>,
    vdj_ref: Option<&VdjReference>,
    inner_enrichment_primers: &Path,
) -> Result<Vec<String>> {
    if !inner_enrichment_primers.is_file() {
        bail!(
            "The file specifying inner enrichment primers ({}) does not exist or is not readable. Please check the path on machine {}.",
            inner_enrichment_primers.display(),
            hostname(),
        );
    }
    let mut rdr = BufReader::new(File::open(inner_enrichment_primers)?);
    let mut buf = String::new();

    let mut primers = vec![];
    let mut lineno = 1;
    while rdr.read_line(&mut buf)? > 0 {
        let primer = buf.trim_end().to_string();
        if primer.is_empty() {
            bail!(
                "Line number {} in the inner enrichment primers file ({}) is empty. \
                 Please specify a newline separated list of primers.",
                lineno,
                inner_enrichment_primers.display(),
            );
        }
        for (col, base) in primer.chars().enumerate() {
            if "ACGT".chars().all(|x| x != base) {
                bail!(
                    "Inner enrichment primers file ({}) contains non-ACGT characters, which are \
                     not supported (Found {} on line {}, character {}). Please specify a newline \
                     separated list of primers.",
                    inner_enrichment_primers.display(),
                    base,
                    lineno,
                    col + 1,
                );
            }
        }
        primers.push(primer);
        buf.clear();
        lineno += 1;
    }

    if let Some(vdj_ref) = vdj_ref {
        let mut invalid_primers = vec![];
        for primer in &primers {
            let primer_rc = revcomp(primer.as_bytes());
            let mut found = false;
            for entry in vdj_ref.iter_region_filtered(VdjRegion::C) {
                // the following is O(n^2), but fast for short things (which we expect here)
                // the rust standard library provides no good way to do this, so we do it ourselves
                if entry
                    .sequence
                    .windows(primer_rc.len())
                    .any(|s| s == primer_rc)
                {
                    found = true;
                    break;
                }
            }
            if !found {
                invalid_primers.push(primer);
            }
        }
        if !invalid_primers.is_empty() {
            bail!(
                "None of the C-REGIONs in the reference ({}) is targeted by the following inner enrichment primer(s): {}",
                vdj_ref_path.unwrap().display(),
                invalid_primers.into_iter().join(", "),
            );
        }
    }
    Ok(primers)
}

const KNOWN_PRIMERS: [&str; 23] = [
    // Human TCR primers
    "AGTCTCTCAGCTGGTACACG",
    "TCTGATGGCTCAAACACAGC",
    // Human IG primers
    "GGGAAGTTTCTGGCGGTCA",
    "GGTGGTACCCAGTTATCAAGCAT",
    "GTGTCCCAGGTCACCATCAC",
    "TCCTGAGGACTGTAGGACAGC",
    "CACGCTGCTCGTATCCGA",
    "TAGCTGCTGGCCGC",
    "GCGTTATCCACCTTCCACTGT",
    // Mouse TCR primers
    "AGTCAAAGTCGGTGAACAGGCA",
    "GGCCAAGCACACGAGGGTA",
    // Mouse IG primers
    "TACACACCAGTGTGGCCTT",
    "CAGGCCACTGTCACACCACT",
    "CAGGTCACATTCATCGTGCCG",
    "GAGGCCAGCACAGTGACCT",
    "GCAGGGAAGTTCACAGTGCT",
    "CTGTTTGAGATCAGTTTGCCATCCT",
    "TGCGAGGTGGCTAGGTACTTG",
    "CCCTTGACCAGGCATCC",
    "AGGTCACGGAGGAACCAGTTG",
    "GGCATCCCAGTGTCACCGA",
    "AGAAGATCCACTTCACCTTGAAC",
    "GAAGCACACGACTGAGGCAC",
];

pub fn check_vdj_known_enrichment_primers(
    vdj_ref_path: &Path,
    vdj_ref: &VdjReference,
) -> Result<()> {
    for &primer in &KNOWN_PRIMERS {
        let primer_rc = revcomp(primer.as_bytes());
        for entry in vdj_ref.iter_region_filtered(VdjRegion::C) {
            // the following is O(n^2), but fast for short things (which we expect here)
            // the rust standard library provides no good way to do this, so we do it ourselves
            if entry
                .sequence
                .windows(primer_rc.len())
                .any(|s| s == primer_rc)
            {
                return Ok(());
            }
        }
    }
    bail!(
        "Inner enrichment primers are required for species other than human or mouse for which \
         primers are not provided by 10x Genomics. None of the constant regions in the reference \
         ({}) is targeted by the known primers.",
        vdj_ref_path.display(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use glob::glob;
    use insta::assert_debug_snapshot;
    use martian::Resource;
    use std::path::PathBuf;
    use tempfile::tempdir;

    fn rover() -> MartianRover {
        let tmpdir = tempdir().unwrap();
        MartianRover::new(
            tmpdir.path(),
            Resource::new().mem_gb(1).vmem_gb(1).threads(1),
        )
    }

    #[test]
    fn test_invalid_target_panels() -> Result<()> {
        let target_panels: Vec<_> = glob("test/target_panels/invalid_csvs/*.csv")?.try_collect()?;
        let outs: Vec<_> = target_panels
            .iter()
            .sorted()
            .map(|target_panel| {
                (
                    target_panel,
                    validate_target_panel(
                        &rover(),
                        None,
                        None,
                        &TargetSetFile::from(target_panel),
                        TargetingMethod::HybridCapture,
                        true,
                    ),
                )
            })
            .collect();
        assert_debug_snapshot!(outs);
        Ok(())
    }

    #[test]
    fn test_valid_target_panels_two_genomes() -> Result<()> {
        let target_panels = [
            "test/probe_sets/GRCh38-fmt3-refv24.csv",
            "test/probe_sets/mm10-fmt3-refv20.csv",
        ]
        .map(PathBuf::from);
        let ref_info = ReferenceInfo::from_reference_path(Path::new(
            "test/reference/GRCh38-and-mm10_ref_tiny",
        ))?;
        let outs: Vec<_> = target_panels
            .iter()
            .sorted()
            .map(|target_panel| {
                (
                    target_panel,
                    validate_target_panel(
                        &rover(),
                        None,
                        Some(&ref_info),
                        &TargetSetFile::from(target_panel),
                        TargetingMethod::TemplatedLigation,
                        true,
                    ),
                )
            })
            .collect();
        assert_debug_snapshot!(outs);
        Ok(())
    }

    #[test]
    fn test_vdj_enrichment_primers() {
        let vdj_ref_path =
            Path::new("../dui_tests/test_resources/reference/vdj_GRCh38_alts_ensembl-4.0.0");
        let vdj_ref = &VdjReference::from_reference_folder(vdj_ref_path).unwrap();
        assert!(
            check_vdj_inner_enrichment_primers(
                Some(vdj_ref_path),
                Some(vdj_ref),
                Path::new("test/preflight/vdj_human_t_inner_primers.txt")
            )
            .is_ok()
        );
        assert!(
            check_vdj_inner_enrichment_primers(
                Some(vdj_ref_path),
                Some(vdj_ref),
                Path::new("test/preflight/vdj_invalid_primers.txt")
            )
            .is_err()
        );
    }
}
