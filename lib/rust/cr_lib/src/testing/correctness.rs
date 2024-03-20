#![allow(dead_code, unused_variables)]

use crate::aligner::BarcodeSummary;
use crate::testing::diff_metrics;
use anyhow::{Context, Result};
use cr_bam::bam::BamPosSort;
use cr_types::{FeatureBarcodeType, MetricsFile};
use itertools::{zip_eq, Itertools};
use martian_filetypes::tabular_file::CsvFile;
use martian_filetypes::{FileTypeRead, LazyFileTypeIO};
use metric::join_metric_name;
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::{self, Read, Record};
use serde_json::value::Value;
use shardio::SortKey;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::path::Path;
use std::process::Command;

// ##################################################################################
// Correctness checks for the metric summary json
// ##################################################################################
pub fn check_metrics_correctness(actual: &MetricsFile, expected: &MetricsFile) -> Result<()> {
    let mut actual_metrics: HashMap<_, _> = actual.read()?;
    let mut expected_metrics: HashMap<_, _> = expected.read()?;

    // Remap these absurdly named metrics to their new equivalents
    let remapped_metrics = [
        (
            "multi_transcriptome_conf_mapped_reads_frac",
            Some("recognized_feature_bc_frac"),
        ),
        ("multi_transcriptome_conf_mapped_barcoded_reads_frac", None),
    ];
    for feature_barcode_type in [
        FeatureBarcodeType::Antibody,
        FeatureBarcodeType::Crispr,
        FeatureBarcodeType::Custom,
        FeatureBarcodeType::Multiplexing,
    ] {
        for (old, new) in remapped_metrics {
            if let Some(v) = expected_metrics.remove(&join_metric_name(feature_barcode_type, old)) {
                if let Some(new) = new {
                    expected_metrics.insert(join_metric_name(feature_barcode_type, new), v);
                }
            }
        }
    }

    // This metric was only computed for 3' on master
    if expected_metrics.get("tso_frac") == Some(&Value::Null) {
        actual_metrics.remove("tso_frac");
    }

    // This metric is not always present
    if !expected_metrics.contains_key("top_read_prefixes") {
        actual_metrics.remove("top_read_prefixes");
    }

    // Retain only keys present in actual_metrics
    expected_metrics.retain(|k, _| actual_metrics.contains_key(k));

    // Filter keys in actual_metrics which are not there in expected_metrics
    // and equal to "NaN"
    actual_metrics
        .retain(|k, v| expected_metrics.contains_key(k) || v != &Value::String("NaN".into()));

    diff_metrics(actual_metrics, expected_metrics);
    Ok(())
}

// ##################################################################################
// Correctness checks for the barcode_correction.csv
// ##################################################################################
pub fn check_barcode_csv_correctness(
    actual: &CsvFile<BarcodeSummary>,
    expected: &CsvFile<BarcodeSummary>,
) -> Result<()> {
    let actual_rows: HashSet<_> = actual.read_all()?;
    let expected_rows: HashSet<_> = expected.read_all()?;
    assert!(
        actual_rows == expected_rows,
        "Barcode csv correctness failed!"
    );
    println!(" > Barcode CSV files pass the correctness test!");
    Ok(())
}

// ##################################################################################
// Correctness checks for the gzipped files storing the mtx
// ##################################################################################
pub fn check_mtx_correctness(actual: &Path, expected: &Path) -> Result<()> {
    use flate2::read::GzDecoder;
    use std::io::{BufRead, BufReader};

    let actual_lines = BufReader::new(GzDecoder::new(File::open(actual)?)).lines();
    let expected_lines = BufReader::new(GzDecoder::new(File::open(expected)?)).lines();

    for (a, e) in actual_lines.zip_eq(expected_lines) {
        let a = a?;
        let e = e?;
        if a.starts_with("%metadata_json") {
            assert!(
                e.starts_with("%metadata_json"),
                "Mtx check failed at metadata"
            );
        } else {
            assert_eq!(a, e);
        }
    }
    println!(" > Mtx files pass the correctness test!");
    Ok(())
}

// ##################################################################################
// Correctness checks for the h5 files using h5diff
// TODO: Replace this with crate::h5::compare::H5Compare
// ##################################################################################
pub fn check_h5_correctness(actual: &Path, expected: &Path) -> Result<()> {
    let output = Command::new("h5diff")
        .arg("-cr")
        .arg(actual.display().to_string())
        .arg(expected.display().to_string())
        .output()
        .with_context(|| "Error: h5diff")?;
    assert!(output.status.code().is_some());
    let stdout = String::from_utf8(output.stdout).unwrap();
    let stderr = String::from_utf8(output.stderr).unwrap();
    if output.status.code().unwrap() > 1 {
        eprintln!("{stdout}{stderr}Error: h5diff failed");
        assert_eq!(output.status.code().unwrap(), 0);
    }
    for line in stdout.lines() {
        if line.ends_with(" differences found") && line != "0 differences found" {
            eprintln!(
                "{}{}Files {} and {} are not identical.",
                stdout,
                stderr,
                actual.display(),
                expected.display()
            );
            assert_eq!(line, "0 differences found");
        }
    }
    println!(
        "Files {} and {} are identical.",
        actual.display(),
        expected.display()
    );
    Ok(())
}

// ##################################################################################
// BAM file correctness is broken down into multiple smaller functions
// ##################################################################################

fn collect_tags(rec: &Record) -> Vec<(String, Option<Aux<'_>>)> {
    use cr_bam::bam_tags::{
        ANTISENSE_TAG, EXTRA_FLAGS_TAG, FEATURE_IDS_TAG, FEATURE_QUAL_TAG, FEATURE_RAW_TAG,
        FEATURE_SEQ_TAG, GENE_ID_TAG, GENE_NAME_TAG, MULTIMAPPER_TAG, PROC_BC_SEQ_TAG,
        PROC_UMI_SEQ_TAG, RAW_BARCODE_QUAL_TAG, RAW_BARCODE_SEQ_TAG, RAW_UMI_QUAL_TAG,
        RAW_UMI_SEQ_TAG, READ_GROUP_TAG, REGION_TAG, TRANSCRIPT_TAG, UNPAIRED_GENE_ID_TAG,
        UNPAIRED_GENE_NAME_TAG,
    };

    let tags_to_check = [
        READ_GROUP_TAG,
        PROC_BC_SEQ_TAG,
        PROC_UMI_SEQ_TAG,
        FEATURE_RAW_TAG,
        FEATURE_QUAL_TAG,
        FEATURE_SEQ_TAG,
        FEATURE_IDS_TAG,
        EXTRA_FLAGS_TAG,
        RAW_UMI_SEQ_TAG,
        RAW_UMI_QUAL_TAG,
        RAW_BARCODE_SEQ_TAG,
        RAW_BARCODE_QUAL_TAG,
        TRANSCRIPT_TAG,
        GENE_ID_TAG,
        GENE_NAME_TAG,
        REGION_TAG,
        MULTIMAPPER_TAG,
        ANTISENSE_TAG,
        UNPAIRED_GENE_ID_TAG,
        UNPAIRED_GENE_NAME_TAG,
    ];
    let mut result = Vec::new();
    for tag in tags_to_check {
        let tag_string = std::str::from_utf8(tag).unwrap().to_string();
        // From the spec: <https://samtools.github.io/hts-specs/SAMv1.pdf>
        // While all single (i.e., non-array) integer types are stored in SAM
        // as ‘i’, in BAM any of ‘cCsSiI’ may be used together with the
        // correspondingly-sized binary integer value, chosen according to
        // the field value’s magnitude
        //
        // So let's read map any int to Aux::I32 for comparison
        let aux = rec.aux(tag).ok().map(|aux| match aux {
            Aux::I8(s) => Aux::I32(s as i32),
            Aux::U8(s) => Aux::I32(s as i32),
            Aux::I16(s) => Aux::I32(s as i32),
            Aux::U16(s) => Aux::I32(s as i32),
            Aux::I32(s) => Aux::I32(s),
            Aux::U32(s) => Aux::I32(s as i32),
            a => a,
        });

        result.push((tag_string, aux));
    }
    result
}

fn bam_tag_display(tags: &[(String, Option<Aux<'_>>)]) -> Result<String> {
    use std::fmt::Write;
    let mut result = String::new();
    for (k, v) in tags {
        writeln!(&mut result, "  - TAG {k:4} VALUE {v:?}")?;
    }
    Ok(result)
}

fn fold_cigar_operations(cigar_view: bam::record::CigarStringView) -> HashMap<char, u32> {
    cigar_view
        .iter()
        .map(|c| (c.char(), c.len()))
        .fold(HashMap::new(), |mut acc, (k, v)| {
            *acc.entry(k).or_insert(0) += v;
            acc
        })
}

fn check_bam_record_correctness(slfe_rec: &Record, master_rec: &Record) -> Result<()> {
    assert_eq!(slfe_rec.qname(), master_rec.qname());
    assert_eq!(slfe_rec.flags(), master_rec.flags());
    assert_eq!(slfe_rec.pos(), master_rec.pos());
    assert_eq!(slfe_rec.mapq(), master_rec.mapq());
    assert_eq!(slfe_rec.mtid(), master_rec.mtid());
    assert_eq!(slfe_rec.mpos(), master_rec.mpos());
    assert_eq!(slfe_rec.insert_size(), master_rec.insert_size());
    if slfe_rec.raw_cigar() != master_rec.raw_cigar() {
        // Could be two alignments with the same alignment score
        assert_eq!(
            fold_cigar_operations(slfe_rec.cigar()),
            fold_cigar_operations(master_rec.cigar())
        );
    }

    assert_eq!(slfe_rec.seq().as_bytes(), master_rec.seq().as_bytes());
    assert_eq!(slfe_rec.qual(), master_rec.qual());
    let slfe_tags = collect_tags(slfe_rec);
    let master_tags = collect_tags(master_rec);
    let slfe_display = format!("SLFE TAGS\n{}", bam_tag_display(&slfe_tags)?);
    let master_display = format!("MASTER TAGS\n{}", bam_tag_display(&master_tags)?);
    for (s, m) in zip_eq(slfe_tags, master_tags) {
        // Master sets inconsistent UMIs for primary and secondary alignments
        // of a read with a low support UMI. So ignore the UMI correctness check for
        // secondary alignments.
        if slfe_rec.is_secondary() && s.0.as_bytes() == cr_bam::bam_tags::PROC_UMI_SEQ_TAG {
            continue;
        }
        assert_eq!(
            s,
            m,
            "Tag check failed on qname {}.\n{slfe_display}\n{master_display}",
            std::str::from_utf8(slfe_rec.qname()).unwrap(),
        );
    }
    Ok(())
}

pub fn check_bam_file_correctness(actual_bam: &Path, expected_bam: &Path) -> Result<()> {
    check_bam_header_correctness(actual_bam, expected_bam)?;

    let mut slfe_records: Vec<_> = bam::Reader::from_path(actual_bam)?
        .records()
        .try_collect()?;
    slfe_records.sort_by_key(|x| BamPosSort::sort_key(x).into_owned());

    let mut master_records: Vec<_> = bam::Reader::from_path(expected_bam)?
        .records()
        .try_collect()?;
    master_records.sort_by_key(|x| BamPosSort::sort_key(x).into_owned());

    for (slfe_rec, master_rec) in slfe_records
        .iter()
        .zip_eq(master_records.iter().filter(|x| !x.is_secondary()))
    {
        check_bam_record_correctness(slfe_rec, master_rec)?;
    }

    println!(" > BAM files pass the correctness test!");
    Ok(())
}

fn check_bam_header_correctness(actual_bam: &Path, expected_bam: &Path) -> Result<()> {
    let tags_to_check = ["@HD", "@SQ", "@RG", "@CO\t10x_bam_to_fastq:R"];
    let pred = |x: &&str| -> bool {
        tags_to_check
            .iter()
            .map(|tag| x.starts_with(tag))
            .any(|x| x)
    };

    let actual_header = String::from_utf8(
        bam::Reader::from_path(actual_bam)?
            .header()
            .as_bytes()
            .to_vec(),
    )?;
    let actual_lines: HashSet<_> = actual_header.lines().filter(pred).collect();

    let expected_header = String::from_utf8(
        bam::Reader::from_path(expected_bam)?
            .header()
            .as_bytes()
            .to_vec(),
    )?;
    let expected_lines: HashSet<_> = expected_header.lines().filter(pred).collect();

    if actual_lines != expected_lines {
        println!("-----------Lines present in actual, but not in expected-----------");
        for line in actual_lines.difference(&expected_lines) {
            println!("{line}");
        }
        println!("\n-----------Lines present in expected, but not in actual-----------");
        for line in expected_lines.difference(&actual_lines) {
            println!("{line}");
        }
        panic!("Bam header check failed!");
    }
    Ok(())
}
