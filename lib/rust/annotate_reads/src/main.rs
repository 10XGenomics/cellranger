//
// Copyright (c) 2017 10x Genomics, Inc. All rights reserved.
//

//#[macro_use] extern crate log;
extern crate histogram;
extern crate csv;
extern crate rust_htslib;
extern crate itertools;
extern crate num;
extern crate rand;
extern crate bincode;
extern crate docopt;
extern crate serde;
#[macro_use] extern crate serde_derive;
#[macro_use] extern crate serde_json;
extern crate flate2;
#[macro_use] extern crate bitflags;
extern crate bio;
extern crate lz4;

mod barcodes;
mod features;
mod metrics;
mod reference;
mod report;
mod transcriptome;
mod utils;

use std::path::Path;
use std::fs::File;
use std::io::{BufReader};
use std::str;
use std::collections::{HashSet, HashMap, BTreeMap};

use docopt::Docopt;

use rust_htslib::bam;
use rust_htslib::bam::Read;
use rust_htslib::bam::header::{Header, HeaderRecord};
use rust_htslib::bam::record::Record;

use serde_json::Value;

use itertools::Itertools;

use transcriptome::{TranscriptAnnotator, AnnotationParams, Strandedness, PairAnnotationData};
use report::{Reporter, Metrics, ReadData, RecordData};
use barcodes::{BarcodeUmiChecker, BarcodeUmiData};
use features::{FeatureChecker, FeatureData};
use metrics::MetricGroup;
use utils::CellrangerFastqHeader;

const USAGE: &'static str = "
Usage:
  annotate_reads main <in-bam> <in-tags> <out-bam> <out-metrics> <reference-path> <gene-index> <bc-counts> <bc-whitelist> <gem-group> <out-metadata> <strandedness> <feature-dist> <library-type> <library-id> <library-info> [--fiveprime] [--bam-comments=F] [--feature-ref=F]
  annotate_reads join <in-chunked-metrics> <out-json> <out-bc-csv>
  annotate_reads (-h | --help)

Options:
  -h --help            Show this screen.
  --neg-strand         Assume negative-strand chemistry.
  --fiveprime          Assume reads originate from 5' end (instead of 3').
  --bam-comments=F     JSON file with array of strings to append as @CO items
  --feature-ref=F      Feature definition file (CSV)
";

#[derive(Debug, Deserialize, Clone)]
struct Args {
    cmd_main:               bool,
    cmd_join:               bool,

    // main args
    arg_in_bam:             Option<String>,
    arg_in_tags:            Option<String>,
    arg_out_bam:            Option<String>,
    arg_out_metrics:        Option<String>,
    arg_reference_path:     Option<String>,
    arg_gene_index:         Option<String>,
    arg_bc_counts:          Option<String>,
    arg_bc_whitelist:       Option<String>,
    arg_gem_group:          Option<u8>,
    arg_out_metadata:       Option<String>,
    arg_strandedness:       Option<String>,
    arg_feature_dist:       Option<String>,
    arg_library_type:       Option<String>,
    arg_library_id:         Option<String>,
    arg_library_info:       Option<String>,

    // main flags
    flag_fiveprime:         bool,
    flag_bam_comments:      Option<String>,
    flag_feature_ref:       Option<String>,

    // join args
    arg_in_chunked_metrics: Option<String>,
    arg_out_json:           Option<String>,
    arg_out_bc_csv:         Option<String>,
}

const LIBRARY_INDEX_TAG: &'static str = "li";

#[derive(Deserialize, Serialize, Debug)]
pub struct LibraryInfo {
    library_id: String,
    library_type: String,
    gem_group: u64,
}

fn main() {
    let args: Args = Docopt::new(USAGE).and_then(|d| d.deserialize()).unwrap_or_else(|e| e.exit());
    if args.cmd_main {
        annotate_reads_main(args);
    } else {
        annotate_reads_join(args);
    }
}

fn annotate_reads_main(args: Args) {

    println!("Loading indices");
    let reference_path = args.arg_reference_path.unwrap();
    let strandedness = Strandedness::from_string(args.arg_strandedness.unwrap());
    let params = AnnotationParams {
        chemistry_strandedness:     strandedness,
        chemistry_fiveprime:        args.flag_fiveprime,
        intergenic_trim_bases:      0,
        intronic_trim_bases:        0,
        junction_trim_bases:        0,
        region_min_overlap:         0.5,
    };
    let annotator = TranscriptAnnotator::new(&reference_path, &args.arg_gene_index.clone().unwrap(), params);


    println!("Loading whitelist");
    let gem_group = &args.arg_gem_group.unwrap();

    // Attempt to translate barcodes only if the library type specifies that we should
    let translate_barcodes = match &args.arg_library_type {
        &Some(ref library_type) =>
            !features::LIBRARY_TYPES_WITHOUT_FEATURES.contains(&library_type.as_str()),
        _ => false,
    };

    let bc_umi_checker = barcodes::BarcodeUmiChecker::new(&args.arg_bc_counts.unwrap(),
                                                          &args.arg_bc_whitelist.unwrap(),
                                                          gem_group,
                                                          translate_barcodes,
                                                          &args.arg_library_type.as_ref().unwrap());


    let feature_checker = match (&args.flag_feature_ref,
                                 &args.arg_feature_dist,
                                 &args.arg_library_type) {
        (&Some(ref fref_path), &Some(ref fdist_path), &Some(ref library_type)) => {
            if features::library_type_requires_feature_ref(library_type) {
                println!("Loading feature reference");

                let fref_file = File::open(fref_path)
                    .expect("Failed to open feature definition file");

                let gene_index_file = File::open(&args.arg_gene_index.unwrap())
                    .expect("Failed to open gene index file");

                let fdist_file = File::open(fdist_path)
                    .expect("Failed to open feature barcode count file");

                Some(features::FeatureChecker::new(fref_file, gene_index_file,
                                                   fdist_file, library_type))
            } else {
                // TODO: Not true anymore
                // Gene Expression library type; don't need feature ref
                None
            }
        },
        _ => None,
    };

    println!("Setting up BAMs");
    let in_bam = bam::Reader::from_path(Path::new(&args.arg_in_bam.unwrap())).unwrap();
    let mut out_header = Header::from_template(in_bam.header());

    // Append PG tag
    let mut pg_rec = HeaderRecord::new(b"PG");
    pg_rec.push_tag(b"ID", &"annotate_reads");
    out_header.push_record(&pg_rec);

    // Append CO header lines
    match &args.flag_bam_comments {
        &Some(ref path_str) => {
            let file = File::open(&path_str).expect("Failed to open BAM comments JSON file");
            let reader = BufReader::new(file);
            let comments: Vec<String> = serde_json::from_reader(reader)
                .expect("Failed to parse BAM comments JSON file");
            for comment in comments {
                out_header.push_comment(comment.as_bytes());
            }
        },
        &None => {}
    }

    // Load library info
    let library_json_file = File::open(args.arg_library_info.unwrap())
        .expect("Failed to open library info JSON");
    let library_info: Vec<LibraryInfo> = serde_json::from_reader(library_json_file)
        .expect("Failed to read library info JSON from file");

    // Append CO header lines describing the libraries
    for lib_info in &library_info {
        let lib_str = serde_json::to_string(lib_info)
            .expect(&format!("Failed to serialize library_info to JSON: {:?}", lib_info));
        out_header.push_comment(format!("library_info:{}", lib_str).as_bytes());
    }

    let mut out_bam = bam::Writer::from_path(Path::new(&args.arg_out_bam.unwrap()), &out_header).unwrap();

    let chroms = in_bam.header().target_names().iter().map(|s| str::from_utf8(s).unwrap().to_string()).collect();
    let genomes = utils::get_reference_genomes(&reference_path);

    let mut reporter = Reporter::new(chroms, genomes, &annotator.get_transcript_index(), &annotator.get_params());

    println!("Processing reads");

    let mut num_alignments = 0;

    let qname_iter = in_bam.records()
        .map(|res| res.unwrap())
        .group_by(|rec| str::from_utf8(rec.qname()).unwrap().to_string());

    let tags_file = utils::open_maybe_compressed(args.arg_in_tags.unwrap());
    let tags_iter = bio::io::fastq::Reader::new(BufReader::new(tags_file)).records();

    let library_id = args.arg_library_id.unwrap();
    let library_index = library_info.iter().position(|lib| lib.library_id == library_id)
        .expect(&format!("Could not find library id {} in library info JSON", library_id));

    for ((_qname, genome_alignments), tags_rec) in qname_iter.zip(tags_iter) {
        num_alignments += genome_alignments.len();
        process_qname(tags_rec.unwrap().id(),
                      genome_alignments, &mut reporter,
                      &annotator, &bc_umi_checker, &feature_checker,
                      &mut out_bam, &gem_group, library_index);
    }

    println!("Writing metrics");
    reporter.get_metrics().write_binary(&args.arg_out_metrics.unwrap());

    println!("Writing metadata");
    let mut meta: HashMap<String, Value> = HashMap::new();
    meta.insert("num_alignments".into(), json!(num_alignments));
    utils::write_json(json!(meta), &args.arg_out_metadata.unwrap());
}

fn annotate_reads_join(args: Args) {
    println!("Merging metrics");

    let chunked_metrics = utils::read_json(&args.arg_in_chunked_metrics.unwrap());

    // Group metrics files by library type
    let mut metric_groups: HashMap<String, Vec<String>> = HashMap::new();
    for chunk in chunked_metrics.as_array().unwrap() {
        let filename = chunk.get("metrics").unwrap().as_str().unwrap();
        let library_type = chunk.get("library_type").unwrap().as_str().unwrap();
        let filenames = metric_groups.entry(library_type.to_owned()).or_insert(vec![]);
        filenames.push(filename.to_owned());
    }

    // Merge metrics by library type
    let mut summary = BTreeMap::new();
    for (lib_type, filenames) in &metric_groups {
        let mut merged_metrics = Metrics::new();

        for filename in filenames {
            merged_metrics.merge(&Metrics::read_binary(&filename));
        }

        // Prefix the metrics keys
        for (k, v) in merged_metrics.report() {
            let prefix = report::get_library_type_metric_prefix(lib_type);
            summary.insert(format!("{}{}", &prefix, k), v);
        }
    }

    let writer = File::create(&args.arg_out_json.unwrap()).unwrap();
    serde_json::to_writer_pretty(writer, &summary).expect("Failed to write JSON");

    //merged_metrics.write_barcodes_csv(&args.arg_out_bc_csv.unwrap());


}

/// Use transcriptome alignments to promote a single genome alignment
/// when none are confidently mapped to the genome.
/// Returns true if rescue took place.
fn rescue_alignments(r1_data: &mut Vec<RecordData>, mut maybe_r2_data: Option<&mut Vec<RecordData>>) -> bool {
    // Check if rescue is appropriate and determine which record to promote
    let mut seen_genes = HashSet::new();
    let mut promote_index: Option<usize> = None;

    for i in 0..r1_data.len() {
        let ref r1 = r1_data[i];
        let ref r2 = maybe_r2_data.as_ref().map(|ref x| &x[i]);

        // Abort if any of the records mapped uniquely to the genome
        if r1.rec.mapq() >= report::HIGH_CONF_MAPQ ||
            r2.map_or(0, |ref x| x.rec.mapq()) >= report::HIGH_CONF_MAPQ {
            return false;
        }

        let genes = report::get_alignment_gene_intersection(&r1.anno, r2.map(|x| &x.anno));

        // Track which record/record-pair we should promote;
        // Take the first record/pair with 1 gene
        if genes.len() == 1 {
            promote_index = promote_index.or(Some(i));
        }

        // Track number of distinct genes we're aligned to
        seen_genes.extend(genes.into_iter());
    }
    // There are >1 candidate genes to align to
    // or there are no candidate genes to align to.
    if seen_genes.len() > 1 || promote_index.is_none() {
        return false;
    }

    // Promote a single alignment
    for i in 0..r1_data.len() {
        let ref mut r1 = r1_data[i];

        if promote_index.unwrap() == i {
            // Promote one alignment
            r1.anno.rescued = true;
            r1.rec.set_mapq(report::HIGH_CONF_MAPQ);
            utils::set_primary(&mut r1.rec);
            if let Some(ref mut r2_data) = maybe_r2_data.as_mut() {
                r2_data[i].anno.rescued = true;
                r2_data[i].rec.set_mapq(report::HIGH_CONF_MAPQ);
                utils::set_primary(&mut r2_data[i].rec);
            }
        } else {
            // Demote the rest
            r1.rec.set_mapq(0);
            r1.rec.set_secondary();
            if let Some(ref mut r2_data) = maybe_r2_data.as_mut() {
                r2_data[i].rec.set_mapq(0);
                r2_data[i].rec.set_secondary();
            }
        }
    }
    true
}

fn process_qname(tag_string: &str,
                 genome_alignments: Vec<Record>,
                 reporter: &mut Reporter,
                 annotator: &TranscriptAnnotator,
                 bc_umi_checker: &BarcodeUmiChecker,
                 feature_checker: &Option<FeatureChecker>,
                 out_bam: &mut bam::Writer,
                 gem_group: &u8,
                 library_index: usize) {
    let fastq_header = CellrangerFastqHeader::new(tag_string);

    let bc_umi_data = bc_umi_checker.process_barcodes_and_umis(&fastq_header.tags);

    let feature_data = feature_checker.as_ref()
        .and_then(|fc| fc.process_feature_data(&fastq_header.tags));

    let stripped_qname = fastq_header.qname.into_bytes();

    let mut r1_data = Vec::new();
    let mut r2_data = Vec::new();

    // Annotate individual alignments
    for genome_alignment in genome_alignments {
        let anno = annotator.annotate_genomic_alignment(&genome_alignment);
        let rec_data = RecordData { rec: genome_alignment, anno: anno };
        if !rec_data.rec.is_last_in_template() {
            r1_data.push(rec_data);
        } else {
            r2_data.push(rec_data);
        }
    }

    // Annotate pairs
    let pair_data = match r2_data.len() > 0 {
        false => Vec::new(),
        true => r1_data.iter().zip(r2_data.iter())
            .map(|(r1, r2)| PairAnnotationData::from_pair(&r1.anno, &r2.anno)).collect(),
    };

    let mut read_data = ReadData {
        r1_data:        r1_data,
        r2_data:        r2_data,
        bc_umi_data:    bc_umi_data,
        feature_data:   feature_data,
        pair_data:      pair_data,
    };

    // STAR does not generate partially-mapped pairs. Double-check.
    let is_paired = read_data.is_paired_end();
    assert!(!is_paired || read_data.is_properly_paired());

    rescue_alignments(&mut read_data.r1_data, match is_paired {
        true => Some(&mut read_data.r2_data),
        false => None,
    });

    let is_conf_mapped = read_data.is_conf_mapped_to_transcriptome();
    let is_gene_discordant = read_data.is_gene_discordant();

    // Interleave pairs of aligned records (R1,R2; R1,R2; ...)
    for i in 0 .. read_data.r1_data.len() {
        {
            let pair_data = match read_data.is_paired_end() {
                true => Some(&read_data.pair_data[i]),
                false => None,
            };

            let ref mut r1 = read_data.r1_data[i];
            process_record(&stripped_qname, r1, pair_data,
                           &read_data.bc_umi_data, &read_data.feature_data, is_conf_mapped,
                           is_gene_discordant, gem_group, library_index);
            out_bam.write(&r1.rec).unwrap();
        }
        if read_data.is_paired_end() {
            let ref mut r2 = read_data.r2_data[i];
            process_record(&stripped_qname, r2, Some(&read_data.pair_data[i]),
                           &read_data.bc_umi_data, &read_data.feature_data, is_conf_mapped,
                           is_gene_discordant, gem_group, library_index);
            out_bam.write(&r2.rec).unwrap();
        }
    }

    reporter.update_metrics(&read_data);
}

/// Attach BAM tags to a single record.
/// is_conf_mapped: the qname is confidently mapped to the transcriptome
/// is_gene_discordant: the qname is gene-discordant
fn process_record(qname: &[u8],
                  record_data: &mut RecordData,
                  pair_anno: Option<&PairAnnotationData>,
                  bc_umi_data: &BarcodeUmiData,
                  feature_data: &Option<FeatureData>,
                  is_conf_mapped: bool,
                  is_gene_discordant: bool,
                  gem_group: &u8,
                  library_index: usize) {
    // Strip auxiliary tags from the qname
    record_data.rec.set_qname(qname);

    // Attach annotation tags
    record_data.anno.attach_tags(&mut record_data.rec, pair_anno);

    // Attach extra flags
    if !record_data.rec.is_secondary() {
        transcriptome::add_extra_flags(&mut record_data.rec,
                                       is_conf_mapped,
                                       is_gene_discordant,
                                       feature_data);
    }

    // Attach library id
    record_data.rec.push_aux(&LIBRARY_INDEX_TAG.as_bytes(),
                             &rust_htslib::bam::record::Aux::Integer(library_index as i32));

    // Attach cell barcode and UMI
    bc_umi_data.attach_tags(&mut record_data.rec, gem_group);

    // Attach feature tags
    if let &Some(ref feature_data) = feature_data {
        feature_data.attach_tags(&mut record_data.rec);
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::path::PathBuf;

    // TODO validate outputs?

    fn test_base(reference: &str, sample: &str) {
        // Setup test output dir
        let out_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("test_output").join("annotate_reads").join(&reference).join(&sample);
        fs::create_dir_all(&out_dir)
            .expect("Failed to create test output dir");

        let args_main = Args {
            cmd_main:               true,
            cmd_join:               false,
            arg_in_bam:             Some(format!("test/{}/{}/align_reads.bam", reference, sample).into()),
            arg_in_tags:            Some(format!("test/{}/{}/tags.fastq.lz4", reference, sample).into()),
            arg_out_bam:            Some(format!("{}/annotate_reads.bam", out_dir.to_str().unwrap()).into()),
            arg_out_metrics:        Some(format!("{}/metrics.bincode", out_dir.to_str().unwrap()).into()),
            arg_reference_path:     Some(format!("test/{}", reference).into()),
            arg_gene_index:         Some(format!("test/{}/gene_index.tab", reference).into()),
            arg_bc_counts:          Some(format!("test/{}/{}/barcode_counts.json", reference, sample).into()),
            arg_bc_whitelist:       Some("../../python/cellranger/barcodes/737K-august-2016.txt".into()),
            arg_gem_group:          Some(1),
            arg_out_metadata:       Some(format!("{}/metadata.json", out_dir.to_str().unwrap()).into()),
            arg_strandedness:       Some("+".into()),
            arg_in_chunked_metrics: None,
            arg_out_json:           None,
            arg_out_bc_csv:         None,
            arg_feature_dist:       None,
            arg_library_type:       Some("Gene Expression".to_owned()),
            arg_library_id:         Some("0".to_owned()),
            arg_library_info:       Some("test/library_info.json".to_owned()),
            flag_fiveprime:         false,
            flag_bam_comments:      None,
            flag_feature_ref:       None,
        };
        annotate_reads_main(args_main);

        let args_join = Args {
            cmd_main:               false,
            cmd_join:               true,
            arg_in_bam:             None,
            arg_in_tags:            None,
            arg_out_bam:            None,
            arg_out_metrics:        None,
            arg_reference_path:     None,
            arg_gene_index:         None,
            arg_bc_counts:          None,
            arg_bc_whitelist:       None,
            arg_gem_group:          None,
            arg_out_metadata:       None,
            arg_strandedness:       None,
            arg_in_chunked_metrics: Some(format!("test/{}/{}/chunked_metrics.json", reference, sample).into()),
            arg_out_json:           Some(format!("{}/summary.json", out_dir.to_str().unwrap()).into()),
            arg_out_bc_csv:         Some(format!("{}/barcodes_detected.csv", out_dir.to_str().unwrap()).into()),
            arg_feature_dist:       None,
            arg_library_type:       Some("Gene Expression".to_owned()),
            arg_library_id:         Some("0".to_owned()),
            arg_library_info:       Some("test/library_info.json".to_owned()),
            flag_fiveprime:         false,
            flag_bam_comments:      None,
            flag_feature_ref:       None,
        };
        annotate_reads_join(args_join);
    }

    #[test]
    fn test_single_end() {
        test_base("GRCh38", "30647_3p")
    }

    #[test]
    fn test_paired_end() {
        test_base("GRCh38", "37358_5p_pe")
    }

    #[test]
    fn test_barnyard() {
        test_base("hg19_and_mm10", "30570_3p")
    }

    #[test]
    fn test_edge_cases() {
        test_base("GRCh38", "edge_cases")
    }
}
