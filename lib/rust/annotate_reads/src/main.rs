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

mod barcodes;
mod metrics;
mod reference;
mod report;
mod transcriptome;
mod utils;

use std::path::Path;
use std::fs::File;
use std::io::{BufReader};
use std::str;
use std::collections::{HashSet, HashMap};

use docopt::Docopt;

use rust_htslib::bam;
use rust_htslib::bam::Read;
use rust_htslib::bam::header::{Header, HeaderRecord};
use rust_htslib::bam::record::Record;

use serde_json::Value;

use itertools::Itertools;

use transcriptome::{TranscriptAnnotator, AnnotationParams, Strandedness};
use report::{Reporter, Metrics, ReadData, RecordData};
use barcodes::{BarcodeUmiChecker, BarcodeUmiData};
use metrics::MetricGroup;
use utils::CellrangerFastqHeader;

const USAGE: &'static str = "
Usage:
  annotate_reads main <in-bam> <out-bam> <out-metrics> <reference-path> <gene-index> <bc-counts> <bc-whitelist> <gem-group> <out-metadata> <strandedness> [--fiveprime] [--bam-comments=F]
  annotate_reads join <in-metrics-list> <out-json> <out-bc-csv>
  annotate_reads (-h | --help)

Options:
  -h --help            Show this screen.
  --neg-strand         Assume negative-strand chemistry.
  --fiveprime          Assume reads originate from 5' end (instead of 3').
  --bam-comments=F     JSON file with array of strings to append as @CO items
";

#[derive(Debug, Deserialize, Clone)]
struct Args {
    cmd_main:               bool,
    cmd_join:               bool,

    // main args
    arg_in_bam:             Option<String>,
    arg_out_bam:            Option<String>,
    arg_out_metrics:        Option<String>,
    arg_reference_path:     Option<String>,
    arg_gene_index:         Option<String>,
    arg_bc_counts:          Option<String>,
    arg_bc_whitelist:       Option<String>,
    arg_gem_group:          Option<u8>,
    arg_out_metadata:       Option<String>,
    arg_strandedness:       Option<String>,

    // main flags
    flag_fiveprime:         bool,
    flag_bam_comments:      Option<String>,

    // join args
    arg_in_metrics_list:    Option<String>,
    arg_out_json:           Option<String>,
    arg_out_bc_csv:         Option<String>,
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
    let annotator = TranscriptAnnotator::new(&reference_path, &args.arg_gene_index.unwrap(), params);

    println!("Loading whitelist");
    let gem_group = &args.arg_gem_group.unwrap();
    let bc_umi_checker = barcodes::BarcodeUmiChecker::new(&args.arg_bc_counts.unwrap(), &args.arg_bc_whitelist.unwrap(), gem_group);

    println!("Setting up BAMs");
    let in_bam = bam::Reader::from_path(Path::new(&args.arg_in_bam.unwrap())).unwrap();
    let mut out_header = Header::from_template(in_bam.header());

    // Append PG tag
    let mut pg_rec = HeaderRecord::new(b"PG");
    pg_rec.push_tag(b"ID", &"annotate_reads");
    out_header.push_record(&pg_rec);

    // Append CO tags
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

    let mut out_bam = bam::Writer::from_path(Path::new(&args.arg_out_bam.unwrap()), &out_header).unwrap();

    let chroms = in_bam.header().target_names().iter().map(|s| str::from_utf8(s).unwrap().to_string()).collect();
    let genomes = utils::get_reference_genomes(&reference_path);

    let mut reporter = Reporter::new(chroms, genomes, &annotator.get_transcript_index(), &annotator.get_params());

    println!("Processing reads");
    // TODO untangle this mess and handle errors
    let mut num_alignments = 0;
    for (qname, genome_alignments) in in_bam.records().map(|res| res.unwrap()).group_by(|rec| str::from_utf8(rec.qname()).unwrap().to_string()) {
        num_alignments += genome_alignments.len();
        process_qname(qname, genome_alignments, &mut reporter, &annotator, &bc_umi_checker, &mut out_bam, &gem_group);
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
    let chunked_metrics = utils::load_txt(&args.arg_in_metrics_list.unwrap());
    let mut merged_metrics = Metrics::new();
    for metrics in chunked_metrics {
        merged_metrics.merge(&Metrics::read_binary(&metrics));
    }
    merged_metrics.write_json_summary(&args.arg_out_json.unwrap());
    merged_metrics.write_barcodes_csv(&args.arg_out_bc_csv.unwrap());
}

fn process_qname(qname: String, genome_alignments: Vec<Record>, reporter: &mut Reporter,
                 annotator: &TranscriptAnnotator, bc_umi_checker: &BarcodeUmiChecker, out_bam: &mut bam::Writer,
                 gem_group: &u8) {
    let fastq_header = CellrangerFastqHeader::new(qname);
    let stripped_qname = fastq_header.qname.into_bytes();
    let bc_umi_data = bc_umi_checker.process_barcodes_and_umis(fastq_header.tags);

    let mut r1_data = Vec::new();
    let mut r2_data = Vec::new();

    let mut seen_genes = HashSet::new();
    let mut seen_primary = false;
    let mut do_rescue = false;

    for genome_alignment in genome_alignments {
        let mut anno = annotator.annotate_genomic_alignment(&genome_alignment);
        for gene in &anno.genes { seen_genes.insert(gene.to_owned()); }

        if genome_alignment.mapq() == 255 {
            seen_primary = true;
        }
        if (!seen_primary) && (!do_rescue) && (anno.genes.len() > 0) {
            do_rescue = true;
            anno.rescued = true;
        }

        let rec_data = RecordData { rec: genome_alignment, anno: anno };
        if !rec_data.rec.is_last_in_template() {
            r1_data.push(rec_data);
        } else {
            r2_data.push(rec_data);
        }
    }

    do_rescue = (seen_genes.len() == 1) && !seen_primary && do_rescue;

    let mut read_data = ReadData {
        r1_data:        r1_data,
        r2_data:        r2_data,
        bc_umi_data:    bc_umi_data,
        num_genes:      seen_genes.len() as i32,
    };

    if read_data.is_properly_paired() {
        // interleave pairs + rescue together
        for i in 0..read_data.r1_data.len() {
            let ref mut r1_data = read_data.r1_data[i];
            let ref mut r2_data = read_data.r2_data[i];
            if do_rescue && (r1_data.anno.rescued || r2_data.anno.rescued) {
                r1_data.anno.rescued = true;
                r2_data.anno.rescued = true;
            } else {
                r1_data.anno.rescued = false;
                r2_data.anno.rescued = false;
            }
            process_record(&stripped_qname, r1_data, &read_data.bc_umi_data, do_rescue, gem_group);
            process_record(&stripped_qname, r2_data, &read_data.bc_umi_data, do_rescue, gem_group);
            out_bam.write(&r1_data.rec).unwrap();
            out_bam.write(&r2_data.rec).unwrap();
        }
    } else {
        for i in 0..read_data.r1_data.len() {
            let ref mut r1_data = read_data.r1_data[i];
            if !do_rescue { r1_data.anno.rescued = false; }
            process_record(&stripped_qname, r1_data, &read_data.bc_umi_data, do_rescue, gem_group);
            out_bam.write(&r1_data.rec).unwrap();
        }
        for i in 0..read_data.r2_data.len() {
            let ref mut r2_data = read_data.r2_data[i];
            if !do_rescue { r2_data.anno.rescued = false; }
            process_record(&stripped_qname, r2_data, &read_data.bc_umi_data, do_rescue, gem_group);
            out_bam.write(&r2_data.rec).unwrap();
        }
    }

    reporter.update_metrics(&read_data);
}

fn process_record(qname: &[u8], record_data: &mut RecordData, bc_umi_data: &BarcodeUmiData,
                  do_rescue: bool, gem_group: &u8) {
    // Strip auxiliary tags from the qname
    record_data.rec.set_qname(qname);

    // Set mapq and flags
    if do_rescue && !record_data.rec.is_unmapped() {
        if record_data.anno.rescued {
            record_data.rec.set_mapq(255);
            utils::set_primary(&mut record_data.rec);
        } else {
            record_data.rec.set_mapq(0);
            record_data.rec.set_secondary();
        }
    }

    // Attach annotation tags
    record_data.anno.attach_tags(&mut record_data.rec);

    // Attach BC and UMI
    bc_umi_data.attach_tags(&mut record_data.rec, gem_group);

}


#[cfg(test)]
mod tests {
    use super::*;

    // TODO validate outputs?

    fn test_base(reference: &str, sample: &str) {
        let args_main = Args {
            cmd_main:               true,
            cmd_join:               false,
            arg_in_bam:             Some(format!("test/{}/{}/align_reads.bam", reference, sample).into()),
            arg_out_bam:            Some(format!("test/{}/{}/annotate_reads.bam", reference, sample).into()),
            arg_out_metrics:        Some(format!("test/{}/{}/metrics.bincode", reference, sample).into()),
            arg_reference_path:     Some(format!("test/{}", reference).into()),
            arg_gene_index:         Some(format!("test/{}/gene_index.tab", reference).into()),
            arg_bc_counts:          Some(format!("test/{}/{}/barcode_counts.json", reference, sample).into()),
            arg_bc_whitelist:       Some("../../python/cellranger/barcodes/737K-august-2016.txt".into()),
            arg_gem_group:          Some(1),
            arg_out_metadata:       Some(format!("test/{}/{}/metadata.json", reference, sample).into()),
            arg_strandedness:       Some("+".into()),
            arg_in_metrics_list:    None,
            arg_out_json:           None,
            arg_out_bc_csv:         None,
            flag_fiveprime:         false,
        };
        annotate_reads_main(args_main);

        let args_join = Args {
            cmd_main:               false,
            cmd_join:               true,
            arg_in_bam:             None,
            arg_out_bam:            None,
            arg_out_metrics:        None,
            arg_reference_path:     None,
            arg_gene_index:         None,
            arg_bc_counts:          None,
            arg_bc_whitelist:       None,
            arg_gem_group:          None,
            arg_out_metadata:       None,
            arg_strandedness:       None,
            arg_in_metrics_list:    Some(format!("test/{}/{}/chunked_metrics.txt", reference, sample).into()),
            arg_out_json:           Some(format!("test/{}/{}/summary.json", reference, sample).into()),
            arg_out_bc_csv:         Some(format!("test/{}/{}/barcodes_detected.csv", reference, sample).into()),
            flag_fiveprime:         false,
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
