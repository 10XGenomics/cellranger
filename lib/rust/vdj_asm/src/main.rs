//
// Copyright (c) 2017 10x Genomics, Inc. All rights reserved.
//

#![crate_name = "vdj_asm"]
#![allow(dead_code)]

extern crate bio;
extern crate tada;
extern crate flate2;
extern crate itertools;
extern crate docopt;
extern crate rustc_serialize;
extern crate time;
extern crate num;
extern crate rand;
extern crate rust_htslib;

pub mod graph;
mod constants;
mod fastq;
mod utils;
mod sw;
mod asm;
mod bam_utils;
mod graph_read;
mod metrics;

use rust_htslib::bam;
use rust_htslib::bam::Read;

use tada::bitenc;
use std::path::{Path, PathBuf};
use std::io::{BufWriter, Write};
use std::collections::{HashSet, HashMap};
use docopt::Docopt;
use constants::{READ_LEN, UmiType, PROCESSED_UMI_TAG,
                SEED_LEN, MATCH_SCORE, MAX_NUM_READPAIRS,
                MISMATCH_SCORE, GAP_OPEN, GAP_EXTEND, CLIP, DUMMY_CONTIG_NAME,
                QUAL_OFFSET};
use std::fs;
use std::fs::File;

use time::PreciseTime;

const USAGE: &'static str = "
Usage:
vdj_asm asm <inbam> <outdir> [--plot][--kmers=<NUM>][--min-contig=<NUM>][--npaths=<NUN>][--frac-reads=<NUM>][--score-factor=<NUM>][--qual-factor=<NUM>][--rt-error=<NUM>][--min-qual=<NUM>][--nx=<NUM>][--match-score=<NUM>][--miss-score=<NUM>][--gap-open=<NUM>][--gap-extend=<NUM>][--min-sw-score=<NUM>][--min-umi-reads=<NUM>][--cons][--single-end][--subsample-rate=<NUM>][--fast-align]
vdj_asm base-quals <inpref> <outdir> [--single-end][--rev-strand][--rt-error=<NUM>][--global][--match-score=<NUM>][--miss-score=<NUM>][--gap-open=<NUM>][--gap-extend=<NUM>][--seed=<NUM>]
vdj_asm read-match <fasta> <fqpref> <outbam> [--rev-strand][--match-score=<NUM>][--miss-score=<NUM>][--gap-open=<NUM>][--gap-extend=<NUM>][--seed=<NUM>][--min-sw-score=<NUM>]
vdj_asm (-h | --help)

Options:
   --plot                 Create dot file with assembly graph.
   --kmers=<NUM>          Minimum number of occurrences to consider a kmer
   --min-contig=<NUM>     Minimum output contig length
   --npaths=<NUM>         Number of paths per component.
   --rev_strand           Data are RF (default is FR).
   --min-sw-score=<NUM>   Minimum SW score to consider an alignment.
   --rt-error=<NUM>       RT error rate.
   --trace=<STR>
   --nx=<NUM>
   --qual-factor=<NUM>
   --score-factor=<NUM>
   --frac-reads=<NUM>
   --min-umi-reads=<NUM>
   --read2
   --fast-align
   --match-score=<NUM>
   --miss-score=<NUM>
   --gap-open=<NUM>
   --gap-extend=<NUM>
   --seed=<NUM>
   --subsample-rate=<NUM>
   --cons                 Compute single consensus sequence
   --single-end           Input data are single end
   --global               Try global alignemnt if local alignment doesn't return a good hit
";

#[derive(Debug, RustcDecodable, Clone)]
pub struct Args {
    cmd_asm: bool,
    cmd_base_quals: bool,
    cmd_read_match: bool,

    arg_inbam: Option<String>,
    arg_inpref: Option<String>,
    arg_outdir: Option<String>,
    arg_fasta: Option<String>,
    arg_fqpref: Option<String>,
    arg_outbam: Option<String>,

    flag_cons: bool,
    flag_kmers: Option<u32>,
    flag_min_contig: Option<usize>,
    flag_npaths: Option<usize>,
    flag_rev_strand: bool,
    flag_plot: bool,
    flag_min_sw_score: Option<f64>,
    flag_rt_error: Option<f64>,
    flag_min_qual: Option<u8>,
    flag_nx: Option<f64>,
    flag_qual_factor: Option<f64>,
    flag_score_factor: Option<f64>,
    flag_frac_reads: Option<f64>,
    flag_min_umi_reads: Option<usize>,
    flag_subsample_rate: Option<f64>,
    flag_fast_align: bool,
    flag_match_score: Option<f32>,
    flag_miss_score: Option<f32>,
    flag_gap_open: Option<f32>,
    flag_gap_extend: Option<f32>,
    flag_clip: Option<f32>,
    flag_seed: Option<usize>,
    flag_single_end: bool,
    flag_global: bool,
}

fn main() {

    let mut args: Args = Docopt::new(USAGE).and_then(|d| d.decode()).unwrap_or_else(|e| e.exit());

    if args.cmd_asm {
        if args.flag_cons {
            args.flag_min_umi_reads = Some(1);
            args.flag_min_contig = Some(1);
        }
        vdj_asm(args, MAX_NUM_READPAIRS);
    } else if args.cmd_base_quals {
        get_base_qualities(args);
    } else {
        get_matches(args);
    }
}

/// Number of reads we're currently storing
/// reads: reads we're going to use
/// filtered_pairs: didn't map to any V(D)J segment sequence
/// dropped_pairs: dropped due to random subsampling
fn reads_stored(reads: &Vec<graph_read::Read>,
                filtered_pairs: &Vec<(graph_read::Read, Option<graph_read::Read>)>,
                dropped_pairs: &Vec<(graph_read::Read, Option<graph_read::Read>)>,
                single_end: bool) -> usize {
    reads.len() + (!single_end as usize + 1) * (filtered_pairs.len() + dropped_pairs.len())
}

/// max_readpairs_per_bc: Maximum number of pairs (if paired-end) or reads (if single-end)
/// per barcode to read before we start assembling. All read pairs (or reads) after the first
/// that many will be skipped. This max limit includes filtered pairs (reads) and dropped pairs (reads).
pub fn vdj_asm(args: Args, max_readpairs_per_bc: usize) {

    let aligner_params = get_align_params(&args);
    let input_bam_filename = args.arg_inbam.unwrap();
    let out_dirname = args.arg_outdir.unwrap();

    let min_kmer_count = match args.flag_kmers {
        Some(c) => c,
        None => 2,
    };
    let min_contig_len = match args.flag_min_contig {
        Some(c) => c,
        None => READ_LEN,
    };
    let paths_per_component = match args.flag_npaths {
        Some(c) => c,
        None => 3,
    };
    let min_align_score = match args.flag_min_sw_score {
        Some(c) => c,
        None => 50.0,
    };
    let rt_error = match args.flag_rt_error {
        Some(c) => c,
        None => 0.0001,
    };
    let min_qual = match args.flag_min_qual {
        Some(c) => c,
        None => 10,
    };
    let nx = match args.flag_nx {
        Some(c) => c,
        None => 99.9,
    };
    let qual_factor = match args.flag_qual_factor {
        Some(c) => c,
        None => 0.5,
    };
    let score_factor = match args.flag_score_factor {
        Some(c) => c,
        None => 0.9,
    };
    let frac_path_reads = match args.flag_frac_reads {
        Some(c) => c,
        None => 0.3,
    };
    let min_umi_reads = args.flag_min_umi_reads;
    let subsample_rate = match args.flag_subsample_rate {
        Some(c) => c,
        None => 1.0,
    };

    let plot = args.flag_plot;

    let fast_align = args.flag_fast_align;
    let single_end = args.flag_single_end;
    let max_reads_per_bc = match single_end {
        true => max_readpairs_per_bc,
        false => 2 * max_readpairs_per_bc,
    };

    let out_pref = Path::new(&out_dirname).join(Path::new(&(input_bam_filename)).file_stem().unwrap())
                                          .into_os_string().into_string().unwrap();

    let fasta_file = match File::create(out_pref.to_string() + ".fasta") {
        Ok(file) => file,
        Err(..) => panic!("Could not create fasta file"),
    };
    let fasta_writer = BufWriter::new(&fasta_file);

    let fastq_file = match File::create(out_pref.to_string() + ".fastq") {
        Ok(file) => file,
        Err(..) => panic!("Could not create fastq file"),
    };
    let fastq_writer = BufWriter::new(&fastq_file);

    let summary_file = match File::create(out_pref.to_string() + "_summary.tsv") {
        Ok(file) => file,
        Err(..) => panic!("Could not create summary file"),
    };
    let mut summary_writer = BufWriter::new(&summary_file);
    summary_writer.write_fmt(format_args!("{}\t{}\t{}\t{}\t{}\t{}\n",
                                          "barcode", "contig_name",
                                          "num_reads", "num_pairs",
                                          "num_umis", "umi_list")).unwrap();

    let umi_summary_file = match File::create(out_pref.to_string() + "_umi_summary.tsv") {
        Ok(file) => file,
        Err(..) => panic!("Could not create UMI summary file"),
    };
    let mut umi_summary_writer = BufWriter::new(&umi_summary_file);
    umi_summary_writer.write_fmt(format_args!("{}\t{}\t{}\t{}\t{}\t{}\t{}\n", "barcode", "umi_id", "umi",
                                              "reads", "min_umi_reads", "good_umi", "contigs")).unwrap();

    let mut assembly_outs = asm::AssemblyOuts::new(fasta_writer, fastq_writer,
                                                   summary_writer, umi_summary_writer);

    let in_bam = bam::Reader::from_path(&Path::new(&input_bam_filename)).ok().expect("Error opening input bam.");

    // Will create one bam per barcode and merge incrementally
    // Can't append to a single bam, because header (contig names) is not known in advance.
    let bam_merge_frequency = 50; // merge every N BAMs / barcodes to reduce number of files.
    let mut out_bams = Vec::new();
    //let all_contig_names = Vec::new();

    let mut last_barcode = "".to_string();
    let mut reads : Vec<graph_read::Read> = Vec::new();
    let mut filtered_pairs : Vec<(graph_read::Read, Option<graph_read::Read>)> = Vec::new();
    let mut dropped_pairs : Vec<(graph_read::Read, Option<graph_read::Read>)> = Vec::new();

    let mut umi_map = HashMap::new(); // UMI string -> UMI id
    let mut umi_counts = HashMap::new();

    // Create a dummy Writer so the compiler doesn't complain that this is uninitialized.
    let (header, _) = get_assembly_contig_header(&Vec::new(), &last_barcode);
    let mut tmp_bam = bam::Writer::from_path(&Path::new(&format!("{}_{}.bam", out_pref, "dummy")), &header).unwrap(); // This should never be used.

    let mut bam_iter = in_bam.records();
    let mut reading_through = false;
    let mut npairs = 0;

    let mut metrics: metrics::AssemblyMetrics = Default::default();

    loop {
        let wrapped_record = bam_iter.next();
        let mut barcode = "".to_string();
        let mut record = bam::record::Record::new();

        let mut done_reading = false;
        match wrapped_record {
            Some(r) => {
                record = r.ok().expect("");
                // We haven't reached the end of the records yet.
                let header = fastq::CellrangerFastqHeader::new(String::from_utf8_lossy(record.qname()).into_owned());
                let wrapped_barcode = header.get_barcode();
                if wrapped_barcode.is_none() {
                    continue;
                }
                barcode = wrapped_barcode.unwrap();
            },
            None => { done_reading = true; },
        };

        if !reads.is_empty() && !reading_through && (barcode != last_barcode ||
                                                     done_reading ||
                                                     reads_stored(&reads, &filtered_pairs, &dropped_pairs, single_end) >= max_reads_per_bc) {

            // if we've reached the BAM limit, squash the current set before proceeding
            if out_bams.len() >= bam_merge_frequency {
                drop(tmp_bam); // make sure last BAM is flushed
                let squashed_bam_filename = format!("{}_{}_{}.bam", out_pref, last_barcode, out_bams.len());
                bam_utils::concatenate_bams(&out_bams, &squashed_bam_filename);
                for bam in out_bams.iter() {
                    match fs::remove_file(bam) {
                       Ok(()) => (),
                       Err(e) => println!("Error removing temporary BAM file: {}. Likely a filesystem hiccup. Continuing", e),
                    }
                }
                out_bams.clear();
                out_bams.push(squashed_bam_filename);
            }

            if reads_stored(&reads, &filtered_pairs, &dropped_pairs, single_end) >= max_reads_per_bc {
                println!("Reached maximum number of reads for barcode {}. Will skip the rest of the reads.", last_barcode);
            }
            println!("Barcode {}: Assembling {:?} reads, {:?} distinct UMIs", last_barcode, reads.len(), umi_counts.len());

            // If min_umi_reads is provided, then use this as the reads/umi cutoff.
            // Otherwise, use the NX reads/UMI. Note that reads without valid UMIs will be included in this calculation.
            let nx_read_count = match min_umi_reads {
                Some(c) => c,
                None => utils::nx_count(&umi_counts.values().map(|x| *x as usize).collect(), nx),
            };
            println!("Reads/UMI cutoff: {}", nx_read_count);

            let good_umis : HashSet<UmiType> = umi_counts.iter().filter(|&(_, c)| *c >= nx_read_count).map(|(u, _)| *u).collect();
            println!("Good UMIs: {:?}", good_umis);

            // Only keep sequences with good UMIs
            let mut assembled_reads = Vec::new();
            for read in reads.iter() {
                if good_umis.contains(&(read.umi)) {
                    assembled_reads.push(read);
                }
            }

            println!("Reads to use in assembly: {}", assembled_reads.len());
            let assemblable_read_pairs = match single_end {
                true => assembled_reads.len() as u64,
                false => assembled_reads.len() as u64 / 2,
            };
            metrics.assemblable_read_pairs_by_bc.insert(last_barcode.clone(), assemblable_read_pairs);

            let cc_start = PreciseTime::now();
            let (graph, contigs) = asm::assemble_reads(&mut assembled_reads,
                                                       &umi_counts,
                                                       min_kmer_count as usize,
                                                       qual_factor,
                                                       paths_per_component,
                                                       rt_error, min_qual, min_contig_len,
                                                       frac_path_reads,
                                                       min_align_score,
                                                       score_factor,
                                                       fast_align,
                                                       args.flag_cons,
                                                       &aligner_params,
                                                       !single_end);

            println!("Assembled {} contigs in {} sec", contigs.len(), cc_start.to(PreciseTime::now()));

            if plot {
                graph.to_dot(PathBuf::from(out_pref.to_string() + "_full_graph.dot"), false);
            }

            let (header, contig_names) = get_assembly_contig_header(&contigs.iter().map(|x| x.0.clone()).collect(), &last_barcode);
            let out_bam_filename = format!("{}_{}.bam", out_pref, last_barcode);

            tmp_bam = bam::Writer::from_path(&Path::new(&out_bam_filename), &header).unwrap();
            out_bams.push(out_bam_filename);

            asm::write_assembly_results(&reads,
                                        &umi_map,
                                        &umi_counts,
                                        &good_umis,
                                        &mut assembly_outs,
                                        &mut tmp_bam,
                                        &contigs,
                                        &contig_names,
                                        &last_barcode,
                                        nx_read_count,
                                        single_end);

            // Write all dropped reads
            for &(ref dropped_read1, ref dropped_read2) in dropped_pairs.iter() {
                let _ = tmp_bam.write(&mut dropped_read1.to_unmapped_bam_record());
                if !single_end {
                    let _ = tmp_bam.write(&mut dropped_read2.clone().unwrap().to_unmapped_bam_record());
                }
            }

            // Write all filtered reads
            for &(ref filtered_read1, ref filtered_read2) in filtered_pairs.iter() {
                let _ = tmp_bam.write(&mut filtered_read1.to_unmapped_bam_record());
                if !single_end {
                    let _ = tmp_bam.write(&mut filtered_read2.clone().unwrap().to_unmapped_bam_record());
                }
            }

            reading_through = barcode == last_barcode;
        }

        if done_reading {
            break;
        }

        if barcode != last_barcode {
            reads.clear();
            filtered_pairs.clear();
            dropped_pairs.clear();
            umi_map.clear();
            umi_counts.clear();
            npairs = 0;
            reading_through = false;
            last_barcode = barcode;
        }

        let header = fastq::CellrangerFastqHeader::new(String::from_utf8_lossy(record.qname()).into_owned());
        let umi_str = header.get_tag(&(PROCESSED_UMI_TAG.to_string())).unwrap_or("".to_string());
        if !umi_map.contains_key(&umi_str.clone()) {
            let new_id = umi_map.len() as UmiType;
            umi_counts.insert(new_id, 0);
            umi_map.insert(umi_str.clone(), new_id);
        }
        let umi = umi_map.get(&umi_str.clone()).unwrap();

        if single_end {
            let mut read1 = graph_read::Read::from_bam_record(2 * npairs, *umi, &record);
            read1.unset_paired();
            if umi_str == "".to_string() {
                assert!(record.is_unmapped());
            }
            if !reading_through {
                if record.is_unmapped() {
                    filtered_pairs.push((read1, None));

                } else if utils::drop_read(subsample_rate, &read1.name) {
                    dropped_pairs.push((read1, None));

                } else {
                    reads.push(read1);
                    *umi_counts.get_mut(umi).unwrap() += 1; // read count
                }
            } else {
                let _ = tmp_bam.write(&mut read1.to_unmapped_bam_record());
            }

        } else {
            let mut read1 = graph_read::Read::from_bam_record(2 * npairs, *umi, &record);
            let mate_record = bam_iter.next().unwrap().ok().expect("Error while reading input bam");
            let mut read2 = graph_read::Read::from_bam_record(2 * npairs + 1, *umi, &mate_record);
            read1.set_paired();
            read2.set_paired();

            // Make sure that the only reads that are missing a UMI are reads that are
            // going to be ignored anyway.
            if umi_str == "".to_string() {
                assert!(record.is_unmapped() && mate_record.is_unmapped());
            }

            if !reading_through {
                if record.is_unmapped() && mate_record.is_unmapped() {
                    filtered_pairs.push((read1, Some(read2)));

                } else if utils::drop_read(subsample_rate, &read1.name) {
                    dropped_pairs.push((read1, Some(read2)));

                } else {
                    reads.push(read1);
                    reads.push(read2);
                    *umi_counts.get_mut(umi).unwrap() += 2; // read count
                }
            } else {
                // Just output the rest of the reads as unmapped.
                let _ = tmp_bam.write(&mut read1.to_unmapped_bam_record());
                let _ = tmp_bam.write(&mut read2.to_unmapped_bam_record());
            }
        }

        npairs += 1; // pairs (or reads if single-end)
    }

    // Force tmp_bam to get flushed
    drop(tmp_bam);

    // Write metrics summary
    {
        let mut metrics_summary_file = File::create(out_pref.to_string() + "_metrics_summary.json")
            .expect("Could not create metrics summary file");
        metrics::write_summary(&mut metrics_summary_file, &metrics)
            .expect("Could not write metrics summary");
    }

    assembly_outs.close();

    bam_utils::concatenate_bams(&out_bams, &(out_pref.to_string() + ".bam"));
    bam_utils::sort_and_index(&(out_pref.to_string() + ".bam"), &(out_pref.to_string() + "_sorted.bam"));
}

fn get_assembly_contig_header(seqs: &Vec<String>, pref: &String) -> (bam::header::Header, Vec<String>) {
    let mut header = bam::header::Header::new();
    let mut header_rec = bam::header::HeaderRecord::new(b"PG");
    let mut contig_names = Vec::with_capacity(seqs.len());
    header_rec.push_tag(b"ID", &"vdj_asm asm");
    header.push_record(&header_rec);

    if seqs.is_empty() {
        println!("isempty");
        bam_utils::add_ref_to_bam_header(&mut header, &DUMMY_CONTIG_NAME.to_string(), 0);
    } else {
        println!("notisempty");
        for (idx, ref seq) in seqs.iter().enumerate() {
            let contig_name = format!("{}_contig_{}", pref, idx + 1);
            contig_names.push(contig_name.clone());
            bam_utils::add_ref_to_bam_header(&mut header, &contig_name, seq.len());
        }
    }
    (header, contig_names)
}

/// Compute base qualities of a sequence based on a set of reads.
pub fn get_base_qualities(args: Args) {

    let aligner_params = get_align_params(&args);
    let file_pref = args.arg_inpref.unwrap();
    let out_dir = args.arg_outdir.unwrap();

    let fasta = bio::io::fasta::Reader::from_file(Path::new(&(file_pref.clone() + ".fasta"))).ok().expect("Could not open input FASTA");
    let mut record_iter = fasta.records();
    let record = record_iter.next().unwrap().ok().expect("Error reading sequence from fasta file");
    let ref_seq_name = record.id().unwrap().to_owned();
    let ref_seq = utils::replace_ns(&String::from_utf8_lossy(record.seq()).into_owned(), &ref_seq_name);

    let record = record_iter.next();
    if record.is_some() {
        println!("FASTA file has more than one sequence. Only the first one will be used.")
    }

    let rt_error = match args.flag_rt_error {
        Some(c) => c,
        None => 0.0001,
    };
    let seed_len = match args.flag_seed {
        Some(c) => c,
        None => SEED_LEN,
    };
    let single_end = args.flag_single_end;

    let mut fq_iter = match single_end {
        true => {
            let name1 = file_pref.clone() + ".fastq";
            let fq_path1 = Path::new(&name1);
            fastq::CellrangerPairedFastqIter::new(fq_path1, None, false)
        },
        false => {
            let name1 = file_pref.clone() + "_1.fastq";
            let name2 = file_pref.clone() + "_2.fastq";
            let fq_path1 = Path::new(&name1);
            let fq_path2 = Path::new(&name2);
            fastq::CellrangerPairedFastqIter::new(fq_path1, Some(fq_path2), false)
        },
    };

    let base_path = Path::new(&(file_pref));
    // components returns a Components, last returns a Component,
    // which we convert to OsStr, which we convert to str. phew!
    let contig_pref = base_path.components().last().unwrap().as_os_str().to_str().unwrap();
    println!("Using contig prefix {}", contig_pref);

    let mut out_path = PathBuf::from(out_dir);
    out_path.push(contig_pref);
    let out_pref = out_path.to_str().unwrap();

    let fastq_file = match File::create(out_pref.to_string() + ".fastq") {
        Ok(file) => file,
        Err(..) => panic!("Could not create fastq file {}", out_pref.to_string() + ".fastq"),
    };
    let mut fastq_writer = BufWriter::new(&fastq_file);

    let mut header = bam::header::Header::new();
    let mut header_rec = bam::header::HeaderRecord::new(b"PG");
    header_rec.push_tag(b"ID", &"vdj_asm base-quals");
    header.push_record(&header_rec);
    bam_utils::add_ref_to_bam_header(&mut header, &ref_seq_name, ref_seq.len());

    let mut out_bam = bam::Writer::from_path(&(out_pref.to_string() + ".bam"), &header).unwrap();

    let aligner = sw::Aligner::new(&vec![ref_seq.clone()], seed_len);

    let mut good_alignments = Vec::new();

    loop {
        match fq_iter.next() {
            Some(pair) => {
                let r1 = pair.r1;

                let bitenc1 = bitenc::BitEnc::from_bytes(&r1.seq);
                let mut read1 = graph_read::Read::new(r1.id, pair.umi, pair.header.clone(), bitenc1, r1.quals.clone());
                let alignment1 = local_or_global_align(&aligner, &read1, &aligner_params, args.flag_global);

                let (read2, alignment2) = match single_end {
                    true => (None, None),
                    false => {
                        let r2 = pair.r2.unwrap();
                        let bitenc2 = bitenc::BitEnc::from_bytes(&r2.seq);
                        let mut read2 = graph_read::Read::new(r2.id, pair.umi, pair.header.clone(), bitenc2, r2.quals.clone());
                        read2.set_paired();
                        read1.set_paired();
                        (Some(read2.clone()),
                         local_or_global_align(&aligner, &read2, &aligner_params, args.flag_global))
                    }
                };

                match alignment1.clone() {
                    Some(alignment) => {
                        good_alignments.push((alignment.clone(), read1.clone()));
                        let _ = out_bam.write(&mut read1.to_bam_record(Some(alignment.clone()), alignment2.clone()));
                    },
                    None => {
                        let _ = out_bam.write(&mut read1.to_unmapped_bam_record());
                    }
                };

                match (read2, alignment2) {
                    (Some(read), Some(alignment)) => {
                        good_alignments.push((alignment.clone(), read.clone()));
                        let _ = out_bam.write(&mut read.to_bam_record(alignment1.clone(),
                                                                      Some(alignment.clone())));
                    },
                    (Some(read), None) => {
                        let _ = out_bam.write(&mut read.to_unmapped_bam_record());
                    },
                    _ => {}
                }
            },
            None => break,
        }
    }

    let quals : Vec<u8>;
    if good_alignments.is_empty() {
        quals = vec![0; ref_seq.len()];
    } else {
        good_alignments.sort_by_key(|x| x.1.umi);
        let pileup = aligner.pileup(0, &good_alignments);

        quals = aligner.base_quals(0, &pileup, rt_error);
        assert!(quals.len() == ref_seq.len());
    }

    fastq_writer.write_fmt(format_args!("@{}\n{}\n+\n{}\n", ref_seq_name, ref_seq,
                           fastq::get_qual_string(&quals, QUAL_OFFSET))).unwrap();
}

/// Try aligning locally, if that doesn't work (eg. no good seed), align globally.
fn local_or_global_align(aligner: &sw::Aligner, read: &graph_read::Read,
                         aligner_params: &sw::AlignParams,
                         try_global: bool) -> Option<sw::Alignment> {

    let alignment = aligner.find_read_matches(read, aligner_params, 0.0, false);

    match alignment {
        Some(al) => {
            Some(al)
        },
        None => {
            if try_global {
                let global_alignment = aligner.global_align(&read.seq, &aligner_params);
                Some(global_alignment)
            } else {
                None
            }
        }
    }
}

/// Extract aligner parameters from input arguments
fn get_align_params(args: &Args) -> sw::AlignParams {
    let match_score = match args.flag_match_score {
        Some(c) => c,
        None => MATCH_SCORE,
    };
    let miss_score = match args.flag_miss_score {
        Some(c) => c,
        None => MISMATCH_SCORE,
    };
    let gap_open = match args.flag_gap_open {
        Some(c) => c,
        None => GAP_OPEN,
    };
    let gap_extend = match args.flag_gap_extend {
        Some(c) => c,
        None => GAP_EXTEND,
    };
    let clip = match args.flag_clip {
        Some(c) => c,
        None => CLIP,
    };

    let aligner_params = sw::AlignParams::new(match_score, miss_score, gap_open, gap_extend, clip);
    aligner_params
}

/// Match input reads against a set of reference sequences.
///
/// This does not perform a proper read alignment. It is meant for cases when you need to
/// check whether there are matches between the reads and the reference but don't care about
/// finding exactly where each read originated from (i.e. don't care about properly dealing with
/// multi-mapping).
///
/// # Usage
/// vdj_asm read-match <fasta> <fq-pref> <out-bam> <options>
///
/// Reads are read from <fq-pref>_1.fastq and <fq-pref>_2.fastq. Reads are reverse complemented according to the
/// rev-strand flag and then matches against the reference sequences in the input fasta.
/// This does NOT perform a read alignment. It just checks whether there exists a match of each read against
/// at least one of the reference sequences that exceeds the specified alignment score.
/// If such an alignment exists, it is written to the output bam. This alignment is NOT guaranteed to be
/// the best alignemnt of the read against the reference, just AN alignment that exceeds the minimum score.
/// If there is no such alignment, the read is written to the output as unmapped.
///
/// Notes on the bam output:
/// - The mate-mapped flag bit is not properly set.
/// - The MAPQ is always 0 because we don't check whether the read maps to multiple references.
pub fn get_matches(args: Args) {

    let mut aligner_params = get_align_params(&args);
    aligner_params.clip = 0.0; // always do local alignment

    let min_align_score = match args.flag_min_sw_score {
        Some(c) => c,
        None => 50.0,
    };
    let seed_len = match args.flag_seed {
        Some(c) => c,
        None => SEED_LEN,
    };

    let fq_pref = args.arg_fqpref.unwrap();

    let fq_path1 = fq_pref.clone() + "_1.fastq";
    let fq_path2 = fq_pref.clone() + "_2.fastq";

    let fasta = bio::io::fasta::Reader::from_file(Path::new(&args.arg_fasta.unwrap())).unwrap();
    let mut ref_seq_names = Vec::new();
    let mut ref_seqs = Vec::new();
    for record_option in fasta.records() {
        let record = record_option.unwrap();
        let seq_name = record.id().unwrap().to_owned();
        ref_seq_names.push(seq_name.clone());
        ref_seqs.push(utils::replace_ns(&String::from_utf8_lossy(record.seq()).into_owned(), &seq_name));
    }

    let mut header = bam::header::Header::new();
    let mut header_rec = bam::header::HeaderRecord::new(b"PG");
    header_rec.push_tag(b"ID", &"read-match");
    header.push_record(&header_rec);
    for (seq_name, seq) in ref_seq_names.iter().zip(ref_seqs.iter()) {
        bam_utils::add_ref_to_bam_header(&mut header, &seq_name, seq.len());
    }

    let aligner = sw::Aligner::new(&ref_seqs, seed_len);

    let mut out_bam = bam::Writer::from_path(&args.arg_outbam.unwrap(), &header).unwrap();

    let mut num_unmapped = 0;
    let mut npairs = 0;
    let mut fq_iter = fastq::CellrangerPairedFastqIter::new(Path::new(&fq_path1), Some(Path::new(&fq_path2)), args.flag_rev_strand);

    let cc_start = PreciseTime::now();
    loop {
        match fq_iter.next() {
            Some(pair) => {
                let r1 = pair.r1;
                let r2 = pair.r2.unwrap();

                let bitenc1 = bitenc::BitEnc::from_bytes(&r1.seq);
                let bitenc2 = bitenc::BitEnc::from_bytes(&r2.seq);
                let read_match = aligner.find_read_matches_sparse(&bitenc1, &aligner_params, min_align_score as f32);
                let mate_match = aligner.find_read_matches_sparse(&bitenc2, &aligner_params, min_align_score as f32);

                let mut read1 = graph_read::Read::new(r1.id, pair.umi, pair.header.clone(), bitenc1, r1.quals.clone());
                let mut read2 = graph_read::Read::new(r2.id, pair.umi, pair.header.clone(), bitenc2, r2.quals.clone());

                npairs += 1;

                read1.set_first_in_template();
                read2.set_last_in_template();
                read1.set_paired();
                read2.set_paired();
                if args.flag_rev_strand {
                    read1.set_reverse();
                } else {
                    read2.set_reverse();
                }
                let mut rec = read1.to_bam_record(read_match.clone(), mate_match.clone());
                let mut mate_rec = read2.to_bam_record(mate_match.clone(), read_match.clone());

                let _ = out_bam.write(&mut rec);
                let _ = out_bam.write(&mut mate_rec);

                num_unmapped += match read_match.clone() {
                    Some(_) => 0,
                    None => 1
                };
                num_unmapped += match mate_match.clone() {
                    Some(_) => 0,
                    None => 1
                };
            },
            None => break,
        }
    }
    println!("Got {} matches in {} sec", 2 * npairs - num_unmapped, cc_start.to(PreciseTime::now()));
    println!("Unmapped reads {}", num_unmapped);
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::io::{BufRead, BufReader};
    use rust_htslib::bam;
    use rust_htslib::bam::Read;
    use std::path::Path;
    use std::fs;
    use bio;
    use std::collections::{HashSet, HashMap};

    fn init_test(outdir: &str) {
        if Path::new(outdir).is_dir() {
            let _ = fs::remove_dir_all(outdir);
        }
        let _ = fs::create_dir_all(outdir);
    }

    fn get_seqs_by_barcode(fasta: &str) -> HashMap<String, HashSet<String>> {

        let mut seqs = HashMap::new();

        let reader = bio::io::fasta::Reader::new(fs::File::open(fasta).unwrap());
        for record in reader.records() {
            let rec = record.unwrap();
            let id = rec.id().unwrap().to_string().split("_").collect::<Vec<&str>>()[0].to_string();
            let seq = String::from_utf8(rec.seq().to_vec()).unwrap();
            if !seqs.contains_key(&id) {
                seqs.insert(id.clone(), HashSet::new());
            }
            (*seqs.get_mut(&id).unwrap()).insert(seq);
        }
        seqs
    }

    fn count_records(bam: &bam::Reader) -> (usize, usize) {
        let mut mapped_records = 0;
        let mut unmapped_records = 0;
        for record in bam.records() {
            let rec = record.unwrap();
            if rec.is_unmapped() {
                unmapped_records += 1;
            } else {
                mapped_records += 1;
            }
        }
        (mapped_records, unmapped_records)
    }

    #[test]
    /// Input test_asm.bam has 3 barcodes, 2 UMIs per barcode, one TRA, one TRB.
    /// The 3 barcodes have 5, 4, and 3 read-pairs from each chain respectively.
    /// The readpairs are inter-leaved (one from TRA, one from TRB).
    /// Each barcode starts with an unmapped pair, so we have
    /// unmapped pair first barcode
    /// TRB pair
    /// TRA pair
    /// (5 such TRA-TRB pairs)
    /// unmapped pair second barcode
    /// 4 TRA-TRB pairs
    /// unmapped pair third barcode
    /// 3 TRA-TRB pairs
    ///
    /// The single end file has a different structure:
    /// TRB read first barcode
    /// unmapped read
    /// TRA read
    /// umapped read
    /// (5 such 4-read blocks)
    /// And similarly for the other 2 barcodes which have 4 and 3 4-read blocks
    /// respectively.
    fn test_vdj_asm() {

        let outdir = "test/outputs/asm";
        let true_fasta_name = "test/inputs/test_asm.fasta";
        let true_seqs = get_seqs_by_barcode(&true_fasta_name);

        let mut args = Args {cmd_asm: true,
                        cmd_base_quals: false,
                        cmd_read_match: false,
                        arg_inbam: Some("test/inputs/test_asm.bam".to_string()),
                        arg_inpref: None,
                        arg_outdir: Some(outdir.to_string()),
                        arg_fasta: None,
                        arg_fqpref: None,
                        arg_outbam: None,
                        flag_cons: false,
                        flag_kmers: Some(0),
                        flag_min_contig: Some(150),
                        flag_npaths: Some(100),
                        flag_rev_strand: false,
                        flag_plot: false,
                        flag_min_sw_score: Some(50.0),
                        flag_rt_error: Some(0.0001),
                        flag_min_qual: Some(2),
                        flag_nx: Some(0.99),
                        flag_qual_factor: Some(0.9),
                        flag_score_factor: Some(0.6),
                        flag_frac_reads: Some(0.3),
                        flag_min_umi_reads: Some(1),
                        flag_subsample_rate: Some(1.0),
                        flag_fast_align: false,
                        flag_match_score: None,
                        flag_miss_score: None,
                        flag_gap_open: None,
                        flag_gap_extend: None,
                        flag_clip: None,
                        flag_seed: Some(20),
                        flag_single_end: false,
                        flag_global: false,
                        };

        init_test(outdir);

        {
            vdj_asm(args.clone(), 11); // use all read-pairs

            let assembled_seqs = get_seqs_by_barcode(&"test/outputs/asm/test_asm.fasta");

            for (bc, seqs) in true_seqs.iter() {
                assert_eq!(*seqs, *assembled_seqs.get(bc).unwrap());
            }

            let bam = bam::Reader::from_path(Path::new(&"test/outputs/asm/test_asm.bam")).ok().expect("Error reading test_asm.bam");
            let (mapped, unmapped) = count_records(&bam);
            assert_eq!(unmapped, 6);
            assert_eq!(mapped, 48);

            let header = bam.header();
            assert_eq!(header.target_count(), 6);
            let header_names : HashSet<String> = header.target_names().iter().map(|x| String::from_utf8_lossy(x).into_owned()).collect();
            let summary = BufReader::new(File::open("test/outputs/asm/test_asm_summary.tsv").unwrap());
            // Skip the header and test that all other contig names are present in the header
            for wrapped_line in summary.lines().skip(1) {
                let line = wrapped_line.unwrap();
                let contig_name = line.split('\t').nth(1).unwrap();
                assert!(header_names.contains(contig_name));
            }

            let summary = BufReader::new(File::open("test/outputs/asm/test_asm_umi_summary.tsv").unwrap());
            for wrapped_line in summary.lines().skip(1) {
                let line = wrapped_line.unwrap();
                let contig_name = line.split('\t').nth(6).unwrap();
                assert!(contig_name == "" || header_names.contains(contig_name));
            }
        }
        {
            // Use only the first 3 read-pairs per barcode per chain.
            // So in all barcodes we'll use as much data as in the last barcode.
            // 4 + 2 readpairs will be left unmapped.
            vdj_asm(args.clone(), 7);

            let assembled_seqs = get_seqs_by_barcode(&"test/outputs/asm/test_asm.fasta");

            for (_, seqs) in assembled_seqs.iter() {
                assert_eq!(*seqs, *true_seqs.get(&"GAACGGGT".to_string()).unwrap());
            }

            let bam = bam::Reader::from_path(Path::new(&"test/outputs/asm/test_asm.bam")).ok().expect("Error reading test_asm.bam");
            let (_, unmapped) = count_records(&bam);
            assert_eq!(unmapped, 18);

            let header = bam.header();
            assert_eq!(header.target_count(), 6);
        }
        {
            args.flag_min_contig = Some(250);

            vdj_asm(args.clone(), 11);

            let bam = bam::Reader::from_path(Path::new(&"test/outputs/asm/test_asm.bam")).ok().expect("Error reading test_asm.bam");
            let (_, unmapped) = count_records(&bam);
            assert_eq!(unmapped, 18);

            // We won't assemble anything from the last barcode.
            let header = bam.header();
            assert_eq!(header.target_count(), 4);
        }
        {
            // Now pretend the data are single-end (but still reading the paired end file)
            // The names and flags of the BAM are ignored, so this should still work.
            args.flag_single_end = true;
            args.flag_min_contig = Some(150);

            vdj_asm(args.clone(), 22); // use all reads

            let assembled_seqs = get_seqs_by_barcode(&"test/outputs/asm/test_asm.fasta");

            for (bc, seqs) in true_seqs.iter() {
                assert_eq!(*seqs, *assembled_seqs.get(bc).unwrap());
            }

            let bam = bam::Reader::from_path(Path::new(&"test/outputs/asm/test_asm.bam")).ok().expect("Error reading test_asm.bam");
            let (mapped, unmapped) = count_records(&bam);
            assert_eq!(unmapped, 6);
            assert_eq!(mapped, 48);

            let header = bam.header();
            assert_eq!(header.target_count(), 6);
        }
        {
            // Truly single end BAM
            args.arg_inbam = Some("test/inputs/test_asm_se.bam".to_string());
            args.flag_single_end = true;
            args.flag_min_contig = Some(150);

            vdj_asm(args.clone(), 20); // use all reads

            let assembled_seqs = get_seqs_by_barcode(&"test/outputs/asm/test_asm_se.fasta");

            for (bc, seqs) in true_seqs.iter() {
                assert_eq!(*seqs, *assembled_seqs.get(bc).unwrap());
            }

            let bam = bam::Reader::from_path(Path::new(&"test/outputs/asm/test_asm_se.bam")).ok().expect("Error reading test_asm.bam");
            let (mapped, unmapped) = count_records(&bam);
            // every other read in the input is unmapped
            assert_eq!(unmapped, 24);
            assert_eq!(mapped, 24);

            let header = bam.header();
            assert_eq!(header.target_count(), 6);
        }
        {
            vdj_asm(args.clone(), 11);

            let assembled_seqs = get_seqs_by_barcode(&"test/outputs/asm/test_asm_se.fasta");

            let bam = bam::Reader::from_path(Path::new(&"test/outputs/asm/test_asm_se.bam")).ok().expect("Error reading test_asm.bam");
            let (mapped, unmapped) = count_records(&bam);
            assert_eq!(unmapped, 30);
            assert_eq!(mapped, 18);

            for (_, seqs) in assembled_seqs.iter() {
                assert_eq!(*seqs, *true_seqs.get(&"GAACGGGT".to_string()).unwrap());
            }
        }
        {
            vdj_asm(args.clone(), 2); // second chain not assembled, the second read we read is filtered

            let bam = bam::Reader::from_path(Path::new(&"test/outputs/asm/test_asm_se.bam")).ok().expect("Error reading test_asm.bam");
            let header = bam.header();
            assert_eq!(header.target_count(), 3);
        }
        {
            args.flag_cons = true;
            args.flag_min_umi_reads = Some(1);
            args.flag_seed = Some(1);
            args.flag_min_sw_score = Some(50.0);
            vdj_asm(args.clone(), 20);

            let bam = bam::Reader::from_path(Path::new(&"test/outputs/asm/test_asm_se.bam")).ok().expect("Error reading test_asm.bam");
            let header = bam.header();
            assert_eq!(header.target_count(), 3);
        }
    }

    fn read_one_fq_record(fq_file: &str) -> (String, String, Vec<u8>) {
        let mut fastq = bio::io::fastq::Reader::from_file(fq_file).ok().expect("Could not open FASTQ file");
        let mut record = bio::io::fastq::Record::new();
        let _ = fastq.read(&mut record);
        let seq_name = record.id().unwrap().to_owned();
        let seq = String::from_utf8_lossy(record.seq()).into_owned();
        let qual = record.qual().to_vec();
        (seq_name, seq, qual)
    }

    #[test]
    fn test_base_qualities() {

        let outdir = "test/outputs/base_quals";
        let out_bam = "test/outputs/base_quals/test_base_quals.bam";
        let out_fastq = "test/outputs/base_quals/test_base_quals.fastq";
        let true_fasta_name = "test/inputs/base_quals/test_base_quals.fasta";

        let fasta = bio::io::fasta::Reader::from_file(true_fasta_name).ok().expect("Could not open input FASTA");
        let mut record_iter = fasta.records();
        let record = record_iter.next().unwrap().ok().expect("Error reading base qual true FASTA");
        let ref_seq_name = record.id().unwrap().to_owned();
        let ref_seq = String::from_utf8_lossy(record.seq()).into_owned();

        let args = Args {cmd_asm: false,
                        cmd_base_quals: true,
                        cmd_read_match: false,
                        arg_inbam: None,
                        arg_inpref: Some("test/inputs/base_quals/test_base_quals".to_string()),
                        arg_outdir: Some(outdir.to_string()),
                        arg_fasta: None,
                        arg_fqpref: None,
                        arg_outbam: None,
                        flag_cons: false,
                        flag_kmers: Some(0),
                        flag_min_contig: Some(150),
                        flag_npaths: Some(100),
                        flag_rev_strand: false,
                        flag_plot: false,
                        flag_min_sw_score: Some(50.0),
                        flag_rt_error: Some(0.0001),
                        flag_min_qual: Some(2),
                        flag_nx: Some(0.99),
                        flag_qual_factor: Some(0.9),
                        flag_score_factor: Some(0.6),
                        flag_frac_reads: Some(0.3),
                        flag_min_umi_reads: Some(1),
                        flag_subsample_rate: Some(1.0),
                        flag_fast_align: false,
                        flag_match_score: None,
                        flag_miss_score: None,
                        flag_gap_open: None,
                        flag_gap_extend: None,
                        flag_clip: None,
                        flag_seed: Some(20),
                        flag_single_end: false,
                        flag_global: false,
                        };

        init_test(outdir);

        {
            get_base_qualities(args.clone());

            let bam = bam::Reader::from_path(Path::new(&out_bam)).ok().expect("Error reading test_base_quals.bam");
            let (mapped, unmapped) = count_records(&bam);
            // Two records left unmapped because they didn't have a good seed for mapping.
            assert_eq!(unmapped, 2);
            assert_eq!(mapped, 6);

            let (seq_name, seq, qual) = read_one_fq_record(&out_fastq);
            // Make sure we left the sequence unchanged
            assert_eq!(seq_name, ref_seq_name);
            assert_eq!(seq, ref_seq);
            assert!(qual[0] > qual[201]);
            assert!(qual[201] > qual[401]);
            // zero qual after the first 400 bases
            assert_eq!(qual[401], 33);
        }
        {
            let mut args2 = args.clone();
            args2.flag_global = true;
            get_base_qualities(args2);

            let bam = bam::Reader::from_path(Path::new(&out_bam)).ok().expect("Error reading test_base_quals.bam");
            let (mapped, unmapped) = count_records(&bam);
            assert_eq!(unmapped, 0);
            assert_eq!(mapped, 8);
        }
        {
            let mut args3 = args.clone();
            args3.flag_global = true;
            args3.flag_single_end = true;
            get_base_qualities(args3);

            let bam = bam::Reader::from_path(Path::new(&out_bam)).ok().expect("Error reading test_base_quals.bam");
            let (mapped, unmapped) = count_records(&bam);
            assert_eq!(unmapped, 0);
            assert_eq!(mapped, 8);
        }
        {
            let mut args4 = args.clone();
            args4.arg_inpref = Some("test/inputs/base_quals/test_base_quals_invalid_bases".to_string());
            get_base_qualities(args4);

            let new_out_bam = "test/outputs/base_quals/test_base_quals_invalid_bases.bam";
            let bam = bam::Reader::from_path(Path::new(&new_out_bam)).ok().expect("Error reading output bam");
            let (mapped, _) = count_records(&bam);
            assert_eq!(mapped, 0);
        }
    }

    fn revcomp_fastq(infile: &str, outfile: &str) {
        let fq_reader = bio::io::fastq::Reader::from_file(infile).unwrap();
        let mut fq_writer = bio::io::fastq::Writer::to_file(outfile).unwrap();
        for ur in fq_reader.records() {
            let rec = ur.unwrap();
            let mut new_qual = rec.qual().to_vec();
            new_qual.reverse();
            let new_seq = String::from_utf8(bio::alphabets::dna::revcomp(rec.seq())).unwrap();
            let _ = fq_writer.write(rec.id().unwrap(), None, new_seq.as_bytes(), &new_qual);
        }
    }

    #[test]
    fn test_get_matches() {

        let outdir = "test/outputs/base_quals";
        let out_bam = "test/outputs/base_quals/test_match.bam";

        let args = Args {cmd_asm: false,
                        cmd_base_quals: false,
                        cmd_read_match: true,
                        arg_inbam: None,
                        arg_inpref: None,
                        arg_outdir: Some(outdir.to_string()),
                        arg_fasta: Some("test/inputs/base_quals/test_base_quals.fasta".to_string()),
                        arg_fqpref: Some("test/inputs/base_quals/test_base_quals".to_string()),
                        arg_outbam: Some(out_bam.to_string()),
                        flag_cons: false,
                        flag_kmers: Some(0),
                        flag_min_contig: Some(150),
                        flag_npaths: Some(100),
                        flag_rev_strand: false,
                        flag_plot: false,
                        flag_min_sw_score: Some(50.0),
                        flag_rt_error: Some(0.0001),
                        flag_min_qual: Some(2),
                        flag_nx: Some(0.99),
                        flag_qual_factor: Some(0.9),
                        flag_score_factor: Some(0.6),
                        flag_frac_reads: Some(0.3),
                        flag_min_umi_reads: Some(1),
                        flag_subsample_rate: Some(1.0),
                        flag_fast_align: false,
                        flag_match_score: None,
                        flag_miss_score: None,
                        flag_gap_open: None,
                        flag_gap_extend: None,
                        flag_clip: None,
                        flag_seed: Some(20),
                        flag_single_end: false,
                        flag_global: false,
                        };

        {
            get_matches(args.clone());

            let bam = bam::Reader::from_path(Path::new(&out_bam)).ok().expect("Error reading test_match.bam");
            let (mapped, unmapped) = count_records(&bam);
            // Two records left unmapped because they didn't have a good seed for mapping.
            assert_eq!(unmapped, 2);
            assert_eq!(mapped, 6);
        }
        {
            // Test reverse complemented inputs
            revcomp_fastq("test/inputs/base_quals/test_base_quals_1.fastq",
                          "test/outputs/base_quals/test_base_quals_rc_1.fastq");
            revcomp_fastq("test/inputs/base_quals/test_base_quals_2.fastq",
                          "test/outputs/base_quals/test_base_quals_rc_2.fastq");
            let mut args5 = args.clone();
            args5.flag_rev_strand = true;
            args5.arg_fqpref = Some("test/outputs/base_quals/test_base_quals_rc".to_string());

            let bam = bam::Reader::from_path(Path::new(&out_bam)).ok().expect("Error reading test_match.bam");
            let (mapped, unmapped) = count_records(&bam);
            // Two records left unmapped because they didn't have a good seed for mapping.
            assert_eq!(unmapped, 2);
            assert_eq!(mapped, 6);
        }
        {
            let mut args2 = args.clone();
            // Input FASTA is all Ns.
            args2.arg_fasta = Some("test/inputs/base_quals/test_base_quals_invalid_bases_1.fasta".to_string());
            get_matches(args2);

            let bam = bam::Reader::from_path(Path::new(&out_bam)).ok().expect("Error reading test_match.bam");
            let (_, unmapped) = count_records(&bam);
            assert_eq!(unmapped, 8);
        }
        {
            // One all Ns, one good sequence. The presence of the invalid one shouldn't prevent us
            // from getting mathces on the good one.
            let mut args3 = args.clone();
            args3.arg_fasta = Some("test/inputs/base_quals/test_base_quals_invalid_bases_2.fasta".to_string());
            get_matches(args3);

            let bam = bam::Reader::from_path(Path::new(&out_bam)).ok().expect("Error reading test_match.bam");
            let (mapped, _) = count_records(&bam);
            assert_eq!(mapped, 6);
        }
        {
            let mut args4 = args.clone();
            args4.arg_fasta = Some("test/inputs/base_quals/test_base_quals_invalid_bases_2.fasta".to_string());
            // First pair has two reads with Ns. Second and third have only one.
            // 4th has both valid reads but no good seed. Last one should be good.
            args4.arg_fqpref = Some("test/inputs/base_quals/test_base_quals_invalid_bases".to_string());
            get_matches(args4);

            let bam = bam::Reader::from_path(Path::new(&out_bam)).ok().expect("Error reading test_match.bam");
            let (mapped, _) = count_records(&bam);
            assert_eq!(mapped, 8);
        }
    }
}
