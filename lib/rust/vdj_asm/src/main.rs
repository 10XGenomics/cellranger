//
// Copyright (c) 2017 10x Genomics, Inc. All rights reserved.
//

#![crate_name = "vdj_asm"]
#![allow(dead_code)]

extern crate bio;
extern crate debruijn;
extern crate flate2;
extern crate itertools;
extern crate docopt;
extern crate time;
extern crate num;
extern crate rand;
extern crate rust_htslib;
extern crate lz4;

#[macro_use]
extern crate serde_derive;
extern crate serde;
#[macro_use]
extern crate serde_json;
#[macro_use]
extern crate lazy_static;


pub mod graph;
mod constants;
mod fastq;
mod utils;
mod sw;
mod asm;
mod asm_helper;
mod bam_utils;
mod graph_read;
mod metrics;
mod perf;

use rust_htslib::bam;
use rust_htslib::bam::Read;

use debruijn::dna_string;
use metrics::AssemblyMetrics;
use asm::UmiCounter;

use itertools::Itertools;
use std::path::{Path, PathBuf};
use std::io::{BufWriter, Write};
use docopt::Docopt;
use constants::{MAX_NUM_READPAIRS, MATCH_SCORE, MISMATCH_SCORE, GAP_OPEN, GAP_EXTEND, CLIP, DUMMY_CONTIG_NAME, QUAL_OFFSET, KMER_LEN_BANDED_ALIGN, WINDOW_SIZE_BANDED_ALIGN};
use std::fs;
use std::fs::File;

use time::PreciseTime;

const USAGE: &'static str = "
Usage:
vdj_asm asm <inbam> <outdir> [--barcode=<BC>][--plot][--plot-json][--kmers=<NUM>][--min-contig=<NUM>][--frac-reads=<NUM>][--reads-per-barcode=<NUM>][--score-factor=<NUM>][--qual-factor=<NUM>][--rt-error=<NUM>][--min-qual=<NUM>][--match-score=<NUM>][--miss-score=<NUM>][--gap-open=<NUM>][--gap-extend=<NUM>][--min-sw-score=<NUM>][--min-umi-reads=<NUM>][--cons][--single-end][--subsample-rate=<NUM>][--use-unmapped][--mixture-filter]
vdj_asm base-quals <inpref> <outdir> [--single-end][--rev-strand][--rt-error=<NUM>][--match-score=<NUM>][--miss-score=<NUM>][--gap-open=<NUM>][--gap-extend=<NUM>][--seed=<NUM>]
vdj_asm read-match [--ref=FASTA][--r1=FASTQ][--r2=FASTQ][--outbam=BAM][--rev-strand][--match-score=<NUM>][--miss-score=<NUM>][--gap-open=<NUM>][--gap-extend=<NUM>][--seed=<NUM>][--min-sw-score=<NUM>]
vdj_asm (-h | --help)

Options:
   --plot                 Create gfa file with assembly graph.
   --kmers=<NUM>          Minimum number of occurrences to consider a kmer
   --min-contig=<NUM>     Minimum output contig length
   --rev_strand           Data are RF (default is FR).
   --min-sw-score=<NUM>   Minimum SW score to consider an alignment.
   --rt-error=<NUM>       RT error rate.
   --barcode=<BC>         Only process selected barcode.
   --trace=<STR>
   --qual-factor=<NUM>
   --score-factor=<NUM>
   --frac-reads=<NUM>
   --reads-per-barcode=<NUM>
   --min-umi-reads=<NUM>
   --ref=FASTA
   --r1=FASTQ
   --r2=FASTQ
   --outbam=BAM
   --match-score=<NUM>
   --miss-score=<NUM>
   --gap-open=<NUM>
   --gap-extend=<NUM>
   --seed=<NUM>
   --subsample-rate=<NUM>
   --cons                 Compute single consensus sequence
   --single-end           Input data are single end
   --use-unmapped         Use unmapped reads in assembly (ignore filter)
   --mixture-filter       Enable filtering of mixed contigs (probable chimeras)
";

#[derive(Debug, Deserialize, Clone, Default)]
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
    flag_rev_strand: bool,
    flag_plot: bool,
    flag_plot_json: bool,
    flag_min_sw_score: Option<f64>,
    flag_rt_error: Option<f64>,
    flag_min_qual: Option<u8>,
    flag_qual_factor: Option<f64>,
    flag_score_factor: Option<f64>,
    flag_frac_reads: Option<f64>,
    flag_reads_per_barcode: Option<usize>,
    flag_min_umi_reads: Option<usize>,
    flag_subsample_rate: Option<f64>,
    flag_ref: String,
    flag_r1: String,
    flag_r2: Option<String>,
    flag_outbam: String,
    flag_match_score: Option<i32>,
    flag_miss_score: Option<i32>,
    flag_gap_open: Option<i32>,
    flag_gap_extend: Option<i32>,
    flag_clip: Option<i32>,
    flag_seed: Option<usize>,
    flag_barcode: Option<String>,
    flag_single_end: bool,
    flag_use_unmapped: bool,
    flag_mixture_filter: bool
}

fn main() {

    let mut args: Args = Docopt::new(USAGE).and_then(|d| d.deserialize()).unwrap_or_else(|e| e.exit());

    if args.cmd_asm {
        if args.flag_cons {
            args.flag_min_umi_reads = Some(1);
            args.flag_min_contig = Some(1);
        }
        vdj_asm(args);
    } else if args.cmd_base_quals {
        get_base_qualities(args);
    } else {
        get_matches(args);
    }
}

pub fn vdj_asm(args: Args) {

    let input_bam_filename = args.arg_inbam.clone().unwrap();
    let out_dirname = args.arg_outdir.clone().unwrap();

    // setup shared output files
    let out_prefix = Path::new(&out_dirname).join(Path::new(&(input_bam_filename)).file_stem().unwrap())
                                          .into_os_string().into_string().unwrap();

    let fasta_file = File::create(out_prefix.to_string() + ".fasta").expect("Could not create fasta file");
    let fasta_writer = BufWriter::new(&fasta_file);

    let fastq_file = File::create(out_prefix.to_string() + ".fastq").expect("Could not create fastq file");
    let fastq_writer = BufWriter::new(&fastq_file);

    let summary_file = File::create(out_prefix.to_string() + "_summary.tsv").expect("Could not create summary file");
    let mut summary_writer = BufWriter::new(&summary_file);
    summary_writer.write_fmt(format_args!("{}\t{}\t{}\t{}\t{}\t{}\n",
                                          "barcode", "contig_name",
                                          "num_reads", "num_pairs",
                                          "num_umis", "umi_list")).unwrap();

    let umi_summary_file = File::create(out_prefix.to_string() + "_umi_summary.tsv").expect("Could not create UMI summary file");
    let mut umi_summary_writer = BufWriter::new(&umi_summary_file);
    umi_summary_writer.write_fmt(format_args!("{}\t{}\t{}\t{}\t{}\t{}\t{}\n", "barcode", "umi_id", "umi",
                                              "reads", "min_umi_reads", "good_umi", "contigs")).unwrap();

    let mut assembly_outs = asm::AssemblyOuts::new(fasta_writer, fastq_writer, summary_writer, umi_summary_writer);


    let in_bam = bam::Reader::from_path(&Path::new(&input_bam_filename)).ok().expect("Error opening input bam.");

    // Will create one bam per barcode and merge incrementally
    // Can't append to a single bam, because header (contig names) is not known in advance.
    let bam_merge_frequency = 50; // merge every N BAMs / barcodes to reduce number of files.
    let mut out_bams = Vec::new();

    let bam_iter = in_bam.records();

    let mut metrics = AssemblyMetrics::default();

    let get_bc = |rec: &bam::Record| {
        let header = fastq::CellrangerFastqHeader::new(String::from_utf8_lossy(rec.qname()).into_owned());
        header.get_barcode().cloned()
    };

    for (barcode, records) in &bam_iter.map(|x| x.expect("trouble when reading")).group_by(get_bc) {

        // Handle non barcoded & barcoded reads separately
        match barcode {
            Some(bc) => {

                // If a single barcode was requested, only process that bc.
                if args.flag_barcode.as_ref().map_or(false, |req_bc| &bc != req_bc) {
                    continue
                }

                let out_bam = asm_bc(records, &bc, &out_prefix, &mut metrics, &mut assembly_outs, &args);
                out_bams.push(out_bam);

                // if we've reached the BAM limit, squash the current set
                if out_bams.len() >= bam_merge_frequency {
                    let squashed_bam_filename = format!("{}_{}_{}.bam", out_prefix, bc.clone(), out_bams.len());
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

            },
            None => {
                // FIXME - what to with unbarcoded reads?
            },
        }
    }

    // Write metrics summary
    {
        let mut metrics_summary_file = File::create(out_prefix.to_string() + "_metrics_summary.json")
            .expect("Could not create metrics summary file");
        metrics::write_summary(&mut metrics_summary_file, &metrics)
            .expect("Could not write metrics summary");
    }

    assembly_outs.close();

    bam_utils::concatenate_bams(&out_bams, &(out_prefix.to_string() + ".bam"));
    bam_utils::sort_and_index(&(out_prefix.to_string() + ".bam"), &(out_prefix.to_string() + "_sorted.bam"));
}

struct ReservoirSampler<T> {
    max_items: usize,
    items_seen: usize,
    items: Vec<T>,
    rand: rand::StdRng,
}

impl<T> ReservoirSampler<T> {
    pub fn new(max_items: usize, seed: usize) -> ReservoirSampler<T> {

        let seed: &[_] = &[seed, 1, 1, 1];
        let rng: rand::StdRng = rand::SeedableRng::from_seed(seed);

        ReservoirSampler {
            max_items: max_items,
            items_seen: 0,
            items: Vec::new(),
            rand: rng,
        }
    }

    /// Add a new item to the sampler
    pub fn add(&mut self, item: T) {
        use rand::Rng;

        self.items_seen += 1;

        if self.items.len() < self.max_items {
            self.items.push(item);
        } else {
            // Swap out an existing item according
            let idx = self.rand.gen_range(0, self.items_seen);
            if idx < self.items.len() {
                self.items[idx] = item;
            }
        }
    }

    pub fn done(self) -> Vec<T> {
        self.items
    }
}


fn asm_bc<I: Iterator<Item=bam::Record>, T: Write>(mut bam_iter: I, barcode: &str, out_prefix: &str,
                                         metrics: &mut AssemblyMetrics, assembly_outs: &mut asm::AssemblyOuts<T>,
                                         args: &Args) -> String {

    let single_end = args.flag_single_end;

    let max_readpairs_per_bc = match args.flag_reads_per_barcode {
        Some(v) => v,
        None => MAX_NUM_READPAIRS,
    };

    let min_umi_reads = match args.flag_min_umi_reads {
        Some(c) => c,
        None => 1,
    };

    let mut sampled_reads = ReservoirSampler::<(graph_read::Read, Option<graph_read::Read>)>::new(max_readpairs_per_bc, 1);
    let mut umi_counts = UmiCounter::new();
    let mut npairs = 0;

    // Iterate until we have enough reads to assemble
    loop {
        // get next read
        let record = match bam_iter.next() {
            Some(r) => r,
            _ => break,
        };

        let header = fastq::CellrangerFastqHeader::new(String::from_utf8_lossy(record.qname()).into_owned());
        let umi = umi_counts.get_umi(&header);

        let rec =
            if single_end {
                let mut read1 = graph_read::Read::from_bam_record(2 * npairs, umi, &record);
                read1.unset_paired();
                umi_counts.count_umi(umi);
                (read1, None)
            } else {
                let mut read1 = graph_read::Read::from_bam_record(2 * npairs, umi, &record);
                let mate_record = bam_iter.next().expect("Error while reading input bam");
                let mut read2 = graph_read::Read::from_bam_record(2 * npairs + 1, umi, &mate_record);
                read1.set_paired();
                read2.set_paired();
                umi_counts.count_umi(umi);
                umi_counts.count_umi(umi);
                (read1, Some(read2))
            };

        sampled_reads.add(rec);
        npairs += 1;
    }

    fn rec_unmapped(rec: &(graph_read::Read, Option<graph_read::Read>)) -> bool {
        match rec {
            &(ref r1, None) => r1.is_unmapped(),
            &(ref r1, Some(ref r2)) => r1.is_unmapped() && r2.is_unmapped(),
        }
    }

    // Single array of selected reads
    let mut reads = Vec::new();
    for rec in sampled_reads.done() {

        if !args.flag_use_unmapped && rec_unmapped(&rec) {
            continue;
        }

        let (read1, read2) = rec;
        reads.push(read1);
        match read2 {
            Some(r2) => reads.push(r2),
            _ => ()
        }
    }

    // Only keep sequences with good UMIs - these will be assembled
    println!("Reads/UMI cutoff: {}", min_umi_reads);
    let good_umis = umi_counts.get_good_umis(min_umi_reads);

    println!("Observed {} read pairs with {} distinct UMIs prior to subsampling", npairs, umi_counts.len());
    let recalculate_umi_counts = (npairs as usize) > max_readpairs_per_bc;
    if recalculate_umi_counts {
        umi_counts.reset_counts();
    }
    let mut assembled_reads = Vec::new();
    for read in reads.iter() {
        if good_umis.contains(&(read.umi)) && read.umi > 0 {
            assembled_reads.push(read);
        }
        if recalculate_umi_counts {
            umi_counts.count_umi(read.umi);
        }
    }

    println!("Barcode {}: Assembling {:?} reads, {:?} distinct UMIs", barcode, reads.len(), umi_counts.count_good_umis(1));
    println!("Good UMIs: {:?}", good_umis);

    // Report number of reads used
    println!("Reads to use in assembly: {}", assembled_reads.len());
    let assemblable_read_pairs = match single_end {
        true => assembled_reads.len() as u64,
        false => assembled_reads.len() as u64 / 2,
    };

    metrics.assemblable_read_pairs_by_bc.insert(barcode.to_string(), assemblable_read_pairs);


    // Args
    let scoring = get_scoring(&args);

    let min_kmer_count = match args.flag_kmers {
        Some(c) => c,
        None => 2,
    };
    let min_contig_len = match args.flag_min_contig {
        Some(c) => c,
        None => 150,
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

    let frac_path_starts = 0.99f64; // Terminate path enumeration
    // after starting from nodes containing this fraction of reads

    // Assemble
    let cc_start = PreciseTime::now();
    let graph;
    let contigs;
    if args.flag_mixture_filter {
        assert!(args.flag_cons == false);
        let result = asm::assemble_reads_with_mixture_filter(&mut assembled_reads,
                                                              &umi_counts,
                                                              min_kmer_count as usize,
                                                              qual_factor,
                                                              rt_error,
                                                              min_qual,
                                                              min_contig_len,
                                                              min_align_score,
                                                              score_factor,
                                                              scoring,
                                                              frac_path_starts,
                                                              frac_path_reads,
                                                              !single_end);
        graph = result.0;
        contigs = result.1;
    } else {
        let result = asm::assemble_reads(&mut assembled_reads,
                                                &umi_counts,
                                                min_kmer_count as usize,
                                                qual_factor,
                                                rt_error,
                                                min_qual,
                                                min_contig_len,
                                                frac_path_reads,
                                                min_align_score,
                                                score_factor,
                                                args.flag_cons,
                                                scoring,
                                                !single_end,
                                                frac_path_starts);
        graph = result.0;
        contigs = result.1;
    }

    println!("Assembled {} contigs in {} sec", contigs.len(), cc_start.to(PreciseTime::now()));
    println!("maxrss: {}", perf::getrusage().ru_maxrss as u64);


    if args.flag_plot {
        println!("Writing gfa");
        let path = out_prefix.to_string() + "_" + &barcode + "_graph.gfa";
        graph.graph.to_gfa(path);
    }

    if args.flag_plot_json {

        let f = |umis: &Vec<u32>| {
            let mut uu = umis.clone();
            uu.sort();
            uu.dedup();

            json!({
                "NU": uu.len(),
                "NR": umis.len(),
            })
        };

        let contig_seqs: Vec<String> = contigs.iter().map(|&(ref s,_,_,_)| s.clone()).collect();
        let contig_json = json!({"contigs": contig_seqs});

        let path = out_prefix.to_string() + "_" + &barcode + "_graph.json";
        let mut file = std::fs::File::create(path).unwrap();
        graph.graph.to_json_rest(f, &mut file, Some(contig_json));
    }

    let (header, contig_names) = get_assembly_contig_header(&contigs.iter().map(|x| x.0.clone()).collect(), &barcode);

    let out_bam_filename = format!("{}_{}.bam", out_prefix, barcode);
    let mut out_bam = bam::Writer::from_path(&Path::new(&out_bam_filename), &header).unwrap();

    asm::write_assembly_results(&reads,
                                &umi_counts,
                                &good_umis,
                                assembly_outs,
                                &mut out_bam,
                                &contigs,
                                &contig_names,
                                &barcode,
                                min_umi_reads,
                                single_end);

    out_bam_filename
}

fn get_assembly_contig_header(seqs: &Vec<String>, pref: &str) -> (bam::header::Header, Vec<String>) {
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

    let scoring = get_scoring(&args);
    let file_pref = args.arg_inpref.unwrap();
    let out_dir = args.arg_outdir.unwrap();

    let fasta = bio::io::fasta::Reader::from_file(Path::new(&(file_pref.clone() + ".fasta"))).ok()
        .expect("Could not open input FASTA");
    let mut record_iter = fasta.records();
    let record = record_iter.next().unwrap().ok()
        .expect("Error reading sequence from fasta file");
    let ref_seq_name = record.id().to_owned();
    let ref_seq = utils::replace_ns(&String::from_utf8_lossy(record.seq()).into_owned(), &ref_seq_name);

    let record = record_iter.next();
    if record.is_some() {
        println!("FASTA file has more than one sequence. Only the first one will be used.")
    }

    let rt_error = match args.flag_rt_error {
        Some(c) => c,
        None => 0.0001,
    };

    let single_end = args.flag_single_end;

    let mut fq_iter = match single_end {
        true => {
            let name1 = utils::find_file_maybe_compressed(&(file_pref.clone() + ".fastq"))
                .expect("Couldn't find FASTQ file");
            fastq::CellrangerPairedFastqIter::new(&name1, None, false)
        },
        false => {
            let name1 = utils::find_file_maybe_compressed(&(file_pref.clone() + "_1.fastq"))
                .expect("Couldn't find R1 FASTQ file");
            let name2 = utils::find_file_maybe_compressed(&(file_pref.clone() + "_2.fastq"))
                .expect("Coldn't find R2 FASTQ file");
            fastq::CellrangerPairedFastqIter::new(&name1, Some(&name2), false)
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

    let refs = vec![ref_seq.clone()];
    let align_helper = sw::AlignHelper::new(&refs, scoring.clone(), KMER_LEN_BANDED_ALIGN, WINDOW_SIZE_BANDED_ALIGN);
    let min_align_score = match args.flag_min_sw_score {
        Some(c) => c,
        None => 50.0,
    };

    let mut good_alignments = Vec::new();

    loop {
        match fq_iter.next() {
            Some(pair) => {
                let r1 = pair.r1;

                let bitenc1 = dna_string::DnaString::from_bytes(&r1.seq);
                let mut read1 = graph_read::Read::new(r1.id, pair.umi, pair.header.clone(), bitenc1, r1.quals.clone());
                let al_pack1 = align_helper.find_read_matches(&read1, min_align_score as i32);

                let (read2, al_pack2) = match single_end {
                    true => (None, None),
                    false => {
                        let r2 = pair.r2.unwrap();
                        let bitenc2 = dna_string::DnaString::from_bytes(&r2.seq);
                        let mut read2 = graph_read::Read::new(r2.id, pair.umi, pair.header.clone(), bitenc2, r2.quals.clone());
                        read2.set_paired();
                        read1.set_paired();
                        (Some(read2.clone()),
                        align_helper.find_read_matches(&read2, min_align_score as i32))
                    }
                };

                match al_pack1.clone() {
                    Some(al_pack) => {
                        good_alignments.push((al_pack.clone(), read1.clone()));
                        let _ = out_bam.write(&mut read1.to_bam_record(&Some(al_pack), &al_pack2));
                    },
                    None => {
                        let _ = out_bam.write(&mut read1.to_unmapped_bam_record());
                    }
                };

                match (read2, al_pack2) {
                    (Some(read), Some(al_pack)) => {
                        good_alignments.push((al_pack.clone(), read.clone()));
                        let _ = out_bam.write(&mut read.to_bam_record(&al_pack1, &Some(al_pack)));
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
        let pileup = align_helper.pileup(0, &good_alignments);

        quals = align_helper.base_quals(0, &pileup, rt_error);
        assert!(quals.len() == ref_seq.len());
    }

    fastq_writer.write_fmt(format_args!("@{}\n{}\n+\n{}\n", ref_seq_name, ref_seq,
                           fastq::get_qual_string(&quals, QUAL_OFFSET))).unwrap();
}

/// Extract aligner parameters from input arguments
fn get_scoring(args: &Args) -> sw::Scoring {
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
    // By convention all these scores are passed as positive values
    // The solver adds the appropriate signs
    assert!(gap_open >= 0);
    assert!(gap_extend >= 0);
    assert!(match_score >= 0);
    assert!(miss_score >= 0);
    assert!(clip >= 0);
    let scoring = sw::Scoring::from_scores(-gap_open,
                                           -gap_extend,
                                           match_score,
                                           -miss_score).
                                           xclip(-clip).
                                           yclip(0);
    scoring
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
    let scoring = get_scoring(&args).
        xclip(0); // always do local alignment

    let min_align_score = match args.flag_min_sw_score {
        Some(c) => c,
        None => 50.0,
    };

    let ref_filename = args.flag_ref;
    let r1_filename = args.flag_r1;
    let r2_filename = args.flag_r2;
    let outbam_filename = args.flag_outbam;
    let single_end = r2_filename.is_none();

    let fasta = bio::io::fasta::Reader::from_file(ref_filename)
        .expect("Failed to load reference FASTA");
    let mut ref_seq_names = Vec::new();
    let mut ref_seqs = Vec::new();
    for record_option in fasta.records() {
        let record = record_option.unwrap();
        let seq_name = record.id().to_owned();
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

    let align_helper = sw::AlignHelper::new(&ref_seqs, scoring.clone(),
                                            KMER_LEN_BANDED_ALIGN, WINDOW_SIZE_BANDED_ALIGN);

    let mut out_bam = bam::Writer::from_path(outbam_filename, &header)
        .expect("Failed to open output BAM for writing");

    let mut num_unmapped = 0;
    let mut npairs = 0;
    let mut fq_iter = fastq::CellrangerPairedFastqIter::new(&r1_filename,
                                                            r2_filename.as_ref().map(String::as_str),
                                                            args.flag_rev_strand);

    let cc_start = PreciseTime::now();
    loop {
        match fq_iter.next() {
            Some(pair) => {
                let bitenc1 = dna_string::DnaString::from_bytes(&pair.r1.seq);
                let bitenc2 = pair.r2.as_ref().map(|r| dna_string::DnaString::from_bytes(&r.seq));

                let mut read1 = graph_read::Read::new(pair.r1.id, pair.umi, pair.header.clone(), bitenc1, pair.r1.quals.clone());
                let mut read2 = bitenc2.map(|b| graph_read::Read::new(pair.r2.as_ref().unwrap().id,
                                                                      pair.umi,
                                                                      pair.header.clone(),
                                                                      b,
                                                                      pair.r2.as_ref().unwrap().quals.clone()));

                let read_match = align_helper.find_read_matches_sparse(&read1, min_align_score as i32);
                let mate_match = read2.as_ref().map(|r| align_helper.find_read_matches_sparse(&r, min_align_score as i32));
                npairs += 1;

                read1.set_first_in_template();
                if !single_end {
                    read1.set_paired();
                    read2.as_mut().unwrap().set_paired();
                }

                if args.flag_rev_strand {
                    read1.set_reverse();
                } else {
                    read2.as_mut().map(|r| r.set_reverse());
                }
                // Note: clone() below just clones a reference
                let mut rec = read1.to_bam_record(&read_match,
                                                  &mate_match.as_ref().and_then(|x| x.clone()));
                let mate_rec = read2.map(|r| r.to_bam_record(&mate_match.as_ref().and_then(|x| x.clone()),
                                                             &read_match));

                let _ = out_bam.write(&mut rec);
                if !single_end {
                    let _ = out_bam.write(&mut mate_rec.unwrap());
                }

                num_unmapped += match read_match {
                    Some(_) => 0,
                    None => 1
                };
                if !single_end {
                    num_unmapped += match mate_match {
                        Some(_) => 0,
                        None => 1
                    };
                }
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
            return;
        }
        let _ = fs::create_dir_all(outdir);
    }

    fn get_seqs_by_barcode(fasta: &str) -> HashMap<String, HashSet<String>> {

        let mut seqs = HashMap::new();
        println!("Opening {}", fasta);
        let reader = bio::io::fasta::Reader::new(fs::File::open(fasta).unwrap());
        for record in reader.records() {
            let rec = record.unwrap();
            let id = rec.id().to_string().split("_").collect::<Vec<&str>>()[0].to_string();
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
    fn test_vdj_asm_se() {

        let outdir = "test/outputs/asm";
        let true_fasta_name = "test/inputs/test_asm.fasta";
        let true_seqs = get_seqs_by_barcode(&true_fasta_name);

        let mut args = Args {
                        cmd_asm: true,
                        arg_inbam: Some("test/inputs/test_asm.bam".to_string()),
                        arg_outdir: Some(outdir.to_string()),
                        flag_kmers: Some(0),
                        flag_min_contig: Some(150),
                        flag_min_sw_score: Some(50.0),
                        flag_rt_error: Some(0.0001),
                        flag_min_qual: Some(2),
                        flag_qual_factor: Some(0.9),
                        flag_score_factor: Some(0.6),
                        flag_frac_reads: Some(0.3),
                        flag_min_umi_reads: Some(0),
                        flag_subsample_rate: Some(1.0),
                        .. Args::default()
                        };

        init_test(outdir);

        // Truly single end BAM
        args.arg_inbam = Some("test/inputs/test_asm_se.bam".to_string());
        args.flag_single_end = true;
        args.flag_min_contig = Some(150);

        vdj_asm(args.clone()); // use all reads

        let assembled_seqs = get_seqs_by_barcode(&"test/outputs/asm/test_asm_se.fasta");

        for (bc, seqs) in true_seqs.iter() {
            assert_eq!(*seqs, *assembled_seqs.get(bc).unwrap());
        }

        let bam = bam::Reader::from_path(Path::new(&"test/outputs/asm/test_asm_se.bam")).ok().expect("Error reading test_asm.bam");
        let (mapped, unmapped) = count_records(&bam);
        assert_eq!(mapped, 24);

        let header = bam.header();
        assert_eq!(header.target_count(), 6);
    }

    #[test]
    fn test_vdj_asm() {

        let outdir = "test/outputs/asm";
        let true_fasta_name = "test/inputs/test_asm.fasta";
        let true_seqs = get_seqs_by_barcode(&true_fasta_name);

        let mut args = Args {
                        cmd_asm: true,
                        arg_inbam: Some("test/inputs/test_asm.bam".to_string()),
                        arg_outdir: Some(outdir.to_string()),
                        flag_kmers: Some(0),
                        flag_min_contig: Some(150),
                        flag_min_sw_score: Some(50.0),
                        flag_rt_error: Some(0.0001),
                        flag_min_qual: Some(2),
                        flag_qual_factor: Some(0.9),
                        flag_score_factor: Some(0.6),
                        flag_frac_reads: Some(0.3),
                        flag_min_umi_reads: Some(1),
                        flag_use_unmapped: true,
                        .. Args::default()
                        };

        init_test(outdir);

        {
            vdj_asm(args.clone()); // use all read-pairs

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
            vdj_asm(args.clone());

            let assembled_seqs = get_seqs_by_barcode(&"test/outputs/asm/test_asm.fasta");

            for (bc, seqs) in assembled_seqs.iter() {
                assert_eq!(*seqs, *true_seqs.get(bc).unwrap());
            }

            let bam = bam::Reader::from_path(Path::new(&"test/outputs/asm/test_asm.bam")).ok().expect("Error reading test_asm.bam");

            let header = bam.header();
            assert_eq!(header.target_count(), 6);
        }
        {
            args.flag_min_contig = Some(250);

            vdj_asm(args.clone());

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

            vdj_asm(args.clone()); // use all reads

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

            vdj_asm(args.clone());

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
            args.flag_cons = true;
            args.flag_min_umi_reads = Some(1);
            args.flag_seed = Some(1);
            args.flag_min_sw_score = Some(50.0);
            vdj_asm(args.clone());

            let bam = bam::Reader::from_path(Path::new(&"test/outputs/asm/test_asm_se.bam")).ok().expect("Error reading test_asm.bam");
            let header = bam.header();
            assert_eq!(header.target_count(), 3);
        }
    }

    fn read_one_fq_record(fq_file: &str) -> (String, String, Vec<u8>) {
        let mut fastq = bio::io::fastq::Reader::from_file(fq_file).ok().expect("Could not open FASTQ file");
        let mut record = bio::io::fastq::Record::new();
        let _ = fastq.read(&mut record);
        let seq_name = record.id().to_owned();
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
        let ref_seq_name = record.id().to_owned();
        let ref_seq = String::from_utf8_lossy(record.seq()).into_owned();

        let args = Args {
                        cmd_base_quals: true,
                        arg_inpref: Some("test/inputs/base_quals/test_base_quals".to_string()),
                        arg_outdir: Some(outdir.to_string()),
                        flag_min_contig: Some(150),
                        flag_min_sw_score: Some(50.0),
                        flag_rt_error: Some(0.0001),
                        flag_min_qual: Some(2),
                        flag_qual_factor: Some(0.9),
                        flag_score_factor: Some(0.6),
                        flag_frac_reads: Some(0.3),
                        .. Args::default()
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
            let mut args3 = args.clone();
            args3.flag_single_end = true;
            get_base_qualities(args3);

            let bam = bam::Reader::from_path(Path::new(&out_bam)).ok().expect("Error reading test_base_quals.bam");
            let (mapped, unmapped) = count_records(&bam);
            assert_eq!(unmapped, 2);
            assert_eq!(mapped, 6);
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
            let _ = fq_writer.write(rec.id(), None, new_seq.as_bytes(), &new_qual);
        }
    }

    #[test]
    fn test_get_matches() {
        ()
        //let outdir = "test/outputs/base_quals";
        //let out_bam = "test/outputs/base_quals/test_match.bam";

        /*
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
                        flag_rev_strand: false,
                        flag_plot: false,
                        flag_plot_json: false,
                        flag_barcode: None,
                        flag_min_sw_score: Some(50.0),
                        flag_rt_error: Some(0.0001),
                        flag_min_qual: Some(2),
                        flag_qual_factor: Some(0.9),
                        flag_score_factor: Some(0.6),
                        flag_frac_reads: Some(0.3),
                        flag_min_umi_reads: Some(1),
                        flag_subsample_rate: Some(1.0),
                        flag_match_score: None,
                        flag_miss_score: None,
                        flag_gap_open: None,
                        flag_gap_extend: None,
                        flag_clip: None,
                        flag_seed: Some(20),
                        flag_single_end: false,
                        flag_use_unmapped: true,
                        };
        */

        // {
        //     get_matches(args.clone());

        //     let bam = bam::Reader::from_path(Path::new(&out_bam)).ok().expect("Error reading test_match.bam");
        //     let (mapped, unmapped) = count_records(&bam);
        //     // Two records left unmapped because they didn't have a good seed for mapping.
        //     assert_eq!(unmapped, 2);
        //     assert_eq!(mapped, 6);
        // }
    //     {
    //         // Test reverse complemented inputs
    //         revcomp_fastq("test/inputs/base_quals/test_base_quals_1.fastq",
    //                       "test/outputs/base_quals/test_base_quals_rc_1.fastq");
    //         revcomp_fastq("test/inputs/base_quals/test_base_quals_2.fastq",
    //                       "test/outputs/base_quals/test_base_quals_rc_2.fastq");
    //         let mut args5 = args.clone();
    //         args5.flag_rev_strand = true;
    //         args5.arg_fqpref = Some("test/outputs/base_quals/test_base_quals_rc".to_string());

    //         let bam = bam::Reader::from_path(Path::new(&out_bam)).ok().expect("Error reading test_match.bam");
    //         let (mapped, unmapped) = count_records(&bam);
    //         // Two records left unmapped because they didn't have a good seed for mapping.
    //         assert_eq!(unmapped, 2);
    //         assert_eq!(mapped, 6);
    //     }
    //     {
    //         let mut args2 = args.clone();
    //         // Input FASTA is all Ns.
    //         args2.arg_fasta = Some("test/inputs/base_quals/test_base_quals_invalid_bases_1.fasta".to_string());
    //         get_matches(args2);

    //         let bam = bam::Reader::from_path(Path::new(&out_bam)).ok().expect("Error reading test_match.bam");
    //         let (_, unmapped) = count_records(&bam);
    //         assert_eq!(unmapped, 8);
    //     }
    //     {
    //         // One all Ns, one good sequence. The presence of the invalid one shouldn't prevent us
    //         // from getting mathces on the good one.
    //         let mut args3 = args.clone();
    //         args3.arg_fasta = Some("test/inputs/base_quals/test_base_quals_invalid_bases_2.fasta".to_string());
    //         get_matches(args3);

    //         let bam = bam::Reader::from_path(Path::new(&out_bam)).ok().expect("Error reading test_match.bam");
    //         let (mapped, _) = count_records(&bam);
    //         assert_eq!(mapped, 6);
    //     }
    //     {
    //         let mut args4 = args.clone();
    //         args4.arg_fasta = Some("test/inputs/base_quals/test_base_quals_invalid_bases_2.fasta".to_string());
    //         // First pair has two reads with Ns. Second and third have only one.
    //         // 4th has both valid reads but no good seed. Last one should be good.
    //         args4.arg_fqpref = Some("test/inputs/base_quals/test_base_quals_invalid_bases".to_string());
    //         get_matches(args4);

    //         let bam = bam::Reader::from_path(Path::new(&out_bam)).ok().expect("Error reading test_match.bam");
    //         let (mapped, _) = count_records(&bam);
    //         assert_eq!(mapped, 8);
    //     }
     }
}
