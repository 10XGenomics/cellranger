//
// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
//

/// This stage does the following:
/// 1) Correct UMI sequences
/// 2) Mark BAM records as low-support (putatively chimeric) UMIs
/// 3) Mark BAM records as PCR duplicates
/// 4) Compute metrics related to the above

use std::fs;
use std::io;
use std::fs::File;
use std::path::Path;
use std::collections::{HashMap, HashSet};
use std::io::{BufWriter, Write};
use std::str;
use std::cmp::{min, max};
use std::default::Default;
use bincode;
use itertools::Itertools;
use rust_htslib::bam;
use tenkit::chunk_bam::{bam_block_offsets, chunk_bam_records, BamChunkIter, GzipHeader, GzipExtra, GzipISIZE, BamTellIter};
use martian::*;
use serde_json;

use bincode::internal::deserialize_from;
use bincode::Infinite;
use byteorder;

use utils;

pub struct MarkDuplicatesStage;

#[derive(Serialize, Deserialize, Default, Clone)]
/// Read accounting
struct Metrics {
    total_reads: u64,
    low_support_umi_reads: u64,
    umi_corrected_reads: u64,
    candidate_dup_reads: u64,
    dup_reads: u64,
    umis: u64,
}

#[derive(Serialize, Deserialize, Default, Clone)]
/// Per-barcode information
struct BarcodeSummary {
    barcodes: HashMap<Vec<u8>, BarcodeSummaryEntry>,
}
#[derive(Serialize, Deserialize, Default, Clone)]
struct BarcodeSummaryEntry {
    reads: u64,
    umis: u64,
    candidate_dup_reads: u64,
    umi_corrected_reads: u64,
}

#[derive(Deserialize, Serialize, Debug)]
pub struct LibraryInfo {
    library_id: String,
    library_type: String,
    gem_group: u64,
}

/// Within each gene, correct Hamming-distance-one UMIs
fn correct_umis(umigene_counts: HashMap<(Vec<u8>, String), u64>)
                -> HashMap<(Vec<u8>, String), Vec<u8>> {
    let nucs = b"ACGT";

    let mut corrections = HashMap::new();

    for (&(ref umi, ref gene_id), orig_count) in &umigene_counts {
        let mut test_umi = umi.clone();

        let mut best_dest_count = *orig_count;
        let mut best_dest_umi = umi.clone();

        for pos in 0..umi.len() {
            // Try each nucleotide at this position
            for test_char in nucs {
                if *test_char == umi[pos] {
                    // Skip the identitical nucleotide
                    continue;
                }
                test_umi[pos] = *test_char;

                // Test for the existence of this mutated UMI
                let test_count = *umigene_counts.get(&(test_umi.clone(), gene_id.clone())).unwrap_or(&0u64);

                // If there's a 1-HD UMI w/ greater count, move to that UMI.
                // If there's a 1-HD UMI w/ equal count, move to the lexicographically larger UMI.
                if test_count > best_dest_count || (test_count == best_dest_count && test_umi > best_dest_umi) {
                    best_dest_umi = test_umi.clone();
                    best_dest_count = test_count;
                }
            }
            // Reset this position to the unmutated sequence
            test_umi[pos] = umi[pos];
        }
        if *umi != best_dest_umi {
            corrections.insert((umi.clone(), gene_id.clone()), best_dest_umi);
        }
    }
    corrections
}

#[derive(Deserialize, Debug)]
struct ChunkArgs {
    input: String,
    chunk_start: i64,
    chunk_end: Option<i64>,
    filter_umis: bool,
    library_info: Vec<LibraryInfo>,
}

#[derive(Deserialize, Debug)]
struct StageArgs {
    library_info: Vec<LibraryInfo>,
}

/// Denotes whether a record is unpaired, first of a read-pair, or second of a read-pair
#[derive(Debug, PartialEq, Hash, Eq)]
enum MateType {
    Unpaired,
    First,
    Second,
}

fn get_mate_type(record: &bam::Record) -> MateType {
    match (record.is_paired(), record.is_first_in_template(), record.is_last_in_template()) {
        (false, false, false) => Some(MateType::Unpaired),
        (true, true, false) => Some(MateType::First),
        (true, false, true) => Some(MateType::Second),
        _ => None,
    }.expect("Invalid combination of is_paired, is_first, and is_second.")
}

/// For all the records w/ the same qname (a read or a read-pair),
/// find the single gene/feature confidently mapped to if any.
fn get_qname_conf_mapped_feature<'a, I>(records: I) -> Option<String>
    where I: Iterator<Item = &'a bam::record::Record> {
    let mut features = Vec::new();
    for record in records.filter(|&r| !r.is_secondary()) {
        let flags = utils::get_read_extra_flags(record);

        if flags.intersects(utils::ExtraFlags::CONF_MAPPED) ||
            flags.intersects(utils::ExtraFlags::CONF_FEATURE) {
                features.extend(utils::get_read_feature_ids(record).into_iter());
            }
    }
    features.sort_unstable();
    features.dedup();
    match features.len() {
        1 => features.into_iter().nth(0),
        _ => None,
    }
}

/// Estimate BAM file compression ratio.
/// Iterate first N BGZF blocks, get the BSIZE (block SIZE) and ISIZE (size of uncompressed
/// data) of block, estimate the compression ratio as the average of (ISIZE/BSIZE).
///
/// Note: this function shares significant amount of code with parse_bgzf_header function.
fn bam_compression_ratio<P: AsRef<Path>>(bam_path : &P) -> f64{
    use self::io::Seek;

    let bam_bytes = fs::metadata(bam_path)
        .expect("Failed to get filesize of input BAM").len();
    let mut bgzf = File::open(bam_path).expect("Failed to open input BAM");
    bgzf.seek(io::SeekFrom::Start(0)).expect("Failed to get the beginning of the BAM file ");

    let mut i = 0;
    let mut compression_ratios : Vec<f64> = Vec::new();
    while i < 1000 {
        i += 1;

        let cur_pos = bgzf.seek(io::SeekFrom::Current(0))
            .expect("Failed to get current position in BAM file");

        if cur_pos == bam_bytes {
            break;
        }

        let header: GzipHeader = deserialize_from::<_, _, _, byteorder::LittleEndian>(&mut bgzf, Infinite)
            .expect("Failed to deserialize gzip header - invalid BAM file");

        if header.id1 != 31 || header.id2 != 139 {
            panic!("Invalid gzip header in BAM file");
        }

        // Determine the BGZF block size
        let extra: GzipExtra = deserialize_from::<_, _, _, byteorder::LittleEndian>(&mut bgzf, Infinite)
            .expect("Failed to deserialize gzip extra field - invalid BAM file");
        if extra.id1 != 66 || extra.id2 != 67 || extra.slen != 2 {
            panic!("BSIZE field not found in BGZF header - invalid BAM file");
        }

        // Determine the size of uncompressed data
        let isize_pos = cur_pos + (extra.bsize as u64) + 1 - 4;
        bgzf.seek(io::SeekFrom::Start(isize_pos))
            .expect("Failed to get isize position in BAM file");
        let isize: GzipISIZE = deserialize_from::<_, _, _, byteorder::LittleEndian>(&mut bgzf, Infinite)
            .expect("Failed to deserialize gzip isize field - invalid BAM file");

        if extra.bsize > 0 && isize.isize > extra.bsize as u32 {
            compression_ratios.push(extra.bsize as f64 * 1.0 / (min(isize.isize, 2u32.pow(16)) as f64));
        }
    }

    if compression_ratios.len() == 0 {
        0.2
    } else {
        compression_ratios.iter().cloned().fold(0./0., f64::min)
    }
}

fn major_barcode_proportion_within_chunk(bam: &mut bam::Reader,
                                         block_offsets : &Vec<u64>,
                                         left_block_idx : usize,
                                         right_block_idx : usize,
                                         num_sample_per_chunk : i32) -> f64 {
    use self::bam::Read;
    use std::collections::btree_map::BTreeMap;

    let step = max(((right_block_idx - left_block_idx) as i32 / num_sample_per_chunk) as usize, 1);
    let mut idx = left_block_idx;
    let mut barcode_counter : BTreeMap<Vec<u8>, u32> = BTreeMap::new();

    while idx < right_block_idx {
        bam.seek((block_offsets[idx] << 16) as i64)
        .expect("Failed to get block position in BAM file");

        let mut bam_iter = BamTellIter { bam: bam };
        if let Some(Ok((_, record))) = bam_iter.next() {
            if let Some(bc) = utils::get_read_barcode(&record) {
                *barcode_counter.entry(bc).or_insert(0) += 1;
            }
        }

        idx += step;
    };

    // calculate major barcode proportion
    let (mut total_bc_count, mut major_bc_count) = (0u32, 0u32);
    for (_, count) in &barcode_counter {
        total_bc_count += count;
        major_bc_count = max(*count, major_bc_count);
    };

    if total_bc_count == 0 {
        1f64
    } else {
        (major_bc_count as f64) / (total_bc_count as f64)
    }
}

/// Do duplicate marking on a single barcode's worth of data
fn process_barcode(bc_group: &mut Vec<bam::Record>,
                   metrics: &mut Vec<Metrics>,
                   bc_summaries: &mut Vec<BarcodeSummary>,
                   out_bam: &mut bam::Writer,
                   lib_idx: usize,
                   bc: &Vec<u8>,
                   filter_umis: bool) {
    // HACK: Rust's entry api still requires key-ownership to query,
    //       so do it the ugly way.
    if bc_summaries[lib_idx].barcodes.get(bc).is_none() {
        // Only clone the key when it's not present
        bc_summaries[lib_idx].barcodes.insert(bc.clone(), Default::default());
    }
    let barcode_summary = bc_summaries[lib_idx].barcodes.get_mut(bc).unwrap();

    // Get raw UMI frequencies
    let mut raw_umigene_counts: HashMap<(Vec<u8>, String), u64> = HashMap::new();
    let mut low_support_umigenes: HashSet<(Vec<u8>, String)> = HashSet::new();

    for (maybe_umi, umi_group) in bc_group.iter().group_by(|x| utils::get_read_umi(&x)) {
        let mut gene_counts: HashMap<String, u64> = HashMap::new();

        if let Some(umi) = maybe_umi {
            // Count (raw UMI, feature) pairs to prepare for UMI correction.
            // Track (raw UMI, feature) pairs with submaximal count per (raw UMI)
            //   to prepare for marking of low-support UMIs (putative chimeras).

            // Assumes records are qname-contiguous in the input.
            for (qname, qname_records) in umi_group.into_iter()
                .filter(|&x| utils::is_read_dup_candidate(x))
                .group_by(|x| x.qname()) {
                    match get_qname_conf_mapped_feature(qname_records.into_iter()) {
                        None => { panic!(format!("Found 0 or >1 features for confidently mapped read/pair {}", str::from_utf8(qname).unwrap())) },
                        Some(gene) => {
                            *raw_umigene_counts.entry((umi.clone(), gene.clone())).or_insert(0) += 1;
                            *gene_counts.entry(gene.clone()).or_insert(0) += 1;
                        },
                    }
                } // for each qname

            // Mark (UMI, gene) pairs w/ frequency below the max for the UMI as low support.
            if let Some((_max_gene, max_count)) = gene_counts.iter().max_by_key(|x| x.1) {
                for (gene, count) in gene_counts.iter() {
                    if filter_umis && count < max_count {
                        low_support_umigenes.insert((umi.clone(), gene.clone()));
                    }
                }
            }
        }
    }

    // Determine which UMIs need to be corrected
    let umi_corrections = correct_umis(raw_umigene_counts);

    // Correct UMIs and mark PCR duplicates
    let mut wrote_umigenes = HashSet::new();

    for (_qname, qname_records) in bc_group.iter().group_by(|x| x.qname()) {
        // Take UMI from first record
        let maybe_umi = utils::get_read_umi(&qname_records.iter().nth(0).unwrap());

        let maybe_gene = get_qname_conf_mapped_feature(qname_records.iter()
                                                       .map(|x| *x));

        for record in qname_records {
            // Clone here because holding a &mut to the record
            // while simultaneously grouping by qname in the outer loop
            // makes rustc very unhappy.
            let mut new_record = record.clone();
            let mut new_flags = utils::get_read_extra_flags(&new_record);

            let is_secondary = record.is_secondary();
            let is_dup_candidate = utils::is_read_dup_candidate(&record);

            if !is_secondary {
                metrics[lib_idx].total_reads += 1;
                barcode_summary.reads += 1
            }

            metrics[lib_idx].candidate_dup_reads += is_dup_candidate as u64;

            if let (Some(umi), Some(gene)) = (maybe_umi.as_ref(), maybe_gene.as_ref()) {
                let key = (umi.clone(), gene.clone());

                if !record.is_secondary() && low_support_umigenes.contains(&key) {
                    // Low support (UMI, gene). Mark as low support.
                    // - Only consider primary alignments for this flag.
                    // - Do not correct the UMI.
                    // - Do not mark duplicates w/ for this (UMI, gene).
                    metrics[lib_idx].low_support_umi_reads += 1;
                    new_flags |= utils::ExtraFlags::LOW_SUPPORT_UMI;

                } else {
                    // Correct UMIs in all records
                    let (corrected_umi, is_corrected) = match umi_corrections.get(&key) {
                        Some(new_umi) => (new_umi.clone(), true),
                        None => (umi.clone(), false),
                    };

                    // Correct the UMI tag
                    if is_corrected {
                        metrics[lib_idx].umi_corrected_reads += 1;
                        barcode_summary.umi_corrected_reads += 1;

                        new_record.remove_aux(utils::PROC_UMI_SEQ_TAG);
                        new_record.push_aux(utils::PROC_UMI_SEQ_TAG,
                                            &bam::record::Aux::String(&corrected_umi));
                    }

                    // Don't try to dup mark secondary alignments.
                    if is_dup_candidate {
                        let dup_key = (corrected_umi, key.1, get_mate_type(&record));

                        barcode_summary.candidate_dup_reads += 1;

                        if wrote_umigenes.contains(&dup_key) {
                            // Duplicate
                            metrics[lib_idx].dup_reads += 1;
                            let flags = record.flags();
                            new_record.set_flags(flags | 1024u16);
                        } else {
                            // Non-duplicate
                            wrote_umigenes.insert(dup_key);

                            // Flag read1 as countable
                            if !record.is_last_in_template() {
                                metrics[lib_idx].umis += 1;
                                barcode_summary.umis += 1;
                                new_flags |= utils::ExtraFlags::UMI_COUNT;
                            }
                        }
                    }
                }
            }
            if new_flags.bits() > 0 {
                new_record.remove_aux(utils::EXTRA_FLAGS_TAG);
                new_record.push_aux(utils::EXTRA_FLAGS_TAG,
                                    &bam::record::Aux::Integer(new_flags.bits() as i32));
            }
            out_bam.write(&new_record).expect("Failed to write BAM record");
        }
    }
}

fn cmd_mark_dups(args: &ChunkArgs, outs: JsonDict) -> JsonDict {
    use bam::Read;

    // Load library info
    let library_info = &args.library_info;

    // Partition metrics by library
    let mut metrics: Vec<Metrics> = vec![Default::default(); library_info.len()];
    let mut bc_summaries: Vec<BarcodeSummary> = vec![Default::default(); library_info.len()];

    let mut bam = bam::Reader::from_path(&args.input)
        .expect("Failed to open input BAM file");

    let mut out_bam = bam::Writer::from_path(outs["alignments"].as_str().unwrap(),
                                             &bam::Header::from_template(bam.header()))
        .expect("Failed to open output BAM file");

    let range = (args.chunk_start, args.chunk_end);
    let chunk_iter = BamChunkIter::new(&mut bam, range);

    let mut maybe_prev_bc: Option<Vec<u8>> = None;
    let mut maybe_prev_lib_idx: Option<usize> = None;

    // Records for a single group
    let mut group_records: Vec<bam::Record> = vec![];

    // Group records by (barcode, library) and process each barcode.
    // The following code previously used itertools group_by().
    // For the pathological case of a large dataset
    // with a typical number of invalid barcodes,
    // it would load an entire chunk of invalid-barcode records
    // into memory unnecessarily. Chaining filter() upstream did not
    // work due to the required side effects of metric-tallying and BAM writing.
    // So here we do filter(...).group_by(...) from scratch.
    for record in chunk_iter.map(|r| r.unwrap()) {
        let maybe_bc = utils::get_read_barcode(&record);
        let lib_idx = utils::get_read_library_index(&record);

        // Skip records without a valid barcode.
        if maybe_bc.is_none() {
            metrics[lib_idx].total_reads += 1;
            out_bam.write(&record).expect("Failed to write BAM record");
            continue;
        }

        let bc = maybe_bc.unwrap();

        if let (Some(ref prev_bc), Some(ref prev_lib_idx)) = (maybe_prev_bc, maybe_prev_lib_idx) {
            // Verify sort order (bc, library_idx)
            // Within a barcode, the library index is non-decreasing
            assert!(bc != *prev_bc || lib_idx >= *prev_lib_idx);

            if bc != *prev_bc || lib_idx != *prev_lib_idx {
                // Hit a group key boundary
                process_barcode(&mut group_records,
                                &mut metrics, &mut bc_summaries, &mut out_bam,
                                *prev_lib_idx, prev_bc, args.filter_umis);
                group_records.clear()
            }
        }

        // Accumulate the records in this group
        group_records.push(record);

        maybe_prev_bc = Some(bc);
        maybe_prev_lib_idx = Some(lib_idx);
    } // for each record

    // Process the final group
    if !group_records.is_empty() {
        process_barcode(&mut group_records,
                        &mut metrics, &mut bc_summaries, &mut out_bam,
                        maybe_prev_lib_idx.unwrap(), &maybe_prev_bc.unwrap(),
                        args.filter_umis);
    }

    // Write metrics
    utils::write_json_file(outs["metrics"].as_str().unwrap(), &metrics)
        .expect("Failed to write metrics to JSON file");

    // Write barcode summary
    let chunk_bc_file = File::create(outs["chunk_barcode_summary"].as_str().unwrap())
        .expect("Failed to open barcode summary file for writing");
    let mut writer = BufWriter::new(chunk_bc_file);
    bincode::serialize_into(&mut writer, &bc_summaries, bincode::Infinite)
        .expect("Failed to serialize barcode summary to file");

    outs
}

impl MartianStage for MarkDuplicatesStage {
    fn split(&self, args: JsonDict) -> JsonDict {
        let bam = args["input"].as_str().unwrap();

        let bam_comp_ratio = bam_compression_ratio(&bam);
        let mut bam_reader = bam::Reader::from_path(bam).expect("Failed to open input BAM");

        let block_offsets = bam_block_offsets(&bam);
        let chunk_intervals =
            chunk_bam_records(&bam, &block_offsets, &|x| utils::get_read_barcode(&x), 0.5, 256);

        let mut left_idx = 0usize;
        let mut chunks : Vec<serde_json::Value> = Vec::new();

        for offset in &chunk_intervals {
            let right_idx = match offset.1 {
                Some(end) => block_offsets.binary_search(&((end >> 16) as u64))
                .expect("BAM file error - chunk offset is not block offset"),
                _ => block_offsets.len() as usize
            };

            let major_bc_prop = major_barcode_proportion_within_chunk(&mut bam_reader, &block_offsets, left_idx, right_idx, 1000);
            let chunk_size_gb = ((block_offsets[right_idx - 1usize] - block_offsets[left_idx]) as f64) / bam_comp_ratio / (1024u64.pow(3) as f64);
            // emprically learned model for memory usage estimation
            let mem_gb = 2f64 + chunk_size_gb + chunk_size_gb * major_bc_prop;
            let mem_gb_round = min(64i32, max(2i32, ((mem_gb / 2.0f64).ceil() * 2.0f64) as i32)); // round to next even

            chunks.push(json!({
                "chunk_start": offset.0,
                "chunk_end": offset.1,
                "chunk_size_gb" : chunk_size_gb,
                "bam_compression_ratio" : bam_comp_ratio,
                "major_barcode_proportion" : major_bc_prop,
                "__mem_gb": mem_gb_round,
            }));

            left_idx = right_idx;
        }

        serde_json::from_value(json!({
            "chunks": json!(chunks),
        })).expect("Failed to serialize split outs")
    }

    fn main(&self, json_args: JsonDict, outs: JsonDict) -> JsonDict {
        let args: ChunkArgs = obj_decode(json_args);
        cmd_mark_dups(&args, outs)
    }

    fn join(&self, json_args: JsonDict, outs: JsonDict, _chunk_defs: Vec<JsonDict>, chunk_outs: Vec<JsonDict>) -> JsonDict {
        let args: StageArgs = obj_decode(json_args);

        let mut final_outs = outs.clone();

        // Collect output BAMs
        final_outs["output"] = json!(
            chunk_outs.iter()
            .map(|x| x["alignments"].clone())
                .collect::<Vec<serde_json::Value>>());

        // Load library info
        let library_info = &args.library_info;

        // Merge summary metrics, grouping by library_type
        let mut metrics: HashMap<String, Metrics> = HashMap::new();
        for lib in library_info {
            metrics.insert(lib.library_type.clone(), Default::default());
        }

        for chunk_out in chunk_outs.iter() {
            let chunk: Vec<Metrics> = utils::read_json_file(chunk_out["metrics"].as_str().unwrap())
                .expect("Failed to read from metrics JSON file");
            for (lib_idx, lib) in chunk.iter().enumerate() {
                let lt = &library_info[lib_idx].library_type;
                let mut m = metrics.get_mut(lt).unwrap();

                m.total_reads += lib.total_reads;
                m.umi_corrected_reads += lib.umi_corrected_reads;
                m.dup_reads += lib.dup_reads;
                m.umis += lib.umis;
                m.candidate_dup_reads += lib.candidate_dup_reads;
                m.low_support_umi_reads += lib.low_support_umi_reads;
            }
        }

        // Merge per-barcode data, grouping by library_type
        let mut bc_summaries: HashMap<String, BarcodeSummary> = HashMap::new();
        for lib in library_info {
            bc_summaries.insert(lib.library_type.clone(), Default::default());
        }

        for chunk_out in chunk_outs.iter() {
            let mut reader = File::open(chunk_out["chunk_barcode_summary"].as_str().unwrap())
                .expect("Failed to open barcode summary binary file for reading");

            let chunk: Vec<BarcodeSummary> =
                bincode::deserialize_from(&mut reader, bincode::Infinite)
                .expect("Failed to deserialize from barcode summary bincode file");

            for (lib_idx, lib) in chunk.iter().enumerate() {
                let lt = &library_info[lib_idx].library_type;
                let mut bc_summaries_lt = bc_summaries.get_mut(lt).unwrap();

                // Merge the per-chunk barcode summaries into a single barcode summary
                // for this library type.
                for (bc, bc_entry) in &lib.barcodes {
                    if bc_summaries_lt.barcodes.get(bc).is_none() {
                        bc_summaries_lt.barcodes.insert(bc.clone(),
                                                        Default::default());
                    }
                    let new_entry = bc_summaries_lt.barcodes.get_mut(bc).unwrap();

                    new_entry.reads += bc_entry.reads;
                    new_entry.umis += bc_entry.umis;
                    new_entry.candidate_dup_reads += bc_entry.candidate_dup_reads;
                    new_entry.umi_corrected_reads += bc_entry.umi_corrected_reads;
                }
            }
        }

        // Write per-barcode data
        let bc_summary_file = File::create(final_outs["barcode_summary"].as_str().unwrap())
            .expect("Failed to open barcode summary file for writing");
        let mut writer = BufWriter::new(bc_summary_file);
        write!(writer,
               "library_type,barcode,reads,umis,candidate_dup_reads,umi_corrected_reads\n")
            .expect("Failed to write to barcode summary file");
        for (lib_type, bc_summary) in bc_summaries {
            for (bc, entry) in bc_summary.barcodes {
                write!(writer, "\"{}\",{},{},{},{},{}\n",
                       lib_type, str::from_utf8(&bc).unwrap(),
                       entry.reads,
                       entry.umis, entry.candidate_dup_reads,
                       entry.umi_corrected_reads)
                    .expect("Failed to write to barcode summary file");
            }
        }
        drop(writer);


        // Write summary metrics
        let mut summary = serde_json::Map::new();
        for (lib_type, lt_metrics) in &metrics {
            let prefix = utils::get_library_type_metric_prefix(lib_type);
            summary.insert(format!("{}corrected_umi_frac", prefix), json!(
                           lt_metrics.umi_corrected_reads as f64 /
                           lt_metrics.total_reads as f64));
            summary.insert(format!("{}low_support_umi_reads_frac", prefix), json!(
                           lt_metrics.low_support_umi_reads as f64 /
                           lt_metrics.candidate_dup_reads as f64));
            summary.insert(format!("{}multi_cdna_pcr_dupe_reads_frac", prefix), json!(
                           lt_metrics.dup_reads as f64 /
                           lt_metrics.candidate_dup_reads as f64));
        }

        utils::write_json_file(final_outs["summary"].as_str().unwrap(), &summary)
            .expect("Failed to write to summary JSON file");

        final_outs
    }
}

#[cfg(test)]
mod test {
    use std::fs;
    use std::collections::HashMap;
    use std::path::{Path, PathBuf};
    use std::str;
    use rust_htslib::bam;
    use rust_htslib::bam::HeaderView;
    use rust_htslib::bam::header::{Header, HeaderRecord};
    use rust_htslib::bam::Read;
    use tenkit::chunk_bam::{bam_block_offsets, chunk_bam_records};
    use serde_json;

    use cmd_mark_dups;
    use cmd_mark_dups::correct_umis;
    use utils;

    const RAW_UMI_SEQ_TAG: &'static [u8]     = b"UR";

    #[test]
    fn test_correct_umis() {
        let mut umis: HashMap<(Vec<u8>, String), u64> = HashMap::new();
        umis.insert((b"AAAA".to_vec(), "G0".to_owned()), 3u64);
        umis.insert((b"AAAT".to_vec(), "G0".to_owned()), 2u64);
        umis.insert((b"AAAA".to_vec(), "G1".to_owned()), 1u64);
        umis.insert((b"AATT".to_vec(), "G1".to_owned()), 1u64);

        let corr = correct_umis(umis);
        assert!(corr.len() == 1);
        assert!(corr.get(&(b"AAAT".to_vec(), "G0".to_owned())).unwrap() == b"AAAA");

        let mut umis: HashMap<(Vec<u8>, String), u64> = HashMap::new();
        umis.insert((b"CCCC".to_vec(), "G0".to_owned()), 1u64);
        umis.insert((b"CGCC".to_vec(), "G0".to_owned()), 1u64);

        let corr = correct_umis(umis);
        assert!(corr.len() == 1);
        assert!(corr.get(&(b"CCCC".to_vec(), "G0".to_owned())).unwrap() == b"CGCC");
    }

    fn sam_to_records(header: &HeaderView, sam: &[u8]) -> Vec<bam::Record> {
        sam
            .split(|x| *x == b'\n')
            .filter(|x| x.len() > 0 && x[0] != b'@' && x[0] != b'#') // strip headers and test comments
            .map(|line| {
                // Split qname from read_id:test_flags
                let qname = str::from_utf8(line).unwrap().split('\t').nth(0).unwrap();
                let parts: Vec<&str> = qname.split(':').collect();
                let new_line = match parts[1].len() > 0 {
                    true => format!("{}\ttf:Z:{}", str::from_utf8(line).unwrap(), parts[1]),
                    false => str::from_utf8(line.clone()).unwrap().to_owned(),
                };
                let mut rec = bam::Record::from_sam(header, new_line.as_bytes()).unwrap();
                rec.set_qname(parts[0].as_bytes());
                rec
            })
            .collect()
    }

    fn get_test_flags(record: &bam::Record) -> String {
        str::from_utf8(record.aux(b"tf").map_or(b"", |x| x.string())).unwrap().to_owned()
    }

    /// Verify a dup-marked BAM and metrics using text in the QNAME that encodes
    /// the expected properties of each dup-marked record
    fn check_outputs<P: AsRef<Path>>(test_records: &Vec<bam::Record>,
                                     bam_filename: P, metrics_filename: P) -> () {
        let mut bam = bam::Reader::from_path(&bam_filename).expect("Failed to read BAM file");

        // Compute expected metric values
        let mut expected_metrics: cmd_mark_dups::Metrics = Default::default();
        for record in test_records {
            let is_secondary = (record.flags() & 256u16) > 0;
            if !is_secondary {
                expected_metrics.total_reads += 1;
            }

            let test_flags = get_test_flags(&record);
            // d => duplicate
            expected_metrics.dup_reads += test_flags.contains("d") as u64;

            // L => low-support UMI
            expected_metrics.low_support_umi_reads += test_flags.contains("L") as u64;

            // u => UMI corrected
            expected_metrics.umi_corrected_reads += test_flags.contains("u") as u64;

            // c => candidate dup read (has all the tags and is confidently mapped)
            expected_metrics.candidate_dup_reads += test_flags.contains("c") as u64;

            // C => countable read
            expected_metrics.umis += test_flags.contains("C") as u64;
        }

        let lib_metrics: Vec<cmd_mark_dups::Metrics> = utils::read_json_file(&metrics_filename)
            .expect("Failed to read metrics file");
        let metrics = &lib_metrics[0];
        assert_eq!(metrics.total_reads,           expected_metrics.total_reads);
        assert_eq!(metrics.low_support_umi_reads, expected_metrics.low_support_umi_reads);
        assert_eq!(metrics.umi_corrected_reads,   expected_metrics.umi_corrected_reads);
        assert_eq!(metrics.candidate_dup_reads,   expected_metrics.candidate_dup_reads);
        assert_eq!(metrics.dup_reads,             expected_metrics.dup_reads);
        assert_eq!(metrics.umis,                  expected_metrics.umis);


        // Check BAM records
        for record in bam.records().map(|x| x.unwrap()) {
            let test_flags = get_test_flags(&record);
            // d => duplicate
            let expect_dup = test_flags.contains("d");

            // L => low-support UMI
            let expect_low_supp = test_flags.contains("L");

            // u => UMI corrected
            let expect_umi_corr = test_flags.contains("u");

            let expect_countable = test_flags.contains("C");

            assert_eq!(expect_dup, record.is_duplicate());

            let is_low_supp = utils::get_read_extra_flags(&record).intersects(utils::ExtraFlags::LOW_SUPPORT_UMI);
            assert_eq!(expect_low_supp, is_low_supp);

            let is_countable = utils::get_read_extra_flags(&record).intersects(utils::ExtraFlags::UMI_COUNT);
            assert_eq!(expect_countable, is_countable);

            let maybe_raw_umi = record.aux(RAW_UMI_SEQ_TAG);
            let maybe_proc_umi = record.aux(utils::PROC_UMI_SEQ_TAG);
            let umi_corrected = match (maybe_raw_umi, maybe_proc_umi) {
                (Some(raw_umi), Some(proc_umi)) => raw_umi.string() != proc_umi.string(),
                _ => false,
            };
            assert_eq!(expect_umi_corr, umi_corrected);
        }
    }

    fn test_mark_dups<P: AsRef<Path>>(sam_text: &[u8], test_subdir: P) {
        use super::LibraryInfo;
        // Setup test dir
        let test_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("test_output").join("mark_dups").join(test_subdir);
        fs::create_dir_all(&test_dir)
            .expect("Failed to create test output dir");

        // Write input SAM to BAM file
        // Generate input SAM
        let mut header = Header::new();
        header.push_record(HeaderRecord::new(b"SQ").push_tag(b"SN", &"1").push_tag(b"LN", &10000000));

        let in_bam_filename = test_dir.join("input.bam");

        let test_records = sam_to_records(&HeaderView::from_header(&header), sam_text);
        {
            let mut writer = bam::Writer::from_path(&in_bam_filename, &header)
                .expect("Failed to create BAM writer");
            for record in &test_records {
                writer.write(&record).expect("Failed to write BAM record");
            }
        }

        // Setup input
        let block_offsets = bam_block_offsets(&in_bam_filename);
        let chunk_intervals =
            chunk_bam_records(&in_bam_filename, &block_offsets, &|x| utils::get_read_barcode(&x), 0.5, 256);

        let lib_info: Vec<LibraryInfo> = serde_json::from_value(
            json!([
                {
                    "library_type": "Gene Expression",
                    "library_id": "0",
                    "gem_group": 1,
                }])).unwrap();

        let chunk_args = cmd_mark_dups::ChunkArgs {
            input: in_bam_filename.to_str().unwrap().to_owned(),
            chunk_start: chunk_intervals[0].0,
            chunk_end: chunk_intervals[0].1,
            filter_umis: true,
            library_info: lib_info,
        };

        // Setup output
        let metrics_filename = test_dir.join("metrics.json");
        let out_bam_filename = test_dir.join("output.bam");
        let bc_summary_filename = test_dir.join("chunk_barcode_summary.bin");
        let chunk_outs = match json!({
            "metrics": metrics_filename.clone(),
            "alignments": out_bam_filename.clone(),
            "chunk_barcode_summary": bc_summary_filename.clone(),
        }) {
            serde_json::Value::Object(x) => Some(x),
            _ => None,
        }.unwrap();

        let outs = cmd_mark_dups::cmd_mark_dups(&chunk_args, chunk_outs);

        // Check outs dict
        assert_eq!(outs["metrics"].as_str().unwrap(), metrics_filename.to_str().unwrap());
        assert_eq!(outs["alignments"].as_str().unwrap(), out_bam_filename.to_str().unwrap());

        // Check output files
        check_outputs(&test_records, &out_bam_filename, &metrics_filename);
    }


    /// In the below tests, the qnames are read_id:<test-flags> where test-flags are the expected properties of the records.
    /// See test code above for test-flag defs.
    #[test]
    fn test_mark_dups_se() {
        let sam_text = b"
r0:cC	0	1	1	255	1M	*	0	0	A	F	CB:Z:AAAA-1	UR:Z:AAAA	UB:Z:AAAA	GX:Z:G1	xf:i:1	fx:Z:G1	li:i:0
r1:dc	0	1	1	255	1M	*	0	0	A	F	CB:Z:AAAA-1	UR:Z:AAAA	UB:Z:AAAA	GX:Z:G1	xf:i:1	fx:Z:G1	li:i:0
r1:	256	1	2	0	1M	*	0	0	A	F	CB:Z:AAAA-1	UR:Z:AAAA	UB:Z:AAAA	GX:Z:G1	fx:Z:G1	li:i:0
r2:Lc	0	1	1	255	1M	*	0	0	A	F	CB:Z:AAAA-1	UR:Z:AAAA	UB:Z:AAAA	GX:Z:G2	xf:i:1	fx:Z:G2	li:i:0
r3:udc	0	1	1	255	1M	*	0	0	A	F	CB:Z:AAAA-1	UR:Z:AATA	UB:Z:AATA	GX:Z:G1	xf:i:1	fx:Z:G1	li:i:0
r4:	0	1	1	255	1M	*	0	0	A	F	CB:Z:AAAA-1	UR:Z:AATA	UB:Z:AATA	li:i:0
r5:	0	1	1	255	1M	*	0	0	A	F	CB:Z:AAAA-1	UR:Z:AATA	UB:Z:AATA	GX:Z:G1;G2	fx:Z:G1;G2	li:i:0
r6:	0	1	1	0	1M	*	0	0	A	F	CB:Z:AAAA-1	UR:Z:AATA	UB:Z:AATA	GX:Z:G1	fx:Z:G1	li:i:0
r12:	0	1	1	255	1M	*	0	0	A	F	CB:Z:AAAA-1	UR:Z:AATA	GX:Z:G1	fx:Z:G1	li:i:0
r7:cC	0	1	1	255	1M	*	0	0	A	F	CB:Z:AAAA-1	UR:Z:CCCC	UB:Z:CCCC	GX:Z:G1	xf:i:1	fx:Z:G1	li:i:0
r8:dc	0	1	1	255	1M	*	0	0	A	F	CB:Z:AAAA-1	UR:Z:CCCC	UB:Z:CCCC	GX:Z:G1	xf:i:1	fx:Z:G1	li:i:0
r9:cC	0	1	1	255	1M	*	0	0	A	F	CB:Z:AATT-1	UR:Z:AAAA	UB:Z:AAAA	GX:Z:G1	xf:i:1	fx:Z:G1	li:i:0
r10:dc	0	1	1	255	1M	*	0	0	A	F	CB:Z:AATT-1	UR:Z:AAAA	UB:Z:AAAA	GX:Z:G1	xf:i:1	fx:Z:G1	li:i:0
r11:	0	1	1	255	1M	*	0	0	A	F	UR:Z:AAAA	UB:Z:AAAA	GX:Z:G3	xf:i:1	fx:Z:G3	li:i:0
";
        test_mark_dups(sam_text, &"se");
    }

    #[test]
    fn test_mark_dups_pe() {
        let sam_text = b"
# Non-duplicate
r0:cC	67	1	1	255	1M	*	0	0	A	F	CB:Z:AAAA-1	UR:Z:AAAA	UB:Z:AAAA	GX:Z:G1	xf:i:1	fx:Z:G1	li:i:0
r0:c	131	1	1	255	1M	*	0	0	T	F	CB:Z:AAAA-1	UR:Z:AAAA	UB:Z:AAAA	GX:Z:G1	xf:i:1	fx:Z:G1	li:i:0
# Duplicate
r1:dc	67	1	1	255	1M	*	0	0	A	F	CB:Z:AAAA-1	UR:Z:AAAA	UB:Z:AAAA	GX:Z:G1	xf:i:1	fx:Z:G1	li:i:0
r1:dc	131	1	1	255	1M	*	0	0	T	F	CB:Z:AAAA-1	UR:Z:AAAA	UB:Z:AAAA	GX:Z:G1	xf:i:1	fx:Z:G1	li:i:0
# Secondary alignment for duplicate
r1:	323	1	2	0	1M	*	0	0	A	F	CB:Z:AAAA-1	UR:Z:AAAA	UB:Z:AAAA	GX:Z:G1	fx:Z:G1	li:i:0
r1:	387	1	2	0	1M	*	0	0	T	F	CB:Z:AAAA-1	UR:Z:AAAA	UB:Z:AAAA	GX:Z:G1	fx:Z:G1	li:i:0
# Low-support gene for this UMI
r2:Lc	67	1	1	255	1M	*	0	0	A	F	CB:Z:AAAA-1	UR:Z:AAAA	UB:Z:AAAA	GX:Z:G2	xf:i:1	fx:Z:G2	li:i:0
r2:Lc	131	1	1	255	1M	*	0	0	T	F	CB:Z:AAAA-1	UR:Z:AAAA	UB:Z:AAAA	GX:Z:G2	xf:i:1	fx:Z:G2	li:i:0
# UMI that should be corrected
r3:udc	67	1	1	255	1M	*	0	0	A	F	CB:Z:AAAA-1	UR:Z:AATA	UB:Z:AATA	GX:Z:G1	xf:i:1	fx:Z:G1	li:i:0
r3:udc	131	1	1	255	1M	*	0	0	T	F	CB:Z:AAAA-1	UR:Z:AATA	UB:Z:AATA	GX:Z:G1	xf:i:1	fx:Z:G1	li:i:0
# Pair maps to no genes
r4:	67	1	1	255	1M	*	0	0	A	F	CB:Z:AAAA-1	UR:Z:AATA	UB:Z:AATA	li:i:0
r4:	131	1	1	255	1M	*	0	0	T	F	CB:Z:AAAA-1	UR:Z:AATA	UB:Z:AATA	li:i:0
# Pair maps to multiple genes
r5:	67	1	1	255	1M	*	0	0	A	F	CB:Z:AAAA-1	UR:Z:AATA	UB:Z:AATA	GX:Z:G1;G2	fx:Z:G1;G2	li:i:0
r5:	131	1	1	255	1M	*	0	0	T	F	CB:Z:AAAA-1	UR:Z:AATA	UB:Z:AATA	GX:Z:G1;G2	fx:Z:G1;G2	li:i:0
# Pair maps to discordant genes
r15:	67	1	1	255	1M	*	0	0	A	F	CB:Z:AAAA-1	UR:Z:AATA	UB:Z:AATA	gX:Z:G1	li:i:0
r15:	131	1	1	255	1M	*	0	0	T	F	CB:Z:AAAA-1	UR:Z:AATA	UB:Z:AATA	gX:Z:G2	li:i:0
# Low MAPQ
r6:	67	1	1	0	1M	*	0	0	A	F	CB:Z:AAAA-1	UR:Z:AATA	UB:Z:AATA	GX:Z:G1	fx:Z:G1	li:i:0
r6:	131	1	1	0	1M	*	0	0	T	F	CB:Z:AAAA-1	UR:Z:AATA	UB:Z:AATA	GX:Z:G1	fx:Z:G1	li:i:0
# No UB
r12:	67	1	1	255	1M	*	0	0	A	F	CB:Z:AAAA-1	UR:Z:AATA	GX:Z:G1	fx:Z:G1	li:i:0
r12:	131	1	1	255	1M	*	0	0	T	F	CB:Z:AAAA-1	UR:Z:AATA	GX:Z:G1	fx:Z:G1	li:i:0
# Repeated UB
r7:cC	67	1	1	255	1M	*	0	0	A	F	CB:Z:AAAA-1	UR:Z:CCCC	UB:Z:CCCC	GX:Z:G1	xf:i:1	fx:Z:G1	li:i:0
r7:c	131	1	1	255	1M	*	0	0	T	F	CB:Z:AAAA-1	UR:Z:CCCC	UB:Z:CCCC	GX:Z:G1	xf:i:1	fx:Z:G1	li:i:0
r8:dc	67	1	1	255	1M	*	0	0	A	F	CB:Z:AAAA-1	UR:Z:CCCC	UB:Z:CCCC	GX:Z:G1	xf:i:1	fx:Z:G1	li:i:0
r8:dc	131	1	1	255	1M	*	0	0	T	F	CB:Z:AAAA-1	UR:Z:CCCC	UB:Z:CCCC	GX:Z:G1	xf:i:1	fx:Z:G1	li:i:0
# Repeated CB
r9:cC	67	1	1	255	1M	*	0	0	A	F	CB:Z:AATT-1	UR:Z:AAAA	UB:Z:AAAA	GX:Z:G1	xf:i:1	fx:Z:G1	li:i:0
r9:c	131	1	1	255	1M	*	0	0	T	F	CB:Z:AATT-1	UR:Z:AAAA	UB:Z:AAAA	GX:Z:G1	xf:i:1	fx:Z:G1	li:i:0
r10:dc	67	1	1	255	1M	*	0	0	A	F	CB:Z:AATT-1	UR:Z:AAAA	UB:Z:AAAA	GX:Z:G1	xf:i:1	fx:Z:G1	li:i:0
r10:dc	131	1	1	255	1M	*	0	0	T	F	CB:Z:AATT-1	UR:Z:AAAA	UB:Z:AAAA	GX:Z:G1	xf:i:1	fx:Z:G1	li:i:0
# Pair half-maps to gene
r13:cC	67	1	1	255	1M	*	0	0	A	F	CB:Z:AAGG-1	UR:Z:AAAA	UB:Z:AAAA	GX:Z:G1	xf:i:1	fx:Z:G1	li:i:0
r13:c	131	1	1	255	1M	*	0	0	T	F	CB:Z:AAGG-1	UR:Z:AAAA	UB:Z:AAAA	GX:Z:G1	xf:i:1	fx:Z:G1	li:i:0
# No CB
r11:	67	1	1	255	1M	*	0	0	A	F	UR:Z:AAAA	UB:Z:AAAA	GX:Z:G3	xf:i:1	fx:Z:G3	li:i:0
r11:	131	1	1	255	1M	*	0	0	T	F	UR:Z:AAAA	UB:Z:AAAA	GX:Z:G3	xf:i:1	fx:Z:G3	li:i:0
# Secondary alignments, primary comes after
r14:	323	1	1	0	1M	*	0	0	A	F	CB:Z:AATT-1	UR:Z:AAAA	UB:Z:AAAA	xf:i:1	li:i:0
r14:	387	1	1	0	1M	*	0	0	T	F	CB:Z:AATT-1	UR:Z:AAAA	UB:Z:AAAA	xf:i:1	li:i:0
r14:cC	67	1	1	255	1M	*	0	0	A	F	CB:Z:AATT-1	UR:Z:AAAA	UB:Z:AAAA	GX:Z:G3	xf:i:1	fx:Z:G3	li:i:0
r14:c	131	1	1	255	1M	*	0	0	T	F	CB:Z:AATT-1	UR:Z:AAAA	UB:Z:AAAA	GX:Z:G3	xf:i:1	fx:Z:G3	li:i:0
";
        test_mark_dups(sam_text, &"pe");
    }


    extern crate rand;
    use self::rand::Rng;

    /// Generate a random nucleotide sequence of length len
    fn random_sequence<R: Rng>(len: usize, rng: &mut R) -> Vec<u8> {
        let nucs = b"ACGT";
        (0..len).map(|_| nucs[rng.gen_range(0, 4)]).collect()
    }

    fn create_sorted_bam_with_bc<P: AsRef<Path>, R: Rng>(filename: P,
                                                         barcodes: &Vec<Vec<u8>>,
                                                         num_tids: u64,
                                                         reads_per_tid: u64,
                                                         rng: &mut R) {
        use rust_htslib::bam::header::{Header, HeaderRecord};

        let read_len = 100;
        let qual = vec![b'I'; read_len];
        let cigar = bam::record::CigarString(vec![bam::record::Cigar::Match(read_len as u32)]);

        let mut header = Header::new();
        for i in 0..num_tids {
            header.push_record(
                HeaderRecord::new(b"SQ")
                    .push_tag(b"SN", &format!("chr{}", 1 + i))
                    .push_tag(b"LN", &10000000),
            );
      	}
        let mut writer =
            bam::Writer::from_path(filename, &header).expect("Failed to create bam writer");

        for tid in 0..num_tids {
            for pos in 0..reads_per_tid {
                let mut rec = bam::Record::new();
                let seq = random_sequence(read_len, rng);
                rec.set(b"1", &cigar, &seq, &qual);
                rec.set_tid(tid as i32);
                rec.set_pos(pos as i32);
                let bc = rng.choose(&barcodes).unwrap();
                rec.push_aux(b"CB", &bam::record::Aux::String(&bc[..]));
                writer.write(&rec).expect("Failed to write BAM record");
     		}
       	}
        drop(writer);
    }

    #[test]
    fn test_bam_comp_ratio() {
		let test_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("test_output").join("comp_ratio");
		fs::create_dir_all(&test_dir).expect("Failed to create test output dir");
		let test_bam = test_dir.join("test.bam");
        let mut rng = rand::thread_rng();

        let barcodes : Vec<Vec<u8>> = (0..10).map(|_| random_sequence(12, &mut rng)).collect();
        create_sorted_bam_with_bc(&test_bam, &barcodes, 6, 5000, &mut rng);

        let comp_ratio = cmd_mark_dups::bam_compression_ratio(&test_bam);
        println!("{:?}", comp_ratio);
        assert!(comp_ratio >= 0.1 && comp_ratio <= 0.5);
    }

    #[test]
    fn test_chunk_sample_bc() {
		let test_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("test_output").join("sample_bc");
        fs::create_dir_all(&test_dir).expect("Failed to create test output dir");
        let mut rng = rand::thread_rng();

        // completely random barcode
        let barcodes : Vec<Vec<u8>> = (0..1000).map(|_| random_sequence(12, &mut rng)).collect();
        let test_bam = test_dir.join("test_random.bam");
        create_sorted_bam_with_bc(&test_bam, &barcodes, 6, 5000, &mut rng);

        let mut bam_reader = bam::Reader::from_path(&test_bam).expect("Failed to open input BAM");
        let block_offsets = bam_block_offsets(&test_bam);
        let major_bc_prop = cmd_mark_dups::major_barcode_proportion_within_chunk(&mut bam_reader, &block_offsets, 0, block_offsets.len(), 1000);
        println!("{:?}", major_bc_prop);
        assert!(major_bc_prop >= 0.0 && major_bc_prop <= 0.1);

        // same random barcode
        let barcodes : Vec<Vec<u8>> = (0..1).map(|_| random_sequence(12, &mut rng)).collect();
        let test_bam = test_dir.join("test_same_bc.bam");
        create_sorted_bam_with_bc(&test_bam, &barcodes, 6, 5000, &mut rng);

        let mut bam_reader = bam::Reader::from_path(&test_bam).expect("Failed to open input BAM");
        let block_offsets = bam_block_offsets(&test_bam);
        let major_bc_prop = cmd_mark_dups::major_barcode_proportion_within_chunk(&mut bam_reader, &block_offsets, 0, block_offsets.len(), 1000);
        println!("{:?}", major_bc_prop);
        assert!(major_bc_prop == 1.0);
    }
}
