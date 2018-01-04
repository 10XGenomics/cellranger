//
// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
//

use std::collections::{HashMap, HashSet};
use std::str;
use std::default::Default;
use itertools::Itertools;
use rust_htslib::bam;
use tenkit::chunk_bam::{chunk_bam_records, BamChunkIter};
use martian::*;
use serde_json;

use utils;

pub struct MarkDuplicatesStage;

#[derive(Serialize, Deserialize, Default)]
struct Metrics {
    total_reads: u64,
    low_support_umi_reads: u64,
    umi_corrected_reads: u64,
    candidate_dup_reads: u64,
    dup_reads: u64,
    umis: u64,
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
/// find the single gene confidently mapped to if any.
fn get_qname_conf_mapped_gene<'a, I>(records: I) -> Option<String>
    where I: Iterator<Item = &'a bam::record::Record> {
    let mut genes = Vec::new();
    for record in records.filter(|&r| !r.is_secondary() &&
                                 utils::get_read_extra_flags(r).intersects(utils::ExtraFlags::CONF_MAPPED)) {
        genes.extend(utils::get_read_gene_ids(record).into_iter());
    }
    genes.sort_unstable();
    genes.dedup();
    match genes.len() {
        1 => genes.into_iter().nth(0),
        _ => None,
    }
}

fn cmd_mark_dups(args: &ChunkArgs, outs: JsonDict) -> JsonDict {
        use bam::Read;

        let mut metrics: Metrics = Default::default();

        let mut bam = bam::Reader::from_path(&args.input)
            .expect("Failed to open input BAM file");

        let mut out_bam = bam::Writer::from_path(outs["alignments"].as_str().unwrap(),
                                                 &bam::Header::from_template(bam.header()))
            .expect("Failed to open output BAM file");

        let range = (args.chunk_start, args.chunk_end);
        let chunk_iter = BamChunkIter::new(&mut bam, range);

        for (_bc, bc_group) in chunk_iter.map(|x| x.unwrap()).group_by(utils::get_read_barcode) {
            if _bc.is_none() {
                // Invalid barcode. Write all records as-is.
                for record in bc_group {
                    metrics.total_reads += 1;
                    out_bam.write(&record).expect("Failed to write BAM record");
                }
                // Skip this barcode.
                continue;
            }

            // Get raw UMI frequencies
            let mut raw_umigene_counts: HashMap<(Vec<u8>, String), u64> = HashMap::new();
            let mut low_support_umigenes: HashSet<(Vec<u8>, String)> = HashSet::new();

            for (maybe_umi, umi_group) in bc_group.iter().group_by(|x| utils::get_read_umi(&x)) {
                let mut gene_counts: HashMap<String, u64> = HashMap::new();

                if let Some(umi) = maybe_umi {
                    // Count (UMI, gene) pairs.

                    // Assumes records are qname-contiguous in the input.
                    for (qname, qname_records) in umi_group.into_iter()
                        .filter(|&x| utils::is_read_dup_candidate(x))
                        .group_by(|x| x.qname()) {
                            match get_qname_conf_mapped_gene(qname_records.into_iter()) {
                                None => { panic!(format!("Found 0 or >1 genes for confidently mapped read/pair {}", str::from_utf8(qname).unwrap())) },
                                Some(gene) => {
                                    *raw_umigene_counts.entry((umi.clone(), gene.clone())).or_insert(0) += 1;
                                    *gene_counts.entry(gene.clone()).or_insert(0) += 1;
                                },
                            }
                    } // for each qname

                    // Mark (UMI, gene) pairs w/ frequency below the max for the UMI as low support.
                    if let Some((_max_gene, max_count)) = gene_counts.iter().max_by_key(|x| x.1) {
                        for (gene, count) in gene_counts.iter() {
                            if args.filter_umis && count < max_count {
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

                let maybe_gene = get_qname_conf_mapped_gene(qname_records.iter()
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
                        metrics.total_reads += 1;
                    }
                    metrics.candidate_dup_reads += is_dup_candidate as u64;

                    if let (Some(umi), Some(gene)) = (maybe_umi.as_ref(), maybe_gene.as_ref()) {
                        let key = (umi.clone(), gene.clone());

                        if !record.is_secondary() && low_support_umigenes.contains(&key) {
                            // Low support (UMI, gene). Mark as low support.
                            // - Only consider primary alignments for this flag.
                            // - Do not correct the UMI.
                            // - Do not mark duplicates w/ for this (UMI, gene).
                            metrics.low_support_umi_reads += 1;
                            new_flags |= utils::ExtraFlags::LOW_SUPPORT_UMI;

                        } else {
                            // Correct UMIs in all records
                            let (corrected_umi, is_corrected) = match umi_corrections.get(&key) {
                                Some(new_umi) => (new_umi.clone(), true),
                                None => (umi.clone(), false),
                            };

                            // Correct the UMI tag
                            if is_corrected {
                                metrics.umi_corrected_reads += 1;
                                new_record.remove_aux(utils::PROC_UMI_SEQ_TAG);
                                new_record.push_aux(utils::PROC_UMI_SEQ_TAG,
                                                    &bam::record::Aux::String(&corrected_umi));
                            }

                            // Don't try to dup mark secondary alignments.
                            if is_dup_candidate {
                                let dup_key = (corrected_umi, key.1, get_mate_type(&record));

                                if wrote_umigenes.contains(&dup_key) {
                                    // Duplicate
                                    metrics.dup_reads += 1;
                                    let flags = record.flags();
                                    new_record.set_flags(flags | 1024u16);
                                } else {
                                    // Non-duplicate
                                    wrote_umigenes.insert(dup_key);

                                    // Flag read1 as countable
                                    if !record.is_last_in_template() {
                                        metrics.umis += 1;
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

            } // for each record
        } // for each barcode

        // Write metrics
        utils::write_json_file(outs["metrics"].as_str().unwrap(), &metrics)
            .expect("Failed to write metrics to JSON file");

        outs
}


impl MartianStage for MarkDuplicatesStage {
    fn split(&self, args: JsonDict) -> JsonDict {
        let bam = args["input"].as_str().unwrap();

        let chunk_intervals =
            chunk_bam_records(&bam, &|x| utils::get_read_barcode(&x), 0.5, 256);

        let chunks: Vec<serde_json::Value> = chunk_intervals.into_iter().map(|ci|
            json!({
                "chunk_start": ci.0,
                "chunk_end": ci.1,
            })
        ).collect();

        serde_json::from_value(json!({
            "chunks": json!(chunks),
        })).expect("Failed to serialize split outs")
    }

    fn main(&self, json_args: JsonDict, outs: JsonDict) -> JsonDict {
        let args: ChunkArgs = obj_decode(json_args);
        cmd_mark_dups(&args, outs)
    }

    fn join(&self, _: JsonDict, outs: JsonDict, _chunk_defs: Vec<JsonDict>, chunk_outs: Vec<JsonDict>) -> JsonDict {
        let mut final_outs = outs.clone();

        // Collect output BAMs
        final_outs["output"] = json!(
            chunk_outs.iter()
            .map(|x| x["alignments"].clone())
                .collect::<Vec<serde_json::Value>>());

        // Merge summary metrics
        let mut metrics: Metrics = Default::default();
        for chunk_out in chunk_outs.iter() {
            let chunk: Metrics = utils::read_json_file(chunk_out["metrics"].as_str().unwrap())
                .expect("Failed to read from metrics JSON file");
            metrics.total_reads += chunk.total_reads;
            metrics.umi_corrected_reads += chunk.umi_corrected_reads;
            metrics.dup_reads += chunk.dup_reads;
            metrics.umis += chunk.umis;
            metrics.candidate_dup_reads += chunk.candidate_dup_reads;
            metrics.low_support_umi_reads += chunk.low_support_umi_reads;
        }

        // Write summary metrics
        utils::write_json_file(final_outs["summary"].as_str().unwrap(), &json!({
            "corrected_umi_frac": metrics.umi_corrected_reads as f64 / metrics.total_reads as f64,
            "low_support_umi_reads_frac": metrics.low_support_umi_reads as f64 / metrics.candidate_dup_reads as f64,
            "multi_cdna_pcr_dupe_reads_frac": metrics.dup_reads as f64 / metrics.candidate_dup_reads as f64,
        }))
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
    use tenkit::chunk_bam::chunk_bam_records;
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

        let metrics: cmd_mark_dups::Metrics = utils::read_json_file(&metrics_filename).expect("Failed to read metrics file");
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
        let chunk_intervals =
            chunk_bam_records(&in_bam_filename, &|x| utils::get_read_barcode(&x), 0.5, 256);

        let chunk_args = cmd_mark_dups::ChunkArgs {
            input: in_bam_filename.to_str().unwrap().to_owned(),
            chunk_start: chunk_intervals[0].0,
            chunk_end: chunk_intervals[0].1,
            filter_umis: true,
        };

        // Setup output
        let metrics_filename = test_dir.join("metrics.json");
        let out_bam_filename = test_dir.join("output.bam");
        let chunk_outs = match json!({
            "metrics": metrics_filename.clone(),
            "alignments": out_bam_filename.clone(),
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
r0:cC	0	1	1	255	1M	*	0	0	A	F	CB:Z:AAAA-1	UR:Z:AAAA	UB:Z:AAAA	GX:Z:G1	xf:i:1
r1:dc	0	1	1	255	1M	*	0	0	A	F	CB:Z:AAAA-1	UR:Z:AAAA	UB:Z:AAAA	GX:Z:G1	xf:i:1
r1:	256	1	2	0	1M	*	0	0	A	F	CB:Z:AAAA-1	UR:Z:AAAA	UB:Z:AAAA	GX:Z:G1
r2:Lc	0	1	1	255	1M	*	0	0	A	F	CB:Z:AAAA-1	UR:Z:AAAA	UB:Z:AAAA	GX:Z:G2	xf:i:1
r3:udc	0	1	1	255	1M	*	0	0	A	F	CB:Z:AAAA-1	UR:Z:AATA	UB:Z:AATA	GX:Z:G1	xf:i:1
r4:	0	1	1	255	1M	*	0	0	A	F	CB:Z:AAAA-1	UR:Z:AATA	UB:Z:AATA
r5:	0	1	1	255	1M	*	0	0	A	F	CB:Z:AAAA-1	UR:Z:AATA	UB:Z:AATA	GX:Z:G1;G2
r6:	0	1	1	0	1M	*	0	0	A	F	CB:Z:AAAA-1	UR:Z:AATA	UB:Z:AATA	GX:Z:G1
r12:	0	1	1	255	1M	*	0	0	A	F	CB:Z:AAAA-1	UR:Z:AATA	GX:Z:G1
r7:cC	0	1	1	255	1M	*	0	0	A	F	CB:Z:AAAA-1	UR:Z:CCCC	UB:Z:CCCC	GX:Z:G1	xf:i:1
r8:dc	0	1	1	255	1M	*	0	0	A	F	CB:Z:AAAA-1	UR:Z:CCCC	UB:Z:CCCC	GX:Z:G1	xf:i:1
r9:cC	0	1	1	255	1M	*	0	0	A	F	CB:Z:AATT-1	UR:Z:AAAA	UB:Z:AAAA	GX:Z:G1	xf:i:1
r10:dc	0	1	1	255	1M	*	0	0	A	F	CB:Z:AATT-1	UR:Z:AAAA	UB:Z:AAAA	GX:Z:G1	xf:i:1
r11:	0	1	1	255	1M	*	0	0	A	F	UR:Z:AAAA	UB:Z:AAAA	GX:Z:G3	xf:i:1
";
        test_mark_dups(sam_text, &"se");
    }

    #[test]
    fn test_mark_dups_pe() {
        let sam_text = b"
# Non-duplicate
r0:cC	67	1	1	255	1M	*	0	0	A	F	CB:Z:AAAA-1	UR:Z:AAAA	UB:Z:AAAA	GX:Z:G1	xf:i:1
r0:c	131	1	1	255	1M	*	0	0	T	F	CB:Z:AAAA-1	UR:Z:AAAA	UB:Z:AAAA	GX:Z:G1	xf:i:1
# Duplicate
r1:dc	67	1	1	255	1M	*	0	0	A	F	CB:Z:AAAA-1	UR:Z:AAAA	UB:Z:AAAA	GX:Z:G1	xf:i:1
r1:dc	131	1	1	255	1M	*	0	0	T	F	CB:Z:AAAA-1	UR:Z:AAAA	UB:Z:AAAA	GX:Z:G1	xf:i:1
# Secondary alignment for duplicate
r1:	323	1	2	0	1M	*	0	0	A	F	CB:Z:AAAA-1	UR:Z:AAAA	UB:Z:AAAA	GX:Z:G1
r1:	387	1	2	0	1M	*	0	0	T	F	CB:Z:AAAA-1	UR:Z:AAAA	UB:Z:AAAA	GX:Z:G1
# Low-support gene for this UMI
r2:Lc	67	1	1	255	1M	*	0	0	A	F	CB:Z:AAAA-1	UR:Z:AAAA	UB:Z:AAAA	GX:Z:G2	xf:i:1
r2:Lc	131	1	1	255	1M	*	0	0	T	F	CB:Z:AAAA-1	UR:Z:AAAA	UB:Z:AAAA	GX:Z:G2	xf:i:1
# UMI that should be corrected
r3:udc	67	1	1	255	1M	*	0	0	A	F	CB:Z:AAAA-1	UR:Z:AATA	UB:Z:AATA	GX:Z:G1	xf:i:1
r3:udc	131	1	1	255	1M	*	0	0	T	F	CB:Z:AAAA-1	UR:Z:AATA	UB:Z:AATA	GX:Z:G1	xf:i:1
# Pair maps to no genes
r4:	67	1	1	255	1M	*	0	0	A	F	CB:Z:AAAA-1	UR:Z:AATA	UB:Z:AATA
r4:	131	1	1	255	1M	*	0	0	T	F	CB:Z:AAAA-1	UR:Z:AATA	UB:Z:AATA
# Pair maps to multiple genes
r5:	67	1	1	255	1M	*	0	0	A	F	CB:Z:AAAA-1	UR:Z:AATA	UB:Z:AATA	GX:Z:G1;G2
r5:	131	1	1	255	1M	*	0	0	T	F	CB:Z:AAAA-1	UR:Z:AATA	UB:Z:AATA	GX:Z:G1;G2
# Pair maps to discordant genes
r15:	67	1	1	255	1M	*	0	0	A	F	CB:Z:AAAA-1	UR:Z:AATA	UB:Z:AATA	gX:Z:G1
r15:	131	1	1	255	1M	*	0	0	T	F	CB:Z:AAAA-1	UR:Z:AATA	UB:Z:AATA	gX:Z:G2
# Low MAPQ
r6:	67	1	1	0	1M	*	0	0	A	F	CB:Z:AAAA-1	UR:Z:AATA	UB:Z:AATA	GX:Z:G1
r6:	131	1	1	0	1M	*	0	0	T	F	CB:Z:AAAA-1	UR:Z:AATA	UB:Z:AATA	GX:Z:G1
# No UB
r12:	67	1	1	255	1M	*	0	0	A	F	CB:Z:AAAA-1	UR:Z:AATA	GX:Z:G1
r12:	131	1	1	255	1M	*	0	0	T	F	CB:Z:AAAA-1	UR:Z:AATA	GX:Z:G1
# Repeated UB
r7:cC	67	1	1	255	1M	*	0	0	A	F	CB:Z:AAAA-1	UR:Z:CCCC	UB:Z:CCCC	GX:Z:G1	xf:i:1
r7:c	131	1	1	255	1M	*	0	0	T	F	CB:Z:AAAA-1	UR:Z:CCCC	UB:Z:CCCC	GX:Z:G1	xf:i:1
r8:dc	67	1	1	255	1M	*	0	0	A	F	CB:Z:AAAA-1	UR:Z:CCCC	UB:Z:CCCC	GX:Z:G1	xf:i:1
r8:dc	131	1	1	255	1M	*	0	0	T	F	CB:Z:AAAA-1	UR:Z:CCCC	UB:Z:CCCC	GX:Z:G1	xf:i:1
# Repeated CB
r9:cC	67	1	1	255	1M	*	0	0	A	F	CB:Z:AATT-1	UR:Z:AAAA	UB:Z:AAAA	GX:Z:G1	xf:i:1
r9:c	131	1	1	255	1M	*	0	0	T	F	CB:Z:AATT-1	UR:Z:AAAA	UB:Z:AAAA	GX:Z:G1	xf:i:1
r10:dc	67	1	1	255	1M	*	0	0	A	F	CB:Z:AATT-1	UR:Z:AAAA	UB:Z:AAAA	GX:Z:G1	xf:i:1
r10:dc	131	1	1	255	1M	*	0	0	T	F	CB:Z:AATT-1	UR:Z:AAAA	UB:Z:AAAA	GX:Z:G1	xf:i:1
# Pair half-maps to gene
r13:cC	67	1	1	255	1M	*	0	0	A	F	CB:Z:AAGG-1	UR:Z:AAAA	UB:Z:AAAA	GX:Z:G1	xf:i:1
r13:c	131	1	1	255	1M	*	0	0	T	F	CB:Z:AAGG-1	UR:Z:AAAA	UB:Z:AAAA	GX:Z:G1	xf:i:1
# No CB
r11:	67	1	1	255	1M	*	0	0	A	F	UR:Z:AAAA	UB:Z:AAAA	GX:Z:G3	xf:i:1
r11:	131	1	1	255	1M	*	0	0	T	F	UR:Z:AAAA	UB:Z:AAAA	GX:Z:G3	xf:i:1
# Secondary alignments, primary comes after
r14:	323	1	1	0	1M	*	0	0	A	F	CB:Z:AATT-1	UR:Z:AAAA	UB:Z:AAAA	xf:i:1
r14:	387	1	1	0	1M	*	0	0	T	F	CB:Z:AATT-1	UR:Z:AAAA	UB:Z:AAAA	xf:i:1
r14:cC	67	1	1	255	1M	*	0	0	A	F	CB:Z:AATT-1	UR:Z:AAAA	UB:Z:AAAA	GX:Z:G3	xf:i:1
r14:c	131	1	1	255	1M	*	0	0	T	F	CB:Z:AATT-1	UR:Z:AAAA	UB:Z:AAAA	GX:Z:G3	xf:i:1
";
        test_mark_dups(sam_text, &"pe");
    }
}
