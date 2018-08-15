//
// Copyright (c) 2017 10x Genomics, Inc. All rights reserved.
//

use debruijn::bits_to_ascii;
use sw;
use utils;
use graph_read;
use fastq;

use std;
use std::process::Command;
use constants::{DUMMY_CONTIG_NAME};
use std::path::{Path};
use std::collections::{HashMap};
use failure::Error;

use rust_htslib::bam::{self, Read};
use rust_htslib::bam::record::CigarString; 

pub fn concatenate_bams(bam_filenames: &Vec<String>, merged_bam: &String) -> Result<(), Error> {
    let mut concat_header = bam::header::Header::new();
    let mut header_rec = bam::header::HeaderRecord::new(b"PG");
    header_rec.push_tag(b"ID", &"vdj_asm concatenate_bams");
    concat_header.push_record(&header_rec);

    // First collect all the headers
    let mut all_contigs = HashMap::new(); // contig -> length
    let mut bam_tid_maps = Vec::new(); // (tid -> contig) for each input bam

    for bam in bam_filenames.iter() {
        let in_bam = bam::Reader::from_path(&Path::new(bam))?;
        let header = in_bam.header();

        let mut bam_tids = HashMap::new();

        for (contig_idx, name_bytes) in header.target_names().iter().enumerate() {
            let name = String::from_utf8_lossy(name_bytes).into_owned();
            bam_tids.insert(contig_idx as i32, name.clone());
            let contig_len = header.target_len(contig_idx as u32).unwrap() as usize;
            if all_contigs.contains_key(&name) {
                assert_eq!(*all_contigs.get(&name).unwrap(), contig_len);
            }
            all_contigs.insert(name, contig_len);
        }
        bam_tid_maps.push(bam_tids);
    }

    let mut sorted_contigs : Vec<String> = all_contigs.keys().map(|x| (*x).clone()).collect();
    sorted_contigs.sort();

    if sorted_contigs.len() > 1 {
        // If the input bam only contains unmapped reads, it might have a dummy contig.
        // This is because we can't have a bam without any reference sequences.
        // Remove the dummy contig if there is at least one non-dummy contig.
        sorted_contigs = sorted_contigs.iter().filter(|&x| *x != DUMMY_CONTIG_NAME.to_string()).map(|x| (*x).clone()).collect();
    }
    let mut header_to_tid = HashMap::new(); // contig -> final_tid


    for (idx, name) in sorted_contigs.iter().enumerate() {
        header_to_tid.insert(name, idx as i32);
        add_ref_to_bam_header(&mut concat_header, name, *all_contigs.get(name).unwrap());
    }

    let mut out_bam = bam::Writer::from_path(Path::new(&merged_bam), &concat_header)?;

    for (bam, ref tid_map) in bam_filenames.iter().zip(bam_tid_maps.iter()) {
        println!("Merging BAM: {}", bam);
        let mut in_bam = bam::Reader::from_path(&Path::new(bam))?;

        for record in in_bam.records() {
            let mut rec = record?;  
            let curr_tid = rec.tid();
            if curr_tid >= 0 {
                rec.set_tid(*header_to_tid.get(tid_map.get(&curr_tid).unwrap()).unwrap());
            }
            let curr_mtid = rec.mtid();
            if curr_mtid >= 0 {
                rec.set_mtid(*header_to_tid.get(tid_map.get(&curr_mtid).unwrap()).unwrap());
            }

            let _ = out_bam.write(&mut rec);
        }
    }

    Ok(())
}

pub fn sort_and_index(unsorted_bam: &String, output_bam: &String) {

    // In older samtools version, the "--version" doesn't work.
    // However, the exit status is still 0.
    let output = Command::new("samtools").arg("--version").output().expect("samtools not installed");

    let output_str = String::from_utf8(output.stdout).unwrap();
    let new_format = output_str.len() > 0;

    if !new_format {
        println!("sort {} {}", unsorted_bam, utils::rm_extension(&output_bam.clone()));
        Command::new("samtools").arg("sort")
                                .arg(unsorted_bam)
                                .arg(utils::rm_extension(&output_bam.clone()))
                                .status().expect("samtools sort failed");
    } else {
        Command::new("samtools").arg("sort")
                                .arg("-o")
                                .arg(output_bam.clone())
                                .arg(unsorted_bam)
                                .status().expect("samtools sort failed");
    }

    Command::new("samtools").arg("index").arg(output_bam.clone()).status().expect("samtools index failed");
}

/// Populate the input Header object with reference sequence info.
pub fn add_ref_to_bam_header(header: &mut bam::header::Header,
                             seq_name: &str, seq_len: usize) {
    let mut header_rec = bam::header::HeaderRecord::new(b"SQ");
    header_rec.push_tag(b"SN", &seq_name);
    header_rec.push_tag(b"LN", &seq_len);
    header.push_record(&header_rec);
}

/// Convert an internal `Read` object into a BAM record.
/// From the SAM spec:
/// 1. For a unmapped paired-end or mate-pair read whose mate is mapped, the unmapped read should
/// have RNAME and POS identical to its mate.
/// 2. If all segments in a template are unmapped, their RNAME should be set as ‘*’ and POS as 0.
/// ...
/// 4. Unmapped reads should be stored in the orientation in which they came off the sequencing machine
/// and have their reverse flag bit (0x10) correspondingly unset.
pub fn read_to_bam_record_opts(read: &graph_read::Read,
                          alignment: &Option<sw::AlignmentPacket>,
                          mate_alignment: &Option<sw::AlignmentPacket>,
                          strip_qname: bool,
                          set_augmented_tags: bool) -> bam::record::Record {


    let mut rec = bam::record::Record::new();

    // Keep these from the input. Alignments never reverse complement, so if a read
    // started as reverse complemented then it remained so after alignment.
    let is_rc = read.is_reverse();
    let is_mate_rc = read.is_mate_reverse();
    let is_paired = read.is_paired();
    assert!(is_paired || mate_alignment.is_none());

    match alignment.as_ref() {
        Some(al) => read.validate_alignment(al),
        None => true,
    };

    // We do NOT ensure that rule 4 is satisfied.
    let adj_qual = read.quals.clone();
    let adj_seq: Vec<u8> = read.seq.iter().map(|x| bits_to_ascii(x)).collect();

    let (unmapped_bit, cigar) = match alignment.as_ref() {
        Some(al) => (0, CigarString::from_alignment(&al.alignment, false)), // bool to specify soft clipping
        None => (4, CigarString(vec![]))
    };

    let new_header = fastq::CellrangerFastqHeader::new(read.name.clone());
    let qname =
        if strip_qname {
            new_header.qname.as_bytes()
        } else {
            new_header.header.as_bytes()
        };

    rec.set(qname, &cigar, &adj_seq, &adj_qual);

    if set_augmented_tags {
        for (tag, value) in new_header.tags.iter() {
            if value.len() > 0 {
                rec.push_aux(tag.as_bytes(), &bam::record::Aux::String(value.as_bytes())).unwrap();
            }
        }
    }

    let (tid, pos, mapq) = match (alignment.as_ref(), mate_alignment.as_ref()) {
        (Some(al), Some(mate_al)) => {
            if is_paired && mate_al.ref_idx == al.ref_idx {
                let max_end = std::cmp::max(al.alignment.yend, mate_al.alignment.yend) as i32;
                let min_start = std::cmp::min(al.alignment.ystart, mate_al.alignment.ystart) as i32;
                if al.alignment.ystart <= mate_al.alignment.ystart {
                    rec.set_insert_size(max_end - min_start);
                } else {
                    rec.set_insert_size(mate_al.alignment.ystart as i32 - al.alignment.yend as i32);
                }
            }
            (al.ref_idx as i32, al.alignment.ystart as i32, 255)
        },
        (Some(al), None) => (al.ref_idx as i32, al.alignment.ystart as i32, 255),
        (None, Some(mate_al)) => (mate_al.ref_idx as i32, mate_al.alignment.ystart as i32, 0),
        (None, None) => (-1, 0, 0),
    };
    rec.set_tid(tid);
    rec.set_pos(pos);

    rec.set_mapq(mapq); // this is not really a MAPQ...

    let mate_unmapped_bit;
    match mate_alignment.as_ref() {
        Some(al) => {
            rec.set_mtid(al.ref_idx as i32);
            rec.set_mpos(al.alignment.ystart as i32);
            mate_unmapped_bit = 0;
        },
        None => {
            rec.set_mtid(tid);
            rec.set_mpos(pos);
            mate_unmapped_bit = 8;
        }
    }

    rec.set_flags(
        is_paired as u16 + // assume always paired
        unmapped_bit +
        mate_unmapped_bit +
        (is_rc as u16) * 16 + // is reverse complemented
        (is_mate_rc as u16) * 32 + // assume if one is RC then the other is not
        {if read.is_first_in_template() {64} else {128}});

    match alignment.as_ref() {
        Some(al) => {
            rec.push_aux(b"AS", &bam::record::Aux::Integer(al.alignment.score as i64)).unwrap();
            rec.push_aux(b"NM", &bam::record::Aux::Integer(al.edit_distance() as i64)).unwrap();
        },
        None => {},
    };

    rec
}

pub fn read_to_bam_record(read: &graph_read::Read,
                          alignment: &Option<sw::AlignmentPacket>,
                          mate_alignment: &Option<sw::AlignmentPacket>) -> bam::record::Record {

    read_to_bam_record_opts(read, alignment, mate_alignment, false, false)
                          }


#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use rust_htslib::bam;
    use rust_htslib::bam::Read;
    use std::path::Path;
    
    const OUTDIR : &'static str = "test/outputs/bam";

    fn init_test() {
        if Path::new(OUTDIR).is_dir() {
            return;
        }
        let _ = fs::create_dir_all(OUTDIR);
    }

    #[test]
    fn test_sort_and_index() {
        init_test();

        let test_bam_name = "test/inputs/test_index.bam";
        let true_sorted_bam_name = "test/inputs/test_index_sorted.bam";
        let out_bam_name = "test/outputs/bam/test_index_sorted.bam";

        sort_and_index(&test_bam_name.to_string(), &out_bam_name.to_string());

        {
            let mut true_sorted_bam = bam::Reader::from_path(&Path::new(&(true_sorted_bam_name.to_string()))).unwrap();
            let mut sorted_bam = bam::Reader::from_path(&out_bam_name).unwrap();

            // Don't test the byte arrays. Some versions of samtools seem to omit the SO line.
            //assert_eq!(true_sorted_bam.header().as_bytes(), sorted_bam.header().as_bytes());
            {
                let true_names = true_sorted_bam.header().target_names();
                let sorted_names = sorted_bam.header().target_names();
                assert_eq!(true_names, sorted_names);
            }

            for (read1, read2) in true_sorted_bam.records().zip(sorted_bam.records()) {
                let r1 = read1.unwrap();
                let r2 = read2.unwrap();
                assert_eq!(r1.qname(), r2.qname());
                assert_eq!(r1.seq().as_bytes(), r2.seq().as_bytes());
            }
        }

        {
            let mut true_sorted_bam = bam::IndexedReader::from_path(&Path::new(&(true_sorted_bam_name.to_string()))).unwrap();
            let mut sorted_bam = bam::IndexedReader::from_path(&out_bam_name).unwrap();

            let _ = true_sorted_bam.fetch(1, 0, 1000);
            let _ = sorted_bam.fetch(1, 0, 1000);
            assert_eq!(sorted_bam.records().count(), 4);

            let _ = true_sorted_bam.fetch(1, 0, 1000);
            let _ = sorted_bam.fetch(1, 0, 1000);
            assert_eq!(true_sorted_bam.records().count(), sorted_bam.records().count());
        }
    }

    // FIXME -- the test_split1.bam has invalid CIGAR / seq length
    //#[test]
    fn test_concatenate_bams() {
        init_test();
        let merged_bam_name = "test/inputs/test_index.bam";
        let out_bam_name = "test/outputs/bam/test_merged.bam";

        let split_bam_names = vec!["test/inputs/test_split1.bam".to_string(),
                                   "test/inputs/test_split2.bam".to_string()];

        if let Err(v) = concatenate_bams(&split_bam_names, &out_bam_name.to_string()) {
            println!("{}, {}", v.cause(), v.backtrace());
            assert!(false);
        };

        let mut true_bam = bam::Reader::from_path(&Path::new(&(merged_bam_name.to_string()))).unwrap();
        let mut merged_bam = bam::Reader::from_path(&out_bam_name).unwrap();

        {
            let true_header = true_bam.header();
            let merged_header = merged_bam.header();

            for (true_name, other_name) in true_header.target_names().iter().zip(merged_header.target_names().iter()) {
                assert_eq!(true_name, other_name);
            }
        }

        for (read1, read2) in true_bam.records().zip(merged_bam.records()) {
            let r1 = read1.unwrap();
            let r2 = read2.unwrap();
            assert_eq!(r1.tid(), r2.tid());
            assert_eq!(r1.qname(), r2.qname());
            assert_eq!(r1.seq().as_bytes(), r2.seq().as_bytes());
        }
    }
}
