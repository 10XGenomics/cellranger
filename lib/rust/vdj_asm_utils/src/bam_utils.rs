// ◼ There are some fundamental problems with rust-htslib that affect our ability
// ◼ to debug our code:
// ◼ (a) the lack of code to print a bam record
// ◼ (b) lack of code to validate a bam record and say what’s wrong with it
// ◼ (c) bad exit handling by rust-htslib -- it doesn’t panic upon failure so you
// ◼     don’t get a traceback
// ◼ (d) inability to test a change to rust-htslib without doing cargo clean on
// ◼     the project that uses it.

use crate::constants::DUMMY_CONTIG_NAME;
use crate::{fastq, graph_read, sw, utils};
use anyhow::Result;
use debruijn::bits_to_ascii;
use rust_htslib::bam::record::CigarString;
use rust_htslib::bam::{self, Read};
use std;
use std::cmp::{max, min};
use std::collections::HashMap;
use std::iter::zip;
use std::path::Path;
use std::process::Command;

pub fn concatenate_bams(bam_filenames: &[&Path], merged_bam: &Path) -> Result<()> {
    let mut concat_header = bam::header::Header::new();
    let mut header_rec = bam::header::HeaderRecord::new(b"PG");
    header_rec.push_tag(b"ID", "vdj_asm concatenate_bams");
    concat_header.push_record(&header_rec);

    // First collect all the headers
    let mut all_contigs = HashMap::new(); // contig -> length
    let mut bam_tid_maps = Vec::new(); // (tid -> contig) for each input bam

    for bam in bam_filenames {
        let in_bam = bam::Reader::from_path(bam)?;
        let header = in_bam.header();

        let mut bam_tids = HashMap::new();

        for (contig_idx, name_bytes) in header.target_names().iter().enumerate() {
            let name = String::from_utf8_lossy(name_bytes).into_owned();
            bam_tids.insert(contig_idx as i32, name.clone());
            let contig_len = header.target_len(contig_idx as u32).unwrap() as usize;
            if let Some(&len) = all_contigs.get(&name) {
                assert_eq!(len, contig_len);
            }
            all_contigs.insert(name, contig_len);
        }
        bam_tid_maps.push(bam_tids);
    }

    let mut sorted_contigs: Vec<String> = all_contigs.keys().map(|x| (*x).clone()).collect();
    sorted_contigs.sort();

    if sorted_contigs.len() > 1 {
        // If the input bam only contains unmapped reads, it might have a dummy contig.
        // This is because we can't have a bam without any reference sequences.
        // Remove the dummy contig if there is at least one non-dummy contig.
        sorted_contigs = sorted_contigs
            .iter()
            .filter(|&x| x != DUMMY_CONTIG_NAME)
            .map(|x| (*x).clone())
            .collect();
    }
    let mut header_to_tid = HashMap::new(); // contig -> final_tid

    for (idx, name) in sorted_contigs.iter().enumerate() {
        header_to_tid.insert(name, idx as i32);
        add_ref_to_bam_header(&mut concat_header, name, all_contigs[name]);
    }

    let mut out_bam = bam::Writer::from_path(merged_bam, &concat_header, bam::Format::Bam)?;

    for (bam, tid_map) in zip(bam_filenames, bam_tid_maps) {
        println!("Merging BAM: {}", bam.display());
        let mut in_bam = bam::Reader::from_path(bam)?;

        for record in in_bam.records() {
            let mut rec = record?;
            let curr_tid = rec.tid();
            if curr_tid >= 0 {
                rec.set_tid(header_to_tid[&tid_map[&curr_tid]]);
            }
            let curr_mtid = rec.mtid();
            if curr_mtid >= 0 {
                rec.set_mtid(header_to_tid[&tid_map[&curr_mtid]]);
            }

            let _ = out_bam.write(&rec);
        }
    }

    Ok(())
}

pub fn sort_and_index(unsorted_bam: &Path, output_bam: &Path) {
    // In older samtools version, the "--version" doesn't work.
    // However, the exit status is still 0.
    let output = Command::new("samtools")
        .arg("--version")
        .output()
        .expect("samtools not installed");

    let output_str = String::from_utf8(output.stdout).unwrap();
    let new_format = !output_str.is_empty();

    if !new_format {
        println!(
            "sort {} {}",
            unsorted_bam.display(),
            utils::rm_extension(output_bam).display()
        );
        Command::new("samtools")
            .arg("sort")
            .arg(unsorted_bam)
            .arg(utils::rm_extension(output_bam))
            .status()
            .expect("samtools sort failed");
    } else {
        Command::new("samtools")
            .arg("sort")
            .arg("-o")
            .arg(output_bam)
            .arg(unsorted_bam)
            .status()
            .expect("samtools sort failed");
    }

    Command::new("samtools")
        .arg("index")
        .arg(output_bam)
        .status()
        .expect("samtools index failed");
}

/// Populate the input Header object with reference sequence info.
pub fn add_ref_to_bam_header(header: &mut bam::header::Header, seq_name: &str, seq_len: usize) {
    let mut header_rec = bam::header::HeaderRecord::new(b"SQ");
    header_rec.push_tag(b"SN", seq_name);
    header_rec.push_tag(b"LN", seq_len);
    header.push_record(&header_rec);
}

/// Convert an internal `Read` object into a BAM record.
/// From the SAM spec:
/// 1. For a unmapped paired-end or mate-pair read whose mate is mapped, the unmapped read should
///    have RNAME and POS identical to its mate.
/// 2. If all segments in a template are unmapped, their RNAME should be set as ‘*’ and POS as 0.
///    ...
/// 4. Unmapped reads should be stored in the orientation in which they came off the sequencing machine
///    and have their reverse flag bit (0x10) correspondingly unset.
pub fn read_to_bam_record_opts(
    read: &graph_read::Read,
    alignment: &Option<sw::AlignmentPacket>,
    mate_alignment: &Option<sw::AlignmentPacket>,
    strip_qname: bool,
    set_augmented_tags: bool,
) -> bam::record::Record {
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
    let adj_seq: Vec<u8> = read.seq.iter().map(bits_to_ascii).collect();

    let (unmapped_bit, cigar) = match alignment.as_ref() {
        Some(al) => (0, CigarString::from_alignment(&al.alignment, false)), // bool to specify soft clipping
        None => (4, CigarString(vec![])),
    };

    let new_header = fastq::CellrangerFastqHeader::new(read.name.clone());
    let qname = if strip_qname {
        new_header.qname.as_bytes()
    } else {
        new_header.header.as_bytes()
    };

    rec.set(qname, Some(&cigar), &adj_seq, &adj_qual);

    if set_augmented_tags {
        for (tag, value) in new_header.tags {
            if !value.is_empty() {
                rec.push_aux(tag.as_bytes(), bam::record::Aux::String(&value))
                    .unwrap();
            }
        }
    }

    let (tid, pos, mapq) = match (alignment.as_ref(), mate_alignment.as_ref()) {
        (Some(al), Some(mate_al)) => {
            if is_paired && mate_al.ref_idx == al.ref_idx {
                let max_end = max(al.alignment.yend, mate_al.alignment.yend) as i32;
                let min_start = min(al.alignment.ystart, mate_al.alignment.ystart) as i32;
                if al.alignment.ystart <= mate_al.alignment.ystart {
                    rec.set_insert_size((max_end - min_start) as i64);
                } else {
                    rec.set_insert_size(mate_al.alignment.ystart as i64 - al.alignment.yend as i64);
                }
            }
            (al.ref_idx as i32, al.alignment.ystart as i32, 255)
        }
        (Some(al), None) => (al.ref_idx as i32, al.alignment.ystart as i32, 255),
        (None, Some(mate_al)) => (mate_al.ref_idx as i32, mate_al.alignment.ystart as i32, 0),
        (None, None) => (-1, 0, 0),
    };
    rec.set_tid(tid);
    rec.set_pos(pos.into());

    rec.set_mapq(mapq); // this is not really a MAPQ...

    let mate_unmapped_bit = match mate_alignment.as_ref() {
        Some(al) => {
            rec.set_mtid(al.ref_idx as i32);
            rec.set_mpos(al.alignment.ystart as i64);
            0
        }
        None => {
            rec.set_mtid(tid);
            rec.set_mpos(pos.into());
            8
        }
    };

    rec.set_flags(
        is_paired as u16 + // assume always paired
        unmapped_bit +
        mate_unmapped_bit +
        (is_rc as u16) * 16 + // is reverse complemented
        (is_mate_rc as u16) * 32 + // assume if one is RC then the other is not
        {if read.is_first_in_template() {64} else {128}},
    );

    if let Some(al) = alignment.as_ref() {
        rec.push_aux(b"AS", bam::record::Aux::I32(al.alignment.score))
            .unwrap();
        rec.push_aux(b"NM", bam::record::Aux::I32(al.edit_distance() as i32))
            .unwrap();
    }

    rec
}

pub fn read_to_bam_record(
    read: &graph_read::Read,
    alignment: &Option<sw::AlignmentPacket>,
    mate_alignment: &Option<sw::AlignmentPacket>,
) -> bam::record::Record {
    read_to_bam_record_opts(read, alignment, mate_alignment, false, false)
}

/// Add a directory that contains samtools to the PATH.
pub fn check_and_setup_path() {
    use std::env::{join_paths, set_var, split_paths, var_os};
    use std::process::Stdio;

    let mut paths = split_paths(var_os("PATH").unwrap().as_os_str()).collect::<Vec<_>>();
    if Path::new("../../../lib/bin/samtools").exists() {
        // When running under bazel test
        paths.insert(0, "../../../lib/bin".into());
        set_var("PATH", join_paths(paths).unwrap());
    } else if Path::new("../../../bazel-bin/lib/bin/samtools").exists() {
        // When running under cargo test in a tree where bazel build has
        // produced a samtools binary at some point.
        paths.insert(0, "../../../bazel-bin/lib/bin".into());
        set_var("PATH", join_paths(paths).unwrap());
    }
    if let Err(e) = Command::new("samtools")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .spawn()
    {
        if let std::io::ErrorKind::NotFound = e.kind() {
            panic!("`samtools` was not found! Check your PATH!");
        } else {
            panic!("Error {e} while executing samtools");
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam;
    use rust_htslib::bam::Read;
    use std::path::PathBuf;
    use std::{env, fs};

    /// Return the directory for writing test output.
    fn test_outdir() -> PathBuf {
        check_and_setup_path();

        if let Ok(dir) = env::var("TEST_TMPDIR") {
            dir.into()
        } else {
            // TODO: Use TMPDIR.
            const OUTDIR: &str = "test/outputs/bam";
            fs::create_dir_all(OUTDIR).unwrap();
            OUTDIR.into()
        }
    }

    #[test]
    fn test_sort_and_index() {
        let outdir = &test_outdir();
        let test_bam_name = Path::new("test/inputs/test_index.bam");
        let true_sorted_bam_name = Path::new("test/inputs/test_index_sorted.bam");
        let out_bam_name = &outdir.join("test_index_sorted.bam");
        sort_and_index(test_bam_name, out_bam_name);

        {
            let mut true_sorted_bam = bam::Reader::from_path(true_sorted_bam_name).unwrap();
            let mut sorted_bam = bam::Reader::from_path(out_bam_name).unwrap();

            // Don't test the byte arrays. Some versions of samtools seem to omit the SO line.
            //assert_eq!(true_sorted_bam.header().as_bytes(), sorted_bam.header().as_bytes());
            {
                let true_names = true_sorted_bam.header().target_names();
                let sorted_names = sorted_bam.header().target_names();
                assert_eq!(true_names, sorted_names);
            }

            for (read1, read2) in zip(true_sorted_bam.records(), sorted_bam.records()) {
                let r1 = read1.unwrap();
                let r2 = read2.unwrap();
                assert_eq!(r1.qname(), r2.qname());
                assert_eq!(r1.seq().as_bytes(), r2.seq().as_bytes());
            }
        }

        {
            let mut true_sorted_bam = bam::IndexedReader::from_path(true_sorted_bam_name).unwrap();
            let mut sorted_bam = bam::IndexedReader::from_path(out_bam_name).unwrap();

            let _ = true_sorted_bam.fetch((1, 0, 1000));
            let _ = sorted_bam.fetch((1, 0, 1000));
            assert_eq!(sorted_bam.records().count(), 4);

            let _ = true_sorted_bam.fetch((1, 0, 1000));
            let _ = sorted_bam.fetch((1, 0, 1000));
            assert_eq!(
                true_sorted_bam.records().count(),
                sorted_bam.records().count()
            );
        }
    }

    // FIXME -- the test_split1.bam has invalid CIGAR / seq length
    //#[test]
    fn _test_concatenate_bams() {
        let outdir = &test_outdir();
        let merged_bam_name = "test/inputs/test_index.bam";
        let out_bam_name = &outdir.join("test_merged.bam");
        let split_bam_names = [
            Path::new("test/inputs/test_split1.bam"),
            Path::new("test/inputs/test_split2.bam"),
        ];
        if let Err(v) = concatenate_bams(&split_bam_names, out_bam_name) {
            println!("{v}, {}", v.backtrace());
            panic!();
        };

        let mut true_bam = bam::Reader::from_path(merged_bam_name).unwrap();
        let mut merged_bam = bam::Reader::from_path(out_bam_name).unwrap();

        {
            let true_header = true_bam.header();
            let merged_header = merged_bam.header();

            for (true_name, other_name) in
                zip(true_header.target_names(), merged_header.target_names())
            {
                assert_eq!(true_name, other_name);
            }
        }

        for (read1, read2) in zip(true_bam.records(), merged_bam.records()) {
            let r1 = read1.unwrap();
            let r2 = read2.unwrap();
            assert_eq!(r1.tid(), r2.tid());
            assert_eq!(r1.qname(), r2.qname());
            assert_eq!(r1.seq().as_bytes(), r2.seq().as_bytes());
        }
    }
}
