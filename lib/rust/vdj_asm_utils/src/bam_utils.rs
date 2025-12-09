// ◼ There are some fundamental problems with rust-htslib that affect our ability
// ◼ to debug our code:
// ◼ (a) the lack of code to print a bam record
// ◼ (b) lack of code to validate a bam record and say what’s wrong with it
// ◼ (c) bad exit handling by rust-htslib -- it doesn’t panic upon failure so you
// ◼     don’t get a traceback
// ◼ (d) inability to test a change to rust-htslib without doing cargo clean on
// ◼     the project that uses it.
#![expect(missing_docs)]

use crate::{fastq, graph_read, sw};
use debruijn::bits_to_ascii;
use rust_htslib::bam::record::CigarString;
use rust_htslib::bam::{self};
use std::cmp::{max, min};

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

#[cfg(test)]
pub(super) fn read_to_bam_record(
    read: &graph_read::Read,
    alignment: &Option<sw::AlignmentPacket>,
    mate_alignment: &Option<sw::AlignmentPacket>,
) -> bam::record::Record {
    read_to_bam_record_opts(read, alignment, mate_alignment, false, false)
}
