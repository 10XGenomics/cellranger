// TODO: fix these.
#![allow(clippy::needless_range_loop)]

use lz4;
use rust_htslib::bam;
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::Record;
use std::fs::File;
use std::io::{BufWriter, Write};
use string_utils::TextUtils;

pub fn write_sam_record_simple(
    x: &bam::Record,
    u: &str, // corrected umi
    tignames: &[String],
    writer: &mut BufWriter<&mut lz4::Encoder<BufWriter<File>>>,
) {
    // Determine if mate is mapped.

    let mate_mapped = (x.flags() & 0x8) == 0;

    // Set reference sequence name field.

    let mut tigname = "*".to_string();
    if x.tid() >= 0 {
        tigname.clone_from(&tignames[x.tid() as usize]);
    }

    // Set mate reference sequence name field.

    let mut mtid = "*".to_string();
    if x.tid() != x.mtid() {
        mtid.clone_from(&tignames[x.mtid() as usize]);
    } else if mate_mapped {
        mtid = "=".to_string();
    }

    // Punt if no NM tag, but read is aligned.
    // ◼ This seems flaky.

    if x.aux(b"NM").is_err() && tigname != "*" {
        return;
    }

    // Fetch tags.

    let extract_aux_string = |r: &Record, tag: &str| -> String {
        match r.aux(tag.as_bytes()) {
            Ok(Aux::String(s)) => format!("\t{tag}:Z:{s}"),
            Ok(a) => panic!("Unexpected aux {a:?} for {tag}"),
            Err(_) => String::default(),
        }
    };

    let bc = extract_aux_string(x, "BC");

    let cb = match x.aux(b"CB") {
        Ok(Aux::String(s)) => {
            format!("\tCB:Z:{}", if s.contains(' ') { s.before(" ") } else { s })
        }
        Ok(a) => panic!("Unexpected aux {a:?} for CB"),
        Err(_) => String::default(),
    };

    let cr = extract_aux_string(x, "CR");

    let cy = extract_aux_string(x, "CY");

    let qt = extract_aux_string(x, "QT");

    let ub = format!("\tUB:Z:{u}");

    let ur = extract_aux_string(x, "UR");

    let uy = extract_aux_string(x, "UY");

    let extract_aux_i32 = |r: &Record, tag: &str| -> String {
        match r.aux(tag.as_bytes()) {
            Ok(Aux::I32(i)) => format!("\t{tag}:i:{i}"),
            Ok(a) => panic!("Unexpected aux {a:?} for {tag}"),
            Err(_) => String::default(),
        }
    };

    let nm = extract_aux_i32(x, "NM");

    let as_tag = extract_aux_i32(x, "AS");

    let mut quals = x.qual().to_vec();
    for i in 0..quals.len() {
        quals[i] += 33;
    }

    let cigar = x.cigar();
    let cigar = if cigar.0.is_empty() {
        "*".to_string()
    } else {
        cigar.to_string()
    };

    // Get read position.  This is super annoying.  The field x.pos() is
    // zero-based, but in a sam record, pos is one-based, and zero is reserved for
    // the case where the read is not aligned.

    let mut pos = x.pos();
    if tigname != "*" {
        pos += 1;
    }

    // Get mate read position.  Ditto above comments.

    let mut mpos = x.mpos();
    if mate_mapped {
        mpos += 1;
    }

    // Write sam record.

    writer
        .write_fmt(format_args!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\
             {}{}{}{}{}{}{}{}{}{}\n",
            String::from_utf8(x.qname().to_vec()).unwrap(),
            x.flags(),
            tigname,
            pos,
            x.mapq(),
            cigar,
            &mtid,
            mpos,
            x.insert_size(),
            String::from_utf8(x.seq().as_bytes()).unwrap(),
            String::from_utf8(quals).unwrap(),
            bc,
            cb,
            cr,
            cy,
            qt,
            ub,
            ur,
            uy,
            as_tag,
            nm
        ))
        .unwrap();
}
