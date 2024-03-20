// Code for reading and testing bam and fastq files.
//
// Argument max_reads: stop after loading this many reads.  If -1, do not cap.
// As used in VDJ assembly, this is somewhat dubious, for one thing because we cap
// the number of reads before sorting by UMI.

use debruijn::dna_string::DnaString;
use io_utils::open_maybe_compressed;
use std::io::{BufRead, BufReader};
use std::path::Path;
use string_utils::TextUtils;

pub struct BamReadOptions {
    pub keep_seqs: bool,     // keep sequences
    pub keep_quals: bool,    // keep quality scores
    pub keep_readname: bool, // keep readname
    pub n_free_tail: usize,  // reads must have this many N-free bases at
                             // end before reversing else we trim
}

impl BamReadOptions {
    pub fn new() -> BamReadOptions {
        BamReadOptions {
            keep_seqs: true,
            keep_quals: true,
            keep_readname: true,
            n_free_tail: 0,
        }
    }
}

impl Default for BamReadOptions {
    fn default() -> Self {
        Self::new()
    }
}

// Extract barcode and umi from a 10x readname.
// This might be specific to certain pipelines.

fn fetch_bc_and_umi(s: &[u8], bc: &mut String, umi: &mut String) {
    let mut field = 0;
    bc.clear();
    umi.clear();
    // ◼ could be more efficient to go backwards
    for c in s {
        if *c == b'|' {
            field += 1;
            if field > 42 {
                break;
            }
        } else if field == 30 {
            umi.push(*c as char);
        } else if field == 42 {
            if *c as char == ' ' {
                break;
            }
            bc.push(*c as char);
        }
    }
}

// Gather read data from a 10x bam file, storing as
// { ( barcode, { ( umi, seq, qual, readname, flags ) } }.

// Read one read from a 10x fastq file.  Reverse complement if requested.
// Result is (seq,qual,bc,umi,readname).  There is an extra argument "line"
// to avoid reallocating it.
//
// ◼ Need to document what is going on here with the many calls to read_line
// ◼ and pop.  At first glance it does not make sense.  And you can see from
// ◼ "ends_with("+\n")" that we're doing something wonky.  If you add temporary
// ◼ print lines after each readline you can see what is happening.
//
// ◼ Ideally this code would do a full check for fastq validity.
#[allow(clippy::too_many_arguments)]
fn read_one_read_from_10x_fastq<R: BufRead>(
    r: &mut R,
    rc: bool,
    seq: &mut Vec<u8>,
    qual: &mut Vec<u8>,
    bc: &mut String,
    umi: &mut String,
    readname: &mut String,
    line: &mut String,
    opt: &BamReadOptions,
) -> Option<()> {
    line.clear();
    r.read_line(line).ok()?;
    if line.is_empty() {
        return None;
    }
    line.pop();
    assert!(
        line.starts_with('@'),
        "Invalid readname line in fastq file: \"{line}\""
    );
    if opt.keep_readname {
        *readname = line.after("@").to_string();
    }
    fetch_bc_and_umi(line.as_bytes(), bc, umi);
    line.clear();
    r.read_line(line).ok()?;
    if line.is_empty() {
        return None;
    }
    r.read_line(line).ok()?;
    assert!(
        line.ends_with("+\n"),
        "Failed to find plus line in fastq file."
    );
    if line.is_empty() {
        return None;
    }
    line.pop();
    line.pop();
    line.pop();
    if opt.keep_seqs {
        seq.clear();
        seq.extend(line.as_bytes());
    }

    // Do some minimal filtering to remove Ns in the sequence.  The motivation is
    // if there are enough Ns, downstream code that kmerizes reads will use too
    // much memory.
    // ◼ This is clumsy.  Find a better filter.

    if opt.keep_seqs && opt.n_free_tail > 0 {
        while seq.len() >= opt.n_free_tail {
            let n = seq.len();
            let mut trimmed = false;
            for j in n - opt.n_free_tail..n {
                if seq[j] == b'N' {
                    seq.truncate(j);
                    trimmed = true;
                    break;
                }
            }
            if !trimmed {
                break;
            }
        }
    }

    // Get the quality scores.

    line.clear();
    r.read_line(line).ok()?;
    if line.is_empty() {
        return None;
    }
    line.pop();
    if opt.keep_quals {
        qual.clear();
        qual.extend(line.as_bytes().iter().take(seq.len()).map(|&i| i - 33));
    }
    if opt.keep_seqs && rc {
        seq.reverse();
        for v in seq {
            if *v == b'A' {
                *v = b'T';
            } else if *v == b'C' {
                *v = b'G';
            } else if *v == b'G' {
                *v = b'C';
            } else if *v == b'T' {
                *v = b'A';
            }
        }
        qual.reverse();
    }
    Some(())
}

/*
// Read a barcode's worth of reads from a 10x fastq file, for single-ended reads.
// This keeps a bunch of data as an "internal state", because the code has to read
// into the next barcode before knowing it's done.  This code is slower than
// expected, for unknown reasons.  See alt implementation of fastq1, below.

fn read_one_barcode_of_reads_from_10x_fastq1<R: BufRead>(
    // where data are read from
    r1: &mut R,
    // max reads to keep, or -1 for all
    max_reads: i64,
    opt: &BamReadOptions,
    // output data: (barcode, { (umi,   seq,      qual,   readname,flags) }, actual reads )
    read_data: &mut (
        String,
        Vec<(String, DnaString, Vec<u8>, String, u16)>,
        usize,
    ),
    // internal state
    seq: &mut Vec<u8>,
    qual: &mut Vec<u8>,
    bc: &mut String,
    umi: &mut String,
    readname: &mut String,
    line: &mut String,
) -> Option<()> {
    if *bc == "done".to_string() {
        return None;
    }
    read_data.0.clear();
    read_data.1.clear();
    read_data.2 = 0;
    loop {
        if bc.is_empty() || !read_data.1.is_empty() {
            if read_one_read_from_10x_fastq(r1, true, seq, qual, bc, umi, readname, line, opt)
                .is_none()
            {
                if read_data.1.is_empty() {
                    return None;
                } else {
                    *bc = "done".to_string();
                    return Some(());
                }
            }
            if read_data.0.is_empty() {
                read_data.0 = bc.clone();
            } else if *bc != read_data.0 {
                return Some(());
            }
        }
        let flags: u16 = 16_u16  // sequence is reverse complemented
            + 128_u16; // the last segment in the template
        if max_reads < 0 || (read_data.1.len() as i64) < max_reads {
            if read_data.1.is_empty() {
                read_data.0 = bc.clone();
            }
            (read_data.1).push((
                umi.to_string(),
                DnaString::from_acgt_bytes(seq),
                qual.clone(),
                readname.clone(),
                flags,
            ));
        }
        read_data.2 += 1;
    }
}
*/

/*
// NOT DOING THIS BECAUSE IT'S SLOWER
// ◼ Need to understand why.
pub fn parse_fastq1( fastq1: &str,
    read_data: &mut Vec<(String,Vec<(String,DnaString,Vec<u8>,String,u16)>)>,
    max_reads: i64, opt: &BamReadOptions ) {
    read_data.clear();
    let mut r1 = BufReader::new(open_maybe_compressed(fastq1));
    let mut readname = String::new();
    let (mut bc, mut umi) = ( String::new(), String::new() );
    let (mut seq, mut qual) = ( Vec::<u8>::new(), Vec::<u8>::new() );
    let mut line = String::new();
    loop {
        let mut read_datax
            = ( String::new(), Vec::<(String,DnaString,Vec<u8>,String,u16)>::new() );
        if read_one_barcode_of_reads_from_10x_fastq1( &mut r1, max_reads,
            &mut read_datax, &mut seq, &mut qual, &mut bc, &mut umi, &mut readname,
            &mut line, opt ).is_none() {
            break;
        }
        read_data.push(read_datax);
    }
}
*/

pub type ReadData = (
    String,
    Vec<(String, DnaString, Vec<u8>, String, u16)>,
    usize,
);

fn parse_fastq1(
    fastq1: &Path,
    read_data: &mut Vec<ReadData>,
    max_reads: i64,
    opt: &BamReadOptions,
) {
    read_data.clear();
    let mut r1 = BufReader::new(open_maybe_compressed(fastq1));
    let mut readname = String::new();
    let (mut bc, mut umi) = (String::new(), String::new());
    let (mut seq, mut qual) = (Vec::<u8>::new(), Vec::<u8>::new());
    let mut this_bc = String::new();
    let mut count = 0_i64;
    let mut line = String::new();
    loop {
        if read_one_read_from_10x_fastq(
            &mut r1,
            true,
            &mut seq,
            &mut qual,
            &mut bc,
            &mut umi,
            &mut readname,
            &mut line,
            opt,
        )
        .is_none()
        {
            break;
        }
        if bc != this_bc {
            let x = Vec::<(String, DnaString, Vec<u8>, String, u16)>::new();
            read_data.push((bc.to_string(), x, 0));
            this_bc = bc.to_string();
            count = 0_i64;
        }
        let n = read_data.len();
        let flags: u16 = 16_u16  // sequence is reverse complemented
            + 128_u16; // the last segment in the template
        if max_reads < 0 || count < max_reads {
            (read_data[n - 1].1).push((
                umi.to_string(),
                DnaString::from_acgt_bytes(&seq),
                qual.clone(),
                readname.clone(),
                flags,
            ));
            count += 1;
        }
        read_data[n - 1].2 += 1;
    }
}

/*
// Read a barcode's worth of reads from a 10x fastq file, for paired-end reads.
// See comments for the single-ended version.

fn read_one_barcode_of_reads_from_10x_fastq2<R: BufRead>(
    // where data are read from
    r1: &mut R,
    r2: &mut R,
    // max reads to keep, or -1 for all
    max_reads: i64,
    opt: &BamReadOptions,
    // output data: (barcode, { (umi,   seq,      qual,   readname,flags) }, actual reads )
    read_data: &mut (
        String,
        Vec<(String, DnaString, Vec<u8>, String, u16)>,
        usize,
    ),
    // internal state
    seq: &mut Vec<u8>,
    qual: &mut Vec<u8>,
    bc: &mut String,
    umi: &mut String,
    readname: &mut String,
    line: &mut String,
) -> Option<()> {
    if *bc == "done".to_string() {
        return None;
    }
    read_data.0.clear();
    read_data.1.clear();
    read_data.2 = 0;
    loop {
        if bc.is_empty() || !read_data.1.is_empty() {
            if read_one_read_from_10x_fastq(r1, false, seq, qual, bc, umi, readname, line, opt)
                .is_none()
            {
                if read_data.1.is_empty() {
                    return None;
                } else {
                    *bc = "done".to_string();
                    return Some(());
                }
            }
            if read_data.0.is_empty() {
                read_data.0 = bc.clone();
            } else if *bc != read_data.0 {
                return Some(());
            }
        }
        let flags: u16 = 1_u16   // have multiple segments
            +  32_u16   // seg of next segment reverse complemented
            +  64_u16; // the first segment in the template
        if max_reads < 0 || (read_data.1.len() as i64) < max_reads {
            if read_data.1.is_empty() {
                read_data.0 = bc.clone();
            }
            (read_data.1).push((
                umi.to_string(),
                DnaString::from_acgt_bytes(seq),
                qual.clone(),
                readname.clone(),
                flags,
            ));
        }
        read_data.2 += 1;
        if read_one_read_from_10x_fastq(r2, true, seq, qual, bc, umi, readname, line, opt).is_none()
        {
            panic!("inconsistent reading from fastq files 1 and 2");
        }
        let flags: u16 = 1_u16  // have multiple segments
            +  16_u16  // sequence is reverse complemented
            + 128_u16; // the last segment in the template
        if max_reads < 0 || (read_data.1.len() as i64) < max_reads {
            (read_data.1).push((
                umi.to_string(),
                DnaString::from_acgt_bytes(seq),
                qual.clone(),
                readname.clone(),
                flags,
            ));
        }
        read_data.2 += 1;
    }
}
*/

/*
// NOT DOING THIS BECAUSE IT'S SLOWER
// ◼ Need to understand why.  Test case is pogo.
pub fn parse_fastq2( fastq1: &str, fastq2: &str,
    read_data: &mut Vec<(String,Vec<(String,DnaString,Vec<u8>,String,u16)>)>,
    max_reads: i64, opt: &BamReadOptions ) {
    read_data.clear();
    let mut r1 = BufReader::new(open_maybe_compressed(fastq1));
    let mut r2 = BufReader::new(open_maybe_compressed(fastq2));
    let mut readname = String::new();
    let (mut bc, mut umi) = ( String::new(), String::new() );
    let (mut seq, mut qual) = ( Vec::<u8>::new(), Vec::<u8>::new() );
    let mut line = String::new();
    loop {
        let mut read_datax
            = ( String::new(), Vec::<(String,DnaString,Vec<u8>,String,u16)>::new() );
        if read_one_barcode_of_reads_from_10x_fastq2( &mut r1, &mut r2, max_reads,
            opt, &mut read_datax, &mut seq, &mut qual, &mut bc, &mut umi,
            &mut readname, &mut line ).is_none() {
            break;
        }
        read_data.push(read_datax);
    }
}
*/

fn parse_fastq2(
    fastq1: &Path,
    fastq2: &Path,
    read_data: &mut Vec<ReadData>,
    max_reads: i64,
    opt: &BamReadOptions,
) {
    read_data.clear();
    let mut r1 = BufReader::new(open_maybe_compressed(fastq1));
    let mut r2 = BufReader::new(open_maybe_compressed(fastq2));
    let mut readname = String::new();
    let (mut bc, mut umi) = (String::new(), String::new());
    let (mut seq, mut qual) = (Vec::<u8>::new(), Vec::<u8>::new());
    let mut this_bc = String::new();
    let mut count = 0_i64;
    let mut line = String::new();
    loop {
        if read_one_read_from_10x_fastq(
            &mut r1,
            false,
            &mut seq,
            &mut qual,
            &mut bc,
            &mut umi,
            &mut readname,
            &mut line,
            opt,
        )
        .is_none()
        {
            break;
        }
        if bc != this_bc {
            let x = Vec::<(String, DnaString, Vec<u8>, String, u16)>::new();
            read_data.push((bc.to_string(), x, 0));
            this_bc = bc.to_string();
            count = 0_i64;
        }
        let n = read_data.len();
        let flags: u16 = 1_u16   // have multiple segments
            +  32_u16   // seg of next segment reverse complemented
            +  64_u16; // the first segment in the template
        if max_reads < 0 || count < max_reads {
            (read_data[n - 1].1).push((
                umi.to_string(),
                DnaString::from_acgt_bytes(&seq),
                qual.clone(),
                readname.clone(),
                flags,
            ));
        }
        read_data[n - 1].2 += 1;
        if read_one_read_from_10x_fastq(
            &mut r2,
            true,
            &mut seq,
            &mut qual,
            &mut bc,
            &mut umi,
            &mut readname,
            &mut line,
            opt,
        )
        .is_none()
        {
            break;
        }
        let n = read_data.len();
        let flags: u16 = 1_u16  // have multiple segments
            +  16_u16  // sequence is reverse complemented
            + 128_u16; // the last segment in the template
        if max_reads < 0 || count < max_reads {
            (read_data[n - 1].1).push((
                umi.to_string(),
                DnaString::from_acgt_bytes(&seq),
                qual.clone(),
                readname.clone(),
                flags,
            ));
            count += 2;
        }
        read_data[n - 1].2 += 1;
    }
}

#[allow(clippy::type_complexity)]
pub fn parse_fastq(
    fastq1: &Path,
    fastq2_option: Option<&Path>,
    read_data: &mut Vec<(
        String,
        Vec<(String, DnaString, Vec<u8>, String, u16)>,
        usize,
    )>,
    max_reads: i64,
    opt: &BamReadOptions,
) {
    if let Some(fastq2) = fastq2_option {
        parse_fastq2(fastq1, fastq2, read_data, max_reads, opt);
    } else {
        parse_fastq1(fastq1, read_data, max_reads, opt);
    }
}

pub fn parse_fastq_as_bam(
    bam_name: &Path,
    read_data: &mut Vec<ReadData>,
    max_reads: i64,
    opt: &BamReadOptions,
) {
    let bam_dir = bam_name.parent().unwrap();
    let fastq = bam_dir.parent().unwrap().join("fastq");
    let ch = bam_dir.file_name().unwrap().to_str().unwrap().before(".");
    let name1 = fastq.join(format!("{ch}.reads_1.fastq.lz4"));
    let name2 = fastq.join(format!("{ch}.reads_2.fastq.lz4"));
    parse_fastq(&name1, Some(&name2), read_data, max_reads, opt);
}
