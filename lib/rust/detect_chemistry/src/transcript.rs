//
// Copyright (c) 2017 10x Genomics, Inc. All rights reserved.
//

use std::collections::{BinaryHeap,HashMap,HashSet};
use std::io::{Read, Seek, BufRead};

use regex::Regex;

use bio::io::{fasta};
use bio::alphabets::dna::revcomp;

pub struct Transcript {
    pub chr: String,             // reference id / chromosome
    pub reverse: bool,           // reverse strand
    pub exons: BinaryHeap<(i64, i64)>,  // exon intervals
}

lazy_static! {
    static ref TRANSCRIPT_ID_REGEX: Regex = Regex::new(r#"transcript_id\s+"(.*?)""#).expect("Failed to compile regex");
    static ref GENE_ID_REGEX: Regex = Regex::new(r#"gene_id\s+"(.*?)""#).expect("Failed to compile regex");
}

fn get_gtf_ids<'a>(text: &'a str) -> (&'a str, &'a str) {
    (&TRANSCRIPT_ID_REGEX.captures(text).expect("Missing transcript_id").get(1).unwrap().as_str(),
    &GENE_ID_REGEX.captures(text).expect("Missing gene_id").get(1).unwrap().as_str())
}

pub fn parse_gtf<R: BufRead>(reader: R, take_first: bool) -> HashMap<String, Transcript> {
    let mut transcripts: HashMap<String, Transcript> = HashMap::new();
    let mut genes: HashSet<String> = HashSet::new();

    // Note: rust-bio gff/gtf parsing still breaks on comments; do it ourselves
    for line in reader.lines().map(|x| x.expect("Failed to read line"))
        .filter(|x| !x.starts_with("#")) {
            let mut row = line.split("\t");

            let chr = row.nth(0).expect("GTF missing 'chrom' column");

            let feature = row.nth(1).expect("GTF missing 'feature_type' column");
            if feature != "exon" {
                continue;
            }

            let start = row.next().expect("GTF missing 'start' column")
                .parse::<i64>().expect("Couldn't interpret 'start' as integer'");
            let end = row.next().expect("GTF missing 'end' column")
                .parse::<i64>().expect("Couldn't interpret 'end' as integer'");
            let strand = row.nth(1).expect("GTF missing 'strand' column");

            let (tx_id, gene_id) = get_gtf_ids(row.nth(1).expect("GTF missing 'info' column"));

            // Only take the 1st encountered transcript per gene
            if take_first && !transcripts.contains_key(tx_id) && genes.contains(gene_id) {
                continue;
            }

            let mut tx = transcripts.entry(tx_id.to_owned()).or_insert(Transcript {
                chr: chr.to_owned(),
                reverse: strand == "-",
                exons: BinaryHeap::new(),
            });

            tx.exons.push((start, end));
            genes.insert(gene_id.to_owned());
        }
    transcripts
}

pub fn fetch_transcript<R: Read + Seek>(reader: &mut fasta::IndexedReader<R>,
                                        transcript: Transcript) -> Vec<u8> {
    // Get transcript sequence from genome FASTA
    let mut exons = transcript.exons.into_sorted_vec();

    if transcript.reverse {
        exons.reverse();
    }

    let tx_len = exons.iter().fold(0, |x, &(a,b)| x + (b - a));
    let mut seq = Vec::with_capacity(tx_len as usize);

    for (start, end) in exons {
        let mut exon_seq = Vec::with_capacity((end-start) as usize);

	let _ = reader.fetch(&transcript.chr, start as u64, end as u64);
        match reader.read(&mut exon_seq) {
            Ok(_) => {},
            Err(e) => {
                println!("WARNING: Failed to fetch sequence from genome FASTA file. \nIgnoring {:?}",e);
                return Vec::new();
            }
        }

        if transcript.reverse {
            seq.append(&mut revcomp(&exon_seq));
        } else {
            seq.append(&mut exon_seq);
        }
    }
    seq
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_parse_gtf() {
        use transcript::parse_gtf;
        let gtf_str = b"#test_comment\n1\tfoo\texon\t1\t100\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\";" as &[u8];
        let transcripts = parse_gtf(gtf_str, true);
        assert!(transcripts.contains_key("t1"));
        assert_eq!(transcripts.get("t1").unwrap().chr, "1");
        assert_eq!(transcripts.get("t1").unwrap().reverse, false);
        assert_eq!(transcripts.get("t1").unwrap().exons.len(), 1);
    }

    #[test]
    fn test_fetch() {
        use transcript::{fetch_transcript, Transcript};
        use bio::io::fasta;
        use std::io::Cursor;
        use std::collections::BinaryHeap;
        use bio::alphabets::dna::revcomp;
        use tests::random_seq;

        let contig_len = 1000usize;
        let seq = random_seq(contig_len);
        let fa_str = format!(">1\n{}", String::from_utf8(seq.clone()).unwrap());
        let faidx_str = format!("1\t{}\t3\t{}\t{}", &contig_len, &contig_len, &(1+contig_len));

        let mut reader = fasta::IndexedReader::new(Cursor::new(fa_str.as_bytes()),
                                                   Cursor::new(faidx_str.as_bytes())).unwrap();
        let mut exons = BinaryHeap::new();
        exons.push((100i64, 200i64));
        exons.push((300i64, 400i64));

        // Test boring transcript
        let tx = Transcript { chr: "1".to_owned(), reverse: false, exons: exons.clone() };
        let tx_seq = fetch_transcript(&mut reader, tx);

        let mut exp_seq = seq[100..200].to_owned();
        exp_seq.extend(seq[300..400].iter());
        assert_eq!(tx_seq, exp_seq);

        // Test rev strand
        let tx = Transcript { chr: "1".to_owned(), reverse: true, exons: exons };
        let tx_seq = fetch_transcript(&mut reader, tx);

        assert_eq!(tx_seq, revcomp(&exp_seq));
    }
}
