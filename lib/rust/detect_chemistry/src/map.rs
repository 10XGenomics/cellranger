//
// Copyright (c) 2017 10x Genomics, Inc. All rights reserved.
//

use std::io::{Read};
use itertools::Itertools;

use bio::io::fastq;

use debruijn::{Kmer, Vmer};
use debruijn::dna_string::DnaString;

use index;

#[derive(Copy, Clone, Debug, Default, Serialize)]
pub struct Metrics {
    total_reads: u64,
    sense_reads: u64,
    antisense_reads: u64,
    ambiguous_reads: u64,
    mapped_reads: u64,
}

fn count_read(metrics: &mut Metrics, sense: bool, anti: bool) {
    metrics.antisense_reads += (anti && !sense) as u64;
    metrics.sense_reads += (sense && !anti) as u64;
    metrics.ambiguous_reads += (sense && anti) as u64;
    metrics.mapped_reads += (sense || anti) as u64;
    metrics.total_reads += 1;
}

pub fn add_metrics(x: &mut Metrics, y: &Metrics) {
    x.total_reads += y.total_reads;
    x.sense_reads += y.sense_reads;
    x.antisense_reads += y.antisense_reads;
    x.ambiguous_reads += y.ambiguous_reads;
    x.mapped_reads += y.mapped_reads;
}

pub fn map_reads<K: Kmer, I: index::KmerPresenceQuery<K> + ?Sized, R: Read> (index: &I, fq_reader: fastq::Reader<R>, min_kmers: usize, skip_bases: usize, initial_reads: usize, interleaved: bool, read_type: Option<String>) -> Metrics {
    let mut metrics: Metrics = Default::default();

    let mut start_at = 0;
    if interleaved {
        match read_type.expect("read_type must be provided if interleaved == true").as_str() {
            "R1" => { start_at = 0;},
            "R2" => { start_at = 1;},
            _ => { panic!("Unexpected read type"); },
        }
    };

    let every = match interleaved {
        true => 2,
        false => 1,
    };

    for rec in fq_reader.records()
        .map(|x| x.expect("Failed to parse FASTQ record"))
        .take(initial_reads)
        .skip(start_at)
        .step(every) {

            let dna_string = DnaString::from_acgt_bytes(&rec.seq());

            let mut sense = 0usize;
            let mut anti = 0usize;

            for kmer in dna_string.iter_kmers().step(1+skip_bases) {
                sense += index.contains(&kmer) as usize;
                anti += index.contains(&kmer.rc()) as usize;
                if (sense >= min_kmers) || (anti >= min_kmers) {
                    break;
                }
            }
            count_read(&mut metrics, sense >= min_kmers, anti >= min_kmers);
        }
    metrics
}


#[cfg(test)]
mod tests {
    fn make_test_fastq(read_seqs: &Vec<String>) -> String {
        let mut fq_lines = Vec::new();
        for read_seq in read_seqs {
            fq_lines.push(String::from("@x"));
            fq_lines.push(read_seq.clone());
            fq_lines.push(String::from("+"));
            fq_lines.push("I".repeat(50));
        }
        fq_lines.join("\n")
    }

    #[test]
    fn test_map_reads() {
        use index;
        use std::io::{Cursor};
        use bio::io::{fasta,fastq};
        use tests::random_seq_rng;
        use debruijn::kmer::{IntKmer};
        use rand::{SeedableRng, StdRng};
        use map::map_reads;
        use bio::alphabets::dna::revcomp;
        use std;

        type MyKmer = IntKmer<u64>;

        let seed: &[_] = &[1, 2, 3, 4];
        let mut rng: StdRng = SeedableRng::from_seed(seed);

        // Generate transcript sequence and index
        let tx_seq = random_seq_rng(100, &mut rng);
        let fa_str = format!(">1\n{}", String::from_utf8(tx_seq.clone()).unwrap());

        let reader = fasta::Reader::new(Cursor::new(fa_str.as_bytes()));

        let index: index::BBHashKmerIndex<MyKmer> = index::index_transcripts_mphf(reader, 1000, 0, 2.0);

        // Generate fastq
        let mut sense_fq_reads = Vec::new();
        for i in 0..50 {
            let read_seq = String::from_utf8(tx_seq.clone()[i..(i+50)].to_owned()).unwrap();
            sense_fq_reads.push(read_seq);
        }
        let sense_fq = make_test_fastq(&sense_fq_reads);
        let sense_fq_reader = fastq::Reader::new(Cursor::new(sense_fq.as_bytes()));

        let metrics = map_reads(&index, sense_fq_reader, 1, 0, std::usize::MAX, false, None);
        assert_eq!(metrics.total_reads, 50);
        assert_eq!(metrics.sense_reads, 50);
        assert_eq!(metrics.mapped_reads, 50);
        assert_eq!(metrics.antisense_reads, 0);
        assert_eq!(metrics.ambiguous_reads, 0);

        // Generate fastq
        let mut anti_fq_reads = Vec::new();
        for i in 0..50 {
            let piece = tx_seq.clone()[i..(i+50)].to_owned();
            let read_seq = String::from_utf8(revcomp(&piece).to_vec()).unwrap();
            anti_fq_reads.push(read_seq);
        }
        let anti_fq = make_test_fastq(&anti_fq_reads);
        let anti_fq_reader = fastq::Reader::new(Cursor::new(anti_fq.as_bytes()));

        let metrics = map_reads(&index, anti_fq_reader, 1, 0, std::usize::MAX, false, None);
        assert_eq!(metrics.total_reads, 50);
        assert_eq!(metrics.sense_reads, 0);
        assert_eq!(metrics.mapped_reads, 50);
        assert_eq!(metrics.antisense_reads, 50);
        assert_eq!(metrics.ambiguous_reads, 0);
    }
}
