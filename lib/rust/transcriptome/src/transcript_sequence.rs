use crate::transcriptome::Exon;
use anyhow::Result;
use bio::io::fasta::IndexedReader;
use bio_types::strand::ReqStrand;
use itertools::process_results;
use std::io::{Read, Seek};

/// Construct a transcript sequence from a FASTA file, a strand, and a set of intervals
/// The intervals must be ordered according to the order of the intervals on the genome
pub(crate) fn get_transcript_sequence(
    reader: &mut IndexedReader<impl Read + Seek>,
    contig: &str,
    exons: &[Exon],
    strand: ReqStrand,
) -> Result<Vec<u8>> {
    let size: usize = exons.iter().map(Exon::len).sum::<u64>() as usize;
    let mut seq = Vec::with_capacity(size);
    for exon in exons {
        reader.fetch(contig, exon.start, exon.end)?;
        process_results(reader.read_iter()?, |iter| seq.extend(iter))?;
    }
    assert_eq!(seq.len(), size);
    Ok(match strand {
        ReqStrand::Forward => seq,
        ReqStrand::Reverse => revcomp_vec(seq),
    })
}

/// Reverse complement DNA.
fn revcomp_slice(dna: &mut [u8]) {
    for i in 0..(dna.len() + 1) / 2 {
        let j = dna.len() - 1 - i;
        let fst = bio::alphabets::dna::complement(dna[i]);
        let lst = bio::alphabets::dna::complement(dna[j]);
        dna[i] = lst;
        dna[j] = fst;
    }
}

/// Reverse complement DNA.
fn revcomp_vec(mut dna: Vec<u8>) -> Vec<u8> {
    revcomp_slice(&mut dna[..]);
    dna
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::Transcriptome;
    use std::collections::HashMap;
    use std::fs::File;
    use std::io::BufReader;

    fn check_rc(fwd: &[u8], rev: &[u8]) {
        assert_eq!(revcomp_vec(Vec::from(fwd)), rev);
    }

    #[test]
    fn test_rc_slice() {
        check_rc(b"ACGT", b"ACGT");
        check_rc(b"A", b"T");
        check_rc(b"GC", b"GC");
        check_rc(b"GCA", b"TGC");
        check_rc(b"TTT", b"AAA");
        check_rc(b"GTGTT", b"AACAC");
    }

    /// Run this test manually to check that transcript sequences produced by Rust from GTF + genome fasta
    /// match what gencode produces.
    #[ignore]
    #[test]
    fn loader() -> Result<()> {
        let rdr = File::open("/Users/patrick/refdata_cellranger/GRCh38-3.0.0/genes/genes.gtf")?;
        let txome = Transcriptome::from_reader(BufReader::new(rdr))?;

        let mut fa = bio::io::fasta::IndexedReader::from_file(
            &"/Users/patrick/refdata_cellranger/GRCh38-3.0.0/fasta/genome.fa",
        )
        .unwrap();

        let txs = "/Users/patrick/code/rust-pseudoaligner/test/gencode.v28.transcripts.fa";
        let gencode_txs = bio::io::fasta::Reader::from_file(txs).unwrap();

        let mut tx_seqs = HashMap::new();

        for r in gencode_txs.records() {
            let r = r?;

            let id = r.id();
            let tx_id = id.split('.').next().unwrap();
            tx_seqs.insert(tx_id.to_string(), r.seq().to_vec());
        }

        let mut not_found = 0;
        for tx in &txome.transcripts {
            if let Some(real_seq) = tx_seqs.get(&tx.id) {
                let loaded_tx_seq =
                    get_transcript_sequence(&mut fa, &tx.chrom, &tx.exons, tx.strand)?;
                assert!(
                    &loaded_tx_seq == real_seq,
                    "tx_id: {} gencode: {} me: {}",
                    &tx.id,
                    std::str::from_utf8(real_seq).unwrap(),
                    std::str::from_utf8(&loaded_tx_seq).unwrap()
                );
            } else {
                not_found += 1;
            }
        }

        println!(
            "tried to validate {}, not found: {}",
            &txome.transcripts.len(),
            not_found
        );

        Ok(())
    }
}
