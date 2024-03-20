//!
//! Create read data object that the assembler needs from a
//! vector of rna reads. This is likely a patch until the more structured
//! form of data makes its way deep into the assembler
//!
use bitflags::bitflags;
use cr_types::chemistry::ChemistryDef;
use cr_types::rna_read::RnaRead;
use debruijn::dna_string::DnaString;
use fastq_set::metric_utils::ILLUMINA_QUAL_OFFSET;
use fastq_set::sseq::HammingIterOpt;
use fastq_set::WhichEnd;
use itertools::Itertools;
use metric::SimpleHistogram;

bitflags! {
    struct ReadFlags: u16 {
        const MULTI_SEGMENT = 0b0000_0001; // Have multiple segments i.e read1 and read2 in our case
        const REV_COMP      = 0b0001_0000; // Sequence in reverse complemented
        const REV_COMP_MATE = 0b0010_0000; // Sequence of the next segment is reverse complemented
        const FIRST_SEGMENT = 0b0100_0000;
        const LAST_SEGMENT  = 0b1000_0000;
        // Derived compond read flags
        const SINGLE_END_READ_FLAG  = Self::REV_COMP.bits() | Self::LAST_SEGMENT.bits(); // This will be the illumina read 2
        const PAIRED_END_READ1_FLAG = Self::MULTI_SEGMENT.bits() | Self::REV_COMP_MATE.bits() | Self::FIRST_SEGMENT.bits();
        const PAIRED_END_READ2_FLAG = Self::MULTI_SEGMENT.bits() | Self::REV_COMP.bits() | Self::LAST_SEGMENT.bits();
    }
}

// TODO: We should ultimately get rid of this complex output type and use a struct instead.
// Output: (barcode, Vec<(umi, seq, qual, readname, flags)>, actual number of reads)
pub(crate) fn make_read_data(
    rna_reads: &[RnaRead],
    n_free_tail: usize,
    chemistry_def: &ChemistryDef,
) -> Vec<(String, DnaString, Vec<u8>, String, u16)> {
    assert!(!rna_reads.is_empty());
    let mut reads = Vec::new();

    for rna_read in rna_reads {
        // Determine if rna read X (not necessarily illumina read X) needs to
        // be reversed and complemented based on endedness.
        //
        // RNA read X: rx       Illumina read X: Rx
        //
        //          RNA read 1      RNA read 2      RC RNA read 1?      RC RNA read 2?
        // 5P-PE    R1              R2              N                   Y
        // 5P-R2    R2                              Y
        // 3P-PE    R1              R2              Y                   N
        // 3P-R2    R2                              N

        let (r1_revcomp, r2_revcomp) = {
            if rna_read.r2_exists() {
                if chemistry_def.endedness == Some(WhichEnd::FivePrime) {
                    (false, Some(true))
                } else {
                    (true, Some(false))
                }
            } else if chemistry_def.endedness == Some(WhichEnd::FivePrime) {
                (true, None)
            } else {
                (false, None)
            }
        };
        let (r1_seq, r1_qual) = {
            let r1_seq = rna_read.r1_seq();
            let r1_trimmed_len = get_seq_len_with_n_free_tail(r1_seq, n_free_tail);
            let mut r1_seq = r1_seq.to_vec();
            r1_seq.truncate(r1_trimmed_len);
            let mut r1_qual: Vec<_> = rna_read
                .r1_qual()
                .iter()
                .map(|q| q - ILLUMINA_QUAL_OFFSET)
                .collect();
            r1_qual.truncate(r1_trimmed_len);

            if r1_revcomp {
                (
                    byteseq::revcomp(r1_seq),
                    r1_qual.into_iter().rev().collect(),
                )
            } else {
                (r1_seq, r1_qual)
            }
        };

        if rna_read.r2_exists() {
            // Paired end

            // Read 1
            reads.push((
                rna_read.umi().to_string(),
                DnaString::from_acgt_bytes(&r1_seq),
                r1_qual,
                String::from_utf8(augmented_fastq_header(rna_read)).unwrap(),
                ReadFlags::PAIRED_END_READ1_FLAG.bits(),
            ));
            // Read 2
            let r2_seq = rna_read.r2_seq().unwrap();
            let r2_trimmed_len = get_seq_len_with_n_free_tail(r2_seq, n_free_tail);
            let mut r2_seq = r2_seq.to_vec();
            r2_seq.truncate(r2_trimmed_len);

            let mut r2_qual: Vec<_> = rna_read
                .r2_qual()
                .unwrap()
                .iter()
                .map(|q| q - ILLUMINA_QUAL_OFFSET)
                .collect();
            r2_qual.truncate(r2_trimmed_len);

            let (r2_seq, r2_qual) = {
                if r2_revcomp == Some(true) {
                    (
                        byteseq::revcomp(r2_seq),
                        r2_qual.into_iter().rev().collect(),
                    )
                } else {
                    (r2_seq, r2_qual)
                }
            };

            reads.push((
                rna_read.umi().to_string(),
                DnaString::from_acgt_bytes(&r2_seq),
                r2_qual,
                String::from_utf8(augmented_fastq_header(rna_read)).unwrap(),
                ReadFlags::PAIRED_END_READ2_FLAG.bits(),
            ));
        } else {
            // Single end
            reads.push((
                rna_read.umi().to_string(),
                DnaString::from_acgt_bytes(&r1_seq),
                r1_qual,
                String::from_utf8(augmented_fastq_header(rna_read)).unwrap(),
                ReadFlags::SINGLE_END_READ_FLAG.bits(),
            ));
        }
    }
    reads
}

#[derive(Default)]
pub struct UmiSortedReads {
    pub reads: Vec<DnaString>,
    pub quals: Vec<Vec<u8>>,
    pub readnames: Vec<String>,
    pub umi_id: Vec<i32>,         // Umi Id of each read
    pub unique_umis: Vec<String>, // List of sorted UMIs
    pub flags: Vec<u16>,
}

pub fn correct_umis(
    rna_reads: &mut [RnaRead],
    chemistry_def: &ChemistryDef,
) -> (i32, UmiSortedReads) {
    let mut ncorrected = 0;
    let umi_counts = SimpleHistogram::from_iter_owned(rna_reads.iter().map(RnaRead::umi));
    for read in &mut *rna_reads {
        let umi = read.umi();
        let self_count = umi_counts.get(&umi);
        if let Some(other) = umi
            .one_hamming_iter(HammingIterOpt::MutateNBase) // Among UMIs 1 HD away from this UMI
            .find(|u| umi_counts.get(u) >= 10 * self_count)
        {
            ncorrected += 1;
            read.correct_umi(other.seq());
        }
    }
    rna_reads.sort_by(|r1, r2| (r1.umi(), r1.header()).cmp(&(r2.umi(), r2.header())));
    let read_data = make_read_data(rna_reads, 8, chemistry_def);

    let mut umi_sorted_reads = UmiSortedReads::default();
    for (i, (umi, group_reads)) in read_data
        .into_iter()
        .group_by(|r| r.0.clone())
        .into_iter()
        .enumerate()
    {
        umi_sorted_reads.unique_umis.push(umi);
        for (_, seq, qual, readname, flags) in group_reads {
            umi_sorted_reads.reads.push(seq);
            umi_sorted_reads.quals.push(qual);
            umi_sorted_reads.readnames.push(readname);
            umi_sorted_reads.umi_id.push(i as i32);
            umi_sorted_reads.flags.push(flags);
        }
    }
    (ncorrected, umi_sorted_reads)
}

fn get_seq_len_with_n_free_tail(seq: &[u8], n_free_tail: usize) -> usize {
    let mut end = seq.len();
    while end > 0 {
        match (end.saturating_sub(n_free_tail)..end).find(|&j| seq[j] == b'N') {
            Some(j) => end = j,
            None => break,
        }
    }
    end
}

fn augmented_fastq_header(read: &RnaRead) -> Vec<u8> {
    const SAMPLE_INDEX_TAG: &[u8] = b"BC";
    const SAMPLE_INDEX_QUAL_TAG: &[u8] = b"QT";
    const RAW_BARCODE_TAG: &[u8] = b"CR";
    const RAW_BARCODE_QUAL_TAG: &[u8] = b"CY";
    const RAW_UMI_TAG: &[u8] = b"UR";
    const RAW_UMI_QUAL_TAG: &[u8] = b"UY";
    const PROCESSED_BARCODE_TAG: &[u8] = b"CB";
    const SEP: &[u8] = b"|||";

    let header = read.header().to_vec();
    let mut iter = header.split(|&m| m == b' ');
    let mut header_augment: Vec<u8> = Vec::new();
    header_augment.extend(iter.next().unwrap());

    header_augment.extend(SEP);
    header_augment.extend(SAMPLE_INDEX_TAG);
    header_augment.extend(SEP);
    // header_augment.extend([]);
    header_augment.extend(SEP);
    header_augment.extend(SAMPLE_INDEX_QUAL_TAG);
    header_augment.extend(SEP);
    // header_augment.extend([]);
    header_augment.extend(SEP);
    header_augment.extend(RAW_BARCODE_TAG);
    header_augment.extend(SEP);
    header_augment.extend(read.raw_bc_seq());
    header_augment.extend(SEP);
    header_augment.extend(RAW_BARCODE_QUAL_TAG);
    header_augment.extend(SEP);
    header_augment.extend(read.raw_bc_qual());
    header_augment.extend(SEP);
    header_augment.extend(RAW_UMI_TAG);
    header_augment.extend(SEP);
    header_augment.extend(read.raw_umi_seq());
    header_augment.extend(SEP);
    header_augment.extend(RAW_UMI_QUAL_TAG);
    header_augment.extend(SEP);
    header_augment.extend(read.raw_umi_qual());
    header_augment.extend(SEP);
    header_augment.extend(PROCESSED_BARCODE_TAG);
    header_augment.extend(SEP);
    header_augment.extend(read.barcode().to_string().as_bytes());
    header_augment.push(b' ');
    header_augment.extend(iter.collect::<Vec<_>>().join(&b' '));

    header_augment
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_n_free_tail() {
        let seq = b"AAAAAGACTTCACAAGAATGAGTAGTATGCCTAATGTAGT";
        assert_eq!(get_seq_len_with_n_free_tail(seq, 8), seq.len());

        let seq = b"AAAAAGACTTCACAAGAATGAGTAGTATGCCTAATGTAGN";
        assert_eq!(get_seq_len_with_n_free_tail(seq, 8), seq.len() - 1);

        let seq = b"AAAAAGACTTCACAAGAATGAGTAGTATGCCTNATGTAGT";
        assert_eq!(get_seq_len_with_n_free_tail(seq, 8), seq.len() - 8);

        let seq = b"AANAAGANTTCANAAGAATNAGTAGNATGCCTNATGTAGT";
        assert_eq!(get_seq_len_with_n_free_tail(seq, 8), 2);

        let seq = b"AAG";
        assert_eq!(get_seq_len_with_n_free_tail(seq, 8), 3);
    }
}
