//
// Copyright (c) 2017 10x Genomics, Inc. All rights reserved.
//

extern crate tada;
use tada::bitenc;

use fastq::CellrangerFastqHeader;
use sw;
use bam_utils;
use constants::{UmiType, ReadType};
use std;
use utils;

use rust_htslib::bam;

/// A macro creating methods for flag access.
/// Copied from rust_htslib.
macro_rules! flag {
    ($get:ident, $set:ident, $unset:ident, $bit:expr) => (
        pub fn $get(&self) -> bool {
            self.flags & $bit != 0
        }

        pub fn $set(&mut self) {
            self.flags |= $bit;
        }

        pub fn $unset(&mut self) {
            self.flags &= std::u16::MAX - $bit;
        }
    )
}

#[derive(Clone, Debug)]
pub struct Read {
    pub name: String,
    pub seq: bitenc::BitEnc,
    pub quals: Vec<u8>,
    pub id: ReadType,
    pub umi: UmiType,
    pub flags: u16,
}

impl Read {
    pub fn new(id: ReadType, umi: UmiType, name: String, seq: bitenc::BitEnc, quals: Vec<u8>) -> Read {
        Read {
            name: name,
            id: id,
            umi: umi,
            seq: seq,
            quals: quals,
            flags: 0,
        }
    }

    /// Create Read object from a BAM Record object.
    ///
    /// NOTE: If the read sequence has any characters other than A, C, G, T (uppercase)
    /// then the flag "Quality check failed" will be set. This is done because
    /// we rely on a 2-bit encoding of bases and other characters cannot be handled properly.
    pub fn from_bam_record(id: ReadType, umi: UmiType, record: &bam::record::Record) -> Read {
        let read_name = String::from_utf8_lossy(record.qname()).into_owned();
        let seq_bytes = record.seq().as_bytes();

        let read = Read {
            name: read_name.clone(),
            id: id,
            umi: umi,
            // NOTE: rust_htslib and bitencs use different representations of the bases!!
            // as bytes returns a vector, one element per base. The elements are u8
            // representations of 'A', 'C', 'G', 'T', or 'N'
            seq: bitenc::BitEnc::from_bytes(&(seq_bytes.iter().enumerate().map(|(i, x)| utils::base_to_bits_hash(*x, &read_name, i)).collect())),
            quals: record.qual().to_vec(),
            flags: record.flags(),
        };
        read
    }

    pub fn get_tag(&self, tag: &String) -> Option<String> {
        CellrangerFastqHeader::new(self.name.clone()).get_tag(tag)
    }

    pub fn barcode(&self) -> Option<String> {
        CellrangerFastqHeader::new(self.name.clone()).get_barcode()
    }

    pub fn len(&self) -> usize {
        self.seq.len()
    }

    pub fn to_bam_record(&self, alignment: Option<sw::Alignment>,
                         mate_alignment: Option<sw::Alignment>) -> bam::record::Record {
        bam_utils::read_to_bam_record(&self, alignment, mate_alignment)
    }

    pub fn to_bam_record_strip_augmented(&self, alignment: Option<sw::Alignment>,
                         mate_alignment: Option<sw::Alignment>) -> bam::record::Record {
        bam_utils::read_to_bam_record_opts(&self, alignment, mate_alignment, true, true)
                         }

    pub fn to_unmapped_bam_record(&self) -> bam::record::Record {
        self.to_bam_record(None, None)
    }

    /// True if the provided alignment has an implied length equal to this read's length.
    pub fn validate_alignment(&self, alignment: &sw::Alignment) -> bool {
        let mut cigar_len = 0;
        for step in alignment.alignment.iter() {
            cigar_len += match *step {
                sw::AlignmentStep::Del => 0,
                _ => 1,
            };
        }
        cigar_len == self.seq.len()
    }

    flag!(is_paired, set_paired, unset_paired, 1u16);
    flag!(is_proper_pair, set_proper_pair, unset_proper_pair, 2u16);
    flag!(is_unmapped, set_unmapped, unset_unmapped, 4u16);
    flag!(is_mate_unmapped, set_mate_unmapped, unset_mate_unmapped, 8u16);
    flag!(is_reverse, set_reverse, unset_reverse, 16u16);
    flag!(is_mate_reverse, set_mate_reverse, unset_mate_reverse, 32u16);
    flag!(is_first_in_template, set_first_in_template, unset_first_in_template, 64u16);
    flag!(is_last_in_template, set_last_in_template, unset_last_in_template, 128u16);
    flag!(is_secondary, set_secondary, unset_secondary, 256u16);
    flag!(is_quality_check_failed, set_quality_check_failed, unset_quality_check_failed, 512u16);
    flag!(is_duplicate, set_duplicate, unset_duplicate, 1024u16);
    flag!(is_supplementary, set_supplementary, unset_supplementary, 2048u16);
}

#[cfg(test)]
mod tests {
    use rust_htslib::bam;
    use rust_htslib::bam::Read;
    use std::path::Path;
    use sw;
    use constants::ReadType;

    /// Test the conversion between bam::record::Record and graph_read::Read as well
    /// as the generation of bam records from Read objects.
    #[test]
    fn test_read_to_bam_record() {

        let test_bam_name = "test/inputs/test_index.bam";

        let bam = bam::Reader::from_path(&Path::new(&(test_bam_name.to_string()))).unwrap();
        let mut records = bam.records();

        for i in 0..4 {
            let rec1 = records.next().unwrap().ok().expect("error reading test file");
            let rec2 = records.next().unwrap().ok().expect("error reading test file");
            let read1 = super::Read::from_bam_record(i * 2 as ReadType, 0, &rec1);
            let read2 = super::Read::from_bam_record(i * 2 + 1 as ReadType, 0, &rec2);

            // read1 is 124, read2 iw 150 bases
            let al1 = sw::Alignment::new(10.0, 0, 0, vec![sw::AlignmentStep::Match; read1.len()]);
            let al2 = sw::Alignment::new(10.0, 0, 100, vec![sw::AlignmentStep::Match; read2.len()]);

            {
                // Test 1: both reads mapped
                let new_rec1 = read1.to_bam_record(Some(al1.clone()), Some(al2.clone()));
                let new_rec2 = read2.to_bam_record(Some(al2.clone()), Some(al1.clone()));
                assert_eq!(new_rec1.insert_size(), 250);
                assert_eq!(new_rec2.insert_size(), -250);
                assert!(!new_rec1.is_unmapped());
                assert!(!new_rec2.is_unmapped());
                assert!(!new_rec1.is_mate_unmapped());
                // We should have kept these as they were.
                assert!(new_rec1.is_first_in_template());
                assert!(new_rec2.is_last_in_template());
                assert_eq!(new_rec1.is_reverse(), rec1.is_reverse());
                assert_eq!(new_rec2.is_reverse(), rec2.is_reverse());
                assert_eq!(new_rec1.is_mate_reverse(), rec1.is_mate_reverse());
                assert_eq!(new_rec1.tid(), new_rec2.mtid());
                assert_eq!(new_rec1.pos(), new_rec2.mpos());
                assert_eq!(new_rec1.tid(), new_rec1.tid());

                assert_eq!(new_rec1.qname(), rec1.qname());
                assert_eq!(new_rec1.qual(), rec1.qual());
                assert_eq!(new_rec1.seq().as_bytes(), rec1.seq().as_bytes());

                assert_eq!(new_rec2.qname(), rec2.qname());
                assert_eq!(new_rec2.qual(), rec2.qual());
                assert_eq!(new_rec2.seq().as_bytes(), rec2.seq().as_bytes());

                assert_eq!(new_rec1.cigar(), vec![bam::record::Cigar::Equal(read1.len() as u32)]);
                assert_eq!(new_rec2.cigar(), vec![bam::record::Cigar::Equal(read2.len() as u32)]);
            }

            {
                // Test 2: First read mapped, second unmapped. Make sure tids are set correctly.
                let new_rec1 = read1.to_bam_record(Some(al1.clone()), None);
                assert!(!new_rec1.is_unmapped());
                assert!(new_rec1.is_mate_unmapped());
                assert_eq!(new_rec1.insert_size(), 0);
                // unmapped mate should have the coordinates of the mapped one
                assert_eq!(new_rec1.tid(), new_rec1.mtid());
                assert_eq!(new_rec1.pos(), new_rec1.mpos());

                assert_eq!(new_rec1.qname(), rec1.qname());
                assert_eq!(new_rec1.qual(), rec1.qual());
                assert_eq!(new_rec1.seq().as_bytes(), rec1.seq().as_bytes());
            }

            {
                // Test 3: First read unmapped, second mapped. Make sure tids are set correctly.
                let new_rec1 = read1.to_bam_record(None, Some(al2.clone()));
                assert!(new_rec1.is_unmapped());
                assert!(!new_rec1.is_mate_unmapped());
                // unmapped mate should have the coordinates of the mapped one
                assert_eq!(new_rec1.tid(), new_rec1.mtid());
                assert_eq!(new_rec1.pos(), new_rec1.mpos());

                assert_eq!(new_rec1.qname(), rec1.qname());
                assert_eq!(new_rec1.qual(), rec1.qual());

                let new_rec2 = read2.to_bam_record(None, Some(al2.clone()));
                assert_eq!(new_rec2.seq().as_bytes(), rec2.seq().as_bytes());
            }

            {
                // Test 4: Both reads unmapped. Make sure tids and flags are set correctly.
                // Make sure insert size is set to 0.
                let new_rec1 = read1.to_bam_record(None, None);
                assert!(new_rec1.is_unmapped());
                assert!(new_rec1.is_mate_unmapped());
                assert_eq!(new_rec1.tid(), -1);
                assert_eq!(new_rec1.mtid(), -1);

                let new_rec2 = read2.to_bam_record(None, None);
                assert_eq!(new_rec2.seq().as_bytes(), rec2.seq().as_bytes());

                assert_eq!(new_rec1.insert_size(), 0);
                assert_eq!(new_rec2.insert_size(), 0);
            }
            {
                // Test 5: Make sure that insert size is set correctly.

                // Read1 cigar should be 3S 101M 1I 1D 19S
                // It extends from base 0 (inclusive) to base 102 (exclusive).
                let mut cigar1 = vec![sw::AlignmentStep::Clip, sw::AlignmentStep::Clip, sw::AlignmentStep::Clip];
                for _ in 0..100 { cigar1.push(sw::AlignmentStep::Match); }
                cigar1.push(sw::AlignmentStep::Mismatch);
                cigar1.push(sw::AlignmentStep::Ins);
                cigar1.push(sw::AlignmentStep::Del);
                for _ in 0..(read1.len() - cigar1.len() + 1) {
                    cigar1.push(sw::AlignmentStep::Clip);
                }
                let al1 = sw::Alignment::new(10.0, 0, 0, cigar1);

                // Read2 cigar should be 3S 100M 47S
                // It extends from base 100 to base 200 (exclusive) (so after read1)
                let mut cigar2 = vec![sw::AlignmentStep::Clip, sw::AlignmentStep::Clip, sw::AlignmentStep::Clip];
                for _ in 0..100 { cigar2.push(sw::AlignmentStep::Match); }
                for _ in 0..(read2.len() - cigar2.len()) {
                    cigar2.push(sw::AlignmentStep::Clip);
                }
                let al2 = sw::Alignment::new(10.0, 0, 100, cigar2);

                let new_rec1 = read1.to_bam_record(Some(al1.clone()), Some(al2.clone()));
                let new_rec2 = read2.to_bam_record(Some(al2.clone()), Some(al1.clone()));
                assert_eq!(new_rec1.insert_size(), 200);
                assert_eq!(new_rec2.insert_size(), -200);

                // Read2 cigar should be 3S 100M 47S
                // It extends from base 10 to base 110 (exclusive) (so partially overlapping read1)
                cigar2 = vec![sw::AlignmentStep::Clip, sw::AlignmentStep::Clip, sw::AlignmentStep::Clip];
                for _ in 0..100 { cigar2.push(sw::AlignmentStep::Match); }
                for _ in 0..(read2.len() - cigar2.len()) {
                    cigar2.push(sw::AlignmentStep::Clip);
                }
                let al2 = sw::Alignment::new(10.0, 0, 10, cigar2);

                let new_rec1 = read1.to_bam_record(Some(al1.clone()), Some(al2.clone()));
                let new_rec2 = read2.to_bam_record(Some(al2.clone()), Some(al1.clone()));
                assert_eq!(new_rec1.insert_size(), 110);
                assert_eq!(new_rec2.insert_size(), -110);

                // Read2 cigar should be 3S 10M 147S
                // It extends from base 0 to base 10 (exclusive) (so completely contained within read1)
                cigar2 = vec![sw::AlignmentStep::Clip, sw::AlignmentStep::Clip, sw::AlignmentStep::Clip];
                for _ in 0..10 { cigar2.push(sw::AlignmentStep::Match); }
                for _ in 0..(read2.len() - cigar2.len()) {
                    cigar2.push(sw::AlignmentStep::Clip);
                }
                let al2 = sw::Alignment::new(10.0, 0, 0, cigar2);

                let new_rec1 = read1.to_bam_record(Some(al1.clone()), Some(al2.clone()));
                let new_rec2 = read2.to_bam_record(Some(al2.clone()), Some(al1.clone()));
                assert_eq!(new_rec1.insert_size(), 102);
                assert_eq!(new_rec2.insert_size(), 102); // Is that the right behavior??
            }
        }
    }
}
