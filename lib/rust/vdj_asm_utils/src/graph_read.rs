use crate::constants::{ReadType, UmiType};
use crate::{bam_utils, sw, utils};
use debruijn::dna_string::DnaString;
use rust_htslib::bam;
use serde::{Deserialize, Serialize};

/// A macro creating methods for flag access.
/// Copied from rust_htslib.
macro_rules! flag {
    ($get:ident, $set:ident, $unset:ident, $bit:expr) => {
        pub fn $get(&self) -> bool {
            self.flags & $bit != 0
        }

        pub fn $set(&mut self) {
            self.flags |= $bit;
        }

        pub fn $unset(&mut self) {
            self.flags &= u16::MAX - $bit;
        }
    };
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Read {
    pub name: String,
    pub seq: DnaString,
    pub quals: Vec<u8>,
    pub id: ReadType,
    pub umi: UmiType,
    pub flags: u16,
}

impl Read {
    pub fn new(id: ReadType, umi: UmiType, name: String, seq: DnaString, quals: Vec<u8>) -> Read {
        Read {
            name,
            id,
            umi,
            seq,
            quals,
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

        let v: Vec<u8> = seq_bytes
            .iter()
            .enumerate()
            .map(|(i, x)| utils::base_to_bits_hash(*x, &read_name, i))
            .collect();
        let seq_dna_string = DnaString::from_bytes(&v);

        let read = Read {
            name: read_name,
            id,
            umi,
            // NOTE: rust_htslib and DnaString use different representations of the bases!!
            // as bytes returns a vector, one element per base. The elements are u8
            // representations of 'A', 'C', 'G', 'T', or 'N'
            seq: seq_dna_string,
            quals: record.qual().to_vec(),
            flags: record.flags(),
        };
        read
    }

    pub fn len(&self) -> usize {
        self.seq.len()
    }
    pub fn is_empty(&self) -> bool {
        self.seq.is_empty()
    }

    pub fn to_bam_record(
        &self,
        alignment: &Option<sw::AlignmentPacket>,
        mate_alignment: &Option<sw::AlignmentPacket>,
    ) -> bam::record::Record {
        bam_utils::read_to_bam_record(self, alignment, mate_alignment)
    }

    pub fn to_bam_record_strip_augmented(
        &self,
        alignment: &Option<sw::AlignmentPacket>,
        mate_alignment: &Option<sw::AlignmentPacket>,
    ) -> bam::record::Record {
        bam_utils::read_to_bam_record_opts(self, alignment, mate_alignment, true, true)
    }

    pub fn to_unmapped_bam_record(&self) -> bam::record::Record {
        self.to_bam_record_strip_augmented(&None, &None)
    }

    /// True if the provided alignment has an implied length equal to this read's length.
    pub fn validate_alignment(&self, al_pack: &sw::AlignmentPacket) -> bool {
        use bio::alignment::AlignmentOperation::{Ins, Match, Subst};
        let alignment = &al_pack.alignment;
        let mut cigar_len = alignment.xstart;
        for &step in &alignment.operations {
            cigar_len += match step {
                Match | Subst | Ins => 1,
                _ => 0,
            };
        }

        cigar_len += alignment.xlen - alignment.xend;
        cigar_len == self.seq.len()
    }

    flag!(is_paired, set_paired, unset_paired, 1u16);
    flag!(is_proper_pair, set_proper_pair, unset_proper_pair, 2u16);
    flag!(is_unmapped, set_unmapped, unset_unmapped, 4u16);
    flag!(
        is_mate_unmapped,
        set_mate_unmapped,
        unset_mate_unmapped,
        8u16
    );
    flag!(is_reverse, set_reverse, unset_reverse, 16u16);
    flag!(is_mate_reverse, set_mate_reverse, unset_mate_reverse, 32u16);
    flag!(
        is_first_in_template,
        set_first_in_template,
        unset_first_in_template,
        64u16
    );
    flag!(
        is_last_in_template,
        set_last_in_template,
        unset_last_in_template,
        128u16
    );
    flag!(is_secondary, set_secondary, unset_secondary, 256u16);
    flag!(
        is_quality_check_failed,
        set_quality_check_failed,
        unset_quality_check_failed,
        512u16
    );
    flag!(is_duplicate, set_duplicate, unset_duplicate, 1024u16);
    flag!(
        is_supplementary,
        set_supplementary,
        unset_supplementary,
        2048u16
    );
}

#[cfg(test)]
mod tests {

    use crate::sw;
    use bio::alignment::{Alignment, AlignmentMode, AlignmentOperation};
    use rust_htslib::bam;
    use rust_htslib::bam::Read;
    use std::path::Path;

    /// Test the conversion between bam::record::Record and graph_read::Read as well
    /// as the generation of bam records from Read objects.
    #[test]
    fn test_read_to_bam_record() {
        let test_bam_name = "test/inputs/test_index.bam";

        let mut bam = bam::Reader::from_path(Path::new(&test_bam_name)).unwrap();
        let mut records = bam.records();

        for i in 0..4 {
            let rec1 = records.next().unwrap().expect("error reading test file");
            let rec2 = records.next().unwrap().expect("error reading test file");
            let read1 = super::Read::from_bam_record(i * 2_u32, 0, &rec1);
            let read2 = super::Read::from_bam_record(i * 2 + 1_u32, 0, &rec2);

            // read1 is 124, read2 iw 150 bases
            let al1 = Alignment {
                score: read1.len() as i32,
                xstart: 0,
                ystart: 0,
                xend: read1.len(),
                yend: read1.len(),
                ylen: 300,
                xlen: read1.len(),
                operations: vec![AlignmentOperation::Match; read1.len()],
                mode: AlignmentMode::Semiglobal,
            };
            let al_pack1 = Some(sw::AlignmentPacket::new(0, al1));
            let al2 = Alignment {
                score: read2.len() as i32,
                xstart: 0,
                ystart: 100,
                xend: read2.len(),
                yend: 100 + read2.len(),
                ylen: 300,
                xlen: read2.len(),
                operations: vec![AlignmentOperation::Match; read2.len()],
                mode: AlignmentMode::Semiglobal,
            };
            let al_pack2 = Some(sw::AlignmentPacket::new(0, al2));

            {
                // Test 1: both reads mapped
                let new_rec1 = read1.to_bam_record(&al_pack1, &al_pack2);
                let new_rec2 = read2.to_bam_record(&al_pack2, &al_pack1);
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

                assert_eq!(
                    new_rec1.cigar().iter().copied().collect::<Vec<_>>(),
                    vec![bam::record::Cigar::Equal(read1.len() as u32)]
                );

                assert_eq!(
                    new_rec2.cigar().iter().copied().collect::<Vec<_>>(),
                    vec![bam::record::Cigar::Equal(read2.len() as u32)]
                );
            }

            {
                // Test 2: First read mapped, second unmapped. Make sure tids are set correctly.
                let new_rec1 = read1.to_bam_record(&al_pack1, &None);
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
                let new_rec1 = read1.to_bam_record(&None, &al_pack2);
                assert!(new_rec1.is_unmapped());
                assert!(!new_rec1.is_mate_unmapped());
                // unmapped mate should have the coordinates of the mapped one
                assert_eq!(new_rec1.tid(), new_rec1.mtid());
                assert_eq!(new_rec1.pos(), new_rec1.mpos());

                assert_eq!(new_rec1.qname(), rec1.qname());
                assert_eq!(new_rec1.qual(), rec1.qual());

                let new_rec2 = read2.to_bam_record(&None, &al_pack2);
                assert_eq!(new_rec2.seq().as_bytes(), rec2.seq().as_bytes());
            }

            {
                // Test 4: Both reads unmapped. Make sure tids and flags are set correctly.
                // Make sure insert size is set to 0.
                let new_rec1 = read1.to_bam_record(&None, &None);
                assert!(new_rec1.is_unmapped());
                assert!(new_rec1.is_mate_unmapped());
                assert_eq!(new_rec1.tid(), -1);
                assert_eq!(new_rec1.mtid(), -1);

                let new_rec2 = read2.to_bam_record(&None, &None);
                assert_eq!(new_rec2.seq().as_bytes(), rec2.seq().as_bytes());

                assert_eq!(new_rec1.insert_size(), 0);
                assert_eq!(new_rec2.insert_size(), 0);
            }
            {
                // Test 5: Make sure that insert size is set correctly.

                // Read1 cigar is 3S 100M 1X 1I 1D 19S
                // It extends from base 0 (inclusive) to base 102 (exclusive).
                let mut op = vec![AlignmentOperation::Xclip(3)];
                for _ in 0..100 {
                    op.push(AlignmentOperation::Match);
                }
                op.push(AlignmentOperation::Subst);
                op.push(AlignmentOperation::Ins);
                op.push(AlignmentOperation::Del);
                op.push(AlignmentOperation::Xclip(19));
                op.push(AlignmentOperation::Yclip(198));

                let al1 = Alignment {
                    score: 95,
                    xstart: 3,
                    ystart: 0,
                    xend: 105,
                    yend: 102,
                    ylen: 300,
                    xlen: read1.len(),
                    operations: op,
                    mode: AlignmentMode::Custom,
                };
                let al_pack1 = Some(sw::AlignmentPacket::new(0, al1));

                // Read2 cigar is 3S 100M 47S
                // It extends from base 100 to base 200 (exclusive) (so after read1)
                let mut op = vec![AlignmentOperation::Xclip(3), AlignmentOperation::Yclip(100)];
                for _ in 0..100 {
                    op.push(AlignmentOperation::Match);
                }
                op.push(AlignmentOperation::Xclip(47));
                op.push(AlignmentOperation::Yclip(100));

                let al2 = Alignment {
                    score: 98,
                    xstart: 3,
                    ystart: 100,
                    xend: 103,
                    yend: 200,
                    ylen: 300,
                    xlen: read2.len(),
                    operations: op,
                    mode: AlignmentMode::Custom,
                };
                let al_pack2 = Some(sw::AlignmentPacket::new(0, al2));

                let new_rec1 = read1.to_bam_record(&al_pack1, &al_pack2);
                let new_rec2 = read2.to_bam_record(&al_pack2, &al_pack1);
                assert_eq!(new_rec1.insert_size(), 200);
                assert_eq!(new_rec2.insert_size(), -200);

                // Read2 cigar is 3S 100M 47S
                // It extends from base 10 to base 110 (exclusive) (so partially overlapping read1)
                let mut op = vec![AlignmentOperation::Xclip(3), AlignmentOperation::Yclip(10)];
                for _ in 0..100 {
                    op.push(AlignmentOperation::Match);
                }
                op.push(AlignmentOperation::Xclip(47));
                op.push(AlignmentOperation::Yclip(190));
                let al2 = Alignment {
                    score: 98,
                    xstart: 3,
                    ystart: 10,
                    xend: 103,
                    yend: 110,
                    ylen: 300,
                    xlen: read2.len(),
                    operations: op,
                    mode: AlignmentMode::Custom,
                };
                let al_pack2 = Some(sw::AlignmentPacket::new(0, al2));

                let new_rec1 = read1.to_bam_record(&al_pack1, &al_pack2);
                let new_rec2 = read2.to_bam_record(&al_pack2, &al_pack1);
                assert_eq!(new_rec1.insert_size(), 110);
                assert_eq!(new_rec2.insert_size(), -110);

                // Read2 cigar should be 3S 10M 147S
                // It extends from base 0 to base 10 (exclusive) (so completely contained within read1)
                let mut op = vec![AlignmentOperation::Xclip(3)];
                for _ in 0..10 {
                    op.push(AlignmentOperation::Match);
                }
                op.push(AlignmentOperation::Xclip(147));
                op.push(AlignmentOperation::Yclip(290));
                let al2 = Alignment {
                    score: 98,
                    xstart: 3,
                    ystart: 0,
                    xend: 13,
                    yend: 10,
                    ylen: 300,
                    xlen: read2.len(),
                    operations: op,
                    mode: AlignmentMode::Custom,
                };
                let al_pack2 = Some(sw::AlignmentPacket::new(0, al2));

                let new_rec1 = read1.to_bam_record(&al_pack1, &al_pack2);
                let new_rec2 = read2.to_bam_record(&al_pack2, &al_pack1);
                assert_eq!(new_rec1.insert_size(), 102);
                assert_eq!(new_rec2.insert_size(), 102); // Is that the right behavior??
            }
        }
    }
}
