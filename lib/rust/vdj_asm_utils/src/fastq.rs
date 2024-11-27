use std::collections::HashMap;

pub struct CellrangerFastqHeader {
    pub qname: String,
    pub header: String,
    pub tags: HashMap<String, String>,
}

impl CellrangerFastqHeader {
    pub fn new(header: String) -> CellrangerFastqHeader {
        let parts: Vec<_> = header.split("|||").map(String::from).collect();
        let mut tags = HashMap::new();
        if parts.len() > 1 {
            for idx in 1..parts.len() {
                if idx % 2 == 0 {
                    continue;
                }
                tags.insert(parts[idx].to_string(), parts[idx + 1].to_string());
            }
        }

        CellrangerFastqHeader {
            qname: parts[0].to_string(),
            header,
            tags,
        }
    }

    #[cfg(test)]
    fn get_tag(&self, tag: &str) -> Option<&str> {
        self.tags.get(tag).map(String::as_str)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::{ReadType, UmiType, PROCESSED_UMI_TAG, QUAL_OFFSET, RAW_UMI_TAG};
    use crate::utils;
    use io_utils::open_maybe_compressed;
    use std::io::{BufRead, BufReader, Lines, Read};
    use std::path::Path;

    #[derive(Clone, Debug)]
    pub struct Record {
        pub id: ReadType,
        pub seq: Vec<u8>,
        #[allow(dead_code)]
        pub quals: Vec<u8>,
    }

    fn extract_head(header: &str) -> String {
        // Remove the @ at the beginning and then take stuff before whitespace
        let no_at = (header.split_at(1).1).split_whitespace().next().unwrap();
        no_at.to_string()
    }
    pub fn get_header_tag(header: &str, tag: &str) -> Option<String> {
        let parts = header.split("|||").collect::<Vec<&str>>();
        for (idx, part) in parts.iter().enumerate() {
            if *part == tag {
                return Some(parts[idx + 1].to_string());
            }
        }
        Some(String::default())
    }
    struct PairedInputRead {
        pub r1: Record,
        pub r2: Option<Record>,
        pub umi: UmiType,
    }
    struct CellrangerPairedFastqIter {
        lines1: Lines<BufReader<Box<dyn Read>>>,
        lines2: Option<Lines<BufReader<Box<dyn Read>>>>,
        rev_strand: bool,
        read_count: usize,
        umis: HashMap<String, UmiType>,
    }
    impl CellrangerPairedFastqIter {
        pub fn new(
            name1: &Path,
            name2: Option<&Path>,
            rev_strand: bool,
        ) -> CellrangerPairedFastqIter {
            let r1 = BufReader::new(open_maybe_compressed(name1));
            let r2 = name2.map(|x| BufReader::new(open_maybe_compressed(x)));
            CellrangerPairedFastqIter {
                lines1: r1.lines(),
                lines2: r2.map(BufRead::lines),
                rev_strand,
                read_count: 0,
                umis: HashMap::new(),
            }
        }
    }
    impl Iterator for CellrangerPairedFastqIter {
        type Item = PairedInputRead;
        fn next(&mut self) -> Option<PairedInputRead> {
            let lines1 = &mut self.lines1;
            match lines1.next() {
                Some(head) => {
                    let read_seq = lines1.next().unwrap().unwrap();
                    let _ = lines1.next().unwrap().unwrap();
                    let q1 = lines1.next().unwrap().unwrap();
                    let header = head.unwrap();
                    let read_name = extract_head(&header);
                    self.read_count += 1;
                    let mut umi = get_header_tag(&header, PROCESSED_UMI_TAG).unwrap();
                    if umi.is_empty() {
                        umi = get_header_tag(&header, RAW_UMI_TAG).unwrap();
                    }
                    let umi_id;
                    if self.umis.contains_key(&umi) {
                        umi_id = self.umis[&umi];
                    } else {
                        umi_id = self.umis.len() as u32;
                        self.umis.insert(umi, umi_id);
                    }
                    let fw_seq: Vec<u8> = match self.rev_strand {
                        true => rc(&read_seq)
                            .into_bytes()
                            .into_iter()
                            .enumerate()
                            .map(|(i, x)| utils::base_to_bits_hash(x, &read_name, i))
                            .collect(),
                        false => read_seq
                            .into_bytes()
                            .into_iter()
                            .enumerate()
                            .map(|(i, x)| utils::base_to_bits_hash(x, &read_name, i))
                            .collect(),
                    };
                    let mut qual: Vec<u8> =
                        q1.into_bytes().iter().map(|x| x - QUAL_OFFSET).collect();
                    if self.rev_strand {
                        qual.reverse();
                    }
                    let (read_id, rec2) = match self.lines2 {
                        Some(ref mut lines2) => {
                            let _ = lines2.next().unwrap().unwrap();
                            let read_seq2 = lines2.next().unwrap().unwrap();
                            let _ = lines2.next().unwrap().unwrap();
                            let q2 = lines2.next().unwrap().unwrap();
                            let fw_seq2: Vec<u8> = match self.rev_strand {
                                true => read_seq2
                                    .into_bytes()
                                    .into_iter()
                                    .enumerate()
                                    .map(|(i, x)| utils::base_to_bits_hash(x, &read_name, i))
                                    .collect(),
                                false => rc(&read_seq2)
                                    .into_bytes()
                                    .into_iter()
                                    .enumerate()
                                    .map(|(i, x)| utils::base_to_bits_hash(x, &read_name, i))
                                    .collect(),
                            };
                            let mut qual2: Vec<u8> =
                                q2.into_bytes().iter().map(|x| x - QUAL_OFFSET).collect();
                            if !self.rev_strand {
                                qual2.reverse();
                            }
                            self.read_count += 1;
                            (
                                (self.read_count as ReadType - 2),
                                Some(Record {
                                    id: self.read_count as ReadType - 1,
                                    seq: fw_seq2,
                                    quals: qual2,
                                }),
                            )
                        }
                        None => (2 * (self.read_count as ReadType) - 2, None),
                    };
                    Some(PairedInputRead {
                        r1: Record {
                            id: read_id as ReadType,
                            seq: fw_seq,
                            quals: qual,
                        },
                        r2: rec2,
                        umi: umi_id,
                    })
                }
                None => None,
            }
        }
    }

    fn rc(seq: &str) -> String {
        // result vector
        let mut rc_seq: String = String::with_capacity(seq.len());
        // iterate through the input &str
        for c in seq.chars().rev() {
            if c == 'A' {
                rc_seq.push('T');
            } else if c == 'C' {
                rc_seq.push('G');
            } else if c == 'G' {
                rc_seq.push('C');
            } else if c == 'T' {
                rc_seq.push('A');
            } else {
                rc_seq.push('X');
            }
        }
        rc_seq
    }
    #[test]
    fn test_fastq_header() {
        let name = "D000684:775:H3YYTBCXY:1:1101:10541:57223|||BC|||TGGTAAAC|||UR|||AGAGCTGCCA|||UY|||IIIIIIIIII|||TR||||||CB|||TCAGATGCAGGCTCAC-1|||UB|||AGAGCTGCCA";
        let header = CellrangerFastqHeader::new(name.to_string());
        assert_eq!(header.get_tag("BC"), Some("TGGTAAAC"));
        assert_eq!(header.get_tag("TR"), Some(""));
        assert_eq!(header.get_tag("FOO"), None);
    }

    #[test]
    fn test_fastq_iter() {
        let fastq1 = Path::new("test/inputs/base_quals/test_base_quals_1.fastq");
        let fastq2 = Path::new("test/inputs/base_quals/test_base_quals_2.fastq");

        let true_umis = [0, 0, 1, 1];

        {
            // Paired, FR reads
            let mut fq_iter = CellrangerPairedFastqIter::new(fastq1, Some(fastq2), false);

            // Four pairs in the input files.
            for i in 0..4 {
                let pair = fq_iter.next().unwrap();
                let r1 = pair.r1;
                let r2 = pair.r2.unwrap();
                assert_eq!(r1.id, 2 * i);
                assert_eq!(r2.id, 2 * i + 1);
                assert_eq!(pair.umi, true_umis[i as usize]);
                if i == 0 {
                    assert_eq!(r1.seq[0..5].to_vec(), vec![1, 3, 0, 0, 3]);
                    assert_eq!(r2.seq[0..5].to_vec(), vec![1, 1, 3, 2, 0]);
                }
            }
            assert_eq!(fq_iter.read_count, 8);
        }
        {
            // Paired, RF reads
            let mut fq_iter = CellrangerPairedFastqIter::new(fastq1, Some(fastq2), true);

            for i in 0..4 {
                let pair = fq_iter.next().unwrap();
                let r1 = pair.r1;
                let r2 = pair.r2.unwrap();
                assert_eq!(r1.id, 2 * i);
                assert_eq!(r2.id, 2 * i + 1);
                assert_eq!(pair.umi, true_umis[i as usize]);
                if i == 0 {
                    assert_eq!(r1.seq[0..5].to_vec(), vec![0, 3, 2, 0, 1]);
                    assert_eq!(r2.seq[0..5].to_vec(), vec![1, 2, 2, 0, 1]);
                }
            }
            assert_eq!(fq_iter.read_count, 8);
        }
        {
            // Single-end, F reads
            let mut fq_iter = CellrangerPairedFastqIter::new(fastq1, None, false);

            for i in 0..4 {
                let pair = fq_iter.next().unwrap();
                let r1 = pair.r1;
                assert!(pair.r2.is_none());
                assert_eq!(r1.id, 2 * i);
                assert_eq!(pair.umi, true_umis[i as usize]);
                if i == 0 {
                    assert_eq!(r1.seq[0..5].to_vec(), vec![1, 3, 0, 0, 3]);
                } else if i == 1 {
                    assert_eq!(r1.seq[0..5].to_vec(), vec![0, 2, 0, 2, 1]);
                }
            }
            assert_eq!(fq_iter.read_count, 4);
        }
        {
            // Single-end, R reads
            let mut fq_iter = CellrangerPairedFastqIter::new(fastq1, None, true);

            for i in 0..4 {
                let pair = fq_iter.next().unwrap();
                let r1 = pair.r1;
                assert!(pair.r2.is_none());
                assert_eq!(r1.id, 2 * i);
                assert_eq!(pair.umi, true_umis[i as usize]);
                if i == 0 {
                    assert_eq!(r1.seq[0..5].to_vec(), vec![0, 3, 2, 0, 1]); // RC of GTCAT
                } else if i == 1 {
                    assert_eq!(r1.seq[0..5].to_vec(), vec![0, 3, 3, 2, 0]); // RC of TCAAT
                }
            }
            assert_eq!(fq_iter.read_count, 4);
        }
    }
}
