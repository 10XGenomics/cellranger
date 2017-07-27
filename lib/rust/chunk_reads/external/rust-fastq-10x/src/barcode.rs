// Copyright (c) 2016 10x Genomics, Inc. All rights reserved.

use std::collections::HashMap;
use std::cmp::max;
use ordered_float::NotNaN;

pub struct BarcodeValidator<'a> {
	pub whitelist: &'a HashMap<Vec<u8>, u32>,
    pub bc_counts: &'a HashMap<(u8, u32), u32>,
	pub max_expected_barcode_errors: f64,
	pub bc_confidence_threshold: f64,
}

const BASE_OPTS: [u8; 4] = [b'A', b'C', b'G', b'T'];

type of64 = NotNaN<f64>;

impl<'a> BarcodeValidator<'a> {

    fn do_correct_barcode(&self, gem_group: u8, barcode: &[u8], qual: &[u8]) -> Option<u32> {
        let mut a = Vec::from(barcode);

        let mut candidates: Vec<(of64, Vec<u8>)> = Vec::new();
        let mut total_likelihood = of64::from(0.0);
         
	    for pos in 0 .. barcode.len() {
            let qv = qual[pos];
            let existing = a[pos];
            for val in BASE_OPTS.iter().cloned() {

                if val == existing { continue; }
                a[pos] = val;

                match self.whitelist.get(&a) {
                    Some(bc_id) => {
                        let bc_count = self.bc_counts.get(&(gem_group, *bc_id)).cloned().unwrap_or(0);
                        let prob_edit = max(of64::from(0.0005), of64::from(probability(qv)));
                        let likelihood = prob_edit * max(of64::from(bc_count as f64), of64::from(0.5));
                        candidates.push((likelihood, a.clone()));
                        total_likelihood += likelihood;
                    },
                    None => (),
                }
            }
            a[pos] = existing;
        }
	
        println!("{:?}", candidates);

        let thresh = of64::from(self.bc_confidence_threshold);

        let best_option = candidates.into_iter().max();
        
        match best_option {
            Some((best_like, bc)) => {
                if best_like / total_likelihood > thresh {
                    self.whitelist.get(&bc).cloned()
                } else {
                    None
                }
            },
            _ => None
        }
    }

    pub fn correct_barcode(&self, gem_group: u8, barcode: &[u8], qual: &[u8]) -> Option<u32> {
        let expected_errors: f64 = qual.iter().cloned().map(probability).sum();

        let bc_id = 
            match self.whitelist.get(barcode) {
                Some(id) => Some(*id),
                None => self.do_correct_barcode(gem_group, barcode, qual),
            };

        if bc_id.is_some() && expected_errors < self.max_expected_barcode_errors {
            bc_id
        } else {
            None
        }
    }
}


pub fn probability(qual: u8) -> f64 {
    //33 is the illumina qual offset
    let q = qual as f64;
    (10_f64).powf(-(q - 33.0) / 10.0) 
}

#[cfg(test)]
mod test {
    use super::*;
    use std::collections::HashMap;

    #[test]
    pub fn test_bc_correct()
    {
        let mut wl = HashMap::new();
        wl.insert(Vec::from(b"AAAAA".as_ref()), 1);
        wl.insert(Vec::from(b"AAGAC".as_ref()), 2);
        wl.insert(Vec::from(b"ACGAA".as_ref()), 3);
        wl.insert(Vec::from(b"ACGTT".as_ref()), 4);

        let mut counts = HashMap::new();
        counts.insert((0,1), 100);
        counts.insert((0,2), 11);
        counts.insert((0,3), 2);

        let val = BarcodeValidator {
            max_expected_barcode_errors: 1.0,
            bc_confidence_threshold: 0.95,
            whitelist: &wl,
            bc_counts: &counts,
        };

        // Easy
        assert_eq!(val.correct_barcode(0, b"AAAAA", &vec![66,66,66,66,66]), Some(1));

        // Low quality
        assert_eq!(val.correct_barcode(0, b"AAAAA", &vec![34,34,34,66,66]), None);

        // Trivial correction
        assert_eq!(val.correct_barcode(0, b"AAAAT", &vec![66,66,66,66,40]), Some(1));

        // Pseudo-count kills you
        assert_eq!(val.correct_barcode(0, b"ACGAT", &vec![66,66,66,66,66]), None);

        // Quality help you
        assert_eq!(val.correct_barcode(0, b"ACGAT", &vec![66,66,66,66,40]), Some(3));

        // Counts help you
        assert_eq!(val.correct_barcode(0, b"ACAAA", &vec![66,66,66,66,40]), Some(1));
    }
}