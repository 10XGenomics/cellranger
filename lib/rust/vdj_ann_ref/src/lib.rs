//! vdj_ann_ref
// Copyright (c) 2018 10X Genomics, Inc. All rights reserved.
#![deny(missing_docs)]

// This file contains code to make reference data.
//
// ◼ Experiment with building a reference from scratch.  This would be a better
// ◼ solution than ad hoc editing of a flawed reference.
//
// ◼ Document reference sequence requirements so that a customer who wishes to
// ◼ create a reference for a new species will know the conventions used by the
// ◼ code.

/// Return the path of the human VDJ reference
pub fn human_ref() -> &'static str {
    include_str!["../vdj_refs/human/fasta/regions.fa"]
}

/// Open a file for reading
#[macro_export]
macro_rules! open_for_read {
    ($filename:expr) => {
        ::std::io::BufReader::new(
            ::std::fs::File::open(::core::convert::AsRef::<::std::path::Path>::as_ref(
                $filename,
            ))
            .unwrap_or_else(|_| {
                panic!(
                    "Could not open file \"{}\"",
                    ::core::convert::AsRef::<::std::path::Path>::as_ref($filename)
                        .to_string_lossy(),
                )
            }),
        )
    };
}

/// Open a file for writing
#[macro_export]
macro_rules! open_for_write_new {
    ($filename:expr) => {
        ::std::io::BufWriter::new(
            ::std::fs::File::create(::core::convert::AsRef::<::std::path::Path>::as_ref(
                $filename,
            ))
            .unwrap_or_else(|_| {
                panic!(
                    "Could not create file \"{}\"",
                    ::core::convert::AsRef::<::std::path::Path>::as_ref($filename)
                        .to_string_lossy()
                )
            }),
        )
    };
}

/// Capitalize first letter
pub fn cap1(s: &str) -> String {
    let mut x = s.as_bytes().to_vec();
    let c = x[0].to_ascii_uppercase();
    x[0] = c;
    String::from_utf8(x.clone()).unwrap()
}

#[cfg(test)]
mod tests {
    use super::*;
    use vdj_ann::refx::RefData;

    // The following test checks for alignment of a D region.  This example was fixed by code
    // changes in March 2020.
    #[test]
    fn test_d_region_alignment() {
        use debruijn::dna_string::DnaString;
        use vdj_ann::annotate::annotate_seq;
        use vdj_ann::refx::make_vdj_ref_data_core;
        let seq = DnaString::from_acgt_bytes(
            b"GGAGGTGCGAATGACTCTGCTCTCTGTCCTGTCTCCTCATCTGCAAAATTAGGAAGCCTGTCTTGATTATCTCCAGGAA\
            CCTCCCACCTCTTCATTCCAGCCTCTGACAAACTCTGCACATTAGGCCAGGAGAAGCCCCCGAGCCAAGTCTCTTTTCTCATTCTC\
            TTCCAACAAGTGCTTGGAGCTCCAAGAAGGCCCCCTTTGCACTATGAGCAACCAGGTGCTCTGCTGTGTGGTCCTTTGTCTCCTGG\
            GAGCAAACACCGTGGATGGTGGAATCACTCAGTCCCCAAAGTACCTGTTCAGAAAGGAAGGACAGAATGTGACCCTGAGTTGTGAA\
            CAGAATTTGAACCACGATGCCATGTACTGGTACCGACAGGACCCAGGGCAAGGGCTGAGATTGATCTACTACTCACAGATAGTAAA\
            TGACTTTCAGAAAGGAGATATAGCTGAAGGGTACAGCGTCTCTCGGGAGAAGAAGGAATCCTTTCCTCTCACTGTGACATCGGCCC\
            AAAAGAACCCGACAGCTTTCTATCTCTGTGCCAGTAGTATTTTTCTTGCCGGGACAGGGGGCTGGAGCGGCACTGAAGCTTTCTTT\
            GGACAAGGCACCAGACTCACAGTTGTAGAGGACCTGAACAAGGTGTTCCCACCCGAGGTCGCTGTGTTTGAGCCATCAGA",
        );
        let (refx, ext_refx) = (human_ref(), String::new());
        let (is_tcr, is_bcr) = (true, false);
        let mut refdata = RefData::new();
        make_vdj_ref_data_core(&mut refdata, refx, &ext_refx, is_tcr, is_bcr, None);
        let ann = annotate_seq(&seq, &refdata, true, false, true);
        let mut have_d = false;
        for a in ann {
            if refdata.is_d(a.ref_id as usize) {
                have_d = true;
            }
        }
        assert!(have_d, "\nFailed to find alignment of D region.\n");
    }
}
