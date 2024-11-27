// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

//!
//! This crate defines the data structure that would represent the clonotypes
//! computed by enclone.
//!

use vdj_ann::annotate::Region as AnnRegion;

include!("./enclone.types.rs");

impl From<&bio_edit::alignment::Alignment> for Alignment {
    fn from(al: &bio_edit::alignment::Alignment) -> Self {
        Alignment {
            ref_start: al.ystart as u32,
            cigar: al.cigar(false),
        }
    }
}

fn make_ann_region(
    start: Option<u32>,
    end: Option<u32>,
    nt_seq: &[u8],
    v_start: u32,
    aa_seq: &[u8],
) -> Option<AnnRegion> {
    if start.is_some() && end.is_some() && start.unwrap() >= end.unwrap() {
        return None;
    }
    let v = v_start as usize;
    match (start, end) {
        (Some(s), Some(e)) => {
            let s = s as usize;
            let e = e as usize;
            Some(AnnRegion {
                start: s,
                stop: e,
                nt_seq: std::str::from_utf8(&nt_seq[s..e]).unwrap().to_string(),
                aa_seq: std::str::from_utf8(&aa_seq[(s - v) / 3..(e - v) / 3])
                    .unwrap()
                    .to_string(),
            })
        }
        _ => None,
    }
}

impl ClonotypeChain {
    pub fn fwr1_region(&self) -> Option<AnnRegion> {
        make_ann_region(
            self.fwr1_start,
            self.cdr1_start,
            &self.nt_sequence,
            self.v_start,
            &self.aa_sequence,
        )
    }
    pub fn cdr1_region(&self) -> Option<AnnRegion> {
        make_ann_region(
            self.cdr1_start,
            self.fwr2_start,
            &self.nt_sequence,
            self.v_start,
            &self.aa_sequence,
        )
    }
    pub fn fwr2_region(&self) -> Option<AnnRegion> {
        make_ann_region(
            self.fwr2_start,
            self.cdr2_start,
            &self.nt_sequence,
            self.v_start,
            &self.aa_sequence,
        )
    }
    pub fn cdr2_region(&self) -> Option<AnnRegion> {
        make_ann_region(
            self.cdr2_start,
            self.fwr3_start,
            &self.nt_sequence,
            self.v_start,
            &self.aa_sequence,
        )
    }
    pub fn fwr3_region(&self) -> Option<AnnRegion> {
        make_ann_region(
            self.fwr3_start,
            Some(self.cdr3_start),
            &self.nt_sequence,
            self.v_start,
            &self.aa_sequence,
        )
    }
    pub fn fwr4_region(&self) -> Option<AnnRegion> {
        make_ann_region(
            Some(self.cdr3_end),
            self.fwr4_end,
            &self.nt_sequence,
            self.v_start,
            &self.aa_sequence,
        )
    }

    pub fn cdr3_nt(&self) -> &[u8] {
        &self.nt_sequence[self.cdr3_start as usize..self.cdr3_end as usize]
    }
    pub fn cdr3_nt_string(&self) -> String {
        std::str::from_utf8(self.cdr3_nt()).unwrap().to_string()
    }
    pub fn cdr3_aa(&self) -> &[u8] {
        let start = (self.cdr3_start - self.v_start) / 3;
        let end = (self.cdr3_end - self.v_start) / 3;
        &self.aa_sequence[start as usize..end as usize]
    }
    pub fn cdr3_aa_string(&self) -> String {
        std::str::from_utf8(self.cdr3_aa()).unwrap().to_string()
    }
}

impl ExactSubClonotypeChain {
    pub fn fwr1_region(&self) -> Option<AnnRegion> {
        make_ann_region(
            self.fwr1_start,
            self.cdr1_start,
            &self.nt_sequence,
            self.v_start,
            &self.aa_sequence,
        )
    }
    pub fn cdr1_region(&self) -> Option<AnnRegion> {
        make_ann_region(
            self.cdr1_start,
            self.fwr2_start,
            &self.nt_sequence,
            self.v_start,
            &self.aa_sequence,
        )
    }
    pub fn fwr2_region(&self) -> Option<AnnRegion> {
        make_ann_region(
            self.fwr2_start,
            self.cdr2_start,
            &self.nt_sequence,
            self.v_start,
            &self.aa_sequence,
        )
    }
    pub fn cdr2_region(&self) -> Option<AnnRegion> {
        make_ann_region(
            self.cdr2_start,
            self.fwr3_start,
            &self.nt_sequence,
            self.v_start,
            &self.aa_sequence,
        )
    }
    pub fn fwr3_region(&self) -> Option<AnnRegion> {
        make_ann_region(
            self.fwr3_start,
            Some(self.cdr3_start),
            &self.nt_sequence,
            self.v_start,
            &self.aa_sequence,
        )
    }
    pub fn fwr4_region(&self) -> Option<AnnRegion> {
        make_ann_region(
            Some(self.cdr3_end),
            self.fwr4_end,
            &self.nt_sequence,
            self.v_start,
            &self.aa_sequence,
        )
    }
    pub fn cdr3_nt(&self) -> &[u8] {
        &self.nt_sequence[self.cdr3_start as usize..self.cdr3_end as usize]
    }
    pub fn cdr3_nt_string(&self) -> String {
        std::str::from_utf8(self.cdr3_nt()).unwrap().to_string()
    }
    pub fn cdr3_aa(&self) -> &[u8] {
        let start = (self.cdr3_start - self.v_start) / 3;
        let end = (self.cdr3_end - self.v_start) / 3;
        &self.aa_sequence[start as usize..end as usize]
    }
    pub fn cdr3_aa_string(&self) -> String {
        std::str::from_utf8(self.cdr3_aa()).unwrap().to_string()
    }
}
