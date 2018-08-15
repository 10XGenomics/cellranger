// Copyright (c) 2017 10x Genomics, Inc. All rights reserved.


//! This is documentation for the `mark_duplicates` crate.
//!
//! This crate is useful in marking PCR duplicates in a sorted BAM file.
//! PCR duplicates are reads that arise from sequencing multiple copies
//! of the same DNA fragment.
//!
//! # Algorithm
//!
//! Reads which have no mapping info
//! Both reads in a pair must be mapped, otherwise we drop the pair
//!
//! # Wishlist
//!
//! * Make the reporting generic. Allow the user to set up a custom
//! reporting function which accepts the duplicate grounps and creates
//! a report. 
//! * Instead of just accepting a barcode tag, allow user to specify a custom
//! function to group read by for dedup. This seems low priority as I do not
//! have a good use case right now

#![doc(html_logo_url =
    "https://upload.wikimedia.org/wikipedia/commons/thumb/4/4c/10x_Genomics_logo.svg/2000px-10x_Genomics_logo.svg.png")]
// #![deny(missing_docs,
//         missing_copy_implementations,
//         non_upper_case_globals,
//         trivial_casts,
//         trivial_numeric_casts,
//         unsafe_code,
//         unstable_features,
//         unused_extern_crates,
//         unused_import_braces,
//         unused_qualifications)]

extern crate rust_htslib;
extern crate itertools;
extern crate serde;
#[macro_use]
extern crate serde_derive;
extern crate serde_json;
extern crate fxhash;
extern crate tenkit;
extern crate num;

use std::path::Path;
use std::fs::File;
use itertools::Itertools;
use rust_htslib::bam;
use std::io::Write;

pub mod metrics;
use metrics::Histogram;
use metrics::Metric;

/// If the dup group gets extremely large we can run out of memory.
/// Process things in groups of 500K to prevent memory blow-up
/// May cause us to miss a few dups, but it doesn't really matter in these crazy regions
pub const MAX_DEDUP_SIZE: usize = 500_000;

pub struct MarkDup<'a> {
    barcode_tag: Option<&'a str>,
    predup_filter: fn(&bam::record::Record) -> bool,
    out_bam: Option<&'a mut bam::Writer>,
    report: MarkDupReport,
}

#[derive(Serialize, Deserialize, Default, Debug)]
struct MarkDupReport {
    description: String,
    all_duplicates: Histogram,
    optical_duplicates: Histogram,
    diffusion_duplicates: Histogram,
    non_proximal_duplicates: Histogram,
    inter_duplicate_quadrance: Histogram,
}

impl MarkDupReport {
    fn with_description<T>(description: T) -> Self where String: From<T> {
        let mut report = MarkDupReport::default();
        report.description = description.into();
        report
    }

    fn update(&mut self, dup_stat: &DupStat) {
        self.all_duplicates.increment(dup_stat.all);
        self.optical_duplicates.increment(dup_stat.optical);
        self.diffusion_duplicates.increment(dup_stat.diffusion);
        self.non_proximal_duplicates.increment(dup_stat.non_proximal);
        self.inter_duplicate_quadrance.merge(&dup_stat.quadrance);
    }
}

impl Metric for MarkDupReport {
    fn new() -> Self {
        MarkDupReport::default()
    }

    fn merge(&mut self, other: &Self) {
        assert!(self.description==other.description, 
            "Only reports(MarkDupReport) with same description can be merged");
        self.all_duplicates.merge(&other.all_duplicates);
        self.optical_duplicates.merge(&other.optical_duplicates);
        self.diffusion_duplicates.merge(&other.diffusion_duplicates);
        self.non_proximal_duplicates.merge(&other.non_proximal_duplicates);
        self.inter_duplicate_quadrance.merge(&other.inter_duplicate_quadrance);
    }
}


impl<'a> std::default::Default for MarkDup<'a> {
    #[inline]
    fn default() -> MarkDup<'a> {
        MarkDup {
            barcode_tag: None,
            predup_filter: |_| true,
            out_bam: None,
            report: MarkDupReport::with_description("no_filter_ignore_bcs")
        }
    }
}

impl<'a> MarkDup<'a> {
    fn dedup_records(&mut self, records: &mut Vec<bam::record::Record>,
        lane_coordinate_system: &LaneCoordinateSystem) {

        let dup_groups: Vec<Vec<usize>> = {
            let mut groups = Vec::new();
            let mut read_tuple: Vec<_> = records
                .iter()
                .enumerate()
                .filter(|&(_, x)| !(x.is_unmapped() || x.is_mate_unmapped()))
                .filter(|&(_, x)| (self.predup_filter)(x))
                .map(|(i, x)| {
                    ((if let Some(tag) = self.barcode_tag {
                          if let Some(bc) = x.aux(tag.as_bytes()) {
                              Some(bc.string())
                          } else {
                              None
                          }
                      } else {
                          None
                      },
                      x.is_first_in_template(),
                      x.is_reverse(),
                      x.tid(),
                      x.pos(),
                      x.mtid(),
                      x.mpos()),
                     i)
                })
                .collect();
            read_tuple.sort_by_key(|x| x.0);

            for (_, data) in &read_tuple.into_iter().group_by(|x| x.0) {
                groups.push(data.into_iter().map(|x| x.1).collect::<Vec<_>>());
            }
            groups
        };

        for group in &dup_groups {
            // If we are splitting on bcs, then only counts
            // stats for read groups with BCs?
            if let Some(tag) = self.barcode_tag {
                if records[group[0]].aux(tag.as_bytes()).is_none() {
                    continue;
                }
            }
            let dup_records: Vec<_> = group.iter().map(|&x| records.get(x).unwrap()).collect();
            let stat = DupStat::new(dup_records, lane_coordinate_system);
            self.report.update(&stat);
        }

        if let Some(ref mut out_bam) = self.out_bam {
            for group in &dup_groups {
                assert!(!records[group[0]].is_duplicate());
                for &g in group.iter().skip(1) {
                    records[g].set_duplicate();
                }
            }
            for r in records.iter_mut() {
                out_bam.write(r).expect("Failed to write record");
            }
        }
    }
}

pub struct MarkDupBuilder<'a> {
    barcode_tag: Option<&'a str>,
    predup_filter: fn(&bam::record::Record) -> bool,
    out_bam: Option<&'a mut bam::Writer>,
    description: Option<String>,
}

impl<'a> std::default::Default for MarkDupBuilder<'a> {
    #[inline]
    fn default() -> MarkDupBuilder<'a> {
        MarkDupBuilder {
            barcode_tag: None,
            predup_filter: |_| true,
            out_bam: None,
            description: None,
        }
    }
}

impl<'a> MarkDupBuilder<'a> {
    pub fn barcode_tag(mut self, tag: &'a str) -> Self {
        self.barcode_tag = Some(tag);
        self
    }
    pub fn predup_filter(mut self, filter: fn(&bam::record::Record) -> bool) -> Self {
        self.predup_filter = filter;
        self
    }
    pub fn out_bam(mut self, out: &'a mut bam::Writer) -> Self {
        self.out_bam = Some(out);
        self
    }

    pub fn description(mut self, desc: String) -> Self {
        self.description = Some(desc);
        self
    }

    pub fn build(self) -> MarkDup<'a> {
        let mut md = MarkDup::default();
        md.barcode_tag = self.barcode_tag;
        md.predup_filter = self.predup_filter;
        md.out_bam = self.out_bam;
        if let Some(desc) = self.description {
            md.report.description = desc;
        } else {
            panic!("ERROR: description is not initialized. Call description() before build()");
        }
        md
    }
}

// Marks duplicates
pub fn mark_duplicates<W: Write>(markdup_flavors: &mut Vec<MarkDup>,
                                 lane_coordinate_system: &LaneCoordinateSystem,
                                 bam: &mut bam::Reader,
                                 start: Option<i64>,
                                 end: Option<i64>,
                                 writer: &mut W) {

    let mut records_to_process = Vec::with_capacity(MAX_DEDUP_SIZE);
    let mut last_bam_key = (-1, -1);

    for r in bam.iter_chunk(start, end) {
        let mut record = r.unwrap();
        let curr_bam_key = (record.tid(), record.pos());
        if (curr_bam_key != last_bam_key) || (records_to_process.len() >= MAX_DEDUP_SIZE) {
            for md in markdup_flavors.iter_mut() {
                md.dedup_records(&mut records_to_process, lane_coordinate_system);
            }
            records_to_process.clear();
        }
        last_bam_key = curr_bam_key;
        records_to_process.push(record);
        if curr_bam_key == (-1, -1) || curr_bam_key == (-1, 0) {
            for md in markdup_flavors.iter_mut() {
                md.dedup_records(&mut records_to_process, lane_coordinate_system);
            }
            records_to_process.clear();
        }

    }

    let mut reports = Vec::new();
    for md in markdup_flavors.iter_mut() {
        md.dedup_records(&mut records_to_process, lane_coordinate_system);
        reports.push(&md.report)
    }
    serde_json::to_writer_pretty(writer, &reports)
        .expect("Could not write JSON in mark_duplicates");

}

pub fn join_reports<P: AsRef<Path>>(chunked_report_paths: Vec<P>, joined_reports_path: P) {

    let joined_report: Vec<MarkDupReport> = chunked_report_paths
        .iter()
        .map(|x| File::open(x).expect("Could not open file in join reports"))
        .map(|x| serde_json::from_reader(x).unwrap())
        .fold(Vec::new(), |mut summary, report: Vec<MarkDupReport>| {
            if summary.is_empty() {
                summary.extend(report);
            } else {
                for (s, r) in summary.iter_mut().zip(report.into_iter()) {
                    s.merge(&r);
                }
            }
            summary
        });

    let writer = File::create(joined_reports_path).expect("Could not open file");
    serde_json::to_writer_pretty(writer, &joined_report)
        .expect("Could not write JSON in mark_duplicates");
}




pub const MAX_DIFFUSION_DUP_DISTANCE: u32 = 25_000;
const SQ_MAX_DIFFUSION_DUP_DISTANCE: u32 = MAX_DIFFUSION_DUP_DISTANCE * MAX_DIFFUSION_DUP_DISTANCE;
pub const OPTICAL_DUPLICATE_DISTANCE: u32 = 100;
const SQ_OPTICAL_DUPLICATE_DISTANCE: u32 = OPTICAL_DUPLICATE_DISTANCE * OPTICAL_DUPLICATE_DISTANCE;

use tenkit::lane::{ReadLoc, LaneCoordinateSystem};
use tenkit::collections::FxHashSet;
use std::cmp::{min, max};

struct DupStat {
    all: i64,
    optical: i64,
    diffusion: i64,
    non_proximal: i64,
    quadrance: Histogram
}

impl DupStat {
    pub fn new(dup_records: Vec<&bam::record::Record>,
            lane_coordinate_system: &LaneCoordinateSystem) -> Self {

        let mut read_locs: Vec<_> = dup_records
            .iter()
            .filter_map(|r| ReadLoc::from_bam_record(r).map(|l| (l, r)))
            .collect();

        read_locs.sort_by_key(|&(ref l, _)| l.to_flowcell_lane_key());

        let mut diff_dups_count = 0i64;
        let mut optical_dups_count = 0i64;
        let mut inter_dup_quadrance = Histogram::new();

        for (_, data) in
            &read_locs
                .into_iter()
                .group_by(|&(ref l, _)| l.to_flowcell_lane_key()) {

            let grouped_read_locs: Vec<_> = data.into_iter().collect();
            let layout = lane_coordinate_system.get_layout_for_read_loc(&grouped_read_locs[0].0);
            let test_diff_dups = layout.has_diffusion_duplicates(MAX_DIFFUSION_DUP_DISTANCE);

            let lane_pos: Vec<_> = grouped_read_locs
                .iter()
                .filter_map(|&(ref l, _)| lane_coordinate_system.convert_to_lane_coords(&l))
                .collect();
                
            let cmp_reads = min(200, lane_pos.len());

            let mut diff_dups = FxHashSet::default();
            let mut optical_dups = FxHashSet::default();

            for i in 0..cmp_reads {
                let posi = &lane_pos[i];
                for j in (i+1)..lane_pos.len() {
                    let posj = &lane_pos[j];
                    let quadrance = posi.quadrance(posj);
                    if test_diff_dups &&  quadrance < SQ_MAX_DIFFUSION_DUP_DISTANCE {
                        diff_dups.insert(j);
                    }
                    if quadrance < SQ_OPTICAL_DUPLICATE_DISTANCE {
                        optical_dups.insert(j);
                    }
                    inter_dup_quadrance.increment(quadrance);
                }
            }
            diff_dups_count += diff_dups.len() as i64;
            optical_dups_count += optical_dups.len() as i64;
        }

        let total_dups_count = dup_records.len() as i64;
        // diffusion dups encompass optical dups, if
        let non_proximal_dups_count = total_dups_count - max(optical_dups_count, diff_dups_count);
        DupStat {
            all: total_dups_count,
            optical: optical_dups_count,
            diffusion: diff_dups_count,
            non_proximal: non_proximal_dups_count,
            quadrance: inter_dup_quadrance
        }

    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {}
}
