// This captures and analyzes data about barcodes.

#![allow(clippy::many_single_char_names)]
// TODO: fix these.
#![allow(clippy::needless_range_loop)]

use crate::analyze_vs::analyze_vs;
use crate::constants::{CHAIN_TYPES, CHAIN_TYPESX, PRIMER_EXT_LEN};
use crate::filter_barcodes::BarcodeCellInfo;
use crate::filter_log::{AsmCellFilter, FilterLogEntry, FilterLogger, FilterSwitch};
use debruijn::dna_string::DnaString;
use itertools::Itertools;
use num_traits::ToPrimitive;
use ordered_float::NotNan;
use serde::Deserialize;
use stats_utils::{len_weighted_mean, mean, n50, n90};
use std::cmp::{max, min};
use std::fs::File;
use std::io::{BufWriter, Write};
use string_utils::strme;
use tables::print_tabular;
use tenkit2::pack_dna::{reverse_complement, unpack_bases_80};
use vdj_ann::refx::RefData;
use vector_utils::{
    bin_member, bin_position, bin_position1_3, contains_at, erase_if, lower_bound, next_diff,
    next_diff12_3, next_diff1_2, next_diff1_5, reverse_sort, unique_sort, upper_bound,
};

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
// BARCODE DATA STRUCTURES
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

#[derive(PartialEq, PartialOrd, Clone, Serialize, Deserialize)]
pub struct BarcodeData {
    // the barcode
    pub barcode: String,

    // number of reads in the barcode
    pub nreads: i32,

    // number of reads hitting TRAV, TRBV, TRGV, and TRDV
    pub trav_reads: i32,
    pub trbv_reads: i32,
    pub trgv_reads: i32,
    pub trdv_reads: i32,

    // number of reads hitting TRAV..TRAJ, TRBV..TRBJ, TRGV..TRGJ, and TRDV..TRDJ
    pub travj_reads: i32,
    pub trbvj_reads: i32,
    pub trgvj_reads: i32,
    pub trdvj_reads: i32,

    // total number of umis
    pub total_ucounts: usize,

    // sorted list of counts, one per umi, showing the number of read pairs
    // assigned to the umi, and excluding umis having only one read pair
    pub ucounts: Vec<i32>,

    // surviving nonsolo ucounts, and their ids
    pub xucounts: Vec<i32>,
    pub xuids: Vec<i32>,

    // one entry (npairs,chain-type,umi) for each umi
    pub umi_info: Vec<(i32, i8, String)>,

    // number of reads whose umi was error corrected
    pub nreads_umi_corrected: i32,

    // chain type of a sample of reads:
    // 15 entries, giving counts for
    // fw.IGH, fw.IGK, fw.IGL, fw.TRA, fw.TRB, fw.TRD, fw.TRG,
    // rc.IGH, rc.IGK, rc.IGL, rc.TRA, rc.TRB, rc.TRD, rc.TRG, None.
    // For a sample of reads, each is assigned to one of these categories,
    // and the total in each category is given here.  Because the orientation
    // of reads is normalized upstream of this, ideal behavior is for all to be fw.
    pub chain_sample: Vec<i32>,

    // number of contigs
    pub ncontigs: i32,

    // indices of reference segments appearing in good contigs
    pub good_refseqs: Vec<i32>,

    // stats for tracking primer counts in reads
    pub primer_sample_count: u16,
    pub inner_hit: Vec<u16>, // one entry for each inner primer
    pub inner_hit_good: Vec<u16>,
    pub outer_hit: Vec<u16>, // one entry for each inner primer
    pub outer_hit_good: Vec<u16>,

    // stats for tracking primer counts in good contigs
    pub inner_hit_good_contigs: Vec<u16>,
    pub outer_hit_good_contigs: Vec<u16>,
    // Fraction of reads in barcode used for assembly
    pub frac: f64,
}

impl BarcodeData {
    pub fn new() -> BarcodeData {
        BarcodeData {
            barcode: String::new(),
            nreads: 0_i32,
            trav_reads: 0_i32,
            trbv_reads: 0_i32,
            trgv_reads: 0_i32,
            trdv_reads: 0_i32,
            travj_reads: 0_i32,
            trbvj_reads: 0_i32,
            trgvj_reads: 0_i32,
            trdvj_reads: 0_i32,
            total_ucounts: 0_usize,
            ucounts: Vec::<i32>::new(),
            xucounts: Vec::<i32>::new(),
            xuids: Vec::<i32>::new(),
            umi_info: Vec::<(i32, i8, String)>::new(),
            nreads_umi_corrected: 0_i32,
            chain_sample: vec![0_i32; 15],
            ncontigs: 0_i32,
            good_refseqs: Vec::<i32>::new(),
            primer_sample_count: 0,
            inner_hit: Vec::<u16>::new(),
            inner_hit_good: Vec::<u16>::new(),
            outer_hit: Vec::<u16>::new(),
            outer_hit_good: Vec::<u16>::new(),
            inner_hit_good_contigs: Vec::<u16>::new(),
            outer_hit_good_contigs: Vec::<u16>::new(),
            frac: 1.0_f64,
        }
    }
}

impl Default for BarcodeData {
    fn default() -> Self {
        Self::new()
    }
}

#[derive(PartialEq, Eq, PartialOrd, Ord, Clone, Serialize, Deserialize, Debug)]
pub struct ContigJunctionData {
    // compressed 80 bp of junction sequence, stopping at the end of the J segment
    pub jxn_seq: [u8; 20],

    // umi support for junction (capped at 65535)
    pub umis: u16,

    // high confidence barcode?
    pub high_confidence: bool,

    // igh contig?
    pub is_igh: bool,
}

#[derive(PartialEq, Eq, PartialOrd, Ord, Clone, Serialize, Deserialize, Debug)]
pub struct ContigChimeraData {
    // cdr3_nt, assuming one was found
    pub cdr3: Vec<u8>,

    // v segment id, assuming only one
    pub v_ref_id: usize,

    // num umis assigned to this contig
    pub umi_count: usize,

    // is this a productive contig
    pub productive: bool,

    // the barcode
    pub barcode: String,
}

#[derive(PartialEq, PartialOrd, Clone, Serialize, Deserialize, Debug)]
pub struct BarcodeDataBrief {
    // the barcode
    pub barcode: String,

    // number of read pairs
    pub read_pairs: u64,

    // total number of umis
    pub total_ucounts: usize,

    // surviving nonsolo ucounts,
    pub xucounts: Vec<i32>,

    // number of contigs
    pub ncontigs: usize,

    // Fraction of reads in barcode used for assembly
    pub frac_reads_used: f64,
}

impl BarcodeDataBrief {
    pub fn new() -> BarcodeDataBrief {
        BarcodeDataBrief {
            barcode: String::new(),
            read_pairs: 0_u64,
            total_ucounts: 0,
            xucounts: Vec::<i32>::new(),
            ncontigs: 0,
            frac_reads_used: 1.0_f64,
        }
    }
}

impl Default for BarcodeDataBrief {
    fn default() -> Self {
        Self::new()
    }
}

#[derive(PartialEq, PartialOrd, Clone, Serialize, Deserialize)]
pub struct BarcodeDataSum {
    pub nreads: i64,                   // number of reads
    pub nreads_umi_corrected: i64,     // number of reads whose umi was corrected
    pub ncontigs: i64,                 // number of contigs
    pub xucount_sum: i64,              // sum over xucounts
    pub chain_sample: Vec<i64>,        // chain sample
    pub adj_chain_sample: Vec<f64>, // adjusted chain_sample (summation done after scaling up based on frac)
    pub umi_dist: Vec<(i8, i32, i32)>, // (chain-type,npairs,multiplicity)
    pub cs_good: Vec<usize>,        // for each constant region, count in good
    // contigs
    pub inner_hit_total: Vec<NotNan<f64>>,
    pub inner_hit_good_total: Vec<NotNan<f64>>,
    pub outer_hit_total: Vec<NotNan<f64>>,
    pub outer_hit_good_total: Vec<NotNan<f64>>,
    pub inner_hit_good_contigs_total: Vec<usize>,
    pub outer_hit_good_contigs_total: Vec<usize>,
}

impl BarcodeDataSum {
    pub fn new() -> BarcodeDataSum {
        BarcodeDataSum {
            nreads: 0_i64,
            nreads_umi_corrected: 0_i64,
            ncontigs: 0_i64,
            xucount_sum: 0_i64,
            chain_sample: vec![0_i64; 15],
            adj_chain_sample: vec![0_f64; 15],
            umi_dist: Vec::<(i8, i32, i32)>::new(),
            cs_good: Vec::<usize>::new(),
            inner_hit_total: Vec::<NotNan<f64>>::new(),
            inner_hit_good_total: Vec::<NotNan<f64>>::new(),
            outer_hit_total: Vec::<NotNan<f64>>::new(),
            outer_hit_good_total: Vec::<NotNan<f64>>::new(),
            inner_hit_good_contigs_total: Vec::<usize>::new(),
            outer_hit_good_contigs_total: Vec::<usize>::new(),
        }
    }
    pub fn sum(d: &[BarcodeData], refdata: &RefData) -> BarcodeDataSum {
        BarcodeDataSum {
            nreads: {
                let mut s = 0_i64;
                for x in d {
                    s += x.nreads as i64;
                }
                s
            },
            nreads_umi_corrected: {
                let mut s = 0_i64;
                for x in d {
                    s += x.nreads_umi_corrected as i64;
                }
                s
            },
            ncontigs: {
                let mut s = 0_i64;
                for x in d {
                    s += x.ncontigs as i64;
                }
                s
            },
            xucount_sum: {
                let mut s = 0_i64;
                for x in d {
                    s += x.xucounts.iter().sum::<i32>() as i64;
                }
                s
            },
            chain_sample: {
                let mut s = vec![0_i64; 15];
                for x in d {
                    for i in 0..15 {
                        s[i] += x.chain_sample[i] as i64;
                    }
                }
                s
            },
            // When summing barcode objects, adjust each chain_sample value based on frac
            adj_chain_sample: {
                let mut s = vec![0_f64; 15];
                for x in d {
                    for i in 0..15 {
                        if x.frac != 0.0 {
                            s[i] += x.chain_sample[i] as f64 / x.frac;
                        } else {
                            s[i] += 0.0;
                        }
                    }
                }
                s
            },
            umi_dist: {
                let mut t = Vec::<(i32, i8)>::new();
                for x in d {
                    for i in 0..x.umi_info.len() {
                        t.push((x.umi_info[i].0, x.umi_info[i].1));
                    }
                }
                t.sort_unstable();
                let mut s = Vec::<(i8, i32, i32)>::new();
                let mut i = 0;
                while i < t.len() {
                    let j = next_diff(&t, i);
                    s.push((t[i].1, t[i].0, (j - i) as i32));
                    i = j;
                }
                s
            },
            cs_good: {
                let mut s = vec![0; refdata.cs.len()];
                for x in d {
                    for y in &x.good_refseqs {
                        let p = bin_position(&refdata.cs, &(*y as usize));
                        if p >= 0 {
                            s[p as usize] += 1;
                        }
                    }
                }
                s
            },
            inner_hit_total: {
                // in single-end case, remember to divide by two later, and below
                if d.is_empty() {
                    Vec::<NotNan<f64>>::new()
                } else {
                    let mut s = vec![0.0; d[0].inner_hit.len()];
                    for x in d {
                        for i in 0..x.inner_hit.len() {
                            s[i] += x.inner_hit[i] as f64 * x.nreads as f64
                                / x.primer_sample_count as f64;
                        }
                    }
                    let mut t = Vec::<NotNan<f64>>::new();
                    for i in 0..s.len() {
                        t.push(NotNan::new(s[i]).unwrap());
                    }
                    t
                }
            },
            inner_hit_good_total: {
                if d.is_empty() {
                    Vec::<NotNan<f64>>::new()
                } else {
                    let mut s = vec![0.0; d[0].inner_hit_good.len()];
                    for x in d {
                        for i in 0..x.inner_hit_good.len() {
                            s[i] += x.inner_hit_good[i] as f64 * x.nreads as f64
                                / x.primer_sample_count as f64;
                        }
                    }
                    let mut t = Vec::<NotNan<f64>>::new();
                    for i in 0..s.len() {
                        t.push(NotNan::new(s[i]).unwrap());
                    }
                    t
                }
            },
            outer_hit_total: {
                if d.is_empty() {
                    Vec::<NotNan<f64>>::new()
                } else {
                    let mut s = vec![0.0; d[0].outer_hit.len()];
                    for x in d {
                        for i in 0..x.outer_hit.len() {
                            s[i] += x.outer_hit[i] as f64 * x.nreads as f64
                                / x.primer_sample_count as f64;
                        }
                    }
                    let mut t = Vec::<NotNan<f64>>::new();
                    for i in 0..s.len() {
                        t.push(NotNan::new(s[i]).unwrap());
                    }
                    t
                }
            },
            outer_hit_good_total: {
                if d.is_empty() {
                    Vec::<NotNan<f64>>::new()
                } else {
                    let mut s = vec![0.0; d[0].outer_hit_good.len()];
                    for x in d {
                        for i in 0..x.outer_hit_good.len() {
                            s[i] += x.outer_hit_good[i] as f64 * x.nreads as f64
                                / x.primer_sample_count as f64;
                        }
                    }
                    let mut t = Vec::<NotNan<f64>>::new();
                    for i in 0..s.len() {
                        t.push(NotNan::new(s[i]).unwrap());
                    }
                    t
                }
            },
            inner_hit_good_contigs_total: {
                if d.is_empty() {
                    Vec::<usize>::new()
                } else {
                    let mut s = vec![0; d[0].inner_hit_good_contigs.len()];
                    for x in d {
                        for i in 0..x.inner_hit_good_contigs.len() {
                            s[i] += x.inner_hit_good_contigs[i] as usize;
                        }
                    }
                    s
                }
            },
            outer_hit_good_contigs_total: {
                if d.is_empty() {
                    Vec::<usize>::new()
                } else {
                    let mut s = vec![0; d[0].outer_hit_good_contigs.len()];
                    for x in d {
                        for i in 0..x.outer_hit_good_contigs.len() {
                            s[i] += x.outer_hit_good_contigs[i] as usize;
                        }
                    }
                    s
                }
            },
        }
    }
    pub fn sumsum(d: &[BarcodeDataSum], refdata: &RefData) -> BarcodeDataSum {
        BarcodeDataSum {
            nreads: {
                let mut s = 0_i64;
                for x in d {
                    s += x.nreads;
                }
                s
            },
            nreads_umi_corrected: {
                let mut s = 0_i64;
                for x in d {
                    s += x.nreads_umi_corrected;
                }
                s
            },
            ncontigs: {
                let mut s = 0_i64;
                for x in d {
                    s += x.ncontigs;
                }
                s
            },
            xucount_sum: {
                let mut s = 0_i64;
                for x in d {
                    s += x.xucount_sum;
                }
                s
            },
            chain_sample: {
                let mut s = vec![0_i64; 15];
                for x in d {
                    for i in 0..15 {
                        s[i] += x.chain_sample[i];
                    }
                }
                s
            },
            // When summing barcodesum objects, simply add adjusted chain samples
            adj_chain_sample: {
                let mut s = vec![0_f64; 15];
                for x in d {
                    for i in 0..15 {
                        s[i] += x.adj_chain_sample[i];
                    }
                }
                s
            },
            umi_dist: {
                let mut t = Vec::<(i8, i32, i32)>::new();
                for x in d {
                    t.append(&mut x.umi_dist.clone());
                }
                t.sort_unstable();
                let mut s = Vec::<(i8, i32, i32)>::new();
                let mut i = 0;
                while i < t.len() {
                    let j = next_diff12_3(&t, i as i32) as usize;
                    let mut mult = 0;
                    for k in i..j {
                        mult += t[k].2;
                    }
                    s.push((t[i].0, t[i].1, mult));
                    i = j;
                }
                s
            },
            cs_good: {
                let mut s = vec![0; refdata.cs.len()];
                for x in d {
                    for i in 0..s.len() {
                        s[i] += x.cs_good[i];
                    }
                }
                s
            },
            inner_hit_total: {
                let mut s = Vec::<NotNan<f64>>::new();
                for x in d {
                    if !x.inner_hit_total.is_empty() {
                        if s.is_empty() {
                            s = x.inner_hit_total.clone();
                        } else {
                            for i in 0..x.inner_hit_total.len() {
                                s[i] += x.inner_hit_total[i];
                            }
                        }
                    }
                }
                s
            },
            inner_hit_good_total: {
                let mut s = Vec::<NotNan<f64>>::new();
                for x in d {
                    if !x.inner_hit_good_total.is_empty() {
                        if s.is_empty() {
                            s = x.inner_hit_good_total.clone();
                        } else {
                            for i in 0..x.inner_hit_good_total.len() {
                                s[i] += x.inner_hit_good_total[i];
                            }
                        }
                    }
                }
                s
            },
            outer_hit_total: {
                let mut s = Vec::<NotNan<f64>>::new();
                for x in d {
                    if !x.outer_hit_total.is_empty() {
                        if s.is_empty() {
                            s = x.outer_hit_total.clone();
                        } else {
                            for i in 0..x.outer_hit_total.len() {
                                s[i] += x.outer_hit_total[i];
                            }
                        }
                    }
                }
                s
            },
            outer_hit_good_total: {
                let mut s = Vec::<NotNan<f64>>::new();
                for x in d {
                    if !x.outer_hit_good_total.is_empty() {
                        if s.is_empty() {
                            s = x.outer_hit_good_total.clone();
                        } else {
                            for i in 0..x.outer_hit_good_total.len() {
                                s[i] += x.outer_hit_good_total[i];
                            }
                        }
                    }
                }
                s
            },
            inner_hit_good_contigs_total: {
                let mut s = Vec::<usize>::new();
                for x in d {
                    if !x.inner_hit_good_contigs_total.is_empty() {
                        if s.is_empty() {
                            s = x.inner_hit_good_contigs_total.clone();
                        } else {
                            for i in 0..x.inner_hit_good_contigs_total.len() {
                                s[i] += x.inner_hit_good_contigs_total[i];
                            }
                        }
                    }
                }
                s
            },
            outer_hit_good_contigs_total: {
                let mut s = Vec::<usize>::new();
                for x in d {
                    if !x.outer_hit_good_contigs_total.is_empty() {
                        if s.is_empty() {
                            s = x.outer_hit_good_contigs_total.clone();
                        } else {
                            for i in 0..x.outer_hit_good_contigs_total.len() {
                                s[i] += x.outer_hit_good_contigs_total[i];
                            }
                        }
                    }
                }
                s
            },
        }
    }
}

impl Default for BarcodeDataSum {
    fn default() -> Self {
        Self::new()
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
// COMPUTE METRICS YIELDING JSON STRING
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

fn write_json_metric_str(metric_name: &str, value: &str, json: &mut Vec<u8>) {
    fwriteln!(json, "    \"{}\": \"{}\",", metric_name, value);
}

pub fn write_json_metric_f64(metric_name: &str, value: f64, json: &mut Vec<u8>) {
    if value.is_nan() {
        fwriteln!(json, "    \"{}\": null,", metric_name);
    } else {
        fwriteln!(json, "    \"{}\": {},", metric_name, value);
    }
}

pub fn write_json_metric_ratio(metric_name: &str, num: usize, den: usize, json: &mut Vec<u8>) {
    let metric = if den != 0 {
        num as f64 / den as f64
    } else {
        0.0
    };
    fwriteln!(json, "    \"{}\": {},", metric_name, metric);
}

#[allow(clippy::too_many_arguments)]
pub fn metrics_json(
    dir: &str,
    is_tcr: bool,
    d: &[BarcodeDataSum],
    json: &mut Vec<u8>,
    single_end: bool,
    refdata: &RefData,
    inner_primers: &[Vec<u8>],
    outer_primers: &[Vec<u8>],
    npairs: usize,
    report: &mut BufWriter<File>,
    is_gd: Option<bool>,
) {
    // Get number of reads.

    let mut nreads = npairs;
    if !single_end {
        nreads *= 2;
    }

    // Sum.

    let x = BarcodeDataSum::sumsum(d, refdata);

    // Add opening bracket.

    json.append(&mut b"{\n".to_vec());

    // Analyze primer hits.

    let mut log = Vec::<u8>::new();
    analyze_primer_hits(
        &x,
        refdata,
        inner_primers,
        outer_primers,
        single_end,
        json,
        &mut log,
    );

    // Write metrics productive_contigs_with_x, where x is a constant region name.

    let mut cnames = Vec::<String>::new();
    for i in 0..refdata.cs.len() {
        let j = refdata.cs[i];
        cnames.push(refdata.name[j].clone());
    }
    unique_sort(&mut cnames);
    let mut cc = vec![0_i64; cnames.len()];
    for i in 0..x.cs_good.len() {
        let j = refdata.cs[i];
        let p = bin_position(&cnames, &refdata.name[j]);
        cc[p as usize] += x.cs_good[i] as i64;
    }
    for i in 0..cc.len() {
        let metric_name = format!("productive_contigs_with_{}", cnames[i]);
        fwriteln!(json, "    \"{}\": {},", metric_name, cc[i]);
    }

    // Compute "vdj_corrected_umi_frac", which is the fraction of reads whose UMI
    // was corrected.

    write_json_metric_ratio(
        "vdj_corrected_umi_frac",
        x.nreads_umi_corrected as usize,
        nreads,
        json,
    );

    // Generate {IGH,IGK,IGL,TRA,TRB,TRD,TRG,multi}
    // _vdj_recombinome_readpairs_per_umi_distribution metrics.  Note that we put
    // the distribution in order by numerical count, as opposed to ordering
    // alphabetically, which is how it was done before.

    let mut chains = CHAIN_TYPES.to_vec();
    chains.push("multi");
    let mut all = Vec::<(i32, i32)>::new();

    let mut ud = Vec::<(i32, i32)>::new();
    for i in 0..x.umi_dist.len() {
        if x.umi_dist[i].0 != 1 {
            ud.push((x.umi_dist[i].1, x.umi_dist[i].2));
        }
    }
    ud.sort_unstable();
    let mut i = 0;
    while i < ud.len() {
        let j = next_diff1_2(&ud, i as i32) as usize;
        let mut s = 0;
        for k in i..j {
            s += ud[k].1;
        }
        all.push((ud[i].0, s));
        i = j;
    }
    for c in 0..8 {
        if c == 1 {
            continue;
        }
        fwriteln!(
            json,
            "    \"{}_vdj_recombinome_readpairs_per_umi_distribution\": {{\n",
            chains[c]
        );
        for i in 0..x.umi_dist.len() {
            if x.umi_dist[i].0 != c as i8 {
                continue;
            }
            fwrite!(json, "        \"{}\": {}", x.umi_dist[i].1, x.umi_dist[i].2);
            if i < x.umi_dist.len() - 1 && x.umi_dist[i].0 == x.umi_dist[i + 1].0 {
                fwrite!(json, ",");
            }
            fwriteln!(json, "");
        }
        fwriteln!(json, "}},");
    }
    fwriteln!(
        json,
        "    \"multi_vdj_recombinome_readpairs_per_umi_distribution\": {{\n"
    );
    for i in 0..all.len() {
        fwrite!(json, "        \"{}\": {}", all[i].0, all[i].1);
        if i < all.len() - 1 {
            fwrite!(json, ",");
        }
        fwriteln!(json, "");
    }
    fwriteln!(json, "}},");

    // Generate {IGH,IGK,IGL,TRA,TRB,TRD,TRG,multi}
    // _vdj_recombinome_antisense_reads_frac metrics.

    let (mut chain_total_fw, mut chain_total_rc) = (0, 0);
    for i in 0..CHAIN_TYPESX.len() {
        let (fw, rc) = (x.chain_sample[i] as usize, x.chain_sample[i + 7] as usize);
        write_json_metric_ratio(
            &format!("{}_vdj_recombinome_antisense_reads_frac", CHAIN_TYPESX[i]),
            rc,
            fw + rc,
            json,
        );
        chain_total_fw += fw;
        chain_total_rc += rc;
    }
    write_json_metric_ratio(
        "multi_vdj_recombinome_antisense_reads_frac",
        chain_total_rc,
        chain_total_fw + chain_total_rc,
        json,
    );

    // Generate {IGH,IGK,IGL,TRA,TRB,TRD,TRG,multi}
    // _vdj_recombinome_mapped_reads_frac metrics.

    // Using adjusted chain_sample to account for capping number of reads used in assembly
    let adj_chain_total = x.adj_chain_sample.iter().sum::<f64>();
    for i in 0..CHAIN_TYPESX.len() {
        write_json_metric_ratio(
            &format!("{}_vdj_recombinome_mapped_reads_frac", CHAIN_TYPESX[i]),
            (x.adj_chain_sample[i] + x.adj_chain_sample[i + 7]) as usize,
            adj_chain_total as usize,
            json,
        );
    }
    write_json_metric_ratio(
        "multi_vdj_recombinome_mapped_reads_frac",
        (adj_chain_total - x.adj_chain_sample[14]) as usize,
        adj_chain_total as usize,
        json,
    );

    // Compute vdj_sequencing_efficiency.  For this, define the number of
    // assemblable read pairs to be the total number of read pairs appearing in
    // xucounts.  Then divide by the total number of read pairs.

    write_json_metric_ratio(
        "vdj_sequencing_efficiency",
        x.xucount_sum as usize,
        npairs,
        json,
    );

    // Add error rate metrics.

    if is_tcr {
        let o = [dir].to_vec();
        analyze_vs(
            o, "unknown", refdata, true, false, false, false, true, true, json, &mut log, is_gd,
        );
    }
    fwrite!(report, "{}", strme(&log));

    // Remove trailing comma.  Add closing bracket.

    let n = json.len();
    json.truncate(n - 2);
    json.extend(b"\n}\n");
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
// ANALYZE BRIEF BARCODE DATA
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Identify productive pairs that should be untrusted based on their relationship
// to large clones.

pub fn analyze_barcode_data_brief(
    d: &[BarcodeCellInfo],
    FilterSwitch {
        asm_shared_contig, ..
    }: FilterSwitch,
    kills: &mut Vec<String>,
    killsc: &mut Vec<String>,
    mut filter_logger: Option<&mut FilterLogger>,
) {
    println!("\nUNTRUSTED CONTIGS");

    chimeric_filters(d, kills, &mut filter_logger);

    if asm_shared_contig {
        junction_filters(d, killsc, kills, &mut filter_logger);
        common_clone_filters(d, killsc, kills, &mut filter_logger);
    }

    unique_sort(kills);
    unique_sort(killsc);
    println!();
}

fn junction_filters(
    d: &[BarcodeCellInfo],
    killsc: &mut Vec<String>,
    kills: &mut Vec<String>,
    filter_logger: &mut Option<&mut FilterLogger>,
) {
    // Kill contigs that appear to arise from some sort of leakage from plasma cells.  We identify
    // these by seeing UMI counts that are much smaller than a dominant UMI count, and for which
    // the median UMI count is very small.  One reason we think this filtering is correct is that
    // we observed large clusters arising in a library, and completely absent from a parallel
    // library made from the same cell lot.  (And this happened multiple times.)  Note that
    // extreme differences in UMI count appear to be biologically possible (and we observe such),
    // so that alone is not characteristic.  Note that there is no mechanism here to address the
    // case of fake plasma clonal expansions arising where *only* background is present.
    //
    // This is a messy problem as we have diverse and relatively limited data.
    //
    // Cases that are fixed now (all GEM-55):
    // fake expansion  good from same cell lot   fake clone
    // 124960_100_per  124952_100_per            IGH:CAREYPTSYGSGTYYVSPAPFDSW;IGK:CQRYTMSPFISF
    // 124552_100_per  124548_100_per            IGH:CAKSGAGEIGEYYFGYW;IGL:CQVWDSTSDHRWVF
    // 124553_100_per  124549_100_per.
    //
    // Cases that appear not to be fixable because UMI counts are too similar (all GEM-U, CR 3.0):
    // fake expansion  good from same cell lot   fake clone
    // 79209           79210                     IGH:CAKHDYSNPQW;IGK:CFQGSHVPFTF
    // 86356           86355                     IGH:CVRVVEGTSAYDIW;IGL:CTSYTSSSTYVF
    // 86228           86227.                    IGH:CARQSDTGYFEFW;IGL:CQVWDSSTDHPIF
    //
    // Other cases that might hypothetically be fixable:
    // fake expansion  good from same cell lot   fake clone
    // 74447           74446
    // 79211           79212                     IGH:CARWGGSSNDYW;IGK:CQQHYSTPYTF
    //
    // Negative control (case where damage might be done by this algorithm):
    // 77225.  Compare replicate 77226;
    // also matched PBMC 77223,77224 and matched splenocytes 77221,77222.
    let mut all = Vec::<([u8; 20], u16, bool, usize, usize)>::new();
    for i in 0..d.len() {
        for (j, jundata) in d[i].jundata.iter().enumerate() {
            // (junction segment of contig, #umis, confident, index)
            all.push((jundata.jxn_seq, jundata.umis, jundata.high_confidence, i, j));
        }
    }
    all.sort_unstable();
    let mut i = 0;
    const MIN_RATIO_UMI: usize = 40;
    const MAX_MEDIAN: u16 = 1;
    const MIN_CLUSTER: usize = 10;
    // was 20 in experiment
    while i < all.len() {
        let j = next_diff1_5(&all, i as i32) as usize;

        // Now i..j is a group of entries in all, each with the same junction segment.

        let median = i + (j - i) / 2;
        if j - i >= MIN_CLUSTER && all[median].1 <= MAX_MEDIAN {
            let mut alts = Vec::<[u8; 20]>::new();
            for k in i..j {
                let u = all[k].3;
                for l in 0..d[u].jundata.len() {
                    if d[u].jundata[l].jxn_seq != all[i].0 {
                        alts.push(d[u].jundata[l].jxn_seq);
                    }
                }
            }
            alts.sort_unstable();
            let mut max_alt = 0;
            let mut r = 0;
            let mut alt_counts = Vec::<usize>::new();
            while r < alts.len() {
                let s = next_diff(&alts, r);
                max_alt = max(max_alt, s - r);
                alt_counts.push(s - r);
                r = s;
            }
            reverse_sort(&mut alt_counts);
            //if
            /* max_alt >= MIN_CLUSTER */
            //0 == 0
            {
                for k in i..j {
                    if all[j - 1].1 as usize >= MIN_RATIO_UMI * max(1, all[k].1 as usize) {
                        let m = all[k].3;
                        let x = &d[m];
                        for j in 0..x.jundata.len() {
                            if x.jundata[j].high_confidence {
                                println!(
                                    "{}: contig {} = possible plasma cell leakage",
                                    x.barcode,
                                    j + 1
                                );
                            }
                            killsc.push(format!("{}_contig_{}", x.barcode, j + 1));
                        }
                        println!("{} = possible plasma cell leakage", x.barcode);
                        println!("alt counts = {}", alt_counts.iter().format(","));
                        kills.push(x.barcode.clone());
                        if let Some(ref mut logger) = filter_logger {
                            logger.log(&FilterLogEntry::cell_calling(
                                x.barcode.clone(),
                                AsmCellFilter::NonDominantJunction {
                                    contig: format!("{}_contig_{}", x.barcode, all[k].4 + 1),
                                    junction_umis: all[k].1 as usize,
                                    param_min_umi_ratio: MIN_RATIO_UMI,
                                    dominant_contig: format!(
                                        "{}_contig_{}",
                                        d[all[j - 1].3].barcode,
                                        all[j - 1].4 + 1
                                    ),
                                    dominant_junction_umis: all[j - 1].1 as usize,
                                    cluster_size: j - i,
                                    param_min_cluster_size: MIN_CLUSTER,
                                    cluster_median_junction_umis: all[median].1,
                                    param_max_median_junction_umis: MAX_MEDIAN,
                                },
                            ));
                        }
                    }
                }
            }
        }
        i = j;
    }
    // Address a different version of plasma cell leakage.  In this case, we observe a junction
    // segment having a very high UMI count in one cell, we observe a very low count in another
    // cell, these two cells share only one chain (allowing for some mutation), and the "weak"
    // cell has at least three chains.
    const ALLOWED_DIFFS: i32 = 10;
    let mut i = 0;
    while i < all.len() {
        let j = next_diff1_5(&all, i as i32) as usize;
        for k1 in i..j {
            let i1 = all[k1].3;
            if all[k1].2 && all[k1].1 >= MIN_RATIO_UMI as u16 && d[i1].jundata.len() >= 2 {
                'couter: for k2 in i..j {
                    let i2 = all[k2].3;
                    if all[k2].2 && all[k2].1 == 1 && d[i2].jundata.len() >= 3 {
                        let mut commons = 0;
                        for c1 in 0..d[i1].jundata.len() {
                            for c2 in 0..d[i2].jundata.len() {
                                let y1 = &&d[i1].jundata[c1].jxn_seq; // first junction segment
                                let y2 = &d[i2].jundata[c2].jxn_seq; // second junction segment
                                if y1 == &y2 {
                                    commons += 1;
                                } else {
                                    let (mut u1, mut u2) = ([0_u8; 80], [0_u8; 80]);
                                    unpack_bases_80(y1, &mut u1);
                                    unpack_bases_80(y2, &mut u2);
                                    // Should have an ndiffs trait and use it here on u1 and u2,
                                    // same as for DnaStrings.  Although, interestingly, the
                                    // structures y1 and y2 are close to DnaStrings and could
                                    // be compared more efficiently as we do for DnaStrings.
                                    let mut dist = 0;
                                    for l in 0..80 {
                                        if u1[l] != u2[l] {
                                            dist += 1;
                                        }
                                    }
                                    if dist <= ALLOWED_DIFFS {
                                        commons += 1;
                                    }
                                }
                                if commons > 1 {
                                    continue 'couter;
                                }
                            }
                        }
                        let x = &d[i2];
                        println!("{} = possible type two plasma cell leakage", x.barcode);
                        kills.push(x.barcode.clone());
                        if let Some(ref mut logger) = filter_logger {
                            logger.log(&FilterLogEntry::cell_calling(
                                x.barcode.clone(),
                                AsmCellFilter::WeakJunction {
                                    contig: format!("{}_contig_{}", x.barcode, all[k2].4 + 1),
                                    param_min_dominant_umis: MIN_RATIO_UMI,
                                    dominant_contig: format!(
                                        "{}_contig_{}",
                                        d[all[k1].3].barcode,
                                        all[k1].4 + 1
                                    ),
                                    dominant_junction_umis: all[k1].1 as usize,
                                },
                            ));
                        }
                    }
                }
            }
        }
        i = j;
    }
}

fn chimeric_filters(
    d: &[BarcodeCellInfo],
    kills: &mut Vec<String>,
    filter_logger: &mut Option<&mut FilterLogger>,
) {
    // Look for chimeric contigs.  For a given cdr3_nt, consider the V segments that appear in
    // contigs.  If one has collective support at least 100 times greater than another, then
    // untrust the weaker contigs.
    const CHIM_RATIO: usize = 100;
    let mut all_chimdata = d
        .iter()
        .flat_map(|bc| bc.chimdata.clone())
        .collect::<Vec<ContigChimeraData>>();
    all_chimdata.sort();
    let mut i = 0;
    while i < all_chimdata.len() {
        let mut j = i + 1;
        let j = loop {
            if j == all_chimdata.len() || all_chimdata[j].cdr3 != all_chimdata[i].cdr3 {
                break j;
            }
            j += 1;
        };
        let mut vu = Vec::<(usize, usize)>::new();
        for k in i..j {
            vu.push((all_chimdata[k].v_ref_id, all_chimdata[k].umi_count));
        }
        let mut uv = Vec::<(usize, usize)>::new();
        let mut r = 0;
        while r < vu.len() {
            let s = next_diff1_2(&vu, r as i32) as usize;
            let mut numi = 0;
            for m in r..s {
                numi += vu[m].1;
            }
            uv.push((numi, vu[r].0));
            r = s;
        }
        reverse_sort(&mut uv);
        let mut bads = Vec::<usize>::new();
        for m in 1..uv.len() {
            if uv[0].0 >= 1 && uv[0].0 >= CHIM_RATIO * uv[m].0 {
                bads.push(uv[m].1);
            }
        }
        bads.sort_unstable();

        for k in i..j {
            if all_chimdata[k].productive {
                let t = all_chimdata[k].v_ref_id;
                if bin_member(&bads, &t) {
                    kills.push(all_chimdata[k].barcode.clone());
                    println!("{} = possible chimera", all_chimdata[k].barcode);
                    if let Some(ref mut logger) = filter_logger {
                        logger.log(&FilterLogEntry::cell_calling(
                            all_chimdata[k].barcode.clone(),
                            AsmCellFilter::ChimericContig {
                                cdr3_nt: DnaString::from_bytes(&all_chimdata[i].cdr3).to_string(),
                                param_chimera_ratio: CHIM_RATIO,
                                contig_v_region_id: t,
                                dominant_v_region_id: uv[0].1,
                                dominant_v_region_umis: uv[0].0,
                            },
                        ));
                    }
                }
            }
        }
        i = j;
    }
}

fn common_clone_filters(
    d: &[BarcodeCellInfo],
    killsc: &mut Vec<String>,
    kills: &mut Vec<String>,
    filter_logger: &mut Option<&mut FilterLogger>,
) {
    // Find two-chain productive pairs and their frequency.
    let mut pairs = Vec::<([u8; 20], [u8; 20])>::new();
    for i in 0..d.len() {
        if d[i].jundata.len() != 2 || !(d[i].paired && d[i].now_a_cell) {
            continue;
        }
        if d[i].jundata[0].jxn_seq <= d[i].jundata[1].jxn_seq {
            pairs.push((d[i].jundata[0].jxn_seq, d[i].jundata[1].jxn_seq));
        } else {
            pairs.push((d[i].jundata[1].jxn_seq, d[i].jundata[0].jxn_seq));
        }
    }
    pairs.sort_unstable();
    let mut pairsu = Vec::<([u8; 20], [u8; 20])>::new();
    let mut pairsf = Vec::<i32>::new();
    let mut j = 0;
    while j < pairs.len() {
        let mut k = j + 1;
        loop {
            if k == pairs.len() || pairs[k] != pairs[j] {
                break;
            }
            k += 1;
        }
        pairsu.push((pairs[j].0, pairs[j].1));
        pairsf.push((k - j) as i32);
        j = k;
    }
    // Find transcripts that appear in these two-chain productive pairs, and the max
    // frequency that was observed for each.  Also track the partner.
    let mut u = Vec::<([u8; 20], i32, [u8; 20])>::new();
    for i in 0..pairsu.len() {
        u.push((pairsu[i].0, pairsf[i], pairsu[i].1));
        u.push((pairsu[i].1, pairsf[i], pairsu[i].0));
    }
    u.sort_unstable();
    let mut to_delete = vec![false; u.len()];
    let mut i = 0;
    while i < u.len() {
        let mut j = i + 1;
        while j < u.len() {
            if u[j].0 != u[i].0 {
                break;
            }
            j += 1;
        }
        for k in i..j - 1 {
            to_delete[k] = true;
        }
        i = j;
    }
    erase_if(&mut u, &to_delete);
    // Make list of the junction segments for the case where there are two or more
    // confident contigs.
    let mut bigs = Vec::<Vec<[u8; 20]>>::new();
    for i in 0..d.len() {
        let x = &d[i];
        let mut jundata = x.jundata.clone();
        let mut to_delete = vec![false; jundata.len()];
        for i in 0..jundata.len() {
            if !jundata[i].high_confidence {
                to_delete[i] = true;
            }
        }
        erase_if(&mut jundata, &to_delete);
        if jundata.len() >= 2 {
            let mut big = Vec::<[u8; 20]>::new();
            for i in 0..jundata.len() {
                big.push(jundata[i].jxn_seq);
            }
            big.sort_unstable();
            bigs.push(big);
        }
    }
    bigs.sort();
    // Identify contigs that should now be labeled low confidence.
    const MAX_KILL: i32 = 3;
    const MIN_RATIO: i32 = 10;
    const MIN_RATIO_BIG: i32 = 50;
    const ALLOWED_DIFFS: i32 = 10;
    for i in 0..d.len() {
        // Remove the low-confidence contigs.  Ignore after that if only one contig is left.

        let x = &d[i];
        let mut jundata = x.jundata.clone();
        let mut to_delete = vec![false; jundata.len()];
        for i in 0..jundata.len() {
            if !jundata[i].high_confidence {
                to_delete[i] = true;
            }
        }
        erase_if(&mut jundata, &to_delete);
        if jundata.len() <= 1 {
            continue;
        }

        // Anything seen rarely and involving a very common clone is deemed
        // dubious, probably a doublet.

        if jundata.len() >= 2 {
            let mut big = Vec::<[u8; 20]>::new();
            for i in 0..jundata.len() {
                big.push(jundata[i].jxn_seq);
            }
            big.sort_unstable();
            let low = lower_bound(&bigs, &big);
            let high = upper_bound(&bigs, &big);
            let mut max_freq = 0_i32;
            let mut best_l = -1_i32;
            let mut best_j = -1_i32;
            for j in 0..jundata.len() {
                let l = bin_position1_3(&u, &jundata[j].jxn_seq);
                if l >= 0 && u[l as usize].1 > max_freq {
                    max_freq = u[l as usize].1;
                    best_l = l;
                    best_j = j as i32;
                }
            }
            let mult = (high - low) as i32;
            if mult <= MAX_KILL && max_freq >= MIN_RATIO_BIG * mult {
                //
                // Try to avoid being tricked by somatic hypermutation.

                let mut protected = false;
                if jundata.len() == 2 {
                    let p1 = &&jundata[(1 - best_j) as usize].jxn_seq;
                    let p2 = &u[best_l as usize].2;
                    let (mut u1, mut u2) = ([0_u8; 80], [0_u8; 80]);
                    unpack_bases_80(p1, &mut u1);
                    unpack_bases_80(p2, &mut u2);
                    let mut dist = 0;
                    for l in 0..80 {
                        if u1[l] != u2[l] {
                            dist += 1;
                        }
                    }
                    if dist <= ALLOWED_DIFFS {
                        protected = true;
                    }
                }
                if !protected {
                    for j in 0..x.jundata.len() {
                        if x.jundata[j].high_confidence {
                            println!("{}: contig {}", x.barcode, j + 1);
                        }
                        killsc.push(format!("{}_contig_{}", x.barcode, j + 1));
                    }
                    kills.push(x.barcode.clone());
                    println!("{}", x.barcode);
                    if let Some(ref mut logger) = *filter_logger {
                        logger.log(&FilterLogEntry::cell_calling(
                            x.barcode.clone(),
                            AsmCellFilter::CommonCloneShadow {
                                multiplicity: mult as usize,
                                max_multiplicity: max_freq as usize,
                                param_max_kill: MAX_KILL as usize,
                                param_min_ratio_big: MIN_RATIO_BIG as usize,
                            },
                        ));
                    }
                    continue;
                }
            }
        }

        // Now assume just two contigs.

        if jundata.len() != 2 {
            continue;
        }
        let min_umis = min(jundata[0].umis, jundata[1].umis);
        let p = if jundata[0].jxn_seq <= jundata[1].jxn_seq {
            (jundata[0].jxn_seq, jundata[1].jxn_seq)
        } else {
            (jundata[1].jxn_seq, jundata[0].jxn_seq)
        };
        let l = bin_position(&pairsu, &p);
        let freq = if l >= 0 { pairsf[l as usize] } else { 0_i32 };
        if freq > MAX_KILL {
            continue;
        }
        let (mut max_alt_freq, mut min_alt_freq) = (0_i32, 1000000000_i32);
        for j in 0..2 {
            let l = bin_position1_3(&u, &jundata[j].jxn_seq);
            if l >= 0 {
                max_alt_freq = max(max_alt_freq, u[l as usize].1);
                min_alt_freq = min(min_alt_freq, u[l as usize].1);
            }
        }

        // The model here is that a single stray UMI from a common clone
        // floats into a GEM.

        if max_alt_freq >= MIN_RATIO * max(1, freq) && min_umis == 1 {
            for j in 0..x.jundata.len() {
                if x.jundata[j].umis <= 1 && x.jundata[j].high_confidence {
                    println!("{}: contig {}", x.barcode, j + 1);
                    killsc.push(format!("{}_contig_{}", x.barcode, j + 1));
                }
            }
            kills.push(x.barcode.clone());
            if let Some(ref mut logger) = filter_logger {
                logger.log(&FilterLogEntry::cell_calling(
                    x.barcode.clone(),
                    AsmCellFilter::CommonCloneShadowSingleUmi {
                        multiplicity: freq as usize,
                        max_multiplicity: max_alt_freq as usize,
                        param_max_kill: MAX_KILL as usize,
                        param_min_ratio_big: MIN_RATIO_BIG as usize,
                    },
                ));
            }
            println!("{}", x.barcode);
        }
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
// ANALYZE BARCODE DATA
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn analyze_barcode_data(
    d: &[BarcodeData],
    refdata: &RefData,
    inner_primers: &[Vec<u8>],
    outer_primers: &[Vec<u8>],
    single_end: bool,
) {
    println!("\nBARCODE SUMMARY STATS\n");

    let dsum: BarcodeDataSum = BarcodeDataSum::sum(d, refdata);

    // Compute mean pairs per umi having more than one pair.

    let mut ucounts_all = Vec::<i32>::new();
    for x in d {
        for i in 0..x.ucounts.len() {
            ucounts_all.push(x.ucounts[i]);
        }
    }
    println!("{:<6.1} = mean pairs per nonsolo umi", mean(&ucounts_all));
    println!("{:<6.1} = N50 pairs per nonsolo umi", n50(&ucounts_all));
    println!("{:<6.1} = N90 pairs per nonsolo umi", n90(&ucounts_all));
    println!(
        "{:<6.1} = length-weighted mean pairs per nonsolo umi",
        len_weighted_mean(&ucounts_all)
    );

    // Compute fraction of reads that survive.

    let mut nreads_all_all = 0_i64;
    for x in d {
        nreads_all_all += x.nreads as i64;
    }
    let mut surviving_pairs = 0_i64;
    for x in d {
        surviving_pairs += x.xucounts.iter().sum::<i32>() as i64;
    }
    let surper = 100_f64 * (surviving_pairs * 2) as f64 / nreads_all_all as f64;
    println!("{surper:<6.1} = percent surviving reads");

    // Compute total number of surviving umis.

    let mut sumi = 0_i64;
    for x in d {
        sumi += x.xucounts.len() as i64;
    }
    println!("{sumi:<6.0} = total surviving umis");

    // Compute TRBV:TRAV read ratio.

    let (mut trav_reads, mut trbv_reads) = (0, 0);
    for x in d {
        trav_reads += x.trav_reads as usize;
        trbv_reads += x.trbv_reads as usize;
    }
    if trav_reads > 0 {
        let ba = trbv_reads as f64 / trav_reads as f64;
        println!("{ba:<6.2} = TRBV/TRAV read ratio");
    }

    let (mut travj_reads, mut trbvj_reads) = (0, 0);
    for x in d {
        travj_reads += x.travj_reads as usize;
        trbvj_reads += x.trbvj_reads as usize;
    }
    if travj_reads > 0 {
        let ba = trbvj_reads as f64 / travj_reads as f64;
        println!("{ba:<6.2} = TRBVJ/TRAVJ read ratio");
    }

    // Report number of contigs.

    let mut ntigs = 0;
    for x in d {
        ntigs += x.ncontigs;
    }
    println!("{ntigs:<6.0} = total contigs");

    // Find abundance of C segments in productive contigs.

    let mut cs = Vec::<String>::new();
    for x in d {
        for id in &x.good_refseqs {
            if refdata.is_c(*id as usize) {
                cs.push(refdata.name[*id as usize].clone());
            }
        }
    }
    cs.sort();
    let mut i = 0;
    while i < cs.len() {
        let j = next_diff(&cs, i);
        println!(
            "{:<6.0} = total occurrences of {} in productive contigs",
            j - i,
            cs[i]
        );
        i = j;
    }

    // Give up if no barcodes.

    if d.is_empty() {
        return;
    }

    // Analyze primer hits.

    let mut json = Vec::<u8>::new();
    let mut log = Vec::<u8>::new();
    analyze_primer_hits(
        &dsum,
        refdata,
        inner_primers,
        outer_primers,
        single_end,
        &mut json,
        &mut log,
    );
    print!("{}", strme(&log));

    // Done.

    println!();
}

fn analyze_primer_hits(
    dsum: &BarcodeDataSum,
    refdata: &RefData,
    inner_primers: &[Vec<u8>],
    outer_primers: &[Vec<u8>],
    single_end: bool,
    json: &mut Vec<u8>,
    log: &mut Vec<u8>,
) {
    // Analyze primer hits.  The math here is tricky because when we sample, we don't get a true
    // random sample, and so we have to normalize.
    // ◼ Ugly inner/outer code duplication.

    // ◼ The naming here for good is very confusing.
    let mut inner_hit_total = vec![0 as f64; inner_primers.len()];
    let mut inner_hit_good_total = vec![0 as f64; inner_primers.len()];
    let mut outer_hit_total = vec![0 as f64; outer_primers.len()];
    let mut outer_hit_good_total = vec![0 as f64; outer_primers.len()];
    let mut all_pairs = dsum.nreads as f64;
    if !single_end {
        all_pairs /= 2.0;
    }
    let inner_hit_good_contigs = &dsum.inner_hit_good_contigs_total;
    let outer_hit_good_contigs = &dsum.outer_hit_good_contigs_total;
    inner_hit_total.clear();
    for i in 0..dsum.inner_hit_total.len() {
        inner_hit_total.push(dsum.inner_hit_total[i].to_f64().unwrap());
    }
    inner_hit_good_total.clear();
    for i in 0..dsum.inner_hit_good_total.len() {
        inner_hit_good_total.push(dsum.inner_hit_good_total[i].to_f64().unwrap());
    }
    outer_hit_total.clear();
    for i in 0..dsum.outer_hit_total.len() {
        outer_hit_total.push(dsum.outer_hit_total[i].to_f64().unwrap());
    }
    outer_hit_good_total.clear();
    for i in 0..dsum.outer_hit_good_total.len() {
        outer_hit_good_total.push(dsum.outer_hit_good_total[i].to_f64().unwrap());
    }
    if !single_end {
        inner_hit_total.iter_mut().for_each(|x| *x /= 2.0);
        inner_hit_good_total.iter_mut().for_each(|x| *x /= 2.0);
        outer_hit_total.iter_mut().for_each(|x| *x /= 2.0);
        outer_hit_good_total.iter_mut().for_each(|x| *x /= 2.0);
    }
    let primed_pairs = inner_hit_total.iter().sum::<f64>() + outer_hit_total.iter().sum::<f64>();
    for i in 0..inner_hit_total.len() {
        if inner_hit_total[i] > 0.0 {
            inner_hit_good_total[i] /= inner_hit_total[i];
            inner_hit_good_total[i] = 1.0 - inner_hit_good_total[i];
        }
        inner_hit_total[i] /= all_pairs;
    }
    for i in 0..outer_hit_total.len() {
        if outer_hit_total[i] > 0.0 {
            outer_hit_good_total[i] /= outer_hit_total[i];
            outer_hit_good_total[i] = 1.0 - outer_hit_good_total[i];
        }
        outer_hit_total[i] /= all_pairs;
    }
    let total_primed = inner_hit_total.iter().sum::<f64>() + outer_hit_total.iter().sum::<f64>();
    inner_hit_total.iter_mut().for_each(|x| *x /= total_primed);
    outer_hit_total.iter_mut().for_each(|x| *x /= total_primed);
    let mut locs_all = Vec::<Vec<String>>::new();
    for i in 0..inner_primers.len() {
        let mut p = inner_primers[i].clone();
        reverse_complement(&mut p);
        let mut locs = Vec::<String>::new();
        for j in 0..refdata.refs.len() {
            if refdata.is_c(j) {
                let c = refdata.refs[j].to_ascii_vec();
                for l in 0..c.len() {
                    if contains_at(&c, &p, l) && l + p.len() >= PRIMER_EXT_LEN {
                        locs.push(refdata.name[j].clone());
                    }
                }
            }
            unique_sort(&mut locs);
        }
        locs_all.push(locs);
    }
    let mut total_off = 0.0;
    for i in 0..inner_primers.len() {
        total_off += inner_hit_total[i] * inner_hit_good_total[i];
    }
    for i in 0..outer_primers.len() {
        total_off += outer_hit_total[i] * outer_hit_good_total[i];
    }
    fwriteln!(
        log,
        "\n▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓"
    );
    fwriteln!(log, "PRIMER ANALYSIS");
    fwriteln!(
        log,
        "▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓\n"
    );
    fwriteln!(log, "OVERALL PRIMER STATS");
    write_json_metric_f64(
        "frac_of_read_pairs_containing_a_primer",
        primed_pairs / all_pairs,
        json,
    );
    fwriteln!(
        log,
        "fraction of read pairs containing primer = {:.1}%",
        100.0 * primed_pairs / all_pairs
    );
    write_json_metric_f64(
        "frac_of_priming_events_that_are_off_target",
        total_off,
        json,
    );
    fwriteln!(
        log,
        "fraction of priming events that are off target = {:.1}%",
        100.0 * total_off
    );
    let mut rows = Vec::<Vec<String>>::new();
    fwriteln!(
        log,
        "--------------------------------------------------------------------"
    );
    fwriteln!(log, "INNER PRIMER STATS");
    fwriteln!(log, "#    = inner primer number");
    fwriteln!(log, "len  = length of this primer");
    fwriteln!(log, "loc  = constant regions that primer aligns to");
    fwriteln!(
        log,
        "prod = number of productive contigs containing rc of primer"
    );
    fwriteln!(
        log,
        "freq = fraction of priming events attributable to this primer"
    );
    fwriteln!(
        log,
        "off  = fraction of priming events by this primer that are off target"
    );
    fwriteln!(
        log,
        "mul  = fraction of priming events that this primer makes off target"
    );
    fwriteln!(
        log,
        "--------------------------------------------------------------------"
    );
    rows.push(vec![
        "#".to_string(),
        "len".to_string(),
        "loc".to_string(),
        "prod".to_string(),
        "freq".to_string(),
        "off".to_string(),
        "mul".to_string(),
    ]);
    for i in 0..inner_primers.len() {
        let mut row = Vec::<String>::new();
        row.push(format!("{}", i + 1));
        row.push(format!("{}", inner_primers[i].len()));
        let metric = format!("inner_primer_{}", i + 1);
        write_json_metric_str(&metric, strme(&inner_primers[i]), json);
        let locs = format!("{}", locs_all[i].iter().format("+"));
        row.push(locs.clone());
        let metric = format!("binding_sites_of_inner_primer_{}", i + 1);
        write_json_metric_str(&metric, &locs, json);
        row.push(format!("{}", inner_hit_good_contigs[i]));
        let metric = format!("productive_contigs_containing_inner_primer_{}", i + 1);
        write_json_metric_f64(&metric, inner_hit_good_contigs[i] as f64, json);
        row.push(format!("{:.1}%", 100.0 * inner_hit_total[i]));
        let metric = format!("frac_of_priming_events_from_inner_primer_{}", i + 1);
        write_json_metric_f64(&metric, inner_hit_total[i], json);
        row.push(format!("{:.1}%", 100.0 * inner_hit_good_total[i]));
        let metric = format!(
            "frac_of_priming_events_from_inner_primer_{}_that_are_off_target",
            i + 1
        );
        write_json_metric_f64(&metric, inner_hit_good_total[i], json);
        row.push(format!(
            "{:.1}%",
            100.0 * inner_hit_total[i] * inner_hit_good_total[i]
        ));
        let metric = format!(
            "frac_of_priming_events_that_inner_primer_{}_makes_off_target",
            i + 1
        );
        write_json_metric_f64(&metric, inner_hit_total[i] * inner_hit_good_total[i], json);
        rows.push(row);
    }
    print_tabular(log, &rows, 2, Some(b"lrlrrrr".to_vec()));
    fwriteln!(
        log,
        "--------------------------------------------------------------------"
    );
    let mut locs_all = Vec::<Vec<String>>::new();
    for i in 0..outer_primers.len() {
        let mut p = outer_primers[i].clone();
        reverse_complement(&mut p);
        let mut locs = Vec::<String>::new();
        for j in 0..refdata.refs.len() {
            if refdata.is_c(j) {
                let c = refdata.refs[j].to_ascii_vec();
                for l in 0..c.len() {
                    if contains_at(&c, &p, l) && l + p.len() >= PRIMER_EXT_LEN {
                        locs.push(refdata.name[j].clone());
                    }
                }
            }
            unique_sort(&mut locs);
        }
        locs_all.push(locs);
    }
    let mut rows = Vec::<Vec<String>>::new();
    fwriteln!(log, "OUTER PRIMER STATS");
    fwriteln!(
        log,
        "--------------------------------------------------------------------"
    );
    rows.push(vec![
        "#".to_string(),
        "len".to_string(),
        "loc".to_string(),
        "prod".to_string(),
        "freq".to_string(),
        "off".to_string(),
        "mul".to_string(),
    ]);
    for i in 0..outer_primers.len() {
        let mut row = Vec::<String>::new();
        row.push(format!("{}", i + 1));
        row.push(format!("{}", outer_primers[i].len()));
        let metric = format!("outer_primer_{}", i + 1);
        write_json_metric_str(&metric, strme(&outer_primers[i]), json);
        let locs = format!("{}", locs_all[i].iter().format("+"));
        row.push(locs.clone());
        let metric = format!("binding_sites_of_outer_primer_{}", i + 1);
        write_json_metric_str(&metric, &locs, json);
        row.push(format!("{}", outer_hit_good_contigs[i]));
        let metric = format!("productive_contigs_containing_outer_primer_{}", i + 1);
        write_json_metric_f64(&metric, outer_hit_good_contigs[i] as f64, json);
        row.push(format!("{:.1}%", 100.0 * outer_hit_total[i]));
        let metric = format!("frac_of_priming_events_from_outer_primer_{}", i + 1);
        write_json_metric_f64(&metric, outer_hit_total[i], json);
        row.push(format!("{:.1}%", 100.0 * outer_hit_good_total[i]));
        let metric = format!(
            "frac_of_priming_events_from_outer_primer_{}_that_are_off_target",
            i + 1
        );
        write_json_metric_f64(&metric, outer_hit_good_total[i], json);
        row.push(format!(
            "{:.1}%",
            100.0 * outer_hit_total[i] * outer_hit_good_total[i]
        ));
        let metric = format!(
            "frac_of_priming_events_that_outer_primer_{}_makes_off_target",
            i + 1
        );
        write_json_metric_f64(&metric, outer_hit_total[i] * outer_hit_good_total[i], json);
        rows.push(row);
    }
    print_tabular(log, &rows, 2, Some(b"lrlrrrr".to_vec()));
    fwriteln!(
        log,
        "--------------------------------------------------------------------"
    );
    fwriteln!(log, "METRICS USED ABOVE");
    fwriteln!(log, "frac_of_read_pairs_containing_a_primer");
    fwriteln!(log, "frac_of_priming_events_that_are_off_target");
    fwriteln!(log, "inner_primer_<n>");
    fwriteln!(log, "binding_sites_of_inner_primer_<n>");
    fwriteln!(log, "productive_contigs_containing_inner_primer_<n>");
    fwriteln!(log, "frac_of_priming_events_from_inner_primer_<n>");
    fwriteln!(
        log,
        "frac_of_priming_events_from_inner_primer_<n>_that_are_off_target"
    );
    fwriteln!(
        log,
        "frac_of_priming_events_that_inner_primer_<n>_makes_off_target"
    );
    fwriteln!(
        log,
        "[and matching outer primer stats for each inner primer stat]"
    );
}
