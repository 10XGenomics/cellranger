// This captures and analyzes data about barcodes.

#![allow(clippy::many_single_char_names)]
// TODO: fix these.
#![allow(clippy::needless_range_loop)]

use crate::constants::{CHAIN_TYPES, CHAIN_TYPESX, PRIMER_EXT_LEN};
use itertools::Itertools;
use metric::{JsonReporter, PercentMetric, SimpleHistogram, TxHashMap};
use ordered_float::NotNan;
use serde::{Deserialize, Serialize};
use stats_utils::{len_weighted_mean, mean, n50, n90};
use std::cmp::max;
use std::fs::File;
use std::io::{BufWriter, Write};
use string_utils::strme;
use tenkit2::pack_dna::reverse_complement;
use vdj_ann::refx::RefData;
use vector_utils::{contains_at, next_diff, unique_sort};

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
// BARCODE DATA STRUCTURES
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

#[derive(PartialEq, PartialOrd, Clone, Serialize, Deserialize)]
pub struct BarcodeData {
    // the barcode
    pub barcode: String,

    // number of reads in the barcode
    pub nreads: i32,
    // Number of subsampled reads used for assembly in barcode
    // if nreads < vdj_max_reads_per_barcode, nreads = nreads_used_for_assembly
    pub nreads_used_for_assembly: i32,

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
    pub outer_hit: Vec<u16>, // one entry for each outer primer
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
            nreads_used_for_assembly: 0_i32,
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

    /// Function returns 1 / subsample fraction for primer search.
    /// Used to mutliply number of hits in sample to estimate the total number of hits.
    fn primer_inv_subsample_frac(&self) -> f64 {
        self.nreads as f64 / self.primer_sample_count as f64
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
    pub is_cell_and_productive: bool,

    // the barcode
    pub barcode: String,
}

impl ContigChimeraData {
    pub fn update_productive(mut self, is_cell: bool) -> Self {
        self.is_cell_and_productive = self.is_cell_and_productive && is_cell;
        self
    }
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

// The number of read pairs and Vdj chain-assignment of a single UMI
#[derive(Debug, Serialize, Deserialize, Hash, PartialEq, Clone, PartialOrd, Ord, Eq)]
pub struct VdjUmi {
    chain_type: String,
    npairs: usize,
}

#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct BarcodeDataSum {
    pub nreads: i64,                       // number of reads
    pub nreads_used_for_assembly: i64,     // number of reads used for assembly
    pub nreads_umi_corrected: i64,         // number of reads whose umi was corrected
    pub ncontigs: i64,                     // number of contigs
    pub xucount_sum: i64,                  // sum over xucounts
    pub chain_sample: Vec<i64>,            // chain sample
    pub adj_chain_sample: Vec<f64>, // adjusted chain_sample (summation done after scaling up based on frac)
    pub umi_dist: SimpleHistogram<VdjUmi>, // (chain-type,npairs,multiplicity)
    pub cs_good: SimpleHistogram<String>, // for each constant region, count in good
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
            nreads_used_for_assembly: 0_i64,
            nreads_umi_corrected: 0_i64,
            ncontigs: 0_i64,
            xucount_sum: 0_i64,
            chain_sample: vec![0_i64; 15],
            adj_chain_sample: vec![0_f64; 15],
            umi_dist: SimpleHistogram::<VdjUmi>::default(),
            cs_good: SimpleHistogram::<String>::default(),
            inner_hit_total: Vec::<NotNan<f64>>::new(),
            inner_hit_good_total: Vec::<NotNan<f64>>::new(),
            outer_hit_total: Vec::<NotNan<f64>>::new(),
            outer_hit_good_total: Vec::<NotNan<f64>>::new(),
            inner_hit_good_contigs_total: Vec::<usize>::new(),
            outer_hit_good_contigs_total: Vec::<usize>::new(),
        }
    }
    pub fn sum(barcode_data_list: &[BarcodeData], refdata: &RefData) -> BarcodeDataSum {
        BarcodeDataSum {
            nreads: barcode_data_list
                .iter()
                .map(|bc_data| bc_data.nreads as i64)
                .sum(),
            nreads_used_for_assembly: barcode_data_list
                .iter()
                .map(|bc_data| bc_data.nreads_used_for_assembly as i64)
                .sum(),
            nreads_umi_corrected: barcode_data_list
                .iter()
                .map(|bc_data| bc_data.nreads_umi_corrected as i64)
                .sum(),
            ncontigs: barcode_data_list
                .iter()
                .map(|bc_data| bc_data.ncontigs as i64)
                .sum(),
            xucount_sum: barcode_data_list
                .iter()
                .map(|bc_data| bc_data.xucounts.iter().sum::<i32>() as i64)
                .sum(),
            chain_sample: {
                let mut s = vec![0_i64; 15];
                for bc_data in barcode_data_list {
                    for i in 0..15 {
                        s[i] += bc_data.chain_sample[i] as i64;
                    }
                }
                s
            },
            // When summing barcode objects, adjust each chain_sample value based on frac
            adj_chain_sample: {
                let mut s = vec![0_f64; 15];
                for bc_data in barcode_data_list {
                    for i in 0..15 {
                        if bc_data.frac != 0.0 {
                            s[i] += bc_data.chain_sample[i] as f64 / bc_data.frac;
                        } else {
                            s[i] += 0.0;
                        }
                    }
                }
                s
            },
            umi_dist: {
                let mut hist = SimpleHistogram::<VdjUmi>::default();
                for bc_data in barcode_data_list {
                    for (npairs, chain_type, _) in &bc_data.umi_info {
                        hist.observe(&VdjUmi {
                            chain_type: CHAIN_TYPES[*chain_type as usize].to_owned(),
                            npairs: npairs.to_owned() as usize,
                        });
                    }
                }

                hist
            },
            cs_good: SimpleHistogram::<String>::from_iter_owned(barcode_data_list.iter().flat_map(
                |bc_data| {
                    bc_data
                        .good_refseqs
                        .iter()
                        // filter out any non constant segment refdata indexes
                        .filter(|&&refseq_idx| {
                            refdata
                                .cs
                                .iter()
                                .any(|constant_idx| constant_idx == &(refseq_idx as usize))
                        })
                        // map refdata index to name
                        .map(|&refseq_idx| refdata.name[refseq_idx as usize].clone())
                },
            )),
            inner_hit_total: {
                // in single-end case, remember to divide by two later, and below
                if barcode_data_list.is_empty() {
                    Vec::<NotNan<f64>>::new()
                } else {
                    let mut s = vec![0.0; barcode_data_list[0].inner_hit.len()];
                    for bc_data in barcode_data_list {
                        for i in 0..bc_data.inner_hit.len() {
                            s[i] +=
                                bc_data.inner_hit[i] as f64 * bc_data.primer_inv_subsample_frac();
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
                if barcode_data_list.is_empty() {
                    Vec::<NotNan<f64>>::new()
                } else {
                    let mut sum_vec = vec![0.0; barcode_data_list[0].inner_hit_good.len()];
                    for bc_data in barcode_data_list {
                        for i in 0..bc_data.inner_hit_good.len() {
                            sum_vec[i] += bc_data.inner_hit_good[i] as f64
                                * bc_data.primer_inv_subsample_frac();
                        }
                    }
                    let mut t = Vec::<NotNan<f64>>::new();
                    for i in 0..sum_vec.len() {
                        t.push(NotNan::new(sum_vec[i]).unwrap());
                    }
                    t
                }
            },
            outer_hit_total: {
                if barcode_data_list.is_empty() {
                    Vec::<NotNan<f64>>::new()
                } else {
                    let mut s = vec![0.0; barcode_data_list[0].outer_hit.len()];
                    for bc_data in barcode_data_list {
                        for i in 0..bc_data.outer_hit.len() {
                            s[i] +=
                                bc_data.outer_hit[i] as f64 * bc_data.primer_inv_subsample_frac();
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
                if barcode_data_list.is_empty() {
                    Vec::<NotNan<f64>>::new()
                } else {
                    let mut s = vec![0.0; barcode_data_list[0].outer_hit_good.len()];
                    for bc_data in barcode_data_list {
                        for i in 0..bc_data.outer_hit_good.len() {
                            s[i] += bc_data.outer_hit_good[i] as f64
                                * bc_data.primer_inv_subsample_frac();
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
                if barcode_data_list.is_empty() {
                    Vec::<usize>::new()
                } else {
                    let mut s = vec![0; barcode_data_list[0].inner_hit_good_contigs.len()];
                    for bc_data in barcode_data_list {
                        for i in 0..bc_data.inner_hit_good_contigs.len() {
                            s[i] += bc_data.inner_hit_good_contigs[i] as usize;
                        }
                    }
                    s
                }
            },
            outer_hit_good_contigs_total: {
                if barcode_data_list.is_empty() {
                    Vec::<usize>::new()
                } else {
                    let mut s = vec![0; barcode_data_list[0].outer_hit_good_contigs.len()];
                    for bc_data in barcode_data_list {
                        for i in 0..bc_data.outer_hit_good_contigs.len() {
                            s[i] += bc_data.outer_hit_good_contigs[i] as usize;
                        }
                    }
                    s
                }
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

pub fn write_json_metric_f64(metric_name: &str, value: f64, json: &mut Vec<u8>) {
    if value.is_nan() {
        fwriteln!(json, "    \"{}\": null,", metric_name);
    } else {
        fwriteln!(json, "    \"{}\": {},", metric_name, value);
    }
}

#[allow(clippy::too_many_arguments)]
pub fn metrics_json(
    x: &BarcodeDataSum,
    single_end: bool,
    refdata: &RefData,
    inner_primers: &[Vec<u8>],
    outer_primers: &[Vec<u8>],
    npairs: Option<i64>,
    report: &mut BufWriter<File>,
) -> JsonReporter {
    let mut json = JsonReporter::default();

    // Analyze primer hits.
    let mut log = Vec::<u8>::new();
    analyze_primer_hits(
        x,
        refdata,
        inner_primers,
        outer_primers,
        single_end,
        &mut json,
        &mut log,
    );

    // Write metrics productive_contigs_with_x, where x is a constant region name.
    // do this for all ref constant segment names i.e add a "0" count metric for any not observed
    for name in refdata
        .cs
        .iter()
        .map(|&const_seg_idx| refdata.name[const_seg_idx].clone())
        .unique()
    {
        json.insert(
            format!("productive_contigs_with_{}", &name),
            x.cs_good.get(&name),
        );
    }

    // Compute "vdj_corrected_umi_frac", which is the fraction of reads whose UMI
    // was corrected.

    let corrected_umi_frac =
        PercentMetric::from_parts(x.nreads_umi_corrected, x.nreads_used_for_assembly)
            .fraction()
            .unwrap_or(0.0);
    json.insert("vdj_corrected_umi_frac", corrected_umi_frac);

    // Generate {IGH,IGK,IGL,TRA,TRB,TRD,TRG,multi}
    // _vdj_recombinome_readpairs_per_umi_distribution metrics.  Note that we put
    // the distribution in order by numerical count, as opposed to ordering
    // alphabetically, which is how it was done before.

    let mut chains = CHAIN_TYPES.to_vec();
    chains.push("multi");
    let mut umi_dist = x.umi_dist.clone();
    for (vdj_umi, value) in x.umi_dist.distribution() {
        if &vdj_umi.chain_type != "None" {
            umi_dist.observe_by_owned(
                VdjUmi {
                    chain_type: "multi".to_owned(),
                    npairs: vdj_umi.npairs,
                },
                value.count(),
            );
        }
    }
    for chain in chains {
        if chain == "None" {
            continue;
        }

        let metric_key = format!("{chain}_vdj_recombinome_readpairs_per_umi_distribution");
        let mut metric_val = TxHashMap::<usize, usize>::default();
        for (vdj_umi, value) in umi_dist.distribution() {
            if vdj_umi.chain_type == chain {
                metric_val.insert(vdj_umi.npairs, value.count().try_into().unwrap());
            }
        }
        json.insert(metric_key, serde_json::to_value(metric_val).unwrap());
    }

    // Generate {IGH,IGK,IGL,TRA,TRB,TRD,TRG,multi}
    // _vdj_recombinome_antisense_reads_frac metrics.

    let (mut chain_total_fw, mut chain_total_rc) = (0, 0);
    for i in 0..CHAIN_TYPESX.len() {
        let (fw, rc) = (x.chain_sample[i] as usize, x.chain_sample[i + 7] as usize);
        let metric_key = format!("{}_vdj_recombinome_antisense_reads_frac", CHAIN_TYPESX[i]);
        let metric_val = PercentMetric::from_parts(rc as i64, (fw + rc).try_into().unwrap());
        json.insert(metric_key, metric_val.fraction());
        chain_total_fw += fw;
        chain_total_rc += rc;
    }
    let multi_recombinome_antisense_total = PercentMetric::from_parts(
        chain_total_rc as i64,
        (chain_total_fw + chain_total_rc).try_into().unwrap(),
    );
    json.insert(
        "multi_vdj_recombinome_antisense_reads_frac",
        &multi_recombinome_antisense_total,
    );

    // Generate {IGH,IGK,IGL,TRA,TRB,TRD,TRG,multi}
    // _vdj_recombinome_mapped_reads_frac metrics.

    // Using adjusted chain_sample to account for capping number of reads used in assembly
    let adj_chain_total = x.adj_chain_sample.iter().sum::<f64>();
    for i in 0..CHAIN_TYPESX.len() {
        let metric_key = format!("{}_vdj_recombinome_mapped_reads_frac", CHAIN_TYPESX[i]);
        let metric_val = PercentMetric::from_parts(
            (x.adj_chain_sample[i] + x.adj_chain_sample[i + 7]) as i64,
            adj_chain_total as i64,
        );
        json.insert(metric_key, &metric_val);
    }
    let recombinome_mapped_reads = PercentMetric::from_parts(
        (adj_chain_total - x.adj_chain_sample[14]) as i64,
        adj_chain_total as i64,
    );
    json.insert(
        "multi_vdj_recombinome_mapped_reads_frac",
        &recombinome_mapped_reads,
    );

    // Compute vdj_sequencing_efficiency.  For this, define the number of
    // assemblable read pairs to be the total number of read pairs appearing in
    // xucounts.  Then divide by the total number of read pairs.
    let sequencing_efficiency = if let Some(npairs) = npairs {
        PercentMetric::from_parts(x.xucount_sum, npairs).fraction()
    } else {
        None
    };
    json.insert("vdj_sequencing_efficiency", sequencing_efficiency);

    fwrite!(report, "{}", strme(&log));
    json
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

    let mut json = JsonReporter::default();
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
    json: &mut JsonReporter,
    log: &mut Vec<u8>,
) {
    let all_pairs = match single_end {
        true => dsum.nreads as f64,
        false => PercentMetric::from_parts(dsum.nreads, 2)
            .fraction()
            .unwrap_or(0.0),
    };
    if all_pairs == 0.0 {
        return;
    }
    // Analyze primer hits.  The math here is tricky because when we sample, we don't get a true
    // random sample, and so we have to normalize.
    // ◼ Ugly inner/outer code duplication.

    let mut inner_hit_total: Vec<f64> = dsum
        .inner_hit_total
        .iter()
        .map(|n| n.into_inner())
        .collect();
    let mut inner_hit_good_total: Vec<f64> = dsum
        .inner_hit_good_total
        .iter()
        .map(|n| n.into_inner())
        .collect();
    let mut outer_hit_total: Vec<f64> = dsum
        .outer_hit_total
        .iter()
        .map(|n| n.into_inner())
        .collect();
    let mut outer_hit_good_total: Vec<f64> = dsum
        .outer_hit_good_total
        .iter()
        .map(|n| n.into_inner())
        .collect();

    if !single_end {
        inner_hit_total.iter_mut().for_each(|x| *x /= 2.0);
        inner_hit_good_total.iter_mut().for_each(|x| *x /= 2.0);
        outer_hit_total.iter_mut().for_each(|x| *x /= 2.0);
        outer_hit_good_total.iter_mut().for_each(|x| *x /= 2.0);
    }

    let inner_offtarget_rate: Vec<f64> = inner_hit_good_total
        .iter()
        .zip(inner_hit_total.iter())
        .map(|(n, d)| if d > &0.0 { 1.0 - (n / d) } else { *d })
        .collect();
    let outer_offtarget_rate: Vec<f64> = outer_hit_good_total
        .iter()
        .zip(outer_hit_total.iter())
        .map(|(n, d)| if d > &0.0 { 1.0 - (n / d) } else { *d })
        .collect();

    let primed_pairs = inner_hit_total.iter().sum::<f64>() + outer_hit_total.iter().sum::<f64>();
    let inner_hit_primer_frac: Vec<f64> =
        inner_hit_total.iter().map(|i| i / primed_pairs).collect();
    let outer_hit_primer_frac: Vec<f64> =
        outer_hit_total.iter().map(|i| i / primed_pairs).collect();

    let inner_hit_good_contigs = &dsum.inner_hit_good_contigs_total;
    let outer_hit_good_contigs = &dsum.outer_hit_good_contigs_total;

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
    let total_offtarget_frac = inner_hit_primer_frac
        .iter()
        .zip(inner_offtarget_rate.iter())
        .map(|(t, o)| t * o)
        .sum::<f64>()
        + outer_hit_primer_frac
            .iter()
            .zip(outer_offtarget_rate.iter())
            .map(|(t, o)| t * o)
            .sum::<f64>();
    let total_primer_frac: f64 = primed_pairs / all_pairs;

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
    json.insert("frac_of_read_pairs_containing_a_primer", total_primer_frac);
    fwriteln!(
        log,
        "fraction of read pairs containing primer = {:.1}%",
        100.0 * total_primer_frac,
    );
    json.insert(
        "frac_of_priming_events_that_are_off_target",
        total_offtarget_frac,
    );
    fwriteln!(
        log,
        "fraction of priming events that are off target = {:.1}%",
        100.0 * total_offtarget_frac,
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
        json.insert(metric, strme(&inner_primers[i]));
        let locs = format!("{}", locs_all[i].iter().format("+"));
        row.push(locs.clone());
        let metric = format!("binding_sites_of_inner_primer_{}", i + 1);
        json.insert(metric, locs);
        row.push(format!("{}", inner_hit_good_contigs[i]));
        let metric = format!("productive_contigs_containing_inner_primer_{}", i + 1);
        json.insert(metric, inner_hit_good_contigs[i]);
        row.push(format!("{:.1}%", 100.0 * inner_hit_primer_frac[i]));
        let metric = format!("frac_of_priming_events_from_inner_primer_{}", i + 1);
        json.insert(metric, inner_hit_primer_frac[i]);
        row.push(format!("{:.1}%", 100.0 * inner_offtarget_rate[i]));
        let metric = format!(
            "frac_of_priming_events_from_inner_primer_{}_that_are_off_target",
            i + 1
        );
        json.insert(metric, inner_offtarget_rate[i]);
        row.push(format!(
            "{:.1}%",
            100.0 * inner_hit_primer_frac[i] * inner_offtarget_rate[i]
        ));
        let metric = format!(
            "frac_of_priming_events_that_inner_primer_{}_makes_off_target",
            i + 1
        );
        json.insert(metric, inner_hit_primer_frac[i] * inner_offtarget_rate[i]);
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
        json.insert(metric, strme(&outer_primers[i]));
        let locs = format!("{}", locs_all[i].iter().format("+"));
        row.push(locs.clone());
        let metric = format!("binding_sites_of_outer_primer_{}", i + 1);
        json.insert(metric, locs);
        row.push(format!("{}", outer_hit_good_contigs[i]));
        let metric = format!("productive_contigs_containing_outer_primer_{}", i + 1);
        json.insert(metric, outer_hit_good_contigs[i]);
        row.push(format!("{:.1}%", 100.0 * outer_hit_primer_frac[i]));
        let metric = format!("frac_of_priming_events_from_outer_primer_{}", i + 1);
        json.insert(metric, outer_hit_primer_frac[i]);
        row.push(format!("{:.1}%", 100.0 * outer_offtarget_rate[i]));
        let metric = format!(
            "frac_of_priming_events_from_outer_primer_{}_that_are_off_target",
            i + 1
        );
        json.insert(metric, outer_offtarget_rate[i]);
        row.push(format!(
            "{:.1}%",
            100.0 * outer_hit_primer_frac[i] * outer_offtarget_rate[i]
        ));
        let metric = format!(
            "frac_of_priming_events_that_outer_primer_{}_makes_off_target",
            i + 1
        );
        json.insert(metric, outer_hit_primer_frac[i] * outer_offtarget_rate[i]);
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

/// Print out a matrix, with left-justified entries, and given separation between
/// columns.  (Justification may be changed by supplying an optional argument
/// consisting of a string of l's and r's.)
pub fn print_tabular(
    log: &mut Vec<u8>,
    rows: &[Vec<String>],
    sep: usize,
    justify: Option<Vec<u8>>,
) {
    let just = justify.unwrap_or_default();
    let nrows = rows.len();
    let mut ncols = 0;
    for row in &rows[0..nrows] {
        ncols = max(ncols, row.len());
    }
    let mut maxcol = vec![0; ncols];
    for row in rows {
        for (j, item) in row.iter().enumerate() {
            maxcol[j] = max(maxcol[j], item.chars().count());
        }
    }
    for row in rows {
        for (j, x) in row.iter().enumerate() {
            if j < just.len() && just[j] == b'r' {
                log.append(&mut vec![b' '; maxcol[j] - x.chars().count()]);
                log.append(&mut x.as_bytes().to_vec());
                if j < row.len() - 1 {
                    log.append(&mut vec![b' '; sep]);
                }
            } else {
                log.append(&mut x.as_bytes().to_vec());
                if j < row.len() - 1 {
                    log.append(&mut vec![b' '; maxcol[j] - x.chars().count() + sep]);
                }
            }
        }
        log.push(b'\n');
    }
}
