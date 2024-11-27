// TODO: fix these.
#![allow(clippy::needless_range_loop)]

use crate::filter_barcode_level::BarcodeCellInfo;
use crate::filter_log::{
    AsmCellFilter, FilterLogEntry, FilterLogger, FilterSwitch, IndelErrorMode,
};
use barcode::{Barcode, BcSeq};
use debruijn::dna_string::DnaString;
use fastq_set::sseq::{HammingIterOpt, InsertionIterOpt};
use itertools::Itertools;
use std::cmp::{max, min};
use std::collections::{HashMap, HashSet};
use tenkit2::pack_dna::unpack_bases_80;
use vdj_asm_utils::barcode_data::ContigChimeraData;
use vdj_asm_utils::exact_clonotyping::{BarcodeContigUMI, ProductiveContig};
use vector_utils::{
    bin_member, bin_position, bin_position1_3, erase_if, lower_bound, next_diff, next_diff1_2,
    next_diff1_5, reverse_sort, unique_sort, upper_bound,
};

pub fn whitelist_indel_filter(
    exact_contigs: HashMap<ProductiveContig, Vec<BarcodeContigUMI>>,
    mut filter_logger: Option<&mut FilterLogger>,
) -> Vec<String> {
    const GB_UMI_MULT: usize = 10;

    struct Indels {
        insertions: Vec<BcSeq>,
        insertions_1hd: Vec<BcSeq>,
        deletions: Vec<BcSeq>,
        deletions_1hd: Vec<BcSeq>,
    }

    impl Indels {
        fn generate(source_bc: &BcSeq) -> Self {
            let insertions: Vec<BcSeq> = source_bc
                .one_insertion_iter(InsertionIterOpt::ExcludeNBase)
                .map(|bc| BcSeq::from_bytes(&bc.as_bytes()[..bc.len() - 1]))
                .collect();
            let insertions_1hd: Vec<BcSeq> = source_bc
                .one_insertion_iter(InsertionIterOpt::ExcludeNBase)
                .flat_map(|i| i.one_hamming_iter(HammingIterOpt::SkipNBase))
                .map(|bc| BcSeq::from_bytes(&bc.as_bytes()[..bc.len() - 1]))
                .collect();
            let deletions: Vec<BcSeq> = source_bc
                .one_deletion_iter()
                .map(|bc| BcSeq::from_bytes(bc.as_bytes()))
                .collect();
            let deletions_1hd: Vec<BcSeq> = source_bc
                .one_deletion_iter()
                .flat_map(|i| i.one_hamming_iter(HammingIterOpt::SkipNBase))
                .map(|bc| BcSeq::from_bytes(bc.as_bytes()))
                .collect();
            Indels {
                insertions,
                insertions_1hd,
                deletions,
                deletions_1hd,
            }
        }
    }

    fn find_source_bcs(exact_contig: &[BarcodeContigUMI]) -> Vec<&BarcodeContigUMI> {
        let mut source_bcs: Vec<&BarcodeContigUMI> = vec![];
        for this_bc in exact_contig {
            if exact_contig
                .iter()
                .any(|bc| this_bc.umi >= bc.umi * GB_UMI_MULT)
            {
                source_bcs.push(this_bc.to_owned());
            }
        }
        source_bcs
    }

    struct GBIndels {
        sink_barcode: Barcode,
        sink_umis: usize,
        source_barcode: Barcode,
        source_umis: usize,
        error_mode: IndelErrorMode,
    }

    fn find_matching_prefix(fulllen_bcs: &HashSet<&[u8]>, prefix: &[u8]) -> Option<Barcode> {
        fulllen_bcs
            .iter()
            .find(|bc| &bc[..15] == prefix)
            .map(|bc| Barcode::with_seq(1, BcSeq::from_bytes(bc), true))
    }

    fn is_contam(exact_contig: &[BarcodeContigUMI]) -> Vec<GBIndels> {
        let mut gb_indels: Vec<GBIndels> = Vec::new();
        let mut contam_bcs: HashSet<Barcode> = HashSet::new();
        for source in find_source_bcs(exact_contig) {
            if contam_bcs.contains(&source.barcode) {
                continue;
            }
            let indels = Indels::generate(source.barcode.content.sequence());
            let sink_umi_map: HashMap<Barcode, usize> = exact_contig
                .iter()
                .filter(|bc| !contam_bcs.contains(&bc.barcode))
                .filter(|bc| source.umi >= bc.umi * GB_UMI_MULT)
                .map(|bc| (bc.barcode, bc.umi))
                .collect();
            let mut sink_fulllen: HashSet<&[u8]> = sink_umi_map
                .keys()
                .map(barcode::Barcode::sequence_bytes)
                .collect();
            let mut sink_truncated: HashSet<&[u8]> = sink_umi_map
                .keys()
                .map(|bc| &bc.sequence_bytes()[..bc.sequence_bytes().len() - 1])
                .collect();

            let mut check = |variants: &[BcSeq], error_mode: &IndelErrorMode| {
                for var in variants {
                    if let Some(sink_bc) = match error_mode {
                        IndelErrorMode::Insertion | IndelErrorMode::InsertionWithCorrection => {
                            if sink_fulllen.contains(var.as_bytes()) {
                                let sink_bc = Barcode::with_seq(1, *var, true);
                                sink_fulllen.remove(var.as_bytes());
                                sink_truncated.remove(&var.as_bytes()[..var.len() - 1]);
                                Some(sink_bc)
                            } else {
                                None
                            }
                        }
                        IndelErrorMode::Deletion | IndelErrorMode::DeletionWithCorrection => {
                            if sink_truncated.contains(var.as_bytes()) {
                                let sink_bc =
                                    find_matching_prefix(&sink_fulllen, var.as_bytes()).unwrap();
                                sink_truncated.remove(var.as_bytes());
                                sink_fulllen.remove(sink_bc.sequence_bytes());
                                Some(sink_bc)
                            } else {
                                None
                            }
                        }
                    } {
                        contam_bcs.insert(sink_bc);
                        gb_indels.push(GBIndels {
                            sink_barcode: sink_bc,
                            sink_umis: sink_umi_map[&sink_bc],
                            source_barcode: source.barcode,
                            source_umis: source.umi,
                            error_mode: *error_mode,
                        });
                    }
                }
            };

            check(&indels.insertions, &IndelErrorMode::Insertion);
            check(
                &indels.insertions_1hd,
                &IndelErrorMode::InsertionWithCorrection,
            );
            check(&indels.deletions, &IndelErrorMode::Deletion);
            check(
                &indels.deletions_1hd,
                &IndelErrorMode::DeletionWithCorrection,
            );
        }
        gb_indels
    }

    exact_contigs
        .into_iter()
        .filter(|(_, bcs)| bcs.len() > 1)
        .filter(|(_, bcs)| {
            let max_umis = bcs.first().unwrap().umi;
            bcs.iter().any(|bc| max_umis >= GB_UMI_MULT * bc.umi)
        })
        .flat_map(|(_, bcs)| is_contam(&bcs))
        .map(|indel| {
            if let Some(ref mut logger) = filter_logger {
                logger.log(&FilterLogEntry::cell_calling(
                    indel.sink_barcode.to_string(),
                    AsmCellFilter::GelBeadIndel {
                        source_sequence: indel.source_barcode.to_string(),
                        error_mode: indel.error_mode,
                        sink_contig_umis: indel.sink_umis,
                        source_contig_umis: indel.source_umis,
                    },
                ));
            }
            indel.sink_barcode.to_string()
        })
        .collect()
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
        let j = next_diff1_5(&all, i);

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
        let j = next_diff1_5(&all, i);
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
            let s = next_diff1_2(&vu, r);
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
            if all_chimdata[k].is_cell_and_productive {
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::str::FromStr;
    use vdj_asm_utils::exact_clonotyping::ProductiveContig;
    use vdj_types::VdjChain;

    fn generate_data(
        source_bc: &str,
        sink_bc: &str,
        source_umi: Option<usize>,
        sink_umi: Option<usize>,
    ) -> HashMap<ProductiveContig, Vec<BarcodeContigUMI>> {
        let prod_contig = ProductiveContig {
            chain: VdjChain::IGH,
            cdr3_nt: String::from("ATGC"),
            vdj_nt: String::from("ATGCATGC"),
            c_ref_name: Some("IGHA1".to_owned()),
            jc_delta: Some(0),
        };
        let this_test = vec![
            BarcodeContigUMI {
                barcode: Barcode::from_str(source_bc).unwrap(),
                umi: source_umi.unwrap_or(100),
            },
            BarcodeContigUMI {
                barcode: Barcode::from_str(sink_bc).unwrap(),
                umi: sink_umi.unwrap_or(10),
            },
        ];
        HashMap::from([(prod_contig, this_test)])
    }

    #[test]
    fn test_1bp_insertion() {
        // AAACCATTCACTCATC 16-bp sink sequence
        // |||||||||||||| |
        // AAACCATTCACTCA-C 15-bp source sequence
        let exact_contigs = generate_data("AAACCATTCACTCACT-1", "AAACCATTCACTCATC-1", None, None);
        let filtered_bcs = whitelist_indel_filter(exact_contigs, None);
        assert_eq!(filtered_bcs, vec![String::from("AAACCATTCACTCATC-1")]);
    }

    #[test]
    fn test_1bp_insertion_with_correction() {
        // AAACCATTCACTCTCG 16-bp sink sequence
        // |||||||||||||.|
        // AAACCATTCACTCAC- 15-bp source sequence
        let exact_contigs = generate_data("AAACCATTCACTCACT-1", "AAACCATTCACTCTCG-1", None, None);
        let filtered_bcs = whitelist_indel_filter(exact_contigs, None);
        assert_eq!(filtered_bcs, vec![String::from("AAACCATTCACTCTCG-1")]);
    }

    #[test]
    fn test_1bp_deletion() {
        // CCTACATTC-ATAAGT 15-bp sink sequence
        // ||||||||| ||||||
        // CCTACATTCCATAAGT 16-bp source sequence
        let exact_contigs = generate_data("CCTACATTCCATAAGT-1", "CCTACATTCATAAGTC-1", None, None);
        let filtered_bcs = whitelist_indel_filter(exact_contigs, None);
        assert_eq!(filtered_bcs, vec![String::from("CCTACATTCATAAGTC-1")]);
    }

    #[test]
    fn test_1bp_deletion_with_correction() {
        // CC-ACATTTCATAAGT 15-bp sink sequence
        // || |||||.|||||||
        // CCTACATTCCATAAGT 16-bp source sequence
        let exact_contigs = generate_data("CCTACATTCCATAAGT-1", "CCACATTTCATAAGTC-1", None, None);
        let filtered_bcs = whitelist_indel_filter(exact_contigs, None);
        assert_eq!(filtered_bcs, vec![String::from("CCACATTTCATAAGTC-1")]);
    }
}
