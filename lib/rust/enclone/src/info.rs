// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// This file provides the single function build_info.

use self::refx::RefData;
use amino::{codon_to_aa, nucleotide_to_aminoacid_sequence};
use debruijn::dna_string::DnaString;
use debruijn::Mer;
use enclone_core::barcode_fate::BarcodeFate;
use enclone_core::defs::{CloneInfo, ExactClonotype};
use enclone_core::enclone_structs::BarcodeFates;
use rayon::prelude::*;
use std::convert::TryInto;
use string_utils::strme;
use vdj_ann::refx;
use vector_utils::unique_sort;

struct Result {
    clone_info: Vec<CloneInfo>,
    exact_clonotype: ExactClonotype,
    fate: Vec<ImproperFateResult>,
}

struct ImproperFateResult {
    dataset_id: usize,
    barcode: String,
}

/// Build info about clonotypes.  We create a data structure info.
/// An entry in info is a clonotype having appropriate properties.
///
/// Much of the information in a CloneInfo object is redundant.  So we could probably
/// improve both time and space computational performance by reducing that redundancy.
pub fn build_info(
    refdata: &RefData,
    exact_clonotypes: &mut [ExactClonotype],
    fate: &mut [BarcodeFates],
) -> Vec<CloneInfo> {
    let results: Vec<_> = exact_clonotypes
        .par_iter()
        .enumerate()
        .map(|(i, exact_clonotype)| {
            let mut lens = Vec::<usize>::new();
            let mut tigs = Vec::<Vec<u8>>::new();
            let mut tigs_amino = Vec::<Vec<u8>>::new();
            let mut aa_mod_indel = Vec::<Vec<u8>>::new();
            let mut tigs_ins = Vec::<Vec<(usize, Vec<u8>)>>::new();
            let mut tigsp = Vec::<DnaString>::new();
            let mut has_del = Vec::<bool>::new();
            let (mut vs, mut js) = (Vec::<DnaString>::new(), Vec::<DnaString>::new());
            let mut vsids = Vec::<usize>::new();
            let mut jsids = Vec::<usize>::new();
            let mut cdr3s = Vec::<String>::new();
            let mut updated_exact_clonotype = exact_clonotype.clone();
            for j in 0..updated_exact_clonotype.share.len() {
                let x = &mut updated_exact_clonotype.share[j];
                tigsp.push(DnaString::from_acgt_bytes(&x.seq));
                let jid = x.j_ref_id;
                js.push(refdata.refs[jid].clone());

                // If there is a deletion in a V segment region, edit the contig sequence,
                // inserting hyphens where the deletion is, and if there is an insertion, delete it.

                let vid = x.v_ref_id;
                let jid = x.j_ref_id;
                let mut annv = x.annv.clone();
                vsids.push(vid);
                jsids.push(jid);
                // DELETION
                if annv.len() == 2 && annv[1].tig_start == annv[0].tig_start + annv[0].match_len {
                    let mut t = Vec::<u8>::new();
                    let (mut del_start, mut del_stop) = (annv[0].match_len, annv[1].ref_start);
                    t.extend(&x.seq[..del_start.try_into().unwrap()]);
                    t.resize(del_stop.try_into().unwrap(), b'-');
                    t.extend(&x.seq[annv[1].tig_start as usize..]);
                    lens.push(t.len());
                    tigs.push(t.clone());
                    if del_start % 3 != 0 {
                        // Bad solution here, should pick optimal choice.
                        let offset = del_start % 3 - 3;
                        del_start -= offset;
                        del_stop -= offset;
                        t.clear();
                        t.extend(&x.seq[..del_start.try_into().unwrap()]);
                        t.resize(del_stop.try_into().unwrap(), b'-');
                        t.extend(&x.seq[((annv[1].tig_start - offset) as usize)..]);
                    }
                    annv[0].match_len += (del_stop - del_start) + annv[1].match_len;
                    annv.truncate(1);
                    tigs_amino.push(t.clone());
                    let mut aa = Vec::<u8>::new();
                    for p in (0..=t.len() - 3).step_by(3) {
                        if t[p] == b'-' {
                            aa.push(b'-');
                        } else {
                            aa.push(codon_to_aa(&t[p..p + 3]));
                        }
                    }

                    aa_mod_indel.push(aa);
                    tigs_ins.push(Vec::new());
                    has_del.push(true);
                // INSERTION
                } else if annv.len() == 2
                    && annv[1].ref_start == annv[0].ref_start + annv[0].match_len
                {
                    let ins_len =
                        (annv[1].tig_start - annv[0].tig_start - annv[0].match_len) as usize;
                    let ins_pos = (annv[0].tig_start + annv[0].match_len) as usize;
                    let mut t = Vec::<u8>::new();
                    let mut nt = Vec::<u8>::new();
                    for i in 0..x.seq.len() {
                        if i < ins_pos || i >= ins_pos + ins_len {
                            t.push(x.seq[i]);
                        } else {
                            nt.push(x.seq[i]);
                        }
                    }
                    has_del.push(true); // DOES NOT MAKE SENSE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    lens.push(t.len());
                    tigs.push(t.clone());
                    let ins = vec![(ins_pos, nt)];
                    tigs_ins.push(ins);
                    tigs_amino.push(t);

                    // Optimize to compute entry in aa_mod_indel and the inserted aa sequence.

                    let aa_full = nucleotide_to_aminoacid_sequence(&x.seq, 0);
                    let ref_aa =
                        nucleotide_to_aminoacid_sequence(&refdata.refs[vid].to_ascii_vec(), 0);
                    let ins_len_aa = ins_len / 3;
                    const EXT: usize = 10;
                    let ins_pos_low = if ins_pos / 3 < EXT {
                        0
                    } else {
                        ins_pos / 3 - EXT
                    };
                    let mut ins_pos_high =
                        std::cmp::min(ins_pos / 3 + EXT, aa_full.len() - ins_len_aa + 1);
                    ins_pos_high = std::cmp::min(ins_pos_high, ref_aa.len() - ins_len_aa + 1);
                    let mut mis = Vec::<(usize, usize, Vec<u8>)>::new();
                    for j in ins_pos_low..ins_pos_high {
                        let mut y = Vec::<u8>::new();
                        for (k, &aa) in aa_full.iter().enumerate() {
                            if k < j || k >= j + ins_len_aa {
                                y.push(aa);
                            }
                        }
                        let mut m = 0;
                        for l in 0..ref_aa.len() {
                            if l < y.len() && ref_aa[l] != y[l] {
                                m += 1;
                            }
                        }
                        mis.push((m, j, y.clone()));
                    }
                    mis.sort();
                    aa_mod_indel.push(mis[0].2.clone());
                    let ins_aa_pos = mis[0].1;
                    let mut aax = Vec::<u8>::new();
                    let b = 3 * ins_aa_pos;
                    for p in 0..ins_len_aa {
                        emit_codon_color_escape(&x.seq[b + 3 * p..b + 3 * p + 3], &mut aax);
                        let aa = codon_to_aa(&x.seq[b + 3 * p..b + 3 * p + 3]);
                        aax.push(aa);
                        aax.append(&mut b"[0m".to_vec());
                    }

                    // Finish up ann.

                    annv[0].match_len += annv[1].match_len;
                    annv.truncate(1);
                } else {
                    has_del.push(false);
                    lens.push(x.seq.len());
                    tigs.push(x.seq.clone());
                    tigs_amino.push(x.seq.clone());
                    aa_mod_indel.push(nucleotide_to_aminoacid_sequence(&x.seq, 0));
                    tigs_ins.push(Vec::new());
                }

                // Save reference V segment.  However in the case where there is a
                // single indel between the contig and the reference sequence, edit the
                // reference V segment accordingly.

                let rt = &refdata.refs[vid];
                if x.annv.len() == 2 {
                    let mut r = rt.slice(0, x.annv[0].match_len as usize).to_owned();
                    // deletion
                    if x.annv[1].tig_start == x.annv[0].tig_start + x.annv[0].match_len {
                        // DEAD CODE
                        for m in x.annv[1].ref_start as usize..rt.len() {
                            r.push(rt.get(m));
                        }
                        vs.push(r.clone());
                    // insertion
                    } else if x.annv[1].ref_start == x.annv[0].ref_start + x.annv[0].match_len {
                        for m in x.annv[1].ref_start as usize..rt.len() {
                            r.push(rt.get(m));
                        }
                        vs.push(r.clone());
                    } else {
                        // maybe can't happen
                        vs.push(rt.clone());
                    }
                } else {
                    vs.push(rt.clone());
                }
                cdr3s.push(x.cdr3_dna.clone());

                // Modify the exact subclonotype to fill in some members.
                // This is the only place where build_info modifies the exact subclonotype.

                x.seq_del.clone_from(&tigs[tigs.len() - 1]);
                x.seq_del_amino
                    .clone_from(&tigs_amino[tigs_amino.len() - 1]);
                x.aa_mod_indel
                    .clone_from(&aa_mod_indel[aa_mod_indel.len() - 1]);
                x.ins.clone_from(&tigs_ins[tigs_ins.len() - 1]);
                x.vs = vs[vs.len() - 1].clone();
                x.js = js[js.len() - 1].clone();
            }
            let mut origin = Vec::<usize>::new();
            for j in 0..exact_clonotype.clones.len() {
                origin.push(exact_clonotype.clones[j][0].dataset_index);
            }
            unique_sort(&mut origin);
            let shares = &exact_clonotype.share;
            let mut placed = false;
            let mut res = Result {
                clone_info: Default::default(),
                exact_clonotype: updated_exact_clonotype,
                fate: Default::default(),
            };
            for i1 in 0..shares.len() {
                if shares[i1].left {
                    for i2 in 0..shares.len() {
                        if !shares[i2].left {
                            placed = true;
                            let lensx = [lens[i1], lens[i2]].to_vec();
                            let tigsx = [tigs[i1].clone(), tigs[i2].clone()].to_vec();
                            let tigs_aminox =
                                [tigs_amino[i1].clone(), tigs_amino[i2].clone()].to_vec();
                            let tigspx = [tigsp[i1].clone(), tigsp[i2].clone()].to_vec();
                            let has_delx = [has_del[i1], has_del[i2]].to_vec();
                            let vsx = [vs[i1].clone(), vs[i2].clone()].to_vec();
                            let jsx = [js[i1].clone(), js[i2].clone()].to_vec();
                            let cdr3sx = [cdr3s[i1].clone(), cdr3s[i2].clone()].to_vec();
                            let vsidsx = [vsids[i1], vsids[i2]].to_vec();
                            let jsidsx = [jsids[i1], jsids[i2]].to_vec();
                            let exact_cols = vec![i1, i2];
                            res.clone_info.push(CloneInfo {
                                lens: lensx,
                                tigs: tigsx,
                                tigs_amino: tigs_aminox,
                                tigsp: tigspx,
                                has_del: has_delx,
                                clonotype_index: i,
                                exact_cols,
                                origin: origin.clone(),
                                vs: vsx.clone(),
                                js: jsx,
                                vsids: vsidsx,
                                jsids: jsidsx,
                                cdr3s: cdr3sx,
                            });
                        }
                    }
                }
            }

            // Incorporate improper cells if they are onesies.  Note that we're dropping the
            // improper cells having two or more chains.

            if !placed && shares.len() > 1 {
                let ex = &exact_clonotypes[i];
                for j in 0..ex.clones.len() {
                    res.fate.push(ImproperFateResult {
                        dataset_id: ex.clones[j][0].dataset_index,
                        barcode: ex.clones[j][0].barcode.clone(),
                    });
                }
            }
            if !placed && shares.len() == 1 {
                let mut exact_cols = Vec::<usize>::new();
                for i in 0..tigs.len() {
                    exact_cols.push(i);
                }
                res.clone_info.push(CloneInfo {
                    lens,
                    tigs,
                    tigs_amino,
                    tigsp,
                    has_del,
                    clonotype_index: i,
                    exact_cols,
                    origin: origin.clone(),
                    vs: vs.clone(),
                    js,
                    vsids,
                    jsids,
                    cdr3s,
                });
            }
            res
        })
        .collect();

    // Cumulate info.  This is single threaded and could probably be speeded up.
    let mut info = Vec::<CloneInfo>::new();
    for (i, r) in results.into_iter().enumerate() {
        info.extend(r.clone_info);
        exact_clonotypes[i] = r.exact_clonotype;
        for f in r.fate {
            fate[f.dataset_id].insert(f.barcode, BarcodeFate::Improper);
        }
    }

    // Sort info.

    info.par_sort();

    // Done.

    info
}

fn emit_codon_color_escape(c: &[u8], log: &mut Vec<u8>) {
    let mut s = 0;
    if c == b"CTG" {
        s = 3;
    } else if c == b"AGG" {
        s = 1;
    } else if c == b"AGT" {
        s = 2;
    } else {
        for i in 0..3 {
            if c[i] == b'A' {
            } else if c[i] == b'C' {
                s += 1;
            } else if c[i] == b'G' {
                s += 2;
            } else if c[i] == b'T' {
                s += 3;
            } else {
                panic!("Illegal codon: \"{}\".", strme(c));
            }
        }
    }
    let s = s % 6;
    print_color(s, log);
}

fn print_color(s: usize, log: &mut Vec<u8>) {
    assert!(s < 7);
    if s == 0 {
        log.append(&mut b"[01m[38;5;75m".to_vec());
    } else if s == 1 {
        log.append(&mut b"[01m[38;5;166m".to_vec());
    } else if s == 2 {
        // At one point this was made bold, which makes it more readable when printed, but
        // it's uglier if bold and overall contrast is reduced.
        log.append(&mut b"[38;5;178m".to_vec());
    } else if s == 3 {
        log.append(&mut b"[01m[38;5;25m".to_vec());
    } else if s == 4 {
        log.append(&mut b"[01m[38;5;175m".to_vec());
    } else if s == 5 {
        log.append(&mut b"[38;5;36m".to_vec());
    } else {
        log.append(&mut b"[01m[38;5;11m".to_vec());
    }
}
