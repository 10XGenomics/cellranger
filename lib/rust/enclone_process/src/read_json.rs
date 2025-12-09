// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
#![deny(missing_docs)]

use self::annotate::Annotation;
use self::refx::RefData;
use crate::core::barcode_fate::BarcodeFate;
use crate::core::defs::{BarcodeContigs, Contig, EncloneControl};
use crate::core::enclone_structs::BarcodeFates;
use crate::enclone::misc3::sort_tig_bc;
use debruijn::dna_string::DnaString;
use martian_filetypes::LazyFileTypeIO;
use martian_filetypes::json_file::JsonFile;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fmt::Write;
use string_utils::TextUtils;
use vdj_ann::annotate::{ContigAnnotation, get_cdr3_using_ann};
use vdj_ann::{annotate, refx};
use vdj_types::{VdjChain, VdjReceptor, VdjRegion};
use vector_utils::{bin_position, unique_sort};

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

fn json_error(json: Option<&str>, msg: &str) -> anyhow::Error {
    let mut msgx =
        "There is something wrong with the contig annotations in the cellranger output file"
            .to_string();
    if let Some(json) = json {
        write!(msgx, "\n{json}.").unwrap();
    } else {
        msgx += ".";
    }
    msgx += "\n\n";
    msgx += msg;
    msgx += "\n\nHere is what you should do:\n\n\
        1. If you used cellranger version ≥ 4.0, the problem is very likely\n\
        that the directory outs/vdj_reference was not retained, so enclone\n\
        didn't see it, and had to guess what the reference sequence was.\n\
        Fix this and everything should be fine.\n\n\
        2. If you used cellranger version 3.1, then you need to add a command-line\n\
        argument REF=<vdj_reference_fasta_file_name>, or if you already did that,\n\
        make sure it is the *same* as that which you gave cellranger.\n\n\
        3. If you used cellranger version < 3.1 (the only other possibility), then\n\
        you have options:\n\
        • rerun cellranger using the current version\n\
        • or provide an argument REF= as above and RE to force reannotation\n\
        • or provide the argument BUILT_IN to use the current reference and force\n  \
            reannotation (and MOUSE if you used mouse); only works with human and mouse.\n\n\
        Note that one way to get the error is to specify TCR when you meant BCR, or the\n\
        other way.\n\n\
        If you're stuck, please write to us at enclone@10xgenomics.com.\n";

    anyhow::format_err!("{msgx}")
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

#[derive(Default)]
struct JsonParseResult {
    vdj_cell: Option<String>,
    gex_cell: Option<String>,
    gex_cells_specified: bool,
    tig: Option<Contig>,
}

#[allow(clippy::too_many_arguments)]
fn process_json_annotation(
    ann: ContigAnnotation,
    json_fpath: &str,
    ctl: &EncloneControl,
    dataset_index: usize,
    refdata: &RefData,
    to_ref_index: &HashMap<usize, usize>,
    donor_list: &[String],
) -> anyhow::Result<JsonParseResult> {
    let mut res: JsonParseResult = Default::default();

    // Get cell status.  Sometime after CR 4.0 was released, and before 4.1 was released,
    // we added new fields is_asm_cell and is_gex_cell to the json file.  The value of
    // is_asm_cell is the original determination of "cell" in the VDJ pipeline, whereas the
    // value of is_gex_cell is that for the GEX pipeline.
    let mut is_cell = ann.is_cell;
    if ann.is_asm_cell.is_some_and(|is_asm_cell| is_asm_cell) {
        is_cell = true;
    }

    if let Some(is_gex_cell) = ann.is_gex_cell {
        res.gex_cells_specified = true;
        if is_gex_cell {
            res.gex_cell = Some(ann.barcode.clone());
        }
    }

    if !is_cell {
        return Ok(res);
    }

    res.vdj_cell = Some(ann.barcode.clone());

    // Proceed.

    if !ann.is_productive() {
        return Ok(res);
    }
    if !ann.high_confidence {
        return Ok(res);
    }
    let mut left = false;
    let (mut v_ref_id, mut j_ref_id) = (1000000, 0);
    let mut d_ref_id: Option<usize> = None;
    let mut c_ref_id = None;
    let mut c_ref_name = None;
    let mut chain_type = String::new();
    let mut u_ref_id = None;
    let (mut tig_start, mut tig_stop) = (-1_isize, -1_isize);
    let mut v_stop = 0;
    let mut v_stop_ref = 0;
    let mut d_start = None;
    let mut j_start = 0;
    let mut j_start_ref = 0;
    let mut c_start = None;
    let mut annv = Vec::<Annotation>::new();
    let mut cdr3_aa: String;
    let mut cdr3_dna: String;
    let mut cdr3_start: usize;

    {
        // Use annotations from json file.

        cdr3_aa = ann.cdr3.unwrap();
        cdr3_dna = ann.cdr3_seq.unwrap();
        cdr3_start = ann.cdr3_start.unwrap();
        let annotations = ann.annotations;
        assert!(!annotations.is_empty());

        let mut cigarv = String::new(); // cigar for V segment
        for a in annotations {
            let region_type = a.feature.region_type;
            let feature_id = a.feature.feature_id;
            if !to_ref_index.contains_key(&feature_id) {
                continue;
            }
            let feature_idx = to_ref_index[&feature_id];
            let ref_start = a.annotation_match_start;
            if region_type == VdjRegion::V {
                v_stop = a.contig_match_end;
                v_stop_ref = a.annotation_match_end;
            }
            let gene_name = a.feature.gene_name;
            if refdata.name[feature_idx] != gene_name {
                return Err(anyhow::format_err!(
                    "\nThere is an inconsistency between the reference \
                     file used to create the Cell Ranger output files in\n{}\nand the \
                     reference that enclone is using.\n\nFor example, the feature \
                     numbered {} is\nthe gene {} in one and the gene {} in the other.\n\n\
                     As far as we know, this type of error can only occur with Cell Ranger \
                     versions before 4.0.\n\n\
                     If this is mouse data, please use the argument MOUSE, and that may \
                     solve the problem.\n\n\
                     If this is human or mouse data, and you are OK with using the current \
                     built-in reference that\nenclone has, \
                     you can instead add the argument BUILT_IN to the command line.  This \
                     forces\nrecomputation of annotations and may be somewhat slower.\n\n\
                     A solution that should always work is to supply\n\
                     REF=vdj_reference_fasta_filename as an argument to enclone.\n",
                    json_fpath.rev_before("/"),
                    feature_id,
                    gene_name,
                    refdata.name[feature_idx]
                ));
            }
            if region_type == VdjRegion::V && ref_start == 0 {
                let chain = a.feature.chain;
                chain_type = chain.to_string();
                tig_start = a.contig_match_start as isize;
                cdr3_start -= tig_start as usize;
                if chain == VdjChain::IGH
                    || chain == VdjChain::TRB
                    || (chain == VdjChain::TRD && ctl.cr_opt.input.receptor == VdjReceptor::TRGD)
                {
                    left = true;
                }
                v_ref_id = feature_idx;
                cigarv = a.cigar;
            } else {
                // also check for IG chain?????????????????????????????????????????
                let ref_stop = a.annotation_match_end;
                let ref_len = a.annotation_length;
                if region_type == VdjRegion::J && ref_stop == ref_len {
                    tig_stop = a.contig_match_end as isize;
                    j_ref_id = feature_idx;
                    j_start = a.contig_match_start;
                    j_start_ref = a.annotation_match_start;
                }
                if region_type == VdjRegion::UTR {
                    u_ref_id = Some(feature_idx);
                }
                if region_type == VdjRegion::D {
                    d_start = Some(a.contig_match_start);
                    d_ref_id = Some(feature_idx);
                }
                if region_type == VdjRegion::C {
                    c_ref_id = Some(feature_idx);
                    c_ref_name = Some(refdata.name[feature_idx].clone());
                    c_start = Some(a.contig_match_start);
                }
            }
        }
        if v_ref_id == 1000000 {
            return Ok(res);
        }

        // Compute annv from cigarv.  We don't compute the mismatch entry.

        let mut cg = Vec::<Vec<u8>>::new(); // pieces of cigar string
        let mut piece = Vec::<u8>::new();
        for c in cigarv.chars() {
            piece.push(c as u8);
            if c.is_ascii_alphabetic() {
                cg.push(piece.clone());
                piece.clear();
            }
        }
        let t = v_ref_id as i32;
        let (mut len1, mut len2) = (0, 0);
        let (mut ins, mut del) = (0, 0);
        for cgi in cg {
            let x = std::str::from_utf8(&cgi[0..cgi.len() - 1])
                .unwrap()
                .force_i32();
            if cgi[cgi.len() - 1] == b'M' {
                if len1 == 0 {
                    len1 = x;
                } else if len2 == 0 {
                    len2 = x;
                } else {
                    // probably can't happen
                    len1 = 0;
                    len2 = 0;
                    break;
                }
            }
            if cgi[cgi.len() - 1] == b'I' {
                ins = x;
            }
            if cgi[cgi.len() - 1] == b'D' {
                del = x;
            }
        }
        annv.push(Annotation {
            tig_start: 0,
            match_len: len1,
            ref_id: t,
            ref_start: 0,
            mismatches: 0,
        });
        if ins > 0 && ins % 3 == 0 && del == 0 && len2 > 0 {
            let start = len1 + ins;
            annv.push(Annotation {
                tig_start: start,
                match_len: len2,
                ref_id: t,
                ref_start: len1,
                mismatches: 0,
            });
        } else if del > 0 && del % 3 == 0 && ins == 0 && len2 > 0 {
            annv.push(Annotation {
                tig_start: len1,
                match_len: len2,
                ref_id: t,
                ref_start: len1 + del,
                mismatches: 0,
            });
        }
        let rt = &refdata.refs[v_ref_id];
        if annv.len() == 2 && annv[0].match_len as usize > rt.len() {
            let msg = format!("annv[0].1 = {}, rt.len() = {}", annv[0].match_len, rt.len());
            return Err(json_error(None, &msg));
        }

        // Check to see if the CDR3 sequence has changed.  This could happen if the cellranger
        // version for all_contig_annotations.json used an older version of the CDR3 calculation
        // than is used in the current version of enclone.  This could result in internal
        // inconsistencies, leading to an assert somewhere downstream.

        let x = DnaString::from_dna_string(&ann.sequence);
        let found_cdr3s = get_cdr3_using_ann(&x, refdata, &annv);
        if found_cdr3s.is_empty() {
            return Ok(res);
        }
        let cdr3 = found_cdr3s.first().unwrap();
        let cdr3_aa_alt = std::str::from_utf8(&cdr3.aa_seq).unwrap();
        if cdr3_aa != cdr3_aa_alt {
            // This is particularly pathological and rare:

            if tig_start as usize > cdr3.start_position_on_contig {
                return Ok(res);
            }

            // Define start.

            cdr3_start = cdr3.start_position_on_contig - tig_start as usize;

            // Define cdr3.

            cdr3_aa = cdr3_aa_alt.to_string();
            cdr3_dna = x
                .slice(cdr3_start, cdr3_start + 3 * cdr3_aa.len())
                .to_string();
        }
    }

    // Test for two very rare conditions where the CDR3 is busted.  This could be confusing to
    // users if they hit one of these.
    // Case 1: seen on 47680, barcode CGCCAAGTCCATGAAC-1.
    // Case 2: seen on 99640, barcode CAGTAACCATGTCGAT-1.
    // It is not known if these correspond to bugs in cellranger that were subsequently fixed.

    if cdr3_aa.contains('*') {
        return Ok(res);
    }
    if cdr3_start + 3 * cdr3_aa.len() > tig_stop as usize - tig_start as usize {
        return Ok(res);
    }

    // Keep going.

    if tig_start < 0 || tig_stop < 0 {
        let msg = format!("tig_start = {tig_start}, tig_stop = {tig_stop}");
        return Err(json_error(Some(json_fpath), &msg));
    }
    let (tig_start, tig_stop) = (tig_start as usize, tig_stop as usize);
    let mut quals = ann.quals.as_bytes().to_vec();
    assert_eq!(ann.sequence.len(), ann.quals.len());
    let seq = &ann.sequence[tig_start..tig_stop].to_string();
    for qual in &mut quals {
        *qual -= 33_u8;
    }
    let full_quals = quals;
    let quals = full_quals[tig_start..tig_stop].to_vec();
    let umi_count = ann.umi_count;
    let read_count = ann.read_count;
    let origin = if !ctl.origin_info[dataset_index].origin_id.is_empty() {
        Some(&ctl.origin_info[dataset_index].origin_id)
    } else {
        None
    };
    let donor = if !ctl.origin_info[dataset_index].origin_id.is_empty() {
        Some(&ctl.origin_info[dataset_index].donor_id)
    } else {
        None
    };
    let mut donor_index = None;
    if origin.is_some()
        && let Some(donor) = donor
    {
        donor_index = Some(bin_position(donor_list, donor) as usize);
    }

    res.tig = Some(Contig {
        cdr3_dna,
        len: seq.len(),
        v_start: tig_start,
        v_stop,
        v_stop_ref,
        d_start,
        j_start,
        j_start_ref,
        j_stop: tig_stop,
        c_start,
        full_seq: ann.sequence.as_bytes().to_vec(),
        v_ref_id,
        d_ref_id,
        j_ref_id,
        c_ref_id,
        c_ref_name,
        u_ref_id,
        fr1_start: 0,
        cdr1_start: None,
        fr2_start: None,
        cdr2_start: None,
        fr3_start: None,
        cdr3_aa,
        cdr3_start,
        quals,
        full_quals,
        barcode: ann.barcode,
        tigname: ann.contig_name,
        left,
        dataset_index,
        donor_index,
        umi_count,
        read_count,
        chain_type,
        annv,
        mix_donors: ctl.cr_opt.mix_donors,
    });
    Ok(res)
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Parse the JSON annotations file.
//
// Tracking contigs using bc_cdr3_aa; could improve later.
//
// This section requires 3.1.  If you want to avoid that, do something to make tig_start
// and tig_stop always nonnegative.  Or use the RE option.

#[derive(Default)]
struct ReadJsonResult {
    vdj_cells: Vec<String>,
    gex_cells: Vec<String>,
    gex_cells_specified: bool,
    contig_data_by_barcode: Vec<BarcodeContigs>,
}

fn read_json(
    ctl: &EncloneControl,
    dataset_index: usize,
    json: &JsonFile<Vec<ContigAnnotation>>,
    refdata: &RefData,
    to_ref_index: &HashMap<usize, usize>,
) -> anyhow::Result<ReadJsonResult> {
    let mut tigs = Vec::new();
    let mut vdj_cells = Vec::new();
    let mut gex_cells = Vec::new();
    let mut gex_cells_specified = false;

    let reader = json.lazy_reader()?;
    let donor_list = ctl.donor_list();

    for entry in reader {
        let result = process_json_annotation(
            entry?,
            json.as_ref().to_str().unwrap(),
            ctl,
            dataset_index,
            refdata,
            to_ref_index,
            &donor_list,
        )?;
        if let Some(tig) = result.tig {
            tigs.push(tig);
        }
        if let Some(vdj_cell) = result.vdj_cell {
            vdj_cells.push(vdj_cell);
        }
        if let Some(gex_cell) = result.gex_cell {
            gex_cells.push(gex_cell);
        }
        if result.gex_cells_specified {
            gex_cells_specified = true;
        }
    }
    unique_sort(&mut gex_cells);
    let mut contig_data_by_barcode = Vec::<BarcodeContigs>::new();
    let mut r = 0;
    while r < tigs.len() {
        let mut s = r + 1;
        while s < tigs.len() {
            if tigs[s].barcode != tigs[r].barcode {
                break;
            }
            s += 1;
        }

        // For now we require at most four contigs (but we don't yet merge foursies).

        if s - r <= 4 {
            let mut bc_tigs: BarcodeContigs = tigs[r..s].to_vec();
            bc_tigs.sort();
            contig_data_by_barcode.push(bc_tigs);
        }
        r = s;
    }
    sort_tig_bc(&mut contig_data_by_barcode);
    unique_sort(&mut vdj_cells);

    // Done.

    Ok(ReadJsonResult {
        vdj_cells,
        gex_cells,
        gex_cells_specified,
        contig_data_by_barcode,
    })
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub(crate) struct Annotations {
    pub(crate) gex_cells: Vec<Vec<String>>,
    pub(crate) gex_cells_specified: Vec<bool>,
    pub(crate) tig_bc: Vec<BarcodeContigs>,
    pub(crate) fate: Vec<BarcodeFates>,
}

pub(crate) fn parse_json_annotations_files(
    ctl: &EncloneControl,
    refdata: &RefData,
) -> anyhow::Result<Annotations> {
    // Note: only tracking truncated seq and quals initially
    let to_ref_index = refdata
        .id
        .iter()
        .take(refdata.refs.len())
        .enumerate()
        .map(|(i, &id)| (id as usize, i))
        .collect();
    let results = ctl
        .origin_info
        .par_iter()
        .enumerate()
        .map(|(li, ds)| read_json(ctl, li, &ds.file, refdata, &to_ref_index).map(|r| (li, r)))
        .collect::<anyhow::Result<Vec<_>>>()?;

    let mut ann = Annotations {
        tig_bc: Default::default(),
        gex_cells: Default::default(),
        gex_cells_specified: Default::default(),
        fate: vec![HashMap::<String, BarcodeFate>::new(); ctl.origin_info.len()],
    };

    for (_, result) in results {
        let cells = &result.vdj_cells;
        let mut found = vec![false; cells.len()];
        let tigs = &result.contig_data_by_barcode;
        for bc_contigs in tigs {
            let p = bin_position(cells, &bc_contigs[0].barcode);
            if p >= 0 {
                found[p as usize] = true;
            }
        }

        ann.tig_bc.extend(result.contig_data_by_barcode.into_iter());
        ann.gex_cells.push(result.gex_cells);
        ann.gex_cells_specified.push(result.gex_cells_specified);
    }
    Ok(ann)
}
