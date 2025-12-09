// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.
#![deny(missing_docs)]

// Populate features.

use amino::nucleotide_to_aminoacid_sequence;
use vdj_ann::refx::RefData;
use vdj_ann::vdj_features::{cdr1_start, cdr2_start, fr1_start, fr2_start, fr3_start};

#[allow(clippy::too_many_arguments)]
pub(crate) fn populate_features(
    refdata: &RefData,
    fr1_starts: &mut Vec<usize>,
    fr2_starts: &mut Vec<Option<usize>>,
    fr3_starts: &mut Vec<Option<usize>>,
    cdr1_starts: &mut Vec<Option<usize>>,
    cdr2_starts: &mut Vec<Option<usize>>,
) {
    *fr1_starts = vec![0; refdata.refs.len()];
    *fr2_starts = vec![None; refdata.refs.len()];
    *fr3_starts = vec![None; refdata.refs.len()];
    *cdr1_starts = vec![None; refdata.refs.len()];
    *cdr2_starts = vec![None; refdata.refs.len()];
    for i in 0..refdata.refs.len() {
        if refdata.is_v(i) {
            let aa = nucleotide_to_aminoacid_sequence(&refdata.refs[i].to_ascii_vec(), 0);
            let rtype = refdata.rtype[i];
            let chain_type = if rtype == 0 {
                "IGH"
            } else if rtype == 1 {
                "IGK"
            } else if rtype == 2 {
                "IGL"
            } else if rtype == 3 {
                "TRA"
            } else if rtype == 4 {
                "TRB"
            } else {
                continue;
            };
            let Some(fs1) = fr1_start(&aa, chain_type) else {
                continue;
            };
            fr1_starts[i] = 3 * fs1;
            let fs2 = fr2_start(&aa, chain_type, false);
            if let Some(fs2) = fs2 {
                fr2_starts[i] = Some(3 * fs2);
            }

            let fs3 = fr3_start(&aa, chain_type, false);
            if let Some(fs3) = fs3 {
                fr3_starts[i] = Some(3 * fs3);
            }

            let cs1 = cdr1_start(&aa, chain_type, false);
            if let Some(cs1) = cs1 {
                cdr1_starts[i] = Some(3 * cs1);
            }

            let cs2 = cdr2_start(&aa, chain_type, false);
            if let Some(cs2) = cs2 {
                cdr2_starts[i] = Some(3 * cs2);
            }
        }
    }
}
