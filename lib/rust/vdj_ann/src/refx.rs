// Copyright (c) 2018 10X Genomics, Inc. All rights reserved.

// This file contains code to make reference data.
//
// ◼ Experiment with building a reference from scratch.  This would be a better
// ◼ solution than ad hoc editing of a flawed reference.
//
// ◼ Document reference sequence requirements so that a customer who wishes to
// ◼ create a reference for a new species will know the conventions used by the
// ◼ code.

use debruijn::dna_string::DnaString;
use debruijn::kmer::Kmer12;
use io_utils::read_to_string_safe;
use kmer_lookup::make_kmer_lookup_12_single;
use std::collections::{HashMap, HashSet};
use std::path::Path;
use string_utils::TextUtils;
use vector_utils::erase_if;

// RefData: this is a packaging of reference data appropriate for VDJ analysis.

#[derive(Default)]
pub struct RefData {
    pub refs: Vec<DnaString>,
    pub rheaders: Vec<String>,
    pub rheaders_orig: Vec<String>,
    pub rkmers_plus: Vec<(Kmer12, i32, i32)>,
    // which V segments have matched UTRs in the reference:
    pub has_utr: HashMap<String, bool>,
    pub name: Vec<String>,
    pub segtype: Vec<&'static str>, // U, V, D, J or C
    pub rtype: Vec<i32>,            // index in "IGH","IGK","IGL","TRA","TRB","TRD","TRG" or -1
    pub igjs: Vec<usize>,           // index of all IGJ segments
    pub cs: Vec<usize>,             // index of all C segments
    pub ds: Vec<usize>,             // index of all D segments
    pub id: Vec<i32>,               // the number after ">" on the header line
    pub transcript: Vec<String>,    // transcript name from header line
}

impl RefData {
    pub fn new() -> RefData {
        RefData {
            refs: Vec::<DnaString>::new(),
            rheaders: Vec::<String>::new(),
            rheaders_orig: Vec::<String>::new(),
            rkmers_plus: Vec::<(Kmer12, i32, i32)>::new(),
            has_utr: HashMap::<String, bool>::new(),
            name: Vec::<String>::new(),
            segtype: Vec::<&'static str>::new(),
            rtype: Vec::<i32>::new(),
            igjs: Vec::<usize>::new(),
            cs: Vec::<usize>::new(),
            ds: Vec::<usize>::new(),
            id: Vec::<i32>::new(),
            transcript: Vec::<String>::new(),
        }
    }
    pub fn is_u(&self, i: usize) -> bool {
        self.segtype[i] == "U"
    }
    pub fn is_v(&self, i: usize) -> bool {
        self.segtype[i] == "V"
    }
    pub fn is_d(&self, i: usize) -> bool {
        self.segtype[i] == "D"
    }
    pub fn is_j(&self, i: usize) -> bool {
        self.segtype[i] == "J"
    }
    pub fn is_c(&self, i: usize) -> bool {
        self.segtype[i] == "C"
    }

    pub fn from_fasta(path: impl AsRef<Path>) -> Self {
        let path = path.as_ref();
        let mut refdata = RefData::new();
        let path_contents = read_to_string_safe(path);
        assert!(
            !path_contents.is_empty(),
            "Reference file at {} has zero length.",
            path.to_string_lossy()
        );
        make_vdj_ref_data_core(
            &mut refdata,
            &path_contents,
            "",
            true, // is_tcr
            true, // is_bcr
            None,
        );
        refdata
    }

    pub fn from_fasta_with_filter(path: impl AsRef<Path>, ids_to_use: &HashSet<i32>) -> Self {
        let mut refdata = RefData::new();
        let path = path.as_ref();
        let path_contents = read_to_string_safe(path);
        assert!(
            !path_contents.is_empty(),
            "Reference file at {} has zero length.",
            path.to_string_lossy()
        );
        make_vdj_ref_data_core(
            &mut refdata,
            &path_contents,
            "",
            true, // is_tcr
            true, // is_bcr
            Some(ids_to_use),
        );
        refdata
    }
}

// ids_to_use_opt: Optional hashSet of ids. If specified only reference
// entries with id in the HashSet is used to construct RefData

pub fn make_vdj_ref_data_core(
    refdata: &mut RefData,
    ref_fasta: &str,
    extended_ref_fasta: &str,
    is_tcr: bool,
    is_bcr: bool,
    ids_to_use_opt: Option<&HashSet<i32>>,
) {
    // Define convenient abbreviations.

    let refs = &mut refdata.refs;
    let rheaders = &mut refdata.rheaders;
    let rkmers_plus = &mut refdata.rkmers_plus;

    // Parse the fasta file.

    refs.clear();
    rheaders.clear();
    read_fasta_contents_into_vec_dna_string_plus_headers(ref_fasta, refs, rheaders);

    // Filter by ids if requested.

    if let Some(ids_to_use) = ids_to_use_opt {
        let mut to_delete = vec![false; refs.len()];
        for i in 0..refs.len() {
            let id = rheaders[i].before("|").force_i32();
            to_delete[i] = !ids_to_use.contains(&id);
        }
        erase_if(refs, &to_delete);
        erase_if(rheaders, &to_delete);
    }
    refdata.rheaders_orig.clone_from(rheaders);

    // Now build stuff.

    let mut rheaders2 = Vec::<String>::new();
    let types = vdj_types::VDJ_CHAINS;
    refdata.rtype = vec![-1_i32; refs.len()];
    for i in 0..rheaders.len() {
        let v: Vec<&str> = rheaders[i].split_terminator('|').collect();
        refdata.name.push(v[2].to_string());
        rheaders2.push(format!("|{}|{}|{}|", v[0], v[2], v[3]));
        match v[3] {
            "5'UTR" => {
                refdata.segtype.push("U");
            }
            "L-REGION+V-REGION" => {
                refdata.segtype.push("V");
            }
            "D-REGION" => {
                refdata.segtype.push("D");
            }
            "J-REGION" => {
                refdata.segtype.push("J");
            }
            "C-REGION" => {
                refdata.segtype.push("C");
            }
            _ => {
                refdata.segtype.push("?");
            }
        }
        for (j, type_) in types.iter().enumerate() {
            if rheaders2[i].contains(type_) {
                refdata.rtype[i] = j as i32;
            }
        }
        refdata.transcript.push(v[1].after(" ").to_string());
    }
    *rheaders = rheaders2;

    // Filter by TCR/BCR.

    if !is_tcr || !is_bcr {
        let mut to_delete = vec![false; refs.len()];
        for i in 0..refs.len() {
            if !is_tcr && (rheaders[i].contains("|TR") || rheaders[i].starts_with("TR")) {
                to_delete[i] = true;
            }
            if !is_bcr && (rheaders[i].contains("|IG") || rheaders[i].starts_with("IG")) {
                to_delete[i] = true;
            }
        }
        erase_if(refs, &to_delete);
        erase_if(rheaders, &to_delete);
        erase_if(&mut refdata.name, &to_delete);
        erase_if(&mut refdata.segtype, &to_delete);
        erase_if(&mut refdata.transcript, &to_delete);
        erase_if(&mut refdata.rtype, &to_delete);
        erase_if(&mut refdata.rheaders_orig, &to_delete);
    }

    // Fill in igjs and cs and ds.

    for i in 0..rheaders.len() {
        if refdata.segtype[i] == "J" && refdata.rtype[i] >= 0 && refdata.rtype[i] < 3 {
            refdata.igjs.push(i);
        }
        if refdata.segtype[i] == "C" {
            refdata.cs.push(i);
        }
        if refdata.segtype[i] == "D" {
            refdata.ds.push(i);
        }
    }

    // Fill in id.

    for header in rheaders.iter() {
        refdata.id.push(header.between("|", "|").force_i32());
    }

    // Extend the reference.

    if !extended_ref_fasta.is_empty() {
        let mut refs2 = Vec::<DnaString>::new();
        let mut rheaders2 = Vec::<String>::new();
        read_fasta_contents_into_vec_dna_string_plus_headers(
            extended_ref_fasta,
            &mut refs2,
            &mut rheaders2,
        );
        refs.append(&mut refs2);
        rheaders.append(&mut rheaders2);
        // ◼ Note not appending to refdata.name.  This may be a bug.
    }

    // Make lookup table for reference.

    make_kmer_lookup_12_single(refs, rkmers_plus);

    // Determine which V segments have matching UTRs in the reference.

    for header in rheaders.iter() {
        if !header.contains("segment") {
            let name = header.after("|").between("|", "|");
            if header.contains("UTR") {
                refdata.has_utr.insert(name.to_string(), true);
            }
        }
    }
    for header in rheaders.iter() {
        if !header.contains("segment") {
            let name = header.after("|").between("|", "|");
            if header.contains("V-REGION") {
                refdata.has_utr.entry(name.to_string()).or_insert(false);
            }
        }
    }
}

pub fn read_fasta_contents_into_vec_dna_string_plus_headers(
    f: &str,
    dv: &mut Vec<DnaString>,
    headers: &mut Vec<String>,
) {
    let mut last: String = String::new();
    let mut first = true;
    for s in f.split('\n') {
        if first {
            assert!(s.starts_with('>'), "fasta format failure reading {f}");
            first = false;
            headers.push(s.get(1..).unwrap().to_string());
        } else if s.starts_with('>') {
            dv.push(DnaString::from_dna_string(&last));
            last.clear();
            headers.push(s.get(1..).unwrap().to_string());
        } else {
            last += s;
        }
    }
    dv.push(DnaString::from_dna_string(&last));
}
