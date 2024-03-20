use super::feature_reference::FeatureType;
use crate::constants::ILLUMINA_QUAL_OFFSET;
use crate::reference::feature_reference::{
    FeatureDef, FeatureReference, LIBRARY_TYPES_WITHOUT_FEATURES,
};
use crate::rna_read::RnaRead;
use anyhow::{anyhow, bail, Result};
use fastq_set::read_pair::{ReadPart, WhichRead};
use itertools::Itertools;
use regex::{Regex, RegexSet};
use serde::{Deserialize, Serialize};
use std::cmp::{self, Reverse};
use std::collections::hash_map::Entry;
use std::collections::{HashMap, HashSet};
use std::str;
use std::string::String;
use std::sync::Arc;

pub const NUCLEOTIDES: &[u8] = b"ACGT";

const FEATURE_CONF_THRESHOLD: f64 = 0.975;
const FEATURE_MAX_QV: u8 = 33; // cap observed QVs at this value

const REGEX_TOO_BIG_ERR: &str =
    "Patterns in the Feature Reference file cannot be compiled as they exceed the \
memory limit used for constructing the regular expression.  Please check your feature reference \
 file and contact support@10xgenomics.com for help.";

/// Correct a feature barcode against a feature reference.
/// * `seq` - ASCII nucleotide sequence
/// * `qual` - ASCII ILMN Phred qual string
/// * `feat_dist` - Feature barcode distrubution as proportions (within-feature-type)
/// * `feat_type_idx` - Integer indicating the feature type
pub fn correct_feature_barcode<'a>(
    bcs: &[(&'a [u8], &'a [u8])],
    feat_pat: &FeaturePattern,
    feat_dist: &[f64],
) -> Option<(&'a [u8], &'a [u8], Vec<u8>)> {
    let check_seq = |seq: &[u8], qual: &[u8], i: Option<usize>| -> Option<f64> {
        if let Some(feature) = feat_pat.get_feature(seq) {
            let p_whitelist = feat_dist[feature.index];
            let likelihood = if let Some(i) = i {
                let qv = cmp::min(qual[i] - ILLUMINA_QUAL_OFFSET, FEATURE_MAX_QV) as f64;
                let p_edit = 10.0_f64.powf(-qv / 10.0);
                p_whitelist * p_edit
            } else {
                p_whitelist
            };
            return Some(likelihood);
        }
        None
    };

    // compute likelihoods for Hamming distance 1 sequences
    #[allow(clippy::type_complexity)]
    let mut whitelist_likelihoods: HashMap<Vec<u8>, (f64, &'a [u8], &'a [u8])> = HashMap::new();
    let mut likelihood_sum = 0.0;

    // this determines the strategy for handling multiple hits, we keep the one
    //   with greatest likelihood
    let mut insert_hit = |likelihood: f64, test_seq: &[u8], seq: &'a [u8], qual: &'a [u8]| {
        match whitelist_likelihoods.entry(test_seq.to_owned()) {
            Entry::Occupied(mut entry) => {
                // if our new likelihood is greater, replace the old hit
                if likelihood > entry.get().0 {
                    let (old_likelihood, _, _) = entry.insert((likelihood, seq, qual));
                    likelihood_sum += likelihood - old_likelihood;
                }
            }
            Entry::Vacant(entry) => {
                entry.insert((likelihood, seq, qual));
                likelihood_sum += likelihood;
            }
        }
    };

    for &(seq, qual) in bcs {
        // fast-path if seq already correct, this edit is "100%"
        if let Some(likelihood) = check_seq(seq, qual, None) {
            insert_hit(likelihood, seq, seq, qual);
            continue;
        }

        let mut test_seq: Vec<u8> = seq.to_owned();
        for i in 0..test_seq.len() {
            let orig_base = test_seq[i];
            for base in NUCLEOTIDES {
                if *base != orig_base {
                    test_seq[i] = *base;
                    if let Some(likelihood) = check_seq(&test_seq, qual, Some(i)) {
                        insert_hit(likelihood, test_seq.as_slice(), seq, qual);
                    }
                }
            }
            test_seq[i] = orig_base;
        }
    }

    // find maximum likelihood
    let mut max_likelihood = -1.0;
    let mut best_whitelist_bc = Vec::new();
    let mut best_seq = None;
    let mut best_qual = None;
    for (whitelist_bc, (likelihood, seq, qual)) in whitelist_likelihoods {
        if likelihood > max_likelihood {
            max_likelihood = likelihood;
            best_whitelist_bc = whitelist_bc;
            best_seq = Some(seq);
            best_qual = Some(qual);
        }
    }

    if (max_likelihood / likelihood_sum) >= FEATURE_CONF_THRESHOLD {
        return Some((best_seq.unwrap(), best_qual.unwrap(), best_whitelist_bc));
    }
    None
}

#[derive(Eq, Copy, Clone, PartialEq)]
enum PatternType {
    Tethered,
    Untethered,
}

pub struct FeaturePattern {
    read: WhichRead,
    regex_str: String,
    regex: Regex,
    feature_type: FeatureType,
    features: HashMap<Vec<u8>, FeatureDef>,
    pattern_type: PatternType,
}

impl FeaturePattern {
    fn least_feature_index(&self) -> usize {
        self.features.values().map(|x| x.index).min().unwrap()
    }

    fn get_feature(&self, seq: &[u8]) -> Option<&FeatureDef> {
        self.features.get(seq)
    }

    fn is_tethered(&self) -> bool {
        self.pattern_type == PatternType::Tethered
    }
}

impl FeaturePattern {
    fn insert(&mut self, bc_seq: Vec<u8>, feature: FeatureDef) -> Result<()> {
        if self.features.contains_key(&bc_seq) {
            bail!(
                "Found two feature definitions with the same read ('{:?}'), pattern ('{}') and \
                 barcode sequence ('{}'). This combination of values must be unique for each row \
                 in the feature definition file.",
                self.read,
                feature.pattern,
                String::from_utf8_lossy(&bc_seq)
            );
        }

        self.features.insert(bc_seq, feature);
        Ok(())
    }
}

pub struct FeatureExtractor {
    reference: Arc<FeatureReference>,
    patterns: HashMap<(FeatureType, WhichRead), (RegexSet, Vec<FeaturePattern>)>,
    feature_dist: Option<Vec<f64>>,
}

#[derive(Serialize, Deserialize)]
pub struct FeatureData {
    #[serde(with = "serde_bytes")]
    pub barcode: Vec<u8>,
    #[serde(with = "serde_bytes")]
    pub corrected_barcode: Option<Vec<u8>>,
    #[serde(with = "serde_bytes")]
    pub qual: Vec<u8>,
    pub ids: Vec<(usize, String)>,
}

impl FeatureExtractor {
    pub fn new(
        feature_ref: Arc<FeatureReference>,
        use_feature_types: Option<&HashSet<FeatureType>>,
        feature_dist: Option<Vec<f64>>,
    ) -> Result<FeatureExtractor> {
        let mut patterns: HashMap<(FeatureType, WhichRead, String), FeaturePattern> =
            HashMap::new();

        let mut bare_patterns = HashMap::new();
        for fd in &feature_ref.feature_defs {
            if fd.feature_type == FeatureType::Gene
                || use_feature_types.is_some_and(|m| !m.contains(&fd.feature_type))
            {
                continue;
            }

            FeatureExtractor::validate_sequence(fd.sequence.as_bytes())?;

            // we prefer longer patterns, but bundle all bare barcodes of the
            //   same length together
            if fd.pattern == "(BC)" {
                bare_patterns
                    .entry((fd.feature_type, fd.read, fd.sequence.len()))
                    .or_insert_with(Vec::new)
                    .push(fd);
                continue;
            }

            let (regex, regex_str, _) =
                FeatureExtractor::compile_pattern(&fd.pattern, fd.sequence.len())?;
            FeatureExtractor::insert(
                &mut patterns,
                fd.feature_type,
                fd.read,
                regex_str,
                regex,
                fd.sequence.as_bytes().to_owned(),
                fd.clone(),
                PatternType::Tethered,
            )?;
        }

        for ((feature_type, read, _), fds) in bare_patterns {
            let (regex, regex_str, _) = FeatureExtractor::compile_bare_patterns(fds.as_slice())?;
            for &fd in &fds {
                FeatureExtractor::insert(
                    &mut patterns,
                    feature_type,
                    read,
                    regex_str.clone(),
                    regex.clone(),
                    fd.sequence.as_bytes().to_owned(),
                    fd.clone(),
                    PatternType::Untethered,
                )?;
            }
        }

        let patterns = patterns
            .into_iter()
            .fold(HashMap::new(), |mut acc: HashMap<_, Vec<_>>, (_, pat)| {
                acc.entry((pat.feature_type, pat.read))
                    .or_default()
                    .push(pat);
                acc
            })
            .into_iter()
            .map(|(r, pats)| {
                let regex_set: Result<RegexSet> = RegexSet::new(pats.iter().map(|x| &x.regex_str))
                    .map_err(|e| match e {
                        regex::Error::CompiledTooBig(_) => anyhow!(REGEX_TOO_BIG_ERR),
                        _ => anyhow::Error::from(e),
                    });
                anyhow::Ok((r, (regex_set?, pats)))
            })
            .try_collect()?;

        Ok(FeatureExtractor {
            reference: feature_ref,
            patterns,
            feature_dist,
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn insert(
        pats: &mut HashMap<(FeatureType, WhichRead, String), FeaturePattern>,
        feature_type: FeatureType,
        read: WhichRead,
        regex_str: String,
        regex: Regex,
        bc_seq: Vec<u8>,
        feature: FeatureDef,
        pattern_type: PatternType,
    ) -> Result<()> {
        pats.entry((feature_type, read, regex_str.clone()))
            .or_insert(FeaturePattern {
                read,
                regex_str,
                regex,
                feature_type,
                features: HashMap::new(),
                pattern_type,
            })
            .insert(bc_seq, feature)
    }

    fn compile_bare_patterns(pats: &[&FeatureDef]) -> Result<(Regex, String, String)> {
        // here we support bare "(BC)" definitions, using a group that matches N versions of the
        // barcode, one each for each position to allow for hamming-1 errors
        let parts: Vec<_> = pats
            .iter()
            .flat_map(|&fd| {
                let length = fd.sequence.len();
                (0..length)
                    .map(move |i| format!("{}.{}", &fd.sequence[0..i], &fd.sequence[(i + 1)..]))
            })
            .collect();
        let pat = format!("({})", parts.join("|"));
        Ok((Regex::new(&pat)?, pat, "(BC)".to_string()))
    }

    fn compile_pattern(orig_pat: &str, length: usize) -> Result<(Regex, String, String)> {
        // bare patterns are handled elsewhere
        assert_ne!(orig_pat, "(BC)");

        // We previously only supported the regex ^/$ markers for start and end,
        // but in order to be more human-friendly, we now also support 5p and 3p too
        // replace them here to make a valid regex
        let re5p = Regex::new("^5[Pp]?[-_]?").unwrap();
        let pat = re5p.replace(orig_pat, "^");

        let re3p = Regex::new("[-_]?3[Pp]?$").unwrap();
        let pat = re3p.replace(&pat, "$");

        let bc_regex = Regex::new(r"\(BC\)").unwrap();
        let check_pattern = bc_regex.replace(&pat, "");

        let valid_regex = Regex::new(r"^\^?[ACGTN]*\$?$").unwrap();
        if !valid_regex.is_match(&check_pattern) || orig_pat.matches("(BC)").count() != 1 {
            bail!(
                "Invalid pattern: '{check_pattern}'. The pattern must optionally start with '5P', \
                 optionally end with '3P', contain exactly one instance of the string '(BC)' and \
                 otherwise contain only the characters A, C, G, T, and N."
            );
        }

        // TODO: should be [ACGTN]?
        let pat = Regex::new("N").unwrap().replace_all(&pat, ".");

        let bc_match = format!("(.{{{length},{length}}})");

        let pat = Regex::new(r"\(BC\)")
            .unwrap()
            .replace::<&str>(&pat, &bc_match)
            .to_string();

        Ok((Regex::new(&pat)?, pat, orig_pat.to_string()))
    }

    fn validate_sequence(seq: &[u8]) -> Result<()> {
        let re = regex::bytes::Regex::new("^[ACGTN]+$").unwrap();

        if re.is_match(seq) {
            Ok(())
        } else {
            bail!(
                "Invalid sequence: '{}'. The only allowed characters are A, C, G, T, and N.",
                std::str::from_utf8(seq)?
            )
        }
    }

    /// Extract the feature bc sequence from the read.
    pub fn match_read(&self, read: &RnaRead) -> Option<FeatureData> {
        // this corresponds to the
        // lib/python/cellranger/rna/feature_ref.py:FeatureExtractor._extract_from_seq method
        // the behavior there is to take the longest whitelist hits, if any
        // or just the longest hit
        // or nothing
        // This is less than ideal, but how it was done.

        let mut pattern_matches = vec![];
        let mut whitelist_matches = vec![];
        let mut barcode_captures = vec![];

        for (&(feature_type, which_read), (regset, patterns)) in &self.patterns {
            if read.library_type != feature_type.into() {
                continue;
            }

            let seq = read.read.get(which_read, ReadPart::Seq).unwrap();
            let qual = read.read.get(which_read, ReadPart::Qual).unwrap();

            let s = std::str::from_utf8(seq).unwrap();

            // these are guaranteed to always be in ascending order of index
            for pat in regset.matches(s).iter().map(|i| &patterns[i]) {
                barcode_captures.clear();
                let mut offset = 0;
                // captures_iter is non-overlapping, so do this instead
                while let Some(cap) = pat.regex.captures(&s[offset..]) {
                    let mat = cap.get(1).unwrap();
                    let bc = mat.as_str().as_bytes();
                    let start = mat.start() + offset;
                    let end = mat.end() + offset;
                    let q = &qual[start..end];
                    barcode_captures.push((bc, q));
                    // pushing the next iteration to start just after this match
                    //   gets us around the non-overlapping limitation
                    offset += mat.start() + 1;
                    // tethered features get only the first chance
                    if pat.is_tethered() {
                        break;
                    }
                }
                if let Some((bc, q, corrected_bc, feat_idx)) =
                    self.find_closest(pat, barcode_captures.as_slice())
                {
                    whitelist_matches.push((bc, q, corrected_bc, feat_idx));
                } else {
                    for (bc, q) in &barcode_captures {
                        pattern_matches.push((*bc, *q, pat.least_feature_index()));
                    }
                }
            }
        }

        if !whitelist_matches.is_empty() {
            let (barcode, qual, corrected_barcode, _feat_idx): (&[u8], &[u8], _, _) =
                whitelist_matches
                    .iter()
                    .max_by_key(|x| (x.2.len(), Reverse(x.3)))
                    .unwrap()
                    .clone();
            let ids = whitelist_matches
                .into_iter()
                .map(|x| (x.3, self.reference.feature_defs[x.3].id.clone()))
                .collect();
            return Some(FeatureData {
                barcode: barcode.to_owned(),
                qual: qual.to_owned(),
                corrected_barcode: Some(corrected_barcode),
                ids,
            });
        }

        if !pattern_matches.is_empty() {
            let (barcode, qual, _feat_idx): (&[u8], &[u8], _) = pattern_matches
                .into_iter()
                .max_by_key(|x| (x.0.len(), Reverse(x.2)))
                .unwrap();
            return Some(FeatureData {
                barcode: barcode.to_owned(),
                qual: qual.to_owned(),
                corrected_barcode: None,
                ids: vec![],
            });
        }

        None
    }

    #[allow(clippy::type_complexity)]
    fn find_closest<'a>(
        &self,
        pat: &FeaturePattern,
        bcs: &[(&'a [u8], &'a [u8])],
    ) -> Option<(&'a [u8], &'a [u8], Vec<u8>, usize)> {
        if bcs.is_empty() {
            return None;
        }
        // if we have but 1 potential hit, fast-path check if it's in the list
        if bcs.len() == 1 {
            let (bc, q) = bcs[0];
            if let Some(bc_def) = pat.features.get(bc) {
                return Some((bc, q, Vec::from(bc), bc_def.index));
            }
        }
        if let Some(feat_dist) = self.feature_dist.as_ref() {
            if let Some((bc, q, hit)) = correct_feature_barcode(bcs, pat, feat_dist) {
                let bc_def = &pat.features[&hit];
                return Some((bc, q, hit, bc_def.index));
            }
        }
        None
    }
}

/// Return true if a library type requires a feature reference
pub fn library_type_requires_feature_ref(library_type: &str) -> bool {
    !LIBRARY_TYPES_WITHOUT_FEATURES
        .iter()
        .any(|x| *x == library_type)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reference::feature_checker::compute_feature_dist;
    use crate::reference::reference_info::ReferenceInfo;
    use crate::types::{FeatureBarcodeType, GenomeName, LibraryType};
    use arrayvec::ArrayVec;
    use barcode::BarcodeConstruct::GelBeadOnly;
    use barcode::BarcodeSegmentState::NotChecked;
    use barcode::SegmentedBarcode;
    use fastq_set::read_pair::{ReadPair, RpRange, WhichRead};
    use fastq_set::Record;
    use std::fs::File;
    use std::io::{BufReader, Cursor, Write};
    use transcriptome::Transcriptome;
    use umi::Umi;

    // helper function for checking corrections
    fn correct_feature_barcode(
        fext: &FeatureExtractor,
        seq: &[u8],
        qual: &[u8],
        ftype: FeatureType,
    ) -> Option<Vec<u8>> {
        // implement a private type to hold the data for ReadPair
        struct Rec<'a> {
            seq: &'a [u8],
            qual: &'a [u8],
        }
        impl<'a> Record for Rec<'a> {
            fn seq(&self) -> &[u8] {
                self.seq
            }
            fn head(&self) -> &[u8] {
                &b"h"[..]
            }
            fn qual(&self) -> &[u8] {
                self.qual
            }
            fn write<W: Write>(&self, _: &mut W) -> std::io::Result<usize> {
                unimplemented!()
            }
        }
        let read = RnaRead {
            read: ReadPair::new([Some(Rec { seq, qual }), Some(Rec { seq, qual }), None, None]),
            barcode: SegmentedBarcode::gel_bead_only(0, b"A", NotChecked),
            umi: Umi::new(b"A"),
            bc_range: GelBeadOnly(RpRange::new(WhichRead::R1, 0, None)),
            umi_parts: ArrayVec::new(),
            r1_range: RpRange::new(WhichRead::R1, 0, Some(seq.len())),
            r2_range: Some(RpRange::new(WhichRead::R2, 0, Some(seq.len()))),
            library_type: LibraryType::from(ftype),
            chunk_id: 0,
        };
        fext.match_read(&read).and_then(|x| x.corrected_barcode)
    }

    #[test]
    fn test_load_feature_ref() -> Result<()> {
        let fdf_path = "test/feature/citeseq.csv";
        let rdr = BufReader::new(File::open(fdf_path).unwrap());
        let fref = FeatureReference::new(
            &ReferenceInfo::default(),
            &Transcriptome::dummy(),
            Some(rdr),
            None,
            None,
            None,
            None,
        )
        .unwrap();
        let _fextr = FeatureExtractor::new(
            Arc::new(fref),
            Some(&[FeatureType::Barcode(FeatureBarcodeType::Antibody)].into()),
            None,
        )
        .unwrap();
        Ok(())
    }

    #[test]
    fn test_bad_feature_ref() -> Result<()> {
        let fdf_path = "test/feature/CRISPR_lib.v5.500.csv";
        let rdr = BufReader::new(File::open(fdf_path)?);
        let fref = FeatureReference::new(
            &ReferenceInfo::default(),
            &Transcriptome::dummy(),
            Some(rdr),
            None,
            None,
            None,
            None,
        )?;
        let fextr = FeatureExtractor::new(
            Arc::new(fref),
            Some(&[FeatureType::Barcode(FeatureBarcodeType::Crispr)].into()),
            None,
        );
        match fextr {
            Err(e) => assert_eq!(e.to_string(), REGEX_TOO_BIG_ERR),
            Ok(_) => panic!("Expected Error not observed"),
        }
        Ok(())
    }

    #[test]
    fn test_compile_pattern() {
        let fd = FeatureDef {
            index: 0,
            id: "ID1".to_string(),
            name: "Name1".to_string(),
            genome: GenomeName::default(),
            sequence: "ACGT".to_string(),
            pattern: "(BC)".to_string(),
            read: WhichRead::R1,
            feature_type: FeatureType::Barcode(FeatureBarcodeType::Antibody),
            tags: HashMap::new(),
        };
        let r = FeatureExtractor::compile_bare_patterns(&[&fd]).unwrap();
        assert_eq!("(.CGT|A.GT|AC.T|ACG.)", r.1);

        let r = FeatureExtractor::compile_pattern("AGTCN(BC)TTT", 5).unwrap();
        assert_eq!("AGTC.(.{5,5})TTT", r.1);

        let r = FeatureExtractor::compile_pattern("5PAGTCN(BC)TTT", 5).unwrap();
        assert_eq!("^AGTC.(.{5,5})TTT", r.1);

        let r = FeatureExtractor::compile_pattern("5PAGTCN(BC)TTT-3p", 5).unwrap();
        assert_eq!("^AGTC.(.{5,5})TTT$", r.1);

        let r = FeatureExtractor::compile_pattern("5P-AGTCN(BC)TTT3p", 5).unwrap();
        assert_eq!("^AGTC.(.{5,5})TTT$", r.1);

        let r = FeatureExtractor::compile_pattern("^AGTCN(BC)TTT$", 5).unwrap();
        assert_eq!("^AGTC.(.{5,5})TTT$", r.1);

        let r = FeatureExtractor::compile_pattern("^AGTCN(BC)TTT3$", 5);
        assert!(r.is_err());

        let r = FeatureExtractor::compile_pattern("5PAGTCN(BCTTT", 5);
        assert!(r.is_err());

        let r = FeatureExtractor::compile_pattern("5PAGTCNTTT", 5);
        assert!(r.is_err());

        let r = FeatureExtractor::compile_pattern("5PAGT(BC)CNTQTT", 5);
        assert!(r.is_err());

        let r = FeatureExtractor::compile_pattern("3PAGT(BC)CNTATT", 5);
        assert!(r.is_err());

        let r = FeatureExtractor::compile_pattern("AGT(BC)CNTATT5P", 5);
        assert!(r.is_err());

        let r = FeatureExtractor::compile_pattern("AGT(BC)CNTATT^", 5);
        assert!(r.is_err());
    }

    #[test]
    fn test_correct_bare_feature() {
        let fdf_csv = r#"
id,name,read,pattern,sequence,feature_type
ID1,Name1,R1,(BC),ACGT,Antibody Capture
ID2,Name2,R1,(BC),ACCT,Antibody Capture
ID3,Name3,R1,(BC),TTTT,Antibody Capture
"#;
        let ref_info = ReferenceInfo::default();
        let txome = Transcriptome::dummy();
        let fref = FeatureReference::new(
            &ref_info,
            &txome,
            Some(Cursor::new(fdf_csv.as_bytes())),
            None,
            None,
            None,
            None,
        )
        .unwrap();
        let fdist = compute_feature_dist(vec![1i64, 10, 10], &fref).unwrap();
        let fext = FeatureExtractor::new(Arc::new(fref), None, Some(fdist)).unwrap();
        let correct_feature_barcode =
            |seq, qual, ftype| correct_feature_barcode(&fext, seq, qual, ftype);

        // matches two patterns, but cannot choose b/c of 97.5% threshold
        let r = correct_feature_barcode(
            b"ACTT",
            b"IIII",
            FeatureType::Barcode(FeatureBarcodeType::Antibody),
        );
        assert_eq!(r, None);

        // also _perfectly_ matches two patterns, but still no dice
        let r = correct_feature_barcode(
            b"ACCTTTT",
            b"IIIIIII",
            FeatureType::Barcode(FeatureBarcodeType::Antibody),
        );
        assert_eq!(r, None);

        let r = correct_feature_barcode(
            b"ACGT",
            b"IIII",
            FeatureType::Barcode(FeatureBarcodeType::Antibody),
        );
        assert_eq!(r, Some(b"ACGT".to_vec()));

        let r = correct_feature_barcode(
            b"ACCT",
            b"IIII",
            FeatureType::Barcode(FeatureBarcodeType::Antibody),
        );
        assert_eq!(r, Some(b"ACCT".to_vec()));

        let r = correct_feature_barcode(
            b"TTTT",
            b"IIII",
            FeatureType::Barcode(FeatureBarcodeType::Antibody),
        );
        assert_eq!(r, Some(b"TTTT".to_vec()));

        let r = correct_feature_barcode(
            b"TTTA",
            b"IIII",
            FeatureType::Barcode(FeatureBarcodeType::Antibody),
        );
        assert_eq!(r, Some(b"TTTT".to_vec()));
    }

    #[test]
    fn test_load_feature_dist() {
        let fdf_csv = r#"
id,name,read,pattern,sequence,feature_type
ID1,Name1,R1,^(BC),AA,Antibody Capture
ID2,Name1,R1,^(BC),AT,Antibody Capture
ID3,Name1,R1,^(BC),TA,CRISPR Guide Capture
ID4,Name1,R1,^(BC),TT,Custom
ID5,Name1,R1,^(BC),TT,Custom
ID6,Name1,R1,^(BC),TT,Custom
"#;
        let ref_info = ReferenceInfo::default();
        let txome = Transcriptome::dummy();
        let fref = FeatureReference::new(
            &ref_info,
            &txome,
            Some(Cursor::new(fdf_csv.as_bytes())),
            None,
            None,
            None,
            None,
        )
        .unwrap();

        let x = compute_feature_dist(vec![1i64, 1, 0, 9, 1, 0], &fref).unwrap();

        assert_eq!(x, vec![0.5, 0.5, 0.0, 9.0 / 10.0, 1.0 / 10.0, 0.0]);
    }

    #[test]
    fn test_correct_feature() {
        let gi_path = "test/feature/gene_index.tab";
        let _gi_file = File::open(gi_path).unwrap();

        let fdf_csv = r#"
id,name,read,pattern,sequence,feature_type
ID1,N,R1,^(BC),AAAA,Antibody Capture
ID2,N,R1,^(BC),CCCC,Antibody Capture
ID3,N,R1,^(BC),GGGG,Antibody Capture
ID4,N,R1,^(BC),TTTT,CRISPR Guide Capture
ID5,N,R1,^(BC),TTTA,CRISPR Guide Capture
ID6,N,R1,^(BC),AAAT,Custom
ID7,N,R1,^(BC),AATA,Custom
ID8,N,R1,^(BC),ATAA,Custom
ID9,N,R1,^(BC),TAAA,Custom
"#;

        let ref_info = ReferenceInfo::default();
        let txome = Transcriptome::dummy();
        let fref = FeatureReference::new(
            &ref_info,
            &txome,
            Some(Cursor::new(fdf_csv.as_bytes())),
            None,
            None,
            None,
            None,
        )
        .unwrap();
        let fdist = compute_feature_dist(vec![0i64, 10, 1, 10, 10, 10, 10, 10, 10], &fref).unwrap();
        let fext = FeatureExtractor::new(Arc::new(fref), None, Some(fdist)).unwrap();
        let correct_feature_barcode =
            |seq, qual, ftype| correct_feature_barcode(&fext, seq, qual, ftype);

        let (seq, qual, ftype) = (
            b"AAAT",
            b"IIII",
            FeatureType::Barcode(FeatureBarcodeType::Antibody),
        );
        let r = correct_feature_barcode(seq, qual, ftype);
        assert_eq!(r, None);

        let (seq, qual, ftype) = (
            b"CGCC",
            b"IIII",
            FeatureType::Barcode(FeatureBarcodeType::Antibody),
        );
        let r = correct_feature_barcode(seq, qual, ftype);
        assert_eq!(r, Some(b"CCCC".to_vec()));

        let (seq, qual, ftype) = (
            b"TTTA",
            b"IIII",
            FeatureType::Barcode(FeatureBarcodeType::Antibody),
        );
        let r = correct_feature_barcode(seq, qual, ftype);
        assert_eq!(r, None);

        let (seq, qual, ftype) = (
            b"TTTC",
            b"IIII",
            FeatureType::Barcode(FeatureBarcodeType::Crispr),
        );
        let r = correct_feature_barcode(seq, qual, ftype);
        assert_eq!(r, None);

        let (seq, qual, ftype) = (
            b"AAAA",
            b"III!",
            FeatureType::Barcode(FeatureBarcodeType::Custom),
        );
        let r = correct_feature_barcode(seq, qual, ftype);
        assert_eq!(r, Some(b"AAAT".to_vec()));

        let (seq, qual, ftype) = (
            b"AAAA",
            b"I!II",
            FeatureType::Barcode(FeatureBarcodeType::Custom),
        );
        let r = correct_feature_barcode(seq, qual, ftype);
        assert_eq!(r, Some(b"ATAA".to_vec()));

        let (seq, qual, ftype) = (
            b"AAAA",
            b"IIII",
            FeatureType::Barcode(FeatureBarcodeType::Custom),
        );
        let r = correct_feature_barcode(seq, qual, ftype);
        assert_eq!(r, None);
    }
}
