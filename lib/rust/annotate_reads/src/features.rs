use csv;
use rust_htslib::bam;
use serde_json;
use std::collections::{HashMap, HashSet};
use std::cmp;
use std::io::{BufReader, Read};
use std::str;
use std::string::String;

use utils;

pub const RAW_FEATURE_BARCODE_TAG: &'static [u8]   = b"fr";
pub const PROC_FEATURE_BARCODE_TAG: &'static [u8]  = b"fb";
pub const FEATURE_BARCODE_QUAL_TAG: &'static [u8]  = b"fq";
pub const FEATURE_IDS_TAG: &'static [u8]           = b"fx";

const FEATURE_CONF_THRESHOLD: f64 = 0.975;
const FEATURE_MAX_QV: u8 = 33; // cap observed QVs at this value

pub const LIBRARY_TYPES_WITHOUT_FEATURES: &'static [&'static str] = &["Gene Expression"];

#[derive(Debug)]
/// Processes raw feature barcode data
pub struct FeatureChecker {
    pub feature_reference: FeatureReference,
    pub feature_dist: Vec<f64>,
    pub feat_type_idx: usize,
}

impl FeatureChecker {
    pub fn new<R: Read>(fref_file: R, gene_index_file: R,
                        fdist_file: R, library_type: &str) -> FeatureChecker {
        let fref = FeatureReference::new(fref_file, gene_index_file);
        let fdist = load_feature_dist(fdist_file, &fref);

        let ftype_idx = fref.feature_types.iter().position(|x| x == library_type)
            .expect(&format!("Couldn't find the library type '{}' as a feature type in the feature definition file.", library_type));

        FeatureChecker {
            feature_reference: fref,
            feature_dist: fdist,
            feat_type_idx: ftype_idx,
        }
    }

    /// Extract and potentially correct the feature barcode
    pub fn process_feature_data(&self, tags: &HashMap<String, String>) -> Option<FeatureData> {
        let raw_seq = tags.get(str::from_utf8(RAW_FEATURE_BARCODE_TAG).unwrap())
            .map(|x| x.as_bytes());
        let qual = tags.get(str::from_utf8(FEATURE_BARCODE_QUAL_TAG).unwrap())
                            .map(|x| x.as_bytes());

        // Try to correct the raw sequence if there is no processed sequence
        let (corrected_seq, was_corrected) =
            tags.get(str::from_utf8(PROC_FEATURE_BARCODE_TAG).unwrap())
            .map_or_else(||
                         // The processed seq doesn't exist; try to correct the raw seq
                         match (raw_seq.as_ref(), qual.as_ref()) {
                             (Some(seq), Some(qual)) =>
                                 correct_feature_barcode(seq, qual,
                                                         &self.feature_reference,
                                                         &self.feature_dist,
                                                         self.feat_type_idx)
                                 .map_or(
                                     (None, false),
                                     |new_seq| (Some(new_seq), true)),
                             // No sequence to correct
                             _ =>  (None, false),
                         },
                         // The processed seq exists, pass it through uncorrected
                         |x| (Some(x.as_bytes().to_owned()), false),
            );

        // Take the given feature(s), or try to find the feature for the corrected sequence
        let ids = tags.get(str::from_utf8(FEATURE_IDS_TAG).unwrap())
            .or(corrected_seq.as_ref()
                .and_then(|seq| self.feature_reference.get_feature_by_barcode(self.feat_type_idx, seq)
                          .map(|feat| &feat.id)));

        match (raw_seq, qual) {
            (Some(raw_seq), Some(qual)) => Some(FeatureData{
                raw_seq: raw_seq.to_owned(),
                corrected_seq: corrected_seq,
                qual: qual.to_owned(),
                ids: ids.cloned(),
                was_corrected: was_corrected,
            }),
            _ => None,
        }
    }
}

#[derive(Debug)]
/// Feature barcode data for a sequence read
pub struct FeatureData {
    pub raw_seq: Vec<u8>,
    pub corrected_seq: Option<Vec<u8>>,
    pub qual: Vec<u8>,
    pub ids: Option<String>, // semicolon-delimited list of feature IDs
    pub was_corrected: bool, // this feature barcode was corrected
}

impl FeatureData {
    /// Attach feature-related BAM tags to a BAM record.
    pub fn attach_tags(&self, record: &mut bam::Record) {
        record.push_aux(&RAW_FEATURE_BARCODE_TAG,
                        &bam::record::Aux::String(&self.raw_seq));

        if let Some(ref seq) = self.corrected_seq {
            record.push_aux(&PROC_FEATURE_BARCODE_TAG,
                            &bam::record::Aux::String(seq));
        }

        record.push_aux(&FEATURE_BARCODE_QUAL_TAG,
                        &bam::record::Aux::String(&self.qual));

        if let Some(ref feature_ids) = self.ids {
            record.push_aux(&FEATURE_IDS_TAG,
                            &bam::record::Aux::String(feature_ids.as_bytes()));
        }
    }
}

/// Gene index TSV row
#[derive(Debug, Deserialize)]
struct GeneIndexRow {
    pub transcript_id: String,
    pub gene_id: String,
    pub gene_name: String,
    pub transcript_length: i64,
}

/// Feature definition CSV row
#[derive(Debug, Deserialize)]
struct FeatureRow {
    pub id: String,
    pub sequence: String,
    pub feature_type: String,
}

/// Feature reference entry
#[derive(Debug)]
pub struct FeatureDef {
    pub index: usize,
    pub id: String,
    pub sequence: Vec<u8>,
    pub feature_type_idx: usize,
}

#[derive(Debug)]
pub struct FeatureReference {
    pub feature_defs: Vec<FeatureDef>,
    pub feature_types: Vec<String>,
    pub feature_maps: Vec<HashMap<Vec<u8>, usize>>,
}

impl FeatureReference {
    pub fn new<R: Read, R2: Read>(csv_stream: R, gene_index_stream: R2) -> FeatureReference {
        let reader = BufReader::new(csv_stream);
        let mut csv_reader = csv::Reader::from_reader(reader);

        let gene_reader = BufReader::new(gene_index_stream);
        let mut gene_csv_reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_reader(gene_reader);

        let mut type_map = HashMap::new();
        let mut type_vec = Vec::new();

        let mut fdefs = Vec::new();
        let mut fmaps = Vec::new();

        // Create gene expression features
        let mut seen_genes = HashSet::new();

        type_map.insert("Gene Expression".to_owned(), 0);
        type_vec.push("Gene Expression".to_owned());
        fmaps.push(HashMap::new());

        for record in gene_csv_reader.deserialize() {
            let row: GeneIndexRow = record.expect("Failed to parse gene index TSV row");

            // The rows are per-transcript but we only care about genes here
            if seen_genes.contains(&row.gene_id) {
                continue;
            }
            let num_fdefs = fdefs.len();
            fdefs.push(FeatureDef {
                index: num_fdefs,
                id: row.gene_id.clone(),
                sequence: vec![],
                feature_type_idx: type_map.len() - 1,
            });
            seen_genes.insert(row.gene_id);
        }

        // Create feature barcode (fBC) features
        for record in csv_reader.deserialize() {
            let frow: FeatureRow = record.expect("Failed to parse feature definition row");

            // Map feature types to integers
            if !type_map.contains_key(&frow.feature_type) {
                let num_types = type_map.len();
                type_map.insert(frow.feature_type.clone(), num_types);
                type_vec.push(frow.feature_type.clone());
                fmaps.push(HashMap::new());
            }

            let num_fdefs = fdefs.len();
            let ftype = *type_map.get(&frow.feature_type).unwrap();
            fdefs.push(FeatureDef {
                index: num_fdefs,
                id: frow.id,
                sequence: frow.sequence.as_bytes().to_owned(),
                feature_type_idx: ftype,
            });

            fmaps[ftype].insert(frow.sequence.as_bytes().to_owned(), num_fdefs);
        }

        FeatureReference {
            feature_defs: fdefs,
            feature_types: type_vec,
            feature_maps: fmaps,
        }
    }

    /// Get a feature by its barcode sequence
    /// * `barcode` - ASCII nucleotide sequence
    fn get_feature_by_barcode(&self, feat_type_idx: usize, barcode: &[u8])
    -> Option<&FeatureDef> {
        self.feature_maps[feat_type_idx].get(barcode)
            .map(|i| &self.feature_defs[*i])
    }
}

/// Load a JSON file containing an array of integer counts, one per feature.
/// Return the proportions, normalized within each feature type.
fn load_feature_dist<R: Read>(counts_json: R, feat_ref: &FeatureReference) -> Vec<f64> {
    let reader = BufReader::new(counts_json);

    let counts: Vec<i64> = serde_json::from_reader(reader)
        .expect("Failed to parse feature counts JSON file");

    let feature_types: Vec<usize> = feat_ref.feature_defs.iter()
        .map(|fd| fd.feature_type_idx).collect();
    if counts.len() != feature_types.len() {
        panic!("Mismatch between number of features in counts JSON ({}) and in feature reference ({})",
               counts.len(), feature_types.len());
    }

    let mut feature_type_sums = vec![0i64; feat_ref.feature_types.len()];

    for (count, feature_type) in counts.iter().zip(feature_types.iter()) {
        feature_type_sums[*feature_type] += *count;
    }

    let mut proportions = vec![0f64; counts.len()];
    for (i, (count, feature_type)) in counts.iter().zip(feature_types.iter()).enumerate() {
        let sum = feature_type_sums[*feature_type];
        if sum > 0 {
            proportions[i] = *count as f64 / sum as f64;
        }
    }

    return proportions;
}

/// Correct a feature barcode against a feature reference.
/// * `seq` - ASCII nucleotide sequence
/// * `qual` - ASCII ILMN Phred qual string
/// * `feat_dist` - Feature barcode distrubution as proportions (within-feature-type)
/// * `feat_type_idx` - Integer indicating the feature type
fn correct_feature_barcode(seq: &[u8], qual: &[u8],
                           feat_ref: &FeatureReference,
                           feat_dist: &[f64],
                           feat_type_idx: usize) -> Option<Vec<u8>> {
    let mut test_seq = seq.to_owned();
    let mut whitelist_likelihoods = HashMap::new();
    let mut likelihood_sum = 0.0;

    // compute likelihoods for Hamming distance 1 sequences
    for i in 0..test_seq.len() {
        let orig_base = test_seq[i];
        for base in utils::NUCLEOTIDES {
            if *base != orig_base {
                test_seq[i] = *base;

                if let Some(ref feature) = feat_ref.get_feature_by_barcode(feat_type_idx, &test_seq) {
                    let p_whitelist = feat_dist[feature.index];
                    let qv = cmp::min(qual[i] - utils::ILLUMINA_QUAL_OFFSET, FEATURE_MAX_QV) as f64;
                    let p_edit = 10.0_f64.powf(-qv / 10.0);
                    let likelihood = p_whitelist * p_edit;
                    whitelist_likelihoods.insert(test_seq.clone(), likelihood);
                    likelihood_sum += likelihood;
                }
            }
        }
        test_seq[i] = orig_base;
    }

    // find maximum likelihood
    let mut max_likelihood = -1.0;
    let mut best_whitelist_bc = Vec::new();
    for (whitelist_bc, likelihood) in whitelist_likelihoods {
        if likelihood > max_likelihood {
            max_likelihood = likelihood;
            best_whitelist_bc = whitelist_bc;
        }
    }

    match max_likelihood / likelihood_sum >= FEATURE_CONF_THRESHOLD {
        true => Some(best_whitelist_bc),
        false => None,
    }
}

/// Return true if a library type requires a feature reference
pub fn library_type_requires_feature_ref(library_type: &str) -> bool{
    LIBRARY_TYPES_WITHOUT_FEATURES.iter()
        .position(|x| *x == library_type).is_none()
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::io::Cursor;

    #[test]
    fn test_load_feature_ref() {
        let fdf_path = "test/feature/citeseq.csv";
        let fdf_file = File::open(fdf_path).unwrap();
        let gi_path = "test/feature/gene_index.tab";
        let gi_file = File::open(gi_path).unwrap();
        let _fref = FeatureReference::new(fdf_file, gi_file);
    }

    #[test]
    fn test_load_feature_dist() {
        let gi_path = "test/feature/gene_index.tab";
        let gi_file = File::open(gi_path).unwrap();

        let fdf_csv = r#"
id,name,read,pattern,sequence,feature_type
ID1,Name1,R1,^(BC),AA,Type1
ID2,Name1,R1,^(BC),AT,Type1
ID3,Name1,R1,^(BC),TA,Type2
ID4,Name1,R1,^(BC),TT,Type3
ID5,Name1,R1,^(BC),TT,Type3
ID6,Name1,R1,^(BC),TT,Type3
"#;
        let fref = FeatureReference::new(Cursor::new(fdf_csv.as_bytes()), gi_file);

        let count_json = "[0, 0, 1, 1, 0, 9, 1, 0]";

        let x = load_feature_dist(Cursor::new(count_json.as_bytes()), &fref);

        assert_eq!(x, vec![0.0, 0.0, 0.5, 0.5, 0.0, 9.0/10.0, 1.0/10.0, 0.0]);
    }

    #[test]
    fn test_correct_feature() {
        let gi_path = "test/feature/gene_index.tab";
        let gi_file = File::open(gi_path).unwrap();

        let fdf_csv = r#"
id,name,read,pattern,sequence,feature_type
ID1,N,R1,^(BC),AAAA,type1
ID2,N,R1,^(BC),CCCC,type1
ID3,N,R1,^(BC),GGGG,type1
ID4,N,R1,^(BC),TTTT,type2
ID5,N,R1,^(BC),TTTA,type2
ID6,N,R1,^(BC),AAAT,type3
ID7,N,R1,^(BC),AATA,type3
ID8,N,R1,^(BC),ATAA,type3
ID9,N,R1,^(BC),TAAA,type3
"#;
        let fref = FeatureReference::new(Cursor::new(fdf_csv.as_bytes()), gi_file);
        let s = "[0, 0, 0, 10, 1, 10, 10, 10, 10, 10, 10]";
        let fdist = load_feature_dist(Cursor::new(s.as_bytes()), &fref);

        let (seq, qual, ftype) = (b"AAAT", b"IIII", 1);
        let r = correct_feature_barcode(seq, qual, &fref, &fdist, ftype);
        assert_eq!(r, None);

        let (seq, qual, ftype) = (b"CGCC", b"IIII", 1);
        let r = correct_feature_barcode(seq, qual, &fref, &fdist, ftype);
        assert_eq!(r, Some(b"CCCC".to_vec()));

        let (seq, qual, ftype) = (b"TTTA", b"IIII", 1);
        let r = correct_feature_barcode(seq, qual, &fref, &fdist, ftype);
        assert_eq!(r, None);

        let (seq, qual, ftype) = (b"TTTC", b"IIII", 2);
        let r = correct_feature_barcode(seq, qual, &fref, &fdist, ftype);
        assert_eq!(r, None);

        let (seq, qual, ftype) = (b"AAAA", b"III!", 3);
        let r = correct_feature_barcode(seq, qual, &fref, &fdist, ftype);
        assert_eq!(r, Some(b"AAAT".to_vec()));

        let (seq, qual, ftype) = (b"AAAA", b"I!II", 3);
        let r = correct_feature_barcode(seq, qual, &fref, &fdist, ftype);
        assert_eq!(r, Some(b"ATAA".to_vec()));

        let (seq, qual, ftype) = (b"AAAA", b"IIII", 3);
        let r = correct_feature_barcode(seq, qual, &fref, &fdist, ftype);
        assert_eq!(r, None);
    }
}
