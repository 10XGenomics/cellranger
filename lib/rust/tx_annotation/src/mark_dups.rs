//! Identify duplicate reads originating from the same molecule.

use crate::read::{AnnotationInfo, ReadAnnotations, RecordAnnotation};
use cr_types::probe_set::ProbeSetReference;
use cr_types::reference::feature_reference::FeatureReference;
use cr_types::rna_read::RnaChunk;
use cr_types::types::UmiCount;
use itertools::Itertools;
use rand::Rng;
use rand_chacha::ChaCha20Rng;
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::sync::Arc;
use umi::{UmiSeq, UmiType};

pub type Gene = usize;

/// Within each gene, correct Hamming-distance-one UMIs
fn correct_umis(umigene_counts: &HashMap<(UmiSeq, Gene), u64>) -> HashMap<(UmiSeq, Gene), UmiSeq> {
    let nucs = b"ACGT";

    let mut corrections = HashMap::new();

    for ((umi, gene), orig_count) in umigene_counts {
        let mut test_umi = *umi;

        let mut best_dest_count = *orig_count;
        let mut best_dest_umi = *umi;

        for pos in 0..umi.len() {
            // Try each nucleotide at this position
            for test_char in nucs {
                if *test_char == umi[pos] {
                    // Skip the identitical nucleotide
                    continue;
                }
                test_umi[pos] = *test_char;

                // Test for the existence of this mutated UMI
                let test_count = *umigene_counts.get(&(test_umi, *gene)).unwrap_or(&0u64);

                // If there's a 1-HD UMI w/ greater count, move to that UMI.
                // If there's a 1-HD UMI w/ equal count, move to the lexicographically larger UMI.
                if test_count > best_dest_count
                    || (test_count == best_dest_count && test_umi > best_dest_umi)
                {
                    best_dest_umi = test_umi;
                    best_dest_count = test_count;
                }
            }
            // Reset this position to the unmutated sequence
            test_umi[pos] = umi[pos];
        }
        if *umi != best_dest_umi {
            corrections.insert((*umi, *gene), best_dest_umi);
        }
    }
    corrections
}

#[derive(Clone, Copy, Serialize, Deserialize)]
pub struct DupInfo {
    pub processed_umi: UmiSeq,
    pub is_corrected: bool,
    umi_count: Option<UmiCount>,
    is_umi_count: bool,
    pub is_low_support_umi: bool,
    /// Was this read filtered because the number of reads
    /// mapping to the (UMI, gene) within this barcode is
    /// lower than the specified threshold.
    pub is_filtered_target_umi: bool,
}

impl DupInfo {
    /// Return true if this read is representative of its UMI.
    pub fn is_umi_count(&self) -> bool {
        self.is_umi_count
    }

    /// Return the UmiCount structure for this molecule.
    pub fn umi_count(&self) -> Option<UmiCount> {
        self.umi_count
    }
}

/// Mark as low support umigenes with frequency below the maximum for that UMI.
fn determine_low_support_umigenes(
    umigene_counts: &HashMap<(UmiSeq, Gene), u64>,
) -> HashSet<(UmiSeq, Gene)> {
    let mut low_support_umigenes: HashSet<(UmiSeq, Gene)> = HashSet::new();
    let mut umigene_count_vec: Vec<_> = umigene_counts
        .iter()
        .map(|x| ((x.0).0, &(x.0).1, x.1))
        .collect();
    umigene_count_vec.sort();
    for (umi, grouped) in &umigene_count_vec.iter().group_by(|x| x.0) {
        let gene_counts: Vec<_> = grouped.collect();
        if let Some((_umi, _gene, max_count)) = gene_counts.iter().max_by_key(|x| x.2) {
            let max_is_tied = gene_counts.iter().filter(|x| x.2 == *max_count).count() >= 2;
            for (_umi, gene, count) in gene_counts {
                if max_is_tied || count < max_count {
                    low_support_umigenes.insert((umi, *(*gene)));
                }
            }
        }
    }
    low_support_umigenes
}

#[derive(PartialEq, PartialOrd, Eq, Ord, Clone)]
pub struct UmiSelectKey {
    utype: UmiType,
    qname: Vec<u8>,
}

pub struct DupBuilder {
    umigene_counts: HashMap<(UmiSeq, Gene), u64>,
    umigene_min_key: HashMap<(UmiSeq, Gene), UmiSelectKey>,
}

impl DupBuilder {
    pub fn new() -> Self {
        DupBuilder {
            umigene_counts: HashMap::new(),
            umigene_min_key: HashMap::new(),
        }
    }
    pub fn observe(&mut self, annotation: &ReadAnnotations, feature_reference: &FeatureReference) {
        if annotation.umi_info.is_valid {
            // Count (raw UMI, feature) pairs to prepare for UMI correction.
            // Track (raw UMI, feature) pairs with submaximal count per (raw UMI)
            //   to prepare for marking of low-support UMIs (putative chimeras).
            if let Some(gene) = annotation.conf_mapped_feature(feature_reference) {
                let key = (annotation.umi_info.seq, gene);
                *self.umigene_counts.entry(key).or_insert(0) += 1;
                let header = annotation.read.header().to_vec(); // TODO: Use a buffer to avoid repeated allocations
                let ann_key = match annotation.is_conf_mapped_unique_txomic() {
                    true => UmiSelectKey {
                        utype: UmiType::Txomic,
                        qname: header,
                    },
                    false => UmiSelectKey {
                        utype: UmiType::NonTxomic,
                        qname: header,
                    },
                };

                let min_key = match self.umigene_min_key.remove(&key) {
                    Some(old_min) => old_min.min(ann_key),
                    None => ann_key,
                };
                self.umigene_min_key.insert(key, min_key);
            }
        }
    }
    pub fn build(
        self,
        filter_umis: bool,
        umi_correction: UmiCorrection,
        targeted_umi_min_read_count: Option<u64>,
    ) -> BarcodeDupMarker {
        BarcodeDupMarker::new(
            self.umigene_counts,
            self.umigene_min_key,
            filter_umis,
            umi_correction,
            targeted_umi_min_read_count,
        )
    }
}

impl Default for DupBuilder {
    fn default() -> Self {
        DupBuilder::new()
    }
}

/// Do duplicate marking on a single barcode's worth of data
/// 1) Correct UMI sequences
/// 2) Mark records as low-support (putatively chimeric) UMIs
/// 3) Mark records as PCR duplicates
/// 4) update ReadAnnotation with DupInfo that will be used for metrics and BAM tags
pub struct BarcodeDupMarker {
    umigene_counts: HashMap<(UmiSeq, Gene), u64>,
    low_support_umigenes: HashSet<(UmiSeq, Gene)>,
    umi_corrections: HashMap<(UmiSeq, Gene), UmiSeq>,
    umigene_min_key: HashMap<(UmiSeq, Gene), UmiSelectKey>,
    /// If this is Some(r), we filter out (UMI, genes) pairs
    /// with **less than** r reads and not include them in UMI counts.
    targeted_umi_min_read_count: Option<u64>,
}

/// Technically, this is just a bool, but a code that looks like
/// `BarcodeDupMarker::new(.., UmiCorrection::Disable)` is a lot more
/// readable than `BarcodeDupMarker::new(.., false)`
pub enum UmiCorrection {
    Enable,
    Disable,
}

impl BarcodeDupMarker {
    pub fn new(
        mut umigene_counts: HashMap<(UmiSeq, Gene), u64>,
        mut umigene_min_key: HashMap<(UmiSeq, Gene), UmiSelectKey>,
        filter_umis: bool,
        umi_correction: UmiCorrection,
        targeted_umi_min_read_count: Option<u64>,
    ) -> Self {
        // Determine which UMIs need to be corrected
        let umi_corrections: HashMap<(UmiSeq, Gene), UmiSeq> = match umi_correction {
            UmiCorrection::Enable => correct_umis(&umigene_counts),
            UmiCorrection::Disable => HashMap::new(),
        };

        let umi_correction_counts = umi_corrections
            .iter()
            .map(|(raw_key, corrected_umi)| {
                (
                    raw_key,
                    (*corrected_umi, raw_key.1),
                    umigene_counts[raw_key],
                )
            })
            .collect::<Vec<_>>();

        // Before determining low-support UMI-genes, count one read of each corrected UMI
        // to match the behaviour of Cell Ranger 3.
        for (raw_key, corrected_key, _raw_count) in &umi_correction_counts {
            // One read has been counted before determining low-support UMI-genes.
            *umigene_counts.get_mut(raw_key).unwrap() -= 1;
            *umigene_counts.get_mut(corrected_key).unwrap() += 1;
        }

        // Determine low-support UMI-genes.
        let low_support_umigenes = if filter_umis {
            determine_low_support_umigenes(&umigene_counts)
        } else {
            HashSet::new()
        };

        // After determining low-support UMI-genes, count the remaining reads of each corrected UMI.
        for (raw_key, corrected_key, raw_count) in &umi_correction_counts {
            // One read has already been counted before determining low-support UMI-genes.
            *umigene_counts.get_mut(raw_key).unwrap() -= raw_count - 1;
            *umigene_counts.get_mut(corrected_key).unwrap() += raw_count - 1;
        }

        // Which is the lowest raw UMI lexicographicaly that would be corrected to (UmiSeq, Gene)
        // and would potentially be a "UMI count".
        let mut min_raw_umis: HashMap<(UmiSeq, Gene), UmiSeq> = HashMap::new();
        for ((raw_seq, gene), corr_seq) in &umi_corrections {
            if raw_seq < corr_seq || umi_corrections.contains_key(&(*corr_seq, *gene)) {
                let min_umi = match min_raw_umis.remove(&(*corr_seq, *gene)) {
                    Some(prev_min) => prev_min.min(*raw_seq),
                    None => *raw_seq,
                };
                min_raw_umis.insert((*corr_seq, *gene), min_umi);
            }
        }
        // The min UmiSelectKey after correction
        let mut min_umi_key_corrections = HashMap::new();
        for ((corr_seq, gene), raw_seq) in min_raw_umis {
            min_umi_key_corrections
                .insert((corr_seq, gene), umigene_min_key[&(raw_seq, gene)].clone());
        }
        for (key, umi_key) in min_umi_key_corrections {
            umigene_min_key.insert(key, umi_key.clone());
        }

        BarcodeDupMarker {
            umigene_counts,
            umi_corrections,
            low_support_umigenes,
            umigene_min_key,
            targeted_umi_min_read_count,
        }
    }

    /// Process a single annotation and compute the `DupInfo`
    pub fn process(
        &mut self,
        annotation: &ReadAnnotations,
        feature_reference: &FeatureReference,
        probe_set_reference: Option<&Arc<ProbeSetReference>>,
        read_chunks: &[RnaChunk],
        barcode_subsample_rate: f64,
        rng: &mut ChaCha20Rng,
    ) -> Option<DupInfo> {
        if !annotation.umi_info.is_valid || !annotation.is_conf_mapped_to_feature() {
            return None;
        }

        // Did we correct the UMI?
        let gene = annotation.conf_mapped_feature(feature_reference).unwrap();
        let raw_key = (annotation.umi_info.seq, gene);
        let (corrected_umi, is_corrected) = match self.umi_corrections.get(&raw_key) {
            Some(umi) => (*umi, true),
            None => (raw_key.0, false),
        };

        let corrected_key = (corrected_umi, raw_key.1);
        let is_low_support_umi = self.low_support_umigenes.contains(&corrected_key);

        let is_min_qname = match self.umigene_min_key.get(&corrected_key) {
            Some(min_key) => annotation.read.header() == min_key.qname.as_slice(),
            None => unreachable!(),
        };

        let read_count = self.umigene_counts[&corrected_key];

        let is_filtered_target_umi = match self.targeted_umi_min_read_count {
            // low support UMIs and low read count UMIs are disjoint
            Some(threshold) => {
                let target_set = feature_reference.target_set.as_ref().unwrap();
                target_set.is_on_target(corrected_key.1 as u32)
                    && read_count < threshold
                    && !is_low_support_umi
            }
            None => false,
        };

        let sampling_factor = rng.gen_bool(barcode_subsample_rate);

        let is_umi_count =
            !is_low_support_umi && is_min_qname && !is_filtered_target_umi && sampling_factor;

        let umi_count = if is_umi_count {
            let umi_type = match annotation.is_conf_mapped_unique_txomic() {
                true => UmiType::Txomic,
                false => UmiType::NonTxomic,
            };
            let probe_id: Option<i32> = match &annotation.primary {
                RecordAnnotation::Probe(_, data) => {
                    if data.is_conf_mapped() {
                        let probe_id = &data.lhs_probe()?.probe_id;
                        Some(probe_set_reference?.probe_id_to_index(probe_id) as i32)
                    } else {
                        None
                    }
                }
                _ => None,
            };
            Some(UmiCount {
                library_idx: annotation.read.read_chunk(read_chunks).library_id(),
                feature_idx: corrected_key.1 as u32,
                probe_idx: probe_id,
                umi: corrected_umi.encode_2bit_u32(),
                read_count: read_count as u32,
                utype: umi_type,
            })
        } else {
            None
        };

        Some(DupInfo {
            processed_umi: corrected_umi,
            is_corrected,
            umi_count,
            is_umi_count,
            is_low_support_umi,
            is_filtered_target_umi,
        })
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use std::collections::HashMap;

    #[test]
    fn test_correct_umis() {
        let (g0, g1) = (0usize, 1usize);

        let mut umis = HashMap::new();
        umis.insert((UmiSeq::from_bytes(b"AAAA"), g0), 3u64);
        umis.insert((UmiSeq::from_bytes(b"AAAT"), g0), 2u64);
        umis.insert((UmiSeq::from_bytes(b"AAAA"), g1), 1u64);
        umis.insert((UmiSeq::from_bytes(b"AATT"), g1), 1u64);

        let corr = correct_umis(&umis);
        assert!(corr.len() == 1);
        assert!(corr[&(UmiSeq::from_bytes(b"AAAT"), g0)] == UmiSeq::from_bytes(b"AAAA"));

        let mut umis = HashMap::new();
        umis.insert((UmiSeq::from_bytes(b"CCCC"), g0), 1u64);
        umis.insert((UmiSeq::from_bytes(b"CGCC"), g0), 1u64);

        let corr = correct_umis(&umis);
        assert!(corr.len() == 1);
        assert!(corr[&(UmiSeq::from_bytes(b"CCCC"), g0)] == UmiSeq::from_bytes(b"CGCC"));
    }

    #[test]
    fn test_umi_type() {
        let key1 = UmiSelectKey {
            utype: UmiType::Txomic,
            qname: vec![1_u8],
        };
        let key2 = UmiSelectKey {
            utype: UmiType::NonTxomic,
            qname: vec![0_u8],
        };
        assert!(key1 < key2);
    }
}
