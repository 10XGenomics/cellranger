//
// Copyright (c) 2017 10x Genomics, Inc. All rights reserved.
//

use std::collections::{HashMap, HashSet, BTreeMap};
use std::fs::File;
use std::iter::FromIterator;
use std::cmp::{min, max};

use rust_htslib::bam::record::Record;
use serde_json;
use serde_json::Value;
use bincode;
use bincode::Infinite;
use rand;
use rand::{Rng, SeedableRng};

use transcriptome::{AnnotationData, PairAnnotationData, AnnotationRegion, AnnotationParams, Strand};
use reference::{TranscriptIndex, Gene};
use barcodes::BarcodeUmiData;
use features::FeatureData;
use metrics::{SequenceDistMetric, PercentMetric, InsertSizeMetric, PrefixGroup, MetricGroup, Metric};
use utils;

const MULTI_REFS_PREFIX: &'static str = "multi";
pub const HIGH_CONF_MAPQ: u8 = 255;

const TOP_RAW_SEQ_SAMPLE_RATE: f64 = 0.01;
const TOP_CORRECTED_SEQ_SAMPLE_RATE: f64 = 1.0;

const REGION_GENOME: &'static str = "genome";
const REGION_TRANSCRIPTOME: &'static str = "transcriptome";
const REGION_EXONIC: &'static str = "exonic";
const REGION_INTRONIC: &'static str = "intronic";
const REGION_INTERGENIC: &'static str = "intergenic";

const MAX_INSERT_SIZE: i64 = 1000;

#[derive(Eq, PartialEq, Ord, PartialOrd, Debug, Clone, Copy, Hash)]
enum InsertType {
    SingleEndThreePrime,
    SingleEndFivePrime,
    PairedEnd,
}

/// All data associated with a qname (single-end read or read pair)
pub struct ReadData {
    /// Data for each Read1 BAM record
    pub r1_data:        Vec<RecordData>,
    /// Data for each Read2 BAM record
    pub r2_data:        Vec<RecordData>,
    /// Cell barcode and UMI data
    pub bc_umi_data:    BarcodeUmiData,
    /// Feature barcode data
    pub feature_data:   Option<FeatureData>,
    /// Data for each pair of BAM records
    pub pair_data:      Vec<PairAnnotationData>,
}

/// All data associated with a single BAM record
pub struct RecordData {
    pub rec:    Record,
    pub anno:   AnnotationData,
}

/// Get the metric prefix for a given library type.
/// Some are hardcoded for historical reasons.
pub fn get_library_type_metric_prefix(lib_type: &str) -> String {
    match lib_type {
        "Gene Expression" => "".to_owned(),
        "CRISPR Guide Capture" => "CRISPR_".to_owned(),
        "Antibody Capture" => "ANTIBODY_".to_owned(),
        _ => format!("{}_", lib_type),
    }
}

/// Get the gene(s) that an alignment or pair of alignments aligned to.
/// For a single-end alignment, this is the genes from the alignment.
/// For a pair of alignments, this is the intersection of the genes
///   among mates with a non-zero number of genes.
pub fn get_alignment_gene_intersection(r1_anno: &AnnotationData, maybe_r2_anno: Option<&AnnotationData>) -> Vec<Gene> {
    let r1_genes = HashSet::from_iter(r1_anno.genes.iter());
    let r2_genes = match maybe_r2_anno {
        Some(r2_anno) => HashSet::from_iter(r2_anno.genes.iter()),
        None => HashSet::new(),
    };
    if r1_genes.len() > 0 && r2_genes.len() == 0 {
        r1_genes.into_iter().map(|x| x.to_owned()).collect()
    } else if r1_genes.len() == 0 && r2_genes.len() > 0 {
        r2_genes.into_iter().map(|x| x.to_owned()).collect()
    } else {
        r1_genes.intersection(&r2_genes).into_iter().map(|x| (*x).to_owned()).collect()
    }
}

impl ReadData {
    pub fn is_paired_end(&self) -> bool {
        return self.r1_data.len() > 0 && self.r2_data.len() > 0
    }

    pub fn is_properly_paired(&self) -> bool {
        return self.r1_data.len() == self.r2_data.len()
    }

    /// Mapped to genome
    pub fn is_mapped(&self) -> bool {
        return (self.r1_data.len() > 0 && !self.r1_data[0].rec.is_unmapped()) ||
               (self.r2_data.len() > 0 && !self.r2_data[0].rec.is_unmapped())
    }

    /// Get the index of the primary alignment record.
    /// Should always exist.
    pub fn get_primary_index(&self) -> Option<usize> {
        self.r1_data.iter().position(|x| !x.rec.is_secondary())
    }

    /// Mapped to a single non-gene-expression feature
    pub fn is_conf_mapped_to_non_gex_feature(&self) -> bool {
        if let Some(ref feature_data) = self.feature_data {
            if let Some(ref id_string) = feature_data.ids {
                // The presence of a semicolon in the feature IDs list
                // indicates mapping to multiple features.
                return id_string.chars().all(|c| c != ';')
            }
        }
        return false;
    }

    /// Is there a primary alignment that mapped uniquely to the
    /// genome and is compatible with a single gene
    pub fn is_conf_mapped_to_transcriptome(&self) -> bool {
        // Get primary record
        let primary = self.get_primary_index();
        if primary.is_none() {
            return false;
        }
        let primary = primary.unwrap();

        let r1 = &self.r1_data[primary].rec;

        if self.is_paired_end() {
            assert!(!self.r2_data[primary].rec.is_secondary());
            let r2 = &self.r2_data[primary].rec;

            // Both reads uniquely map and their gene intersection is unique
            return !r1.is_unmapped() && r1.mapq() >= HIGH_CONF_MAPQ &&
                !r2.is_unmapped() && r2.mapq() >= HIGH_CONF_MAPQ &&
                self.pair_data[primary].genes.len() == 1
        } else {

            // Single read uniquely maps and has a single gene
            return !r1.is_unmapped() && r1.mapq() >= HIGH_CONF_MAPQ &&
                self.r1_data[primary].anno.genes.len() == 1
        }
    }

    /// Do the records in the primary uniquely-mapped alignment pair
    /// that are compatible with at least one gene share
    /// no genes?
    pub fn is_gene_discordant(&self) -> bool {
        if !self.is_paired_end() {
            return false;
        }

        // Get primary record
        let primary = self.get_primary_index();
        if primary.is_none() {
            return false;
        }
        let primary = primary.unwrap();

        let r1 = &self.r1_data[primary].rec;
        let r2 = &self.r2_data[primary].rec;

        return !r1.is_unmapped() && r1.mapq() >= HIGH_CONF_MAPQ &&
            !r2.is_unmapped() && r2.mapq() >= HIGH_CONF_MAPQ &&
            self.pair_data[primary].genes.len() == 0
    }

    // TODO add iterators for zipping R1, R2, etc.
}

struct PairedEndStats {
    insert_sizes:    Vec<i64>,
    discordant_pair: bool, // inconsistent gene annotations
    improper_pair:   bool, // inconsistent number of alignments (e.g. only one end mapped)
}

fn get_mapping_region(ann: &AnnotationData) -> Option<String> {
    match ann.region {
        AnnotationRegion::Exonic        => Some(REGION_EXONIC.into()),
        AnnotationRegion::Intronic      => Some(REGION_INTRONIC.into()),
        AnnotationRegion::Intergenic    => Some(REGION_INTERGENIC.into()),
        AnnotationRegion::Unmapped      => None,
    }
}

#[derive(Serialize, Deserialize)]
struct RegionMetrics {
    mapped_reads:           PercentMetric,
    conf_mapped_reads:      PercentMetric,
    conf_mapped_bc_reads:   PercentMetric,
}

impl MetricGroup for RegionMetrics {
    fn new() -> Self {
        RegionMetrics {
            mapped_reads: PercentMetric::new(),
            conf_mapped_reads: PercentMetric::new(),
            conf_mapped_bc_reads: PercentMetric::new(),
        }
    }

    fn merge(&mut self, other: &Self) {
        self.mapped_reads.merge(&other.mapped_reads);
        self.conf_mapped_reads.merge(&other.conf_mapped_reads);
        self.conf_mapped_bc_reads.merge(&other.conf_mapped_bc_reads);
    }

    fn report(&self) -> BTreeMap<String, Value> {
        let mut results = BTreeMap::new();
        results.insert("mapped_reads_frac".into(), self.mapped_reads.report());
        results.insert("conf_mapped_reads_frac".into(), self.conf_mapped_reads.report());
        results.insert("conf_mapped_barcoded_reads_frac".into(), self.conf_mapped_bc_reads.report());
        return results
    }
}

#[derive(Serialize, Deserialize)]
struct MappingMetrics {
    antisense_reads:    PercentMetric,
    insert_sizes:       InsertSizeMetric,
    region_metrics:     PrefixGroup<RegionMetrics>,
    conf_mapped_bcs:    SequenceDistMetric,
    conf_mapped_umis:   SequenceDistMetric,
    improper_pairs:     PercentMetric,
    discordant_pairs:   PercentMetric,
}

impl MetricGroup for MappingMetrics {
    fn new() -> Self {
        let regions = vec![REGION_GENOME.into(), REGION_TRANSCRIPTOME.into(), REGION_EXONIC.into(),
                           REGION_INTRONIC.into(), REGION_INTERGENIC.into()];
        return MappingMetrics {
            antisense_reads:    PercentMetric::new(),
            insert_sizes:       InsertSizeMetric::new(),
            region_metrics:     PrefixGroup::new(&regions),
            conf_mapped_bcs:    SequenceDistMetric::new(),
            conf_mapped_umis:   SequenceDistMetric::new(),
            improper_pairs:     PercentMetric::new(),
            discordant_pairs:   PercentMetric::new(),
        }
    }

    fn merge(&mut self, other: &Self) {
        self.antisense_reads.merge(&other.antisense_reads);
        self.insert_sizes.merge(&other.insert_sizes);
        self.region_metrics.merge(&other.region_metrics);
        self.conf_mapped_bcs.merge(&other.conf_mapped_bcs);
        self.conf_mapped_umis.merge(&other.conf_mapped_umis);
        self.improper_pairs.merge(&other.improper_pairs);
        self.discordant_pairs.merge(&other.discordant_pairs);
    }

    fn report(&self) -> BTreeMap<String, Value> {
        let mut results = BTreeMap::new();
        results.insert("antisense_reads_frac".into(), self.antisense_reads.report());
        results.insert("median_insert_size".into(), self.insert_sizes.report_median());
        results.insert("iqr_insert_size".into(), self.insert_sizes.report_iqr());
        results.insert("insert_size_histogram".into(), self.insert_sizes.report_binned());
        results.insert("conf_mapped_effective_barcode_diversity".into(), self.conf_mapped_bcs.report_effective_diversity());
        results.insert("conf_mapped_effective_umi_diversity".into(), self.conf_mapped_umis.report_effective_diversity());
        results.insert("conf_mapped_top_processed_umis".into(), self.conf_mapped_umis.report_top_n());
        results.insert("discordant_pairs_frac".into(), self.discordant_pairs.report());
        results.insert("improper_pairs_frac".into(), self.improper_pairs.report());
        results.extend(self.region_metrics.report());
        return results
    }
}

#[derive(Serialize, Deserialize)]
struct ReadMetrics {
    unmapped_reads:     PercentMetric,      // Reads not mapped to genome
    good_umi_reads:     PercentMetric,      // UMI passes the various filters
    good_bc_reads:      PercentMetric,      // barcode on whitelist post-correction
    corrected_bc_reads: PercentMetric,      // reads where cell barcode was corrected
    mapping_metrics:    PrefixGroup<MappingMetrics>, // per-genome
    feature_bc_extracted_reads: PercentMetric, // reads were feature barcode was extracted
    corrected_feature_bc_reads: PercentMetric, // reads where feature barcode was corrected, given extracted
    unrecognized_feature_bc_reads: PercentMetric, // reads feature barcode was unrecognized, given extracted
}

impl ReadMetrics {
    fn set_genomes(&mut self, genomes: &[String]) {
        self.mapping_metrics = PrefixGroup::new(&genomes);
    }
}

impl MetricGroup for ReadMetrics {
    fn new() -> Self {
        let genomes = Vec::new();
        ReadMetrics {
            unmapped_reads:     PercentMetric::new(),
            good_umi_reads:     PercentMetric::new(),
            good_bc_reads:      PercentMetric::new(),
            corrected_bc_reads: PercentMetric::new(),
            mapping_metrics:    PrefixGroup::new(&genomes),
            feature_bc_extracted_reads: PercentMetric::new(),
            corrected_feature_bc_reads: PercentMetric::new(),
            unrecognized_feature_bc_reads: PercentMetric::new(),
        }
    }

    fn merge(&mut self, other: &Self) {
        self.unmapped_reads.merge(&other.unmapped_reads);
        self.good_umi_reads.merge(&other.good_umi_reads);
        self.good_bc_reads.merge(&other.good_bc_reads);
        self.corrected_bc_reads.merge(&other.corrected_bc_reads);
        self.feature_bc_extracted_reads.merge(&other.feature_bc_extracted_reads);
        self.corrected_feature_bc_reads.merge(&other.corrected_feature_bc_reads);
        self.unrecognized_feature_bc_reads.merge(&other.unrecognized_feature_bc_reads);
        self.mapping_metrics.merge(&other.mapping_metrics);
    }

    fn report(&self) -> BTreeMap<String, Value> {
        let mut results = BTreeMap::new();
        results.insert("unmapped_reads_frac".into(), self.unmapped_reads.report());
        results.insert("good_umi_frac".into(), self.good_umi_reads.report());
        results.insert("good_bc_frac".into(), self.good_bc_reads.report());
        results.insert("corrected_bc_frac".into(), self.corrected_bc_reads.report());
        results.insert("feature_bc_extracted_frac".into(), self.feature_bc_extracted_reads.report());
        results.insert("corrected_feature_bc_frac".into(), self.corrected_feature_bc_reads.report());
        results.insert("unrecognized_feature_bc_frac".into(), self.unrecognized_feature_bc_reads.report());
        results.extend(self.mapping_metrics.report());
        return results
    }
}

#[derive(Serialize, Deserialize)]
struct UmiMetrics {
    has_n_umis:                 PercentMetric,
    homopolymer_umis:           PercentMetric,
    low_min_qual_umis:          PercentMetric,
    primer_umis:                PercentMetric,
    polyt_umis:                 PercentMetric,
    filtered_has_n_umis:        PercentMetric,
    filtered_homopolymer_umis:  PercentMetric,
    filtered_low_min_qual_umis: PercentMetric,
    top_raw_umis:               SequenceDistMetric,
    top_processed_umis:         SequenceDistMetric,
}

impl MetricGroup for UmiMetrics {
    fn new() -> Self {
        UmiMetrics {
            has_n_umis:                 PercentMetric::new(),
            homopolymer_umis:           PercentMetric::new(),
            low_min_qual_umis:          PercentMetric::new(),
            primer_umis:                PercentMetric::new(),
            polyt_umis:                 PercentMetric::new(),
            filtered_has_n_umis:        PercentMetric::new(),
            filtered_homopolymer_umis:  PercentMetric::new(),
            filtered_low_min_qual_umis: PercentMetric::new(),
            top_raw_umis:               SequenceDistMetric::new(),
            top_processed_umis:         SequenceDistMetric::new(),
        }
    }

    fn merge(&mut self, other: &Self) {
        self.has_n_umis.merge(&other.has_n_umis);
        self.homopolymer_umis.merge(&other.homopolymer_umis);
        self.low_min_qual_umis.merge(&other.low_min_qual_umis);
        self.primer_umis.merge(&other.primer_umis);
        self.polyt_umis.merge(&other.polyt_umis);
        self.filtered_has_n_umis.merge(&other.filtered_has_n_umis);
        self.filtered_homopolymer_umis.merge(&other.filtered_homopolymer_umis);
        self.filtered_low_min_qual_umis.merge(&other.filtered_low_min_qual_umis);
        self.top_raw_umis.merge(&other.top_raw_umis);
        self.top_processed_umis.merge(&other.top_processed_umis);
    }

    fn report(&self) -> BTreeMap<String, Value> {
        let mut results = BTreeMap::new();
        results.insert("has_n_umi_property_frac".into(), self.has_n_umis.report());
        results.insert("homopolymer_umi_property_frac".into(), self.homopolymer_umis.report());
        results.insert("low_min_qual_umi_property_frac".into(), self.low_min_qual_umis.report());
        results.insert("primer_umi_property_frac".into(), self.primer_umis.report());
        results.insert("polyt_suffix_umi_property_frac".into(), self.polyt_umis.report());
        results.insert("has_n_umi_filter_frac".into(), self.filtered_has_n_umis.report());
        results.insert("homopolymer_umi_filter_frac".into(), self.filtered_homopolymer_umis.report());
        results.insert("low_min_qual_umi_filter_frac".into(), self.filtered_low_min_qual_umis.report());
        results.insert("top_raw_umis".into(), self.top_raw_umis.report_top_n());
        results.insert("top_processed_umis".into(), self.top_processed_umis.report_top_n());
        results.insert("effective_umi_diversity".into(), self.top_processed_umis.report_effective_diversity());
        return results
    }
}

#[derive(Serialize, Deserialize)]
struct BarcodeMetrics {
    miss_whitelist_bcs:             PercentMetric,
    has_n_bcs:                      PercentMetric,
    homopolymer_bcs:                PercentMetric,
    low_min_qual_bcs:               PercentMetric,
    filtered_miss_whitelist_bcs:    PercentMetric,
    top_raw_bcs:                    SequenceDistMetric,
    top_processed_bcs:              SequenceDistMetric,
}

impl MetricGroup for BarcodeMetrics {
    fn new() -> Self {
        BarcodeMetrics {
            miss_whitelist_bcs:             PercentMetric::new(),
            has_n_bcs:                      PercentMetric::new(),
            homopolymer_bcs:                PercentMetric::new(),
            low_min_qual_bcs:               PercentMetric::new(),
            filtered_miss_whitelist_bcs:    PercentMetric::new(),
            top_raw_bcs:                    SequenceDistMetric::new(),
            top_processed_bcs:              SequenceDistMetric::new(),
        }
    }

    fn merge(&mut self, other: &Self) {
        self.miss_whitelist_bcs.merge(&other.miss_whitelist_bcs);
        self.has_n_bcs.merge(&other.has_n_bcs);
        self.homopolymer_bcs.merge(&other.homopolymer_bcs);
        self.low_min_qual_bcs.merge(&other.low_min_qual_bcs);
        self.filtered_miss_whitelist_bcs.merge(&other.filtered_miss_whitelist_bcs);
        self.top_raw_bcs.merge(&other.top_raw_bcs);
        self.top_processed_bcs.merge(&other.top_processed_bcs);
    }

    fn report(&self) -> BTreeMap<String, Value> {
        let mut results = BTreeMap::new();
        results.insert("has_n_barcode_property_frac".into(), self.has_n_bcs.report());
        results.insert("homopolymer_barcode_property_frac".into(), self.homopolymer_bcs.report());
        results.insert("low_min_qual_barcode_property_frac".into(), self.low_min_qual_bcs.report());
        results.insert("miss_whitelist_barcode_property_frac".into(), self.miss_whitelist_bcs.report());
        results.insert("miss_whitelist_barcode_filter_frac".into(), self.filtered_miss_whitelist_bcs.report());
        results.insert("top_raw_barcodes".into(), self.top_raw_bcs.report_top_n());
        results.insert("effective_barcode_diversity".into(), self.top_processed_bcs.report_effective_diversity());
        results.insert("barcodes_detected".into(), self.top_processed_bcs.report_num_seqs());
        return results
    }
}

#[derive(Serialize, Deserialize)]
pub struct Metrics {
    read_metrics:       ReadMetrics,
    umi_metrics:        UmiMetrics,
    barcode_metrics:    BarcodeMetrics,
}

impl MetricGroup for Metrics {
    fn new() -> Self {
        Metrics {
            read_metrics:       ReadMetrics::new(),
            umi_metrics:        UmiMetrics::new(),
            barcode_metrics:    BarcodeMetrics::new(),

        }
    }

    fn merge(&mut self, other: &Self) {
        self.read_metrics.merge(&other.read_metrics);
        self.umi_metrics.merge(&other.umi_metrics);
        self.barcode_metrics.merge(&other.barcode_metrics);
    }

    fn report(&self) -> BTreeMap<String, Value> {
        let mut results = BTreeMap::new();
        results.extend(self.read_metrics.report());
        results.extend(self.umi_metrics.report());
        results.extend(self.barcode_metrics.report());
        return results
    }
}

impl Metrics {
    pub fn write_json_summary(&self, filename: &String) {
        let summary = self.report();
        let writer = File::create(filename).unwrap();
        serde_json::to_writer_pretty(writer, &summary).expect("Failed to write JSON");
    }

    pub fn write_barcodes_csv(&self, filename: &String) {
        let barcodes = self.barcode_metrics.top_processed_bcs.get_seqs();
        utils::save_txt(barcodes, filename);
    }

    pub fn write_binary(&self, filename: &String) {
        let mut writer = File::create(filename).unwrap();
        bincode::serialize_into(&mut writer, self, Infinite).expect("Failed to serialize binary data");
    }

    pub fn read_binary(filename: &String) -> Self {
        let mut reader = File::open(filename).unwrap();
        return bincode::deserialize_from(&mut reader, Infinite).unwrap()
    }
}

pub struct Reporter<'a> {
    genomes:            Vec<String>,
    chroms:             Vec<String>,
    chroms_to_genomes:  HashMap<String, String>,
    transcript_idx:     &'a TranscriptIndex,
    rng:                rand::XorShiftRng,
    metrics:            Metrics,
    params:             &'a AnnotationParams,
}

fn map_chroms_to_genomes(chroms: &[String], genomes: &[String]) -> HashMap<String, String> {
    let mut mapping = HashMap::new();
    for chrom in chroms {
        if genomes.len() == 1 {
            mapping.insert(chrom.clone(), genomes[0].clone());
        } else {
            let genome = chrom.split("_").next().unwrap();
            mapping.insert(chrom.clone(), genome.to_owned());
        }
    }
    return mapping
}

impl<'a> Reporter<'a> {
    pub fn new(chroms: Vec<String>, genomes: Vec<String>, transcript_idx: &'a TranscriptIndex, params: &'a AnnotationParams) -> Reporter<'a> {
        let chroms_to_genomes = map_chroms_to_genomes(&chroms, &genomes);
        let seed = [1, 2, 3, 4];
        let mut reporter = Reporter {
            genomes:            genomes,
            chroms:             chroms,
            chroms_to_genomes:  chroms_to_genomes,
            transcript_idx:     transcript_idx,
            rng:                rand::XorShiftRng::from_seed(seed),
            metrics:            Metrics::new(),
            params:             params,
        };
        // TODO - hack to get around the fact that metric constructors take no args
        let mut genome_prefixes = vec![MULTI_REFS_PREFIX.into()];
        for genome in &reporter.genomes { genome_prefixes.push(genome.clone()); }
        reporter.metrics.read_metrics.set_genomes(&genome_prefixes);
        return reporter
    }

    pub fn get_genome(&self, read: &Record) -> Option<&String> {
        if read.is_unmapped() {
            return None
        } else if self.genomes.len() == 1 {
            return Some(&self.genomes[0])
        } else {
            let chrom = &self.chroms[read.tid() as usize];
            return self.chroms_to_genomes.get(chrom)
        }
    }

    fn get_paired_metrics(&mut self, data: &ReadData, insert_type: &InsertType) -> PairedEndStats {
        let mut insert_sizes = Vec::new();
        let mut discordant_pair = false;
        let mut improper_pair = false;
        match insert_type {
            &InsertType::PairedEnd => {
                if data.is_properly_paired() {
                    for i in 0..data.r1_data.len() {
                        let tx_ids1: HashSet<_> = HashSet::from_iter(data.r1_data[i].anno.transcripts.keys());
                        let tx_ids2: HashSet<_> = HashSet::from_iter(data.r2_data[i].anno.transcripts.keys());
                        let tx_intersection: HashSet<_> = tx_ids1.intersection(&tx_ids2).collect();
                        let tx_union: HashSet<_> = tx_ids1.union(&tx_ids2).collect();
                        for shared_tx_id in &tx_intersection {
                            let tx1 = data.r1_data[i].anno.transcripts.get(**shared_tx_id).unwrap();
                            let tx2 = data.r2_data[i].anno.transcripts.get(**shared_tx_id).unwrap();
                            let start = min(tx1.pos, tx2.pos);
                            let end = max(tx1.pos + tx1.alen, tx2.pos + tx2.alen);
                            let insert_size = end - start;
                            if insert_size <= MAX_INSERT_SIZE {
                                insert_sizes.push(insert_size);
                            }
                        }
                        if tx_intersection.len() == 0 && tx_union.len() > 0 {
                            discordant_pair = true;
                        }
                    }
                } else {
                    improper_pair = true;
                }
            },
            _ => {
                for read in &data.r1_data {
                    for tx in read.anno.transcripts.values() {
                        let insert_size = match insert_type {
                            &InsertType::SingleEndThreePrime => self.transcript_idx.get_transcript_length(&tx.id) - tx.pos,
                            &InsertType::SingleEndFivePrime => tx.pos + tx.alen,
                            &InsertType::PairedEnd => panic!("Should never happen"),
                        };
                        if insert_size <= MAX_INSERT_SIZE && tx.strand == Strand::Forward {
                            insert_sizes.push(insert_size);
                        }
                    }
                }
            },
        }
        return PairedEndStats { insert_sizes: insert_sizes, discordant_pair: discordant_pair, improper_pair: improper_pair }
    }

    pub fn update_read_metrics(&mut self, read_data: &ReadData) {
        match &read_data.bc_umi_data.barcode_data {
            &Some(ref bc_data) => {
                let mut good_bc = false;
                let mut corrected_bc = false;
                if bc_data.processed_seq.is_some() {
                    good_bc = true;
                    // If original bc not on whitelist, it was corrected
                    if !bc_data.on_whitelist {
                        corrected_bc = true;
                    }
                }
                self.metrics.read_metrics.good_bc_reads.add(&good_bc);
                self.metrics.read_metrics.corrected_bc_reads.add(&corrected_bc);
            },
            &None => {},
        }

        match &read_data.bc_umi_data.umi_data {
            &Some(ref umi_data) => {
                let good_umi = umi_data.is_valid;
                self.metrics.read_metrics.good_umi_reads.add(&good_umi);
            },
            &None => {},
        }

        self.metrics.read_metrics.feature_bc_extracted_reads.add(&read_data.feature_data.is_some());

        if let &Some(ref feature_data) = &read_data.feature_data {
            // In order to have feature_data, there must have been a raw feature sequence
            // So these metrics are conditional on "a feature barcode was extracted."
            self.metrics.read_metrics.corrected_feature_bc_reads.add(&feature_data.was_corrected);
            self.metrics.read_metrics.unrecognized_feature_bc_reads.add(&feature_data.corrected_seq.is_none());
        }

        let unmapped = !read_data.is_mapped();
        self.metrics.read_metrics.unmapped_reads.add(&unmapped);

        self.update_mapping_metrics(read_data);
    }

    fn count_regions(&self, genome_region_set: &mut HashSet<(String, String)>, genome: &String, region: &String, include_txome: bool) {
        // helper function for update_mapping_metrics
        for prefix in vec![genome.clone(), MULTI_REFS_PREFIX.into()] {
            genome_region_set.insert((prefix.clone(), REGION_GENOME.into()));
            genome_region_set.insert((prefix.clone(), region.clone()));
            if include_txome { genome_region_set.insert((prefix.clone(), REGION_TRANSCRIPTOME.into())); }
        }
    }

    pub fn update_mapping_metrics(&mut self, read_data: &ReadData) {
        let corrected_bc_seq = match &read_data.bc_umi_data.barcode_data {
            &Some(ref bc_data) => match &bc_data.processed_seq {
                &Some(ref seq) => Some(seq),
                &None => None,
            },
            &None => None,
        };
        let valid_umi = match &read_data.bc_umi_data.umi_data {
            &Some(ref umi_data) => if umi_data.is_valid { Some(&umi_data.raw_seq) } else { None },
            &None => None,
        };

        // precompute genomes, genes, and regions
        // TODO make this less ugly
        let mut read_mapped: HashSet<(String, String)> = HashSet::new();
        let mut read_conf_mapped: HashSet<(String, String)> = HashSet::new();
        let mut read_conf_mapped_barcoded: HashSet<(String, String)> = HashSet::new();
        let mut read_antisense: HashSet<String> = HashSet::new();
        let mut pair_discordant: HashSet<String> = HashSet::new();
        let mut pair_improper: HashSet<String> = HashSet::new();

        let insert_type = if read_data.is_paired_end() {
            InsertType::PairedEnd
        } else if self.params.chemistry_fiveprime {
            InsertType::SingleEndFivePrime
        } else {
            InsertType::SingleEndThreePrime
        };

        let mut paired_stats = self.get_paired_metrics(read_data, &insert_type);
        paired_stats.insert_sizes.sort();
        let median_insert_size = utils::median(&paired_stats.insert_sizes) as u64;

        let is_mapped = read_data.is_mapped();
        let is_conf_mapped_to_non_gex_feature = read_data.is_conf_mapped_to_non_gex_feature();
        let is_conf_mapped_to_transcriptome = read_data.is_conf_mapped_to_transcriptome();

        // Assert that "mapped to genome" and "conf mapped to non-gene-expression feature"
        //   are mutually exclusive. I.e., non-GEX feature reads are always unmapped.
        assert!(is_mapped && !is_conf_mapped_to_non_gex_feature || !is_mapped);

        if is_conf_mapped_to_non_gex_feature {
            // Confidently mapped
            read_conf_mapped.insert((MULTI_REFS_PREFIX.to_string(), // genome is meaningless
                                     REGION_TRANSCRIPTOME.to_string())); // pretend it's the transcriptome

            // and barcoded
            if corrected_bc_seq.is_some() {
                read_conf_mapped_barcoded.insert((MULTI_REFS_PREFIX.to_string(),
                                                  REGION_TRANSCRIPTOME.to_string()));
            }
        }

        if is_mapped {
            for rdata in read_data.r1_data.iter().chain(read_data.r2_data.iter()) {
                let alignment = &rdata.rec;
                let annotation = &rdata.anno;
                let (genome, region) = match (self.get_genome(alignment),
                                              get_mapping_region(annotation)) {
                    (Some(g), Some(r)) => (g, r),
                    _ => continue, // should only happen for improper / halfmapped pairs
                };
                self.count_regions(&mut read_mapped, &genome, &region,
                                   annotation.genes.len() > 0);

                if paired_stats.discordant_pair {
                    pair_discordant.insert(genome.clone());
                    pair_discordant.insert(MULTI_REFS_PREFIX.into());
                }

                if paired_stats.improper_pair {
                    pair_improper.insert(genome.clone());
                    pair_improper.insert(MULTI_REFS_PREFIX.into());
                }

                if alignment.mapq() >= HIGH_CONF_MAPQ {
                    self.count_regions(&mut read_conf_mapped, &genome, &region,
                                       is_conf_mapped_to_transcriptome);

                    if annotation.is_antisense() {
                        read_antisense.insert(genome.clone());
                        read_antisense.insert(MULTI_REFS_PREFIX.into());
                    }

                    if corrected_bc_seq.is_some() { // conf mapped and barcoded
                        self.count_regions(&mut read_conf_mapped_barcoded, &genome, &region,
                                           is_conf_mapped_to_transcriptome);
                    }
                }
            }
        }

        // update metrics
        let ref mut mapping_metrics = self.metrics.read_metrics.mapping_metrics;
        for (genome, genome_metrics) in mapping_metrics.iter_mut() {
            // NOTE: in the Python code, we only track insert size for confidently mapped reads, so do that here as well
            if read_data.is_conf_mapped_to_transcriptome() && paired_stats.insert_sizes.len() > 0 { genome_metrics.insert_sizes.add(&median_insert_size); }
            genome_metrics.antisense_reads.add(&read_antisense.contains(genome));
            genome_metrics.improper_pairs.add(&pair_improper.contains(genome));
            genome_metrics.discordant_pairs.add(&pair_discordant.contains(genome));
            for (region, region_metrics) in genome_metrics.region_metrics.iter_mut() {
                let tuple = &(genome.clone(), region.clone());
                region_metrics.mapped_reads.add(&read_mapped.contains(tuple));
                region_metrics.conf_mapped_reads.add(&read_conf_mapped.contains(tuple));
                region_metrics.conf_mapped_bc_reads.add(&read_conf_mapped_barcoded.contains(tuple));
            }
            if read_conf_mapped.contains(&(genome.clone(), REGION_GENOME.into())) {
                match corrected_bc_seq {
                    Some(seq) => genome_metrics.conf_mapped_bcs.add(seq),
                    None => {},
                };
                match valid_umi {
                    Some(seq) => genome_metrics.conf_mapped_umis.add(seq),
                    None => {},
                };
            }
        }
    }

    pub fn update_bc_umi_metrics(&mut self, data: &BarcodeUmiData) {
        // barcodes
        match &data.barcode_data {
            &Some(ref bc_data) => {
                let ref mut bc_metrics = self.metrics.barcode_metrics;

                if self.rng.next_f64() < TOP_RAW_SEQ_SAMPLE_RATE {
                    bc_metrics.top_raw_bcs.add(&bc_data.raw_seq);
                }

                match &bc_data.processed_seq {
                    &Some(ref seq) => bc_metrics.top_processed_bcs.add(&seq),
                    &None => {},
                }

                bc_metrics.has_n_bcs.add(&bc_data.has_n);
                bc_metrics.homopolymer_bcs.add(&bc_data.is_homopolymer);
                bc_metrics.low_min_qual_bcs.add(&bc_data.low_min_qual);
                bc_metrics.miss_whitelist_bcs.add(&!bc_data.on_whitelist);

                bc_metrics.filtered_miss_whitelist_bcs.add(&!bc_data.on_whitelist);
            },
            &None => {},
        }

        // UMIs
        match &data.umi_data {
            &Some(ref umi_data) => {
                let ref mut umi_metrics = self.metrics.umi_metrics;

                if self.rng.next_f64() < TOP_RAW_SEQ_SAMPLE_RATE {
                    umi_metrics.top_raw_umis.add(&umi_data.raw_seq);
                }
                if umi_data.is_valid && self.rng.next_f64() < TOP_CORRECTED_SEQ_SAMPLE_RATE {
                    umi_metrics.top_processed_umis.add(&umi_data.raw_seq);
                }

                umi_metrics.low_min_qual_umis.add(&umi_data.low_min_qual);
                umi_metrics.has_n_umis.add(&umi_data.has_n);
                umi_metrics.homopolymer_umis.add(&umi_data.is_homopolymer);
                umi_metrics.polyt_umis.add(&umi_data.has_polyt);
                umi_metrics.primer_umis.add(&umi_data.has_primer);

                let mut umi_filter = 0;
                if umi_data.low_min_qual { umi_filter = 1; }
                else if umi_data.has_n { umi_filter = 2; }
                else if umi_data.is_homopolymer { umi_filter = 3; }

                umi_metrics.filtered_low_min_qual_umis.add(&(umi_filter == 1));
                umi_metrics.filtered_has_n_umis.add(&(umi_filter == 2));
                umi_metrics.filtered_homopolymer_umis.add(&(umi_filter == 3));
            },
            &None => {},
        }
    }

    pub fn update_metrics(&mut self, read_data: &ReadData) {
        self.update_read_metrics(read_data);
        self.update_bc_umi_metrics(&read_data.bc_umi_data);
    }

    pub fn get_metrics(&self) -> &Metrics {
        return &self.metrics
    }
}
