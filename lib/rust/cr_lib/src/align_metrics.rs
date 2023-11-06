use crate::aligner::BarcodeSummary;
use anyhow::Result;
use barcode::{Barcode, HasBarcode};
use cr_types::probe_set::{MAPQ_HALF_MAPPED, MAPQ_SPLIT_MAPPED};
use cr_types::reference::genome_of_chrom::GenomeName;
use cr_types::rna_read::RnaRead;
use cr_types::types::LibraryFeatures;
use fxhash::FxHashMap;
use itertools::Itertools;
use json_report_derive::JsonReport;
use metric::{CountMetric, JsonReport, JsonReporter, MeanMetric, Metric, PercentMetric};
use metric_derive::Metric;
use serde::{Deserialize, Serialize};
use shardio::ShardSender;
use std::borrow::Cow;
use std::collections::HashSet;
use std::convert::From;
use std::fmt;
use std::io::Write;
use strum::IntoEnumIterator;
use strum_macros::{Display, EnumIter};
use transcriptome::Gene;
use tx_annotation::mark_dups::DupInfo;
use tx_annotation::read::{AnnotationInfo, ReadAnnotations};
use tx_annotation::transcript::AnnotationRegion;
use tx_annotation::visitor::AnnotatedReadVisitor;
use MappingRegion::{Genome, Transcriptome};
pub const MULTI_GENOME: &str = "multi";

#[derive(Copy, Clone, PartialEq, Eq, Serialize, Deserialize, Display, EnumIter, Hash)]
#[strum(serialize_all = "snake_case")]
pub enum MappingRegion {
    Genome,
    Transcriptome,
    Exonic,
    Intronic,
    Intergenic,
}

impl From<AnnotationRegion> for MappingRegion {
    fn from(item: AnnotationRegion) -> Self {
        match item {
            AnnotationRegion::Exonic => MappingRegion::Exonic,
            AnnotationRegion::Intronic => MappingRegion::Intronic,
            AnnotationRegion::Intergenic => MappingRegion::Intergenic,
        }
    }
}

#[derive(Serialize, Deserialize, Clone, Hash, PartialEq, Eq)]
pub struct GenomeMapping {
    pub genome: GenomeName,
    pub region: MappingRegion,
}

impl GenomeMapping {
    fn new(genome: impl Into<GenomeName>, region: impl Into<MappingRegion>) -> Self {
        GenomeMapping {
            genome: genome.into(),
            region: region.into(),
        }
    }
    fn multi(region: impl Into<MappingRegion>) -> Self {
        GenomeMapping::new(MULTI_GENOME, region)
    }
}

impl fmt::Display for GenomeMapping {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}_{}", self.genome, self.region)
    }
}

// When there is a target set, a read is targeted if it is annotated
// with at least one targeted gene, and untargeted otherwise
#[derive(Copy, Clone, PartialEq, Eq, Serialize, Deserialize, Display, EnumIter, Hash)]
#[strum(serialize_all = "snake_case")]
pub enum TargetingStatus {
    Targeted,
    Untargeted,
}

// When there is a target set, each read is bucketed for metrics
// based on genomic/transcriptomic alignment and targeting status
#[derive(Serialize, Deserialize, Clone, Hash, PartialEq, Eq)]
struct TargetedMapping {
    genome: GenomeName,
    region: MappingRegion,
    targeting: TargetingStatus,
}

impl TargetedMapping {
    fn new(
        genome: impl Into<GenomeName>,
        region: impl Into<MappingRegion>,
        targeting: impl Into<TargetingStatus>,
    ) -> Self {
        TargetedMapping {
            genome: genome.into(),
            region: region.into(),
            targeting: targeting.into(),
        }
    }
    fn multi(region: impl Into<MappingRegion>, targeting: TargetingStatus) -> Self {
        TargetedMapping::new(MULTI_GENOME, region, targeting)
    }
}

impl fmt::Display for TargetedMapping {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}_{}_{}", self.genome, self.region, self.targeting)
    }
}

/// Each AlignAndCountVisitor processes one library type
pub struct AlignAndCountVisitor {
    /// Buffer to store metrics for the current barcode
    metrics: VisitorMetrics,
    /// The visitor is expected to process a single library type
    library_feats: LibraryFeatures,
    #[allow(dead_code)]
    /// If targeted, the set of genes which are on target
    target_genes: Option<HashSet<Gene>>,
    pub barcode_summaries: Vec<BarcodeSummary>,
    pub last_barcode: Option<Barcode>,
    metrics_sender: ShardSender<BarcodeMetrics, LibFeatThenBarcodeOrder>,
}

impl AlignAndCountVisitor {
    pub fn new(
        library_feats: LibraryFeatures,
        metrics_sender: ShardSender<BarcodeMetrics, LibFeatThenBarcodeOrder>,
        target_genes: Option<HashSet<Gene>>,
    ) -> Self {
        AlignAndCountVisitor {
            metrics: VisitorMetrics::new(),
            library_feats,
            target_genes,
            barcode_summaries: Vec::new(),
            last_barcode: None,
            metrics_sender,
        }
    }
    fn final_send(&mut self) -> Result<()> {
        if let Some(bc) = self.last_barcode {
            self.metrics_sender.send(BarcodeMetrics {
                barcode: bc.into(),
                library_feats: self.library_feats,
                metrics: self.metrics.clone(),
            })?;
        }
        self.metrics_sender.finished()?;
        self.last_barcode = None;
        Ok(())
    }
    pub fn finish(mut self) -> Result<()> {
        self.final_send()
    }
}

impl Drop for AlignAndCountVisitor {
    fn drop(&mut self) {
        self.final_send().unwrap();
    }
}

#[derive(Serialize, Deserialize, Clone)]
pub struct BarcodeMetrics {
    pub barcode: BarcodeKind,
    pub library_feats: LibraryFeatures,
    pub metrics: VisitorMetrics,
}

impl BarcodeMetrics {
    pub fn to_csv_header<W: Write>(genomes: &[GenomeName], stream: &mut W) {
        let mut per_genome_metrics = [
            "mapped_reads",
            "conf_reads",
            "ontarget_reads",
            "offtarget_reads",
            "conf_ontarget_reads",
            "conf_offtarget_reads",
            "conf_intergenic_reads",
            "conf_exonic_reads",
            "conf_intronic_reads",
        ]
        .map(String::from)
        .to_vec();
        for region in ["exonic", "intronic"] {
            per_genome_metrics.push(format!("conf_{region}_unique_reads"));
            per_genome_metrics.push(format!("conf_{region}_antisense_reads"));
            per_genome_metrics.push(format!("conf_{region}_low_support_umi_reads"));
            per_genome_metrics.push(format!("conf_{region}_dup_reads"));
            per_genome_metrics.push(format!("{region}_umis"));
        }
        per_genome_metrics.push("conf_txomic_unique_reads".to_string());
        let mut header = [
            "barcode",
            "sequenced_reads",
            "barcode_corrected_sequenced_reads",
            "raw_reads",
        ]
        .map(String::from)
        .to_vec();
        header.reserve_exact(1 + per_genome_metrics.len() * genomes.len());
        for genome in genomes {
            for metric in &per_genome_metrics {
                header.push(format!("{metric}_{genome}"));
            }
        }
        header.push("mean_insert_size".to_string());
        stream
            .write_all(format!("{}\n", header.join(",")).as_bytes())
            .unwrap();
    }
    pub fn to_csv_row<W: Write>(&self, genomes: &[GenomeName], stream: &mut W) {
        if self.library_feats != LibraryFeatures::gex() {
            return;
        }
        let bc_string = match self.barcode {
            BarcodeKind::Valid(bc) => bc.to_string(),
            BarcodeKind::Invalid => String::from("NO_BARCODE"),
        };
        let mut csv_row = vec![
            bc_string,
            self.metrics.sequenced_reads.count().to_string(),
            self.metrics
                .barcode_corrected_sequenced_reads
                .count()
                .to_string(),
            self.metrics.total_aligner_reads.count().to_string(),
        ];
        let zero_count = CountMetric::new();
        for genome in genomes {
            let per_genome_mapping_metrics = self
                .metrics
                .per_genome_mapping
                .get(&GenomeMapping {
                    genome: genome.to_string(),
                    region: MappingRegion::Genome,
                })
                .cloned()
                .unwrap_or_else(Metric::new);
            // mapped_reads
            csv_row.push(per_genome_mapping_metrics.mapped.count().to_string());

            // conf reads
            csv_row.push(per_genome_mapping_metrics.conf_mapped.count().to_string());

            // on and off target mapped and conf_mapped reads
            for which_metric in [
                |m: &PerMappingStatusMetrics| m.mapped,
                |m: &PerMappingStatusMetrics| m.conf_mapped,
            ] {
                for targeting in [TargetingStatus::Targeted, TargetingStatus::Untargeted] {
                    csv_row.push(
                        self.metrics
                            .per_targeted_mapping
                            .get(&TargetedMapping {
                                genome: genome.to_string(),
                                region: MappingRegion::Transcriptome,
                                targeting,
                            })
                            .map_or(zero_count, which_metric)
                            .count()
                            .to_string(),
                    );
                }
            }

            // conf <region> reads
            for region in [
                MappingRegion::Intergenic,
                MappingRegion::Exonic,
                MappingRegion::Intronic,
            ] {
                csv_row.push(
                    self.metrics
                        .per_genome_mapping
                        .get(&GenomeMapping {
                            genome: genome.to_string(),
                            region,
                        })
                        .map_or(zero_count, |m| m.conf_mapped)
                        .count()
                        .to_string(),
                );
            }

            // region metrics
            for region in [AnnotationRegion::Exonic, AnnotationRegion::Intronic] {
                let PerGenomeAnnotationRegionMetrics {
                    conf_mapped_unique_region,
                    low_support_umi_region,
                    dup_reads_region,
                    antisense_region,
                    umis_region,
                    ..
                } = self
                    .metrics
                    .per_genome_annotation_region
                    .get(&GenomeAnnotationRegion {
                        genome: genome.to_string(),
                        region,
                    })
                    .cloned()
                    .unwrap_or_else(Metric::new);

                csv_row.push(conf_mapped_unique_region.count().to_string());
                csv_row.push(antisense_region.count().to_string());
                csv_row.push(low_support_umi_region.count().to_string());
                csv_row.push(dup_reads_region.count().to_string());
                csv_row.push(umis_region.count().to_string());
            }
            csv_row.push(
                self.metrics
                    .per_genome_name
                    .get(genome)
                    .map_or(zero_count, |f| f.txomic_unique)
                    .count()
                    .to_string(),
            );
        }
        // insert sizes
        let insert_size_string = format!("{:?}", self.metrics.insert_sizes.mean());
        csv_row.push(insert_size_string);
        stream
            .write_all(format!("{}\n", csv_row.join(",")).as_bytes())
            .unwrap();
    }
}
/// Metrics are sorted by library features and then by barcode.
pub struct LibFeatThenBarcodeOrder;
impl shardio::SortKey<BarcodeMetrics> for LibFeatThenBarcodeOrder {
    type Key = (LibraryFeatures, BarcodeKind);
    fn sort_key(barcode_metrics: &BarcodeMetrics) -> Cow<'_, Self::Key> {
        Cow::Owned((barcode_metrics.library_feats, barcode_metrics.barcode))
    }
}

#[derive(Serialize, Deserialize, Clone, Copy, PartialOrd, Ord, PartialEq, Eq, Hash)]
pub enum BarcodeKind {
    Valid(Barcode),
    Invalid,
}

impl From<Barcode> for BarcodeKind {
    fn from(barcode: Barcode) -> BarcodeKind {
        if barcode.is_valid() {
            BarcodeKind::Valid(barcode)
        } else {
            BarcodeKind::Invalid
        }
    }
}

#[derive(Serialize, Deserialize, Clone, Hash, PartialEq, Eq)]
pub struct GenomeAnnotationRegion {
    pub genome: GenomeName,
    pub region: AnnotationRegion,
}

impl GenomeAnnotationRegion {
    fn new(genome: impl Into<GenomeName>, region: impl Into<AnnotationRegion>) -> Self {
        GenomeAnnotationRegion {
            genome: genome.into(),
            region: region.into(),
        }
    }
}

#[derive(Serialize, Deserialize, Clone, Metric)]
pub struct PerGenomeNameMetrics {
    antisense: CountMetric,
    txomic_unique: CountMetric,
    /// Total UMI counts
    pub umi_counts: CountMetric,
    /// Total reads that contributed to UMI counts.
    pub usable_reads: CountMetric,
}

#[derive(Serialize, Deserialize, Clone, Metric)]
struct PerGenomeAnnotationRegionMetrics {
    conf_mapped_unique_region: CountMetric,
    low_support_umi_region: CountMetric,
    filtered_target_umi_region: CountMetric,
    dup_reads_region: CountMetric,
    antisense_region: CountMetric,
    umis_region: CountMetric,
}

#[derive(Serialize, Deserialize, Clone, Metric)]
pub struct PerMappingStatusMetrics {
    /// Mapped reads
    mapped: CountMetric,

    /// Confidently mapped reads
    conf_mapped: CountMetric,

    /// Confidently mapped reads with a valid barcode
    conf_mapped_bc: CountMetric,

    /// Half-mapped probe reads (MAPQ = 1)
    pub half_mapped: CountMetric,

    /// Split-mapped probe reads (MAPQ = 3)
    pub split_mapped: CountMetric,
}

#[derive(Serialize, Deserialize, Clone, Metric)]
pub struct VisitorMetrics {
    /// Number of sequenced reads where the raw barcode was corrected
    pub barcode_corrected_sequenced_reads: CountMetric,
    /// Number of reads including both discarded reads (discarded_reads) and
    /// the reads processed by the aligner (total_aligner_reads). See comments
    /// above MakeShardMetrics.sequenced_reads for more information
    pub sequenced_reads: CountMetric,
    /// Reads discarded by subsampling.
    discarded_reads: CountMetric,
    /// Reads not mapped to any of the genomes
    pub unmapped_reads: CountMetric,
    /// Total number of reads processed by the aligner
    total_aligner_reads: CountMetric,
    /// Total number or reads with a valid barcode where the UMI was corrected
    /// TODO: Not sure if we are intentionally conditioning on reads with valid barcodes
    corrected_umi_reads: CountMetric,
    /// Total number of reads with a valid barcode where the UMI was marked as low support
    low_support_umi_reads: CountMetric,
    /// Total number of reads with a valid barcode, conf. mapped to an on-target gene,
    /// but not counted because the UMI did not have high enough read support. Applies
    /// only to Targeted Gene Expression data
    filtered_target_umi_reads: CountMetric,
    candidate_dup_reads: CountMetric,
    /// Candidate dup reads which are not used for UMI counts
    dup_reads: CountMetric,
    /// Feature barcode metrics
    corrected_feature_bc_reads: CountMetric,
    extracted_feature_bc_reads: CountMetric,
    invalid_feature_bc_reads: CountMetric,
    total_feature_bc_reads: CountMetric,
    // Per-barcode insert sizes
    insert_sizes: MeanMetric,
    /// Number of reads that match the template-switching oligo (TSO) sequence.
    tso_reads: CountMetric,
    pub per_genome_name: FxHashMap<GenomeName, PerGenomeNameMetrics>,
    // Additional metrics split by region needed for per barcode metrics
    per_genome_annotation_region:
        FxHashMap<GenomeAnnotationRegion, PerGenomeAnnotationRegionMetrics>,
    pub per_genome_mapping: FxHashMap<GenomeMapping, PerMappingStatusMetrics>,
    /// Additional metrics split by targeting status
    per_targeted_mapping: FxHashMap<TargetedMapping, PerMappingStatusMetrics>,
    /// Fraction of reads with a valid UMI among the reads that have a valid barcode
    good_umi_with_good_bc: PercentMetric,
    /// Total UMI counts
    pub umi_counts: i64,
    /// Total reads that contributed to UMI counts.
    pub usable_reads: i64,
}

impl VisitorMetrics {
    fn increment_mapped(&mut self, genome_mapping: GenomeMapping, mapq: u8) {
        let metric = self
            .per_genome_mapping
            .entry(genome_mapping)
            .or_insert_with(Metric::new);
        metric.mapped.increment();

        match mapq {
            MAPQ_HALF_MAPPED => metric.half_mapped.increment(),
            MAPQ_SPLIT_MAPPED => metric.split_mapped.increment(),
            _ => (),
        };
    }

    fn increment_conf_mapped(&mut self, barcode: Barcode, genome_mapping: GenomeMapping) {
        let entry = self
            .per_genome_mapping
            .entry(genome_mapping)
            .or_insert_with(Metric::new);
        if barcode.is_valid() {
            entry.conf_mapped_bc.increment();
        }
        entry.conf_mapped.increment();
    }

    fn increment_mapped_targeting(&mut self, targeted_mapping: TargetedMapping, mapq: u8) {
        let metric = self
            .per_targeted_mapping
            .entry(targeted_mapping)
            .or_insert_with(Metric::new);
        metric.mapped.increment();

        match mapq {
            MAPQ_HALF_MAPPED => metric.half_mapped.increment(),
            MAPQ_SPLIT_MAPPED => metric.split_mapped.increment(),
            _ => (),
        }
    }

    fn increment_conf_mapped_targeting(
        &mut self,
        barcode: Barcode,
        targeted_mapping: TargetedMapping,
    ) {
        let entry = self
            .per_targeted_mapping
            .entry(targeted_mapping)
            .or_insert_with(Metric::new);
        if barcode.is_valid() {
            entry.conf_mapped_bc.increment();
        }
        entry.conf_mapped.increment();
    }

    pub fn make_report(
        self,
        library_feats: LibraryFeatures,
        genomes: &[GenomeName],
        has_targeted_gex: bool,
    ) -> JsonReporter {
        match library_feats {
            LibraryFeatures::GeneExpression(_) => {
                AlignCountReport::Gex(GexReport::new(self, genomes, has_targeted_gex))
                    .to_json_reporter()
            }
            LibraryFeatures::FeatureBarcodes(_) => {
                AlignCountReport::Feature(FeatureReport::from(self)).to_json_reporter()
            }
            _ => unreachable!(),
        }
    }
}

impl AnnotatedReadVisitor for AlignAndCountVisitor {
    fn visit_mapped_read_genomes(&mut self, ann: &ReadAnnotations, genomes: HashSet<&GenomeName>) {
        if !genomes.is_empty() {
            self.metrics
                .increment_mapped(GenomeMapping::multi(Genome), ann.mapq());
        }
        for genome in genomes {
            self.metrics
                .increment_mapped(GenomeMapping::new(genome, Genome), ann.mapq());
        }
    }

    fn visit_mapped_read_genes(
        &mut self,
        ann: &ReadAnnotations,
        genes: HashSet<(&GenomeName, &Gene)>,
    ) {
        let genomes: HashSet<_> = genes.iter().copied().map(|(g, _)| g).collect();
        if !genomes.is_empty() {
            self.metrics
                .increment_mapped(GenomeMapping::multi(Transcriptome), ann.mapq())
        }
        for genome in &genomes {
            self.metrics
                .increment_mapped(GenomeMapping::new(*genome, Transcriptome), ann.mapq())
        }
        if let Some(target_genes) = &self.target_genes {
            // If mapped to more than one gene, count on-target if at least one gene is targeted
            let is_targeted = genes.iter().any(|(_, g)| target_genes.contains(g));
            let targeting = if is_targeted {
                TargetingStatus::Targeted
            } else {
                TargetingStatus::Untargeted
            };
            if !genomes.is_empty() {
                self.metrics.increment_mapped_targeting(
                    TargetedMapping::multi(Transcriptome, targeting),
                    ann.mapq(),
                );
            };
            for genome in genomes {
                self.metrics.increment_mapped_targeting(
                    TargetedMapping::new(genome, Transcriptome, targeting),
                    ann.mapq(),
                );
            }
        }
    }

    fn visit_mapped_read_regions(
        &mut self,
        ann: &ReadAnnotations,
        genome_regions: HashSet<(&GenomeName, AnnotationRegion)>,
    ) {
        for region in genome_regions.iter().map(|&(_g, r)| r).unique() {
            self.metrics
                .increment_mapped(GenomeMapping::multi(region), ann.mapq());
        }
        for (genome, region) in genome_regions {
            self.metrics
                .increment_mapped(GenomeMapping::new(genome, region), ann.mapq());
        }
    }

    fn visit_conf_mapped_read_genome(
        &mut self,
        ReadAnnotations { read, .. }: &ReadAnnotations,
        genome: &GenomeName,
    ) {
        self.metrics
            .increment_conf_mapped(read.barcode(), GenomeMapping::multi(Genome));
        self.metrics
            .increment_conf_mapped(read.barcode(), GenomeMapping::new(genome, Genome));
    }

    fn visit_conf_mapped_read_gene(
        &mut self,
        annotation: &ReadAnnotations,
        genome: &GenomeName,
        gene: &Gene,
    ) {
        self.metrics.increment_conf_mapped(
            annotation.read.barcode(),
            GenomeMapping::multi(Transcriptome),
        );
        self.metrics.increment_conf_mapped(
            annotation.read.barcode(),
            GenomeMapping::new(genome, Transcriptome),
        );
        if let Some(target_genes) = &self.target_genes {
            let targeting = if target_genes.contains(gene) {
                TargetingStatus::Targeted
            } else {
                TargetingStatus::Untargeted
            };
            self.metrics.increment_conf_mapped_targeting(
                annotation.read.barcode(),
                TargetedMapping::multi(Transcriptome, targeting),
            );
            self.metrics.increment_conf_mapped_targeting(
                annotation.read.barcode(),
                TargetedMapping::new(genome, Transcriptome, targeting),
            );
        }
        if let Some((_, region)) = annotation.primary.conf_mapped_region() {
            self.metrics
                .per_genome_annotation_region
                .entry(GenomeAnnotationRegion::new(genome, region))
                .or_insert_with(Metric::new)
                .conf_mapped_unique_region
                .increment();
        }
        if annotation.primary.is_conf_mapped_unique_txomic() {
            for g in [genome, MULTI_GENOME] {
                self.metrics
                    .per_genome_name
                    .entry(g.to_string())
                    .or_insert_with(Metric::new)
                    .txomic_unique
                    .increment();
            }
        }
        if let Some(insert_size) = annotation.insert_size() {
            self.metrics.insert_sizes.record(insert_size as f64);
        }
    }

    fn visit_conf_mapped_read_region(
        &mut self,
        ReadAnnotations { read, .. }: &ReadAnnotations,
        genome: &GenomeName,
        region: AnnotationRegion,
    ) {
        self.metrics
            .increment_conf_mapped(read.barcode(), GenomeMapping::multi(region));
        self.metrics
            .increment_conf_mapped(read.barcode(), GenomeMapping::new(genome, region));
    }

    fn visit_antisense_conf_mapped_read(&mut self, annotation: &ReadAnnotations) {
        let genome = annotation.conf_mapped_genome().unwrap();
        let region = match annotation.conf_mapped_region() {
            Some(reg) => reg.1,
            // it is possible for conf mapped paired end reads to map to different genomes
            // conf_mapped_region() may return None in that case
            None => {
                return;
            }
        };
        for g in [genome, MULTI_GENOME] {
            self.metrics
                .per_genome_name
                .entry(g.to_string())
                .or_insert_with(Metric::new)
                .antisense
                .increment();
        }
        self.metrics
            .per_genome_annotation_region
            .entry(GenomeAnnotationRegion::new(genome, region))
            .or_insert_with(Metric::new)
            .antisense_region
            .increment();
    }

    fn visit_every_read(&mut self, read: &RnaRead) {
        self.metrics.sequenced_reads.increment();
        if read.barcode.is_valid_and_at_least_one_segment_corrected() {
            self.metrics.barcode_corrected_sequenced_reads.increment();
        }
    }

    fn visit_discarded_read(&mut self, _: &RnaRead) {
        self.metrics.discarded_reads.increment()
    }

    fn visit_non_discarded_read(&mut self, _read: &RnaRead) {
        self.metrics.total_aligner_reads.increment();
    }

    fn visit_read_annotation(&mut self, annotation: &ReadAnnotations) {
        let barcode = annotation.read.barcode();
        if barcode.is_valid() {
            // last_barcode is initialized as None
            if self.last_barcode.map_or(true, |bc| bc != barcode) {
                // Add a BarcodeSummary element for the new barcode
                self.barcode_summaries.push(BarcodeSummary::new(
                    barcode,
                    self.library_feats.legacy_library_type(),
                ));
            }
            self.barcode_summaries
                .last_mut()
                .unwrap()
                .observe(annotation);
            self.metrics
                .good_umi_with_good_bc
                .increment(annotation.umi_info.is_valid);
        }

        match self.last_barcode {
            Some(bc) if BarcodeKind::from(bc) != BarcodeKind::from(barcode) => {
                self.metrics_sender
                    .send(BarcodeMetrics {
                        barcode: bc.into(),
                        library_feats: self.library_feats,
                        metrics: self.metrics.clone(),
                    })
                    .expect("Error while sending metrics");
                self.metrics = VisitorMetrics::new();
            }
            _ => {}
        }
        self.last_barcode = Some(barcode);

        debug_assert!(annotation.read.library_feats() == self.library_feats);
        tx_annotation::visitor::walk_read_annotation(self, annotation);
        self.metrics.tso_reads += annotation.matched_tso.into();
    }

    fn visit_feature_read(&mut self, annotation: &ReadAnnotations) {
        self.metrics.total_feature_bc_reads.increment();
        if annotation.is_feature_extracted() {
            self.metrics.extracted_feature_bc_reads.increment();
            if annotation.is_feature_corrected() {
                self.metrics.corrected_feature_bc_reads.increment();
            } else if annotation.is_feature_invalid() {
                self.metrics.invalid_feature_bc_reads.increment();
            }
        }
    }

    fn visit_unmapped_read(&mut self, _: &ReadAnnotations) {
        self.metrics.unmapped_reads.increment()
    }

    fn visit_dup_info(&mut self, annotation: &ReadAnnotations, dup_info: &DupInfo) {
        self.metrics.corrected_umi_reads += dup_info.is_corrected.into();
        let genome = annotation.conf_mapped_genome();

        if let Some(umi_count) = dup_info.umi_count() {
            if let Some(genome) = genome {
                for g in [genome, MULTI_GENOME] {
                    let entry = self
                        .metrics
                        .per_genome_name
                        .entry(g.to_string())
                        .or_insert_with(Metric::new);
                    entry.umi_counts.increment();
                    entry.usable_reads.increment_by(umi_count.read_count);
                }
            }
            self.metrics.umi_counts += 1;
            self.metrics.usable_reads += umi_count.read_count as i64;
        }

        let region = annotation.conf_mapped_region();
        let is_mapped_to_region = genome.is_some() && region.is_some();
        if dup_info.is_low_support_umi {
            self.metrics.low_support_umi_reads.increment();
            if is_mapped_to_region {
                self.metrics
                    .per_genome_annotation_region
                    .entry(GenomeAnnotationRegion::new(
                        genome.unwrap(),
                        region.unwrap().1,
                    ))
                    .or_insert_with(Metric::new)
                    .low_support_umi_region
                    .increment();
            }
        } else if dup_info.is_filtered_target_umi {
            self.metrics.filtered_target_umi_reads.increment();
            if is_mapped_to_region {
                self.metrics
                    .per_genome_annotation_region
                    .entry(GenomeAnnotationRegion::new(
                        genome.unwrap(),
                        region.unwrap().1,
                    ))
                    .or_insert_with(Metric::new)
                    .filtered_target_umi_region
                    .increment();
            }
        } else {
            self.metrics.candidate_dup_reads.increment();
            self.metrics.dup_reads += (!dup_info.is_umi_count()).into();
            if is_mapped_to_region {
                if dup_info.is_umi_count() {
                    self.metrics
                        .per_genome_annotation_region
                        .entry(GenomeAnnotationRegion::new(
                            genome.unwrap(),
                            region.unwrap().1,
                        ))
                        .or_insert_with(Metric::new)
                        .umis_region
                        .increment();
                } else {
                    self.metrics
                        .per_genome_annotation_region
                        .entry(GenomeAnnotationRegion::new(
                            genome.unwrap(),
                            region.unwrap().1,
                        ))
                        .or_insert_with(Metric::new)
                        .dup_reads_region
                        .increment();
                }
            }
        }
    }
}

#[derive(Serialize, Deserialize, Clone, Metric, JsonReport)]
struct MappingMetrics {
    antisense_reads: PercentMetric,
}

#[derive(Serialize, Deserialize, Clone, Metric, JsonReport)]
struct RegionMetrics {
    /// Mapped reads
    mapped_reads: PercentMetric,

    // Confidently mapped reads
    conf_mapped_reads: PercentMetric,

    // Confidently mapped reads with a valid barcode
    conf_mapped_barcoded_reads: PercentMetric,

    /// Half-mapped probe reads (MAPQ = 1)
    half_mapped_reads: PercentMetric,

    /// Split-mapped probe reads (MAPQ = 3)
    split_mapped_reads: PercentMetric,
}

impl RegionMetrics {
    fn with_total_and_mapping(
        total: CountMetric,
        PerMappingStatusMetrics {
            mapped,
            conf_mapped,
            conf_mapped_bc,
            half_mapped,
            split_mapped,
        }: PerMappingStatusMetrics,
    ) -> Self {
        RegionMetrics {
            mapped_reads: (mapped, total).into(),
            conf_mapped_reads: (conf_mapped, total).into(),
            conf_mapped_barcoded_reads: (conf_mapped_bc, total).into(),
            half_mapped_reads: (half_mapped, total).into(),
            split_mapped_reads: (split_mapped, total).into(),
        }
    }

    fn with_total(total: CountMetric) -> RegionMetrics {
        RegionMetrics::with_total_and_mapping(total, Metric::new())
    }
}

/// Metrics which need to be reported for both gene expression
/// and feature barcoding
#[derive(Serialize, Deserialize, Clone, JsonReport)]
pub struct CommonReport {
    /// Number of analyzed reads, not including reads discarded by spaceranger specific subsampling
    /// (See `aligner_subsample_rate`). See the comment above
    /// `make_shard_metrics::MakeShardMetrics.sequenced_reads` for more info
    total_read_pairs: i64,

    /// - Numerator: Total number of reads with a valid barcode where the UMI was corrected
    /// - Denominator: Total number of reads in the gene expression library
    corrected_umi: PercentMetric,

    /// - Numerator: Total number of reads with a valid barcode where the UMI was marked as low support
    /// - Denominator: candidate_dup_reads
    low_support_umi_reads: PercentMetric,

    /// - Numerator: Total number of dup reads
    /// - Denominator: candidate_dup_reads
    multi_cdna_pcr_dupe_reads: PercentMetric,

    /// - Numerator: Number of read with a valid UMI among the reads with a valid barcode
    /// - Denominator: Number of reads with a valid barcode
    good_umi_with_good_bc: PercentMetric,
}

impl From<&VisitorMetrics> for CommonReport {
    fn from(visitor: &VisitorMetrics) -> Self {
        CommonReport {
            total_read_pairs: visitor.total_aligner_reads.count(),
            corrected_umi: (visitor.corrected_umi_reads, visitor.total_aligner_reads).into(),
            low_support_umi_reads: (
                visitor.low_support_umi_reads,
                visitor.candidate_dup_reads + visitor.low_support_umi_reads,
            )
                .into(),
            multi_cdna_pcr_dupe_reads: (visitor.dup_reads, visitor.candidate_dup_reads).into(),
            good_umi_with_good_bc: visitor.good_umi_with_good_bc,
        }
    }
}

/// Metrics which are reported for gene expression
#[derive(Serialize, Deserialize, Clone, JsonReport)]
pub struct GexReport {
    #[json_report(inline)]
    common: CommonReport,
    #[json_report(inline)]
    targeted: Option<TargetedGexReport>,
    #[json_report(inline)]
    antisense: FxHashMap<GenomeName, MappingMetrics>,
    #[json_report(inline)]
    mapping: FxHashMap<GenomeMapping, RegionMetrics>,
    unmapped_reads: PercentMetric,
    /// - Numerator: Number of reads that match the template-switching oligo (TSO) sequence.
    /// - Denominator: Total number of reads in the gene expression library.
    tso: PercentMetric,
    /// Number of reads discarded by subsampling.
    discarded_reads: CountMetric,
    /// Fraction of subsampled reads.
    /// - Numerator: total_reads
    /// - Denominator: total_reads + discarded_reads
    subsampled: PercentMetric,
}

impl GexReport {
    fn new(mut visitor: VisitorMetrics, genomes: &[GenomeName], is_targeted: bool) -> Self {
        let common = CommonReport::from(&visitor);
        let targeted = if is_targeted {
            Some(TargetedGexReport::from_visitor(&mut visitor, genomes))
        } else {
            None
        };

        let den = visitor.total_aligner_reads; // denominator

        let mut antisense: FxHashMap<GenomeName, MappingMetrics> = visitor
            .per_genome_name
            .into_iter()
            .map(|(k, PerGenomeNameMetrics { antisense, .. })| {
                (
                    k,
                    MappingMetrics {
                        antisense_reads: (antisense, den).into(),
                    },
                )
            })
            .collect();

        let mut mapping: FxHashMap<_, _> = visitor
            .per_genome_mapping
            .into_iter()
            .map(|(genome_mapping, per_mapping_status_metrics)| {
                (
                    genome_mapping,
                    RegionMetrics::with_total_and_mapping(den, per_mapping_status_metrics),
                )
            })
            .collect();

        // It is possible that visitor do not see all the genomes for all
        // metrics. Explicitly set the missing metrics to zero.
        for g in genomes
            .iter()
            .chain(std::iter::once(&MULTI_GENOME.to_string()))
        {
            antisense.entry(g.to_string()).or_insert(MappingMetrics {
                antisense_reads: (0, den.count()).into(),
            });
            for m in MappingRegion::iter() {
                mapping
                    .entry(GenomeMapping::new(g, m))
                    .or_insert_with(|| RegionMetrics::with_total(den));
            }
        }

        GexReport {
            common,
            targeted,
            antisense,
            mapping,
            discarded_reads: visitor.discarded_reads,
            subsampled: (
                visitor.total_aligner_reads,
                visitor.total_aligner_reads + visitor.discarded_reads,
            )
                .into(),
            unmapped_reads: (visitor.unmapped_reads, visitor.total_aligner_reads).into(),
            tso: (visitor.tso_reads, visitor.total_aligner_reads).into(),
        }
    }
}

/// Metrics which are reported for targeted gene expression
#[derive(Serialize, Deserialize, Clone, JsonReport)]
pub struct TargetedGexReport {
    #[json_report(inline)]
    mapping: FxHashMap<TargetedMapping, RegionMetrics>,
    /// - Numerator: Total number of reads with a valid barcode, conf. mapped to an on-target gene,
    ///              but not counted because the UMI did not have high enough read support. Applies
    ///              only to Targeted Gene Expression data
    /// - Denominator: total_reads
    filtered_target_umi_reads: PercentMetric,
}

impl TargetedGexReport {
    /// Construct a TargetedGexReport from a VisitorMetrics.
    /// This method consumes (drains) its *_targeting member variables.
    fn from_visitor(visitor: &mut VisitorMetrics, genomes: &[GenomeName]) -> Self {
        let den = visitor.total_aligner_reads; // denominator

        let mut mapping: FxHashMap<_, _> = visitor
            .per_targeted_mapping
            .drain()
            .map(|(genome_mapping, per_mapping_status_metrics)| {
                (
                    genome_mapping,
                    RegionMetrics::with_total_and_mapping(den, per_mapping_status_metrics),
                )
            })
            .collect();

        // Same as in GexReport, set the missing metrics to zero in case
        // we do not observe all genomes.
        for g in genomes
            .iter()
            .chain(std::iter::once(&MULTI_GENOME.to_string()))
        {
            for ts in TargetingStatus::iter() {
                mapping
                    .entry(TargetedMapping::new(g, MappingRegion::Transcriptome, ts))
                    .or_insert_with(|| RegionMetrics::with_total(den));
            }
        }

        TargetedGexReport {
            mapping,
            filtered_target_umi_reads: (
                visitor.filtered_target_umi_reads,
                visitor.total_aligner_reads,
            )
                .into(),
        }
    }
}

/// Things which are reported when we have Feature Barcoding
#[derive(Serialize, Deserialize, Clone, JsonReport)]
pub struct FeatureReport {
    corrected_feature_bc: PercentMetric,
    feature_bc_extracted: PercentMetric,
    recognized_feature_bc: PercentMetric,
    unrecognized_feature_bc: PercentMetric,
    #[json_report(inline)]
    common: CommonReport,
}

impl From<VisitorMetrics> for FeatureReport {
    fn from(visitor: VisitorMetrics) -> Self {
        let corrected = visitor.corrected_feature_bc_reads.count();
        let extracted = visitor.extracted_feature_bc_reads.count();
        let invalid = visitor.invalid_feature_bc_reads.count();
        let valid = extracted - invalid;
        let total = visitor.total_feature_bc_reads.count();
        FeatureReport {
            corrected_feature_bc: (corrected, extracted).into(),
            feature_bc_extracted: (extracted, total).into(),
            recognized_feature_bc: (valid, total).into(),
            unrecognized_feature_bc: (invalid, extracted).into(),
            common: CommonReport::from(&visitor),
        }
    }
}

#[derive(Serialize, Deserialize, Clone)]
enum AlignCountReport {
    Gex(GexReport),
    Feature(FeatureReport),
}

impl JsonReport for AlignCountReport {
    fn to_json_reporter(&self) -> JsonReporter {
        match self {
            AlignCountReport::Gex(ref gex_report) => gex_report.to_json_reporter(),
            AlignCountReport::Feature(ref feature_report) => feature_report.to_json_reporter(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use pretty_assertions::assert_eq;

    #[test]
    fn test_default_keys() {
        let mut visitor_metrics = VisitorMetrics::new();
        visitor_metrics.total_aligner_reads = 100.into();
        visitor_metrics.candidate_dup_reads = 1.into();
        visitor_metrics.good_umi_with_good_bc = (45, 90).into();
        let genomes = ["multi", "GRCh38", "mm10"].map(String::from).to_vec();
        // Note: barnyard + targeting not fully supported, but VisitorMetrics works
        let has_targeted_gex = true;

        let report =
            visitor_metrics.make_report(LibraryFeatures::gex(), &genomes, has_targeted_gex);
        let actual = report.to_json_reporter();
        let mut expected = JsonReporter::new();
        expected.insert("corrected_umi_frac", 0.0f64);
        expected.insert("discarded_reads_count", 0);
        expected.insert("low_support_umi_reads_frac", 0.0f64);
        expected.insert("filtered_target_umi_reads_frac", 0.0f64);
        expected.insert("multi_cdna_pcr_dupe_reads_frac", 0.0f64);
        expected.insert("subsampled_frac", 1.0f64);
        expected.insert("total_read_pairs", 100);
        expected.insert("tso_frac", 0.0f64);
        expected.insert("unmapped_reads_frac", 0.0f64);
        expected.insert("good_umi_with_good_bc_frac", 0.5f64);
        for g in genomes {
            expected.insert(format!("{g}_antisense_reads_frac"), 0.0f64);
            for m in MappingRegion::iter() {
                let gm = GenomeMapping::new(g.to_string(), m);
                expected.insert(format!("{gm}_mapped_reads_frac"), 0.0f64);
                expected.insert(format!("{gm}_conf_mapped_reads_frac"), 0.0f64);
                expected.insert(format!("{gm}_conf_mapped_barcoded_reads_frac"), 0.0f64);
                expected.insert(format!("{gm}_half_mapped_reads_frac"), 0.0f64);
                expected.insert(format!("{gm}_split_mapped_reads_frac"), 0.0f64);
            }
            for ts in TargetingStatus::iter() {
                let tm = TargetedMapping::new(g.to_string(), MappingRegion::Transcriptome, ts);
                expected.insert(format!("{tm}_mapped_reads_frac"), 0.0f64);
                expected.insert(format!("{tm}_conf_mapped_reads_frac"), 0.0f64);
                expected.insert(format!("{tm}_conf_mapped_barcoded_reads_frac"), 0.0f64);
                expected.insert(format!("{tm}_half_mapped_reads_frac"), 0.0f64);
                expected.insert(format!("{tm}_split_mapped_reads_frac"), 0.0f64);
            }
        }
        assert_eq!(actual, expected);
    }
}
