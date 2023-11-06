// There are a **lot** of metrics which need to be imported here.
#![allow(clippy::wildcard_imports)]
use super::tables::*;
use crate::alert::AlertLevel;
use crate::multi::svg::SvgGraph;
use crate::{
    Alert, AlertContext, AlertSpec, CardWithMetric, ChartWithHelp, GenericTable, MakePretty,
    MetricCard, PlotlyChart, RawChartWithHelp, Tab, TableRow, TitleWithHelp, WsSample,
};
use anyhow::Result;
use cr_types::rna_read::LegacyLibraryType;
use cr_types::{AlignerParam, TargetingMethod};
use csv::Writer;
use itertools::Itertools;
use metric::TxHashMap;
use serde::{Deserialize, Serialize};
use serde_json::value::Value;
use std::path::Path;
use websummary_derive::{Alert, ToCsvRows};

/// The threshold to trigger a web summary alert when an unexpected probe barcode is observed or
/// an expected probe barcode is not observed.
pub const UMI_PER_PROBE_BARCODE_BACKGROUND_THRESHOLD: f64 = 0.005;

const METRICS_SUMMARY_CSV_HEADER: [&str; 6] = [
    "Category",
    "Library Type",
    "Grouped By",
    "Group Name",
    "Metric Name",
    "Metric Value",
];

// websummary structs implementing this trait know how to convert their contents into a metrics CSV row
// the metrics CSV is represented as a Vec of CSV rows, each a Vec of String
pub trait ToCsvRows {
    fn to_csv_rows(self) -> Vec<Vec<String>>
    where
        Self: Sized,
    {
        vec![]
    }
}

impl<T> ToCsvRows for Option<T>
where
    T: ToCsvRows,
{
    fn to_csv_rows(self) -> Vec<Vec<String>> {
        match self {
            Some(t) => t.to_csv_rows(),
            None => vec![],
        }
    }
}

#[derive(Debug, Serialize, PartialEq, Eq, Clone, Default)]
pub struct SampleHeaderInfo {
    #[serde(rename = "Sample ID")]
    pub sample_id: String,
    #[serde(rename = "Sample Description")]
    pub sample_desc: String,
    #[serde(rename = "Pipeline Version")]
    pub pipeline_version: String,
}

#[derive(Debug, Serialize, PartialEq, Eq, Clone, Default)]
pub struct LibraryHeaderInfo {
    #[serde(rename = "Run ID")]
    pub run_id: String,
    #[serde(rename = "Run Description")]
    pub run_desc: String,
    #[serde(rename = "Pipeline Version")]
    pub pipeline_version: String,
}

#[derive(Debug, Serialize, PartialEq, Eq, Clone)]
#[serde(into = "GenericTable")]
pub struct CountParametersTable {
    pub chemistry: String,
    pub introns_included: bool,
    pub reference_path: String,
    pub transcriptome: String,
    pub feature_ref_path: Option<String>,
    pub cmo_set_path: Option<String>,
    pub target_set_name: Option<String>,
    pub targeting_method: Option<TargetingMethod>,
    pub filter_probes: Option<bool>,
    pub num_genes_on_target: Option<usize>,
    pub library_type: LegacyLibraryType,
    pub throughput: Option<String>,
    pub tenx_cmos: Option<bool>,
    pub aligner: AlignerParam,
    pub antigen_negative_control: bool,
    // Things we may need to alert on.
    pub dropped_tags: Vec<String>,
    pub probe_barcodes_high_gem_overlap: Vec<String>,
    pub unspecified_probe_barcodes_detected: Vec<String>,
    pub specified_probe_barcodes_missing: Vec<String>,
    pub mismatched_probe_barcode_pairings: Option<MismatchedProbeBarcodePairings>,
}

impl From<CountParametersTable> for GenericTable {
    fn from(info: CountParametersTable) -> GenericTable {
        let CountParametersTable {
            chemistry,
            introns_included,
            reference_path,
            transcriptome,
            feature_ref_path,
            cmo_set_path,
            target_set_name,
            targeting_method,
            filter_probes,
            num_genes_on_target,
            library_type,
            throughput,
            aligner,
            antigen_negative_control,
            // Supporting data that should not appear in the table itself.
            tenx_cmos: _,
            dropped_tags: _,
            probe_barcodes_high_gem_overlap: _,
            unspecified_probe_barcodes_detected: _,
            specified_probe_barcodes_missing: _,
            mismatched_probe_barcode_pairings: _,
        } = info;

        let chemistry_with_throughput = match throughput {
            Some(ref tp) => {
                if tp == "HT" && !["HT", "MT", "LT"].iter().any(|x| chemistry.contains(x)) {
                    format!("{chemistry} HT")
                } else {
                    chemistry
                }
            }
            None => chemistry,
        };

        let mut rows = vec![TableRow::two_col("Chemistry", chemistry_with_throughput)];

        // GEX tab
        if library_type == LegacyLibraryType::GeneExpression {
            rows.extend([
                TableRow::two_col("Reference Path", reference_path),
                TableRow::two_col("Transcriptome", transcriptome),
            ]);
            if aligner != AlignerParam::Hurtle {
                rows.push(TableRow::two_col("Include Introns", introns_included));
            }
            if let Some(targeting_method) = targeting_method {
                rows.push(TableRow::two_col(
                    match targeting_method {
                        TargetingMethod::HybridCapture => "Target Panel Name",
                        TargetingMethod::TemplatedLigation => "Probe Set Name",
                    },
                    target_set_name.unwrap(),
                ));
            }
            if let Some(n_genes) = num_genes_on_target {
                rows.push(TableRow::two_col(
                    "Number of Genes Targeted",
                    n_genes.make_pretty(),
                ));
            }
            if let Some(filter_probes) = filter_probes {
                rows.push(TableRow::two_col(
                    "Filter Probes",
                    if filter_probes { "On" } else { "Off" },
                ));
            }
        }
        // Feature Tabs
        match (library_type, feature_ref_path) {
            (x, Some(feature_ref)) if x != LegacyLibraryType::GeneExpression => {
                rows.push(TableRow::two_col("Feature Reference", feature_ref));
            }
            _ => {}
        }

        if let (Some(cmo_set_path), LegacyLibraryType::Multiplexing) = (cmo_set_path, library_type)
        {
            rows.push(TableRow::two_col("CMO Set", cmo_set_path));
        }

        if library_type == LegacyLibraryType::AntigenCapture {
            rows.push(TableRow::two_col(
                "Control Specified",
                antigen_negative_control,
            ));
        }
        GenericTable { header: None, rows }
    }
}

impl Alert for CountParametersTable {
    fn alerts(&self, ctx: &AlertContext) -> Vec<AlertSpec> {
        let mut alerts = vec![];

        if !ctx.is_hybrid_capture
            && !ctx.is_rtl
            && (self.introns_included || ctx.include_introns)
            && self.library_type == LegacyLibraryType::GeneExpression
        {
            alerts.push(AlertSpec {
                level: AlertLevel::Info,
                title: "Intron mode used".into(),
                formatted_value: String::default(),
                message: r#"This data has been analyzed with intronic reads included in the count matrix. This behavior is different from previous Cell Ranger versions. If you would not like to count intronic reads, please rerun with the "include-introns" option set to "false". Please contact support@10xgenomics.com for any further questions."#.into(),
            });
        }
        if ctx.is_hybrid_capture
            && ctx.include_introns
            && self.library_type == LegacyLibraryType::GeneExpression
        {
            alerts.push(AlertSpec {
                level: AlertLevel::Warn,
                title: "Unsupported workflow used".to_string(),
                formatted_value: String::default(),
                message: "Your data has been analyzed with targeted panel (target-panel) and included introns (include-introns). This is not recommended, because in the Targeted Gene Expression assay, 10x Genomics supported baits are designed to capture exons only. Results cannot be guaranteed.".to_string(),
            });
        }
        if ctx.is_antigen
            && self.library_type == LegacyLibraryType::AntigenCapture
            && !self.antigen_negative_control
        {
            alerts.push(AlertSpec {
                level: AlertLevel::Info,
                title: "No antigen negative control specified".to_string(),
                formatted_value: String::default(),
                message: r#"Your data has been analyzed without a negative control. Antigen specificity scores and related outputs will not be generated."#.into(),
            });
        }
        if let Some(tp) = &self.throughput {
            if self.chemistry.ends_with("HT") && tp != "HT" {
                alerts.push(AlertSpec {
                    level: AlertLevel::Warn,
                    title: "Inconsistent throughput detected".to_string(),
                    formatted_value: "HT".to_string(),
                    message: "High throughput (HT) chemistry was specified, but Cell Ranger detected that data was not HT. This may be due to low sequencing depth, sample quality issues, or incorrect chemistry input to Cell Ranger. If applicable, cell multiplexing results may not be accurate.".to_string(),
                });
            }
        }
        if self.tenx_cmos == Some(false) {
            alerts.push(AlertSpec {
                level: AlertLevel::Warn,
                title: "Unsupported multiplexing tag used".to_string(),
                formatted_value: String::default(),
                message: "Multiplexing performance cannot be guaranteed".to_string(),
            });
        }
        if ctx.is_fiveprime && ctx.is_multiplexing {
            alerts.push(AlertSpec {
                level: AlertLevel::Warn,
                title: "Unsupported combination of 5' chemistry with multiplexing".to_string(),
                formatted_value: String::default(),
                message: "Multiplexing performance cannot be guaranteed".to_string(),
            });
        }
        if ctx.is_lt_chemistry && ctx.is_multiplexing {
            alerts.push(AlertSpec {
                level: AlertLevel::Warn,
                title: "Unsupported combination of 3' v3 LT chemistry with multiplexing"
                    .to_string(),
                formatted_value: String::default(),
                message: "Multiplexing performance cannot be guaranteed".to_string(),
            });
        }
        if ctx.is_arc_chemistry {
            alerts.push(AlertSpec {
                level: AlertLevel::Warn,
                title: "Unsupported workflow used".to_string(),
                formatted_value: String::default(),
                message: "Multiome Gene Expression only analysis is not a supported workflow. Results may vary.".to_string(),
            })
        }

        if !self.dropped_tags.is_empty() {
            let filtered_tags = self.dropped_tags.join(", ");
            let message = format!(
                "Tag(s) {filtered_tags} were identified as contaminants.  This will happen if \
            any of the following conditions are met.  The total UMI count for that tag is < 1,000, \
            the sum of UMIs for a tag is less than 2% of the maximum UMI count for all other tags, \
            or the pearson correlation with a tag and another was >0.5 and it had fewer UMIs than the \
            correlated tag. \
            You should check if your input CSV contained the right list of tags used."
            );
            alerts.push(AlertSpec {
                level: AlertLevel::Warn,
                title: "Contaminant tags detected".to_string(),
                formatted_value: String::default(),
                message,
            });
        }
        if !self.probe_barcodes_high_gem_overlap.is_empty() {
            let total_high_overlap_pairs = self.probe_barcodes_high_gem_overlap.len();

            const MAX_DISPLAY_PAIRS: usize = 5;

            let message = if total_high_overlap_pairs <= MAX_DISPLAY_PAIRS {
                let probe_pairs = self.probe_barcodes_high_gem_overlap.join(", ");
                format!(
                        "The following pair(s) of probe barcode IDs were identified as having much \
                        higher than expected GEM barcode overlap: {probe_pairs}. This may be due to cell \
                        clumping or indicate that probes with two or more different probe barcodes \
                        were used in the same hybridization reaction. \
                        See frp_gem_barcode_overlap.csv for a complete set of overlap metrics."
                    )
            } else {
                // Provide a modified alert if the number of pairs is high
                let probe_pairs =
                    self.probe_barcodes_high_gem_overlap[0..MAX_DISPLAY_PAIRS].join(", ");
                let remaining_pairs = total_high_overlap_pairs - MAX_DISPLAY_PAIRS;
                format!(
                        "The following pair(s) of probe barcode IDs were identified as having much \
                        higher than expected GEM barcode overlap: {probe_pairs}... [{remaining_pairs} additional pair(s)]. \
                        This may be due to cell clumping or indicate that probes with two or more \
                        different barcodes were used in the same hybridization reaction. \
                        See frp_gem_barcode_overlap.csv for a complete set of overlap metrics."
                    )
            };

            alerts.push(AlertSpec {
                level: AlertLevel::Error,
                title: "Probe barcodes with high GEM barcode overlap detected".to_string(),
                formatted_value: format!("{total_high_overlap_pairs} probe barcode pair(s)"),
                message,
            });
        }

        if !self.unspecified_probe_barcodes_detected.is_empty() {
            let barcodes_string = self.unspecified_probe_barcodes_detected.join(", ");
            let message = format!(
                    "The following probe barcodes ID(s) were not specified in config CSV but account for at least {:.2}% of UMIs: {}. \
                This could result from omitting one or more probe barcode IDs used in the experiment \
                or from accidental contamination coming from a barcode that was not used in this experiment. \
                Please check your config CSV file.",
                    UMI_PER_PROBE_BARCODE_BACKGROUND_THRESHOLD * 100.0,
                    barcodes_string
                );
            alerts.push(AlertSpec {
                level: AlertLevel::Error,
                title: "Probe barcode ID(s) missing from config CSV".to_string(),
                formatted_value: String::default(),
                message,
            });
        }

        if !self.specified_probe_barcodes_missing.is_empty() {
            let barcodes_string = self.specified_probe_barcodes_missing.join(", ");
            let message = format!(
                    "The following probe barcodes ID(s) were specified in config CSV but account for less than {:.2}% of UMIs: {}. \
            This could result from adding a probe barcode to your config CSV when it was not actually used in the experiment. \
            Please check your config CSV file.",
                    UMI_PER_PROBE_BARCODE_BACKGROUND_THRESHOLD * 100.0,
                    barcodes_string
                );
            alerts.push(AlertSpec {
                level: AlertLevel::Warn,
                title: "Specified probe barcode ID(s) not detected".to_string(),
                formatted_value: String::default(),
                message,
            });
        }

        if let Some(mismatched_probe_bcs) = &self.mismatched_probe_barcode_pairings {
            alerts.push(AlertSpec {
                level: AlertLevel::Warn,
                title: "Mismatch between declared and detected probe barcode pairings".to_string(),
                formatted_value: mismatched_probe_bcs.formatted_value(),
                message: mismatched_probe_bcs.to_string(),
            });
        }

        alerts
    }
}

impl ToCsvRows for CountParametersTable {}

#[derive(Serialize, Deserialize, Clone)]
pub enum GexOrRtl<G, R> {
    Gex(G),
    Rtl(R),
}

impl<G, R> From<GexOrRtl<G, R>> for CardWithMetric
where
    CardWithMetric: From<G>,
    CardWithMetric: From<R>,
{
    fn from(src: GexOrRtl<G, R>) -> CardWithMetric {
        match src {
            GexOrRtl::Gex(g) => g.into(),
            GexOrRtl::Rtl(r) => r.into(),
        }
    }
}

impl<G, R> Alert for GexOrRtl<G, R>
where
    G: Alert,
    R: Alert,
{
    fn alerts(&self, ctx: &AlertContext) -> Vec<AlertSpec> {
        match self {
            GexOrRtl::Gex(g) => g.alerts(ctx),
            GexOrRtl::Rtl(r) => r.alerts(ctx),
        }
    }
}

pub type GexOrRtlLibraryMappingMetricsTable =
    GexOrRtl<GexLibraryMappingMetricsTable, RtlLibraryMappingMetricsTable>;

pub type GexOrRtlSampleMappingMetricsTable =
    GexOrRtl<GexSampleMappingMetricsTable, RtlSampleMappingMetricsTable>;

pub type SampleCellMetricsTable = GexOrRtl<GexSampleCellMetricsTable, RtlSampleCellMetricsTable>;

#[derive(Serialize, Deserialize, Clone)]
pub enum AntibodyOrAntigen<B, G> {
    Antibody(B),
    Antigen(G),
}

impl<B, G> From<AntibodyOrAntigen<B, G>> for CardWithMetric
where
    CardWithMetric: From<B>,
    CardWithMetric: From<G>,
{
    fn from(src: AntibodyOrAntigen<B, G>) -> CardWithMetric {
        match src {
            AntibodyOrAntigen::Antibody(b) => b.into(),
            AntibodyOrAntigen::Antigen(g) => g.into(),
        }
    }
}

impl<B, G> Alert for AntibodyOrAntigen<B, G>
where
    B: Alert,
    G: Alert,
{
    fn alerts(&self, ctx: &AlertContext) -> Vec<AlertSpec> {
        match self {
            AntibodyOrAntigen::Antibody(b) => b.alerts(ctx),
            AntibodyOrAntigen::Antigen(g) => g.alerts(ctx),
        }
    }
}

pub type AntibodyOrAntigenPhysicalLibraryMetricsTable =
    AntibodyOrAntigen<AntibodyPhysicalLibraryMetricsTable, AntigenPhysicalLibraryMetricsTable>;

// Websummary data structures may have shared _resources that they access by key
// If a websummary struct has any _resources they are emptied and bubbled up and stored in the top-level _resources
pub type MultiSharedResource = TxHashMap<String, Value>;

impl Alert for MultiSharedResource {}
impl ToCsvRows for MultiSharedResource {}

#[derive(Debug, Serialize, PartialEq, Eq, Clone, Default)]
#[serde(into = "GenericTable")]
pub struct VdjParametersTable {
    pub chemistry: String,
    pub vdj_reference: String,
    pub vdj_reference_path: String,
    pub gamma_delta: bool,
}

impl From<VdjParametersTable> for GenericTable {
    fn from(info: VdjParametersTable) -> GenericTable {
        let VdjParametersTable {
            chemistry,
            vdj_reference,
            vdj_reference_path,
            gamma_delta: _, // Only used for alert
        } = info;
        let rows = vec![
            TableRow::two_col("Chemistry", chemistry),
            TableRow::two_col("V(D)J Reference", vdj_reference),
            TableRow::two_col("V(D)J Reference Path", vdj_reference_path),
        ];
        GenericTable { header: None, rows }
    }
}

impl Alert for VdjParametersTable {
    fn alerts(&self, _: &AlertContext) -> Vec<AlertSpec> {
        let mut alerts = vec![];
        if self.gamma_delta {
            alerts.push(AlertSpec {
                level: AlertLevel::Warn,
                title: "Unsupported workflow used".to_string(),
                formatted_value: String::default(),
                message: "Gamma Delta TCR analysis is not a supported workflow. Algorithm performance cannot be guaranteed.".to_string(),
            });
        }

        alerts
    }
}
impl ToCsvRows for VdjParametersTable {}

/// A web summary element that is specific for
/// each chain type.
#[derive(Serialize, Deserialize, Clone)]
pub enum VdjChainTypeSpecific<T, B, GD> {
    VdjT(T),    // TRA, TRB
    VdjB(B),    // IGH, IGL, IGK
    VdjTgd(GD), // TRG, TRD
}

impl<T, B, GD> From<VdjChainTypeSpecific<T, B, GD>> for CardWithMetric
where
    CardWithMetric: From<T>,
    CardWithMetric: From<B>,
    CardWithMetric: From<GD>,
{
    fn from(src: VdjChainTypeSpecific<T, B, GD>) -> CardWithMetric {
        match src {
            VdjChainTypeSpecific::VdjT(t) => t.into(),
            VdjChainTypeSpecific::VdjB(b) => b.into(),
            VdjChainTypeSpecific::VdjTgd(t_gd) => t_gd.into(),
        }
    }
}

impl<T, B, GD> Alert for VdjChainTypeSpecific<T, B, GD>
where
    T: Alert,
    B: Alert,
    GD: Alert,
{
    fn alerts(&self, ctx: &AlertContext) -> Vec<AlertSpec> {
        match self {
            VdjChainTypeSpecific::VdjT(t) => t.alerts(ctx),
            VdjChainTypeSpecific::VdjB(b) => b.alerts(ctx),
            VdjChainTypeSpecific::VdjTgd(t_gd) => t_gd.alerts(ctx),
        }
    }
}

pub type VdjEnrichmentMetricsTable = VdjChainTypeSpecific<
    VdjTEnrichmentMetricsTable,
    VdjBEnrichmentMetricsTable,
    VdjTgdEnrichmentMetricsTable,
>;
pub type VdjSampleHeroMetricsTable = VdjChainTypeSpecific<
    VdjTSampleHeroMetricsTable,
    VdjBSampleHeroMetricsTable,
    VdjTgdSampleHeroMetricsTable,
>;
pub type VdjSampleAnnotationMetricsTable = VdjChainTypeSpecific<
    VdjTSampleAnnotationMetricsTable,
    VdjBSampleAnnotationMetricsTable,
    VdjTgdSampleAnnotationMetricsTable,
>;
impl Alert for Value {}
impl ToCsvRows for Value {}

#[derive(Serialize, Clone)]
pub struct MultiWebSummary {
    pub sample: WsSample,
    pub data: MultiWebSummaryData,
    pub diagnostics: MultiDiagnostics,
    #[serde(rename = "_resources")]
    pub resources: MultiSharedResource,
}

#[derive(Serialize, Clone, Default)]
pub struct VdjDiagnostics {
    pub filter_metrics: Option<Value>,
}

#[derive(Serialize, Clone, Default)]
#[allow(non_snake_case)]
pub struct MultiDiagnostics {
    pub corrected_bc_frac: Option<Value>,
    pub corrected_bc_in_gel_bead_frac: Option<Value>,
    pub corrected_bc_in_probe_frac: Option<Value>,
    pub ANTIBODY_corrected_bc_frac: Option<Value>,
    pub ANTIBODY_corrected_bc_in_gel_bead_frac: Option<Value>,
    pub ANTIBODY_corrected_bc_in_probe_frac: Option<Value>,
    pub i1_bases_with_q30_frac: Option<Value>,
    pub i2_bases_with_q30_frac: Option<Value>,
    pub low_support_umi_reads_frac: Option<Value>,
    pub tag_contaminant_info: Option<Value>,
    pub tso_frac: Option<Value>,
    pub probe_barcode_overlap_coefficients: Option<Value>,
    pub fraction_reads_high_occupancy_gems: Option<Value>,
    pub high_occupancy_probe_barcode_count_threshold: Option<Value>,
    pub vdj_t: Option<VdjDiagnostics>,
    pub vdj_b: Option<VdjDiagnostics>,
    pub vdj_t_gd: Option<VdjDiagnostics>,
}

impl MultiWebSummary {
    pub fn to_csv(self, filename: &Path) -> Result<()> {
        let mut writer = Writer::from_path(filename)?;
        writer.write_record(METRICS_SUMMARY_CSV_HEADER.iter())?;
        for row in self.data.to_csv_rows().iter().sorted().dedup() {
            writer.write_record(row)?;
        }

        writer.flush()?;
        Ok(())
    }
}

impl ToCsvRows for MultiWebSummary {
    fn to_csv_rows(self) -> Vec<Vec<String>> {
        self.data.to_csv_rows()
    }
}

#[derive(Serialize, Clone)]
pub struct MultiWebSummaryData {
    pub library_websummary: LibraryWebSummary,
    pub sample_websummary: SampleWebSummary,
    pub experimental_design: ExperimentalDesign,
}

impl ToCsvRows for MultiWebSummaryData {
    fn to_csv_rows(self) -> Vec<Vec<String>> {
        let mut rows = self.library_websummary.to_csv_rows();
        rows.append(&mut self.sample_websummary.to_csv_rows());
        rows
    }
}

#[derive(Serialize, Clone)]
pub struct ExperimentalDesign {
    pub svg: SvgGraph,
    pub csv: String,
}

fn to_csv_rows_helper<T: ToCsvRows>(
    tab: Option<T>,
    label: &str,
    sample_or_library: &str,
) -> Vec<Vec<String>> {
    let mut rows = vec![];
    if let Some(t) = tab {
        let mut tab_rows = t
            .to_csv_rows()
            .iter()
            .cloned()
            .map(|mut r| {
                r.insert(0, label.to_string());
                r.insert(0, sample_or_library.to_string());
                r
            })
            .collect();
        rows.append(&mut tab_rows);
    }
    rows
}

#[derive(Serialize, Clone, Default)]
pub struct LibraryWebSummary {
    pub header_info: LibraryHeaderInfo,
    pub gex_tab: Option<Tab<LibraryGexWebSummary>>,
    pub vdj_t_tab: Option<Tab<LibraryVdjWebSummary>>,
    pub vdj_t_gd_tab: Option<Tab<LibraryVdjWebSummary>>,
    pub vdj_b_tab: Option<Tab<LibraryVdjWebSummary>>,
    pub antibody_tab: Option<Tab<LibraryAntibodyOrAntigenWebSummary>>,
    pub antigen_tab: Option<Tab<LibraryAntibodyOrAntigenWebSummary>>,
    pub crispr_tab: Option<Tab<LibraryCrisprWebSummary>>,
    pub custom_feature_tab: Option<Tab<LibraryCustomFeatureWebSummary>>,
    pub cmo_tab: Option<Tab<LibraryCmoWebSummary>>,
    #[serde(rename = "_resources")]
    pub resources: MultiSharedResource,
}

impl ToCsvRows for LibraryWebSummary {
    fn to_csv_rows(self) -> Vec<Vec<String>> {
        [
            to_csv_rows_helper(self.gex_tab, "Gene Expression", "Library"),
            to_csv_rows_helper(self.vdj_t_tab, "VDJ T", "Library"),
            to_csv_rows_helper(self.vdj_t_gd_tab, "VDJ T GD", "Library"),
            to_csv_rows_helper(self.vdj_b_tab, "VDJ B", "Library"),
            to_csv_rows_helper(self.antibody_tab, "Antibody Capture", "Library"),
            to_csv_rows_helper(self.antigen_tab, "Antigen Capture", "Library"),
            to_csv_rows_helper(self.crispr_tab, "CRISPR Guide Capture", "Library"),
            to_csv_rows_helper(self.custom_feature_tab, "Custom Feature", "Library"),
            to_csv_rows_helper(self.cmo_tab, "Multiplexing Capture", "Library"),
        ]
        .concat()
    }
}

#[derive(Serialize, Clone, Default)]
pub struct SampleWebSummary {
    pub header_info: SampleHeaderInfo,
    pub gex_tab: Option<Tab<SampleGexWebSummary>>,
    pub vdj_t_tab: Option<Tab<SampleVdjWebSummary>>,
    pub vdj_t_gd_tab: Option<Tab<SampleVdjWebSummary>>,
    pub vdj_b_tab: Option<Tab<SampleVdjWebSummary>>,
    pub antibody_tab: Option<Tab<SampleAntibodyWebSummary>>,
    pub antigen_tab: Option<Tab<SampleAntigenWebSummary>>,
    pub crispr_tab: Option<Tab<SampleCrisprWebSummary>>,
    pub custom_feature_tab: Option<Tab<SampleCustomFeatureWebSummary>>,
}

impl ToCsvRows for SampleWebSummary {
    fn to_csv_rows(self) -> Vec<Vec<String>> {
        let mut rows = Vec::new();
        for mut v in [
            to_csv_rows_helper(self.gex_tab, "Gene Expression", "Cells"),
            to_csv_rows_helper(self.vdj_t_tab, "VDJ T", "Cells"),
            to_csv_rows_helper(self.vdj_t_gd_tab, "VDJ T GD", "Cells"),
            to_csv_rows_helper(self.vdj_b_tab, "VDJ B", "Cells"),
            to_csv_rows_helper(self.antibody_tab, "Antibody Capture", "Cells"),
            to_csv_rows_helper(self.antigen_tab, "Antigen Capture", "Cells"),
            to_csv_rows_helper(self.crispr_tab, "CRISPR Guide Capture", "Cells"),
            to_csv_rows_helper(self.custom_feature_tab, "Custom Feature", "Cells"),
        ] {
            rows.append(&mut v);
        }
        rows
    }
}

impl<T: Alert> Alert for Option<T> {
    fn alerts(&self, ctx: &AlertContext) -> Vec<AlertSpec> {
        match self {
            Some(t) => t.alerts(ctx),
            None => vec![],
        }
    }
}

#[derive(Serialize, Clone, ToCsvRows, Alert)]
pub struct LibraryGexWebSummary {
    pub parameters_table: CountParametersTable,
    pub cell_metrics_table: MetricCard<LibraryCellMetricsTable>,
    pub sequencing_metrics_table: MetricCard<SequencingMetricsTable>,
    pub mapping_metrics_table: MetricCard<GexOrRtlLibraryMappingMetricsTable>,
    pub physical_library_metrics_table: MetricCard<GexPhysicalLibraryMetricsTable>,
    pub rtl_probe_barcode_metrics_table: Option<MetricCard<RtlProbeBarcodeMetricsTable>>,
    pub gdna_table: Option<MetricCard<GdnaMetricsTable>>,
    pub targeted_plot: Option<ChartWithHelp>,
    pub targeted_table: Option<MetricCard<GexLibraryTargetedEnrichmentMetricsTable>>,
    pub targeted_alerts: Option<MetricCard<GexLibraryTargetedEnrichmentAlertsTable>>,
    pub sequencing_saturation_plot: ChartWithHelp,
    pub median_genes_per_cell_plot: ChartWithHelp,
    pub barcode_rank_plot: ChartWithHelp,
}

#[derive(Serialize, Clone, ToCsvRows, Alert)]
pub struct LibraryVdjWebSummary {
    pub parameters_table: VdjParametersTable,
    pub cell_metrics_table: MetricCard<VdjLibraryCellMetricsTable>,
    pub enrichment_metrics_table: MetricCard<VdjEnrichmentMetricsTable>,
    pub sequencing_metrics_table: MetricCard<SequencingMetricsTable>,
    pub physical_library_metrics_table: MetricCard<VdjPhysicalLibraryMetricsTable>,
    pub barcode_rank_plot: Option<ChartWithHelp>, // None if there are 0 cells
}

#[derive(Serialize, Clone, ToCsvRows, Alert)]
pub struct LibraryAntibodyOrAntigenWebSummary {
    pub parameters_table: CountParametersTable,
    pub cell_metrics_table: MetricCard<LibraryCellMetricsTable>,
    pub sequencing_metrics_table: MetricCard<SequencingMetricsTable>,
    pub mapping_metrics_table: Option<MetricCard<AntibodyLibraryMappingMetricsTable>>,
    pub physical_library_metrics_table: MetricCard<AntibodyOrAntigenPhysicalLibraryMetricsTable>,
    pub rtl_probe_barcode_metrics_table: Option<MetricCard<RtlProbeBarcodeMetricsTable>>,
    pub barcode_rank_plot: ChartWithHelp,
    pub feature_histogram: Option<RawChartWithHelp>,
}

#[derive(Serialize, Clone, ToCsvRows, Alert)]
pub struct LibraryCrisprWebSummary {
    pub parameters_table: CountParametersTable,
    pub cell_metrics_table: MetricCard<LibraryCellMetricsTable>,
    pub sequencing_metrics_table: MetricCard<SequencingMetricsTable>,
    pub physical_library_metrics_table: MetricCard<CrisprPhysicalLibraryMetricsTable>,
    pub barcode_rank_plot: ChartWithHelp,
}

#[derive(Serialize, Clone, ToCsvRows, Alert)]
pub struct LibraryCustomFeatureWebSummary {
    pub parameters_table: CountParametersTable,
    pub cell_metrics_table: MetricCard<LibraryCellMetricsTable>,
    pub sequencing_metrics_table: MetricCard<SequencingMetricsTable>,
    pub physical_library_metrics_table: MetricCard<CustomFeaturePhysicalLibraryMetricsTable>,
    pub barcode_rank_plot: ChartWithHelp,
}

#[derive(Serialize, Clone, ToCsvRows, Alert)]
pub struct LibraryCmoWebSummary {
    pub parameters_table: CountParametersTable,
    pub multiplexing_metrics_table: MetricCard<MultiplexingLibraryCellMetricsTable>,
    pub sample_assignments_table: MetricCard<MultiplexingSampleAssignmentsTable>,
    pub sequencing_metrics_table: MetricCard<SequencingMetricsTable>,
    pub physical_library_metrics_table: MetricCard<MultiplexingPhysicalLibraryMetricsTable>,
    pub cmo_metrics_table: MetricCard<MultiplexingCmoMetricsTable>,
    pub barcode_rank_plot: ChartWithHelp,
    pub jibes_biplot: Option<RawChartWithHelp>,
    pub jibes_histogram: Option<RawChartWithHelp>,
    pub cmo_umi_tsne_plot: Option<RawChartWithHelp>,
    pub cmo_tags_tsne_plot: Option<RawChartWithHelp>,
    #[serde(rename = "_resources")]
    pub resources: MultiSharedResource,
}

#[derive(Serialize, Clone, ToCsvRows, Alert)]
pub struct SampleGexWebSummary {
    pub parameters_table: CountParametersTable,
    pub hero_metrics: MetricCard<GexSampleHeroMetricsTable>,
    /// Only for multiplexed runs
    pub cell_metrics_table: Option<MetricCard<SampleCellMetricsTable>>,
    pub mapping_metrics_table: Option<MetricCard<GexOrRtlSampleMappingMetricsTable>>,
    pub gdna_table: Option<MetricCard<GdnaMetricsTable>>,
    pub barcode_rank_plot: Option<ChartWithHelp>,
    pub median_genes_per_cell_plot: Option<ChartWithHelp>,
    pub clustering_and_diffexp_plots: Value,
}

#[derive(Serialize, Deserialize, Clone)]
pub struct ClonotypeInfo {
    pub table: GenericTable,
    pub plot: PlotlyChart,
    pub help: TitleWithHelp,
}

impl Alert for ClonotypeInfo {}
impl ToCsvRows for ClonotypeInfo {}

#[derive(Serialize, Clone, ToCsvRows, Alert)]
pub struct SampleVdjWebSummary {
    pub parameters_table: VdjParametersTable,
    pub hero_metrics: MetricCard<VdjSampleHeroMetricsTable>,
    pub annotation_metrics_table: MetricCard<VdjSampleAnnotationMetricsTable>,
    pub clonotype_info: ClonotypeInfo,
    pub barcode_rank_plot: Option<ChartWithHelp>,
}

#[derive(Serialize, Clone, ToCsvRows, Alert)]
pub struct SampleAntibodyWebSummary {
    pub parameters_table: CountParametersTable,
    pub hero_metrics: MetricCard<AntibodySampleHeroMetricsTable>,
    pub antibody_treemap: Option<RawChartWithHelp>,
    /// Only for multiplexed runs
    pub cell_metrics_table: Option<MetricCard<SampleCellMetricsTable>>,
    pub mapping_metrics_table: MetricCard<AntibodySampleMappingMetricsTable>,
    pub barcode_rank_plot: Option<ChartWithHelp>,
    pub clustering_and_diffexp_plots: Option<Value>,
    pub tsne_plot: Option<RawChartWithHelp>,
    pub feature_histogram: Option<RawChartWithHelp>,
}

#[derive(Serialize, Clone, ToCsvRows, Alert)]
pub struct SampleAntigenWebSummary {
    pub parameters_table: CountParametersTable,
    pub hero_metrics: MetricCard<AntigenSampleHeroMetricsTable>,
    pub antigen_treemap: Option<RawChartWithHelp>,
    // Heatmap of clonotypes x antigen specificity hierarchically clustered
    pub clonotype_clustermap: Option<RawChartWithHelp>,
}

#[derive(Serialize, Clone, ToCsvRows, Alert)]
pub struct SampleCrisprWebSummary {
    pub parameters_table: CountParametersTable,
    pub hero_metrics: MetricCard<CrisprSampleHeroMetricsTable>,
    /// Only for multiplexed runs
    pub cell_metrics_table: Option<MetricCard<SampleCellMetricsTable>>,
    pub barcode_rank_plot: Option<ChartWithHelp>,
    pub tsne_plot: Option<RawChartWithHelp>,
}

#[derive(Serialize, Clone, ToCsvRows, Alert)]
pub struct SampleCustomFeatureWebSummary {
    pub parameters_table: CountParametersTable,
    pub hero_metrics: MetricCard<CustomFeatureSampleHeroMetricsTable>,
    /// Only for multiplexed runs
    pub cell_metrics_table: Option<MetricCard<SampleCellMetricsTable>>,
    pub barcode_rank_plot: Option<ChartWithHelp>,
    pub tsne_plot: Option<RawChartWithHelp>,
}

#[derive(Default, Debug, Clone, PartialEq, Eq, Serialize)]
pub struct MismatchedProbeBarcodePairings {
    /// Probe barcode pairings declared in the config but not detected.
    pub configured_not_detected: Vec<String>,
    /// Probe barcode pairings detected during detect chemistry but not declared
    /// in the config.
    pub detected_not_configured: Vec<String>,
}

impl MismatchedProbeBarcodePairings {
    fn formatted_value(&self) -> String {
        let count = self
            .configured_not_detected
            .len()
            .max(self.detected_not_configured.len());
        format!("{count} probe barcode(s)")
    }
}

impl ToString for MismatchedProbeBarcodePairings {
    fn to_string(&self) -> String {
        let format_pairings = |pairings: &[String]| {
            if pairings.is_empty() {
                return "none".to_string();
            }
            pairings.join(", ")
        };
        format!("Mismatch found between probe barcode pairing specified in config CSV file and chemistry detection. \
            Pairings specified in the config CSV but not detected: {}. \
            Pairings detected but not specified in the config CSV: {}. \
            This may be due to an error in the config CSV, or \
            a workflow error resulting in cross-contamination between samples or \
            some hybridization reactions having the wrong probe barcodes.",
            format_pairings(&self.configured_not_detected),
            format_pairings(&self.detected_not_configured),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::multi::plots::{format_histogram, format_jibes_biplots};
    use crate::multi::tables::GexSampleHeroMetricsTable;
    use cr_websummary::*;
    use insta;
    use metric::PercentMetric;
    use std::fs;

    const SAMPLE_ID: &str = "manual_5k_pbmc_NGSC3_ch1";
    const SAMPLE_DESC: &str = "Peripheral blood mononuclear cells (PBMCs) from a healthy donor (the same cells were used to generate connect_5k_pbmc_NGSC3_ch1)";

    pub(crate) fn sequencing_saturation_plot(x_data: Vec<f64>, y_data: Vec<f64>) -> ChartWithHelp {
        ChartWithHelp {
            plot: PlotlyChart::new_line_plot(
                x_data,
                y_data,
                "Mean Reads per Cell".to_string(),
                "Sequencing Saturation".to_string(),
            ),
            help: TitleWithHelp {
                help: "TODO".to_string(),
                // this help text is "TODO" but it doesn't show up in websummary
                title: "Sequencing saturation".to_string(),
            },
        }
    }
    pub(crate) fn median_genes_per_cell_plot(x_data: Vec<f64>, y_data: Vec<f64>) -> ChartWithHelp {
        ChartWithHelp {
            plot: PlotlyChart::new_line_plot(
                x_data,
                y_data,
                "Mean reads per cell".to_string(),
                "Median genes per cell".to_string(),
            ),
            help: TitleWithHelp {
                help: "TODO".to_string(),
                // this help text is "TODO" but it doesn't show up in websummary
                title: "Median genes per cell".to_string(),
            },
        }
    }

    pub(crate) fn barcode_rank_plot(x_data: Vec<f64>, y_data: Vec<f64>) -> ChartWithHelp {
        ChartWithHelp {
            plot: PlotlyChart::new_line_plot(
                x_data,
                y_data,
                "Barcodes".to_string(),
                "UMI counts".to_string(),
            ),
            help: TitleWithHelp {
                help: "TODO".to_string(),
                // this help text is "TODO" but it doesn't show up in websummary
                title: "Barcode rank plot".to_string(),
            },
        }
    }

    fn make_percent(percent: f64) -> Percent {
        Percent::Float(percent / 100.0)
    }

    fn make_count_percent(num: usize, percent: f64) -> CountAndPercent {
        let denom = (num as f64 / (percent / 100.0)) as usize;
        CountAndPercent(PercentMetric {
            numerator: (num as i64).into(),
            denominator: (denom as i64).into(),
        })
    }

    fn gen_count_param_table(library_type: LegacyLibraryType) -> CountParametersTable {
        CountParametersTable {
            chemistry: "Single Cell 5' R2-only".into(),
            introns_included: false,
            reference_path: "refdata-gex-GRCh38-2020-A".into(),
            transcriptome: "GRCh38-2020-A".into(),
            feature_ref_path: Some("some/feature/ref".into()),
            cmo_set_path: None,
            target_set_name: None,
            targeting_method: None,
            filter_probes: None,
            num_genes_on_target: None,
            library_type,
            throughput: None,
            tenx_cmos: Some(false),
            aligner: AlignerParam::Star,
            antigen_negative_control: false,
            dropped_tags: Default::default(),
            probe_barcodes_high_gem_overlap: Default::default(),
            mismatched_probe_barcode_pairings: None,
            unspecified_probe_barcodes_detected: Default::default(),
            specified_probe_barcodes_missing: Default::default(),
        }
    }

    fn gen_library_gex_tab(count_param_table: &CountParametersTable) -> LibraryGexWebSummary {
        LibraryGexWebSummary {
            parameters_table: count_param_table.clone(),
            cell_metrics_table: LibraryCellMetricsTable(vec![LibraryCellMetricsRow {
                physical_library_id: Some("GEX_1".to_string()),
                cell_associated_partitions: Some(973),
                mean_reads_per_cell_associated_partition: Some(FloatAsInt(47_716.0)),
                singlets_assigned_sample: Some(make_count_percent(404, 41.6)),
                partitions_with_no_cmos: Some(make_count_percent(250, 29.9)),
                partitions_called_multiplets: Some(make_count_percent(230, 24.5)),
                fraction_cells_passing_high_occupancy_filtering: None,
            }])
            .into(),
            sequencing_metrics_table: SequencingMetricsTable(vec![
                SequencingMetricsRow {
                    fastq_id: Some("seq_GEX_1_1".to_string()),
                    number_of_reads: Some(29_374_662),
                    unprocessed_reads: Some(0),
                    q30_barcode: Some(PercentF1(make_percent(90.65))),
                    q30_gem_barcode: None,
                    q30_probe_barcode: None,
                    q30_umi: Some(PercentF1(make_percent(90.65))),
                    q30_read1: Some(PercentF1(make_percent(91.7))),
                    q30_read2: Some(PercentF1(make_percent(94.65))),
                },
                SequencingMetricsRow {
                    fastq_id: Some("seq_GEX_1_2".to_string()),
                    number_of_reads: Some(29_374_662),
                    unprocessed_reads: Some(0),
                    q30_barcode: Some(PercentF1(make_percent(90.65))),
                    q30_gem_barcode: None,
                    q30_probe_barcode: None,
                    q30_umi: Some(PercentF1(make_percent(90.65))),
                    q30_read1: Some(PercentF1(make_percent(91.7))),
                    q30_read2: Some(PercentF1(make_percent(94.65))),
                },
            ])
            .into(),
            mapping_metrics_table: GexOrRtl::Gex(GexLibraryMappingMetricsTable(vec![
                GexLibraryMappingMetricsRow {
                    physical_library_id: Some("GEX_1".to_string()),
                    reads_in_library: Some(1_000_000),
                    mapped_to_genome: Some(make_percent(95.0)),
                    confidently_mapped_to_genome: Some(make_percent(91.65)),
                    confidently_mapped_to_transcriptome: Some(make_percent(61.85)),
                    confidently_mapped_to_targeted_transcriptome: Some(make_percent(61.85)),
                    confidently_mapped_to_intronic_regions: Some(make_percent(13.51)),
                    confidently_mapped_to_exonic_regions: Some(make_percent(65.59)),
                    confidently_mapped_to_intergenic_regions: Some(make_percent(3.47)),
                    confidently_mapped_antisense: Some(make_percent(0.72)),
                },
            ]))
            .into(),
            physical_library_metrics_table: GexPhysicalLibraryMetricsTable(vec![
                GexPhysicalLibraryMetricsRow {
                    physical_library_id: Some("GEX_1".to_string()),
                    number_of_reads: Some(29_374_662),
                    valid_barcodes: Some(make_percent(99.65)),
                    valid_gem_barcodes: None,
                    valid_probe_barcodes: None,
                    valid_umis: Some(make_percent(92.80)),
                    targeted_sequencing_saturation: None,
                    sequencing_saturation: Some(make_percent(15.21)),
                    reads_in_cell_associated_partitions: Some(make_percent(91.70)),
                    mean_reads_per_cell_associated_partition: Some(FloatAsInt(47_716.0)),
                    mean_targeted_reads_per_cell_associated_partition: None,
                },
            ])
            .into(),
            rtl_probe_barcode_metrics_table: None,
            gdna_table: None,
            targeted_table: None,
            targeted_plot: None,
            targeted_alerts: None,
            sequencing_saturation_plot: sequencing_saturation_plot(
                vec![0.0, 20_000.0, 40_000.0, 60_000.0],
                vec![0.0, 0.43, 0.63, 0.75],
            ),
            median_genes_per_cell_plot: median_genes_per_cell_plot(
                vec![0.0, 20_000.0, 40_000.0, 60_000.0],
                vec![0.0, 1_500.0, 1_800.0, 2_000.0],
            ),
            barcode_rank_plot: barcode_rank_plot(
                vec![1.0, 100.0, 10_000.0, 10_000.0, 500_000.0, 1_000_000.0],
                vec![10_500.0, 10_000.0, 1000.0, 10.0, 5.0, 1.0],
            ),
        }
    }

    fn gen_sample_gex_tab(count_param_table: &CountParametersTable) -> SampleGexWebSummary {
        SampleGexWebSummary {
            parameters_table: count_param_table.clone(),
            hero_metrics: GexSampleHeroMetricsTable(vec![GexSampleHeroMetricsRow {
                genome: None,
                total_singlets: Some(1_023),
                median_reads_per_singlet: Some(FloatAsInt(63_575.0)),
                median_reads_per_cell_on_target: None,
                median_genes_per_singlet: Some(FloatAsInt(2_149.0)),
                total_genes_detected: Some(16_313),
                median_umi_per_singlet: Some(FloatAsInt(19_209.0)),
                median_genes_per_cell_on_target: None,
                num_genes_detected_on_target: None,
                median_umis_per_cell_on_target: None,
                confidently_mapped_reads_in_cells: Some(make_percent(99.0)),
            }])
            .into(),
            cell_metrics_table: Some(
                GexOrRtl::Gex(GexSampleCellMetricsTable(vec![GexSampleCellMetricsRow {
                    physical_library_id: Some("GEX_1".to_string()),
                    singlets_assigned_to_this_sample: Some(make_count_percent(404, 41.6)),
                    singlets_assigned_to_other_samples: Some(make_count_percent(450, 48.6)),
                    cell_associated_partitions_not_assigned_any_samples: Some(make_count_percent(
                        20, 2.0,
                    )),
                    cell_associated_partitions_identified_as_multiplets: Some(make_count_percent(
                        200, 20.5,
                    )),
                }]))
                .into(),
            ),
            mapping_metrics_table: Some(
                GexOrRtl::Gex(GexSampleMappingMetricsTable(vec![
                    GexSampleMappingMetricsRow {
                        reads_from_cells_assigned_to_sample: Some(1_000_000),
                        mapped_to_genome: Some(make_percent(95.0)),
                        confidently_mapped_to_genome: Some(make_percent(91.65)),
                        confidently_mapped_to_transcriptome: Some(make_percent(61.85)),
                        confidently_mapped_to_targeted_transcriptome: Some(make_percent(61.85)),
                        confidently_mapped_to_intronic_regions: Some(make_percent(13.51)),
                        confidently_mapped_to_exonic_regions: Some(make_percent(65.59)),
                        confidently_mapped_to_intergenic_regions: Some(make_percent(3.47)),
                        confidently_mapped_antisense: Some(make_percent(0.72)),
                    },
                ]))
                .into(),
            ),
            gdna_table: None,
            median_genes_per_cell_plot: Some(median_genes_per_cell_plot(
                vec![0.0, 20_000.0, 40_000.0, 60_000.0],
                vec![0.0, 1_500.0, 1_800.0, 2_000.0],
            )),
            clustering_and_diffexp_plots: Value::String("CLUSTERING_PLOTS_GO_HERE".to_string()),
            barcode_rank_plot: None,
        }
    }

    fn gen_sample_ab_tab(count_param_table: &CountParametersTable) -> SampleAntibodyWebSummary {
        SampleAntibodyWebSummary {
            parameters_table: count_param_table.clone(),
            hero_metrics: AntibodySampleHeroMetricsTable(vec![AntibodySampleHeroMetricsRow {
                total_singlets: Some(1013),
                median_umis_per_singlet: Some(FloatAsInt(68_379.0)),
                antibody_reads_usable_per_cell: Some(FloatAsInt(4_608.0)),
                reads_in_cell_associated_partitions: Some(make_percent(96.7)),
            }])
            .into(),
            mapping_metrics_table: AntibodySampleMappingMetricsTable(vec![
                AntibodySampleMappingMetricsRow {
                    reads_from_cells_assigned_to_sample: Some(12345),
                    fraction_antibody_reads: Some(make_percent(95.15)),
                    fraction_reads_in_aggregate_barcodes: Some(make_percent(5.0)),
                },
            ])
            .into(),
            antibody_treemap: None,
            cell_metrics_table: Some(
                GexOrRtl::Gex(GexSampleCellMetricsTable(vec![GexSampleCellMetricsRow {
                    physical_library_id: Some("CC_1".to_string()),
                    singlets_assigned_to_this_sample: Some(make_count_percent(404, 41.6)),
                    singlets_assigned_to_other_samples: Some(make_count_percent(450, 48.6)),
                    cell_associated_partitions_not_assigned_any_samples: Some(make_count_percent(
                        20, 2.0,
                    )),
                    cell_associated_partitions_identified_as_multiplets: Some(make_count_percent(
                        200, 20.5,
                    )),
                }]))
                .into(),
            ),
            clustering_and_diffexp_plots: None,
            tsne_plot: None,
            barcode_rank_plot: None,
            feature_histogram: None,
        }
    }

    #[cfg(test)]
    fn generate_test_websummary() -> MultiWebSummary {
        let mut count_param_table = gen_count_param_table(LegacyLibraryType::GeneExpression);
        let library_gex_tab = gen_library_gex_tab(&count_param_table);

        let vdj_param_table = VdjParametersTable {
            chemistry: "Single Cell V(D)J R2-only".into(),
            vdj_reference: "vdj_GRCh38_alts_ensembl-5.0.0".into(),
            vdj_reference_path: test_refdata::refdata_path("vdj/vdj_GRCh38_alts_ensembl")
                .to_str()
                .map(String::from)
                .unwrap(),
            gamma_delta: false,
        };

        let library_vdj_tab = LibraryVdjWebSummary {
            parameters_table: vdj_param_table.clone(),
            cell_metrics_table: VdjLibraryCellMetricsTable(vec![VdjLibraryCellMetricsRow {
                physical_library_id: Some("VDJ_1".to_string()),
                vdj_filtered_bcs: Some(200),
                vdj_total_raw_read_pairs_per_filtered_bc: Some(FloatAsInt(1200.0)),
            }])
            .into(),
            enrichment_metrics_table: VdjChainTypeSpecific::VdjT(VdjTEnrichmentMetricsTable(vec![
                VdjTEnrichmentMetricsRow {
                    physical_library_id: Some("VDJ_1".to_string()),
                    multi_vdj_recombinome_mapped_reads_frac: Some(make_percent(0.6)),
                    TRA_vdj_recombinome_mapped_reads_frac: Some(make_percent(0.2)),
                    TRB_vdj_recombinome_mapped_reads_frac: Some(make_percent(0.4)),
                },
            ]))
            .into(),
            sequencing_metrics_table: SequencingMetricsTable(vec![
                SequencingMetricsRow {
                    fastq_id: Some("seq_VDJ_1_1".to_string()),
                    number_of_reads: Some(838_586),
                    unprocessed_reads: Some(0),
                    q30_barcode: Some(PercentF1(make_percent(90.65))),
                    q30_gem_barcode: None,
                    q30_probe_barcode: None,
                    q30_umi: Some(PercentF1(make_percent(90.65))),
                    q30_read1: Some(PercentF1(make_percent(91.70))),
                    q30_read2: Some(PercentF1(make_percent(94.65))),
                },
                SequencingMetricsRow {
                    fastq_id: Some("seq_VDJ_1_2".to_string()),
                    number_of_reads: Some(720_004),
                    unprocessed_reads: Some(0),
                    q30_barcode: Some(PercentF1(make_percent(90.65))),
                    q30_gem_barcode: None,
                    q30_probe_barcode: None,
                    q30_umi: Some(PercentF1(make_percent(90.65))),
                    q30_read1: Some(PercentF1(make_percent(91.70))),
                    q30_read2: Some(PercentF1(make_percent(94.65))),
                },
            ])
            .into(),
            physical_library_metrics_table: VdjPhysicalLibraryMetricsTable(vec![
                VdjPhysicalLibraryMetricsRow {
                    physical_library_id: Some("VDJ_1".to_string()),
                    VDJ_total_read_pairs: Some(291330),
                    vdj_good_bc_frac: Some(make_percent(93.20)),
                    vdj_total_raw_read_pairs_per_filtered_bc: Some(FloatAsInt(5410.0)),
                    vdj_assemblable_read_pairs_per_filtered_bc: Some(FloatAsInt(2896.0)),
                    vdj_filtered_bcs_cum_frac: Some(make_percent(84.82)),
                },
            ])
            .into(),
            barcode_rank_plot: Some(barcode_rank_plot(
                vec![1.0, 100.0, 10_000.0, 10_000.0, 500_000.0, 1_000_000.0],
                vec![10_500.0, 10_000.0, 1000.0, 10.0, 5.0, 1.0],
            )),
        };

        let library_antibody_tab = LibraryAntibodyOrAntigenWebSummary {
            parameters_table: gen_count_param_table(LegacyLibraryType::AntibodyCapture),
            cell_metrics_table: LibraryCellMetricsTable(vec![LibraryCellMetricsRow {
                physical_library_id: Some("AB_1".to_string()),
                cell_associated_partitions: Some(973),
                mean_reads_per_cell_associated_partition: Some(FloatAsInt(5471.0)),
                singlets_assigned_sample: Some(make_count_percent(404, 41.6)),
                partitions_with_no_cmos: Some(make_count_percent(250, 29.9)),
                partitions_called_multiplets: Some(make_count_percent(230, 24.5)),
                fraction_cells_passing_high_occupancy_filtering: None,
            }])
            .into(),
            sequencing_metrics_table: SequencingMetricsTable(vec![
                SequencingMetricsRow {
                    fastq_id: Some("seq_AB_1_1".to_string()),
                    number_of_reads: Some(29_374_662),
                    unprocessed_reads: Some(0),
                    q30_barcode: Some(PercentF1(make_percent(90.65))),
                    q30_gem_barcode: None,
                    q30_probe_barcode: None,
                    q30_umi: Some(PercentF1(make_percent(90.65))),
                    q30_read1: Some(PercentF1(make_percent(91.7))),
                    q30_read2: Some(PercentF1(make_percent(94.65))),
                },
                SequencingMetricsRow {
                    fastq_id: Some("seq_AB_1_2".to_string()),
                    number_of_reads: Some(29_374_662),
                    unprocessed_reads: Some(0),
                    q30_barcode: Some(PercentF1(make_percent(90.65))),
                    q30_gem_barcode: None,
                    q30_probe_barcode: None,
                    q30_umi: Some(PercentF1(make_percent(90.65))),
                    q30_read1: Some(PercentF1(make_percent(91.7))),
                    q30_read2: Some(PercentF1(make_percent(94.65))),
                },
            ])
            .into(),
            physical_library_metrics_table: AntibodyOrAntigen::Antibody(
                AntibodyPhysicalLibraryMetricsTable(vec![AntibodyPhysicalLibraryMetricsRow {
                    physical_library_id: Some("AB_1".to_string()),
                    number_of_reads: Some(8_385_586),
                    valid_barcodes: Some(make_percent(99.65)),
                    valid_gem_barcodes: Some(make_percent(99.7)),
                    valid_probe_barcodes: Some(make_percent(99.8)),
                    valid_umis: Some(make_percent(92.80)),
                    sequencing_saturation: Some(make_percent(15.21)),
                    reads_in_cell_associated_partitions: Some(make_percent(90.65)),
                    mean_reads_per_cell_associated_partition: Some(FloatAsInt(5471.0)),
                }]),
            )
            .into(),
            mapping_metrics_table: Some(
                AntibodyLibraryMappingMetricsTable(vec![AntibodyLibraryMappingMetricsRow {
                    physical_library_id: Some("AB_1".to_string()),
                    reads_in_library: Some(8_385_586),
                    fraction_antibody_reads: Some(make_percent(90.65)),
                    fraction_antibody_reads_usable: Some(make_percent(66.42)),
                    fraction_reads_in_aggregate_barcodes: Some(make_percent(3.09)),
                }])
                .into(),
            ),
            rtl_probe_barcode_metrics_table: None,
            barcode_rank_plot: barcode_rank_plot(
                vec![1.0, 100.0, 10_000.0, 10_000.0, 500_000.0, 1_000_000.0],
                vec![10_500.0, 10_000.0, 1000.0, 10.0, 5.0, 1.0],
            ),
            feature_histogram: Some(format_histogram(
                &serde_json::Value::String("HISTOGRAM_GOES_HERE".to_string()),
                "",
            )),
        };

        let library_custom_feature_tab = LibraryCustomFeatureWebSummary {
            parameters_table: gen_count_param_table(LegacyLibraryType::Custom),
            cell_metrics_table: LibraryCellMetricsTable(vec![LibraryCellMetricsRow {
                physical_library_id: Some("CC_1".to_string()),
                cell_associated_partitions: Some(973),
                mean_reads_per_cell_associated_partition: Some(FloatAsInt(5471.0)),
                singlets_assigned_sample: Some(make_count_percent(404, 41.6)),
                partitions_with_no_cmos: Some(make_count_percent(250, 29.9)),
                partitions_called_multiplets: Some(make_count_percent(230, 24.5)),
                fraction_cells_passing_high_occupancy_filtering: None,
            }])
            .into(),
            sequencing_metrics_table: SequencingMetricsTable(vec![
                SequencingMetricsRow {
                    fastq_id: Some("seq_CC_1_1".to_string()),
                    number_of_reads: Some(29_374_662),
                    unprocessed_reads: Some(0),
                    q30_barcode: Some(PercentF1(make_percent(90.65))),
                    q30_gem_barcode: None,
                    q30_probe_barcode: None,
                    q30_umi: Some(PercentF1(make_percent(90.65))),
                    q30_read1: Some(PercentF1(make_percent(91.7))),
                    q30_read2: Some(PercentF1(make_percent(94.65))),
                },
                SequencingMetricsRow {
                    fastq_id: Some("seq_CC_1_2".to_string()),
                    number_of_reads: Some(29_374_662),
                    unprocessed_reads: Some(0),
                    q30_barcode: Some(PercentF1(make_percent(90.65))),
                    q30_gem_barcode: None,
                    q30_probe_barcode: None,
                    q30_umi: Some(PercentF1(make_percent(90.65))),
                    q30_read1: Some(PercentF1(make_percent(91.7))),
                    q30_read2: Some(PercentF1(make_percent(94.65))),
                },
            ])
            .into(),
            physical_library_metrics_table: CustomFeaturePhysicalLibraryMetricsTable(vec![
                CustomFeaturePhysicalLibraryMetricsRow {
                    physical_library_id: Some("CC_1".to_string()),
                    number_of_reads: Some(8_385_586),
                    valid_barcodes: Some(make_percent(99.65)),
                    valid_umis: Some(make_percent(92.80)),
                    sequencing_saturation: Some(make_percent(15.21)),
                    reads_in_cell_associated_partitions: Some(make_percent(90.65)),
                    mean_reads_per_cell_associated_partition: Some(FloatAsInt(5471.0)),
                    fraction_feature_reads: Some(make_percent(90.65)),
                    fraction_feature_reads_usable: Some(make_percent(66.42)),
                    fraction_unknown_feature: Some(make_percent(67.01)),
                },
            ])
            .into(),
            barcode_rank_plot: barcode_rank_plot(
                vec![1.0, 100.0, 10_000.0, 10_000.0, 500_000.0, 1_000_000.0],
                vec![10_500.0, 10_000.0, 1000.0, 10.0, 5.0, 1.0],
            ),
        };

        let library_crispr_tab = LibraryCrisprWebSummary {
            parameters_table: gen_count_param_table(LegacyLibraryType::CrisprGuideCapture),
            cell_metrics_table: LibraryCellMetricsTable(vec![LibraryCellMetricsRow {
                physical_library_id: Some("GC_1".to_string()),
                cell_associated_partitions: Some(973),
                mean_reads_per_cell_associated_partition: Some(FloatAsInt(5471.0)),
                singlets_assigned_sample: Some(make_count_percent(404, 41.6)),
                partitions_with_no_cmos: Some(make_count_percent(250, 29.9)),
                partitions_called_multiplets: Some(make_count_percent(230, 24.5)),
                fraction_cells_passing_high_occupancy_filtering: None,
            }])
            .into(),
            sequencing_metrics_table: SequencingMetricsTable(vec![
                SequencingMetricsRow {
                    fastq_id: Some("seq_GC_1_1".to_string()),
                    number_of_reads: Some(29_374_662),
                    unprocessed_reads: Some(0),
                    q30_barcode: Some(PercentF1(make_percent(90.65))),
                    q30_gem_barcode: None,
                    q30_probe_barcode: None,
                    q30_umi: Some(PercentF1(make_percent(90.65))),
                    q30_read1: Some(PercentF1(make_percent(91.7))),
                    q30_read2: Some(PercentF1(make_percent(94.65))),
                },
                SequencingMetricsRow {
                    fastq_id: Some("seq_GC_1_2".to_string()),
                    number_of_reads: Some(29_374_662),
                    unprocessed_reads: Some(0),
                    q30_barcode: Some(PercentF1(make_percent(90.65))),
                    q30_gem_barcode: None,
                    q30_probe_barcode: None,
                    q30_umi: Some(PercentF1(make_percent(90.65))),
                    q30_read1: Some(PercentF1(make_percent(91.7))),
                    q30_read2: Some(PercentF1(make_percent(94.65))),
                },
            ])
            .into(),
            physical_library_metrics_table: CrisprPhysicalLibraryMetricsTable(vec![
                CrisprPhysicalLibraryMetricsRow {
                    physical_library_id: Some("GC_1".to_string()),
                    number_of_reads: Some(8_385_586),
                    valid_barcodes: Some(make_percent(99.65)),
                    valid_umis: Some(make_percent(92.80)),
                    sequencing_saturation: Some(make_percent(15.21)),
                    reads_in_cell_associated_partitions: Some(make_percent(90.65)),
                    mean_reads_per_cell_associated_partition: Some(FloatAsInt(5471.0)),
                    fraction_reads_with_putative_protospacer: Some(make_percent(90.65)),
                    fraction_guide_reads: Some(make_percent(90.65)),
                    fraction_guide_reads_usable: Some(make_percent(66.42)),
                    fraction_protospacer_not_recognized: Some(make_percent(67.01)),
                },
            ])
            .into(),
            barcode_rank_plot: barcode_rank_plot(
                vec![1.0, 100.0, 10_000.0, 10_000.0, 500_000.0, 1_000_000.0],
                vec![10_500.0, 10_000.0, 1000.0, 10.0, 5.0, 1.0],
            ),
        };

        let library_antigen_tab = LibraryAntibodyOrAntigenWebSummary {
            parameters_table: count_param_table.clone(),
            cell_metrics_table: LibraryCellMetricsTable(vec![]).into(),
            sequencing_metrics_table: SequencingMetricsTable(vec![SequencingMetricsRow {
                fastq_id: Some("seq_AG_1_1".to_string()),
                number_of_reads: Some(29_374_662),
                unprocessed_reads: Some(0),
                q30_barcode: Some(PercentF1(make_percent(90.65))),
                q30_gem_barcode: None,
                q30_probe_barcode: None,
                q30_umi: Some(PercentF1(make_percent(90.65))),
                q30_read1: Some(PercentF1(make_percent(91.7))),
                q30_read2: Some(PercentF1(make_percent(94.65))),
            }])
            .into(),
            physical_library_metrics_table: AntibodyOrAntigen::Antigen(
                AntigenPhysicalLibraryMetricsTable(vec![AntigenPhysicalLibraryMetricsRow {
                    physical_library_id: Some("AB_1".to_string()),
                    number_of_reads: Some(8_385_586),
                    valid_barcodes: Some(make_percent(99.65)),
                    valid_umis: Some(make_percent(92.80)),
                    sequencing_saturation: Some(make_percent(15.21)),
                    reads_in_cell_associated_partitions: Some(make_percent(90.65)),
                    mean_reads_per_cell_associated_partition: Some(FloatAsInt(5471.0)),
                    fraction_antigen_reads: Some(make_percent(90.65)),
                    fraction_antigen_reads_usable: Some(make_percent(66.42)),
                    fraction_unknown_antigen: Some(make_percent(67.01)),
                    fraction_reads_in_aggregate_barcodes: Some(make_percent(3.09)),
                }]),
            )
            .into(),
            mapping_metrics_table: None,
            rtl_probe_barcode_metrics_table: None,
            barcode_rank_plot: barcode_rank_plot(
                vec![1.0, 100.0, 10_000.0, 10_000.0, 500_000.0, 1_000_000.0],
                vec![10_500.0, 10_000.0, 1000.0, 10.0, 5.0, 1.0],
            ),
            feature_histogram: Some(format_histogram(
                &serde_json::Value::String("HISTOGRAM_GOES_HERE".to_string()),
                "",
            )),
        };

        let library_cmo_tab = LibraryCmoWebSummary {
            cmo_umi_tsne_plot: None,
            cmo_tags_tsne_plot: None,
            parameters_table: gen_count_param_table(LegacyLibraryType::Multiplexing),
            multiplexing_metrics_table: MultiplexingLibraryCellMetricsTable(vec![
                MultiplexingLibraryCellMetricsRow {
                    cell_associated_partitions: Some(973),
                    samples_assigned_at_least_one_singlet: Some(9),
                    singlets_assigned_to_sample: Some(make_count_percent(2700, 92.5)),
                    singlet_capture_ratio: Some(0.5),
                    cell_associated_partitions_identified_as_multiplet: Some(make_count_percent(
                        220, 20.1,
                    )),
                    median_cmo_umis_per_singlet: Some(FloatAsInt(458.0)),
                },
            ])
            .into(),
            sample_assignments_table: MultiplexingSampleAssignmentsTable(vec![
                MultiplexingSampleAssignmentsRow {
                    physical_library_id: Some("MC_1".to_string()),
                    cell_associated_partitions: Some(1_012),
                    mean_reads_per_cell: Some(FloatAsInt(966.0)),
                    samples_assigned_at_least_one_singlet: Some(1),
                    singlets_assigned_to_a_sample: Some(make_count_percent(900, 90.8)),
                    cell_associated_partitions_identified_as_multiplets: Some(make_count_percent(
                        60, 19.2,
                    )),
                    cell_associated_partitions_not_assigned_any_cmos: Some(make_count_percent(
                        6, 1.9,
                    )),
                    median_cmo_umis_per_cell_associated_partition: Some(FloatAsInt(479.0)),
                },
            ])
            .into(),
            sequencing_metrics_table: SequencingMetricsTable(vec![
                SequencingMetricsRow {
                    fastq_id: Some("seq_MC_1_1".to_string()),
                    number_of_reads: Some(29_374_662),
                    unprocessed_reads: Some(0),
                    q30_barcode: Some(PercentF1(make_percent(90.65))),
                    q30_gem_barcode: None,
                    q30_probe_barcode: None,
                    q30_umi: Some(PercentF1(make_percent(90.65))),
                    q30_read1: Some(PercentF1(make_percent(91.7))),
                    q30_read2: Some(PercentF1(make_percent(94.65))),
                },
                SequencingMetricsRow {
                    fastq_id: Some("seq_MC_1_2".to_string()),
                    number_of_reads: Some(29_374_662),
                    unprocessed_reads: Some(0),
                    q30_barcode: Some(PercentF1(make_percent(90.65))),
                    q30_gem_barcode: None,
                    q30_probe_barcode: None,
                    q30_umi: Some(PercentF1(make_percent(90.65))),
                    q30_read1: Some(PercentF1(make_percent(91.7))),
                    q30_read2: Some(PercentF1(make_percent(94.65))),
                },
            ])
            .into(),
            physical_library_metrics_table: MultiplexingPhysicalLibraryMetricsTable(vec![
                MultiplexingPhysicalLibraryMetricsRow {
                    physical_library_id: Some("MC_1".to_string()),
                    number_of_reads: Some(29_374_662),
                    valid_barcodes: Some(make_percent(99.65)),
                    valid_umis: Some(make_percent(92.80)),
                    sequencing_saturation: Some(make_percent(15.21)),
                    reads_in_cell_associated_partitions: Some(make_percent(90.65)),
                    mean_reads_per_cell_associated_partition: Some(FloatAsInt(5471.0)),
                    fraction_cmo_reads: Some(make_percent(90.65)),
                    fraction_cmo_reads_usable: Some(make_percent(66.42)),
                    fraction_unknown_cmo: Some(make_percent(94.65)),
                    fraction_reads_from_multiplets: Some(make_percent(90.77)),
                },
            ])
            .into(),
            cmo_metrics_table: MultiplexingCmoMetricsTable(vec![
                MultiplexingCmoMetricsRow {
                    gem_well_cmo: Some("CMO_2_1".to_string()),
                    reads_in_cell_associated_partitions: Some(make_percent(83.39)),
                    singlets_assigned_to_cmo: Some(make_percent(10.27)),
                    cmo_signal_to_background_ratio: Some(3.0),
                },
                MultiplexingCmoMetricsRow {
                    gem_well_cmo: Some("CMO_2_2".to_string()),
                    reads_in_cell_associated_partitions: Some(make_percent(82.2)),
                    singlets_assigned_to_cmo: Some(make_percent(11.63)),
                    cmo_signal_to_background_ratio: Some(3.0),
                },
                MultiplexingCmoMetricsRow {
                    gem_well_cmo: Some("CMO_2_3".to_string()),
                    reads_in_cell_associated_partitions: Some(make_percent(84.41)),
                    singlets_assigned_to_cmo: Some(make_percent(10.46)),
                    cmo_signal_to_background_ratio: Some(2.0),
                },
            ])
            .into(),
            barcode_rank_plot: barcode_rank_plot(
                vec![1.0, 100.0, 10_000.0, 10_000.0, 500_000.0, 1_000_000.0],
                vec![10_500.0, 10_000.0, 1000.0, 10.0, 5.0, 1.0],
            ),
            jibes_biplot: Some(format_jibes_biplots(&serde_json::Value::String(
                "BIPLOTS_GO_HERE".to_string(),
            ))),
            jibes_histogram: Some(format_histogram(
                &serde_json::Value::String("HISTOGRAM_GOES_HERE".to_string()),
                "",
            )),
            resources: TxHashMap::default(),
        };

        let ctx = AlertContext::default();
        let library_websummary = LibraryWebSummary {
            header_info: LibraryHeaderInfo {
                run_id: "Run id".to_string(),
                run_desc: "Run Description".to_string(),
                pipeline_version: "6.0.0".to_string(),
            },
            gex_tab: Some(Tab::new(library_gex_tab, &ctx)),
            vdj_t_tab: Some(Tab::new(library_vdj_tab, &ctx)),
            vdj_t_gd_tab: None,
            vdj_b_tab: None,
            antibody_tab: Some(Tab::new(library_antibody_tab, &ctx)),
            antigen_tab: Some(Tab::new(library_antigen_tab, &ctx)),
            crispr_tab: Some(Tab::new(library_crispr_tab, &ctx)),
            custom_feature_tab: Some(Tab::new(library_custom_feature_tab, &ctx)),
            cmo_tab: Some(Tab::new(library_cmo_tab, &ctx)),
            resources: TxHashMap::default(),
        };

        let sample_gex_tab = gen_sample_gex_tab(&count_param_table);
        count_param_table.library_type = LegacyLibraryType::AntibodyCapture;
        let sample_antibody_tab = gen_sample_ab_tab(&count_param_table);

        let sample_antigen_tab = SampleAntigenWebSummary {
            parameters_table: count_param_table,
            hero_metrics: AntigenSampleHeroMetricsTable(vec![
                AntigenSampleHeroMetricsRow {
                    total_singlets: Some(971),
                    median_umis_per_singlet: Some(FloatAsInt(100.0)),
                    antigen_reads_usable_per_cell: Some(FloatAsInt(200.0)),
                    feature_type: Some("Gene Expression".into()),
                },
                AntigenSampleHeroMetricsRow {
                    total_singlets: Some(921),
                    median_umis_per_singlet: Some(FloatAsInt(110.0)),
                    antigen_reads_usable_per_cell: Some(FloatAsInt(180.0)),
                    feature_type: Some("VDJ-T".into()),
                },
            ])
            .into(),
            antigen_treemap: None,
            clonotype_clustermap: None,
        };

        let sample_vdj_tab = SampleVdjWebSummary {
            parameters_table: vdj_param_table,
            hero_metrics: VdjChainTypeSpecific::VdjT(VdjTSampleHeroMetricsTable(vec![
                VdjTSampleHeroMetricsRow {
                    vdj_filtered_bcs: Some(438),
                    multi_vdj_assembly_contig_pair_productive_full_len_bc_count: Some(380),
                    TRA_vdj_assembly_umis_per_cell_median: Some(FloatAsInt(2.0)),
                    TRB_vdj_assembly_umis_per_cell_median: Some(FloatAsInt(12.0)),
                },
            ]))
            .into(),
            annotation_metrics_table: VdjChainTypeSpecific::VdjT(VdjTSampleAnnotationMetricsTable(
                vec![VdjTSampleAnnotationMetricsRow {
                    multi_vdj_assembly_contig_pair_productive_full_len_bc_frac: Some(make_percent(
                        84.93,
                    )),
                    TRA_TRB_vdj_assembly_contig_pair_productive_full_len_bc_frac: Some(
                        make_percent(69.98),
                    ),
                    TRA_vdj_assembly_prod_cdr_bc_frac: Some(make_percent(69.98)),
                    TRB_vdj_assembly_prod_cdr_bc_frac: Some(make_percent(69.98)),
                    multi_raw_vdj_paired_clonotype_diversity: Some(143.0),
                }],
            ))
            .into(),
            clonotype_info: ClonotypeInfo {
                table: GenericTable {
                    header: None,
                    rows: Vec::new(),
                },
                plot: PlotlyChart::default(),
                help: TitleWithHelp {
                    help: "Clonotype Table/Plot help".into(),
                    title: "Top 10 Clonotypes".into(),
                },
            },
            barcode_rank_plot: None,
        };

        let sample_crispr_tab = SampleCrisprWebSummary {
            parameters_table: gen_count_param_table(LegacyLibraryType::CrisprGuideCapture),
            hero_metrics: CrisprSampleHeroMetricsTable(vec![CrisprSampleHeroMetricsRow {
                total_singlets: Some(1013),
                median_umis_per_singlet: Some(FloatAsInt(68_379.0)),
                guide_reads_usable_per_cell: Some(FloatAsInt(4_608.0)),
                cells_with_one_or_more_protospacers_detected: Some(make_percent(91.0)), //make_count_percent(910, 91.0)),
                cells_with_two_or_more_protospacers_detected: Some(make_percent(0.5)), //make_count_percent(6, 0.5)),
            }])
            .into(),
            cell_metrics_table: Some(
                GexOrRtl::Gex(GexSampleCellMetricsTable(vec![GexSampleCellMetricsRow {
                    physical_library_id: Some("GC_1".to_string()),
                    singlets_assigned_to_this_sample: Some(make_count_percent(404, 41.6)),
                    singlets_assigned_to_other_samples: Some(make_count_percent(450, 48.6)),
                    cell_associated_partitions_not_assigned_any_samples: Some(make_count_percent(
                        20, 2.0,
                    )),
                    cell_associated_partitions_identified_as_multiplets: Some(make_count_percent(
                        200, 20.5,
                    )),
                }]))
                .into(),
            ),
            tsne_plot: None,
            barcode_rank_plot: None,
        };

        let sample_custom_tab = SampleCustomFeatureWebSummary {
            parameters_table: gen_count_param_table(LegacyLibraryType::Custom),
            hero_metrics: CustomFeatureSampleHeroMetricsTable(vec![
                CustomFeatureSampleHeroMetricsRow {
                    total_singlets: Some(1013),
                    median_umis_per_singlet: Some(FloatAsInt(68_379.0)),
                    feature_reads_usable_per_cell: Some(FloatAsInt(4_608.0)),
                },
            ])
            .into(),
            cell_metrics_table: Some(
                GexOrRtl::Gex(GexSampleCellMetricsTable(vec![GexSampleCellMetricsRow {
                    physical_library_id: Some("CC_1".to_string()),
                    singlets_assigned_to_this_sample: Some(make_count_percent(404, 41.6)),
                    singlets_assigned_to_other_samples: Some(make_count_percent(450, 48.6)),
                    cell_associated_partitions_not_assigned_any_samples: Some(make_count_percent(
                        20, 2.0,
                    )),
                    cell_associated_partitions_identified_as_multiplets: Some(make_count_percent(
                        200, 20.5,
                    )),
                }]))
                .into(),
            ),
            tsne_plot: None,
            barcode_rank_plot: None,
        };

        let sample_websummary = SampleWebSummary {
            header_info: SampleHeaderInfo {
                sample_id: SAMPLE_ID.to_string(),
                sample_desc: SAMPLE_DESC.to_string(),
                pipeline_version: "6.0.0".to_string(),
            },
            gex_tab: Some(Tab::new(sample_gex_tab, &ctx)),
            vdj_t_tab: Some(Tab::new(sample_vdj_tab, &ctx)),
            vdj_t_gd_tab: None,
            vdj_b_tab: None,
            antibody_tab: Some(Tab::new(sample_antibody_tab, &ctx)),
            antigen_tab: Some(Tab::new(sample_antigen_tab, &ctx)),
            crispr_tab: Some(Tab::new(sample_crispr_tab, &ctx)),
            custom_feature_tab: Some(Tab::new(sample_custom_tab, &ctx)),
        };

        MultiWebSummary {
            sample: WsSample::multi(SAMPLE_ID.to_string(), SAMPLE_DESC.to_string()),
            data: MultiWebSummaryData {
                library_websummary,
                sample_websummary,
                experimental_design: ExperimentalDesign {
                    svg: SvgGraph::new(
                        "svg string".into(),
                        "sample_1".into(),
                        Some(cr_types::CellMultiplexingType::CMO),
                    ),
                    csv: "csv goes here".to_string(),
                },
            },
            diagnostics: MultiDiagnostics::default(),
            resources: TxHashMap::default(),
        }
    }

    #[test]
    fn test_multi_websummary() {
        let multi_websummary = generate_test_websummary();
        insta::assert_json_snapshot!(&multi_websummary);
    }

    #[test]
    fn test_multi_metrics_csv() {
        let csv_path = tempfile::Builder::new()
            .prefix("test_multi_metrics")
            .suffix(".csv")
            .tempfile()
            .unwrap()
            .into_temp_path();
        generate_test_websummary()
            .to_csv(&csv_path)
            .expect("Error generating test metrics CSV.");
        let csv_str: String =
            fs::read_to_string(&csv_path).expect("Error reading test metrics CSV.");
        insta::assert_snapshot!(&csv_str);
        csv_path.close().expect("failed to delete temp csv");
    }

    #[test]
    fn multi_websummary_include_introns_info() {
        let mut count_param_table = gen_count_param_table(LegacyLibraryType::GeneExpression);
        count_param_table.introns_included = true;

        let library_gex_tab = gen_library_gex_tab(&count_param_table);

        let ctx = AlertContext::default();

        let library_websummary = LibraryWebSummary {
            header_info: LibraryHeaderInfo {
                run_id: "Run id".to_string(),
                run_desc: "Run Description".to_string(),
                pipeline_version: "7.0.0".to_string(),
            },
            gex_tab: Some(Tab::new(library_gex_tab, &ctx)),
            ..Default::default()
        };

        let sample_gex_tab = gen_sample_gex_tab(&count_param_table);

        let sample_websummary = SampleWebSummary {
            header_info: SampleHeaderInfo {
                sample_id: SAMPLE_ID.to_string(),
                sample_desc: SAMPLE_DESC.to_string(),
                pipeline_version: "7.0.0".to_string(),
            },
            gex_tab: Some(Tab::new(sample_gex_tab, &ctx)),
            ..Default::default()
        };

        let multi_websummary_introns_info = MultiWebSummary {
            sample: WsSample::multi(SAMPLE_ID.to_string(), SAMPLE_DESC.to_string()),
            data: MultiWebSummaryData {
                library_websummary,
                sample_websummary,
                experimental_design: ExperimentalDesign {
                    svg: SvgGraph::new(
                        "svg string".into(),
                        "sample_1".into(),
                        Some(cr_types::CellMultiplexingType::CMO),
                    ),
                    csv: "csv goes here".to_string(),
                },
            },
            diagnostics: Default::default(),
            resources: Default::default(),
        };

        let json =
            serde_json::to_value(&multi_websummary_introns_info).expect("Error generating JSON");
        for ws in ["library_websummary", "sample_websummary"] {
            assert!(json["data"][ws]["gex_tab"]["alerts"]
                .as_array()
                .expect("No alerts generated")
                .iter()
                .any(|v| {
                    v["level"].to_string().contains("INFO")
                        && v["title"].to_string().contains("Intron mode used")
                }));
        }
    }

    #[test]
    fn multi_websummary_ab_aggregate_rtl_alert_trigger() {
        let mut count_param_table = gen_count_param_table(LegacyLibraryType::GeneExpression);
        count_param_table.targeting_method = Some(TargetingMethod::TemplatedLigation);
        count_param_table.target_set_name = Some("Turtle".into());

        let library_gex_tab = gen_library_gex_tab(&count_param_table);
        count_param_table.library_type = LegacyLibraryType::AntibodyCapture;
        let library_ab_tab = LibraryAntibodyOrAntigenWebSummary {
            parameters_table: count_param_table.clone(),
            cell_metrics_table: LibraryCellMetricsTable(vec![]).into(),
            sequencing_metrics_table: SequencingMetricsTable(vec![]).into(),
            physical_library_metrics_table: AntibodyOrAntigen::Antibody(
                AntibodyPhysicalLibraryMetricsTable(vec![Default::default()]),
            )
            .into(),
            mapping_metrics_table: Some(
                AntibodyLibraryMappingMetricsTable(vec![AntibodyLibraryMappingMetricsRow {
                    fraction_reads_in_aggregate_barcodes: Some(make_percent(20.09)),
                    ..Default::default()
                }])
                .into(),
            ),
            rtl_probe_barcode_metrics_table: None,
            barcode_rank_plot: barcode_rank_plot(
                vec![1.0, 100.0, 10_000.0, 10_000.0, 500_000.0, 1_000_000.0],
                vec![10_500.0, 10_000.0, 1000.0, 10.0, 5.0, 1.0],
            ),
            feature_histogram: None,
        };

        let ctx = AlertContext {
            is_rtl: true,
            ..Default::default()
        };

        let library_websummary = LibraryWebSummary {
            header_info: LibraryHeaderInfo {
                run_id: "Run id".to_string(),
                run_desc: "Run Description".to_string(),
                pipeline_version: "7.0.0".to_string(),
            },
            gex_tab: Some(Tab::new(library_gex_tab, &ctx)),
            antibody_tab: Some(Tab::new(library_ab_tab, &ctx)),
            ..Default::default()
        };

        let sample_gex_tab = gen_sample_gex_tab(&count_param_table);
        count_param_table.library_type = LegacyLibraryType::AntibodyCapture;
        let sample_ab_tab = gen_sample_ab_tab(&count_param_table);

        let sample_websummary = SampleWebSummary {
            header_info: SampleHeaderInfo {
                sample_id: SAMPLE_ID.to_string(),
                sample_desc: SAMPLE_DESC.to_string(),
                pipeline_version: "7.0.0".to_string(),
            },
            gex_tab: Some(Tab::new(sample_gex_tab, &ctx)),
            antibody_tab: Some(Tab::new(sample_ab_tab, &ctx)),
            ..Default::default()
        };

        let multi_websummary = MultiWebSummary {
            sample: WsSample::multi(SAMPLE_ID.to_string(), SAMPLE_DESC.to_string()),
            data: MultiWebSummaryData {
                library_websummary,
                sample_websummary,
                experimental_design: ExperimentalDesign {
                    svg: SvgGraph::new(
                        "svg string".into(),
                        "sample_1".into(),
                        Some(cr_types::CellMultiplexingType::CMO),
                    ),
                    csv: "csv goes here".to_string(),
                },
            },
            diagnostics: Default::default(),
            resources: Default::default(),
        };

        let json = serde_json::to_value(multi_websummary).expect("Error generating JSON");
        assert!(json["data"]["library_websummary"]["antibody_tab"]["alerts"]
            .as_array()
            .expect("No alerts generated")
            .iter()
            .any(|v| {
                v["level"].to_string().contains("WARN")
                    && v["title"]
                        .to_string()
                        .contains("High Fraction of Antibody Reads in Aggregate Barcodes")
            }));
    }

    #[test]
    fn multi_websummary_ab_aggregate_rtl_alert_notrigger() {
        let mut count_param_table = gen_count_param_table(LegacyLibraryType::GeneExpression);
        count_param_table.targeting_method = Some(TargetingMethod::TemplatedLigation);
        count_param_table.target_set_name = Some("Turtle".into());

        let library_gex_tab = gen_library_gex_tab(&count_param_table);
        count_param_table.library_type = LegacyLibraryType::AntibodyCapture;
        let library_ab_tab = LibraryAntibodyOrAntigenWebSummary {
            parameters_table: count_param_table.clone(),
            cell_metrics_table: LibraryCellMetricsTable(vec![]).into(),
            sequencing_metrics_table: SequencingMetricsTable(vec![]).into(),
            physical_library_metrics_table: AntibodyOrAntigen::Antibody(
                AntibodyPhysicalLibraryMetricsTable(vec![Default::default()]),
            )
            .into(),
            mapping_metrics_table: Some(
                AntibodyLibraryMappingMetricsTable(vec![AntibodyLibraryMappingMetricsRow {
                    fraction_reads_in_aggregate_barcodes: Some(make_percent(7.09)),
                    ..Default::default()
                }])
                .into(),
            ),
            rtl_probe_barcode_metrics_table: None,
            barcode_rank_plot: barcode_rank_plot(
                vec![1.0, 100.0, 10_000.0, 10_000.0, 500_000.0, 1_000_000.0],
                vec![10_500.0, 10_000.0, 1000.0, 10.0, 5.0, 1.0],
            ),
            feature_histogram: None,
        };

        let ctx = AlertContext {
            is_rtl: true,
            ..Default::default()
        };

        let library_websummary = LibraryWebSummary {
            header_info: LibraryHeaderInfo {
                run_id: "Run id".to_string(),
                run_desc: "Run Description".to_string(),
                pipeline_version: "7.0.0".to_string(),
            },
            gex_tab: Some(Tab::new(library_gex_tab, &ctx)),
            antibody_tab: Some(Tab::new(library_ab_tab, &ctx)),
            ..Default::default()
        };

        let sample_gex_tab = gen_sample_gex_tab(&count_param_table);
        count_param_table.library_type = LegacyLibraryType::AntibodyCapture;
        let sample_ab_tab = gen_sample_ab_tab(&count_param_table);

        let sample_websummary = SampleWebSummary {
            header_info: SampleHeaderInfo {
                sample_id: SAMPLE_ID.to_string(),
                sample_desc: SAMPLE_DESC.to_string(),
                pipeline_version: "7.0.0".to_string(),
            },
            gex_tab: Some(Tab::new(sample_gex_tab, &ctx)),
            antibody_tab: Some(Tab::new(sample_ab_tab, &ctx)),
            ..Default::default()
        };

        let multi_websummary = MultiWebSummary {
            sample: WsSample::multi(SAMPLE_ID.to_string(), SAMPLE_DESC.to_string()),
            data: MultiWebSummaryData {
                library_websummary,
                sample_websummary,
                experimental_design: ExperimentalDesign {
                    svg: SvgGraph::new(
                        "svg string".into(),
                        "sample_1".into(),
                        Some(cr_types::CellMultiplexingType::CMO),
                    ),
                    csv: "csv goes here".to_string(),
                },
            },
            diagnostics: Default::default(),
            resources: Default::default(),
        };

        let json = serde_json::to_value(multi_websummary).expect("Error generating JSON");
        assert!(
            !json["data"]["library_websummary"]["antibody_tab"]["alerts"]
                .as_array()
                .expect("No alerts generated")
                .iter()
                .any(|v| {
                    v["level"].to_string().contains("WARN")
                        && v["title"]
                            .to_string()
                            .contains("High Fraction of Antibody Reads in Aggregate Barcodes")
                })
        );
    }
}
