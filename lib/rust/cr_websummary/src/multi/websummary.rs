// There are a **lot** of metrics which need to be imported here.
#![allow(clippy::wildcard_imports)]
use super::tables::*;
use crate::alert::AlertLevel;
use crate::multi::svg::SvgGraph;
use crate::{
    Alert, AlertContext, AlertSpec, CardWithMetric, ChartWithHelp, GenericTable, MakePretty,
    MetricCard, Percent, PlotlyChart, RawChartWithHelp, Tab, TableRow, TitleWithHelp, WsSample,
};
use anyhow::Result;
use cr_types::websummary::{AlertConfig, AlertIfMetricIs, MetricConfig};
use cr_types::{
    AlignerParam, BarcodeMultiplexingType, FeatureBarcodeType, LibraryType, TargetingMethod,
};
use cr_websummary::CountAndPercent;
use csv::Writer;
use itertools::Itertools;
use metric::TxHashMap;
use serde::{Deserialize, Serialize};
use serde_json::value::Value;
use std::borrow::Cow;
use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};
use std::fmt::Display;
use std::path::Path;
use tenx_websummary::components::{
    DifferentialExpressionTable, TableMetric, VegaLitePlot, WithTitle,
};
use websummary_derive::{Alert, ToCsvRows, ToJsonSummary};

const CELL_ANNOTATION_FAILURE_TITLE: &str = "No cell type annotations produced!";

const CELL_ANNOTATION_FAILURE_MESSAGE: &str = r#"<p>Please check your cellranger logs.
If you wish to attempt cell type annotation again please use
<a href="https://www.10xgenomics.com/support/software/cloud-analysis/latest/tutorials/CA-cell-annotation-pipeline">cellranger annotate</a>.
</p>"#;
const CELL_ANNOTATION_DE_WARN_TITLE: &str = "Cell type differential expression not run";
const CELL_ANNOTATION_DE_WARN_MESSAGE: &str = "Too few cell types to run differential expression.";
pub const CELL_ANNOTATION_ADVERTISEMENT_STRING: &str = r#"<p><b>Automated cell type annotation is now available for Cell Ranger!</b><br>
For details on how to run cell type annotation and species we have it available for, visit our
<a href='https://www.10xgenomics.com/support/software/cell-ranger/latest/getting-started/cr-what-is-cell-ranger'
target='_blank' title='Cell Type Annotation' rel='noopener noreferrer'>support page</a>.</p>"#;

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
    fn to_csv_rows(&self) -> Vec<Vec<String>>
    where
        Self: Sized,
    {
        vec![]
    }
}

impl ToCsvRows for String {}

#[derive(Clone, Serialize)]
pub struct JsonMetricSummary {
    pub key: String,
    pub value: Value,
    pub category: String,
    pub library_type: String,
    pub config: MetricConfig,
    /// For metrics of the same key that may have multiple values,
    /// this optional grouping key can be used to distinguish the different
    /// instances from each other, and associate them with related data.
    /// For example, in a table of probe barcode metrics, this would be
    /// the probe barcode associated with a particular metric.
    ///
    /// This is a stopgap solution - improving the namespacing of metrics in the
    /// first place would be a better long-term solution.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub grouping_key: Option<Value>,
    /// The collection of alerts triggered for this metric.
    pub alerts: Vec<AlertSpec>,
}

impl JsonMetricSummary {
    /// Construct a JSON metric summary from input data.
    ///
    /// This includes computing any alerts that were triggered.
    pub fn new<T: Serialize + MakePretty>(
        key: String,
        value: Option<&T>,
        category: String,
        library_type: String,
        config: MetricConfig,
        grouping_key: Option<Value>,
        ctx: &AlertContext,
    ) -> Self {
        Self {
            alerts: Self::construct_alerts(&key, value, &config.alerts, ctx),
            value: serde_json::to_value(value).unwrap(),
            key,
            category,
            library_type,
            config,
            grouping_key,
        }
    }

    /// Construct alerts for the provided data.
    ///
    /// This oddball static method is used as a shim to avoid code duplication
    /// in the make_tables macro until it goes away.
    pub fn construct_alerts<T: MakePretty>(
        key: &str,
        val: Option<&T>,
        alerts: &[AlertConfig],
        ctx: &AlertContext,
    ) -> Vec<AlertSpec> {
        if alerts.is_empty() {
            return vec![];
        }
        let Some(val) = val else {
            return vec![];
        };
        let (as_f64, pretty) = (val.as_f64(), val.make_pretty());
        let triggered = |thresh: Option<f64>, cond: AlertIfMetricIs| {
            let Some(thresh) = thresh else {
                return false;
            };
            match cond {
                AlertIfMetricIs::GreaterThanOrEqual => as_f64 >= thresh,
                AlertIfMetricIs::LessThanOrEqual => as_f64 <= thresh,
            }
        };
        alerts
            .iter()
            .filter(|config| {
                let mut matches_conditions = true;
                if let Some(required_include_introns_setting) = config.conditions.include_introns {
                    matches_conditions = ctx.include_introns == required_include_introns_setting;
                }
                if let Some(required_is_rtl_setting) = config.conditions.is_rtl {
                    matches_conditions = ctx.is_rtl == required_is_rtl_setting;
                }
                matches_conditions
            })
            .filter_map(|config| {
                let cond = config.alert_if().unwrap();
                if triggered(config.error_threshold, cond) {
                    return Some(AlertSpec {
                        level: AlertLevel::Error,
                        title: config.error_title(key).to_string(),
                        formatted_value: pretty.clone(),
                        message: config.detail.clone(),
                    });
                }
                if triggered(config.warn_threshold, cond) {
                    return Some(AlertSpec {
                        level: AlertLevel::Warn,
                        title: config.warn_title(key).to_string(),
                        formatted_value: pretty.clone(),
                        message: config.detail.clone(),
                    });
                }
                None
            })
            .collect()
    }
}

/// Websummary structs implementing this trait know how to convert their contents
/// into a JSON summary containing all metrics data and metadata, and triggered
/// alerts.
pub trait ToJsonSummary {
    fn to_json_summary(&self, _ctx: &AlertContext) -> Vec<JsonMetricSummary> {
        vec![]
    }
}

impl ToJsonSummary for String {}

impl<T> ToCsvRows for Option<T>
where
    T: ToCsvRows,
{
    fn to_csv_rows(&self) -> Vec<Vec<String>> {
        match self {
            Some(t) => t.to_csv_rows(),
            None => vec![],
        }
    }
}

impl<T> ToJsonSummary for Option<T>
where
    T: ToJsonSummary,
{
    fn to_json_summary(&self, ctx: &AlertContext) -> Vec<JsonMetricSummary> {
        match self {
            Some(t) => t.to_json_summary(ctx),
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
    pub reference_path: Option<String>,
    pub transcriptome: String,
    pub feature_ref_path: Option<String>,
    pub cmo_set_path: Option<String>,
    pub target_set_name: Option<String>,
    pub targeting_method: Option<TargetingMethod>,
    pub filter_probes: Option<bool>,
    pub disable_ab_aggregate_detection: bool,
    pub disable_high_occupancy_gem_detection: bool,
    pub num_genes_on_target: Option<usize>,
    pub library_type: LibraryType,
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
            disable_ab_aggregate_detection,
            disable_high_occupancy_gem_detection,
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
        if library_type.is_gex() {
            rows.extend([
                TableRow::two_col(
                    "Reference Path",
                    reference_path.unwrap_or_else(|| "None".to_string()),
                ),
                TableRow::two_col("Transcriptome", transcriptome),
            ]);
            if aligner != AlignerParam::Hurtle {
                rows.push(TableRow::two_col("Include Introns", introns_included));
            }
            if targeting_method.is_some() {
                rows.push(TableRow::two_col(
                    "Probe Set Name",
                    target_set_name.unwrap(),
                ));
            }
            if let Some(n_genes) = num_genes_on_target {
                rows.push(TableRow::two_col(
                    "Number of Genes Targeted",
                    n_genes.make_pretty(),
                ));
            }
            if filter_probes == Some(false) {
                rows.push(TableRow::two_col("Filter Probes", "Disabled"));
            }
        }
        // Feature Tabs
        if let (LibraryType::FeatureBarcodes(_), Some(feature_ref)) =
            (library_type, feature_ref_path)
        {
            rows.push(TableRow::two_col("Feature Reference", feature_ref));
        }

        let fb_type = library_type.feature_barcode_type();
        if let (Some(cmo_set_path), Some(FeatureBarcodeType::Multiplexing)) =
            (cmo_set_path, fb_type)
        {
            rows.push(TableRow::two_col("CMO Set", cmo_set_path));
        }

        if fb_type == Some(FeatureBarcodeType::Antigen) {
            rows.push(TableRow::two_col(
                "Control Specified",
                antigen_negative_control,
            ));
        }
        if disable_ab_aggregate_detection {
            rows.push(TableRow::two_col(
                "Antibody Aggregate Filtering",
                "Disabled",
            ));
        }
        if disable_high_occupancy_gem_detection {
            rows.push(TableRow::two_col(
                "High-occupancy GEM Filtering",
                "Disabled",
            ));
        }

        GenericTable {
            header: None,
            rows,
            grouping_header: None,
        }
    }
}

impl Alert for CountParametersTable {
    fn alerts(&self, ctx: &AlertContext) -> Vec<AlertSpec> {
        let mut alerts = vec![];

        if self.library_type.is_fb_type(FeatureBarcodeType::Antigen)
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
        if ctx.is_fiveprime
            && ctx.multiplexing_method
                == Some(BarcodeMultiplexingType::CellLevel(cr_types::CellLevel::CMO))
        {
            alerts.push(AlertSpec {
                level: AlertLevel::Warn,
                title: "Unsupported combination of 5' chemistry with multiplexing".to_string(),
                formatted_value: String::default(),
                message: "Multiplexing performance cannot be guaranteed".to_string(),
            });
        }
        if ctx.is_hashtag_multiplexed {
            alerts.push(AlertSpec {
                level: AlertLevel::Warn,
                title: "Unsupported multiplexing tag used".to_string(),
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
            });
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
impl ToJsonSummary for CountParametersTable {}

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

impl<G, R> ToJsonSummary for GexOrRtl<G, R>
where
    G: ToJsonSummary,
    R: ToJsonSummary,
{
    fn to_json_summary(&self, ctx: &AlertContext) -> Vec<JsonMetricSummary> {
        match self {
            GexOrRtl::Gex(g) => g.to_json_summary(ctx),
            GexOrRtl::Rtl(r) => r.to_json_summary(ctx),
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

impl<B, G> ToJsonSummary for AntibodyOrAntigen<B, G>
where
    B: ToJsonSummary,
    G: ToJsonSummary,
{
    fn to_json_summary(&self, ctx: &AlertContext) -> Vec<JsonMetricSummary> {
        match self {
            AntibodyOrAntigen::Antibody(b) => b.to_json_summary(ctx),
            AntibodyOrAntigen::Antigen(g) => g.to_json_summary(ctx),
        }
    }
}

pub type AntibodyOrAntigenPhysicalLibraryMetricsTable =
    AntibodyOrAntigen<AntibodyPhysicalLibraryMetricsTable, AntigenPhysicalLibraryMetricsTable>;

#[derive(Serialize, Deserialize, Clone)]
pub enum CmoOrHashtag<C, H> {
    Cmo(C),
    Hashtag(H),
}

impl<C, H> From<CmoOrHashtag<C, H>> for CardWithMetric
where
    CardWithMetric: From<C>,
    CardWithMetric: From<H>,
{
    fn from(src: CmoOrHashtag<C, H>) -> CardWithMetric {
        match src {
            CmoOrHashtag::Cmo(c) => c.into(),
            CmoOrHashtag::Hashtag(h) => h.into(),
        }
    }
}

impl<C, H> Alert for CmoOrHashtag<C, H>
where
    C: Alert,
    H: Alert,
{
    fn alerts(&self, ctx: &AlertContext) -> Vec<AlertSpec> {
        match self {
            CmoOrHashtag::Cmo(c) => c.alerts(ctx),
            CmoOrHashtag::Hashtag(h) => h.alerts(ctx),
        }
    }
}

impl<C, H> ToJsonSummary for CmoOrHashtag<C, H>
where
    C: ToJsonSummary,
    H: ToJsonSummary,
{
    fn to_json_summary(&self, ctx: &AlertContext) -> Vec<JsonMetricSummary> {
        match self {
            CmoOrHashtag::Cmo(c) => c.to_json_summary(ctx),
            CmoOrHashtag::Hashtag(h) => h.to_json_summary(ctx),
        }
    }
}

pub type CmoOrHashtagPerTagMetricsTable =
    CmoOrHashtag<CmoPerTagMetricsTable, HashtagPerTagMetricsTable>;

pub type CmoOrHashtagMultiplexingQualityTable =
    CmoOrHashtag<CmoMultiplexingQualityTable, HashtagMultiplexingQualityTable>;

pub struct MultiplexingTagMetricsRow {
    pub gem_well_tag: Option<String>,
    pub sample_id: Option<String>,
    pub tag_reads_in_cell_associated_partitions: Option<Percent>,
    pub singlets_assigned_to_tag: Option<CountAndPercent>,
    pub tag_signal_to_background_ratio: Option<f64>,
}

#[derive(Default)]
pub struct MultiplexingTagMetricsRows(pub Vec<MultiplexingTagMetricsRow>);

impl MultiplexingTagMetricsRows {
    pub fn sort_tag_rows(&mut self) {
        self.0
            .sort_by(|a, b| match (&a.gem_well_tag, &b.gem_well_tag) {
                (Some(ref _x), None) => Ordering::Greater,
                (None, Some(ref _y)) => Ordering::Less,
                (Some(ref x), Some(ref y)) => x.cmp(y),
                (None, None) => Ordering::Equal,
            });
    }
}

impl From<MultiplexingTagMetricsRows> for crate::multi::tables::CmoPerTagMetricsTable {
    fn from(src: MultiplexingTagMetricsRows) -> Self {
        crate::multi::tables::CmoPerTagMetricsTable(
            src.0
                .into_iter()
                .map(|row| CmoPerTagMetricsRow {
                    gem_well_cmo: row.gem_well_tag,
                    sample_id: row.sample_id,
                    cmo_reads_in_cell_associated_partitions: row
                        .tag_reads_in_cell_associated_partitions,
                    singlets_assigned_to_cmo: row.singlets_assigned_to_tag,
                    cmo_signal_to_background_ratio: row.tag_signal_to_background_ratio,
                })
                .collect::<Vec<CmoPerTagMetricsRow>>(),
        )
    }
}

impl From<MultiplexingTagMetricsRows> for crate::multi::tables::HashtagPerTagMetricsTable {
    fn from(src: MultiplexingTagMetricsRows) -> Self {
        crate::multi::tables::HashtagPerTagMetricsTable(
            src.0
                .into_iter()
                .map(|row| HashtagPerTagMetricsRow {
                    gem_well_hashtag: row.gem_well_tag,
                    sample_id: row.sample_id,
                    hashtag_reads_in_cell_associated_partitions: row
                        .tag_reads_in_cell_associated_partitions,
                    singlets_assigned_to_hashtag: row.singlets_assigned_to_tag,
                    hashtag_signal_to_background_ratio: row.tag_signal_to_background_ratio,
                })
                .collect::<Vec<HashtagPerTagMetricsRow>>(),
        )
    }
}

// Websummary data structures may have shared _resources that they access by key
// If a websummary struct has any _resources they are emptied and bubbled up and stored in the top-level _resources
pub type MultiSharedResource = TxHashMap<String, Value>;

impl Alert for MultiSharedResource {}
impl ToCsvRows for MultiSharedResource {}
impl ToJsonSummary for MultiSharedResource {}

#[derive(Debug, Serialize, Deserialize, PartialEq, Eq, Clone, Default)]
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
        GenericTable {
            header: None,
            rows,
            grouping_header: None,
        }
    }
}

impl Alert for VdjParametersTable {
    fn alerts(&self, ctx: &AlertContext) -> Vec<AlertSpec> {
        let mut alerts = vec![];
        if self.gamma_delta {
            alerts.push(AlertSpec {
                level: AlertLevel::Warn,
                title: "Unsupported workflow used".to_string(),
                formatted_value: String::default(),
                message: "Gamma Delta TCR analysis is not a supported workflow. Algorithm performance cannot be guaranteed.".to_string(),
            });
        }
        if ctx.multiplexing_method
            == Some(BarcodeMultiplexingType::ReadLevel(cr_types::ReadLevel::OH))
            && !ctx.library_types.contains(&LibraryType::Gex)
            && ctx.library_types.iter().any(|t: &LibraryType| t.is_vdj())
        {
            alerts.push(AlertSpec {
            level: AlertLevel::Info,
            title: "GEX and VDJ libraries are recommended to be analyzed together for optimal results.".to_string(),
            formatted_value: String::default(),
            message: r#"Multiplexing performance cannot be guaranteed"#.into(),
            });
        }

        alerts
    }
}

impl ToCsvRows for VdjParametersTable {}
impl ToJsonSummary for VdjParametersTable {}

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

impl<T, B, GD> ToJsonSummary for VdjChainTypeSpecific<T, B, GD>
where
    T: ToJsonSummary,
    B: ToJsonSummary,
    GD: ToJsonSummary,
{
    fn to_json_summary(&self, ctx: &AlertContext) -> Vec<JsonMetricSummary> {
        match self {
            VdjChainTypeSpecific::VdjT(t) => t.to_json_summary(ctx),
            VdjChainTypeSpecific::VdjB(b) => b.to_json_summary(ctx),
            VdjChainTypeSpecific::VdjTgd(t_gd) => t_gd.to_json_summary(ctx),
        }
    }
}

pub type VdjLibraryCellMetricsTable = VdjChainTypeSpecific<
    VdjTLibraryCellMetricsTable,
    VdjBLibraryCellMetricsTable,
    VdjTgdLibraryCellMetricsTable,
>;
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
impl ToJsonSummary for Value {}

#[derive(Serialize, Clone)]
pub struct MultiWebSummary {
    // FIXME delete this CELLRANGER-8423
    pub sample: WsSample,
    pub library: MultiWebSummaryLibraryData,
    pub per_sample: Vec<MultiWebSummarySampleData>,
    pub experimental_design: ExperimentalDesign,
    pub diagnostics: MultiDiagnostics,
    pub sample_diagnostics: Vec<SampleDiagnostics>,
    // FIXME delete this CELLRANGER-8423
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
    pub unknown_feature_barcode_seqs: HashMap<String, Value>,
}

#[derive(Serialize, Clone, Default)]
#[allow(non_snake_case)]
pub struct SampleDiagnostics {
    pub vdj_t: Option<VdjDiagnostics>,
    pub vdj_b: Option<VdjDiagnostics>,
    pub vdj_t_gd: Option<VdjDiagnostics>,
}

impl MultiWebSummary {
    pub fn to_csv(&self, filename: &Path, sample_index: usize) -> Result<()> {
        let mut rows = self.library.data.to_csv_rows();
        rows.extend(self.per_sample[sample_index].data.to_csv_rows());
        rows.sort();
        rows.dedup();

        let mut writer = Writer::from_path(filename)?;
        writer.write_record(METRICS_SUMMARY_CSV_HEADER)?;
        for row in rows {
            writer.write_record(row)?;
        }

        writer.flush()?;
        Ok(())
    }
}

#[derive(Serialize, Clone)]
pub struct MultiWebSummaryLibraryData {
    pub data: LibraryWebSummary,
    /// All unique library types used in the analysis.
    pub types: Vec<LibraryType>,
    pub metrics: Vec<JsonMetricSummary>,
}

#[derive(Serialize, Clone)]
pub struct MultiWebSummarySampleData {
    pub data: SampleWebSummary,
    pub metrics: Vec<JsonMetricSummary>,
}

#[derive(Serialize, Clone)]
pub struct ExperimentalDesign {
    pub svg: SvgGraph,
    pub csv: String,
    pub multiplexing_method: Option<BarcodeMultiplexingType>,
    /// True if this experiment is using a Flex kit.
    pub is_rtl: bool,
    /// A flag indicating if this is a barnyard analysis.
    pub is_barnyard: bool,
}

fn to_csv_rows_helper<T: ToCsvRows>(
    tab: &Option<T>,
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

fn to_json_summary_helper<T: ToJsonSummary>(
    tab: Option<&T>,
    library_type: &str,
    sample_or_library: &str,
    ctx: &AlertContext,
) -> Vec<JsonMetricSummary> {
    let Some(tab) = tab else { return vec![] };
    let mut rows = tab.to_json_summary(ctx);
    for row in &mut rows {
        row.category = sample_or_library.to_string();
        row.library_type = library_type.to_string();
    }
    rows
}

mod section {
    pub const GEX: &str = "Gene Expression";
    pub const VDJ_T: &str = "VDJ T";
    pub const VDJ_T_GD: &str = "VDJ T GD";
    pub const VDJ_B: &str = "VDJ B";
    pub const AB: &str = "Antibody Capture";
    pub const AG: &str = "Antigen Capture";
    pub const CRISPR: &str = "CRISPR Guide Capture";
    pub const CUSTOM: &str = "Custom Feature";
    pub const CMO: &str = "Multiplexing Capture";
    pub const CA: &str = "Cell Annotation";
    pub const HASHTAG: &str = "Hashtag";
}

const TAB_CELLS: &str = "Cells";
const TAB_LIBRARY: &str = "Library";

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
    pub hashtag_tab: Option<Tab<LibraryHashtagWebSummary>>,
    #[serde(rename = "_resources")]
    pub resources: MultiSharedResource,
}

impl ToCsvRows for LibraryWebSummary {
    fn to_csv_rows(&self) -> Vec<Vec<String>> {
        [
            to_csv_rows_helper(&self.gex_tab, section::GEX, TAB_LIBRARY),
            to_csv_rows_helper(&self.vdj_t_tab, section::VDJ_T, TAB_LIBRARY),
            to_csv_rows_helper(&self.vdj_t_gd_tab, section::VDJ_T_GD, TAB_LIBRARY),
            to_csv_rows_helper(&self.vdj_b_tab, section::VDJ_B, TAB_LIBRARY),
            to_csv_rows_helper(&self.antibody_tab, section::AB, TAB_LIBRARY),
            to_csv_rows_helper(&self.antigen_tab, section::AG, TAB_LIBRARY),
            to_csv_rows_helper(&self.crispr_tab, section::CRISPR, TAB_LIBRARY),
            to_csv_rows_helper(&self.custom_feature_tab, section::CUSTOM, TAB_LIBRARY),
            to_csv_rows_helper(&self.cmo_tab, section::CMO, TAB_LIBRARY),
            to_csv_rows_helper(&self.hashtag_tab, section::HASHTAG, TAB_LIBRARY),
        ]
        .concat()
    }
}

impl ToJsonSummary for LibraryWebSummary {
    fn to_json_summary(&self, ctx: &AlertContext) -> Vec<JsonMetricSummary> {
        [
            to_json_summary_helper(self.gex_tab.as_ref(), section::GEX, TAB_LIBRARY, ctx),
            to_json_summary_helper(self.vdj_t_tab.as_ref(), section::VDJ_T, TAB_LIBRARY, ctx),
            to_json_summary_helper(
                self.vdj_t_gd_tab.as_ref(),
                section::VDJ_T_GD,
                TAB_LIBRARY,
                ctx,
            ),
            to_json_summary_helper(self.vdj_b_tab.as_ref(), section::VDJ_B, TAB_LIBRARY, ctx),
            to_json_summary_helper(self.antibody_tab.as_ref(), section::AB, TAB_LIBRARY, ctx),
            to_json_summary_helper(self.antigen_tab.as_ref(), section::AG, TAB_LIBRARY, ctx),
            to_json_summary_helper(self.crispr_tab.as_ref(), section::CRISPR, TAB_LIBRARY, ctx),
            to_json_summary_helper(
                self.custom_feature_tab.as_ref(),
                section::CUSTOM,
                TAB_LIBRARY,
                ctx,
            ),
            to_json_summary_helper(self.cmo_tab.as_ref(), section::CMO, TAB_LIBRARY, ctx),
            to_json_summary_helper(
                self.hashtag_tab.as_ref(),
                section::HASHTAG,
                TAB_LIBRARY,
                ctx,
            ),
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
    pub cell_annotation_tab: Option<Tab<SampleCellAnnotationWebSummary>>,
}

impl ToCsvRows for SampleWebSummary {
    fn to_csv_rows(&self) -> Vec<Vec<String>> {
        let mut rows = Vec::new();
        for mut v in [
            to_csv_rows_helper(&self.gex_tab, section::GEX, TAB_CELLS),
            to_csv_rows_helper(&self.vdj_t_tab, section::VDJ_T, TAB_CELLS),
            to_csv_rows_helper(&self.vdj_t_gd_tab, section::VDJ_T_GD, TAB_CELLS),
            to_csv_rows_helper(&self.vdj_b_tab, section::VDJ_B, TAB_CELLS),
            to_csv_rows_helper(&self.antibody_tab, section::AB, TAB_CELLS),
            to_csv_rows_helper(&self.antigen_tab, section::AG, TAB_CELLS),
            to_csv_rows_helper(&self.crispr_tab, section::CRISPR, TAB_CELLS),
            to_csv_rows_helper(&self.custom_feature_tab, section::CUSTOM, TAB_CELLS),
            to_csv_rows_helper(&self.cell_annotation_tab, section::CA, TAB_CELLS),
        ] {
            rows.append(&mut v);
        }
        rows
    }
}

impl ToJsonSummary for SampleWebSummary {
    fn to_json_summary(&self, ctx: &AlertContext) -> Vec<JsonMetricSummary> {
        [
            to_json_summary_helper(self.gex_tab.as_ref(), section::GEX, TAB_CELLS, ctx),
            to_json_summary_helper(self.vdj_t_tab.as_ref(), section::VDJ_T, TAB_CELLS, ctx),
            to_json_summary_helper(
                self.vdj_t_gd_tab.as_ref(),
                section::VDJ_T_GD,
                TAB_CELLS,
                ctx,
            ),
            to_json_summary_helper(self.vdj_b_tab.as_ref(), section::VDJ_B, TAB_CELLS, ctx),
            to_json_summary_helper(self.antibody_tab.as_ref(), section::AB, TAB_CELLS, ctx),
            to_json_summary_helper(self.antigen_tab.as_ref(), section::AG, TAB_CELLS, ctx),
            to_json_summary_helper(self.crispr_tab.as_ref(), section::CRISPR, TAB_CELLS, ctx),
            to_json_summary_helper(
                self.custom_feature_tab.as_ref(),
                section::CUSTOM,
                TAB_CELLS,
                ctx,
            ),
            to_json_summary_helper(
                self.cell_annotation_tab.as_ref(),
                section::CA,
                TAB_CELLS,
                ctx,
            ),
        ]
        .concat()
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

#[derive(Serialize, Clone, ToCsvRows, ToJsonSummary, Alert)]
pub struct LibraryGexWebSummary {
    pub parameters_table: CountParametersTable,
    #[serde(skip)]
    pub cell_metrics_table: MetricCard<LibraryCellMetricsTable>,
    #[serde(skip)]
    pub sequencing_metrics_table: MetricCard<SequencingMetricsTable>,
    #[serde(skip)]
    pub mapping_metrics_table: MetricCard<GexOrRtlLibraryMappingMetricsTable>,
    #[serde(skip)]
    pub physical_library_metrics_table: MetricCard<GexPhysicalLibraryMetricsTable>,
    #[serde(skip)]
    pub rtl_probe_barcode_metrics_table: Option<MetricCard<RtlProbeBarcodeMetricsTable>>,
    #[serde(skip)]
    pub ocm_per_overhang_metrics_table: Option<MetricCard<OcmPerOverhangMetricsTable>>,
    #[serde(skip)]
    pub gdna_table: Option<MetricCard<GdnaMetricsTable>>,
    pub sequencing_saturation_plot: ChartWithHelp,
    pub median_genes_per_cell_plot: ChartWithHelp,
    pub barcode_rank_plot: ChartWithHelp,
}

#[derive(Serialize, Clone, ToCsvRows, ToJsonSummary, Alert)]
pub struct LibraryVdjWebSummary {
    pub parameters_table: VdjParametersTable,
    #[serde(skip)]
    pub cell_metrics_table: MetricCard<VdjLibraryCellMetricsTable>,
    #[serde(skip)]
    pub enrichment_metrics_table: MetricCard<VdjEnrichmentMetricsTable>,
    #[serde(skip)]
    pub sequencing_metrics_table: MetricCard<SequencingMetricsTable>,
    #[serde(skip)]
    pub physical_library_metrics_table: MetricCard<VdjPhysicalLibraryMetricsTable>,
    pub barcode_rank_plot: Option<ChartWithHelp>, // None if there are 0 cells
    pub metrics_per_ocm_barcode_table: Option<MetricCard<VdjLibraryMetricsPerOcmBarcodeTable>>,
    pub metrics_per_hashtag_id_table: Option<MetricCard<VdjLibraryMetricsPerHashtagIdTable>>,
}

#[derive(Serialize, Clone, ToCsvRows, ToJsonSummary, Alert)]
pub struct LibraryAntibodyOrAntigenWebSummary {
    pub parameters_table: CountParametersTable,
    #[serde(skip)]
    pub cell_metrics_table: MetricCard<LibraryCellMetricsTable>,
    #[serde(skip)]
    pub sequencing_metrics_table: MetricCard<SequencingMetricsTable>,
    #[serde(skip)]
    pub mapping_metrics_table: Option<MetricCard<AntibodyLibraryMappingMetricsTable>>,
    #[serde(skip)]
    pub physical_library_metrics_table: MetricCard<AntibodyOrAntigenPhysicalLibraryMetricsTable>,
    #[serde(skip)]
    pub rtl_probe_barcode_metrics_table: Option<MetricCard<RtlProbeBarcodeMetricsTable>>,
    pub barcode_rank_plot: ChartWithHelp,
    pub feature_histogram: Option<RawChartWithHelp>,
}

#[derive(Serialize, Clone, ToCsvRows, ToJsonSummary, Alert)]
pub struct LibraryCrisprWebSummary {
    pub parameters_table: CountParametersTable,
    #[serde(skip)]
    pub cell_metrics_table: MetricCard<LibraryCellMetricsTable>,
    #[serde(skip)]
    pub sequencing_metrics_table: MetricCard<SequencingMetricsTable>,
    #[serde(skip)]
    pub mapping_metrics_table: MetricCard<CrisprLibraryMappingMetricsTable>,
    #[serde(skip)]
    pub physical_library_metrics_table: MetricCard<CrisprPhysicalLibraryMetricsTable>,
    #[serde(skip)]
    pub rtl_probe_barcode_metrics_table: Option<MetricCard<RtlProbeBarcodeMetricsTable>>,
    pub barcode_rank_plot: ChartWithHelp,
}

#[derive(Serialize, Clone, ToCsvRows, ToJsonSummary, Alert)]
pub struct LibraryCustomFeatureWebSummary {
    pub parameters_table: CountParametersTable,
    #[serde(skip)]
    pub cell_metrics_table: MetricCard<LibraryCellMetricsTable>,
    #[serde(skip)]
    pub sequencing_metrics_table: MetricCard<SequencingMetricsTable>,
    #[serde(skip)]
    pub physical_library_metrics_table: MetricCard<CustomFeaturePhysicalLibraryMetricsTable>,
    pub barcode_rank_plot: ChartWithHelp,
}

#[derive(Serialize, Clone, ToCsvRows, ToJsonSummary, Alert)]
pub struct LibraryCmoWebSummary {
    pub parameters_table: CountParametersTable,
    #[serde(skip)]
    pub multiplexing_quality_table: MetricCard<CmoOrHashtagMultiplexingQualityTable>,
    #[serde(skip)]
    pub sequencing_metrics_table: MetricCard<SequencingMetricsTable>,
    #[serde(skip)]
    pub physical_library_metrics_table: MetricCard<MultiplexingPhysicalLibraryMetricsTable>,
    #[serde(skip)]
    pub cmo_metrics_table: MetricCard<CmoOrHashtagPerTagMetricsTable>,
    pub barcode_rank_plot: ChartWithHelp,
    pub jibes_biplot: Option<RawChartWithHelp>,
    pub jibes_histogram: Option<RawChartWithHelp>,
    pub cmo_umi_projection_plot: Option<RawChartWithHelp>,
    pub cmo_tags_projection_plot: Option<RawChartWithHelp>,
    #[serde(rename = "_resources")]
    pub resources: MultiSharedResource,
}

#[derive(Serialize, Clone, ToCsvRows, ToJsonSummary, Alert)]
pub struct LibraryHashtagWebSummary {
    pub parameters_table: CountParametersTable,
    #[serde(skip)]
    pub multiplexing_quality_table: MetricCard<CmoOrHashtagMultiplexingQualityTable>,
    #[serde(skip)]
    pub hashtag_metrics_table: MetricCard<CmoOrHashtagPerTagMetricsTable>,
    pub jibes_biplot: Option<RawChartWithHelp>,
    pub jibes_histogram: Option<RawChartWithHelp>,
    pub hashtag_umi_projection_plot: Option<RawChartWithHelp>,
    pub hashtag_tags_projection_plot: Option<RawChartWithHelp>,
    #[serde(rename = "_resources")]
    pub resources: MultiSharedResource,
}

impl Alert for Option<String> {}
#[derive(Serialize, Clone, ToCsvRows, ToJsonSummary, Alert)]
pub struct SampleGexWebSummary {
    pub parameters_table: CountParametersTable,
    #[serde(skip)]
    pub hero_metrics: MetricCard<GexSampleHeroMetricsTable>,
    /// Only for multiplexed runs
    #[serde(skip)]
    pub cell_metrics_table: Option<MetricCard<SampleCellMetricsTable>>,
    #[serde(skip)]
    pub mapping_metrics_table: Option<MetricCard<GexOrRtlSampleMappingMetricsTable>>,
    #[serde(skip)]
    pub gdna_table: Option<MetricCard<GdnaMetricsTable>>,
    pub barcode_rank_plot: Option<ChartWithHelp>,
    pub median_genes_per_cell_plot: Option<ChartWithHelp>,
    pub clustering_and_diffexp_plots: Value,
    pub disclaimer: Option<String>,
}

#[derive(Serialize, Deserialize, Clone)]
pub struct ClonotypeInfo {
    pub table: GenericTable,
    pub plot: PlotlyChart,
    pub help: TitleWithHelp,
}

impl Alert for ClonotypeInfo {}
impl ToCsvRows for ClonotypeInfo {}
impl ToJsonSummary for ClonotypeInfo {}

#[derive(Serialize, Clone, ToCsvRows, ToJsonSummary, Alert)]
pub struct SampleVdjWebSummary {
    pub parameters_table: VdjParametersTable,
    #[serde(skip)]
    pub hero_metrics: MetricCard<VdjSampleHeroMetricsTable>,
    #[serde(skip)]
    pub annotation_metrics_table: MetricCard<VdjSampleAnnotationMetricsTable>,
    #[serde(skip)]
    pub enrichment_metrics_table: Option<MetricCard<VdjEnrichmentMetricsTable>>,
    pub clonotype_info: Option<ClonotypeInfo>,
    pub barcode_rank_plot: Option<ChartWithHelp>,
}

#[derive(Serialize, Clone, ToCsvRows, ToJsonSummary, Alert)]
pub struct SampleAntibodyWebSummary {
    pub parameters_table: CountParametersTable,
    #[serde(skip)]
    pub hero_metrics: MetricCard<AntibodySampleHeroMetricsTable>,
    pub antibody_treemap: Option<RawChartWithHelp>,
    /// Only for multiplexed runs
    #[serde(skip)]
    pub cell_metrics_table: Option<MetricCard<SampleCellMetricsTable>>,
    #[serde(skip)]
    pub mapping_metrics_table: MetricCard<AntibodySampleMappingMetricsTable>,
    pub barcode_rank_plot: Option<ChartWithHelp>,
    pub clustering_and_diffexp_plots: Option<Value>,
    pub projection_plot: Option<RawChartWithHelp>,
    pub feature_histogram: Option<RawChartWithHelp>,
}

#[derive(Serialize, Clone, ToCsvRows, ToJsonSummary, Alert)]
pub struct SampleAntigenWebSummary {
    pub parameters_table: CountParametersTable,
    #[serde(skip)]
    pub hero_metrics: MetricCard<AntigenSampleHeroMetricsTable>,
    pub antigen_treemap: Option<RawChartWithHelp>,
    // Heatmap of clonotypes x antigen specificity hierarchically clustered
    pub clonotype_clustermap: Option<RawChartWithHelp>,
}

#[derive(Serialize, Clone, ToCsvRows, ToJsonSummary, Alert)]
pub struct SampleCrisprWebSummary {
    pub parameters_table: CountParametersTable,
    #[serde(skip)]
    pub hero_metrics: MetricCard<CrisprSampleHeroMetricsTable>,
    /// Only for multiplexed runs
    #[serde(skip)]
    pub cell_metrics_table: Option<MetricCard<SampleCellMetricsTable>>,
    #[serde(skip)]
    pub mapping_metrics_table: MetricCard<CrisprSampleMappingMetricsTable>,
    pub barcode_rank_plot: Option<ChartWithHelp>,
    pub projection_plot: Option<RawChartWithHelp>,
}

#[derive(Serialize, Clone, ToCsvRows, ToJsonSummary, Alert)]
pub struct SampleCustomFeatureWebSummary {
    pub parameters_table: CountParametersTable,
    #[serde(skip)]
    pub hero_metrics: MetricCard<CustomFeatureSampleHeroMetricsTable>,
    /// Only for multiplexed runs
    #[serde(skip)]
    pub cell_metrics_table: Option<MetricCard<SampleCellMetricsTable>>,
    pub barcode_rank_plot: Option<ChartWithHelp>,
    pub projection_plot: Option<RawChartWithHelp>,
}

#[derive(Serialize, Clone)]
pub struct SampleCellAnnotationWebSummary {
    pub parameters_table: Option<TableMetric>,
    pub disclaimer: Option<String>,
    pub cas_success: Option<bool>,
    pub cell_annotation_disable_differential_expression: Option<bool>,
    pub cell_annotation_cell_types_chart: Option<WithTitle<VegaLitePlot>>,
    pub cell_annotation_violin_plot_chart: Option<WithTitle<VegaLitePlot>>,
    pub cell_annotation_umap_plot_chart: Option<WithTitle<VegaLitePlot>>,
    pub cell_annotation_diffexp_table: Option<WithTitle<DifferentialExpressionTable>>,
}

impl Alert for SampleCellAnnotationWebSummary {
    fn alerts(&self, _ctx: &AlertContext) -> Vec<AlertSpec> {
        match (
            self.cas_success,
            self.cell_annotation_disable_differential_expression,
        ) {
            (Some(false), _) => vec![AlertSpec {
                level: cr_websummary::alert::AlertLevel::Error,
                title: CELL_ANNOTATION_FAILURE_TITLE.to_string(),
                formatted_value: String::default(),
                message: CELL_ANNOTATION_FAILURE_MESSAGE.to_string(),
            }],
            (_, Some(true)) => vec![AlertSpec {
                level: AlertLevel::Warn,
                title: CELL_ANNOTATION_DE_WARN_TITLE.to_string(),
                formatted_value: String::default(),
                message: CELL_ANNOTATION_DE_WARN_MESSAGE.to_string(),
            }],
            _ => Vec::new(),
        }
    }
}

impl ToCsvRows for SampleCellAnnotationWebSummary {}
impl ToJsonSummary for SampleCellAnnotationWebSummary {}

#[derive(Default, Debug, Clone, PartialEq, Eq, Serialize)]
pub struct MismatchedProbeBarcodePairings {
    /// Probe barcode pairings declared in the config but not detected.
    configured_not_detected: Vec<String>,
    /// Probe barcode pairings detected during detect chemistry but not declared
    /// in the config.
    detected_not_configured: Vec<String>,
}

impl MismatchedProbeBarcodePairings {
    pub fn new(
        configured_probe_barcode_pairing: &HashSet<String>,
        detected_probe_barcode_pairing: &HashSet<String>,
    ) -> Self {
        Self {
            configured_not_detected: configured_probe_barcode_pairing
                .difference(detected_probe_barcode_pairing)
                .cloned()
                .sorted()
                .collect(),
            detected_not_configured: detected_probe_barcode_pairing
                .difference(configured_probe_barcode_pairing)
                .cloned()
                .sorted()
                .collect(),
        }
    }

    fn formatted_value(&self) -> String {
        let count = self
            .configured_not_detected
            .len()
            .max(self.detected_not_configured.len());
        format!("{count} probe barcode(s)")
    }
}

impl Display for MismatchedProbeBarcodePairings {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let format_pairings = |pairings: &[String]| {
            if pairings.is_empty() {
                return Cow::Borrowed("none");
            }
            Cow::Owned(pairings.join(", "))
        };
        write!(f,
            "Mismatch found between probe barcode pairing specified in config CSV file and chemistry detection. \
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
    use cr_types::{CellLevel, ReadLevel};
    use cr_websummary::*;
    use insta::assert_snapshot;
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

    fn gen_count_param_table(library_type: LibraryType) -> CountParametersTable {
        CountParametersTable {
            chemistry: "Single Cell 5' R2-only".to_string(),
            introns_included: false,
            reference_path: Some("refdata-gex-GRCh38-2020-A".to_string()),
            transcriptome: "GRCh38-2020-A".to_string(),
            feature_ref_path: Some("some/feature/ref".to_string()),
            cmo_set_path: None,
            target_set_name: None,
            targeting_method: None,
            filter_probes: None,
            disable_ab_aggregate_detection: false,
            disable_high_occupancy_gem_detection: false,
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
                    q30_barcode: Some(make_percent(90.65)),
                    q30_gem_barcode: None,
                    q30_probe_barcode: None,
                    q30_umi: Some(make_percent(90.65)),
                    q30_read1: Some(make_percent(91.7)),
                    q30_read2: Some(make_percent(94.65)),
                },
                SequencingMetricsRow {
                    fastq_id: Some("seq_GEX_1_2".to_string()),
                    number_of_reads: Some(29_374_662),
                    unprocessed_reads: Some(0),
                    q30_barcode: Some(make_percent(90.65)),
                    q30_gem_barcode: None,
                    q30_probe_barcode: None,
                    q30_umi: Some(make_percent(90.65)),
                    q30_read1: Some(make_percent(91.7)),
                    q30_read2: Some(make_percent(94.65)),
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
                    sequencing_saturation: Some(make_percent(15.21)),
                    reads_in_cell_associated_partitions: Some(make_percent(91.70)),
                    mean_reads_per_cell_associated_partition: Some(FloatAsInt(47_716.0)),
                },
            ])
            .into(),
            rtl_probe_barcode_metrics_table: None,
            ocm_per_overhang_metrics_table: None,
            gdna_table: None,
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
                mean_reads_per_cell: Some(FloatAsInt(63_575.0)),
                median_genes_per_singlet: Some(FloatAsInt(2_149.0)),
                total_genes_detected: Some(16_313),
                median_umi_per_singlet: Some(FloatAsInt(19_209.0)),
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
            disclaimer: Some("This is a disclaimer".to_string()),
        }
    }

    fn gen_sample_ab_tab(count_param_table: &CountParametersTable) -> SampleAntibodyWebSummary {
        SampleAntibodyWebSummary {
            parameters_table: count_param_table.clone(),
            hero_metrics: AntibodySampleHeroMetricsTable(vec![AntibodySampleHeroMetricsRow {
                total_singlets: Some(1013),
                median_umis_per_singlet: Some(FloatAsInt(68_379.0)),
                antibody_reads_usable_per_cell: Some(FloatAsInt(4_608.0)),
                reads_in_cells: Some(make_percent(96.7)),
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
            projection_plot: None,
            barcode_rank_plot: None,
            feature_histogram: None,
        }
    }

    #[cfg(test)]
    fn generate_test_websummary() -> MultiWebSummary {
        let mut count_param_table = gen_count_param_table(LibraryType::Gex);
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
            cell_metrics_table: VdjChainTypeSpecific::VdjT(VdjTLibraryCellMetricsTable(vec![
                VdjTLibraryCellMetricsRow {
                    physical_library_id: Some("VDJ_1".to_string()),
                    vdj_filtered_bcs: Some(200),
                    multi_vdj_assembly_contig_pair_productive_full_len_bc_frac: Some(make_percent(
                        81.4,
                    )),
                    TRA_TRB_vdj_assembly_contig_pair_productive_full_len_bc_frac: Some(
                        make_percent(92.3),
                    ),
                    TRA_vdj_assembly_prod_cdr_bc_frac: Some(make_percent(87.6)),
                    TRB_vdj_assembly_prod_cdr_bc_frac: Some(make_percent(91.3)),
                },
            ]))
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
                    q30_barcode: Some(make_percent(90.65)),
                    q30_gem_barcode: None,
                    q30_probe_barcode: None,
                    q30_umi: Some(make_percent(90.65)),
                    q30_read1: Some(make_percent(91.70)),
                    q30_read2: Some(make_percent(94.65)),
                },
                SequencingMetricsRow {
                    fastq_id: Some("seq_VDJ_1_2".to_string()),
                    number_of_reads: Some(720_004),
                    unprocessed_reads: Some(0),
                    q30_barcode: Some(make_percent(90.65)),
                    q30_gem_barcode: None,
                    q30_probe_barcode: None,
                    q30_umi: Some(make_percent(90.65)),
                    q30_read1: Some(make_percent(91.70)),
                    q30_read2: Some(make_percent(94.65)),
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
            metrics_per_ocm_barcode_table: None,
            metrics_per_hashtag_id_table: None,
        };

        let library_antibody_tab = LibraryAntibodyOrAntigenWebSummary {
            parameters_table: gen_count_param_table(LibraryType::FeatureBarcodes(
                FeatureBarcodeType::Antibody,
            )),
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
                    q30_barcode: Some(make_percent(90.65)),
                    q30_gem_barcode: None,
                    q30_probe_barcode: None,
                    q30_umi: Some(make_percent(90.65)),
                    q30_read1: Some(make_percent(91.7)),
                    q30_read2: Some(make_percent(94.65)),
                },
                SequencingMetricsRow {
                    fastq_id: Some("seq_AB_1_2".to_string()),
                    number_of_reads: Some(29_374_662),
                    unprocessed_reads: Some(0),
                    q30_barcode: Some(make_percent(90.65)),
                    q30_gem_barcode: None,
                    q30_probe_barcode: None,
                    q30_umi: Some(make_percent(90.65)),
                    q30_read1: Some(make_percent(91.7)),
                    q30_read2: Some(make_percent(94.65)),
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
            parameters_table: gen_count_param_table(LibraryType::FeatureBarcodes(
                FeatureBarcodeType::Custom,
            )),
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
                    q30_barcode: Some(make_percent(90.65)),
                    q30_gem_barcode: None,
                    q30_probe_barcode: None,
                    q30_umi: Some(make_percent(90.65)),
                    q30_read1: Some(make_percent(91.7)),
                    q30_read2: Some(make_percent(94.65)),
                },
                SequencingMetricsRow {
                    fastq_id: Some("seq_CC_1_2".to_string()),
                    number_of_reads: Some(29_374_662),
                    unprocessed_reads: Some(0),
                    q30_barcode: Some(make_percent(90.65)),
                    q30_gem_barcode: None,
                    q30_probe_barcode: None,
                    q30_umi: Some(make_percent(90.65)),
                    q30_read1: Some(make_percent(91.7)),
                    q30_read2: Some(make_percent(94.65)),
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
            parameters_table: gen_count_param_table(LibraryType::FeatureBarcodes(
                FeatureBarcodeType::Crispr,
            )),
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
                    q30_barcode: Some(make_percent(90.65)),
                    q30_gem_barcode: None,
                    q30_probe_barcode: None,
                    q30_umi: Some(make_percent(90.65)),
                    q30_read1: Some(make_percent(91.7)),
                    q30_read2: Some(make_percent(94.65)),
                },
                SequencingMetricsRow {
                    fastq_id: Some("seq_GC_1_2".to_string()),
                    number_of_reads: Some(29_374_662),
                    unprocessed_reads: Some(0),
                    q30_barcode: Some(make_percent(90.65)),
                    q30_gem_barcode: None,
                    q30_probe_barcode: None,
                    q30_umi: Some(make_percent(90.65)),
                    q30_read1: Some(make_percent(91.7)),
                    q30_read2: Some(make_percent(94.65)),
                },
            ])
            .into(),
            mapping_metrics_table: CrisprLibraryMappingMetricsTable(vec![
                CrisprLibraryMappingMetricsRow {
                    physical_library_id: Some("GC_1".to_string()),
                    number_of_reads: Some(8_385_586),
                    fraction_reads_with_putative_protospacer: Some(make_percent(90.65)),
                    fraction_guide_reads: Some(make_percent(90.65)),
                    fraction_guide_reads_usable: Some(make_percent(66.42)),
                    fraction_protospacer_not_recognized: Some(make_percent(67.01)),
                },
            ])
            .into(),
            physical_library_metrics_table: CrisprPhysicalLibraryMetricsTable(vec![
                CrisprPhysicalLibraryMetricsRow {
                    physical_library_id: Some("GC_1".to_string()),
                    number_of_reads: Some(8_385_586),
                    valid_barcodes: Some(make_percent(99.65)),
                    valid_gem_barcodes: Some(make_percent(99.8)),
                    valid_probe_barcodes: Some(make_percent(99.7)),
                    valid_umis: Some(make_percent(92.80)),
                    sequencing_saturation: Some(make_percent(15.21)),
                    reads_in_cell_associated_partitions: Some(make_percent(90.65)),
                    mean_reads_per_cell_associated_partition: Some(FloatAsInt(5471.0)),
                },
            ])
            .into(),
            rtl_probe_barcode_metrics_table: None,
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
                q30_barcode: Some(make_percent(90.65)),
                q30_gem_barcode: None,
                q30_probe_barcode: None,
                q30_umi: Some(make_percent(90.65)),
                q30_read1: Some(make_percent(91.7)),
                q30_read2: Some(make_percent(94.65)),
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
            cmo_umi_projection_plot: None,
            cmo_tags_projection_plot: None,
            parameters_table: gen_count_param_table(LibraryType::FeatureBarcodes(
                FeatureBarcodeType::Multiplexing,
            )),
            multiplexing_quality_table: CmoOrHashtag::Cmo(CmoMultiplexingQualityTable(vec![
                CmoMultiplexingQualityRow {
                    cell_associated_partitions: Some(973),
                    samples_assigned_at_least_one_singlet: Some(9),
                    singlets_assigned_to_a_sample: Some(make_count_percent(2700, 92.5)),
                    singlet_capture_ratio: Some(0.5),
                    median_cmo_umis_per_singlet: Some(FloatAsInt(458.0)),
                    cell_associated_partitions_identified_as_multiplets: Some(make_count_percent(
                        60, 19.2,
                    )),
                    cell_associated_partitions_not_assigned_any_tags: Some(make_count_percent(
                        6, 1.9,
                    )),
                },
            ]))
            .into(),
            sequencing_metrics_table: SequencingMetricsTable(vec![
                SequencingMetricsRow {
                    fastq_id: Some("seq_MC_1_1".to_string()),
                    number_of_reads: Some(29_374_662),
                    unprocessed_reads: Some(0),
                    q30_barcode: Some(make_percent(90.65)),
                    q30_gem_barcode: None,
                    q30_probe_barcode: None,
                    q30_umi: Some(make_percent(90.65)),
                    q30_read1: Some(make_percent(91.7)),
                    q30_read2: Some(make_percent(94.65)),
                },
                SequencingMetricsRow {
                    fastq_id: Some("seq_MC_1_2".to_string()),
                    number_of_reads: Some(29_374_662),
                    unprocessed_reads: Some(0),
                    q30_barcode: Some(make_percent(90.65)),
                    q30_gem_barcode: None,
                    q30_probe_barcode: None,
                    q30_umi: Some(make_percent(90.65)),
                    q30_read1: Some(make_percent(91.7)),
                    q30_read2: Some(make_percent(94.65)),
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
            cmo_metrics_table: CmoOrHashtag::Cmo(CmoPerTagMetricsTable(vec![
                CmoPerTagMetricsRow {
                    gem_well_cmo: Some("CMO_2_1".to_string()),
                    sample_id: Some("sample_1".to_string()),
                    cmo_reads_in_cell_associated_partitions: Some(make_percent(83.39)),
                    singlets_assigned_to_cmo: Some(make_count_percent(1027, 10.27)),
                    cmo_signal_to_background_ratio: Some(3.0),
                },
                CmoPerTagMetricsRow {
                    gem_well_cmo: Some("CMO_2_2".to_string()),
                    sample_id: Some("sample_2".to_string()),
                    cmo_reads_in_cell_associated_partitions: Some(make_percent(82.2)),
                    singlets_assigned_to_cmo: Some(make_count_percent(1163, 11.63)),
                    cmo_signal_to_background_ratio: Some(3.0),
                },
                CmoPerTagMetricsRow {
                    gem_well_cmo: Some("CMO_2_3".to_string()),
                    sample_id: Some("sample_2".to_string()),
                    cmo_reads_in_cell_associated_partitions: Some(make_percent(84.41)),
                    singlets_assigned_to_cmo: Some(make_count_percent(1046, 10.46)),
                    cmo_signal_to_background_ratio: Some(2.0),
                },
            ]))
            .into(),
            barcode_rank_plot: barcode_rank_plot(
                vec![1.0, 100.0, 10_000.0, 10_000.0, 500_000.0, 1_000_000.0],
                vec![10_500.0, 10_000.0, 1000.0, 10.0, 5.0, 1.0],
            ),
            jibes_biplot: Some(format_jibes_biplots(
                &serde_json::Value::String("BIPLOTS_GO_HERE".to_string()),
                "CMO",
            )),
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
            hashtag_tab: None,
            resources: TxHashMap::default(),
        };

        let sample_gex_tab = gen_sample_gex_tab(&count_param_table);
        count_param_table.library_type = LibraryType::Antibody;
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
                    vdj_filtered_bcs_cum_frac: None,
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
            enrichment_metrics_table: None,
            clonotype_info: Some(ClonotypeInfo {
                table: GenericTable {
                    header: None,
                    rows: Vec::new(),
                    grouping_header: None,
                },
                plot: PlotlyChart::default(),
                help: TitleWithHelp {
                    help: "Clonotype Table/Plot help".into(),
                    title: "Top 10 Clonotypes".into(),
                },
            }),
            barcode_rank_plot: None,
        };

        let sample_crispr_tab = SampleCrisprWebSummary {
            parameters_table: gen_count_param_table(LibraryType::FeatureBarcodes(
                FeatureBarcodeType::Crispr,
            )),
            hero_metrics: CrisprSampleHeroMetricsTable(vec![CrisprSampleHeroMetricsRow {
                total_singlets: Some(1013),
                median_umis_per_singlet: Some(FloatAsInt(68_379.0)),
                guide_reads_usable_per_cell: Some(FloatAsInt(4_608.0)),
                reads_in_cells: Some(make_percent(96.7)),
                cells_with_one_or_more_protospacers_detected: Some(make_percent(91.0)),
                cells_with_two_or_more_protospacers_detected: Some(make_percent(0.5)),
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
            mapping_metrics_table: CrisprSampleMappingMetricsTable(vec![
                CrisprSampleMappingMetricsRow {
                    number_of_reads: Some(8_385_586),
                    fraction_reads_with_putative_protospacer: Some(make_percent(90.65)),
                    fraction_guide_reads: Some(make_percent(90.65)),
                    fraction_protospacer_not_recognized: Some(make_percent(67.01)),
                },
            ])
            .into(),
            projection_plot: None,
            barcode_rank_plot: None,
        };

        let sample_custom_tab = SampleCustomFeatureWebSummary {
            parameters_table: gen_count_param_table(LibraryType::FeatureBarcodes(
                FeatureBarcodeType::Custom,
            )),
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
            projection_plot: None,
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
            cell_annotation_tab: None,
        };

        MultiWebSummary {
            sample: WsSample::multi(SAMPLE_ID.to_string(), SAMPLE_DESC.to_string()),
            library: MultiWebSummaryLibraryData {
                metrics: library_websummary.to_json_summary(&ctx),
                types: vec![
                    LibraryType::Gex,
                    LibraryType::Antibody,
                    LibraryType::Antigen,
                    LibraryType::Crispr,
                    LibraryType::Custom,
                    LibraryType::Vdj(cr_types::VdjChainType::VdjT),
                ],
                data: library_websummary,
            },
            per_sample: vec![MultiWebSummarySampleData {
                metrics: sample_websummary.to_json_summary(&ctx),
                data: sample_websummary,
            }],
            experimental_design: ExperimentalDesign {
                svg: SvgGraph::new(
                    "svg string".into(),
                    "sample_1".into(),
                    Some(BarcodeMultiplexingType::CellLevel(CellLevel::CMO)),
                ),
                csv: "csv goes here".to_string(),
                multiplexing_method: Some(BarcodeMultiplexingType::CellLevel(CellLevel::CMO)),
                is_rtl: false,
                is_barnyard: false,
            },
            diagnostics: MultiDiagnostics::default(),
            sample_diagnostics: Default::default(),
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
            .to_csv(&csv_path, 0)
            .expect("Error generating test metrics CSV.");
        let csv_str: String =
            fs::read_to_string(&csv_path).expect("Error reading test metrics CSV.");
        assert_snapshot!(&csv_str);
        csv_path.close().expect("failed to delete temp csv");
    }

    #[test]
    fn multi_websummary_ab_aggregate_rtl_alert_trigger() {
        let mut count_param_table = gen_count_param_table(LibraryType::Gex);
        count_param_table.targeting_method = Some(TargetingMethod::TemplatedLigation);
        count_param_table.target_set_name = Some("Turtle".into());

        let library_gex_tab = gen_library_gex_tab(&count_param_table);
        count_param_table.library_type = LibraryType::Antibody;
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
        count_param_table.library_type = LibraryType::Antibody;
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
            library: MultiWebSummaryLibraryData {
                metrics: library_websummary.to_json_summary(&ctx),
                types: vec![LibraryType::Gex, LibraryType::Antibody],
                data: library_websummary,
            },
            per_sample: vec![MultiWebSummarySampleData {
                metrics: sample_websummary.to_json_summary(&ctx),
                data: sample_websummary,
            }],
            experimental_design: ExperimentalDesign {
                svg: SvgGraph::new(
                    "svg string".into(),
                    "sample_1".into(),
                    Some(BarcodeMultiplexingType::ReadLevel(ReadLevel::RTL)),
                ),
                csv: "csv goes here".to_string(),
                multiplexing_method: Some(BarcodeMultiplexingType::ReadLevel(ReadLevel::RTL)),
                is_rtl: true,
                is_barnyard: false,
            },
            diagnostics: Default::default(),
            sample_diagnostics: Default::default(),
            resources: Default::default(),
        };

        let json = serde_json::to_value(multi_websummary).expect("Error generating JSON");
        assert!(json["library"]["data"]["antibody_tab"]["alerts"]
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
        let mut count_param_table = gen_count_param_table(LibraryType::Gex);
        count_param_table.targeting_method = Some(TargetingMethod::TemplatedLigation);
        count_param_table.target_set_name = Some("Turtle".into());

        let library_gex_tab = gen_library_gex_tab(&count_param_table);
        count_param_table.library_type = LibraryType::Antibody;
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
        count_param_table.library_type = LibraryType::Antibody;
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
            library: MultiWebSummaryLibraryData {
                metrics: library_websummary.to_json_summary(&ctx),
                types: vec![LibraryType::Gex, LibraryType::Antibody],
                data: library_websummary,
            },
            per_sample: vec![MultiWebSummarySampleData {
                metrics: sample_websummary.to_json_summary(&ctx),
                data: sample_websummary,
            }],
            experimental_design: ExperimentalDesign {
                svg: SvgGraph::new(
                    "svg string".into(),
                    "sample_1".into(),
                    Some(BarcodeMultiplexingType::ReadLevel(ReadLevel::RTL)),
                ),
                csv: "csv goes here".to_string(),
                multiplexing_method: Some(BarcodeMultiplexingType::ReadLevel(ReadLevel::RTL)),
                is_rtl: true,
                is_barnyard: false,
            },
            diagnostics: Default::default(),
            sample_diagnostics: Default::default(),
            resources: Default::default(),
        };

        let json = serde_json::to_value(multi_websummary).expect("Error generating JSON");
        assert!(!json["library"]["data"]["antibody_tab"]["alerts"]
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
}
