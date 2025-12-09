//! Process a specification file into a metrics ETL operation.
//!
//! ETL = Extract, Transform, Load
//! In this case, we extract metrics from the pipestance scope, transform their
//! data types and process possible alerts, then load them into output data
//! structures.
#![deny(missing_docs)]

use super::websummary::{JsonMetricSummary, Section};
use crate::alert::{AlertContext, AlertSpec};
use crate::value::TryFromValue;
use crate::{CountAndPercent, FloatAsInt, MakePretty, Percent};
use anyhow::{Context, Result, anyhow, bail, ensure};
use cr_types::websummary::{AlertConfig, ExtractMode, MetricEtlConfig};
use itertools::Itertools;
use metric::{JsonReport, PercentMetric, TxHashMap};
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::collections::{BTreeMap, HashMap};
use std::marker::PhantomData;
use toml;
use vdj_reference::VdjReceptor;

/// Load the metric groups.
/// Return metric processors for library and sample level.
pub fn load_metrics_etl() -> Result<(MetricsProcessor, MetricsProcessor)> {
    load_and_divide_groups(include_str!("metrics_etl.toml"))
}

/// Load the VDJ metric groups.
/// Return metric processors for library and sample level.
pub fn load_metrics_etl_vdj() -> Result<(MetricsProcessor, MetricsProcessor)> {
    load_and_divide_groups(include_str!("metrics_etl_vdj.toml"))
}

/// Load the special metric groups.
/// Return a single metric processor.
pub fn load_metrics_etl_special() -> Result<MetricsProcessor> {
    Ok(MetricsProcessor::new(load_metric_configs_from_str(
        include_str!("metrics_etl_special.toml"),
    )?))
}

fn load_metric_configs_from_str(contents: &str) -> Result<BTreeMap<String, MetricEtlGroup>> {
    let parsed: BTreeMap<String, MetricEtlGroup> = toml::from_str(contents)?;
    for (key, t) in &parsed {
        t.validate().with_context(|| key.clone())?;
    }
    Ok(parsed)
}

/// Load metric ETLs from the contents string.
/// Divide them up into two processors based on the tier of metrics being extracted.
fn load_and_divide_groups(contents: &str) -> Result<(MetricsProcessor, MetricsProcessor)> {
    let (lib_metrics, sample_metrics) = load_metric_configs_from_str(contents)?
        .into_iter()
        .partition::<BTreeMap<String, _>, _>(|(_key, table)| table.tier == MetricTier::Library);

    Ok((
        MetricsProcessor::new(lib_metrics),
        MetricsProcessor::new(sample_metrics),
    ))
}

/// Contains metric ETL definitions as well as a set of transformers.
pub struct MetricsProcessor {
    groups: BTreeMap<String, MetricEtlGroup>,
    transformers: Transformers,
}

impl MetricsProcessor {
    fn new(groups: BTreeMap<String, MetricEtlGroup>) -> Self {
        Self {
            groups,
            transformers: Transformers::default(),
        }
    }

    /// Return a copy of this processor with a fresh collection of transformers.
    pub fn with_default_transformers(&self) -> Self {
        Self::new(self.groups.clone())
    }

    /// Register a new transformer with this processor with the provided type key.
    pub fn add_transformer<T: Transformer + 'static, K: Into<String>>(&mut self, key: K, t: T) {
        self.transformers.add(key, t);
    }

    /// Process the provided metrics.
    pub fn process(
        &self,
        section: Section,
        active_conditions: &ActiveConditions,
        alert_ctx: &AlertContext,
        metrics: &TxHashMap<String, Value>,
    ) -> Result<Vec<JsonMetricSummary>> {
        let mut out = Vec::new();
        for (key, group) in &self.groups {
            out.extend(
                process_metric_group(
                    section,
                    group,
                    active_conditions,
                    alert_ctx,
                    metrics,
                    &self.transformers,
                )
                .with_context(|| key.clone())?
                .into_iter(),
            );
        }
        Ok(out)
    }

    /// Process all metrics for the specified named group.
    ///
    /// This is used for handling individual metric groups that are going to be
    /// extracted from a type that represents a collection of metrics.
    pub fn process_group<T: JsonReport>(
        &self,
        group: &str,
        section: Section,
        active_conditions: &ActiveConditions,
        alert_ctx: &AlertContext,
        metrics: &T,
    ) -> Result<Vec<JsonMetricSummary>> {
        let reporter = metrics.to_json_reporter();
        let table = self.groups.get(group).ok_or_else(|| {
            anyhow!(
                "metrics group {group} not in processor; groups: {}",
                self.groups.keys().format(", ")
            )
        })?;
        process_metric_group(
            section,
            table,
            active_conditions,
            alert_ctx,
            reporter.inner(),
            &self.transformers,
        )
        .with_context(|| group.to_string())
    }

    /// Process a specific provided metric config.
    /// This is intended to be used as a shim to process metrics whose keys
    /// or other config properties must be determined dynamically.
    ///
    /// TODO CELLRANGER-8444 this is janky, consider replacing this with a
    /// metrics pre-processor that remaps metric keys to a single expected value.
    pub fn process_one(
        &self,
        alert_ctx: &AlertContext,
        metrics: &TxHashMap<String, Value>,
        config: &MetricEtlConfig,
        tier: MetricTier,
        section: Section,
    ) -> Result<JsonMetricSummary> {
        assert!(matches!(
            config.extract,
            ExtractMode::Required | ExtractMode::Placeholder
        ));
        Ok(extract_and_transform(
            config.json_key.as_ref().unwrap(),
            tier,
            section,
            config,
            alert_ctx,
            metrics,
            &self.transformers,
        )?
        .unwrap())
    }
}

/// Whether a metric should be extracted from library-level or sample-level metrics.
#[derive(Serialize, Deserialize, Debug, PartialEq, Eq, Clone, Copy, strum::Display)]
pub enum MetricTier {
    /// Metrics should be extracted from library level.
    Library,
    /// Metrics should be extracted from sample level.
    Cells,
}

/// Define the conditions under which a metric should be extracted.
#[derive(Default, Debug, Deserialize, Clone)]
pub struct MetricConditions {
    /// Which sections of the websummary is this metric extracted for?
    /// FIXME this is basically just LibraryType with the addition of Hashtag.
    #[serde(default)]
    pub sections: Vec<Section>,
    /// Is this metric associated with a specific VDJ receptor?
    pub vdj_receptor: Option<VdjReceptor>,
    /// Should this metric be extracted if we do/don't have a VDJ reference?
    pub has_vdj_reference: Option<bool>,
    /// Covers all forms of multiplexing.
    pub is_multiplexed: Option<bool>,
    /// Covers both CMO and Hashtag multiplexing.
    pub is_cell_multiplexed: Option<bool>,
    /// Covers both RTL and OCM multiplexing.
    pub is_read_multiplexed: Option<bool>,
    /// Is/is not RTL data.
    pub is_rtl: Option<bool>,
    /// Does/doesn't have GDNA.
    pub has_gdna: Option<bool>,
    /// Include introns mode is/isn't set.
    pub include_introns: Option<bool>,
}

fn matches_option<T: PartialEq>(opt: Option<T>, val: T) -> bool {
    opt.is_none_or(|v| v == val)
}

impl MetricConditions {
    fn matches(&self, cond: &ActiveConditions) -> bool {
        let matches_receptor = cond
            .vdj_receptor
            .is_none_or(|receptor| matches_option(self.vdj_receptor, receptor));
        self.matches_section(cond.section)
            && matches_receptor
            && matches_option(self.has_vdj_reference, cond.has_vdj_reference)
            && matches_option(self.is_multiplexed, cond.is_multiplexed)
            && matches_option(self.is_cell_multiplexed, cond.is_cell_multiplexed)
            && matches_option(self.is_read_multiplexed, cond.is_read_multiplexed)
            && matches_option(self.is_rtl, cond.is_rtl)
            && matches_option(self.has_gdna, cond.has_gdna)
    }

    fn matches_section(&self, this_section: Section) -> bool {
        if self.sections.is_empty() {
            return true;
        }
        self.sections.contains(&this_section)
    }
}

/// The current conditions under which a metrics ETL is taking place.
#[derive(Default, Debug)]
pub struct ActiveConditions {
    /// The section of the websummary.
    pub section: Section,
    /// The VDJ library type.
    pub vdj_receptor: Option<VdjReceptor>,
    /// Has a VDJ reference, otherwise de novo.
    pub has_vdj_reference: bool,
    /// Is multiplexed, otherwise singleplex.
    pub is_multiplexed: bool,
    /// Cells are multiplexed, as with Hashtag.
    pub is_cell_multiplexed: bool,
    /// Reads are multiplexed, as with Flex and OCM.
    pub is_read_multiplexed: bool,
    /// Is probe-based RTL, otherwise Universal.
    pub is_rtl: bool,
    /// Has gDNA metrics.
    pub has_gdna: bool,
    /// The parameter include-introns is enabled.
    pub include_introns: bool,
}

#[derive(Deserialize, Debug, Clone)]
struct MetricEtlGroup {
    tier: MetricTier,
    #[serde(default)]
    conditions: MetricConditions,
    /// table that can be used to group the metrics from the individual rows
    /// together. This is a stopgap solution until we refactor tables to not be
    /// table-structured.
    ///
    /// The value associated with the key column will be extracted and added to
    /// every metric in the row when we output flattened JSON metrics.
    ///
    /// This mechanism is also used to populate the first two columns in the
    /// metrics CSV output, so most tables should specify this key even if
    /// they are only going to have one row in the websummary.
    #[serde(default)]
    group_by_key: Option<String>,
    #[serde(flatten)]
    entry_info: BTreeMap<String, MetricEtlConfig>,
}

impl MetricEtlGroup {
    fn validate(&self) -> Result<()> {
        if let Some(group_by_key) = &self.group_by_key {
            ensure!(
                self.entry_info.contains_key(group_by_key),
                "group_by_key {group_by_key} is not present among the keys",
            );
        }
        for cfg in self.entry_info.values() {
            cfg.validate()?;
        }
        Ok(())
    }
}

/// The result of running a metric transformation.
pub struct TransformResult {
    value: Value,
    string_value: String,
    alerts: Vec<AlertSpec>,
}

/// Behavior required to act as a metric transformer.
pub trait Transformer {
    /// Transform a particular metric value.
    ///
    /// In the process, generate alerts based on the provided configuration.
    /// Processing alerts is done here since the value only takes a strong type
    /// inside the context of this method.
    fn transform(
        &self,
        key: &str,
        val: &Value,
        alerts: &[AlertConfig],
        ctx: &AlertContext,
    ) -> Result<TransformResult>;
}

/// A transformer for any value that is parsed as a strong type and then passed on.
#[derive(Default)]
struct PassThroughTransformer<T: MakePretty + TryFromValue + Serialize>(PhantomData<T>);

impl<T: MakePretty + TryFromValue + Serialize> Transformer for PassThroughTransformer<T> {
    fn transform(
        &self,
        key: &str,
        val: &Value,
        alerts: &[AlertConfig],
        ctx: &AlertContext,
    ) -> Result<TransformResult> {
        let val = T::try_from_value(val)?;
        let alerts = JsonMetricSummary::construct_alerts(key, val.as_ref(), alerts, ctx);
        let string_value = val
            .as_ref()
            .map_or_else(|| "---".to_string(), MakePretty::to_string_for_csv);
        Ok(TransformResult {
            value: serde_json::to_value(val).unwrap(),
            string_value,
            alerts,
        })
    }
}

/// A metric transformer that converts the value into a fraction using a fixed denominator.
#[derive(Clone)]
pub struct CountAndPercentTransformer {
    denominator: usize,
}

impl CountAndPercentTransformer {
    /// Divide the incoming metric by the provided denominator.
    pub fn new(denominator: usize) -> Self {
        Self { denominator }
    }
}

impl Transformer for CountAndPercentTransformer {
    fn transform(
        &self,
        key: &str,
        val: &Value,
        alerts: &[AlertConfig],
        ctx: &AlertContext,
    ) -> Result<TransformResult> {
        let val = usize::try_from_value(val)?;
        let Some(val) = val else {
            return Ok(TransformResult {
                value: Value::Null,
                string_value: String::new(),
                alerts: vec![],
            });
        };
        let val = CountAndPercent(PercentMetric::from_parts(
            val as i64,
            self.denominator as i64,
        ));
        let alerts = JsonMetricSummary::construct_alerts(key, Some(&val), alerts, ctx);
        Ok(TransformResult {
            value: serde_json::to_value(val).unwrap(),
            string_value: val.to_string_for_csv(),
            alerts,
        })
    }
}

/// A percent metric transformer that provides the complement of the given percentage.
///
/// In other words, if the input metric is 1%, the output will be 99%.
struct ComplementPercentTransformer;

impl Transformer for ComplementPercentTransformer {
    fn transform(
        &self,
        key: &str,
        val: &Value,
        alerts: &[AlertConfig],
        ctx: &AlertContext,
    ) -> Result<TransformResult> {
        let val = f64::try_from_value(val)?.map(|p| Percent::Float(1.0 - p));
        let alerts = JsonMetricSummary::construct_alerts(key, val.as_ref(), alerts, ctx);
        let string_value = val
            .as_ref()
            .map_or_else(|| "---".to_string(), MakePretty::to_string_for_csv);
        Ok(TransformResult {
            value: serde_json::to_value(val).unwrap(),
            string_value,
            alerts,
        })
    }
}

/// Key lookup for metric transformers.
pub struct Transformers(HashMap<String, Box<dyn Transformer>>);

impl Transformers {
    /// Add a transformer to this collection with the provided key.
    ///
    /// Panic if the key is already present.
    pub fn add<T: Transformer + 'static, K: Into<String>>(&mut self, key: K, f: T) {
        let key = key.into();
        assert!(
            !self.0.contains_key(&key),
            "duplicate transformer key: {key}"
        );
        self.0.insert(key, Box::new(f) as Box<dyn Transformer>);
    }

    /// Get a transformer by key.
    pub fn get(&self, key: &str) -> Option<&dyn Transformer> {
        self.0.get(key).map(|v| &**v)
    }
}

impl Default for Transformers {
    fn default() -> Self {
        let mut trans = Self(HashMap::new());
        trans.add("usize", PassThroughTransformer::<usize>::default());
        trans.add("f64", PassThroughTransformer::<f64>::default());
        trans.add("Percent", PassThroughTransformer::<Percent>::default());
        trans.add("String", PassThroughTransformer::<String>::default());
        trans.add(
            "FloatAsInt",
            PassThroughTransformer::<FloatAsInt>::default(),
        );
        trans.add("ComplementPercent", ComplementPercentTransformer);
        trans
    }
}

/// Extract a value for the provided metric and transform it.
///
/// The metric config's extract mode determines how we handle optionality.
fn extract_and_transform(
    output_key: &str,
    tier: MetricTier,
    section: Section,
    cfg: &MetricEtlConfig,
    alert_ctx: &AlertContext,
    metrics: &TxHashMap<String, Value>,
    transformers: &Transformers,
) -> Result<Option<JsonMetricSummary>> {
    // Use a specified transformer if one is provided; otherwise use the
    // transformer with the same key as the specified output type.
    let transformer = cfg.transformer.as_ref().unwrap_or(&cfg.ty);
    let Some(transformer) = transformers.get(transformer) else {
        bail!("transformer {transformer} not found for metric {output_key}");
    };

    let key = cfg.json_key.as_deref().unwrap_or(output_key);

    let TransformResult {
        value,
        string_value,
        alerts,
    } = match (metrics.get(key), cfg.extract) {
        (None, ExtractMode::Required) => bail!("required metric {key} was not found"),
        (None, ExtractMode::Optional) => {
            return Ok(None);
        }
        (None, ExtractMode::Placeholder) => {
            // For a missing optional metric, fill in a null in the output.
            TransformResult {
                value: Value::Null,
                string_value: "—".to_string(),
                alerts: vec![],
            }
        }
        (Some(v), _) => transformer.transform(key, v, &cfg.alerts, alert_ctx)?,
    };

    Ok(Some(JsonMetricSummary {
        key: output_key.to_string(),
        value,
        string_value,
        category: tier,
        section,
        config: cfg.clone(),
        // Will populate these after all metrics are extracted.
        grouping_key: None,
        grouping_header: None,
        alerts,
    }))
}

fn process_metric_group(
    section: Section,
    group: &MetricEtlGroup,
    active_conditions: &ActiveConditions,
    alert_ctx: &AlertContext,
    metrics: &TxHashMap<String, Value>,
    transformers: &Transformers,
) -> Result<Vec<JsonMetricSummary>> {
    if !group.conditions.matches(active_conditions) {
        return Ok(vec![]);
    }
    let mut grouping = None;
    let mut out = Vec::new();
    for (key, cfg) in &group.entry_info {
        let val = extract_and_transform(
            key,
            group.tier,
            section,
            cfg,
            alert_ctx,
            metrics,
            transformers,
        )?;
        let Some(val) = val else {
            continue;
        };
        if let Some(group_by_key) = &group.group_by_key
            && group_by_key == key
        {
            let grouping_key = val.value.as_str().unwrap_or_else(|| {
                panic!(
                    "invalid metric grouping key value for metric {key}; must be a string, got {}",
                    val.value
                )
            });
            grouping = Some((grouping_key.to_string(), val.config.header.clone()));
        }
        out.push(val);
    }

    if let Some((grouping_val, grouping_header)) = grouping {
        for metric in &mut out {
            metric.grouping_key = Some(grouping_val.clone());
            metric.grouping_header = Some(grouping_header.clone());
        }
    }
    Ok(out)
}

#[cfg(test)]
mod test {
    use super::*;
    use anyhow::Result;
    use insta::assert_json_snapshot;
    use serde_json::json;

    const TEST_GROUPS: &str = include_str!("../../tests/metrics_etl_test.toml");

    #[test]
    fn test_conditions_match() {
        let default_cond = MetricConditions::default();
        let active_cond = ActiveConditions::default();
        assert!(default_cond.matches(&active_cond));
    }

    #[test]
    fn test_load_etls() -> Result<()> {
        load_metrics_etl()?;
        load_metrics_etl_vdj()?;
        load_metrics_etl_special()?;
        Ok(())
    }

    #[test]
    fn test_extract_and_transform() -> Result<()> {
        let (_lib_metrics_proc, mut metrics_proc) = load_and_divide_groups(TEST_GROUPS)?;
        metrics_proc.add_transformer(
            "CountAndPercent",
            CountAndPercentTransformer { denominator: 100 },
        );

        let mut metrics: TxHashMap<String, Value> = TxHashMap::default();
        metrics.insert("num_filtered_bcs".to_string(), json!(10));
        metrics.insert("physical_library_id".to_string(), json!("an ID"));
        metrics.insert("cell_associated_partitions".to_string(), json!(50));
        metrics.insert("singlets_assigned_sample".to_string(), json!(10));
        metrics.insert("partitions_with_no_cmos".to_string(), json!(0));
        metrics.insert("partitions_called_multiplets".to_string(), json!(0.5));
        metrics.insert("multiplet_rate".to_string(), json!(0.05));
        metrics.insert("reads_in_cells".to_string(), json!(0.3));

        let extracted = metrics_proc.process(
            Section::Gex,
            &Default::default(),
            &Default::default(),
            &metrics,
        )?;
        assert_json_snapshot!(extracted);
        Ok(())
    }
}
