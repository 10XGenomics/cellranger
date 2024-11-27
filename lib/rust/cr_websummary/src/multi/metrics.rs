//! Process a specification file into a metrics ETL operation.
//! NOTE: this currently parses the tables.toml format.
//! This will be refactored along with the format of that file to be intrinsically
//! flat, without any explicit table structure.

use super::websummary::JsonMetricSummary;
use crate::value::TryFromValue;
use crate::{FloatAsInt, Percent};
use anyhow::{bail, Result};
use cr_types::websummary::MetricConfig;
use cr_types::LibraryType;
use serde::Deserialize;
use serde_json::Value;
use std::collections::{BTreeMap, HashMap, HashSet};
use toml;

#[derive(Deserialize, Debug)]
#[serde(rename_all = "lowercase")]
pub enum MetricTier {
    Library,
    Sample,
}
#[derive(Default, Debug, Deserialize)]
pub struct MetricConditions {
    pub library_types: Vec<LibraryType>,
    pub is_cmo_multiplexed: Option<bool>,
    /// Covers both RTL and OCM multiplexing.
    pub is_read_multiplexed: Option<bool>,
    pub is_rtl: Option<bool>,
    pub has_gdna: Option<bool>,
    pub include_introns: Option<bool>,
}

fn matches_option(opt: Option<bool>, val: bool) -> bool {
    opt.map_or(true, |v| v == val)
}

impl MetricConditions {
    pub fn matches(&self, cond: &ActiveConditions) -> bool {
        self.matches_library_types(&cond.library_types)
            && matches_option(self.is_cmo_multiplexed, cond.is_cmo_multiplexed)
            && matches_option(self.is_read_multiplexed, cond.is_read_multiplexed)
            && matches_option(self.is_rtl, cond.is_rtl)
            && matches_option(self.has_gdna, cond.has_gdna)
    }

    fn matches_library_types(&self, types: &HashSet<LibraryType>) -> bool {
        if self.library_types.is_empty() {
            return true;
        }
        self.library_types
            .iter()
            .any(|lib_type| types.contains(lib_type))
    }
}

#[derive(Default)]
pub struct ActiveConditions {
    pub library_types: HashSet<LibraryType>,
    pub is_cmo_multiplexed: bool,
    pub is_read_multiplexed: bool,
    pub is_rtl: bool,
    pub has_gdna: bool,
    pub include_introns: bool,
}

#[allow(unused)]
#[derive(Deserialize, Debug)]
pub struct CardWithTableToml {
    title: String,
    help: Option<String>,
    tier: MetricTier,
    #[serde(default)]
    conditions: MetricConditions,
    entries: Vec<String>,
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
    entry_info: BTreeMap<String, MetricConfig>,
}

impl CardWithTableToml {
    fn validate(&self) {
        for e in &self.entries {
            assert!(
                self.entry_info.contains_key(e),
                "Information for {} does not exist in table '{}'",
                e,
                self.title
            );
        }
        assert!(
            self.entries.len() == self.entry_info.len(),
            "Duplicate entries listed under table '{}'",
            &self.title
        );
        if let Some(group_by_key) = &self.group_by_key {
            assert!(
                self.entries.contains(group_by_key),
                "group_by_key {group_by_key} is not present among the entries for table '{}'",
                self.title
            );
        }
        for cfg in self.entry_info.values() {
            cfg.validate().unwrap();
        }
    }
}

pub fn load_tables() {
    let content = include_str!("tables.toml");

    let parsed: BTreeMap<String, CardWithTableToml> = toml::from_str(content).unwrap();

    for table in parsed.values() {
        table.validate();
        process_table_metrics(
            table,
            &Default::default(),
            &HashMap::new(),
            &Transformers::default(),
        )
        .unwrap();
    }

    // let mut keys: HashMap<String, String> = HashMap::new();

    // for (table_key, table) in parsed {
    //     table.validate();
    //     for (key, _entry) in table.entry_info {
    //         if keys.contains_key(&key) {
    //             println!(
    //                 "duplicate key: {key} from tables: {table_key}, {}",
    //                 keys[&key]
    //             );
    //         }
    //         keys.insert(key.clone(), table_key.clone());
    //     }
    // }
}

pub type Transformer = Box<dyn Fn(&Value) -> Result<Value>>;

pub struct Transformers(HashMap<String, Transformer>);

impl Transformers {
    pub fn add(&mut self, key: &str, f: impl Fn(&Value) -> Result<Value> + 'static) {
        self.0.insert(key.to_string(), Box::new(f));
    }

    pub fn get(&self, key: &str) -> Option<&Transformer> {
        self.0.get(key)
    }
}

impl Default for Transformers {
    fn default() -> Self {
        let mut trans = Self(HashMap::new());
        trans.add("usize", |v| {
            Ok(serde_json::json!(usize::try_from_value(v)?))
        });
        trans.add("f64", |v| Ok(serde_json::json!(f64::try_from_value(v)?)));
        trans.add("Percent", |v| {
            Ok(serde_json::json!(Percent::try_from_value(v)?))
        });
        trans.add("String", |v| {
            Ok(serde_json::json!(String::try_from_value(v)?))
        });
        trans.add("FloatAsInt", |v| {
            Ok(serde_json::json!(FloatAsInt::try_from_value(v)?))
        });
        trans.add("CountAndPercent", |_| {
            panic!("This extractor must be customized for each use case.");
        });
        trans
    }
}

pub fn extract_and_transform(
    key: &str,
    cfg: &MetricConfig,
    _conditions: &MetricConditions,
    _active_conditions: &ActiveConditions,
    metrics: &HashMap<String, Value>,
    transformers: &Transformers,
) -> Result<Option<JsonMetricSummary>> {
    // if !conditions.matches(active_conditions) {
    //     return None;
    // }

    let Some(transformer) = transformers.get(&cfg.ty) else {
        bail!("transformer {} not found", cfg.ty);
    };

    if !transformers.0.is_empty() {
        return Ok(None);
    }

    let key = cfg.json_key.clone().unwrap_or_else(|| key.to_string());
    let val = match metrics.get(&key) {
        None if !cfg.optional => bail!("required metric {key} was not found"),
        None => {
            return Ok(None);
        }
        Some(v) => v,
    };

    Ok(Some(JsonMetricSummary {
        key,
        value: transformer(val)?,
        category: String::new(),
        library_type: String::new(),
        config: cfg.clone(),
        // Will populate this after all metrics are extracted.
        grouping_key: None,
        // TODO: populate alerts
        alerts: Default::default(),
    }))
}

pub fn process_table_metrics(
    table: &CardWithTableToml,
    active_conditions: &ActiveConditions,
    metrics: &HashMap<String, Value>,
    transformers: &Transformers,
) -> Result<Vec<JsonMetricSummary>> {
    let mut grouping_val = None;
    let mut out = Vec::new();
    for (key, cfg) in &table.entry_info {
        let val = extract_and_transform(
            key,
            cfg,
            &table.conditions,
            active_conditions,
            metrics,
            transformers,
        )?;
        let Some(val) = val else {
            continue;
        };
        if let Some(group_by_key) = &table.group_by_key {
            if group_by_key == key {
                grouping_val = Some(val.value.clone());
            }
        }
        out.push(val);
    }

    if let Some(grouping_val) = grouping_val {
        for metric in &mut out {
            metric.grouping_key = Some(grouping_val.clone());
        }
    }
    Ok(out)
}
