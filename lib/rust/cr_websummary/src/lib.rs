pub mod alert;
pub mod multi;
pub mod value;

extern crate self as cr_websummary;

use alert::{Alert, AlertContext, AlertLevel, AlertSpec};
use anyhow::Result;
use metric::{PercentMetric, TxHashMap};
use multi::websummary::ToCsvRows;
use plotly::common::{Anchor, Title};
use plotly::layout::{Axis, AxisType, HoverMode, Legend};
use plotly::Layout;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::default::Default;
use std::fmt::{Display, Formatter};
use std::iter::{zip, FromIterator};
use std::str::FromStr;
use thousands::Separable;

const TABLE_HEADER_METRIC_KEYS: [&str; 7] = [
    "Physical library ID",
    "Fastq ID",
    "CMO Name",
    "Probe barcode ID",
    "Genome",
    "Targeting Status",
    "Feature Type",
];

#[derive(Debug, Serialize, Deserialize, Clone)]
pub enum WsCommand {
    #[serde(rename = "Cell Ranger")]
    Cellranger,
}

#[derive(Debug, Serialize, Deserialize, Clone)]
#[serde(rename_all = "lowercase")]
pub enum WsSubCommand {
    Aggr,
    Multi,
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct WsSample {
    pub subcommand: WsSubCommand,
    pub command: WsCommand,
    pub id: String,
    pub description: String,
}

impl WsSample {
    pub fn aggr(id: String, description: String) -> Self {
        WsSample {
            subcommand: WsSubCommand::Aggr,
            command: WsCommand::Cellranger,
            id,
            description,
        }
    }
    pub fn multi(id: String, description: String) -> Self {
        WsSample {
            subcommand: WsSubCommand::Multi,
            command: WsCommand::Cellranger,
            id,
            description,
        }
    }
}

#[derive(Debug, Clone, Serialize)]
pub struct Tab<C: Alert + ToCsvRows> {
    pub content: C,
    alerts: Vec<AlertSpec>,
}

impl<C: Alert + ToCsvRows> Tab<C> {
    pub fn new(content: C, context: &AlertContext) -> Self {
        let mut alerts = content.alerts(context);
        // All tabs should have an alert if preflight was skipped.
        if context.no_preflight {
            alerts.push(AlertSpec {
                level: AlertLevel::Warn,
                title: "Analysis preflight checks were skipped".to_string(),
                formatted_value: String::default(),
                message: "Your analysis was run without preflight checks and may have used non-standard settings. Please carefully review your results.".to_string(),
            })
        }
        Tab { content, alerts }
    }
}

impl<C: Alert + ToCsvRows> ToCsvRows for Tab<C> {
    fn to_csv_rows(self) -> Vec<Vec<String>> {
        self.content.to_csv_rows()
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "lowercase")]
pub enum Threshold {
    Pass,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
#[serde(untagged)]
pub enum Percent {
    Metric(PercentMetric),
    Float(f64),
}

impl Default for Percent {
    fn default() -> Percent {
        Percent::Metric(PercentMetric::default())
    }
}

impl FromStr for Percent {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Percent> {
        let f = s.parse::<f64>()?;
        Ok(Percent::Float(f))
    }
}

impl From<PercentMetric> for Percent {
    fn from(src: PercentMetric) -> Self {
        Percent::Metric(src)
    }
}

impl From<Percent> for f64 {
    fn from(metric: Percent) -> f64 {
        match metric {
            Percent::Metric(metric) => metric.fraction().unwrap_or(f64::NAN),
            Percent::Float(f) => f,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct PrettyMetric(pub String);
impl Display for PrettyMetric {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}
impl PrettyMetric {
    pub fn integer(src: impl Display) -> Self {
        PrettyMetric(src.separate_with_commas())
    }
    pub fn percent(m: f64) -> Self {
        PrettyMetric(format!("{:.2}%", 100.0 * m))
    }
    pub fn percent_f1(m: f64) -> Self {
        PrettyMetric(format!("{:.1}%", 100.0 * m))
    }
    pub fn count_and_percent(m: PercentMetric) -> Self {
        match m.fraction() {
            Some(f) => PrettyMetric(format!(
                "{} ({})",
                m.numerator.count().separate_with_commas(),
                PrettyMetric::percent(f)
            )),
            None => PrettyMetric(format!("{} (---%)", m.numerator.count())),
        }
    }
    pub fn decimal(m: f64) -> Self {
        let without_comma = format!("{m:.2}");
        let parts: Vec<_> = without_comma.split('.').collect();
        PrettyMetric(format!(
            "{}.{}",
            PrettyMetric::integer(parts[0].parse::<i64>().unwrap()),
            parts[1]
        ))
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct HeroMetric {
    pub threshold: Threshold,
    pub name: String,
    pub metric: PrettyMetric,
}

/// Usually used to attach heading to a card with a help snippet
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct TitleWithHelp {
    #[serde(rename = "helpText")]
    pub help: String,
    pub title: String,
}

/// Description of a specific term or metric under the collapsible help
/// First element of the tuple is the term and the second element is the description
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct TermDesc(pub String, pub Vec<String>);
impl TermDesc {
    pub fn with_one_desc(term: impl ToString, desc: impl ToString) -> Self {
        TermDesc(term.to_string(), vec![desc.to_string()])
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct TitleWithTermDesc {
    pub title: String,
    pub data: Vec<TermDesc>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct TableRow(pub Vec<String>);

impl TableRow {
    pub fn two_col(c1: impl ToString, c2: impl ToString) -> Self {
        TableRow(vec![c1.to_string(), c2.to_string()])
    }
}

impl From<Vec<String>> for TableRow {
    fn from(item: Vec<String>) -> Self {
        TableRow(item)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct GenericTable {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub header: Option<Vec<String>>,
    pub rows: Vec<TableRow>,
}

impl ToCsvRows for GenericTable {
    fn to_csv_rows(self) -> Vec<Vec<String>> {
        let mut csv_rows: Vec<Vec<String>> = vec![];

        let table_header = self
            .header
            .expect("No header for table being converted to CSV.");
        if table_header.is_empty() {
            return csv_rows;
        }
        if TABLE_HEADER_METRIC_KEYS
            .iter()
            .any(|x| *x == table_header[0])
        {
            // in this case, the metric table has a key in the first column such as "Physical Library ID"
            // we write one CSV column that says what type of ID/key it is (Physical library ID, Fastq ID, etc)
            // and another CSV column that specifies that ID/key itself
            let metric_key_type: String = table_header[0].clone(); // e.g. "Physical library ID"
            for table_row in self.rows {
                let metric_key: String = table_row.0[0].clone(); // e.g. "GEX_1"
                                                                 // header value and metric value for rest of columns (non-key columns)
                for (metric_name, metric_val) in zip(&table_header[1..], &table_row.0[1..]) {
                    let csv_row: Vec<String> = vec![
                        metric_key_type.clone(),
                        metric_key.clone(),
                        metric_name.clone(),
                        metric_val.clone(),
                    ];
                    csv_rows.push(csv_row);
                }
            }
        } else {
            // in this case, the table has no "key" in the first column and the metrics are free-standing.
            for table_row in self.rows {
                // header value and metric value for rest of columns (non-key columns)
                for (metric_name, metric_val) in zip(&table_header, &table_row.0) {
                    let csv_row: Vec<String> = vec![
                        String::default(),
                        String::default(),
                        metric_name.clone(),
                        metric_val.clone(),
                    ];
                    csv_rows.push(csv_row);
                }
            }
        }

        csv_rows
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct NamedTable {
    pub name: String,
    pub table: GenericTable,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct CardWithTable {
    pub table: GenericTable,
    pub help: TitleWithHelp,
}

impl Alert for CardWithTable {}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct CardWithMetric {
    pub table: GenericTable,
    pub help: TitleWithTermDesc,
}

#[derive(Debug, Clone, Serialize, PartialEq, Eq)]
#[serde(transparent)]
pub struct MetricCard<C: Into<CardWithMetric>> {
    card: CardWithMetric,
    #[serde(skip)]
    raw: C,
}

impl<C: Into<CardWithMetric>> ToCsvRows for MetricCard<C> {
    fn to_csv_rows(self) -> Vec<Vec<String>> {
        self.card.table.to_csv_rows()
    }
}

impl<C: Clone + Into<CardWithMetric>> From<C> for MetricCard<C> {
    fn from(src: C) -> Self {
        MetricCard {
            card: src.clone().into(),
            raw: src,
        }
    }
}

impl<C: Alert + Into<CardWithMetric>> Alert for MetricCard<C> {
    fn alerts(&self, ctx: &AlertContext) -> Vec<AlertSpec> {
        self.raw.alerts(ctx)
    }
}

const DEFAULT_PLOTLY_CONFIG: &str = r#"{
    "displayModeBar": true,
    "staticPlot": false,
    "dragmode": "zoom",
    "modeBarButtons": [
        [
            "toImage"
        ]
    ]
}"#;

pub fn default_plotly_config() -> Value {
    serde_json::from_str::<Value>(DEFAULT_PLOTLY_CONFIG).unwrap()
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct ChartWithHelp {
    pub plot: PlotlyChart,
    pub help: TitleWithHelp,
}

impl Alert for ChartWithHelp {}
impl ToCsvRows for ChartWithHelp {}

// a titled plot where the plot can have
// whatever data structure we want
// this is useful for e.g. the JIBES biplot
#[derive(Debug, Clone, Deserialize, Serialize, PartialEq, Eq)]
pub struct RawChartWithHelp {
    pub plot: Value,
    pub help: TitleWithHelp,
}

impl Alert for RawChartWithHelp {}
impl ToCsvRows for RawChartWithHelp {}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, Default)]
pub struct PlotlyChart {
    pub config: Value,
    pub data: Vec<Value>,
    pub layout: Value,
}

impl PlotlyChart {
    pub fn with_layout_and_data<D: Serialize>(layout: Layout, data: Vec<D>) -> Self {
        PlotlyChart {
            config: default_plotly_config(),
            data: data
                .into_iter()
                .map(|d| serde_json::to_value(&d).unwrap())
                .collect(),
            layout: serde_json::to_value(&layout).unwrap(),
        }
    }

    pub fn new_line_plot(
        x_data: Vec<f64>,
        y_data: Vec<f64>,
        x_label: String,
        y_label: String,
    ) -> Self {
        let data: TxHashMap<String, Vec<f64>> = TxHashMap::from_iter(
            [("x".to_string(), x_data), ("y".to_string(), y_data)]
                .iter()
                .cloned(),
        );
        let layout = Layout::new()
            .show_legend(true)
            .margin(plotly::layout::Margin::new().right(40).left(60).top(0))
            .hover_mode(HoverMode::Closest)
            .x_axis(
                Axis::new()
                    .type_(AxisType::Category)
                    .title(Title::new(&x_label)),
            )
            .y_axis(Axis::new().title(Title::new(&y_label)))
            .legend(
                Legend::new()
                    .y_anchor(Anchor::Top)
                    .y(0.99)
                    .x_anchor(Anchor::Left)
                    .x(0.01)
                    .background_color("#ffffff"),
            );

        Self::with_layout_and_data(layout, vec![data])
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
pub struct CountAndPercent(pub PercentMetric);

// a wrapper type for Percents with only 1 decimal place precision
#[derive(Debug, Clone, Copy, Serialize, Deserialize, Default, PartialEq)]
pub struct PercentF1(pub Percent);

impl From<PercentMetric> for PercentF1 {
    fn from(src: PercentMetric) -> Self {
        PercentF1(Percent::Metric(src))
    }
}

pub trait MakePretty {
    fn make_pretty(&self) -> String;
    fn as_f64(&self) -> f64;
}

impl MakePretty for Vec<String> {
    fn make_pretty(&self) -> String {
        self.join("<br>")
    }
    fn as_f64(&self) -> f64 {
        unreachable!("Vec<String> cannot be cast into f64")
    }
}

impl MakePretty for String {
    fn make_pretty(&self) -> String {
        self.clone()
    }
    fn as_f64(&self) -> f64 {
        unreachable!("String cannot be cast into f64")
    }
}

impl MakePretty for usize {
    fn make_pretty(&self) -> String {
        PrettyMetric::integer(*self).0
    }
    fn as_f64(&self) -> f64 {
        *self as f64
    }
}

impl MakePretty for Percent {
    fn make_pretty(&self) -> String {
        match &self {
            Percent::Metric(m) => match m.fraction() {
                Some(f) => PrettyMetric::percent(f).0,
                None => "---%".to_string(),
            },
            Percent::Float(f) => PrettyMetric::percent(*f).0,
        }
    }
    fn as_f64(&self) -> f64 {
        match &self {
            Percent::Metric(m) => match m.fraction() {
                Some(f) => f,
                None => f64::NAN,
            },
            Percent::Float(f) => *f,
        }
    }
}

impl MakePretty for CountAndPercent {
    fn make_pretty(&self) -> String {
        PrettyMetric::count_and_percent(self.0).0
    }
    fn as_f64(&self) -> f64 {
        self.0.numerator.count() as f64
    }
}

impl MakePretty for PercentF1 {
    fn make_pretty(&self) -> String {
        match &self.0 {
            Percent::Metric(m) => match m.fraction() {
                Some(f) => PrettyMetric::percent_f1(f).0,
                None => "---%".to_string(),
            },
            Percent::Float(f) => PrettyMetric::percent(*f).0,
        }
    }
    fn as_f64(&self) -> f64 {
        match &self.0 {
            Percent::Metric(m) => match m.fraction() {
                Some(f) => f,
                None => f64::NAN,
            },
            Percent::Float(f) => *f,
        }
    }
}

impl MakePretty for f64 {
    fn make_pretty(&self) -> String {
        PrettyMetric::decimal(*self).0
    }
    fn as_f64(&self) -> f64 {
        *self
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
pub struct FloatAsInt(pub f64);

impl MakePretty for FloatAsInt {
    fn make_pretty(&self) -> String {
        (self.0.round() as usize).make_pretty()
    }
    fn as_f64(&self) -> f64 {
        self.0
    }
}

///
/// The macro below allows us to declare a new card with a web summary table inside
/// in terms of a Table struct and a Row struct
/// The table serializes to JSON in the correct format for the web summary.
/// The Row can be converted into a generic TableRow (holding strings)
/// The table struct is defined as a vector of Rows. Conversion into generic CardWithTable is implemented
/// After conversion, the generic table type can be serialized to JSON
/// The MakePretty trait is how each Row element knows how to convert itself into a neatly formatted String
///
/// Example invocation:
///
///card_with_table!(
///    table_name: VdjAggrClonotypeTable,
///    row_name: VdjAggrClonotypeRow,
///    title: "Top 10 Clonotype CDR3 Sequences",
///    help_text: "TODO",
///    columns: {
///        id: usize : "Clonotype ID",
///        cdr3_aas : Vec<String> : "CDR3s",
///        num_cells: usize : "Frequency",
///        proportion: f64 : "Proportion"
///    }
///);
///
/// The implementation of the macro is sort of gnarly and that's because of the variable number of
/// columns and the fact that Rust macros can't return "partial" syntax elements (e.g. part of a struct).
/// Therefore things like the Row struct definition, the header formatting, and the Row printing
/// must be done with a "muncher" pattern where the desired things are built up recursively in the
/// macro inputs and then turned into the desired struct or vec all at once.
///
#[macro_export]
macro_rules! card_with_table {
    (
        table_name: $table_name:ident,
        row_name: $row_name:ident,
        title: $card_title:literal,
        help_text: $help_text:literal,
        columns: {
            $($col_id:ident : $col_ty:ty : $col_header:literal,)*
        }
    ) => {
        #[derive(Debug, Clone, Serialize, PartialEq)]
        #[serde(into = "CardWithTable")]
        pub struct $table_name(pub Vec<$row_name>);

        /*impl Default for $table_name {
            let mut v = Vec::new();
            v.push($row_name::default());
            fn default() -> $table_name{0: v}
        }*/

        impl From<$table_name> for CardWithTable {
            fn from(src: $table_name) -> CardWithTable {
                let table = GenericTable {
                    header: Some(vec![$($col_header.to_string(),)*]),
                    rows: src.0.into_iter().map(|row| row.into()).collect(),
                };
                let help = TitleWithHelp {
                    title: $card_title.to_string(),
                    help: $help_text.to_string(),
                };

                CardWithTable { table, help }
            }
        }

        // One row within the table
        #[derive(Debug, Clone, Serialize, PartialEq, Default)]
        #[serde(into = "TableRow")]
        pub struct $row_name {
            $(pub $col_id: $col_ty,)*
        }

        impl From<$row_name> for TableRow {
            fn from(src: $row_name) -> TableRow {
                TableRow(vec![$(src.$col_id.make_pretty(),)*])
            }
        }
    };
    (
        table_name: $table_name:ident,
        row_name: $row_name:ident,
        title: $card_title:literal,
        help_text: $help_text:literal,
        columns: {
            $($col_id:ident : $col_ty:ty : $col_header:literal),*
        }
    ) => {
        card_with_table!(
            table_name: $table_name,
            row_name: $row_name,
            title: $card_title,
            help_text: $help_text,
            columns: {
                $($col_id : $col_ty : $col_header,)*
            }
        );
    };
}

#[cfg(test)]
#[ctor::ctor]
fn init() {
    // this ensures insta knows where to find its snap tests
    let cwd = std::env::current_dir().unwrap();
    let workspace_root = cwd.parent().unwrap();
    std::env::set_var("INSTA_WORKSPACE_ROOT", workspace_root);
}
