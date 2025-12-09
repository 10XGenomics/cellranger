#![expect(missing_docs)]
use super::websummary::MetricsTraitWrapper;
use crate::alert::AlertLevel;
use crate::{
    Alert, AlertContext, AlertSpec, ChartWithHelp, GenericTable, PlotlyChart, TableRow,
    TitleWithHelp,
};
use cr_types::{BarcodeMultiplexingType, LibraryType};
use serde::{Deserialize, Serialize};
use serde_json::value::Value;
use websummary_derive::Alert;

#[derive(Debug, Serialize, Deserialize, PartialEq, Eq, Clone, Default)]
#[serde(into = "GenericTable")]
pub struct VdjParametersTable {
    pub chemistry: String,
    pub vdj_reference: Option<String>,
    pub vdj_reference_path: Option<String>,
    pub gamma_delta: bool,
    pub denovo: bool,
}

impl From<VdjParametersTable> for GenericTable {
    fn from(info: VdjParametersTable) -> GenericTable {
        let VdjParametersTable {
            chemistry,
            vdj_reference,
            vdj_reference_path,
            gamma_delta: _, // Only used for alert
            denovo: _,      // Only used for alert
        } = info;
        let rows = vec![
            TableRow::two_col("Chemistry", chemistry),
            TableRow::two_col("V(D)J Reference", vdj_reference.unwrap_or_default()),
            TableRow::two_col(
                "V(D)J Reference Path",
                vdj_reference_path.unwrap_or_default(),
            ),
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
            message: r"Multiplexing performance cannot be guaranteed".into(),
            });
        }
        if self.denovo {
            alerts.push(AlertSpec {
                level: AlertLevel::Warn,
                title: "Denovo workflow used".to_string(),
                formatted_value: String::default(),
                message: "Denovo is a non-standard VDJ analysis workflow, if run without a reference only a subset of outputs and metrics are generated.".to_string(),
            });
        }

        alerts
    }
}

#[derive(Serialize, Clone, Default)]
pub struct VdjDiagnostics {
    pub filter_metrics: Option<Value>,
}

#[derive(Serialize, Deserialize, Clone)]
pub struct ClonotypeInfo {
    pub table: GenericTable,
    pub plot: PlotlyChart,
    pub help: TitleWithHelp,
}

impl Alert for ClonotypeInfo {}

#[derive(Serialize, Clone, Alert)]
pub struct SampleVdjWebSummary {
    pub clonotype_info: Option<ClonotypeInfo>,
    pub barcode_rank_plot: Option<ChartWithHelp>,
    #[serde(skip)]
    pub metrics: MetricsTraitWrapper,
}

impl SampleVdjWebSummary {
    pub fn remove_large_data(&mut self) {
        self.clonotype_info = None;
    }
}

#[derive(Serialize, Clone, Alert)]
pub struct LibraryVdjWebSummary {
    pub parameters_table: VdjParametersTable,
    pub barcode_rank_plot: Option<ChartWithHelp>, // None if there are 0 cells
    #[serde(skip)]
    pub metrics: MetricsTraitWrapper,
}
