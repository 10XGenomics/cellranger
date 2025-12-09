//! cr_aggr::websummary
#![deny(missing_docs)]

use cr_vdj::clonotype_table::VdjAggrClonotypeTable;
use cr_websummary::{
    CardWithTable, GenericTable, HeroMetric, NamedTable, PrettyMetric, TableRow, TermDesc,
    Threshold, TitleWithHelp, TitleWithTermDesc, WsSample,
};
use serde::Serialize;

pub(crate) mod annotation_card;
pub(crate) mod cdr3_table;
pub(crate) mod cells_card;
pub(crate) mod hero_metrics;

#[derive(Debug, Serialize, PartialEq, Eq, Clone)]
#[serde(into = "GenericTable")]
pub struct VdjAggrPipelineInfo {
    pub run_id: String,
    pub run_desc: String,
    pub vdj_ref: String,
    pub pipeline_version: String,
    pub num_samples: usize,
    pub num_donors: usize,
    pub num_origins: usize,
}

impl From<VdjAggrPipelineInfo> for GenericTable {
    fn from(info: VdjAggrPipelineInfo) -> GenericTable {
        GenericTable {
            header: None,
            rows: vec![
                TableRow::two_col("Run ID", info.run_id),
                TableRow::two_col("Run Description", info.run_desc),
                TableRow::two_col("V(D)J Reference", info.vdj_ref),
                TableRow::two_col("Pipeline Version", info.pipeline_version),
                TableRow::two_col("Number of Samples", info.num_samples),
                TableRow::two_col("Number of Donors", info.num_donors),
                TableRow::two_col("Number of Origins", info.num_origins),
            ],
            grouping_header: None,
        }
    }
}

#[derive(Debug, Serialize)]
pub struct VdjAggrWsContent {
    pub sample: WsSample,
    pub summary_tab: VdjAggrWsSummaryTab,
}

#[derive(Debug, Serialize)]
pub struct VdjAggrWsSummaryTab {
    pub hero_metrics: hero_metrics::VdjAggrHeroMetrics,
    pub pipeline_info_table: VdjAggrPipelineInfo,
    pub vdj_clonotype: VdjAggrClonotypeTable,
    pub vdj_annotation: annotation_card::VdjAggrAnnotationTable,
    pub vdj_aggr_cells: cells_card::VdjAggrCellsTable,
    pub vdj_shared_cdr3: cdr3_table::VdjAggrSharedCdr3,
    pub vdj_clonotype_hist: cr_vdj::clonotype_hist::ClonotypeHist,
}

#[cfg(test)]
mod tests {
    use super::*;
    use cr_vdj::{check_eq_json, test_json_roundtrip};
    use pretty_assertions::assert_eq;

    #[test]
    fn test_ws_sample_parse() {
        test_json_roundtrip::<WsSample>(
            r#"{
                "subcommand": "aggr",
                "command": "Cell Ranger",
                "id": "Human PBMC BCR",
                "description": "Pre vs post vaccination"
            }"#,
        );
    }

    #[test]
    fn test_metrics_parse() {
        let m = test_json_roundtrip::<HeroMetric>(
            r#"{
                "threshold": "pass",
                "metric": "18,000",
                "name": "Total number of Cells assigned to a clonotype"
            }"#,
        );
        assert_eq!(
            m,
            HeroMetric {
                threshold: Threshold::Pass,
                name: "Total number of Cells assigned to a clonotype".into(),
                metric: PrettyMetric::integer(18000)
            }
        );
    }

    #[test]
    fn test_header_with_help() {
        test_json_roundtrip::<TitleWithHelp>(
            r#"{
            "helpText": "This is the help text",
            "title": "This is the title"
        }"#,
        );
    }

    #[test]
    fn test_term_descriptions() {
        test_json_roundtrip::<TermDesc>(
            r#"[
                "Inner title",
                [
                    "Inner help description 1",
                    "Inner help description 2"
                ]
            ]"#,
        );
    }

    #[test]
    fn test_dyn_help() {
        test_json_roundtrip::<TitleWithTermDesc>(
            r#"{
                "data": [
                    [
                        "Metric 1",
                        [
                            "Help 1"
                        ]
                    ],
                    [
                        "Metric 2",
                        [
                            "Help 2"
                        ]
                    ]
                ],
                "title": "Title text"
            }"#,
        );
    }

    #[test]
    fn test_generic_table() {
        test_json_roundtrip::<GenericTable>(
            r#"{
                "rows": [
                    [
                        "Sample ID",
                        "Human PBMC BCR"
                    ],
                    [
                        "Sample Description",
                        "Pre vs post vaccination"
                    ]
                ]
            }
            "#,
        );
    }

    #[test]
    fn test_generic_table_with_header() {
        test_json_roundtrip::<GenericTable>(
            r#"{
                "header": [
                    "Donor",
                    "Origin",
                    "Cells",
                    "Clonotypes"
                ],
                "rows": [
                    [
                        "Donor1",
                        "PreVac",
                        "10,000",
                        "7,000"
                    ],
                    [
                        "Donor1",
                        "PostVac",
                        "8,000",
                        "2,000"
                    ]
                ]
            }"#,
        );
    }

    #[test]
    fn test_pipeline_info() {
        let info = VdjAggrPipelineInfo {
            run_id: "Human PBMC BCR".into(),
            run_desc: "Pre vs post vaccination".into(),
            vdj_ref: "vdj_GRCh38_alts_ensembl-4.0.0".into(),
            pipeline_version: "4.1.0".into(),
            num_samples: 3,
            num_donors: 1,
            num_origins: 2,
        };

        check_eq_json(
            &serde_json::to_string(&info).unwrap(),
            r#"{
                "rows": [
                    [
                        "Run ID",
                        "Human PBMC BCR"
                    ],
                    [
                        "Run Description",
                        "Pre vs post vaccination"
                    ],
                    [
                        "V(D)J Reference",
                        "vdj_GRCh38_alts_ensembl-4.0.0"
                    ],
                    [
                        "Pipeline Version",
                        "4.1.0"
                    ],
                    [
                        "Number of Samples",
                        "3"
                    ],
                    [
                        "Number of Donors",
                        "1"
                    ],
                    [
                        "Number of Origins",
                        "2"
                    ]
                ]
            }"#,
        );
    }

    #[test]
    fn test_pretty_metric() {
        assert_eq!(PrettyMetric::integer(100).0, "100");
        assert_eq!(PrettyMetric::integer(1000).0, "1,000");
        assert_eq!(PrettyMetric::percent(0.10).0, "10.0%");
        assert_eq!(PrettyMetric::percent(0.114567).0, "11.5%");
        assert_eq!(PrettyMetric::decimal(0.114567).0, "0.11");
        assert_eq!(PrettyMetric::decimal(11.4567).0, "11.46");
        assert_eq!(PrettyMetric::decimal(10011.4567).0, "10,011.46");
        assert_eq!(PrettyMetric::decimal(2.99).0, "2.99");
        assert_eq!(PrettyMetric::decimal(2.999).0, "3.00");
        assert_eq!(PrettyMetric::decimal(-2123.999).0, "-2,124.00");
    }

    #[test]
    fn test_config_valid_json() {
        let _ = cr_websummary::default_plotly_config();
    }
}
