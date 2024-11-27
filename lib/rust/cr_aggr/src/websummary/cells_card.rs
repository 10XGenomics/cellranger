use crate::websummary::{CardWithTable, GenericTable, PrettyMetric, TableRow, TitleWithHelp};
use serde::Serialize;
use std::convert::Into;

const HELP_TEXT: &str = "This table displays summary statistics for the aggregated samples on the basis of donor and origin. Each line lists the number of cells, number of clonotypes, and clonotypic diversity of the donor/origin combination. Note that the cell counts provided here may be slightly different from those in the individual pipelines due to the nature of the filters applied.";
const CARD_TITLE: &str = "Cells";

/// Hold all the rows. Where serialized to json, this will produce the appropriate structure
/// required for the web summary. See the tests below for the structure
#[derive(Debug, Serialize, PartialEq, Clone)]
#[serde(into = "CardWithTable")]
pub struct VdjAggrCellsTable(pub Vec<VdjAggrCellsRow>);

impl From<VdjAggrCellsTable> for CardWithTable {
    fn from(src: VdjAggrCellsTable) -> CardWithTable {
        let table = GenericTable {
            header: Some(vec![
                "Sample ID".into(),
                "Donor".into(),
                "Origin".into(),
                "Cells".into(),
                "Clonotypes".into(),
                "Diversity".into(),
            ]),
            rows: src.0.into_iter().map(Into::into).collect(),
            grouping_header: None,
        };
        let help = TitleWithHelp {
            title: CARD_TITLE.into(),
            help: HELP_TEXT.into(),
        };

        CardWithTable { table, help }
    }
}

/// One row within the table
#[derive(Debug, Serialize, PartialEq, Clone)]
#[serde(into = "TableRow")]
pub struct VdjAggrCellsRow {
    pub sample_id: String,
    pub donor: String,
    pub origin: String,
    pub cells: usize,
    pub clonotypes: usize,
    pub diversity: f64,
}

impl From<VdjAggrCellsRow> for TableRow {
    fn from(src: VdjAggrCellsRow) -> TableRow {
        TableRow(vec![
            src.sample_id,
            src.donor,
            src.origin,
            PrettyMetric::integer(src.cells as i64).0,
            PrettyMetric::integer(src.clonotypes as i64).0,
            PrettyMetric::decimal(src.diversity).0,
        ])
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::websummary::check_eq_json;

    #[test]
    fn test_vdj_aggr_cells_table_row() {
        let row = VdjAggrCellsRow {
            sample_id: "Sample1".into(),
            donor: "Donor1".into(),
            origin: "PreVac".into(),
            cells: 10000,
            clonotypes: 7000,
            diversity: 6123.4512f64,
        };
        check_eq_json(
            &serde_json::to_string(&row).unwrap(),
            r#"[
                "Sample1",
                "Donor1",
                "PreVac",
                "10,000",
                "7,000",
                "6,123.45"
            ]"#,
        );
    }

    #[test]
    fn test_vdj_aggr_cells_card() {
        let card = VdjAggrCellsTable(vec![
            VdjAggrCellsRow {
                sample_id: "Sample1".into(),
                donor: "Donor1".into(),
                origin: "PreVac".into(),
                cells: 10000,
                clonotypes: 7000,
                diversity: 6123.4512f64,
            },
            VdjAggrCellsRow {
                sample_id: "Sample1".into(),
                donor: "Donor1".into(),
                origin: "PostVac".into(),
                cells: 8000,
                clonotypes: 2000,
                diversity: 110.4167f64,
            },
        ]);
        check_eq_json(
            &serde_json::to_string(&card).unwrap(),
            r#"{
                "table": {
                    "header": [
                        "Sample ID",
                        "Donor",
                        "Origin",
                        "Cells",
                        "Clonotypes",
                        "Diversity"
                    ],
                    "rows": [
                        [
                            "Sample1",
                            "Donor1",
                            "PreVac",
                            "10,000",
                            "7,000",
                            "6,123.45"
                        ],
                        [
                            "Sample1",
                            "Donor1",
                            "PostVac",
                            "8,000",
                            "2,000",
                            "110.42"
                        ]
                    ]
                },
                "help": {
                    "helpText": "This table displays summary statistics for the aggregated samples on the basis of donor and origin. Each line lists the number of cells, number of clonotypes, and clonotypic diversity of the donor/origin combination. Note that the cell counts provided here may be slightly different from those in the individual pipelines due to the nature of the filters applied.",
                    "title": "Cells"
                }
            }"#,
        );
    }
}
