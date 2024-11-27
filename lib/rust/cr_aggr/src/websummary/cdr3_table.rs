use crate::websummary::{GenericTable, NamedTable, PrettyMetric, TableRow, TitleWithHelp};
use serde::Serialize;
use std::convert::Into;

const HELP_TEXT: &str = "This table displays the 10 most shared CDR3 amino acid sequences. It is ordered by the number of shared donor/origin/libraries and the number of cells. The \"Categorize by:\" button allows you to display the 10 most prevalent CDR3 amino acid sequences by donor (number of donors a CDR3 is shared between), origin (number of origins a CDR3 is shared between), and libraries (number of libraries a CDR3 is shared between). This table also displays the number of cells sharing a given CDR3 amino acid sequence.";
const CARD_TITLE: &str = "Top Shared CDR3s";

#[derive(Debug, Serialize, PartialEq, Eq, Clone)]
#[serde(into = "VdjAggrSharedCdr3View")]
pub struct VdjAggrSharedCdr3 {
    pub by_donor: Vec<VdjAggrSharedCdr3Row>,
    pub by_origin: Vec<VdjAggrSharedCdr3Row>,
    pub by_library: Vec<VdjAggrSharedCdr3Row>,
}

#[derive(Debug, Serialize, PartialEq)]
struct VdjAggrSharedCdr3View {
    tables: Vec<NamedTable>,
    help: TitleWithHelp,
}

impl From<VdjAggrSharedCdr3> for VdjAggrSharedCdr3View {
    fn from(src: VdjAggrSharedCdr3) -> Self {
        VdjAggrSharedCdr3View {
            tables: vec![
                to_named_table(src.by_donor, "Donors"),
                to_named_table(src.by_origin, "Origins"),
                to_named_table(src.by_library, "Libraries"),
            ],
            help: TitleWithHelp {
                help: HELP_TEXT.into(),
                title: CARD_TITLE.into(),
            },
        }
    }
}

/// One row within the table
#[derive(Debug, Serialize, PartialEq, Clone, Ord, PartialOrd, Eq)]
#[serde(into = "TableRow")]
pub struct VdjAggrSharedCdr3Row {
    pub num_categories: usize,
    pub num_cells: usize,
    pub cdr3_aa: String,
}

fn to_named_table(rows: Vec<VdjAggrSharedCdr3Row>, category: &str) -> NamedTable {
    NamedTable {
        name: category.into(),
        table: GenericTable {
            header: Some(vec![
                "CDR3".into(),
                format!("Number of {category}"),
                "Total Number of Cells".into(),
            ]),
            rows: rows.into_iter().map(Into::into).collect(),
            grouping_header: None,
        },
    }
}

impl From<VdjAggrSharedCdr3Row> for TableRow {
    fn from(src: VdjAggrSharedCdr3Row) -> TableRow {
        TableRow(vec![
            src.cdr3_aa,
            PrettyMetric::integer(src.num_categories as i64).0,
            PrettyMetric::integer(src.num_cells as i64).0,
        ])
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::websummary::check_eq_json;

    #[test]
    fn test_row() {
        let row = VdjAggrSharedCdr3Row {
            cdr3_aa: "IGH:CASSFSTCSANYGYTF".into(),
            num_categories: 1,
            num_cells: 5,
        };
        check_eq_json(
            &serde_json::to_string(&row).unwrap(),
            r#"[
                "IGH:CASSFSTCSANYGYTF",
                "1",
                "5"
            ]"#,
        );
    }

    #[test]
    fn test_table() {
        let table = VdjAggrSharedCdr3 {
            by_donor: vec![
                VdjAggrSharedCdr3Row {
                    cdr3_aa: "IGH:CASSFSTCSANYGYTF".into(),
                    num_categories: 2,
                    num_cells: 5,
                },
                VdjAggrSharedCdr3Row {
                    cdr3_aa: "IGH:CASSFSTCSANYGYTF".into(),
                    num_categories: 1,
                    num_cells: 5,
                },
            ],
            by_origin: vec![VdjAggrSharedCdr3Row {
                cdr3_aa: "IGK:CAVRDPDTGRRALTF".into(),
                num_categories: 2,
                num_cells: 300,
            }],
            by_library: vec![],
        };
        check_eq_json(
            &serde_json::to_string(&table).unwrap(),
            r#"{
                "tables": [
                    {
                        "name": "Donors",
                        "table": {
                            "header": [
                                "CDR3",
                                "Number of Donors",
                                "Total Number of Cells"
                            ],
                            "rows": [
                                [
                                    "IGH:CASSFSTCSANYGYTF",
                                    "2",
                                    "5"
                                ],
                                [
                                    "IGH:CASSFSTCSANYGYTF",
                                    "1",
                                    "5"
                                ]
                            ]
                        }
                    },
                    {
                        "name": "Origins",
                        "table": {
                            "header": [
                                "CDR3",
                                "Number of Origins",
                                "Total Number of Cells"
                            ],
                            "rows": [
                                [
                                    "IGK:CAVRDPDTGRRALTF",
                                    "2",
                                    "300"
                                ]
                            ]
                        }
                    },
                    {
                        "name": "Libraries",
                        "table": {
                            "header": [
                                "CDR3",
                                "Number of Libraries",
                                "Total Number of Cells"
                            ],
                            "rows": []
                        }
                    }
                ],
                "help": {
                    "helpText": "This table displays the 10 most shared CDR3 amino acid sequences. It is ordered by the number of shared donor/origin/libraries and the number of cells. The \"Categorize by:\" button allows you to display the 10 most prevalent CDR3 amino acid sequences by donor (number of donors a CDR3 is shared between), origin (number of origins a CDR3 is shared between), and libraries (number of libraries a CDR3 is shared between). This table also displays the number of cells sharing a given CDR3 amino acid sequence.",
                    "title": "Top Shared CDR3s"
                }
            }"#,
        );
    }

    #[test]
    fn test_ord() {
        assert!(
            VdjAggrSharedCdr3Row {
                cdr3_aa: "IGL:CAAWDDSLNGWVF".into(),
                num_categories: 2,
                num_cells: 21
            } > VdjAggrSharedCdr3Row {
                cdr3_aa: "IGK:CQQSYSTPRTF".into(),
                num_categories: 1,
                num_cells: 100
            }
        );
        assert!(
            VdjAggrSharedCdr3Row {
                cdr3_aa: "IGL:CAAWDDSLNGWVF".into(),
                num_categories: 1,
                num_cells: 21
            } > VdjAggrSharedCdr3Row {
                cdr3_aa: "IGK:CQQSYSTPRTF".into(),
                num_categories: 1,
                num_cells: 19
            }
        );
        assert!(
            VdjAggrSharedCdr3Row {
                cdr3_aa: "IGL:CAAWDDSLNGWVF".into(),
                num_categories: 2,
                num_cells: 21
            } > VdjAggrSharedCdr3Row {
                cdr3_aa: "IGK:CQQSYSTPRTF".into(),
                num_categories: 2,
                num_cells: 21
            }
        );
    }
}
