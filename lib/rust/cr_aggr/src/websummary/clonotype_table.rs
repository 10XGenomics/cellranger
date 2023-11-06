use cr_websummary::*;
use serde::Serialize;

card_with_table!(
    table_name: VdjAggrClonotypeTable,
    row_name: VdjAggrClonotypeRow,
    title: "Top 10 Clonotype CDR3 Sequences",
    help_text: "This table lists the CDR3 sequence of the first exact subclonotype of the 10 most abundant clonotypes in the aggregated datasets. For each of the top 10 clonotypes, the constant region, number of cells (frequency), and what percentage of the dataset those cells occupy (proportion) are also displayed. For the full table and more details, please refer to the clonotypes.csv and consensus_annotations.csv files produced by the pipeline.",
    columns: {
        id: usize : "Clonotype ID",
        cdr3_aas : Vec<String> : "CDR3s",
        num_cells: usize : "Frequency",
        proportion: Percent : "Proportion"
    }
);

#[cfg(test)]
mod tests {
    use super::*;
    use crate::websummary::check_eq_json;

    fn make_percent(percent: f64) -> Percent {
        Percent::Float(percent / 100.0)
    }

    #[test]
    fn test_table_row() {
        let row = VdjAggrClonotypeRow {
            id: 3,
            cdr3_aas: vec![
                "IGH:CAVSDLEPNSSASKIIF".into(),
                "IGK:CASSFSTCSANYGYTF".into(),
            ],
            num_cells: 4,
            proportion: make_percent(4.3971),
        };
        check_eq_json(
            &serde_json::to_string(&row).unwrap(),
            r#"[
                "3",
                "IGH:CAVSDLEPNSSASKIIF<br>IGK:CASSFSTCSANYGYTF",
                "4",
                "4.40%"
            ]"#,
        );
    }

    #[test]
    fn test_full_table() {
        let table = VdjAggrClonotypeTable(vec![
            VdjAggrClonotypeRow {
                id: 1,
                cdr3_aas: vec!["IGH:CASSFSTCSANYGYTF".into()],
                num_cells: 43,
                proportion: make_percent(47.2511),
            },
            VdjAggrClonotypeRow {
                id: 3,
                cdr3_aas: vec![
                    "IGH:CAVSDLEPNSSASKIIF".into(),
                    "IGK:CASSFSTCSANYGYTF".into(),
                ],
                num_cells: 4,
                proportion: make_percent(4.3971),
            },
        ]);
        check_eq_json(
            &serde_json::to_string(&table).unwrap(),
            r#"{
                "table": {
                    "header": [
                        "Clonotype ID",
                        "CDR3s",
                        "Frequency",
                        "Proportion"
                    ],
                    "rows": [
                        [
                            "1",
                            "IGH:CASSFSTCSANYGYTF",
                            "43",
                            "47.25%"
                        ],
                        [
                            "3",
                            "IGH:CAVSDLEPNSSASKIIF<br>IGK:CASSFSTCSANYGYTF",
                            "4",
                            "4.40%"
                        ]
                    ]
                },
                "help": {
                    "helpText": "This table lists the CDR3 sequence of the first exact subclonotype of the 10 most abundant clonotypes in the aggregated datasets. For each of the top 10 clonotypes, the constant region, number of cells (frequency), and what percentage of the dataset those cells occupy (proportion) are also displayed. For the full table and more details, please refer to the clonotypes.csv and consensus_annotations.csv files produced by the pipeline.",
                    "title": "Top 10 Clonotype CDR3 Sequences"
                }
            }"#,
        );
    }
}
