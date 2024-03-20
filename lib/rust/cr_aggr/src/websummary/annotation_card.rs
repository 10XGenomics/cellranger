//
use crate::websummary::{
    CardWithMetric, GenericTable, PrettyMetric, TableRow, TermDesc, TitleWithTermDesc,
};
use serde::Serialize;

const PRODUCTIVE_CONTIG: (&str, &str) = (
    "Productive Contig",
    "A productive contig satisfies the following conditions: the contig annotations span \
    the 5' end of the V region to the 3' end of the J region of the chain, a start codon was \
    found in the expected part of the V sequence, an in-frame CDR3 amino acid motif was found, \
    and no stop codons were found in the aligned V-J region.",
);
const IGHK: &str = "(IGK, IGH) "; // NOTE: The terminal space is needed. See the fmt string below!
const IGHL: &str = "(IGL, IGH) "; // NOTE: The terminal space is needed. See the fmt string below!
const TRAB: &str = "(TRA, TRB) "; // NOTE: The terminal space is needed. See the fmt string below!
const TRGD: &str = "(TRG, TRD) "; // NOTE: The terminal space is needed. See the fmt string below!
const ANNOTATION_TABLE_TITLE: &str = "Annotation";

fn paired_cells_name(pair: Option<&str>) -> String {
    format!(
        "Cells With Productive V-J Spanning {}Pair",
        pair.unwrap_or("")
    )
}
fn paired_cells_help(pair: Option<&str>) -> String {
    format!(
        "Fraction of cell-associated barcodes with at least one productive contig for \
        each chain of the {}receptor pair.",
        pair.unwrap_or("")
    )
}

#[derive(Debug, Serialize, PartialEq, Clone)]
#[serde(into = "CardWithMetric")]
pub struct VdjAggrAnnotationTable {
    pub paired_cells_frac: f64,
    pub trab_paired_cells_frac: Option<f64>,
    pub ighk_paired_cells_frac: Option<f64>,
    pub ighl_paired_cells_frac: Option<f64>,
    pub trgd_paired_cells_frac: Option<f64>,
}

impl From<VdjAggrAnnotationTable> for CardWithMetric {
    fn from(src: VdjAggrAnnotationTable) -> CardWithMetric {
        let mut term_descs = vec![
            TermDesc::with_one_desc(PRODUCTIVE_CONTIG.0, PRODUCTIVE_CONTIG.1),
            TermDesc::with_one_desc(paired_cells_name(None), paired_cells_help(None)),
        ];

        let mut rows = vec![TableRow::two_col(
            paired_cells_name(None),
            PrettyMetric::percent(src.paired_cells_frac),
        )];

        if let Some(f) = src.trab_paired_cells_frac {
            rows.push(TableRow::two_col(
                paired_cells_name(Some(TRAB)),
                PrettyMetric::percent(f),
            ));
            term_descs.push(TermDesc::with_one_desc(
                paired_cells_name(Some(TRAB)),
                paired_cells_help(Some(TRAB)),
            ));
        }
        if let Some(f) = src.trgd_paired_cells_frac {
            rows.push(TableRow::two_col(
                paired_cells_name(Some(TRGD)),
                PrettyMetric::percent(f),
            ));
            term_descs.push(TermDesc::with_one_desc(
                paired_cells_name(Some(TRGD)),
                paired_cells_help(Some(TRGD)),
            ));
        }
        if let Some(f) = src.ighk_paired_cells_frac {
            rows.push(TableRow::two_col(
                paired_cells_name(Some(IGHK)),
                PrettyMetric::percent(f),
            ));
            term_descs.push(TermDesc::with_one_desc(
                paired_cells_name(Some(IGHK)),
                paired_cells_help(Some(IGHK)),
            ));
        }
        if let Some(f) = src.ighl_paired_cells_frac {
            rows.push(TableRow::two_col(
                paired_cells_name(Some(IGHL)),
                PrettyMetric::percent(f),
            ));
            term_descs.push(TermDesc::with_one_desc(
                paired_cells_name(Some(IGHL)),
                paired_cells_help(Some(IGHL)),
            ));
        }

        CardWithMetric {
            table: GenericTable { header: None, rows },
            help: TitleWithTermDesc {
                title: ANNOTATION_TABLE_TITLE.into(),
                data: term_descs,
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::websummary::check_eq_json;

    #[test]
    fn test_vdj_annotation_table_tcr() {
        let table = VdjAggrAnnotationTable {
            paired_cells_frac: 2.0f64 / 3.0f64,
            trab_paired_cells_frac: Some(2.0f64 / 3.0f64),
            ighk_paired_cells_frac: None,
            ighl_paired_cells_frac: None,
            trgd_paired_cells_frac: None,
        };
        check_eq_json(
            &serde_json::to_string(&table).unwrap(),
            r#"{
                "table": {
                    "rows": [
                        [
                            "Cells With Productive V-J Spanning Pair",
                            "66.67%"
                        ],
                        [
                            "Cells With Productive V-J Spanning (TRA, TRB) Pair",
                            "66.67%"
                        ]
                    ]
                },
                "help": {
                    "data": [
                        [
                            "Productive Contig",
                            [
                                "A productive contig satisfies the following conditions: the contig annotations span the 5' end of the V region to the 3' end of the J region of the chain, a start codon was found in the expected part of the V sequence, an in-frame CDR3 amino acid motif was found, and no stop codons were found in the aligned V-J region."
                            ]
                        ],
                        [
                            "Cells With Productive V-J Spanning Pair",
                            [
                                "Fraction of cell-associated barcodes with at least one productive contig for each chain of the receptor pair."
                            ]
                        ],
                        [
                            "Cells With Productive V-J Spanning (TRA, TRB) Pair",
                            [
                                "Fraction of cell-associated barcodes with at least one productive contig for each chain of the (TRA, TRB) receptor pair."
                            ]
                        ]
                    ],
                    "title": "Annotation"
                }
            }"#,
        );
    }

    #[test]
    fn test_vdj_annotation_table_ig() {
        let table = VdjAggrAnnotationTable {
            paired_cells_frac: 2.0f64 / 3.0f64,
            trab_paired_cells_frac: None,
            ighk_paired_cells_frac: Some(1.0f64 / 2.0f64),
            ighl_paired_cells_frac: Some(1.0f64 / 6.0f64),
            trgd_paired_cells_frac: None,
        };
        check_eq_json(
            &serde_json::to_string(&table).unwrap(),
            r#"{
                "table": {
                    "rows": [
                        [
                            "Cells With Productive V-J Spanning Pair",
                            "66.67%"
                        ],
                        [
                            "Cells With Productive V-J Spanning (IGK, IGH) Pair",
                            "50.00%"
                        ],
                        [
                            "Cells With Productive V-J Spanning (IGL, IGH) Pair",
                            "16.67%"
                        ]
                    ]
                },
                "help": {
                    "data": [
                        [
                            "Productive Contig",
                            [
                                "A productive contig satisfies the following conditions: the contig annotations span the 5' end of the V region to the 3' end of the J region of the chain, a start codon was found in the expected part of the V sequence, an in-frame CDR3 amino acid motif was found, and no stop codons were found in the aligned V-J region."
                            ]
                        ],
                        [
                            "Cells With Productive V-J Spanning Pair",
                            [
                                "Fraction of cell-associated barcodes with at least one productive contig for each chain of the receptor pair."
                            ]
                        ],
                        [
                            "Cells With Productive V-J Spanning (IGK, IGH) Pair",
                            [
                                "Fraction of cell-associated barcodes with at least one productive contig for each chain of the (IGK, IGH) receptor pair."
                            ]
                        ],
                        [
                            "Cells With Productive V-J Spanning (IGL, IGH) Pair",
                            [
                                "Fraction of cell-associated barcodes with at least one productive contig for each chain of the (IGL, IGH) receptor pair."
                            ]
                        ]
                    ],
                    "title": "Annotation"
                }
            }"#,
        );
    }
}
