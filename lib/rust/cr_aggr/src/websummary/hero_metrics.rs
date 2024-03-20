use crate::websummary::{HeroMetric, PrettyMetric, Threshold};
use serde::{Deserialize, Serialize};

const NUMBER_OF_CELLS: &str = "Total Number of Cells";
const NUMBER_OF_CLONOTYPES: &str = "Total Number of Clonotypes";
const NUMBER_OF_PAIRED_CELLS: &str = "Number of Cells With Productive V-J Spanning Pair";

#[derive(Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct VdjAggrHeroMetrics {
    total_num_cells: HeroMetric,
    total_num_clonotypes: HeroMetric,
    total_paired_cells: HeroMetric,
}

impl VdjAggrHeroMetrics {
    pub fn new(num_cells: usize, num_clonotypes: usize, num_paired_cells: usize) -> Self {
        VdjAggrHeroMetrics {
            total_num_cells: HeroMetric {
                threshold: Threshold::Pass,
                name: NUMBER_OF_CELLS.into(),
                metric: PrettyMetric::integer(num_cells as i64),
            },
            total_num_clonotypes: HeroMetric {
                threshold: Threshold::Pass,
                name: NUMBER_OF_CLONOTYPES.into(),
                metric: PrettyMetric::integer(num_clonotypes as i64),
            },
            total_paired_cells: HeroMetric {
                threshold: Threshold::Pass,
                name: NUMBER_OF_PAIRED_CELLS.into(),
                metric: PrettyMetric::integer(num_paired_cells as i64),
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::websummary::{check_eq_json, test_json_roundtrip};

    #[test]
    fn test_hero_metrics_parse() {
        let json = r#"{
            "total_num_cells": {
                "threshold": "pass",
                "metric": "18,000",
                "name": "Total Number of Cells"
            },
            "total_num_clonotypes": {
                "threshold": "pass",
                "metric": "9,000",
                "name": "Total Number of Clonotypes"
            },
            "total_paired_cells": {
                "threshold": "pass",
                "metric": "12,000",
                "name": "Number of Cells With Productive V-J Spanning Pair"
            }
        }"#;
        test_json_roundtrip::<VdjAggrHeroMetrics>(json);
        check_eq_json(
            &serde_json::to_string(&VdjAggrHeroMetrics::new(18000, 9000, 12000)).unwrap(),
            json,
        );
    }
}
