//! Metric implementation for f64 and i64

use crate::{JsonReport, JsonReporter, Metric};

macro_rules! impl_metric_and_json_report {
    ($type:ty) => {
        impl Metric for $type {
            fn merge(&mut self, other: $type) {
                *self += other
            }
        }

        impl JsonReport for $type {
            fn to_json_reporter(&self) -> JsonReporter {
                let mut reporter = JsonReporter::default();
                reporter.insert("", *self);
                reporter
            }
        }
    };
}

// Not implementing for other types
impl_metric_and_json_report! { f64 }
impl_metric_and_json_report! { i64 }
impl_metric_and_json_report! { usize }
