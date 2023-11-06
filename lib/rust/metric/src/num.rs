//! Metric implementation for f64 and i64

use crate::{JsonReport, JsonReporter, Metric};
use num_traits::Zero;

macro_rules! impl_metric_and_json_report {
    ($type:ty) => {
        impl Metric for $type {
            fn new() -> $type {
                <$type>::zero()
            }
            fn merge(&mut self, other: $type) {
                *self += other
            }
        }
        impl JsonReport for $type {
            fn to_json_reporter(&self) -> JsonReporter {
                let mut reporter = JsonReporter::new();
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
