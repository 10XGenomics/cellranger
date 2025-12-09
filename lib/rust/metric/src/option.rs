//!
//! This module implements `Metric` and `JsonReport` for `Option<T>`
//! where T is a `Metric` or `JsonReport`
//!
#![deny(missing_docs)]

use crate::{JsonReport, JsonReporter, Metric};

/// `JsonReport` for `Option<J>` where `J: JsonReport` is straightforward.
/// We will report J if is it Some and not report it if it is None.
impl<J> JsonReport for Option<J>
where
    J: JsonReport,
{
    fn to_json_reporter(&self) -> JsonReporter {
        match self {
            Some(j) => j.to_json_reporter(),
            None => JsonReporter::default(),
        }
    }
}

/// Implementing `Metric` for `Option<M>` where `M: Metric` is more tricky.
/// We have different choices on how new() and merge() would work, whether
/// "zero" is `Some(0)` or `None` and whether it's legal to merge `Some(t)`
/// and `None`. We are choosing a simple implementation here, where we would
/// define None as the "zero" and we will follow simple rules for merging:
/// `Some(m1) + Some(m2) = Some(m1 + m2)` and `Some(m) + None = Some(m)`
impl<M> Metric for Option<M>
where
    M: Metric,
{
    fn merge(&mut self, other: Self) {
        match (self.as_mut(), other) {
            (Some(t1), Some(t2)) => t1.merge(t2),
            (None, Some(t2)) => {
                self.replace(t2);
            }
            _ => {}
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::count_metric::CountMetric;

    #[test]
    fn test_option_metric() {
        let mut count1: Option<CountMetric> = Option::default();
        assert_eq!(count1, None);
        // None + None = None
        count1.merge(None);
        assert_eq!(count1, None);
        // None + Some(m) = Some(m)
        count1.merge(Some(20.into()));
        assert_eq!(count1, Some(CountMetric::from(20)));

        let mut count2 = Some(CountMetric::from(10));
        // Some(m) + None = Some(m)
        count2.merge(None);
        assert_eq!(count2, Some(CountMetric::from(10)));
        // Some(m1) + Some(m2) = Some(m1 + m2)
        count2.merge(Some(40.into()));
        assert_eq!(count2, Some(CountMetric::from(50)));
    }

    #[test]
    fn test_option_report() {
        let count: Option<CountMetric> = None;
        assert_eq!(count.to_json_reporter(), JsonReporter::default());

        let inner = CountMetric::from(100);
        let inner_report = inner.to_json_reporter();
        let count = Some(inner);
        assert_eq!(count.to_json_reporter(), inner_report);
    }
}
