//!
//! This module defines the `MeanMetric` struct. This struct
//! helps accumulate a mean value over observations.

use crate::{JsonReport, JsonReporter, Metric};
use serde_json::Value;
use std::iter::{FromIterator, Sum};

/// Use this struct to accumulate mean metrics
///
/// # Example
/// ```rust
/// use metric::{Metric, MeanMetric};
/// let mut c1 = MeanMetric::new(); // Initialize to zero
/// c1.record(1.0); // Now it is 1
/// let mut c2 = MeanMetric::new();
/// c2.record(3.0);
/// c1.merge(c2);
/// assert!((c1.mean() - 2.0).abs() < 0.00001)
/// ```
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Default, PartialOrd, Metric)]
pub struct MeanMetric {
    pub(crate) numerator: f64,
    pub(crate) denominator: usize,
}

impl MeanMetric {
    /// Record a new observation
    ///
    /// # Example
    /// ```rust
    /// use metric::{Metric, MeanMetric};
    /// let mut c1 = MeanMetric::new();
    /// c1.record(5.0);
    /// c1.record(7);
    /// assert!((c1.mean() - 6.0).abs() < 0.0001)
    /// ```
    pub fn record<T>(&mut self, value: T)
    where
        T: Into<f64>,
    {
        self.denominator += 1;
        self.numerator += value.into();
    }

    /// Return the mean
    pub fn mean(self) -> f64 {
        self.numerator / (self.denominator as f64)
    }
}

/// Calculate a mean metric from an iterator over integers.
impl<T> FromIterator<T> for MeanMetric
where
    i64: Sum<T>,
{
    /// Calculate a mean metric from an iterator over integers.
    fn from_iter<I: IntoIterator<Item = T>>(iter: I) -> Self {
        let mut denominator = 0;
        let numerator: i64 = iter.into_iter().inspect(|_| denominator += 1).sum();
        MeanMetric {
            // This conversion from i64 to f64 is lossy when it exceeds 2**52 or ~4.5e15.
            numerator: numerator as f64,
            denominator,
        }
    }
}

impl From<&MeanMetric> for Value {
    fn from(metric: &MeanMetric) -> Value {
        let mean = metric.mean();
        if mean.is_nan() {
            Value::from("NaN")
        } else {
            Value::from(mean)
        }
    }
}

impl JsonReport for MeanMetric {
    /// Use the key "mean" for reporting the mean metric.
    /// Use "NaN" if the denominator is zero.
    fn to_json_reporter(&self) -> JsonReporter {
        JsonReporter::from(("mean", self))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::arbitrary::any;

    #[test]
    fn test_merge_mean_metric() {
        let means: Vec<_> = (1..10)
            .map(|f| {
                let mut m = MeanMetric::new();
                m.record(f);
                m
            })
            .collect();
        let merged = means.iter().fold(MeanMetric::new(), |mut acc, &c| {
            acc.merge(c);
            acc
        });
        println!("{}", merged.mean());
        assert!((merged.mean() - 5.0).abs() < 0.000001f64);
    }

    #[test]
    fn test_from_iter() {
        assert_eq!(MeanMetric::from_iter(0..10).mean(), 4.5);
        let xs: Vec<_> = (0..10).collect();
        assert_eq!(MeanMetric::from_iter(&xs).mean(), 4.5);
        assert_eq!(MeanMetric::from_iter(xs).mean(), 4.5);
        let xs: Vec<i64> = vec![];
        assert!(MeanMetric::from_iter(xs).mean().is_nan());
    }

    #[test]
    fn test_collect() {
        assert_eq!((0..10).collect::<MeanMetric>().mean(), 4.5);
        let xs: Vec<_> = (0..10).collect();
        assert_eq!(xs.iter().collect::<MeanMetric>().mean(), 4.5);
        assert_eq!(xs.into_iter().collect::<MeanMetric>().mean(), 4.5);
    }

    proptest! {
        #[test]
        fn prop_test_mean_metric_add_and_addassign(
            x in any::<f64>(),
            y in any::<f64>(),
            z in any::<f64>()
        ) {
            let mean = (x + y + z) / 3.0;
            let tolerance = 0.0001 * [x, y, z].into_iter().map(f64::abs).max_by(f64::total_cmp).unwrap();

            let mut m_x = MeanMetric::new();
            m_x.record(x);
            let mut m_yz = MeanMetric::new();
            m_yz.record(y);
            m_yz.record(z);

            assert!(((m_yz + m_x).mean() - mean).abs() <= tolerance);
            m_x += m_yz;
            assert!((m_x.mean() - mean).abs() <= tolerance);
        }
    }
}
