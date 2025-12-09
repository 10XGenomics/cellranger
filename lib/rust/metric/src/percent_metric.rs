//! This module defines the `PercentMetric` struct. As the name implies, this struct
//! is geared toward tracking percentages. Internally it is represented using a
//! numerator and denominator counter.
#![deny(missing_docs)]

use crate::{CountMetric, JsonReport, JsonReporter, Metric};
use metric_derive::Metric;
use serde::{Deserialize, Serialize};
use serde_json::Value;

/// Use this struct to keep track of metrics which can be represented
/// as a fraction with numerator and denominator of type `i64`
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default, Metric)]
pub struct PercentMetric {
    /// Numerator
    pub numerator: CountMetric,
    /// Denominator
    pub denominator: CountMetric,
}

impl PercentMetric {
    /// Generate a `PercentMetric` from a numerator and a denominator
    ///
    /// # Example
    /// ```rust
    /// use metric::{Metric, PercentMetric};
    /// let c1 = PercentMetric::from_parts(50, 100);
    /// let mut c2 = PercentMetric::default();
    /// c2.increment_by(50, true);
    /// c2.increment_by(50, false);
    /// assert!(c1==c2);
    /// let c3 = PercentMetric::from_parts(150, 200);
    /// c2 += c3; // Can also use c2.merge(c3)
    /// assert!(c2==(200, 300).into());
    /// ```
    pub fn from_parts<T>(num: T, den: T) -> Self
    where
        T: Into<CountMetric>,
    {
        PercentMetric {
            numerator: num.into(),
            denominator: den.into(),
        }
    }
    /// Add a single observation. The denominator is incremented by 1
    /// and the numerator is incremented if the filter is true
    ///
    pub fn increment(&mut self, filter: bool) {
        self.denominator.increment();
        if filter {
            self.numerator.increment();
        }
    }

    /// Increment the denominator by the specified value. Increment the numerator
    /// if the filter is true
    ///
    /// # Example
    /// ```rust
    /// use metric::{Metric, PercentMetric};
    /// let mut c1 = PercentMetric::default(); // Initialize to zero
    /// c1.increment_by(10, true); // Now it is 10
    /// assert!(c1 == PercentMetric::from_parts(10, 10))
    /// ```
    pub fn increment_by<T>(&mut self, val: T, filter: bool)
    where
        T: Copy + Into<i64>,
    {
        self.denominator.increment_by(val);
        if filter {
            self.numerator.increment_by(val);
        }
    }

    /// Return the fraction as an `Option`. It is `None` if the denominator is zero
    ///
    /// # Example
    /// ```rust
    /// use metric::{Metric, PercentMetric};
    /// let mut c = PercentMetric::default(); // Initialize to zero
    /// assert!(c.fraction().is_none());
    /// c.increment_by(10, false);
    /// c.increment_by(10, true);
    /// assert!(c.fraction() == Some(0.5f64));
    /// ```
    pub fn fraction(&self) -> Option<f64> {
        if self.denominator.count == 0 {
            None
        } else {
            Some((self.numerator.count as f64) / (self.denominator.count as f64))
        }
    }
}

/// A convenience function to create a new instance of a `PercentMetric`
/// from a boolean.
///
/// * `True` will be converted to a `PercentMetric` with numerator and denominator
///   equal to one
/// * `False` will be converted to a `PercentMetric` with numerator equal to zero
///   and denominator equal to one
///
/// # Example
///
/// ```rust
/// use metric::{Metric, PercentMetric};
/// let c1: PercentMetric = true.into();
/// let c2 = PercentMetric::from_parts(1, 1);
/// assert!(c1==c2);
/// let mut c1: PercentMetric = false.into();
/// c1.increment(true);
/// let c2 = PercentMetric::from_parts(1, 2);
/// assert!(c1==c2);
/// ```
impl From<bool> for PercentMetric {
    fn from(b: bool) -> Self {
        let num = i32::from(b);
        PercentMetric::from_parts(num, 1)
    }
}

/// A convenience function to create a new instance of a `PercentMetric`
/// from a tuple with two elements of a type which can be converter to a
/// `CountMetric`.
///
/// # Example
///
/// ```rust
/// use metric::{Metric, PercentMetric};
/// let c1: PercentMetric = (5, 37).into();
/// let c2 = PercentMetric::from_parts(5, 37);
/// assert!(c1==c2);
/// ```
impl<T> From<(T, T)> for PercentMetric
where
    T: Into<CountMetric>,
{
    fn from((num, den): (T, T)) -> Self {
        PercentMetric::from_parts(num, den)
    }
}

impl From<&PercentMetric> for Value {
    fn from(metric: &PercentMetric) -> Value {
        if let Some(frac) = metric.fraction() {
            Value::from(frac)
        } else {
            Value::from("NaN")
        }
    }
}

impl JsonReport for PercentMetric {
    /// Use the key "frac" for reporting the percent metric.
    /// Use "NaN" if the denominator is zero.
    fn to_json_reporter(&self) -> JsonReporter {
        JsonReporter::from(("frac", self))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::arbitrary::any;
    use std::cmp::{max, min};

    #[test]
    fn test_merge_percent_metric() {
        let mut p1 = PercentMetric {
            numerator: CountMetric::from(5),
            denominator: CountMetric::from(10),
        };

        let p2 = PercentMetric {
            numerator: CountMetric::from(7),
            denominator: CountMetric::from(20),
        };
        p1.merge(p2);

        let p = PercentMetric {
            numerator: CountMetric::from(12),
            denominator: CountMetric::from(30),
        };

        assert!(p1 == p);
    }

    proptest::proptest! {
        #[test]
        fn prop_test_percent_metric_add_and_addassign(
            x1 in any::<u32>(),
            x2 in any::<u32>(),
            y1 in any::<u32>(),
            y2 in any::<u32>()
        ) {
            let x_num = min(x1, x2);
            let x_den = max(x1, x2);
            let mut p_x: PercentMetric = (x_num, x_den).into();
            let y_num = min(y1, y2);
            let y_den = max(y1, y2);
            let p_y: PercentMetric = (y_num, y_den).into();

            let p_total: PercentMetric = (x_num as i64 + y_num as i64, x_den as i64 + y_den as i64).into();

            assert_eq!(p_total, p_x + p_y);
            p_x += p_y;
            assert_eq!(p_total, p_x);
        }
    }
}
