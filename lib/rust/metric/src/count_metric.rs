//!
//! This module defines the `CountMetric` struct. As the name implies, this struct
//! is geared toward tracking counts, typically read counts. Internally it is
//! represented using an `i64` variable.

use crate::{JsonReport, JsonReporter, Metric};

/// Use this struct to keep track of metrics which can be represented using
/// a single `i64` variable
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default, Ord, PartialOrd)]
#[serde(transparent)]
pub struct CountMetric {
    pub(crate) count: i64,
}

impl Metric for CountMetric {
    /// Merging two `CountMetric` objects is just adding up the two counts
    ///
    /// # Example
    /// ```rust
    /// use metric::{Metric, CountMetric};
    /// let mut c1 = CountMetric::default(); // Initialize to zero
    /// c1.increment(); // Now it is 1
    /// let c2 = CountMetric::from(50);
    /// c1.merge(c2);
    /// assert!(c1 == CountMetric::from(51))
    /// ```
    fn merge(&mut self, other: Self) {
        self.count += other.count;
    }
}

impl CountMetric {
    /// Increment the counter by 1
    ///
    /// # Example
    /// ```rust
    /// use metric::{Metric, CountMetric};
    /// let mut c1 = CountMetric::default(); // Initialize to zero
    /// c1.increment(); // Now it is 1
    /// assert!(c1 == CountMetric::from(1))
    /// ```
    pub fn increment(&mut self) {
        self.count += 1;
    }

    /// Increment the counter by the specified value. The value can be anything
    /// which can be cast into `i64`
    ///
    /// # Example
    /// ```rust
    /// use metric::{Metric, CountMetric};
    /// let mut c1 = CountMetric::default(); // Initialize to zero
    /// c1.increment_by(10); // Now it is 10
    /// assert!(c1 == CountMetric::from(10))
    /// ```
    pub fn increment_by<T>(&mut self, val: T)
    where
        T: Into<i64>,
    {
        self.count += val.into();
    }

    /// Return the count
    pub fn count(self) -> i64 {
        self.count
    }
}

/// A convenience function to create a new instance of a `CountMetric`
/// from any type which can be cast into an `i64`
///
/// # Example
///
/// ```rust
/// use metric::{Metric, CountMetric};
/// let c1 = CountMetric::from(100);
/// let mut c2 = CountMetric::default();
/// c2.increment_by(100);
/// assert!(c1==c2);
/// ```
impl<T> From<T> for CountMetric
where
    T: Into<i64>,
{
    fn from(val: T) -> Self {
        CountMetric { count: val.into() }
    }
}

/// Use the key "count" for reporting the count metric
///
impl JsonReport for CountMetric {
    fn to_json_reporter(&self) -> JsonReporter {
        let mut reporter = JsonReporter::default();
        reporter.insert("count", self.count);
        reporter
    }
}

/// Impl `Add` trait to overload the + operator.
/// Maybe useful to have a macro if this pattern needs to be repeated >3 times
///
/// ## Tests
///
/// * prop_test_count_metric_add_and_addassign() - Make sure that the basic operation
///   is correct for assorted inputs
impl ::std::ops::Add for CountMetric {
    type Output = CountMetric;
    fn add(mut self, other: CountMetric) -> CountMetric {
        self.merge(other);
        self
    }
}

/// Impl `AddAssign` trait to overload the += operator.
/// Maybe useful to have a macro if this pattern needs to be repeated >3 times
///
/// ## Tests
///
/// * prop_test_count_metric_add_and_addassign() - Make sure that the basic operation
///   is correct for assorted inputs
impl ::std::ops::AddAssign for CountMetric {
    fn add_assign(&mut self, other: CountMetric) {
        self.merge(other);
    }
}

/// Impl `Sum` trait for an iterator which yields a `CountMetric`
/// Maybe useful to have a macro if this pattern needs to be repeated >3 times
///
/// ## Tests
///
/// * test_merge_count_metric() - Merge a vector of CountMetrics
impl ::std::iter::Sum for CountMetric {
    fn sum<I: Iterator<Item = CountMetric>>(iter: I) -> CountMetric {
        iter.fold(CountMetric::default(), |a, b| a + b)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::arbitrary::any;

    #[test]
    fn test_merge_count_metric() {
        let counts: Vec<_> = (1..10).map(CountMetric::from).collect();
        let merged = counts.iter().fold(CountMetric::default(), |mut acc, &c| {
            acc.merge(c);
            acc
        });
        assert!(merged == CountMetric::from(45));

        let merged: CountMetric = (1..10).map(Into::into).sum();
        assert!(merged == CountMetric::from(45));
    }

    proptest! {
        #[test]
        fn prop_test_count_metric_add_and_addassign(
            x in any::<u32>(),
            y in any::<u32>()
        ) {
            let mut c_x = CountMetric::from(x);
            let c_y = CountMetric::from(y);
            let c_total = CountMetric::from(x as i64 + y as i64);
            assert_eq!(c_total, c_x + c_y);
            c_x += c_y;
            assert_eq!(c_total, c_x);
        }
    }
}
