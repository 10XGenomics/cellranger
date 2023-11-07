//!
//! This module implements the `Metric` and `JsonReport` trait
//! for common collections (`Vec`, `HashMap`)
//!

use crate::{JsonReport, JsonReporter, Metric};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap};
use std::hash::{BuildHasher, Hash};

/// `Metric` trait for `BTreeMap` where the value implements the `Metric` trait
impl<K, V> Metric for BTreeMap<K, V>
where
    K: Ord + Eq + Hash + Serialize + for<'de> Deserialize<'de>,
    V: Metric,
{
    fn new() -> Self {
        BTreeMap::default()
    }

    fn merge(&mut self, other: Self) {
        use std::collections::btree_map::Entry::{Occupied, Vacant};
        for (key, value) in other {
            match self.entry(key) {
                Vacant(e) => {
                    e.insert(value);
                }
                Occupied(ref mut e) => {
                    e.get_mut().merge(value);
                }
            }
        }
    }
}

/// `Metric` trait for `HashMap` where the value implements the `Metric` trait
impl<K, V, S> Metric for HashMap<K, V, S>
where
    K: Eq + Hash + Serialize + for<'de> Deserialize<'de>,
    V: Metric,
    S: BuildHasher + Default,
{
    /// Returns the default HashMap
    fn new() -> Self {
        HashMap::default()
    }

    /// Merging two hashmaps correspond to merging the values for keys
    /// which are common to both hashmaps and copying into `self` any new (key, value)
    /// pair from the `other` hashmap
    fn merge(&mut self, other: Self) {
        use std::collections::hash_map::Entry::{Occupied, Vacant};
        for (key, value) in other {
            match self.entry(key) {
                Vacant(e) => {
                    e.insert(value);
                }
                Occupied(ref mut e) => {
                    e.get_mut().merge(value);
                }
            }
        }
    }
}

/// `JsonReport` trait for `HashMap` where the value implements the `JsonReport` trait
/// and the key can be converted into a string
impl<K: ToString, V: JsonReport, S> JsonReport for HashMap<K, V, S> {
    fn to_json_reporter(&self) -> JsonReporter {
        let mut reporter = JsonReporter::new();
        for (key, value) in self {
            let mut this_report = value.to_json_reporter();
            this_report.add_prefix(key);
            reporter.merge(this_report);
        }
        reporter
    }
}

/// `Metric` trait for `Vec<T>` where T implements the `Metric` trait.
impl<T> Metric for Vec<T>
where
    T: Metric,
{
    fn new() -> Self {
        Vec::new()
    }

    /// `self` needs to be either empty or the same size as `other`.
    /// Elements at the same positions are merged
    fn merge(&mut self, mut other: Self) {
        if self.is_empty() {
            self.append(&mut other);
            return;
        }
        assert!(self.len() == other.len());
        for (v1, v2) in self.iter_mut().zip(other.into_iter()) {
            v1.merge(v2);
        }
    }
}

/// `JsonReport` trait for `Vec<T>` where T implements the `JsonReport` trait.
/// Create a reporter from each element and merge them. This means each element
/// in the vector must have unique key(s) when converted to a reporter.
impl<T: JsonReport> JsonReport for Vec<T> {
    fn to_json_reporter(&self) -> JsonReporter {
        let mut reporter = JsonReporter::new();
        for value in self {
            let this_report = value.to_json_reporter();
            reporter.merge(this_report);
        }
        reporter
    }
}

/// Implement `JsonReport` trait for a tuple.
impl<K: ToString, V: Serialize> JsonReport for (K, V) {
    fn to_json_reporter(&self) -> JsonReporter {
        let mut reporter = JsonReporter::new();
        reporter.insert(self.0.to_string(), serde_json::to_value(&self.1).unwrap());
        reporter
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::count_metric::CountMetric;
    use crate::TxHashMap;
    use std::collections::BTreeMap;

    #[test]
    fn check_btree_merge() {
        let mut m1: BTreeMap<String, CountMetric> = BTreeMap::new();
        m1.insert("Hello".into(), CountMetric::from(10));
        m1.insert("World".into(), CountMetric::from(5));

        let mut m2: BTreeMap<String, CountMetric> = BTreeMap::new();
        m2.insert("Hello".into(), CountMetric::from(20));
        m2.insert("Blah".into(), CountMetric::from(40));
        m1.merge(m2);

        let mut m: BTreeMap<String, CountMetric> = BTreeMap::new();
        m.insert("Hello".into(), CountMetric::from(30));
        m.insert("World".into(), CountMetric::from(5));
        m.insert("Blah".into(), CountMetric::from(40));

        assert_eq!(m, m1);
    }

    #[test]
    fn check_hashmap_merge() {
        let mut m1: TxHashMap<String, CountMetric> = TxHashMap::new();
        m1.insert("Hello".into(), CountMetric::from(10));
        m1.insert("World".into(), CountMetric::from(5));

        let mut m2: TxHashMap<String, CountMetric> = TxHashMap::new();
        m2.insert("Hello".into(), CountMetric::from(20));
        m2.insert("Blah".into(), CountMetric::from(40));
        m1.merge(m2);

        let mut m: TxHashMap<String, CountMetric> = TxHashMap::new();
        m.insert("Hello".into(), CountMetric::from(30));
        m.insert("World".into(), CountMetric::from(5));
        m.insert("Blah".into(), CountMetric::from(40));

        assert_eq!(m, m1);
    }

    #[test]
    fn check_hashmap_report() {
        let mut m1: TxHashMap<String, CountMetric> = TxHashMap::new();
        m1.insert("Hello".into(), CountMetric::from(10));
        m1.insert("World".into(), CountMetric::from(5));

        let report = m1.to_json_reporter();

        let mut expected_report = JsonReporter::new();
        expected_report.insert("Hello_count", 10);
        expected_report.insert("World_count", 5);

        assert_eq!(report, expected_report);
    }

    #[test]
    fn check_vec_merge() {
        let mut m1 = vec![CountMetric::from(10), CountMetric::from(5)];

        let m2 = vec![CountMetric::from(20), CountMetric::from(40)];
        m1.merge(m2);

        let m = vec![CountMetric::from(30), CountMetric::from(45)];

        assert_eq!(m, m1);
    }

    #[test]
    #[should_panic]
    fn check_invalid_vec_merge() {
        let mut m1 = vec![CountMetric::from(10)];

        let m2 = vec![CountMetric::from(20), CountMetric::from(40)];
        m1.merge(m2);
    }

    #[test]
    fn check_empty_vec_merge() {
        let mut m1 = Vec::new();

        let m2 = vec![CountMetric::from(20), CountMetric::from(40)];

        let m = m2.clone();

        m1.merge(m2);

        assert_eq!(m, m1);
    }
}
