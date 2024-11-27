#![deny(
    missing_docs,
    missing_copy_implementations,
    non_upper_case_globals,
    trivial_casts,
    trivial_numeric_casts,
    unsafe_code,
    unstable_features,
    unused_extern_crates,
    unused_import_braces,
    unused_qualifications
)]

//!
//! This is documentation for the `metric` crate.
//!
//! This crate defines a number of data structures used for tracking
//! various metrics in a pipeline. `Metric` trait forms the core of this
//! crate and all the data structure which represents a metric implements
//! this trait. The crate `metric_derive` allows one to derive the `Metric`
//! trait. The data structures and their intended uses are noted below:
//!
//!

#[macro_use]
extern crate serde_derive;

#[macro_use]
extern crate metric_derive;

#[cfg(test)]
#[macro_use]
extern crate proptest;

use ahash::AHasher;
use anyhow::{Context, Error};
use serde::ser::SerializeMap;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::borrow::Borrow;
use std::collections::{hash_map, HashMap, HashSet};
use std::fs::File;
use std::hash::{BuildHasher, Hash};
use std::io::{Read, Write};
use std::iter::FromIterator;
use std::path::Path;

pub mod count_metric;
pub use crate::count_metric::CountMetric;
pub mod percent_metric;
pub use crate::percent_metric::PercentMetric;
pub mod histogram;
pub use crate::histogram::SimpleHistogram;
pub mod mean_metric;
pub use crate::mean_metric::MeanMetric;
pub mod collections;
pub mod num;
pub mod option;

/// A deterministic and fast hasher.
#[derive(Clone, Copy, Default)]
pub struct TxHasher;

impl TxHasher {
    fn random_state() -> ahash::RandomState {
        ahash::RandomState::with_seeds(0, 0, 0, 0)
    }

    /// Return a new hasher.
    pub fn hasher() -> AHasher {
        Self::random_state().build_hasher()
    }

    /// Calculate the hash of a single value.
    pub fn hash(x: impl Hash) -> u64 {
        Self::random_state().hash_one(x)
    }
}

impl BuildHasher for TxHasher {
    type Hasher = AHasher;

    fn build_hasher(&self) -> Self::Hasher {
        Self::hasher()
    }
}

/// A default HashMap using some faster hashing scheme
pub type TxHashMap<K, V> = HashMap<K, V, TxHasher>;

/// A default HashSet using some faster hashing scheme
pub type TxHashSet<K> = HashSet<K, TxHasher>;

/// Create a TxHashSet with the given elements. Similar to a `vec![]`
#[macro_export]
macro_rules! set {
    ( $( $x:expr ),* ) => {
        {
            let mut result = TxHashSet::default();
            $(result.insert($x);)*
            result
        }
    };
    ( $( $x:expr, )* ) => {set!($( $x:expr ),*)};
}

/// The string used in a metric name.
pub trait AsMetricPrefix {
    /// Return the metric prefix as a borrowed string, or None if it has none.
    fn as_metric_prefix(&self) -> Option<&str>;
}

impl AsMetricPrefix for Option<&str> {
    fn as_metric_prefix(&self) -> Option<&str> {
        *self
    }
}

impl AsMetricPrefix for String {
    fn as_metric_prefix(&self) -> Option<&str> {
        Some(self)
    }
}

/// The string used in a metric name.
pub trait ToMetricPrefix {
    /// Return the metric prefix as an owned string, or None if it has none.
    fn to_metric_prefix(&self) -> Option<String>;
}

impl<T: AsMetricPrefix> ToMetricPrefix for T {
    fn to_metric_prefix(&self) -> Option<String> {
        self.as_metric_prefix().map(String::from)
    }
}

impl ToMetricPrefix for u16 {
    fn to_metric_prefix(&self) -> Option<String> {
        Some(self.to_string())
    }
}

/// Join the metric prefix and the metric name.
pub fn join_metric_name(prefix: impl AsMetricPrefix, name: &str) -> String {
    let Some(prefix) = prefix.as_metric_prefix() else {
        return name.to_string();
    };
    if name.is_empty() {
        prefix.to_string()
    } else {
        format!("{prefix}_{name}")
    }
}

/// This enum defines the format for serializing/deserializing a metric
/// Supported formats are
/// * `Json` using the `serde_json` crate
/// * `Binary` using the `bincode` crate
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
pub enum SerdeFormat {
    /// Use serde_json to write the output in json format
    Json,
    /// Use `bincode` to write the output in binary format
    Binary,
}

/// The core trait in this crate. Any data structure representing a metric
/// will implement this trait. It may be auto derived using the macros in
/// `metric_derive`
pub trait Metric: Serialize + for<'de> Deserialize<'de> {
    /// Combined two Metric objects, modifying self in place,
    /// consuming the `other`
    fn merge(&mut self, other: Self);

    /// Write(Serialize) the `Metric` to a file using serde
    ///
    /// # Arguments
    ///
    /// * `filename` - Output will be written to this file. It may be of any type
    ///   which can be interpreted as a `Path`, typically a `String`
    /// * `serde_format` - Defines the format for serialization, typically
    ///   `SerdeFormat::Json`
    ///
    /// # Example
    ///
    /// ```rust
    /// use metric::{Metric, CountMetric};
    /// use metric::SerdeFormat;
    /// let file = "example1.json";
    /// let c1 = CountMetric::from(10); // Initialize to 10
    /// c1.to_file(file, SerdeFormat::Json);
    /// let c2 = CountMetric::from_file(file, SerdeFormat::Json);
    /// assert!(c1==c2);
    /// std::fs::remove_file(file);
    /// ```
    fn to_file<P>(&self, filename: P, serde_format: SerdeFormat) -> Result<(), Error>
    where
        P: AsRef<Path>,
    {
        let writer = File::create(filename.as_ref()).with_context(|| {
            format!(
                "Could not open file '{:?}' for writing in Metric::to_file()",
                filename.as_ref()
            )
        })?;
        self.to_writer(writer, serde_format)
    }

    /// Write(Serialize) the `Metric` to an implementor of `Write`
    ///
    /// # Arguments
    ///
    /// * `writer` - An object which implements `Write` trait. `File`, `Stdout`, `Stderr` etc
    /// * `serde_format` - Defines the format for serialization, typically `SerdeFormat::Json`
    fn to_writer<W>(&self, writer: W, serde_format: SerdeFormat) -> Result<(), Error>
    where
        W: Write,
    {
        match serde_format {
            SerdeFormat::Json => serde_json::to_writer_pretty(writer, &self)
                .with_context(|| "Could not write JSON in Metric::to_writer()")?,
            SerdeFormat::Binary => bincode::serialize_into(writer, &self)
                .with_context(|| "Could not write Binary in Metric::to_writer()")?,
        }

        Ok(())
    }

    /// Read(Deserialize) the `Metric` from a file using serde
    ///
    /// # Arguments
    ///
    /// * `filename` - Read from this file. It may be of any type
    ///   which can be interpreted as a `Path`, typically a `String`
    /// * `serde_format` - Defines the format for deserialization, typically
    ///   `SerdeFormat::Json`
    ///
    /// # Example
    ///
    /// ```rust
    /// use metric::{Metric, CountMetric};
    /// use metric::SerdeFormat;
    /// let file = "example2.json";
    /// let c1 = CountMetric::from(10);
    /// c1.to_file(file, SerdeFormat::Json);
    /// let c2 = CountMetric::from_file(file, SerdeFormat::Json);
    /// assert!(c1==c2);
    /// std::fs::remove_file(file);
    /// ```
    fn from_file<P>(filename: P, serde_format: SerdeFormat) -> Self
    where
        P: AsRef<Path>,
    {
        let reader =
            File::open(filename).expect("Could not open file for reading in Metric::from_file()");
        Metric::from_reader(reader, serde_format)
    }

    /// Read(Deserialize) the `Metric` from a set of files using serde
    ///
    /// # Arguments
    ///
    /// * `filenames` - Read from this set of file. It may be a slice/vector of any type
    ///   which can be interpreted as a `Path`, typically a `String`
    /// * `serde_format` - Defines the format for deserialization, typically
    ///   `SerdeFormat::Json` or `SerdeFormat::Binary`
    ///
    /// # Example
    /// ```rust
    /// use metric::{Metric, CountMetric};
    /// use metric::SerdeFormat;
    /// let files = ["example3.json", "example4.json"];
    /// let c0 = CountMetric::from(10);
    /// c0.to_file(files[0], SerdeFormat::Json);
    /// let c1 = CountMetric::from(6);
    /// c1.to_file(files[1], SerdeFormat::Json);
    ///
    /// let c2 = CountMetric::from_files(&files, SerdeFormat::Json);
    /// assert!(c2 == 16.into());
    /// for file in files {
    ///     std::fs::remove_file(file);
    /// }
    /// ```
    fn from_files<I, P>(filenames: I, serde_format: SerdeFormat) -> Self
    where
        I: IntoIterator<Item = P>,
        P: AsRef<Path>,
        Self: Default,
    {
        filenames
            .into_iter()
            .map(|f| Metric::from_file(f, serde_format))
            .fold(Default::default(), |mut merged, this| {
                merged.merge(this);
                merged
            })
    }

    /// Read(Deserialize) the `Metric` from an implementor of `Read`
    ///
    /// # Arguments
    ///
    /// * `reader` - An object which implements `Read` trait. `File`, `Stdin` etc
    /// * `serde_format` - Defines the format for serialization, typically `SerdeFormat::Json`
    fn from_reader<R>(reader: R, serde_format: SerdeFormat) -> Self
    where
        R: Read,
    {
        match serde_format {
            SerdeFormat::Json => serde_json::from_reader(reader)
                .expect("Could not deserialize JSON in Metric::from_reader()"),
            SerdeFormat::Binary => bincode::deserialize_from(reader)
                .expect("Could not deserialize Binary in Metric::from_reader()"),
        }
    }

    /// Merge metrics from an iterator
    fn from_chunks<I>(chunks: I) -> Self
    where
        I: IntoIterator<Item = Self>,
        Self: Default,
    {
        chunks
            .into_iter()
            .fold(Default::default(), |mut merged, this| {
                merged.merge(this);
                merged
            })
    }
}

/// `JsonReporter` is just a hashmap from `String` to `serde_json::Value`.
/// It represents a json summary file. Concretely, a serialized json output
/// of this data structure would correspond to the json report file. All
/// metrics which can be reported will implement the `JsonReport` trait, which
/// converts the Metric into a `JsonReporter` object.
///
#[derive(Debug, Clone, Deserialize, Default, PartialEq, Eq)]
#[serde(transparent)]
pub struct JsonReporter {
    hashmap: TxHashMap<String, Value>,
}

impl Serialize for JsonReporter {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        let mut xs: Vec<_> = self.iter().collect();
        xs.sort_by_key(|&(k, _v)| k);
        let mut map = serializer.serialize_map(Some(xs.len()))?;
        for (k, v) in xs {
            map.serialize_entry(k, v)?;
        }
        map.end()
    }
}

/// `JsonReport` trait defines how a structure would be written out to a
/// metric summary json file. It has a single function `to_json_reporter` which
/// creats a `JsonReporter` object from the data structure
pub trait JsonReport {
    /// Convert `self` into a `JsonReporter` object
    fn to_json_reporter(&self) -> JsonReporter;

    /// Generate a report by converting `self` into a `JsonReporter` and
    /// serializing it
    ///
    /// # Arguments
    /// * `filename` - Output will be written to this file. It may be of any type
    ///   which can be interpreted as a `Path`, typically a `String`
    fn report(&self, filename: &dyn AsRef<Path>) -> Result<(), Error> {
        let reporter = self.to_json_reporter();
        reporter.to_file(filename, SerdeFormat::Json)
    }
}

impl JsonReport for JsonReporter {
    fn to_json_reporter(&self) -> JsonReporter {
        self.clone()
    }
}

impl<K, V> From<(K, V)> for JsonReporter
where
    K: Eq + Hash + ToString,
    Value: From<V>,
{
    /// Convert a (key, value) tuple into a JsonReporter.
    fn from((key, value): (K, V)) -> JsonReporter {
        let mut reporter = JsonReporter::default();
        reporter.insert(key, Value::from(value));
        reporter
    }
}

impl<K: ToString, V: Serialize> Extend<(K, V)> for JsonReporter {
    fn extend<T: IntoIterator<Item = (K, V)>>(&mut self, iter: T) {
        self.hashmap.extend(
            iter.into_iter()
                .map(|(k, v)| (k.to_string(), serde_json::to_value(v).unwrap())),
        );
    }
}

impl<K: ToString, V: Serialize> FromIterator<(K, V)> for JsonReporter {
    fn from_iter<T: IntoIterator<Item = (K, V)>>(iter: T) -> JsonReporter {
        let mut reporter = JsonReporter::default();
        reporter.extend(iter);
        reporter
    }
}

impl<'a> IntoIterator for &'a JsonReporter {
    type Item = (&'a String, &'a Value);
    type IntoIter = hash_map::Iter<'a, String, Value>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

impl<'a> IntoIterator for &'a mut JsonReporter {
    type Item = (&'a String, &'a mut Value);
    type IntoIter = hash_map::IterMut<'a, String, Value>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter_mut()
    }
}

impl IntoIterator for JsonReporter {
    type Item = (String, Value);
    type IntoIter = hash_map::IntoIter<String, Value>;

    fn into_iter(self) -> Self::IntoIter {
        self.hashmap.into_iter()
    }
}

impl JsonReporter {
    /// Insert a new (key, value) pair into the JsonReporter. The key could be
    /// any type which can be casted into a `String` and value can be any type
    /// which can be casted into `serde_json::Value`
    pub fn insert(&mut self, key: impl ToString, value: impl Into<Value>) {
        let key_str = key.to_string();
        assert!(!self.hashmap.contains_key(&key_str));
        self.hashmap.insert(key_str, value.into());
    }

    /// Insert a new (key, value) pair or update the value if the key already
    /// exists in the `JsonReporter`. The key could be any type which can be
    /// casted into a `String` and value can be any type which can be casted
    /// into `serde_json::Value`
    pub fn insert_or_update(&mut self, key: impl ToString, value: impl Into<Value>) {
        let key_str = key.to_string();
        self.hashmap.insert(key_str, value.into());
    }

    /// Removes any (key, value) pair from the `JsonReporter`. The value is
    /// returned as an `Option<Value>` if the key existed in the `JsonReporter`.
    pub fn remove<Q>(&mut self, key: &Q) -> Option<Value>
    where
        String: Borrow<Q>,
        Q: Hash + Eq + ?Sized,
    {
        self.hashmap.remove(key)
    }

    /// Add `prefix` to the keys in the `JsonReporter`.
    /// Replace each "key" with "prefix_key", separated by an underscore.
    pub fn add_prefix(self, prefix: &str) -> Self {
        self.into_iter()
            .map(|(k, v)| {
                (
                    if k.is_empty() {
                        prefix.to_string()
                    } else {
                        format!("{prefix}_{k}")
                    },
                    v,
                )
            })
            .collect()
    }

    /// Get a reference to the value corresponding to the key from the reporter
    pub fn get<Q>(&self, key: &Q) -> Option<&Value>
    where
        String: Borrow<Q>,
        Q: Hash + Eq + ?Sized,
    {
        self.hashmap.get(key)
    }

    /// Get a mutable reference to the value corresponding to the key from the reporter
    pub fn get_mut<Q>(&mut self, key: &Q) -> Option<&mut Value>
    where
        String: Borrow<Q>,
        Q: Hash + Eq + ?Sized,
    {
        self.hashmap.get_mut(key)
    }

    /// Iterates over the `JsonReporter`, draining it along the way
    pub fn drain(&mut self) -> hash_map::Drain<'_, String, Value> {
        self.hashmap.drain()
    }

    /// Iterates over pairs of references to the keys and values of the `JsonReporter`
    pub fn iter(&self) -> hash_map::Iter<'_, String, Value> {
        self.hashmap.iter()
    }

    /// Iterates over pairs of references to the keys and (mutable) values of the `JsonReporter`
    pub fn iter_mut(&mut self) -> hash_map::IterMut<'_, String, Value> {
        self.hashmap.iter_mut()
    }
}

impl Metric for JsonReporter {
    /// Merges two `JsonReporter` objects, consuming the `other` reporter.
    ///
    /// # Panics
    /// Both reporters need to be disjoint. If there is any key shared
    /// between the two reporters, this function panics with an assertion failed.
    fn merge(&mut self, other: Self) {
        // The `merge` operation for the underlying `HashMap` is different
        // from the `Metric::merge()` implemented for `HashMap`
        for (key, val) in other {
            assert!(
                !self.hashmap.contains_key(&key),
                "Common key ({key}) found when merging reporters"
            );
            self.insert(key, val);
        }
    }

    /// Write(Serialize) the underlying `HashMap` to an implementor of `Write`
    ///
    /// # Arguments
    ///
    /// * `writer` - An object which implements `Write` trait. `File`, `Stdout`, `Stderr` etc
    /// * `serde_format` - Defines the format for serialization, typically `SerdeFormat::Json`
    fn to_writer<W>(&self, writer: W, serde_format: SerdeFormat) -> Result<(), Error>
    where
        W: Write,
    {
        match serde_format {
            SerdeFormat::Json => serde_json::to_writer_pretty(writer, &self)
                .context("Could not write JSON in JsonReporter::to_writer()")?,
            SerdeFormat::Binary => bincode::serialize_into(writer, &self)
                .context("Could not write Binary in Metric::to_writer()")?,
        }

        Ok(())
    }

    /// Read(Deserialize) the underlying `Hashmap` from an implementor of `Read` and
    /// create an instance of `JsonReporter`
    ///
    /// # Arguments
    ///
    /// * `reader` - An object which implements `Read` trait. `File`, `Stdin` etc
    /// * `serde_format` - Defines the format for serialization, typically `SerdeFormat::Json`
    fn from_reader<R>(reader: R, serde_format: SerdeFormat) -> Self
    where
        R: Read,
    {
        let hashmap = match serde_format {
            SerdeFormat::Json => serde_json::from_reader(reader)
                .expect("Could not deserialize JSON in Metric::from_reader()"),
            SerdeFormat::Binary => bincode::deserialize_from(reader)
                .expect("Could not deserialize Binary in Metric::from_reader()"),
        };
        JsonReporter { hashmap }
    }
}

/// Impl Add trait to overload the + operator.
/// Maybe useful to have a macro if this pattern needs to be repeated >3 times
///
/// ## Tests
///
/// * test_add_simple() - Adds two simple JsonReporters and make sure that the outputs match
impl ::std::ops::Add for JsonReporter {
    type Output = JsonReporter;
    fn add(mut self, other: JsonReporter) -> JsonReporter {
        self.merge(other);
        self
    }
}

/// Impl AddAssign trait to overload the += operator.
/// Maybe useful to have a macro if this pattern needs to be repeated >3 times
///
/// ## Tests
///
/// * test_addassign_simple() - Adds two simple JsonReporters and make sure that the outputs match
impl ::std::ops::AddAssign for JsonReporter {
    fn add_assign(&mut self, other: JsonReporter) {
        self.merge(other);
    }
}

/// Impl `Sum` trait for an iterator which yields a `JsonReporter`
/// Maybe useful to have a macro if this pattern needs to be repeated >3 times
///
/// ## Tests
///
/// * test_json_reporter_sum() - Merge a vector of JsonReporters
impl ::std::iter::Sum for JsonReporter {
    fn sum<I: Iterator<Item = JsonReporter>>(iter: I) -> JsonReporter {
        iter.fold(JsonReporter::default(), |a, b| a + b)
    }
}

impl ::std::fmt::Display for JsonReporter {
    fn fmt(&self, f: &mut ::std::fmt::Formatter<'_>) -> ::std::fmt::Result {
        write!(f, "{}", serde_json::to_string_pretty(&self).unwrap())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::hash::Hasher;

    #[test]
    fn test_txhasher() {
        assert_eq!(TxHasher::hash(1234567890u64), 11377898369596099812);
        assert_eq!(TxHasher::hash(b"Hello, world!"), 909132802233200014);
        assert_eq!(TxHasher::hash("Hello, world!"), 1677595690243835615);

        let mut hasher = TxHasher::hasher();
        hasher.write(b"Hello, world!");
        assert_eq!(hasher.finish(), 17599674632180055185);
    }

    #[derive(Default, Serialize, Deserialize, Metric)]
    struct NamedMetric {
        named: CountMetric,
    }

    impl NamedMetric {
        pub fn increment(&mut self) {
            self.named.increment();
        }
        pub fn increment_by(&mut self, v: i64) {
            self.named.increment_by(v);
        }
        pub fn count(&self) -> i64 {
            self.named.count()
        }
    }

    #[test]
    fn test_derive_named_struct() {
        let mut m1 = NamedMetric::default();
        let mut m2 = NamedMetric::default();
        assert_eq!(m1.count(), 0);
        assert_eq!(m2.count(), 0);
        m1.increment();
        m2.increment_by(2);
        assert_eq!(m1.count(), 1);
        assert_eq!(m2.count(), 2);
        m1.merge(m2);
        assert_eq!(m1.count(), 3);
        let m3 = NamedMetric { named: 20.into() };
        m1 += m3;
        assert_eq!(m1.count(), 23);
        m1 += NamedMetric { named: 10.into() };
        assert_eq!(m1.count(), 33);
    }

    #[derive(Default, Serialize, Deserialize, Metric)]
    struct TupleMetric(CountMetric);

    impl TupleMetric {
        pub fn increment(&mut self) {
            self.0.increment();
        }
        pub fn increment_by(&mut self, v: i64) {
            self.0.increment_by(v);
        }
        pub fn count(&self) -> i64 {
            self.0.count()
        }
    }

    #[test]
    fn test_derive_tuple_struct() {
        let mut m1 = TupleMetric::default();
        let mut m2 = TupleMetric::default();
        assert_eq!(m1.count(), 0);
        assert_eq!(m2.count(), 0);
        m1.increment();
        m2.increment_by(2);
        assert_eq!(m1.count(), 1);
        assert_eq!(m2.count(), 2);
        m1.merge(m2);
        assert_eq!(m1.count(), 3);
        let m3 = TupleMetric(20.into());
        m1 += m3;
        assert_eq!(m1.count(), 23);
        m1 += TupleMetric(10.into());
        assert_eq!(m1.count(), 33);
    }

    #[test]
    fn test_add_prefix() {
        let mut report = JsonReporter::default();
        report.insert("Hello", 10);
        report.insert("World", "Foo");

        let mut expected_report = JsonReporter::default();
        expected_report.insert("Bar_Hello", 10);
        expected_report.insert("Bar_World", "Foo");

        assert_eq!(report.add_prefix("Bar"), expected_report);
    }

    #[test]
    fn test_add_prefix_blank_key() {
        let mut report = JsonReporter::default();
        report.insert("", 10);
        report.insert("World", "Foo");

        let mut expected_report = JsonReporter::default();
        expected_report.insert("Bar", 10); // No trailing `_`
        expected_report.insert("Bar_World", "Foo");

        assert_eq!(report.add_prefix("Bar"), expected_report);
    }

    #[test]
    fn test_combine_simple() {
        let mut report1 = JsonReporter::default();
        report1.insert("Hello", 10);
        report1.insert("World", "Foo");

        let mut report2 = JsonReporter::default();
        report2.insert("Bar", 20);

        report1.merge(report2);

        let mut expected_report = JsonReporter::default();
        expected_report.insert("Hello", 10);
        expected_report.insert("World", "Foo");
        expected_report.insert("Bar", 20);

        assert_eq!(report1, expected_report);
    }

    #[test]
    fn test_add_simple() {
        let mut report1 = JsonReporter::default();
        report1.insert("Hello", 10);
        report1.insert("World", "Foo");

        let mut report2 = JsonReporter::default();
        report2.insert("Bar", 20);

        let report = report1 + report2;

        let mut expected_report = JsonReporter::default();
        expected_report.insert("Hello", 10);
        expected_report.insert("World", "Foo");
        expected_report.insert("Bar", 20);

        assert_eq!(report, expected_report);
    }

    #[test]
    fn test_addassign_simple() {
        let mut report1 = JsonReporter::default();
        report1.insert("Hello", 10);
        report1.insert("World", "Foo");

        let mut report2 = JsonReporter::default();
        report2.insert("Bar", 20);

        report1 += report2;

        let mut expected_report = JsonReporter::default();
        expected_report.insert("Hello", 10);
        expected_report.insert("World", "Foo");
        expected_report.insert("Bar", 20);

        assert_eq!(report1, expected_report);
    }

    #[test]
    fn test_json_reporter_sum() {
        let report = {
            let mut report1 = JsonReporter::default();
            report1.insert("Hello", 10);
            report1.insert("World", "Foo");

            let mut report2 = JsonReporter::default();
            report2.insert("Bar", 20);
            vec![report1, report2]
        };

        let merged: JsonReporter = report.into_iter().sum();
        let mut expected_report = JsonReporter::default();
        expected_report.insert("Hello", 10);
        expected_report.insert("World", "Foo");
        expected_report.insert("Bar", 20);

        assert_eq!(merged, expected_report);
    }

    #[test]
    #[should_panic]
    fn test_combine_panic() {
        let mut report1 = JsonReporter::default();
        report1.insert("Hello", 10);
        report1.insert("World", "Foo");

        let mut report2 = JsonReporter::default();
        report2.insert("Hello", 20);

        report1.merge(report2);
    }
    #[test]
    fn test_json_report_serialize_order() {
        let mut reporter = JsonReporter::default();
        reporter.insert("umi_bases_with_q30_frac", 0.97);
        reporter.insert("bc_bases_with_q30_frac", 0.92);
        reporter.insert("A_homopolymer_frac", 0.01);
        reporter.insert("reference_genomes", "GRCh38");
        reporter.insert("total_read_pairs", 1000);
        let mut writer = Vec::new();
        reporter.to_writer(&mut writer, SerdeFormat::Json).unwrap();
        let expected = r#"{
  "A_homopolymer_frac": 0.01,
  "bc_bases_with_q30_frac": 0.92,
  "reference_genomes": "GRCh38",
  "total_read_pairs": 1000,
  "umi_bases_with_q30_frac": 0.97
}"#;
        assert_eq!(std::str::from_utf8(&writer).unwrap(), expected);
        assert_eq!(serde_json::to_string_pretty(&reporter).unwrap(), expected);
        let deserialized: JsonReporter = serde_json::from_str(expected).unwrap();
        assert_eq!(deserialized, reporter);
    }

    #[test]
    fn test_set_macro() {
        assert_eq!(set![1, 2, 3, 4], vec![1, 2, 3, 4].into_iter().collect());
    }
}
