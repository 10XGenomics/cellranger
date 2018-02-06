//
// Copyright (c) 2017 10x Genomics, Inc. All rights reserved.
//

use std::collections::BTreeMap;
use std::collections::btree_map::IterMut;
use std::iter::FromIterator;
use std::cmp;
use std::f64;

use serde_json::Value;

const TOP_N_SEQS: usize = 5;
const INSERT_SIZE_BIN_MAX: u64 = 1550;
const INSERT_SIZE_BIN_STEP: u64 = 50;

pub trait Metric {
    fn new() -> Self;
    fn merge(&mut self, other: &Self);
    fn report(&self) -> Value;
}

#[derive(Serialize, Deserialize)]
pub struct PercentMetric {
    numerator: i32,
    denominator: i32,
}

impl PercentMetric {
    pub fn add(&mut self, data: &bool) {
        if *data { self.numerator += 1; }
        self.denominator += 1;
    }
}

impl Metric for PercentMetric {
    fn new() -> Self {
        PercentMetric {
            numerator: 0,
            denominator: 0,
        }
    }

    fn merge(&mut self, other: &Self) {
        self.numerator += other.numerator;
        self.denominator += other.denominator;
    }

    fn report(&self) -> Value {
        return json!(self.numerator as f64 / self.denominator as f64)
    }
}

#[derive(Serialize, Deserialize)]
pub struct SequenceDistMetric {
    histogram: SequenceHistogram,
}

impl SequenceDistMetric {
    pub fn add(&mut self, data: &String) {
        self.histogram.increment(data);
    }

    pub fn get_seqs(&self) -> Vec<&String> {
        let mut seqs = self.histogram.keys();
        seqs.sort();
        return seqs
    }

    pub fn report_num_seqs(&self) -> Value {
        return json!(self.histogram.keys().len())
    }

    pub fn report_top_n(&self) -> Value {
        return json!(self.histogram.top_n(TOP_N_SEQS))
    }

    pub fn report_effective_diversity(&self) -> Value {
        return json!(self.histogram.effective_diversity())
    }
}

impl Metric for SequenceDistMetric {
    fn new() -> Self {
        SequenceDistMetric {
            histogram: SequenceHistogram::new(),
        }
    }

    fn merge(&mut self, other: &Self) {
        self.histogram.merge(&other.histogram);
    }

    fn report(&self) -> Value {
        return json!("Unimplemented")
    }
}

#[derive(Serialize, Deserialize)]
pub struct InsertSizeMetric {
    histogram: LengthHistogram,
}

impl InsertSizeMetric {
    pub fn add(&mut self, data: &u64) {
        self.histogram.increment(*data);
    }

    pub fn get_bins(&self) -> Vec<u64> {
        let mut bins = Vec::new();
        for i in 0..INSERT_SIZE_BIN_MAX/INSERT_SIZE_BIN_STEP { bins.push(i*INSERT_SIZE_BIN_STEP); }
        return bins
    }

    pub fn report_median(&self) -> Value {
        return json!(self.histogram.get_median())
    }

    pub fn report_iqr(&self) -> Value {
        return json!(self.histogram.get_iqr())
    }

    pub fn report_binned(&self) -> Value {
        let cutoffs = self.get_bins();
        return json!(self.histogram.get_binned(&cutoffs))
    }
}

impl Metric for InsertSizeMetric {
    fn new() -> Self {
        InsertSizeMetric {
            histogram: LengthHistogram::new(),
        }
    }

    fn merge(&mut self, other: &Self) {
        self.histogram.merge(&other.histogram);
    }

    fn report(&self) -> Value {
        return json!("Unimplemented")
    }
}

pub trait MetricGroup {
    fn new() -> Self;
    fn merge(&mut self, other: &Self);
    fn report(&self) -> BTreeMap<String, Value>;
}

#[derive(Serialize, Deserialize)]
pub struct PrefixGroup<T: MetricGroup> {
    mapping: BTreeMap<String, T>,
}

impl<T: MetricGroup> PrefixGroup<T> {
    pub fn new(prefixes: &[String]) -> Self {
        let mut mapping = BTreeMap::new();
        for prefix in prefixes {
            mapping.insert(prefix.clone(), T::new());
        }
        return PrefixGroup { mapping: mapping }
    }

    pub fn merge(&mut self, other: &Self) {
        for (prefix, other_group) in &other.mapping {
            let mut group = self.mapping.entry(prefix.clone()).or_insert(T::new());
            group.merge(&other_group);
        }
    }

    pub fn iter_mut(&mut self) -> IterMut<String, T> {
        self.mapping.iter_mut()
    }

    pub fn report(&self) -> BTreeMap<String, Value> {
        let mut results = BTreeMap::new();
        for (prefix, group) in &self.mapping {
            for (name, metric) in group.report() {
                results.insert(format!("{}_{}", prefix, name), metric);
            }
        }
        return results
    }
}

/* Additional Data Structures */

#[derive(Serialize, Deserialize)]
pub struct SequenceHistogram {
    hist: BTreeMap<String, u64>,
}

impl SequenceHistogram {
    pub fn new() -> SequenceHistogram {
        let hist = BTreeMap::new();
        return SequenceHistogram {
            hist: hist,
        }
    }

    pub fn increment(&mut self, seq: &String) {
        let ref mut hist = self.hist;
        let count = hist.entry(seq.to_owned()).or_insert(0);
        *count += 1;
    }

    //pub fn get_count(&self, seq: &String) -> u64 {
    //    return self.hist.get(seq).unwrap_or(&0).to_owned()
    //}

    pub fn merge(&mut self, other: &Self) {
        let ref mut hist = self.hist;
        for (seq, other_count) in &other.hist {
            let count = hist.entry(seq.to_owned()).or_insert(0);
            *count += *other_count;
        }
    }

    pub fn top_n(&self, n: usize) -> BTreeMap<&String, &u64> {
        let mut seqs = Vec::new();
        for k in self.hist.keys() {
            seqs.push(k);
        }
        seqs.sort_by_key(|s| self.hist.get(*s).unwrap());
        seqs.reverse();
        let mut top_n_seqs = BTreeMap::new();
        for seq in &seqs[0..cmp::min(n, seqs.len())] {
            top_n_seqs.insert(*seq, self.hist.get(*seq).unwrap());
        }
        return top_n_seqs
    }

    pub fn effective_diversity(&self) -> f64 {
        // inverse Simpson index
        let mut s = 0_f64;
        let mut s2 = 0_f64;
        for (_, count) in &self.hist {
            s += *count as f64;
            s2 += (*count as f64).powi(2);
        }
        return s.powi(2) / s2
    }

    pub fn keys(&self) -> Vec<&String> {
        return Vec::from_iter(self.hist.keys())
    }
}

#[derive(Serialize, Deserialize)]
pub struct LengthHistogram {
    hist: BTreeMap<u64, u64>,
}

impl LengthHistogram {
    pub fn new() -> LengthHistogram {
        let hist = BTreeMap::new();
        LengthHistogram {
            hist: hist,
        }
    }

    pub fn increment(&mut self, length: u64) {
        let ref mut hist = self.hist;
        let count = hist.entry(length).or_insert(0);
        *count += 1;
    }

    //pub fn get_count(&self, length: &u64) -> u64 {
    //    return self.hist.get(length).unwrap_or(&0).to_owned()
    //}

    pub fn merge(&mut self, other: &Self) {
        let ref mut hist = self.hist;
        for (seq, other_count) in &other.hist {
            let count = hist.entry(seq.to_owned()).or_insert(0);
            *count += *other_count;
        }
    }

    pub fn percentile(&self, percentile: f64) -> f64 {
        assert!(percentile >= 0.0 && percentile <= 100.0);

        let mut n = 0;
        for count in self.hist.values() { n += *count; }
        if n == 0 {
            return f64::NAN
        }

        let h = (n-1) as f64 * (percentile / 100.0);
        let mut lower_value = 0;
        let mut saw_floor = false;
        let mut cum_sum = 0;

        for (value, count) in self.hist.iter() {
            cum_sum += *count;
            if cum_sum as f64 > h.floor() && !saw_floor {
                lower_value = *value;
                saw_floor = true;
            }
            if cum_sum as f64 > h.ceil() {
                assert!(saw_floor);
                return lower_value as f64 + (h - h.floor()) * (value - lower_value) as f64
            }
        }

        return f64::NAN
    }

    pub fn get_median(&self) -> f64 {
        return self.percentile(50.0)
    }

    pub fn get_iqr(&self) -> f64 {
        return self.percentile(75.0) - self.percentile(25.0)
    }

    pub fn get_binned(&self, cutoffs: &[u64]) -> BTreeMap<u64, u64> {
        let mut binned_histogram = BTreeMap::new();
        for cutoff in cutoffs { binned_histogram.insert(cutoff.clone(), 0); }
        // assumes histogram is sorted
        let mut curr_cutoff_idx = 0;
        for (value, count) in self.hist.iter() {
            if (curr_cutoff_idx < cutoffs.len() - 1) && (value >= &cutoffs[curr_cutoff_idx+1]) {
                curr_cutoff_idx += 1;
            }
            let mut binned_count = binned_histogram.get_mut(&cutoffs[curr_cutoff_idx]).unwrap();
            *binned_count += *count;
        }

        return binned_histogram
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_percent_metric() {
        let mut m0 = PercentMetric::new();
        m0.add(&false);
        assert_eq!(m0.report(), 0.0);
        m0.add(&true);
        assert_eq!(m0.report(), 0.5);
        let mut m1 = PercentMetric::new();
        m1.add(&false);
        m1.add(&false);
        m0.merge(&m1);
        assert_eq!(m0.report(), 0.25);
    }

    #[test]
    fn test_sequence_histogram() {
        let mut m0 = SequenceHistogram::new();
        let foo = "foo".into();
        let bar = "bar".into();
        let baz = "baz".into();
        m0.increment(&foo);
        m0.increment(&foo);
        m0.increment(&bar);
        assert_eq!(m0.get_count(&foo), 2);
        assert_eq!(m0.get_count(&bar), 1);
        assert_eq!(m0.top_n(1).keys().next().unwrap().to_owned(), &foo);
        assert_eq!(m0.effective_diversity(), 1.8);
        let mut m1 = SequenceHistogram::new();
        m1.increment(&bar);
        m1.increment(&baz);
        m0.merge(&m1);
        assert_eq!(m0.get_count(&bar), 2);
        assert_eq!(m0.get_count(&baz), 1);
    }

    #[test]
    fn test_length_histogram() {
        let mut m0 = LengthHistogram::new();
        m0.increment(0);
        m0.increment(10);
        m0.increment(10);
        m0.increment(20);
        assert_eq!(m0.get_count(&0), 1);
        assert_eq!(m0.get_count(&10), 2);
        assert_eq!(m0.get_count(&20), 1);
        assert_eq!(m0.percentile(0.0), 0.0);
        assert_eq!(m0.percentile(25.0), 7.5);
        assert_eq!(m0.percentile(50.0), 10.0);
        assert_eq!(m0.percentile(75.0), 12.5);
        assert_eq!(m0.percentile(100.0), 20.0);
        let mut m1 = LengthHistogram::new();
        m1.increment(0);
        m0.merge(&m1);
        assert_eq!(m0.get_count(&0), 2);
    }
}
