
use std::collections::HashMap;
use std::hash::BuildHasherDefault;
use fxhash::FxHasher;

pub type FxHashMap<K, V> = HashMap<K, V, BuildHasherDefault<FxHasher>>;
// use serde_json::Value;

pub trait Metric {
    fn new() -> Self;
    fn merge(&mut self, other: &Self);
    // fn report(&self) -> Value;
}

#[derive(Serialize, Deserialize, Default, Debug)]
pub struct Histogram {
    distribution: FxHashMap<i64, u64>,
    zeroth_moment: i64,
    first_moment: i64,
}

impl Histogram {
    pub fn increment<T>(&mut self, v: T) where i64: From<T> {
        let val = i64::from(v);
        self.zeroth_moment += 1;
        self.first_moment += val;
        *self.distribution.entry(val).or_insert(0) += 1;
    }
}


impl Metric for Histogram {
    fn new() -> Self {
        Histogram::default()
    }

    fn merge(&mut self, other: &Self) {
        for (&key, &val) in &other.distribution {
            *self.distribution.entry(key).or_insert(0) += val;
        }
        self.zeroth_moment += other.zeroth_moment;
        self.first_moment += other.first_moment;
    }
}
