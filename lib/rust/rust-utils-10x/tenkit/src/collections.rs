
use std::collections::{HashMap, HashSet};
use std::hash::BuildHasherDefault;
use fxhash::FxHasher;

pub type FxHashMap<K, V> = HashMap<K, V, BuildHasherDefault<FxHasher>>;
pub type FxHashSet<V> = HashSet<V, BuildHasherDefault<FxHasher>>;