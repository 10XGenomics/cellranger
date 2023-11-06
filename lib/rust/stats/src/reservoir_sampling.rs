use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoshiro256StarStar;
use std::iter::IntoIterator;

/// Randomly choose a sample of k items from a list containing n items,
/// where n is often a very large or unknown number
/// See https://en.wikipedia.org/wiki/Reservoir_sampling
pub struct ReservoirSampler<T> {
    // `capacity` is the maximum number of items
    capacity: usize,
    // Number of items seen so far
    items_seen: usize,
    // Sampled items
    items: Vec<T>,
    // Fast and high quality psuedo random number generator
    rng: Xoshiro256StarStar,
}

impl<T> ReservoirSampler<T> {
    /// Create a new `ReservoirSampler` with the specified capacity and random seed
    pub fn new(capacity: usize, seed: u64) -> ReservoirSampler<T> {
        ReservoirSampler {
            capacity,
            items_seen: 0,
            items: Vec::with_capacity(capacity),
            rng: Xoshiro256StarStar::seed_from_u64(seed),
        }
    }

    /// Add a new item to the sampler
    pub fn add(&mut self, item: T) {
        self.items_seen += 1;

        if self.items.len() < self.capacity {
            self.items.push(item);
        } else {
            // Swap out an existing item according
            let idx = self.rng.gen_range(0..self.items_seen);
            if idx < self.items.len() {
                self.items[idx] = item;
            }
        }
    }

    /// Consume the `ReservoirSampler` and return the sampled items
    pub fn done(mut self) -> Vec<T> {
        self.items.shrink_to_fit();
        self.items
    }

    /// Number of items seen so far by the reservoir samples
    pub fn num_items_seen(&self) -> usize {
        self.items_seen
    }

    /// Reservoir sample from an iterator and return the sampled items
    pub fn sample_from_iter<I: IntoIterator<Item = T>>(
        iter: I,
        capacity: usize,
        seed: u64,
    ) -> Vec<T> {
        let mut sampler = ReservoirSampler::new(capacity, seed);
        for item in iter {
            sampler.add(item);
        }
        sampler.done()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::arbitrary::any;
    use proptest::{prop_assert, proptest};
    use std::cmp::min;
    use std::collections::HashSet;

    proptest! {
        #[test]
        fn prop_test_reservoir_sampling(
            num_items in 0usize..1000usize,
            capacity in 0usize..500usize,
            seed in any::<u64>(),
        ) {
            let sampled_items1 = ReservoirSampler::sample_from_iter(0..num_items, capacity, seed);
            prop_assert!(sampled_items1.len()==min(capacity, num_items));
            let sampled_items2 = ReservoirSampler::sample_from_iter(0..num_items, capacity, seed);
            // Repeatability
            assert_eq!(sampled_items1, sampled_items2);
            let mut seen = HashSet::new();
            for item in sampled_items1 {
                prop_assert!(!seen.contains(&item));
                prop_assert!(item < num_items);
                seen.insert(item);
            }

        }
    }
}
