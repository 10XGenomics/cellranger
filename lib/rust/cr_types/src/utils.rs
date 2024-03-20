use anyhow::{Context, Result};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::ops::{Add, Div};
use std::path::Path;
use std::str::FromStr;

/// Read a text file with one record per line into a vector of records.
pub fn load_txt<T: FromStr>(path: &Path) -> Result<Vec<T>>
where
    <T as FromStr>::Err: Sync + Send + std::error::Error + 'static,
{
    BufReader::new(File::open(path).with_context(|| path.display().to_string())?)
        .lines()
        .map(|line| Ok(line?.parse::<T>()?))
        .collect()
}

#[derive(Eq, PartialEq, Debug, Clone, Copy)]
pub enum BisectDirection {
    Left,
    Right,
}

pub fn bisect<T, U, F>(haystack: &[T], needle: U, key_func: &F, direction: BisectDirection) -> usize
where
    U: Ord,
    F: Fn(&T) -> U,
{
    // NOTE: assumes that input is sorted by key_func
    match haystack.binary_search_by_key(&needle, key_func) {
        Ok(init_idx) => {
            // found an exact match -> seek to left / right if necessary
            let mut idx = init_idx as i32;
            let inc = match direction {
                BisectDirection::Right => 1,
                BisectDirection::Left => -1,
            };
            while (idx >= 0) && (idx < haystack.len() as i32) {
                if key_func(&haystack[idx as usize]) != needle {
                    break;
                }
                idx += inc;
            }
            if direction == BisectDirection::Left {
                idx += 1;
            }
            idx as usize
        }
        Err(idx) => idx, // no exact match -> this is the only correct insertion point
    }
}

/// Calculate the median of a sorted slice.
/// Returns the mean (rounding down) of the two middle values if the slice has an even length.
/// Returns None if the slice is empty.
pub fn calculate_median_of_sorted<T>(xs: &[T]) -> Option<T>
where
    T: Copy + Add<Output = T> + Div<Output = T> + From<u8>,
{
    let n = xs.len();
    if n == 0 {
        None
    } else if n % 2 == 0 {
        let i = n / 2;
        Some((xs[i - 1] + xs[i]) / T::from(2))
    } else {
        Some(xs[n / 2])
    }
}

/// Find intersection of two sorted, deduped vectors
pub fn intersect_vecs<T: PartialOrd + Eq + Clone>(x: &[T], y: &[T]) -> Vec<T> {
    let mut i = 0;
    let mut j = 0;
    let mut result = Vec::new();

    while i < x.len() && j < y.len() {
        if x[i] < y[j] {
            i += 1;
        } else if y[j] < x[i] {
            j += 1;
        } else {
            result.push(x[i].clone());
            i += 1;
            j += 1;
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bisect() {
        let v = [-1, -1, 0, 2, 4];
        let id = |x: &i32| *x;
        assert_eq!(bisect(&v, -2, &id, BisectDirection::Left), 0);
        assert_eq!(bisect(&v, -2, &id, BisectDirection::Right), 0);
        assert_eq!(bisect(&v, -1, &id, BisectDirection::Left), 0);
        assert_eq!(bisect(&v, -1, &id, BisectDirection::Right), 2);
        assert_eq!(bisect(&v, 0, &id, BisectDirection::Left), 2);
        assert_eq!(bisect(&v, 0, &id, BisectDirection::Right), 3);
        assert_eq!(bisect(&v, 3, &id, BisectDirection::Left), 4);
        assert_eq!(bisect(&v, 3, &id, BisectDirection::Right), 4);
        assert_eq!(bisect(&v, 5, &id, BisectDirection::Left), 5);
        assert_eq!(bisect(&v, 5, &id, BisectDirection::Right), 5);
        // NOTE: input is assumed to be sorted by key_func, so only apply functions that maintain sort
        let minus = |x: &i32| (*x) - 1;
        assert_eq!(bisect(&v, 0, &minus, BisectDirection::Left), 3);
        assert_eq!(bisect(&v, 0, &minus, BisectDirection::Right), 3);
    }

    #[test]
    fn test_calculate_median_of_sorted() {
        assert_eq!(calculate_median_of_sorted::<usize>(&[]), None);
        assert_eq!(calculate_median_of_sorted(&[1]), Some(1));
        assert_eq!(calculate_median_of_sorted(&[1, 2]), Some(1));
        assert_eq!(calculate_median_of_sorted(&[2, 2]), Some(2));
        assert_eq!(calculate_median_of_sorted(&[1, 2, 3]), Some(2));
    }

    #[test]
    fn test_intersect() {
        assert_eq!(intersect_vecs(&[0u64, 1], &[2, 3]), Vec::<u64>::new());
        assert_eq!(intersect_vecs(&[0, 2, 4], &[1, 2, 5]), vec![2]);
    }
}
