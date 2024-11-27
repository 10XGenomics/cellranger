use num_traits::PrimInt;
use std::fmt::Display;

// Compute nx of a list of positive integers
// Helper function to compute n50 and n90 scores.
pub fn nx<'a, T, I>(items: I, fraction: f64) -> Option<T>
where
    T: 'a + PrimInt + Display,
    I: IntoIterator<Item = &'a T>,
{
    assert!(fraction > 0f64 && fraction < 1f64);

    let (rev_sorted_items, sum) = {
        let mut items_owned = Vec::new();
        let mut s = 0f64;
        for item in items {
            assert!(
                *item > T::zero(),
                "Found a number that is not positive while computing Nx : {}",
                *item
            );
            s += item.to_f64().unwrap();
            items_owned.push(*item);
        }
        items_owned.sort_unstable_by(|x, y| x.cmp(y).reverse());
        (items_owned, s)
    };

    let mut cumulative_sum = 0f64;
    let cutoff_sum = sum * fraction;
    for item in rev_sorted_items {
        cumulative_sum += item.to_f64().unwrap();
        if cumulative_sum >= cutoff_sum {
            return Some(item);
        }
    }
    None // Empty Iterator
}

/// Compute the N50 from a list of positive numbers. The numbers
/// can be any of the primary integer type.
///
/// # Inputs
/// - `items`: Any type which can be converted to an iterator over a primary integer
///            type T. E.g &Vec<_>
///
/// # Outputs
/// - `Option<T>`: `None` variant is returned only when the input is an empty iterator.
///                The caller should handle this case explicitly.
///
/// # Panics
/// - If any of the numbers is <=0. This likely points to a bug in the caller code,
///   since one should not be attempting to compute N50 of a list of non-positive
///   numbers. Hence decided on a panic over returning a Result.
///
/// # Warning
/// - It is implicitly assumed that the sum of all the elements will
///   fit within an f64.
///
/// # Example
/// ```rust
/// use stats::n50;
/// let items: Vec<i32> = vec![2, 3, 4, 5, 6, 7, 8, 9, 10];
/// assert_eq!(n50(&items), Some(8));
/// let empty_items: Vec<i32> = vec![];
/// assert_eq!(n50(&empty_items), None);
/// ```
pub fn n50<'a, T, I>(items: I) -> Option<T>
where
    T: 'a + PrimInt + Display,
    I: IntoIterator<Item = &'a T>,
{
    nx(items, 0.5f64)
}

/// Compute the N90 from a list of positive numbers. The numbers
/// can be any of the primary integer type.
///
/// # Inputs
/// - `items`: Any type which can be converted to an iterator over a primary integer
///            type T. E.g &Vec<_>
///
/// # Outputs
/// - `Option<T>`: `None` variant is returned only when the input is an empty iterator.
///                The caller should handle this case explicitly.
///
/// # Panics
/// - If any of the numbers is <=0. This likely points to a bug in the caller code,
///   since one should not be attempting to compute N90 of a list of non-positive
///   numbers. Hence decided on a panic over returning a Result.
///
/// # Warning
/// - It is implicitly assumed that the sum of all the elements will
///   fit within an f64.
///
/// # Example
/// ```rust
/// use stats::n90;
/// let items: Vec<i32> = vec![2, 3, 4, 5, 6, 7, 8, 9, 10];
/// assert_eq!(n90(&items), Some(4));
/// let empty_items: Vec<i32> = vec![];
/// assert_eq!(n90(&empty_items), None);
/// ```
pub fn n90<'a, T, I>(items: I) -> Option<T>
where
    T: 'a + PrimInt + Display,
    I: IntoIterator<Item = &'a T>,
{
    nx(items, 0.9f64)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_n50_1() {
        let items = vec![100, 70, 60, 50, 50, 40, 30];
        assert_eq!(n50(&items), Some(60));
    }
    #[test]
    fn test_n50_2() {
        let items = vec![70, 60, 50, 40, 30, 100, 50];
        assert_eq!(n50(&items), Some(60));
    }
    #[test]
    fn test_n50_3() {
        let items = vec![68, 90, 11, 50, 15, 57, 27, 67, 24, 45];
        assert_eq!(n50(&items), Some(57));
    }
    #[test]
    fn test_n50_empty() {
        let items: Vec<i32> = Vec::new();
        assert_eq!(n50(&items), None);
    }
    #[test]
    #[should_panic]
    fn test_n50_panic_zero() {
        let items = vec![68, 90, 0, 50];
        let _ = n50(&items);
    }
    #[test]
    #[should_panic]
    fn test_n50_panic_negative() {
        let items = vec![68, 90, -2, 50];
        let _ = n50(&items);
    }
    #[test]
    fn test_n90_1() {
        let items = vec![100, 70, 60, 50, 50, 40, 30];
        assert_eq!(n90(&items), Some(40));
    }
    #[test]
    fn test_n90_2() {
        let items = vec![70, 60, 50, 40, 30, 100, 50];
        assert_eq!(n90(&items), Some(40));
    }
}
