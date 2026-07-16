//! In-flight counter pool: a single-producer semaphore with RAII reservation tokens.
//!
//! The producer calls [`InFlightPool::acquire`] which blocks until the total in-flight
//! count drops below the caller-supplied threshold, then wraps a payload in a
//! [`Reservation`] and increments the counter. When the [`Reservation`] is dropped (or
//! consumed via [`IntoIterator`]) the counter is decremented and the producer is woken.

use anyhow::{Result, ensure};
use std::ops::{Deref, DerefMut};
use std::sync::{Arc, Condvar, Mutex};

struct Inner {
    count: Mutex<usize>,
    condvar: Condvar,
}

/// Single-producer in-flight counter pool.
///
/// Tracks the total "weight" of outstanding [`Reservation`]s. Call [`InFlightPool::acquire`]
/// to block until capacity is available, then ship the returned [`Reservation`] downstream.
pub struct InFlightPool {
    inner: Arc<Inner>,
}

impl InFlightPool {
    /// Create a new pool with an initial in-flight count of zero.
    pub fn new() -> Self {
        Self {
            inner: Arc::new(Inner {
                count: Mutex::new(0),
                condvar: Condvar::new(),
            }),
        }
    }

    /// Block until the in-flight count is below `threshold`, then reserve `amount` units
    /// and wrap `value` in a [`Reservation`].
    ///
    /// Returns an error if `amount > threshold` (the reservation could never succeed).
    pub fn acquire<T>(&self, value: T, amount: usize, threshold: usize) -> Result<Reservation<T>> {
        ensure!(
            amount <= threshold,
            "requested amount {amount} exceeds threshold {threshold}"
        );
        let mut count = self.inner.count.lock().unwrap();
        loop {
            if *count + amount <= threshold {
                *count += amount;
                return Ok(Reservation {
                    guard: ReservationGuard {
                        inner: Arc::clone(&self.inner),
                        amount,
                    },
                    value,
                });
            }
            count = self.inner.condvar.wait(count).unwrap();
        }
    }

    /// Return the current total in-flight count.
    pub fn in_flight(&self) -> usize {
        *self.inner.count.lock().unwrap()
    }
}

impl Default for InFlightPool {
    fn default() -> Self {
        Self::new()
    }
}

/// Private guard that owns the decrement logic. By keeping `Drop` here rather than on
/// `Reservation<T>`, `Reservation<T>` can be freely destructured in `into_iter`.
struct ReservationGuard {
    inner: Arc<Inner>,
    amount: usize,
}

impl Drop for ReservationGuard {
    fn drop(&mut self) {
        let mut count = self.inner.count.lock().unwrap();
        *count = count.saturating_sub(self.amount);
        drop(count);
        self.inner.condvar.notify_one();
    }
}

/// RAII reservation returned by [`InFlightPool::acquire`].
///
/// Derefs to the wrapped value. The in-flight count is held until this reservation is
/// dropped or consumed via [`TryIntoIterator`].
pub struct Reservation<T> {
    guard: ReservationGuard,
    value: T,
}

impl<T> Reservation<T> {
    /// The weight that was reserved when this token was created.
    pub fn amount(&self) -> usize {
        self.guard.amount
    }
}

impl<T> Deref for Reservation<T> {
    type Target = T;

    fn deref(&self) -> &T {
        &self.value
    }
}

impl<T> DerefMut for Reservation<T> {
    fn deref_mut(&mut self) -> &mut T {
        &mut self.value
    }
}

/// Iterator wrapper that keeps the [`ReservationGuard`] alive until iteration completes.
pub struct ReservationIter<Iter> {
    iter: Iter,
    _guard: ReservationGuard,
}

impl<Iter: Iterator> Iterator for ReservationIter<Iter> {
    type Item = Iter::Item;

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next()
    }
}

impl<T: IntoIterator> IntoIterator for Reservation<T> {
    type Item = T::Item;
    type IntoIter = ReservationIter<T::IntoIter>;

    fn into_iter(self) -> Self::IntoIter {
        let Reservation { guard, value } = self;
        ReservationIter {
            iter: value.into_iter(),
            _guard: guard,
        }
    }
}

impl<T: TryIntoIterator> TryIntoIterator for Reservation<T> {
    type Item = T::Item;
    type IntoIter = ReservationIter<T::IntoIter>;

    fn try_into_iter(self) -> Result<Self::IntoIter> {
        let Reservation { guard, value } = self;
        Ok(ReservationIter {
            iter: value.try_into_iter()?,
            _guard: guard,
        })
    }
}

/// A fallible version of IntoIterator.
pub trait TryIntoIterator {
    /// The type of item returned by the iterator.
    type Item;
    /// The concrete type that will serve as the iterator.
    type IntoIter: Iterator<Item = Self::Item>;

    /// Fallibly try to convert this value into an iterator.
    fn try_into_iter(self) -> Result<Self::IntoIter>;
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::thread;
    use std::time::Duration;

    #[test]
    fn acquire_and_drop_tracks_count() {
        let pool = InFlightPool::new();
        assert_eq!(pool.in_flight(), 0);

        let res = pool.acquire("hello", 10, 100).unwrap();
        assert_eq!(pool.in_flight(), 10);
        assert_eq!(res.amount(), 10);
        assert_eq!(*res, "hello");

        drop(res);
        assert_eq!(pool.in_flight(), 0);
    }

    #[test]
    fn deref_mut_works() {
        let pool = InFlightPool::new();
        let mut res = pool.acquire(vec![1u32, 2, 3], 30, 100).unwrap();
        res.push(4);
        assert_eq!(*res, vec![1, 2, 3, 4]);
    }

    #[test]
    fn into_iter_releases_count_when_exhausted() {
        let pool = InFlightPool::new();
        let res = pool.acquire(vec![1u32, 2, 3], 30, 100).unwrap();
        assert_eq!(pool.in_flight(), 30);

        let items: Vec<_> = res.into_iter().collect();
        assert_eq!(items, vec![1, 2, 3]);
        assert_eq!(pool.in_flight(), 0);
    }

    #[test]
    fn into_iter_releases_count_when_dropped_early() {
        let pool = InFlightPool::new();
        let res = pool.acquire(vec![1u32, 2, 3, 4, 5], 50, 100).unwrap();
        assert_eq!(pool.in_flight(), 50);

        let mut iter = res.into_iter();
        assert_eq!(iter.next(), Some(1));
        drop(iter); // dropped before exhausted
        assert_eq!(pool.in_flight(), 0);
    }

    #[test]
    fn exceeds_threshold_returns_error() {
        let pool = InFlightPool::new();
        assert!(pool.acquire("big", 200, 100).is_err());
    }

    #[test]
    fn blocks_until_reservation_dropped() {
        let pool = InFlightPool::new();
        let res = pool.acquire((), 10, 10).unwrap();

        thread::spawn(move || {
            thread::sleep(Duration::from_millis(50));
            drop(res);
        });

        let _res2 = pool.acquire((), 5, 10).unwrap();
        assert_eq!(pool.in_flight(), 5);
    }

    #[test]
    fn blocks_until_enough_space_for_full_amount() {
        let pool = InFlightPool::new();
        let res = pool.acquire((), 7, 10).unwrap();

        thread::spawn(move || {
            thread::sleep(Duration::from_millis(50));
            drop(res);
        });

        let _res2 = pool.acquire((), 5, 10).unwrap();
        assert_eq!(pool.in_flight(), 5);
    }
}
