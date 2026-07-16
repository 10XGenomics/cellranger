//! Measure fraction of local wall time consumed by an action.

use std::io::{Write, stderr};
use std::time::{Duration, Instant};

/// Quick one-liner function to time a single action.
///
/// Writes the elapsed time to stderr.
pub fn timed<T>(name: &str, action: impl FnOnce() -> T) -> T {
    let start = Instant::now();
    let r = action();
    let _ = writeln!(
        stderr(),
        "action time for {name}: {:.2}s",
        start.elapsed().as_secs_f64()
    );
    r
}

/// Measure fraction of lifetime spent taking some action.
///
/// Useful for measuring things like how much processing time was dedicated to
/// processing an iterator, sending items over channels, etc.
///
/// When this struct is dropped, it will report its statistics to stderr.
pub struct ActionTimer {
    start: Instant,
    accum: Duration,
    name: String,
}

impl ActionTimer {
    /// Initialize an action timer with the provided name.
    pub fn new(name: &str) -> Self {
        Self {
            start: Instant::now(),
            accum: Duration::ZERO,
            name: name.to_string(),
        }
    }

    /// Accumlate the time needed to perform action.
    pub fn time<T>(&mut self, action: impl FnOnce() -> T) -> T {
        let start = Instant::now();
        let result = action();
        self.accum += start.elapsed();
        result
    }
}

impl Drop for ActionTimer {
    /// Report fraction of elapsed wall time accumulated.
    fn drop(&mut self) {
        let total = self.start.elapsed();
        let frac = self.accum.as_secs_f64() / total.as_secs_f64();
        let _ = writeln!(
            stderr(),
            "action time for {}: {:.2}s / {:.2}s ({:.2}%)",
            self.name,
            self.accum.as_secs_f64(),
            total.as_secs_f64(),
            100. * frac
        );
    }
}

/// Wrap an iterator to measure fraction of wall time used in iterating it.
pub struct IteratorTimer<T, I: Iterator<Item = T>> {
    iter: I,
    timer: ActionTimer,
}

impl<T, I: Iterator<Item = T>> Iterator for IteratorTimer<T, I> {
    type Item = T;
    fn next(&mut self) -> Option<Self::Item> {
        self.timer.time(|| self.iter.next())
    }
}

/// Wrap an iterator to measure fraction of wall time used in iterating it.
pub fn timed_iter<T, I: Iterator<Item = T>>(iter: I, name: &str) -> IteratorTimer<T, I> {
    IteratorTimer {
        iter,
        timer: ActionTimer::new(name),
    }
}

/// Extension trait to time an iterator.
pub trait TimedIterExt<T>
where
    Self: Iterator<Item = T> + Sized,
{
    /// Iterator adapter that reports action timing for the input iterator.
    fn timed(self, name: &str) -> IteratorTimer<T, Self>;
}

impl<T, I: Iterator<Item = T>> TimedIterExt<T> for I {
    fn timed(self, name: &str) -> IteratorTimer<T, Self> {
        timed_iter(self, name)
    }
}
