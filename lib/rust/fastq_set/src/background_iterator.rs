#![deny(missing_docs)]

use itertools::Itertools;
use std::sync::mpsc::{Receiver, sync_channel};
use std::thread::{JoinHandle, Scope, ScopedJoinHandle, spawn};

/// Abstract over scoped or regular thread join handles.
enum AnyJoinHandle<'scope> {
    Regular(JoinHandle<()>),
    Scoped(ScopedJoinHandle<'scope, ()>),
}

impl<'scope> AnyJoinHandle<'scope> {
    /// Join on this handle. If the thread panicked, return Err with the
    /// panic message.
    fn join(self) -> Result<(), String> {
        let Err(e) = (match self {
            Self::Regular(h) => h.join(),
            Self::Scoped(h) => h.join(),
        }) else {
            return Ok(());
        };
        // need to do a little dance here to get the panic message
        let msg = match e.downcast_ref::<&'static str>() {
            Some(s) => (*s).to_string(),
            None => match e.downcast_ref::<String>() {
                Some(s) => s.clone(),
                None => "unknown panic on BackgroundIterator worker thread".to_string(),
            },
        };
        Err(msg)
    }
}

/// Execute an iterator on a worker thread.
///
/// The worker will buffer a finite number of processed items, then send that
/// buffer on for further processing.
///
/// To unpack the items one by one, just call .flatten() on the output.
pub struct BackgroundIterator<'scope, T> {
    rx: Receiver<Option<Vec<T>>>,
    done: bool,
    handle: Option<AnyJoinHandle<'scope>>,
}

impl<'scope, T> BackgroundIterator<'scope, T> {
    // Join on the producer thread dies, and propagate the panic
    fn join_and_propagate_panic(&mut self) {
        match self.handle.take().unwrap().join() {
            Ok(()) => (),
            Err(msg) => {
                // sender panicked - we will re-panic the panic.
                panic!("{}", msg);
            }
        }
    }
}

impl<'scope, T: Send + 'static> Iterator for BackgroundIterator<'scope, T> {
    type Item = Vec<T>;

    fn next(&mut self) -> Option<Vec<T>> {
        if self.done {
            return None;
        }

        match self.rx.recv() {
            Ok(Some(v)) => Some(v),
            Ok(None) | Err(_) => {
                // if the producer thread dies, then the sender thread must have panicked
                // propagate the panic, otherwise we will silently continue with likely incomplete results
                self.done = true;
                self.join_and_propagate_panic();
                None
            }
        }
    }
}

impl<'scope, T: 'static + Send> BackgroundIterator<'scope, T> {
    /// Iterate through `itr` on a newly created thread, and send items back to the returned
    /// `BackgroundIterator` for consumption on the calling thread. The worker thread will
    /// buffer max_read_ahead items before sending to the consuming thread.
    /// The worker will pre-buffer at most one buffer.
    pub fn new<I: 'static + Send + Iterator<Item = T>>(
        itr: I,
        max_read_ahead: usize,
    ) -> BackgroundIterator<'scope, T> {
        let (tx, rx) = sync_channel::<Option<Vec<T>>>(1);
        let handle = spawn(move || {
            for chunk in &itr.chunks(max_read_ahead) {
                let mut buf = Vec::with_capacity(max_read_ahead);
                buf.extend(chunk);
                if tx.send(Some(buf)).is_err() {
                    return;
                };
            }

            tx.send(None).unwrap();
        });

        BackgroundIterator {
            rx,
            handle: Some(AnyJoinHandle::Regular(handle)),
            done: false,
        }
    }
}
impl<'scope, T: 'static + Send> BackgroundIterator<'scope, T> {
    /// Iterate through `itr` on a scoped thread, and send items back to the returned
    /// `BackgroundIterator` for consumption on the calling thread. The worker thread will
    /// buffer max_read_ahead items before sending to the consuming thread.
    /// The worker will pre-buffer at most one buffer.
    pub fn new_scoped<I: 'scope + Send + Iterator<Item = T>>(
        scope: &'scope Scope<'scope, '_>,
        itr: I,
        max_read_ahead: usize,
    ) -> BackgroundIterator<'scope, T> {
        let (tx, rx) = sync_channel::<Option<Vec<T>>>(1);
        let handle = scope.spawn(move || {
            for chunk in &itr.chunks(max_read_ahead) {
                let mut buf = Vec::with_capacity(max_read_ahead);
                buf.extend(chunk);
                if tx.send(Some(buf)).is_err() {
                    return;
                };
            }

            tx.send(None).unwrap();
        });

        BackgroundIterator {
            rx,
            handle: Some(AnyJoinHandle::Scoped(handle)),
            done: false,
        }
    }
}

/// Extension trait to allow backgrounding any iterator.
pub trait BackgroundIteratorExt<I, T>
where
    I: Iterator<Item = T> + Send + 'static,
    T: Send + 'static,
{
    /// Run this iter on a background thread, accumulating buffers of results.
    fn background(self, buffer_size: usize) -> BackgroundIterator<'static, T>;
}

/// Extension trait to allow backgrounding any iterator on a scoped thread.
pub trait BackgroundIteratorScopedExt<'scope, I, T>
where
    I: Iterator<Item = T> + Send + 'scope,
    T: Send + 'static,
{
    /// Run this iter on a scoped thread, accumulating buffers of results.
    fn background_scoped(
        self,
        scope: &'scope Scope<'scope, '_>,
        buffer_size: usize,
    ) -> BackgroundIterator<'scope, T>;
}

impl<I, T> BackgroundIteratorExt<I, T> for I
where
    I: Iterator<Item = T> + Send + 'static,
    T: Send + 'static,
{
    fn background(self, buffer_size: usize) -> BackgroundIterator<'static, T> {
        BackgroundIterator::new(self, buffer_size)
    }
}
impl<'scope, I, T> BackgroundIteratorScopedExt<'scope, I, T> for I
where
    I: Iterator<Item = T> + Send + 'scope,
    T: Send + 'static,
{
    fn background_scoped(
        self,
        scope: &'scope Scope<'scope, '_>,
        buffer_size: usize,
    ) -> BackgroundIterator<'scope, T> {
        BackgroundIterator::new_scoped(scope, self, buffer_size)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use std::time::Duration;

    #[test]
    #[should_panic]
    fn panic_bg_iter() {
        let v = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9];

        let iter = (0..11usize).map(move |i| {
            assert!(v[i] != 5, "simulated panic");
            v[i]
        });

        let bg_iter = BackgroundIterator::new(iter, 5);

        let mut sum = 0;
        for v in bg_iter.flatten() {
            sum += v;
        }
        println!("sum :{sum}");
    }

    #[test]
    fn bg_iter() {
        let v = [5; 10];

        let iter = (0..10usize).map(move |i| v[i]);

        let bg_iter = BackgroundIterator::new(iter, 5);

        let mut sum = 0;
        for v in bg_iter.flatten() {
            sum += v;
        }

        assert_eq!(sum, 50);
    }

    #[test]
    fn slow_reader_buffer() {
        let n_send = 1_000_000;

        // big iterator with lots of buffering
        let iter = (0..n_send).map(|i| i * i);
        let bg_iter = BackgroundIterator::new(iter, n_send >> 2);

        // slow reader -- the sender thread will complete before the reading is complete
        // make sure we get all the items

        let mut n_read = 0;
        for v in bg_iter.flatten() {
            assert_eq!(v, n_read * n_read);
            n_read += 1;

            // make sure the reader goes slowly
            if n_read % 50_000 == 0 {
                std::thread::sleep(Duration::from_millis(1));
            }
        }

        assert_eq!(n_send, n_read);
    }
}
