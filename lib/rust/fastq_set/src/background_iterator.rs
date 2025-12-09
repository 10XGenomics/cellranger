#![deny(missing_docs)]

use std::sync::mpsc::{Receiver, sync_channel};
use std::thread::{JoinHandle, spawn};

/// Execute an iterator on a worker thread, which can work ahead a configurable number of items.
pub struct BackgroundIterator<T> {
    rx: Receiver<Option<T>>,
    done: bool,
    handle: Option<JoinHandle<()>>,
}

impl<T> BackgroundIterator<T> {
    // Join on the producer thread dies, and propagate the panic
    fn join_and_propagate_panic(&mut self) {
        match self.handle.take().unwrap().join() {
            Ok(()) => (),
            Err(e) => {
                // sender panicked - we will re-panic the panic.
                // need to do a little dance here to get the panic message
                let msg = match e.downcast_ref::<&'static str>() {
                    Some(s) => (*s).to_string(),
                    None => match e.downcast_ref::<String>() {
                        Some(s) => s.clone(),
                        None => "unknown panic on BackgroundIterator worker thread".to_string(),
                    },
                };
                panic!("{}", msg);
            }
        }
    }
}

impl<T: Send> Iterator for BackgroundIterator<T> {
    type Item = T;

    fn next(&mut self) -> Option<T> {
        if self.done {
            return None;
        }

        match self.rx.recv() {
            Ok(Some(v)) => Some(v),
            Ok(None) => {
                self.done = true;
                self.join_and_propagate_panic();
                None
            }
            // if the producer thread dies, then the sender thread must have panicked
            // propagate the panic, otherwise we will silently continue with likely incomplete results
            Err(_) => {
                self.done = true;
                self.join_and_propagate_panic();
                None
            }
        }
    }
}

impl<T: 'static + Send> BackgroundIterator<T> {
    /// Iterate through `itr` on a newly created thread, and send items back to the returned
    /// `BackgroundIterator` for consumption on the calling thread. The worker thread will
    /// continue to produce elements until it is `max_read_ahead` items ahead of the consumer iterator.
    pub fn new<I: 'static + Send + Iterator<Item = T>>(
        itr: I,
        max_read_ahead: usize,
    ) -> BackgroundIterator<T> {
        let (tx, rx) = sync_channel::<Option<T>>(max_read_ahead);
        let handle = spawn(move || {
            for item in itr {
                if tx.send(Some(item)).is_err() {
                    return;
                };
            }

            tx.send(None).unwrap();
        });

        BackgroundIterator {
            rx,
            handle: Some(handle),
            done: false,
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

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
        for v in bg_iter {
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
        for v in bg_iter {
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
        for v in bg_iter {
            assert_eq!(v, n_read * n_read);
            n_read += 1;

            // make sure the reader goes slowly
            if n_read % 50_000 == 0 {
                std::thread::sleep(std::time::Duration::from_millis(1));
            }
        }

        assert_eq!(n_send, n_read);
    }
}
