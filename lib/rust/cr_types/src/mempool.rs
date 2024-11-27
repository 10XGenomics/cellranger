//! A thread-safe, fixed-size pool of T.
//! This type can be used to provide reusable allocations that can be passed
//! between threads, such as in-memory buffers.
//!
//! The pool will never allocate more than the specified number of members,
//! and requesting a member from the pool will block until one is available.

use std::ops::{Deref, DerefMut};
use std::sync::mpsc::{sync_channel, Receiver, SyncSender};

pub struct MemPool<T: Send> {
    recv: Receiver<T>,
    load: SyncSender<T>,
}

impl<T: Send> MemPool<T> {
    /// Initialize a new pool of size count.
    /// Calls init repeatedly to produce the items.
    pub fn new<F: FnMut() -> T>(count: usize, mut init: F) -> Self {
        let (load, recv) = sync_channel(count);
        for _ in 0..count {
            load.send(init()).unwrap();
        }
        Self { recv, load }
    }

    /// Fetch a member from the pool, blocking until one is available.
    pub fn fetch(&self) -> PoolMember<T> {
        PoolMember {
            val: Some(self.recv.recv().unwrap()),
            return_to_pool: self.load.clone(),
        }
    }
}

/// A value of type T that will be returned to the source memory pool when dropped.
pub struct PoolMember<T: Send> {
    val: Option<T>,
    return_to_pool: SyncSender<T>,
}

impl<T: Send> Drop for PoolMember<T> {
    fn drop(&mut self) {
        // This can only fail if the entire pool itself has been dropped.
        let _ = self.return_to_pool.send(self.val.take().unwrap());
    }
}

impl<T: Send> Deref for PoolMember<T> {
    type Target = T;
    fn deref(&self) -> &Self::Target {
        self.val.as_ref().unwrap()
    }
}

impl<T: Send> DerefMut for PoolMember<T> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        self.val.as_mut().unwrap()
    }
}
