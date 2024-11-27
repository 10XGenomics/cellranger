use anyhow::{anyhow, Result};
use core::iter::IntoIterator;
use cr_types::spill_vec::SpillVec;
use crossbeam_channel::{bounded, Receiver};
use crossbeam_utils::thread::{Scope, ScopedJoinHandle};
use itertools::Itertools;
use martian_filetypes::LazyFileTypeIO;
use std::any::Any;
use std::path::Path;

// FIXME -- decide where this code should live.  If it's going to stay here & not be in shardio, fix the tests below.

/// A generic processor trait that accepts items and can return an error
/// if the processing of the item fails.
pub trait Proc {
    /// The type of the items to be processed
    type Item;

    /// The type of the error produced on a processing failure
    type Err;

    /// Process one item
    fn process(&mut self, item: Self::Item) -> Result<(), Self::Err>;
}

// Given a receiver (`rec`) that can receive items of type `T` through a channel
// and a `processor` that knows how to process an item of type `T`, this function
// processes all the items that the receiver can receive.
fn process<P: Proc<Item = T, Err = E>, T: Send, E>(
    mut processor: P,
    rec: Receiver<T>,
) -> Result<P, E> {
    loop {
        let item = rec.recv(); // Block the thread until the next item is received.
        match item {
            Ok(item) => processor.process(item)?, // We received an item, process it.
            Err(_) => break,                      // The channel is disconnected
        }
    }

    Ok(processor)
}

pub const MAX_ITEMS_IN_MEM: usize = 500_000;

#[allow(clippy::type_complexity)]
fn start_processors<'a, 'b, P, T, K>(
    processors: Vec<P>,
    scope: &'a Scope<'b>,
) -> (
    crossbeam_channel::Sender<(K, T)>,
    Vec<ScopedJoinHandle<'a, Result<P>>>,
)
where
    T: 'b, // Type of the iterable, wrapped in a result.
    // The values are ordered by the key `K'
    K: 'static + Eq + Send + Sync + Clone,
    P: 'static + Send + Proc<Item = (K, T), Err = anyhow::Error>,
    // The first element of (K, Vec<T>) is the key and the second element are the items
    // associated with the key after grouping.
    <P as Proc>::Item: Send,
{
    // allow a little bit of read-ahead
    // Create a channel that can hold at most 2 messages at a time
    let (send, recv) = bounded(2);

    // Make one thread dedicated to each processor object
    let mut handles = Vec::with_capacity(processors.len());
    for processor in processors {
        // Each thread gets a copy of the receiver
        let recv = recv.clone();
        // All that this thread does is repeatedly look for any new items in the
        // receive channel and process any item it receives until the channel is
        // disconnected.
        let thread = scope.spawn(move |_| process(processor, recv));
        handles.push(thread);
    }
    (send, handles)
}

fn collect_jobs<P, T, K>(handles: Vec<ScopedJoinHandle<'_, Result<P>>>) -> Result<Vec<P>>
where
    P: Proc<Item = (K, T), Err = anyhow::Error>,
{
    // let the threads finish up & return the processors
    let mut results = Vec::with_capacity(handles.len());
    for h in handles {
        // join() Returns Result<Result<P, E>, E>
        match h.join() {
            // return the inner processor value, bubbling up an Err
            Ok(v) => results.push(v?),

            // if a thread panicked, capture the stack trace
            Err(e) => return Err(anyhow!(decipher_panic(e))),
        }
    }

    Ok(results)
}

/// Iterate over `iterable` of type `Result<T>` returning if any errors are encounterd.
/// Group items according to the `key` function.  For each group, process the data with one of the
/// processors of type `P`. Return the processor objects once all the groups have been processed.
/// Processor objects should carry any output IO handles or metrics/results accumulators as
/// state.
pub fn group_by_processor<I, F, P, T, K, S>(
    iterable: I,
    processors: Vec<P>,
    key: F,
    spill_folder: &Path,
) -> Result<Vec<P>>
where
    T: 'static + Send + Sync, // Type returned by the iterable, wrapped in a result.
    S: LazyFileTypeIO<T> + Send + Sync,
    // The values are ordered by the key `K'
    K: 'static + Eq + Send + Sync + Clone + ToString,
    P: 'static + Send + Proc<Item = (K, SpillVec<T, S>), Err = anyhow::Error>,
    // The first element of (K, Vec<T>) is the key and the second element are the items
    // associated with the key after grouping.
    <P as Proc>::Item: Send,
    F: Fn(&T) -> K,
    I: IntoIterator<Item = Result<T>>,
{
    // Create a scope where the processor threads will run
    // Note that scoped threads are a mechanism to guarantee to the compiler that spawned
    // threads will be joined before the scope ends, thereby allowing the threads to borrow
    // variables in the stack
    let r = crossbeam_utils::thread::scope(move |s| -> Result<Vec<P>> {
        let (send, handles) = start_processors(processors, s);

        for (key, group_iter) in &iterable
            .into_iter()
            .group_by(move |item_result: &Result<T>| {
                // Map an Error to None and Ok(t) to Some(key(&t)). This is needed because
                // `failure::Error` doesn't satisfy all the traits to be a valid key to group by
                // We will bubble up the error via the `group_iter`
                key(item_result.as_ref().unwrap())
            })
        {
            let spill_file = S::new(spill_folder, key.to_string());
            let mut spill_vec = SpillVec::new(MAX_ITEMS_IN_MEM, spill_file);
            for item in group_iter {
                spill_vec.push(item?)?;
            }
            let send_res = send.send((key, spill_vec));
            // stop sending if all the recievers have hung up
            // go down and get the Error/panic from the dead
            // worker
            if send_res.is_err() {
                break;
            }
        }

        // Close the send channel - this will cause the threads to exit
        drop(send);

        collect_jobs(handles)
    });

    match r {
        Ok(v) => v,
        Err(e) => Err(anyhow!(decipher_panic(e))),
    }
}

/// Note: this duplicates the logic in `group_by_processor` above with two differences:
/// - there is no SpillVec and all objects are in memory
/// - if a groupby key has a very long list of items to process that list is capped at `max_items`.
///
/// The remaining items are processed in a separate group with the same key.
///
/// Iterate over `iterable` of type `Result<T>` returning if any errors are encounterd.
/// Group items according to the `key` function.  For each group, process the data with one of the
/// processors of type `P`. Return the processor objects once all the groups have been processed.
/// Processor objects should carry any output IO handles or metrics/results accumulators as
/// state. Note that the group by operation creates an in-memory Vec<T>, which can be huge. The
/// parameter `max_items` controls this size.
pub fn group_by_processor_in_memory<I, F, P, T, K>(
    iterable: I,
    processors: Vec<P>,
    key: F,
    max_items: usize,
) -> Result<Vec<P>>
where
    T: 'static + Send + Sync, // Type returned by the iterable, wrapped in a result.
    // The values are ordered by the key `K'
    K: 'static + Eq + Send + Sync + Clone,
    P: 'static + Send + Proc<Item = (K, Vec<T>), Err = anyhow::Error>,
    // The first element of (K, Vec<T>) is the key and the second element are the items
    // associated with the key after grouping.
    <P as Proc>::Item: Send,
    F: Fn(&T) -> &K,
    I: IntoIterator<Item = Result<T>>,
{
    // Create a scope where the processor threads will run
    // Note that scoped threads are a mechanism to guarantee to the compiler that spawned
    // threads will be joined before the scope ends, thereby allowing the threads to borrow
    // variables in the stack
    let r = crossbeam_utils::thread::scope(move |s| -> Result<Vec<P>> {
        let (send, handles) = start_processors(processors, s);

        for (key, group_iter) in &iterable
            .into_iter()
            .group_by(move |item_result: &Result<T>| {
                // Map an Error to None and Ok(t) to Some(key(&t)). This is needed because
                // `failure::Error` doesn't satisfy all the traits to be a valid key to group by
                // We will bubble up the error via the `group_iter`
                key(item_result.as_ref().unwrap()).clone()
            })
        {
            // Split up the groupby into sets of `max_items` items to bound the in-memory size
            // of the Vec<T>
            for (_, piece) in &group_iter
                .map(Result::unwrap)
                .enumerate()
                .group_by(|(i, _)| *i / max_items)
            {
                let items: Vec<T> = piece.map(|(_, item)| item).collect();
                let send_res = send.send((key.clone(), items));
                // stop sending if all the recievers have hung up
                // go down and get the Error/panic from the dead
                // worker
                if send_res.is_err() {
                    break;
                }
            }
        }

        // Close the send channel - this will cause the threads to exit
        drop(send);

        collect_jobs(handles)
    });

    match r {
        Ok(v) => v,
        Err(e) => Err(anyhow!(decipher_panic(e))),
    }
}

fn decipher_panic(p: Box<dyn Any + 'static + Send>) -> String {
    if let Some(&s) = p.downcast_ref::<&'static str>() {
        s.to_string()
    } else if let Ok(s) = p.downcast::<String>() {
        *s
    } else {
        "thread panicked with unrecognized type".to_string()
    }
}

#[cfg(test)]
mod proc_tests {
    use super::*;
    use martian_filetypes::bin_file::BincodeFile;
    use quickcheck::{Arbitrary, Gen};
    use serde::{Deserialize, Serialize};
    use shardio::{Range, ShardReader, ShardWriter, SortKey};
    use std::borrow::Cow;
    use std::marker::PhantomData;

    // ----------------------------------------------------------------------------------
    // ----------------------------------------------------------------------------------
    // Copied from shardio::shard_tests

    #[derive(Copy, Clone, Eq, PartialEq, Serialize, Deserialize, Debug, PartialOrd, Ord, Hash)]
    struct T1 {
        a: u64,
        b: u32,
        c: u16,
        d: u8,
    }

    impl Arbitrary for T1 {
        fn arbitrary(g: &mut Gen) -> T1 {
            T1 {
                a: u64::arbitrary(g),
                b: u32::arbitrary(g),
                c: u16::arbitrary(g),
                d: u8::arbitrary(g),
            }
        }
    }

    struct FieldDSort;
    impl SortKey<T1> for FieldDSort {
        type Key = u8;
        fn sort_key(item: &T1) -> Cow<'_, u8> {
            Cow::Borrowed(&item.d)
        }
    }

    fn rand_items(n: usize) -> Vec<T1> {
        let mut items = Vec::new();

        for i in 0..n {
            let tt = T1 {
                a: (((i / 2) + (i * 10)) % 3 + (i * 5) % 2) as u64,
                b: i as u32,
                c: (i * 2) as u16,
                d: ((i % 5) * 10 + (if i % 10 > 7 { i / 10 } else { 0 })) as u8,
            };
            items.push(tt);
        }

        items
    }

    // End of copy from shardio
    // ----------------------------------------------------------------------------------
    // ----------------------------------------------------------------------------------

    struct ProcIntoVec<K, T>
    where
        BincodeFile<Vec<T>>: LazyFileTypeIO<T>,
    {
        #[allow(clippy::type_complexity)]
        items: Vec<(K, SpillVec<T, BincodeFile<Vec<T>>>)>,
        phantom: PhantomData<K>,
    }

    impl<K, T> ProcIntoVec<K, T>
    where
        BincodeFile<Vec<T>>: LazyFileTypeIO<T>,
    {
        pub fn new() -> ProcIntoVec<K, T> {
            ProcIntoVec {
                items: Vec::new(),
                phantom: PhantomData,
            }
        }
    }

    impl<K, T> Proc for ProcIntoVec<K, T>
    where
        BincodeFile<Vec<T>>: LazyFileTypeIO<T>,
    {
        type Item = (K, SpillVec<T, BincodeFile<Vec<T>>>);
        type Err = anyhow::Error;

        fn process(&mut self, item: Self::Item) -> Result<()> {
            self.items.push(item);
            Ok(())
        }
    }

    #[test]
    fn run_test_proc() -> Result<()> {
        test_proc(12)?;
        test_proc(1000)?;
        test_proc(100_000)
    }

    fn test_proc(n_items: usize) -> Result<()> {
        let disk_chunk_size = 10;
        let producer_chunk_size = 2;
        let buffer_size = 20;

        println!(
            "test round trip: disk_chunk_size: {disk_chunk_size}, producer_chunk_size: {producer_chunk_size}, n_items: {n_items}"
        );

        let tmp = tempfile::NamedTempFile::new()?;

        // Write and close file
        let true_items = {
            let manager: ShardWriter<T1, FieldDSort> = ShardWriter::new(
                tmp.path(),
                producer_chunk_size,
                disk_chunk_size,
                buffer_size,
            )?;
            let mut true_items = rand_items(n_items);

            // Sender must be closed
            {
                for chunk in true_items.chunks(n_items / 8) {
                    let mut sender = manager.get_sender();
                    for item in chunk {
                        sender.send(*item)?;
                    }
                }
            }
            true_items.sort();
            true_items
        };

        // Open finished file
        let reader = ShardReader::<T1, FieldDSort>::open(tmp.path())?;
        let iter = reader.iter_range(&Range::all())?;
        let tmp_dir = tempfile::tempdir()?;
        let processors = vec![ProcIntoVec::new(), ProcIntoVec::new(), ProcIntoVec::new()];
        let r = group_by_processor(iter, processors, |item| item.d, tmp_dir.path())?;

        let mut all_items = Vec::new();

        for p in r {
            for (k, items) in p.items {
                for i in items.iter()? {
                    let i = i?;
                    assert_eq!(k, i.d);
                    all_items.push(i);
                }
            }
        }

        all_items.sort();
        assert_eq!(true_items, all_items);
        Ok(())
    }

    // Make sure we don't deadlock & propagate the error if all threads panic
    // NOTE: this test puts a bunch of scary text in the log because it
    // intentionally creates threads that panic.  Disabling for routine testing.
    /*
    #[derive(Debug)]
    struct ProcPanic;

    impl Proc for ProcPanic {
        type Item = (usize, Vec<usize>);
        type Err = failure::Error;

        fn process(&mut self, _item: (usize, Vec<usize>)) -> Result<()> {
            panic!("ProcPanic panicked");
        }
    }

    #[test]
    fn test_proc_panic() {
        let iter = (0..10000).map(Ok);
        let processors = vec![ProcPanic, ProcPanic, ProcPanic, ProcPanic];
        let r = group_by_processor(iter, processors, |item| item);
        assert!(r.is_err());
        println!("{:?}", r);
    }
    */

    #[derive(Debug)]
    struct ProcErr;

    impl Proc for ProcErr {
        type Item = (usize, SpillVec<usize, BincodeFile<Vec<usize>>>);
        type Err = anyhow::Error;

        fn process(&mut self, _item: Self::Item) -> Result<()> {
            Err(anyhow!("ProcErr errored"))
        }
    }

    // Make sure we don't deadlock & propagate the error if all threads error
    #[test]
    fn test_proc_error() {
        let iter = (0..10000).map(Ok);
        let processors = vec![ProcErr, ProcErr, ProcErr, ProcErr];
        let tmp_dir = tempfile::tempdir().unwrap();
        let r = group_by_processor(iter, processors, |item| *item, tmp_dir.path());
        assert!(r.is_err());
        println!("{r:?}");
    }
}
