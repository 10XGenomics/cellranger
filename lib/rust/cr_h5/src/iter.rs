//! A generalized approach for streaming 1D data from HDF5 using chunked reads.
//!
//! This could be further optimized to pre-fetch the next buffer asynchronously
//! in a worker thread.
//!
//! TODO: refactor MoleculeInfoIterator to use this as the underlying backing.

use anyhow::Result;
use hdf5::{Dataset, H5Type};
use ndarray::Array1;
use std::cmp::min;
use std::sync::mpsc::{sync_channel, Receiver};
use std::thread;

struct H5Buffer<T: H5Type> {
    /// Index of the end of the buffer in the origin dataset.
    end: usize,
    /// Loaded data.
    data: Array1<T>,
}

impl<T: H5Type> H5Buffer<T> {
    fn new(dataset: &Dataset, start: usize, end: usize) -> Result<Self> {
        Ok(Self {
            end,
            data: dataset.read_slice_1d(start..end)?,
        })
    }
}

/// A chunked-loading iterator over a single HDF5 dataset.
pub struct H5Iterator<T: H5Type> {
    buf: Option<H5Buffer<T>>,
    /// Current iteration index into the dataset.
    index: usize,
    /// Total elements in the dataset.
    size: usize,
    /// The size of each chunk to read.
    chunk_size: usize,
    /// Have we returned a read error?
    read_error_returned: bool,
    dataset: Dataset,
}

impl<T: H5Type> H5Iterator<T> {
    /// Create an iterator over the provided dataset, loading chunk_size elements at a time.
    pub fn new(dataset: Dataset, chunk_size: usize) -> Self {
        let size = dataset.size();
        Self {
            buf: None,
            index: 0,
            size,
            chunk_size,
            read_error_returned: false,
            dataset,
        }
    }
}

impl<T: H5Type + Clone> Iterator for H5Iterator<T> {
    type Item = Result<T>;
    fn next(&mut self) -> Option<Self::Item> {
        if self.read_error_returned || self.index >= self.size {
            return None;
        }

        if self.buf.is_none() {
            match H5Buffer::new(
                &self.dataset,
                self.index,
                min(self.index + self.chunk_size, self.size),
            ) {
                Ok(buf) => {
                    self.buf = Some(buf);
                }
                Err(err) => {
                    self.read_error_returned = true;
                    return Some(Err(err));
                }
            }
        }

        let (clear_buf, cur_data) = {
            let buf = self.buf.as_ref().unwrap();
            let cur_index = self.index % self.chunk_size;
            self.index += 1;
            (self.index >= buf.end, buf.data[cur_index].clone())
        };
        if clear_buf {
            self.buf = None;
        }

        Some(Ok(cur_data))
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = if self.index < self.size {
            self.size - self.index
        } else {
            0
        };
        (remaining, Some(remaining))
    }
}

/// Spawn a thread to read chunks asynchronously.
/// Chunks are sent back to the calling thread on a channel bounded to hold no
/// more than max_chunks.
fn read_chunks_async<T: H5Type + Send>(
    dataset: Dataset,
    chunk_size: usize,
    max_chunks: usize,
) -> Receiver<Result<H5Buffer<T>>> {
    let (sender, receiver) = sync_channel(max_chunks);
    thread::spawn(move || {
        let size = dataset.size();
        let mut start = 0;
        loop {
            if start >= size {
                return;
            }
            let end = (start + chunk_size).min(size);
            let read_result = H5Buffer::new(&dataset, start, end);
            if sender.send(read_result).is_err() {
                // Receiver hung up, stop reading.
                return;
            }
            start = end;
        }
    });
    receiver
}

/// A threaded chunked-loading iterator over a single HDF5 1D dataset.
/// Spawns a thread to load chunks in the background, buffering a fixed number.
pub struct ThreadedH5Iterator<T: H5Type + Send> {
    dataset: Dataset,
    chunk_size: usize,
    max_chunks: usize,
    /// Channel receiver from reader thread.
    /// Lazily initialized on the first call to the iterator.
    receiver: Option<Receiver<Result<H5Buffer<T>>>>,
    /// Iterator over the last chunk receieved from the reader.
    cur_iter: Option<<Array1<T> as std::iter::IntoIterator>::IntoIter>,
    /// Current iteration index into the dataset.
    index: usize,
    /// Total elements in the dataset.
    size: usize,
    /// Have we returned a read error?
    read_error_returned: bool,
}

#[allow(unused)]
impl<T: H5Type + Send> ThreadedH5Iterator<T> {
    fn new(dataset: Dataset, chunk_size: usize, max_chunks: usize) -> Self {
        let size = dataset.size();
        Self {
            dataset,
            chunk_size,
            max_chunks,
            receiver: None,
            cur_iter: None,
            index: 0,
            size,
            read_error_returned: false,
        }
    }
}

impl<T: H5Type + Clone + Send> Iterator for ThreadedH5Iterator<T> {
    type Item = Result<T>;
    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if self.read_error_returned {
                return None;
            }

            // Lazily initialize receiver thread.
            if self.receiver.is_none() {
                self.receiver = Some(read_chunks_async(
                    self.dataset.clone(),
                    self.chunk_size,
                    self.max_chunks,
                ));
            }

            if self.cur_iter.is_none() {
                match self.receiver.as_mut().unwrap().recv() {
                    Ok(Ok(buf)) => {
                        self.cur_iter = Some(buf.data.into_iter());
                    }
                    Ok(Err(err)) => {
                        self.read_error_returned = true;
                        return Some(Err(err));
                    }
                    Err(_) => {
                        // Sending thread terminated - no more data to read.
                        return None;
                    }
                }
            }

            match self.cur_iter.as_mut().unwrap().next() {
                Some(item) => {
                    self.index += 1;
                    return Some(Ok(item));
                }
                None => {
                    // We've exhausted this iterator; load another chunk.
                    self.cur_iter = None;
                }
            }
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = if self.index < self.size {
            self.size - self.index
        } else {
            0
        };
        (remaining, Some(remaining))
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use hdf5::File;
    use std::path::Path;

    #[test]
    fn test_h5_iter() -> Result<()> {
        let path = Path::new("test/h5/pbmc_1k_v2_molecule_info.h5");

        let file = File::open(path)?;

        let items: Vec<_> =
            H5Iterator::<u64>::new(file.dataset("barcode_idx")?, 1024).collect::<Result<_>>()?;

        assert_eq!(items.len(), 4303635);

        // Ensure chunked loading produces same results as loading entire collection.
        assert_eq!(items, file.dataset("barcode_idx")?.read_raw::<u64>()?);

        Ok(())
    }

    #[test]
    fn test_threaded_h5_iter() -> Result<()> {
        let path = Path::new("test/h5/pbmc_1k_v2_molecule_info.h5");

        let file = File::open(path)?;

        let items: Vec<_> = ThreadedH5Iterator::<u64>::new(file.dataset("barcode_idx")?, 1024, 1)
            .collect::<Result<_>>()?;

        assert_eq!(items.len(), 4303635);

        // Ensure chunked loading produces same results as loading entire collection.
        assert_eq!(items, file.dataset("barcode_idx")?.read_raw::<u64>()?);

        Ok(())
    }
}
