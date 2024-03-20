use anyhow::Result;
use martian_filetypes::{LazyFileTypeIO, LazyWrite};
use std::vec;

enum SpillableStorage<M, D> {
    MemOnly(M),
    MemAndDisk(M, D),
}

pub struct SpillVec<T, F>
where
    F: LazyFileTypeIO<T>,
{
    storage: Option<SpillableStorage<Vec<T>, <F as LazyFileTypeIO<T>>::LazyWriter>>,
    spill_file: F,
    max_items: usize,
    num_items: usize,
}

impl<T, F> SpillVec<T, F>
where
    F: LazyFileTypeIO<T>,
{
    pub fn new(max_items: usize, spill_file: F) -> Self {
        let buffer = Vec::with_capacity(16);
        SpillVec {
            storage: Some(SpillableStorage::MemOnly(buffer)),
            spill_file,
            max_items,
            num_items: 0,
        }
    }

    pub fn push(&mut self, item: T) -> Result<()> {
        self.push_weighted(item, 1)
    }

    pub fn push_weighted(&mut self, item: T, weight: usize) -> Result<()> {
        assert!(weight > 0);
        self.num_items += weight;
        // State transition from MemOnly to MemAndDisk if the buffer is full.
        self.storage = match self.storage.take().unwrap() {
            SpillableStorage::MemOnly(buffer) if self.num_items >= self.max_items => {
                let writer = self.spill_file.lazy_writer()?;
                Some(SpillableStorage::MemAndDisk(buffer, writer))
            }
            s => Some(s),
        };

        // Add the item to the mem or disk appropriately
        match self.storage {
            Some(SpillableStorage::MemOnly(ref mut buffer)) => buffer.push(item),
            Some(SpillableStorage::MemAndDisk(_, ref mut writer)) => {
                writer.write_item(&item)?;
            }
            _ => unreachable!(),
        }
        Ok(())
    }

    pub fn iter(self) -> Result<SpillVecReader<T, F>> {
        match self.storage.unwrap() {
            SpillableStorage::MemOnly(buffer) => Ok(SpillVecReader {
                inner: Some(SpillableStorage::MemOnly(buffer.into_iter())),
                spill_file: self.spill_file,
            }),
            SpillableStorage::MemAndDisk(buffer, writer) => {
                writer.finish()?;
                let reader = self.spill_file.lazy_reader()?;
                Ok(SpillVecReader {
                    inner: Some(SpillableStorage::MemAndDisk(buffer.into_iter(), reader)),
                    spill_file: self.spill_file,
                })
            }
        }
    }

    #[cfg(test)]
    fn buffer_len(&self) -> usize {
        match self.storage.as_ref().unwrap() {
            SpillableStorage::MemOnly(ref buf) => buf.len(),
            SpillableStorage::MemAndDisk(ref buf, _) => buf.len(),
        }
    }
}

pub struct SpillVecReader<T, F>
where
    F: LazyFileTypeIO<T>,
{
    inner: Option<SpillableStorage<vec::IntoIter<T>, <F as LazyFileTypeIO<T>>::LazyReader>>,
    spill_file: F,
}

impl<T, F> Iterator for SpillVecReader<T, F>
where
    F: LazyFileTypeIO<T>,
{
    type Item = Result<T>;
    fn next(&mut self) -> Option<Self::Item> {
        let inner = self.inner.take().unwrap();
        match inner {
            SpillableStorage::MemOnly(mut vec_iter) => {
                let next = vec_iter.next().map(Ok);
                self.inner = Some(SpillableStorage::MemOnly(vec_iter));
                next
            }
            SpillableStorage::MemAndDisk(vec_iter, mut reader) => match reader.next() {
                Some(item) => {
                    self.inner = Some(SpillableStorage::MemAndDisk(vec_iter, reader));
                    Some(item)
                }
                None => {
                    self.inner = Some(SpillableStorage::MemOnly(vec_iter));
                    match std::fs::remove_file(&self.spill_file) {
                        Ok(_) => self.next(),
                        Err(e) => Some(Err(e.into())),
                    }
                }
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use itertools::Itertools;
    use martian::MartianFileType;
    use martian_filetypes::bin_file::BincodeFile;
    use proptest::arbitrary::any;
    use proptest::collection::vec;
    use proptest::proptest;
    use std::path::PathBuf;

    fn spill_vec_roundtrip<T: Ord + Clone, F: LazyFileTypeIO<T>>(
        items: Vec<T>,
        capacity: usize,
        spill_file: F,
    ) -> Result<()> {
        let spill_pathbuf = PathBuf::from(spill_file.as_ref());
        let mut spill_vec = SpillVec::new(capacity, spill_file);
        let expected: Vec<_> = items.iter().cloned().sorted().collect();
        for item in items {
            spill_vec.push(item)?;
        }
        assert!(spill_vec.buffer_len() <= capacity);
        assert_eq!(spill_vec.num_items, expected.len());
        if expected.len() > capacity {
            assert!(spill_pathbuf.exists());
        }
        let mut actual: Vec<_> = spill_vec.iter()?.try_collect()?;
        assert!(!spill_pathbuf.exists()); // Should be cleaned up
        actual.sort();
        assert!(expected == actual);
        Ok(())
    }

    proptest! {
        #[test]
        fn prop_test_spill_vec_roundtrip(
            items in vec(any::<u8>(), 0usize..10000usize),
            capacity in 0..1000usize,
        ) {
            let tmp_dir = tempfile::tempdir().unwrap();
            let spill_file = BincodeFile::new(&tmp_dir, "prop_test_spill_vec_roundtrip");
            spill_vec_roundtrip(items, capacity, spill_file).unwrap();
            drop(tmp_dir);
        }
    }
}
