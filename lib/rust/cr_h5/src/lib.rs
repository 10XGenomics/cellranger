//! cr_h5
#![deny(missing_docs)]
use anyhow::{Result, ensure};
use hdf5::{Dataset, Extents, Group, H5Type};

mod compare;
pub mod count_matrix;
pub mod feature_reference_io;
mod iter;
pub mod molecule_info;
pub mod probe_reference_io;

/// Write a scalar attribute to the given group, with the given attribute name.
pub fn scalar_attribute<T: H5Type>(group: &Group, name: &str, value: T) -> Result<()> {
    let attr = group.new_attr::<T>().shape(()).create(name)?;

    let v = ndarray::Array0::from_shape_vec((), vec![value])?;
    attr.as_writer().write(&v)?;
    Ok(())
}

/// Extend the dataset `ds` with the data in slice `data`. `ds` must have type `T` and be
/// 1-dimensional.
pub fn extend_dataset<T: H5Type>(ds: &Dataset, data: &[T]) -> Result<()> {
    let sz = ds.shape()[0];
    let new_sz = sz + data.len();

    ds.resize([new_sz])?;
    ds.as_writer().write_slice(data, sz..new_sz)?;

    Ok(())
}

/// Create a HDF5 dataset for one column of the molecule info file, with the appropriate chunking and
/// compression settings.
pub fn make_column_ds<T: H5Type>(group: &Group, name: &str) -> Result<Dataset> {
    let ds = group
        .new_dataset::<T>()
        .shuffle()
        .deflate(1)
        .chunk(1 << 16)
        .shape(Extents::resizable((0).into()))
        .create(name)?;

    Ok(ds)
}

/// Create a HDF5 dataset for one column of the molecule info file, with the appropriate chunking and
/// compression settings.
pub fn write_column_ds<T: H5Type>(group: &Group, name: &str, data: &[T]) -> Result<()> {
    let dataset = group.new_dataset::<T>().shape(data.len());
    if data.is_empty() {
        // Filters are not permitted for an empty dataset.
        dataset.create(name)?;
    } else {
        dataset.shuffle().deflate(1).create(name)?.write(data)?;
    }
    Ok(())
}

/// Action to take when creating a new column.
#[derive(Clone, Copy, Debug)]
pub enum ColumnAction {
    /// Create a new column, ensuring it does not already exist.
    CreateNew,
    /// Replace an existing column. If the column does not exist, it is created.
    ReplaceExisting,
}

/// Handle buffered writing to a HDF5 dataset.
pub(crate) struct ChunkedWriter<T: H5Type> {
    dataset: Dataset,
    buf: Vec<T>,
    buf_size: usize,
}

impl<T: H5Type> ChunkedWriter<T> {
    /// Use chunked writing to write out an entire iterator of items.
    ///
    /// See the `in_group` constructor for parameter reference.
    pub fn write_all(
        group: &Group,
        name: &str,
        buf_size: usize,
        action: ColumnAction,
        items: impl Iterator<Item = T>,
    ) -> Result<()> {
        let mut writer = Self::in_group(group, name, buf_size, action)?;
        for item in items {
            writer.write(item)?;
        }
        writer.flush()?;
        Ok(())
    }

    /// Initialize a chunked writer in the given group with column name.
    ///
    /// Use the provided buffer size.
    /// Returns an error if the dataset already exists in the group and `action` is `CreateNew`.
    pub fn in_group(
        group: &Group,
        name: &str,
        buf_size: usize,
        action: ColumnAction,
    ) -> Result<Self> {
        match action {
            ColumnAction::CreateNew => {
                ensure!(
                    !group.link_exists(name),
                    "Dataset {} already exists in group {}",
                    name,
                    group.name()
                );
            }
            ColumnAction::ReplaceExisting => {
                if group.link_exists(name) {
                    group.unlink(name)?;
                }
            }
        }
        Ok(Self::new(make_column_ds::<T>(group, name)?, buf_size))
    }

    /// Initialize a chunked writer for the provided dataset.
    ///
    /// Use the provided buffer size.
    pub fn new(dataset: Dataset, buf_size: usize) -> Self {
        Self {
            dataset,
            buf_size,
            buf: Vec::with_capacity(buf_size),
        }
    }

    /// Write an item.
    ///
    /// Flushes to disk if the buffer is full.
    pub fn write(&mut self, item: T) -> Result<()> {
        self.buf.push(item);
        if self.buf.len() >= self.buf_size {
            self.flush()?;
        }
        Ok(())
    }

    /// Flush all current buffer contents to disk.
    ///
    /// This method should be called manually after all items have been written.
    pub fn flush(&mut self) -> Result<()> {
        if self.buf.is_empty() {
            return Ok(());
        }
        extend_dataset(&self.dataset, &self.buf)?;
        self.buf.clear();
        Ok(())
    }
}

impl<T: H5Type> Drop for ChunkedWriter<T> {
    fn drop(&mut self) {
        // Follow the pattern of BufWriter. IO errors are lost if they occur
        let _ = self.flush();
    }
}
