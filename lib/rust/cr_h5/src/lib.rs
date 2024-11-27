use anyhow::Result;
use hdf5::{Dataset, Extents, Group, H5Type};

pub mod compare;
pub mod count_matrix;
pub mod feature_reference_io;
pub mod iter;
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
        .chunk((1 << 16,))
        .shape(Extents::resizable((0,).into()))
        .create(name)?;

    Ok(ds)
}

#[allow(clippy::deref_addrof)]
/// Create a HDF5 dataset for one column of the molecule info file, with the appropriate chunking and
/// compression settings.
pub fn write_column_ds<T: H5Type>(group: &Group, name: &str, data: &[T]) -> Result<()> {
    group
        .new_dataset::<T>()
        .shuffle()
        .deflate(1)
        .shape((data.len(),))
        .create(name)?
        .write(data)?;
    Ok(())
}
