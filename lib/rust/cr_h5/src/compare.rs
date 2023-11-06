use anyhow::{Context, Result};
use hdf5::types::{FixedAscii, FloatSize, IntSize, TypeDescriptor, VarLenAscii, VarLenUnicode};
use hdf5::{Attribute, Container, Dataset, File, Group, H5Type, Location};
use itertools::{EitherOrBoth, Itertools};
use std::cmp::Ord;
use std::path::Path;

type Cmp<T> = EitherOrBoth<T, T>;

/// Visitor for comparing a pair of HDF5 files
pub trait CompareVisitor: Sized {
    /// entry point -- call this to compare two files
    fn visit_files(&mut self, f1: &File, f2: &File) -> Result<()> {
        walk_group(self, "", Cmp::Both(f1, f2))
    }

    /// Visit a group
    fn visit_group(&mut self, name: &str, values: Cmp<&Group>) -> Result<()> {
        walk_group(self, name, values)
    }

    /// Visit an attribute & compare the values
    fn visit_attribute(&mut self, name: &str, values: Cmp<&Attribute>) -> Result<()> {
        walk_container(
            self,
            name,
            values.map_any(|l| l as &Container, |r| r as &Container),
        )
    }

    /// Visit a dataset & compare the values
    fn visit_dataset(&mut self, name: &str, values: Cmp<&Dataset>) -> Result<()> {
        walk_container(
            self,
            name,
            values.map_any(|l| l as &Container, |r| r as &Container),
        )
    }

    /// Test whether to skip the comparison of given item
    fn skip_comparison(&self, _name: &str) -> bool {
        false
    }

    /// Customize how to report on group vs dataset mismatch
    fn report_conflict(&mut self, name: &str);

    /// Customize how to report on a group or dataset that is missing
    /// from one file
    fn report_missing(&mut self, name: &str);

    /// Customize how to report on a difference between datasets or attributes
    fn report_diff(&mut self, name: &str);
}

/// H5 file comparison tool. Can be customized with a list of
/// h5 objects to ignore, for example "/group/dataset", "/group/attribute",
// or "/group/dataset/attribute" will all be recognized.
// Implement CompareVisitor to walk over h5 objects.
pub struct H5Compare {
    ignore_items: Vec<String>,
    correct: bool,
    errors: Vec<String>,
}

impl H5Compare {
    pub fn new(ignore_items: Vec<String>) -> H5Compare {
        H5Compare {
            ignore_items,
            correct: true,
            errors: Vec::new(),
        }
    }

    /// Compare the contents of two hdf5 files and report detected difference.
    /// Ignore items listed in `ignore_items`.
    pub fn compare_files(
        f1: &Path,
        f2: &Path,
        ignore_items: Vec<String>,
    ) -> Result<(bool, Vec<String>)> {
        let f1 = File::open(f1).with_context(|| f1.display().to_string())?;
        let f2 = File::open(f2).with_context(|| f2.display().to_string())?;
        let mut cmp = H5Compare::new(ignore_items);
        cmp.visit_files(&f1, &f2)?;
        Ok((cmp.correct, cmp.errors))
    }
}

impl CompareVisitor for H5Compare {
    fn skip_comparison(&self, name: &str) -> bool {
        self.ignore_items.contains(&name.to_string())
    }

    fn report_missing(&mut self, name: &str) {
        self.correct = false;
        let m = format!("missing item: {name}");
        self.errors.push(m);
    }

    fn report_diff(&mut self, name: &str) {
        self.correct = false;
        let m = format!("difference detected in item: {name}");
        self.errors.push(m);
    }

    fn report_conflict(&mut self, name: &str) {
        self.correct = false;
        let m = format!("conflicting type in item: {name}");
        self.errors.push(m);
    }
}

/// utility for disentangling groups and datasets
/// when iterating over members of a group.
#[derive(Debug)]
enum GroupOrDataset {
    Group(Group),
    Dataset(Dataset),
}

impl GroupOrDataset {
    fn load(obj: &Group, name: &str) -> Result<Self> {
        match obj.group(name) {
            Ok(g) => Ok(Self::Group(g)),
            Err(_) => Ok(Self::Dataset(obj.dataset(name)?)),
        }
    }
}

/// The default traversal from the root of the barcode annotation graph that visits all the
/// read nodes.
pub fn walk_group(
    visitor: &mut impl CompareVisitor,
    name: &str,
    values: Cmp<&Group>,
) -> Result<()> {
    let (left, right) = match values {
        Cmp::Left(_) | Cmp::Right(_) => {
            visitor.report_missing(name);
            return Ok(());
        }
        Cmp::Both(l, r) => (l, r),
    };

    // compare attributes
    walk_attributes(visitor, name, left, right)?;

    let mut left_members = left.member_names()?;
    left_members.sort();

    let mut right_members = right.member_names()?;
    right_members.sort();

    let merged_members = left_members.iter().merge_join_by(&right_members, Ord::cmp);

    for m in merged_members {
        let (v, name_suffix) = match m {
            Cmp::Left(n) => (Cmp::Left(GroupOrDataset::load(left, n)?), n),
            Cmp::Right(n) => (Cmp::Right(GroupOrDataset::load(right, n)?), n),
            Cmp::Both(n, _) => (
                Cmp::Both(
                    GroupOrDataset::load(left, n)?,
                    GroupOrDataset::load(right, n)?,
                ),
                n,
            ),
        };

        let name = &format!("{name}/{name_suffix}");
        match v {
            Cmp::Left(GroupOrDataset::Group(l)) => {
                visitor.visit_group(name, Cmp::Left(&l))?;
            }
            Cmp::Right(GroupOrDataset::Group(r)) => {
                visitor.visit_group(name, Cmp::Right(&r))?;
            }
            Cmp::Left(GroupOrDataset::Dataset(l)) => {
                visitor.visit_dataset(name, Cmp::Left(&l))?;
            }
            Cmp::Right(GroupOrDataset::Dataset(r)) => {
                visitor.visit_dataset(name, Cmp::Right(&r))?;
            }
            Cmp::Both(GroupOrDataset::Group(l), GroupOrDataset::Group(r)) => {
                visitor.visit_group(name, Cmp::Both(&l, &r))?;
            }
            Cmp::Both(GroupOrDataset::Dataset(l), GroupOrDataset::Dataset(r)) => {
                walk_attributes(visitor, name, &l, &r)?;
                visitor.visit_dataset(name, Cmp::Both(&l, &r))?;
            }
            _ => visitor.report_conflict(name),
        }
    }

    Ok(())
}

/// Compare a pair of HDF5 datasets or attributes.
pub fn walk_container(
    visitor: &mut impl CompareVisitor,
    name: &str,
    values: Cmp<&Container>,
) -> Result<()> {
    if visitor.skip_comparison(name) {
        return Ok(());
    }

    match values {
        Cmp::Left(_) | Cmp::Right(_) => visitor.report_missing(name),
        Cmp::Both(l, r) => {
            let td = l.dtype().unwrap().to_descriptor()?;
            match dispatch_h5_type(&td, DataCompare, (l, r)) {
                Some(Ok(false)) => {
                    visitor.report_diff(name);
                }
                Some(Err(e)) => {
                    println!("couldn't compare: {name} {e}");
                    visitor.report_diff(name);
                }
                _ => (),
            }
        }
    }

    Ok(())
}

/// The default traversal from the root of the barcode annotation graph that visits all the
/// read nodes.
pub fn walk_attributes(
    visitor: &mut impl CompareVisitor,
    prefix: &str,
    left: &Location,
    right: &Location,
) -> Result<()> {
    let mut left_members = left.attr_names()?;
    left_members.sort();

    let mut right_members = right.attr_names()?;
    right_members.sort();

    let merged_members = left_members.iter().merge_join_by(&right_members, Ord::cmp);

    for m in merged_members {
        let (v, name) = match m {
            Cmp::Left(n) => (Cmp::Left(left.attr(n)?), n),
            Cmp::Right(n) => (Cmp::Right(right.attr(n)?), n),
            Cmp::Both(n, _) => (Cmp::Both(left.attr(n)?, right.attr(n)?), n),
        };
        visitor.visit_attribute(&format!("{prefix}/{name}"), v.as_ref())?;
    }

    Ok(())
}

pub trait Action<T: H5Type, A, R> {
    fn action(&self, args: A) -> R;
}

pub struct DataCompare;

/// Load a dataset from each container with type T and compare the results.
impl<T: H5Type + PartialEq> Action<T, (&Container, &Container), Result<bool>> for DataCompare {
    fn action(&self, d: (&Container, &Container)) -> Result<bool> {
        let ds1 = d.0.as_reader().read::<T, ndarray::IxDyn>()?;
        let ds2 = d.1.as_reader().read::<T, ndarray::IxDyn>()?;
        Ok(ds1 == ds2)
    }
}

pub fn dispatch_h5_type<T, A, R>(td: &TypeDescriptor, action: T, args: A) -> Option<R>
where
    T: Action<u8, A, R>
        + Action<u16, A, R>
        + Action<u32, A, R>
        + Action<u64, A, R>
        + Action<i8, A, R>
        + Action<i16, A, R>
        + Action<i32, A, R>
        + Action<i64, A, R>
        + Action<f32, A, R>
        + Action<f64, A, R>
        + Action<VarLenAscii, A, R>
        + Action<VarLenUnicode, A, R>
        + Action<FixedAscii<256>, A, R>,
{
    match td {
        TypeDescriptor::Unsigned(IntSize::U1) => Some(Action::<u8, _, _>::action(&action, args)),
        TypeDescriptor::Unsigned(IntSize::U2) => Some(Action::<u16, _, _>::action(&action, args)),
        TypeDescriptor::Unsigned(IntSize::U4) => Some(Action::<u32, _, _>::action(&action, args)),
        TypeDescriptor::Unsigned(IntSize::U8) => Some(Action::<u64, _, _>::action(&action, args)),
        TypeDescriptor::Integer(IntSize::U1) => Some(Action::<i8, _, _>::action(&action, args)),
        TypeDescriptor::Integer(IntSize::U2) => Some(Action::<i16, _, _>::action(&action, args)),
        TypeDescriptor::Integer(IntSize::U4) => Some(Action::<i32, _, _>::action(&action, args)),
        TypeDescriptor::Integer(IntSize::U8) => Some(Action::<i64, _, _>::action(&action, args)),
        TypeDescriptor::Float(FloatSize::U4) => Some(Action::<f32, _, _>::action(&action, args)),
        TypeDescriptor::Float(FloatSize::U8) => Some(Action::<f64, _, _>::action(&action, args)),
        TypeDescriptor::VarLenAscii => Some(Action::<VarLenAscii, _, _>::action(&action, args)),
        TypeDescriptor::VarLenUnicode => Some(Action::<VarLenUnicode, _, _>::action(&action, args)),
        TypeDescriptor::FixedAscii(n) if *n <= 256 => {
            Some(Action::<FixedAscii<256>, _, _>::action(&action, args))
        }
        _ => None,
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_hdf5_compare() -> Result<()> {
        hdf5::silence_errors(true);
        let f1 = File::open("test/h5/diff_test_fbm_a.h5")?;
        let f2 = File::open("test/h5/diff_test_fbm_b.h5")?;
        let cmp = {
            // H5 still dumps a bunch of annoying messages. Gag them if possible
            let _err_gag = gag::Gag::stderr().ok();
            let mut cmp = H5Compare::new(Vec::new());
            cmp.visit_files(&f1, &f2)?;
            cmp
        };

        assert!(!cmp.correct);
        assert_eq!(
            cmp.errors,
            vec!["difference detected in item: /software_version".to_string()]
        );

        println!("{:#?}", cmp.errors);

        Ok(())
    }

    #[test]
    fn test_hdf5_compare_many_diffs() -> Result<()> {
        hdf5::silence_errors(true);
        let f1 = File::open("test/h5/diff_test_fbm_a.h5")?;
        let f2 = File::open("test/h5/diff_test_fbm_c.h5")?;
        let cmp = {
            let _err_gag = gag::Gag::stderr().ok();
            let mut cmp = H5Compare::new(Vec::new());
            cmp.visit_files(&f1, &f2)?;
            cmp
        };

        assert!(!cmp.correct);
        assert_eq!(
            cmp.errors,
            vec![
                "difference detected in item: /chemistry_description".to_string(),
                "difference detected in item: /library_ids".to_string(),
                "difference detected in item: /matrix/data".to_string(),
                "difference detected in item: /matrix/indices".to_string(),
                "difference detected in item: /matrix/indptr".to_string()
            ]
        );

        println!("{:#?}", cmp.errors);

        Ok(())
    }
}
