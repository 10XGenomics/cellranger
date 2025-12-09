#![deny(missing_docs)]

use serde::de::{self, Visitor};
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use std::borrow::Borrow;
use std::fmt;
use std::hash::{Hash, Hasher};
use std::iter::FromIterator;
use std::marker::PhantomData;
use std::ops::{Index, IndexMut};

pub trait ArrayContent {
    fn validate_bytes(bytes: &[u8]);
    fn expected_contents() -> &'static str;
}

/// Fixed-sized container for a short DNA sequence or quality.
/// The capacity is determined by the type `N` and the contents are validates based on type `T`
/// Typically used as a convenient container for barcode or UMI sequences or quality.
#[derive(Clone, Copy, PartialOrd, Ord, Eq)]
pub struct ByteArray<T, const N: usize>
where
    T: ArrayContent,
{
    bytes: [u8; N],
    length: u8,
    phantom: PhantomData<T>,
}

impl<T, const N: usize> ByteArray<T, N>
where
    T: ArrayContent,
{
    pub fn new() -> Self {
        ByteArray {
            length: 0,
            bytes: [0; N],
            phantom: PhantomData,
        }
    }

    /// Caller needs to ensure that the bytes are valid
    pub fn push_unchecked(&mut self, src: &[u8]) {
        let len = self.length as usize;
        assert!(
            src.len() <= (N - len),
            "Input slice has length {} which exceeds the remaining capacity of {} bytes in the ByteArray",
            src.len(),
            N - len
        );
        self.bytes[len..len + src.len()].copy_from_slice(src);
        self.length += src.len() as u8;
    }

    pub fn push(&mut self, src: &[u8]) {
        T::validate_bytes(src);
        self.push_unchecked(src);
    }

    /// Remove the element at pos and return it
    pub fn remove(&mut self, pos: usize) -> u8 {
        let val = self[pos];
        self.as_mut_bytes().split_at_mut(pos).1.rotate_left(1);
        self.length -= 1;
        self.bytes[self.length as usize] = 0;
        val
    }

    ///Insert the element at pos without validating
    ///
    /// Panics if the capacity is full or if the pos is > self.len()
    pub fn insert_unchecked(&mut self, pos: usize, byte: u8) {
        let curr_len = self.len();
        assert!(
            curr_len <= (N - 1),
            "No remaining capacity to insert a byte into ByteArray (N={N}).",
        );
        assert!(
            pos <= curr_len,
            "Cannot insert at position {pos} that is more that the array length {curr_len}.",
        );
        self.length += 1;
        self.as_mut_bytes().split_at_mut(pos).1.rotate_right(1);
        self[pos] = byte;
    }

    /// Create a new ByteArray from the given byte slice
    /// The byte slice should contain only valid alphabets as defined by ArrayContent trait
    /// otherwise this function will panic
    pub fn from_bytes(src: &[u8]) -> Self {
        let mut arr = Self::new();
        arr.push(src);
        arr
    }

    /// Create a new ByteArray from the given byte slice
    /// Caller needs to ensure that the byte slice contains only valid alphabets as defined by ArrayContent trait
    pub fn from_bytes_unchecked(src: &[u8]) -> Self {
        let mut arr = Self::new();
        arr.push_unchecked(src);
        arr
    }

    pub fn from_iter_unchecked<C, D>(src: D) -> Self
    where
        C: Borrow<u8>,
        D: IntoIterator<Item = C>,
    {
        let mut src = src.into_iter().fuse();
        let mut bytes = [0; N];
        let mut len = 0;
        for (l, r) in bytes.iter_mut().zip(&mut src) {
            *l = *r.borrow();
            len += 1;
        }
        assert!(
            src.next().is_none(),
            "Error: Input iter exceeds capacity of {} bytes.",
            bytes.len()
        );

        ByteArray {
            length: len,
            bytes,
            phantom: PhantomData,
        }
    }

    /// Returns a byte slice of the contents.
    pub fn as_bytes(&self) -> &[u8] {
        &self.bytes[0..self.length as usize]
    }

    /// Returns a str of the contents.
    pub fn as_str(&self) -> &str {
        std::str::from_utf8(self.as_bytes()).unwrap()
    }

    /// Returns a mutable byte slice of the contents.
    pub fn as_mut_bytes(&mut self) -> &mut [u8] {
        &mut self.bytes[0..self.length as usize]
    }

    /// Returns the length of this sequence, in bytes.
    pub fn len(&self) -> usize {
        self.length as usize
    }

    /// Returns true if self has a length of zero bytes.
    pub fn is_empty(&self) -> bool {
        self.length == 0
    }

    /// Returns an iterator over the bytes.
    pub fn iter(&self) -> std::slice::Iter<'_, u8> {
        self.as_bytes().iter()
    }
}

impl<T, const N: usize> Default for ByteArray<T, N>
where
    T: ArrayContent,
{
    fn default() -> Self {
        Self::new()
    }
}

impl<T, C, const N: usize> FromIterator<C> for ByteArray<T, N>
where
    T: ArrayContent,
    C: Borrow<u8>,
{
    fn from_iter<D: IntoIterator<Item = C>>(src: D) -> Self
    where
        C: Borrow<u8>,
    {
        let array = ByteArray::from_iter_unchecked(src);
        T::validate_bytes(array.as_bytes());
        array
    }
}

impl<T, const N: usize> fmt::Display for ByteArray<T, N>
where
    T: ArrayContent,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.as_str())
    }
}

impl<T, const N: usize> fmt::Debug for ByteArray<T, N>
where
    T: ArrayContent,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Display::fmt(&self, f)
    }
}

impl<T, const N: usize> Index<usize> for ByteArray<T, N>
where
    T: ArrayContent,
{
    type Output = u8;

    fn index(&self, index: usize) -> &Self::Output {
        assert!(index < self.length as usize, "index out of bounds");

        &self.bytes[index]
    }
}

impl<T, const N: usize> IndexMut<usize> for ByteArray<T, N>
where
    T: ArrayContent,
{
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        assert!(index < self.length as usize, "index out of bounds");
        &mut self.bytes[index]
    }
}

impl<T, const N: usize> AsRef<[u8]> for ByteArray<T, N>
where
    T: ArrayContent,
{
    fn as_ref(&self) -> &[u8] {
        self.as_bytes()
    }
}

impl<T, const N: usize> From<ByteArray<T, N>> for String
where
    T: ArrayContent,
{
    fn from(v: ByteArray<T, N>) -> String {
        String::from(v.as_str())
    }
}

impl<T, const N: usize> Borrow<[u8]> for ByteArray<T, N>
where
    T: ArrayContent,
{
    fn borrow(&self) -> &[u8] {
        self.as_bytes()
    }
}

impl<T, const N: usize> Hash for ByteArray<T, N>
where
    T: ArrayContent,
{
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.as_bytes().hash(state);
    }
}

impl<T, const N: usize> PartialEq for ByteArray<T, N>
where
    T: ArrayContent,
{
    fn eq(&self, other: &Self) -> bool {
        self.as_bytes() == other.as_bytes()
    }
}

impl<T, const N: usize> IntoIterator for &ByteArray<T, N>
where
    T: ArrayContent,
{
    type Item = u8;
    type IntoIter = std::iter::Take<std::array::IntoIter<u8, N>>;

    fn into_iter(self) -> Self::IntoIter {
        self.bytes.into_iter().take(self.length as usize)
    }
}

impl<T, const N: usize> IntoIterator for ByteArray<T, N>
where
    T: ArrayContent,
{
    type Item = u8;
    type IntoIter = std::iter::Take<std::array::IntoIter<u8, N>>;

    fn into_iter(self) -> Self::IntoIter {
        self.bytes.into_iter().take(self.length as usize)
    }
}

impl<T, const N: usize> Serialize for ByteArray<T, N>
where
    T: ArrayContent,
{
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        serializer.serialize_str(self.as_str())
    }
}

impl<'de, T, const N: usize> Deserialize<'de> for ByteArray<T, N>
where
    T: ArrayContent,
{
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        deserializer.deserialize_str(ByteArrayVisitor {
            phantom_t: PhantomData,
        })
    }
}

struct ByteArrayVisitor<T, const N: usize> {
    phantom_t: PhantomData<[T; N]>,
}

impl<T, const N: usize> Visitor<'_> for ByteArrayVisitor<T, N>
where
    T: ArrayContent,
{
    type Value = ByteArray<T, N>;

    fn expecting(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(T::expected_contents())
    }

    fn visit_str<E>(self, value: &str) -> Result<Self::Value, E>
    where
        E: de::Error,
    {
        Ok(ByteArray::from_bytes(value.as_bytes()))
    }
}
