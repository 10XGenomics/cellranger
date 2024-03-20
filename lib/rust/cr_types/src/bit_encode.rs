use serde::{Deserialize, Serialize};
use std::marker::PhantomData;

// Gives a `bit_pos`, check that it is in the range 0..8
// and return 1 << bit_pos
fn bit_flag(bit_pos: u8) -> u8 {
    assert!(
        bit_pos < 8,
        "Bit position {bit_pos} cannot be 8-bit encoded"
    );
    1u8 << bit_pos
}

/// A single byte to store a set of items (upto 8) of type `T` with one bit per item
/// The type `T` should itself be something that can be converted to u8 and created from
/// a u8. This is typically used to store multiple variants of an enum in a compact manner
#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialOrd, Ord, PartialEq, Eq, Hash)]
#[serde(transparent)]
pub struct BitEncoded<T: From<u8> + Into<u8>> {
    #[serde(skip)]
    phantom: PhantomData<T>,
    bits: u8,
}

impl<T: From<u8> + Into<u8>> BitEncoded<T> {
    pub fn new() -> Self {
        BitEncoded {
            phantom: PhantomData,
            bits: 0,
        }
    }
    pub fn push(&mut self, val: T) {
        self.bits |= bit_flag(val.into());
    }
    #[cfg(test)]
    #[allow(unused)]
    pub fn inspect_bits(self) -> u8 {
        self.bits
    }
    pub fn iter(self) -> BitEncodedIter<T> {
        BitEncodedIter {
            phantom: PhantomData,
            pos: 0,
            bits: self.bits,
        }
    }
}

impl<T: From<u8> + Into<u8>> Default for BitEncoded<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: From<u8> + Into<u8>> From<T> for BitEncoded<T> {
    fn from(val: T) -> Self {
        BitEncoded {
            phantom: PhantomData,
            bits: bit_flag(val.into()),
        }
    }
}

impl<T: From<u8> + Into<u8>> From<Vec<T>> for BitEncoded<T> {
    fn from(vals: Vec<T>) -> Self {
        let mut encoded = BitEncoded::new();
        for val in vals {
            encoded.push(val);
        }
        encoded
    }
}

impl<T: From<u8> + Into<u8>> IntoIterator for BitEncoded<T> {
    type IntoIter = BitEncodedIter<T>;
    type Item = T;

    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

#[derive(Debug)]
pub struct BitEncodedIter<T: From<u8> + Into<u8>> {
    phantom: PhantomData<T>,
    pos: u8,
    bits: u8,
}

impl<T: From<u8> + Into<u8>> Iterator for BitEncodedIter<T> {
    type Item = T;
    fn next(&mut self) -> Option<T> {
        if self.pos > 7 {
            return None;
        }
        let flag = bit_flag(self.pos);
        self.pos += 1;
        if (self.bits & flag) > 0 {
            return Some(T::from(self.pos - 1));
        }
        self.next()
    }
}
