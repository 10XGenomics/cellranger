//! A stack-allocated immutable string.

use anyhow::{bail, Result};
use metric::AsMetricPrefix;
use serde::de::{self, Visitor};
use serde::{Deserialize, Deserializer, Serialize};
use std::borrow::Borrow;
use std::fmt::{Debug, Display};
use std::hash::Hash;
use std::ops::Deref;

/// A stack-allocated string occupying the same space as a 64-bit machine word.
pub type ShortString7 = ShortString<7>;

/// A stack-allocated string occupying the same space as two 64-bit machine words.
pub type ShortString15 = ShortString<15>;

/// A stack-allocated string occupying the same space as String.
pub type ShortString23 = ShortString<23>;

/// A stack-allocated immutable string.
/// The contents are guaranteed to be utf-8 as the only way to construct this
/// type is from an instance of &str, either via TryFrom or indirectly via
/// deserialization.
#[derive(Clone, Copy, PartialEq, Eq)]
pub struct ShortString<const N: usize> {
    contents: [u8; N],
    len: u8,
}

impl<const N: usize> ShortString<N> {
    pub fn as_str(&self) -> &str {
        unsafe { std::str::from_utf8_unchecked(&self.contents[0..self.len as usize]) }
    }

    /// Attempts to pack a str into this size of short string.
    /// Panics if the input string is too long.
    pub fn pack(s: &str) -> Self {
        Self::try_from(s).unwrap()
    }

    /// Return the length of this string, in bytes.
    pub fn len(&self) -> usize {
        self.len as usize
    }

    /// Return true if this string is empty.
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }
}

impl<const N: usize> Default for ShortString<N> {
    fn default() -> Self {
        Self {
            contents: [0; N],
            len: 0,
        }
    }
}

impl<const N: usize> Hash for ShortString<N> {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.as_str().hash(state);
    }
}

impl<const N: usize> PartialOrd for ShortString<N> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl<const N: usize> Ord for ShortString<N> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.as_str().cmp(other.as_str())
    }
}

impl<const N: usize> TryFrom<&str> for ShortString<N> {
    type Error = anyhow::Error;

    fn try_from(value: &str) -> Result<Self> {
        if value.len() > N {
            // if the input string is very long, truncate it.
            let (truncated_input, more) = if value.len() > 40 {
                (&value[..37], "...")
            } else {
                (value, "")
            };
            bail!(
                "the string \"{truncated_input}{more}\" (length {}) is too long to represent as a ShortString<{N}>",
                value.len()
            );
        }
        let mut contents = [0; N];
        contents[..value.len()].copy_from_slice(value.as_bytes());
        Ok(Self {
            contents,
            len: value.len() as u8,
        })
    }
}

impl<const N: usize> Deref for ShortString<N> {
    type Target = str;

    fn deref(&self) -> &Self::Target {
        self.as_str()
    }
}

impl<const N: usize> AsRef<str> for ShortString<N> {
    fn as_ref(&self) -> &str {
        self.as_str()
    }
}

impl<const N: usize> Borrow<str> for ShortString<N> {
    fn borrow(&self) -> &str {
        self.as_str()
    }
}

impl<const N: usize> Debug for ShortString<N> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        Debug::fmt(self.as_str(), f)
    }
}

impl<const N: usize> Display for ShortString<N> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        Display::fmt(self.as_str(), f)
    }
}

impl<const N: usize> AsMetricPrefix for ShortString<N> {
    fn as_metric_prefix(&self) -> Option<&str> {
        Some(self.as_str())
    }
}

// serde implementations

impl<const N: usize> Serialize for ShortString<N> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        self.as_str().serialize(serializer)
    }
}

impl<'de, const N: usize> Deserialize<'de> for ShortString<N> {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        deserializer.deserialize_str(ShortStringVisitor)
    }
}

struct ShortStringVisitor<const N: usize>;

impl<'de, const N: usize> Visitor<'de> for ShortStringVisitor<N> {
    type Value = ShortString<N>;

    fn expecting(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(formatter, "a string shorter than {N} bytes")
    }

    fn visit_str<E>(self, value: &str) -> Result<Self::Value, E>
    where
        E: de::Error,
    {
        ShortString::try_from(value).map_err(E::custom)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use std::mem::size_of;
    const TEST_SIZE: usize = 7;

    #[test]
    fn test_too_long() {
        type SS = ShortString<TEST_SIZE>;

        let too_long = "foobarba";
        assert_eq!(TEST_SIZE + 1, too_long.as_bytes().len());
        assert!(SS::try_from(&too_long[..TEST_SIZE]).is_ok());
        assert!(SS::try_from(too_long).is_err());

        let dangerous_snakes = "🐍🐍";
        assert_eq!(TEST_SIZE + 1, dangerous_snakes.len());
        assert!(SS::try_from(dangerous_snakes).is_err());
    }

    #[test]
    fn test_alignment() {
        assert_eq!(size_of::<u64>(), size_of::<ShortString7>());
        assert_eq!(size_of::<u128>(), size_of::<ShortString15>());
        assert_eq!(size_of::<String>(), size_of::<ShortString23>());
    }
}
