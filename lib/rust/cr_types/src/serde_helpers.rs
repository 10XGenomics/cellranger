use anyhow::Result;
use serde::de::Visitor;
use serde::{de, Deserialize, Deserializer};
use std::convert::TryFrom;
use std::fmt::Formatter;

#[derive(Debug)]
pub enum NumberOrStr {
    Number(u64),
    Str(String),
}

impl TryFrom<NumberOrStr> for usize {
    type Error = de::value::Error;

    fn try_from(val: NumberOrStr) -> Result<usize, Self::Error> {
        match val {
            NumberOrStr::Number(n) => usize::try_from(n).map_err(de::Error::custom),
            NumberOrStr::Str(s) => s.parse().map_err(de::Error::custom),
        }
    }
}

impl TryFrom<NumberOrStr> for u16 {
    type Error = de::value::Error;

    fn try_from(val: NumberOrStr) -> Result<u16, Self::Error> {
        match val {
            NumberOrStr::Number(n) => u16::try_from(n).map_err(de::Error::custom),
            NumberOrStr::Str(s) => s.parse().map_err(de::Error::custom),
        }
    }
}

impl<'de> Deserialize<'de> for NumberOrStr {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        struct NumberOrStrVisitor;
        impl<'de> Visitor<'de> for NumberOrStrVisitor {
            type Value = NumberOrStr;
            fn expecting(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
                f.write_str("A number or string")
            }

            fn visit_u64<E>(self, id: u64) -> Result<Self::Value, E>
            where
                E: de::Error,
            {
                Ok(NumberOrStr::Number(id))
            }

            fn visit_str<E>(self, id: &str) -> Result<Self::Value, E>
            where
                E: de::Error,
            {
                Ok(NumberOrStr::Str(id.into()))
            }
        }
        deserializer.deserialize_any(NumberOrStrVisitor)
    }
}
