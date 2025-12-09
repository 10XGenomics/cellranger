#![deny(missing_docs)]

use crate::{FloatAsInt, Percent};
use anyhow::{Result, anyhow, bail};
use serde_json::{Number, Value};

pub trait TryFromValue {
    fn try_from_value(value: &Value) -> Result<Option<Self>>
    where
        Self: Sized;
}

impl TryFromValue for Number {
    fn try_from_value(value: &Value) -> Result<Option<Self>> {
        match value {
            Value::Number(n) => Ok(Some(n.clone())),
            Value::String(s) if s.to_lowercase() == "nan" => Ok(None),
            Value::Null => Ok(None),
            _ => {
                bail!("Metric had unexpected type {value:?}, expected numeric");
            }
        }
    }
}

impl TryFromValue for f64 {
    fn try_from_value(value: &Value) -> Result<Option<Self>> {
        Number::try_from_value(value)?
            .map(|num| {
                num.as_f64()
                    .ok_or_else(|| anyhow!("Converting {num} to float in get_metric_f64."))
            })
            .transpose()
    }
}

impl TryFromValue for usize {
    fn try_from_value(value: &Value) -> Result<Option<Self>> {
        Ok(f64::try_from_value(value)?.map(|num| num.round() as usize))
    }
}

impl TryFromValue for Percent {
    fn try_from_value(value: &Value) -> Result<Option<Self>> {
        Ok(f64::try_from_value(value)?.map(Percent::Float))
    }
}

impl TryFromValue for FloatAsInt {
    fn try_from_value(value: &Value) -> Result<Option<Self>> {
        Ok(f64::try_from_value(value)?.map(FloatAsInt))
    }
}

impl TryFromValue for String {
    fn try_from_value(value: &Value) -> Result<Option<Self>> {
        Ok(value.as_str().map(String::from))
    }
}
