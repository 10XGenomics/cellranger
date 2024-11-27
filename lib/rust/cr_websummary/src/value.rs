use crate::{CountAndPercent, FloatAsInt, Percent};
use anyhow::{anyhow, bail, Result};
use log::warn;
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
                bail!("Metric had unexpected type {:?}, expected numeric", value);
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

impl TryFromValue for CountAndPercent {
    fn try_from_value(_: &Value) -> Result<Option<Self>> {
        // TEMP: log a warning and return None
        // FIXME CELLRANGER-8444
        // We should never end up using the result of this call in output,
        // since all CountAndPercent fields need to be manually constructed.
        // However, we cannot use struct update syntax for populating metrics
        // tables for any tables that contain CountAndPercent if we panic here.
        // unimplemented!("TryFromValue is not implemented for CountAndPercent")
        warn!("Attempted to call try_from_value for CountAndPercent.");
        Ok(None)
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
