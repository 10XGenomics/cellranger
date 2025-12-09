//! multiconfig_converter
#![deny(missing_docs)]
use anyhow::Result;
use multi::config::MultiConfigCsv;
use serde_json::to_vec;
use std::slice::{from_raw_parts, from_raw_parts_mut};

fn safe_convert_multi_csv(input: &[u8]) -> Result<Vec<u8>> {
    let config = MultiConfigCsv::from_reader(input, "api")?;
    Ok(to_vec(&config)?)
}

/// convert_multi_config is a C-ABI function for converting a multi CSV to json.
///
/// # Arguments
///
/// * `buffer`: a destination buffer to write the resulting content to.
/// * `buf_len`: the size of the buffer.  Content will only be written to the
///   buffer if this is greater than or equal to the return value.
/// * `input`: The input data to process.
/// * `input_len`: The number of bytes in `input`.
///
/// # Returns
///
/// The return value is the number of bytes in the result json.
/// If (and only if) this value is less than or equal to the `buf_len` argument,
/// the content will have been written back to `buffer`.
/// Otherwise, `buffer` will not be touched, so it's ok to call this function
/// once with NULL buffer to get the appropriate length and then make a second
/// call with the appropriately-allocated buffer.
///
/// A negative return value indicates an error.  In that case, the buffer is
/// populated with the (possibly truncated) error message.
///
/// # Safety
///
/// `buffer` must point to an allocation with space for at least `buf_len`
/// bytes, unless `buf_len` is zero.
///
/// The `input` array must contain at least `input_len` bytes.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn convert_multi_config(
    buffer: *mut u8,
    buf_len: usize,
    input: *const u8,
    input_len: usize,
) -> i64 {
    match safe_convert_multi_csv(unsafe { from_raw_parts(input, input_len) }) {
        Ok(v) => {
            if v.len() <= buf_len && buf_len > 0 {
                unsafe { from_raw_parts_mut(buffer, v.len()) }.copy_from_slice(&v);
            }
            v.len().try_into().unwrap_or(i64::MAX)
        }
        Err(e) => {
            let e_str = format!("{e}");
            let e_bytes = e_str.as_bytes();
            if buf_len > 0 {
                if e_bytes.len() <= buf_len {
                    unsafe { from_raw_parts_mut(buffer, e_bytes.len()) }.copy_from_slice(e_bytes);
                } else {
                    unsafe { from_raw_parts_mut(buffer, buf_len) }
                        .copy_from_slice(&e_bytes[..buf_len]);
                }
            }
            -i64::try_from(e_bytes.len()).unwrap_or(i64::MAX)
        }
    }
}
