//! slide_design
#![deny(missing_docs)]
#[cfg(feature = "tenx_source_available")]
mod stubs;
#[cfg(feature = "tenx_source_available")]
pub use stubs::*;
#[cfg(feature = "tenx_internal")]
mod internal;
#[cfg(feature = "tenx_internal")]
pub use internal::*;
