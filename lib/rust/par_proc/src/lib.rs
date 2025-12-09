//! Parallel processing
#![deny(missing_docs)]

mod par_proc;

pub use par_proc::{MAX_ITEMS_IN_MEM, Proc, group_by_processor, group_by_processor_in_memory};
