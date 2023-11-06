// Bug: we report VDJ_total_reads = VDJ_total_reads_pairs.

use serde::{Deserialize, Serialize};
use std::collections::HashMap;

#[derive(Serialize, Deserialize, Default)]
pub struct AssemblyMetrics {
    pub assemblable_read_pairs_by_bc: HashMap<String, u64>,
}
