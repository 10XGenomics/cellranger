// Constants used for ShardIO operations
// in align_and_count and barcode_sort
// we only need these here because we need DISK_CHUNK_SZ to estimate memory usage. Once ShardReader
// has API to provide that, these can return to the stages where they are used
// TODO: Set the parameters optimally?
pub const ALN_BC_SEND_BUFFER_SZ: usize = 256;
pub const ALN_BC_ITEM_BUFFER_SZ: usize = 1_048_576;
pub const ALN_BC_DISK_CHUNK_SZ: usize = 8_192;
// 8MiB if ~1kB per record
pub const ALN_BC_GIB: f64 = 1_073_741_824.0;
