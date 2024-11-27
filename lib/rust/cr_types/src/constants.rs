/// The quality offset for Illumina FASTQ files.
pub const ILLUMINA_QUAL_OFFSET: u8 = 33;

/// Default minimum CRISPR UMI threshold.
pub const DEFAULT_MIN_CRISPR_UMI_THRESHOLD: usize = 3;

/// Name of the environment variable command line inputs are stashed in
pub const COMMAND_LINE_ENV_VARIABLE_NAME: &str = "CMDLINE";

/// Default value for non command line runs
pub const COMMAND_LINE_ENV_DEFAULT_VALUE: &str = "NA";
