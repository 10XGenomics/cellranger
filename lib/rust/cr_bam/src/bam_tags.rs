use bitflags::bitflags;

pub const READ_GROUP_TAG: &[u8] = b"RG";
pub const PROC_BC_SEQ_TAG: &[u8] = b"CB";
pub const PROC_UMI_SEQ_TAG: &[u8] = b"UB";
pub const FEATURE_RAW_TAG: &[u8] = b"fr";
pub const FEATURE_QUAL_TAG: &[u8] = b"fq";
pub const FEATURE_SEQ_TAG: &[u8] = b"fb";
pub const FEATURE_IDS_TAG: &[u8] = b"fx";
pub const PROBE_TAG: &[u8] = b"pr";
pub const GAP_COORDINATES_TAG: &[u8] = b"gc";
pub const EXTRA_FLAGS_TAG: &[u8] = b"xf";

pub const RAW_UMI_SEQ_TAG: &[u8] = b"UR";
pub const RAW_UMI_QUAL_TAG: &[u8] = b"UY";

pub const RAW_BARCODE_SEQ_TAG: &[u8] = b"CR";
pub const RAW_BARCODE_QUAL_TAG: &[u8] = b"CY";

pub const RAW_GEL_BEAD_BARCODE_SEQ_TAG: &[u8] = b"GR";
pub const RAW_GEL_BEAD_BARCODE_QUAL_TAG: &[u8] = b"GY";

pub const REST_R1_SEQ_TAG: &[u8] = b"1R";
pub const REST_R1_QUAL_TAG: &[u8] = b"1Y";

pub const REST_R2_SEQ_TAG: &[u8] = b"2R";
pub const REST_R2_QUAL_TAG: &[u8] = b"2Y";

// These only exist in BAM land; just use bytes.
pub const TRANSCRIPT_TAG: &[u8] = b"TX";
pub const GENE_ID_TAG: &[u8] = b"GX";
pub const GENE_NAME_TAG: &[u8] = b"GN";
pub const REGION_TAG: &[u8] = b"RE";
pub const MULTIMAPPER_TAG: &[u8] = b"mm";
pub const ANTISENSE_TAG: &[u8] = b"AN"; // equivalent to TX, but for antisense alignments

// These are set to the original single-read annotations
// if a read-pair had gene disagreement.
pub const UNPAIRED_GENE_ID_TAG: &[u8] = b"gX";
pub const UNPAIRED_GENE_NAME_TAG: &[u8] = b"gN";

bitflags! {
    #[derive(Default)]
    /// Extra alignment bit flags, stored in the BAM tag xf:i.
    pub struct ExtraFlags: u32 {
        /// Confidently mapped to transcriptome
        const CONF_MAPPED = 1u32;
        /// This read's (BC, UMI, feature) combination was discarded in favor of a different feature
        /// with higher read support.
        const LOW_SUPPORT_UMI = 2u32;
        /// Mates mapped to incompatible sets of genes
        const GENE_DISCORDANT = 4u32;
        /// This read is representative for a molecule and can be treated as a UMI count.
        const UMI_COUNT = 8u32;
        /// Confidently assigned feature barcode
        const CONF_FEATURE = 16u32;
        /// This read's (BC, UMI, feature) combination was not counted only because the read
        /// count was too low. Implies not a LOW_SUPPORT_UMI
        const FILTERED_TARGET_UMI = 32u32;
    }
}
