//!
//! Alarms in the websummary (aka alerts)
//!

use serde::Serialize;

#[derive(Debug, Clone, Serialize)]
#[serde(rename_all = "UPPERCASE")]
pub enum AlertLevel {
    Error,
    Warn,
    Info,
}

/// Specification of Alarm/Alert that appear in the web summary
///
///
/// NOTE: Currently in cellranger, an alert specification looks like this:
/// ```text
/// {
///     "raw_value": 0.15510204081632653,
///     "raised": true,
///     "parent": "multi_vdj_assembly_contig_pair_productive_full_len_bc_frac",
///     "level": "ERROR",
///     "title": "Low Cells With Productive V-J Spanning Pair",
///     "formatted_value": "15.5%",
///     "test": "",
///     "message": "Ideal > 30%. This can indicate poor cell quality, low yield from the RT reaction, poor specificity of the V(D)J enrichment, poor sequencing quality, or the use of an unsupported chemistry type (e.g., using Single Cell 3' for V(D)J assembly). Application performance is likely to be affected",
///     "id": "multi_vdj_assembly_contig_pair_productive_full_len_bc_frac"
/// }
/// ```
/// Only the fields `level`, `title`, `formatted_value` and `message` are used in the web summary
/// react code.
#[derive(Debug, Clone, Serialize)]
pub struct AlertSpec {
    pub level: AlertLevel,
    pub title: String,
    pub formatted_value: String,
    pub message: String,
}

/// Common alert configuration data that may be relevant in many contexts.
#[derive(Debug, Clone, Default)]
pub struct AlertContext {
    pub is_hybrid_capture: bool,
    pub is_rtl: bool,
    pub is_lt_chemistry: bool,
    pub is_arc_chemistry: bool,
    pub is_fiveprime: bool,
    pub is_multiplexing: bool,
    pub is_antigen: bool,
    pub include_introns: bool,
    pub no_preflight: bool,
}

pub trait Alert {
    fn alerts(&self, _ctx: &AlertContext) -> Vec<AlertSpec> {
        // Default impl has no alerts
        vec![]
    }
}
