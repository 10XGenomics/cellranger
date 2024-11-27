use crate::reference::probe_set_reference::TargetSetFile;
use crate::types::TargetingMethod;
use martian_derive::martian_filetype;
use martian_filetypes::json_file::JsonFormat;
use serde::{Deserialize, Serialize};

martian_filetype!(_TargetPanelSummaryFile, "tps");
pub type TargetPanelSummaryFormat = JsonFormat<_TargetPanelSummaryFile, TargetPanelSummary>;

#[derive(Serialize, Deserialize)]
pub struct TargetPanelSummary {
    pub target_panel_hash: String,
    pub target_panel_name: String,
    pub target_panel_path: TargetSetFile,
    pub target_panel_gene_count: i32,
    pub target_panel_type: String,
    pub targeting_method: TargetingMethod,
}

#[cfg(test)]
mod tests {
    use super::*;
    use martian_filetypes::FileTypeRead;
    use std::io::Write;

    #[test]
    fn test_deserialize_python_target_panel_summary() {
        let data = r#"{
    "target_panel_hash": "6dcce063f83fdd9e3ebf352a3e45656be9b6cf79",
    "target_panel_name": "Spatial Neuroscience",
    "target_panel_path": "Spatial_neuroscience_genes.csv",
    "target_panel_gene_count": 1186,
    "target_panel_type": "custom",
    "targeting_method": "hybrid_capture"
    }"#;
        let mut file = tempfile::Builder::new()
            .suffix(".tps.json")
            .tempfile()
            .unwrap();
        write!(file.as_file_mut(), "{data}").expect("Could not write data");
        file.as_file().flush().expect("Could not flush file");
        let tps = TargetPanelSummaryFormat::from(file.as_ref())
            .read()
            .unwrap();
        assert_eq!(tps.targeting_method, TargetingMethod::HybridCapture);
    }
}
