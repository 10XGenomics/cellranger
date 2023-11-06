//! Martian stage CREATE_MULTI_GRAPH
//! Create the multiplexing sample graph.
//! This stage requires information from both the multi config CSV and chemistry detection.
use super::detect_chemistry::DetectedProbeBarcodePairingFile;
use anyhow::Result;
use cr_types::CrMultiGraph;
use martian::prelude::*;
use martian_derive::{make_mro, MartianStruct};
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::{FileTypeRead, FileTypeWrite};
use multi::config::MultiConfigCsvFile;
use serde::{Deserialize, Serialize};

pub type CrMultiGraphFile = JsonFile<CrMultiGraph>;

#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct CreateMultiGraphStageInputs {
    pub sample_id: String,
    pub sample_desc: String,
    pub multi_config: MultiConfigCsvFile,
    pub detected_probe_barcode_pairing: Option<DetectedProbeBarcodePairingFile>,
}

#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct CreateMultiGraphStageOutputs {
    pub multi_graph: CrMultiGraphFile,
}

pub struct CreateMultiGraph;

#[make_mro(volatile = strict)]
impl MartianMain for CreateMultiGraph {
    type StageInputs = CreateMultiGraphStageInputs;
    type StageOutputs = CreateMultiGraphStageOutputs;

    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        let detected_probe_barcode_pairing = args
            .detected_probe_barcode_pairing
            .map(|f| f.read())
            .transpose()?;
        let multi_graph: JsonFile<_> = rover.make_path("multi_graph");
        multi_graph.write(&args.multi_config.read()?.to_multi_graph(
            &args.sample_id,
            &args.sample_desc,
            detected_probe_barcode_pairing.as_ref(),
        )?)?;
        Ok(Self::StageOutputs { multi_graph })
    }
}
