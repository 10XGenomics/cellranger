//! WriteContigProto stage code

use anyhow::Result;
use cr_types::MetricsFile;
use martian::prelude::*;
use martian_derive::{make_mro, martian_filetype, MartianStruct};
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::{FileTypeRead, LazyFileTypeIO};
use serde::{Deserialize, Serialize};
use std::path::PathBuf;
use vdj_ann::annotate::ContigAnnotation;
use vdj_asm_asm::BarcodeDataBriefFile;
use vdj_proto::io::{VdjProtoWriter, PROTOBUF_VERSION};
use vdj_proto::types::{MetricsSummary, VdjMetadata, VdjReferenceRaw};
use vdj_reference::VdjReceptor;

martian_filetype! {ProtoFile, "pb"}

#[derive(Debug, Clone, Deserialize, MartianStruct)]
pub struct WriteContigProtoStageInputs {
    vdj_reference_path: PathBuf,
    contig_annotations_json: JsonFile<Vec<ContigAnnotation>>,
    metrics_summary_json: MetricsFile,
    receptor: VdjReceptor,
    gem_wells: Vec<u32>,
    cell_barcodes: JsonFile<Vec<String>>,
    sample_id: String,
    sample_desc: String,
    multi_config_sha: Option<String>,
    barcode_brief: BarcodeDataBriefFile,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct WriteContigProtoStageOutputs {
    vdj_contig_info: ProtoFile,
}

// This is our stage struct
pub struct WriteContigProto;

#[make_mro]
impl MartianMain for WriteContigProto {
    type StageInputs = WriteContigProtoStageInputs;
    type StageOutputs = WriteContigProtoStageOutputs;
    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        let proto_file: ProtoFile = rover.make_path("vdj_contig_info");

        let reference = VdjReferenceRaw::new(&args.vdj_reference_path)?;

        let cell_barcodes = args.cell_barcodes.read()?;
        let metadata = VdjMetadata {
            reference_fasta_hash: reference.fasta_hash(),
            pipeline_version: rover.pipelines_version(),
            receptor: vdj_proto::types::Receptor::from(args.receptor).into(),
            gem_wells: args.gem_wells,
            number_of_cells: cell_barcodes.len() as u32,
            sample_id: args.sample_id,
            sample_desc: args.sample_desc,
            multi_config_sha: args.multi_config_sha.unwrap_or_default(),
            protobuf_version: String::from(PROTOBUF_VERSION),
        };

        let metrics = MetricsSummary::from_metrics_json(&args.metrics_summary_json)?;

        let mut writer = VdjProtoWriter::new(&proto_file, metadata, reference, metrics)?;
        for ann in args.contig_annotations_json.lazy_reader()? {
            writer.write_annotation(ann?)?;
        }
        for brief in args.barcode_brief.lazy_reader()? {
            writer.write_barcode_data(brief?)?;
        }
        writer.finish()?;
        Ok(WriteContigProtoStageOutputs {
            vdj_contig_info: proto_file,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use itertools::Itertools;
    use martian_filetypes::FileTypeWrite;
    use metric::JsonReporter;
    use vdj_ann::annotate::ContigAnnotation;
    use vdj_asm_utils::barcode_data::BarcodeDataBrief;
    use vdj_proto::io::VdjProtoReader;

    #[test]
    fn test_write_proto() -> Result<()> {
        let contig_annotations_json = JsonFile::from("../vdj_asm_asm/test.json");
        let anns: Vec<ContigAnnotation> = contig_annotations_json.read()?;
        let cells: Vec<_> = anns
            .iter()
            .filter_map(|ann| {
                if ann.is_cell {
                    Some(ann.barcode.clone())
                } else {
                    None
                }
            })
            .collect();
        let dir = tempfile::tempdir()?;
        let cell_barcodes_json = JsonFile::new(&dir, "cell_barcodes");
        cell_barcodes_json.write(&cells)?;
        let metrics_summary_json = MetricsFile::new(&dir, "metrics");
        metrics_summary_json.write(&JsonReporter::default())?;
        let barcode_data_brief_file = BarcodeDataBriefFile::new(&dir, "barcode_brief");
        let barcode_brief = vec![BarcodeDataBrief::new()];
        barcode_data_brief_file.write(&barcode_brief)?;

        let args = WriteContigProtoStageInputs {
            vdj_reference_path: PathBuf::from(
                "../dui_tests/test_resources/reference/vdj_GRCh38_alts_ensembl-4.0.0",
            ),
            contig_annotations_json,
            receptor: VdjReceptor::TR,
            gem_wells: vec![1],
            cell_barcodes: cell_barcodes_json,
            sample_id: "sample_id".into(),
            sample_desc: "sample_desc".into(),
            multi_config_sha: None,
            metrics_summary_json,
            barcode_brief: barcode_data_brief_file,
        };

        let outs = WriteContigProto.test_run(&dir, args)?;
        let anns_proto: Vec<_> =
            VdjProtoReader::read_annotations(&outs.vdj_contig_info)?.try_collect()?;
        assert_eq!(anns, anns_proto);
        Ok(())
    }
}
