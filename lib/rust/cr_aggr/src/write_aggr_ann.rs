//! WriteAggrAnn stage code

use crate::parse_aggr_csv::{DONOR_HEADER, ORIGIN_HEADER};
use crate::setup_vdj_aggr::EncloneProtoMetaFormat;
use anyhow::Result;
use enclone_proto::types::Metadata;
use martian::prelude::*;
use martian_derive::{make_mro, MartianStruct};
use martian_filetypes::tabular_file::CsvFile;
use martian_filetypes::FileTypeRead;
use serde::{Deserialize, Serialize};

const META_HEADER_PREFIX: &str = "Meta:";

fn meta_headers(meta: &Metadata) -> Vec<String> {
    let mut result = vec![DONOR_HEADER.to_string(), ORIGIN_HEADER.to_string()];
    for col in &meta.additional_columns {
        result.push(format!("{META_HEADER_PREFIX}{col}"));
    }
    result
}

fn meta_records(meta: &Metadata, gem_well: u32) -> Vec<String> {
    let info = &meta.per_gem_well_info[&gem_well];
    let mut result = vec![info.donor.clone(), info.origin.clone()];
    for col in &meta.additional_columns {
        result.push(info.additional_data[col].clone());
    }
    result
}

#[derive(Debug, Clone, Deserialize, MartianStruct)]
pub struct WriteAggrAnnStageInputs {
    enclone_gem_well_meta: EncloneProtoMetaFormat,
    annotation_csv: CsvFile<()>,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct WriteAggrAnnStageOutputs {
    augmented_annotation_csv: CsvFile<()>,
}

pub struct WriteAggrAnn;

#[make_mro]
impl MartianMain for WriteAggrAnn {
    type StageInputs = WriteAggrAnnStageInputs;
    type StageOutputs = WriteAggrAnnStageOutputs;
    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        let metadata = args.enclone_gem_well_meta.read()?;
        let augmented_annotation_csv: CsvFile<()> = rover.make_path("augmented_annotation");

        let mut rdr = csv::Reader::from_reader(args.annotation_csv.buf_reader()?);
        let mut headers = rdr.headers()?.clone();
        for h in meta_headers(&metadata) {
            headers.push_field(&h);
        }
        let bc_col_num = headers.iter().position(|h| h == "barcode").unwrap();

        let mut wtr = csv::Writer::from_writer(augmented_annotation_csv.buf_writer()?);
        wtr.write_record(&headers)?;

        for r in rdr.records() {
            let mut records: Vec<_> = r?.iter().map(String::from).collect();
            let gem_well: u32 = records[bc_col_num].split('-').last().unwrap().parse()?;
            records.extend(meta_records(&metadata, gem_well));
            wtr.write_record(&records)?;
        }

        Ok(WriteAggrAnnStageOutputs {
            augmented_annotation_csv,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use pretty_assertions::assert_eq;

    #[test]
    fn test_meta_headers() -> Result<()> {
        let meta =
            EncloneProtoMetaFormat::from("test_resources/write_aggr_ann/1d2o.em.json").read()?;
        assert_eq!(meta_headers(&meta), ["donor", "origin", "Meta:AMLStatus"]);
        Ok(())
    }

    #[test]
    fn test_meta_column() -> Result<()> {
        let meta =
            EncloneProtoMetaFormat::from("test_resources/write_aggr_ann/1d2o.em.json").read()?;
        assert_eq!(meta_records(&meta, 1), ["Donor 1", "PBMC", "Normal"]);
        Ok(())
    }

    #[test]
    fn test_ann_csv() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let outs = WriteAggrAnn.test_run(
            &dir,
            WriteAggrAnnStageInputs {
                enclone_gem_well_meta: "test_resources/write_aggr_ann/1d2o.em.json".into(),
                annotation_csv: "test_resources/write_aggr_ann/filtered_ann.csv".into(),
            },
        )?;
        assert_eq!(
            std::fs::read_to_string(outs.augmented_annotation_csv)?,
            std::fs::read_to_string("test_resources/write_aggr_ann/expected_ann.csv")?
        );
        Ok(())
    }
}
