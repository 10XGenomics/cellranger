//! WriteAggrAnn stage code
#![expect(missing_docs)]

use crate::parse_aggr_csv::{DONOR_HEADER, ORIGIN_HEADER};
use crate::setup_vdj_aggr::EncloneProtoMetaFormat;
use anyhow::Result;
use enclone_proto::types::Metadata;
use martian::prelude::*;
use martian_derive::{MartianStruct, make_mro};
use martian_filetypes::FileTypeRead;
use martian_filetypes::tabular_file::CsvFile;
use serde::{Deserialize, Serialize};
use std::iter::chain;

const META_HEADER_PREFIX: &str = "Meta:";

fn meta_headers(meta: &Metadata) -> impl Iterator<Item = String> {
    chain(
        [DONOR_HEADER, ORIGIN_HEADER].map(String::from),
        meta.additional_columns
            .iter()
            .map(|x| format!("{META_HEADER_PREFIX}{x}")),
    )
}

fn meta_records(meta: &Metadata, gem_well: u32) -> impl Iterator<Item = &str> {
    let info = &meta.per_gem_well_info[&gem_well];
    chain(
        [&info.donor, &info.origin],
        meta.additional_columns
            .iter()
            .map(|x| &info.additional_data[x]),
    )
    .map(String::as_str)
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

        let mut reader = csv::Reader::from_reader(args.annotation_csv.buf_reader()?);
        let headers: Vec<String> = chain(
            reader.headers()?.iter().map(String::from),
            meta_headers(&metadata),
        )
        .collect();
        let bc_col_num = headers.iter().position(|x| x == "barcode").unwrap();

        let mut writer = csv::Writer::from_writer(augmented_annotation_csv.buf_writer()?);
        writer.write_record(&headers)?;
        for record in reader.into_records() {
            let record: csv::StringRecord = record?;
            let gem_well = record[bc_col_num].rsplit_once('-').unwrap().1.parse()?;
            writer.write_record(chain(&record, meta_records(&metadata, gem_well)))?;
        }
        writer.flush()?;

        Ok(WriteAggrAnnStageOutputs {
            augmented_annotation_csv,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use itertools::assert_equal;
    use pretty_assertions::assert_eq;

    #[test]
    fn test_meta_headers() -> Result<()> {
        let meta =
            EncloneProtoMetaFormat::from("test_resources/write_aggr_ann/1d2o.em.json").read()?;
        assert_equal(meta_headers(&meta), ["donor", "origin", "Meta:AMLStatus"]);
        Ok(())
    }

    #[test]
    fn test_meta_column() -> Result<()> {
        let meta =
            EncloneProtoMetaFormat::from("test_resources/write_aggr_ann/1d2o.em.json").read()?;
        assert_equal(meta_records(&meta, 1), ["Donor 1", "PBMC", "Normal"]);
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
