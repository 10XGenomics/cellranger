//!
//! Compute per_barcode.csv for BEAM analysis. Note that this is different from the
//! per_barcode_metrics.csv computed in BASIC_SC_RNA_COUNTER.
//!

use crate::matrix::{load_barcodes_from_matrix, H5File};
use anyhow::{ensure, Result};
use martian::prelude::*;
use martian_derive::{make_mro, MartianStruct, MartianType};
use martian_filetypes::tabular_file::CsvFile;
use martian_filetypes::{LazyFileTypeIO, LazyWrite};
use metric::{TxHashMap, TxHashSet};
use serde::{Deserialize, Serialize};
use std::collections::hash_map::Entry;
use std::collections::{BTreeSet, HashMap};
use vdj_asm_asm::write_ann_csv::ContigAnnotationCsvRow;

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, MartianType)]
#[serde(try_from = "GemInfoHelper", into = "GemInfoHelper")]
pub struct GemInfo {
    sample_id: String,
    index: i64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct GemInfoHelper([serde_json::Value; 2]);

impl From<GemInfo> for GemInfoHelper {
    fn from(value: GemInfo) -> Self {
        GemInfoHelper([value.sample_id.into(), value.index.into()])
    }
}

impl TryFrom<GemInfoHelper> for GemInfo {
    type Error = anyhow::Error;
    fn try_from(value: GemInfoHelper) -> Result<Self, Error> {
        ensure!(value.0.len() == 2);
        Ok(GemInfo {
            sample_id: value.0[0].as_str().unwrap().to_string(),
            index: value.0[1].to_string().parse::<i64>()?,
        })
    }
}

// We require both GEX and VDJ libraries to be present in a beam analysis, so
// all the inputs are non-optional
#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct StageInputs {
    gex_filtered_matrix: H5File,
    vdj_filtered_annotations: CsvFile<ContigAnnotationCsvRow>,
    count_gem_well_map: Option<HashMap<String, GemInfo>>,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct StageOutputs {
    per_barcode_csv: CsvFile<()>,
}

// BarcodeCsvRowGeneral includes optional sample_id, only used in aggr
#[derive(Serialize, Deserialize)]
struct BarcodeCsvRowGeneral {
    barcode: String,
    sample_id: Option<String>,
    is_gex_cell: bool,
    is_vdj_cell: bool,
    raw_clonotype_id: Option<String>,
    exact_subclonotype_id: Option<String>,
}

impl From<BarcodeCsvRowGeneral> for BarcodeCsvRow {
    fn from(value: BarcodeCsvRowGeneral) -> Self {
        BarcodeCsvRow {
            barcode: value.barcode,
            is_gex_cell: value.is_gex_cell,
            is_vdj_cell: value.is_vdj_cell,
            raw_clonotype_id: value.raw_clonotype_id,
            exact_subclonotype_id: value.exact_subclonotype_id,
        }
    }
}

#[derive(Serialize, Deserialize)]
struct BarcodeCsvRow {
    barcode: String,
    is_gex_cell: bool,
    is_vdj_cell: bool,
    raw_clonotype_id: Option<String>,
    exact_subclonotype_id: Option<String>,
}

#[derive(PartialEq, Eq, Debug)]
struct VdjCellInfo {
    raw_clonotype_id: String,
    exact_subclonotype_id: String,
}

pub struct CreateBarcodeCsv;

#[make_mro]
impl MartianMain for CreateBarcodeCsv {
    type StageInputs = StageInputs;
    type StageOutputs = StageOutputs;
    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        let per_barcode_csv: CsvFile<()> = rover.make_path("per_barcode");

        let gex_cells: TxHashSet<_> = load_barcodes_from_matrix(&args.gex_filtered_matrix)?
            .into_iter()
            .collect();
        let mut vdj_cell_info = TxHashMap::default();

        let ann_reader = args.vdj_filtered_annotations.lazy_reader()?;

        for ann in ann_reader {
            let ContigAnnotationCsvRow {
                barcode,
                is_cell,
                raw_clonotype_id,
                exact_subclonotype_id,
                ..
            } = ann?;
            assert!(is_cell); // Sanity check
            if raw_clonotype_id.is_none() {
                continue; // CELLRANGER-6602
            }
            let info = VdjCellInfo {
                raw_clonotype_id: raw_clonotype_id.unwrap(),
                exact_subclonotype_id: exact_subclonotype_id.unwrap(),
            };
            match vdj_cell_info.entry(barcode) {
                Entry::Occupied(val) => assert_eq!(val.get(), &info), // Sanity check
                Entry::Vacant(v) => {
                    v.insert(info);
                }
            }
        }

        let all_barcodes: BTreeSet<_> = gex_cells
            .iter()
            .chain(vdj_cell_info.keys())
            .cloned()
            .collect();

        let mut writer = None;
        let mut writer_aggr = None;
        if args.count_gem_well_map.is_none() {
            writer = Some(CsvFile::from_path(&per_barcode_csv).lazy_writer()?);
        } else {
            writer_aggr = Some(CsvFile::from_path(&per_barcode_csv).lazy_writer()?);
        }

        for barcode in all_barcodes {
            let is_gex_cell = gex_cells.contains(&barcode);

            let sample_id = args
                .count_gem_well_map
                .as_ref()
                .and_then(|gem_well_map| gem_well_map.get(barcode.split('-').last().unwrap()))
                .map(|gw_info| gw_info.sample_id.clone());

            let row = match vdj_cell_info.remove(&barcode) {
                Some(VdjCellInfo {
                    raw_clonotype_id,
                    exact_subclonotype_id,
                }) => BarcodeCsvRowGeneral {
                    barcode,
                    sample_id,
                    is_gex_cell,
                    is_vdj_cell: true,
                    raw_clonotype_id: Some(raw_clonotype_id),
                    exact_subclonotype_id: Some(exact_subclonotype_id),
                },
                None => BarcodeCsvRowGeneral {
                    barcode,
                    sample_id,
                    is_gex_cell,
                    is_vdj_cell: false,
                    raw_clonotype_id: None,
                    exact_subclonotype_id: None,
                },
            };
            if let Some(wr_aggr) = &mut writer_aggr {
                wr_aggr.write_item(&row)?;
            }
            if let Some(wr) = &mut writer {
                wr.write_item(&BarcodeCsvRow::from(row))?;
            }
        }
        if let Some(wr) = writer {
            wr.finish()?;
        }
        if let Some(wr_aggr) = writer_aggr {
            wr_aggr.finish()?;
        }

        Ok(StageOutputs { per_barcode_csv })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use pretty_assertions::assert_eq;

    #[test]
    fn test_gem_well_info() -> Result<()> {
        let gw_info = GemInfo {
            sample_id: "s1".to_string(),
            index: 1,
        };
        let serialized = serde_json::to_string(&gw_info)?;
        assert_eq!(serialized, "[\"s1\",1]");
        let deserialized: GemInfo = serde_json::from_str(&serialized)?;
        assert_eq!(deserialized, gw_info);

        Ok(())
    }
}
