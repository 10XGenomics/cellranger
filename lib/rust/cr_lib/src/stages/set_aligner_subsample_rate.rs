//! Martian stage SET_ALIGNER_SUBSAMPLE_RATE

use anyhow::Result;
use barcode::Barcode;
use cr_types::types::LibraryFeatures;
use cr_types::BcCountFormat;
use itertools::Itertools;
use martian::{MartianMain, MartianRover};
use martian_derive::{make_mro, MartianStruct};
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::FileTypeRead;
use serde::{Deserialize, Serialize};

/// Compute susbsample rate for gene expression data.
pub struct SetAlignerSubsampleRate;

#[derive(Deserialize, Clone, MartianStruct)]
pub struct StageInputs {
    pub barcodes_under_tissue: Option<JsonFile<Vec<String>>>,
    pub corrected_barcode_counts: BcCountFormat,
    pub rps_limit: Option<u64>,
}

#[derive(Serialize, Deserialize, Clone, MartianStruct)]
pub struct StageOutputs {
    pub aligner_subsample_rate: Option<f64>,
}

#[make_mro(mem_gb = 8, volatile = strict)]
impl MartianMain for SetAlignerSubsampleRate {
    type StageInputs = StageInputs;
    type StageOutputs = StageOutputs;

    fn main(&self, args: Self::StageInputs, _rover: MartianRover) -> Result<Self::StageOutputs> {
        if args.rps_limit.is_none() {
            return Ok(StageOutputs {
                aligner_subsample_rate: None,
            });
        }
        let rps_limit = args.rps_limit.unwrap();
        assert!(
            rps_limit > 0,
            "Error: rps_limit must be non-zero. Use rps_limit=null to disable subsampling."
        );

        let barcodes_under_tissue_json = args
            .barcodes_under_tissue
            .expect("Error: barcodes_under_tissue is required for subsampling.");

        let barcode_counts = args
            .corrected_barcode_counts
            .read()?
            .into_iter()
            .filter_map(|(k, v)| {
                if matches!(k, LibraryFeatures::GeneExpression(_)) {
                    Some(v)
                } else {
                    None
                }
            })
            .exactly_one()
            .unwrap();

        let barcodes_under_tissue = barcodes_under_tissue_json.read()?;
        let barcodes_under_tissue: Vec<Barcode> = barcodes_under_tissue
            .iter()
            .map(|x| x.parse().unwrap())
            .collect();
        let reads_under_tissue = barcodes_under_tissue
            .iter()
            .map(|x| barcode_counts.get(x))
            .sum::<i64>() as u64;

        let desired_reads_under_tissue = rps_limit * barcodes_under_tissue.len() as u64;
        let aligner_subsample_rate = Some(if desired_reads_under_tissue < reads_under_tissue {
            desired_reads_under_tissue as f64 / reads_under_tissue as f64
        } else {
            1.0
        });

        Ok(StageOutputs {
            aligner_subsample_rate,
        })
    }
}
