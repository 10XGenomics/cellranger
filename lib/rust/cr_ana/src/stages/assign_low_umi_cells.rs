//! Martian stage ASSIGN_LOW_UMI_CELLS

use crate::io::h5;
use crate::pca::PCA_COMPONENTS;
use crate::types::{H5File, ReductionType};
use anyhow::{Result, bail};
use barcode::{Barcode, BarcodeContent};
use cr_types::reference::feature_reference::FeatureType;
use hdf5_io::matrix::read_adaptive_csr_matrix;
use itertools::Itertools;
use martian::{MartianRover, MartianStage, MartianVoid, Resource, StageDef};
use martian_derive::{MartianStruct, make_mro};
use martian_filetypes::json_file::{Json, JsonFile, JsonFormat};
use martian_filetypes::{LazyFileTypeIO, LazyWrite};
use ndarray::Axis;
use rayon::ThreadPoolBuilder;
use scan_rs::nn::knn;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::path::PathBuf;
use std::str::FromStr;

#[derive(Debug, Clone, Serialize, Deserialize)]
/// Barcode mapping.
pub struct BarcodeMapping {
    new_barcode: String,
    old_barcodes: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
/// Stage input
pub struct AssignLowUmiCellsStageInputs {
    analysis: Option<PathBuf>,
    umi_upper_threshold: Option<u32>,
    umi_lower_threshold: Option<u32>,
    aggregate_and_annotate: Option<bool>,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
/// Stage output
pub struct AssignLowUmiCellsStageOutputs {
    metacell_json: Option<JsonFile<Vec<BarcodeMapping>>>,
}

/// Martian stage ASSIGN_LOW_UMI_CELLS
pub struct AssignLowUmiCells;

const GENE_EXPRESSION_STR: &str = "Gene Expression";
const UMI_UPPER_THRESH: u32 = 500;
const UMI_LOWER_THRESH: u32 = 30;
const NUM_NEIGHBOURS: usize = 11;
const ANALYSIS_H5_DEFAULT: &str = "analysis.h5";

#[make_mro(volatile = strict)]
impl MartianStage for AssignLowUmiCells {
    type StageInputs = AssignLowUmiCellsStageInputs;
    type StageOutputs = AssignLowUmiCellsStageOutputs;
    type ChunkInputs = MartianVoid;
    type ChunkOutputs = MartianVoid;

    fn split(
        &self,
        args: Self::StageInputs,
        _rover: MartianRover,
    ) -> Result<StageDef<Self::ChunkInputs>> {
        if let Some(analysis_path) = args.analysis
            && args.aggregate_and_annotate.unwrap_or(false)
        {
            let analysis_h5: H5File = analysis_path.join(ANALYSIS_H5_DEFAULT).into();
            if !analysis_h5.exists() {
                bail!("Analysis H5 file {analysis_h5:?} does not exist");
            }

            let mem_gib = 2
                + h5::estimate_mem_gib_for_pca_outs(&analysis_h5, PCA_COMPONENTS)?.ceil() as isize
                + (1.2 * h5::estimate_mem_gib_from_nnz(&analysis_h5)?.ceil()) as isize;
            Ok(StageDef::with_join_resource(Resource::with_mem_gb(mem_gib)))
        } else {
            Ok(StageDef::with_join_resource(Resource::with_mem_gb(1)))
        }
    }

    fn main(
        &self,
        _args: Self::StageInputs,
        _chunk_args: Self::ChunkInputs,
        _rover: MartianRover,
    ) -> Result<Self::ChunkOutputs> {
        unreachable!()
    }

    fn join(
        &self,
        args: Self::StageInputs,
        _chunk_defs: Vec<Self::ChunkInputs>,
        _chunk_outs: Vec<Self::ChunkOutputs>,
        rover: MartianRover,
    ) -> Result<Self::StageOutputs> {
        ThreadPoolBuilder::new()
            .num_threads(rover.get_threads())
            .build_global()?;
        if !args.aggregate_and_annotate.unwrap_or(false) || args.analysis.is_none() {
            println!("Not running aggregation and annotation as part of low UMI cell assignment");
            Ok(AssignLowUmiCellsStageOutputs {
                metacell_json: None,
            })
        } else {
            let analysis_h5 = H5File::from(args.analysis.unwrap().join(ANALYSIS_H5_DEFAULT));
            let matrix = read_adaptive_csr_matrix(&analysis_h5, Some(GENE_EXPRESSION_STR), None)?.0;
            let mut per_bc_umis = matrix.matrix.sum_axis::<u32>(Axis(0));
            let barcodes = matrix.barcodes;
            let max_barcode_id = barcodes
                .iter()
                .map(|x| {
                    if x.starts_with("cellid")
                        && let Ok(Barcode {
                            gem_group: _,
                            valid: _,
                            content: BarcodeContent::CellName(cellid),
                        }) = Barcode::from_str(x)
                    {
                        return cellid.id as usize;
                    }
                    0
                })
                .max()
                .unwrap_or(0);

            let umi_upper_thresh = args.umi_upper_threshold.unwrap_or(UMI_UPPER_THRESH);
            let umi_lower_thresh = args.umi_lower_threshold.unwrap_or(UMI_LOWER_THRESH);

            let bcs_to_keep: Vec<_> = per_bc_umis
                .iter()
                .enumerate()
                .filter(|x| *x.1 <= umi_upper_thresh && *x.1 >= umi_lower_thresh)
                // descending order
                .sorted_by(|a, b| Ord::cmp(&b.1, &a.1))
                .map(|x| x.0)
                .collect();
            println!(
                "Reassigning {}/{} barcodes",
                bcs_to_keep.len(),
                barcodes.len()
            );

            let projection_matrix = h5::load_transformed_matrix(
                &analysis_h5,
                FeatureType::Gene,
                None,
                ReductionType::pca,
            )?;
            let (kneighbors, _ball_tree) = knn::<u32>(&projection_matrix.view(), NUM_NEIGHBOURS);
            println!("projection matrix shape {:?}", projection_matrix.dim());
            println!("knn shape {:?}", kneighbors.dim());
            assert!(projection_matrix.dim().0 == per_bc_umis.len());
            println!("per bc umis shape {:?}", per_bc_umis.dim());

            let mut old_to_new_map: HashMap<_, _> = HashMap::new();
            let mut last_meta_cell_index = max_barcode_id + 1;
            println!("Naming meta-cells with IDs starting from {last_meta_cell_index}");
            for bc in bcs_to_keep {
                let mut new_bc_umis = per_bc_umis[bc];

                // have to do this check, because the bc may initially have fewer UMIs than
                // umi_upper_thresh but got combined into a metacell in a previous
                // iteration the loop.
                if new_bc_umis > umi_upper_thresh {
                    continue;
                }
                let mut neighbour_under_considration = 0;
                while new_bc_umis <= umi_upper_thresh
                    && neighbour_under_considration < NUM_NEIGHBOURS
                {
                    // Always get into this loop if we get here as we are only iterating
                    // on barcodes with less barcodes than the umi_upper_thresh.
                    new_bc_umis +=
                        per_bc_umis[kneighbors[[bc, neighbour_under_considration]] as usize];
                    neighbour_under_considration += 1;
                    // Note that we always stop if neighbour_under_considration is a meta cell
                }

                // If we exceed the threshold; then we need to update per_bc_umis
                // and cell indices
                if new_bc_umis > umi_upper_thresh {
                    // We always stop if we hit a barcode which is a metacell because it already
                    // has more UMIs than the threshold.
                    // Hence check if the last neighbour_under_consideration was a meta cell.
                    // if so we use that meta cell ID. else generate a new meta cell ID.
                    let current_metacell_id = if let Some(cell_index) = old_to_new_map
                        .get(&(kneighbors[[bc, neighbour_under_considration - 1]] as usize))
                    {
                        *cell_index
                    } else {
                        last_meta_cell_index += 1;
                        last_meta_cell_index
                    };

                    // update metacell ID of all neighbours; and UMIs of the barcode
                    for neighbour_to_update in 0..neighbour_under_considration {
                        old_to_new_map.insert(
                            kneighbors[[bc, neighbour_to_update]] as usize,
                            current_metacell_id,
                        );
                        per_bc_umis[kneighbors[[bc, neighbour_to_update]] as usize] = new_bc_umis;
                    }
                    // The original barcode was not updated. Hence update it here.
                    per_bc_umis[bc] = new_bc_umis;
                    old_to_new_map.insert(bc, current_metacell_id);
                } else {
                    println!(
                        "Could not get {} to have more than {} UMIs despite going through {} k-NN",
                        barcodes[bc], umi_upper_thresh, NUM_NEIGHBOURS
                    );
                }
            }

            let metacell_to_cells_map = old_to_new_map
                .into_iter()
                .map(|(a, b)| (b, a))
                .into_group_map();

            println!("{} new metacells created!", metacell_to_cells_map.len());
            let lens_of_vecs = metacell_to_cells_map.iter().map(|x| x.1.len()).counts();
            println!("Statistics of new metacells created: {lens_of_vecs:#?}");

            let barcode_maps = metacell_to_cells_map
                .into_iter()
                .map(|(x, y)| BarcodeMapping {
                    new_barcode: Barcode::with_cell_id(1, x as u32).to_string(),
                    old_barcodes: y.into_iter().map(|z| barcodes[z].clone()).collect(),
                });

            let metacell_json_file: JsonFormat<Json, Vec<BarcodeMapping>> =
                rover.make_path("metacell_mappings.json");
            let mut metacell_json_writer = metacell_json_file.lazy_writer()?;
            for barcode_map in barcode_maps {
                metacell_json_writer.write_item(&barcode_map)?;
            }
            metacell_json_writer.finish()?;

            Ok(AssignLowUmiCellsStageOutputs {
                metacell_json: Some(metacell_json_file),
            })
        }
    }
}
