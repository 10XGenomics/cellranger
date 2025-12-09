//! MergePerSampleAnnotations stage code
#![expect(missing_docs)]
use anyhow::Result;
use cr_h5::count_matrix::BarcodeWithGemGroup;
use cr_lib::map_multiplexing_seq_to_id;
use cr_types::chemistry::ChemistryDef;
use cr_types::{BarcodeMultiplexingType, ReadLevel};
use hdf5::types::FixedAscii;
use itertools::Itertools;
use martian::prelude::*;
use martian_derive::{MartianStruct, make_mro};
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::{FileTypeRead, FileTypeWrite, LazyFileTypeIO, LazyWrite};
use metric::TxHashMap;
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use vdj_ann::annotate::ContigAnnotation;
#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct MergePerSampleAnnotationsStageInputs {
    pub per_sample_annotations: HashMap<String, JsonFile<Vec<ContigAnnotation>>>,
    pub asm_contig_annotations: JsonFile<Vec<ContigAnnotation>>,
    pub multiplexing_method: Option<BarcodeMultiplexingType>,
    pub gex_cells_per_tag: Option<JsonFile<TxHashMap<String, HashSet<String>>>>,
    pub vdj_chemistry_def: ChemistryDef,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct MergePerSampleAnnotationsStageOutputs {
    pub contig_annotations: JsonFile<Vec<ContigAnnotation>>,
    pub vdj_cells_per_tag_json: Option<JsonFile<HashMap<String, usize>>>,
}

pub struct MergePerSampleAnnotations;

fn get_vdj_cells_per_tag_cell_level_multiplexing(
    gex_cells_per_tag: Option<JsonFile<TxHashMap<String, HashSet<String>>>>,
    cell_bcs: &HashSet<String>,
) -> HashMap<String, usize> {
    gex_cells_per_tag
        .unwrap()
        .read()
        .unwrap()
        .into_iter()
        .map(|(tag_id, gex_cells_set)| {
            (
                tag_id.to_string(),
                gex_cells_set.intersection(cell_bcs).count(),
            )
        })
        .collect()
}

fn get_vdj_cells_per_tag_overhang_multiplexing(
    vdj_chemistry_def: ChemistryDef,
    cell_bcs: &HashSet<String>,
) -> HashMap<String, usize> {
    let overhang_barcodes = vdj_chemistry_def.overhang_read_barcode().unwrap();

    let overhang_offset = overhang_barcodes.offset();
    let overhang_length = overhang_barcodes.length();
    let overhang_range = overhang_offset..(overhang_offset + overhang_length);

    // this contains a path to a text file with information about the overhang sequence to overhang ID map.
    let whitelist_sources = overhang_barcodes.whitelist().as_source().unwrap();

    let overhang_seq_to_id: TxHashMap<_, _> = whitelist_sources.as_translation_seq_to_id().unwrap();

    // Do this initialization. that way even if no barcodes are assigned to a particular overhang, we will still populate it in the table.
    let vdj_cells_per_tag_initial: HashMap<String, usize> = overhang_seq_to_id
        .values()
        .map(|v| (v.to_string().clone(), 0))
        .collect();

    cell_bcs
        .iter()
        .map(|barcode| {
            map_multiplexing_seq_to_id(
                &BarcodeWithGemGroup::from(FixedAscii::from_ascii(barcode.as_bytes()).unwrap()),
                &overhang_seq_to_id,
                &overhang_range,
            )
        })
        .fold(vdj_cells_per_tag_initial, |mut acc, overhang_id| {
            *acc.entry(overhang_id.to_string()).or_insert(0) += 1;
            acc
        })
}

#[make_mro]
impl MartianMain for MergePerSampleAnnotations {
    type StageInputs = MergePerSampleAnnotationsStageInputs;
    type StageOutputs = MergePerSampleAnnotationsStageOutputs;
    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        // In case of CMO & HASHTAG multiplexing sample-level assignments can only be made for cell bcs
        let only_cell_bcs_in_samples = matches!(
            args.multiplexing_method,
            Some(BarcodeMultiplexingType::CellLevel(_))
        );

        let contig_annotations: JsonFile<Vec<ContigAnnotation>> =
            rover.make_path("contig_annotations");
        let mut ann_writer = contig_annotations.lazy_writer()?;
        let mut cell_bcs: HashSet<String> = HashSet::new();
        for (sample, annot) in args
            .per_sample_annotations
            .into_iter()
            .sorted_by_key(|(sample_id, _)| sample_id.clone())
        {
            let contig_reader = annot.lazy_reader()?;
            for ann in contig_reader {
                let mut ann: ContigAnnotation = ann?;
                // ann.clonotype = None;
                // ann.info = Default::default();
                ann.sample = Some(sample.clone());
                ann_writer.write_item(&ann)?;
                if ann.is_cell {
                    cell_bcs.insert(ann.barcode);
                }
            }
        }

        let vdj_cells_per_tag = if let Some(multiplexing_method) = args.multiplexing_method {
            match multiplexing_method {
                BarcodeMultiplexingType::CellLevel(_) => {
                    Some(get_vdj_cells_per_tag_cell_level_multiplexing(
                        args.gex_cells_per_tag,
                        &cell_bcs,
                    ))
                }
                BarcodeMultiplexingType::ReadLevel(ReadLevel::OH) => Some(
                    get_vdj_cells_per_tag_overhang_multiplexing(args.vdj_chemistry_def, &cell_bcs),
                ),
                BarcodeMultiplexingType::ReadLevel(ReadLevel::RTL) => {
                    panic!("Unsupported vdj multiplexing method!")
                }
            }
        } else {
            None
        };

        let vdj_cells_per_tag_json: Option<JsonFile<HashMap<String, usize>>> =
            if let Some(vdj_cells_per_tag) = vdj_cells_per_tag {
                Some(
                    rover
                        .make_path::<JsonFile<HashMap<String, usize>>>("vdj_cells_per_tag")
                        .with_content(&vdj_cells_per_tag)?,
                )
            } else {
                None
            };

        // Add non-cell bcs
        if only_cell_bcs_in_samples {
            for ann in args.asm_contig_annotations.lazy_reader()? {
                let mut ann: ContigAnnotation = ann?;
                if !cell_bcs.contains(&ann.barcode) {
                    ann.sample = None;
                    ann_writer.write_item(&ann)?;
                }
            }
        }
        ann_writer.finish()?;

        Ok(MergePerSampleAnnotationsStageOutputs {
            contig_annotations,
            vdj_cells_per_tag_json,
        })
    }
}
