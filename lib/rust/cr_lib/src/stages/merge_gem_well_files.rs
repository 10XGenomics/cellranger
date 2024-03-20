//! Martian stage MERGE_GEM_WELL_FILES
//! Take in lists of matrix computer intermediary files for each GEM well
//! Merge each into a single list for multi-GEM well

use crate::types::{BarcodeMetricsShardFile, FeatureReferenceFormat};
use crate::utils::hard_link_martianfile;
use crate::{AlignShardFile, BcUmiInfoShardFile};
use anyhow::{bail, Result};
use cr_types::rna_read::RnaChunk;
use martian::prelude::*;
use martian_derive::{make_mro, MartianStruct};
use serde::{Deserialize, Serialize};
use std::path::{Path, PathBuf};
use tx_annotation::read::AnnotationFiles;

pub struct MergeGemWellFiles;

// intermediary matrix_computer files for a single gem well
#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct GemWellFiles {
    // things that we need to combine
    pub gem_groups: Vec<usize>,
    pub alignments: Vec<AlignShardFile>,
    pub read_chunks: Vec<RnaChunk>,
    pub bc_umi_info: Vec<BcUmiInfoShardFile>,
    pub per_barcode_metrics_shard: Vec<BarcodeMetricsShardFile>,
    pub annotation_files: Option<AnnotationFiles>,
    // things that we only need one of
    pub target_set_name: Option<String>,
    pub bam_header: PathBuf,
    pub slfe_feature_reference: FeatureReferenceFormat,
}

#[derive(Clone, Deserialize, MartianStruct)]
pub struct StageInputs {
    pub unmerged_gem_well_files: Vec<GemWellFiles>,
}

#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct StageOutputs {
    pub merged_gem_well_files: GemWellFiles,
}

pub fn hard_link_files<P: Clone + MartianFileType + AsRef<Path>>(
    v: &mut Vec<P>,
    rover: &MartianRover,
) {
    *v = v
        .iter()
        .cloned()
        .map(|m| {
            hard_link_martianfile(m, rover).expect("Error hard linking file in merge_gem_wells")
        })
        .collect();
}

#[make_mro(volatile = strict)]
impl MartianMain for MergeGemWellFiles {
    type StageInputs = StageInputs;
    type StageOutputs = StageOutputs;

    fn main(&self, mut args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        // header and target_set_name can be left the same
        // read_chunks and alignments are to be combined over the gem wells

        // Merge the lists of files over GEM wells
        let mut merged = args.unmerged_gem_well_files[0].clone();
        for unmerged in &mut args.unmerged_gem_well_files[1..] {
            merged.alignments.append(&mut unmerged.alignments);
            merged.read_chunks.append(&mut unmerged.read_chunks);
            merged.bc_umi_info.append(&mut unmerged.bc_umi_info);
            merged.gem_groups.append(&mut unmerged.gem_groups);

            match (&mut merged.annotation_files, &mut unmerged.annotation_files) {
                (Some(m), Some(u)) => {
                    m.num_reads += u.num_reads;
                    m.files.append(&mut u.files);
                }
                (None, None) => {}
                _ => {
                    bail!("PD annotation_files were None for some gem wells and Some for others.");
                }
            }

            merged
                .per_barcode_metrics_shard
                .append(&mut unmerged.per_barcode_metrics_shard);
        }

        // Replace the files in the merged lists with a new hard link to play nice with Martian
        if let Some(m) = &mut merged.annotation_files {
            hard_link_files(&mut m.files, &rover);
        }
        hard_link_files(&mut merged.alignments, &rover);
        hard_link_files(&mut merged.bc_umi_info, &rover);
        hard_link_files(&mut merged.per_barcode_metrics_shard, &rover);

        // hard link bam and feature ref
        merged.slfe_feature_reference =
            hard_link_martianfile(merged.slfe_feature_reference, &rover)?;
        // PathBuf doesn't implement MartianFiletype...
        let bam_link: PathBuf = rover.make_path("bam_header");
        std::fs::hard_link(&merged.bam_header, &bam_link)
            .or_else(|_| std::fs::copy(&merged.bam_header, &bam_link).map(|_| ()))?;
        merged.bam_header = bam_link;

        Ok(StageOutputs {
            merged_gem_well_files: merged,
        })
    }
}
