//! MergePerSampleAnnotations stage code

use anyhow::Result;
use martian::prelude::*;
use martian_derive::{make_mro, MartianStruct};
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::{LazyFileTypeIO, LazyWrite};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use vdj_ann::annotate::ContigAnnotation;

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct MergePerSampleAnnotationsStageInputs {
    pub per_sample_annotations: HashMap<String, JsonFile<Vec<ContigAnnotation>>>,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct MergePerSampleAnnotationsStageOutputs {
    pub contig_annotations: JsonFile<Vec<ContigAnnotation>>,
}

pub struct MergePerSampleAnnotations;

#[make_mro]
impl MartianMain for MergePerSampleAnnotations {
    type StageInputs = MergePerSampleAnnotationsStageInputs;
    type StageOutputs = MergePerSampleAnnotationsStageOutputs;
    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        let contig_annotations: JsonFile<Vec<ContigAnnotation>> =
            rover.make_path("contig_annotations");
        let mut ann_writer = contig_annotations.lazy_writer()?;
        for (sample, annot) in args.per_sample_annotations {
            let contig_reader = annot.lazy_reader()?;
            for ann in contig_reader {
                let mut ann: ContigAnnotation = ann?;
                // ann.clonotype = None;
                // ann.info = Default::default();
                ann.sample = Some(sample.clone());
                ann_writer.write_item(&ann)?;
            }
        }
        ann_writer.finish()?;

        Ok(MergePerSampleAnnotationsStageOutputs { contig_annotations })
    }
}
