//! WriteConsensusBam stage code

use crate::assigner::{
    bam_record_from_align, replace_sam_by_indexed_bam, BamBaiFile, BamFile, ProtoBinFile, SamFile,
    QUALITY_OFFSET,
};
use anyhow::Result;
use bio::alignment::pairwise::Aligner;
use cr_types::clonotype::ClonotypeId;
use enclone_proto::proto_io::read_proto;
use martian::prelude::*;
use martian_derive::{make_mro, MartianStruct};
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::LazyFileTypeIO;
use rayon::prelude::*;
use rust_htslib::bam;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::time::Instant;
use vdj_ann::annotate::ContigAnnotation;
use vdj_asm_utils::bam_utils::add_ref_to_bam_header;

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct WriteConsensusBamStageInputs {
    pub sample_id: Option<String>,
    pub enclone_output: ProtoBinFile,
    pub all_contig_annotations_json: JsonFile<Vec<ContigAnnotation>>,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct WriteConsensusBamStageOutputs {
    pub consensus_bam: BamFile,
    pub consensus_bam_bai: BamBaiFile,
}

// This is our stage struct
pub struct WriteConsensusBam;

#[make_mro(mem_gb = 4, threads = 4)]
impl MartianMain for WriteConsensusBam {
    type StageInputs = WriteConsensusBamStageInputs;
    type StageOutputs = WriteConsensusBamStageOutputs;
    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        // Get the contig seqs and quals from all_contig_annotations.json.
        let mut to_seq = HashMap::<String, String>::new(); // Key: Contig name, Value: Contig seq
        let mut to_quals = HashMap::<String, Vec<u8>>::new(); // Key: Contig name, Value: Contig qual

        for ann in args.all_contig_annotations_json.lazy_reader()? {
            let ann: ContigAnnotation = ann?;

            if ann.is_cell && ann.is_productive() && ann.high_confidence {
                to_seq.insert(ann.contig_name.clone(), ann.sequence.clone());
                to_quals.insert(
                    ann.contig_name.clone(),
                    ann.quals.bytes().map(|q| q - QUALITY_OFFSET).collect(),
                );
            }
        }

        let enclone_outs = read_proto(&args.enclone_output)?;
        // Create consensus.bam and consensus.bam.bai.

        let mut to_tid = HashMap::<(usize, usize), usize>::new();
        let mut count = 0;
        for (i, x) in enclone_outs.clonotypes.iter().enumerate() {
            for j in 0..x.chains.len() {
                to_tid.insert((i, j), count);
                count += 1;
            }
        }

        let t = Instant::now();
        let mut header = bam::header::Header::new();
        for (i, x) in enclone_outs.clonotypes.iter().enumerate() {
            for j in 0..x.chains.len() {
                let seq = &x.chains[j].nt_sequence;
                add_ref_to_bam_header(
                    &mut header,
                    &ClonotypeId {
                        id: i + 1,
                        sample_id: args.sample_id.as_deref(),
                    }
                    .consensus_name(j + 1),
                    seq.len(),
                );
            }
        }
        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let mut header_rec = bam::header::HeaderRecord::new(b"PG");
        header_rec.push_tag(b"ID", "cellranger");
        header.push_record(&header_rec);
        let consensus_bam_file: BamFile = rover.make_path("consensus.bam");
        let consensus_bam_bai_file: BamBaiFile = rover.make_path("consensus.bam.bai");
        let consensus_sam_file: &SamFile = &rover.make_path("consensus.sam");
        {
            let mut out_sam =
                bam::Writer::from_path(consensus_sam_file, &header, bam::Format::Sam)?;
            let mut results = Vec::<(usize, Vec<bam::Record>)>::new();
            for i in 0..enclone_outs.clonotypes.len() {
                results.push((i, Vec::<bam::Record>::new()));
            }
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(rover.get_threads())
                .build()
                .unwrap();
            pool.install(|| {
                results.par_iter_mut().for_each(|res| {
                    let i = res.0;
                    let x = &enclone_outs.clonotypes[i];
                    let clonotype_chains = &x.chains;
                    for j in 0..x.exact_clonotypes.len() {
                        let barcodes: &Vec<String> = &x.exact_clonotypes[j].cell_barcodes;
                        for k in 0..x.exact_clonotypes[j].chains.len() {
                            let index = x.exact_clonotypes[j].chains[k].index as usize;
                            let refid = to_tid[&(i, index)];
                            let seq = &clonotype_chains[index].nt_sequence;
                            let contig_ids = &x.exact_clonotypes[j].chains[k].chain.contig_ids;
                            for (l, id) in contig_ids.iter().enumerate() {
                                let tig_seq = &to_seq[id];
                                let tig_quals = &to_quals[id];
                                let mut aligner = Aligner::new(-6, -1, &score);
                                let al = aligner.semiglobal(tig_seq.as_bytes(), seq);
                                let rec = bam_record_from_align(
                                    &al,
                                    id,
                                    refid,
                                    tig_seq.as_bytes(),
                                    tig_quals,
                                    &barcodes[l],
                                );
                                res.1.push(rec);
                            }
                        }
                    }
                });
            });
            for result in results {
                for rec in result.1 {
                    let _ = out_sam.write(&rec);
                }
            }
        }
        println!(
            "used {:.2} seconds building consensus sam",
            t.elapsed().as_secs_f64()
        );

        // Make the indexed bam.

        replace_sam_by_indexed_bam(consensus_sam_file);

        Ok(WriteConsensusBamStageOutputs {
            consensus_bam: consensus_bam_file,
            consensus_bam_bai: consensus_bam_bai_file,
        })
    }
}
