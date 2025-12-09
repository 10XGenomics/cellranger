//! Martian stage WRITE_CONCAT_REF_OUTS creates concat_ref.bam and concat_ref.bam.bai.
#![expect(missing_docs)]

use crate::assigner::{
    BamBaiFile, BamFile, ProtoBinFile, QUALITY_OFFSET, SamFile, bam_record_from_align,
    replace_sam_by_indexed_bam,
};
use anyhow::Result;
use bio::alignment::pairwise::Aligner;
use cr_types::clonotype::ClonotypeId;
use enclone_proto::proto_io::read_proto;
use martian::prelude::*;
use martian_derive::{MartianStruct, make_mro, martian_filetype};
use martian_filetypes::LazyFileTypeIO;
use martian_filetypes::json_file::JsonFile;
use rayon::prelude::*;
use rust_htslib::bam;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::io::Write;
use vdj_ann::annotate::ContigAnnotation;
use vdj_asm_utils::bam_utils::add_ref_to_bam_header;

martian_filetype!(FastaFile, "fasta");
martian_filetype!(FastaFaiFile, "fasta.fai");

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct WriteConcatRefOutsStageInputs {
    pub sample_id: Option<String>,
    pub enclone_output: ProtoBinFile,
    pub all_contig_annotations_json: JsonFile<Vec<ContigAnnotation>>,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct WriteConcatRefOutsStageOutputs {
    pub concat_ref_bam: BamFile,
    pub concat_ref_bam_bai: BamBaiFile,
    pub concat_ref_fasta: FastaFile,
    pub concat_ref_fasta_fai: FastaFaiFile,
}

// This is our stage struct
pub struct WriteConcatRefOuts;

const CONCAT_REF: &str = "concat_ref";

#[make_mro(mem_gb = 12, threads = 4)]
impl MartianMain for WriteConcatRefOuts {
    type StageInputs = WriteConcatRefOutsStageInputs;
    type StageOutputs = WriteConcatRefOutsStageOutputs;
    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        // Get the contig seqs and quals from all_contig_annotations.json.
        let mut to_seq = HashMap::<String, String>::new(); // Key: Contig name, Value: Contig seq
        let mut to_quals = HashMap::<String, Vec<u8>>::new(); // Key: Contig name, Value: Contig qual

        // Define contig_ann reader
        let contig_reader = args.all_contig_annotations_json.lazy_reader()?;
        for ann in contig_reader {
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

        let mut bytes_written = 0;
        let concat_ref_fasta_file: FastaFile = rover.make_path(CONCAT_REF);
        let concat_ref_fasta_fai_file: FastaFaiFile = rover.make_path(CONCAT_REF);
        let mut writer = concat_ref_fasta_file.buf_writer()?;
        let mut writer_fai = concat_ref_fasta_fai_file.buf_writer()?;
        let mut header = bam::header::Header::new();
        let mut to_tid = HashMap::<(usize, usize), usize>::new();
        let mut count = 0;
        for (i, x) in enclone_outs.clonotypes.iter().enumerate() {
            let clonotype_id = ClonotypeId {
                id: i + 1,
                sample_id: args.sample_id.as_deref(),
            };
            for j in 0..x.chains.len() {
                let record_name = clonotype_id.concat_ref_name(j + 1);
                let seq = &x.chains[j].universal_reference;
                add_ref_to_bam_header(&mut header, &record_name, seq.len());
                to_tid.insert((i, j), count);
                count += 1;

                writeln!(writer, ">{record_name}")?;
                bytes_written += record_name.len() + 2;
                writeln!(writer, "{}", std::str::from_utf8(seq).unwrap())?;
                writeln!(
                    writer_fai,
                    "{}\t{}\t{}\t{}\t{}",
                    record_name,
                    seq.len(),
                    bytes_written,
                    seq.len(),
                    seq.len() + 1
                )?;
                bytes_written += seq.len() + 1;
            }
        }
        let mut header_rec = bam::header::HeaderRecord::new(b"PG");
        header_rec.push_tag(b"ID", "cellranger");
        header.push_record(&header_rec);
        let concat_ref_bam_file: BamFile = rover.make_path(CONCAT_REF);
        let concat_ref_bam_bai_file: BamBaiFile = rover.make_path(CONCAT_REF);
        let concat_ref_sam_file: &SamFile = &rover.make_path(CONCAT_REF);
        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };

        // Using "let pool = ..." as below because the following commented out line did not work,
        // for unknown reasons:
        // let _ = rayon::ThreadPoolBuilder::new().num_threads(nthreads).build_global();

        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(rover.get_threads())
            .build()
            .unwrap();
        pool.install(|| {
            let mut out_sam =
                bam::Writer::from_path(concat_ref_sam_file, &header, bam::Format::Sam).unwrap();

            let mut results = Vec::<(usize, Vec<bam::Record>)>::new();
            for i in 0..enclone_outs.clonotypes.len() {
                results.push((i, Vec::<bam::Record>::new()));
            }
            results.par_iter_mut().for_each(|res| {
                let i = res.0;
                let x = &enclone_outs.clonotypes[i];
                let clonotype_id = ClonotypeId {
                    id: i + 1,
                    sample_id: args.sample_id.as_deref(),
                };
                for j in 0..x.exact_clonotypes.len() {
                    let barcodes: &Vec<String> = &x.exact_clonotypes[j].cell_barcodes;
                    for k in 0..x.exact_clonotypes[j].chains.len() {
                        let index = x.exact_clonotypes[j].chains[k].index as usize;
                        let refid = to_tid[&(i, index)];
                        let chain = &x.exact_clonotypes[j].chains[k].chain;
                        let seq = &chain.nt_sequence;
                        let concat = &x.chains[index].universal_reference;

                        // The following is the line that one needs to declare that no quality
                        // scores are provided.  This can be empirically demonstrated by
                        // uncommenting the line that deletes the sam file.  Then you'll see a * in
                        // in the sam file where the quality scores would be, as per the sam spec.

                        let quals = vec![255_u8; seq.len()];

                        // And now for the rest.

                        let id = clonotype_id.consensus_name(index + 1);
                        let mut aligner = Aligner::new(-6, -1, &score);
                        let al = aligner.semiglobal(seq, concat);
                        let rec = bam_record_from_align(&al, &id, refid, seq, &quals, "");
                        res.1.push(rec);
                        let contig_ids: &Vec<String> = &chain.contig_ids;
                        for (l, id) in contig_ids.iter().enumerate() {
                            let tig_seq = &to_seq[id];
                            let tig_quals = &to_quals[id];
                            let mut aligner = Aligner::new(-6, -1, &score);
                            let al = aligner.semiglobal(tig_seq.as_bytes(), concat);
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
            for result in results {
                for rec in result.1 {
                    let _ = out_sam.write(&rec);
                }
            }
        });

        // Make the indexed bam.

        replace_sam_by_indexed_bam(concat_ref_sam_file);

        Ok(WriteConcatRefOutsStageOutputs {
            concat_ref_bam: concat_ref_bam_file,
            concat_ref_bam_bai: concat_ref_bam_bai_file,
            concat_ref_fasta: concat_ref_fasta_file,
            concat_ref_fasta_fai: concat_ref_fasta_fai_file,
        })
    }
}
