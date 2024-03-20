//! WriteConsensusTxt stage code

use crate::assigner::ProtoBinFile;
use crate::write_concat_ref_outs::{FastaFaiFile, FastaFile};
use amino::aa_seq;
use anyhow::Result;
use cr_types::clonotype::ClonotypeId;
use enclone_proto::proto_io::read_proto;
use martian::prelude::*;
use martian_derive::{make_mro, MartianStruct};
use martian_filetypes::tabular_file::CsvFile;
use serde::{Deserialize, Serialize};
use std::io::Write;
use string_utils::stringme;

#[derive(Clone, Debug, Serialize)]
pub struct ConsensusAnnotationCsvRow {
    clonotype_id: String,
    consensus_id: String,
    length: usize,
    chain: String,
    v_gene: String,
    d_gene: Option<String>,
    j_gene: String,
    c_gene: Option<String>,
    full_length: bool,
    productive: bool,
    fwr1: Option<String>,
    fwr1_nt: Option<String>,
    cdr1: Option<String>,
    cdr1_nt: Option<String>,
    fwr2: Option<String>,
    fwr2_nt: Option<String>,
    cdr2: Option<String>,
    cdr2_nt: Option<String>,
    fwr3: Option<String>,
    fwr3_nt: Option<String>,
    cdr3: String,
    cdr3_nt: String,
    fwr4: Option<String>,
    fwr4_nt: Option<String>,
    reads: usize,
    umis: usize,
    v_start: usize,
    v_end: usize,
    v_end_ref: usize,
    j_start: usize,
    j_start_ref: usize,
    j_end: usize,
    fwr1_start: Option<usize>,
    fwr1_end: Option<usize>,
    cdr1_start: Option<usize>,
    cdr1_end: Option<usize>,
    fwr2_start: Option<usize>,
    fwr2_end: Option<usize>,
    cdr2_start: Option<usize>,
    cdr2_end: Option<usize>,
    fwr3_start: Option<usize>,
    fwr3_end: Option<usize>,
    cdr3_start: usize,
    cdr3_end: usize,
    fwr4_start: Option<usize>,
    fwr4_end: Option<usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct WriteConsensusTxtStageInputs {
    pub sample_number: Option<usize>,
    pub enclone_output: ProtoBinFile,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct WriteConsensusTxtStageOutputs {
    pub consensus_fasta: FastaFile,
    pub consensus_fasta_fai: FastaFaiFile,
    pub consensus_annotations_csv: CsvFile<ConsensusAnnotationCsvRow>,
}

// This is our stage struct
pub struct WriteConsensusTxt;

#[make_mro(mem_gb = 4, threads = 1)]
impl MartianMain for WriteConsensusTxt {
    type StageInputs = WriteConsensusTxtStageInputs;
    type StageOutputs = WriteConsensusTxtStageOutputs;
    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        let enclone_outs = read_proto(args.enclone_output)?;

        // Generate consensus_annotations.csv.
        let consensus_ann_csv_file = rover.make_path("consensus_annotations.csv");
        let mut all_writer = csv::Writer::from_path(&consensus_ann_csv_file)?;
        let uref_items = &enclone_outs.universal_reference.items;
        for (i, x) in enclone_outs.clonotypes.iter().enumerate() {
            for (j, cl_chain) in x.chains.iter().enumerate() {
                let seq = &cl_chain.nt_sequence;
                let cdr3_start = cl_chain.cdr3_start as usize;
                let cdr3_end = cl_chain.cdr3_end as usize;
                let cdr3_nt = seq[cdr3_start..cdr3_end].to_vec();
                let mut reads: u32 = 0;
                let mut umis: u32 = 0;
                for ex_cl in &x.exact_clonotypes {
                    for chain_info in &ex_cl.chains {
                        if chain_info.index as usize == j {
                            let exchain = &chain_info.chain;
                            reads += exchain.read_counts.iter().sum::<u32>();
                            umis += exchain.umi_counts.iter().sum::<u32>();
                        }
                    }
                }
                let fwr1_region = cl_chain.fwr1_region();
                let cdr1_region = cl_chain.cdr1_region();
                let fwr2_region = cl_chain.fwr2_region();
                let cdr2_region = cl_chain.cdr2_region();
                let fwr3_region = cl_chain.fwr3_region();
                let fwr4_region = cl_chain.fwr4_region();

                let clonotype_id = ClonotypeId {
                    id: i + 1,
                    sample_number: args.sample_number,
                };
                let row = ConsensusAnnotationCsvRow {
                    clonotype_id: clonotype_id.to_string(),
                    consensus_id: clonotype_id.consensus_name(j + 1),
                    length: seq.len(),
                    chain: cl_chain.chain_type.clone(),
                    v_gene: uref_items[cl_chain.v_idx as usize].display_name.clone(),
                    d_gene: cl_chain
                        .d_idx
                        .map(|i| uref_items[i as usize].display_name.clone()),
                    j_gene: uref_items[cl_chain.j_idx as usize].display_name.clone(),
                    c_gene: cl_chain
                        .c_idx
                        .map(|i| uref_items[i as usize].display_name.clone()),
                    full_length: true,
                    productive: true,
                    fwr1: fwr1_region.as_ref().map(|r| &r.aa_seq).cloned(),
                    fwr1_nt: fwr1_region.as_ref().map(|r| &r.nt_seq).cloned(),
                    cdr1: cdr1_region.as_ref().map(|r| &r.aa_seq).cloned(),
                    cdr1_nt: cdr1_region.as_ref().map(|r| &r.nt_seq).cloned(),
                    fwr2: fwr2_region.as_ref().map(|r| &r.aa_seq).cloned(),
                    fwr2_nt: fwr2_region.as_ref().map(|r| &r.nt_seq).cloned(),
                    cdr2: cdr2_region.as_ref().map(|r| &r.aa_seq).cloned(),
                    cdr2_nt: cdr2_region.as_ref().map(|r| &r.nt_seq).cloned(),
                    fwr3: fwr3_region.as_ref().map(|r| &r.aa_seq).cloned(),
                    fwr3_nt: fwr3_region.as_ref().map(|r| &r.nt_seq).cloned(),
                    cdr3: stringme(&aa_seq(&cdr3_nt, 0)),
                    cdr3_nt: stringme(&cdr3_nt),
                    fwr4: fwr4_region.as_ref().map(|r| &r.aa_seq).cloned(),
                    fwr4_nt: fwr4_region.as_ref().map(|r| &r.nt_seq).cloned(),
                    reads: reads as usize,
                    umis: umis as usize,
                    v_start: cl_chain.v_start as usize,
                    v_end: cl_chain.v_end as usize,
                    v_end_ref: cl_chain.v_end_ref as usize,
                    j_start: cl_chain.j_start as usize,
                    j_start_ref: cl_chain.j_start_ref as usize,
                    j_end: cl_chain.j_end as usize,
                    fwr1_start: fwr1_region.as_ref().map(|r| r.start),
                    fwr1_end: fwr1_region.as_ref().map(|r| r.stop),
                    cdr1_start: cdr1_region.as_ref().map(|r| r.start),
                    cdr1_end: cdr1_region.as_ref().map(|r| r.stop),
                    fwr2_start: fwr2_region.as_ref().map(|r| r.start),
                    fwr2_end: fwr2_region.as_ref().map(|r| r.stop),
                    cdr2_start: cdr2_region.as_ref().map(|r| r.start),
                    cdr2_end: cdr2_region.as_ref().map(|r| r.stop),
                    fwr3_start: fwr3_region.as_ref().map(|r| r.start),
                    fwr3_end: fwr3_region.as_ref().map(|r| r.stop),
                    cdr3_start,
                    cdr3_end,
                    fwr4_start: fwr4_region.as_ref().map(|r| r.start),
                    fwr4_end: fwr4_region.as_ref().map(|r| r.stop),
                };
                all_writer.serialize(row)?;
            }
        }

        // Write consensus.fasta and consensus.fasta.fai.

        let consensus_fasta: FastaFile = rover.make_path("consensus.fasta");
        let consensus_fasta_fai: FastaFaiFile = rover.make_path("consensus.fasta.fai");
        let mut writer = consensus_fasta.buf_writer()?;
        let mut writer_fai = consensus_fasta_fai.buf_writer()?;
        let mut bytes_written = 0;
        for (i, clonotype) in enclone_outs.clonotypes.iter().enumerate() {
            let clonotype_id = ClonotypeId {
                id: i + 1,
                sample_number: args.sample_number,
            };
            for (j, chain) in clonotype.chains.iter().enumerate() {
                let record_name = clonotype_id.consensus_name(j + 1);
                writeln!(writer, ">{record_name}")?;
                bytes_written += record_name.len() + 2;
                let seq = &chain.nt_sequence;
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
        Ok(WriteConsensusTxtStageOutputs {
            consensus_fasta,
            consensus_fasta_fai,
            consensus_annotations_csv: consensus_ann_csv_file,
        })
    }
}
