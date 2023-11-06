//! WriteClonotypeOuts stage code

use crate::assigner::ProtoBinFile;
use amino::aa_seq;
use anyhow::Result;
use enclone_proto::proto_io::read_proto;
use enclone_proto::types::InvariantTCellAnnotation;
use martian::prelude::*;
use martian_derive::{make_mro, MartianStruct};
use martian_filetypes::tabular_file::CsvFile;
use martian_filetypes::{LazyFileTypeIO, LazyWrite};
use serde::{Deserialize, Serialize};
use std::str;
use vdj_reference::VdjReceptor;

#[derive(Clone, Debug, Serialize, Deserialize, Default)]
pub struct ClonotypesCsvRow {
    clonotype_id: String,
    frequency: usize,
    proportion: f64,
    cdr3s_aa: String,
    cdr3s_nt: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    inkt_evidence: Option<String>, // Will be None for B cells
    #[serde(skip_serializing_if = "Option::is_none")]
    mait_evidence: Option<String>, // Will be None for B cells
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct WriteClonotypeOutsStageInputs {
    pub receptor: VdjReceptor,
    pub enclone_output: ProtoBinFile,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct WriteClonotypeOutsStageOutputs {
    pub clonotypes_csv: CsvFile<ClonotypesCsvRow>,
}

// This is our stage struct
pub struct WriteClonotypeOuts;

pub fn invariant_evidence_display(evidences: &[InvariantTCellAnnotation]) -> String {
    let chain_evidence = |chain, gene, junction| -> Option<String> {
        match (gene, junction) {
            (true, true) => Some(format!("{chain}:gene+junction")),
            (true, false) => Some(format!("{chain}:gene")),
            (false, true) => Some(format!("{chain}:junction")),
            (false, false) => None,
        }
    };

    let mut result = Vec::new();

    let alpha_gene = evidences.iter().any(|e| e.alpha_chain_gene_match);
    let alpha_junction = evidences.iter().any(|e| e.alpha_chain_junction_match);
    if let Some(e) = chain_evidence("TRA", alpha_gene, alpha_junction) {
        result.push(e)
    }

    let beta_gene = evidences.iter().any(|e| e.beta_chain_gene_match);
    let beta_junction = evidences.iter().any(|e| e.beta_chain_junction_match);
    if let Some(e) = chain_evidence("TRB", beta_gene, beta_junction) {
        result.push(e)
    }

    result.join(";")
}

#[make_mro(mem_gb = 8)]
impl MartianMain for WriteClonotypeOuts {
    type StageInputs = WriteClonotypeOutsStageInputs;
    type StageOutputs = WriteClonotypeOutsStageOutputs;
    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        let enclone_outs = read_proto(&args.enclone_output)?;

        // Write clonotypes.csv.  If there are no clonotypes, csv::Writer does not write the
        // header line.  This is not the desired behavior, so we have separate code for that case.
        let outputs = WriteClonotypeOutsStageOutputs {
            clonotypes_csv: rover.make_path("clonotypes.csv"),
        };
        let mut clonotypes_csv_writer = outputs.clonotypes_csv.lazy_writer()?;

        if enclone_outs.clonotypes.is_empty() {
            clonotypes_csv_writer.write_header()?;
        } else {
            let ncells: u32 = enclone_outs.clonotypes.iter().map(|cl| cl.frequency).sum();
            for (i, clonotype) in enclone_outs.clonotypes.iter().enumerate() {
                let mut cdr3s_aa = Vec::new();
                let mut cdr3s_nt = Vec::new();
                for chain in &clonotype.chains {
                    let cdr3_nt =
                        &chain.nt_sequence[chain.cdr3_start as usize..chain.cdr3_end as usize];
                    cdr3s_nt.push(format!(
                        "{}:{}",
                        chain.chain_type,
                        str::from_utf8(cdr3_nt).unwrap()
                    ));
                    cdr3s_aa.push(format!(
                        "{}:{}",
                        chain.chain_type,
                        str::from_utf8(&aa_seq(cdr3_nt, 0)).unwrap()
                    ));
                }
                let (inkt_evidence, mait_evidence) = match args.receptor {
                    VdjReceptor::TR => {
                        let mut inkt_evidences = Vec::new();
                        let mut mait_evidences = Vec::new();
                        for ex_cl in &clonotype.exact_clonotypes {
                            inkt_evidences.push(ex_cl.inkt_evidence.clone());
                            mait_evidences.push(ex_cl.mait_evidence.clone());
                        }
                        (
                            Some(invariant_evidence_display(&inkt_evidences)),
                            Some(invariant_evidence_display(&mait_evidences)),
                        )
                    }
                    VdjReceptor::TRGD => (None, None),
                    VdjReceptor::IG => (None, None),
                };

                let row = ClonotypesCsvRow {
                    clonotype_id: format!("clonotype{}", i + 1),
                    frequency: clonotype.frequency as usize,
                    proportion: clonotype.frequency as f64 / ncells as f64,
                    cdr3s_aa: cdr3s_aa.join(";"),
                    cdr3s_nt: cdr3s_nt.join(";"),
                    inkt_evidence,
                    mait_evidence,
                };
                clonotypes_csv_writer.write_item(&row)?;
            }
        }

        Ok(outputs)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_invariant_evidence() {
        assert_eq!(
            invariant_evidence_display(&[InvariantTCellAnnotation {
                alpha_chain_gene_match: false,
                alpha_chain_junction_match: false,
                beta_chain_gene_match: false,
                beta_chain_junction_match: false,
            }]),
            ""
        );
        assert_eq!(
            invariant_evidence_display(&[InvariantTCellAnnotation {
                alpha_chain_gene_match: false,
                alpha_chain_junction_match: false,
                beta_chain_gene_match: true,
                beta_chain_junction_match: false,
            }]),
            "TRB:gene"
        );

        assert_eq!(
            invariant_evidence_display(&[InvariantTCellAnnotation {
                alpha_chain_gene_match: true,
                alpha_chain_junction_match: false,
                beta_chain_gene_match: false,
                beta_chain_junction_match: false,
            }]),
            "TRA:gene"
        );

        assert_eq!(
            invariant_evidence_display(&[
                InvariantTCellAnnotation {
                    alpha_chain_gene_match: true,
                    alpha_chain_junction_match: false,
                    beta_chain_gene_match: false,
                    beta_chain_junction_match: false,
                },
                InvariantTCellAnnotation {
                    alpha_chain_gene_match: true,
                    alpha_chain_junction_match: true,
                    beta_chain_gene_match: false,
                    beta_chain_junction_match: false,
                }
            ]),
            "TRA:gene+junction"
        );

        assert_eq!(
            invariant_evidence_display(&[
                InvariantTCellAnnotation {
                    alpha_chain_gene_match: true,
                    alpha_chain_junction_match: false,
                    beta_chain_gene_match: false,
                    beta_chain_junction_match: false,
                },
                InvariantTCellAnnotation {
                    alpha_chain_gene_match: true,
                    alpha_chain_junction_match: false,
                    beta_chain_gene_match: true,
                    beta_chain_junction_match: true,
                }
            ]),
            "TRA:gene;TRB:gene+junction"
        );
    }
}
