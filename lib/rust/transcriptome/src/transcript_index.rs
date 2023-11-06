use crate::transcriptome::Transcriptome;
use anyhow::Result;
use martian_filetypes::tabular_file::TsvFile;
use martian_filetypes::LazyFileTypeIO;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TranscriptIndex {
    pub transcript_genes: HashMap<String, Gene>,
    pub transcript_lengths: HashMap<String, i64>,
}

#[derive(Hash, Eq, PartialEq, Debug, Clone, Ord, PartialOrd, Serialize, Deserialize)]
pub struct Gene {
    pub id: String,
    pub name: String,
}

// A summary representation of the Gene/Transcript table. Used in previous versions of CR.
// retained here for testing purposes
impl TranscriptIndex {
    pub fn new(tab_file: &TsvFile<GeneTranscriptRow>) -> Result<TranscriptIndex> {
        let mut transcript_lengths = HashMap::new();
        let mut transcript_genes = HashMap::new();
        for tx in tab_file.lazy_reader()? {
            let tx = tx?;
            transcript_lengths.insert(tx.transcript_id.clone(), tx.transcript_length);
            transcript_genes.insert(
                tx.transcript_id.clone(),
                Gene {
                    id: tx.gene_id,
                    name: tx.gene_name,
                },
            );
        }

        Ok(TranscriptIndex {
            transcript_genes,
            transcript_lengths,
        })
    }

    pub fn get_gene_from_transcript(&self, tx: &str) -> &Gene {
        self.transcript_genes
            .get(tx)
            .unwrap_or_else(|| panic!("Error: No corresponding gene for transcript: {tx}"))
    }

    pub fn get_transcript_length(&self, tx: &str) -> i64 {
        self.transcript_lengths[tx]
    }

    pub fn from_transcriptome(txome: &Transcriptome) -> TranscriptIndex {
        let mut transcript_lengths = HashMap::new();
        let mut transcript_genes = HashMap::new();

        for tx in &txome.transcripts {
            transcript_lengths.insert(tx.id.clone(), tx.len() as i64);

            let g = &txome.genes[tx.gene_idx.0 as usize];
            let g = Gene {
                id: g.id.clone(),
                name: g.name.clone(),
            };

            transcript_genes.insert(tx.id.clone(), g);
        }

        TranscriptIndex {
            transcript_genes,
            transcript_lengths,
        }
    }
}

#[derive(Deserialize, Serialize)]
pub struct GeneTranscriptRow {
    transcript_id: String,
    gene_id: String,
    gene_name: String,
    transcript_length: i64,
}
