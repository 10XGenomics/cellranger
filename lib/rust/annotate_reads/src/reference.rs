//
// Copyright (c) 2017 10x Genomics, Inc. All rights reserved.
//

use std::collections::HashMap;

use utils;

pub struct TranscriptIndex {
    pub transcript_genes:   HashMap<String, Gene>,
    pub transcript_lengths: HashMap<String, i64>,
}

impl TranscriptIndex {
    pub fn new(tab_file: &str) -> TranscriptIndex {
        let mut transcript_lengths = HashMap::new();
        let mut transcript_genes = HashMap::new();
        let tx_index: Vec<GeneTranscriptRow> = utils::load_tabular(tab_file, true);
        for tx in tx_index {
            transcript_lengths.insert(tx.transcript_id.clone(), tx.transcript_len);
            transcript_genes.insert(tx.transcript_id.clone(), Gene { id: tx.gene_id, name: tx.gene_name } );
        }
        return TranscriptIndex {
            transcript_genes:   transcript_genes,
            transcript_lengths: transcript_lengths,
        }
    }

    pub fn get_gene_from_transcript(&self, tx: &String) -> &Gene {
        self.transcript_genes.get(tx).unwrap()
    }

    pub fn get_transcript_length(&self, tx: &String) -> &i64 {
        self.transcript_lengths.get(tx).unwrap()
    }
}

#[derive(Deserialize)]
struct GeneTranscriptRow {
    transcript_id:  String,
    gene_id:        String,
    gene_name:      String,
    transcript_len: i64,
}

#[derive(Hash, Eq, PartialEq, Debug, Clone, Ord, PartialOrd)]
pub struct Gene {
    pub id:     String,
    pub name:   String,
}
