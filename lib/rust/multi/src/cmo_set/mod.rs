use anyhow::{bail, Result};
use cr_types::reference::feature_reference::{FeatureDef, FeatureType};
use cr_types::types::FeatureBarcodeType;
use cr_types::GenomeName;
use fastq_set::read_pair::WhichRead;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::io::Read;

pub const SC3P_CELLPLEX_SET_A: &str = include_str!("SC3P_CellPlex_SetA.csv");
pub const DEFAULT_CMO_SET: &str = SC3P_CELLPLEX_SET_A;

#[derive(Debug, Deserialize, Eq, PartialEq, Serialize)]
pub struct CmoDef {
    pub id: String,
    pub name: String,
    pub read: WhichRead,
    pub pattern: String,
    pub sequence: String,
    pub feature_type: FeatureBarcodeType,
}

impl CmoDef {
    pub fn into_feature_def(self, index: usize) -> FeatureDef {
        let CmoDef {
            id,
            name,
            read,
            pattern,
            sequence,
            feature_type,
        } = self;
        FeatureDef {
            index,
            id,
            name,
            genome: GenomeName::default(),
            sequence,
            pattern,
            read,
            feature_type: FeatureType::Barcode(feature_type),
            tags: HashMap::new(),
        }
    }

    pub fn is_equivalent(&self, fdef: &FeatureDef) -> bool {
        self.read == fdef.read && self.pattern == fdef.pattern && self.sequence == fdef.sequence
    }
}

pub fn load_cmo_set<R: Read>(reader: R) -> Result<Vec<CmoDef>> {
    let mut rdr = csv::ReaderBuilder::new()
        .trim(csv::Trim::All)
        .from_reader(reader);
    rdr.deserialize()
        .enumerate()
        .map(|(i, record)| {
            let record: CmoDef = record?;
            if record.feature_type != FeatureBarcodeType::Multiplexing {
                bail!(
                    "CMO definition {} must have feature_type \"Multiplexing Capture\"",
                    i + 1
                );
            }
            Ok(record)
        })
        .collect()
}

pub fn load_default_cmo_set() -> Result<Vec<CmoDef>> {
    load_cmo_set(DEFAULT_CMO_SET.as_bytes())
}
