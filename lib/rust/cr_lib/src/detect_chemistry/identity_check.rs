//!
//! A module to check that R1 != R2 (detect accidental duplication of fastqs)
//!

use super::chemistry_filter::DetectChemistryUnit;
use anyhow::{bail, Result};
use fastq_set::read_pair::{ReadPair, ReadPart, WhichRead};
use metric::{TxHashMap, TxHasher};
use std::collections::hash_map::Entry;
use std::hash::Hasher;

pub(crate) fn check_read_identity(
    unit: &DetectChemistryUnit,
    read_pairs: &[ReadPair],
) -> Result<(u64, u64)> {
    use ReadPart::{Header, Qual, Seq};
    use WhichRead::{R1, R2};
    let mut r1_hasher = TxHasher::default();
    let mut r2_hasher = TxHasher::default();
    for read_pair in read_pairs {
        let _ = read_pair.get(R1, Header).map(|x| r1_hasher.write(x));
        let _ = read_pair.get(R1, Seq).map(|x| r1_hasher.write(x));
        let _ = read_pair.get(R1, Qual).map(|x| r1_hasher.write(x));
        let _ = read_pair.get(R2, Header).map(|x| r2_hasher.write(x));
        let _ = read_pair.get(R2, Seq).map(|x| r2_hasher.write(x));
        let _ = read_pair.get(R2, Qual).map(|x| r2_hasher.write(x));
    }
    let r1_hash = r1_hasher.finish();
    let r2_hash = r2_hasher.finish();
    if r1_hash == r2_hash {
        match unit.group.as_ref() {
            Some(sample_name) => bail!(
                "R1 and R2 reads identical in sample \"{}\" at \"{}\"",
                sample_name,
                unit.read_path.display()
            ),
            None => bail!(
                "R1 and R2 reads identical at \"{}\"",
                unit.read_path.display()
            ),
        };
    }
    Ok((r1_hash, r2_hash))
}

pub(crate) fn check_fastq_identity(units: &[DetectChemistryUnit]) -> Result<()> {
    let mut hashes = TxHashMap::default();
    for unit in units {
        for pair_hash in unit.check_read_identity()? {
            match hashes.entry(pair_hash) {
                Entry::Occupied(o) => {
                    bail!("Duplicate FASTQs found between {} and {}", unit, o.get());
                }
                Entry::Vacant(v) => {
                    v.insert(unit);
                }
            }
        }
    }
    Ok(())
}
