//! Detect accidental duplication of FASTQ files.
#![deny(missing_docs)]

use super::chemistry_filter::DetectChemistryUnit;
use anyhow::{Result, bail};
use fastq_set::read_pair::{ReadPair, ReadPart, WhichRead};
use metric::{TxHashMap, TxHasher};
use std::collections::hash_map::Entry;
use std::hash::Hasher;
use std::path::Path;

pub(crate) fn check_read_identity(
    unit: &DetectChemistryUnit,
    read_pairs: &[ReadPair],
) -> Result<(u64, u64)> {
    use ReadPart::{Header, Qual, Seq};
    use WhichRead::{R1, R2};
    let mut r1_hasher = TxHasher::hasher();
    let mut r2_hasher = TxHasher::hasher();
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
    fn format_unit_fastq(
        unit: &DetectChemistryUnit,
        fastq: &fastq_set::read_pair_iter::InputFastqs,
    ) -> String {
        format!(
            "{}\"{}\"",
            match unit.group {
                Some(ref g) => format!("{g} in "),
                None => String::new(),
            },
            Path::new(&fastq.r1).parent().unwrap().display(),
        )
    }
    let mut hashes = TxHashMap::default();
    for unit in units {
        for (fastq, pair_hash) in unit.check_read_identity()? {
            match hashes.entry(pair_hash) {
                Entry::Vacant(v) => {
                    v.insert((unit, fastq));
                }
                Entry::Occupied(o) => {
                    let (o_unit, o_fastq) = o.get();
                    bail!(
                        "duplicate FASTQs found between {} and {}",
                        format_unit_fastq(unit, fastq),
                        format_unit_fastq(o_unit, o_fastq),
                    );
                }
            }
        }
    }
    Ok(())
}
