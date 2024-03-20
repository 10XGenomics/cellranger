//! Martian stage SUBSAMPLE_BARCODES

use anyhow::Result;
use cr_types::types::{BarcodeSetFormat, LibraryType};
use cr_types::BcCountFormat;
use itertools::Itertools;
use martian::{MartianMain, MartianRover};
use martian_derive::{make_mro, MartianStruct};
use martian_filetypes::{FileTypeRead, FileTypeWrite};
use metric::TxHashSet;
use rand::seq::SliceRandom;
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use serde::{Deserialize, Serialize};
use std::cmp::max;

const MIN_FRAC_BARCODES: f64 = 0.10;
const MIN_FRAC_READS: f64 = 0.10;
const MIN_TOTAL_READS: u32 = 10_000_000;
const RNG_SEED: [u8; 32] = [0; 32];

/// Select a random sample of barcodes.
/// Receive counts of reads from corrected barcodes and select a small random
/// subset of those barcodes that covers some minimum threshold amount of data.
pub struct SubsampleBarcodes;

#[derive(Deserialize, Clone, MartianStruct)]
pub struct StageInputs {
    pub corrected_barcode_counts: BcCountFormat,
}

#[derive(Serialize, Deserialize, Clone, MartianStruct)]
pub struct StageOutputs {
    pub barcode_subset: BarcodeSetFormat,
}

#[make_mro(mem_gb = 4, volatile = strict)]
impl MartianMain for SubsampleBarcodes {
    type StageInputs = StageInputs;
    type StageOutputs = StageOutputs;

    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        let gex_barcode_counts = args
            .corrected_barcode_counts
            .read()?
            .into_iter()
            .filter_map(|(k, v)| match k {
                LibraryType::Gex => Some(v),
                _ => None,
            })
            .exactly_one()
            .unwrap();

        let total_reads = gex_barcode_counts.raw_counts().sum::<i64>() as u32;
        let min_sampled_reads = max(
            MIN_TOTAL_READS,
            (MIN_FRAC_READS * total_reads as f64).ceil() as u32,
        );

        let total_barcodes = gex_barcode_counts.distribution().len();
        let min_sampled_barcodes = (MIN_FRAC_BARCODES * total_barcodes as f64).ceil() as usize;

        let mut rng: ChaCha20Rng = ChaCha20Rng::from_seed(RNG_SEED);
        let mut shuffled_barcodes: Vec<_> =
            gex_barcode_counts.distribution().keys().copied().collect();
        shuffled_barcodes.shuffle(&mut rng);

        let mut sampled_barcodes = TxHashSet::<_>::default();
        let mut num_sampled_reads = 0_u32;
        for bc in shuffled_barcodes {
            num_sampled_reads += gex_barcode_counts.get(&bc) as u32;
            sampled_barcodes.insert(bc);
            if sampled_barcodes.len() >= min_sampled_barcodes
                && num_sampled_reads >= min_sampled_reads
            {
                break;
            }
        }

        println!(
            "sampled {} corrected barcodes comprised of {} reads",
            sampled_barcodes.len(),
            num_sampled_reads
        );
        let barcode_subset_file: BarcodeSetFormat = rover.make_path("bc_subset");
        barcode_subset_file.write(&sampled_barcodes)?;

        Ok(StageOutputs {
            barcode_subset: barcode_subset_file,
        })
    }
}

#[cfg(test)]
pub mod tests {
    use super::*;
    use barcode::{Barcode, BcSeq};
    use martian::prelude::{MartianFileType, MartianStage};
    use metric::{SimpleHistogram, TxHashMap};
    use std::path::{Path, PathBuf};
    use tempfile;

    const ACGT: &[u8; 4] = b"ACGT";

    fn prep_path(path: &Path, subdir: &str) -> Result<PathBuf> {
        let sub_path = path.join(subdir);
        std::fs::create_dir(&sub_path)?;
        Ok(sub_path)
    }

    fn u32_to_seq(seq: u32) -> BcSeq {
        let mut bc = [0_u8; 8];
        for idx in 0..8 {
            bc[idx] = ACGT[((seq >> idx) & 3) as usize];
        }
        BcSeq::from_bytes(&bc)
    }

    fn test_stage(num_barcodes: u32) -> Result<bool> {
        let mut barcode_counts = TxHashMap::default();
        let mut histogram: SimpleHistogram<Barcode> = SimpleHistogram::default();
        for i in 0..num_barcodes {
            histogram.observe_by_owned(Barcode::with_seq(0, u32_to_seq(i), true), i + 1);
        }
        let total_num_reads = histogram.raw_counts().sum::<i64>() as u32;
        barcode_counts.insert(LibraryType::Gex, histogram.clone());

        let tmp_dir = tempfile::tempdir().unwrap();
        let tmp_dir = tmp_dir.path();
        let barcode_counts_file: BcCountFormat = BcCountFormat::new(tmp_dir, "barcode_counts");
        barcode_counts_file.write(&barcode_counts)?;

        let stage_args = StageInputs {
            corrected_barcode_counts: barcode_counts_file,
        };
        let stage_path = prep_path(tmp_dir, "subsample_barcodes")?;
        let stage_outs = SubsampleBarcodes.test_run(stage_path, stage_args)?;

        let barcode_subset = stage_outs.barcode_subset.read()?;

        let total_frac_barcodes =
            barcode_subset.len() as f64 / histogram.distribution().len() as f64;
        assert!(total_frac_barcodes >= MIN_FRAC_BARCODES);

        let total_reads_sampled = barcode_subset
            .iter()
            .map(|bc| histogram.get(bc) as u32)
            .sum::<u32>();
        let low_num_reads = total_num_reads < MIN_TOTAL_READS;
        assert!(total_reads_sampled >= MIN_TOTAL_READS || total_num_reads < MIN_TOTAL_READS);

        let total_frac_reads = total_reads_sampled as f64 / total_num_reads as f64;
        assert!(total_frac_reads >= MIN_FRAC_READS);

        Ok(low_num_reads)
    }

    #[test]
    fn test_high_read_count() -> Result<()> {
        // num_counts = len * (len + 1)/2 ==> sqrt(2 * num_counts) - 1 > len
        let num_barcodes = ((2.0 * MIN_TOTAL_READS as f64).sqrt() - 1.0).ceil() as u32;
        let was_low_num_reads = test_stage(num_barcodes)?;
        assert!(!was_low_num_reads);
        Ok(())
    }

    #[test]
    fn test_low_read_count() -> Result<()> {
        let num_barcodes = 10;
        let was_low_num_reads = test_stage(num_barcodes)?;
        assert!(was_low_num_reads);
        Ok(())
    }
}
