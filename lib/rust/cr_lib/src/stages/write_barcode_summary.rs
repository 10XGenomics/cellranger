//! Martian stage WRITE_BARCODE_SUMMARY
//! Write the file barcode_summary.h5.

use crate::align_metrics::{
    BarcodeKind, BarcodeMetrics, GenomeMapping, LibFeatThenBarcodeOrder, MappingRegion,
    MULTI_GENOME,
};
use crate::types::FeatureReferenceFormat;
use crate::BarcodeMetricsShardFile;
use anyhow::{bail, Result};
use cr_h5::count_matrix::write_barcodes_column;
use cr_types::barcode_index::BarcodeIndex;
use cr_types::reference::feature_reference::{FeatureReference, FeatureType};
use cr_types::{BarcodeIndexFormat, GenomeName, H5File, LibraryType};
use itertools::zip_eq;
use martian::prelude::{MartianMain, MartianRover};
use martian_derive::{make_mro, MartianStruct};
use martian_filetypes::FileTypeRead;
use metric::join_metric_name;
use serde::{Deserialize, Serialize};
use shardio::UnsortedShardReader;
use std::collections::HashMap;
use std::ops::AddAssign;
use std::path::Path;

/// The Martian stage inputs.
#[derive(Clone, Deserialize, MartianStruct)]
pub struct WriteBarcodeSummaryStageInputs {
    /// Metrics computed per barcode sorted by (library type, barcode)
    pub per_barcode_metrics: Vec<BarcodeMetricsShardFile>,
    pub feature_reference: FeatureReferenceFormat,
    pub barcode_index: BarcodeIndexFormat,
}

/// The Martian stage outputs.
#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct WriteBarcodeSummaryStageOutputs {
    pub barcode_summary: H5File,
}

/// Output the barcode summary HDF5 file `barcode_summary.h5`.
pub struct WriteBarcodeSummary;

/// A count of reads and molecules per barcode for each genome.
struct BarcodeCounts {
    /// The number of reads per barcode.
    conf_mapped_read_counts: Vec<u32>,
    /// The number of molecules per barcode.
    umi_counts: Vec<u32>,
    // Half-mapped probe reads (MAPQ = 1)
    half_mapped_counts: Vec<u32>,
    // Split-mapped probe reads (MAPQ = 3)
    split_mapped_counts: Vec<u32>,
}

impl BarcodeCounts {
    /// Return a new BarcodeCounts for the specified number of barcodes.
    fn with_size(num_barcodes: usize) -> Self {
        BarcodeCounts {
            conf_mapped_read_counts: vec![0; num_barcodes],
            half_mapped_counts: vec![0; num_barcodes],
            split_mapped_counts: vec![0; num_barcodes],
            umi_counts: vec![0; num_barcodes],
        }
    }
}

/// Add one vector to another.
fn add_assign(lhs: &mut [u32], rhs: &[u32]) {
    for (x, y) in zip_eq(lhs, rhs) {
        *x += y;
    }
}

impl AddAssign<&Self> for BarcodeCounts {
    /// Add one BarcodeCounts to another.
    fn add_assign(&mut self, rhs: &Self) {
        add_assign(
            &mut self.conf_mapped_read_counts,
            &rhs.conf_mapped_read_counts,
        );
        add_assign(&mut self.umi_counts, &rhs.umi_counts);
        add_assign(&mut self.half_mapped_counts, &rhs.half_mapped_counts);
        add_assign(&mut self.split_mapped_counts, &rhs.split_mapped_counts);
    }
}

struct BarcodeReads {
    // Reads that are assigned to this barcode before or after correction
    sequenced_reads: Vec<u32>,
    // Reads that are assigned to this barcode after correction
    // Will be a subset of the reads couted towards `sequenced_reads`
    barcode_corrected_sequenced_reads: Vec<u32>,
    unmapped_reads: Vec<u32>,
}
impl BarcodeReads {
    fn with_size(num_barcodes: usize) -> Self {
        BarcodeReads {
            sequenced_reads: vec![0; num_barcodes],
            barcode_corrected_sequenced_reads: vec![0; num_barcodes],
            unmapped_reads: vec![0; num_barcodes],
        }
    }
}

fn checked_add(sum: &mut u32, add: i64) {
    let new_sum = sum
        .checked_add(u32::try_from(add).unwrap())
        .expect("Overflow when adding counts");
    *sum = new_sum;
}

struct BarcodeSummaryData {
    /// Number of reads/umis per barcode for each genome.
    genome_barcode_counts: HashMap<(FeatureType, GenomeName), BarcodeCounts>,
    read_counts: HashMap<FeatureType, BarcodeReads>,
}

impl BarcodeSummaryData {
    fn new(
        per_barcode_metrics: &[BarcodeMetricsShardFile],
        feature_ref: &FeatureReference,
        barcode_index: &BarcodeIndex,
    ) -> Result<Self> {
        let mut genome_barcode_counts = HashMap::new();
        let mut read_counts = HashMap::new();
        for feature_def in &feature_ref.feature_defs {
            let genome = if feature_def.genome.is_empty() {
                MULTI_GENOME.into()
            } else {
                feature_def.genome.clone()
            };
            genome_barcode_counts
                .entry((feature_def.feature_type, genome))
                .or_insert_with(|| BarcodeCounts::with_size(barcode_index.len()));
            read_counts
                .entry(feature_def.feature_type)
                .or_insert_with(|| BarcodeReads::with_size(barcode_index.len()));
        }

        let metrics_reader =
            UnsortedShardReader::<BarcodeMetrics, LibFeatThenBarcodeOrder>::open_set(
                per_barcode_metrics,
            )?;

        fn update_counts(
            genome_barcode_counts: &mut HashMap<(FeatureType, GenomeName), BarcodeCounts>,
            key: &(FeatureType, GenomeName),
            index: usize,
            umi_counts: i64,
            read_counts: i64,
        ) {
            if let Some(entry) = genome_barcode_counts.get_mut(key) {
                checked_add(&mut entry.umi_counts[index], umi_counts);
                checked_add(&mut entry.conf_mapped_read_counts[index], read_counts);
            }
        }

        for barcode_metrics in metrics_reader {
            let BarcodeMetrics {
                barcode,
                library_type,
                metrics,
            } = barcode_metrics?;
            if let BarcodeKind::Valid(bc) = barcode {
                let index = barcode_index.get_index(&bc);
                let feature_type = match library_type {
                    LibraryType::Gex => {
                        for (genome, genome_metrics) in metrics.per_genome_name {
                            update_counts(
                                &mut genome_barcode_counts,
                                &(FeatureType::Gene, genome),
                                index,
                                genome_metrics.umi_counts.count(),
                                genome_metrics.usable_reads.count(),
                            );
                        }
                        for (GenomeMapping { genome, .. }, genome_mapping_metrics) in metrics
                            .per_genome_mapping
                            .into_iter()
                            .filter(|(k, _)| k.region == MappingRegion::Transcriptome)
                        {
                            // We want to ignore the "multi" genome here
                            if let Some(entry) =
                                genome_barcode_counts.get_mut(&(FeatureType::Gene, genome))
                            {
                                checked_add(
                                    &mut entry.half_mapped_counts[index],
                                    genome_mapping_metrics.half_mapped.count(),
                                );
                                checked_add(
                                    &mut entry.split_mapped_counts[index],
                                    genome_mapping_metrics.split_mapped.count(),
                                );
                            }
                        }
                        FeatureType::Gene
                    }
                    LibraryType::FeatureBarcodes(feature_type) => {
                        // No "genome" in other features but we use the key "multi"
                        update_counts(
                            &mut genome_barcode_counts,
                            &(FeatureType::Barcode(feature_type), MULTI_GENOME.into()),
                            index,
                            metrics.umi_counts,
                            metrics.usable_reads,
                        );
                        FeatureType::Barcode(feature_type)
                    }
                    LibraryType::Vdj(_) | LibraryType::Atac => unreachable!("{library_type}"),
                };

                let read_counts_entry = read_counts.get_mut(&feature_type).unwrap();
                checked_add(
                    &mut read_counts_entry.sequenced_reads[index],
                    metrics.total_reads.count(),
                );
                checked_add(
                    &mut read_counts_entry.barcode_corrected_sequenced_reads[index],
                    metrics.barcode_corrected_sequenced_reads.count(),
                );
                checked_add(
                    &mut read_counts_entry.unmapped_reads[index],
                    metrics.unmapped_reads.count(),
                );
            }
        }
        Ok(BarcodeSummaryData {
            genome_barcode_counts,
            read_counts,
        })
    }
}

/// Write the barcode summary HDF5 file.
fn write_barcode_summary_h5(
    path: &Path,
    per_barcode_metrics: &[BarcodeMetricsShardFile],
    feature_ref: &FeatureReference,
    barcode_index: &BarcodeIndex,
) -> Result<()> {
    let file = hdf5::File::create(path)?;

    // Avoid failure in writing a hdf5 dataset of length 0 with gzip chunk size 1
    if barcode_index.is_empty() {
        bail!(
            "No 10x barcodes were observed in the experiment. This is likely the consequence of a \
            sample mixup or very poor sequencing quality on the barcode bases. Further execution \
            is halted."
        );
    }
    write_barcodes_column(&file, "bc_sequence", barcode_index.sorted_barcodes())?;

    let BarcodeSummaryData {
        genome_barcode_counts,
        read_counts,
    } = BarcodeSummaryData::new(per_barcode_metrics, feature_ref, barcode_index)?;

    for (&(feat, ref genome), barcode_counts) in &genome_barcode_counts {
        // The additional underscore is a bug and retained for backward compatibility.
        // Example metric names are
        // _multi_transcriptome_conf_mapped_barcoded_reads
        // ANTIBODY__multi_transcriptome_conf_mapped_barcoded_reads
        let underscore_genome = format!("_{genome}");
        let library_type_genome = join_metric_name(feat, &underscore_genome);
        let read_counts_name =
            format!("{library_type_genome}_transcriptome_conf_mapped_barcoded_reads");
        file.new_dataset::<u32>()
            .deflate(1)
            .shape((barcode_counts.conf_mapped_read_counts.len(),))
            .create(read_counts_name.as_str())?
            .write(&barcode_counts.conf_mapped_read_counts)?;

        let umi_counts_name =
            format!("{library_type_genome}_transcriptome_conf_mapped_deduped_barcoded_reads");
        file.new_dataset::<u32>()
            .deflate(1)
            .shape((barcode_counts.umi_counts.len(),))
            .create(umi_counts_name.as_str())?
            .write(&barcode_counts.umi_counts)?;
    }

    let BarcodeCounts {
        conf_mapped_read_counts: multi_read_counts,
        umi_counts,
        half_mapped_counts,
        split_mapped_counts,
    } = genome_barcode_counts
        .iter()
        .filter(|((t, _), _)| *t == FeatureType::Gene)
        .fold(
            BarcodeCounts::with_size(barcode_index.len()),
            |mut x, (_, y)| {
                x += y;
                x
            },
        );

    for (name, data) in [
        (
            "_multi_transcriptome_conf_mapped_barcoded_reads",
            multi_read_counts,
        ),
        (
            "_multi_transcriptome_conf_mapped_deduped_barcoded_reads",
            umi_counts,
        ),
        (
            "_multi_transcriptome_split_mapped_barcoded_reads",
            split_mapped_counts,
        ),
        (
            "_multi_transcriptome_half_mapped_barcoded_reads",
            half_mapped_counts,
        ),
    ] {
        if !data.is_empty() {
            file.new_dataset::<u32>()
                .deflate(1)
                .shape((data.len(),))
                .create(name)?
                .write(&data)?;
        }
    }

    for (
        feature_type,
        BarcodeReads {
            sequenced_reads,
            barcode_corrected_sequenced_reads,
            unmapped_reads,
        },
    ) in read_counts
    {
        for (name, data) in [
            ("sequenced_reads", sequenced_reads),
            (
                "barcode_corrected_sequenced_reads",
                barcode_corrected_sequenced_reads,
            ),
            ("unmapped_reads", unmapped_reads),
        ] {
            file.new_dataset::<u32>()
                .deflate(1)
                .shape((data.len(),))
                .create(join_metric_name(feature_type, name).as_str())?
                .write(&data)?;
        }
    }

    Ok(())
}

#[make_mro(mem_gb = 7, volatile = strict)]
impl MartianMain for WriteBarcodeSummary {
    type StageInputs = WriteBarcodeSummaryStageInputs;
    type StageOutputs = WriteBarcodeSummaryStageOutputs;

    /// Run the Martian stage WRITE_BARCODE_SUMMARY.
    fn main(
        &self,
        WriteBarcodeSummaryStageInputs {
            per_barcode_metrics,
            feature_reference,
            barcode_index,
        }: Self::StageInputs,
        rover: MartianRover,
    ) -> Result<Self::StageOutputs> {
        let barcode_summary: H5File = rover.make_path("barcode_summary");
        let barcode_index = barcode_index.read()?;
        write_barcode_summary_h5(
            &barcode_summary,
            &per_barcode_metrics,
            &feature_reference.read()?,
            &barcode_index,
        )?;
        Ok(Self::StageOutputs { barcode_summary })
    }
}
