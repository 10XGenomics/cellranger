//! Martian stage WRITE_CONTIG_OUTS

use crate::assembly_types::{AsmReadsPerBcFormat, FastaFile, FastqFile};
use anyhow::Result;
use cr_types::MetricsFile;
use json_report_derive::JsonReport;
use martian::prelude::*;
use martian_derive::{make_mro, martian_filetype, MartianStruct};
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::{FileTypeRead, FileTypeWrite, LazyFileTypeIO};
use metric::{JsonReport, MeanMetric, PercentMetric, SimpleHistogram};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeSet, HashSet};
use std::io::Write;
use std::path::Path;
use std::process::Command;
use vdj_ann::annotate::ContigAnnotation;
use vdj_reference::{VdjChain, VdjChainCategory};
use vdj_types::VdjRegion;

martian_filetype!(FastaFaiFile, "fasta.fai");
martian_filetype!(BedFile, "bed");

#[derive(JsonReport)]
struct StageMetrics {
    /// The number of barcodes that are reported to be a cell by the VDJ pipeline.
    vdj_filtered_bcs: i64,
    /// This is the customer facing metric "Mean Read Pairs per Cell", defined as:
    /// "Number of input read pairs divided by the estimated number of cells."
    vdj_total_raw_read_pairs_per_filtered_bc: f64,
    /// This is the customer facing metric Fraction Reads in Cells, defined as:
    /// "Number of read pairs with cell-associated barcodes divided by the number of
    /// read pairs with valid barcodes."
    vdj_filtered_bcs_cum: PercentMetric,
    /// The fraction of productive contigs that are in barcodes that were called a cell.
    /// Numerator: Number of productive high confidence contigs in cell barcodes
    /// Denominator: Number of productive contigs in all barcodes
    vdj_high_conf_prod_contig: PercentMetric,
    /// The CS metric Mean Used Read Pairs per Cell.
    vdj_assemblable_read_pairs_per_filtered_bc: f64,
    // the number of contigs
    vdj_total_contig_count: i64,
    // the number of contig annotations
    vdj_total_annotated_regions: i64,
    // number of barcodes that contain a contig
    num_barcodes_with_contig: usize,
    // number of barcodes with a full length contig
    num_barcodes_with_full_length_contig: usize,
    // number of barcodes with a productive contig
    num_barcodes_with_productive_contig: usize,
    // number of barcodes with a cdr3 sequence
    num_barcodes_with_cdr3: usize,
    // average contig length (all contigs)
    mean_contig_length: f64,
    // fraction of barcodes with productive contig (num barcodes with productive contig/num barcodes with contig)
    frac_barcodes_with_productive_contig: f64,
    // fraction of barcodes with with annotated cdr3 (num barcodes with annotated cdr3/num barcodes with contig)
    frac_barcodes_with_cdr3: f64,
    // Number of barcodes with paired annotated cdr3s
    num_barcodes_with_paired_cdr3: usize,
    // Number of barcodes with paired productive contigs
    num_barcodes_with_paired_productive_contig: usize,
    // fraction of barcodes with a contig annotated as TRA
    frac_barcodes_with_tra: f64,
    // fraction of barcodes with a contig annotated as TRB
    frac_barcodes_with_trb: f64,
}

#[derive(Clone, Deserialize, MartianStruct)]
pub struct WriteContigOutsStageInputs {
    pub contig_annotations: JsonFile<Vec<ContigAnnotation>>,
    pub total_read_pairs: i64,
    pub corrected_bc_counts: JsonFile<SimpleHistogram<String>>,
    pub assemblable_reads_per_bc: AsmReadsPerBcFormat,
}

#[derive(Serialize, Deserialize, MartianStruct)]
pub struct WriteContigOutsStageOutputs {
    pub contig_fastq: FastqFile,
    pub filtered_contig_fastq: FastqFile,
    pub contig_fasta: FastaFile,
    pub contig_fasta_fai: Option<FastaFaiFile>,
    pub filtered_contig_fasta: FastaFile,
    pub annotations_bed: BedFile,
    pub cell_barcodes: JsonFile<BTreeSet<String>>,
    pub paired_cell_barcodes: JsonFile<BTreeSet<String>>,
    pub paired_prod_barcodes: JsonFile<BTreeSet<String>>,
    pub paired_cdr3_barcodes: JsonFile<BTreeSet<String>>,
    pub prod_barcodes: JsonFile<BTreeSet<String>>,
    pub cdr3_barcodes: JsonFile<BTreeSet<String>>,
    pub all_contig_barcodes: JsonFile<BTreeSet<String>>,
    pub summary: MetricsFile,
}

fn index_fasta(fasta_file: &Path) {
    Command::new("samtools")
        .arg("faidx")
        .arg(fasta_file)
        .status()
        .expect("samtools faidx failed");
}

pub fn extract_all_chain_types_from_contig(can: &ContigAnnotation) -> Vec<VdjChain> {
    let chain_v = can.get_region(VdjRegion::V).map(|c| c.feature.chain);
    let chain_d = can.get_region(VdjRegion::D).map(|c| c.feature.chain);
    let chain_j = can.get_region(VdjRegion::J).map(|c| c.feature.chain);
    let chain_c = can.get_region(VdjRegion::C).map(|c| c.feature.chain);
    [chain_v, chain_d, chain_j, chain_c]
        .into_iter()
        .flatten()
        .collect()
}

pub fn extract_chain_type_from_contig(can: &ContigAnnotation) -> Option<VdjChain> {
    let all_chain_types: Vec<VdjChain> = extract_all_chain_types_from_contig(can);

    if !all_chain_types.is_empty() {
        let head = all_chain_types[0];
        if all_chain_types.iter().all(|&c| c == head) {
            Some(head)
        } else {
            None
        }
    } else {
        None
    }
}

pub fn extract_chain_category_from_contig(can: &ContigAnnotation) -> Option<VdjChainCategory> {
    let all_chain_cats: Vec<VdjChainCategory> = extract_all_chain_types_from_contig(can)
        .into_iter()
        .map(VdjChainCategory::from)
        .collect();

    if !all_chain_cats.is_empty() {
        let head = all_chain_cats[0];
        if all_chain_cats.iter().all(|&c| c == head) {
            Some(head)
        } else {
            None
        }
    } else {
        None
    }
}

pub struct WriteContigOuts;

#[make_mro]
impl MartianMain for WriteContigOuts {
    type StageInputs = WriteContigOutsStageInputs;
    type StageOutputs = WriteContigOutsStageOutputs;

    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        let contig_fastq_file: FastqFile = rover.make_path("contig_fastq");
        let filtered_contig_fastq_file: FastqFile = rover.make_path("filtered_contig_fastq");
        let contig_fasta_file: FastaFile = rover.make_path("contig_fasta");
        let filtered_contig_fasta_file: FastaFile = rover.make_path("filtered_contig_fasta");
        let annotations_bed_file: BedFile = rover.make_path("annotations_bed");
        let cell_barcodes_file: JsonFile<_> = rover.make_path("cell_barcodes.json");
        let paired_cell_barcodes_file: JsonFile<_> = rover.make_path("paired_cell_barcodes");
        let paired_prod_barcodes_file: JsonFile<_> = rover.make_path("paired_prod_barcodes");
        let paired_cdr3_barcodes_file: JsonFile<_> = rover.make_path("paired_cdr3_barcodes");
        let prod_barcodes_file: JsonFile<_> = rover.make_path("prod_barcodes");
        let cdr3_barcodes_file: JsonFile<_> = rover.make_path("cdr3_barcodes");
        let all_contig_barcodes_file: JsonFile<_> = rover.make_path("all_contig_barcodes");
        let summary_file: MetricsFile = rover.make_path("summary");

        let mut have_fasta = false;

        // We want to drop all the writers so that buffer is flushed out and we can read them
        // later (for indexing etc.). Hence adding a block here.
        {
            let mut fastq = contig_fastq_file.buf_writer()?;
            let mut filt_fastq = filtered_contig_fastq_file.buf_writer()?;
            let mut fasta = contig_fasta_file.buf_writer()?;
            let mut filt_fasta = filtered_contig_fasta_file.buf_writer()?;
            let mut bed_writer = annotations_bed_file.buf_writer()?;
            let mut cell_barcodes = BTreeSet::new();
            let mut all_barcodes = BTreeSet::new();

            // Set of cells with productive light and heavy chains to compute pairing
            let mut cells_prod_light_chain = HashSet::<String>::new();
            let mut cells_prod_heavy_chain = HashSet::<String>::new();
            // Set of barcodes with productive light and heavy chains
            let mut barcodes_prod_light_chain = HashSet::<String>::new();
            let mut barcodes_prod_heavy_chain = HashSet::<String>::new();
            // Set of barcodes with annotated cdrs light and heavy chains
            let mut barcodes_cdr3_light_chain = HashSet::<String>::new();
            let mut barcodes_cdr3_heavy_chain = HashSet::<String>::new();
            // Set of barcodes with TRA/TRB
            let mut barcodes_tra = HashSet::<String>::new();
            let mut barcodes_trb = HashSet::<String>::new();

            let mut num_contigs = 0;
            let mut num_contig_annotations = 0;
            let mut num_productive_contigs = 0;
            let mut num_cell_productive_contigs = 0;

            let mut barcodes_with_contig = HashSet::new();
            let mut barcodes_with_full_length_contig = HashSet::new();
            let mut barcodes_with_productive_contig = HashSet::new();
            let mut barcodes_with_cdr3 = HashSet::new();
            let mut contig_lengths = MeanMetric::default();

            let reader = args.contig_annotations.lazy_reader()?;
            for can in reader {
                let can: ContigAnnotation = can?;
                num_contigs += 1;
                all_barcodes.insert(can.barcode.clone());
                let fastq_string =
                    format!("@{}\n{}\n+\n{}", can.contig_name, can.sequence, can.quals);

                let fasta_string = format!(">{}\n{}", can.contig_name, can.sequence);

                writeln!(&mut fastq, "{fastq_string}")?;
                writeln!(&mut fasta, "{fasta_string}")?;

                // Detecting chain type
                let chain_type = extract_chain_type_from_contig(&can);
                match chain_type {
                    Some(VdjChain::TRA) => {
                        barcodes_tra.insert(can.barcode.clone());
                    }
                    Some(VdjChain::TRB) => {
                        barcodes_trb.insert(can.barcode.clone());
                    }
                    _ => {}
                }
                // Using a different method for detecting chain_type, since non-productives may not have all of the VDJC regions
                let chain_category = extract_chain_category_from_contig(&can);
                if can.is_cell && can.high_confidence && can.is_productive() {
                    writeln!(&mut filt_fastq, "{fastq_string}")?;
                    writeln!(&mut filt_fasta, "{fasta_string}")?;
                    if let Some(chain_category) = chain_category {
                        num_cell_productive_contigs += 1;
                        match chain_category {
                            VdjChainCategory::Light => {
                                cells_prod_light_chain.insert(can.barcode.clone());
                            }
                            VdjChainCategory::Heavy => {
                                cells_prod_heavy_chain.insert(can.barcode.clone());
                            }
                        }
                    }
                }

                if can.is_productive() {
                    if let Some(chain_category) = chain_category {
                        num_productive_contigs += 1;
                        match chain_category {
                            VdjChainCategory::Light => {
                                barcodes_prod_light_chain.insert(can.barcode.clone());
                            }
                            VdjChainCategory::Heavy => {
                                barcodes_prod_heavy_chain.insert(can.barcode.clone());
                            }
                        }
                    }
                }

                if can.is_cell {
                    cell_barcodes.insert(can.barcode.clone());
                }

                barcodes_with_contig.insert(can.barcode.clone());
                if can.is_full_length() {
                    barcodes_with_full_length_contig.insert(can.barcode.clone());
                };
                if can.is_productive() {
                    barcodes_with_productive_contig.insert(can.barcode.clone());
                }
                if can.cdr3.is_some() {
                    barcodes_with_cdr3.insert(can.barcode.clone());

                    if let Some(c) = chain_category {
                        match c {
                            VdjChainCategory::Light => {
                                barcodes_cdr3_light_chain.insert(can.barcode.clone());
                            }
                            VdjChainCategory::Heavy => {
                                barcodes_cdr3_heavy_chain.insert(can.barcode.clone());
                            }
                        }
                    }
                };
                contig_lengths.record(can.sequence.len() as f64);

                for annot in &can.annotations {
                    num_contig_annotations += 1;
                    writeln!(
                        bed_writer,
                        "{}\t{}\t{}\t{}_{}",
                        can.contig_name,
                        annot.contig_match_start,
                        annot.contig_match_end,
                        annot.feature.gene_name,
                        annot.feature.region_type
                    )?;
                }

                have_fasta = true;
            }
            cell_barcodes_file.write(&cell_barcodes)?;
            let paired_cell_barcodes: BTreeSet<_> = cell_barcodes
                .iter()
                .filter(|&bc| {
                    cells_prod_light_chain.contains(bc) && cells_prod_heavy_chain.contains(bc)
                })
                .cloned()
                .collect();
            paired_cell_barcodes_file.write(&paired_cell_barcodes)?;

            let paired_prod_barcodes: BTreeSet<_> = all_barcodes
                .iter()
                .filter(|&bc| {
                    barcodes_prod_light_chain.contains(bc) && barcodes_prod_heavy_chain.contains(bc)
                })
                .cloned()
                .collect();
            paired_prod_barcodes_file.write(&paired_prod_barcodes)?;

            let paired_cdr3_barcodes: BTreeSet<_> = all_barcodes
                .iter()
                .filter(|&bc| {
                    barcodes_cdr3_light_chain.contains(bc) && barcodes_cdr3_heavy_chain.contains(bc)
                })
                .cloned()
                .collect();
            paired_cdr3_barcodes_file.write(&paired_cdr3_barcodes)?;

            let prod_barcodes: BTreeSet<_> = all_barcodes
                .iter()
                .filter(|&bc| {
                    barcodes_prod_light_chain.contains(bc) || barcodes_prod_heavy_chain.contains(bc)
                })
                .cloned()
                .collect();
            prod_barcodes_file.write(&prod_barcodes)?;

            let cdr3_barcodes: BTreeSet<_> = all_barcodes
                .iter()
                .filter(|&bc| {
                    barcodes_cdr3_light_chain.contains(bc) || barcodes_cdr3_heavy_chain.contains(bc)
                })
                .cloned()
                .collect();
            cdr3_barcodes_file.write(&cdr3_barcodes)?;
            all_contig_barcodes_file.write(&all_barcodes)?;
            // -------------------------------------------------------------------------------------
            // Compute the metrics here
            let metrics = StageMetrics {
                vdj_filtered_bcs: cell_barcodes.len() as i64,
                vdj_total_raw_read_pairs_per_filtered_bc: if !cell_barcodes.is_empty() {
                    (args.total_read_pairs as f64) / (cell_barcodes.len() as f64)
                } else {
                    0.0f64
                },
                vdj_filtered_bcs_cum: {
                    let bc_read_counts = args.corrected_bc_counts.read()?;
                    let reads_in_valid_barcodes: i64 = bc_read_counts.raw_counts().sum();
                    let reads_in_cell_barcodes =
                        cell_barcodes.iter().map(|bc| bc_read_counts.get(bc)).sum();
                    (reads_in_cell_barcodes, reads_in_valid_barcodes).into()
                },
                vdj_high_conf_prod_contig: (num_cell_productive_contigs, num_productive_contigs)
                    .into(),
                vdj_assemblable_read_pairs_per_filtered_bc: if !cell_barcodes.is_empty() {
                    let bc_read_counts = args.assemblable_reads_per_bc.read()?;
                    let assemblable_reads_in_cells: i64 =
                        cell_barcodes.iter().map(|bc| bc_read_counts.get(bc)).sum();
                    (assemblable_reads_in_cells as f64) / (cell_barcodes.len() as f64)
                } else {
                    0.0f64
                },
                vdj_total_contig_count: num_contigs,
                vdj_total_annotated_regions: num_contig_annotations,
                num_barcodes_with_contig: barcodes_with_contig.len(),
                num_barcodes_with_full_length_contig: barcodes_with_full_length_contig.len(),
                num_barcodes_with_productive_contig: barcodes_with_productive_contig.len(),
                num_barcodes_with_cdr3: barcodes_with_cdr3.len(),
                mean_contig_length: contig_lengths.mean(),
                frac_barcodes_with_productive_contig: (barcodes_with_productive_contig.len()
                    as f64)
                    / (barcodes_with_contig.len() as f64),
                frac_barcodes_with_cdr3: (barcodes_with_cdr3.len() as f64)
                    / (barcodes_with_contig.len() as f64),
                num_barcodes_with_paired_cdr3: paired_cdr3_barcodes.len(),
                num_barcodes_with_paired_productive_contig: paired_prod_barcodes.len(),
                frac_barcodes_with_tra: (barcodes_tra.len() as f64)
                    / (barcodes_with_contig.len() as f64),
                frac_barcodes_with_trb: (barcodes_trb.len() as f64)
                    / (barcodes_with_contig.len() as f64),
            };

            summary_file.write(&metrics.to_json_reporter())?;
        }

        // Index the fasta file.
        let contig_fasta_fai_file = if have_fasta {
            index_fasta(&contig_fasta_file);
            Some(rover.make_path("contig_fasta"))
        } else {
            None
        };

        Ok(WriteContigOutsStageOutputs {
            contig_fastq: contig_fastq_file,
            filtered_contig_fastq: filtered_contig_fastq_file,
            contig_fasta: contig_fasta_file,
            contig_fasta_fai: contig_fasta_fai_file,
            filtered_contig_fasta: filtered_contig_fasta_file,
            annotations_bed: annotations_bed_file,
            cell_barcodes: cell_barcodes_file,
            paired_cell_barcodes: paired_cell_barcodes_file,
            paired_prod_barcodes: paired_prod_barcodes_file,
            paired_cdr3_barcodes: paired_cdr3_barcodes_file,
            prod_barcodes: prod_barcodes_file,
            cdr3_barcodes: cdr3_barcodes_file,
            all_contig_barcodes: all_contig_barcodes_file,
            summary: summary_file,
        })
    }
}
