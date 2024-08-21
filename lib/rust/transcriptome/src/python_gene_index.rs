//! Datatypes for creating the GeneIndex class used in the
//! python codebase from the reference GTF file.

use crate::Transcriptome;
use anyhow::{bail, Context, Result};
use bio::io::fasta::IndexedReader;
use bio_types::strand::ReqStrand;
use itertools::Itertools;
use ordered_float::NotNan;
use serde::Serialize;
use serde_json;
use statrs::statistics::{Data, OrderStatistics};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, Read, Seek, Write};
use std::path::Path;
use std::process::Command;

#[derive(Serialize, Clone)]
struct Interval {
    chrom: String,
    start: u64,
    end: u64,
    length: u64,
    strand: Option<String>,
}

#[derive(Serialize)]
struct Transcript {
    gene: Gene,
    length: u64,
    gc_content: NotNan<f64>,
    intervals: Vec<Interval>,
}

#[derive(Serialize, Clone)]
struct Gene {
    id: String,
    name: String,
    length: f64,
    gc_content: f64,
    intervals: Vec<Interval>,
}

fn strand_string(s: ReqStrand) -> String {
    match s {
        ReqStrand::Forward => "+".to_string(),
        ReqStrand::Reverse => "-".to_string(),
    }
}

fn python_gene_index<R: Read + Seek>(
    txome: &Transcriptome,
    fasta_reader: &mut IndexedReader<R>,
) -> Result<(HashMap<String, Transcript>, Vec<Gene>)> {
    let chroms: HashSet<_> = fasta_reader
        .index
        .sequences()
        .into_iter()
        .map(|s| s.name)
        .collect();

    let mut transcripts = HashMap::new();
    let mut genes = Vec::new();
    for (gene_idx, txs) in &txome.gene_to_transcripts {
        let txs: Vec<_> = txs
            .iter()
            .map(|t| &txome.transcripts[t.0 as usize])
            // only include transcripts on an chromosome actually in the ref
            .filter(|t| chroms.contains(&t.chrom))
            .collect();

        let mut all_intervals = Vec::new();
        for tx in &txs {
            for exon in &tx.exons {
                let interval = Interval {
                    chrom: tx.chrom.clone(),
                    start: exon.start,
                    end: exon.end,
                    length: exon.end - exon.start,
                    strand: Some(strand_string(tx.strand)),
                };
                all_intervals.push(interval);
            }
        }

        let mut gene_intervals = Vec::new();

        all_intervals.sort_by(|i, j| i.chrom.cmp(&j.chrom));
        for (chrom, chrom_intervals) in &all_intervals.iter().group_by(|i| &i.chrom) {
            let chrom_intervals: Vec<_> = chrom_intervals.collect();
            let start = chrom_intervals.iter().map(|i| i.start).min().unwrap();
            let end = chrom_intervals.iter().map(|i| i.end).max().unwrap();

            let gene_interval = Interval {
                chrom: chrom.clone(),
                start,
                end,
                length: end - start,
                strand: None,
            };
            gene_intervals.push(gene_interval);
        }

        let gene = &txome.genes[gene_idx.0 as usize];

        // stats abou the gene are the median over the transcripts.
        let tx_lens: Vec<f64> = txs.iter().map(|t| t.len() as f64).collect();
        let length = Data::new(tx_lens.clone()).median();

        let mut tx_gcs: Vec<f64> = Vec::new();
        for tx in &txs {
            if chroms.contains(&tx.chrom) {
                let gc = tx.gc_content(fasta_reader)?;
                tx_gcs.push(gc);
            }
        }
        let gc_content = Data::new(tx_gcs.clone()).median();

        let new_gene = Gene {
            id: gene.id.clone(),
            name: gene.name.clone(),
            length,
            gc_content,
            intervals: gene_intervals,
        };
        genes.push(new_gene.clone());

        for (i, tx) in txs.iter().enumerate() {
            let mut intervals = Vec::new();
            for e in &tx.exons {
                let i = Interval {
                    chrom: tx.chrom.clone(),
                    start: e.start,
                    end: e.end,
                    length: e.end - e.start,
                    strand: Some(strand_string(tx.strand)),
                };
                intervals.push(i);
            }

            let new_tx = Transcript {
                gene: new_gene.clone(),
                length: tx_lens[i] as u64,
                gc_content: tx_gcs[i].try_into()?,
                intervals,
            };

            transcripts.insert(tx.id.clone(), new_tx);
        }
    }
    assert_eq!(txome.genes.len(), genes.len());

    Ok((transcripts, genes))
}

/// Open up an indexed fasta reader.  If the .fai isn't available, create it transiently on the fly
/// using samtools.
fn open_fasta_reader(path: &Path) -> Result<bio::io::fasta::IndexedReader<File>> {
    if let Ok(r) = bio::io::fasta::IndexedReader::from_file(&path) {
        return Ok(r);
    }

    // Create a symlink from the fasta to a tmp file in the current dir
    let tmpdir = tempfile::tempdir_in(".")?;
    let fa_symlink = tmpdir.path().join("genome.fa");
    std::os::unix::fs::symlink(path.canonicalize()?, &fa_symlink)?;

    // generate the .fa.fai index file
    let output = Command::new("samtools")
        .args(["faidx", fa_symlink.to_str().unwrap()])
        .output()
        .with_context(|| format!("samtools faidx {}", fa_symlink.display()))?;

    if output.status.success() {
        // load the index from the .fa.fai file
        let fai = fa_symlink.with_extension("fa.fai");
        let index = bio::io::fasta::Index::from_file(&fai).unwrap();
        // load the original fasta file, using the in-memory index
        Ok(bio::io::fasta::IndexedReader::with_index(
            File::open(path)?,
            index,
        ))
    } else {
        bail!(
            "Unable to create fasta index file. {}",
            std::str::from_utf8(&output.stderr).unwrap()
        )
    }
}

pub fn generate_gene_index(gtf_fn: &Path, fasta_fn: &Path) -> Result<String> {
    let txome = Transcriptome::from_gtf_path(gtf_fn)?;

    // throw an error if there are transcripts with no exons
    txome.check_for_transcripts_with_no_exons()?;

    let mut fasta = open_fasta_reader(fasta_fn).with_context(|| fasta_fn.display().to_string())?;

    let gene_index = python_gene_index(&txome, &mut fasta)?;

    Ok(serde_json::json!({
        "transcripts": gene_index.0,
        "genes": gene_index.1,
        "in_gtf_fn": gtf_fn.to_str().unwrap(),
        "in_fasta_fn": fasta_fn.to_str().unwrap(),
    })
    .to_string())
}

/// Load the transcriptome and gene property data from the given `gtf_fn` and `fasta_fn`.
/// Compute the GeneIndex data and convert to a json string.
pub fn write_gene_index(reference_path: &Path, json_fn: &Path) -> Result<()> {
    let gtf_path = reference_path.join("genes/genes.gtf");
    let fasta_path = reference_path.join("fasta/genome.fa");
    if !fasta_path.exists() {
        bail!("FASTA file doesn't exist: {}", fasta_path.display());
    }

    let r = generate_gene_index(&gtf_path, &fasta_path)?;
    let mut w = BufWriter::new(File::create(json_fn)?);
    w.write_all(r.as_bytes())?;
    Ok(())
}
