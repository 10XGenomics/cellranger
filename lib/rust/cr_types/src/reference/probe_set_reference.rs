//! Probe set reference CSV

use crate::probe_set::{Probe, ProbeRegion, ProbeSetReferenceMetadata, ProbeType};
use crate::reference::reference_info::ReferenceInfo;
use crate::{FeatureID, FeatureName, GenomeName};
use anyhow::{bail, ensure, Context, Result};
use itertools::Itertools;
use martian_derive::martian_filetype;
use metric::TxHashMap;
use serde::{Deserialize, Serialize};
use std::path::Path;
use transcriptome::{Gene, Transcriptome};

martian_filetype! { TargetSetFile, "csv" }

impl TargetSetFile {
    /// Read a probe set reference CSV.
    pub fn read(&self, transcriptome_reference_path: Option<&Path>) -> Result<Vec<Probe>> {
        // Read the transcriptome GTF to map gene IDs to gene names.
        let gene_id_to_name: TxHashMap<_, _> =
            if let Some(reference_path) = transcriptome_reference_path {
                Transcriptome::from_reference_path(reference_path)?
                    .genes
                    .into_iter()
                    .map(|x| (x.id, x.name))
                    .collect()
            } else {
                TxHashMap::default()
            };

        // Read the probe set reference CSV file.
        let mut reader = csv::ReaderBuilder::new()
            .comment(Some(b'#'))
            .from_path(self)
            .with_context(|| self.display().to_string())?;

        // Ensure that the headers are correct.
        let header: Vec<_> = reader.headers().unwrap().iter().collect();
        assert_eq!(header[0], "gene_id");
        assert_eq!(header[1], "probe_seq");
        assert_eq!(header[2], "probe_id");
        if let Some(&included_header) = header.get(3) {
            assert_eq!(included_header, "included");
        }
        if let Some(&region_header) = header.get(4) {
            assert_eq!(region_header, "region");
        }
        if let Some(&gene_name_header) = header.get(5) {
            assert_eq!(gene_name_header, "gene_name");
        }
        if let Some(&ref_name_header) = header.get(6) {
            assert_eq!(ref_name_header, "ref_name");
        }
        if let Some(&ref_pos_header) = header.get(7) {
            assert_eq!(ref_pos_header, "ref_pos");
        }
        if let Some(&cigar_header) = header.get(8) {
            assert_eq!(cigar_header, "cigar");
        }

        reader
            .records()
            .map(|record| {
                let record = record?;
                let gene_id = record[0].to_string();
                let probe_seq_str = &record[1];
                let probe_seq = probe_seq_str.as_bytes();
                let probe_id = record[2].to_string();
                let included = record.get(3).map_or(Ok(true), |x| {
                    x.to_lowercase().parse().with_context(|| {
                        format!(r#"The column "included" must be "true" or "false" but saw "{x}""#)
                    })
                })?;
                let region = record
                    .get(4)
                    .map(str::to_string)
                    .map(|r| ProbeRegion::new(&r));

                // The transcriptome gene name, or the probe set gene name, or the gene ID.
                let gene_name = gene_id_to_name
                    .get(&gene_id)
                    .map(String::as_str)
                    .or_else(|| record.get(5))
                    .unwrap_or(&gene_id)
                    .to_string();

                let ref_sequence_name = record.get(6).unwrap_or("").to_string();
                let ref_sequence_pos: Option<usize> = record
                    .get(7)
                    .map(|pos| {
                        pos.parse().with_context(|| {
                            format!(r#"The column "ref_pos" must be an integer but saw "{pos}""#)
                        })
                    })
                    .transpose()?;
                let cigar_string = record.get(8).unwrap_or("").to_string();

                let num_hyphens = probe_seq.iter().filter(|&&x| x == b'-').count();
                let probe_type = if probe_seq.starts_with(b"-") || probe_seq.ends_with(b"-") {
                    ensure!(
                        num_hyphens == 1,
                        "An unpaired probe must have exactly one hyphen \
                         for probe {probe_id}: {probe_seq_str}"
                    );
                    ProbeType::UnpairedGapAlign
                } else {
                    match num_hyphens {
                        0 | 1 => ProbeType::RTL,
                        2 => ProbeType::PairedGapAlign,
                        3.. => {
                            bail!(
                                "Too many hyphens in sequence of probe {probe_id}: {probe_seq_str}"
                            )
                        }
                    }
                };

                Ok(Probe {
                    probe_id,
                    gene: Gene {
                        id: gene_id,
                        name: gene_name,
                    },
                    included,
                    region,
                    probe_type,
                    ref_sequence_name,
                    ref_sequence_pos,
                    cigar_string,
                })
            })
            .try_collect()
            .with_context(|| self.display().to_string())
    }

    /// Return the gene ID, gene name, and their `included` status.
    pub fn read_genes_and_included(
        &self,
        transcriptome_reference_path: Option<&Path>,
    ) -> Result<Vec<(FeatureID, FeatureName, bool)>> {
        let probes = self.read(transcriptome_reference_path)?;
        let gene_id_to_included: TxHashMap<_, _> = probes
            .iter()
            .map(|probe| (&probe.gene.id, probe.included))
            .into_group_map()
            .into_iter()
            .map(|(gene_id, included)| (gene_id.clone(), included.into_iter().all(|x| x)))
            .collect();
        Ok(probes
            .into_iter()
            .map(|x| (x.gene.id, x.gene.name))
            .unique()
            .map(|(gene_id, gene_name)| {
                let included = gene_id_to_included[&gene_id];
                (gene_id, gene_name, included)
            })
            .collect())
    }
}

/// Return the genome names from the reference transcriptome, or target set, or ["NONE"].
pub fn get_reference_genome_names(
    reference_info: Option<&ReferenceInfo>,
    probe_set_metadata: Option<&ProbeSetReferenceMetadata>,
) -> Vec<GenomeName> {
    if let Some(reference_info) = reference_info {
        reference_info.genomes.iter().sorted().cloned().collect()
    } else if let Some(probe_set_metadata) = probe_set_metadata {
        vec![GenomeName::from(probe_set_metadata.reference_genome())]
    } else {
        vec![GenomeName::from("NONE")]
    }
}
