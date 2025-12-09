//! Probe set reference CSV
#![deny(missing_docs)]

use crate::probe_set::{Probe, ProbeRegion, ProbeType};
use crate::reference::feature_reference::get_genome_from_feature_id;
use crate::reference::reference_info::ReferenceInfo;
use crate::{FeatureID, FeatureName};
use anyhow::{Context, Result, bail, ensure};
use itertools::Itertools;
use metric::TxHashMap;
use std::path::Path;
use transcriptome::{Gene, Transcriptome};

/// Workaround martian_filetype not being document-able
pub mod filetypes {
    #![expect(missing_docs)]
    use martian_derive::martian_filetype;
    use serde::{Deserialize, Serialize};
    martian_filetype! { TargetSetFile, "csv" }
}

/// Target probe set CSV file
pub use filetypes::TargetSetFile;

/// The probe set reference CSV file format version.
pub const PROBE_SET_FILE_FORMAT: &str = "3.0";

/// The probe set reference CSV header.
pub const PROBE_SET_HEADER: [&str; 6] = [
    "gene_id",
    "probe_seq",
    "probe_id",
    "included",
    "region",
    "gene_name",
];

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
        let genomes = if let Some(transcriptome_reference_path) = transcriptome_reference_path {
            ReferenceInfo::from_reference_path(transcriptome_reference_path)?.genomes
        } else {
            ReferenceInfo::from_probe_set_csv(self)?.genomes
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
                let genome = get_genome_from_feature_id(&gene_id, &genomes);
                let probe_seq = record[1].to_string();
                let probe_id = record[2].to_string();
                let included = record.get(3).map_or(Ok(true), |x| {
                    x.to_lowercase().parse().with_context(|| {
                        format!(r#"The column "included" must be "true" or "false" but saw "{x}""#)
                    })
                })?;
                let region = record.get(4).and_then(|r| {
                    if r.is_empty() {
                        None
                    } else {
                        Some(ProbeRegion::new(r))
                    }
                });

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

                let num_hyphens = probe_seq.chars().filter(|&x| x == '-').count();
                let probe_type = if probe_seq.starts_with('-') || probe_seq.ends_with('-') {
                    ensure!(
                        num_hyphens == 1,
                        "An unpaired probe must have exactly one hyphen \
                         for probe {probe_id}: {probe_seq}"
                    );
                    ProbeType::UnpairedGapAlign
                } else {
                    match num_hyphens {
                        0 | 1 => ProbeType::RTL,
                        2 => ProbeType::PairedGapAlign,
                        3.. => {
                            bail!("Too many hyphens in sequence of probe {probe_id}: {probe_seq}")
                        }
                    }
                };

                Ok(Probe {
                    probe_id,
                    gene: Gene {
                        id: gene_id,
                        name: gene_name,
                    },
                    probe_seq,
                    included,
                    region,
                    probe_type,
                    ref_sequence_name,
                    ref_sequence_pos,
                    cigar_string,
                    genome: genome.clone(),
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::probe_set::merge_probe_set_csvs;
    use insta::assert_json_snapshot;
    use tempfile::NamedTempFile;

    const HUMAN_PANEL_A: &str = "test/probe_set_merge/GRCh38-fmt3-refv24_a.csv";
    const MOUSE_PANEL_A: &str = "test/probe_set_merge/mm10-fmt3-refv20_a.csv";
    const HUMAN_MOUSE_REF: &str = "test/reference/GRCh38-and-mm10_ref_tiny";
    const HUMAN_REF: &str = "test/reference/GRCh38_ref_tiny";
    #[test]
    fn test_read_one_panel_with_ref() -> Result<()> {
        let human_panel_a = TargetSetFile::from(HUMAN_PANEL_A);
        let human_ref = Path::new(HUMAN_REF);

        assert_json_snapshot!(human_panel_a.read(Some(human_ref))?);
        Ok(())
    }

    #[test]
    fn test_read_one_panel_no_ref() -> Result<()> {
        let human_panel_a = TargetSetFile::from(HUMAN_PANEL_A);

        assert_json_snapshot!(human_panel_a.read(None)?);
        Ok(())
    }

    #[test]
    fn test_read_two_panels_with_ref() -> Result<()> {
        let human_panel_a = TargetSetFile::from(HUMAN_PANEL_A);
        let mouse_panel_a = TargetSetFile::from(MOUSE_PANEL_A);
        let human_mouse_ref = Path::new(HUMAN_MOUSE_REF);

        let mut merged_panel_file = NamedTempFile::with_suffix(".csv")?;
        merge_probe_set_csvs(
            &[human_panel_a, mouse_panel_a],
            &mut merged_panel_file,
            Some(human_mouse_ref),
        )?;
        let merged_panel = TargetSetFile::from(merged_panel_file.path());

        assert_json_snapshot!(merged_panel.read(Some(human_mouse_ref))?);
        Ok(())
    }

    #[test]
    fn test_read_two_panels_no_ref() -> Result<()> {
        let human_panel_a = TargetSetFile::from(HUMAN_PANEL_A);
        let mouse_panel_a = TargetSetFile::from(MOUSE_PANEL_A);

        let mut merged_panel_file = NamedTempFile::with_suffix(".csv")?;
        merge_probe_set_csvs(
            &[human_panel_a, mouse_panel_a],
            &mut merged_panel_file,
            None,
        )?;
        let merged_panel = TargetSetFile::from(merged_panel_file.path());

        assert_json_snapshot!(merged_panel.read(None)?);
        Ok(())
    }

    #[test]
    fn test_read_two_panels_with_ref_at_read() -> Result<()> {
        let human_panel_a = TargetSetFile::from(HUMAN_PANEL_A);
        let mouse_panel_a = TargetSetFile::from(MOUSE_PANEL_A);
        let human_mouse_ref = Path::new(HUMAN_MOUSE_REF);

        let mut merged_panel_file = NamedTempFile::with_suffix(".csv")?;
        merge_probe_set_csvs(
            &[human_panel_a, mouse_panel_a],
            &mut merged_panel_file,
            None,
        )?;
        let merged_panel = TargetSetFile::from(merged_panel_file.path());

        assert_json_snapshot!(merged_panel.read(Some(human_mouse_ref))?);
        Ok(())
    }

    #[test]
    #[should_panic]
    fn test_read_one_panel_wrong_ref_at_read_error() {
        let human_panel_a = TargetSetFile::from(HUMAN_PANEL_A);
        let human_mouse_ref = Path::new(HUMAN_MOUSE_REF);

        let _ = human_panel_a.read(Some(human_mouse_ref));
    }
}
