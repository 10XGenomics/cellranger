//!
//! A thin wrapper around orbit used for chemistry detection to distinguish between 3' and 5'
//!
#![deny(missing_docs)]

use super::chemistry_filter::{ChemistryFilter, DetectChemistryUnit};
use super::errors::DetectChemistryErrors;
use crate::detect_chemistry::chemistry_filter::DetectChemistryCandidates;
use anyhow::Result;
use barcode::whitelist::ReqStrand;
use cr_types::AlignerParam;
use cr_types::chemistry::ChemistryName;
use fastq_set::WhichEnd;
use fastq_set::read_pair::{ReadPair, ReadPart, WhichRead};
use metric::{TxHashSet, set};
use orbit::{StarAligner, StarReference, StarSettings};
use std::fmt;
use std::path::Path;
use tx_annotation::transcript::{AnnotationParams, TranscriptAnnotator};

#[derive(Debug, Default)]
pub(crate) struct MappingStats {
    pub(crate) total_reads: i64,
    pub(crate) sense_reads: i64,
    pub(crate) antisense_reads: i64,
    pub(crate) conf_mapped_reads: i64,
}

impl fmt::Display for MappingStats {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "{:20} = {}", "Total Reads", self.total_reads)?;
        writeln!(f, "{:20} = {}", "Mapped reads", self.conf_mapped_reads)?;
        writeln!(f, "{:20} = {}", "Sense reads", self.sense_reads)?;
        write!(f, "{:20} = {}", "Antisense reads", self.antisense_reads)
    }
}

const MIN_CONF_MAPPED_READS: usize = 1_000;
const MIN_CONF_MAPPED_READS_FRAC: f64 = 0.05;
const MIN_MARGIN: i64 = 2;

impl MappingStats {
    pub(crate) fn compatible_chemistries(&self) -> TxHashSet<ChemistryName> {
        use ChemistryName::{
            FivePrimeHT, FivePrimePE, FivePrimePEV3, FivePrimeR2, ThreePrimeV2, ThreePrimeV3CS1,
            ThreePrimeV3HTCS1, ThreePrimeV3HTPolyA, ThreePrimeV3LT, ThreePrimeV3PolyA,
        };
        if (self.conf_mapped_reads < MIN_CONF_MAPPED_READS as i64)
            || ((self.conf_mapped_reads as f64)
                < MIN_CONF_MAPPED_READS_FRAC * self.total_reads as f64)
        {
            set![]
        } else if self.sense_reads > MIN_MARGIN * self.antisense_reads {
            set![
                ThreePrimeV2,
                ThreePrimeV3PolyA,
                ThreePrimeV3CS1,
                ThreePrimeV3LT,
                ThreePrimeV3HTPolyA,
                ThreePrimeV3HTCS1
            ]
        } else if self.antisense_reads > MIN_MARGIN * self.sense_reads {
            set![FivePrimeR2, FivePrimePE, FivePrimePEV3, FivePrimeHT]
        } else {
            set![]
        }
    }

    pub(crate) fn help_text() -> String {
        format!(
            "In order to distinguish between the 3' vs 5' assay configuration the following conditions \
    		need to be satisfied:\n\
    		- A minimum of {} confidently mapped reads\n\
    		- A minimum of {:.1}% of the total reads considered needs to be confidently mapped\n\
    		- The number of sense reads need to be at least {}x compared to the antisense reads or \
    		vice versa",
            MIN_CONF_MAPPED_READS,
            100.0 * MIN_CONF_MAPPED_READS_FRAC,
            MIN_MARGIN
        )
    }
}

pub(crate) struct ReadMappingFilter {
    aligner: StarAligner,
    annotator: TranscriptAnnotator,
}

impl ChemistryFilter for ReadMappingFilter {
    const CONTEXT: &'static str = "Mapping based filtering";
    fn process_unit(
        &mut self,
        candidates: &DetectChemistryCandidates,
        unit: &DetectChemistryUnit,
        reads: &[ReadPair],
    ) -> Result<DetectChemistryCandidates, DetectChemistryErrors> {
        if !unit.library_type.is_gex() {
            return Ok(candidates.clone());
        }
        let stats = self.map_reads(reads);
        println!("{stats:#?}");
        let compatible_with_stats = stats.compatible_chemistries();
        let compatible_chems: DetectChemistryCandidates = candidates
            .iter()
            .filter_map(|(name, def)| {
                compatible_with_stats
                    .contains(name)
                    .then_some((*name, def.clone()))
            })
            .collect();
        if compatible_chems.is_empty() {
            return Err(DetectChemistryErrors::NotEnoughMapping {
                stats,
                unit: Box::new(unit.clone()),
                chems: candidates.clone(),
            });
        }
        Ok(compatible_chems)
    }
}

impl ReadMappingFilter {
    pub(crate) fn new(reference_path: &Path) -> Result<Self> {
        let star_path = reference_path.join("star");
        let settings = StarSettings::new(star_path.to_str().unwrap());
        let reference = StarReference::load(settings)?;
        let aligner = reference.get_aligner();
        let annotator = TranscriptAnnotator::new(
            reference_path,
            AnnotationParams {
                chemistry_strandedness: ReqStrand::Forward,
                chemistry_endedness: WhichEnd::ThreePrime,
                intergenic_trim_bases: 0,
                intronic_trim_bases: 0,
                junction_trim_bases: 0,
                region_min_overlap: 0.5,
                include_exons: true,
                include_introns: false,
                aligner: AlignerParam::Star,
            },
        )?;

        Ok(ReadMappingFilter { aligner, annotator })
    }

    pub(crate) fn map_reads(&mut self, read_pairs: &[ReadPair]) -> MappingStats {
        let mut mapping_stats = MappingStats::default();
        for read_pair in read_pairs {
            // do not crash if R2 doesn't exist (SC5P-R1)
            if read_pair.get(WhichRead::R2, ReadPart::Header).is_none() {
                continue;
            }
            // this will fail to handle SC5P-R1 data, but that has to be explicitly
            // set via the command-line
            let read_recs = self.aligner.align_read(
                read_pair.get(WhichRead::R2, ReadPart::Header).unwrap(),
                read_pair.get(WhichRead::R2, ReadPart::Seq).unwrap(),
                read_pair.get(WhichRead::R2, ReadPart::Qual).unwrap(),
            );
            mapping_stats.total_reads += 1;

            #[allow(clippy::absurd_extreme_comparisons)]
            // Allow the lint so that the logic remains the same if HIGH_CONF_MAPQ
            // is changed from 255
            if let Some(rec) = read_recs
                .into_iter()
                .find(|rec| rec.mapq() >= self.annotator.params.aligner.high_conf_mapq())
            // Should we .find(|rec| !rec.is_secondary() && !rec.is_unmapped())?
            {
                mapping_stats.conf_mapped_reads += 1;
                let ann = self.annotator.annotate_alignment(&rec);
                if ann.is_antisense() {
                    mapping_stats.antisense_reads += 1;
                }
                if ann.is_sense() {
                    mapping_stats.sense_reads += 1;
                }
            }
        }
        mapping_stats
    }
}
