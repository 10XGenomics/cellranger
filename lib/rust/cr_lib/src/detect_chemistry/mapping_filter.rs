//!
//! A thin wrapper around orbit used for chemistry detection to distinguish between 3' and 5'
//!

use super::chemistry_filter::{ChemistryFilter, DetectChemistryUnit};
use super::errors::DetectChemistryErrors;
use anyhow::Result;
use cr_types::chemistry::ChemistryName;
use cr_types::rna_read::HIGH_CONF_MAPQ;
use cr_types::ReqStrand;
use fastq_set::read_pair::{ReadPair, ReadPart, WhichRead};
use fastq_set::WhichEnd;
use metric::{set, TxHashSet};
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
            FivePrimeHT, FivePrimePE, FivePrimePEV3, FivePrimeR2, ThreePrimeV2, ThreePrimeV3,
            ThreePrimeV3HT, ThreePrimeV3LT,
        };
        if (self.conf_mapped_reads < MIN_CONF_MAPPED_READS as i64)
            || ((self.conf_mapped_reads as f64)
                < MIN_CONF_MAPPED_READS_FRAC * self.total_reads as f64)
        {
            set![]
        } else if self.sense_reads > MIN_MARGIN * self.antisense_reads {
            set![ThreePrimeV2, ThreePrimeV3, ThreePrimeV3LT, ThreePrimeV3HT]
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

pub(crate) struct ReadMappingFilter<'a> {
    aligner: StarAligner,
    annotator: TranscriptAnnotator,
    allowed_chems: &'a TxHashSet<ChemistryName>,
    chems: TxHashSet<ChemistryName>,
}

impl<'a> ChemistryFilter<'a> for ReadMappingFilter<'a> {
    fn context() -> &'static str {
        "Mapping based filtering"
    }
    fn allowed_chemistries(&self) -> &'a TxHashSet<ChemistryName> {
        self.allowed_chems
    }
    fn input_chemistries(&self) -> TxHashSet<ChemistryName> {
        self.chems.clone()
    }
    fn process_unit(
        &mut self,
        unit: &DetectChemistryUnit,
        reads: &[ReadPair],
    ) -> Result<TxHashSet<ChemistryName>, DetectChemistryErrors> {
        if !unit.library_type.is_gex() {
            return Ok(self.input_chemistries());
        }
        let stats = self.map_reads(reads);
        println!("{stats:#?}");
        let compatible_chems: TxHashSet<_> = stats
            .compatible_chemistries()
            .intersection(&self.input_chemistries())
            .copied()
            .collect();
        if compatible_chems.is_empty() {
            return Err(DetectChemistryErrors::NotEnoughMapping {
                stats,
                unit: Box::new(unit.clone()),
                chems: self.input_chemistries(),
            });
        }
        Ok(compatible_chems)
    }
}

impl<'a> ReadMappingFilter<'a> {
    pub(crate) fn new(
        reference_path: &Path,
        allowed_chems: &'a TxHashSet<ChemistryName>,
        chems: TxHashSet<ChemistryName>,
    ) -> Result<Self> {
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
            },
        )?;

        Ok(ReadMappingFilter {
            aligner,
            annotator,
            allowed_chems,
            chems,
        })
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
                .find(|rec| rec.mapq() >= HIGH_CONF_MAPQ)
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
