use crate::mark_dups::DupInfo;
use crate::transcript::{
    AnnotationData, AnnotationParams, AnnotationRegion, PairAnnotationData, TranscriptAnnotator,
};
use anyhow::Result;
use barcode::HasBarcode;
use bio_types::strand::Strand;
use cr_bam::bam_tags::{
    ExtraFlags, ANTISENSE_TAG, EXTRA_FLAGS_TAG, FEATURE_IDS_TAG, FEATURE_QUAL_TAG, FEATURE_RAW_TAG,
    FEATURE_SEQ_TAG, GENE_ID_TAG, GENE_NAME_TAG, MULTIMAPPER_TAG, PROBE_TAG, PROC_BC_SEQ_TAG,
    PROC_UMI_SEQ_TAG, RAW_BARCODE_QUAL_TAG, RAW_BARCODE_SEQ_TAG, RAW_GEL_BEAD_BARCODE_QUAL_TAG,
    RAW_GEL_BEAD_BARCODE_SEQ_TAG, RAW_UMI_QUAL_TAG, RAW_UMI_SEQ_TAG, READ_GROUP_TAG, REGION_TAG,
    REST_R1_QUAL_TAG, REST_R1_SEQ_TAG, REST_R2_QUAL_TAG, REST_R2_SEQ_TAG, TRANSCRIPT_TAG,
    UNPAIRED_GENE_ID_TAG, UNPAIRED_GENE_NAME_TAG,
};
use cr_types::chemistry::ChemistryDef;
use cr_types::probe_set::MappedProbe;
use cr_types::reference::feature_extraction::FeatureData;
use cr_types::reference::feature_reference::FeatureReference;
use cr_types::reference::genome_of_chrom::GenomeName;
use cr_types::rna_read::{RnaChunk, RnaRead, UmiPart, HIGH_CONF_MAPQ};
use cr_types::types::{LibraryFeatures, ReqStrand};
use cr_types::utils::calculate_median_of_sorted;
use cr_types::UmiCount;
use itertools::Itertools;
use martian_derive::{martian_filetype, MartianStruct};
use martian_filetypes::bin_file::BinaryFormat;
use martian_filetypes::lz4_file::Lz4;
use rust_htslib::bam::record::{Aux, Record};
use serde::{Deserialize, Serialize};
use std::cmp::min;
use std::collections::HashSet;
use std::iter::{zip, FromIterator};
use std::path::Path;
use std::slice::Chunks;
use transcriptome::Gene;
use umi::UmiInfo;

pub const MAX_INSERT_SIZE: i64 = 1000;

fn attach_umi_tags(umi: &UmiInfo, record: &mut Record) {
    record
        .push_aux(RAW_UMI_SEQ_TAG, Aux::String(umi.seq.as_str()))
        .unwrap();
    record
        .push_aux(RAW_UMI_QUAL_TAG, Aux::String(umi.qual.as_str()))
        .unwrap();
}

/// TranscriptAnnotator handles annotating individual alignment records.  ReadAnnotator handle additional details
/// relating to paired-end reads and reads that generate multiple alignment records.
pub struct ReadAnnotator {
    annotator: TranscriptAnnotator,
}

impl ReadAnnotator {
    pub fn new(
        reference_path: &Path,
        chemistry_def: &ChemistryDef,
        include_exons: bool,
        include_introns: bool,
    ) -> Result<ReadAnnotator> {
        // NOTE: `Strand` can be Unknown as well, but ReqStrand is Forward or Reverse
        let chemistry_strandedness = match chemistry_def.strandedness {
            ReqStrand::Forward => Strand::Forward,
            ReqStrand::Reverse => Strand::Reverse,
        };
        let chemistry_fiveprime = chemistry_def.endedness == Some(fastq_set::WhichEnd::FivePrime);
        let params = AnnotationParams {
            chemistry_strandedness,
            chemistry_fiveprime,
            intergenic_trim_bases: 0,
            intronic_trim_bases: 0,
            junction_trim_bases: 0,
            region_min_overlap: 0.5,
            include_exons,
            include_introns,
        };
        let annotator = TranscriptAnnotator::new(reference_path, params)?;
        Ok(ReadAnnotator { annotator })
    }

    pub fn annotate_read_se(
        &self,
        read: RnaRead,
        alignments: Vec<Record>,
        umi_info: UmiInfo,
    ) -> ReadAnnotations {
        let mut data = alignments
            .into_iter()
            .map(|aln| RecordAnnotation::new_se(&self.annotator, aln))
            .collect::<Vec<_>>();
        rescue_alignments_se(&mut data);
        ReadAnnotations::from_records(read, data, umi_info, false)
    }

    pub fn annotate_read_pe(
        &self,
        read: RnaRead,
        r1_data: Vec<Record>,
        r2_data: Vec<Record>,
        umi_info: UmiInfo,
    ) -> ReadAnnotations {
        // There should be at least one paired-end alignment
        assert!(!r2_data.is_empty());
        let pair_improper = r1_data.len() != r2_data.len();
        // Annotate pairs
        let mut data: Vec<_> = zip(r1_data, r2_data)
            .map(|(r1, r2)| RecordAnnotation::new_pe(&self.annotator, r1, r2))
            .collect();
        rescue_alignments_pe(data.as_mut_slice());
        ReadAnnotations::from_records(read, data, umi_info, pair_improper)
    }
}

/// Use transcriptome alignments to promote a single genome alignment
/// when none are confidently mapped to the genome.
/// Returns true if rescue took place.
fn rescue_alignments_se(recs: &mut [RecordAnnotation]) -> bool {
    // Check if rescue is appropriate and determine which record to promote
    let mut seen_genes = HashSet::new();
    let mut promote_index: Option<usize> = None;

    for (i, r1) in recs.iter().enumerate() {
        // Abort if any of the records mapped uniquely to the genome
        if is_conf_mapped(r1.rec().0) {
            return false;
        }

        if let Some(anno) = r1.annotation().0 {
            // Only consider transcriptomic alignments for rescue
            if !anno.transcripts.iter().any(|x| x.tx_align.is_some()) {
                continue;
            }
            let genes = &anno.genes;
            // Track which record/record-pair we should promote;
            // Take the first record/pair with 1 gene
            if genes.len() == 1 {
                promote_index = promote_index.or(Some(i));
            }

            // Track number of distinct genes we're aligned to
            seen_genes.extend(genes);
        };
    }
    // There are >1 candidate genes to align to
    // or there are no candidate genes to align to.
    if seen_genes.len() > 1 || promote_index.is_none() {
        return false;
    }

    // Promote a single alignment
    for (i, r1) in recs.iter_mut().enumerate() {
        if promote_index.unwrap() == i {
            // Promote one alignment
            r1.set_rescued();
            let (rec, _) = r1.mut_rec();
            rec.set_mapq(HIGH_CONF_MAPQ);
            cr_bam::bam::set_primary(rec);
        } else {
            let (rec, _) = r1.mut_rec();
            // Demote the rest
            rec.set_mapq(0);
            rec.set_secondary();
        }
    }

    true
}

/// in PeMapped, either read could be mapped, the other unmapped, so
/// is_conf_mapped must test !is_unmapped
fn is_conf_mapped(rec: &Record) -> bool {
    !rec.is_unmapped() && rec.mapq() == HIGH_CONF_MAPQ
}

/// Use transcriptome alignments to promote a single genome alignment
/// when none are confidently mapped to the genome.
/// Returns true if rescue took place.
fn rescue_alignments_pe(pairs: &mut [RecordAnnotation]) -> bool {
    // Check if rescue is appropriate and determine which record to promote
    let mut seen_genes = HashSet::new();
    let mut promote_index: Option<usize> = None;

    for (i, pair) in pairs.iter().enumerate() {
        match pair {
            RecordAnnotation::PeMapped(rec1, _, rec2, _, anno) => {
                // Abort if any of the records mapped uniquely to the genome
                if is_conf_mapped(rec1) || is_conf_mapped(rec2) {
                    return false;
                }

                let genes = &anno.genes;

                // Track which record/record-pair we should promote;
                // Take the first record/pair with 1 gene
                if genes.len() == 1 {
                    promote_index = promote_index.or(Some(i));
                }

                // Track number of distinct genes we're aligned to
                seen_genes.extend(genes);
            }
            RecordAnnotation::Unmapped(_, _) => {}
            _ => unimplemented!(),
        }
    }
    // There are >1 candidate genes to align to
    // or there are no candidate genes to align to.
    if seen_genes.len() > 1 || promote_index.is_none() {
        return false;
    }

    // Promote a single alignment
    for (i, pair) in pairs.iter_mut().enumerate() {
        match pair {
            RecordAnnotation::PeMapped(
                ref mut rec1,
                ref mut anno1,
                ref mut rec2,
                ref mut anno2,
                _,
            ) => {
                if promote_index.unwrap() == i {
                    anno1.rescued = true;
                    rec1.set_mapq(HIGH_CONF_MAPQ);
                    cr_bam::bam::set_primary(rec1);
                    anno2.rescued = true;
                    rec2.set_mapq(HIGH_CONF_MAPQ);
                    cr_bam::bam::set_primary(rec2);
                } else {
                    rec1.set_mapq(0);
                    rec1.set_secondary();
                    rec2.set_mapq(0);
                    rec2.set_secondary();
                }
            }
            _ => {
                if promote_index.unwrap() == i {
                    panic!("attempted to promote invalid records!")
                }
            }
        }
    }

    true
}

pub trait AnnotationInfo {
    /// Return whether the read was discarded by subsampling.
    fn is_discarded(&self) -> bool;

    // Return whether the read is mapped.
    fn is_mapped(&self) -> bool;

    // Return whether the read is confidently mapped.
    fn is_conf_mapped(&self) -> bool;

    /// Return the mapping quality.
    fn mapq(&self) -> u8;

    // Return whether the read is mapped antisense.
    fn is_mapped_antisense(&self) -> bool;

    // Return whether the read is confidently mapped antisense.
    fn is_conf_mapped_antisense(&self) -> bool;

    /// Did we get a concordant mapping
    fn is_gene_discordant(&self) -> bool;

    /// Is a read pair "improper", AKA an imbalanced number of R1/R2 pairs from STAR
    /// this should not happen
    fn is_pair_improper(&self) -> bool;

    /// Is the alignment mapped uniquely to the
    /// genome and compatible with a single gene
    fn is_conf_mapped_to_transcriptome(&self) -> bool;

    // In intron mode on/off is the alignment high mapq, uniquely mapped to one
    // gene and has a TX tag associated with a transcript, i.e.,
    // TX:i:<transcipt id>,<strand><pos>,<cigar>
    fn is_conf_mapped_unique_txomic(&self) -> bool;

    /// Does the read correspond uniquely to a single gene or feature
    fn is_conf_mapped_to_feature(&self) -> bool;

    /// The unique genome this annotation confidently maps to if any
    fn conf_mapped_genome(&self) -> Option<&GenomeName>;

    /// If we confidently mapped to a unique gene, which is it.
    fn conf_mapped_gene(&self) -> Option<(&GenomeName, &Gene)>;

    /// If we confidently mapped to a unique feature, which is it.
    fn conf_mapped_feature(&self, _: &FeatureReference) -> Option<usize>;

    /// The unique region in a genome this annotation confidently maps to if any
    fn conf_mapped_region(&self) -> Option<(&GenomeName, AnnotationRegion)>;

    fn records(self) -> Vec<Record>;

    /// The set of genomes this annotation maps to. The set will
    /// be empty for an unmapped alignment and guaranteed to have at least
    /// one element for a mapped alignment.
    fn mapped_genomes(&self) -> HashSet<&GenomeName>;

    /// The set of (genome, gene) pairs this annotation maps to.
    fn mapped_genes(&self) -> HashSet<(&GenomeName, &Gene)>;

    /// The set of (genome, region) pairs this annotation maps to. The set will
    /// be empty for an unmapped alignment and guaranteed to have at least
    /// one element for a mapped alignment.
    fn mapped_regions(&self) -> HashSet<(&GenomeName, AnnotationRegion)>;

    /// If the read was a feature-barcode read
    fn is_feature_read(&self) -> bool;

    /// If the feature data was extracted
    fn is_feature_extracted(&self) -> bool;

    /// If the feature data was corrected
    fn is_feature_corrected(&self) -> bool;

    /// If extracted features weren't found
    fn is_feature_invalid(&self) -> bool;
}

martian_filetype! { ReadAnnotationsFile, "ann" }
pub type ReadAnnotationsFormat = Lz4<BinaryFormat<ReadAnnotationsFile, Vec<ReadAnnotations>>>;

#[derive(Clone, Debug, Serialize, Deserialize, MartianStruct)]
pub struct AnnotationFiles {
    pub num_reads: usize,
    pub files: Vec<ReadAnnotationsFormat>,
}

impl AnnotationFiles {
    pub fn make_chunks(&self, max_chunks: usize) -> Chunks<'_, ReadAnnotationsFormat> {
        let files_per_chunk = ((self.files.len() + max_chunks - 1) / max_chunks).max(1);
        self.files.chunks(files_per_chunk)
    }
}

/// All data associated with a qname (single-end read or read pair)
/// The desired order of SAM tags is
/// NH:i HI:i AS:i nM:i TX:Z GX:Z GN:Z fr:Z fq:Z fb:Z fx:Z
/// RE:A MM:i xf:i ts:i li:i CR:Z CY:Z CB:Z UR:Z UY:Z UB:Z RG:Z
#[derive(Serialize, Deserialize)]
pub struct ReadAnnotations {
    /// Original read data
    pub read: RnaRead,

    /// Matched TSO sequence.
    pub matched_tso: bool,

    /// when paired-end, if there were imbalanced number of R1/R2 alignments
    pub pair_improper: bool,

    /// Primary alignment
    pub primary: RecordAnnotation,

    pub umi_info: UmiInfo,

    // filled in by mark_dups
    pub dup_info: Option<DupInfo>,

    /// alignment for PD purposes, generally for aligning FB reads to the reference
    pub pd_alignment: Option<RecordAnnotation>,
}

impl ReadAnnotations {
    pub fn from_records(
        read: RnaRead,
        mut recs: Vec<RecordAnnotation>,
        umi_info: UmiInfo,
        pair_improper: bool,
    ) -> ReadAnnotations {
        // there must be a primary alignment or we panic
        let primary = recs.iter().position(|x| !x.rec().0.is_secondary()).unwrap();
        let primary = recs.remove(primary);

        ReadAnnotations {
            read,
            matched_tso: false,
            pair_improper,
            primary,
            umi_info,
            dup_info: None,
            pd_alignment: None,
        }
    }

    /// Construct a new ReadAnnotations from a MappedProbe and SAM Record.
    pub fn from_mapped_probe(
        read: RnaRead,
        recs: Vec<Record>,
        mut mapped_probe: MappedProbe,
        umi_info: UmiInfo,
    ) -> ReadAnnotations {
        // There ought to be precisely one primary (non-secondary) alignment record.
        let mut rec = recs.into_iter().find(|x| !x.is_secondary()).unwrap();
        let star_multimapped = !rec.is_unmapped() && rec.mapq() < HIGH_CONF_MAPQ;

        // Clear the UNMAPPED flag if the read maps to a probe.
        if mapped_probe.is_mapped() {
            rec.unset_unmapped();
        }

        // Set the mapping quality.
        rec.set_mapq(mapped_probe.mapq());

        // Set the multimapper tag MM:i:1 if the read maps confidently to a probe but multimapped with STAR.
        mapped_probe.set_rescued(star_multimapped && mapped_probe.is_conf_mapped());

        ReadAnnotations {
            read,
            matched_tso: false,
            pair_improper: false,
            primary: RecordAnnotation::new_probe(rec, mapped_probe),
            umi_info,
            dup_info: None,
            pd_alignment: None,
        }
    }

    fn for_each_rec<F: Fn(&mut RecordAnnotation)>(&mut self, f: F) {
        f(&mut self.primary);
    }

    /// Attach tags needed by bamtofastq for RTL libraries.
    fn attach_rtl_tags(&mut self) {
        assert!(!self.read.r2_exists());
        let rec = self.primary.mut_rec().0;

        // Attach the gel-bead barcode.
        rec.push_aux(
            RAW_GEL_BEAD_BARCODE_SEQ_TAG,
            Aux::String(self.read.raw_bc_construct_seq().gel_bead().as_str()),
        )
        .unwrap();
        rec.push_aux(
            RAW_GEL_BEAD_BARCODE_QUAL_TAG,
            Aux::String(self.read.raw_bc_construct_qual().gel_bead().as_str()),
        )
        .unwrap();

        // Attach the rest of R1, which follows the UMI.
        let [UmiPart::Untranslated { range: umi }] = self.read.umi_parts.as_slice() else {
                unreachable!()
            };
        assert_eq!(umi.read(), fastq_set::WhichRead::R1);
        let umi_end = umi.offset() + umi.len().unwrap();
        let rest_r1_seq =
            std::str::from_utf8(&self.read.raw_illumina_read1_seq()[umi_end..]).unwrap();
        if !rest_r1_seq.is_empty() {
            let rest_r1_qual =
                std::str::from_utf8(&self.read.raw_illumina_read1_qual()[umi_end..]).unwrap();
            rec.push_aux(REST_R1_SEQ_TAG, Aux::String(rest_r1_seq))
                .unwrap();
            rec.push_aux(REST_R1_QUAL_TAG, Aux::String(rest_r1_qual))
                .unwrap();
        }

        // Attach the rest of R2, which follows the probe.
        assert_eq!(self.read.r1_range().read(), fastq_set::WhichRead::R2);
        assert_eq!(self.read.r1_range().offset(), 0);
        let probe_end = self.read.r1_range().offset() + self.read.r1_range().len().unwrap();
        let rest_r2_seq =
            std::str::from_utf8(&self.read.raw_illumina_read2_seq()[probe_end..]).unwrap();
        if !rest_r2_seq.is_empty() {
            let rest_r2_qual =
                std::str::from_utf8(&self.read.raw_illumina_read2_qual()[probe_end..]).unwrap();
            rec.push_aux(REST_R2_SEQ_TAG, Aux::String(rest_r2_seq))
                .unwrap();
            rec.push_aux(REST_R2_QUAL_TAG, Aux::String(rest_r2_qual))
                .unwrap();
        }
    }

    /// Attach the barcode tags to the SAM records.
    fn attach_barcode_tags(&mut self) {
        let bc_seq = self.read.raw_bc_seq();
        self.primary.for_each_rec(|rec| {
            rec.push_aux(RAW_BARCODE_SEQ_TAG, Aux::String(bc_seq.as_str()))
                .unwrap()
        });

        let bc_qual = self.read.raw_bc_qual();
        self.primary.for_each_rec(|rec| {
            rec.push_aux(RAW_BARCODE_QUAL_TAG, Aux::String(bc_qual.as_str()))
                .unwrap()
        });

        if self.read.barcode().is_valid() {
            let bc_vec = self.read.barcode().to_string();
            self.primary
                .for_each_rec(|rec| rec.push_aux(PROC_BC_SEQ_TAG, Aux::String(&bc_vec)).unwrap());
        }

        if self.read.has_probe_barcode() {
            self.attach_rtl_tags();
        }
    }

    /// Attach the UMI tags to the SAM records.
    fn attach_umi_tags(&mut self) {
        // Attach the raw UMI tags to the SAM records.
        // have to do this the hard way b/c borrowck
        let (rec1, rec2) = self.primary.mut_rec();
        attach_umi_tags(&self.umi_info, rec1);
        if let Some(rec2) = rec2 {
            attach_umi_tags(&self.umi_info, rec2);
        }

        // Attach the corrected UMI tag to the SAM records.
        if self.umi_info.is_valid {
            let umi = if let Some(dup_info) = &self.dup_info {
                dup_info.processed_umi.as_str()
            } else {
                self.umi_info.seq.as_str()
            };
            self.primary
                .for_each_rec(|rec| rec.push_aux(PROC_UMI_SEQ_TAG, Aux::String(umi)).unwrap());
        }
    }

    pub fn attach_tags(&mut self, read_chunks: &[RnaChunk]) {
        if let Some(dup_info) = self.dup_info.as_ref() {
            // Set dup
            if !dup_info.is_umi_count()
                && !dup_info.is_low_support_umi
                && !dup_info.is_filtered_target_umi
            {
                self.primary
                    .for_each_rec(rust_htslib::bam::Record::set_duplicate);
            }
        }

        let (is_low_support_umi, is_filtered_target_umi) =
            if let Some(dup_info) = self.dup_info.as_ref() {
                (dup_info.is_low_support_umi, dup_info.is_filtered_target_umi)
            } else {
                (false, false)
            };

        let is_valid_umi = self.umi_info.is_valid;
        let is_valid_bc = self.read.barcode().is_valid();
        let read_group = &self.read.read_chunk(read_chunks).read_group;
        self.for_each_rec(|ann| {
            ann.attach_tags(
                is_valid_umi,
                is_low_support_umi,
                is_filtered_target_umi,
                is_valid_bc,
                read_group,
            );
        });
        self.attach_barcode_tags();
        self.attach_umi_tags();
    }

    pub fn insert_size(&self) -> Option<i64> {
        match &self.primary {
            RecordAnnotation::SeMapped(_, ref anno) => calculate_median_of_sorted(
                &anno
                    .transcripts
                    .iter()
                    .filter_map(|tx| {
                        tx.tx_align.as_ref().and_then(|aln| {
                            let insert_size = aln.se_insert_size as i64;
                            if insert_size <= MAX_INSERT_SIZE && tx.strand == ReqStrand::Forward {
                                Some(insert_size)
                            } else {
                                None
                            }
                        })
                    })
                    .sorted()
                    .collect_vec(),
            ),
            RecordAnnotation::PeMapped(_, ref anno1, _, ref anno2, annop) => {
                calculate_median_of_sorted(
                    &annop
                        .genes
                        .iter()
                        .filter_map(|gene| {
                            let tx1 = anno1.transcripts.iter().find(|x| &x.gene == gene);
                            let tx2 = anno2.transcripts.iter().find(|x| &x.gene == gene);
                            match (tx1, tx2) {
                                (Some(tx1), Some(tx2)) => {
                                    if let Some(aln1) = tx1.tx_align.as_ref() {
                                        if let Some(aln2) = tx2.tx_align.as_ref() {
                                            let start = aln1.pos.min(aln2.pos);
                                            let end =
                                                (aln1.pos + aln1.alen).max(aln2.pos + aln2.alen);
                                            let insert_size = end - start;
                                            if insert_size <= MAX_INSERT_SIZE {
                                                return Some(insert_size);
                                            }
                                        }
                                    }
                                    None
                                }
                                _ => None,
                            }
                        })
                        .sorted()
                        .collect_vec(),
                )
            }
            _ => None,
        }
    }
}

impl ReadAnnotations {
    pub fn iter_ann_info(&self) -> impl Iterator<Item = &RecordAnnotation> {
        std::iter::once(&self.primary)
    }
    pub fn umi_count(&self) -> Option<UmiCount> {
        self.dup_info.and_then(|info| info.umi_count())
    }
}

impl AnnotationInfo for ReadAnnotations {
    /// Return whether the read was discarded by subsampling.
    fn is_discarded(&self) -> bool {
        self.primary.is_discarded()
    }

    /// Return whether the read is mapped.
    fn is_mapped(&self) -> bool {
        self.primary.is_mapped()
    }

    /// Return whether the read is confidently mapped.
    fn is_conf_mapped(&self) -> bool {
        self.primary.is_conf_mapped()
    }

    /// Return the mapping quality.
    fn mapq(&self) -> u8 {
        self.primary.mapq()
    }

    fn is_mapped_antisense(&self) -> bool {
        self.iter_ann_info()
            .any(AnnotationInfo::is_mapped_antisense)
    }

    fn is_conf_mapped_antisense(&self) -> bool {
        self.primary.is_conf_mapped_antisense()
    }

    fn is_gene_discordant(&self) -> bool {
        self.primary.is_gene_discordant()
    }

    fn is_pair_improper(&self) -> bool {
        self.pair_improper
    }

    fn is_conf_mapped_to_transcriptome(&self) -> bool {
        self.primary.is_conf_mapped_to_transcriptome()
    }

    fn is_conf_mapped_unique_txomic(&self) -> bool {
        self.primary.is_conf_mapped_unique_txomic()
    }

    fn is_conf_mapped_to_feature(&self) -> bool {
        self.primary.is_conf_mapped_to_feature()
    }

    fn conf_mapped_genome(&self) -> Option<&GenomeName> {
        self.primary.conf_mapped_genome()
    }

    fn conf_mapped_gene(&self) -> Option<(&GenomeName, &Gene)> {
        self.primary.conf_mapped_gene()
    }

    fn conf_mapped_feature(&self, feature_reference: &FeatureReference) -> Option<usize> {
        self.primary.conf_mapped_feature(feature_reference)
    }

    fn conf_mapped_region(&self) -> Option<(&GenomeName, AnnotationRegion)> {
        self.primary.conf_mapped_region()
    }

    fn records(self) -> Vec<Record> {
        self.primary.records()
    }

    fn mapped_genomes(&self) -> HashSet<&GenomeName> {
        self.primary.mapped_genomes()
    }

    fn mapped_genes(&self) -> HashSet<(&GenomeName, &Gene)> {
        self.iter_ann_info()
            .flat_map(AnnotationInfo::mapped_genes)
            .collect()
        // let mut result = self.primary.mapped_genes();
        // result
    }

    fn mapped_regions(&self) -> HashSet<(&GenomeName, AnnotationRegion)> {
        self.primary.mapped_regions()
    }

    fn is_feature_read(&self) -> bool {
        matches!(
            self.read.library_feats(),
            LibraryFeatures::FeatureBarcodes(_)
        )
    }

    fn is_feature_extracted(&self) -> bool {
        self.primary.is_feature_extracted()
    }

    fn is_feature_corrected(&self) -> bool {
        self.primary.is_feature_corrected()
    }

    fn is_feature_invalid(&self) -> bool {
        self.primary.is_feature_invalid()
    }
}

/// All data associated with a single BAM record
#[derive(Serialize, Deserialize)]
pub enum RecordAnnotation {
    Discarded(Record),
    Unmapped(Record, Option<Record>),
    SeMapped(Record, AnnotationData),
    PeMapped(
        Record,
        AnnotationData,
        Record,
        AnnotationData,
        PairAnnotationData,
    ),
    Probe(Record, MappedProbe),
    FeatureExtracted(Record, FeatureData, Option<Record>),
}

impl AnnotationInfo for RecordAnnotation {
    fn is_discarded(&self) -> bool {
        matches!(self, RecordAnnotation::Discarded(_))
    }

    fn is_mapped(&self) -> bool {
        match self {
            RecordAnnotation::SeMapped(_, _) => true,
            RecordAnnotation::PeMapped(_, _, _, _, _) => true,
            RecordAnnotation::Probe(_, data) => data.is_mapped(),
            RecordAnnotation::FeatureExtracted(_, _, _) => false,
            RecordAnnotation::Unmapped(_, _) => false,
            RecordAnnotation::Discarded(_) => unreachable!(),
        }
    }

    fn is_conf_mapped(&self) -> bool {
        match self {
            RecordAnnotation::SeMapped(ref rec, _) => is_conf_mapped(rec),
            RecordAnnotation::PeMapped(ref rec1, _, ref rec2, _, _) => {
                is_conf_mapped(rec1) || is_conf_mapped(rec2)
            }
            RecordAnnotation::Probe(_, data) => data.is_conf_mapped(),
            RecordAnnotation::FeatureExtracted(_, _, _) => false,
            RecordAnnotation::Unmapped(_, _) => false,
            RecordAnnotation::Discarded(_) => unreachable!(),
        }
    }

    /// Return the mapping quality of a read mapped to the transcriptome.
    /// Return the minimum mapping quality of paired-end reads.
    /// Return 0 for feature barcoding reads.
    /// Return 0 for unmapped reads.
    fn mapq(&self) -> u8 {
        match self {
            RecordAnnotation::SeMapped(rec, _) => rec.mapq(),
            RecordAnnotation::PeMapped(rec1, _, rec2, _, _) => min(rec1.mapq(), rec2.mapq()),
            RecordAnnotation::Probe(_rec, data) => data.mapq(),
            RecordAnnotation::FeatureExtracted(_, _, _) => 0,
            RecordAnnotation::Unmapped(_, _) => 0,
            RecordAnnotation::Discarded(_) => unreachable!(),
        }
    }

    fn is_mapped_antisense(&self) -> bool {
        match self {
            RecordAnnotation::SeMapped(_, anno) => anno.is_antisense(),
            RecordAnnotation::PeMapped(_, anno1, _, anno2, _) => {
                anno1.is_antisense() && anno2.is_antisense()
            }
            RecordAnnotation::Probe(_, _) => false,
            RecordAnnotation::FeatureExtracted(_, _, _) => false,
            RecordAnnotation::Unmapped(_, _) => false,
            RecordAnnotation::Discarded(_) => unreachable!(),
        }
    }

    fn is_conf_mapped_antisense(&self) -> bool {
        match self {
            RecordAnnotation::SeMapped(rec, anno) => is_conf_mapped(rec) && anno.is_antisense(),
            RecordAnnotation::PeMapped(rec1, anno1, rec2, anno2, _) => {
                let rec1_is_conf_mapped = is_conf_mapped(rec1);
                let rec2_is_conf_mapped = is_conf_mapped(rec2);
                if rec1_is_conf_mapped && rec2_is_conf_mapped {
                    return anno1.is_antisense() && anno2.is_antisense();
                } else if rec1_is_conf_mapped {
                    return anno1.is_antisense();
                } else if rec2_is_conf_mapped {
                    return anno2.is_antisense();
                }
                false
            }
            RecordAnnotation::Probe(_, _) => false,
            RecordAnnotation::FeatureExtracted(_, _, _) => false,
            RecordAnnotation::Unmapped(_, _) => false,
            RecordAnnotation::Discarded(_) => unreachable!(),
        }
    }

    fn is_gene_discordant(&self) -> bool {
        match self {
            RecordAnnotation::SeMapped(_, _) => false,
            RecordAnnotation::PeMapped(rec1, _, rec2, _, annop) => {
                is_conf_mapped(rec1) && is_conf_mapped(rec2) && annop.genes.is_empty()
            }
            RecordAnnotation::Probe(_, _) => false,
            RecordAnnotation::FeatureExtracted(_, _, _) => false,
            RecordAnnotation::Unmapped(_, _) => false,
            RecordAnnotation::Discarded(_) => false,
        }
    }

    // not implemented here, carried by the ReadAnnotations
    fn is_pair_improper(&self) -> bool {
        unimplemented!()
    }

    fn is_conf_mapped_to_feature(&self) -> bool {
        match self {
            RecordAnnotation::FeatureExtracted(_, data, _) => data.ids.len() == 1,
            _ => self.is_conf_mapped_to_transcriptome(),
        }
    }

    /// Is the alignment mapped uniquely to the
    /// genome and compatible with a single gene
    fn is_conf_mapped_to_transcriptome(&self) -> bool {
        self.conf_mapped_gene().is_some()
    }

    // In intron mode on/off is the alignment high mapq, uniquely mapped to one
    // gene and has a TX tag associated with a transcript, i.e.,
    // TX:i:<transcipt id>,<strand><pos>,<cigar>
    fn is_conf_mapped_unique_txomic(&self) -> bool {
        match self {
            RecordAnnotation::SeMapped(rec, anno) => {
                if is_conf_mapped(rec) && anno.genes.len() == 1 {
                    return anno
                        .transcripts
                        .iter()
                        .any(super::transcript::TranscriptAlignment::is_txomic);
                }
                false
            }
            RecordAnnotation::PeMapped(rec1, anno1, rec2, anno2, _) => {
                let rec1_is_conf_mapped = is_conf_mapped(rec1);
                let rec2_is_conf_mapped = is_conf_mapped(rec2);
                if rec1_is_conf_mapped && rec2_is_conf_mapped {
                    return anno1
                        .transcripts
                        .iter()
                        .any(super::transcript::TranscriptAlignment::is_txomic)
                        || anno2
                            .transcripts
                            .iter()
                            .any(super::transcript::TranscriptAlignment::is_txomic);
                } else if rec1_is_conf_mapped {
                    return anno1
                        .transcripts
                        .iter()
                        .any(super::transcript::TranscriptAlignment::is_txomic);
                } else if rec2_is_conf_mapped {
                    return anno2
                        .transcripts
                        .iter()
                        .any(super::transcript::TranscriptAlignment::is_txomic);
                }
                false
            }
            RecordAnnotation::Probe(_, data) => data.is_conf_mapped(),
            RecordAnnotation::FeatureExtracted(_, _, _) => false,
            RecordAnnotation::Unmapped(_, _) => false,
            RecordAnnotation::Discarded(_) => unreachable!(),
        }
    }

    fn conf_mapped_genome(&self) -> Option<&GenomeName> {
        match self {
            RecordAnnotation::SeMapped(rec, anno) => {
                if is_conf_mapped(rec) {
                    Some(&anno.genome)
                } else {
                    None
                }
            }
            RecordAnnotation::PeMapped(rec1, anno1, rec2, anno2, _) => {
                let rec1_is_conf_mapped = is_conf_mapped(rec1);
                let rec2_is_conf_mapped = is_conf_mapped(rec2);
                if rec1_is_conf_mapped && rec2_is_conf_mapped {
                    if anno1.genome == anno2.genome {
                        return Some(&anno1.genome);
                    }
                } else if rec1_is_conf_mapped {
                    return Some(&anno1.genome);
                } else if rec2_is_conf_mapped {
                    return Some(&anno2.genome);
                }
                None
            }
            RecordAnnotation::Probe(_, data) => {
                if data.is_conf_mapped() {
                    data.genome()
                } else {
                    None
                }
            }
            RecordAnnotation::FeatureExtracted(_, _, _) => None,
            RecordAnnotation::Unmapped(_, _) => None,
            RecordAnnotation::Discarded(_) => None,
        }
    }

    fn conf_mapped_gene(&self) -> Option<(&GenomeName, &Gene)> {
        match self {
            RecordAnnotation::SeMapped(rec, anno) => {
                if is_conf_mapped(rec) && anno.genes.len() == 1 {
                    Some((&anno.genome, &anno.genes[0]))
                } else {
                    None
                }
            }
            RecordAnnotation::PeMapped(rec1, anno1, rec2, anno2, annop) => {
                let rec1_is_conf_mapped = is_conf_mapped(rec1);
                let rec2_is_conf_mapped = is_conf_mapped(rec2);
                if rec1_is_conf_mapped && rec2_is_conf_mapped {
                    if annop.genes.len() == 1 {
                        return Some((&anno1.genome, &annop.genes[0]));
                    }
                } else if rec1_is_conf_mapped {
                    if anno1.genes.len() == 1 {
                        return Some((&anno1.genome, &anno1.genes[0]));
                    }
                } else if rec2_is_conf_mapped && anno2.genes.len() == 1 {
                    return Some((&anno2.genome, &anno2.genes[0]));
                }
                None
            }
            RecordAnnotation::Probe(_, data) => {
                data.conf_gene().map(|x| (data.genome().unwrap(), x))
            }
            RecordAnnotation::FeatureExtracted(_, _, _) => None,
            RecordAnnotation::Unmapped(_, _) => None,
            RecordAnnotation::Discarded(_) => None,
        }
    }

    fn conf_mapped_feature(&self, feature_reference: &FeatureReference) -> Option<usize> {
        match self {
            RecordAnnotation::FeatureExtracted(_, data, _) => {
                if data.ids.len() == 1 {
                    Some(data.ids[0].0)
                } else {
                    None
                }
            }
            _ => self
                .conf_mapped_gene()
                .map(|(_, gene)| feature_reference.gene_index(gene)),
        }
    }

    fn conf_mapped_region(&self) -> Option<(&GenomeName, AnnotationRegion)> {
        match self {
            RecordAnnotation::SeMapped(rec, anno) => {
                if is_conf_mapped(rec) {
                    Some((&anno.genome, anno.region))
                } else {
                    None
                }
            }
            RecordAnnotation::PeMapped(rec1, anno1, rec2, anno2, _) => {
                let rec1_is_conf_mapped = is_conf_mapped(rec1);
                let rec2_is_conf_mapped = is_conf_mapped(rec2);
                if rec1_is_conf_mapped && rec2_is_conf_mapped {
                    if anno1.genome == anno2.genome {
                        // doing this to ensure paired-end reads spanning intron-exon get counted
                        // if one of the reads is exonic -> pair is exonic
                        if anno1.region == AnnotationRegion::Exonic
                            || anno2.region == AnnotationRegion::Exonic
                        {
                            return Some((&anno1.genome, AnnotationRegion::Exonic));
                        }
                        // default the region of anno1
                        return Some((&anno1.genome, anno1.region));
                        // this function returns None if read1 read2 genomes don't match
                        // possibility of None must be handled downstream!
                        // ideally is_conf_mapped would take an annotation and ensure genomes are same
                    }
                } else if rec1_is_conf_mapped {
                    return Some((&anno1.genome, anno1.region));
                } else if rec2_is_conf_mapped {
                    return Some((&anno2.genome, anno2.region));
                }
                None
            }
            RecordAnnotation::Probe(_, data) => {
                if data.is_conf_mapped() {
                    Some((data.genome().unwrap(), AnnotationRegion::Exonic))
                } else {
                    None
                }
            }
            RecordAnnotation::FeatureExtracted(_, _, _) => None,
            RecordAnnotation::Unmapped(_, _) => None,
            RecordAnnotation::Discarded(_) => None,
        }
    }

    fn records(self) -> Vec<Record> {
        match self {
            RecordAnnotation::Unmapped(rec1, rec2)
            | RecordAnnotation::FeatureExtracted(rec1, _, rec2) => {
                if let Some(rec2) = rec2 {
                    vec![rec1, rec2]
                } else {
                    vec![rec1]
                }
            }
            RecordAnnotation::SeMapped(rec, _) => vec![rec],
            RecordAnnotation::PeMapped(rec1, _, rec2, _, _) => vec![rec1, rec2],
            RecordAnnotation::Probe(rec, _) => vec![rec],
            RecordAnnotation::Discarded(rec) => vec![rec],
        }
    }

    fn mapped_genomes(&self) -> HashSet<&GenomeName> {
        if let RecordAnnotation::Probe(_, data) = self {
            return data.genome().into_iter().collect();
        }

        let (ann1, ann2) = self.annotation();
        HashSet::from_iter([ann1, ann2].iter().flatten().map(|ann| &ann.genome))
    }

    fn mapped_genes(&self) -> HashSet<(&GenomeName, &Gene)> {
        if let RecordAnnotation::Probe(_, data) = self {
            return data.genes().map(|x| (data.genome().unwrap(), x)).collect();
        }

        let (ann1, ann2) = self.annotation();
        HashSet::from_iter(
            [ann1, ann2]
                .iter()
                .flatten()
                .flat_map(|ann| ann.genes.iter().map(move |gene| (&ann.genome, gene))),
        )
    }

    fn mapped_regions(&self) -> HashSet<(&GenomeName, AnnotationRegion)> {
        if let RecordAnnotation::Probe(_, data) = self {
            return if data.is_mapped() {
                Some((data.genome().unwrap(), AnnotationRegion::Exonic))
            } else {
                None
            }
            .into_iter()
            .collect();
        }

        let (ann1, ann2) = self.annotation();
        HashSet::from_iter(
            [ann1, ann2]
                .iter()
                .flatten()
                .map(|ann| (&ann.genome, ann.region)),
        )
    }

    fn is_feature_read(&self) -> bool {
        unreachable!()
    }

    fn is_feature_extracted(&self) -> bool {
        match self {
            RecordAnnotation::FeatureExtracted(_, _, _) => true,
            RecordAnnotation::Unmapped(_, _) => false,
            _ => unreachable!(),
        }
    }

    fn is_feature_corrected(&self) -> bool {
        match self {
            RecordAnnotation::FeatureExtracted(_, ref data, _) => {
                if let Some(ref corrected) = data.corrected_barcode {
                    corrected != &data.barcode
                } else {
                    false
                }
            }
            _ => unreachable!(),
        }
    }

    fn is_feature_invalid(&self) -> bool {
        match self {
            RecordAnnotation::FeatureExtracted(_, ref data, _) => data.corrected_barcode.is_none(),
            _ => unreachable!(),
        }
    }
}

impl RecordAnnotation {
    pub fn new_se(annotator: &TranscriptAnnotator, rec: Record) -> Self {
        if rec.is_unmapped() {
            RecordAnnotation::Unmapped(rec, None)
        } else {
            let anno = annotator.annotate_alignment(&rec);
            RecordAnnotation::SeMapped(rec, anno)
        }
    }

    pub fn new_pe(annotator: &TranscriptAnnotator, rec1: Record, rec2: Record) -> Self {
        // STAR _shouldn't_ return pairs where only a single end is mapped,
        //   but if it does, consider the pair unmapped
        if rec1.is_unmapped() || rec2.is_unmapped() {
            RecordAnnotation::Unmapped(rec1, Some(rec2))
        } else {
            let anno1 = annotator.annotate_alignment(&rec1);
            let anno2 = annotator.annotate_alignment(&rec2);
            let annop = PairAnnotationData::from_pair(&anno1, &anno2);
            RecordAnnotation::PeMapped(rec1, anno1, rec2, anno2, annop)
        }
    }

    /// Construct a new RecordAnnotation from a SAM Record and MappedProbe.
    /// Return a new Probe if either Record or MappedProbe is mapped.
    /// Return an new Unmapped if both Record and MappedProbe are unmapped.
    pub fn new_probe(rec: Record, mapped_probe: MappedProbe) -> Self {
        if rec.is_unmapped() && !mapped_probe.is_mapped() {
            RecordAnnotation::Unmapped(rec, None)
        } else {
            RecordAnnotation::Probe(rec, mapped_probe)
        }
    }

    pub fn rec(&self) -> (&Record, Option<&Record>) {
        match self {
            RecordAnnotation::Unmapped(ref rec, ref rec2) => (rec, rec2.as_ref()),
            RecordAnnotation::SeMapped(ref rec, _) => (rec, None),
            RecordAnnotation::PeMapped(ref rec1, _, ref rec2, _, _) => (rec1, Some(rec2)),
            RecordAnnotation::FeatureExtracted(ref rec1, _, ref rec2) => (rec1, rec2.as_ref()),
            RecordAnnotation::Probe(ref rec, _) => (rec, None),
            RecordAnnotation::Discarded(rec) => (rec, None),
        }
    }

    pub fn for_each_rec<F: Fn(&mut Record)>(&mut self, f: F) {
        let (rec1, rec2) = self.mut_rec();
        f(rec1);
        rec2.map(f);
    }

    pub fn mut_rec(&mut self) -> (&mut Record, Option<&mut Record>) {
        match self {
            RecordAnnotation::Unmapped(ref mut rec1, ref mut rec2) => (rec1, rec2.as_mut()),
            RecordAnnotation::SeMapped(ref mut rec, _) => (rec, None),
            RecordAnnotation::PeMapped(ref mut rec1, _, ref mut rec2, _, _) => (rec1, Some(rec2)),
            RecordAnnotation::FeatureExtracted(ref mut rec1, _, ref mut rec2) => {
                (rec1, rec2.as_mut())
            }
            RecordAnnotation::Probe(ref mut rec, _) => (rec, None),
            RecordAnnotation::Discarded(rec) => (rec, None),
        }
    }

    pub fn annotation(&self) -> (Option<&AnnotationData>, Option<&AnnotationData>) {
        match self {
            RecordAnnotation::SeMapped(_, ref anno) => (Some(anno), None),
            RecordAnnotation::PeMapped(_, ref anno1, _, ref anno2, _) => (Some(anno1), Some(anno2)),
            RecordAnnotation::Probe(_, _) => unreachable!(),
            RecordAnnotation::FeatureExtracted(_, _, _) => (None, None),
            RecordAnnotation::Unmapped(_, _) => (None, None),
            RecordAnnotation::Discarded(_) => (None, None),
        }
    }

    pub fn set_rescued(&mut self) {
        match self {
            RecordAnnotation::SeMapped(_, ref mut anno) => anno.rescued = true,
            RecordAnnotation::PeMapped(_, ref mut anno1, _, ref mut anno2, _) => {
                anno1.rescued = true;
                anno2.rescued = true;
            }
            RecordAnnotation::Probe(_, _) => unreachable!(),
            RecordAnnotation::FeatureExtracted(_, _, _) => unreachable!(),
            RecordAnnotation::Unmapped(_, _) => panic!("Unmapped annotations cannot be rescued"),
            RecordAnnotation::Discarded(_) => unreachable!(),
        }
    }

    /// Add transcript TX tag to a BAM record.
    pub fn attach_transcript_tag(&mut self) {
        let attach_transcript_tag = |rec: &mut Record, anno: &AnnotationData| {
            if let Some(tag) = anno.make_tx_tag() {
                rec.push_aux(TRANSCRIPT_TAG, Aux::String(&tag)).unwrap();
            }
        };
        match self {
            RecordAnnotation::SeMapped(ref mut rec, ref anno) => attach_transcript_tag(rec, anno),
            RecordAnnotation::PeMapped(ref mut rec1, ref anno1, ref mut rec2, ref anno2, _) => {
                attach_transcript_tag(rec1, anno1);
                attach_transcript_tag(rec2, anno2);
            }
            RecordAnnotation::Probe(_, _) => {} // TODO
            RecordAnnotation::FeatureExtracted(_, _, _) => {}
            RecordAnnotation::Unmapped(_, _) => {}
            RecordAnnotation::Discarded(_) => {}
        }
    }

    /// Add tags to a BAM record.
    /// Set is_conf_mapped to true if the qname is confidently mapped to
    /// the transcriptome.
    pub fn attach_basic_tags(&mut self) {
        let attach_basic_tags = |record: &mut Record, anno: &AnnotationData| {
            if let Some(tag) = anno.make_re_tag() {
                record
                    .push_aux(REGION_TAG, Aux::Char(tag.as_bytes()[0]))
                    .unwrap();
            }

            if let Some(tag) = anno.make_mm_tag() {
                record.push_aux(MULTIMAPPER_TAG, Aux::I32(tag)).unwrap();
            }

            if let Some(tag) = anno.make_an_tag() {
                record.push_aux(ANTISENSE_TAG, Aux::String(&tag)).unwrap();
            }
        };
        match self {
            RecordAnnotation::SeMapped(ref mut rec, ref anno) => attach_basic_tags(rec, anno),
            RecordAnnotation::PeMapped(ref mut rec1, ref anno1, ref mut rec2, ref anno2, _) => {
                attach_basic_tags(rec1, anno1);
                attach_basic_tags(rec2, anno2);
            }
            RecordAnnotation::Probe(ref mut rec, ref data) => {
                if data.is_mapped() {
                    rec.push_aux(REGION_TAG, Aux::Char(b'E')).unwrap();
                }
                if data.is_rescued() {
                    rec.push_aux(MULTIMAPPER_TAG, Aux::I32(1)).unwrap();
                }
            }
            RecordAnnotation::FeatureExtracted(_, _, _) => {}
            RecordAnnotation::Unmapped(_, _) => {}
            RecordAnnotation::Discarded(_) => {}
        }
    }

    /// Add tags to a BAM record.
    /// Set is_conf_mapped to true if the qname is confidently mapped to
    /// the transcriptome.
    pub fn attach_tags(
        &mut self,
        is_valid_umi: bool,
        is_low_support_umi: bool,
        is_filtered_target_umi: bool,
        is_valid_bc: bool,
        read_group: &str,
    ) {
        self.for_each_rec(|rec| {
            rec.push_aux(READ_GROUP_TAG, Aux::String(read_group))
                .unwrap()
        });
        self.attach_transcript_tag();

        match self {
            RecordAnnotation::SeMapped(ref mut rec, ref anno) => {
                if let Some((tag_gx, tag_gn)) = anno.make_gx_gn_tags() {
                    rec.push_aux(GENE_ID_TAG, Aux::String(&tag_gx)).unwrap();
                    rec.push_aux(GENE_NAME_TAG, Aux::String(&tag_gn)).unwrap();
                    rec.push_aux(FEATURE_IDS_TAG, Aux::String(&tag_gx)).unwrap();
                }
            }
            RecordAnnotation::PeMapped(
                ref mut rec1,
                ref anno1,
                ref mut rec2,
                ref anno2,
                ref annop,
            ) => {
                if let Some((tag_gx, tag_gn)) = annop.make_gx_gn_tags() {
                    for (rec, anno) in &mut [(rec1, anno1), (rec2, anno2)] {
                        rec.push_aux(GENE_ID_TAG, Aux::String(&tag_gx)).unwrap();
                        rec.push_aux(GENE_NAME_TAG, Aux::String(&tag_gn)).unwrap();
                        rec.push_aux(FEATURE_IDS_TAG, Aux::String(&tag_gx)).unwrap();
                        if anno.genes != annop.genes {
                            if let Some((tag_gx, tag_gn)) = anno.make_gx_gn_tags() {
                                rec.push_aux(UNPAIRED_GENE_ID_TAG, Aux::String(&tag_gx))
                                    .unwrap();
                                rec.push_aux(UNPAIRED_GENE_NAME_TAG, Aux::String(&tag_gn))
                                    .unwrap();
                            }
                        }
                    }
                }
            }
            RecordAnnotation::Probe(ref mut rec, data) => {
                if data.is_mapped() {
                    let ids = data.genes().map(|x| &x.id).join(";");
                    let names = data.genes().map(|x| &x.name).join(";");
                    rec.push_aux(GENE_ID_TAG, Aux::String(&ids)).unwrap();
                    rec.push_aux(GENE_NAME_TAG, Aux::String(&names)).unwrap();
                    rec.push_aux(FEATURE_IDS_TAG, Aux::String(&ids)).unwrap();

                    let probes = data.probes().join(";");
                    rec.push_aux(PROBE_TAG, Aux::String(&probes)).unwrap();
                }
            }
            RecordAnnotation::FeatureExtracted(ref mut rec1, ref data, ref mut rec2) => {
                for rec in std::iter::once(rec1).chain(rec2.iter_mut()) {
                    rec.push_aux(
                        FEATURE_RAW_TAG,
                        Aux::String(std::str::from_utf8(&data.barcode).unwrap()),
                    )
                    .unwrap();
                    rec.push_aux(
                        FEATURE_QUAL_TAG,
                        Aux::String(std::str::from_utf8(&data.qual).unwrap()),
                    )
                    .unwrap();
                    if let Some(fb) = data.corrected_barcode.as_ref() {
                        rec.push_aux(
                            FEATURE_SEQ_TAG,
                            Aux::String(std::str::from_utf8(fb).unwrap()),
                        )
                        .unwrap();
                    }
                    if !data.ids.is_empty() {
                        let mut ids = data.ids.iter().map(|x| x.1.as_bytes()).collect::<Vec<_>>();
                        ids.sort();
                        let ids = ids.join(&b';');
                        rec.push_aux(
                            FEATURE_IDS_TAG,
                            Aux::String(std::str::from_utf8(&ids).unwrap()),
                        )
                        .unwrap();
                    }
                }
            }
            RecordAnnotation::Unmapped(_, _) => {}
            RecordAnnotation::Discarded(_) => {}
        }

        self.attach_basic_tags();
        self.attach_primary_tags(
            is_valid_umi,
            is_low_support_umi,
            is_filtered_target_umi,
            is_valid_bc,
        );
    }

    /// Add SAM tags to the primary alignment.
    fn attach_primary_tags(
        &mut self,
        is_valid_umi: bool,
        is_low_support_umi: bool,
        is_filtered_target_umi: bool,
        is_valid_bc: bool,
    ) {
        let is_conf_tx = self.is_conf_mapped_to_transcriptome();
        let is_conf_feat = self.is_conf_mapped_to_feature();
        let is_discordant = self.is_gene_discordant();
        self.for_each_rec(|rec| {
            attach_primary_tags(
                rec,
                is_conf_tx,
                is_conf_feat,
                is_discordant,
                is_valid_umi,
                is_low_support_umi,
                is_filtered_target_umi,
                is_valid_bc,
            )
        });
    }
}

#[allow(clippy::too_many_arguments)]
fn attach_primary_tags(
    rec: &mut Record,
    is_conf_mapped_to_transcriptome: bool,
    is_conf_mapped_to_feature: bool,
    is_gene_discordant: bool,
    is_valid_umi: bool,
    is_low_support_umi: bool,
    is_filtered_target_umi: bool,
    is_valid_bc: bool,
) {
    // Note: only attach these flags to primary alignment
    if !rec.is_secondary() {
        let mut flags: ExtraFlags = Default::default();

        if is_conf_mapped_to_transcriptome {
            flags |= ExtraFlags::CONF_MAPPED;
        }

        if is_conf_mapped_to_feature {
            flags |= ExtraFlags::CONF_FEATURE;
        }

        if is_low_support_umi {
            flags |= ExtraFlags::LOW_SUPPORT_UMI;
        }

        if is_gene_discordant {
            flags |= ExtraFlags::GENE_DISCORDANT;
        }

        if is_filtered_target_umi {
            flags |= ExtraFlags::FILTERED_TARGET_UMI;
        }

        if is_conf_mapped_to_feature
            && is_valid_bc
            && is_valid_umi
            && !rec.is_duplicate()
            && !is_low_support_umi
            && !is_filtered_target_umi
        {
            flags |= ExtraFlags::UMI_COUNT;
        }

        rec.push_aux(EXTRA_FLAGS_TAG, Aux::I32(flags.bits() as i32))
            .unwrap();
    }
}
