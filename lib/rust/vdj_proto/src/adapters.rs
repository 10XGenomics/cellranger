//!
//! - Converting from crate::types to vdj_types
//! - Additional functions to auto derived structs

use crate::types::vdj_proto_message::MessageContent;
use crate::types::{
    AsmData, GexData, MetricsSummary, VdjMetadata, VdjProtoMessage, VdjReferenceRaw,
};
use anyhow::Result;
use sha2::{Digest, Sha256};
use std::convert::Into;
use std::fs;
use std::io::Write;
use std::path::Path;

impl From<vdj_types::VdjChain> for crate::types::VdjChain {
    fn from(src: vdj_types::VdjChain) -> Self {
        match src {
            vdj_types::VdjChain::IGH => crate::types::VdjChain::Igh,
            vdj_types::VdjChain::IGK => crate::types::VdjChain::Igk,
            vdj_types::VdjChain::IGL => crate::types::VdjChain::Igl,
            vdj_types::VdjChain::TRA => crate::types::VdjChain::Tra,
            vdj_types::VdjChain::TRB => crate::types::VdjChain::Trb,
            vdj_types::VdjChain::TRG => crate::types::VdjChain::Trg,
            vdj_types::VdjChain::TRD => crate::types::VdjChain::Trd,
        }
    }
}

impl From<crate::types::VdjChain> for vdj_types::VdjChain {
    fn from(src: crate::types::VdjChain) -> Self {
        match src {
            crate::types::VdjChain::Igh => vdj_types::VdjChain::IGH,
            crate::types::VdjChain::Igk => vdj_types::VdjChain::IGK,
            crate::types::VdjChain::Igl => vdj_types::VdjChain::IGL,
            crate::types::VdjChain::Tra => vdj_types::VdjChain::TRA,
            crate::types::VdjChain::Trb => vdj_types::VdjChain::TRB,
            crate::types::VdjChain::Trg => vdj_types::VdjChain::TRG,
            crate::types::VdjChain::Trd => vdj_types::VdjChain::TRD,
        }
    }
}

impl From<vdj_types::VdjRegion> for crate::types::VdjRegion {
    fn from(src: vdj_types::VdjRegion) -> Self {
        match src {
            vdj_types::VdjRegion::UTR => crate::types::VdjRegion::U,
            vdj_types::VdjRegion::V => crate::types::VdjRegion::V,
            vdj_types::VdjRegion::D => crate::types::VdjRegion::D,
            vdj_types::VdjRegion::J => crate::types::VdjRegion::J,
            vdj_types::VdjRegion::C => crate::types::VdjRegion::C,
        }
    }
}

impl From<crate::types::VdjRegion> for vdj_types::VdjRegion {
    fn from(src: crate::types::VdjRegion) -> Self {
        match src {
            crate::types::VdjRegion::U => vdj_types::VdjRegion::UTR,
            crate::types::VdjRegion::V => vdj_types::VdjRegion::V,
            crate::types::VdjRegion::D => vdj_types::VdjRegion::D,
            crate::types::VdjRegion::J => vdj_types::VdjRegion::J,
            crate::types::VdjRegion::C => vdj_types::VdjRegion::C,
        }
    }
}

impl From<crate::types::AnnotationFeature> for vdj_ann::annotate::AnnotationFeature {
    fn from(src: crate::types::AnnotationFeature) -> Self {
        vdj_ann::annotate::AnnotationFeature {
            chain: crate::types::VdjChain::try_from(src.chain).unwrap().into(),
            display_name: src.display_name,
            feature_id: src.feature_id as usize,
            gene_name: src.gene_name,
            region_type: crate::types::VdjRegion::try_from(src.region_type)
                .unwrap()
                .into(),
        }
    }
}

impl From<vdj_ann::annotate::AnnotationFeature> for crate::types::AnnotationFeature {
    fn from(src: vdj_ann::annotate::AnnotationFeature) -> Self {
        crate::types::AnnotationFeature {
            chain: crate::types::VdjChain::from(src.chain).into(),
            display_name: src.display_name,
            feature_id: src.feature_id as u32,
            gene_name: src.gene_name,
            region_type: crate::types::VdjRegion::from(src.region_type).into(),
        }
    }
}

impl From<vdj_ann::annotate::AnnotationUnit> for crate::types::AnnotationUnit {
    fn from(src: vdj_ann::annotate::AnnotationUnit) -> Self {
        crate::types::AnnotationUnit {
            contig_match_start: src.contig_match_start as u32,
            contig_match_end: src.contig_match_end as u32,
            annotation_match_start: src.annotation_match_start as u32,
            annotation_match_end: src.annotation_match_end as u32,
            annotation_length: src.annotation_length as u32,
            cigar: src.cigar,
            score: src.score,
            feature: Some(src.feature.into()),
        }
    }
}

impl From<crate::types::AnnotationUnit> for vdj_ann::annotate::AnnotationUnit {
    fn from(src: crate::types::AnnotationUnit) -> Self {
        vdj_ann::annotate::AnnotationUnit {
            contig_match_start: src.contig_match_start as usize,
            contig_match_end: src.contig_match_end as usize,
            annotation_match_start: src.annotation_match_start as usize,
            annotation_match_end: src.annotation_match_end as usize,
            annotation_length: src.annotation_length as usize,
            cigar: src.cigar,
            score: src.score,
            feature: src.feature.unwrap().into(),
        }
    }
}

impl From<crate::types::JunctionSupport> for vdj_ann::annotate::JunctionSupport {
    fn from(src: crate::types::JunctionSupport) -> Self {
        vdj_ann::annotate::JunctionSupport {
            reads: src.reads,
            umis: src.umis,
        }
    }
}

impl From<vdj_ann::annotate::JunctionSupport> for crate::types::JunctionSupport {
    fn from(src: vdj_ann::annotate::JunctionSupport) -> Self {
        crate::types::JunctionSupport {
            reads: src.umis,
            umis: src.umis,
        }
    }
}

impl From<crate::types::Region> for vdj_ann::annotate::Region {
    fn from(src: crate::types::Region) -> Self {
        vdj_ann::annotate::Region {
            start: src.start as usize,
            stop: src.stop as usize,
            nt_seq: src.nt_seq,
            aa_seq: src.aa_seq,
        }
    }
}

impl From<vdj_ann::annotate::Region> for crate::types::Region {
    fn from(src: vdj_ann::annotate::Region) -> Self {
        crate::types::Region {
            start: src.start as i32,
            stop: src.stop as i32,
            nt_seq: src.nt_seq,
            aa_seq: src.aa_seq,
        }
    }
}

impl From<vdj_asm_utils::barcode_data::BarcodeDataBrief> for crate::types::BarcodeData {
    fn from(src: vdj_asm_utils::barcode_data::BarcodeDataBrief) -> Self {
        crate::types::BarcodeData {
            barcode: src.barcode,
            read_count: src.read_pairs,
            umi_count: src.total_ucounts as u64,
            xucounts: src.xucounts,
            ncontigs: src.ncontigs as u64,
            frac_reads_used: src.frac_reads_used,
        }
    }
}

impl From<crate::types::BarcodeData> for vdj_asm_utils::barcode_data::BarcodeDataBrief {
    fn from(src: crate::types::BarcodeData) -> Self {
        vdj_asm_utils::barcode_data::BarcodeDataBrief {
            barcode: src.barcode,
            read_pairs: src.read_count,
            total_ucounts: src.umi_count as usize,
            xucounts: src.xucounts,
            ncontigs: src.ncontigs as usize,
            frac_reads_used: src.frac_reads_used,
        }
    }
}

impl From<crate::types::ContigStatus> for vdj_ann::transcript::ContigStatus {
    fn from(src: crate::types::ContigStatus) -> Self {
        vdj_ann::transcript::ContigStatus {
            full_length: src.full_length,
            has_v_start: src.has_v_start,
            in_frame: src.in_frame,
            no_premature_stop: src.no_premature_stop,
            has_cdr3: src.has_cdr3,
            has_expected_size: src.has_expected_size,
            correct_ann_order: src.correct_ann_order,
        }
    }
}

impl From<vdj_ann::transcript::ContigStatus> for crate::types::ContigStatus {
    fn from(src: vdj_ann::transcript::ContigStatus) -> Self {
        crate::types::ContigStatus {
            full_length: src.full_length,
            has_v_start: src.has_v_start,
            in_frame: src.in_frame,
            no_premature_stop: src.no_premature_stop,
            has_cdr3: src.has_cdr3,
            has_expected_size: src.has_expected_size,
            correct_ann_order: src.correct_ann_order,
        }
    }
}

fn option_usize_to_i32(src: Option<usize>) -> i32 {
    src.map_or(-1i32, |x| x as i32)
}

fn i32_to_option_usize(src: i32) -> Option<usize> {
    if src < 0 {
        None
    } else {
        Some(src as usize)
    }
}

fn option_string_to_string(src: Option<String>) -> String {
    src.unwrap_or_default()
}

fn string_to_option_string(src: String) -> Option<String> {
    if src.is_empty() {
        None
    } else {
        Some(src)
    }
}

impl From<vdj_ann::annotate::ContigAnnotation> for crate::types::ContigAnnotation {
    fn from(src: vdj_ann::annotate::ContigAnnotation) -> Self {
        let full_length = src.is_full_length();
        crate::types::ContigAnnotation {
            barcode: src.barcode,
            contig_name: src.contig_name,
            sequence: src.sequence,
            quals: src.quals,
            read_count: src.read_count as u64,
            umi_count: src.umi_count as u64,
            start_codon_pos: option_usize_to_i32(src.start_codon_pos),
            stop_codon_pos: option_usize_to_i32(src.stop_codon_pos),
            aa_sequence: option_string_to_string(src.aa_sequence),
            cdr3: option_string_to_string(src.cdr3),
            cdr3_nt: option_string_to_string(src.cdr3_seq),
            cdr3_start: option_usize_to_i32(src.cdr3_start),
            cdr3_stop: option_usize_to_i32(src.cdr3_stop),
            annotations: src.annotations.into_iter().map(Into::into).collect(),
            clonotype: option_string_to_string(src.clonotype),
            raw_clonotype_id: option_string_to_string(src.info.raw_clonotype_id),
            raw_consensus_id: option_string_to_string(src.info.raw_consensus_id),
            exact_subclonotype_id: option_string_to_string(src.info.exact_subclonotype_id),
            high_confidence: src.high_confidence,
            is_cell: src.is_cell,
            productive: src.productive.unwrap_or(false),
            productive_criteria: src.productive_criteria.map(Into::into),
            filtered: src.filtered,
            frame: option_usize_to_i32(src.frame),
            asm_data: src.is_asm_cell.map(|c| AsmData { is_asm_cell: c }),
            gex_data: src.is_gex_cell.map(|c| GexData { is_gex_cell: c }),
            full_length,
            junction_support: src.junction_support.map(Into::into),

            fwr1: src.fwr1.map(Into::into),
            cdr1: src.cdr1.map(Into::into),
            fwr2: src.fwr2.map(Into::into),
            cdr2: src.cdr2.map(Into::into),
            fwr3: src.fwr3.map(Into::into),
            fwr4: src.fwr4.map(Into::into),
            sample: option_string_to_string(src.sample),
        }
    }
}

impl From<crate::types::ContigAnnotation> for vdj_ann::annotate::ContigAnnotation {
    fn from(src: crate::types::ContigAnnotation) -> Self {
        vdj_ann::annotate::ContigAnnotation {
            barcode: src.barcode,
            contig_name: src.contig_name,
            sequence: src.sequence,
            quals: src.quals,
            fraction_of_reads_for_this_barcode_provided_as_input_to_assembly: None,
            read_count: src.read_count as usize,
            umi_count: src.umi_count as usize,
            start_codon_pos: i32_to_option_usize(src.start_codon_pos),
            stop_codon_pos: i32_to_option_usize(src.stop_codon_pos),
            aa_sequence: string_to_option_string(src.aa_sequence),
            cdr3: string_to_option_string(src.cdr3),
            cdr3_seq: string_to_option_string(src.cdr3_nt),
            cdr3_start: i32_to_option_usize(src.cdr3_start),
            cdr3_stop: i32_to_option_usize(src.cdr3_stop),
            annotations: src.annotations.into_iter().map(Into::into).collect(),
            clonotype: string_to_option_string(src.clonotype),
            info: vdj_ann::annotate::ClonotypeInfo {
                raw_clonotype_id: string_to_option_string(src.raw_clonotype_id),
                raw_consensus_id: string_to_option_string(src.raw_consensus_id),
                exact_subclonotype_id: string_to_option_string(src.exact_subclonotype_id),
            },
            high_confidence: src.high_confidence,
            validated_umis: /* src.validated_umis */ None, // ?????????????????????????????????????
            non_validated_umis: /* src.validated_umis */ None, // ?????????????????????????????????
            invalidated_umis: /* src.invalidated_umis */ None, // ?????????????????????????????????
            is_cell: src.is_cell,
            is_asm_cell: src.asm_data.map(|a| a.is_asm_cell),
            is_gex_cell: src.gex_data.map(|g| g.is_gex_cell),
            productive: Some(src.productive),
            productive_criteria: src.productive_criteria.map(Into::into),
            filtered: src.filtered,
            frame: i32_to_option_usize(src.frame),
            full_length: Some(src.full_length),

            fwr1: src.fwr1.map(Into::into),
            cdr1: src.cdr1.map(Into::into),
            fwr2: src.fwr2.map(Into::into),
            cdr2: src.cdr2.map(Into::into),
            fwr3: src.fwr3.map(Into::into),
            fwr4: src.fwr4.map(Into::into),
            junction_support: src.junction_support.map(Into::into),
            sample: string_to_option_string(src.sample),
        }
    }
}

impl From<VdjMetadata> for VdjProtoMessage {
    fn from(metadata: VdjMetadata) -> Self {
        VdjProtoMessage {
            message_content: Some(MessageContent::Metadata(metadata)),
        }
    }
}

impl From<VdjReferenceRaw> for VdjProtoMessage {
    fn from(reference: VdjReferenceRaw) -> Self {
        VdjProtoMessage {
            message_content: Some(MessageContent::Reference(reference)),
        }
    }
}

impl From<MetricsSummary> for VdjProtoMessage {
    fn from(metrics: MetricsSummary) -> Self {
        VdjProtoMessage {
            message_content: Some(MessageContent::Metrics(metrics)),
        }
    }
}

impl From<vdj_ann::annotate::ContigAnnotation> for VdjProtoMessage {
    fn from(ann: vdj_ann::annotate::ContigAnnotation) -> Self {
        VdjProtoMessage {
            message_content: Some(MessageContent::Annotation(ann.into())),
        }
    }
}

impl From<vdj_asm_utils::barcode_data::BarcodeDataBrief> for VdjProtoMessage {
    fn from(brief: vdj_asm_utils::barcode_data::BarcodeDataBrief) -> Self {
        VdjProtoMessage {
            message_content: Some(MessageContent::BarcodeData(brief.into())),
        }
    }
}

impl VdjProtoMessage {
    pub fn content(self) -> MessageContent {
        self.message_content.unwrap()
    }
}

impl VdjReferenceRaw {
    pub fn new(reference_folder: &Path) -> Result<Self> {
        Ok(VdjReferenceRaw {
            regions: std::fs::read_to_string(reference_folder.join("fasta/regions.fa"))?,
            ref_json: std::fs::read_to_string(reference_folder.join("reference.json"))?,
        })
    }

    pub fn write_to_folder(&self, folder: &Path) -> Result<()> {
        fs::create_dir_all(folder.join("fasta"))?;
        fs::File::create(folder.join("fasta/regions.fa"))?.write_all(self.regions.as_bytes())?;
        fs::File::create(folder.join("reference.json"))?.write_all(self.ref_json.as_bytes())?;
        Ok(())
    }

    pub fn fasta_hash(&self) -> String {
        let mut hasher = Sha256::new();
        hasher.update(self.regions.as_bytes());
        let result = hasher.finalize();
        format!("{result:x}")
    }
}

impl From<crate::types::Receptor> for vdj_reference::VdjReceptor {
    fn from(src: crate::types::Receptor) -> vdj_reference::VdjReceptor {
        match src {
            crate::types::Receptor::Tcr => vdj_reference::VdjReceptor::TR,
            crate::types::Receptor::Tcrgd => vdj_reference::VdjReceptor::TRGD,
            crate::types::Receptor::Ig => vdj_reference::VdjReceptor::IG,
        }
    }
}

impl From<vdj_reference::VdjReceptor> for crate::types::Receptor {
    fn from(src: vdj_reference::VdjReceptor) -> crate::types::Receptor {
        match src {
            vdj_reference::VdjReceptor::TR => crate::types::Receptor::Tcr,
            vdj_reference::VdjReceptor::TRGD => crate::types::Receptor::Tcrgd,
            vdj_reference::VdjReceptor::IG => crate::types::Receptor::Ig,
        }
    }
}

impl VdjMetadata {
    pub fn vdj_receptor(&self) -> vdj_reference::VdjReceptor {
        crate::types::Receptor::try_from(self.receptor)
            .unwrap()
            .into()
    }
}

impl MetricsSummary {
    pub fn from_metrics_json(path: &Path) -> Result<Self> {
        let contents = std::fs::read_to_string(path)?;
        // Ensure that this is a valid metric json. Panic if not
        let v: serde_json::Value = serde_json::from_str(&contents).unwrap();
        assert!(
            v.is_object(),
            "Expecting metrics to be an Object(Map<String, Value>)"
        );
        Ok(MetricsSummary { raw_json: contents })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use strum::IntoEnumIterator;

    #[test]
    fn test_vdj_chain_roundtrip() {
        for chain in vdj_types::VdjChain::all() {
            let crate_chain: crate::types::VdjChain = chain.into();
            let back_again: vdj_types::VdjChain = crate_chain.into();
            assert_eq!(chain, back_again);
        }
    }

    #[test]
    fn test_vdj_region_roundtrip() {
        for region in vdj_types::VdjRegion::all() {
            let crate_region: crate::types::VdjRegion = region.into();
            let back_again: vdj_types::VdjRegion = crate_region.into();
            assert_eq!(region, back_again);
        }
    }

    #[test]
    fn test_vdj_receptor_roundtrip() {
        for receptor in vdj_reference::VdjReceptor::iter() {
            let crate_receptor: crate::types::Receptor = receptor.into();
            let back_again: vdj_reference::VdjReceptor = crate_receptor.into();
            assert_eq!(receptor, back_again);
        }
    }

    #[test]
    fn test_annotation_roundtrip() {
        let anns: Vec<vdj_ann::annotate::ContigAnnotation> =
            serde_json::from_reader(File::open("../vdj_asm_asm/test.json").unwrap()).unwrap();
        for ann in anns {
            let crate_anns: crate::types::ContigAnnotation = ann.clone().into();
            let back_again: vdj_ann::annotate::ContigAnnotation = crate_anns.into();
            assert!(ann == back_again);
        }
    }

    #[test]
    fn test_barcode_data_roundtrip() {
        let brief = vdj_asm_utils::barcode_data::BarcodeDataBrief {
            barcode: String::from("AAACCTGGTCCGAAGA-1"),
            read_pairs: 101,
            total_ucounts: 63,
            xucounts: vec![
                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4,
            ],
            ncontigs: 2,
            frac_reads_used: 1.0_f64,
        };
        let crate_barcode_data: crate::types::BarcodeData = brief.clone().into();
        let back_again: vdj_asm_utils::barcode_data::BarcodeDataBrief = crate_barcode_data.into();
        assert_eq!(brief, back_again);
    }

    #[test]
    fn test_vdj_reference_raw() -> Result<()> {
        let vdj_reference = VdjReferenceRaw::new(Path::new(
            "../dui_tests/test_resources/reference/vdj_GRCh38_alts_ensembl-4.0.0",
        ))?;

        let dir = tempfile::tempdir()?;
        let dir = dir.path();
        vdj_reference.write_to_folder(dir)?;

        assert!(file_diff::diff(
            "../dui_tests/test_resources/reference/vdj_GRCh38_alts_ensembl-4.0.0/fasta/regions.fa",
            dir.join("fasta/regions.fa").to_str().unwrap()
        ));

        assert!(file_diff::diff(
            "../dui_tests/test_resources/reference/vdj_GRCh38_alts_ensembl-4.0.0/reference.json",
            dir.join("reference.json").to_str().unwrap()
        ));
        Ok(())
    }

    #[test]
    fn test_fasta_hash() -> Result<()> {
        let vdj_reference = VdjReferenceRaw::new(Path::new(
            "../dui_tests/test_resources/reference/vdj_GRCh38_alts_ensembl-4.0.0",
        ))?;
        assert_eq!(
            vdj_reference.fasta_hash(),
            "17ddc1bb956ab6a33b90c9e4734e3471d8a3249e1f538aa6ce330d3765d60c92" // shasum -a 256
        );
        Ok(())
    }
}
