//! AirrFilter stage code

use anyhow::Result;
use bio::io::fasta::Reader;
// The prelude brings the following items in scope:
// - Traits: MartianMain, MartianStage, RawMartianStage, MartianFileType, MartianMakePath
// - Struct/Enum: MartianRover, Resource, StageDef, MartianVoid,
//                Error (from failure crate), LevelFilter (from log crate)
// - Macros: martian_stages!
// - Functions: martian_main, martian_main_with_log_level, martian_make_mro
use martian::prelude::*;
// Bring the procedural macros in scope:
// #[derive(MartianStruct)], #[derive(MartianType)], #[make_mro], martian_filetype!
use martian_derive::{make_mro, martian_filetype, MartianStruct};
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::tabular_file::TsvFile;
use martian_filetypes::{LazyFileTypeIO, LazyWrite};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap};
use vdj_ann::annotate::ContigAnnotation;
use vdj_types::VdjRegion;

// NOTE: The following two structs will serve as the associated type for the
// trait. The struct fields need to be owned and are limited to
// - Basic int/float/bool/String types, PathBuf, Vec, Option, HashMap, HashSet
// - Structs/Enums implementing "AsMartianPrimaryType" (You can use #[derive(MartianType)])
// - Filetype (see the note below, representing as a filetype in mro)

// If you want to declare a new filetype use the `martian_filetype!` macro:
martian_filetype!(FastaFile, "fasta");

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct AirrFilterStageInputs {
    contig_annotations: JsonFile<Vec<ContigAnnotation>>,
    concat_ref_fasta: Option<FastaFile>,
    #[mro_type = "map"]
    pub gem_well_map: Option<BTreeMap<u32, (String, u32)>>,
}

#[derive(Debug, Clone, Serialize, Deserialize, MartianStruct)]
pub struct AirrFilterStageOutputs {
    airr_annotations: Option<TsvFile<Rearrangement>>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct CustomBool(String); // T for true and F for false

impl From<bool> for CustomBool {
    fn from(v: bool) -> CustomBool {
        if v {
            CustomBool("T".into())
        } else {
            CustomBool("F".into())
        }
    }
}

// This is our stage struct
pub struct AirrFilter;

/// See these pages fro details on the format:
/// - https://docs.airr-community.org/en/stable/datarep/rearrangements.html
/// - https://github.com/airr-community/airr-standards/blob/master/specs/airr-schema.yaml#L1774
///
/// The help text for the fields are copied from the links above
#[derive(Debug, Clone, Serialize, Deserialize)]
struct Rearrangement {
    // Identifier defining the cell of origin for the query sequence.
    // Cell ID (field: barcode)
    cell_id: String,
    // Clonal cluster assignment for the query sequence.
    // Clone ID (field: clonotype_id)
    clone_id: Option<String>,
    // Unique query sequence identifier within the file. Most often this will be the input sequence
    // header or a substring thereof, but may also be a custom identifier defined by the tool in
    // cases where query sequences have been combined in some fashion prior to alignment.
    // Sequence ID (field: contig_name)
    sequence_id: String,
    // The query nucleotide sequence. Usually, this is the unmodified input sequence, which may be
    // reverse complemented if necessary. In some cases, this field may contain consensus sequences
    // or other types of collapsed input sequences if these steps are performed prior to alignment.
    // Full-length NT sequence (field: sequence)
    sequence: String,
    // Amino acid translation of the query nucleotide sequence.
    // Full-length AA sequence (field: aa_sequence)
    sequence_aa: Option<String>,
    // Boolean values must be encoded as T for true and F for false.
    // Productive (field: productive)
    productive: CustomBool,
    // True if the alignment is on the opposite strand (reverse complemented) with respect to the
    // query sequence. If True then all output data, such as alignment coordinates and sequences,
    // are based on the reverse complement of ‘sequence’.
    // Reverse complement (set to False)
    rev_comp: CustomBool,
    // V gene assignment (field: annotations.[].feature.region_type==V-REGION &&
    // annotations.[].feature.gene_name)
    v_call: String,
    // V gene cigar string (field: annotations.[].cigar)
    v_cigar: String,
    // D gene assignment (proposed field: annotations.feature.region_type==D-REGION &&
    // annotations.feature.gene_name
    d_call: Option<String>,
    // D gene cigar string (field: annotations.[].cigar)
    d_cigar: Option<String>,
    // J gene assignment (field: annotations.feature.region_type==J-REGION &&
    // annotations.feature.gene_name
    j_call: String,
    // J gene cigar string (field: annotations.[].cigar)
    j_cigar: String,
    // C gene assignment (field: annotations.[].feature.region_type==C-REGION &&
    // annotations.[].feature.gene_name)
    c_call: Option<String>,
    // C gene cigar string (field: annotations.[].cigar)
    c_cigar: Option<String>,
    // Aligned portion of query sequence, including any indel corrections or numbering spacers,
    // such as IMGT-gaps. Typically, this will include only the V(D)J region, but that is not a
    // requirement.
    // Sequence alignment
    sequence_alignment: String,
    // Assembled, aligned, fully length inferred germline sequence spanning the same region as the
    // sequence_alignment field (typically the V(D)J region) and including the same set of
    // corrections and spacers (if any).
    // Germline alignment (TODO)
    germline_alignment: String,
    // Junction region nucleotide sequence, where the junction is defined as the CDR3 plus the two
    // flanking conserved codons.
    // NT seq of junction (field: cdr3_seq)
    junction: Option<String>,
    // AA seq of junction (field: cdr3)
    junction_aa: Option<String>,
    // How long is the NT seq of the junction (field: cdr3_seq.len)
    junction_length: usize,
    // How long is the AA seq of the junction (field: cdr3.len)
    junction_aa_length: Option<usize>,
    // Start position of the V segment in the query sequence (1-based closed interval).
    // Start of V gene on contig (field: annotations.[].contig_match_start)
    v_sequence_start: usize,
    // End position of the V segment in the query sequence (1-based closed interval).
    // End of V gene on contig (field: annotations.[].contig_match_end)
    v_sequence_end: usize,
    // Start position of the D segment in the query sequence (1-based closed interval).
    // Start of D gene on contig (field: annotations.[].contig_match_start)
    d_sequence_start: Option<usize>,
    // End position of the D segment in the query sequence (1-based closed interval).
    // End of D gene on contig (field: annotations.[].contig_match_end)
    d_sequence_end: Option<usize>,
    // Start position of the J segment in the query sequence (1-based closed interval).
    // Start of J gene on contig (field: annotations.[].contig_match_start)
    j_sequence_start: usize,
    // End position of the J segment in the query sequence (1-based closed interval).
    // End of J gene on contig (field: annotations.[].contig_match_end)
    j_sequence_end: usize,
    // Start of C gene on contig (field: annotations.[].contig_match_start)
    c_sequence_start: Option<usize>,
    // End of C gene on contig (field: annotations.[].contig_match_end)
    c_sequence_end: Option<usize>,
    // How many reads were allocated to this rearrangement (field: read_count)
    consensus_count: usize,
    // How many UMIs were allocated to this rearrangement (field: umi_count)
    duplicate_count: usize,
    // Is this rearrangement cell-associated (field: is_cell)
    is_cell: CustomBool,
    // Unique identifier to disambiguate source of rearrangement when combining multiple repertoires.
    // Populated with library_id in VDJ_AGGR runs.
    #[serde(skip_serializing_if = "Option::is_none")]
    repertoire_id: Option<String>,
}
#[make_mro(stage_name = CREATE_AIRR_TSV)]
impl MartianMain for AirrFilter {
    type StageInputs = AirrFilterStageInputs;
    type StageOutputs = AirrFilterStageOutputs;
    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        // Define FASTA reader for concat_ref
        let ref_reader = match args.concat_ref_fasta {
            Some(ref_fasta) => Reader::from_file(ref_fasta).unwrap(),
            None => {
                // Denovo no reference mode.
                return Ok(AirrFilterStageOutputs {
                    airr_annotations: None,
                });
            }
        };
        let mut concat_ref_map = HashMap::new();
        for record in ref_reader.records() {
            let record = record?;
            concat_ref_map.insert(
                record.id().replace("concat_ref", "consensus"),
                record.seq().to_vec(),
            );
        }

        // Define contig_ann reader
        let contig_reader = args.contig_annotations.lazy_reader()?;

        // Create TSV file in correct path; define path (each chunk has a path, only create files
        // in that path).
        let airr_annotations: TsvFile<_> = rover.make_path("airr_annotations");
        let mut airr_writer = airr_annotations.lazy_writer()?;

        for ann in contig_reader {
            let ann: ContigAnnotation = ann?; // Can you read from the JSON and create object of type ContigAnnotation
            if !ann.is_productive() || !ann.is_cell || ann.info.raw_consensus_id.is_none()
            // We don't generate consensus for "Multi" chains, so need to check if raw_consensus_id
            // is present
            {
                continue;
            }

            let v_region = ann.get_region(VdjRegion::V).unwrap(); // Guaranteed to exist
            let d_region = ann.get_region(VdjRegion::D); // Could be None
            let j_region = ann.get_region(VdjRegion::J).unwrap(); // Guaranteed to exist
            let c_region = ann.get_region(VdjRegion::C); // Could be None

            // Look for raw_consensus_id in hashmapped FASTA
            let germline_char = String::from_utf8(
                concat_ref_map[ann.info.raw_consensus_id.as_ref().unwrap()].clone(),
            )?;

            let repertoire_id = args.gem_well_map.as_ref().map(|gw_map| {
                let gw = ann.barcode.parse::<barcode::Barcode>().unwrap().gem_group() as u32;
                let (library_id, _) = gw_map.get(&gw).unwrap();
                library_id.to_string()
            });

            airr_writer.write_item(&Rearrangement {
                cell_id: ann.barcode.clone(),
                clone_id: ann.info.raw_clonotype_id.clone(),
                rev_comp: false.into(),
                sequence_id: ann.contig_name.clone(), // A barcode usually has >1 rearrangement, so we use contig_name
                sequence: ann.sequence.clone(),
                sequence_aa: ann.aa_sequence.clone(),
                productive: ann.productive.unwrap().into(),

                v_call: v_region.feature.gene_name.clone(),
                v_cigar: v_region.cigar.clone(),
                v_sequence_start: v_region.contig_match_start + 1, // 1-based index
                v_sequence_end: v_region.contig_match_end,

                d_call: d_region.map(|ann| ann.feature.gene_name.clone()),
                d_cigar: d_region.map(|ann| ann.cigar.clone()),
                d_sequence_start: d_region.map(|ann| ann.contig_match_start + 1), // 1-based index
                d_sequence_end: d_region.map(|ann| ann.contig_match_end),

                j_call: j_region.feature.gene_name.clone(),
                j_cigar: j_region.cigar.clone(),
                j_sequence_start: j_region.contig_match_start + 1, // 1-based index
                j_sequence_end: j_region.contig_match_end,

                c_call: c_region.map(|ann| ann.feature.gene_name.clone()),
                c_cigar: c_region.map(|ann| ann.cigar.clone()),
                c_sequence_start: c_region.map(|ann| ann.contig_match_start + 1), // 1-based index
                c_sequence_end: c_region.map(|ann| ann.contig_match_end),

                sequence_alignment: ann.sequence,
                germline_alignment: germline_char,
                junction: ann.cdr3_seq.clone(),
                junction_aa: ann.cdr3.clone(),
                duplicate_count: ann.umi_count,
                consensus_count: ann.read_count,
                junction_length: ann.cdr3_seq.unwrap().len(),
                junction_aa_length: Some(ann.cdr3.unwrap().len()),
                is_cell: ann.is_cell.into(),
                repertoire_id,
            })?;
        }
        airr_writer.finish()?;

        Ok(AirrFilterStageOutputs {
            airr_annotations: Some(airr_annotations),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use insta::{assert_ron_snapshot, assert_snapshot};
    use martian_filetypes::FileTypeRead;

    #[test]
    fn test_micro_tsv_snapshot() {
        let args = AirrFilterStageInputs {
            contig_annotations: JsonFile::from("test.json"),
            concat_ref_fasta: Some(FastaFile::from("test.fasta")),
            gem_well_map: None,
        };
        let run_dir = tempfile::tempdir().unwrap();
        let outs = AirrFilter.test_run(&run_dir, args).unwrap();
        let airr_annotations = outs.airr_annotations.unwrap();
        let tsv_contents: Vec<Rearrangement> = airr_annotations.read().unwrap();
        assert_ron_snapshot!(tsv_contents);
        assert_snapshot!(std::fs::read_to_string(&airr_annotations).unwrap());
    }
}
