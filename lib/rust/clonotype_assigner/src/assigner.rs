//! Martian stage RUN_ENCLONE

use anyhow::Result;
use bio::alignment::{Alignment, AlignmentOperation};
use enclone_proto::proto_io::read_proto;
use enclone_proto::types::EncloneOutputs;
use enclone_ranger::main_enclone::main_enclone_ranger;
use martian::prelude::*;
use martian_derive::{make_mro, martian_filetype, MartianStruct};
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::FileTypeWrite;
use perf_stats::{elapsed, peak_mem_usage_gb};
use rust_htslib::bam;
use rust_htslib::bam::index;
use rust_htslib::bam::record::CigarString;
use serde::{Deserialize, Serialize};
use std::cmp::min;
use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::Instant;
use string_utils::TextUtils;
use vdj_ann::annotate::ContigAnnotation;
use vdj_asm_utils::filter_log::FilterSwitch;
use vdj_reference::VdjReceptor;
use vector_utils::next_diff;

pub(crate) const QUALITY_OFFSET: u8 = 33;

pub struct Assigner;

martian_filetype!(Lz4File, "lz4");
martian_filetype!(FaFile, "fa");
martian_filetype!(BamFile, "bam");
martian_filetype!(BamBaiFile, "bam.bai");
martian_filetype!(SamFile, "sam");
martian_filetype!(BinFile, "bin");
martian_filetype!(ProtoBinFile, "pb");

/// RUN_ENCLONE metrics
#[derive(Serialize)]
pub struct ClonotypeAssignerMetrics {
    /// Raw VDJ paired clonotype diversity
    multi_raw_vdj_paired_clonotype_diversity: f64,
    /// Raw CDRs per barcode histogram
    raw_cdrs_per_bc_histogram: HashMap<String, usize>,
}

#[derive(Debug, Clone, Deserialize, MartianStruct)]
pub struct ClonotypeAssignerStageInputs {
    pub filter_switch: FilterSwitch,
    pub vdj_reference_path: PathBuf,
    pub contig_annotations: JsonFile<Vec<ContigAnnotation>>,
    pub receptor: VdjReceptor, // TODO: denovo needs to be handled
}

#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct ClonotypeAssignerStageOutputs {
    pub summary: JsonFile<ClonotypeAssignerMetrics>,
    pub enclone_output: ProtoBinFile,
    pub donor_ref_fa: FaFile,
    pub barcode_fate: JsonFile<()>,
    // True if there are no cells
    pub disable_vloupe: bool,
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub(crate) fn bam_record_from_align(
    al: &Alignment,
    tigname: &str,
    refid: usize,
    seq: &[u8],
    quals: &[u8],
    barcode: &str,
) -> bam::Record {
    let cigar = CigarString::from_alignment(al, false);
    let mut rec = bam::record::Record::new();
    rec.set(tigname.as_bytes(), Some(&cigar), seq, quals);
    rec.push_aux(b"AS", bam::record::Aux::I32(al.score))
        .unwrap();
    let mut edit_distance = 0;
    for &op in &al.operations {
        edit_distance += match op {
            AlignmentOperation::Match => 0,
            AlignmentOperation::Xclip(_) => 0,
            AlignmentOperation::Yclip(_) => 0,
            _ => 1,
        };
    }
    rec.push_aux(b"NM", bam::record::Aux::I32(edit_distance))
        .unwrap();
    if !barcode.is_empty() {
        rec.push_aux(b"CB", bam::record::Aux::String(barcode))
            .unwrap();
    }
    rec.set_mapq(255_u8);
    rec.set_flags(136_u16);
    rec.set_tid(refid as i32);
    rec.set_pos(al.ystart as i64);
    rec.set_mpos(-1_i64);
    rec
}

pub(crate) fn replace_sam_by_indexed_bam(sam_name: &Path) {
    // Sort the sam file.  At the time of this writing, there is no rust code to sort sam or
    // bam files, so we have to call out to the samtools executable.  This is horrible.
    // See https://github.com/rust-bio/rust-htslib/issues/57.
    let sam_name = sam_name.to_str().unwrap();
    let sorted_sam_name = format!("{sam_name}.sorted");
    Command::new("samtools")
        .args(["sort", sam_name, "-o", &sorted_sam_name])
        .output()
        .expect("failed to sort sam file");

    // Convert the sam to bam, and delete the same file.

    let bam_name = format!("{}.bam", sam_name.rev_before(".sam"));
    Command::new("samtools")
        .args(["view", "-bS", "-o", &bam_name, &sorted_sam_name])
        .output()
        .expect("failed to sort sam file");
    let _ = std::fs::remove_file(sam_name);
    let _ = std::fs::remove_file(&sorted_sam_name);

    // Index the bam file.

    let bam_bai_name = format!("{bam_name}.bai");
    index::build(&bam_name, Some(&bam_bai_name), index::Type::Bai, 1).unwrap();
}

// The allowed mem_gb is based on observed mem 3.59 GB for lena 180517.
// Note that by using -4, in local mode, this forces this process to use the entire server.
// This is not necessarily desirable.
#[make_mro(mem_gb = 5, threads = -4, stage_name = RUN_ENCLONE)]
impl MartianMain for Assigner {
    type StageInputs = ClonotypeAssignerStageInputs;
    type StageOutputs = ClonotypeAssignerStageOutputs;
    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        // Set up logging.

        let t0 = Instant::now();

        // Cap number of threads.
        // let nthreads = max(4, min(rover.get_threads(), 96)); // NOTE: This could panic in local mode.
        let nthreads = min(rover.get_threads(), 96);

        // Call enclone.  We instruct it to generate a binary output file, and print nothing.
        // Need to bound number of threads used by enclone.
        let mut argsx = vec!["enclone".to_string()];
        let contig_annotations = args.contig_annotations.as_ref().to_str().unwrap();

        argsx.push(format!(
            "{}={}",
            match args.receptor {
                VdjReceptor::TR => "TCR",
                VdjReceptor::TRGD => "TCRGD",
                VdjReceptor::IG => "BCR",
            },
            contig_annotations.rev_before("/")
        ));
        let fasta_path = args
            .vdj_reference_path
            .join("fasta/regions.fa")
            .to_str()
            .unwrap()
            .to_string();
        let enclone_output: ProtoBinFile = rover.make_path("enclone_output");
        let donor_ref_fa: FaFile = rover.make_path("donor_ref");
        let barcode_fate: JsonFile<()> = rover.make_path("barcode_fate");
        argsx.push("PRE=".to_string());
        argsx.push(format!("REF={fasta_path}"));
        argsx.push("NOPRINT".to_string());
        argsx.push("CELLRANGER".to_string());
        argsx.push("NOPAGER".to_string());
        argsx.push("FORCE_EXTERNAL".to_string());
        argsx.push(format!("MAX_CORES={nthreads}"));
        argsx.push(format!(
            "PROTO={}",
            enclone_output.as_ref().to_str().unwrap()
        ));
        argsx.push(format!(
            "DONOR_REF_FILE={}",
            donor_ref_fa.as_ref().to_str().unwrap()
        ));
        argsx.push(format!(
            "FATE_FILE={}",
            barcode_fate.as_ref().to_str().unwrap()
        ));
        if args.receptor == VdjReceptor::TRGD {
            argsx.push("GAMMA_DELTA".to_string());
        }

        // Option to split clonotypes that have 4 chains or more.
        // These are most likely false joins due to a shared chain
        argsx.push("SPLIT_MAX_CHAINS=4".to_string());

        let FilterSwitch {
            asm_shared_contig: _,
            enclone_shared_contig,
            enclone_umi,
            enclone_multiplet,
        } = args.filter_switch;

        // ! because we want the filter to be off when the boolean is false
        if !enclone_shared_contig {
            argsx.push("NGRAPH_FILTER".to_string());
        }

        if !enclone_multiplet {
            argsx.extend(["NWEAK_CHAINS", "NFOURSIE_KILL", "NDOUBLET", "NSIG"].map(String::from));
        }

        if !enclone_umi {
            argsx.extend(["NUMI", "NUMI_RATIO"].map(String::from));
        }

        println!("Invoking enclone:\n{}", argsx.join(" "));
        main_enclone_ranger(&argsx).unwrap();
        println!("used {:.2} seconds up through enclone", elapsed(&t0));

        // Read back in the enclone output file.
        let enclone_outs: EncloneOutputs = read_proto(&enclone_output)?;

        // Compute metrics.

        let mut clonotype_sizes = Vec::<usize>::new();
        let mut total_cells = 0;
        let mut chain_counts = Vec::<usize>::new();
        for x in enclone_outs.clonotypes {
            let mut cells = 0;
            for j in 0..x.exact_clonotypes.len() {
                cells += x.exact_clonotypes[j].cell_barcodes.len();
            }
            for _ in 0..cells {
                chain_counts.push(x.chains.len());
            }
            clonotype_sizes.push(cells);
            total_cells += cells;
        }
        chain_counts.sort_unstable();

        let mut multi_raw_vdj_paired_clonotype_diversity = 0.0;
        for x in clonotype_sizes {
            multi_raw_vdj_paired_clonotype_diversity +=
                (x as f64 / total_cells as f64) * (x as f64 / total_cells as f64);
        }
        multi_raw_vdj_paired_clonotype_diversity = 1.0 / multi_raw_vdj_paired_clonotype_diversity;

        let mut raw_cdrs_per_bc_histogram = HashMap::new();
        let mut i = 0;
        while i < chain_counts.len() {
            let j = next_diff(&chain_counts, i);
            let n = chain_counts[i];
            if n <= 10 {
                raw_cdrs_per_bc_histogram.insert(n.to_string(), j - i);
            } else {
                raw_cdrs_per_bc_histogram.insert(">10".to_string(), chain_counts.len() - i);
                break;
            }
            i = j;
        }

        let summary: JsonFile<_> = rover.make_path("summary");
        summary.write(&ClonotypeAssignerMetrics {
            multi_raw_vdj_paired_clonotype_diversity,
            raw_cdrs_per_bc_histogram,
        })?;

        println!(
            "used {:.2} seconds total, peak mem = {:.2} GB",
            elapsed(&t0),
            peak_mem_usage_gb()
        );

        // Return outputs.
        Ok(ClonotypeAssignerStageOutputs {
            summary,
            enclone_output,
            disable_vloupe: total_cells == 0,
            donor_ref_fa,
            barcode_fate,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use enclone_proto::proto_io::ClonotypeIter;
    use std::collections::HashSet;
    use std::path::PathBuf;

    #[test]
    fn test_with_vdj_only_inputs() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let args = ClonotypeAssignerStageInputs {
            filter_switch: Default::default(),
            vdj_reference_path: PathBuf::from(
                "../dui_tests/test_resources/reference/vdj_GRCh38_alts_ensembl-7.1.0/",
            ),
            contig_annotations: JsonFile::from(
                "../dui_tests/test_resources/cellranger-vdj/vdj_tcell_tiny/contig_annotations.json",
            ),
            receptor: vdj_reference::VdjReceptor::TR,
        };
        let out = Assigner.test_run(&dir, args).unwrap();

        let mut cell_barcodes = HashSet::new();
        let mut num_clonotypes = 0;
        for cl in ClonotypeIter::from_file(out.enclone_output)? {
            cell_barcodes.extend(
                cl.exact_clonotypes
                    .clone()
                    .into_iter()
                    .flat_map(|ex_cl| ex_cl.cell_barcodes),
            );
            num_clonotypes += 1;
        }
        assert_eq!(cell_barcodes.len(), 13);
        assert_eq!(num_clonotypes, 13);
        Ok(())
    }
}
