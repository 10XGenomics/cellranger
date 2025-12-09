// This is code for a rust-only implementation of stage ASSEMBLE_VDJ.
#![expect(missing_docs)]
//
// ◼ The way we generate a bam file here doesn't make sense.  In the current
// ◼ implementation, we build sorted bam files, then merge them, which involves
// ◼ interweaving the files.  It seems like we could avoid the interweaving.  For
// ◼ each chunk, we could create two sorted bam files, one for the placed reads and
// ◼ one for the unplaced reads.  Then in the join we should be able to
// ◼ 'samtools reheader' the bam files, then 'samtools cat' all the placed bam
// ◼ files, followed by the unplaced bam files, and that should be sorted without
// ◼ doing any interweaving.

use crate::assembly_types::{
    AssemblyStageInputs, BamBaiFile, BamFile, ContigSummaryRow, FastqFile, UmiList, UmiSummaryRow,
};
use crate::contig_aligner::ContigAligner;
use anyhow::{Context, Error, Result, bail};
use cr_types::LibraryType;
use cr_types::rna_read::RnaRead;
use debruijn::Mer;
use debruijn::dna_string::DnaString;
use debruijn::kmer::Kmer20;
use io_utils::{fwrite, fwriteln};
use itertools::Itertools;
use kmer_lookup::make_kmer_lookup_20_single;
use martian::{MartianFileType, MartianRover, MartianStage, Resource, StageDef};
use martian_derive::{MartianStruct, make_mro, martian_filetype};
use martian_filetypes::bin_file::{BinaryFormat, BincodeFile};
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::lz4_file::Lz4;
use martian_filetypes::tabular_file::{TsvFile, TsvFileNoHeader};
use martian_filetypes::{FileTypeRead, FileTypeWrite, LazyFileTypeIO, LazyWrite};
use metric::Histogram;
use rust_htslib::bam;
use rust_htslib::bam::HeaderView;
use serde::{Deserialize, Serialize};
use std::collections::HashSet;
use std::fs::{File, read_to_string, remove_file, rename};
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::{env, fs};
use string_utils::TextUtils;
use vdj_ann::annotate::{ContigAnnotation, chain_type};
use vdj_ann::refx::{RefData, make_vdj_ref_data_core};
use vdj_asm_utils::asm::write_sam_record_simple;
use vdj_asm_utils::barcode_data::{BarcodeData, BarcodeDataBrief};
use vdj_asm_utils::constants::{
    CHAIN_TYPESX, CLIP, GAP_EXTEND, GAP_OPEN, KMER_LEN_BANDED_ALIGN, MATCH_SCORE, MISMATCH_SCORE,
    OR_CHAIN_TYPES, ReadType, WINDOW_SIZE_BANDED_ALIGN,
};
use vdj_asm_utils::heuristics::Heuristics;
use vdj_asm_utils::log_opts::LogOpts;
use vdj_asm_utils::primers::{get_primer_exts, inner_primers, outer_primers};
use vdj_asm_utils::process::{ProcessBarcodeResult, process_barcode};
use vdj_asm_utils::umi_data::UmiInfo;
use vdj_asm_utils::{bam_utils, graph_read, sw};
use vdj_reference::VdjReceptor;
use vdj_types::{VdjChain, get_max_read_pairs_per_barcode};
pub struct Assembly;
martian_filetype!(Lz4File, "lz4");
martian_filetype!(TxtFile, "txt");

martian_filetype!(_BarcodeDataBriefFile, "bdf");
pub type BarcodeDataBriefFile = BinaryFormat<_BarcodeDataBriefFile, Vec<BarcodeDataBrief>>;

martian_filetype!(_BarcodeDataFile, "bd");
pub type BarcodeDataFile = BinaryFormat<_BarcodeDataFile, Vec<BarcodeData>>;

martian_filetype!(_UmiInfoFile, "uinfo");
pub type UmiInfoFile = BinaryFormat<_UmiInfoFile, Vec<UmiInfo>>;

// martian_filetype!(_ContigAnnFile, "can");
// impl FileStorage<Vec<ContigAnnotation>> for _ContigAnnFile {}
// pub type ContigAnnFormat = BinaryFormat<_ContigAnnFile>;

// =================================================================================
// FUNCTION TO MERGE BAMS
// =================================================================================

// merge coordinate-sorted bam files to yield a new coordinate-sorted bam file
//
// ◼ Duplicated code, with a couple of changes to get it to compile.  If
// ◼ we're going to use it, it shouldn't be in two places.
//
// ◼ sort_bam and index_bam also copied
//
// Note that this uses FOUR threads.  Not tested to determine if this makes it
// faster.

fn merge_bams(paths: &[BamFile], out: &Path) -> Result<()> {
    let mut paths = paths
        .iter()
        .map(|p| PathBuf::from(p.as_ref()))
        .collect::<Vec<_>>();
    // if there's only a single input path, copy the input
    if paths.len() == 1 {
        let _ = fs::copy(&paths[0], out)?;
        return Ok(());
    }
    // get the NOFILE ulimit so we don't ask samtools to open too many files
    let rlim = rustix::process::getrlimit(rustix::process::Resource::Nofile)
        .current
        .unwrap() as usize;
    assert!(rlim >= 102, "soft NOFILE ulimit is unworkably low");
    let rlim = rlim - 100;
    // keep merging files, rlim at a time, until only 1 remains
    let mut created = vec![];
    while paths.len() > 1 {
        let mut fofn = PathBuf::from(out);
        fofn.set_extension("fofn");
        let fofn = fofn.to_str().unwrap();
        {
            let mut handle = File::create(fofn)?;
            let n = paths.len();
            let rest = paths.split_off(rlim.min(n));
            for path in paths {
                handle.write_all(path.to_str().unwrap().as_bytes())?;
                handle.write_all(b"\n")?;
            }
            paths = rest;
        }
        let mut bam = PathBuf::from(out);
        bam.set_extension(format!("{}.bam", created.len()));
        created.push(bam.clone());
        paths.push(bam.clone());
        Command::new("samtools")
            .args([
                "merge",
                "-@",
                "3",
                "-c",
                "-p",
                "-s",
                "0",
                "-b",
                fofn,
                bam.to_str().unwrap(),
            ])
            .output()
            .expect("failed to merge bam files!");
        remove_file(fofn)?;
    }
    // move the final file into position
    rename(&paths[0], out)?;
    // remove all other intermediate bam files
    for bam in created {
        if bam != paths[0] {
            remove_file(bam)?;
        }
    }
    Ok(())
}

// modified to use only one thread, and to reduce memory usage

fn sort_bam(input: &str, output: &str) {
    println!("running samtools sort -l 8G -o {output} {input}");
    Command::new("samtools")
        .args(["sort", "-l", "8G", "-m", "600M", "-o", output, input])
        .output()
        .unwrap_or_else(|_| panic!("failed to sort {}", &input));
}

// Note that index_bam uses FOUR threads.  Not tested to determine if this makes
// it faster.

fn index_bam(bam: &Path) {
    Command::new("samtools")
        .args(["index", "-@", "3", bam.to_str().unwrap()])
        .output()
        .unwrap_or_else(|_| panic!("failed to index {}", bam.display()));
}

pub fn line_by_line_copy<R, W>(reader: &mut R, writer: &mut BufWriter<W>) -> Result<()>
where
    R: BufRead,
    W: Write,
{
    for line in reader.lines() {
        writeln!(writer, "{}", line?)?;
    }
    Ok(())
}

// =================================================================================
// FIND ENRICHMENT PRIMERS
// =================================================================================

fn enrichment_primers(
    primer_file: Option<&Path>,
    refdata: &RefData,
    is_tcr: bool,
    is_bcr: bool,
    inner_primersx: &mut Vec<Vec<u8>>,
    outer_primersx: &mut Vec<Vec<u8>>,
) {
    // Specify inner primers.  If the customer has not specified primers, we use the reference
    // sequence to decide if the species is human or mouse.

    let iprimers = match primer_file {
        Some(path) => BufReader::new(File::open(path).unwrap())
            .lines()
            .map(Result::unwrap)
            .collect(),
        None => Vec::<String>::new(),
    };
    if iprimers.is_empty() {
        let (mut is_human, mut is_mouse) = (false, false);
        let (mut human_count, mut mouse_count) = (0, 0);
        let mut human_inner_primers = inner_primers("human", "tcr");
        human_inner_primers.append(&mut inner_primers("human", "bcr"));
        let mut mouse_inner_primers = inner_primers("mouse", "tcr");
        mouse_inner_primers.append(&mut inner_primers("mouse", "bcr"));
        for i in 0..refdata.refs.len() {
            if !refdata.is_c(i) {
                continue;
            }
            let x = refdata.refs[i].clone().rc().to_string();
            for primer in &human_inner_primers {
                if x.contains(std::str::from_utf8(primer).unwrap()) {
                    human_count += 1;
                }
            }
            for primer in &mouse_inner_primers {
                if x.contains(std::str::from_utf8(primer).unwrap()) {
                    mouse_count += 1;
                }
            }
        }
        if human_count > 0 {
            is_human = true;
        }
        if mouse_count > 0 {
            is_mouse = true;
        }
        if is_human && is_tcr {
            inner_primersx.append(&mut inner_primers("human", "tcr"));
            outer_primersx.append(&mut outer_primers("human", "tcr"));
        }
        if is_human && is_bcr {
            inner_primersx.append(&mut inner_primers("human", "bcr"));
            outer_primersx.append(&mut outer_primers("human", "bcr"));
        }
        if is_mouse && is_tcr {
            inner_primersx.append(&mut inner_primers("mouse", "tcr"));
            outer_primersx.append(&mut outer_primers("mouse", "tcr"));
        }
        if is_mouse && is_bcr {
            inner_primersx.append(&mut inner_primers("mouse", "bcr"));
            outer_primersx.append(&mut outer_primers("mouse", "bcr"));
        }
    } else {
        for x in iprimers {
            let p = x.as_bytes().to_vec();
            inner_primersx.push(p);
        }
    }
}

fn sam_to_bam(out_bam_file: &Path, sam_header: bam::header::Header, out_sam_filenamex: &Path) {
    // Convert sam to bam.
    // ◼ The whole business of first writing sam.lz4, then converting
    // ◼ to bam here seems nuts and inefficient.

    let out_bam_filename_str = out_bam_file.to_str().unwrap();
    let mut out = bam::Writer::from_path(out_bam_file, &sam_header, bam::Format::Bam).unwrap();
    let h = HeaderView::from_header(&sam_header);
    let fin = File::open(out_sam_filenamex).expect("Failed to open file for reading");
    let mut fin = lz4::Decoder::new(fin).expect("Failed to create lz4 decoder");
    let fin = BufReader::new(&mut fin);
    for line in fin.lines() {
        let s = line.unwrap();
        let t = s.as_bytes();
        let rec: bam::Record = bam::Record::from_sam(&h, t).unwrap();
        out.write(&rec).unwrap();
    }
    drop(out);
    drop(sam_header);
    drop(h);

    // Sort the bam file.  This can use a lot of memory.  We delete a temporary file if it
    // exists in case you're debugging this code, as if it's left around it will crash
    // samtools.

    let out_bam_sorted_filename =
        out_bam_filename_str.rev_before("/").to_string() + "/contig_bam_sorted.bam";
    let tmp_filename = out_bam_sorted_filename.clone() + ".tmp.0000.bam";
    if Path::new(&tmp_filename).exists() {
        remove_file(&tmp_filename).unwrap();
    }
    sort_bam(out_bam_filename_str, &out_bam_sorted_filename);
    rename(&out_bam_sorted_filename, out_bam_file).unwrap();
}

// =================================================================================
// DEFINE THE STAGE INPUTS AND OUTPUTS
// =================================================================================

#[derive(Debug, Serialize, Deserialize, MartianStruct)]
pub struct AssemblyStageOutputs {
    pub contig_bam: BamFile,
    pub contig_bam_bai: BamBaiFile,
    pub summary_tsv: TsvFile<ContigSummaryRow>,
    pub umi_summary_tsv: TsvFile<UmiSummaryRow>,
    pub umi_info: UmiInfoFile,
    pub contig_annotations: JsonFile<Vec<ContigAnnotation>>,
    pub barcode_brief: BarcodeDataBriefFile,
    pub barcode_full: BarcodeDataFile,
    pub barcodes_in_chunks: Vec<JsonFile<Vec<String>>>,

    // The outputs below are simply bubbled up to the outs folder
    pub align_info: TxtFile,
    pub unmapped_sample_fastq: FastqFile,
}

// =================================================================================
// DEFINE THE CHUNK INPUTS AND OUTPUTS
// =================================================================================

#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct AssemblyChunkInputs {
    pub chunk_rna_reads: Lz4<BincodeFile<Vec<RnaRead>>>,
    pub chunk_id: usize,
    pub primers: VdjPrimers,
}

#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct AssemblyChunkOutputs {
    pub contig_bam: BamFile,
    pub summary_tsv: TsvFileNoHeader<ContigSummaryRow>,
    pub umi_summary_tsv: TsvFileNoHeader<UmiSummaryRow>,
    pub contig_annotations: JsonFile<Vec<ContigAnnotation>>,
    pub barcodes_in_chunk: JsonFile<Vec<String>>,
    pub align_info: TxtFile,
    pub unmapped_sample_fastq: FastqFile,
    pub barcode_data: BarcodeDataFile,
    pub umi_info: UmiInfoFile,
}

#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct VdjPrimers {
    pub inner_primers: Vec<Vec<u8>>,
    pub outer_primers: Vec<Vec<u8>>,
}

impl VdjPrimers {
    pub fn new(
        ref_path: Option<PathBuf>,
        primer_file: Option<PathBuf>,
        receptor: Option<VdjReceptor>,
    ) -> Result<Self, Error> {
        let is_tcr = receptor == Some(VdjReceptor::TR) || receptor == Some(VdjReceptor::TRGD);
        let is_bcr = receptor == Some(VdjReceptor::IG);
        let mut refdata = RefData::new();
        if let Some(ref ref_path) = ref_path {
            let fasta_path = ref_path.join("fasta/regions.fa");
            let fasta = read_to_string(&fasta_path)
                .with_context(|| fasta_path.to_string_lossy().to_string())?;
            make_vdj_ref_data_core(&mut refdata, &fasta, "", is_tcr, is_bcr, None);
        }
        let mut inner_primersx = Vec::<Vec<u8>>::new();
        let mut outer_primersx = Vec::<Vec<u8>>::new();
        enrichment_primers(
            primer_file.as_deref(),
            &refdata,
            is_tcr,
            is_bcr,
            &mut inner_primersx,
            &mut outer_primersx,
        );

        Ok(VdjPrimers {
            inner_primers: inner_primersx,
            outer_primers: outer_primersx,
        })
    }
}

// =================================================================================
// STAGE CODE BOILERPLATE
// =================================================================================

#[make_mro(stage_name = ASSEMBLE_VDJ)]
impl MartianStage for Assembly {
    type StageInputs = AssemblyStageInputs;
    type StageOutputs = AssemblyStageOutputs;
    type ChunkInputs = AssemblyChunkInputs;
    type ChunkOutputs = AssemblyChunkOutputs;

    // =============================================================================
    // THIS IS THE SPLIT CODE
    // =============================================================================

    fn split(
        &self,
        args: Self::StageInputs,
        _rover: MartianRover,
    ) -> Result<StageDef<Self::ChunkInputs>> {
        // If the customer has not specified primers, we use the
        // reference sequence to decide if the species is human or mouse.
        let primers = VdjPrimers::new(
            args.vdj_reference_path,
            args.inner_enrichment_primers,
            args.receptor,
        )?;
        // Set up chunks.
        // ◼ Join memory highwater mark was 5.4 GB (rounded).
        // ◼ See comments about memory inefficiency in the join step.
        // 4 threads in join for `samtools merge/index`
        Ok(args
            .bc_sorted_rna_reads
            .into_iter()
            .enumerate()
            .map(|(i, chunk_rna_reads)| {
                (
                    AssemblyChunkInputs {
                        chunk_rna_reads,
                        chunk_id: i,
                        primers: primers.clone(),
                    },
                    Resource::with_mem_gb(2),
                )
            })
            .collect::<StageDef<_>>()
            .join_resource(Resource::with_mem_gb(6).threads(4)))
    }

    // =============================================================================
    // THIS IS THE CHUNK CODE
    // =============================================================================

    fn main(
        &self,
        args: Self::StageInputs,
        split_args: Self::ChunkInputs,
        rover: MartianRover,
    ) -> Result<Self::ChunkOutputs> {
        // Print the command.
        println!("{}", env::args().collect::<Vec<_>>().join(" "));

        // Load results from assembly prep stage.

        println!("n50_n50_rpu = {}", args.n50_n50_rpu);
        let is_tcr =
            args.receptor == Some(VdjReceptor::TR) || args.receptor == Some(VdjReceptor::TRGD);
        let is_bcr = args.receptor == Some(VdjReceptor::IG);
        let is_gd = Some(args.receptor == Some(VdjReceptor::TRGD));
        let (refdata, refdatax, refdata_full, rkmers_plus_full_20) =
            load_refdata(args.vdj_reference_path.as_deref(), is_tcr, is_bcr)?;
        let refs = &refdata.refs;

        // Get filenames and set up writers.

        let contig_annotations_file: JsonFile<Vec<ContigAnnotation>> =
            rover.make_path("contig_annotations");

        let umi_summary_file: TsvFileNoHeader<UmiSummaryRow> = rover.make_path("umi_summary");

        let summary_file: TsvFileNoHeader<ContigSummaryRow> = rover.make_path("summary");

        // Set up to write a sam file.

        let out_sam_filenamex: Lz4File = rover.make_path("contig_sam.lz4");

        // Start of new code to process all barcodes.

        let log_opts = LogOpts::default();

        // Set up for alignment.

        let align_info_file: TxtFile = rover.make_path("align_info");

        let unmapped_sample_fastq_file: FastqFile = rover.make_path("unmapped_sample_fastq");

        let (sam_header, barcodes, barcode_data, umi_info) = write_simple_sam(
            &args,
            &split_args,
            &out_sam_filenamex,
            &contig_annotations_file,
            &align_info_file,
            &umi_summary_file,
            &summary_file,
            &unmapped_sample_fastq_file,
            &refdata,
            refdatax,
            refdata_full,
            rkmers_plus_full_20,
            is_tcr,
            is_bcr,
            is_gd,
            !refs.is_empty(),
            &log_opts,
        )?;

        // Write barcode data.

        let barcode_data_file: BarcodeDataFile = rover.make_path("barcode_data");
        barcode_data_file.write(&barcode_data)?;
        drop(barcode_data);
        let umi_info_file: UmiInfoFile = rover.make_path("umi_info");
        umi_info_file.write(&umi_info)?;
        drop(umi_info);

        let out_bam_filename: BamFile = rover.make_path("contig_bam");
        sam_to_bam(&out_bam_filename, sam_header, &out_sam_filenamex);
        remove_file(&out_sam_filenamex).unwrap();

        // Create barcodes_in_chunk.json.

        let barcodes_in_chunk_file: JsonFile<_> = rover.make_path("barcodes_in_chunk");
        barcodes_in_chunk_file.write(&barcodes)?;
        drop(barcodes);

        Ok(AssemblyChunkOutputs {
            contig_bam: out_bam_filename,
            summary_tsv: summary_file,
            umi_summary_tsv: umi_summary_file,
            contig_annotations: contig_annotations_file,
            barcodes_in_chunk: barcodes_in_chunk_file,
            align_info: align_info_file,
            unmapped_sample_fastq: unmapped_sample_fastq_file,
            barcode_data: barcode_data_file,
            umi_info: umi_info_file,
        })
    }

    // =============================================================================
    // THIS IS THE JOIN CODE
    // =============================================================================

    fn join(
        &self,
        args: Self::StageInputs,
        _chunk_defs: Vec<Self::ChunkInputs>,
        chunk_outs: Vec<Self::ChunkOutputs>,
        rover: MartianRover,
    ) -> Result<Self::StageOutputs> {
        // CELLRANGER-7889: "VDJ" is hardcoded in mro injection of chemistry defs map.
        let chemistry_def = &args.chemistry_defs[&LibraryType::VdjAuto];

        // Merge summary_tsv and umi_summary_tsv files.
        let summary_tsv_file: TsvFile<ContigSummaryRow> = rover.make_path("summary_tsv");
        write_contig_summary_tsv(&summary_tsv_file, &chunk_outs)?;

        let umi_summary_tsv_file: TsvFile<UmiSummaryRow> = rover.make_path("umi_summary_tsv");
        write_umi_summary_tsv(&umi_summary_tsv_file, &chunk_outs)?;

        // Merge align info files.

        let align_info_file: TxtFile = rover.make_path("align_info");
        let mut align_info = align_info_file.buf_writer()?;
        for co in &chunk_outs {
            line_by_line_copy(&mut co.align_info.buf_reader()?, &mut align_info)?;
        }
        drop(align_info);

        // Merge unmapped fastq files.

        let unmapped_sample_fastq_file: FastqFile = rover.make_path("unmapped_sample_fastq");
        let mut unmapped_sample_fastq = unmapped_sample_fastq_file.buf_writer()?;
        for co in &chunk_outs {
            line_by_line_copy(
                &mut co.unmapped_sample_fastq.buf_reader()?,
                &mut unmapped_sample_fastq,
            )?;
        }
        drop(unmapped_sample_fastq);

        // Load the number of read pairs assigned to each barcode.

        let bc_counts = args.corrected_bc_counts.as_ref().to_str().unwrap();
        let reader = BufReader::new(File::open(bc_counts).unwrap());
        let npairs: serde_json::Value = serde_json::from_reader(reader)?;

        // Merge contig annotation json files.  This also makes some adjustments:
        // 1. Add read count for barcode.
        // 2. Filter unrelated chains based on Gamma/Delta mode
        //      -> in GD mode, filter all AB (NOT VICE VERSA)
        //      -> in AB mode, no filter is applied to maintain consistency
        let gd_mode = args.receptor == Some(VdjReceptor::TRGD);

        let contig_annotations_file: JsonFile<_> = rover.make_path("contig_annotations");
        let mut ann_writer = contig_annotations_file.lazy_writer()?;
        for co in &chunk_outs {
            // In GD mode, take a pass to identify cells which have at least one productive
            // G/D chain

            let mut gd_barcodes = HashSet::new();
            if gd_mode {
                for can in co.contig_annotations.lazy_reader()? {
                    let can: ContigAnnotation = can?;
                    if can.productive.unwrap_or(false)
                        && can
                            .annotations
                            .iter()
                            .any(|ann| matches!(ann.feature.chain, VdjChain::TRG | VdjChain::TRD))
                    {
                        gd_barcodes.insert(can.barcode);
                    }
                }
            }

            // second pass
            let max_read_pairs_per_barcode = get_max_read_pairs_per_barcode(
                chemistry_def.is_paired_end(),
                args.max_reads_per_barcode,
            );
            let reader2 = co.contig_annotations.lazy_reader()?;
            for can in reader2 {
                let mut can: ContigAnnotation = can?;
                let used = npairs[&can.barcode].as_u64().unwrap() as usize;
                let frac = if used <= max_read_pairs_per_barcode {
                    1.0_f64
                } else {
                    max_read_pairs_per_barcode as f64 / used as f64
                };
                can.fraction_of_reads_for_this_barcode_provided_as_input_to_assembly = Some(frac);
                ann_writer.write_item(&can)?;
            }
        }
        ann_writer.finish()?;

        // Merge bam files, then index.

        let out_bam_filename: BamFile = rover.make_path("contig_bam.bam");

        let bams: Vec<_> = chunk_outs.iter().map(|co| co.contig_bam.clone()).collect();
        let contig_bam_bai_filename: BamBaiFile = rover.make_path("contig_bam.bam.bai");
        if !chunk_outs.is_empty() {
            merge_bams(&bams, out_bam_filename.as_ref()).unwrap();

            index_bam(&out_bam_filename);
        } else {
            out_bam_filename.buf_writer()?;
        }

        let mut refdata = RefData::new();

        let is_tcr =
            args.receptor == Some(VdjReceptor::TR) || args.receptor == Some(VdjReceptor::TRGD);
        let is_bcr = args.receptor == Some(VdjReceptor::IG);

        if let Some(ref ref_path) = args.vdj_reference_path {
            let fasta_path = ref_path.join("fasta/regions.fa");
            let fasta = read_to_string(&fasta_path)
                .with_context(|| fasta_path.to_string_lossy().to_string())?;
            make_vdj_ref_data_core(&mut refdata, &fasta, "", is_tcr, is_bcr, None);
        }

        let barcode_data_brief_file: BarcodeDataBriefFile = rover.make_path("barcode_data_brief");
        write_bc_databrief(&barcode_data_brief_file, &chunk_outs)?;

        let barcode_data_full_file: BarcodeDataFile = rover.make_path("barcode_data_full");
        write_bc_data_full(&barcode_data_full_file, &chunk_outs)?;

        let umi_info_full_file: UmiInfoFile = rover.make_path("umi_info");
        write_umi_info_full(&umi_info_full_file, &chunk_outs)?;

        // Return results.

        let barcodes_in_chunks = chunk_outs
            .iter()
            .map(|co| co.barcodes_in_chunk.clone())
            .collect();

        Ok(AssemblyStageOutputs {
            contig_bam: out_bam_filename,
            contig_bam_bai: contig_bam_bai_filename,
            summary_tsv: summary_tsv_file,
            umi_summary_tsv: umi_summary_tsv_file,
            umi_info: umi_info_full_file,
            contig_annotations: contig_annotations_file,
            barcode_brief: barcode_data_brief_file,
            barcode_full: barcode_data_full_file,
            barcodes_in_chunks,
            align_info: align_info_file,
            unmapped_sample_fastq: unmapped_sample_fastq_file,
        })
    }
}

#[allow(clippy::type_complexity)]
fn load_refdata(
    vdj_reference_path: Option<&Path>,
    is_tcr: bool,
    is_bcr: bool,
) -> Result<(RefData, RefData, RefData, Vec<(Kmer20, i32, i32)>)> {
    // Load reference and make a lookup table for it.  Actually there are three
    // versions:
    // (1) just for TCR or BCR (so long as we know which we have);
    // (2) for TCR and BCR, and with extra k=20 lookup table;
    // (3) just for TCR or BCR but also with extra sequences thrown in.

    let mut refdata = RefData::new();
    let mut refdatax = RefData::new();
    let mut refdata_full = RefData::new();
    let mut rkmers_plus_full_20 = Vec::<(Kmer20, i32, i32)>::new();
    if let Some(ref_path) = vdj_reference_path {
        let ref_path = ref_path.to_str().unwrap();
        let fasta = read_regions(Path::new(ref_path))?;
        let ext_fasta = match read_to_string(format!("{ref_path}/fasta/supp_regions.fa")) {
            Err(err) if err.kind() == io::ErrorKind::NotFound => String::new(),
            Err(err) => bail!(err),
            Ok(contents) => contents,
        };

        make_vdj_ref_data_core(&mut refdata, &fasta, "", is_tcr, is_bcr, None);
        make_vdj_ref_data_core(&mut refdata_full, &fasta, "", true, true, None);
        make_kmer_lookup_20_single(&refdata_full.refs, &mut rkmers_plus_full_20);
        make_vdj_ref_data_core(&mut refdatax, &fasta, &ext_fasta, is_tcr, is_bcr, None);
    }
    Ok((refdata, refdatax, refdata_full, rkmers_plus_full_20))
}

fn make_rtype(refdata_full: &RefData, has_refs: bool) -> Vec<i32> {
    if has_refs {
        refdata_full
            .rheaders
            .iter()
            .take(refdata_full.refs.len())
            .map(|header| {
                CHAIN_TYPESX
                    .iter()
                    .enumerate()
                    .filter_map(|(j, chain_type)| {
                        if header.contains(chain_type) {
                            Some(j as i32)
                        } else {
                            None
                        }
                    })
                    .next_back()
                    .unwrap_or(-1)
            })
            .collect()
    } else {
        Vec::<i32>::new()
    }
}

#[allow(clippy::too_many_arguments)]
#[allow(clippy::type_complexity)]
fn write_simple_sam(
    args: &AssemblyStageInputs,
    split_args: &AssemblyChunkInputs,
    out_sam_filenamex: &Lz4File,
    contig_annotations_file: &JsonFile<Vec<ContigAnnotation>>,
    align_info_file: &TxtFile,
    umi_summary_file: &TsvFileNoHeader<UmiSummaryRow>,
    summary_file: &TsvFileNoHeader<ContigSummaryRow>,
    unmapped_sample_fastq_file: &FastqFile,
    refdata: &RefData,
    refdatax: RefData,
    refdata_full: RefData,
    rkmers_plus_full_20: Vec<(Kmer20, i32, i32)>,
    is_tcr: bool,
    is_bcr: bool,
    is_gd: Option<bool>,
    has_refs: bool,
    log_opts: &LogOpts,
) -> Result<
    (
        bam::header::Header,
        Vec<String>,
        Vec<BarcodeData>,
        Vec<UmiInfo>,
    ),
    Error,
> {
    // CELLRANGER-7889: "VDJ" is hardcoded in mro injection of chemistry defs map.
    let chemistry_def = &args.chemistry_defs[&LibraryType::VdjAuto];
    let rtype = make_rtype(&refdata_full, has_refs);
    // Initialize heuristics.

    let mut heur = Heuristics::new();
    if args.denovo {
        heur.free = true;
    }
    // Null stuff.
    let n50_n50_rpu = args.n50_n50_rpu;
    // Compute primer extensions.
    let inner_primersx = &split_args.primers.inner_primers;
    let outer_primersx = &split_args.primers.outer_primers;
    let inner_primer_exts = get_primer_exts(inner_primersx, refdata);
    let outer_primer_exts = get_primer_exts(outer_primersx, refdata);

    let mut log = Vec::<u8>::new();
    let mut ann_writer = contig_annotations_file.lazy_writer()?;
    let mut align_info = align_info_file.buf_writer()?;
    let mut umi_summary_writer = umi_summary_file.lazy_writer()?;
    let mut summary_writer = summary_file.lazy_writer()?;
    let mut unmapped_sample_fastq = unmapped_sample_fastq_file.buf_writer()?;
    let vdj_adapters = crate::adapter::get_vdj_adapters();

    // Track cell barcodes, barcodes, and barcode data.

    let mut barcodes = Vec::<String>::new();
    let mut barcode_data = Vec::<BarcodeData>::new();
    // Sam header.

    let mut sam_header = bam::header::Header::new();
    // Determine the id of this chunk.

    let ch = split_args.chunk_id;

    // Get number of read pairs.
    let npairs = args.npairs as usize;
    // Determine if single end.
    let single_end = !chemistry_def.is_paired_end();
    // Define fraction of pairs to align to determine chain type.
    // target_frac = number out of 1000 to keep

    const TARGET_READS: usize = 1_000_000;
    let target_pairs = if single_end {
        TARGET_READS
    } else {
        TARGET_READS / 2
    };
    let target_frac = if target_pairs < npairs {
        (1000 * target_pairs) / npairs
    } else {
        1000
    };
    // Create sam writer.

    let mut simple_sam_writer = lz4::EncoderBuilder::new()
        .build(BufWriter::new(
            File::create(out_sam_filenamex).expect("could not open sam file for writing"),
        ))
        .expect("could not start building lz4");

    // Scope to force closure of simple_sam_writer2.

    {
        let barcode_counts_full = args.corrected_bc_counts.read()?;
        let mut trimmer = crate::adapter::VdjTrimmer::new(&vdj_adapters);
        let mut simple_sam_writer2 = BufWriter::new(&mut simple_sam_writer);

        let bc_sorted_lazy_reader = split_args.chunk_rna_reads.lazy_reader()?;

        // Loop over all barcodes in the chunk.

        let (mut bid, mut rid) = (0, 0);
        let mut unmapped = 0;
        for (barcode, read_iter) in &bc_sorted_lazy_reader
            .chunk_by(|read: &Result<RnaRead>| read.as_ref().ok().map(RnaRead::barcode))
        {
            let mut barcode_data_this = BarcodeData::new();

            let mut this_bc_reads = Vec::new();
            for rna_read in read_iter {
                let mut rna_read = rna_read?;
                trimmer.trim(&mut rna_read);
                this_bc_reads.push(rna_read);
            }
            let barcode = barcode.unwrap();
            // Vec<(umi, seq, qual, readname, flags)>
            let read_inner_data =
                crate::translator::make_read_data(&this_bc_reads, 8, chemistry_def);
            // (barcode, Vec<(umi, seq, qual, readname, flags)>, actual reads)
            let barcode_string = barcode.to_string();
            let actual_reads = barcode_counts_full.get(&barcode_string);
            let read_data = (barcode_string, read_inner_data, actual_reads);
            barcode_data_this.nreads = read_data.2 as i32;
            barcode_data_this.nreads_used_for_assembly = this_bc_reads.len() as i32;
            // Assign fraction of reads used for assembly of each barcode
            let max_read_pairs_per_barcode = get_max_read_pairs_per_barcode(
                chemistry_def.is_paired_end(),
                args.max_reads_per_barcode,
            );
            barcode_data_this.frac =
                if barcode_data_this.nreads as usize <= max_read_pairs_per_barcode {
                    1.0_f64
                } else if barcode_data_this.nreads > 0 {
                    max_read_pairs_per_barcode as f64 / barcode_data_this.nreads as f64
                } else {
                    0.0_f64
                };

            // Set thread message.

            let chbid = format!("{ch}.{bid}");

            // Align a fixed fraction of the reads to determine their chain type.
            // This is annoying but it's expensive to align them all.

            if has_refs {
                for i in 0..read_data.1.len() {
                    if rid % 1000 < target_frac {
                        let b = read_data.1[i].1.clone(); // Read sequence
                        let mut best = chain_type(&b, &rkmers_plus_full_20, &rtype);
                        fwrite!(align_info, "{} ==> ", read_data.1[i].0); // UMI
                        if best >= 0 {
                            fwriteln!(align_info, "{}", OR_CHAIN_TYPES[best as usize]);
                        } else {
                            fwriteln!(align_info, "unmapped");
                            if unmapped % 50 == 0 {
                                fwriteln!(unmapped_sample_fastq, "@{}", read_data.1[i].3);
                                fwriteln!(
                                    unmapped_sample_fastq,
                                    "{}",
                                    read_data.1[i].1.to_string()
                                );
                                fwriteln!(unmapped_sample_fastq, "+");
                                let mut qual = read_data.1[i].2.clone();
                                for q in &mut qual {
                                    *q += 33;
                                }
                                fwriteln!(
                                    unmapped_sample_fastq,
                                    "{}",
                                    std::str::from_utf8(&qual).unwrap()
                                );
                            }
                            unmapped += 1;
                        }
                        if best == -1_i8 {
                            best = 14_i8;
                        }
                        barcode_data_this.chain_sample[best as usize] += 1;
                    }
                    rid += 1;
                }
            }
            unmapped_sample_fastq.flush()?;
            // Keep going.

            let barcode = read_data.0.clone();
            barcodes.push(barcode.clone());

            // if track { println!( "START {}", chbid ); }
            let mut log2 = Vec::<u8>::new();

            fwriteln!(
                log2,
                "\n▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓\
                         ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓\n"
            );
            fwriteln!(log2, "BARCODE {} = {}\n", chbid, read_data.0);
            fwriteln!(log2, "using {} reads", read_data.1.len());
            barcode_data_this.barcode.clone_from(&read_data.0);

            drop(read_data);

            let (corrected, umi_sorted_reads) =
                crate::translator::correct_umis(&mut this_bc_reads, chemistry_def);
            let mut reads = umi_sorted_reads.reads;
            let mut quals = umi_sorted_reads.quals;
            let umi_id = umi_sorted_reads.umi_id;
            let uu = umi_sorted_reads.unique_umis;
            let flags = umi_sorted_reads.flags;
            let readnames = umi_sorted_reads.readnames;
            barcode_data_this.nreads_umi_corrected = corrected;

            // Create assemblies.
            let ProcessBarcodeResult {
                conx,
                conxq,
                cids,
                cumi,
                junction_support,
                productive,
                nedges,
                umi_info,
            } = process_barcode(
                &chbid,
                single_end,
                is_tcr,
                is_bcr,
                is_gd,
                inner_primersx,
                outer_primersx,
                &inner_primer_exts,
                &outer_primer_exts,
                &mut reads,
                &mut quals,
                &umi_id,
                &uu,
                n50_n50_rpu as i32,
                refdata,
                &refdata_full,
                &rkmers_plus_full_20,
                &refdatax,
                args.min_contig_length,
                &mut barcode_data_this,
                &mut log2,
                log_opts,
                &heur,
            )?;

            // Write contigs and annotation for them.
            for i in 0..conx.len() {
                let tigname = format!("{barcode}_contig_{}", i + 1);
                // Build annotation.

                let mut can = ContigAnnotation::from_seq(
                    &conx[i],
                    &conxq[i],
                    &tigname,
                    refdata,
                    cids[i].len(),
                    cumi[i].len(),
                    false, // determined in ASM_CALL_CELLS stage
                    false, // determined in ASM_CALL_CELLS stage
                    junction_support[i].clone(),
                );

                can.sample = Some(args.sample_id.clone());

                // confirm that productive labels assigned to contigs
                // match up with assembler con and con2 piles.
                // this is not true for denovo mode
                if !heur.free {
                    assert_eq!(can.productive.unwrap(), productive[i]);
                }

                // Output as json with four space indentation
                ann_writer.write_item(&can)?;
            }

            // Make a vector showing which umis/reads are assigned to which contigs.
            // By design, each UMI is assigned to at most one contig, and each
            // read is assigned to at most one contig.
            let mut contig_of_read = vec![None; reads.len()];
            let mut contig_of_umi = match umi_id.last() {
                Some(&id) => vec![None; id as usize + 1],
                None => Vec::new(),
            };
            for (t, cid) in cids.iter().enumerate() {
                for &read_id in cid {
                    contig_of_read[read_id as usize] = Some(t);
                    contig_of_umi[umi_id[read_id as usize] as usize] = Some(t);
                }
            }
            // Write umi_summary
            const MIN_READS_PER_UMI: i32 = 1;
            for info in umi_info {
                umi_summary_writer.write_item(&UmiSummaryRow {
                    barcode: barcode.clone(),
                    umi_id: info.umi_id,
                    umi: info.umi.to_string(),
                    reads: info.reads,
                    min_umi_reads: MIN_READS_PER_UMI,
                    good_umi: true,
                    contigs: info
                        .contig
                        .map_or_else(String::new, |c| c.tigname(&barcode)),
                })?;
            }

            // The read counts generated using UmiInfo do not match the original counts
            // see: CELLRANGER-9116
            // let mut contig_umi_info: Vec<UmiInfo> = umi_info
            //     .into_iter()
            //     .filter(|info| info.contig.is_some())
            //     .collect();
            // contig_umi_info.sort_by_key(|info| info.contig.clone());
            // for (contig_id, group) in &contig_umi_info
            //     .into_iter()
            //     .chunk_by(|info| info.contig.clone().unwrap())
            // {
            //     let group: Vec<UmiInfo> = group.into_iter().collect();
            //     let barcode = group.first().unwrap().barcode.to_string();
            //     summary_writer.write_item(&ContigSummaryRow {
            //         barcode: barcode.clone(),
            //         contig_name: contig_id.tigname(&barcode),
            //         num_reads: group.iter().map(|info| info.reads).sum(),
            //         num_pairs: group.iter().map(|info| info.read_pairs).sum(),
            //         num_umis: group.len(),
            //         umi_list: UmiList(group.iter().map(|info| info.umi_id).collect()),
            //     });
            // }

            // Write summary.
            for (t, cid) in cids.iter().enumerate() {
                let npairs = if single_end {
                    cid.len()
                } else {
                    // The following code is used to count the read
                    // pairs by checking whether both read1 and read2 are
                    // part of the reads in this contig. In paired end mode,
                    // the even ids are read1 and the adjacent odd ids are read 2.
                    cid.iter()
                        .tuple_windows()
                        .filter(|&(a, b)| *a % 2 == 0 && *b == *a + 1)
                        .count()
                };

                summary_writer.write_item(&ContigSummaryRow {
                    barcode: barcode.clone(),
                    contig_name: format!("{barcode}_contig_{}", t + 1),
                    num_reads: cid.len(),
                    num_pairs: npairs,
                    num_umis: cumi[t].len(),
                    umi_list: UmiList(cumi[t].clone().into_iter().collect()),
                })?;
            }
            // Define scoring for alignments of reads to contigs.

            let scoring =
                sw::Scoring::from_scores(-GAP_OPEN, -GAP_EXTEND, MATCH_SCORE, -MISMATCH_SCORE)
                    .xclip(-CLIP)
                    .yclip(0);
            let min_align_score = 50.0;

            // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

            // Create a vector of Read objects.  We only fill in the entries
            // that are used below.
            // ◼ This is pretty horrible, as we're
            // ◼ just shoveling data into a different form.

            let mut readsx = Vec::<graph_read::Read>::new();
            for i in 0..reads.len() {
                let mut r = graph_read::Read::new(
                    i as ReadType,        // read id
                    0,                    // umi, "leaving blank"
                    readnames[i].clone(), // read name
                    reads[i].clone(),     // read sequence
                    quals[i].clone(),     // quality score for read
                );
                r.flags = flags[i]; // flags
                readsx.push(r);
            }

            // Add to the sam file.  First gather contig names.
            // ◼ The creation of bam headers here is gratuituous, since what
            // ◼ we put in there just gets regurgitated as sam records.  Since
            // ◼ the rust_htslib code is flaky, we should get rid of this
            // ◼ pointless conversion.

            let mut tignames = Vec::<String>::new();
            for (i, contig) in conx.iter().enumerate() {
                let contig_name = format!("{barcode}_contig_{}", i + 1);
                tignames.push(contig_name.clone());
                bam_utils::add_ref_to_bam_header(&mut sam_header, &contig_name, contig.len());
            }

            // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
            // Create a Vec<Vec<u8>> from Vec<DnaString>
            let contig_seqs: Vec<_> = conx.iter().map(DnaString::to_ascii_vec).collect();
            let contig_aligners: Vec<_> = contig_seqs
                .iter()
                .map(|seq| {
                    ContigAligner::new(
                        seq,
                        scoring,
                        KMER_LEN_BANDED_ALIGN,
                        WINDOW_SIZE_BANDED_ALIGN,
                    )
                })
                .collect();
            let mut i = 0;
            while i < reads.len() {
                let read = reads[i].to_ascii_vec();
                let u = &uu[umi_id[i] as usize];
                let aln_packet = contig_of_read[i].and_then(|contig_id| {
                    contig_aligners[contig_id]
                        .align_read(&read, min_align_score as i32)
                        .map(|al| sw::AlignmentPacket {
                            ref_idx: contig_id,
                            alignment: al,
                        })
                });
                if single_end {
                    let rec = bam_utils::read_to_bam_record_opts(
                        &readsx[i],
                        &aln_packet,
                        &None,
                        true,
                        true,
                    );
                    write_sam_record_simple(&rec, u, &tignames, &mut simple_sam_writer2);
                    i += 1;
                } else {
                    let mate_read = reads[i + 1].to_ascii_vec();
                    let mate_aln_packet = contig_of_read[i + 1].and_then(|contig_id| {
                        contig_aligners[contig_id]
                            .align_read(&mate_read, min_align_score as i32)
                            .map(|al| sw::AlignmentPacket {
                                ref_idx: contig_id,
                                alignment: al,
                            })
                    });

                    let rec = bam_utils::read_to_bam_record_opts(
                        &readsx[i],
                        &aln_packet,
                        &mate_aln_packet,
                        true,
                        true,
                    );
                    let mate_rec = bam_utils::read_to_bam_record_opts(
                        &readsx[i + 1],
                        &mate_aln_packet,
                        &aln_packet,
                        true,
                        true,
                    );
                    write_sam_record_simple(&rec, u, &tignames, &mut simple_sam_writer2);
                    write_sam_record_simple(&mate_rec, u, &tignames, &mut simple_sam_writer2);

                    i += 2;
                }
            }
            // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

            // Finish up.

            bid += 1;
            barcode_data.push(barcode_data_this);

            if nedges > 0 {
                log.append(&mut log2);
            }
        }
    }

    // Print log.
    // ◼ Figure out if we want to print this and if so to _stdout or elsewhere.

    print!("{}", String::from_utf8(log).unwrap());

    // Finish writing of the simple_sam_writer.  It is not enough to have it go
    // out of scope, which seems like a bug in lz4.
    // Then end scope.  See commments at beginning of scope.

    let (_, res) = simple_sam_writer.finish();
    res?;
    ann_writer.finish()?;
    umi_summary_writer.finish()?;
    summary_writer.finish()?;
    let umi_info = Vec::<UmiInfo>::new();
    Ok((sam_header, barcodes, barcode_data, umi_info))
}

fn write_contig_summary_tsv(
    summary: &TsvFile<ContigSummaryRow>,
    chunk_outs: &[AssemblyChunkOutputs],
) -> Result<()> {
    let mut writer = summary.lazy_writer()?;
    let mut has_contents = false;
    for co in chunk_outs {
        let summary = co.summary_tsv.read()?;
        has_contents |= !summary.is_empty();
        for row in summary {
            writer.write_item(&row)?;
        }
    }
    if !has_contents {
        writer.write_header()?;
    }
    writer.finish()?;
    Ok(())
}

fn write_umi_summary_tsv(
    summary: &TsvFile<UmiSummaryRow>,
    chunk_outs: &[AssemblyChunkOutputs],
) -> Result<()> {
    let mut writer = summary.lazy_writer()?;
    let mut has_contents = false;
    for co in chunk_outs {
        let summary = co.umi_summary_tsv.read()?;
        has_contents |= !summary.is_empty();
        for row in summary {
            writer.write_item(&row)?;
        }
    }
    if !has_contents {
        writer.write_header()?;
    }
    writer.finish()?;
    Ok(())
}

fn write_bc_data_full(file: &BarcodeDataFile, chunk_outs: &[AssemblyChunkOutputs]) -> Result<()> {
    let mut writer = file.lazy_writer()?;
    for chunk in chunk_outs {
        let reader = chunk.barcode_data.lazy_reader()?;
        for bc_data in reader {
            let bc_data = bc_data?;
            writer.write_item(&bc_data)?;
        }
    }
    writer.finish()?;
    Ok(())
}

fn write_bc_databrief(
    file: &BarcodeDataBriefFile,
    chunk_outs: &[AssemblyChunkOutputs],
) -> Result<()> {
    let mut writer = file.lazy_writer()?;
    for chunk in chunk_outs {
        let reader = chunk.barcode_data.lazy_reader()?;
        for bc_data in reader {
            let bc_data = bc_data?;
            let bc_brief: BarcodeDataBrief = bc_data.into();
            writer.write_item(&bc_brief)?;
        }
    }
    writer.finish()?;
    Ok(())
}

fn write_umi_info_full(file: &UmiInfoFile, chunk_outs: &[AssemblyChunkOutputs]) -> Result<()> {
    let mut writer = file.lazy_writer()?;
    for chunk in chunk_outs {
        let reader = chunk.umi_info.lazy_reader()?;
        for umi_info in reader {
            let umi_info = umi_info?;
            writer.write_item(&umi_info)?;
        }
    }
    writer.finish()?;
    Ok(())
}

/// Read the regions.fa file from the provided reference path.
fn read_regions(ref_path: &Path) -> Result<String> {
    let fasta_path = ref_path.join("fasta/regions.fa");
    let fasta =
        read_to_string(&fasta_path).with_context(|| fasta_path.to_string_lossy().to_string())?;
    assert!(
        !fasta.is_empty(),
        "Reference file at {} has zero length.",
        fasta_path.display(),
    );
    Ok(fasta)
}
