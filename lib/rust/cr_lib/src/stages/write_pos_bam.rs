//! Martian stage WRITE_POS_BAM

use crate::{env, AlignShardFile, BamFile};
use anyhow::{bail, Result};
use barcode::Barcode;
use cr_bam::bam::{estimate_memory_for_range_gib, BamPosSort};
use cr_types::rna_read::{make_library_info, RnaChunk};
use cr_types::types::{SampleAssignment, SampleBarcodes};
use cr_types::{BarcodeToSample, SampleBarcodesFile};
use itertools::Itertools;
use martian::prelude::*;
use martian_derive::{make_mro, martian_filetype, MartianStruct};
use metric::TxHashMap;
use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::{self, Header, Record};
use serde::{Deserialize, Serialize};
use shardio::ShardReader;
use std::cmp::max;
use std::fs::File;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

fn use_csi_instead_of_bai(header: Header) -> bool {
    // The standard BAM index file can only handle chromosomes up to 512 MB in length,
    // this function examines the BAM header and decides if we'll need to use a CSI
    // index file instead by finding the largest reference chromosome

    // Function to get the largest sequence observed.
    let map = header.to_hashmap();
    let seqs_op = map.get("SQ");
    let mut max_sq_len: u64 = 0;
    if let Some(seqs) = seqs_op {
        for seq in seqs {
            if let Some(ln) = seq.get("LN") {
                if let Ok(length) = ln.parse::<u64>() {
                    max_sq_len = max_sq_len.max(length);
                }
            }
        }
    }
    // BAI files can only take up to 512MB as a max index
    let largest_allowed: u64 = 1 << 29;
    max_sq_len > largest_allowed
}

/// Adapted from cellranger-dna
pub fn concat_bams(paths: &[BamFile], out: &Path, threads: usize) -> Result<()> {
    let paths = paths.iter().collect::<Vec<_>>();
    // if there's only a single input path, copy the input
    if paths.len() == 1 {
        std::fs::hard_link(paths[0], out).or_else(|_| std::fs::copy(paths[0], out).map(|_| ()))?;
        return Ok(());
    }
    let fofn = out.with_extension("fofn");
    let fofn = fofn.to_str().unwrap();
    {
        let mut handle = File::create(fofn)?;
        for path in paths {
            handle.write_all(path.as_ref().to_str().unwrap().as_bytes())?;
            handle.write_all(b"\n")?;
        }
    }
    let status = Command::new("samtools")
        .args([
            "cat",
            "--threads",
            &threads.to_string(),
            "-b",
            fofn,
            "-o",
            out.to_str().unwrap(),
        ])
        .spawn()?
        .wait()?;
    if status.success() {
        std::fs::remove_file(fofn)?;
        Ok(())
    } else {
        panic!("failed to `samtools cat {fofn}`");
    }
}

martian_filetype! { BaiIndexFile, "bam.bai" }
martian_filetype! { CsiIndexFile, "bam.csi" }

/// Adapted from assemble_vdj
pub fn index_bam(bam: &Path, bam_index: &Path, threads: usize, use_csi: bool) {
    let index_type = if use_csi { "-c" } else { "-b" };
    Command::new("samtools")
        .args([
            "index",
            index_type,
            "-@",
            &threads.to_string(),
            bam.to_str().unwrap(),
            bam_index.to_str().unwrap(),
        ])
        .output()
        .unwrap_or_else(|_| panic!("failed to `samtools index {}`", bam.display()));
}

const WRITE_BAM_THREADS: isize = 2;

/// Output the sorted and indexed BAM file `possorted_genome_bam.bam`.
pub struct WritePosBam;

#[derive(Deserialize, Clone, MartianStruct)]
pub struct StageInputs {
    pub bam_header: PathBuf,
    pub alignments: Vec<AlignShardFile>,
    pub read_chunks: Vec<RnaChunk>,
    pub target_set_name: Option<String>,
    pub sample_barcodes: Option<SampleBarcodesFile>,
    pub slide_serial_capture_area: Option<String>,
}

#[derive(Serialize, Deserialize, Clone, MartianStruct)]
pub struct ChunkInputs {
    #[mro_type = "map"]
    pub range: shardio::Range<<BamPosSort as shardio::SortKey<Record>>::Key>,
    pub write_header: bool,
}

#[derive(Serialize, Deserialize, Clone, MartianStruct)]
pub struct ChunkOutputs {
    // if there is only one BamFile, that implies single-sample mode
    // if there are two or more BamFiles, that implies multiplexed
    pub sample_pos_sorted_bam_chunks: TxHashMap<SampleAssignment, BamFile>,
}

#[derive(Serialize, Deserialize, Clone, MartianStruct)]
pub struct SampleBamFile {
    pub sample: SampleAssignment,
    pub bam_file: BamFile,
    pub bai_index_file: Option<BaiIndexFile>,
    pub csi_index_file: Option<CsiIndexFile>,
}

#[derive(Serialize, Clone, MartianStruct)]
pub struct StageOutputs {
    pub pos_sorted_bam: Option<SampleBamFile>,
    pub multi_pos_sorted_bam: Option<Vec<SampleBamFile>>,
}

fn attach_read_group_tags<'a>(
    header: &mut bam::Header,
    read_chunks: impl IntoIterator<Item = &'a RnaChunk>,
) {
    for read_group in read_chunks.into_iter().map(|x| &x.read_group).unique() {
        // sample_id, library_id, gem_group, flowcell, lane
        let parts = read_group.split(':').collect::<Vec<_>>();
        let mut record = HeaderRecord::new(b"RG");
        // ID: {read_group}, SM: {sample_id}, LB: {library_id}.{gem_group}, PU: {read_group}
        record.push_tag(b"ID", read_group);
        record.push_tag(b"SM", parts[0]);
        record.push_tag(b"LB", &format!("{}.{}", parts[1], parts[2]));
        record.push_tag(b"PU", read_group);
        record.push_tag(b"PL", "ILLUMINA");
        header.push_record(&record);
    }
}

/// Add BAM header comments that are used by bamtofastq.
fn attach_bamtofastq_comments(header: &mut bam::Header, read_chunks: &[RnaChunk]) {
    let bamtofastq_header_comments = read_chunks
        .iter()
        .map(|x| x.chemistry.name.bamtofastq_headers())
        .unique()
        .exactly_one()
        .unwrap();
    for comment in bamtofastq_header_comments {
        header.push_comment(comment.as_bytes());
    }
}

fn attach_library_info_comments(
    header: &mut bam::Header,
    read_chunks: &[RnaChunk],
    target_set_name: Option<&str>,
) {
    let library_infos = make_library_info(read_chunks, target_set_name);
    for info in library_infos {
        let comment = format!(r#"library_info:{}"#, serde_json::to_string(&info).unwrap());
        header.push_comment(comment.as_bytes());
    }
}

fn attach_slide_serial_capture_area(
    header: &mut bam::Header,
    slide_serial_capture_area: Option<String>,
) {
    if let Some(capture_area) = slide_serial_capture_area {
        let comment = format!(r#"slide_serial_capture_area:{capture_area}"#);
        header.push_comment(comment.as_bytes());
    }
}

pub fn make_header(
    bam_header: Option<&Path>,
    read_chunks: &[RnaChunk],
    target_set_name: Option<&str>,
    write_header: bool,
    slide_serial_capture_area: Option<String>,
    software_version: &str,
) -> Result<Header> {
    let mut header = Vec::new();
    if let Some(bam_header) = bam_header {
        writeln!(&mut header, "@HD\tVN:1.4\tSO:coordinate")?;
        File::open(bam_header)?.read_to_end(&mut header)?;
    } else {
        writeln!(&mut header, "@HD\tVN:1.4\tSO:unsorted")?;
    }
    let product_name = env::get_tenx_product_name();
    let product_id = env::get_tenx_product_id();
    writeln!(
        &mut header,
        "@PG\tID:{product_id}\tPN:{product_name}\tVN:{software_version}"
    )?;

    let header_view = bam::HeaderView::from_bytes(&header);
    let mut header = bam::Header::from_template(&header_view);
    // only the first of these is used in the `samtools cat` later
    if write_header {
        attach_read_group_tags(&mut header, read_chunks);
        attach_bamtofastq_comments(&mut header, read_chunks);
        attach_library_info_comments(&mut header, read_chunks, target_set_name);
        attach_slide_serial_capture_area(&mut header, slide_serial_capture_area);
    }
    Ok(header)
}

#[make_mro(volatile = strict)]
impl MartianStage for WritePosBam {
    type StageInputs = StageInputs;
    type StageOutputs = StageOutputs;
    type ChunkInputs = ChunkInputs;
    type ChunkOutputs = ChunkOutputs;

    fn split(
        &self,
        args: Self::StageInputs,
        _rover: MartianRover,
    ) -> Result<StageDef<Self::ChunkInputs>> {
        let reader = ShardReader::<Record, BamPosSort>::open_set(&args.alignments)?;
        // TODO: no particular reason I use this heuristic
        let nchunks = args.alignments.len().min(150);
        Ok(reader
            .make_chunks(nchunks, &shardio::Range::all())
            .into_iter()
            .enumerate()
            .map(|(chunk_id, range)| {
                let additional_mem_gib = if args.sample_barcodes.is_some() {
                    1.5
                } else {
                    1.0
                };
                let mem_gb = (additional_mem_gib + estimate_memory_for_range_gib(&range, &reader))
                    .ceil() as isize;
                let write_header = chunk_id == 0;
                (
                    ChunkInputs {
                        write_header,
                        range,
                    },
                    Resource::with_mem_gb(mem_gb).threads(WRITE_BAM_THREADS),
                )
            })
            .collect::<StageDef<_>>()
            .join_resource(Resource::with_mem_gb(1).threads(4)))
    }

    fn main(
        &self,
        args: Self::StageInputs,
        chunk_args: Self::ChunkInputs,
        rover: MartianRover,
    ) -> Result<Self::ChunkOutputs> {
        let sample_barcodes = SampleBarcodes::read_from_json(args.sample_barcodes.as_ref())?;
        let barcode_to_sample = BarcodeToSample::construct(&sample_barcodes);

        // Read in the header text & initialize the BAM file with it.
        let header = make_header(
            Some(&args.bam_header),
            &args.read_chunks,
            args.target_set_name.as_deref(),
            chunk_args.write_header,
            args.slide_serial_capture_area,
            &rover.pipelines_version(),
        )?;

        let threads: u32 = rover.get_threads() as u32;
        let threadpool = rust_htslib::tpool::ThreadPool::new(max(1, threads - 1))?;
        let mut writers: TxHashMap<SampleAssignment, bam::Writer> = TxHashMap::default();
        let mut bam_paths: TxHashMap<SampleAssignment, BamFile> = TxHashMap::default();
        for sample in sample_barcodes.get_samples_with_unassigned() {
            let bam_path: BamFile = rover.make_path(format!("pos_sorted{sample}"));
            let mut writer = bam::Writer::from_path(&bam_path, &header, bam::Format::Bam)?;
            writer.set_thread_pool(&threadpool)?;
            writers.insert(sample.clone(), writer);
            bam_paths.insert(sample.clone(), bam_path);
        }

        let reader = ShardReader::<Record, BamPosSort>::open_set(&args.alignments)?;
        for r in reader.iter_range(&chunk_args.range)? {
            let rec = &r?;

            let sample = if sample_barcodes.is_multiplexed() {
                let barcode: Option<Barcode> = match rec.aux(b"CB:Z") {
                    Ok(aux_data) => match aux_data {
                        Aux::String(aux_bc) => Some(Barcode::from_bytes(aux_bc.as_bytes())?),
                        _ => {
                            bail!("CB tag with non-string data!");
                        }
                    },
                    Err(_) => None,
                };

                if let Some(barcode) = barcode {
                    barcode_to_sample.get_sample(&barcode)
                } else {
                    &SampleAssignment::Unassigned
                }
            } else {
                &SampleAssignment::NonMultiplexed
            };
            writers.get_mut(sample).unwrap().write(rec)?;
        }

        Ok(ChunkOutputs {
            sample_pos_sorted_bam_chunks: bam_paths,
        })
    }

    fn join(
        &self,
        args: Self::StageInputs,
        _chunk_defs: Vec<Self::ChunkInputs>,
        chunk_outs: Vec<Self::ChunkOutputs>,
        rover: MartianRover,
    ) -> Result<Self::StageOutputs> {
        //Determine the time of BAM index needed
        let header: Header = make_header(
            Some(&args.bam_header),
            &args.read_chunks,
            args.target_set_name.as_deref(),
            false,
            args.slide_serial_capture_area,
            &rover.pipelines_version(),
        )?;
        let use_csi = use_csi_instead_of_bai(header);

        let samples_with_unassigned =
            SampleBarcodes::read_samples_with_unassigned(args.sample_barcodes.as_ref())?;
        let mut multi_bam_paths = Vec::with_capacity(samples_with_unassigned.len());
        for sample in samples_with_unassigned {
            let bam_path: BamFile = rover.make_path(format!("pos_sorted_{sample}"));
            let bams: Vec<BamFile> = chunk_outs
                .iter()
                .map(|x| {
                    x.sample_pos_sorted_bam_chunks
                        .get(&sample)
                        .expect("BAM writing chunk output was missing a sample id in hashmap.")
                        .clone()
                })
                .collect();
            let threads = rover.get_threads();

            concat_bams(&bams, &bam_path, (threads - 1).max(1))?;
            let index_prefix = format!("pos_sorted_{sample}");
            if use_csi {
                let index_path: CsiIndexFile = rover.make_path(index_prefix);
                index_bam(&bam_path, &index_path, (threads - 1).max(1), use_csi);
                multi_bam_paths.push(SampleBamFile {
                    sample: sample.clone(),
                    bam_file: bam_path,
                    bai_index_file: None,
                    csi_index_file: Some(index_path),
                });
            } else {
                let index_path: BaiIndexFile = rover.make_path(index_prefix);
                index_bam(&bam_path, &index_path, (threads - 1).max(1), use_csi);
                multi_bam_paths.push(SampleBamFile {
                    sample: sample.clone(),
                    bam_file: bam_path,
                    bai_index_file: Some(index_path),
                    csi_index_file: None,
                });
            }
        }

        Ok(if args.sample_barcodes.is_some() {
            // multiplexed output (has at least 1 sample BAM, and unassigned BAM)
            StageOutputs {
                pos_sorted_bam: None,
                multi_pos_sorted_bam: Some(multi_bam_paths),
            }
        } else {
            // single-sample non-multiplexed output (single BAM for all reads) StageOutputs {
            StageOutputs {
                pos_sorted_bam: Some(multi_bam_paths.into_iter().next().unwrap()),
                multi_pos_sorted_bam: None,
            }
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testing::correctness::check_bam_file_correctness;
    use cr_types::chemistry::{ChemistryDef, ChemistryName};
    use cr_types::types::SampleAssignment::{Assigned, Unassigned};
    use cr_types::LibraryType;
    use fastq_set::read_pair_iter::InputFastqs;

    fn get_test_chunk() -> RnaChunk {
        use cr_types::sample_def::SampleDef;
        RnaChunk::new(
            &ChemistryDef::named(ChemistryName::FivePrimePE),
            &SampleDef::default(),
            LibraryType::Gex,
            InputFastqs {
                r1: String::default(),
                r2: None,
                i1: None,
                i2: None,
                r1_interleaved: false,
            },
            "",
            0,
            0,
        )
    }

    #[test]
    #[should_panic]
    fn test_fail_bam_comparison() {
        let chunk = get_test_chunk();
        let args = StageInputs {
            bam_header: PathBuf::from("test/multi/write_bam_multi/test1/bam_header"),
            alignments: vec![AlignShardFile::from(
                "test/multi/write_bam_multi/test1/pos_sorted.asf",
            )],
            read_chunks: vec![
                RnaChunk {
                    read_group: "dev_pbmc1k_5p_micro_42172:0:1:H2NV3DMXX:1".to_string(),
                    ..chunk.clone()
                },
                RnaChunk {
                    read_group: "dev_pbmc1k_5p_micro_42172:0:1:H2NV3DMXX:2".to_string(),
                    ..chunk
                },
            ],
            target_set_name: None,
            sample_barcodes: None,
            slide_serial_capture_area: None,
        };

        let outs = WritePosBam.test_run_tmpdir(args).unwrap();
        check_bam_file_correctness(
            &outs.pos_sorted_bam.unwrap().bam_file,
            Path::new("test/multi/write_bam_multi/test1/wrong_bam.bam"),
        )
        .unwrap();
    }

    #[test]
    fn test_single_sample() {
        let chunk = get_test_chunk();
        let args = StageInputs {
            bam_header: PathBuf::from("test/multi/write_bam_multi/test1/bam_header"),
            alignments: vec![AlignShardFile::from(
                "test/multi/write_bam_multi/test1/pos_sorted.asf",
            )],
            read_chunks: vec![
                RnaChunk {
                    read_group: "dev_pbmc1k_5p_micro_42172:0:1:H2NV3DMXX:1".to_string(),
                    ..chunk.clone()
                },
                RnaChunk {
                    read_group: "dev_pbmc1k_5p_micro_42172:0:1:H2NV3DMXX:2".to_string(),
                    ..chunk
                },
            ],
            target_set_name: None,
            sample_barcodes: None,
            slide_serial_capture_area: None,
        };

        let tmp_dir = tempfile::tempdir().unwrap();
        let outs = WritePosBam.test_run(&tmp_dir, args).unwrap();
        check_bam_file_correctness(
            &outs.pos_sorted_bam.unwrap().bam_file,
            Path::new("test/multi/write_bam_multi/test1/expected_bam.bam"),
        )
        .unwrap();
    }

    #[test]
    fn test_multiple_sample() {
        let chunk = get_test_chunk();

        let args = StageInputs {
            bam_header: PathBuf::from("test/multi/write_bam_multi/test1/bam_header"),
            alignments: vec![AlignShardFile::from(
                "test/multi/write_bam_multi/test1/pos_sorted.asf",
            )],
            read_chunks: vec![
                RnaChunk {
                    read_group: "dev_pbmc1k_5p_micro_42172:0:1:H2NV3DMXX:1".to_string(),
                    ..chunk.clone()
                },
                RnaChunk {
                    read_group: "dev_pbmc1k_5p_micro_42172:0:1:H2NV3DMXX:2".to_string(),
                    ..chunk
                },
            ],
            target_set_name: None,
            sample_barcodes: Some("test/multi/write_bam_multi/test1/sample_barcodes.json".into()),
            slide_serial_capture_area: None,
        };

        let tmp_dir = tempfile::tempdir().unwrap();
        let outs = WritePosBam.test_run(&tmp_dir, args).unwrap();
        check_bam_file_correctness(
            &outs
                .multi_pos_sorted_bam
                .as_ref()
                .unwrap()
                .iter()
                .find(|x| x.sample == Assigned("test_sample1".to_string()))
                .unwrap()
                .bam_file,
            Path::new("test/multi/write_bam_multi/test1/expected_sample1.bam"),
        )
        .unwrap();

        check_bam_file_correctness(
            &outs
                .multi_pos_sorted_bam
                .as_ref()
                .unwrap()
                .iter()
                .find(|x| x.sample == Assigned("test_sample2".to_string()))
                .unwrap()
                .bam_file,
            Path::new("test/multi/write_bam_multi/test1/expected_sample2.bam"),
        )
        .unwrap();

        check_bam_file_correctness(
            &outs
                .multi_pos_sorted_bam
                .unwrap()
                .iter()
                .find(|x| x.sample == Unassigned)
                .unwrap()
                .bam_file,
            Path::new("test/multi/write_bam_multi/test1/expected_unassigned.bam"),
        )
        .unwrap();
    }
}
