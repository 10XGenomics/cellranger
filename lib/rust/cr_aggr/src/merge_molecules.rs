//! Martian stage MERGE_MOLECULES
#![expect(missing_docs)]

use anyhow::Result;
use barcode::BarcodeContent;
use cr_h5::molecule_info::{
    BarcodeIdxType, LibraryIdxType, MoleculeInfoIterator, MoleculeInfoReader, MoleculeInfoWriter,
    bin_barcodes, open_for_modification, trim_barcodes,
};
use cr_types::aggr::SampleDef;
use cr_types::probe_set::Probe;
use cr_types::reference::feature_reference::FeatureReference;
use cr_types::{LibraryInfo, aggr};
use itertools::Itertools;
use martian::prelude::*;
use martian_derive::{MartianStruct, make_mro, martian_filetype};
use metric::{JsonReporter, TxHashSet};
use ndarray::s;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::collections::HashMap;
use std::iter::zip;
use std::path::PathBuf;

martian_filetype!(H5File, "h5");

#[derive(Debug, Deserialize, Clone, MartianStruct)]
pub struct MergeMoleculesStageInputs {
    sample_defs: Vec<SampleDef>,
    libraries: Vec<aggr::LibraryInfo>,
}

#[derive(Debug, Serialize, Deserialize, Clone, MartianStruct)]
pub struct MergeMoleculesChunkInputs {
    sample_def: SampleDef,
}

#[derive(Debug, Serialize, Deserialize, MartianStruct)]
pub struct MergeMoleculesChunkOutputs {
    sample_def: SampleDef,
}

#[derive(Debug, Serialize, Deserialize, MartianStruct)]
pub struct MergeMoleculesStageOutputs {
    merged_molecules: H5File,
    #[mro_type = "map<int[]>"]
    gem_group_barcode_ranges: HashMap<String, Vec<u64>>,
}

// NOTE: Sync with switchboard.py
const HD_BINNING_SCALE: u32 = 8;

pub struct MergeMolecules;

#[make_mro(volatile = strict)]
impl MartianStage for MergeMolecules {
    type StageInputs = MergeMoleculesStageInputs;
    type StageOutputs = MergeMoleculesStageOutputs;
    type ChunkInputs = MergeMoleculesChunkInputs;
    type ChunkOutputs = MergeMoleculesChunkOutputs;

    fn split(
        &self,
        args: Self::StageInputs,
        _rover: MartianRover,
    ) -> Result<StageDef<Self::ChunkInputs>> {
        /// Return the memory required to load n barcodes.
        fn barcode_mem_gib(num_barcodes: usize) -> usize {
            (num_barcodes * size_of::<BarcodeContent>()).div_ceil(1024 * 1024 * 1024)
        }

        // The number of non-zero barcodes per chunk.
        let num_nz_barcodes: Vec<_> = args
            .sample_defs
            .iter()
            .map(|sample_def| MoleculeInfoReader::count_non_zero_bcs(&sample_def.molecule_h5))
            .try_collect()?;
        let num_barcodes: Vec<_> = args
            .sample_defs
            .iter()
            .map(|sample_def| MoleculeInfoReader::read_barcodes_size(&sample_def.molecule_h5))
            .try_collect()?;
        let chunk_mem_gib: Vec<_> = num_barcodes
            .iter()
            .map(|&n| 1 + barcode_mem_gib(n))
            .collect();

        // The maximum number of barcodes is the sum of the non-zero barcodes in each chunk.
        let num_nz_barcodes_sum: usize = num_nz_barcodes.iter().sum();
        let join_mem_gib = 5 + barcode_mem_gib(num_nz_barcodes_sum);

        for (i, (n, mem_gib)) in zip(num_barcodes, &chunk_mem_gib).enumerate() {
            println!("chunk={i},num_barcodes={n},mem_gib={mem_gib}");
        }
        println!("chunk=join,num_barcodes={num_nz_barcodes_sum},mem_gib={join_mem_gib}");

        Ok(zip(args.sample_defs, chunk_mem_gib)
            .map(|(sample_def, mem_gib)| {
                (
                    MergeMoleculesChunkInputs { sample_def },
                    Resource::with_mem_gb(mem_gib as isize),
                )
            })
            .collect::<StageDef<_>>()
            .join_resource(Resource::with_mem_gb(join_mem_gib as isize)))
    }

    fn main(
        &self,
        _args: Self::StageInputs,
        chunk_args: Self::ChunkInputs,
        rover: MartianRover,
    ) -> Result<Self::ChunkOutputs> {
        let molecule_h5: PathBuf = rover.make_path(format!(
            "{}_molecule_info.h5",
            chunk_args.sample_def.library_id
        ));
        std::fs::copy(chunk_args.sample_def.molecule_h5, &molecule_h5)?;

        // Ensure we can write to the file and that it is a compatible version.
        open_for_modification(&molecule_h5)?;

        if MoleculeInfoReader::is_visium_hd(&molecule_h5)? {
            let tmp_shard_path: PathBuf = rover.make_path("tmp_shard.shard");
            bin_barcodes(&molecule_h5, HD_BINNING_SCALE, &tmp_shard_path)?;
        }
        trim_barcodes(&molecule_h5, false)?;

        Ok(MergeMoleculesChunkOutputs {
            sample_def: SampleDef {
                molecule_h5,
                ..chunk_args.sample_def
            },
        })
    }

    fn join(
        &self,
        args: Self::StageInputs,
        _chunk_defs: Vec<Self::ChunkInputs>,
        chunk_outs: Vec<Self::ChunkOutputs>,
        rover: MartianRover,
    ) -> Result<Self::StageOutputs> {
        let new_sample_defs: Vec<_> = chunk_outs.into_iter().map(|x| x.sample_def).collect();

        // TODO: merged_barcodes can have duplicated barcodes,
        // since new barcode_idx is calculated by updating an offset
        // into the concatenated trimmed barcodes from previous step.
        // But this can be collected into a BTreeSet, sorted,
        // and reuse idx. Maybe? Need to check.
        let merged_barcodes: Vec<_> = new_sample_defs
            .iter()
            .map(|x| MoleculeInfoReader::read_barcodes(&x.molecule_h5))
            .flatten_ok()
            .map(|x| x?)
            .try_collect()?;

        let (feature_ref, probes) = get_fref_and_probes(&new_sample_defs)?;

        let mut gem_group_barcode_ranges: HashMap<String, Vec<u64>> = Default::default();
        let mut barcode_idx_offset: usize = 0;

        let mut total_rows = 0;

        struct Lib {
            lib: aggr::LibraryInfo,
            index: usize,
        }
        let libraries: Vec<Lib> = args
            .libraries
            .into_iter()
            .enumerate()
            .map(|(index, lib)| Lib { lib, index })
            .collect();

        // TODO: can probably do all this on the fly, no need to collect and
        // process later like the Python version because we have access to
        // everything already
        let mut in_h5s = Vec::with_capacity(new_sample_defs.len());
        let mut gg_maps = Vec::with_capacity(new_sample_defs.len());
        let mut lib_idx_maps: Vec<Vec<LibraryIdxType>> = Vec::with_capacity(new_sample_defs.len());
        let mut metrics = Vec::with_capacity(new_sample_defs.len());
        let mut bc_idx_offsets: Vec<BarcodeIdxType> = Vec::with_capacity(new_sample_defs.len());
        let mut bc_infos = Vec::with_capacity(new_sample_defs.len());

        for sample_def in new_sample_defs {
            let in_h5 = sample_def.molecule_h5;
            in_h5s.push(in_h5.clone());
            let aggr_id = sample_def.library_id;

            let in_mc_size = MoleculeInfoReader::read_barcodes_size(&in_h5)?;
            let barcode_idx_end = barcode_idx_offset + in_mc_size;

            // Get the gem group and library mappings
            let my_libs: Vec<_> = libraries
                .iter()
                .filter(|l| l.lib.aggr_id == aggr_id)
                .collect();
            let max_old_gg = my_libs
                .iter()
                .map(|l| l.lib.old_gem_group)
                .max()
                .unwrap_or(0) as usize;
            let max_old_lib_idx = my_libs
                .iter()
                .map(|l| l.lib.old_library_index)
                .max()
                .unwrap_or(0) as usize;

            let mut gg_map = vec![0; 1 + max_old_gg];
            let mut lib_idx_map: Vec<LibraryIdxType> = vec![0; 1 + max_old_lib_idx];

            for l in my_libs {
                gg_map[l.lib.old_gem_group as usize] = l.lib.gem_group;
                lib_idx_map[l.lib.old_library_index as usize] = l.index as u16;
            }

            bc_idx_offsets.push(barcode_idx_offset as BarcodeIdxType);

            let mut pass_filter = MoleculeInfoReader::read_barcode_info_pass_filter(&in_h5)?;
            // Offset the barcode index
            for x in pass_filter.slice_mut(s![.., 0]) {
                *x += barcode_idx_offset as u64;
            }
            // Update libraries to new mapping
            for x in pass_filter.slice_mut(s![.., 1]) {
                *x = lib_idx_map[*x as usize] as u64;
            }

            let genomes = MoleculeInfoReader::read_barcode_info_genomes(&in_h5)?;
            bc_infos.push((pass_filter, genomes));

            // Now get the metrics
            let reporter: JsonReporter =
                serde_json::from_str(&MoleculeInfoReader::read_metrics(&in_h5)?)?;

            // Remap the per-gem-group and per-library metrics
            let out_metrics: serde_json::Map<_, _> = reporter
                .into_iter()
                .map(|(k, v)| {
                    if k == "gem_groups" {
                        let new_gg_metrics: serde_json::Map<String, Value> = v
                            .as_object()
                            .unwrap()
                            .into_iter()
                            .map(|(og, m)| {
                                (gg_map[og.parse::<usize>().unwrap()].to_string(), m.clone())
                            })
                            .collect();
                        (k, Value::from(new_gg_metrics))
                    } else if k == "libraries" {
                        let new_lib_metrics: serde_json::Map<String, Value> = v
                            .as_object()
                            .unwrap()
                            .into_iter()
                            .map(|(og, m)| {
                                (
                                    lib_idx_map[og.parse::<LibraryIdxType>().unwrap() as usize]
                                        .to_string(),
                                    m.clone(),
                                )
                            })
                            .collect();
                        (k, Value::from(new_lib_metrics))
                    } else {
                        (k, v)
                    }
                })
                .collect();

            lib_idx_maps.push(lib_idx_map);
            metrics.push(out_metrics);

            let unique_gem_groups: TxHashSet<_> =
                MoleculeInfoReader::iter_gem_groups(&in_h5)?.collect::<Result<_>>()?;
            for x in unique_gem_groups {
                gem_group_barcode_ranges.insert(
                    gg_map[x as usize].to_string(),
                    vec![barcode_idx_offset as u64, barcode_idx_end as u64],
                );
            }

            gg_maps.push(gg_map);

            // Verify feature_refs match
            for (new, old) in zip(
                MoleculeInfoIterator::new(&in_h5)?.feature_ref.feature_defs,
                &feature_ref.feature_defs,
            ) {
                assert_eq!(new.id, old.id);
            }

            total_rows += MoleculeInfoReader::nrows(&in_h5)?;
            barcode_idx_offset = barcode_idx_end;
        }

        let combined_metrics = MoleculeInfoWriter::concatenate_metrics(metrics)?;
        let (pass_filter, genomes) = MoleculeInfoWriter::merge_barcode_infos(bc_infos);

        let merged_molecules: H5File = rover.make_path("merged_molecules.h5");

        let library_info: Vec<_> = libraries
            .into_iter()
            .map(|v| LibraryInfo::Aggr(v.lib))
            .collect();

        // Pass a placeholder vec of true for filtered probes.
        // TODO: it might be best to derive this from the inputs, but at present
        // nothing in the pipeline reads this data.
        let filtered_probes = probes.as_ref().map(|probes| vec![true; probes.len()]);

        let mut merged = MoleculeInfoWriter::new(
            &merged_molecules,
            &feature_ref,
            probes.as_deref(),
            filtered_probes.as_deref(),
            merged_barcodes,
            &library_info,
        )?;
        drop(feature_ref);
        drop(library_info);
        merged.write_metrics(&combined_metrics)?;
        drop(combined_metrics);
        merged.write_barcode_info(&pass_filter, &genomes)?;
        drop(pass_filter);
        drop(genomes);

        // concatenate all the columns
        merged.concatenate_many(in_h5s, bc_idx_offsets, gg_maps, lib_idx_maps)?;

        // Slow validation check here

        assert_eq!(
            total_rows,
            MoleculeInfoReader::nrows(&merged_molecules)?,
            "Concatenation did not produce expected results."
        );

        Ok(MergeMoleculesStageOutputs {
            merged_molecules,
            gem_group_barcode_ranges,
        })
    }
}

enum ProbesState {
    /// Initial state.
    Undetermined,
    /// We are using probes.
    Probes(Vec<Probe>),
    /// One or more inputs do not have probes.
    Skip,
}

impl ProbesState {
    /// Update our probes collection state with an optional collection.
    ///
    /// If we have more than one collection of probes, assert they are identical.
    fn update(self, probes: Option<Vec<Probe>>) -> Self {
        match (self, probes) {
            (Self::Undetermined, Some(p)) => Self::Probes(p),
            (Self::Probes(existing), Some(p)) => {
                assert_eq!(existing, p);
                Self::Probes(existing)
            }
            (_, None) | (Self::Skip, _) => Self::Skip,
        }
    }

    fn into_option(self) -> Option<Vec<Probe>> {
        match self {
            Self::Probes(p) => Some(p),
            _ => None,
        }
    }
}

/// Extract feature reference and probes data.
/// Prefer targeted feature reference if present.
/// Only include probes data if all inputs have probes.
fn get_fref_and_probes(
    new_sample_defs: &[SampleDef],
) -> Result<(FeatureReference, Option<Vec<Probe>>)> {
    let mut feature_ref = None;
    let mut probes = ProbesState::Undetermined;

    for sample_def in new_sample_defs {
        let reader = MoleculeInfoIterator::new(&sample_def.molecule_h5)?;
        probes = probes.update(reader.probes);
        let fref = reader.feature_ref;
        if feature_ref.is_none() || fref.has_target_features() {
            feature_ref = Some(fref);
        }
    }
    // The only way this unwrap can fail is if we have no samples...
    Ok((feature_ref.unwrap(), probes.into_option()))
}

#[cfg(test)]
mod merge_molecules_tests {
    use super::*;
    use cr_types::reference::feature_reference::FeatureReferenceFile;
    use cr_types::{LibraryType, UmiCount};
    use martian::MartianTempFile;
    use std::io::Write;
    use tempfile::TempDir;
    use umi::UmiType;

    const TEMPLATE_METRICS: &str = r#"{
"analysis_parameters": {"include_introns": false},
"cellranger_version": "2021.0428.1",
"chemistry_barcode": [{"kind": "gel_bead",
  "length": 16,
  "offset": 0,
  "read_type": "R1",
  "whitelist": "3M-february-2018"}],
"chemistry_description": "Single Cell 3' v3",
"chemistry_endedness": "three_prime",
"chemistry_name": "SC3Pv3",
"chemistry_rna": {"length": null,
  "min_length": 15,
  "offset": 0,
  "read_type": "R2"},
"chemistry_rna2": null,
"chemistry_strandedness": "+",
"chemistry_umi": {"length": 12,
 "min_length": 10,
 "offset": 16,
 "read_type": "R1"},
"gem_groups": {"0": {"force_cells": null, "recovered_cells": 2000}},
"libraries": {"0": {"feature_read_pairs": 443764534,
  "raw_read_pairs": 453355603,
  "usable_read_pairs": 239308782}},
"molecule_info_type": "count",
"reference_fasta_hash": "b6f131840f9f337e7b858c3d1e89d7ce0321b243",
"reference_gtf_hash": "d9aa710e1eab4e2bdb1a87c25d4cc8a9397db121",
"reference_gtf_hash.gz": "",
"reference_mkref_version": ""
}"#;

    #[test]
    fn test_merge_molecules_1() {
        let feature_def = {
            let fcsv = FeatureReferenceFile::tempfile().expect("Failed to create temp file");
            let mut writer = fcsv.buf_writer().unwrap();
            write!(
                &mut writer,
                "id,name,read,pattern,sequence,feature_type
                CMO301,Cell Multiplexing Oligo 301,R2,5P(BC),ATGAGGAATTCCTGC,Multiplexing Capture
                CMO302,Cell Multiplexing Oligo 302,R2,5P(BC),CATGCCAATAGAGCG,Multiplexing Capture
                CMO303,Cell Multiplexing Oligo 303,R2,5P(BC),CCGTCGTCCAAGCAT,Multiplexing Capture"
            )
            .expect("failed to write features CSV");
            writer.flush().expect("failed to flush file");
            fcsv.read(None).expect("error creating feature reference")
        };

        let barcode_lists = [
            [
                b"AAACCCAAGAAACACT",
                b"AAACCCAAGAAACCAT",
                b"AAACCCAAGAAACCCA",
            ]
            .map(|x| BarcodeContent::from_bytes(x).unwrap()),
            [
                b"AAACCCAAGAAACCCG",
                b"AAACCCAAGAAACGTC",
                b"AAACCCAAGAAACTTC",
            ]
            .map(|x| BarcodeContent::from_bytes(x).unwrap()),
        ];
        let tfiles: Vec<_> = barcode_lists
            .into_iter()
            .map(|barcodes| {
                let tfile = tempfile::NamedTempFile::new().expect("Failed to create temp file");

                let library_info = [LibraryInfo::Count(cr_types::rna_read::LibraryInfo {
                    library_id: 0,
                    gem_group: 0,
                    target_set_name: None,
                    library_type: LibraryType::Gex,
                })];

                let mut molinfo = MoleculeInfoWriter::new(
                    tfile.as_ref(),
                    &feature_def,
                    None,
                    None,
                    barcodes,
                    &library_info,
                )
                .expect("unable to create molecule file");

                // Note: only including the first barcode:
                //  - using the second for barcodes with count > 0 not in the pass_filter
                //  - skipping the third on purpose (no pass_filter, no counts)
                let pass_filter = ndarray::arr2(&[[0u64, 0, 0]]);

                molinfo
                    .write_barcode_info(&pass_filter, &["GRCh38".into()])
                    .expect("error writing barcode info");

                molinfo
                    .write_metrics(TEMPLATE_METRICS)
                    .expect("error writing metrics");

                for j in 0..10 {
                    // There are 3 barcodes per file, but the third one didn't pass_filter
                    molinfo
                        .write(
                            0,
                            j % 2_u64,
                            &UmiCount {
                                library_idx: 0,
                                feature_idx: (j % 3) as _,
                                probe_idx: None, //TODO: no aggring of probe_idx for now.
                                umi: j as _,
                                read_count: (j * 2) as _,
                                utype: UmiType::Txomic,
                            },
                        )
                        .expect("error adding umi to molecule_info");
                }
                molinfo.flush().expect("error flushing molinfo to disk");
                tfile
            })
            .collect();

        let libraries: Vec<_> = tfiles
            .iter()
            .enumerate()
            .map(|(i, _)| aggr::LibraryInfo {
                library_id: i as u16,
                library_type: Some("Gene Expression".into()),
                gem_group: i as u16,
                target_set_name: None,
                old_gem_group: 0,
                old_library_index: 0,
                aggr_id: format!("mi_{i}"),
                batch_name: i.to_string(),
                batch_id: i as u16,
            })
            .collect();

        let sample_defs: Vec<_> = tfiles
            .iter()
            .enumerate()
            .map(|(i, tfile)| SampleDef {
                library_id: format!("mi_{i}"),
                molecule_h5: tfile.as_ref().into(),
            })
            .collect();

        let args = MergeMoleculesStageInputs {
            libraries,
            sample_defs,
        };
        let stage = MergeMolecules;

        // Run the stage
        let testdir = TempDir::new().expect("Error creating tempdir");
        let res = stage.test_run(&testdir, args).unwrap();

        // Check results
        let gem_group_barcode_ranges: HashMap<String, Vec<u64>> = [
            ("0".to_string(), vec![0u64, 2]),
            ("1".to_string(), vec![2u64, 4]),
        ]
        .into_iter()
        .collect();
        assert_eq!(res.gem_group_barcode_ranges, gem_group_barcode_ranges);

        // Verify barcodes length for original and merged files
        let original_bcs: Vec<Vec<BarcodeContent>> = tfiles
            .iter()
            .map(|t| {
                MoleculeInfoReader::read_barcodes(t.path())
                    .unwrap()
                    .try_collect()
                    .unwrap()
            })
            .collect();
        for bcs in &original_bcs {
            assert_eq!(bcs.len(), 3);
        }

        let merged_bcs: Vec<BarcodeContent> =
            MoleculeInfoReader::read_barcodes(&res.merged_molecules)
                .unwrap()
                .try_collect()
                .unwrap();
        assert_eq!(merged_bcs.len(), 4);

        // Iter over all UMIs for both merged and original files
        let original_iter = zip(
            std::iter::repeat(0),
            MoleculeInfoIterator::new(tfiles[0].as_ref())
                .expect("Failed to open original molecule_info"),
        )
        .chain(zip(
            std::iter::repeat(1),
            MoleculeInfoIterator::new(tfiles[1].as_ref())
                .expect("failed to open original molecule_info"),
        ));

        let merged_iter = MoleculeInfoIterator::new(res.merged_molecules.as_ref())
            .expect("Failed to open merged molecule_info");

        for (merged, (fidx, original)) in zip(merged_iter, original_iter) {
            assert_eq!(
                merged_bcs[merged.barcode_idx as usize],
                original_bcs[fidx][original.barcode_idx as usize]
            );
            assert_eq!(merged.gem_group, fidx as u16);

            assert_eq!(merged.umi_data.library_idx, fidx as u16);
            assert_eq!(merged.umi_data.feature_idx, original.umi_data.feature_idx);
            assert_eq!(merged.umi_data.umi, original.umi_data.umi);
            assert_eq!(merged.umi_data.utype, original.umi_data.utype);
            assert_eq!(merged.umi_data.read_count, original.umi_data.read_count);
        }
    }
}
