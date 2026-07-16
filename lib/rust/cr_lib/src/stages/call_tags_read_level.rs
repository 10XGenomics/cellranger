//! Martian stage CALL_TAGS_READ_LEVEL
//! Assign cells to samples using their sample barcode sequence.
#![deny(missing_docs)]

use crate::barcode_overlap::{
    FRPGemBarcodeOverlapRow, GelBeadBarcodesPerProbeBarcode, ProbeBarcodeGelBeadGrouper,
    calculate_frp_gem_barcode_overlap,
};
use crate::detect_chemistry::probe_bc_pairing::get_rtl_and_ab_barcode_from_row;
use crate::read_level_multiplexing::{
    get_barcodes_per_multiplexing_identifier, get_umi_per_multiplexing_identifier,
    get_umi_per_multiplexing_identifier_for_feature_type,
};
use anyhow::Result;
use barcode::whitelist::{
    BarcodeId, RTLMultiplexingBarcodeType, categorize_rtl_multiplexing_barcode_id,
};
use barcode::{BarcodeConstruct, BcSeq, GelBeadAndProbeConstruct};
use cr_h5::count_matrix::{CountMatrix, CountMatrixFile, CountMatrixStreaming};
use cr_types::chemistry::{BarcodeKind, BarcodeReadComponent, ChemistryDefs};
use cr_types::reference::feature_reference::FeatureType;
use cr_types::types::FeatureBarcodeType;
use cr_types::utils::calculate_median_of_sorted;
use cr_types::{
    BarcodeMultiplexingType, CrMultiGraph, Fingerprint, LibraryType, MetricsFile, ReadLevel,
};
use itertools::{Itertools, process_results};
use martian::{MartianRover, MartianStage, MartianVoid, Resource, StageDef};
use martian_derive::{MartianStruct, make_mro};
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::tabular_file::CsvFile;
use martian_filetypes::{FileTypeRead, FileTypeWrite};
use metric::{TxHashMap, TxHashSet, join_metric_name};
use serde::{Deserialize, Serialize};
use serde_json::{Map, Value};
use std::collections::HashMap;
use std::ops::{Range, RangeTo};

/// Partitions a barcode sequence into its gel bead and sample barcode ranges.
struct BarcodePartition {
    gel_bead_range: RangeTo<usize>,
    sample_barcode_range: Range<usize>,
}

// Helper struct for barcode partition calculation
struct SampleBarcodeStartAndGelBeadEnd {
    sample_barcode_start: usize,
    gel_bead_end: usize,
}

impl BarcodePartition {
    /// Extracts the sample barcode start position and gel bead end position from a barcode construct.
    fn get_sample_barcode_start_and_gel_bead_end(
        barcode_construct: &BarcodeConstruct<&BarcodeReadComponent>,
    ) -> SampleBarcodeStartAndGelBeadEnd {
        match barcode_construct {
            BarcodeConstruct::GelBeadAndProbe(x) => match x.probe.kind() {
                BarcodeKind::Overhang => SampleBarcodeStartAndGelBeadEnd {
                    sample_barcode_start: x.probe.offset(),
                    gel_bead_end: x.gel_bead.length(),
                },
                BarcodeKind::RightProbe => SampleBarcodeStartAndGelBeadEnd {
                    sample_barcode_start: x.gel_bead.length(),
                    gel_bead_end: x.gel_bead.length(),
                },
                _ => panic!("Unsupported probe barcode kind."),
            },
            BarcodeConstruct::SegmentedGelBeadAndSample(x) => {
                let gb_len = x.gel_bead_segments.iter().map(|s| s.length()).sum();
                SampleBarcodeStartAndGelBeadEnd {
                    sample_barcode_start: gb_len,
                    gel_bead_end: gb_len,
                }
            }
            _ => {
                panic!("Barcode partitioning is unsupported for this construct")
            }
        }
    }

    /// Calculate barcode partition from chemistry barcode constructs.
    fn from_constructs(
        barcode_constructs: &HashMap<&LibraryType, BarcodeConstruct<&BarcodeReadComponent>>,
    ) -> Self {
        let (gel_bead_range, sample_barcode_range) = barcode_constructs
            .values()
            .map(|construct| {
                let result = Self::get_sample_barcode_start_and_gel_bead_end(construct);
                (
                    ..result.gel_bead_end,
                    result.sample_barcode_start
                        ..result.sample_barcode_start + construct.option_probe().unwrap().length(),
                )
            })
            .dedup()
            .exactly_one()
            .unwrap();

        Self {
            gel_bead_range,
            sample_barcode_range,
        }
    }

    /// Extract gel bead barcode from a full barcode.
    fn extract_gel_bead(&self, barcode: &[u8]) -> BcSeq {
        BcSeq::from_bytes(&barcode[self.gel_bead_range])
    }

    /// Extract sample barcode from a full barcode.
    fn extract_sample(&self, barcode: &[u8]) -> BcSeq {
        BcSeq::from_bytes(&barcode[self.sample_barcode_range.clone()])
    }
}

/// The Martian stage inputs.
#[derive(Clone, Deserialize, MartianStruct)]
pub struct CallTagsReadLevelStageInputs {
    pub chemistry_defs: ChemistryDefs,
    pub raw_feature_bc_matrix: CountMatrixFile,
    pub filtered_feature_bc_matrix: CountMatrixFile,
    pub multi_graph: JsonFile<CrMultiGraph>,
}

/// The Martian stage outputs.
#[derive(Serialize, MartianStruct)]
pub struct CallTagsReadLevelStageOutputs {
    pub barcodes_per_tag: Option<JsonFile<TxHashMap<BarcodeId, Vec<String>>>>,
    pub frp_gem_barcode_overlap: Option<CsvFile<FRPGemBarcodeOverlapRow>>,
    pub summary: Option<MetricsFile>,
}

impl CallTagsReadLevelStageOutputs {
    /// Build metrics and outputs for CALL_TAGS_READ_LEVEL stage.
    /// Handles both OH/CMX (with valid_barcodes) and RTL (with overlap data).
    fn build(
        rover: &MartianRover,
        barcodes_per_tag_file: JsonFile<TxHashMap<BarcodeId, Vec<String>>>,
        gel_bead_barcodes_per_probe_barcode: GelBeadBarcodesPerProbeBarcode<BcSeq>,
        umi_per_sample_barcode: Vec<(String, Value)>,
        valid_barcodes_per_sample_barcode: TxHashMap<BarcodeId, Vec<String>>,
        frp_gem_barcode_overlap: Option<Vec<FRPGemBarcodeOverlapRow>>,
    ) -> Result<Self> {
        // Build metrics
        let metrics: Value = CallTagsReadLevelMetrics::new(
            valid_barcodes_per_sample_barcode,
            gel_bead_barcodes_per_probe_barcode,
            frp_gem_barcode_overlap.as_deref(),
        )
        .to_map()
        .into_iter()
        .chain(umi_per_sample_barcode)
        .collect();

        // Write summary
        let summary: MetricsFile = rover.make_path("summary");
        summary.write(&metrics)?;

        // Conditionally write overlap file if it exists
        let frp_gem_barcode_overlap_file = if let Some(overlap_data) = frp_gem_barcode_overlap {
            let filename: CsvFile<_> = rover.make_path("frp_gem_barcode_overlap");
            filename.write(&overlap_data)?;
            Some(filename)
        } else {
            None
        };

        // Build outputs
        Ok(Self {
            barcodes_per_tag: Some(barcodes_per_tag_file),
            frp_gem_barcode_overlap: frp_gem_barcode_overlap_file,
            summary: Some(summary),
        })
    }
}

/// Martian stage CALL_TAGS_READ_LEVEL
/// Assign valid bcs to samples using their sample barcode sequence.
pub struct CallTagsReadLevel;

/// Read-level multiplexing metrics.
#[derive(Serialize)]
struct CallTagsReadLevelMetrics {
    /// The number of distinct filtered gel bead barcodes.
    filtered_gel_bead_barcodes_count: i64,

    /// The number of valid barcodes per sample barcode
    valid_barcodes_per_sample_barcode: TxHashMap<BarcodeId, i64>,

    /// The number of filtered barcodes for each sample barcode.
    filtered_barcodes_per_sample_barcode: TxHashMap<BarcodeId, i64>,

    /// The overlap coefficients of the sets of gel bead barcodes for each pair of probe barcodes.
    #[serde(skip_serializing_if = "Option::is_none")]
    probe_barcode_overlap_coefficients: Option<TxHashMap<String, f64>>,
}

impl CallTagsReadLevelMetrics {
    /// Calculate the sample barcode metrics.
    fn new(
        valid_barcodes_per_sample_barcode: TxHashMap<BarcodeId, Vec<String>>,
        gel_bead_barcodes_per_probe_barcode: GelBeadBarcodesPerProbeBarcode<BcSeq>,
        frp_gem_barcode_overlap: Option<&[FRPGemBarcodeOverlapRow]>,
    ) -> Self {
        let valid_barcodes_per_sample_barcode = valid_barcodes_per_sample_barcode
            .into_iter()
            .map(|(sample_barcode, gel_bead_barcodes)| {
                (sample_barcode, gel_bead_barcodes.len() as i64)
            })
            .collect();

        let filtered_barcodes_per_sample_barcode = gel_bead_barcodes_per_probe_barcode
            .iter()
            .map(|(sample_barcode_id, gel_bead_barcodes)| {
                (*sample_barcode_id, gel_bead_barcodes.len() as i64)
            })
            .collect();

        let filtered_gel_bead_barcodes_count = gel_bead_barcodes_per_probe_barcode
            .into_values()
            .reduce(|mut acc, gel_bead_barcodes| {
                acc.extend(gel_bead_barcodes);
                acc
            })
            .unwrap_or_default()
            .len() as i64;

        let probe_barcode_overlap_coefficients = frp_gem_barcode_overlap.map(|x| {
            x.iter()
                .map(|x| (format!("{}_{}", x.barcode1_id, x.barcode2_id), x.overlap))
                .collect()
        });

        Self {
            filtered_gel_bead_barcodes_count,
            valid_barcodes_per_sample_barcode,
            filtered_barcodes_per_sample_barcode,
            probe_barcode_overlap_coefficients,
        }
    }

    /// Convert the metrics to a JSON Map.
    fn to_map(&self) -> Map<String, Value> {
        if let Value::Object(map) = serde_json::to_value(self).unwrap() {
            map
        } else {
            unreachable!();
        }
    }
}

#[make_mro(volatile = strict)]
impl MartianStage for CallTagsReadLevel {
    type StageInputs = CallTagsReadLevelStageInputs;
    type StageOutputs = CallTagsReadLevelStageOutputs;
    type ChunkInputs = MartianVoid;
    type ChunkOutputs = MartianVoid;

    fn split(
        &self,
        args: Self::StageInputs,
        _rover: MartianRover,
    ) -> Result<StageDef<Self::ChunkInputs>> {
        let raw_matrix_gib = args
            .raw_feature_bc_matrix
            .estimate_mem_gib_without_counts()?;
        println!("raw_matrix_gib={raw_matrix_gib:.1}");
        // We end up creating a data structure that holds all of the barcodes,
        // so we need about twice as much memory to handle the raw matrix.
        Ok(StageDef::with_join_resource(Resource::with_mem_gb(
            4 + 2 * raw_matrix_gib.ceil() as isize,
        )))
    }

    fn main(
        &self,
        _args: Self::StageInputs,
        _chunk_args: Self::ChunkInputs,
        _rover: MartianRover,
    ) -> Result<Self::ChunkOutputs> {
        unreachable!()
    }

    /// Run the Martian stage CALL_TAGS_READ_LEVEL.
    fn join(
        &self,
        args: Self::StageInputs,
        _chunk_defs: Vec<Self::ChunkInputs>,
        _chunk_outs: Vec<Self::ChunkOutputs>,
        rover: MartianRover,
    ) -> Result<Self::StageOutputs> {
        let multi_graph = args.multi_graph.read()?;
        let Some(BarcodeMultiplexingType::ReadLevel(multiplexing_type)) =
            multi_graph.barcode_multiplexing_type()
        else {
            anyhow::bail!("This stage is only valid for read-level-multiplexed data.");
        };

        let barcode_constructs: HashMap<_, _> = args
            .chemistry_defs
            .iter()
            .map(|(library_type, chemistry_def)| {
                let barcode_construct = chemistry_def.barcode_construct();
                match barcode_construct {
                    BarcodeConstruct::GelBeadAndProbe(_)
                    | BarcodeConstruct::SegmentedGelBeadAndSample(_) => {
                        (library_type, barcode_construct)
                    }
                    BarcodeConstruct::GelBeadOnly(gb) => {
                        let oh = chemistry_def
                            .overhang_read_barcode()
                            .expect("Expected overhang read barcode for GelBeadOnly construct");
                        (
                            library_type,
                            BarcodeConstruct::new_gel_bead_and_probe(gb, oh),
                        )
                    }
                    BarcodeConstruct::SpatialSegmented(_) => {
                        panic!("This stage is not valid for spatial data.")
                    }
                }
            })
            .collect();

        let barcode_partition = BarcodePartition::from_constructs(&barcode_constructs);

        let sample_barcode_seq_to_id: TxHashMap<_, _> = barcode_constructs
            .into_values()
            .map(|x| {
                x.probe()
                    .whitelist()
                    .as_source()?
                    .as_translation_seq_to_id()
            })
            .flatten_ok()
            .map_ok(|(seq, id)| (BcSeq::from_bytes(seq.as_bytes()), id))
            .try_collect()?;

        let raw_matrix = args.raw_feature_bc_matrix.read_streaming()?;

        let valid_barcodes_per_sample_barcode = get_barcodes_per_multiplexing_identifier(
            &raw_matrix,
            &sample_barcode_seq_to_id,
            &barcode_partition.sample_barcode_range,
        );
        let barcodes_per_tag_file: JsonFile<_> = rover.make_path("barcodes_per_tag");
        barcodes_per_tag_file.write(&valid_barcodes_per_sample_barcode)?;

        // Setting it up this way for consistency with previous CALL_TAGS_OH. Can be simplified in the future if we output this metric for all feature types.
        let umi_per_sample_barcode: Vec<_> = (match multiplexing_type {
            ReadLevel::OH => get_umi_per_multiplexing_identifier_for_feature_type(
                &raw_matrix,
                &sample_barcode_seq_to_id,
                &barcode_partition.sample_barcode_range,
                FeatureType::Gene,
            )
            .map(|umis| [(FeatureType::Gene, umis)].into_iter().collect()),
            ReadLevel::RTL | ReadLevel::CMX => get_umi_per_multiplexing_identifier(
                &raw_matrix,
                &sample_barcode_seq_to_id,
                &barcode_partition.sample_barcode_range,
            ),
        })?
        .into_iter()
        .map(|(feature_type, umi_per_sample_barcode)| {
            (
                join_metric_name(feature_type, "umi_per_sample_barcode"),
                serde_json::to_value(umi_per_sample_barcode).unwrap(),
            )
        })
        .collect();
        let filtered_matrix = args.filtered_feature_bc_matrix.read_streaming()?;
        // Gather data for probe barcode overlap computations.
        let gel_bead_barcodes_per_probe_barcode = ProbeBarcodeGelBeadGrouper::group_all(
            filtered_matrix
                .barcodes()
                .iter()
                .map(|barcode| GelBeadAndProbeConstruct {
                    gel_bead: barcode_partition.extract_gel_bead(barcode.as_bytes()),
                    probe: barcode_partition.extract_sample(barcode.as_bytes()),
                }),
            &sample_barcode_seq_to_id,
        );

        if multiplexing_type == ReadLevel::OH || multiplexing_type == ReadLevel::CMX {
            return CallTagsReadLevelStageOutputs::build(
                &rover,
                barcodes_per_tag_file,
                gel_bead_barcodes_per_probe_barcode,
                umi_per_sample_barcode,
                valid_barcodes_per_sample_barcode,
                None,
            );
        }

        // RTL specific computation:

        let mut frp_gem_barcode_overlap =
            calculate_frp_gem_barcode_overlap(&gel_bead_barcodes_per_probe_barcode);
        // Determine if we have barcode pairings; if so, run suspicious pairing
        // detection.
        let ab_barcode_pairings = get_ab_barcode_pairings(&multi_graph);
        let rtl_ab_gem_barcode_overlap = handle_rtl_ab_pairings(
            &raw_matrix,
            &filtered_matrix,
            &barcode_partition,
            &sample_barcode_seq_to_id,
            &ab_barcode_pairings,
            &gel_bead_barcodes_per_probe_barcode,
        )?;
        // Combine the RTL+AB overlap data with the RTL+RTL overlap data.
        frp_gem_barcode_overlap.extend(rtl_ab_gem_barcode_overlap);

        CallTagsReadLevelStageOutputs::build(
            &rover,
            barcodes_per_tag_file,
            gel_bead_barcodes_per_probe_barcode,
            umi_per_sample_barcode,
            valid_barcodes_per_sample_barcode,
            Some(frp_gem_barcode_overlap),
        )
    }
}

fn get_ab_barcode_pairings(multi_graph: &CrMultiGraph) -> TxHashMap<BarcodeId, BarcodeId> {
    multi_graph
        .samples
        .iter()
        .flat_map(|sample| &sample.fingerprints)
        .filter_map(|fingerprint| {
            let Fingerprint::Tagged {
                tag_name,
                translated_tag_names,
                ..
            } = fingerprint
            else {
                return None;
            };

            translated_tag_names
                .iter()
                .filter(|x| {
                    categorize_rtl_multiplexing_barcode_id(x)
                        .expect("Missing Multiplexing Barcode!")
                        == RTLMultiplexingBarcodeType::Antibody
                })
                .at_most_one()
                .unwrap()
                .map(|translated_tag_name| {
                    (
                        BarcodeId::pack(tag_name),
                        BarcodeId::pack(translated_tag_name),
                    )
                })
        })
        .collect()
}
/// Run RTL+AB pairing detection if needed.
fn handle_rtl_ab_pairings(
    raw_matrix: &CountMatrixStreaming,
    filtered_matrix: &CountMatrixStreaming,
    barcode_partition: &BarcodePartition,
    sample_barcode_seq_to_id: &TxHashMap<BcSeq, BarcodeId>,
    barcode_pairings: &TxHashMap<BarcodeId, BarcodeId>,
    // The RTL probe barcode grouping data.
    gex_gel_bead_barcodes_per_probe_barcode: &GelBeadBarcodesPerProbeBarcode<BcSeq>,
) -> Result<Vec<FRPGemBarcodeOverlapRow>> {
    if barcode_pairings.is_empty() {
        return Ok(Vec::default());
    }
    let median_umi_per_cell =
        get_median_umi_per_cell(filtered_matrix, &barcode_partition.sample_barcode_range)?;

    detect_suspicious_rtl_ab_pairings(
        raw_matrix,
        barcode_partition,
        sample_barcode_seq_to_id,
        barcode_pairings,
        &median_umi_per_cell,
        gex_gel_bead_barcodes_per_probe_barcode,
    )
}

/// Compute overlap between RTL and AB probe barcodes.
///
/// The goal is to identify RTL+AB barcode overlap that might represent a
///
/// Use the provided barcode ID mapping to reverse-translate feature barcoding
/// counts back into their pre-translation ID.
/// Remove all barcodes per feature type with counts less than 10% of the median
/// UMI count per cell for that probe barcode and feature type. This excludes
/// all probe barcodes that were not observed in at least one cell, and cleans
/// up the background.
/// Compute the gel bead barcode overlap between all probe barcode pairings.
/// Filter out all pairings besides those between Gex+Ab.
fn detect_suspicious_rtl_ab_pairings(
    raw_matrix: &CountMatrixStreaming,
    barcode_partition: &BarcodePartition,
    sample_barcode_seq_to_id: &TxHashMap<BcSeq, BarcodeId>,
    // Mapping from Gex probe barcode ID back to FB probe barcode ID.
    barcode_pairings: &TxHashMap<BarcodeId, BarcodeId>,
    // The median number of UMI observed per cell per probe barcode/feature type.
    // If a feature type/barcode pair is missing, then there were no cells observed
    // containing it.
    median_umi_per_cell_per_probe_barcode: &TxHashMap<(FeatureType, BcSeq), usize>,
    // The RTL probe barcode grouping data.
    gex_gel_bead_barcodes_per_probe_barcode: &GelBeadBarcodesPerProbeBarcode<BcSeq>,
) -> Result<Vec<FRPGemBarcodeOverlapRow>> {
    // Mapping to reverse-translate antibody probe barcode sequences into IDs.
    let reverse_translation: TxHashMap<_, _> = sample_barcode_seq_to_id
        .iter()
        .map(|(bc, id)| (*bc, *barcode_pairings.get(id).unwrap_or(id)))
        .collect();
    let mut fb_barcode_groups = ProbeBarcodeGelBeadGrouper::new(&reverse_translation);

    for count in raw_matrix.counts()? {
        let count = count?;
        if count.feature.feature_type != FeatureType::Barcode(FeatureBarcodeType::Antibody) {
            continue;
        }
        let construct = GelBeadAndProbeConstruct {
            gel_bead: barcode_partition.extract_gel_bead(count.barcode.as_bytes()),
            probe: barcode_partition.extract_sample(count.barcode.as_bytes()),
        };
        fb_barcode_groups.group(construct, count.count as usize);
    }

    // Reverse-translate the median UMI data for the antibody barcodes.
    // This ensures the barcode IDs match the reverse-translated IDs we collected above.
    let median_umi_per_cell_per_probe_barcode: TxHashMap<_, _> =
        median_umi_per_cell_per_probe_barcode
            .iter()
            .filter_map(|((feature_type, probe_bc), median_umi)| {
                let probe_bc_id = match *feature_type {
                    FeatureType::Gene => sample_barcode_seq_to_id[probe_bc],
                    FeatureType::Barcode(FeatureBarcodeType::Antibody) => {
                        reverse_translation[probe_bc]
                    }
                    // TODO(CELLRANGER-7915): Detect suspicious pairings for Flex CRISPR
                    FeatureType::Barcode(
                        FeatureBarcodeType::Crispr | FeatureBarcodeType::Custom,
                    ) => return None,
                    FeatureType::Barcode(_) => unreachable!(),
                    FeatureType::ProteinExpression
                    | FeatureType::Peaks
                    | FeatureType::DeprecatedCodeword
                    | FeatureType::Unknown
                    | FeatureType::NegativeControlProbe
                    | FeatureType::NegativeControlCodeword
                    | FeatureType::UnassignedCodeword
                    | FeatureType::GenomicControl => unimplemented!(),
                };
                Some(((*feature_type, probe_bc_id), median_umi))
            })
            .collect();

    let mut ab_gel_bead_barcodes_per_probe_barcode = fb_barcode_groups.finish();

    // Filter the antibody data based on the median UMI in any called cell for this feature type.
    // Allow counts down to 10% of the median UMI count.
    // Remove probe barcodes that were never seen in a cell.
    let mut keys_to_remove = Vec::new();
    for (probe_id, gel_bead_counts) in &mut ab_gel_bead_barcodes_per_probe_barcode {
        // If we never saw this probe barcode in a cell, ignore it.
        let Some(median_umi_per_cell) = median_umi_per_cell_per_probe_barcode.get(&(
            FeatureType::Barcode(FeatureBarcodeType::Antibody),
            *probe_id,
        )) else {
            keys_to_remove.push(*probe_id);
            continue;
        };
        // Should never have a non-antibody probe barcode here after CELLRANGER-7501.
        // Note: this check was moved to come after the "seen in cells" check
        // above, because we may be analyzing a sub-set of all the samples in a
        // library. In this case, the probe barcodes associated with the samples
        // that are not participating in the analysis (and thus which do not
        // appear in the multi config graph) do not end up in the barcode
        // pairings that make it into this function. However, since those
        // probe barcodes are not associated with a sample, they should never
        // end up being called as cells.
        assert_eq!(
            categorize_rtl_multiplexing_barcode_id(probe_id).unwrap(),
            RTLMultiplexingBarcodeType::Antibody,
            "probe ID {probe_id} was expected to be an antibody probe barcode",
        );
        // Remove all gel beads with counts below the threshold.
        let min_count = (0.1 * **median_umi_per_cell as f64).round() as usize;
        gel_bead_counts.retain(|_, count| *count >= min_count);
    }
    for key in keys_to_remove {
        ab_gel_bead_barcodes_per_probe_barcode.remove(&key);
    }

    // Combine the Gex and Ab groupings/sum the counts.
    let mut gel_bead_barcodes_per_probe_barcode = ab_gel_bead_barcodes_per_probe_barcode;

    for (probe_bc, gel_bead_bcs) in gex_gel_bead_barcodes_per_probe_barcode {
        let gel_bead_counts = gel_bead_barcodes_per_probe_barcode
            .entry(*probe_bc)
            .or_default();
        for (gel_bead_bc, count) in gel_bead_bcs {
            *gel_bead_counts.entry(*gel_bead_bc).or_default() += count;
        }
    }

    let pairings_to_ignore: TxHashSet<_> = barcode_pairings.iter().map(|(k, v)| (*k, *v)).collect();

    // Filter out all pairings besides Gex+Ab and remove explicitly paired barcodes,
    // and return an ordered collection of pairings.
    Ok(
        calculate_frp_gem_barcode_overlap(&gel_bead_barcodes_per_probe_barcode)
            .into_iter()
            .filter(|row| {
                // Remove pairings besides RTL+AB.
                let Some(pairing) = get_rtl_and_ab_barcode_from_row(row) else {
                    return false;
                };
                !pairings_to_ignore.contains(&pairing)
            })
            .map(|mut row| {
                // Canonicalize so pairings are always ordered as RTL+AB.
                if categorize_rtl_multiplexing_barcode_id(&row.barcode1_id).unwrap()
                    != RTLMultiplexingBarcodeType::Gene
                {
                    row.swap_order();
                }
                row
            })
            .sorted_by_key(|row| (row.barcode1_id, row.barcode2_id))
            .collect(),
    )
}

/// Get median UMI per cell for each probe barcode/feature type pair.
/// Barcodes with 0 UMI for the feature type are not included in the median.
fn get_median_umi_per_cell(
    filtered_matrix: &CountMatrixStreaming,
    sample_barcode_range: &Range<usize>,
) -> Result<TxHashMap<(FeatureType, BcSeq), usize>> {
    process_results(filtered_matrix.counts()?, |counts| {
        counts
            .map(|x| ((x.feature.feature_type, x.barcode), x.count))
            .into_grouping_map()
            .sum()
            .into_iter()
            .filter(|&(_, count)| count > 0)
            .map(|((feature_type, barcode), count)| {
                let probe_bc = BcSeq::from_bytes(&barcode.as_bytes()[sample_barcode_range.clone()]);
                ((feature_type, probe_bc), count)
            })
            .into_group_map()
            // Find medians.
            .into_iter()
            .map(|(id, mut counts)| {
                counts.sort();
                (id, calculate_median_of_sorted(&counts).unwrap() as usize)
            })
            .collect()
    })
}
