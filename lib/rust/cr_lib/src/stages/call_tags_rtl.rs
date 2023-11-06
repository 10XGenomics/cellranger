//! Martian stage CALL_TAGS_RTL
//! Assign cells to samples using their probe barcode sequence.
use crate::barcode_overlap::{
    calculate_frp_gem_barcode_overlap, FRPGemBarcodeOverlapRow, GelBeadBarcodesPerProbeBarcode,
    ProbeBarcodeGelBeadGrouper,
};
use crate::count_matrix::{CountMatrix, CountMatrixFile, LazyCountMatrix};
use crate::detect_chemistry::probe_bc_pairing::get_rtl_and_ab_barcode_from_row;
use crate::read_level_multiplexing::{
    get_barcodes_per_multiplexing_identifier, get_umi_per_multiplexing_identifier,
};
use anyhow::Result;
use barcode::whitelist::{categorize_multiplexing_barcode_id, BarcodeId, MultiplexingBarcodeType};
use barcode::{BarcodeConstruct, BcSegSeq, GelBeadAndProbeConstruct, WhitelistSource};
use cr_types::chemistry::ChemistryDef;
use cr_types::utils::calculate_median_of_sorted;
use cr_types::{CrMultiGraph, FeatureType, Fingerprint, MetricsFile};
use itertools::Itertools;
use json_report_derive::JsonReport;
use martian::prelude::{MartianRover, MartianStage};
use martian::{MartianVoid, Resource, StageDef};
use martian_derive::{make_mro, MartianStruct};
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::tabular_file::CsvFile;
use martian_filetypes::{FileTypeRead, FileTypeWrite};
use metric::{JsonReporter, Metric, TxHashMap, TxHashSet};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::ops::{Range, RangeTo};

/// The Martian stage inputs.
#[derive(Clone, Deserialize, MartianStruct)]
pub struct CallTagsRTLStageInputs {
    pub chemistry_def: ChemistryDef,
    pub raw_feature_bc_matrix: CountMatrixFile,
    pub filtered_feature_bc_matrix: CountMatrixFile,
    pub multi_graph: JsonFile<CrMultiGraph>,
}

/// The Martian stage outputs.
#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct CallTagsRTLStageOutputs {
    pub barcodes_per_tag: Option<JsonFile<TxHashMap<BarcodeId, Vec<String>>>>,
    pub frp_gem_barcode_overlap: Option<CsvFile<FRPGemBarcodeOverlapRow>>,
    pub summary: Option<MetricsFile>,
}

/// Martian stage CALL_TAGS_RTL
/// Assign valid bcs to samples using their probe barcode sequence.
pub struct CallTagsRTL;

/// Probe barcode metrics.
#[derive(JsonReport)]
#[json_report(extend = "extra_reports")]
#[allow(non_snake_case)]
struct CallTagsRTLMetrics {
    /// The number of UMI per probe barcode for each feature type.
    #[json_report(skip)]
    umi_per_probe_barcode: HashMap<FeatureType, TxHashMap<BarcodeId, i64>>,

    /// The number of distinct filtered gel bead barcodes.
    filtered_gel_bead_barcodes_count: i64,

    /// The number of filtered barcodes for each probe barcode.
    #[json_report(block)]
    filtered_barcodes_per_probe_barcode: TxHashMap<BarcodeId, i64>,

    /// The overlap coefficients of the sets of gel bead barcodes for each pair of probe barcodes.
    #[json_report(block)]
    probe_barcode_overlap_coefficients: TxHashMap<String, f64>,
}

impl CallTagsRTLMetrics {
    /// Calculate the probe barcode metrics.
    fn new(
        umi_per_probe_barcode: HashMap<FeatureType, TxHashMap<BarcodeId, i64>>,
        gel_bead_barcodes_per_probe_barcode: GelBeadBarcodesPerProbeBarcode,
        frp_gem_barcode_overlap: &[FRPGemBarcodeOverlapRow],
    ) -> Self {
        let filtered_barcodes_per_probe_barcode = gel_bead_barcodes_per_probe_barcode
            .iter()
            .map(|(probe_barcode_id, gel_bead_barcodes)| {
                (*probe_barcode_id, gel_bead_barcodes.len() as i64)
            })
            .collect();

        let filtered_gel_bead_barcodes_count = gel_bead_barcodes_per_probe_barcode
            .into_values()
            .reduce(|mut acc, gel_bead_barcodes| {
                acc.extend(gel_bead_barcodes.into_iter());
                acc
            })
            .unwrap_or_default()
            .len() as i64;

        let probe_barcode_overlap_coefficients = frp_gem_barcode_overlap
            .iter()
            .map(|x| (format!("{}_{}", x.barcode1_id, x.barcode2_id), x.overlap))
            .collect();

        Self {
            umi_per_probe_barcode,
            filtered_gel_bead_barcodes_count,
            filtered_barcodes_per_probe_barcode,
            probe_barcode_overlap_coefficients,
        }
    }

    /// Emit umi_per_probe_barcode with the feature type prefix.
    fn extra_reports(&self) -> JsonReporter {
        self.umi_per_probe_barcode
            .iter()
            .map(|(feature_type, umi_per_probe_barcode)| {
                (
                    feature_type.join("umi_per_probe_barcode").into_owned(),
                    umi_per_probe_barcode,
                )
            })
            .collect()
    }
}

#[make_mro(volatile = strict)]
impl MartianStage for CallTagsRTL {
    type StageInputs = CallTagsRTLStageInputs;
    type StageOutputs = CallTagsRTLStageOutputs;
    type ChunkInputs = MartianVoid;
    type ChunkOutputs = MartianVoid;

    fn split(
        &self,
        args: Self::StageInputs,
        _rover: MartianRover,
    ) -> Result<StageDef<Self::ChunkInputs>> {
        // The raw matrix is larger than the filtered matrix, and we ensure
        // in the implementation of this stage that we never load both matrices
        // into memory at the same time.
        let raw_matrix_gib = args.raw_feature_bc_matrix.estimate_mem_gib()?;
        println!("raw_matrix_gib={raw_matrix_gib:.1}");
        Ok(StageDef::with_join_resource(Resource::with_mem_gb(
            2 + raw_matrix_gib.ceil() as isize,
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

    /// Run the Martian stage CALL_TAGS_RTL.
    fn join(
        &self,
        args: Self::StageInputs,
        _chunk_defs: Vec<Self::ChunkInputs>,
        _chunk_outs: Vec<Self::ChunkOutputs>,
        rover: MartianRover,
    ) -> Result<Self::StageOutputs> {
        let bc_construct = if let BarcodeConstruct::GelBeadAndProbe(bcc) =
            args.chemistry_def.barcode_construct()
        {
            bcc
        } else {
            return Ok(Self::StageOutputs {
                barcodes_per_tag: None,
                frp_gem_barcode_overlap: None,
                summary: None,
            });
        };

        let probe_barcode_range = bc_construct.gel_bead.length()
            ..(bc_construct.gel_bead.length() + bc_construct.probe.length());
        let gel_bead_barcode_range = ..bc_construct.gel_bead.length();
        let probe_barcode_seq_to_id =
            WhitelistSource::from_spec(bc_construct.probe.whitelist(), true, None)?
                .as_translation_seq_to_id()?;

        let filtered_matrix = args.filtered_feature_bc_matrix.read()?;

        // Gather data for probe barcode overlap computations.
        let gel_bead_barcodes_per_probe_barcode = ProbeBarcodeGelBeadGrouper::group_all(
            filtered_matrix
                .barcodes()
                .iter()
                .map(|barcode| GelBeadAndProbeConstruct {
                    gel_bead: BcSegSeq::from_bytes(&barcode.as_bytes()[gel_bead_barcode_range]),
                    probe: BcSegSeq::from_bytes(&barcode.as_bytes()[probe_barcode_range.clone()]),
                }),
            &probe_barcode_seq_to_id,
        );

        let mut frp_gem_barcode_overlap =
            calculate_frp_gem_barcode_overlap(&gel_bead_barcodes_per_probe_barcode);

        // Determine if we have barcode pairings; if so, run suspicious pairing
        // detection.
        let multi_graph = args.multi_graph.read()?;
        let barcode_pairings: TxHashMap<_, _> = multi_graph
            .samples
            .iter()
            .flat_map(|sample| &sample.fingerprints)
            .filter_map(|fingerprint| {
                let Fingerprint::Tagged {tag_name, translated_tag_names, ..} = fingerprint else {return None;};
                assert!(translated_tag_names.len() < 2);
                Some(translated_tag_names.iter().map(|translated_tag_name| (BarcodeId::pack(tag_name), BarcodeId::pack(translated_tag_name))))
            }).flatten().collect();

        let raw_matrix = LazyCountMatrix::new(args.raw_feature_bc_matrix.clone());

        let rtl_ab_gem_barcode_overlap = handle_rtl_ab_pairings(
            &raw_matrix,
            filtered_matrix,
            gel_bead_barcode_range,
            probe_barcode_range.clone(),
            &probe_barcode_seq_to_id,
            &barcode_pairings,
            &gel_bead_barcodes_per_probe_barcode,
        )?;

        // Combine the RTL+AB overlap data with the RTL+RTL overlap data.
        frp_gem_barcode_overlap.extend(rtl_ab_gem_barcode_overlap);

        let frp_gem_barcode_overlap_filename: CsvFile<_> =
            rover.make_path("frp_gem_barcode_overlap");
        frp_gem_barcode_overlap_filename.write(&frp_gem_barcode_overlap)?;

        let barcodes_per_tag_file: JsonFile<_> = rover.make_path("barcodes_per_tag");
        barcodes_per_tag_file.write(&get_barcodes_per_multiplexing_identifier(
            raw_matrix.loaded()?,
            &probe_barcode_seq_to_id,
            &probe_barcode_range,
        )?)?;

        let umi_per_probe_barcode = get_umi_per_multiplexing_identifier(
            raw_matrix.loaded()?,
            &probe_barcode_seq_to_id,
            &probe_barcode_range,
        );

        let metrics = CallTagsRTLMetrics::new(
            umi_per_probe_barcode,
            gel_bead_barcodes_per_probe_barcode,
            &frp_gem_barcode_overlap,
        );
        let metrics_file = MetricsFile::from_reporter(&rover, "summary", &metrics)?;

        Ok(Self::StageOutputs {
            barcodes_per_tag: Some(barcodes_per_tag_file),
            frp_gem_barcode_overlap: Some(frp_gem_barcode_overlap_filename),
            summary: Some(metrics_file),
        })
    }
}

/// Run RTL+AB pairing detection if needed.
fn handle_rtl_ab_pairings(
    raw_matrix: &LazyCountMatrix,
    // We take ownership of this matrix so we can drop it before loading the raw matrix.
    filtered_matrix: CountMatrix,
    gel_bead_barcode_range: RangeTo<usize>,
    probe_barcode_range: Range<usize>,
    probe_barcode_seq_to_id: &TxHashMap<BcSegSeq, BarcodeId>,
    barcode_pairings: &TxHashMap<BarcodeId, BarcodeId>,
    // The RTL probe barcode grouping data.
    gex_gel_bead_barcodes_per_probe_barcode: &GelBeadBarcodesPerProbeBarcode,
) -> Result<Vec<FRPGemBarcodeOverlapRow>> {
    if barcode_pairings.is_empty() {
        return Ok(Default::default());
    }
    let median_umi_per_cell =
        get_median_umi_per_cell(&filtered_matrix, probe_barcode_range.clone());

    // Drop the filtered matrix to free memory before loading the raw matrix.
    std::mem::drop(filtered_matrix);

    Ok(detect_suspicious_rtl_ab_pairings(
        raw_matrix.loaded()?,
        gel_bead_barcode_range,
        probe_barcode_range,
        probe_barcode_seq_to_id,
        barcode_pairings,
        &median_umi_per_cell,
        gex_gel_bead_barcodes_per_probe_barcode,
    ))
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
    raw_matrix: &CountMatrix,
    gel_bead_barcode_range: RangeTo<usize>,
    probe_barcode_range: Range<usize>,
    probe_barcode_seq_to_id: &TxHashMap<BcSegSeq, BarcodeId>,
    // Mapping from Gex probe barcode ID back to FB probe barcode ID.
    barcode_pairings: &TxHashMap<BarcodeId, BarcodeId>,
    // The median number of UMI observed per cell per probe barcode/feature type.
    // If a feature type/barcode pair is missing, then there were no cells observed
    // containing it.
    median_umi_per_cell_per_probe_barcode: &TxHashMap<(FeatureType, BcSegSeq), usize>,
    // The RTL probe barcode grouping data.
    gex_gel_bead_barcodes_per_probe_barcode: &GelBeadBarcodesPerProbeBarcode,
) -> Vec<FRPGemBarcodeOverlapRow> {
    // Mapping to reverse-translate antibody probe barcode sequences into IDs.
    let reverse_translation: TxHashMap<_, _> = probe_barcode_seq_to_id
        .iter()
        .map(|(bc, id)| (*bc, *barcode_pairings.get(id).unwrap_or(id)))
        .collect();
    let mut fb_barcode_groups = ProbeBarcodeGelBeadGrouper::new(&reverse_translation);

    for count in raw_matrix.counts() {
        if count.feature.feature_type == FeatureType::Gene {
            continue;
        } else {
            assert_eq!(FeatureType::Antibody, count.feature.feature_type);
        }
        let construct = GelBeadAndProbeConstruct {
            gel_bead: BcSegSeq::from_bytes(&count.barcode.as_bytes()[gel_bead_barcode_range]),
            probe: BcSegSeq::from_bytes(&count.barcode.as_bytes()[probe_barcode_range.clone()]),
        };
        fb_barcode_groups.group(construct, count.count as usize);
    }

    // Reverse-translate the median UMI data for the antibody barcodes.
    // This ensures the barcode IDs match the reverse-translated IDs we collected above.
    let median_umi_per_cell_per_probe_barcode: TxHashMap<_, _> =
        median_umi_per_cell_per_probe_barcode
            .iter()
            .map(|((feature_type, probe_bc), median_umi)| {
                let probe_bc_id = if *feature_type == FeatureType::Antibody {
                    reverse_translation[probe_bc]
                } else {
                    assert_eq!(*feature_type, FeatureType::Gene);
                    probe_barcode_seq_to_id[probe_bc]
                };
                ((*feature_type, probe_bc_id), median_umi)
            })
            .collect();

    let mut ab_gel_bead_barcodes_per_probe_barcode = fb_barcode_groups.finish();

    // Filter the antibody data based on the median UMI in any called cell for this feature type.
    // Allow counts down to 10% of the median UMI count.
    // Remove probe barcodes that were never seen in a cell, and also remove any
    // stray GEX probe barcodes that we ended up with (presumably due to barcode
    // correction artifacts - CELLRANGER-7501).
    let mut keys_to_remove = Vec::new();
    for (probe_id, gel_bead_counts) in &mut ab_gel_bead_barcodes_per_probe_barcode {
        // If this is not actually an antibody probe barcode, discard it.
        if categorize_multiplexing_barcode_id(probe_id) != MultiplexingBarcodeType::Antibody {
            keys_to_remove.push(*probe_id);
            continue;
        }
        // If we never saw this probe barcode in a cell, ignore it.
        let Some(median_umi_per_cell) = median_umi_per_cell_per_probe_barcode.get(&(FeatureType::Antibody, *probe_id)) else {
                keys_to_remove.push(*probe_id);
                continue;
            };
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
    calculate_frp_gem_barcode_overlap(&gel_bead_barcodes_per_probe_barcode)
        .into_iter()
        .filter(|row| {
            // Remove pairings besides RTL+AB.
            let Some(pairing) = get_rtl_and_ab_barcode_from_row(row) else { return false; };
            !pairings_to_ignore.contains(&pairing)
        })
        .map(|mut row| {
            // Canonicalize so pairings are always ordered as RTL+AB.
            if categorize_multiplexing_barcode_id(&row.barcode1_id) != MultiplexingBarcodeType::RTL
            {
                row.swap_order();
            }
            row
        })
        .sorted_by_key(|row| (row.barcode1_id, row.barcode2_id))
        .collect()
}

/// Get median UMI per cell for each probe barcode/feature type pair.
/// Barcodes with 0 UMI for the feature type are not included in the median.
fn get_median_umi_per_cell(
    filtered_matrix: &CountMatrix,
    probe_barcode_range: Range<usize>,
) -> TxHashMap<(FeatureType, BcSegSeq), usize> {
    filtered_matrix
        .counts()
        .map(|x| ((x.feature.feature_type, x.barcode), x.count))
        .into_grouping_map()
        .sum()
        .into_iter()
        .filter(|&(_, count)| count > 0)
        .map(|((feature_type, barcode), count)| {
            let probe_bc = BcSegSeq::from_bytes(&barcode.as_bytes()[probe_barcode_range.clone()]);
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
}
