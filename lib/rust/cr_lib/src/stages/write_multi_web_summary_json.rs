//! Martian stage WRITE_MULTI_WEB_SUMMARY_JSON
//! Write a JSON file with the metrics and plots data for the websummary.

use crate::itertools::Itertools;
use crate::stages::build_vdj_ws_contents::{VdjWsContents, VdjWsContentsFormat};
use crate::stages::compute_antigen_vdj_metrics::{AntigenVdjMetrics, AntigenVdjMetricsFormat};
use crate::stages::detect_chemistry::DetectedProbeBarcodePairingFile;
use crate::stages::parse_multi_config::{CommonInputs, CountInputs};
use crate::{SequencingMetricsFormat, SvgFile};
use anyhow::{bail, Result};
use barcode::whitelist::{categorize_multiplexing_barcode_id, MultiplexingBarcodeType};
use cr_types::chemistry::ChemistryName;
use cr_types::reference::feature_reference::{FeatureConfig, SpecificityControls};
use cr_types::reference::reference_info::MULTI_GENOME_SEPARATOR;
use cr_types::rna_read::LegacyLibraryType;
use cr_types::{
    AlignerParam, CellMultiplexingType, CrMultiGraph, Fingerprint, HasMultiplexing, Sample,
    SampleAssignment, TargetingMethod,
};
use cr_websummary::alert::AlertContext;
use cr_websummary::multi::antigen::{clonotype_specificity_heatmap, AntigenSpecificityRow};
use cr_websummary::multi::plots::{
    format_barcode_rank_plot, format_histogram, format_jibes_biplots, format_tags_on_tsne_plot,
    format_umi_on_tsne_plot, library_median_genes_plot_from_metrics,
    library_sequencing_saturation_plot_from_metrics, sample_median_genes_plot_from_metrics,
    targeted_enrichment_plot, PlotType,
};
use cr_websummary::multi::svg::SvgGraph;
use cr_websummary::multi::tables::{
    AntibodyLibraryMappingMetricsRow, AntibodyLibraryMappingMetricsTable,
    AntibodyPhysicalLibraryMetricsRow, AntibodyPhysicalLibraryMetricsTable,
    AntibodySampleHeroMetricsTable, AntibodySampleMappingMetricsTable,
    AntigenPhysicalLibraryMetricsTable, AntigenSampleHeroMetricsRow, AntigenSampleHeroMetricsTable,
    CrisprPhysicalLibraryMetricsRow, CrisprPhysicalLibraryMetricsTable, CrisprSampleHeroMetricsRow,
    CrisprSampleHeroMetricsTable, CustomFeaturePhysicalLibraryMetricsRow,
    CustomFeaturePhysicalLibraryMetricsTable, CustomFeatureSampleHeroMetricsRow,
    CustomFeatureSampleHeroMetricsTable, GdnaMetricsTable, GexLibraryMappingMetricsRow,
    GexLibraryMappingMetricsTable, GexLibraryTargetedEnrichmentAlertsRow,
    GexLibraryTargetedEnrichmentAlertsTable, GexLibraryTargetedEnrichmentMetricsRow,
    GexLibraryTargetedEnrichmentMetricsTable, GexPhysicalLibraryMetricsRow,
    GexPhysicalLibraryMetricsTable, GexSampleCellMetricsRow, GexSampleCellMetricsTable,
    GexSampleHeroMetricsRow, GexSampleHeroMetricsTable, GexSampleMappingMetricsRow,
    GexSampleMappingMetricsTable, LibraryCellMetricsRow, LibraryCellMetricsTable,
    MultiplexingCmoMetricsRow, MultiplexingCmoMetricsTable, MultiplexingLibraryCellMetricsRow,
    MultiplexingLibraryCellMetricsTable, MultiplexingPhysicalLibraryMetricsRow,
    MultiplexingPhysicalLibraryMetricsTable, MultiplexingSampleAssignmentsRow,
    MultiplexingSampleAssignmentsTable, RtlLibraryMappingMetricsRow, RtlLibraryMappingMetricsTable,
    RtlProbeBarcodeMetricsRow, RtlProbeBarcodeMetricsTable, RtlSampleCellMetricsRow,
    RtlSampleCellMetricsTable, RtlSampleMappingMetricsRow, RtlSampleMappingMetricsTable,
    SequencingMetricsTable,
};
use cr_websummary::multi::websummary::{
    AntibodyOrAntigen, CountParametersTable, ExperimentalDesign, GexOrRtl,
    GexOrRtlLibraryMappingMetricsTable, GexOrRtlSampleMappingMetricsTable,
    LibraryAntibodyOrAntigenWebSummary, LibraryCmoWebSummary, LibraryCrisprWebSummary,
    LibraryCustomFeatureWebSummary, LibraryGexWebSummary, LibraryHeaderInfo, LibraryWebSummary,
    MismatchedProbeBarcodePairings, MultiDiagnostics, MultiSharedResource, MultiWebSummary,
    MultiWebSummaryData, SampleAntibodyWebSummary, SampleAntigenWebSummary, SampleCellMetricsTable,
    SampleCrisprWebSummary, SampleCustomFeatureWebSummary, SampleGexWebSummary, SampleHeaderInfo,
    SampleWebSummary, UMI_PER_PROBE_BARCODE_BACKGROUND_THRESHOLD,
};
use cr_websummary::{
    CountAndPercent, FloatAsInt, Percent, PlotlyChart, RawChartWithHelp, Tab, WsSample,
};
use martian::prelude::*;
use martian_derive::{make_mro, MartianStruct};
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::tabular_file::CsvFile;
use martian_filetypes::FileTypeRead;
use metric::{PercentMetric, TxHashMap, TxHashSet};
use multi::config::{
    ChemistryParam, MultiConfigCsv, MultiConfigCsvFile, ProbeBarcodeIterationMode,
    PROBE_BARCODE_ID_GROUPING,
};
use ordered_float::NotNan;
use serde::{Deserialize, Serialize};
use serde_json::value::Value;
use serde_json::{self, Number};
use std::cmp::Ordering;
use std::collections::HashSet;
use std::hash::Hash;
use std::iter::FromIterator;
use std::string::ToString;
use ChemistryName::{ArcV1, FivePrimePE, FivePrimeR1, FivePrimeR2, ThreePrimeV3LT};

pub const GEM_BARCODE_OVERLAP_ALERT_THRESHOLD: f64 = 0.6;

pub struct WriteMultiWebSummaryJson;

#[derive(Clone, Serialize, Deserialize)]
pub struct JibesBiplotHistogramData {
    pub biplot: Value,
    pub histogram: Value,
    #[serde(rename = "_resources")]
    pub resources: MultiSharedResource,
}

#[allow(clippy::type_complexity)]
#[derive(Clone, Deserialize, MartianStruct)]
pub struct StageInputs {
    pub per_sample_metrics: TxHashMap<SampleAssignment, Option<JsonFile<TxHashMap<String, Value>>>>,
    pub library_metrics: Option<JsonFile<TxHashMap<String, Value>>>,
    pub sequencing_metrics: Option<SequencingMetricsFormat>,
    pub multi_config: MultiConfigCsvFile,
    pub multi_graph: JsonFile<CrMultiGraph>,
    pub multi_graph_svg: SvgFile,
    pub common_inputs: CommonInputs,
    pub count_inputs: Option<CountInputs>,
    pub tag_contaminant_info: Option<JsonFile<Value>>,
    pub sample_tsne_plots: TxHashMap<SampleAssignment, Option<JsonFile<SampleTsnePlots>>>,
    pub sample_barcode_rank_plots:
        TxHashMap<SampleAssignment, Option<JsonFile<TxHashMap<LegacyLibraryType, PlotlyChart>>>>,
    pub sample_treemap_plots: Option<
        TxHashMap<
            SampleAssignment,
            Option<JsonFile<TxHashMap<LegacyLibraryType, RawChartWithHelp>>>,
        >,
    >,
    pub barcode_rank_plots: Option<JsonFile<TxHashMap<LegacyLibraryType, PlotlyChart>>>,
    pub jibes_biplot_histogram: Option<JsonFile<Value>>,
    pub antibody_histograms: Option<JsonFile<RawChartWithHelp>>,
    pub sample_antibody_histograms:
        Option<TxHashMap<SampleAssignment, Option<JsonFile<RawChartWithHelp>>>>,
    pub antigen_histograms: Option<JsonFile<RawChartWithHelp>>,
    pub targeted_per_feature_metrics: Option<CsvFile<()>>,
    pub cmo_tsne_plot: Option<JsonFile<MultiplexingTsnePlots>>,
    pub vdj_t_contents: Option<VdjWsContentsFormat>,
    pub vdj_t_gd_contents: Option<VdjWsContentsFormat>,
    pub vdj_b_contents: Option<VdjWsContentsFormat>,
    pub target_set_name: Option<String>,
    pub antigen_vdj_metrics: Option<AntigenVdjMetricsFormat>,
    pub antigen_specificity: Option<CsvFile<AntigenSpecificityRow>>,
    pub feature_config: Option<FeatureConfig>,
    pub detected_probe_barcode_pairing: Option<DetectedProbeBarcodePairingFile>,
    pub no_preflight: bool,
}

#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct StageOutputs {
    #[mro_retain]
    pub web_summary_json: TxHashMap<SampleAssignment, JsonFile<()>>,
    pub metrics_summary_csv: TxHashMap<SampleAssignment, CsvFile<()>>,
}

struct LibWsBuilder {
    common_inputs: CommonInputs,
    multi_graph: CrMultiGraph,
    multi_config: MultiConfigCsv,
    lib_metrics: TxHashMap<String, Value>,
    barcode_rank_plots: TxHashMap<LegacyLibraryType, PlotlyChart>,
    sequencing_metrics: TxHashMap<LegacyLibraryType, SequencingMetricsTable>,
    count_inputs: Option<CountInputs>,
    jibes_biplot_histogram: Option<Value>,
    antibody_histograms: Option<RawChartWithHelp>,
    antigen_histograms: Option<RawChartWithHelp>,
    cmo_tsne_plot: Option<MultiplexingTsnePlots>,
    target_set_name: Option<String>,
    vdj_t_contents: Option<VdjWsContents>,
    vdj_t_gd_contents: Option<VdjWsContents>,
    vdj_b_contents: Option<VdjWsContents>,
    targeted_per_feature_metrics: Option<CsvFile<()>>,
    is_multiplexed: bool,
    specificity_controls: Option<SpecificityControls>,
    dropped_tags: Vec<String>,
    probe_barcodes_high_gem_overlap: Vec<String>,
    mismatched_probe_barcode_pairings: Option<MismatchedProbeBarcodePairings>,
}

impl LibWsBuilder {
    fn vdj_tab_names(&self) -> [Option<&'static str>; 3] {
        let mut tab_names = [None; 3];
        for (i, (name, content)) in [
            ("VDJ-T", &self.vdj_t_contents),
            ("VDJ-B", &self.vdj_b_contents),
            ("VDJ-T-GD", &self.vdj_t_gd_contents),
        ]
        .into_iter()
        .enumerate()
        {
            if content.is_some() {
                tab_names[i] = Some(name);
            }
        }
        tab_names
    }
    /// Return the targeting method.
    fn targeting_method(&self) -> Option<TargetingMethod> {
        self.count_inputs.as_ref().and_then(|x| x.targeting_method)
    }

    /// Return true if the targeting method is hybridization capture.
    fn is_hybrid_capture(&self) -> bool {
        self.targeting_method() == Some(TargetingMethod::HybridCapture)
    }

    /// Return true if the targeting method is RTL.
    fn is_rtl(&self) -> bool {
        self.targeting_method() == Some(TargetingMethod::TemplatedLigation)
    }

    /// Return true if multiplexed using CMO.
    fn is_cmo_multiplexed(&self) -> bool {
        self.is_multiplexed && !self.is_rtl()
    }

    /// Return true if multiplexed using RTL .
    fn is_rtl_multiplexed(&self) -> bool {
        self.is_multiplexed && self.is_rtl()
    }

    /// Return true if antigen specificity controls are specified.
    fn has_antigen_controls(&self) -> bool {
        self.specificity_controls.is_some()
    }

    /// Return the aligner.
    fn aligner(&self) -> AlignerParam {
        get_metric_string(&self.lib_metrics, "alignment_aligner")
            .unwrap()
            .unwrap()
            .parse()
            .unwrap()
    }

    /// Return the feature ref path input by the user
    fn get_feature_ref(&self) -> Option<String> {
        self.multi_config
            .feature
            .as_ref()
            .and_then(|x| x.reference_path.as_ref().map(|r| r.display().to_string()))
    }

    /// Return the cmo-set input by the user
    fn get_cmo_set(&self) -> Option<String> {
        self.multi_config
            .gene_expression
            .as_ref()
            .and_then(|x| x.cmo_set.as_ref().map(|r| r.display().to_string()))
    }

    fn build(&self, pipeline_version: String, context: &AlertContext) -> Result<LibraryWebSummary> {
        let mut lib_ws = LibraryWebSummary {
            header_info: LibraryHeaderInfo {
                run_id: self.common_inputs.sample_id.clone(),
                run_desc: self.common_inputs.sample_desc.clone(),
                pipeline_version,
            },
            vdj_t_tab: self
                .vdj_t_contents
                .clone()
                .map(|content| Tab::new(content.to_library_ws(), context)),
            vdj_t_gd_tab: self
                .vdj_t_gd_contents
                .clone()
                .map(|content| Tab::new(content.to_library_ws(), context)),
            vdj_b_tab: self
                .vdj_b_contents
                .clone()
                .map(|content| Tab::new(content.to_library_ws(), context)),
            ..Default::default()
        };

        for lib in &self.multi_graph.libraries {
            match lib.legacy_library_type() {
                LegacyLibraryType::GeneExpression => {
                    assert!(lib_ws.gex_tab.is_none());
                    lib_ws.gex_tab = Some(Tab::new(
                        self.build_gex_ws(&lib.physical_library_id)?,
                        context,
                    ))
                }
                LegacyLibraryType::AntibodyCapture => {
                    assert!(lib_ws.antibody_tab.is_none());
                    lib_ws.antibody_tab = Some(Tab::new(
                        self.build_antibody_or_antigen_ws(&lib.physical_library_id, true)?,
                        context,
                    ))
                }
                LegacyLibraryType::AntigenCapture => {
                    assert!(lib_ws.antigen_tab.is_none());
                    lib_ws.antigen_tab = Some(Tab::new(
                        self.build_antibody_or_antigen_ws(&lib.physical_library_id, false)?,
                        context,
                    ))
                }
                LegacyLibraryType::CrisprGuideCapture => {
                    assert!(lib_ws.crispr_tab.is_none());
                    lib_ws.crispr_tab = Some(Tab::new(
                        self.build_crispr_ws(&lib.physical_library_id)?,
                        context,
                    ))
                }
                LegacyLibraryType::Vdj => {}
                LegacyLibraryType::Custom => {
                    assert!(lib_ws.custom_feature_tab.is_none());
                    lib_ws.custom_feature_tab = Some(Tab::new(
                        self.build_custom_ws(&lib.physical_library_id)?,
                        context,
                    ))
                }
                LegacyLibraryType::Multiplexing => {
                    assert!(lib_ws.cmo_tab.is_none());
                    let mut cmo_tab: LibraryCmoWebSummary =
                        self.build_cmo_ws(&lib.physical_library_id)?;
                    // bubble up shared resources, currently just cmo counts
                    lib_ws.resources.extend(cmo_tab.resources.drain());
                    lib_ws.cmo_tab = Some(Tab::new(cmo_tab, context));
                }
                LegacyLibraryType::ATAC => unreachable!(),
            }
        }
        Ok(lib_ws)
    }
    fn lib_fraction_cell_partitions(&self, num_key: &str) -> Result<Option<CountAndPercent>> {
        let num = get_metric_usize(&self.lib_metrics, num_key)?;
        let denom = get_metric_usize(&self.lib_metrics, "total_cell_associated_partitions")?;
        if num.is_none() || denom.is_none() {
            return Ok(None);
        }
        Ok(Some(CountAndPercent(PercentMetric::from_parts(
            num.unwrap() as i64,
            denom.unwrap() as i64,
        ))))
    }

    fn get_library_cell_metrics_row(
        &self,
        physical_library_id: &str,
        mean_reads_per_cell_associated_partition: Option<FloatAsInt>,
    ) -> Result<LibraryCellMetricsRow> {
        let has_gex = self
            .multi_graph
            .has_legacy_library_type(LegacyLibraryType::GeneExpression);
        let has_antibody = self
            .multi_graph
            .has_legacy_library_type(LegacyLibraryType::AntibodyCapture);

        let cell_associated_partitions = if has_antibody && !has_gex {
            get_metric_usize(
                &self.lib_metrics,
                "ANTIBODY_filtered_bcs_transcriptome_union",
            )?
        } else {
            get_metric_usize(&self.lib_metrics, "filtered_bcs_transcriptome_union")?
        };

        let get_if_cmo_multiplexed = |key: &str| -> Result<Option<CountAndPercent>> {
            if self.is_cmo_multiplexed() {
                self.lib_fraction_cell_partitions(key)
            } else {
                Ok(None)
            }
        };

        let get_if_rtl_multiplexed = |key: &str| -> Result<Option<f64>> {
            if self.is_rtl_multiplexed() {
                get_metric_f64(&self.lib_metrics, key)
            } else {
                Ok(None)
            }
        };

        Ok(LibraryCellMetricsRow {
            physical_library_id: Some(physical_library_id.to_string()),
            cell_associated_partitions,
            mean_reads_per_cell_associated_partition,
            singlets_assigned_sample: get_if_cmo_multiplexed("total_singlets")?,
            partitions_with_no_cmos: get_if_cmo_multiplexed(
                "cell_associated_partitions_not_assigned_any_samples",
            )?,
            partitions_called_multiplets: get_if_cmo_multiplexed(
                "cell_associated_partitions_identified_as_multiplets",
            )?,
            fraction_cells_passing_high_occupancy_filtering: get_if_rtl_multiplexed(
                "rtl_multiplexing_fraction_cells_in_high_occupancy_gems",
            )?
            .map(|x| Percent::Float(1.0 - x)),
        })
    }

    /// Return the chemistry description.
    fn chemistry_description(&self) -> Option<String> {
        get_metric_string(&self.lib_metrics, "chemistry_description").unwrap()
    }

    /// Return the chemistry description and append "(manual)" if a manual chemistry is specified.
    fn chemistry_description_with_manual(&self) -> String {
        let chemistry_description = self.chemistry_description().unwrap();
        let Some(gex) = &self.multi_config.gene_expression else {
            return chemistry_description;
        };
        if gex.chemistry == ChemistryParam::Auto {
            chemistry_description
        } else {
            format!("{chemistry_description} (manual)")
        }
    }

    /// Return the chemistry name.
    fn chemistry(&self) -> Option<ChemistryName> {
        self.chemistry_description().map(|s| {
            if s.starts_with("custom") {
                ChemistryName::Custom
            } else {
                ChemistryName::from_description(&s).unwrap()
            }
        })
    }

    fn count_param_table(&self, library_type: LegacyLibraryType) -> Result<CountParametersTable> {
        let count_inputs = self.count_inputs.as_ref().unwrap();

        let (unspecified_probe_barcodes_detected, specified_probe_barcodes_missing) =
            self.identify_unexpected_or_missing_probe_barcodes(library_type);

        Ok(CountParametersTable {
            chemistry: self.chemistry_description_with_manual(),
            introns_included: count_inputs.include_introns,
            reference_path: count_inputs.reference_path.display().to_string(),
            transcriptome: format!(
                "{}-{}",
                get_metric_string(&self.lib_metrics, "reference_genomes")?.unwrap(),
                get_metric_string(&self.lib_metrics, "reference_version")?.unwrap()
            ),
            feature_ref_path: self.get_feature_ref(),
            cmo_set_path: self.get_cmo_set(),
            target_set_name: self.target_set_name.clone(),
            targeting_method: self.targeting_method(),
            filter_probes: count_inputs.filter_probes,
            num_genes_on_target: self
                .lib_metrics
                .get("num_genes_on_target")
                .map(|x| x.as_u64().unwrap() as usize),
            library_type,
            throughput: get_metric_string(&self.lib_metrics, "throughput_inferred")?,
            tenx_cmos: count_inputs.tenx_cmos,
            aligner: self.aligner(),
            antigen_negative_control: self.has_antigen_controls(),
            dropped_tags: self.dropped_tags.clone(),
            probe_barcodes_high_gem_overlap: self.probe_barcodes_high_gem_overlap.clone(),
            mismatched_probe_barcode_pairings: self.mismatched_probe_barcode_pairings.clone(),
            unspecified_probe_barcodes_detected,
            specified_probe_barcodes_missing,
        })
    }

    fn genomes(&self) -> Vec<String> {
        let genomes_str = get_metric_string(&self.lib_metrics, "reference_genomes")
            .unwrap()
            .unwrap();
        genomes_str
            .split(MULTI_GENOME_SEPARATOR)
            .map(String::from)
            .collect()
    }

    fn get_mapping_metrics_table(
        &self,
        physical_library_id: &str,
    ) -> Result<GexOrRtlLibraryMappingMetricsTable> {
        let metrics = &self.lib_metrics;
        if self.aligner() == AlignerParam::Hurtle {
            Ok(GexOrRtl::Rtl(RtlLibraryMappingMetricsTable(vec![
                RtlLibraryMappingMetricsRow {
                    physical_library_id: Some(physical_library_id.to_string()),
                    reads_in_library: get_metric_usize(metrics, "total_read_pairs")?,
                    reads_half_mapped_to_probe_set: get_metric_percent(
                        metrics,
                        "multi_transcriptome_half_mapped_reads_frac",
                    )?,
                    reads_split_mapped_to_probe_set: get_metric_percent(
                        metrics,
                        "multi_transcriptome_split_mapped_reads_frac",
                    )?,
                    reads_mapped_to_probe_set: get_metric_percent(
                        metrics,
                        "multi_transcriptome_mapped_reads_frac",
                    )?,
                    reads_confidently_mapped_to_probe_set: get_metric_percent(
                        metrics,
                        "multi_transcriptome_conf_mapped_reads_frac",
                    )?,
                    reads_confidently_mapped_to_filtered_probe_set: get_metric_percent(
                        metrics,
                        "multi_transcriptome_targeted_conf_mapped_reads_frac",
                    )?,
                },
            ])))
        } else {
            Ok(GexOrRtl::Gex(GexLibraryMappingMetricsTable(vec![
                GexLibraryMappingMetricsRow {
                    physical_library_id: Some(physical_library_id.to_string()),
                    reads_in_library: get_metric_usize(metrics, "total_read_pairs")?,
                    mapped_to_genome: get_metric_percent(
                        metrics,
                        "multi_genome_mapped_reads_frac",
                    )?,
                    confidently_mapped_to_genome: get_metric_percent(
                        metrics,
                        "multi_genome_conf_mapped_reads_frac",
                    )?,
                    confidently_mapped_to_transcriptome: get_metric_percent(
                        metrics,
                        "multi_transcriptome_conf_mapped_reads_frac",
                    )?,
                    confidently_mapped_to_targeted_transcriptome: get_metric_percent(
                        metrics,
                        "multi_transcriptome_targeted_conf_mapped_reads_frac",
                    )?,
                    confidently_mapped_to_intronic_regions: get_metric_percent(
                        metrics,
                        "multi_intronic_conf_mapped_reads_frac",
                    )?,
                    confidently_mapped_to_exonic_regions: get_metric_percent(
                        metrics,
                        "multi_exonic_conf_mapped_reads_frac",
                    )?,
                    confidently_mapped_to_intergenic_regions: get_metric_percent(
                        metrics,
                        "multi_intergenic_conf_mapped_reads_frac",
                    )?,
                    confidently_mapped_antisense: get_metric_percent(
                        metrics,
                        "multi_antisense_reads_frac",
                    )?,
                },
            ])))
        }
    }

    /// Return the probe barcode metrics table for the specified library type.
    fn build_rtl_probe_barcode_metrics_table(
        &self,
        library_type: LegacyLibraryType,
    ) -> Option<RtlProbeBarcodeMetricsTable> {
        if !self.is_rtl_multiplexed() {
            return None;
        }

        let expected_multiplexing_barcode_type = match library_type {
            LegacyLibraryType::GeneExpression => MultiplexingBarcodeType::RTL,
            LegacyLibraryType::AntibodyCapture => MultiplexingBarcodeType::Antibody,
            _ => unreachable!("unexpected library type {library_type}"),
        };

        let probe_barcode_id_to_sample_id_and_mapped_probe_barcode_ids: TxHashMap<_, _> = self
            .multi_graph
            .samples
            .iter()
            .flat_map(|sample| {
                sample
                    .fingerprints
                    .iter()
                    .map(|fingerprint| match fingerprint {
                        Fingerprint::Tagged {
                            tag_name,
                            cell_multiplexing_type: CellMultiplexingType::RTL,
                            translated_tag_names,
                            ..
                        } => (
                            tag_name.as_str(),
                            (sample.sample_id.as_str(), translated_tag_names.as_slice()),
                        ),
                        Fingerprint::Tagged {
                            cell_multiplexing_type: CellMultiplexingType::CMO,
                            ..
                        } => unreachable!(),
                        Fingerprint::Tagged {
                            cell_multiplexing_type: CellMultiplexingType::OH,
                            ..
                        } => unreachable!(),
                        Fingerprint::Untagged { .. } => unreachable!(),
                    })
            })
            .collect();

        let umi_per_probe_barcode: Vec<(&str, i64)> = if let Some(umi_per_probe_barcode) = self
            .lib_metrics
            .get(library_type.join("umi_per_probe_barcode").as_ref())
        {
            umi_per_probe_barcode
                .as_object()
                .unwrap()
                .iter()
                .map(|(id, n)| (id.as_str(), n.as_i64().unwrap()))
                .sorted()
                .collect()
        } else {
            Vec::default()
        };

        let umi_sum = umi_per_probe_barcode.iter().map(|(_id, n)| n).sum();

        let filtered_barcodes_per_probe_barcode: TxHashMap<&str, i64> = self.lib_metrics
            ["filtered_barcodes_per_probe_barcode"]
            .as_object()
            .unwrap()
            .iter()
            .map(|(id, n)| (id.as_str(), n.as_i64().unwrap()))
            .collect();

        let filtered_barcodes_sum = filtered_barcodes_per_probe_barcode
            .iter()
            .map(|(_id, &n)| n)
            .sum();

        Some(RtlProbeBarcodeMetricsTable(
            umi_per_probe_barcode
                .into_iter()
                .filter_map(|(probe_barcode_id, umi_count)| {
                    let (sample_id, mapped_barcode_ids) =
                        probe_barcode_id_to_sample_id_and_mapped_probe_barcode_ids
                            .get(probe_barcode_id)
                            .copied()
                            .unzip();
                    let barcode_ids: Vec<_> = std::iter::once(probe_barcode_id)
                        .chain(mapped_barcode_ids.into_iter().flatten().map(String::as_str))
                        .collect();

                    // If this probe barcode isn't of the expected type, and all
                    // of the mapped probe barcodes also aren't of the expected type,
                    // do not include this probe barcode in the report.
                    // See CELLRANGER-7501 for the root of this problem.
                    if !barcode_ids.iter().any(|x| {
                        categorize_multiplexing_barcode_id(x) == expected_multiplexing_barcode_type
                    }) {
                        return None;
                    }

                    let umi_per_probe_barcode =
                        CountAndPercent(PercentMetric::from_parts(umi_count, umi_sum));

                    let cells = *filtered_barcodes_per_probe_barcode
                        .get(probe_barcode_id)
                        .unwrap_or(&0);
                    let cells_per_probe_barcode = if sample_id.is_some() {
                        Some(CountAndPercent(PercentMetric::from_parts(
                            cells,
                            filtered_barcodes_sum,
                        )))
                    } else {
                        assert_eq!(cells, 0);
                        None
                    };

                    // Elide probe barcodes not assigned to a sample and with few UMI.
                    match sample_id {
                        Some(sample_id) => Some(RtlProbeBarcodeMetricsRow {
                            probe_barcode_id: Some(barcode_ids.join(PROBE_BARCODE_ID_GROUPING)),
                            sample_id: Some(sample_id.to_string()),
                            umi_per_probe_barcode: Some(umi_per_probe_barcode),
                            cells_per_probe_barcode,
                        }),
                        None if umi_per_probe_barcode.0.fraction().unwrap_or(0.0)
                            < UMI_PER_PROBE_BARCODE_BACKGROUND_THRESHOLD =>
                        {
                            None
                        }
                        None => Some(RtlProbeBarcodeMetricsRow {
                            probe_barcode_id: Some(probe_barcode_id.to_string()),
                            sample_id: None,
                            umi_per_probe_barcode: Some(umi_per_probe_barcode),
                            cells_per_probe_barcode,
                        }),
                    }
                })
                .sorted_by(|a, b| a.sample_id.cmp(&b.sample_id))
                .collect(),
        ))
    }

    /// Use the RTL probe metrics table to find unexpected or missing probe barcodes.
    /// Return a tuple of (unexpected_but_detected, specified_but_missing).
    fn identify_unexpected_or_missing_probe_barcodes(
        &self,
        library_type: LegacyLibraryType,
    ) -> (Vec<String>, Vec<String>) {
        let Some(bc_metrics_table) = self.build_rtl_probe_barcode_metrics_table(library_type) else {return Default::default()};
        let unspecified_probe_barcodes_detected = bc_metrics_table
            .0
            .iter()
            .filter_map(|row| {
                if row.sample_id.is_none() {
                    row.probe_barcode_id.clone()
                } else {
                    None
                }
            })
            .collect();

        let specified_probe_barcodes_missing = bc_metrics_table
            .0
            .iter()
            .filter_map(|row| {
                let below_threshold = row
                    .umi_per_probe_barcode
                    .and_then(|umi_per_probe_barcode| umi_per_probe_barcode.0.fraction())
                    .map_or(false, |frac| {
                        row.sample_id.is_some() && frac < UMI_PER_PROBE_BARCODE_BACKGROUND_THRESHOLD
                    });

                if below_threshold && row.sample_id.is_some() {
                    row.probe_barcode_id.clone()
                } else {
                    None
                }
            })
            .collect();
        (
            unspecified_probe_barcodes_detected,
            specified_probe_barcodes_missing,
        )
    }

    fn build_gex_ws(&self, physical_library_id: &str) -> Result<LibraryGexWebSummary> {
        let barcode_rank_plot = format_barcode_rank_plot(
            self.barcode_rank_plots
                .get(&LegacyLibraryType::GeneExpression)
                .expect("Gene expression barcode rank plot was missing."),
            "GEX",
        );

        let mean_reads_per_cell_associated_partition = get_metric_f64(
            &self.lib_metrics,
            "multi_transcriptome_total_raw_reads_per_filtered_bc",
        )?
        .map(FloatAsInt);

        // Verify that the gDNA metric is present before trying to add to the WS.
        let gdna_table = if matches!(
            self.lib_metrics.get("estimated_gdna_content"),
            None | Some(Value::Null)
        ) {
            None
        } else {
            assert!(self.is_rtl());
            Some(GdnaMetricsTable::from_metrics(&self.lib_metrics)?.into())
        };

        let log_rpu_threshold = get_metric_f64(&self.lib_metrics, "log_rpu_threshold")?;

        let targeted_table = if self.is_hybrid_capture() {
            let umi_filtering_active =
                get_metric_f64(&self.lib_metrics, "filtered_target_umi_count_threshold")?
                    .unwrap_or(0.0)
                    > 1.0;

            let filtered_target_umi_count_threshold = if umi_filtering_active {
                get_metric_f64(&self.lib_metrics, "filtered_target_umi_count_threshold")?
            } else {
                None
            };

            let filtered_target_umi_reads_frac = if umi_filtering_active {
                get_metric_percent(&self.lib_metrics, "filtered_target_umi_reads_frac")?
            } else {
                None
            };

            Some(
                GexLibraryTargetedEnrichmentMetricsTable(vec![
                    GexLibraryTargetedEnrichmentMetricsRow {
                        targeting_status: Some("Targeted".to_string()),
                        multi_frac_conf_transcriptomic_reads: get_metric_percent(
                            &self.lib_metrics,
                            "multi_frac_conf_transcriptomic_reads_on_target",
                        )?,
                        num_genes: get_metric_usize(&self.lib_metrics, "num_genes_on_target")?,
                        num_genes_quantifiable: get_metric_usize(
                            &self.lib_metrics,
                            "num_genes_quantifiable_on_target",
                        )?,
                        num_rpu_enriched_genes: get_metric_usize(
                            &self.lib_metrics,
                            "num_rpu_enriched_genes_on_target",
                        )?,
                        mean_reads_per_umi_per_gene_cells: get_metric_f64(
                            &self.lib_metrics,
                            "mean_reads_per_umi_per_gene_cells_on_target",
                        )?,
                        filtered_target_umi_count_threshold,
                        filtered_target_umi_reads_frac,
                    },
                    GexLibraryTargetedEnrichmentMetricsRow {
                        targeting_status: Some("Non-Targeted".to_string()),
                        multi_frac_conf_transcriptomic_reads: get_metric_percent(
                            &self.lib_metrics,
                            "multi_frac_conf_transcriptomic_reads_off_target",
                        )?,
                        num_genes: get_metric_usize(&self.lib_metrics, "num_genes_off_target")?,
                        num_genes_quantifiable: get_metric_usize(
                            &self.lib_metrics,
                            "num_genes_quantifiable_off_target",
                        )?,
                        num_rpu_enriched_genes: get_metric_usize(
                            &self.lib_metrics,
                            "num_rpu_enriched_genes_off_target",
                        )?,
                        mean_reads_per_umi_per_gene_cells: get_metric_f64(
                            &self.lib_metrics,
                            "mean_reads_per_umi_per_gene_cells_off_target",
                        )?,
                        filtered_target_umi_count_threshold: None,
                        filtered_target_umi_reads_frac: None,
                    },
                ])
                .into(),
            )
        } else {
            None
        };

        let targeted_alerts = if self.is_hybrid_capture() {
            Some(
                GexLibraryTargetedEnrichmentAlertsTable(vec![
                    GexLibraryTargetedEnrichmentAlertsRow {
                        frac_on_target_genes_enriched: get_metric_percent(
                            &self.lib_metrics,
                            "frac_on_target_genes_enriched",
                        )?,
                        frac_off_target_genes_enriched: get_metric_percent(
                            &self.lib_metrics,
                            "frac_off_target_genes_enriched",
                        )?,
                    },
                ])
                .into(),
            )
        } else {
            None
        };

        let targeted_plot = if self.is_hybrid_capture() {
            targeted_enrichment_plot(
                self.targeted_per_feature_metrics.as_deref(),
                log_rpu_threshold,
            )?
        } else {
            None
        };

        let targeted_sequencing_saturation = if self.is_hybrid_capture() {
            get_metric_percent(
                &self.lib_metrics,
                "multi_cdna_pcr_dupe_reads_frac_on_target",
            )?
        } else {
            None
        };

        let mean_targeted_reads_per_cell_associated_partition = if self.is_hybrid_capture() {
            get_metric_f64(&self.lib_metrics, "total_targeted_reads_per_filtered_bc")?
                .map(FloatAsInt)
        } else {
            None
        };

        let (valid_gem_barcodes, valid_probe_barcodes) = if self.is_rtl_multiplexed() {
            (
                get_metric_percent(&self.lib_metrics, "good_bc_in_gel_bead_frac")?,
                get_metric_percent(&self.lib_metrics, "good_bc_in_probe_frac")?,
            )
        } else {
            (None, None)
        };

        // Construct final websummary
        Ok(LibraryGexWebSummary {
            parameters_table: self.count_param_table(LegacyLibraryType::GeneExpression)?,
            cell_metrics_table: LibraryCellMetricsTable(vec![self.get_library_cell_metrics_row(
                physical_library_id,
                mean_reads_per_cell_associated_partition,
            )?])
            .into(),
            sequencing_metrics_table: self
                .sequencing_metrics
                .get(&LegacyLibraryType::GeneExpression)
                .expect("Sequencing metrics for Gene Expression not found.")
                .clone()
                .into(),
            mapping_metrics_table: self.get_mapping_metrics_table(physical_library_id)?.into(),
            physical_library_metrics_table: GexPhysicalLibraryMetricsTable(vec![
                GexPhysicalLibraryMetricsRow {
                    physical_library_id: Some(physical_library_id.to_string()),
                    number_of_reads: get_metric_usize(
                        &self.lib_metrics,
                        "total_read_pairs", // DOUBLE CHECK THIS METRIC
                    )?,
                    valid_barcodes: get_metric_percent(&self.lib_metrics, "good_bc_frac")?,
                    valid_gem_barcodes,
                    valid_probe_barcodes,
                    valid_umis: get_metric_percent(&self.lib_metrics, "good_umi_frac")?,
                    sequencing_saturation: get_metric_percent(
                        &self.lib_metrics,
                        "multi_cdna_pcr_dupe_reads_frac", // DOUBLE CHECK THIS METRIC exactly the same thing as sequencing saturation?
                    )?,
                    targeted_sequencing_saturation,
                    reads_in_cell_associated_partitions: get_metric_percent(
                        &self.lib_metrics,
                        "multi_filtered_bcs_conf_mapped_barcoded_reads_cum_frac",
                    )?,
                    mean_reads_per_cell_associated_partition,
                    mean_targeted_reads_per_cell_associated_partition,
                },
            ])
            .into(),
            rtl_probe_barcode_metrics_table: self
                .build_rtl_probe_barcode_metrics_table(LegacyLibraryType::GeneExpression)
                .map(Into::into),
            gdna_table,
            targeted_table,
            targeted_plot,
            targeted_alerts,
            sequencing_saturation_plot: library_sequencing_saturation_plot_from_metrics(
                &self.lib_metrics,
                self.is_hybrid_capture(),
            ),
            median_genes_per_cell_plot: library_median_genes_plot_from_metrics(
                &self.lib_metrics,
                TxHashSet::from_iter(self.genomes().into_iter()),
                PlotType::LibraryPlot,
                self.targeting_method(),
            ),
            barcode_rank_plot,
        })
    }

    fn build_antibody_or_antigen_ws(
        &self,
        physical_library_id: &str,
        // true for antibody, false for antigen
        is_antibody: bool,
    ) -> Result<LibraryAntibodyOrAntigenWebSummary> {
        let library_type = if is_antibody {
            LegacyLibraryType::AntibodyCapture
        } else {
            LegacyLibraryType::AntigenCapture
        };

        let barcode_rank_plot = format_barcode_rank_plot(
            self.barcode_rank_plots
                .get(&library_type)
                .unwrap_or_else(|| panic!("{} barcode rank plot was missing.", &library_type)),
            if is_antibody { "AB" } else { "AG" },
        );

        let mean_reads_per_cell_associated_partition =
            get_metric_f64(&self.lib_metrics, &library_type.join("reads_per_cell"))?
                .map(FloatAsInt);

        let mut metrics = self.lib_metrics.clone();
        metrics.insert("physical_library_id".into(), physical_library_id.into());
        let physical_library_metrics_table = if is_antibody {
            let mut metrics_row = AntibodyPhysicalLibraryMetricsRow::from_metrics(&metrics)?;
            // Only show probe/GEM barcode mapping metrics for RTL multiplexing.
            if !self.is_rtl_multiplexed() {
                metrics_row.valid_gem_barcodes = None;
                metrics_row.valid_probe_barcodes = None;
            }
            AntibodyOrAntigen::Antibody(AntibodyPhysicalLibraryMetricsTable(vec![metrics_row]))
        } else {
            AntibodyOrAntigen::Antigen(AntigenPhysicalLibraryMetricsTable::from_metrics(&metrics)?)
        };

        let mapping_metrics_table = if is_antibody {
            Some(AntibodyLibraryMappingMetricsTable(vec![
                AntibodyLibraryMappingMetricsRow::from_metrics(&metrics)?,
            ]))
        } else {
            None
        };

        // Only show the antibody histogram at the library level for CMO.
        let feature_histogram = match (is_antibody, self.is_cmo_multiplexed()) {
            (true, true) => self.antibody_histograms.clone(),
            (true, false) => None,
            (false, _) => self.antigen_histograms.clone(),
        };

        let mut library_cell_metrics_row = self.get_library_cell_metrics_row(
            physical_library_id,
            mean_reads_per_cell_associated_partition,
        )?;

        // Only show high occupancy GEM filtering metric if there's no GEX data,
        // because if there is GEX data then the antibody data was not used to
        // create this filter.
        if self
            .multi_graph
            .has_legacy_library_type(LegacyLibraryType::GeneExpression)
        {
            library_cell_metrics_row.fraction_cells_passing_high_occupancy_filtering = None;
        }

        Ok(LibraryAntibodyOrAntigenWebSummary {
            parameters_table: self.count_param_table(library_type)?,
            cell_metrics_table: LibraryCellMetricsTable(vec![library_cell_metrics_row]).into(),
            sequencing_metrics_table: self
                .sequencing_metrics
                .get(&library_type)
                .expect("Sequencing metrics for Antibody or Antigen Capture not found.")
                .clone()
                .into(),
            mapping_metrics_table: mapping_metrics_table.map(Into::into),
            physical_library_metrics_table: physical_library_metrics_table.into(),
            rtl_probe_barcode_metrics_table: self
                .build_rtl_probe_barcode_metrics_table(library_type)
                .map(Into::into),
            barcode_rank_plot,
            feature_histogram,
        })
    }

    fn build_crispr_ws(&self, physical_library_id: &str) -> Result<LibraryCrisprWebSummary> {
        let mean_reads_per_cell_associated_partition =
            get_metric_f64(&self.lib_metrics, "CRISPR_reads_per_cell")?.map(FloatAsInt);
        Ok(LibraryCrisprWebSummary {
            parameters_table: self.count_param_table(LegacyLibraryType::CrisprGuideCapture)?,
            cell_metrics_table: LibraryCellMetricsTable(vec![self.get_library_cell_metrics_row(
                physical_library_id,
                mean_reads_per_cell_associated_partition,
            )?])
            .into(),
            sequencing_metrics_table: self
                .sequencing_metrics
                .get(&LegacyLibraryType::CrisprGuideCapture)
                .expect("Sequencing metrics for CRISPR not found.")
                .clone()
                .into(),
            physical_library_metrics_table: CrisprPhysicalLibraryMetricsTable(vec![
                CrisprPhysicalLibraryMetricsRow {
                    physical_library_id: Some(physical_library_id.to_string()),
                    number_of_reads: get_metric_usize(
                        &self.lib_metrics,
                        "CRISPR_total_read_pairs", // DOUBLE CHECK THIS METRIC
                    )?,
                    valid_barcodes: get_metric_percent(&self.lib_metrics, "CRISPR_good_bc_frac")?,
                    valid_umis: get_metric_percent(&self.lib_metrics, "CRISPR_good_umi_frac")?,
                    sequencing_saturation: get_metric_percent(
                        &self.lib_metrics,
                        "CRISPR_multi_cdna_pcr_dupe_reads_frac",
                    )?,
                    reads_in_cell_associated_partitions: get_metric_percent(
                        &self.lib_metrics,
                        "CRISPR_feature_reads_in_cells",
                    )?,
                    mean_reads_per_cell_associated_partition,
                    fraction_reads_with_putative_protospacer: get_metric_percent(
                        &self.lib_metrics,
                        "CRISPR_feature_bc_extracted_frac",
                    )?,
                    fraction_guide_reads: get_metric_percent(
                        &self.lib_metrics,
                        "CRISPR_recognized_feature_bc_frac",
                    )?,
                    fraction_guide_reads_usable: get_metric_percent(
                        &self.lib_metrics,
                        "CRISPR_frac_feature_reads_usable",
                    )?,
                    fraction_protospacer_not_recognized: get_metric_percent(
                        &self.lib_metrics,
                        "CRISPR_unrecognized_feature_bc_frac", // DOUBLE CHECK THIS METRIC
                    )?,
                },
            ])
            .into(),
            barcode_rank_plot: format_barcode_rank_plot(
                self.barcode_rank_plots
                    .get(&LegacyLibraryType::CrisprGuideCapture)
                    .expect("CRISPR Guide Capture barcode rank plot was missing."),
                "CRISPR",
            ),
        })
    }

    fn build_custom_ws(&self, physical_library_id: &str) -> Result<LibraryCustomFeatureWebSummary> {
        let mean_reads_per_cell_associated_partition =
            get_metric_f64(&self.lib_metrics, "Custom_reads_per_cell")?.map(FloatAsInt);
        Ok(LibraryCustomFeatureWebSummary {
            parameters_table: self.count_param_table(LegacyLibraryType::Custom)?,
            cell_metrics_table: LibraryCellMetricsTable(vec![self.get_library_cell_metrics_row(
                physical_library_id,
                mean_reads_per_cell_associated_partition,
            )?])
            .into(),
            sequencing_metrics_table: self
                .sequencing_metrics
                .get(&LegacyLibraryType::Custom)
                .expect("Sequencing metrics for Custom Feature not found.")
                .clone()
                .into(),
            physical_library_metrics_table: CustomFeaturePhysicalLibraryMetricsTable(vec![
                CustomFeaturePhysicalLibraryMetricsRow {
                    physical_library_id: Some(physical_library_id.to_string()),
                    number_of_reads: get_metric_usize(
                        &self.lib_metrics,
                        "Custom_total_read_pairs",
                    )?,
                    valid_barcodes: get_metric_percent(&self.lib_metrics, "Custom_good_bc_frac")?,
                    valid_umis: get_metric_percent(&self.lib_metrics, "Custom_good_umi_frac")?,
                    sequencing_saturation: get_metric_percent(
                        &self.lib_metrics,
                        "Custom_multi_cdna_pcr_dupe_reads_frac",
                    )?,
                    reads_in_cell_associated_partitions: get_metric_percent(
                        &self.lib_metrics,
                        "Custom_feature_reads_in_cells",
                    )?,
                    mean_reads_per_cell_associated_partition,
                    fraction_feature_reads: get_metric_percent(
                        &self.lib_metrics,
                        "Custom_recognized_feature_bc_frac",
                    )?,
                    fraction_feature_reads_usable: get_metric_percent(
                        &self.lib_metrics,
                        "Custom_frac_feature_reads_usable",
                    )?,
                    fraction_unknown_feature: get_metric_percent(
                        &self.lib_metrics,
                        "Custom_unrecognized_feature_bc_frac", // DOUBLE CHECK THIS METRIC
                    )?,
                },
            ])
            .into(),
            barcode_rank_plot: format_barcode_rank_plot(
                self.barcode_rank_plots
                    .get(&LegacyLibraryType::Custom)
                    .expect("Custom feature barcode rank plot was missing."),
                "Custom",
            ),
        })
    }

    fn build_cmo_ws(&self, physical_library_id: &str) -> Result<LibraryCmoWebSummary> {
        let (jibes_biplot, jibes_histogram, resources): (
            Option<RawChartWithHelp>,
            Option<RawChartWithHelp>,
            TxHashMap<String, Value>,
        ) = match self.jibes_biplot_histogram.clone() {
            Some(jibes_biplot) => {
                let jibes_biplot_histogram: JibesBiplotHistogramData =
                    serde_json::from_value(jibes_biplot)?;
                (
                    Some(format_jibes_biplots(&jibes_biplot_histogram.biplot)),
                    Some(format_histogram(&jibes_biplot_histogram.histogram, "CMO")),
                    jibes_biplot_histogram.resources.clone(),
                )
            }
            None => (None, None, TxHashMap::default()),
        };

        Ok(LibraryCmoWebSummary {
            parameters_table: self.count_param_table(LegacyLibraryType::Multiplexing)?,
            multiplexing_metrics_table: MultiplexingLibraryCellMetricsTable(vec![
                MultiplexingLibraryCellMetricsRow {
                    cell_associated_partitions: get_metric_usize(
                        &self.lib_metrics,
                        "multi_filtered_bcs",
                    )?,
                    samples_assigned_at_least_one_singlet: get_metric_usize(
                        &self.lib_metrics,
                        "samples_with_any_singlets",
                    )?,
                    singlets_assigned_to_sample: self
                        .lib_fraction_cell_partitions("total_singlets")?,
                    singlet_capture_ratio: get_metric_f64(
                        &self.lib_metrics,
                        "MULTIPLEXING_sc_rec_efficiency_jibes",
                    )?,
                    cell_associated_partitions_identified_as_multiplet: self
                        .lib_fraction_cell_partitions(
                            "cell_associated_partitions_identified_as_multiplets",
                        )?,
                    median_cmo_umis_per_singlet: get_metric_round_float_to_int(
                        &self.lib_metrics,
                        "MULTIPLEXING_median_cmo_umis_per_singlet",
                    )?,
                },
            ])
            .into(),
            sample_assignments_table: MultiplexingSampleAssignmentsTable(vec![
                MultiplexingSampleAssignmentsRow {
                    physical_library_id: Some(physical_library_id.to_string()),
                    cell_associated_partitions: get_metric_usize(
                        &self.lib_metrics,
                        "multi_filtered_bcs",
                    )?,
                    mean_reads_per_cell: get_metric_f64(
                        &self.lib_metrics,
                        "MULTIPLEXING_multi_total_raw_reads_per_filtered_bc",
                    )?
                    .map(FloatAsInt),
                    samples_assigned_at_least_one_singlet: get_metric_usize(
                        &self.lib_metrics,
                        "samples_with_any_singlets",
                    )?,
                    singlets_assigned_to_a_sample: self
                        .lib_fraction_cell_partitions("total_singlets")?,
                    cell_associated_partitions_identified_as_multiplets: self
                        .lib_fraction_cell_partitions(
                            "cell_associated_partitions_identified_as_multiplets",
                        )?,
                    cell_associated_partitions_not_assigned_any_cmos: self
                        .lib_fraction_cell_partitions(
                            "cell_associated_partitions_not_assigned_any_samples",
                        )?,
                    median_cmo_umis_per_cell_associated_partition: get_metric_round_float_to_int(
                        &self.lib_metrics,
                        "MULTIPLEXING_multi_filtered_bcs_median_counts",
                    )?,
                },
            ])
            .into(),
            sequencing_metrics_table: self
                .sequencing_metrics
                .get(&LegacyLibraryType::Multiplexing)
                .expect("Sequencing metrics for Multiplexing not found.")
                .clone()
                .into(),
            physical_library_metrics_table: MultiplexingPhysicalLibraryMetricsTable(vec![
                MultiplexingPhysicalLibraryMetricsRow {
                    physical_library_id: Some(physical_library_id.to_string()),
                    number_of_reads: get_metric_usize(
                        &self.lib_metrics,
                        "MULTIPLEXING_total_read_pairs", // DOUBLE CHECK THIS METRIC
                    )?,
                    valid_barcodes: get_metric_percent(
                        &self.lib_metrics,
                        "MULTIPLEXING_good_bc_frac",
                    )?,
                    valid_umis: get_metric_percent(
                        &self.lib_metrics,
                        "MULTIPLEXING_good_umi_frac",
                    )?,
                    sequencing_saturation: get_metric_percent(
                        &self.lib_metrics,
                        "MULTIPLEXING_multi_cdna_pcr_dupe_reads_frac",
                    )?,
                    reads_in_cell_associated_partitions: get_metric_percent(
                        &self.lib_metrics,
                        "MULTIPLEXING_feature_reads_in_cells",
                    )?,
                    mean_reads_per_cell_associated_partition: get_metric_f64(
                        &self.lib_metrics,
                        "MULTIPLEXING_multi_total_raw_reads_per_filtered_bc",
                    )?
                    .map(FloatAsInt),
                    fraction_cmo_reads: get_metric_percent(
                        &self.lib_metrics,
                        "MULTIPLEXING_recognized_feature_bc_frac",
                    )?,
                    fraction_cmo_reads_usable: get_metric_percent(
                        &self.lib_metrics,
                        "MULTIPLEXING_frac_feature_reads_usable",
                    )?,
                    fraction_unknown_cmo: get_metric_percent(
                        &self.lib_metrics,
                        "MULTIPLEXING_unrecognized_feature_bc_frac", // DOUBLE CHECK THIS METRIC
                    )?,
                    fraction_reads_from_multiplets: get_metric_percent(
                        &self.lib_metrics,
                        "MULTIPLEXING_frac_reads_from_multiplets",
                    )?,
                },
            ])
            .into(),
            cmo_metrics_table: self.build_per_cmo_metrics_table()?.into(),
            barcode_rank_plot: format_barcode_rank_plot(
                self.barcode_rank_plots
                    .get(&LegacyLibraryType::Multiplexing)
                    .expect("CMO barcode rank plot was missing."),
                "CMO",
            ),
            jibes_biplot,
            jibes_histogram,
            cmo_umi_tsne_plot: self.cmo_tsne_plot.as_ref().map(|x| {
                format_umi_on_tsne_plot(
                    &x.cmo_umi_tsne_plot,
                    "Multiplexing Capture",
                    "t-SNE Projection of Cells Colored by UMI Counts",
                )
            }),
            cmo_tags_tsne_plot: self
                .cmo_tsne_plot
                .as_ref()
                .map(|x| format_tags_on_tsne_plot(&x.cmo_tags_tsne_plot)),
            resources,
        })
    }

    fn build_per_cmo_metrics_table(&self) -> Result<MultiplexingCmoMetricsTable> {
        let mut gem_wells = TxHashSet::default();
        let mut rows = Vec::new();
        for sample in &self.multi_graph.samples {
            rows.reserve(sample.fingerprints.len());
            for fingerprint in &sample.fingerprints {
                match fingerprint {
                    Fingerprint::Tagged {
                        gem_well, tag_name, ..
                    } => {
                        gem_wells.insert(*gem_well);
                        // assume input is single-gem-well
                        if gem_wells.len() > 1 {
                            bail!(
                                "Multiple gem well data passed into per-tag metrics table generation, which is not yet supported."
                            );
                        }
                        let singlets_metric_name = format!("tag_{tag_name}_number_of_singlets");
                        let reads_in_cell_associated_partitions_metric_name =
                            format!("tag_{tag_name}_frac_reads_in_cells");
                        let snr_metric_name = format!("snr_{tag_name}_jibes");
                        let total_singlets =
                            get_metric_usize(&self.lib_metrics, "total_singlets")?.unwrap() as i64;

                        let singlets_assigned_to_cmo =
                            get_metric_usize(&self.lib_metrics, &singlets_metric_name)?.map(
                                |singlets_assigned_to_cmo_count| {
                                    Percent::Metric(PercentMetric::from_parts(
                                        singlets_assigned_to_cmo_count as i64,
                                        total_singlets,
                                    ))
                                },
                            );
                        let reads_in_cell_associated_partitions = get_metric_percent(
                            &self.lib_metrics,
                            &reads_in_cell_associated_partitions_metric_name,
                        )?;

                        let cmo_signal_to_background_ratio =
                            get_metric_f64(&self.lib_metrics, &snr_metric_name)?;
                        rows.push(MultiplexingCmoMetricsRow {
                            gem_well_cmo: Some(tag_name.to_string()),
                            reads_in_cell_associated_partitions,
                            singlets_assigned_to_cmo,
                            cmo_signal_to_background_ratio,
                        });
                    }
                    Fingerprint::Untagged { .. } => {
                        bail!("Unable to build CMO metrics table involving untagged fingerprint.");
                    }
                }
            }
        }
        sort_cmo_rows(&mut rows);
        Ok(MultiplexingCmoMetricsTable(rows))
    }
}

// Helper function to sort rows by their CMO name
fn sort_cmo_rows(rows: &mut [MultiplexingCmoMetricsRow]) {
    rows.sort_by(|a, b| match (&a.gem_well_cmo, &b.gem_well_cmo) {
        (Some(ref _x), None) => Ordering::Greater,
        (None, Some(ref _y)) => Ordering::Less,
        (Some(ref x), Some(ref y)) => x.cmp(y),
        (None, None) => Ordering::Equal,
    })
}

// describes the various TSNE plots a single sample may have.
#[derive(Serialize, Deserialize, Clone, Default)]
pub struct SampleTsnePlots {
    // full TSNE/clustering/diffexp plots for gene expression run
    gex_diffexp_clustering_plots: Value,
    // full TSNE/clustering/diffexp plots for antibody-only case only
    antibody_diffexp_clustering_plots: Value,
    // the normal colored-umi TSNEs that feature barcode libraries get otherwise
    crispr_umi_on_tsne: Value,
    antibody_umi_on_tsne: Value,
    custom_umi_on_tsne: Value,
}

// describes the various TSNE plots a multiplexing experiment may have.
#[derive(Serialize, Deserialize, Clone, Default)]
pub struct MultiplexingTsnePlots {
    cmo_umi_tsne_plot: Value,
    cmo_tags_tsne_plot: Value,
}

struct MultiWsBuilder {
    lib_ws_builder: LibWsBuilder,
    per_sample_metrics: TxHashMap<SampleAssignment, TxHashMap<String, Value>>,
    sample_barcode_rank_plots:
        TxHashMap<SampleAssignment, TxHashMap<LegacyLibraryType, PlotlyChart>>,
    sample_treemap_plots:
        Option<TxHashMap<SampleAssignment, TxHashMap<LegacyLibraryType, RawChartWithHelp>>>,
    sample_tsne_plots: TxHashMap<SampleAssignment, SampleTsnePlots>,
    sample_antibody_histograms: Option<TxHashMap<SampleAssignment, RawChartWithHelp>>,
    svg_str: String,
    csv_str: String,
    diagnostics: MultiDiagnostics,
    pipeline_version: String,
    context: AlertContext,
    multiplexing_method: Option<CellMultiplexingType>,
    targeting_method: Option<TargetingMethod>,
    antigen_vdj_metrics: Option<AntigenVdjMetrics>,
    clonotype_clustermap: Option<RawChartWithHelp>,
}

impl MultiWsBuilder {
    /// Return true if the targeting method is hybridization capture.
    fn is_hybrid_capture(&self) -> bool {
        self.targeting_method == Some(TargetingMethod::HybridCapture)
    }

    /// Return true if multiplexed using CMO.
    fn is_cmo_multiplexed(&self) -> bool {
        self.multiplexing_method == Some(CellMultiplexingType::CMO)
    }

    /// Return true if multiplexed using RTL .
    fn is_rtl_multiplexed(&self) -> bool {
        self.multiplexing_method == Some(CellMultiplexingType::RTL)
    }

    /// Return true if multiplexed using OH.
    fn is_overhang_multiplexed(&self) -> bool {
        self.multiplexing_method == Some(CellMultiplexingType::OH)
    }

    /// Return the multi experimental design graph.
    fn multi_graph(&self) -> &CrMultiGraph {
        &self.lib_ws_builder.multi_graph
    }

    fn build(self) -> Result<TxHashMap<SampleAssignment, MultiWebSummary>> {
        let mut result = TxHashMap::default();

        let mut library_websummary = self
            .lib_ws_builder
            .build(self.pipeline_version.clone(), &self.context)?;

        let multi_graph = self.multi_graph();

        let mut resources = TxHashMap::default();
        resources.extend(library_websummary.resources.drain());
        // iterate over samples
        // read in their respective metrics JSONs
        // populate the sample web summary JSON for each one
        for (i, sample) in multi_graph.samples.iter().enumerate() {
            let mut sample_ws = SampleWebSummary {
                header_info: SampleHeaderInfo {
                    sample_id: sample.sample_id.clone(),
                    sample_desc: sample.description.clone(),
                    pipeline_version: self.pipeline_version.clone(),
                },
                vdj_t_tab: self
                    .lib_ws_builder
                    .vdj_t_contents
                    .clone()
                    .map(|content| Tab::new(content.to_sample_ws(), &self.context)),
                vdj_b_tab: self
                    .lib_ws_builder
                    .vdj_b_contents
                    .clone()
                    .map(|content| Tab::new(content.to_sample_ws(), &self.context)),
                vdj_t_gd_tab: self
                    .lib_ws_builder
                    .vdj_t_gd_contents
                    .clone()
                    .map(|content| Tab::new(content.to_sample_ws(), &self.context)),
                ..Default::default()
            };
            let sample_assignment = SampleAssignment::Assigned(sample.sample_id.clone());
            for lib in &multi_graph.libraries {
                match lib.legacy_library_type() {
                    LegacyLibraryType::GeneExpression => {
                        assert!(sample_ws.gex_tab.is_none());

                        sample_ws.gex_tab = Some(Tab::new(
                            self.build_gex_ws(&sample_assignment, &lib.physical_library_id)?,
                            &self.context,
                        ))
                    }
                    LegacyLibraryType::AntibodyCapture => {
                        assert!(sample_ws.antibody_tab.is_none());
                        // For FLEX, only show the antibody tab if we actually
                        // have at least one antibody probe barcode for this sample.
                        if sample.cell_multiplexing_type() == Some(CellMultiplexingType::RTL)
                            && !sample_is_associated_with_flex_library_type(
                                sample,
                                LegacyLibraryType::AntibodyCapture,
                            )
                        {
                            continue;
                        }
                        sample_ws.antibody_tab = Some(Tab::new(
                            self.build_antibody_ws(&sample_assignment, &lib.physical_library_id)?,
                            &self.context,
                        ))
                    }
                    LegacyLibraryType::AntigenCapture => {
                        assert!(sample_ws.antigen_tab.is_none());
                        sample_ws.antigen_tab = Some(Tab::new(
                            self.build_antigen_ws(&sample_assignment)?,
                            &self.context,
                        ))
                    }
                    LegacyLibraryType::CrisprGuideCapture => {
                        assert!(sample_ws.crispr_tab.is_none());
                        sample_ws.crispr_tab = Some(Tab::new(
                            self.build_crispr_ws(&sample_assignment, &lib.physical_library_id)?,
                            &self.context,
                        ))
                    }
                    LegacyLibraryType::Vdj => {}
                    LegacyLibraryType::Custom => {
                        assert!(sample_ws.custom_feature_tab.is_none());
                        sample_ws.custom_feature_tab = Some(Tab::new(
                            self.build_custom_ws(&sample_assignment, &lib.physical_library_id)?,
                            &self.context,
                        ))
                    }
                    LegacyLibraryType::Multiplexing => {}
                    LegacyLibraryType::ATAC => unreachable!(),
                }
            }
            result.insert(
                sample_assignment,
                MultiWebSummary {
                    sample: WsSample::multi(sample.sample_id.clone(), sample.description.clone()),
                    data: MultiWebSummaryData {
                        library_websummary: library_websummary.clone(),
                        sample_websummary: sample_ws,
                        experimental_design: ExperimentalDesign {
                            svg: SvgGraph::new(
                                self.svg_str.clone(),
                                format!("sample_{}", i + 1),
                                self.multiplexing_method,
                            ),
                            csv: self.csv_str.clone(),
                        },
                    },
                    diagnostics: self.diagnostics.clone(),
                    resources: resources.clone(),
                },
            );
        }
        Ok(result)
    }

    fn sample_cell_metrics_table(
        &self,
        metrics: &TxHashMap<String, Value>,
        physical_library_id: &str,
    ) -> Result<SampleCellMetricsTable> {
        let fraction_cell_partitions = |num_key| -> Result<Option<CountAndPercent>> {
            let num = get_metric_usize(metrics, num_key)?;
            let denom = get_metric_usize(metrics, "total_cell_associated_partitions")?;
            if num.is_none() || denom.is_none() {
                return Ok(None);
            }
            Ok(Some(CountAndPercent(PercentMetric::from_parts(
                num.unwrap() as i64,
                denom.unwrap() as i64,
            ))))
        };
        let physical_library_id_str = Some(physical_library_id.to_string());
        let singlets_assigned_to_this_sample =
            fraction_cell_partitions("singlets_assigned_to_this_sample")?;
        let singlets_assigned_to_other_samples =
            fraction_cell_partitions("singlets_assigned_to_other_samples")?;
        if self.is_rtl_multiplexed() | self.is_overhang_multiplexed() {
            Ok(GexOrRtl::Rtl(RtlSampleCellMetricsTable(vec![
                RtlSampleCellMetricsRow {
                    physical_library_id: physical_library_id_str,
                    singlets_assigned_to_this_sample,
                    singlets_assigned_to_other_samples,
                },
            ])))
        } else {
            Ok(GexOrRtl::Gex(GexSampleCellMetricsTable(vec![
                GexSampleCellMetricsRow {
                    physical_library_id: physical_library_id_str,
                    singlets_assigned_to_this_sample,
                    singlets_assigned_to_other_samples,
                    cell_associated_partitions_not_assigned_any_samples: fraction_cell_partitions(
                        "cell_associated_partitions_not_assigned_any_samples",
                    )?,
                    cell_associated_partitions_identified_as_multiplets: fraction_cell_partitions(
                        "cell_associated_partitions_identified_as_multiplets",
                    )?,
                },
            ])))
        }
    }

    fn get_mapping_metrics_table(
        &self,
        metrics: &TxHashMap<String, Value>,
    ) -> Result<GexOrRtlSampleMappingMetricsTable> {
        if self.lib_ws_builder.aligner() == AlignerParam::Hurtle {
            Ok(GexOrRtl::Rtl(RtlSampleMappingMetricsTable(vec![
                RtlSampleMappingMetricsRow {
                    reads_from_cells_assigned_to_sample: get_metric_usize(
                        metrics,
                        "total_read_pairs_in_filtered_barcodes",
                    )?,
                    reads_half_mapped_to_probe_set: get_metric_percent(
                        metrics,
                        "multi_transcriptome_half_mapped_reads_frac_in_filtered_barcodes",
                    )?,
                    reads_split_mapped_to_probe_set: get_metric_percent(
                        metrics,
                        "multi_transcriptome_split_mapped_reads_frac_in_filtered_barcodes",
                    )?,
                    reads_mapped_to_probe_set: get_metric_percent(
                        metrics,
                        "multi_transcriptome_mapped_reads_frac_in_filtered_barcodes",
                    )?,
                    reads_confidently_mapped_to_probe_set: get_metric_percent(
                        metrics,
                        "multi_transcriptome_conf_mapped_reads_frac_in_filtered_barcodes",
                    )?,
                    reads_confidently_mapped_to_filtered_probe_set: get_metric_percent(
                        metrics,
                        "multi_transcriptome_targeted_conf_mapped_reads_frac_in_filtered_barcodes",
                    )?,
                },
            ])))
        } else {
            Ok(GexOrRtl::Gex(GexSampleMappingMetricsTable(vec![
                GexSampleMappingMetricsRow {
                    reads_from_cells_assigned_to_sample: get_metric_usize(
                        metrics,
                        "total_read_pairs_in_filtered_barcodes",
                    )?,
                    mapped_to_genome: get_metric_percent(
                        metrics,
                        "multi_genome_mapped_reads_frac_in_filtered_barcodes",
                    )?,
                    confidently_mapped_to_genome: get_metric_percent(
                        metrics,
                        "multi_genome_conf_mapped_reads_frac_in_filtered_barcodes",
                    )?,
                    confidently_mapped_to_transcriptome: get_metric_percent(
                        metrics,
                        "multi_transcriptome_conf_mapped_reads_frac_in_filtered_barcodes",
                    )?,
                    confidently_mapped_to_targeted_transcriptome: get_metric_percent(
                        metrics,
                        "multi_transcriptome_targeted_conf_mapped_reads_frac_in_filtered_barcodes",
                    )?,
                    confidently_mapped_to_intronic_regions: get_metric_percent(
                        metrics,
                        "multi_intronic_conf_mapped_reads_frac_in_filtered_barcodes",
                    )?,
                    confidently_mapped_to_exonic_regions: get_metric_percent(
                        metrics,
                        "multi_exonic_conf_mapped_reads_frac_in_filtered_barcodes",
                    )?,
                    confidently_mapped_to_intergenic_regions: get_metric_percent(
                        metrics,
                        "multi_intergenic_conf_mapped_reads_frac_in_filtered_barcodes",
                    )?,
                    confidently_mapped_antisense: get_metric_percent(
                        metrics,
                        "multi_antisense_reads_frac_in_filtered_barcodes",
                    )?,
                },
            ])))
        }
    }

    fn build_gex_ws(
        &self,
        sample_assignment: &SampleAssignment,
        physical_library_id: &str,
    ) -> Result<SampleGexWebSummary> {
        let genomes = self.lib_ws_builder.genomes();
        let metrics = self
            .per_sample_metrics
            .get(sample_assignment)
            .expect("Sample metrics file not found.");

        let [median_genes_per_cell_on_target, num_genes_detected_on_target, median_umis_per_cell_on_target] =
            if self.is_hybrid_capture() && genomes.len() == 1 {
                [
                    "median_genes_per_cell_on_target",
                    "num_genes_detected_on_target",
                    "median_umis_per_cell_on_target",
                ]
                .map(|key| get_metric_round_float_to_int(metrics, key).unwrap())
            } else {
                [None, None, None]
            };

        let mut cell_hero_metrics_rows = Vec::with_capacity(genomes.len());
        for genome in &genomes {
            // use median RAW READS metric if this is single genome, but use per-genome mapped reads if it's genome specific
            let median_reads_per_singlet = if genomes.len() == 1 {
                get_metric_round_float_to_int(metrics, "median_total_reads_per_singlet")?
            } else {
                get_metric_round_float_to_int(
                    metrics,
                    &format!("{genome}_median_reads_per_singlet"),
                )?
            };
            let median_reads_per_cell_on_target = if self.is_hybrid_capture() && genomes.len() == 1
            {
                get_metric_round_float_to_int(
                    metrics,
                    &format!("{genome}_median_reads_per_singlet_ontarget"),
                )?
            } else {
                None
            };

            let [median_genes_per_singlet, total_genes_detected, median_umi_per_singlet] =
                match self.targeting_method {
                    Some(TargetingMethod::TemplatedLigation) => {
                        // Use the filtered probe set metrics.
                        [
                            "median_genes_per_cell_on_target",
                            "num_genes_detected_on_target",
                            "median_umis_per_cell_on_target",
                        ]
                        .map(|key| get_metric_round_float_to_int(metrics, key).unwrap())
                    }
                    Some(TargetingMethod::HybridCapture) if genomes.len() == 1 => {
                        // Hide these metrics for targeting with hybridization-capture.
                        [None, None, None]
                    }
                    _ => {
                        // Non-targeted metrics.
                        [
                            "filtered_bcs_median_unique_genes_detected",
                            "filtered_bcs_total_unique_genes_detected",
                            "filtered_bcs_median_counts",
                        ]
                        .map(|key| {
                            get_metric_round_float_to_int(metrics, &format!("{genome}_{key}"))
                                .unwrap()
                        })
                    }
                };

            let confidently_mapped_reads_in_cells = if self.is_cmo_multiplexed() {
                // CMO multiplexing does not demultiplex reads outside of cells.
                None
            } else {
                get_metric_percent(
                    metrics,
                    &format!("{genome}_filtered_bcs_conf_mapped_barcoded_reads_cum_frac"),
                )?
            };

            cell_hero_metrics_rows.push(GexSampleHeroMetricsRow {
                genome: if genomes.len() > 1 {
                    Some(genome.to_string())
                } else {
                    None
                },
                total_singlets: get_metric_usize(
                    metrics,
                    &format!("{genome}_singlets_assigned_to_this_sample"),
                )?,
                median_reads_per_singlet,
                median_reads_per_cell_on_target,
                median_genes_per_singlet,
                total_genes_detected: total_genes_detected.map(|x| x.0 as usize),
                median_umi_per_singlet,
                median_genes_per_cell_on_target,
                num_genes_detected_on_target: num_genes_detected_on_target.map(|x| x.0 as usize),
                median_umis_per_cell_on_target,
                confidently_mapped_reads_in_cells,
            });
        }

        let barcode_rank_plot = if self.is_rtl_multiplexed() | self.is_overhang_multiplexed() {
            Some(format_barcode_rank_plot(
                self.sample_barcode_rank_plots
                    .get(sample_assignment)
                    .unwrap()
                    .get(&LegacyLibraryType::GeneExpression)
                    .expect("GEX sample barcode rank plot not found"),
                "GEX",
            ))
        } else {
            None
        };

        let gdna_table = if self.is_rtl_multiplexed() {
            Some(GdnaMetricsTable::from_metrics(metrics)?.into())
        } else {
            None
        };

        Ok(SampleGexWebSummary {
            parameters_table: self
                .lib_ws_builder
                .count_param_table(LegacyLibraryType::GeneExpression)?,
            hero_metrics: GexSampleHeroMetricsTable(cell_hero_metrics_rows).into(),
            cell_metrics_table: if self.multiplexing_method.is_some() {
                Some(
                    self.sample_cell_metrics_table(metrics, physical_library_id)?
                        .into(),
                )
            } else {
                None
            },
            mapping_metrics_table: if self.multiplexing_method.is_some() {
                Some(self.get_mapping_metrics_table(metrics)?.into())
            } else {
                None
            },
            gdna_table,
            median_genes_per_cell_plot: if self.multiplexing_method.is_some() {
                Some(sample_median_genes_plot_from_metrics(
                    metrics,
                    TxHashSet::from_iter(genomes.into_iter()),
                    PlotType::SamplePlot,
                    self.targeting_method,
                ))
            } else {
                None
            },
            clustering_and_diffexp_plots: self
                .sample_tsne_plots
                .get(sample_assignment)
                .cloned()
                .unwrap_or_default()
                .gex_diffexp_clustering_plots,
            barcode_rank_plot,
        })
    }

    fn build_antibody_ws(
        &self,
        sample_assignment: &SampleAssignment,
        physical_library_id: &str,
    ) -> Result<SampleAntibodyWebSummary> {
        let metrics = self
            .per_sample_metrics
            .get(sample_assignment)
            .expect("Sample metrics file not found.");
        let tsne_plots = self
            .sample_tsne_plots
            .get(sample_assignment)
            .cloned()
            .unwrap_or_default();

        // if this is antibody_only case, we don't show the single umi-on-TSNE plot.
        let umi_on_tsne_plot = match tsne_plots.antibody_diffexp_clustering_plots {
            Value::Null => Some(format_umi_on_tsne_plot(
                &tsne_plots.antibody_umi_on_tsne,
                "Antibody Capture",
                "t-SNE Projection",
            )),
            _ => None,
        };

        let antibody_treemap: Option<RawChartWithHelp> =
            if let Some(ref sample_treemap) = self.sample_treemap_plots {
                sample_treemap[sample_assignment]
                    .get(&LegacyLibraryType::AntibodyCapture)
                    .cloned()
            } else {
                None
            };

        let barcode_rank_plot = if self.is_rtl_multiplexed() {
            Some(format_barcode_rank_plot(
                &self.sample_barcode_rank_plots[sample_assignment]
                    [&LegacyLibraryType::AntibodyCapture],
                "AB",
            ))
        } else {
            None
        };

        Ok(SampleAntibodyWebSummary {
            parameters_table: self
                .lib_ws_builder
                .count_param_table(LegacyLibraryType::AntibodyCapture)?,
            hero_metrics: AntibodySampleHeroMetricsTable::from_metrics(metrics)?.into(),
            cell_metrics_table: if self.multiplexing_method.is_some() {
                Some(
                    self.sample_cell_metrics_table(metrics, physical_library_id)?
                        .into(),
                )
            } else {
                None
            },
            mapping_metrics_table: AntibodySampleMappingMetricsTable::from_metrics(metrics)?.into(),
            antibody_treemap,
            clustering_and_diffexp_plots: Some(tsne_plots.antibody_diffexp_clustering_plots),
            tsne_plot: umi_on_tsne_plot,
            barcode_rank_plot,
            feature_histogram: self
                .sample_antibody_histograms
                .as_ref()
                .map(|histos_per_sample| histos_per_sample[sample_assignment].clone()),
        })
    }

    fn build_antigen_ws(
        &self,
        sample_assignment: &SampleAssignment,
    ) -> Result<SampleAntigenWebSummary> {
        let metrics = self
            .per_sample_metrics
            .get(sample_assignment)
            .expect("Sample metrics file not found.");

        let mut hero_metrics_row = vec![AntigenSampleHeroMetricsRow {
            feature_type: Some(LegacyLibraryType::GeneExpression.to_string()),
            total_singlets: get_metric_usize(metrics, "ANTIGEN_multi_filtered_bcs")?,
            median_umis_per_singlet: get_metric_round_float_to_int(
                metrics,
                "ANTIGEN_multi_filtered_bcs_median_counts",
            )?,
            antigen_reads_usable_per_cell: get_metric_round_float_to_int(
                metrics,
                "ANTIGEN_multi_usable_reads_per_filtered_bc",
            )?,
        }];

        if let Some(AntigenVdjMetrics {
            num_cells,
            mean_usable_reads_per_cell,
            median_umis_per_cell,
        }) = self.antigen_vdj_metrics.as_ref()
        {
            hero_metrics_row.push(AntigenSampleHeroMetricsRow {
                feature_type: Some(
                    self.lib_ws_builder
                        .vdj_tab_names()
                        .into_iter()
                        .flatten()
                        .exactly_one() // Exactly 1 VDJ library guaranteed for antigen runs
                        .unwrap()
                        .to_string(),
                ),
                total_singlets: Some(*num_cells as usize),
                median_umis_per_singlet: Some(FloatAsInt(*median_umis_per_cell)),
                antigen_reads_usable_per_cell: Some(FloatAsInt(*mean_usable_reads_per_cell)),
            })
        }

        let antigen_treemap: Option<RawChartWithHelp> =
            if let Some(ref sample_treemap) = self.sample_treemap_plots {
                sample_treemap
                    .get(sample_assignment)
                    .unwrap()
                    .get(&LegacyLibraryType::AntigenCapture)
                    .cloned()
            } else {
                None
            };

        Ok(SampleAntigenWebSummary {
            parameters_table: self
                .lib_ws_builder
                .count_param_table(LegacyLibraryType::AntigenCapture)?,
            hero_metrics: AntigenSampleHeroMetricsTable(hero_metrics_row).into(),
            antigen_treemap,
            clonotype_clustermap: self.clonotype_clustermap.clone(),
        })
    }

    fn build_crispr_ws(
        &self,
        sample_assignment: &SampleAssignment,
        physical_library_id: &str,
    ) -> Result<SampleCrisprWebSummary> {
        let metrics = self
            .per_sample_metrics
            .get(sample_assignment)
            .expect("Sample metrics file not found.");

        Ok(SampleCrisprWebSummary {
            parameters_table: self
                .lib_ws_builder
                .count_param_table(LegacyLibraryType::CrisprGuideCapture)?,
            hero_metrics: CrisprSampleHeroMetricsTable(vec![CrisprSampleHeroMetricsRow {
                total_singlets: get_metric_usize(metrics, "CRISPR_multi_filtered_bcs")?,
                median_umis_per_singlet: get_metric_round_float_to_int(
                    metrics,
                    "CRISPR_multi_filtered_bcs_median_counts",
                )?,
                guide_reads_usable_per_cell: get_metric_round_float_to_int(
                    metrics,
                    "CRISPR_multi_usable_reads_per_filtered_bc",
                )?,
                cells_with_one_or_more_protospacers_detected: get_metric_percent(
                    metrics,
                    "CRISPR_frac_cells_with_protospacer",
                )?,
                cells_with_two_or_more_protospacers_detected: get_metric_percent(
                    metrics,
                    "CRISPR_frac_cells_with_multiple_protospacer",
                )?,
            }])
            .into(),
            cell_metrics_table: if self.multiplexing_method.is_some() {
                Some(
                    self.sample_cell_metrics_table(metrics, physical_library_id)?
                        .into(),
                )
            } else {
                None
            },
            tsne_plot: self.sample_tsne_plots.get(sample_assignment).map(|x| {
                format_umi_on_tsne_plot(
                    &x.crispr_umi_on_tsne,
                    "CRISPR Guide Capture",
                    "t-SNE Projection",
                )
            }),
            barcode_rank_plot: None,
        })
    }

    fn build_custom_ws(
        &self,
        sample_assignment: &SampleAssignment,
        physical_library_id: &str,
    ) -> Result<SampleCustomFeatureWebSummary> {
        let metrics = self
            .per_sample_metrics
            .get(sample_assignment)
            .expect("Sample metrics file not found.");

        Ok(SampleCustomFeatureWebSummary {
            parameters_table: self
                .lib_ws_builder
                .count_param_table(LegacyLibraryType::Custom)?,
            hero_metrics: CustomFeatureSampleHeroMetricsTable(vec![
                CustomFeatureSampleHeroMetricsRow {
                    total_singlets: get_metric_usize(metrics, "Custom_multi_filtered_bcs")?,
                    median_umis_per_singlet: get_metric_round_float_to_int(
                        metrics,
                        "Custom_multi_filtered_bcs_median_counts",
                    )?,
                    feature_reads_usable_per_cell: get_metric_round_float_to_int(
                        metrics,
                        "Custom_multi_usable_reads_per_filtered_bc",
                    )?,
                },
            ])
            .into(),
            cell_metrics_table: if self.multiplexing_method.is_some() {
                Some(
                    self.sample_cell_metrics_table(metrics, physical_library_id)?
                        .into(),
                )
            } else {
                None
            },
            tsne_plot: Some(format_umi_on_tsne_plot(
                &self
                    .sample_tsne_plots
                    .get(sample_assignment)
                    .cloned()
                    .unwrap_or_default()
                    .custom_umi_on_tsne,
                "Custom Feature",
                "t-SNE Projection",
            )),
            barcode_rank_plot: None,
        })
    }
}

fn get_metric(metrics: &TxHashMap<String, Value>, key: &str) -> Result<Option<Number>> {
    match metrics.get(key) {
        Some(val) => match val {
            Value::Number(n) => Ok(Some(n.clone())),
            Value::String(s) => {
                if s.to_lowercase() == "nan" {
                    Ok(None)
                } else {
                    bail!("JSON metric {} had unexpected type {:?}", key, val);
                }
            }
            _ => {
                bail!("JSON metric {} had unexpected type {:?}", key, val);
            }
        },
        None => Ok(None),
    }
}

fn get_metric_usize(metrics: &TxHashMap<String, Value>, key: &str) -> Result<Option<usize>> {
    Ok(get_metric(metrics, key)?.map(|x| {
        x.as_f64()
            .expect("Error converting metric to float in get_metric_usize.")
            .round() as usize
    }))
}

fn get_metric_f64(metrics: &TxHashMap<String, Value>, key: &str) -> Result<Option<f64>> {
    Ok(get_metric(metrics, key)?.map(|x| {
        x.as_f64()
            .expect("Error converting metric to float in get_metric_f64.")
    }))
}

fn get_metric_percent(metrics: &TxHashMap<String, Value>, key: &str) -> Result<Option<Percent>> {
    Ok(get_metric_f64(metrics, key)?.map(Percent::Float))
}

fn get_metric_round_float_to_int(
    metrics: &TxHashMap<String, Value>,
    key: &str,
) -> Result<Option<FloatAsInt>> {
    Ok(get_metric_f64(metrics, key)?.map(FloatAsInt))
}

fn get_metric_string(metrics: &TxHashMap<String, Value>, key: &str) -> Result<Option<String>> {
    match metrics.get(key) {
        Some(val) => match val {
            Value::String(n) => Ok(Some(n.clone())),
            Value::Bool(b) => Ok(Some(b.to_string())),
            Value::Number(n) => Ok(Some(n.to_string())),
            Value::Null => Ok(None),
            Value::Array(_) | Value::Object(_) => {
                bail!(
                    // probably makes sense to have these error out
                    "JSON metric {} was accessed as a String but was an Array or an Object.",
                    key
                );
            }
        },
        None => Ok(None),
    }
}

fn identify_dropped_tags(
    multi_cfg_file: &MultiConfigCsvFile,
    tag_contaminant_info: &Option<Value>,
) -> Vec<String> {
    let Some(contaminants) = tag_contaminant_info.as_ref() else {return Default::default()};

    // Load expected tags
    let cfg = multi_cfg_file
        .read()
        .expect("Could not load configuration CSV.");
    let used_tags = cfg.sample_barcode_ids_used_in_experiment(ProbeBarcodeIterationMode::Mapped);

    contaminants
        .as_object()
        .expect("Contaminant info not a map.")
        .iter()
        .filter(|t| {
            t.1["is_contaminant"]
                .as_bool()
                .expect("Contaminant flag missing from JSON")
                & used_tags.contains(t.0)
        })
        .map(|z| z.0.to_string())
        .collect()
}

/// Return true if this sample is associated with the provided FLEX-compatible library type.
/// Return false if not multiplexed, or if the multiplexing type is not RTL.
fn sample_is_associated_with_flex_library_type(
    sample: &Sample,
    library_type: LegacyLibraryType,
) -> bool {
    if sample.cell_multiplexing_type() != Some(CellMultiplexingType::RTL) {
        return false;
    }
    match library_type {
        LegacyLibraryType::GeneExpression => sample
            .tag_names()
            .any(|tag| categorize_multiplexing_barcode_id(tag) == MultiplexingBarcodeType::RTL),
        LegacyLibraryType::AntibodyCapture => sample.tag_names().any(|tag| {
            categorize_multiplexing_barcode_id(tag) == MultiplexingBarcodeType::Antibody
        }),
        _ => false,
    }
}

#[make_mro(volatile = strict, mem_gb = 5)]
impl MartianMain for WriteMultiWebSummaryJson {
    type StageInputs = StageInputs;
    type StageOutputs = StageOutputs;

    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        let multi_graph = args.multi_graph.read()?;

        let lib_metrics = args
            .library_metrics
            .as_ref()
            .map(FileTypeRead::read)
            .transpose()?
            .unwrap_or_default();

        // maps sample assignment to websummary
        let mut sample_to_web_summary = TxHashMap::default();
        let mut sample_to_metrics_csv = TxHashMap::default();

        // Determine if any tags specified in the CSV config were dropped as contaminants.
        let contaminant_tags = args.tag_contaminant_info.map(|f| f.read()).transpose()?;
        let dropped_tags = identify_dropped_tags(&args.multi_config, &contaminant_tags);

        // For RTL samples alert on any probe barcodes that have very high GEM barcode overlap
        let probe_barcode_overlap_coefficients = lib_metrics
            .get("probe_barcode_overlap_coefficients")
            .cloned();
        let probe_barcodes_high_gem_overlap = probe_barcode_overlap_coefficients
            .as_ref()
            .map(|coefficients| {
                coefficients
                    .as_object()
                    .unwrap()
                    .into_iter()
                    .sorted_by_key(|&(_, x)| NotNan::new(x.as_f64().unwrap()).unwrap())
                    .rev()
                    .filter_map(|(k, v)| {
                        let v_f64 = v.as_f64().unwrap();
                        if v_f64 >= GEM_BARCODE_OVERLAP_ALERT_THRESHOLD {
                            Some(format!("{k} ({:.2}%)", v_f64 * 100.0))
                        } else {
                            None
                        }
                    })
                    .collect()
            })
            .unwrap_or_default();

        //******************************************************************************************
        // POPULATE LIBRARY WEBSUMMARIES
        //******************************************************************************************
        let lib_ws_builder = LibWsBuilder {
            common_inputs: args.common_inputs.clone(),
            multi_graph: multi_graph.clone(),
            multi_config: args.multi_config.read()?,
            lib_metrics,
            barcode_rank_plots: match args.barcode_rank_plots {
                Some(ref plots) => plots.read()?,
                None => TxHashMap::default(),
            },
            sequencing_metrics: match args.sequencing_metrics {
                Some(f) => f.read()?,
                None => TxHashMap::default(),
            },
            count_inputs: args.count_inputs.clone(),
            jibes_biplot_histogram: args.jibes_biplot_histogram.map(|j| j.read()).transpose()?,
            antibody_histograms: args.antibody_histograms.map(|j| j.read()).transpose()?,
            antigen_histograms: args.antigen_histograms.map(|j| j.read()).transpose()?,
            cmo_tsne_plot: args.cmo_tsne_plot.map(|j| j.read()).transpose()?,
            vdj_t_contents: args.vdj_t_contents.map(|f| f.read()).transpose()?,
            vdj_t_gd_contents: args.vdj_t_gd_contents.map(|f| f.read()).transpose()?,
            vdj_b_contents: args.vdj_b_contents.map(|f| f.read()).transpose()?,
            target_set_name: args.target_set_name,
            targeted_per_feature_metrics: args.targeted_per_feature_metrics.clone(),
            is_multiplexed: multi_graph.has_multiplexing(),
            specificity_controls: args
                .feature_config
                .unwrap_or(FeatureConfig {
                    beam_mode: None,
                    specificity_controls: None,
                    functional_map: None,
                })
                .specificity_controls,
            dropped_tags,
            probe_barcodes_high_gem_overlap,
            mismatched_probe_barcode_pairings: get_mismatched_probe_barcode_pairings(
                &args.multi_config,
                args.detected_probe_barcode_pairing.as_ref(),
            )?,
        };

        //******************************************************************************************
        // POPULATE DIAGNOSTIC DATA
        //******************************************************************************************
        let get_metric = |metric| lib_ws_builder.lib_metrics.get(metric).cloned();

        let diagnostics = MultiDiagnostics {
            corrected_bc_frac: get_metric("corrected_bc_frac"),
            corrected_bc_in_gel_bead_frac: get_metric("corrected_bc_in_gel_bead_frac"),
            corrected_bc_in_probe_frac: get_metric("corrected_bc_in_probe_frac"),
            ANTIBODY_corrected_bc_frac: get_metric("ANTIBODY_corrected_bc_frac"),
            ANTIBODY_corrected_bc_in_gel_bead_frac: get_metric(
                "ANTIBODY_corrected_bc_in_gel_bead_frac",
            ),
            ANTIBODY_corrected_bc_in_probe_frac: get_metric("ANTIBODY_corrected_bc_in_probe_frac"),
            i1_bases_with_q30_frac: get_metric("i1_bases_with_q30_frac"),
            i2_bases_with_q30_frac: get_metric("i2_bases_with_q30_frac"),
            low_support_umi_reads_frac: get_metric("low_support_umi_reads_frac"),
            tag_contaminant_info: contaminant_tags,
            tso_frac: get_metric("tso_frac"),
            probe_barcode_overlap_coefficients,
            fraction_reads_high_occupancy_gems: get_metric(
                "rtl_multiplexing_fraction_reads_high_occupancy_gems",
            ),
            high_occupancy_probe_barcode_count_threshold: get_metric(
                "rtl_multiplexing_high_occupancy_probe_barcode_count_threshold",
            ),
            vdj_t: lib_ws_builder
                .vdj_t_contents
                .as_ref()
                .map(VdjWsContents::diagnostics),
            vdj_t_gd: lib_ws_builder
                .vdj_t_gd_contents
                .as_ref()
                .map(VdjWsContents::diagnostics),
            vdj_b: lib_ws_builder
                .vdj_b_contents
                .as_ref()
                .map(VdjWsContents::diagnostics),
        };

        //******************************************************************************************
        // POPULATE SAMPLE WEBSUMMARIES
        //******************************************************************************************
        let chemistry = lib_ws_builder.chemistry();
        let is_hybrid_capture = lib_ws_builder.is_hybrid_capture();
        let is_rtl = lib_ws_builder.is_rtl();

        let multi_ws_builder = MultiWsBuilder {
            lib_ws_builder,
            per_sample_metrics: read_optional_file_map(&args.per_sample_metrics)?,
            sample_tsne_plots: read_optional_file_map(&args.sample_tsne_plots)?,
            sample_barcode_rank_plots: read_optional_file_map(&args.sample_barcode_rank_plots)?,
            sample_treemap_plots: args
                .sample_treemap_plots
                .as_ref()
                .map(read_optional_file_map)
                .transpose()?,
            sample_antibody_histograms: args
                .sample_antibody_histograms
                .as_ref()
                .map(read_optional_file_map)
                .transpose()?,
            svg_str: std::fs::read_to_string(args.multi_graph_svg)
                .expect("Error reading  multi graph svg"),
            csv_str: std::fs::read_to_string(&args.multi_config)?
                .replace("\r\n", "\n")
                .replace('\r', "\n"),
            diagnostics,
            pipeline_version: rover.pipelines_version(),
            context: AlertContext {
                is_hybrid_capture,
                is_rtl,
                is_lt_chemistry: chemistry == Some(ThreePrimeV3LT),
                is_arc_chemistry: chemistry == Some(ArcV1),
                is_fiveprime: chemistry.map_or(false, |x| {
                    matches!(x, FivePrimeR1 | FivePrimeR2 | FivePrimePE)
                }),
                is_multiplexing: args.count_inputs.as_ref().map_or(false, |inputs| {
                    inputs
                        .sample_def
                        .iter()
                        .any(|sdef| sdef.library_type == Some(LegacyLibraryType::Multiplexing))
                }),
                is_antigen: args.count_inputs.as_ref().map_or(false, |inputs| {
                    inputs
                        .sample_def
                        .iter()
                        .any(|sdef| sdef.library_type == Some(LegacyLibraryType::AntigenCapture))
                }),
                include_introns: args
                    .count_inputs
                    .as_ref()
                    .map_or(false, |inputs| inputs.include_introns),
                no_preflight: args.no_preflight,
            },
            multiplexing_method: multi_graph.get_multiplexing_method(),
            targeting_method: args
                .count_inputs
                .as_ref()
                .and_then(|inputs| inputs.targeting_method),
            antigen_vdj_metrics: args.antigen_vdj_metrics.map(|f| f.read()).transpose()?,
            clonotype_clustermap: args
                .antigen_specificity
                .map(clonotype_specificity_heatmap)
                .transpose()?
                .flatten(),
        };

        let per_sample_ws = multi_ws_builder.build()?;

        for (sample, sample_ws) in per_sample_ws {
            // Write the web summary data to JSON
            // put the path to the JSON into sample_to_web_summary, with key being sample ID
            let json_file: JsonFile<()> =
                rover.make_path(format!("{}_web_summary_data.json", sample.clone()));

            serde_json::to_writer_pretty(json_file.buf_writer()?, &sample_ws)?;
            sample_to_web_summary.insert(sample.clone(), json_file);

            let csv_file: CsvFile<()> =
                rover.make_path(format!("{}_metric_summary_csv", sample.clone()));
            sample_ws.to_csv(&csv_file)?;
            sample_to_metrics_csv.insert(sample, csv_file);
        }

        Ok(StageOutputs {
            web_summary_json: sample_to_web_summary,
            metrics_summary_csv: sample_to_metrics_csv,
        })
    }
}

/// Read a mapping of optional files.
fn read_optional_file_map<F: FileTypeRead<T>, T, K: Clone + Hash + Eq>(
    files: &TxHashMap<K, Option<F>>,
) -> Result<TxHashMap<K, T>> {
    read_files_into_map(
        files
            .iter()
            .filter_map(|(k, opt)| opt.as_ref().map(|v| (k, v))),
    )
}

/// Read an iterator of keys/martian files into a map containing the parsed file contents.
fn read_files_into_map<'a, F: FileTypeRead<T> + 'a, T, K: Clone + Hash + Eq + 'a>(
    files: impl IntoIterator<Item = (&'a K, &'a F)>,
) -> Result<TxHashMap<K, T>> {
    let mut out = TxHashMap::default();
    for (key, file) in files {
        out.insert(key.clone(), file.read()?);
    }
    Ok(out)
}

/// Return any differences between configured and detected probe barcode pairings.
/// Return None if there are no differences, or if there is no detected or
/// configured pairing.
fn get_mismatched_probe_barcode_pairings(
    multi_config: &MultiConfigCsvFile,
    detected_probe_barcode_pairing: Option<&DetectedProbeBarcodePairingFile>,
) -> Result<Option<MismatchedProbeBarcodePairings>> {
    let Some(samples) = multi_config.read()?.samples else {return Ok(None);};
    let Some(detected_probe_barcode_pairing) = detected_probe_barcode_pairing else {return Ok(None);};
    let detected_probe_barcode_pairing: HashSet<_> = detected_probe_barcode_pairing
        .read()?
        .into_iter()
        .map(|(target_bc, source_bc)| format!("{target_bc}{PROBE_BARCODE_ID_GROUPING}{source_bc}"))
        .collect();

    let configured_probe_barcode_pairing: HashSet<_> = samples
        .get_translated_probe_barcodes()
        .into_iter()
        .map(|(source_bc, target_bc)| format!("{target_bc}{PROBE_BARCODE_ID_GROUPING}{source_bc}"))
        .collect();
    if configured_probe_barcode_pairing.is_empty()
        || configured_probe_barcode_pairing == detected_probe_barcode_pairing
    {
        return Ok(None);
    }
    Ok(Some(MismatchedProbeBarcodePairings {
        configured_not_detected: configured_probe_barcode_pairing
            .difference(&detected_probe_barcode_pairing)
            .cloned()
            .collect(),
        detected_not_configured: detected_probe_barcode_pairing
            .difference(&configured_probe_barcode_pairing)
            .cloned()
            .collect(),
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sorted_cmos() {
        let mut rows = vec![
            MultiplexingCmoMetricsRow {
                gem_well_cmo: Some("CMO_2_5".to_string()),
                reads_in_cell_associated_partitions: Some(Percent::Float(0.822)),
                singlets_assigned_to_cmo: Some(Percent::Float(0.822)),
                cmo_signal_to_background_ratio: Some(3.0),
            },
            MultiplexingCmoMetricsRow {
                gem_well_cmo: Some("CMO_2_3".to_string()),
                reads_in_cell_associated_partitions: Some(Percent::Float(0.822)),
                singlets_assigned_to_cmo: Some(Percent::Float(0.822)),
                cmo_signal_to_background_ratio: Some(3.0),
            },
            MultiplexingCmoMetricsRow {
                gem_well_cmo: Some("CMO_2_4".to_string()),
                reads_in_cell_associated_partitions: Some(Percent::Float(0.822)),
                singlets_assigned_to_cmo: Some(Percent::Float(0.822)),
                cmo_signal_to_background_ratio: Some(2.0),
            },
        ];
        sort_cmo_rows(&mut rows);
        let names: Vec<Option<String>> = rows.into_iter().map(|x| x.gem_well_cmo).collect();
        assert_eq!(
            names,
            vec![
                Some("CMO_2_3".to_string()),
                Some("CMO_2_4".to_string()),
                Some("CMO_2_5".to_string()),
            ]
        );
    }
}
