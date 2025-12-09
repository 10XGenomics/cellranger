//! Martian stage WRITE_MULTI_WEB_SUMMARY
//! Write HTML files for websummaries.
#![deny(missing_docs)]

use super::parse_multi_config::CellCalling;
use crate::cell_annotation_ws_parameters::{
    CellAnnotationMetrics, CellTypeWebSummaryBundle, CellTypeWebSummaryBundleValue,
    generate_cell_type_barchart_from_value, generate_cell_type_diffexp_from_value,
    generate_cell_type_parameter_table, generate_cell_type_umap_plot_from_value,
    generate_cell_type_violin_plot_from_value,
};
use crate::stages::DetectedProbeBarcodePairingFile;
use crate::stages::build_per_sample_vdj_ws_contents::{VdjWsContents, VdjWsContentsFormat};
use crate::stages::compute_antigen_vdj_metrics::AntigenVdjMetrics;
use crate::stages::parse_multi_config::{CommonInputs, CountInputs};
use crate::stages::pick_beam_analyzer::BeamAnalyzerOutputs;
use crate::{HtmlFile, PerLibrarySequencingMetrics, SequencingMetricsFormat};
use anyhow::{Context, Result, anyhow, bail};
use barcode::whitelist::{RTLMultiplexingBarcodeType, categorize_rtl_multiplexing_barcode_id};
use cr_types::chemistry::{ChemistryDefs, ChemistryDefsExt, ChemistryName};
use cr_types::reference::feature_reference::{FeatureConfig, SpecificityControls};
use cr_types::reference::reference_info::MULTI_GENOME_SEPARATOR;
use cr_types::websummary::{AlertConfig, ExtractMode, MetricEtlConfig};
use cr_types::{
    AlignerParam, BarcodeMultiplexingType, CellLevel, CrMultiGraph, Fingerprint, GenomeName,
    LibraryType, ReadLevel, Sample, SampleAssignment, TargetingMethod,
};
use cr_websummary::multi::antigen::clonotype_specificity_heatmap;
use cr_websummary::multi::metrics::{
    ActiveConditions, CountAndPercentTransformer, MetricTier, MetricsProcessor, load_metrics_etl,
    load_metrics_etl_special,
};
use cr_websummary::multi::plots::{
    PlotType, format_barcode_rank_plot, format_barnyard_biplot, format_histogram,
    format_jibes_biplots, format_tags_on_umap_plot, format_umi_on_umap_plot,
    library_median_genes_plot_from_metrics, library_sequencing_saturation_plot_from_metrics,
    sample_median_genes_plot_from_metrics,
};
use cr_websummary::multi::websummary::{
    CELL_ANNOTATION_ADVERTISEMENT_STRING, CountParametersTable, ExperimentalDesign,
    JsonMetricSummary, LibraryAntibodyOrAntigenWebSummary, LibraryCmoWebSummary,
    LibraryCrisprWebSummary, LibraryCustomFeatureWebSummary, LibraryGexWebSummary,
    LibraryHashtagWebSummary, LibraryWebSummary, MetricsTraitWrapper,
    MismatchedProbeBarcodePairings, MultiDiagnostics, MultiSharedResource, MultiWebSummary,
    MultiWebSummaryLibraryData, MultiWebSummarySampleData, SampleAntibodyWebSummary,
    SampleAntigenWebSummary, SampleCellAnnotationWebSummary, SampleCrisprWebSummary,
    SampleCustomFeatureWebSummary, SampleDiagnostics, SampleGexWebSummary, SampleWebSummary,
    Section,
};
use cr_websummary::{AlertContext, AlertLevel, AlertSpec, ChartWithHelp, Tab};
use fastq_set::WhichEnd;
use itertools::Itertools;
use json_report_derive::JsonReport;
use martian::prelude::*;
use martian_derive::{MartianStruct, make_mro};
use martian_filetypes::FileTypeRead;
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::tabular_file::CsvFile;
use metric::{TxHashMap, join_metric_name};
use multi::config::{
    MultiConfigCsv, MultiConfigCsvFile, PROBE_BARCODE_ID_GROUPING, ProbeBarcodeIterationMode,
};
use ordered_float::NotNan;
use serde::{Deserialize, Serialize};
use serde_json::value::{RawValue, Value};
use serde_json::{self, Number, json};
use std::collections::{HashMap, HashSet};
use std::fs::read_to_string;
use std::hash::Hash;
use std::io::Write;
use std::string::ToString;

pub const GEM_BARCODE_OVERLAP_ALERT_THRESHOLD: f64 = 0.6;

/// Martian stage WRITE_MULTI_WEB_SUMMARY
pub struct WriteMultiWebSummary;

#[derive(Clone, Serialize, Deserialize)]
pub struct JibesBiplotHistogramData {
    pub biplot: Box<RawValue>,
    pub histogram: Box<RawValue>,
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
    pub common_inputs: CommonInputs,
    pub count_inputs: Option<CountInputs>,
    pub count_cell_calling_config: Option<CellCalling>,
    pub tag_contaminant_info: Option<JsonFile<Value>>,
    pub sample_projection_plots: TxHashMap<SampleAssignment, Option<JsonFile<SampleUmapPlots>>>,
    pub sample_barcode_rank_plots:
        TxHashMap<SampleAssignment, Option<JsonFile<TxHashMap<LibraryType, Box<RawValue>>>>>,
    pub sample_treemap_plots: Option<
        TxHashMap<SampleAssignment, Option<JsonFile<TxHashMap<LibraryType, Box<RawValue>>>>>,
    >,
    pub barcode_rank_plots: Option<JsonFile<TxHashMap<LibraryType, Box<RawValue>>>>,
    pub jibes_biplot_histogram: Option<JsonFile<JibesBiplotHistogramData>>,
    pub barnyard_biplot: Option<JsonFile<Box<RawValue>>>,
    pub antibody_histograms: Option<JsonFile<Box<RawValue>>>,
    pub sample_antibody_histograms:
        Option<TxHashMap<SampleAssignment, Option<JsonFile<Box<RawValue>>>>>,
    pub antigen_histograms: Option<JsonFile<Box<RawValue>>>,
    pub cmo_projection_plot: Option<JsonFile<MultiplexingUmapPlots>>,
    pub vdj_t_contents: Option<TxHashMap<SampleAssignment, VdjWsContentsFormat>>,
    pub vdj_t_gd_contents: Option<TxHashMap<SampleAssignment, VdjWsContentsFormat>>,
    pub vdj_b_contents: Option<TxHashMap<SampleAssignment, VdjWsContentsFormat>>,
    pub target_set_name: Option<String>,
    pub beam_analyzer: Option<TxHashMap<SampleAssignment, Option<BeamAnalyzerOutputs>>>,
    pub cell_type_websummary_bundles:
        Option<TxHashMap<SampleAssignment, Option<CellTypeWebSummaryBundle>>>,
    pub cell_annotation_viable_but_not_requested: Option<TxHashMap<SampleAssignment, Option<bool>>>,
    pub feature_config: Option<FeatureConfig>,
    /// chemistry_defs is None for a VDJ-only analysis.
    pub chemistry_defs: Option<ChemistryDefs>,
    pub detected_probe_barcode_pairing: Option<DetectedProbeBarcodePairingFile>,
    pub no_preflight: bool,
}

#[derive(Clone, Serialize, Deserialize, MartianStruct)]
pub struct StageOutputs {
    pub library_web_summary: HtmlFile,
    pub qc_library_metrics: CsvFile<()>,
    pub qc_sample_metrics: CsvFile<()>,
    pub web_summary: TxHashMap<SampleAssignment, HtmlFile>,
    pub web_summary_json: TxHashMap<SampleAssignment, JsonFile<()>>,
    pub metrics_summary_csv: TxHashMap<SampleAssignment, CsvFile<()>>,
}

/// The data injected into the library websummary.
///
/// This is a "declining-balance account". As we build the different websummary
/// tabs, we *remove* data from this data structure. It should be completely
/// empty at the end of web summary construction.
struct LibWsBuilderContents {
    barcode_rank_plots: TxHashMap<LibraryType, Box<RawValue>>,
    sequencing_metrics: PerLibrarySequencingMetrics,
    jibes_biplot_histogram: Option<JibesBiplotHistogramData>,
    antibody_histograms: Option<Box<RawValue>>,
    antigen_histograms: Option<Box<RawValue>>,
    barnyard_biplot: Option<Box<RawValue>>,
    cmo_projection_plot: Option<MultiplexingUmapPlots>,
    vdj_t_contents: Option<VdjWsContents>,
    vdj_t_gd_contents: Option<VdjWsContents>,
    vdj_b_contents: Option<VdjWsContents>,
}

struct LibWsBuilder {
    alert_context: AlertContext,
    common_inputs: CommonInputs,
    multi_graph: CrMultiGraph,
    multi_config: MultiConfigCsv,
    chemistry_defs: Option<ChemistryDefs>,
    lib_metrics: TxHashMap<String, Value>,
    lib_metrics_proc: MetricsProcessor,
    special_metrics_proc: MetricsProcessor,
    count_inputs: Option<CountInputs>,
    count_cell_calling_config: Option<CellCalling>,
    target_set_name: Option<String>,

    /// True if any form of multiplexing is used in this analysis.
    multiplexing_method: Option<BarcodeMultiplexingType>,
    specificity_controls: Option<SpecificityControls>,

    dropped_tags: Vec<String>,
    probe_barcodes_high_gem_overlap: Vec<String>,
    mismatched_probe_barcode_pairings: Option<MismatchedProbeBarcodePairings>,
    /// Derived parameters.
    is_rtl: bool,
    targeting_method: Option<TargetingMethod>,
}

fn vdj_tab_names(ws: &LibraryWebSummary) -> [Option<&'static str>; 3] {
    [
        ws.vdj_t_tab.is_some().then_some("VDJ-T"),
        ws.vdj_b_tab.is_some().then_some("VDJ-B"),
        ws.vdj_t_gd_tab.is_some().then_some("VDJ-T-GD"),
    ]
}

impl LibWsBuilder {
    /// Return true if multiplexed using CMO.
    fn is_cmo_multiplexed(&self) -> bool {
        self.multiplexing_method == Some(BarcodeMultiplexingType::CellLevel(CellLevel::CMO))
    }

    /// Return true if multiplexed using Hashtag.
    fn is_hashtag_multiplexed(&self) -> bool {
        self.multiplexing_method == Some(BarcodeMultiplexingType::CellLevel(CellLevel::Hashtag))
    }

    /// Return true if multiplexed using OH
    fn is_oh_multiplexed(&self) -> bool {
        self.multiplexing_method == Some(BarcodeMultiplexingType::ReadLevel(ReadLevel::OH))
    }

    /// Return true if multiplexed using RTL .
    fn is_rtl_multiplexed(&self) -> bool {
        self.multiplexing_method == Some(BarcodeMultiplexingType::ReadLevel(ReadLevel::RTL))
    }

    /// Return true if antigen specificity controls are specified.
    fn has_antigen_controls(&self) -> bool {
        self.specificity_controls.is_some()
    }

    /// Return the aligner.
    fn aligner(&self) -> Option<AlignerParam> {
        get_aligner_from_metrics(&self.lib_metrics)
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

    /// Build the library websummary.
    fn build(
        &self,
        context: &AlertContext,
        mut contents: LibWsBuilderContents,
    ) -> Result<(MultiWebSummaryLibraryData, MultiSharedResource)> {
        let take_vdj = |content: &mut Option<VdjWsContents>| {
            content
                .take()
                .map(|content| Tab::new(content.into_library_ws(), context))
        };

        let mut lib_ws = LibraryWebSummary {
            id: self.common_inputs.sample_id.clone(),
            description: self.common_inputs.sample_desc.clone(),
            vdj_t_tab: take_vdj(&mut contents.vdj_t_contents),
            vdj_t_gd_tab: take_vdj(&mut contents.vdj_t_gd_contents),
            vdj_b_tab: take_vdj(&mut contents.vdj_b_contents),
            ..Default::default()
        };
        let mut resources = MultiSharedResource::default();
        for lib in &self.multi_graph.libraries {
            match lib.library_type {
                LibraryType::Gex => {
                    assert!(lib_ws.gex_tab.is_none());
                    lib_ws.gex_tab = Some(Tab::new(
                        self.build_gex_ws(&lib.physical_library_id, &mut contents)?,
                        context,
                    ));
                }
                LibraryType::Antibody => {
                    assert!(lib_ws.antibody_tab.is_none());
                    lib_ws.antibody_tab = Some(Tab::new(
                        self.build_antibody_or_antigen_ws(
                            &lib.physical_library_id,
                            true,
                            &mut contents,
                        )?,
                        context,
                    ));
                    if self.is_hashtag_multiplexed() {
                        let (hashtag_tab, tab_resources) = self.build_hashtag_ws(&mut contents)?;
                        // bubble up shared resources
                        // FIXME destroy this mechanism
                        resources.extend(tab_resources);
                        lib_ws.hashtag_tab = Some(Tab::new(hashtag_tab, context));
                    }
                }
                LibraryType::Antigen => {
                    assert!(lib_ws.antigen_tab.is_none());
                    lib_ws.antigen_tab = Some(Tab::new(
                        self.build_antibody_or_antigen_ws(
                            &lib.physical_library_id,
                            false,
                            &mut contents,
                        )?,
                        context,
                    ));
                }
                LibraryType::Crispr => {
                    assert!(lib_ws.crispr_tab.is_none());
                    lib_ws.crispr_tab = Some(Tab::new(
                        self.build_crispr_ws(&lib.physical_library_id, &mut contents)?,
                        context,
                    ));
                }
                LibraryType::Vdj(_) => {}
                LibraryType::Custom => {
                    assert!(lib_ws.custom_feature_tab.is_none());
                    lib_ws.custom_feature_tab = Some(Tab::new(
                        self.build_custom_ws(&lib.physical_library_id, &mut contents)?,
                        context,
                    ));
                }
                LibraryType::Cellplex => {
                    assert!(lib_ws.cmo_tab.is_none());
                    let (cmo_tab, tab_resources) =
                        self.build_cmo_ws(&lib.physical_library_id, &mut contents)?;
                    // bubble up shared resources
                    // FIXME destroy this mechanism
                    resources.extend(tab_resources);
                    lib_ws.cmo_tab = Some(Tab::new(cmo_tab, context));
                }
                LibraryType::Atac => unreachable!(),
            }
        }

        let data = MultiWebSummaryLibraryData {
            metrics: lib_ws.to_json_summary(),
            types: self
                .multi_graph
                .libraries
                .iter()
                .map(|lib| lib.library_type)
                .sorted()
                .dedup()
                .collect(),
            data: lib_ws,
        };

        Ok((data, resources))
    }

    /// Conditionally get the high-occupancy GEM metric.
    ///
    /// FIXME CELLRANGER-8444 refactor to always get this metric for RTL multiplexing,
    /// but consider moving the display logic to the frontend. This would impact
    /// whether the metric appears in the CSV output.
    fn get_high_occupancy_gem_metric(
        &self,
        section: Section,
        library_type: LibraryType,
        physical_library_id: &str,
    ) -> Result<Option<JsonMetricSummary>> {
        if !self.is_rtl_multiplexed() {
            return Ok(None);
        }
        let has_gex = self.multi_graph.has_library_type(LibraryType::Gex);
        if !(library_type == LibraryType::Gex || !has_gex && library_type == LibraryType::Antibody)
        {
            return Ok(None);
        }
        let alert = AlertConfig {
                error_threshold: Some(0.0),
                warn_threshold: Some(0.9),
                warn_title: Some("Low fraction of initial cell calls pass high occupancy GEM filtering.".to_string()),
                detail: "Numbers under 90% could be due to partial clogs, wetting failures, cell clumping, or significant deviations from the recommended chip loading protocol.".to_string(),
                ..Default::default()
            };
        let config = MetricEtlConfig {
            json_key: Some("rtl_multiplexing_fraction_cells_in_high_occupancy_gems".to_string()),
            ty: "Percent".to_string(),
            transformer: Some("ComplementPercent".to_string()),
            header: "Fraction of initial cell barcodes passing high occupancy GEM filtering"
                .to_string(),
            alerts: vec![alert],
            ..Default::default()
        };
        let mut metric = self.lib_metrics_proc.process_one(
            &self.alert_context,
            &self.lib_metrics,
            &config,
            MetricTier::Library,
            section,
        )?;

        metric.grouping_key = Some(physical_library_id.to_string());
        metric.grouping_header = Some("Physical library ID".to_string());
        Ok(Some(metric))
    }

    /// Return the chemistry description and append "(manual)" if a manual chemistry is specified.
    fn chemistry_description_with_manual(&self, library_type: LibraryType) -> String {
        let chemistry_defs = self.chemistry_defs.as_ref().unwrap();
        let chemistry_description = &chemistry_defs[&library_type].description;
        let chemistry_spec = self.multi_config.chemistry_specs().unwrap()[&library_type];

        if chemistry_spec.is_auto() {
            chemistry_description.to_string()
        } else {
            format!("{chemistry_description} (manual)")
        }
    }

    fn count_param_table(
        &self,
        library_type: LibraryType,
        probe_barcode_data: PerProbeBarcodeData,
    ) -> Result<CountParametersTable> {
        let count_inputs = self.count_inputs.as_ref().unwrap();
        let cell_calling_config = self.count_cell_calling_config.as_ref().unwrap();

        let transcriptome = if let Some(ref_info) = count_inputs.reference_info.as_ref() {
            format!(
                "{}{}",
                ref_info.genomes.iter().join(MULTI_GENOME_SEPARATOR),
                if let Some(version) = ref_info.version.as_ref() {
                    format!("-{version}")
                } else {
                    String::new()
                }
            )
        } else {
            String::new()
        };

        let reference_path = count_inputs
            .reference_info
            .as_ref()
            .and_then(|x| x.get_reference_path());
        Ok(CountParametersTable {
            chemistry: self.chemistry_description_with_manual(library_type),
            introns_included: count_inputs.include_introns,
            reference_path: reference_path.map(|x| x.display().to_string()),
            transcriptome,
            feature_ref_path: self.get_feature_ref(),
            cmo_set_path: self.get_cmo_set(),
            target_set_name: self.target_set_name.clone(),
            targeting_method: self.targeting_method,
            filter_probes: count_inputs.filter_probes,
            disable_ab_aggregate_detection: cell_calling_config.disable_ab_aggregate_detection,
            disable_high_occupancy_gem_detection: cell_calling_config
                .disable_high_occupancy_gem_detection,
            num_genes_on_target: self
                .lib_metrics
                .get("num_genes_on_target")
                .map(|x| x.as_u64().unwrap() as usize),
            library_type,
            throughput: get_metric_string(&self.lib_metrics, "throughput_inferred")?,
            tenx_cmos: count_inputs.tenx_cmos,
            aligner: self.aligner().unwrap(),
            antigen_negative_control: self.has_antigen_controls(),
            dropped_tags: self.dropped_tags.clone(),
            probe_barcodes_high_gem_overlap: self.probe_barcodes_high_gem_overlap.clone(),
            mismatched_probe_barcode_pairings: self.mismatched_probe_barcode_pairings.clone(),
            unspecified_probe_barcodes_detected: probe_barcode_data.unexpected,
            specified_probe_barcodes_missing: probe_barcode_data.missing,
            unexpected_missing_probe_barcode_threshold: probe_barcode_data
                .unexpected_missing_threshold,
        })
    }

    fn genomes(&self) -> &[GenomeName] {
        self.count_inputs
            .as_ref()
            .expect("count_inputs is None")
            .get_genomes()
            .unwrap_or_default()
    }

    fn build_ocm_per_overhang_metrics(
        &self,
        library_type: LibraryType,
        section: Section,
        active_conditions: &ActiveConditions,
    ) -> Result<Vec<JsonMetricSummary>> {
        if !self.is_oh_multiplexed() {
            return Ok(Default::default());
        }

        let ocm_barcode_id_to_sample_id = self.multi_graph.get_tag_name_to_sample_id_map().unwrap();

        let umi_per_ocm_barcode: Vec<(&str, usize)> = if let Some(umi_per_ocm_barcode) = self
            .lib_metrics
            .get(&join_metric_name(library_type, "umi_per_overhang"))
        {
            umi_per_ocm_barcode
                .as_object()
                .unwrap()
                .iter()
                .map(|(id, n)| (id.as_str(), n.as_u64().unwrap() as usize))
                .sorted()
                .collect()
        } else {
            Vec::default()
        };

        let umi_sum = umi_per_ocm_barcode.iter().map(|(_id, n)| n).sum();

        let filtered_barcodes_per_ocm_barcode: TxHashMap<&str, usize> =
            self.lib_metrics["filtered_barcodes_per_overhang"]
                .as_object()
                .unwrap()
                .iter()
                .map(|(id, n)| (id.as_str(), n.as_u64().unwrap() as usize))
                .collect();

        let filtered_barcodes_sum: usize = filtered_barcodes_per_ocm_barcode
            .iter()
            .map(|(_id, &n)| n)
            .sum();

        let mut special_metrics_proc = self.special_metrics_proc.with_default_transformers();
        special_metrics_proc
            .add_transformer("UmiFraction", CountAndPercentTransformer::new(umi_sum));
        special_metrics_proc.add_transformer(
            "CellsFraction",
            CountAndPercentTransformer::new(filtered_barcodes_sum),
        );

        let mut metrics = vec![];
        for (ocm_barcode_id, umi_count) in umi_per_ocm_barcode {
            let sample_id = ocm_barcode_id_to_sample_id
                .get(ocm_barcode_id)
                .map(|(sample_id, _)| *sample_id);
            let cells = *filtered_barcodes_per_ocm_barcode
                .get(ocm_barcode_id)
                .unwrap_or(&0);

            let cells_per_ocm_barcode = if sample_id.is_some() {
                Some(cells)
            } else {
                assert_eq!(cells, 0);
                None
            };

            #[derive(JsonReport)]
            struct Metrics<'a> {
                ocm_barcode_id: &'a str,
                sample_id: Option<&'a str>,
                umi_per_ocm_barcode: usize,
                cells_per_ocm_barcode: Option<usize>,
            }

            metrics.extend(special_metrics_proc.process_group(
                "ocm_per_overhang_metrics",
                section,
                active_conditions,
                &self.alert_context,
                &Metrics {
                    ocm_barcode_id,
                    sample_id,
                    umi_per_ocm_barcode: umi_count,
                    cells_per_ocm_barcode,
                },
            )?);
        }

        Ok(metrics)
    }
    /// Return per probe barcode metrics for the specified library type.
    /// Also compile missing and unexpected probe barcodes.
    fn build_rtl_probe_barcode_metrics(
        &self,
        library_type: LibraryType,
        section: Section,
        active_conditions: &ActiveConditions,
    ) -> Result<PerProbeBarcodeData> {
        if !self.is_rtl_multiplexed() {
            return Ok(Default::default());
        }

        let mut special_metrics_proc = self.special_metrics_proc.with_default_transformers();

        let probe_barcode_id_to_sample_id_and_mapped_probe_barcode_ids =
            self.multi_graph.get_tag_name_to_sample_id_map().unwrap();

        let umi_per_probe_barcode: Vec<(&str, usize)> = if let Some(umi_per_probe_barcode) = self
            .lib_metrics
            .get(&join_metric_name(library_type, "umi_per_probe_barcode"))
        {
            umi_per_probe_barcode
                .as_object()
                .unwrap()
                .iter()
                .map(|(id, n)| (id.as_str(), n.as_u64().unwrap() as usize))
                .sorted()
                .collect()
        } else {
            Vec::default()
        };

        let umi_sum = umi_per_probe_barcode.iter().map(|(_id, n)| n).sum();

        let filtered_barcodes_per_probe_barcode: TxHashMap<&str, usize> =
            self.lib_metrics["filtered_barcodes_per_probe_barcode"]
                .as_object()
                .unwrap()
                .iter()
                .map(|(id, n)| (id.as_str(), n.as_u64().unwrap() as usize))
                .collect();

        let filtered_barcodes_sum = filtered_barcodes_per_probe_barcode
            .iter()
            .map(|(_id, &n)| n)
            .sum();

        special_metrics_proc
            .add_transformer("UmiFraction", CountAndPercentTransformer::new(umi_sum));
        special_metrics_proc.add_transformer(
            "CellsFraction",
            CountAndPercentTransformer::new(filtered_barcodes_sum),
        );

        let mut metrics = vec![];
        let mut unspecified_probe_barcodes_detected = vec![];
        let mut specified_probe_barcodes_missing = vec![];

        // The threshold to trigger a web summary alert when an unexpected probe barcode is observed or
        // an expected probe barcode is not observed.
        //
        // This is computed dynamically since the absolute concentration of any specific
        // probe barcode will scale inversely with the total number.
        //
        // At 16 probe bcs and below, we maintain the previous behavior of using a fixed
        // cutoff. This cutoff represents 8% of what we would expect to find if UMIs
        // were equally distributed across all probe barcodes. Above 16 probe bcs, we
        // maintain the same 8% of equal distribution criteria.
        let background_threshold = {
            let expected_probe_barcode_count =
                probe_barcode_id_to_sample_id_and_mapped_probe_barcode_ids.len();
            if expected_probe_barcode_count <= 16 {
                0.005
            } else {
                0.08 / expected_probe_barcode_count as f64
            }
        };

        for (probe_barcode_id, umi_count) in umi_per_probe_barcode {
            let (sample_id, mapped_barcode_ids) =
                probe_barcode_id_to_sample_id_and_mapped_probe_barcode_ids
                    .get(probe_barcode_id)
                    .copied()
                    .unzip();

            // Elide probe barcodes not assigned to a sample and with few UMI.
            let umi_frac = umi_count as f64 / umi_sum as f64;
            if sample_id.is_none() && umi_frac < background_threshold {
                continue;
            }

            let barcode_ids: Vec<_> = std::iter::once(probe_barcode_id)
                .chain(mapped_barcode_ids.into_iter().flatten().map(String::as_str))
                .collect();

            let cells = *filtered_barcodes_per_probe_barcode
                .get(probe_barcode_id)
                .unwrap_or(&0);
            let cells_per_probe_barcode = if sample_id.is_some() {
                Some(cells)
            } else {
                assert_eq!(cells, 0);
                None
            };

            #[derive(JsonReport)]
            struct Metrics {
                probe_barcode_id: String,
                sample_id: Option<String>,
                umi_per_probe_barcode: usize,
                cells_per_probe_barcode: Option<usize>,
            }

            // Unpack probe barcode and sample IDs, and collect missing or unexpected
            // probe barcodes.

            let (probe_barcode_id, sample_id) = match sample_id {
                Some(sample_id) => {
                    let probe_barcode_ids_str = barcode_ids.join(PROBE_BARCODE_ID_GROUPING);
                    if umi_frac < background_threshold {
                        specified_probe_barcodes_missing.push(probe_barcode_ids_str.clone());
                    }
                    (probe_barcode_ids_str, Some(sample_id.to_string()))
                }
                None => {
                    unspecified_probe_barcodes_detected.push(probe_barcode_id.to_string());
                    (probe_barcode_id.to_string(), None)
                }
            };

            metrics.extend(special_metrics_proc.process_group(
                "rtl_probe_barcode_metrics",
                section,
                active_conditions,
                &self.alert_context,
                &Metrics {
                    probe_barcode_id,
                    sample_id,
                    umi_per_probe_barcode: umi_count,
                    cells_per_probe_barcode,
                },
            )?);
        }
        Ok(PerProbeBarcodeData {
            metrics,
            unexpected: unspecified_probe_barcodes_detected,
            missing: specified_probe_barcodes_missing,
            unexpected_missing_threshold: background_threshold,
        })
    }

    fn active_metric_conditions(&self, section: Section) -> ActiveConditions {
        let has_gdna = !matches!(
            self.lib_metrics.get("estimated_gdna_content"),
            None | Some(Value::Null)
        );
        if has_gdna {
            assert!(self.is_rtl);
        }

        ActiveConditions {
            section,
            vdj_receptor: None,
            is_multiplexed: self.multiplexing_method.is_some(),
            is_cell_multiplexed: self.is_cmo_multiplexed() || self.is_hashtag_multiplexed(),
            is_read_multiplexed: self.is_oh_multiplexed() || self.is_rtl_multiplexed(),
            is_rtl: self.is_rtl,
            has_gdna,
            include_introns: self.alert_context.include_introns,
            has_vdj_reference: false,
        }
    }

    fn add_cell_metrics(
        &self,
        section: Section,
        library_type: LibraryType,
        physical_library_id: &str,
        metrics_out: &mut Vec<JsonMetricSummary>,
    ) -> Result<()> {
        if let Some(metric) =
            self.get_high_occupancy_gem_metric(section, library_type, physical_library_id)?
        {
            metrics_out.push(metric);
        }
        Ok(())
    }

    /// Return a clone of the library metrics, with the provided physical library ID injected.
    fn get_metrics_with_physical_library_id(
        &self,
        physical_library_id: &str,
    ) -> TxHashMap<String, Value> {
        // FIXME CELLRANGER-8444 pick a better way to inject this than copying
        let mut metrics = self.lib_metrics.clone();
        metrics.insert(
            "physical_library_id".to_string(),
            json!(physical_library_id),
        );
        metrics
    }

    fn build_gex_ws(
        &self,
        physical_library_id: &str,
        contents: &mut LibWsBuilderContents,
    ) -> Result<LibraryGexWebSummary> {
        let library_type = LibraryType::Gex;
        let section = Section::Gex;
        let lib_metrics = self.get_metrics_with_physical_library_id(physical_library_id);

        let active_conditions = self.active_metric_conditions(section);

        let mut metrics = self.lib_metrics_proc.process(
            section,
            &active_conditions,
            &self.alert_context,
            &lib_metrics,
        )?;

        self.add_cell_metrics(section, library_type, physical_library_id, &mut metrics)?;

        let mut per_probe_barcode_data =
            self.build_rtl_probe_barcode_metrics(library_type, section, &active_conditions)?;
        metrics.append(&mut per_probe_barcode_data.metrics);

        metrics.extend(self.build_ocm_per_overhang_metrics(
            library_type,
            section,
            &active_conditions,
        )?);

        metrics.extend(self.get_sequencing_metrics(
            section,
            library_type,
            &active_conditions,
            contents,
        )?);

        let barcode_rank_plot = format_barcode_rank_plot(
            contents
                .barcode_rank_plots
                .remove(&LibraryType::Gex)
                .expect("Gene expression barcode rank plot was missing."),
            "GEX",
        );

        let barnyard_biplot = contents
            .barnyard_biplot
            .take()
            .map(|plot| format_barnyard_biplot(plot, "GEX"));

        // Construct final websummary
        Ok(LibraryGexWebSummary {
            parameters_table: self.count_param_table(library_type, per_probe_barcode_data)?,
            sequencing_saturation_plot: library_sequencing_saturation_plot_from_metrics(
                library_type,
                &self.lib_metrics,
            ),
            median_genes_per_cell_plot: library_median_genes_plot_from_metrics(
                &self.lib_metrics,
                self.genomes().iter().cloned().collect(),
                PlotType::LibraryPlot,
                self.targeting_method,
            ),
            barcode_rank_plot,
            barnyard_biplot,
            metrics: MetricsTraitWrapper(metrics),
        })
    }

    fn build_antibody_or_antigen_ws(
        &self,
        physical_library_id: &str,
        // true for antibody, false for antigen
        is_antibody: bool,
        contents: &mut LibWsBuilderContents,
    ) -> Result<LibraryAntibodyOrAntigenWebSummary> {
        let (section, library_type) = if is_antibody {
            (Section::Antibody, LibraryType::Antibody)
        } else {
            (Section::Antigen, LibraryType::Antigen)
        };

        let lib_metrics = self.get_metrics_with_physical_library_id(physical_library_id);

        let active_conditions = self.active_metric_conditions(section);

        let mut metrics = self.lib_metrics_proc.process(
            section,
            &active_conditions,
            &self.alert_context,
            &lib_metrics,
        )?;

        self.add_cell_metrics(section, library_type, physical_library_id, &mut metrics)?;

        let mut per_probe_barcode_data =
            self.build_rtl_probe_barcode_metrics(library_type, section, &active_conditions)?;
        metrics.append(&mut per_probe_barcode_data.metrics);

        metrics.extend(self.get_sequencing_metrics(
            section,
            library_type,
            &active_conditions,
            contents,
        )?);

        let barcode_rank_plot = format_barcode_rank_plot(
            contents
                .barcode_rank_plots
                .remove(&library_type)
                .unwrap_or_else(|| panic!("{library_type} barcode rank plot was missing.")),
            if is_antibody { "Antibody" } else { "Antigen" },
        );

        // Only show the antibody histogram at the library level for CMO and HASHTAG.
        let feature_histogram = match (
            is_antibody,
            (self.is_cmo_multiplexed() || self.is_hashtag_multiplexed()),
        ) {
            (true, true) => contents.antibody_histograms.take(),
            (true, false) => None,
            (false, _) => contents.antigen_histograms.take(),
        };

        Ok(LibraryAntibodyOrAntigenWebSummary {
            parameters_table: self.count_param_table(library_type, per_probe_barcode_data)?,
            barcode_rank_plot,
            sequencing_saturation_plot: library_sequencing_saturation_plot_from_metrics(
                library_type,
                &self.lib_metrics,
            ),
            feature_histogram,
            metrics: MetricsTraitWrapper(metrics),
        })
    }

    fn build_crispr_ws(
        &self,
        physical_library_id: &str,
        contents: &mut LibWsBuilderContents,
    ) -> Result<LibraryCrisprWebSummary> {
        let library_type = LibraryType::Crispr;
        let section = Section::Crispr;
        let lib_metrics = self.get_metrics_with_physical_library_id(physical_library_id);

        let active_conditions = self.active_metric_conditions(section);

        let mut metrics = self.lib_metrics_proc.process(
            section,
            &active_conditions,
            &self.alert_context,
            &lib_metrics,
        )?;

        self.add_cell_metrics(section, library_type, physical_library_id, &mut metrics)?;

        let mut per_probe_barcode_data =
            self.build_rtl_probe_barcode_metrics(library_type, section, &active_conditions)?;
        metrics.append(&mut per_probe_barcode_data.metrics);

        metrics.extend(self.get_sequencing_metrics(
            section,
            library_type,
            &active_conditions,
            contents,
        )?);

        Ok(LibraryCrisprWebSummary {
            parameters_table: self.count_param_table(library_type, per_probe_barcode_data)?,
            barcode_rank_plot: format_barcode_rank_plot(
                contents
                    .barcode_rank_plots
                    .remove(&LibraryType::Crispr)
                    .expect("CRISPR Guide Capture barcode rank plot was missing."),
                "CRISPR",
            ),
            sequencing_saturation_plot: library_sequencing_saturation_plot_from_metrics(
                library_type,
                &self.lib_metrics,
            ),
            metrics: MetricsTraitWrapper(metrics),
        })
    }

    fn build_custom_ws(
        &self,
        physical_library_id: &str,
        contents: &mut LibWsBuilderContents,
    ) -> Result<LibraryCustomFeatureWebSummary> {
        let library_type = LibraryType::Custom;
        let section = Section::Custom;
        let lib_metrics = self.get_metrics_with_physical_library_id(physical_library_id);

        let active_conditions = self.active_metric_conditions(section);

        let mut metrics = self.lib_metrics_proc.process(
            section,
            &active_conditions,
            &self.alert_context,
            &lib_metrics,
        )?;

        self.add_cell_metrics(section, library_type, physical_library_id, &mut metrics)?;

        let per_probe_barcode_data =
            self.build_rtl_probe_barcode_metrics(library_type, section, &active_conditions)?;
        // NOTE: we do not currently show these metrics for Custom, is this an oversight?
        // metrics.extend(per_probe_barcode_data.metrics.drain(..));

        metrics.extend(self.get_sequencing_metrics(
            section,
            library_type,
            &active_conditions,
            contents,
        )?);

        Ok(LibraryCustomFeatureWebSummary {
            parameters_table: self.count_param_table(library_type, per_probe_barcode_data)?,
            barcode_rank_plot: format_barcode_rank_plot(
                contents
                    .barcode_rank_plots
                    .remove(&LibraryType::Custom)
                    .expect("Custom feature barcode rank plot was missing."),
                "Custom",
            ),
            sequencing_saturation_plot: library_sequencing_saturation_plot_from_metrics(
                library_type,
                &self.lib_metrics,
            ),
            metrics: MetricsTraitWrapper(metrics),
        })
    }

    fn build_cmo_ws(
        &self,
        physical_library_id: &str,
        contents: &mut LibWsBuilderContents,
    ) -> Result<(LibraryCmoWebSummary, MultiSharedResource)> {
        let library_type = LibraryType::Cellplex;
        let section = Section::Cmo;
        let lib_metrics = self.get_metrics_with_physical_library_id(physical_library_id);

        let active_conditions = self.active_metric_conditions(section);

        let mut metrics = self.lib_metrics_proc.process(
            section,
            &active_conditions,
            &self.alert_context,
            &lib_metrics,
        )?;

        metrics.extend(self.get_sequencing_metrics(
            section,
            library_type,
            &active_conditions,
            contents,
        )?);

        metrics.extend(self.build_per_tag_metrics(section, &active_conditions)?);

        let (jibes_biplot, jibes_histogram, resources) =
            match contents.jibes_biplot_histogram.take() {
                Some(jibes_biplot_histogram) => (
                    Some(format_jibes_biplots(jibes_biplot_histogram.biplot, "CMO")),
                    Some(format_histogram(jibes_biplot_histogram.histogram, "CMO")),
                    jibes_biplot_histogram.resources,
                ),
                None => (None, None, TxHashMap::default()),
            };

        let (umi_proj, tags_proj) = if let Some(plots) = contents.cmo_projection_plot.take() {
            (
                Some(format_umi_on_umap_plot(
                    plots.cmo_umi_projection_plot,
                    "Multiplexing Capture",
                    "UMAP Projection of Cells Colored by UMI Counts",
                )),
                Some(format_tags_on_umap_plot(
                    plots.cmo_tags_projection_plot,
                    "Multiplexing Capture",
                )),
            )
        } else {
            (None, None)
        };

        Ok((
            LibraryCmoWebSummary {
                parameters_table: self
                    .count_param_table(LibraryType::Cellplex, Default::default())?,
                barcode_rank_plot: format_barcode_rank_plot(
                    contents
                        .barcode_rank_plots
                        .remove(&LibraryType::Cellplex)
                        .expect("CMO barcode rank plot was missing."),
                    "CMO",
                ),
                jibes_biplot,
                jibes_histogram,
                cmo_umi_projection_plot: umi_proj,
                cmo_tags_projection_plot: tags_proj,
                metrics: MetricsTraitWrapper(metrics),
            },
            resources,
        ))
    }

    fn build_hashtag_ws(
        &self,
        contents: &mut LibWsBuilderContents,
    ) -> Result<(LibraryHashtagWebSummary, MultiSharedResource)> {
        let library_type = LibraryType::Antibody;
        let section = Section::Hashtag;
        let active_conditions = self.active_metric_conditions(section);

        let mut metrics = self.lib_metrics_proc.process(
            section,
            &active_conditions,
            &self.alert_context,
            &self.lib_metrics,
        )?;

        metrics.extend(self.build_per_tag_metrics(section, &active_conditions)?);

        let (jibes_biplot, jibes_histogram, resources) =
            match contents.jibes_biplot_histogram.take() {
                Some(jibes_biplot_histogram) => (
                    Some(format_jibes_biplots(
                        jibes_biplot_histogram.biplot,
                        "Hashtag",
                    )),
                    Some(format_histogram(
                        jibes_biplot_histogram.histogram,
                        "Hashtag",
                    )),
                    jibes_biplot_histogram.resources,
                ),
                None => (None, None, TxHashMap::default()),
            };

        let (umi_proj, tags_proj) = if let Some(plots) = contents.cmo_projection_plot.take() {
            (
                Some(format_umi_on_umap_plot(
                    plots.cmo_umi_projection_plot,
                    "Antibody Capture",
                    "UMAP Projection of Cells Colored by Hashtag UMI Counts",
                )),
                Some(format_tags_on_umap_plot(
                    plots.cmo_tags_projection_plot,
                    "Antibody Capture",
                )),
            )
        } else {
            (None, None)
        };

        Ok((
            LibraryHashtagWebSummary {
                parameters_table: self.count_param_table(library_type, Default::default())?,
                jibes_biplot,
                jibes_histogram,
                hashtag_umi_projection_plot: umi_proj,
                hashtag_tags_projection_plot: tags_proj,
                metrics: MetricsTraitWrapper(metrics),
            },
            resources,
        ))
    }

    fn build_per_tag_metrics(
        &self,
        section: Section,
        active_conditions: &ActiveConditions,
    ) -> Result<Vec<JsonMetricSummary>> {
        let mut special_metric_proc = self.special_metrics_proc.with_default_transformers();
        special_metric_proc.add_transformer(
            "CellsFraction",
            CountAndPercentTransformer::new(
                get_metric_usize(&self.lib_metrics, "total_singlets")?.unwrap(),
            ),
        );

        let group = match section {
            Section::Cmo => "cmo_per_tag_metrics",
            Section::Hashtag => "hashtag_per_tag_metrics",
            _ => panic!("invalid section to create per tag metrics: {section}"),
        };

        #[derive(JsonReport)]
        struct Metrics<'a> {
            tag_id: &'a str,
            sample_id: &'a str,
            /// tag_{tag_name}_frac_reads_in_cells
            /// absent if no cells assigned to tag
            tag_reads_in_cell_associated_partitions: Option<f64>,
            /// tag_{tag_name}_number_of_singlets
            singlets_assigned_to_tag: usize,
            /// snr_{tag_name}_jibes
            /// absent if no cells assigned to tag
            tag_signal_to_background_ratio: Option<f64>,
        }

        let mut metrics = vec![];
        for sample in &self.multi_graph.samples {
            for fingerprint in &sample.fingerprints {
                let Fingerprint::Tagged { tag_name, .. } = fingerprint else {
                    panic!("Unable to process Tag metrics group involving untagged fingerprint.");
                };

                let singlets_assigned_to_tag_key = format!("tag_{tag_name}_number_of_singlets");

                metrics.extend(
                    special_metric_proc.process_group(
                        group,
                        section,
                        active_conditions,
                        &self.alert_context,
                        &Metrics {
                            tag_id: tag_name.as_str(),
                            sample_id: sample.sample_id.as_str(),
                            tag_reads_in_cell_associated_partitions: get_metric_f64(
                                &self.lib_metrics,
                                &format!("tag_{tag_name}_frac_reads_in_cells"),
                            )?,
                            singlets_assigned_to_tag: get_metric_usize(
                                &self.lib_metrics,
                                &singlets_assigned_to_tag_key,
                            )?
                            .expect(&singlets_assigned_to_tag_key),
                            tag_signal_to_background_ratio: get_metric_f64(
                                &self.lib_metrics,
                                &format!("snr_{tag_name}_jibes"),
                            )?,
                        },
                    )?,
                );
            }
        }
        Ok(metrics)
    }

    /// Build sequencing metrics for a particular library type.
    fn get_sequencing_metrics(
        &self,
        section: Section,
        library_type: LibraryType,
        active_conditions: &ActiveConditions,
        contents: &mut LibWsBuilderContents,
    ) -> Result<Vec<JsonMetricSummary>> {
        let mut metrics = vec![];
        let Some(all_seq_metrics) = contents.sequencing_metrics.remove(&library_type) else {
            panic!("sequencing metrics for {library_type} not found")
        };
        for seq_metrics in all_seq_metrics {
            metrics.extend(self.special_metrics_proc.process_group(
                "sequencing_metrics",
                section,
                active_conditions,
                &self.alert_context,
                &seq_metrics,
            )?);
        }
        Ok(metrics)
    }
}

fn get_aligner_from_metrics(metrics: &TxHashMap<String, Value>) -> Option<AlignerParam> {
    get_metric_string(metrics, "alignment_aligner")
        .unwrap()
        .map(|v| v.parse())
        .transpose()
        .unwrap()
}

#[derive(Default)]
struct PerProbeBarcodeData {
    metrics: Vec<JsonMetricSummary>,
    unexpected: Vec<String>,
    missing: Vec<String>,
    /// The threshold used to determine if an expected probe barcode was missing,
    /// or if an unexpected probe barcode was present at sufficient level to
    /// merit a warning.
    unexpected_missing_threshold: f64,
}

// describes the various UMAP plots a single sample may have.
#[derive(Serialize, Deserialize, Clone, Default)]
pub struct SampleUmapPlots {
    // full UMAP/clustering/diffexp plots for gene expression run
    gex_diffexp_clustering_plots: Box<RawValue>,
    // full UMAP/clustering/diffexp plots for antibody-only case only
    antibody_diffexp_clustering_plots: Option<Box<RawValue>>,
    // the normal colored-umi UMAPs that feature barcode libraries get otherwise
    crispr_umi_on_umap: Box<RawValue>,
    antibody_umi_on_umap: Box<RawValue>,
    custom_umi_on_umap: Box<RawValue>,
}

// describes the various UMAP plots a multiplexing experiment may have.
#[derive(Serialize, Deserialize, Clone, Default)]
pub struct MultiplexingUmapPlots {
    cmo_umi_projection_plot: Box<RawValue>,
    cmo_tags_projection_plot: Box<RawValue>,
}

struct SampleWsBuilder<'a> {
    lib_ws_builder: &'a LibWsBuilder,
    sample: Sample,
    // Absent if VDJ-only analysis.
    sample_metrics: Option<TxHashMap<String, Value>>,
    sample_metrics_proc: &'a MetricsProcessor,
    sample_barcode_rank_plots: Option<TxHashMap<LibraryType, Box<RawValue>>>,
    sample_treemap_plots: Option<TxHashMap<LibraryType, Box<RawValue>>>,
    sample_projection_plots: SampleUmapPlots,
    sample_antibody_histograms: Option<Box<RawValue>>,
    csv_str: &'a str,
    is_barnyard: bool,
    diagnostics: &'a MultiDiagnostics,
    pipeline_version: String,
    alert_context: &'a AlertContext,
    multiplexing_method: Option<BarcodeMultiplexingType>,
    targeting_method: Option<TargetingMethod>,
    antigen_vdj_metrics: Option<AntigenVdjMetrics>,
    clonotype_clustermap: Option<ChartWithHelp>,
    vdj_t_contents: Option<VdjWsContents>,
    vdj_t_gd_contents: Option<VdjWsContents>,
    vdj_b_contents: Option<VdjWsContents>,
    cell_annotation_websummary_bundle: Option<CellTypeWebSummaryBundleValue>,
    cell_annotation_viable_but_not_requested: bool,
}

impl SampleWsBuilder<'_> {
    fn build<'a>(
        mut self,
        library_websummary: &'a MultiWebSummaryLibraryData,
        resources: &'a MultiSharedResource,
    ) -> Result<MultiWebSummary<'a>> {
        let unpack_vdj = |contents: Option<VdjWsContents>| {
            let contents = contents?;
            let diagnostics = contents.diagnostics();
            let tab = Tab::new(contents.into_sample_ws(), self.alert_context);
            Some((tab, diagnostics))
        };

        let mut sample_ws = SampleWebSummary {
            id: self.sample.sample_id.clone(),
            description: self.sample.description.clone(),
            multiplexing_barcode_ids: self
                .sample
                .fingerprints
                .iter()
                .flat_map(Fingerprint::tag_names)
                .map(ToString::to_string)
                .collect(),
            ..Default::default()
        };

        let mut sample_diagnostics = SampleDiagnostics::default();
        if let Some((tab, diagnostics)) = unpack_vdj(self.vdj_t_contents.take()) {
            sample_diagnostics.vdj_t = Some(diagnostics);
            sample_ws.vdj_t_tab = Some(tab);
        }
        if let Some((tab, diagnostics)) = unpack_vdj(self.vdj_b_contents.take()) {
            sample_diagnostics.vdj_b = Some(diagnostics);
            sample_ws.vdj_b_tab = Some(tab);
        }
        if let Some((tab, diagnostics)) = unpack_vdj(self.vdj_t_gd_contents.take()) {
            sample_diagnostics.vdj_t_gd = Some(diagnostics);
            sample_ws.vdj_t_gd_tab = Some(tab);
        }

        for lib in &self.lib_ws_builder.multi_graph.libraries {
            match lib.library_type {
                LibraryType::Gex => {
                    assert!(sample_ws.gex_tab.is_none());

                    sample_ws.gex_tab = Some(Tab::new(
                        self.build_gex_ws(&lib.physical_library_id)?,
                        self.alert_context,
                    ));
                }
                LibraryType::Antibody => {
                    assert!(sample_ws.antibody_tab.is_none());
                    sample_ws.antibody_tab = Some(Tab::new(
                        self.build_antibody_ws(&lib.physical_library_id)?,
                        self.alert_context,
                    ));
                }
                LibraryType::Antigen => {
                    assert!(sample_ws.antigen_tab.is_none());
                    sample_ws.antigen_tab = Some(Tab::new(
                        self.build_antigen_ws(&lib.physical_library_id, &library_websummary.data)?,
                        self.alert_context,
                    ));
                }
                LibraryType::Crispr => {
                    assert!(sample_ws.crispr_tab.is_none());
                    sample_ws.crispr_tab = Some(Tab::new(
                        self.build_crispr_ws(&lib.physical_library_id)?,
                        self.alert_context,
                    ));
                }
                LibraryType::Vdj(_) => {}
                LibraryType::Custom => {
                    assert!(sample_ws.custom_feature_tab.is_none());
                    sample_ws.custom_feature_tab = Some(Tab::new(
                        self.build_custom_ws(&lib.physical_library_id)?,
                        self.alert_context,
                    ));
                }
                LibraryType::Cellplex => {}
                LibraryType::Atac => unreachable!(),
            }
        }

        let include_cell_annotation_tab = self
            .cell_annotation_websummary_bundle
            .as_ref()
            .is_some_and(CellTypeWebSummaryBundleValue::has_any_value);
        if include_cell_annotation_tab {
            sample_ws.cell_annotation_tab = Some(Tab::new(
                self.build_cell_annotation_ws()?,
                self.alert_context,
            ));
        }

        Ok(MultiWebSummary {
            page_title: format!("{} - {}", self.sample.sample_id, self.sample.description),
            pipeline_version: self.pipeline_version.clone(),
            library: library_websummary,
            per_sample: vec![MultiWebSummarySampleData {
                metrics: sample_ws.to_json_summary(),
                data: sample_ws,
            }],
            experimental_design: ExperimentalDesign {
                csv: self.csv_str.to_string(),
                multiplexing_method: self.multiplexing_method,
                is_rtl: self.lib_ws_builder.is_rtl,
                is_barnyard: self.is_barnyard,
                include_introns: self.alert_context.include_introns,
            },
            diagnostics: self.diagnostics.clone(),
            sample_diagnostics: vec![sample_diagnostics],
            resources,
        })
    }

    fn build_gex_ws(&mut self, physical_library_id: &str) -> Result<SampleGexWebSummary> {
        let section = Section::Gex;

        let sample_metrics = MetricsWithPhysicalLibraryId::wrap(
            self.sample_metrics.as_mut().unwrap(),
            physical_library_id,
        );

        let active_conditions = self.lib_ws_builder.active_metric_conditions(section);

        let mut metrics = self.sample_metrics_proc.process(
            section,
            &active_conditions,
            self.alert_context,
            &sample_metrics,
        )?;

        let genomes = self.lib_ws_builder.genomes();

        let disclaimer_banner = if self.cell_annotation_viable_but_not_requested {
            Some(CELL_ANNOTATION_ADVERTISEMENT_STRING.to_string())
        } else {
            None
        };

        #[derive(JsonReport)]
        struct HeroMetrics<'a> {
            genome: Option<&'a str>,
            total_singlets: Option<usize>,
            median_genes_per_singlet: Option<f64>,
            total_genes_detected: Option<usize>,
            median_umi_per_singlet: Option<f64>,
            confidently_mapped_reads_in_cells: Option<f64>,
        }

        for genome in genomes {
            let (median_genes_per_singlet, total_genes_detected, median_umi_per_singlet) = {
                let (keys, prefix) = match self.targeting_method {
                    Some(TargetingMethod::TemplatedLigation) => {
                        // Use the filtered probe set metrics.
                        // Targeted metrics are prefixed by genome only if there are >=2 genomes.
                        (
                            [
                                "median_genes_per_cell_on_target",
                                "num_genes_detected_on_target",
                                "median_umis_per_cell_on_target",
                            ],
                            (genomes.len() >= 2).then_some(genome.as_str()),
                        )
                    }
                    // Non-targeted metrics are always prefixed by genome, even if there is only one genome.
                    _ => (
                        [
                            "filtered_bcs_median_unique_genes_detected",
                            "filtered_bcs_total_unique_genes_detected",
                            "filtered_bcs_median_counts",
                        ],
                        Some(genome.as_str()),
                    ),
                };
                let keys = keys.map(|key| join_metric_name(prefix, key));
                (
                    get_metric_f64(&sample_metrics, &keys[0]).unwrap(),
                    get_metric_usize(&sample_metrics, &keys[1]).unwrap(),
                    get_metric_f64(&sample_metrics, &keys[2]).unwrap(),
                )
            };

            let confidently_mapped_reads_in_cells = if self
                .multiplexing_method
                .is_some_and(|mm| mm.is_cmo() || mm.is_hashtag())
            {
                // CMO and HASHTAG multiplexing does not demultiplex reads outside of cells.
                None
            } else {
                get_metric_f64(
                    &sample_metrics,
                    &format!("{genome}_filtered_bcs_conf_mapped_barcoded_reads_cum_frac"),
                )?
            };

            metrics.extend(self.lib_ws_builder.special_metrics_proc.process_group(
                "gex_sample_hero_metrics",
                section,
                &active_conditions,
                self.alert_context,
                &HeroMetrics {
                    genome: (genomes.len() >= 2).then_some(genome),
                    total_singlets: get_metric_usize(
                        &sample_metrics,
                        &format!("{genome}_singlets_assigned_to_this_sample"),
                    )?,
                    median_genes_per_singlet,
                    total_genes_detected,
                    median_umi_per_singlet,
                    confidently_mapped_reads_in_cells,
                },
            )?);
        }

        let barcode_rank_plot = if self
            .multiplexing_method
            .is_some_and(|mm| mm.is_rtl() || mm.is_oh())
        {
            Some(format_barcode_rank_plot(
                self.sample_barcode_rank_plots
                    .as_mut()
                    .ok_or_else(|| anyhow!("per-sample barcode rank plots not loaded"))?
                    .remove(&LibraryType::Gex)
                    .expect("GEX sample barcode rank plot not found"),
                "GEX",
            ))
        } else {
            None
        };

        Ok(SampleGexWebSummary {
            median_genes_per_cell_plot: if self.multiplexing_method.is_some() {
                Some(sample_median_genes_plot_from_metrics(
                    &sample_metrics,
                    genomes.iter().cloned().collect(),
                    PlotType::SamplePlot,
                    self.targeting_method,
                ))
            } else {
                None
            },
            clustering_and_diffexp_plots: self
                .sample_projection_plots
                .gex_diffexp_clustering_plots
                .clone(), // FIXME should move instead of cloning
            barcode_rank_plot,
            disclaimer: disclaimer_banner,
            metrics: MetricsTraitWrapper(metrics),
        })
    }

    fn build_antibody_ws(&mut self, physical_library_id: &str) -> Result<SampleAntibodyWebSummary> {
        let section = Section::Antibody;

        // For a Flex sample with no antibody reads, we won't have any plots or metrics.
        // Add an info-level warning but otherwise show an empty tab.
        let antibody_reads = self.sample_metrics.as_ref().unwrap()["ANTIBODY_total_read_pairs"]
            .as_u64()
            .unwrap();

        if antibody_reads == 0
            && self.sample.barcode_multiplexing_type() == Some(BarcodeMultiplexingType::RTL)
        {
            // Issue an informational alert by populating just a single hero metric.
            let reads_metric = JsonMetricSummary {
                key: "number_of_reads_in_cells".to_string(),
                value: antibody_reads.into(),
                string_value: "0".to_string(),
                category: MetricTier::Cells,
                section,
                config: MetricEtlConfig {
                    header: "Number of reads in cells".to_string(),
                    json_key: None,
                    ty: "usize".to_string(),
                    transformer: None,
                    extract: ExtractMode::Required,
                    alerts: vec![],
                },
                grouping_key: None,
                grouping_header: None,
                alerts: vec![AlertSpec {
                    level: AlertLevel::Info,
                    title: "No Antibody Reads Assigned to Sample".to_string(),
                    formatted_value: "0".to_string(),
                    message: "No antibody library reads were assigned to this sample.".to_string(),
                }],
            };
            return Ok(SampleAntibodyWebSummary {
                metrics: MetricsTraitWrapper(vec![reads_metric]),
                antibody_treemap: None,
                barcode_rank_plot: None,
                clustering_and_diffexp_plots: None,
                projection_plot: None,
                feature_histogram: None,
            });
        }

        let sample_metrics = MetricsWithPhysicalLibraryId::wrap(
            self.sample_metrics.as_mut().unwrap(),
            physical_library_id,
        );

        let active_conditions = self.lib_ws_builder.active_metric_conditions(section);

        let metrics = self.sample_metrics_proc.process(
            section,
            &active_conditions,
            self.alert_context,
            &sample_metrics,
        )?;

        // if this is antibody_only case, we don't show the single umi-on-UMAP plot.
        let umi_on_projection_plot = if self
            .sample_projection_plots
            .antibody_diffexp_clustering_plots
            .is_none()
        {
            Some(format_umi_on_umap_plot(
                std::mem::take(&mut self.sample_projection_plots.antibody_umi_on_umap),
                "Antibody Capture",
                "UMAP Projection",
            ))
        } else {
            None
        };

        let antibody_treemap = self
            .sample_treemap_plots
            .as_ref()
            .and_then(|sample_treemap| sample_treemap.get(&LibraryType::Antibody).cloned());

        let barcode_rank_plot = if self
            .multiplexing_method
            .is_some_and(|mm| mm.is_rtl() || mm.is_oh())
        {
            Some(format_barcode_rank_plot(
                self.sample_barcode_rank_plots
                    .as_mut()
                    .ok_or_else(|| anyhow!("per-sample barcode rank plots not loaded"))?
                    .remove(&LibraryType::Antibody)
                    .unwrap(),
                "Antibody",
            ))
        } else {
            None
        };

        Ok(SampleAntibodyWebSummary {
            antibody_treemap,
            clustering_and_diffexp_plots: self
                .sample_projection_plots
                .antibody_diffexp_clustering_plots
                .take(),

            projection_plot: umi_on_projection_plot,
            barcode_rank_plot,
            feature_histogram: self.sample_antibody_histograms.clone(),
            metrics: MetricsTraitWrapper(metrics),
        })
    }

    fn build_antigen_ws(
        &mut self,
        physical_library_id: &str,
        library_websummary: &LibraryWebSummary,
    ) -> Result<SampleAntigenWebSummary> {
        let section = Section::Antigen;

        let mut sample_metrics = MetricsWithPhysicalLibraryId::wrap(
            self.sample_metrics.as_mut().unwrap(),
            physical_library_id,
        );
        sample_metrics.insert_temp(
            "feature_type".to_string(),
            Value::from(LibraryType::Gex.to_string()),
        );

        let active_conditions = self.lib_ws_builder.active_metric_conditions(section);

        let mut metrics = self.sample_metrics_proc.process(
            section,
            &active_conditions,
            self.alert_context,
            &sample_metrics,
        )?;

        if let Some(antigen_vdj_metrics) = self.antigen_vdj_metrics.as_ref() {
            #[allow(non_snake_case)]
            #[derive(JsonReport)]
            struct Metrics<'a> {
                feature_type: &'a str,
                ANTIGEN_multi_filtered_bcs: usize,
                ANTIGEN_multi_filtered_bcs_median_counts: f64,
                ANTIGEN_multi_usable_reads_per_filtered_bc: f64,
            }

            metrics.extend(
                self.sample_metrics_proc.process_group(
                    "antigen_sample_hero_metrics",
                    section,
                    &active_conditions,
                    self.alert_context,
                    &Metrics {
                        feature_type: vdj_tab_names(library_websummary)
                            .into_iter()
                            .flatten()
                            .exactly_one() // Exactly 1 VDJ library guaranteed for antigen runs
                            .unwrap(),
                        ANTIGEN_multi_filtered_bcs: antigen_vdj_metrics.num_cells as usize,
                        ANTIGEN_multi_filtered_bcs_median_counts: antigen_vdj_metrics
                            .median_umis_per_cell,
                        ANTIGEN_multi_usable_reads_per_filtered_bc: antigen_vdj_metrics
                            .mean_usable_reads_per_cell,
                    },
                )?,
            );
        }

        let antigen_treemap = self
            .sample_treemap_plots
            .as_ref()
            .and_then(|sample_treemap| sample_treemap.get(&LibraryType::Antigen).cloned());

        Ok(SampleAntigenWebSummary {
            antigen_treemap,
            clonotype_clustermap: self.clonotype_clustermap.clone(),
            metrics: MetricsTraitWrapper(metrics),
        })
    }

    fn build_crispr_ws(&mut self, physical_library_id: &str) -> Result<SampleCrisprWebSummary> {
        let section = Section::Crispr;

        let sample_metrics = MetricsWithPhysicalLibraryId::wrap(
            self.sample_metrics.as_mut().unwrap(),
            physical_library_id,
        );

        let active_conditions = self.lib_ws_builder.active_metric_conditions(section);

        let metrics = self.sample_metrics_proc.process(
            section,
            &active_conditions,
            self.alert_context,
            &sample_metrics,
        )?;

        let barcode_rank_plot = if self
            .multiplexing_method
            .is_some_and(|mm| mm.is_rtl() || mm.is_oh())
        {
            Some(format_barcode_rank_plot(
                self.sample_barcode_rank_plots
                    .as_mut()
                    .ok_or_else(|| anyhow!("per-sample barcode rank plots not loaded"))?
                    .remove(&LibraryType::Crispr)
                    .unwrap(),
                "CRISPR",
            ))
        } else {
            None
        };

        Ok(SampleCrisprWebSummary {
            projection_plot: Some(format_umi_on_umap_plot(
                std::mem::take(&mut self.sample_projection_plots.crispr_umi_on_umap),
                "CRISPR Guide Capture",
                "UMAP Projection",
            )),
            barcode_rank_plot,
            metrics: MetricsTraitWrapper(metrics),
        })
    }

    fn build_custom_ws(
        &mut self,
        physical_library_id: &str,
    ) -> Result<SampleCustomFeatureWebSummary> {
        let section = Section::Custom;

        let sample_metrics = MetricsWithPhysicalLibraryId::wrap(
            self.sample_metrics.as_mut().unwrap(),
            physical_library_id,
        );

        let active_conditions = self.lib_ws_builder.active_metric_conditions(section);

        let metrics = self.sample_metrics_proc.process(
            section,
            &active_conditions,
            self.alert_context,
            &sample_metrics,
        )?;

        let barcode_rank_plot = if self.multiplexing_method.is_some_and(|mm| mm.is_rtl()) {
            Some(format_barcode_rank_plot(
                self.sample_barcode_rank_plots
                    .as_mut()
                    .ok_or_else(|| anyhow!("per-sample barcode rank plots not loaded"))?
                    .remove(&LibraryType::Custom)
                    .unwrap(),
                "Custom Feature",
            ))
        } else {
            None
        };

        Ok(SampleCustomFeatureWebSummary {
            projection_plot: Some(format_umi_on_umap_plot(
                std::mem::take(&mut self.sample_projection_plots.custom_umi_on_umap),
                "Custom Feature",
                "UMAP Projection",
            )),
            barcode_rank_plot,
            metrics: MetricsTraitWrapper(metrics),
        })
    }

    fn build_cell_annotation_ws(&mut self) -> Result<SampleCellAnnotationWebSummary> {
        let cell_annotation_metric_struct = self
            .cell_annotation_websummary_bundle
            .as_mut()
            .and_then(|x| x.cell_annotation_metrics.take());
        Ok(SampleCellAnnotationWebSummary {
            cas_success: cell_annotation_metric_struct
                .as_ref()
                .and_then(|x| x.cell_annotation_success),
            cell_annotation_disable_differential_expression: cell_annotation_metric_struct
                .as_ref()
                .and_then(|x| x.cell_annotation_differential_expression),
            disclaimer: cell_annotation_metric_struct
                .as_ref()
                .and_then(CellAnnotationMetrics::generate_disclaimer_html_fragment),
            parameters_table: cell_annotation_metric_struct
                .map(generate_cell_type_parameter_table)
                .transpose()?,
            cell_annotation_cell_types_chart: self
                .cell_annotation_websummary_bundle
                .as_mut()
                .and_then(|x| {
                    x.cell_type_interactive_bar_chart
                        .take()
                        .map(generate_cell_type_barchart_from_value)
                })
                .transpose()?,
            cell_annotation_violin_plot_chart: self
                .cell_annotation_websummary_bundle
                .as_mut()
                .and_then(|x| {
                    x.cell_types_box_plot
                        .take()
                        .map(generate_cell_type_violin_plot_from_value)
                })
                .transpose()?,
            cell_annotation_umap_plot_chart: self
                .cell_annotation_websummary_bundle
                .as_mut()
                .and_then(|x| {
                    x.cell_types_umap_plot
                        .take()
                        .map(generate_cell_type_umap_plot_from_value)
                })
                .transpose()?,
            cell_annotation_diffexp_table: self
                .cell_annotation_websummary_bundle
                .as_mut()
                .and_then(|x| x.diffexp.take().map(generate_cell_type_diffexp_from_value))
                .transpose()?,
        })
    }
}

/// Wrap a collection of metrics, injecting the physical library ID provided.
/// When the type is dropped, the added metric is removed from the collection.
///
/// Additional metrics can also be added, and their keys will also be removed
/// when dropped.
struct MetricsWithPhysicalLibraryId<'a> {
    metrics: &'a mut TxHashMap<String, Value>,
    extra_keys: Vec<String>,
}

impl<'a> MetricsWithPhysicalLibraryId<'a> {
    const KEY: &'static str = "physical_library_id";
    fn wrap(metrics: &'a mut TxHashMap<String, Value>, physical_library_id: &str) -> Self {
        metrics.insert(Self::KEY.to_string(), physical_library_id.into());
        Self {
            metrics,
            extra_keys: vec![],
        }
    }

    fn insert_temp(&mut self, k: impl Into<String>, v: impl Into<Value>) {
        let k = k.into();
        self.extra_keys.push(k.clone());
        self.metrics.insert(k, v.into());
    }
}

impl Drop for MetricsWithPhysicalLibraryId<'_> {
    fn drop(&mut self) {
        self.metrics.remove(Self::KEY);
        for k in &self.extra_keys {
            self.metrics.remove(k);
        }
    }
}

impl std::ops::Deref for MetricsWithPhysicalLibraryId<'_> {
    type Target = TxHashMap<String, Value>;

    fn deref(&self) -> &Self::Target {
        self.metrics
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
                    bail!("JSON metric {key} had unexpected type {val:?}");
                }
            }
            _ => {
                bail!("JSON metric {key} had unexpected type {val:?}");
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
                    "JSON metric {key} was accessed as a String but was an Array or an Object."
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
    let Some(contaminants) = tag_contaminant_info.as_ref() else {
        return Default::default();
    };

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
                && used_tags.contains(t.0)
        })
        .map(|z| z.0.to_string())
        .collect()
}

#[make_mro(volatile = strict)]
impl MartianStage for WriteMultiWebSummary {
    type StageInputs = StageInputs;
    type StageOutputs = StageOutputs;
    type ChunkInputs = MartianVoid;
    type ChunkOutputs = MartianVoid;

    fn split(
        &self,
        args: Self::StageInputs,
        _rover: MartianRover,
    ) -> Result<StageDef<Self::ChunkInputs>> {
        let mem_gib = if let Some(library_metrics) = args.library_metrics {
            let library_metrics = library_metrics.read()?;
            let num_samples = library_metrics["samples_with_any_singlets"]
                .as_u64()
                .expect("samples_with_any_singlets");
            let num_cells = library_metrics["total_singlets"]
                .as_u64()
                .expect("total_singlets");
            let samples_gib = (27_000_000 * num_samples).div_ceil(1024 * 1024 * 1024);
            let cells_gib = (38_000 * num_cells).div_ceil(1024 * 1024 * 1024);
            let mem_gib = 8 + samples_gib + cells_gib;
            println!(
                "num_samples={num_samples},num_cells={num_cells},\
                 samples_gib={samples_gib},cells_gib={cells_gib},mem_gib={mem_gib}"
            );
            mem_gib
        } else {
            // VDJ only
            8
        };
        Ok(StageDef::with_join_resource(Resource::with_mem_gb(
            mem_gib as isize,
        )))
    }

    fn main(
        &self,
        _args: Self::StageInputs,
        _chunk_args: Self::ChunkInputs,
        _rover: MartianRover,
    ) -> Result<Self::ChunkOutputs> {
        unreachable!();
    }

    fn join(
        &self,
        args: Self::StageInputs,
        _chunk_defs: Vec<Self::ChunkInputs>,
        _chunk_outs: Vec<Self::ChunkOutputs>,
        rover: MartianRover,
    ) -> Result<Self::StageOutputs> {
        let multi_graph = args.multi_graph.read()?;

        let lib_metrics = args
            .library_metrics
            .as_ref()
            .map(FileTypeRead::read)
            .transpose()?
            .unwrap_or_default();

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

        let chemistries = if let Some(chemistry_defs) = &args.chemistry_defs {
            chemistry_defs.values().map(|x| x.name).collect()
        } else {
            HashSet::new()
        };

        // chemistry_defs is None for a VDJ-only analysis, which is a 5' chemistry.
        assert!(
            args.chemistry_defs.is_some()
                || args.vdj_b_contents.is_some()
                || args.vdj_t_contents.is_some()
                || args.vdj_t_gd_contents.is_some()
        );
        let is_fiveprime = args
            .chemistry_defs
            .as_ref()
            .is_none_or(|x| x.endedness() == Some(WhichEnd::FivePrime));

        let targeting_method = args.count_inputs.as_ref().and_then(|x| x.targeting_method);
        let is_rtl = args
            .chemistry_defs
            .as_ref()
            .and_then(ChemistryDefsExt::is_rtl)
            .or_else(|| {
                // Not enough info in the chemistries for a clear answer; try
                // to rummage in the pipeline inputs.
                // If we're using RTL multiplexing, then we know we're RTL.
                multi_graph
                    .barcode_multiplexing_type()
                    .map(|m| m == BarcodeMultiplexingType::ReadLevel(ReadLevel::RTL))
            })
            .unwrap_or_else(|| {
                // Still can't figure it out; realistically the reason we need
                // to branch on this is because we're going to try to unpack
                // metrics produced by probe alignment, so if we used hurtle as
                // the aligner, then we're definitely running RTL.
                // TODO: we could eliminate the rest of the conditions above
                // and solely branch on this condition, but we'll need to determine
                // if this is totally correct. We only hit this branch for
                // singleplex analyses using custom chemistries.
                get_aligner_from_metrics(&lib_metrics) == Some(AlignerParam::Hurtle)
            });

        let alert_context = AlertContext {
            is_rtl,
            is_arc_chemistry: chemistries.contains(&ChemistryName::ArcV1),
            library_types: multi_graph.library_types().collect(),
            is_fiveprime,
            multiplexing_method: multi_graph.barcode_multiplexing_type(),
            include_introns: args
                .count_inputs
                .as_ref()
                .is_some_and(|inputs| inputs.include_introns),
            no_preflight: args.no_preflight,
        };

        //******************************************************************************************
        // POPULATE LIBRARY WEBSUMMARIES
        //******************************************************************************************
        let (mut lib_metrics_proc, mut sample_metrics_proc) = load_metrics_etl()?;

        let cells_fraction_transformer = CountAndPercentTransformer::new(
            get_metric_usize(&lib_metrics, "total_cell_associated_partitions")?.unwrap_or_default(),
        );
        lib_metrics_proc.add_transformer("CellsFraction", cells_fraction_transformer.clone());
        sample_metrics_proc.add_transformer("CellsFraction", cells_fraction_transformer);

        let lib_ws_builder = LibWsBuilder {
            common_inputs: args.common_inputs.clone(),
            multi_graph: multi_graph.clone(),
            multi_config: args.multi_config.read()?,
            chemistry_defs: args.chemistry_defs.clone(),
            lib_metrics,
            lib_metrics_proc,
            alert_context: alert_context.clone(),
            special_metrics_proc: load_metrics_etl_special()?,

            count_inputs: args.count_inputs.clone(),
            count_cell_calling_config: args.count_cell_calling_config.clone(),

            target_set_name: args.target_set_name,
            multiplexing_method: multi_graph.barcode_multiplexing_type(),
            specificity_controls: args
                .feature_config
                .unwrap_or(FeatureConfig {
                    beam_mode: None,
                    specificity_controls: None,
                    functional_map: None,
                    hashtag_ids: None,
                })
                .specificity_controls,
            dropped_tags,
            probe_barcodes_high_gem_overlap,
            mismatched_probe_barcode_pairings: get_mismatched_probe_barcode_pairings(
                &args.multi_config,
                args.detected_probe_barcode_pairing.as_ref(),
            )?,
            is_rtl,
            targeting_method,
        };

        let lib_ws_builder_contents = LibWsBuilderContents {
            barcode_rank_plots: match args.barcode_rank_plots {
                Some(ref plots) => plots.read()?,
                None => TxHashMap::default(),
            },
            sequencing_metrics: match args.sequencing_metrics {
                Some(f) => f.read()?,
                None => TxHashMap::default(),
            },
            jibes_biplot_histogram: args.jibes_biplot_histogram.map(|j| j.read()).transpose()?,
            barnyard_biplot: args.barnyard_biplot.map(|j| j.read()).transpose()?,
            antibody_histograms: args.antibody_histograms.map(|j| j.read()).transpose()?,
            antigen_histograms: args.antigen_histograms.map(|j| j.read()).transpose()?,
            cmo_projection_plot: args.cmo_projection_plot.map(|j| j.read()).transpose()?,
            vdj_t_contents: args
                .vdj_t_contents
                .clone()
                .map(|f| {
                    f.values()
                        .next()
                        .expect("Missing VDJ-T web-summary contents")
                        .read()
                })
                .transpose()?,
            vdj_t_gd_contents: args
                .vdj_t_gd_contents
                .clone()
                .map(|f| {
                    f.values()
                        .next()
                        .expect("Missing VDJ-T-GD web-summary contents")
                        .read()
                })
                .transpose()?,
            vdj_b_contents: args
                .vdj_b_contents
                .clone()
                .map(|f| {
                    f.values()
                        .next()
                        .expect("Missing VDJ-B web-summary contents")
                        .read()
                })
                .transpose()?,
        };

        //******************************************************************************************
        // POPULATE DIAGNOSTIC DATA
        //******************************************************************************************
        let get_metric = |metric| lib_ws_builder.lib_metrics.get(metric).cloned();
        let unknown_feature_bcs: HashMap<String, Value> = lib_ws_builder
            .lib_metrics
            .iter()
            .filter(|(key, _)| key.contains("unknown_feature_bcs"))
            .map(|(key, value)| (key.clone(), value.clone()))
            .collect();

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
            unknown_feature_barcode_seqs: unknown_feature_bcs,
        };

        //******************************************************************************************
        // POPULATE SAMPLE WEBSUMMARIES
        //******************************************************************************************

        let (library_ws, resources) =
            lib_ws_builder.build(&alert_context, lib_ws_builder_contents)?;

        // maps sample assignment to websummary
        let mut sample_to_json = TxHashMap::default();
        let mut sample_to_web_summary = TxHashMap::default();
        let mut sample_to_metrics_csv = TxHashMap::default();
        let mut per_sample_ws_data_trimmed = TxHashMap::default();

        let dist = load_dist()?;

        let csv_str = read_to_string(&args.multi_config)?
            .replace("\r\n", "\n")
            .replace('\r', "\n");

        for full_sample in &multi_graph.samples {
            let sample = SampleAssignment::Assigned(full_sample.sample_id.clone());
            let builder = SampleWsBuilder {
                sample: full_sample.clone(),
                lib_ws_builder: &lib_ws_builder,
                sample_metrics: read_optional_file_from_map(&args.per_sample_metrics, &sample)
                    .context("sample metrics")?,
                sample_metrics_proc: &sample_metrics_proc,
                sample_projection_plots: read_optional_file_from_map(
                    &args.sample_projection_plots,
                    &sample,
                )?
                .unwrap_or_default(),
                sample_barcode_rank_plots: read_optional_file_from_map(
                    &args.sample_barcode_rank_plots,
                    &sample,
                )
                .context("sample barcode rank plots")?,
                sample_treemap_plots: args
                    .sample_treemap_plots
                    .as_ref()
                    .map(|files| read_optional_file_from_map(files, &sample))
                    .transpose()?
                    .flatten(),
                sample_antibody_histograms: args
                    .sample_antibody_histograms
                    .as_ref()
                    .map(|files| read_optional_file_from_map(files, &sample))
                    .transpose()?
                    .flatten(),
                csv_str: &csv_str,
                is_barnyard: is_barnyard(args.count_inputs.as_ref()),
                diagnostics: &diagnostics,
                pipeline_version: rover.pipelines_version(),
                alert_context: &alert_context,
                multiplexing_method: multi_graph.barcode_multiplexing_type(),
                targeting_method: args
                    .count_inputs
                    .as_ref()
                    .and_then(|inputs| inputs.targeting_method),
                antigen_vdj_metrics: args
                    .beam_analyzer
                    .as_ref()
                    .and_then(|outs| outs.get(&sample))
                    .and_then(|inner| inner.as_ref())
                    .and_then(|beam_analyzer| beam_analyzer.antigen_vdj_metrics_bin.as_ref())
                    .and_then(|file| file.read().ok()),
                clonotype_clustermap: args
                    .beam_analyzer
                    .as_ref()
                    .and_then(|outs| outs.get(&sample))
                    .and_then(|inner| inner.as_ref())
                    .and_then(|beam_analyzer| beam_analyzer.antigen_specificity_scores.as_ref())
                    .and_then(|file| clonotype_specificity_heatmap(file.clone()).unwrap()),
                vdj_t_contents: args
                    .vdj_t_contents
                    .as_ref()
                    .map(|files| files[&sample].read())
                    .transpose()?,
                vdj_t_gd_contents: args
                    .vdj_t_gd_contents
                    .as_ref()
                    .map(|files| files[&sample].read())
                    .transpose()?,
                vdj_b_contents: args
                    .vdj_b_contents
                    .as_ref()
                    .map(|files| files[&sample].read())
                    .transpose()?,
                cell_annotation_websummary_bundle: args
                    .cell_type_websummary_bundles
                    .as_ref()
                    .map(|files| {
                        files
                            .get(&sample)
                            .and_then(|x| x.as_ref().map(CellTypeWebSummaryBundle::read))
                            .transpose()
                    })
                    .transpose()?
                    .flatten(),
                cell_annotation_viable_but_not_requested: args
                    .cell_annotation_viable_but_not_requested
                    .as_ref()
                    .is_some_and(|vals| vals.get(&sample).copied().flatten().unwrap_or_default()),
            };

            let mut sample_ws = builder.build(&library_ws, &resources)?;

            let json_file: JsonFile<()> =
                rover.make_path(format!("{sample}_web_summary_data.json"));

            serde_json::to_writer_pretty(json_file.buf_writer()?, &sample_ws)?;
            sample_to_json.insert(sample.clone(), json_file);

            let html_file: HtmlFile = rover.make_path(format!("{sample}_web_summary.html"));
            write_web_summary_html(&sample_ws, &dist, &mut html_file.buf_writer()?)?;
            sample_to_web_summary.insert(sample.clone(), html_file);

            let csv_file: CsvFile<()> = rover.make_path(format!("{sample}_metrics_summary"));
            sample_ws.write_sample_csv(&csv_file, 0)?;
            sample_to_metrics_csv.insert(sample.clone(), csv_file);

            // Since we've already written out the websummary, we only use this
            // data to write out the library-level websummary, which will never
            // need the large plot data.
            sample_ws.per_sample[0].data.remove_large_data();
            // We do not include any data from cell annotation in the multi-
            // sample websummary.
            sample_ws.per_sample[0].data.cell_annotation_tab = None;
            per_sample_ws_data_trimmed.insert(sample, sample_ws);
        }

        // Write the QC library metrics CSV.
        let qc_library_metrics: CsvFile<()> = rover.make_path("qc_library_metrics");
        per_sample_ws_data_trimmed
            .values()
            .next()
            .unwrap()
            .write_library_metrics_csv(&qc_library_metrics)?;

        // Combine all per-sample websummaries into a library-level websummary.
        // If this is singleplex data or the analysis only has a single sample,
        // use the same file we already generated.
        // Write the QC sample metrics CSV.
        let qc_sample_metrics: CsvFile<()> = rover.make_path("qc_sample_metrics");
        let library_web_summary = if sample_to_web_summary.len() >= 2 {
            let combined_ws_data = combine_ws_data(&multi_graph, per_sample_ws_data_trimmed);
            combined_ws_data.write_wide_csv(&qc_sample_metrics)?;

            let html_file: HtmlFile = rover.make_path("web_summary");
            write_web_summary_html(&combined_ws_data, &dist, &mut html_file.buf_writer()?)?;
            html_file
        } else {
            assert_eq!(per_sample_ws_data_trimmed.len(), 1);
            per_sample_ws_data_trimmed
                .values()
                .next()
                .unwrap()
                .write_wide_csv(&qc_sample_metrics)?;

            sample_to_web_summary
                .values()
                .exactly_one()
                .unwrap()
                .clone()
        };

        Ok(StageOutputs {
            web_summary: sample_to_web_summary,
            web_summary_json: sample_to_json,
            metrics_summary_csv: sample_to_metrics_csv,
            qc_library_metrics,
            qc_sample_metrics,
            library_web_summary,
        })
    }
}

/// Format data into websummary HTML and write into the provided writer.
pub fn write_web_summary_html<T: Serialize>(
    data: &T,
    dist: &WsDist,
    mut w: impl Write,
) -> Result<()> {
    // TODO: stream this directly into the file instead of buffering in memory.
    let json_payload = sanitize_json_str(&serde_json::to_string(data)?);
    let (css_bundle, js_bundle) = (&dist.css, &dist.js);

    write!(
        w,
        r#"<!DOCTYPE html>
<html>
  <head>
    <meta charset="utf-8" />
    <style>{css_bundle}</style>
  </head>

  <body>

    <div id="root">

    <script type="text/javascript">
      const data = {json_payload}
    </script>
    <script>{js_bundle}</script>
  </body>
</html>
    "#
    )?;
    Ok(())
}

/// Escape certain characters in a JSON string to make it safe to embed in HTML.
///
/// We currently escape the same set of characters as golang:
/// `<`, `>`, `&`, U+2028 and U+2029 (line and paragraph separator characters).
fn sanitize_json_str(input: &str) -> String {
    let mut sanitized = String::with_capacity(input.len());

    for char in input.chars() {
        match char {
            '<' => sanitized.push_str(r"\u003c"),
            '>' => sanitized.push_str(r"\u003e"),
            '&' => sanitized.push_str(r"\u0026"),
            '\u{2028}' => sanitized.push_str(r"\u2028"),
            '\u{2029}' => sanitized.push_str(r"\u2029"),
            other => sanitized.push(other),
        }
    }
    sanitized
}

/// The built script and CSS assets.
pub struct WsDist {
    css: String,
    js: String,
}

/// Load the built WS script and CSS assets.
pub fn load_dist() -> Result<WsDist> {
    /// Load a named file from the websummary dist directory.
    fn load_dist_file(name: &str) -> Result<String> {
        let file_path = {
            if let Ok(runfiles_dir) = bazel_utils::runfiles_dir() {
                runfiles_dir.join("cellranger/lib/typescript/websummary/dist/")
            } else {
                bazel_utils::current_exe()?
                    .parent()
                    .unwrap()
                    .parent()
                    .unwrap()
                    .join("typescript/websummary/dist/")
            }
            .join(name)
        };
        read_to_string(&file_path).with_context(move || file_path.to_string_lossy().to_string())
    }

    Ok(WsDist {
        css: load_dist_file("tenx-websummary-styles.min.css")?,
        js: load_dist_file("tenx-websummary-script.min.js")?,
    })
}

/// Combine all websummary data for more than one sample into a single rolled-up entity.
fn combine_ws_data<'a>(
    multi_graph: &CrMultiGraph,
    mut per_sample_ws: TxHashMap<SampleAssignment, MultiWebSummary<'a>>,
) -> MultiWebSummary<'a> {
    assert!(per_sample_ws.len() > 1);
    let mut output = {
        let first_sample =
            &per_sample_ws[&SampleAssignment::Assigned(multi_graph.samples[0].sample_id.clone())];
        MultiWebSummary {
            page_title: "Cell Ranger multi QC report".to_string(),
            pipeline_version: first_sample.pipeline_version.clone(),
            library: first_sample.library,
            experimental_design: first_sample.experimental_design.clone(),
            diagnostics: first_sample.diagnostics.clone(),
            resources: first_sample.resources,
            per_sample: vec![],
            sample_diagnostics: vec![],
        }
    };
    for sample in multi_graph
        .samples
        .iter()
        .map(|s| SampleAssignment::Assigned(s.sample_id.clone()))
    {
        let sample_ws = per_sample_ws.remove(&sample).unwrap();
        assert_eq!(1, sample_ws.per_sample.len());
        output.per_sample.extend(sample_ws.per_sample);
        // TODO: decide if we actually want all of this in the library-level ws.
        output
            .sample_diagnostics
            .extend(sample_ws.sample_diagnostics);
    }
    // We should have consumed all of the data.
    assert!(per_sample_ws.is_empty());

    output
}

/// Read a single file from a mapping of optional files.
fn read_optional_file_from_map<F: FileTypeRead<T>, T, K: Clone + Hash + Eq>(
    files: &TxHashMap<K, Option<F>>,
    key: &K,
) -> Result<Option<T>> {
    files
        .get(key)
        .and_then(Option::as_ref)
        .map(FileTypeRead::read)
        .transpose()
}

/// Return any differences between configured and detected probe barcode pairings.
/// Return None if there are no differences, or if there is no detected or
/// configured pairing.
fn get_mismatched_probe_barcode_pairings(
    multi_config: &MultiConfigCsvFile,
    detected_probe_barcode_pairing: Option<&DetectedProbeBarcodePairingFile>,
) -> Result<Option<MismatchedProbeBarcodePairings>> {
    let Some(samples) = multi_config.read()?.samples else {
        return Ok(None);
    };
    let Some(detected_probe_barcode_pairing) = detected_probe_barcode_pairing else {
        return Ok(None);
    };
    let detected_probe_barcode_pairing: HashSet<_> = detected_probe_barcode_pairing
        .read()?
        .into_iter()
        .filter(|(_, source_bc)| {
            // CRISPR not yet supported in pairing detection
            categorize_rtl_multiplexing_barcode_id(source_bc).unwrap()
                != RTLMultiplexingBarcodeType::Crispr
        })
        .map(|(target_bc, source_bc)| format!("{target_bc}{PROBE_BARCODE_ID_GROUPING}{source_bc}"))
        .collect();

    let configured_probe_barcode_pairing: HashSet<_> = samples
        .get_translated_probe_barcodes()
        .into_iter()
        .filter(|(source_bc, _)| {
            // CRISPR not yet supported in pairing detection
            categorize_rtl_multiplexing_barcode_id(source_bc).unwrap()
                != RTLMultiplexingBarcodeType::Crispr
        })
        .map(|(source_bc, target_bc)| format!("{target_bc}{PROBE_BARCODE_ID_GROUPING}{source_bc}"))
        .collect();
    if configured_probe_barcode_pairing.is_empty()
        || configured_probe_barcode_pairing == detected_probe_barcode_pairing
    {
        return Ok(None);
    }
    Ok(Some(MismatchedProbeBarcodePairings::new(
        &configured_probe_barcode_pairing,
        &detected_probe_barcode_pairing,
    )))
}

/// Determine if this is a barnyard analysis.
fn is_barnyard(count_inputs: Option<&CountInputs>) -> bool {
    count_inputs.is_some_and(|x| {
        x.reference_info
            .as_ref()
            .is_some_and(|x| x.genomes.len() >= 2)
    })
}

#[cfg(test)]
mod test {
    use crate::stages::write_multi_web_summary::sanitize_json_str;

    #[test]
    fn test_sanitize_json_str() {
        let input = "hello<my>funky&chars\u{2028}\u{2029}";
        assert_eq!(
            r"hello\u003cmy\u003efunky\u0026chars\u2028\u2029",
            sanitize_json_str(input),
        );
    }
}
