//! Martian stage SUMMARIZE_VDJ_FILTERS

use crate::matrix::{gex_umi_counts_per_barcode, H5File};
use anyhow::Result;
use cr_lib::parquet_file::{ParquetFile, ParquetWriter, PerBarcodeFilter};
use enclone_core::barcode_fate::BarcodeFate;
use itertools::Itertools;
use martian::prelude::*;
use martian_derive::{make_mro, martian_filetype, MartianStruct};
use martian_filetypes::json_file::JsonFile;
use martian_filetypes::tabular_file::CsvFile;
use martian_filetypes::{FileTypeRead, FileTypeWrite, LazyFileTypeIO};
use metric::{AsMetricPrefix, TxHashMap};
use plotly::color::Rgba;
use plotly::common::{Marker, TextPosition};
use plotly::layout::{Axis, AxisType, BarMode, HoverMode};
use plotly::{Bar, Layout, Scatter};
use rand::rngs::StdRng;
use rand_distr::Distribution;
use serde::{Deserialize, Serialize};
use statrs::statistics::OrderStatistics;
use std::cmp::{Ordering, Reverse};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::path::{Path, PathBuf};
use strum_macros::IntoStaticStr;
use tenx_websummary::components::{HeroMetric, PlotlyChart, TitleWithHelp, WithTitle, WsNavBar};
use tenx_websummary::{HtmlTemplate, SinglePageHtml};
use vdj_ann::annotate::ContigAnnotation;
use vdj_asm_asm::write_ann_csv::ContigAnnotationCsvRow;
use vdj_filter_barcodes::filter_log::{FilterLogEntry, VdjFilterLogFormat};

struct Cdr3Info {
    cdr3_aa: String,
    num_umis: usize,
    chain: String,
}

#[derive(
    Clone, Copy, Hash, Eq, PartialEq, Ord, PartialOrd, IntoStaticStr, Serialize, Deserialize,
)]
#[serde(rename_all = "snake_case")]
#[strum(serialize_all = "snake_case")]
pub enum BarcodeCategory {
    Cell,
    GexFilter,
    AsmFilter,
    AsmGexFilter,
    EncloneFilter,
    EncloneGexFilter,
    UnknownFilter,
}

impl BarcodeCategory {
    pub fn label(self) -> String {
        match self {
            BarcodeCategory::Cell => "Cell barcodes",
            BarcodeCategory::GexFilter => "GEX filter",
            BarcodeCategory::AsmFilter => "Assembly filter",
            BarcodeCategory::EncloneFilter => "Enclone filter",
            BarcodeCategory::AsmGexFilter => "Assembly+GEX filter",
            BarcodeCategory::EncloneGexFilter => "Enclone+GEX filter",
            BarcodeCategory::UnknownFilter => "Unknown Filter",
        }
        .to_string()
    }

    pub fn color(self, alpha: Option<f64>) -> Rgba {
        let a = alpha.unwrap_or(1.0);
        match self {
            BarcodeCategory::Cell => Rgba::new(44, 160, 44, a), // cooked asparagus green
            BarcodeCategory::GexFilter => Rgba::new(255, 127, 14, a), // safety orange
            BarcodeCategory::AsmFilter => Rgba::new(31, 119, 180, a), // muted blue
            BarcodeCategory::EncloneFilter => Rgba::new(214, 39, 40, a), // brick red
            BarcodeCategory::AsmGexFilter => Rgba::new(148, 103, 189, a), // muted purple
            BarcodeCategory::EncloneGexFilter => Rgba::new(140, 86, 75, a), // chestnut brown
            BarcodeCategory::UnknownFilter => Rgba::new(127, 127, 127, a), // middle gray
        }
    }

    fn abbreviation(self) -> &'static str {
        match self {
            BarcodeCategory::Cell => "C",
            BarcodeCategory::GexFilter => "G",
            BarcodeCategory::AsmFilter => "A",
            BarcodeCategory::EncloneFilter => "E",
            BarcodeCategory::AsmGexFilter => "A+G",
            BarcodeCategory::EncloneGexFilter => "E+G",
            BarcodeCategory::UnknownFilter => "U",
        }
    }
}

impl AsMetricPrefix for BarcodeCategory {
    fn as_metric_prefix(&self) -> Option<&'static str> {
        Some(self.into())
    }
}

struct BarcodeInfo {
    is_cell: bool,
    is_gex_cell: Option<bool>,
    is_asm_cell: Option<bool>,
    enclone_fate: Option<BarcodeFate>,
    num_productive_contigs: usize,
    num_umis: usize, // Summed across all productive contigs
    cdr3s: Vec<Cdr3Info>,
}

impl BarcodeInfo {
    fn category(&self) -> BarcodeCategory {
        if self.is_cell {
            BarcodeCategory::Cell
        } else if self.is_asm_cell == Some(false) && self.is_gex_cell == Some(false) {
            BarcodeCategory::AsmGexFilter
        } else if self.is_asm_cell == Some(false) {
            BarcodeCategory::AsmFilter
        } else {
            match &self.enclone_fate {
                Some(fate) => match fate {
                    BarcodeFate::NotGexCell => BarcodeCategory::GexFilter,
                    _ => {
                        if self.is_gex_cell == Some(false) {
                            BarcodeCategory::EncloneGexFilter
                        } else {
                            BarcodeCategory::EncloneFilter
                        }
                    }
                },
                None => {
                    if self.is_gex_cell == Some(false) {
                        BarcodeCategory::GexFilter
                    } else {
                        BarcodeCategory::UnknownFilter
                    }
                }
            }
        }
    }
}

struct TreeMapPlotUnit {
    label: String,
    parent: Option<String>,
    value: usize,
    color: Rgba,
}

const BC_WITH_PROD_CONTIGS: &str = "Barcodes with productive contig";
const PROD_CONTIGS_UMI: &str = "Productive contig UMIs";

fn make_treemap(plot_data: Vec<TreeMapPlotUnit>) -> PlotlyChart {
    let (labels, parents, values, colors): (Vec<_>, Vec<_>, Vec<_>, Vec<_>) =
        itertools::multiunzip(plot_data.into_iter().map(
            |TreeMapPlotUnit {
                 label,
                 parent,
                 value,
                 color,
             }| (label, parent.unwrap_or_default(), value, color),
        ));
    PlotlyChart::with_layout_and_data(
        serde_json::json!({}),
        vec![serde_json::json!({
            "type": "treemap",
            "branchvalues": "total",
            "labels": labels,
            "parents": parents,
            "values": values,
            "textinfo": "label+value+percent root",
            "hovertemplate": "%{label}: %{value} (%{percentRoot:.1%})<extra></extra>",
            "marker": {
                "colors": colors,
            }
        })],
    )
}

#[derive(Serialize)]
pub struct VdjFilterMetrics {
    pub barcodes_with_productive_contig: i64,
    pub productive_contig_umis_in_cells_frac: f64,
    pub barcodes_with_productive_contig_per_filter_category: TxHashMap<BarcodeCategory, i64>,
    pub productive_contig_umis_per_filter_category: TxHashMap<BarcodeCategory, i64>,
    pub barcodes_with_productive_contig_per_asm_filter: Option<TxHashMap<String, i64>>,
    pub barcodes_with_productive_contig_per_enclone_filter: Option<TxHashMap<String, i64>>,
    pub unique_cdr3s_per_chain: HashMap<String, i64>,
}

#[derive(Serialize, Clone, HtmlTemplate)]
pub struct WebSummaryContent {
    umis_in_cells: Option<HeroMetric>,
    barcode_filter_summary: WithTitle<PlotlyChart>,
    umi_filter_summary: WithTitle<PlotlyChart>,
    productive_contigs_per_barcode: WithTitle<PlotlyChart>,
    unique_cdr3s_per_chain: WithTitle<PlotlyChart>,
    umi_histogram_categorized: WithTitle<PlotlyChart>,
    umi_histogram_all: WithTitle<PlotlyChart>,
    cdr3_frequency: WithTitle<PlotlyChart>,
    cdr3_frequency_umi: WithTitle<PlotlyChart>,
    cdr3_umi_distribution: WithTitle<PlotlyChart>,
    gex_vdj_umi_scatter: Option<WithTitle<PlotlyChart>>,
    #[html(row = "barcode_rank")]
    gex_barcode_rank_plot: Option<WithTitle<PlotlyChart>>,
    #[html(row = "barcode_rank")]
    vdj_barcode_rank_plot: WithTitle<PlotlyChart>,
}

pub struct FilterSummary {
    pub html: SinglePageHtml<WebSummaryContent>,
    pub metrics: VdjFilterMetrics,
    pub per_barcode_filters: Vec<PerBarcodeFilter>,
}

pub enum Annotations {
    Json {
        all_contig_json: PathBuf,
        enclone_barcode_fate: PathBuf,
    },
    Csv {
        all_contig_csv: PathBuf,
    },
}

impl Annotations {
    pub fn json(all_contig_json: &Path, enclone_barcode_fate: &Path) -> Self {
        Annotations::Json {
            all_contig_json: all_contig_json.to_path_buf(),
            enclone_barcode_fate: enclone_barcode_fate.to_path_buf(),
        }
    }
    pub fn csv(all_contig_csv: &Path) -> Self {
        Annotations::Csv {
            all_contig_csv: all_contig_csv.to_path_buf(),
        }
    }
    fn _load_enclone_barcode_fate(bc_fate: &Path) -> Result<HashMap<String, BarcodeFate>> {
        let filters_vec: Vec<HashMap<String, BarcodeFate>> = JsonFile::from(bc_fate).read()?;
        let mut filters: HashMap<String, BarcodeFate> = HashMap::new();
        for map in filters_vec {
            for (bc, fate) in map {
                filters.insert(bc, fate);
            }
        }
        Ok(filters)
    }

    fn _load_barcode_info_from_json(
        annotations_path: &Path,
        enclone_filters: &HashMap<String, BarcodeFate>,
    ) -> Result<HashMap<String, BarcodeInfo>> {
        let contig_reader = JsonFile::from(annotations_path).lazy_reader()?;
        let mut barcode_info = HashMap::new();
        for ann in contig_reader {
            let ann: ContigAnnotation = ann?;
            if ann.is_productive() {
                let enclone_fate = enclone_filters.get(&ann.barcode).copied();
                let chain_type = ann.chain_type().unwrap();
                let e = barcode_info.entry(ann.barcode).or_insert(BarcodeInfo {
                    is_cell: ann.is_cell,
                    is_gex_cell: ann.is_gex_cell,
                    is_asm_cell: ann.is_asm_cell,
                    enclone_fate,
                    num_productive_contigs: 0,
                    num_umis: 0,
                    cdr3s: vec![],
                });
                e.num_productive_contigs += 1;
                e.num_umis += ann.umi_count;
                e.cdr3s.push(Cdr3Info {
                    cdr3_aa: format!("{chain_type}:{}", ann.cdr3.unwrap()),
                    num_umis: ann.umi_count,
                    chain: chain_type.to_string(),
                });
            }
        }
        Ok(barcode_info)
    }

    fn _load_barcode_info_from_csv(all_contig_csv: &Path) -> Result<HashMap<String, BarcodeInfo>> {
        let mut barcode_info = HashMap::new();
        let reader = CsvFile::from(all_contig_csv).lazy_reader()?;
        for ann in reader {
            let ann: ContigAnnotationCsvRow = ann?;
            if ann.productive {
                let e = barcode_info.entry(ann.barcode).or_insert(BarcodeInfo {
                    is_cell: ann.is_cell,
                    is_gex_cell: None,
                    is_asm_cell: None,
                    enclone_fate: None,
                    num_productive_contigs: 0,
                    num_umis: 0,
                    cdr3s: vec![],
                });
                e.num_productive_contigs += 1;
                e.num_umis += ann.umis;
                e.cdr3s.push(Cdr3Info {
                    cdr3_aa: format!("{}:{}", ann.chain, ann.cdr3.unwrap()),
                    num_umis: ann.umis,
                    chain: ann.chain.to_string(),
                });
            }
        }
        Ok(barcode_info)
    }
    #[allow(clippy::type_complexity)]
    fn load(
        &self,
    ) -> Result<
        (
            HashMap<String, BarcodeInfo>,
            Option<HashMap<String, BarcodeFate>>,
        ),
        Error,
    > {
        match self {
            Self::Json {
                all_contig_json,
                enclone_barcode_fate,
            } => {
                let enclone_filters = Self::_load_enclone_barcode_fate(enclone_barcode_fate)?;
                let barcode_info =
                    Self::_load_barcode_info_from_json(all_contig_json, &enclone_filters)?;
                Ok((barcode_info, Some(enclone_filters)))
            }
            Self::Csv { all_contig_csv } => {
                Ok((Self::_load_barcode_info_from_csv(all_contig_csv)?, None))
            }
        }
    }
}

pub fn generate_filter_summary(
    analysis_id: &str,
    description: Option<&str>,
    annotations: Annotations,
    filter_diagnostics_path_opt: Option<&VdjFilterLogFormat>,
    raw_matrix_h5: Option<&H5File>,
) -> Result<FilterSummary> {
    // Keys are valid barcodes across the whole library
    let asm_filters = match filter_diagnostics_path_opt {
        Some(filter_diagnostics_path) => {
            let filter_reader = filter_diagnostics_path.lazy_reader()?;
            let mut asm_filters = HashMap::new();
            for log_entry in filter_reader {
                let log_entry = log_entry?;
                let FilterLogEntry::CellCalling { barcode, filter } = log_entry;

                asm_filters
                    .entry(barcode)
                    .or_insert_with(Vec::new)
                    .push(filter.label());
            }
            Some(asm_filters)
        }
        None => None,
    };

    // Keys are valid barcodes restricted to this sample
    let (barcode_info, enclone_filters) = annotations.load()?;

    let vdj_filters_per_barcode = asm_filters
        .as_ref()
        .map_or(Vec::<PerBarcodeFilter>::new(), |asm_filters| {
            make_vdj_filters_per_barcode(asm_filters, &barcode_info)
        });

    // Match definition of BarcodeCategory::AsmFilter i.e. barcodes filtered only by the assembler
    let barcodes_per_asm_filter: Option<TxHashMap<_, _>> = asm_filters.as_ref().map(|filters| {
        filters
            .iter()
            .filter_map(|(bc, f)| {
                if let Some(info) = barcode_info.get(bc) {
                    if info.is_gex_cell.unwrap_or(true) {
                        return Some(f);
                    }
                }
                None
            })
            .flatten()
            .counts()
            .into_iter()
            .map(|(&k, v)| (k.to_string(), v as i64))
            .collect()
    });

    let barcodes_per_enclone_filter: Option<TxHashMap<_, _>> =
        enclone_filters.as_ref().map(|filters| {
            filters
                .values()
                .map(BarcodeFate::label)
                .counts()
                .into_iter()
                .map(|(k, v)| (k.to_string(), v as i64))
                .collect()
        });

    let mut barcode_treemap_data = vec![TreeMapPlotUnit {
        label: BC_WITH_PROD_CONTIGS.into(),
        parent: None,
        value: barcode_info.len(),
        color: Rgba::new(255, 255, 255, 1.0),
    }];

    let mut umi_treemap_data = vec![TreeMapPlotUnit {
        label: PROD_CONTIGS_UMI.into(),
        parent: None,
        value: barcode_info.values().map(|i| i.num_umis).sum(),
        color: Rgba::new(255, 255, 255, 1.0),
    }];

    let mut categorized_barcodes = BTreeMap::new();
    for (bc, info) in &barcode_info {
        categorized_barcodes
            .entry(info.category())
            .or_insert_with(Vec::new)
            .push(bc.as_str());
    }
    for (category, barcodes) in &categorized_barcodes {
        barcode_treemap_data.push(TreeMapPlotUnit {
            label: category.label(),
            parent: Some(BC_WITH_PROD_CONTIGS.into()),
            value: barcodes.len(),
            color: category.color(Some(0.85)),
        });

        umi_treemap_data.push(TreeMapPlotUnit {
            label: category.label(),
            parent: Some(PROD_CONTIGS_UMI.into()),
            value: barcodes.iter().map(|bc| barcode_info[*bc].num_umis).sum(),
            color: category.color(Some(0.85)),
        });
        // Sub-categorize asm-filters
        match category {
            BarcodeCategory::AsmFilter | BarcodeCategory::AsmGexFilter => {
                if let Some(filters) = &asm_filters {
                    for (&filter_label, filter_barcodes) in barcodes
                        .iter()
                        // Take the first filter which the barcode failed
                        .into_group_map_by(|&&barcode| &filters[barcode][0])
                    {
                        barcode_treemap_data.push(TreeMapPlotUnit {
                            label: format!("{filter_label} [{}]", category.abbreviation()),
                            parent: Some(category.label()),
                            value: filter_barcodes.len(),
                            color: category.color(Some(1.0)),
                        });
                        umi_treemap_data.push(TreeMapPlotUnit {
                            label: format!("{filter_label} [{}]", category.abbreviation()),
                            parent: Some(category.label()),
                            value: filter_barcodes
                                .iter()
                                .map(|bc| barcode_info[**bc].num_umis)
                                .sum(),
                            color: category.color(Some(1.0)),
                        });
                    }
                }
            }
            BarcodeCategory::EncloneFilter | BarcodeCategory::EncloneGexFilter => {
                for (filter_label, filter_barcodes) in
                    barcodes.iter().into_group_map_by(|&&barcode| {
                        enclone_filters
                            .as_ref()
                            .unwrap()
                            .get(barcode)
                            .map_or("UNKNOWN", BarcodeFate::label)
                    })
                {
                    barcode_treemap_data.push(TreeMapPlotUnit {
                        label: format!("{filter_label} [{}]", category.abbreviation()),
                        parent: Some(category.label()),
                        value: filter_barcodes.len(),
                        color: category.color(Some(1.0)),
                    });
                    umi_treemap_data.push(TreeMapPlotUnit {
                        label: format!("{filter_label} [{}]", category.abbreviation()),
                        parent: Some(category.label()),
                        value: filter_barcodes
                            .iter()
                            .map(|bc| barcode_info[**bc].num_umis)
                            .sum(),
                        color: category.color(Some(1.0)),
                    });
                }
            }
            BarcodeCategory::Cell | BarcodeCategory::GexFilter | BarcodeCategory::UnknownFilter => {
            }
        }
    }
    let (cdr3_freq_plot, cdr3_freq_umi_plot, cdr3_umi_plot) =
        make_cdr3_plots(&categorized_barcodes, &barcode_info);

    let gex_umis = match raw_matrix_h5 {
        Some(h5_path) => Some(gex_umi_counts_per_barcode(h5_path)?),
        None => None,
    };
    let gex_vdj_umi_scatter = match &gex_umis {
        Some(gex_umis) => {
            let sub_filters: HashMap<_, _> = asm_filters
                .iter()
                .flatten()
                .map(|(k, v)| (k.as_str(), v.iter().unique().join("<br>")))
                .chain(
                    enclone_filters
                        .iter()
                        .flatten()
                        .map(|(k, v)| (k.as_str(), v.label().to_string())),
                )
                .collect();
            Some(WithTitle {
                title: TitleWithHelp {
                    help: String::new(),
                    title: "GEX vs VDJ UMIs".into(),
                }
                .into(),
                inner: make_gex_vdj_umi_scatter(
                    gex_umis,
                    &categorized_barcodes,
                    &barcode_info,
                    &sub_filters,
                ),
            })
        }
        None => None,
    };

    let productive_contig_umis_in_cells_frac = barcode_info
        .values()
        .filter_map(|info| {
            if info.is_cell {
                Some(info.num_umis as i64)
            } else {
                None
            }
        })
        .sum::<i64>() as f64
        / barcode_info
            .values()
            .map(|info| info.num_umis as i64)
            .sum::<i64>() as f64;

    let unique_cdr3s_per_chain: HashMap<String, i64> = {
        let mut cdr3s_per_chain = HashMap::new();
        for Cdr3Info { cdr3_aa, chain, .. } in barcode_info.values().flat_map(|info| &info.cdr3s) {
            cdr3s_per_chain
                .entry(chain)
                .or_insert_with(HashSet::new)
                .insert(cdr3_aa);
        }
        cdr3s_per_chain
            .iter()
            .map(|(&k, v)| (k.clone(), v.len() as i64))
            .collect()
    };

    let content = WebSummaryContent {
        umis_in_cells: (!productive_contig_umis_in_cells_frac.is_nan()).then(|| HeroMetric::new("VDJ UMIs in cells", format!("{:.1}%", 100.0 * productive_contig_umis_in_cells_frac))),
        barcode_filter_summary: WithTitle {
            title: TitleWithHelp {
                help: "The plot shows a summary of where a barcode witha productive contig gets filtered in the pipeline".into(),
                title: "Barcode filter summary".into(),
            }.into(),
            inner: make_treemap(barcode_treemap_data),
        },
        umi_filter_summary: WithTitle {
            title: TitleWithHelp {
                help: "The plot shows a summary of where a productive contig UMIs gets filtered in the pipeline".into(),
                title: "UMI filter summary".into(),
            }.into(),
            inner: make_treemap(umi_treemap_data),
        },
        productive_contigs_per_barcode: WithTitle {
            title: TitleWithHelp {
                help: String::new(),
                title: "Productive contigs per barcode".into(),
            }.into(),
            inner: make_contigs_per_barcode_plot(&categorized_barcodes, &barcode_info),
        },
        unique_cdr3s_per_chain: WithTitle {
            title: TitleWithHelp {
                help: String::new(),
                title: "Unique CDR3s per chain".into(),
            }.into(),
            inner: make_unique_cdr3s_plot(&unique_cdr3s_per_chain),
        },
        umi_histogram_all: WithTitle {
            title: TitleWithHelp {
                help: String::new(),
                title: "VDJ UMI Histogram".into(),
            }.into(),
            inner: make_umi_histogram(&categorized_barcodes, &barcode_info, false),
        },
        umi_histogram_categorized: WithTitle {
            title: TitleWithHelp {
                help: String::new(),
                title: "VDJ UMI Histogram (Categorized)".into(),
            }.into(),
            inner: make_umi_histogram(&categorized_barcodes, &barcode_info, true),
        },
        cdr3_frequency: WithTitle {
            title: TitleWithHelp {
                help: String::new(),
                title: "Top CDR3 AA distribution".into(),
            }.into(),
            inner: cdr3_freq_plot,
        },
        cdr3_frequency_umi: WithTitle {
            title: TitleWithHelp {
                help: String::new(),
                title: "Top CDR3 AA distribution (UMI weighted)".into(),
            }.into(),
            inner: cdr3_freq_umi_plot,
        },
        cdr3_umi_distribution: WithTitle {
            title: TitleWithHelp {
                help: String::new(),
                title: "Top CDR3 UMI distribution".into(),
            }.into(),
            inner: cdr3_umi_plot,
        },
        gex_vdj_umi_scatter,
        gex_barcode_rank_plot: gex_umis.map(|g| {
            WithTitle {
            title: TitleWithHelp {
                help: String::new(),
                title: "Gex Barcode rank".into(),
            }.into(),
            inner: gex_barcode_rank_plot(&g),
        }}),
        vdj_barcode_rank_plot: WithTitle {
            title: TitleWithHelp {
                help: String::new(),
                title: "VDJ Barcode rank".into(),
            }.into(),
            inner: vdj_barcode_rank_plot(&barcode_info),
        },
    };

    Ok(FilterSummary {
        html: SinglePageHtml::new(
            WsNavBar {
                pipeline: "VDJ Barcode Filtering".to_string(),
                id: analysis_id.to_string(),
                description: description.unwrap_or_default().to_string(),
            },
            content,
            None,
        ),
        metrics: VdjFilterMetrics {
            barcodes_with_productive_contig: barcode_info.len() as i64,
            barcodes_with_productive_contig_per_filter_category: categorized_barcodes
                .iter()
                .map(|(category, barcodes)| (*category, barcodes.len() as i64))
                .collect(),
            barcodes_with_productive_contig_per_asm_filter: barcodes_per_asm_filter,
            barcodes_with_productive_contig_per_enclone_filter: barcodes_per_enclone_filter,
            productive_contig_umis_in_cells_frac,
            productive_contig_umis_per_filter_category: categorized_barcodes
                .iter()
                .map(|(category, barcodes)| {
                    (
                        *category,
                        barcodes
                            .iter()
                            .map(|bc| barcode_info[*bc].num_umis as i64)
                            .sum(),
                    )
                })
                .collect(),
            unique_cdr3s_per_chain,
        },
        per_barcode_filters: vdj_filters_per_barcode,
    })
}

fn make_unique_cdr3s_plot(unique_cdr3s_per_chain: &HashMap<String, i64>) -> PlotlyChart {
    let chains: Vec<_> = unique_cdr3s_per_chain.keys().cloned().sorted().collect();
    let unique_cdr3s = chains.iter().map(|c| unique_cdr3s_per_chain[c]).collect();

    PlotlyChart::with_layout_and_data(
        Layout::new()
            .x_axis(Axis::new().title("Chain"))
            .y_axis(Axis::new().title("Unique CDR3s"))
            .hover_mode(HoverMode::X),
        vec![Bar::new(chains, unique_cdr3s)],
    )
}

fn vdj_barcode_rank_plot(barcode_info: &HashMap<String, BarcodeInfo>) -> PlotlyChart {
    barcode_rank_plot(barcode_info.values().map(|info| info.num_umis as u32))
}

fn barcode_rank_plot(umi_counts: impl Iterator<Item = u32>) -> PlotlyChart {
    let (x, y): (Vec<_>, Vec<_>) = umi_counts
        .sorted()
        .rev()
        .enumerate()
        .map(|(r, u)| (r + 1, u))
        .group_by(|(_, umis)| *umis)
        .into_iter()
        .flat_map(|(_, rank_umis)| {
            let rank_umis: Vec<_> = rank_umis.collect();
            let n = rank_umis.len();
            if n > 1 {
                vec![rank_umis[0], rank_umis[n - 1]].into_iter()
            } else {
                vec![rank_umis[0]].into_iter()
            }
        })
        .unzip();

    PlotlyChart::with_layout_and_data(
        Layout::new()
            .x_axis(Axis::new().type_(AxisType::Log))
            .y_axis(Axis::new().type_(AxisType::Log)),
        vec![Scatter::new(x, y)],
    )
}

fn gex_barcode_rank_plot(umis: &HashMap<String, u32>) -> PlotlyChart {
    barcode_rank_plot(umis.values().copied())
}

fn binned_hist(
    values: impl IntoIterator<Item = f64>,
    bin_width: f64,
) -> (Vec<f64>, Vec<usize>, Vec<String>) {
    let get_index = |value: f64| (value / bin_width).floor() as usize;
    let get_mid_x = |index: usize| bin_width * (index as f64 + 0.5);
    let get_umi_range_txt = |index: usize| {
        let start = bin_width * index as f64;
        let end = start + bin_width;
        let start = (10.0f64.powf(start) - 1.0).ceil() as usize;
        let end = (10.0f64.powf(end) - 1.0).floor() as usize;
        match start.cmp(&end) {
            Ordering::Less => format!("{start}-{end} UMIs"),
            Ordering::Equal => format!("{start} UMIs"),
            Ordering::Greater => String::new(),
        }
    };

    let hist = values.into_iter().map(get_index).counts();
    let max_index = *hist.keys().max().unwrap_or(&0);

    let x = (0..=max_index).map(get_mid_x).collect();
    let y = (0..=max_index)
        .map(|i| hist.get(&i).copied().unwrap_or(0))
        .collect();
    let hover_text = (0..=max_index).map(get_umi_range_txt).collect();

    (x, y, hover_text)
}

fn make_umi_histogram(
    categorized_barcodes: &BTreeMap<BarcodeCategory, Vec<&str>>,
    barcode_info: &HashMap<String, BarcodeInfo>,
    per_category: bool,
) -> PlotlyChart {
    // Compute bin width using <https://en.wikipedia.org/wiki/Freedman%E2%80%93Diaconis_rule
    // with all the barcodes
    let log_transform = |umis: usize| -> f64 { (umis as f64 + 1.0).log10() };
    let transformed_counts: Vec<_> = barcode_info
        .values()
        .map(|info| log_transform(info.num_umis))
        .collect();

    let bin_width = {
        let n = transformed_counts.len();
        let iqr = statrs::statistics::Data::new(transformed_counts.clone()).interquartile_range();
        (iqr / (n as f64).cbrt()).max(0.05)
    };

    let layout = Layout::new()
        .x_axis(Axis::new().title("Log10(1+ UMIs)"))
        .y_axis(Axis::new().title("Frequency"))
        .bar_mode(BarMode::Overlay)
        .hover_mode(HoverMode::X);

    let data = if per_category {
        categorized_barcodes
            .iter()
            .map(|(category, barcodes)| {
                // Using `Bar` is more efficient that `Histogram::new` although
                // the latter requires much less code
                // Histogram::new(
                //     barcodes
                //         .iter()
                //         .map(|bc| log_transform(barcode_info[*bc].num_umis)),
                // )
                // .x_bins(Bins::new(0.0, max_count, bin_width))

                let (x, y, text) = binned_hist(
                    barcodes
                        .iter()
                        .map(|bc| log_transform(barcode_info[*bc].num_umis)),
                    bin_width,
                );
                let bar = Bar::new(x, y)
                    .hover_template("%{text}: %{y}")
                    .name(category.label())
                    .marker(Marker::new().color(category.color(Some(0.75))))
                    .text_array(text)
                    .text_position(TextPosition::None);

                width_updated(bar, bin_width)
            })
            .collect()
    } else {
        let (x, y, text) = binned_hist(transformed_counts, bin_width);
        let bar = Bar::new(x, y)
            .hover_template("%{text}: %{y}")
            .text_array(text)
            .text_position(TextPosition::None);

        vec![width_updated(bar, bin_width)]
    };

    PlotlyChart::with_layout_and_data(layout, data)
}

fn width_updated(bar: Box<Bar<f64, usize>>, bin_width: f64) -> serde_json::Value {
    // HACK: plotly rust only supports usize bar widths which is a bug
    let mut bar_json = serde_json::to_value(&bar).unwrap();
    bar_json
        .as_object_mut()
        .unwrap()
        .insert("width".to_string(), bin_width.into());
    bar_json
}

fn make_contigs_per_barcode_plot(
    categorized_barcodes: &BTreeMap<BarcodeCategory, Vec<&str>>,
    barcode_info: &HashMap<String, BarcodeInfo>,
) -> PlotlyChart {
    let per_category_histogram: BTreeMap<_, _> = categorized_barcodes
        .iter()
        .map(|(category, barcodes)| {
            (
                *category,
                barcodes
                    .iter()
                    .map(|bc| barcode_info[*bc].num_productive_contigs)
                    .counts(),
            )
        })
        .collect();

    let max_contigs = *per_category_histogram
        .values()
        .flat_map(HashMap::keys)
        .max()
        .unwrap_or(&0);

    let data = per_category_histogram
        .into_iter()
        .map(|(category, histogram)| {
            let x: Vec<_> = (1..=max_contigs).collect();
            let y = x.iter().map(|i| *histogram.get(i).unwrap_or(&0)).collect();
            Bar::new(x, y)
                .name(category.label())
                .marker(Marker::new().color(category.color(None)))
        })
        .collect();

    PlotlyChart::with_layout_and_data(
        Layout::new()
            .x_axis(Axis::new().title("Number of productive contigs"))
            .y_axis(Axis::new().title("Frequency"))
            .bar_mode(BarMode::Stack)
            .hover_mode(HoverMode::X),
        data,
    )
}

fn make_cdr3_plots(
    categorized_barcodes: &BTreeMap<BarcodeCategory, Vec<&str>>,
    barcode_info: &HashMap<String, BarcodeInfo>,
) -> (PlotlyChart, PlotlyChart, PlotlyChart) {
    const N_CDR3: usize = 15;
    let top_cdr3s: Vec<_> = barcode_info
        .values()
        .flat_map(|i| i.cdr3s.iter().map(|c| &c.cdr3_aa))
        .counts()
        .into_iter()
        .sorted_by_key(|(_, c)| Reverse(*c))
        .take(N_CDR3)
        .map(|(k, _)| k.clone())
        .collect();

    let index_of_cdr3: HashMap<_, _> = top_cdr3s
        .iter()
        .enumerate()
        .map(|(i, cdr3)| (cdr3, i))
        .collect();

    let mut rng: StdRng = rand::SeedableRng::seed_from_u64(0);
    let jitter = rand_distr::Uniform::new(-0.2f32, 0.2f32);

    let umi_scatters = categorized_barcodes
        .iter()
        .map(|(category, barcodes)| {
            let (x, y): (Vec<_>, Vec<_>) = barcodes
                .iter()
                .flat_map(|&bc| &barcode_info[bc].cdr3s)
                .filter_map(
                    |Cdr3Info {
                         cdr3_aa, num_umis, ..
                     }| {
                        index_of_cdr3
                            .get(cdr3_aa)
                            .map(|&i| (i as f32 + jitter.sample(&mut rng), *num_umis))
                    },
                )
                .unzip();
            Scatter::new(x, y)
                .name(category.label())
                .marker(Marker::new().color(category.color(None)))
                .mode(plotly::common::Mode::Markers)
                .text_array(top_cdr3s.clone())
        })
        .collect();

    let umi_scatter_plot = PlotlyChart::with_layout_and_data(
        Layout::new()
            .x_axis(
                Axis::new()
                    .title("CDR3 Amino acid")
                    .tick_values((0..N_CDR3).map(|i| i as f64).collect())
                    .tick_text(top_cdr3s.clone()),
            )
            .y_axis(Axis::new().title("UMI Counts").type_(AxisType::Log)),
        umi_scatters,
    );

    let frequency_plot = cdr3_frequency_plot(categorized_barcodes, barcode_info, &top_cdr3s, false);
    let frequency_plot_umi =
        cdr3_frequency_plot(categorized_barcodes, barcode_info, &top_cdr3s, true);
    (frequency_plot, frequency_plot_umi, umi_scatter_plot)
}

fn cdr3_frequency_plot(
    categorized_barcodes: &BTreeMap<BarcodeCategory, Vec<&str>>,
    barcode_info: &HashMap<String, BarcodeInfo>,
    top_cdr3s: &[String],
    umi_weighted: bool,
) -> PlotlyChart {
    let frequency_bars = categorized_barcodes
        .iter()
        .map(|(category, barcodes)| {
            let cdr3_counts = if umi_weighted {
                let mut counts = HashMap::new();
                for Cdr3Info {
                    cdr3_aa, num_umis, ..
                } in barcodes.iter().flat_map(|&bc| &barcode_info[bc].cdr3s)
                {
                    *counts.entry(cdr3_aa).or_insert(0) += num_umis;
                }
                counts
            } else {
                barcodes
                    .iter()
                    .flat_map(|&bc| barcode_info[bc].cdr3s.iter().map(|c| &c.cdr3_aa))
                    .counts()
            };
            Bar::new(
                top_cdr3s.to_vec(),
                top_cdr3s
                    .iter()
                    .map(|cdr3| *cdr3_counts.get(cdr3).unwrap_or(&0))
                    .collect(),
            )
            .name(category.label())
            .marker(Marker::new().color(category.color(None)))
        })
        .collect();
    PlotlyChart::with_layout_and_data(
        Layout::new()
            .x_axis(Axis::new().title("CDR3 Amino acid"))
            .y_axis(Axis::new().title("Frequency"))
            .bar_mode(BarMode::Stack)
            .hover_mode(HoverMode::X),
        frequency_bars,
    )
}

fn make_gex_vdj_umi_scatter(
    gex_umis: &HashMap<String, u32>,
    categorized_barcodes: &BTreeMap<BarcodeCategory, Vec<&str>>,
    barcode_info: &HashMap<String, BarcodeInfo>,
    sub_filters: &HashMap<&str, String>,
) -> PlotlyChart {
    let data = categorized_barcodes
        .iter()
        .map(|(category, barcodes)| {
            let (vdj_umi_counts, gex_umi_counts, hover): (Vec<_>, Vec<_>, Vec<_>) =
                itertools::multiunzip(barcodes.iter().map(|bc| {
                    let info = &barcode_info[*bc];
                    let g = gex_umis.get(*bc).copied().unwrap_or(0);
                    let hover_suffix = match sub_filters.get(*bc) {
                        Some(filt) => filt.to_string(),
                        None => String::new(),
                    };
                    let hover_text = format!(
                        "{bc}<br>Contigs:{}<br>{hover_suffix}",
                        info.num_productive_contigs
                    );
                    (info.num_umis, g, hover_text)
                }));
            Scatter::new(vdj_umi_counts, gex_umi_counts)
                .name(category.label())
                .marker(Marker::new().color(category.color(None)))
                .mode(plotly::common::Mode::Markers)
                .text_array(hover)
                .hover_template("%{text}<br>VDJ: %{x}<br>GEX: %{y}")
        })
        .collect();

    PlotlyChart::with_layout_and_data(
        Layout::new()
            .x_axis(Axis::new().title("VDJ UMIs").type_(AxisType::Log))
            .y_axis(Axis::new().title("GEX UMIs").type_(AxisType::Log)),
        data,
    )
}

fn make_vdj_filters_per_barcode(
    asm_filters: &HashMap<String, Vec<&str>>,
    barcode_info: &HashMap<String, BarcodeInfo>,
) -> Vec<PerBarcodeFilter> {
    let mut bc_filters: Vec<PerBarcodeFilter> = Vec::new();
    for (bc, info) in barcode_info {
        let mut this_bc = PerBarcodeFilter {
            barcode: bc.clone(),
            is_cell: info.is_cell,
            is_gex_cell: info.is_gex_cell,
            is_asm_cell: info.is_asm_cell,
            low_umi: false,
            no_v_region: false,
            low_junction_support: false,
            no_conf_contig: false,
            low_rpu: false,
            non_dominant_junction: false,
            weak_junction: false,
            chimeric: false,
            common_clone: false,
            gel_bead_contamination: false,
            gel_bead_indel: None,
            enclone_fate: info.enclone_fate.map(|f| f.label().to_string()),
            insert_priming: false,
        };
        // Evaluate state of per barcode assembly filters
        if let Some(filters) = asm_filters.get(bc) {
            // filters scoped at the barcode level
            this_bc.low_umi = filters.contains(&"LOW_UMI");
            this_bc.no_v_region = filters.contains(&"NO_V_REGION");
            this_bc.low_junction_support = filters.contains(&"LOW_JUNCTION_SUPPORT");
            this_bc.no_conf_contig = filters.contains(&"NO_CONF_CONTIG");
            this_bc.low_rpu = filters.contains(&"LOW_RPU");

            // filters scoped at contig level and restricted to bcs that passed
            let was_excluded = this_bc.low_umi
                || this_bc.no_v_region
                || this_bc.low_junction_support
                || this_bc.no_conf_contig
                || this_bc.low_rpu;
            if !was_excluded {
                this_bc.gel_bead_indel = Some(filters.contains(&"GB_INDEL"));
            }

            // filters scoped across the library level
            this_bc.non_dominant_junction = filters.contains(&"NON_DOMINANT_JUNCTION");
            this_bc.weak_junction = filters.contains(&"WEAK_JUNCTION");
            this_bc.chimeric = filters.contains(&"CHIMERIC");
            this_bc.common_clone = filters.contains(&"COMMON_CLONE");
            this_bc.gel_bead_contamination = filters.contains(&"GB_CONTAMINATION");
            this_bc.insert_priming = filters.contains(&"INSERT_PRIMING");
        }

        bc_filters.push(this_bc);
    }
    bc_filters
}
martian_filetype! {HtmlFile, "html"}

#[derive(Clone, Deserialize, MartianStruct)]
pub struct SummarizeVdjFiltersStageInputs {
    sample_id: String,
    sample_description: Option<String>,
    all_contig_annotations: JsonFile<Vec<ContigAnnotation>>,
    asm_filter_diagnostics: Option<VdjFilterLogFormat>,
    enclone_barcode_fate: JsonFile<()>,
    raw_matrix_h5: Option<H5File>,
}

#[derive(Serialize, Deserialize, MartianStruct)]
pub struct SummarizeVdjFiltersStageOutputs {
    filter_summary: HtmlFile,
    metrics_summary: JsonFile<VdjFilterMetrics>,
    per_bc_filters: ParquetFile,
}

pub struct SummarizeVdjFilters;

#[make_mro(mem_gb = 5)]
impl MartianMain for SummarizeVdjFilters {
    type StageInputs = SummarizeVdjFiltersStageInputs;
    type StageOutputs = SummarizeVdjFiltersStageOutputs;

    fn main(&self, args: Self::StageInputs, rover: MartianRover) -> Result<Self::StageOutputs> {
        let FilterSummary {
            html,
            metrics,
            per_barcode_filters,
        } = generate_filter_summary(
            &args.sample_id,
            args.sample_description.as_deref(),
            Annotations::json(&args.all_contig_annotations, &args.enclone_barcode_fate),
            args.asm_filter_diagnostics.as_ref(),
            args.raw_matrix_h5.as_ref(),
        )?;

        let filter_summary_html: HtmlFile = rover.make_path("filter_summary");

        html.generate_html_file_with_build_files(
            &filter_summary_html,
            websummary_build::build_files()?,
        )?;

        let metrics_summary: JsonFile<_> = rover.make_path("metrics_summary");
        metrics_summary.write(&metrics)?;

        let per_bc_filter_pq: ParquetFile = rover.make_path("per_bc_filter");
        let mut pq_writer: ParquetWriter<PerBarcodeFilter> =
            per_bc_filter_pq.writer(PerBarcodeFilter::ROW_GROUP_SIZE)?;
        pq_writer.write_all(&per_barcode_filters)?;

        Ok(SummarizeVdjFiltersStageOutputs {
            filter_summary: filter_summary_html,
            metrics_summary,
            per_bc_filters: per_bc_filter_pq,
        })
    }
}
