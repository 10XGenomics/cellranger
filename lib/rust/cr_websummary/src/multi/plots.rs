use crate::{CardWithMetric, ChartWithHelp, PlotlyChart, RawChartWithHelp, TitleWithHelp};
use anyhow::Result;
use cr_types::TargetingMethod;
use csv::Reader;
use itertools::Itertools;
use metric::{TxHashMap, TxHashSet};
use plotly::common::{Anchor, DashType, Line, Marker, Mode, Orientation, Title};
use plotly::layout::{Axis, HoverMode, Legend, Margin, Shape, ShapeLine, ShapeType};
use plotly::traces::Scatter;
use plotly::Layout;
use regex::Regex;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::path::Path;

// CONSTANTS / LABELS
const SATURATION_LINE: f64 = 0.9;
const ON_TARGET_COLOR: &str = "#0071D9";
const ON_TARGET_LABEL: &str = "Targeted";
const OFF_TARGET_COLOR: &str = "#555555";
const OFF_TARGET_LABEL: &str = "Non-Targeted";

const TARGETED_ENRICHMENT_PLOT_TITLE: &str = "Reads vs UMIs per Gene by Targeting Status";
const TARGETED_ENRICHMENT_PLOT_X_LABEL: &str = "UMIs per gene, log10";
const TARGETED_ENRICHMENT_PLOT_Y_LABEL: &str = "Reads per gene, log10";
const TARGETED_ENRICHMENT_PLOT_HELP_TEXT: &str = "This plot shows the number of reads per gene (log10) in cell-associated barcodes as a function of the number of UMIs per gene (log10) in cell-associated barcodes, with genes colored by whether or not they were targeted. The yellow line represents the optimal threshold (specifically, its y-intercept = log10(Reads per UMI threshold)) in a mixture model that classifies genes into two classes based on Reads per UMI values. If sequencing saturation is low, this line will not be shown. Ideally, targeted genes will lie above the dotted line while non-targeted genes will be below it.";
const TARGETED_ENRICHMENT_PLOT_SEPARATION_BOUNDARY_COLOR: &str = "#E6B72B";

const CMO_JIBES_BIPLOT_TITLE: &str = "Biplots of CMO Count";
const CMO_JIBES_BIPLOT_HELP_TEXT: &str = "The plot shows the relationships between Cell Multiplexing Oligo (CMO) UMI counts for cells. Each point is a cell and the X and Y axes are UMI counts for a given CMO in the log10 scale. The CMOs on the axes can be changed with the selector. The cells which are not confidently assigned to CMO Tags are indicated. The number of cells has been downsampled to a maximum count of 100,000.";

const CMO_TAGS_ON_TSNE_PLOT_HELP_TEXT: &str = "Shown here are the CMO tag assignments for each cell-barcode. The axes correspond to the 2-dimensional embedding produced by the t-SNE algorithm over multiplexing features. In this space, pairs of cells that are close to each other have more similar Multiplexing Capture profiles than cells that are distant from each other. The display is limited to a random subset of cells.";
const CMO_TAGS_ON_TSNE_PLOT_TITLE: &str = "t-SNE Projection of Cells by CMO";

const BARCODE_RANK_PLOT_TITLE: &str = "Barcode Rank Plot";
const BARCODE_RANK_PLOT_HELP_TEXT: &str = "The plot shows filtered UMI counts mapped to each GEM barcode. Barcode-cell associations can be determined by UMI count or expression profile, or removed by Protein Aggregate Detection and Filtering and/or High Occupancy GEM Filtering steps. Therefore, some regions of the graph contain both cell-associated and background-associated barcodes. When present, Gene Expression data is used to identify these barcode populations. The color of the graph is based on the local density of barcodes that are cell-associated in these regions.";

const SEQUENCING_SATURATION_PLOT_X_LABEL: &str = "Mean Reads per Cell";
const SEQUENCING_SATURATION_PLOT_Y_LABEL: &str = "Sequencing Saturation";
const SEQUENCING_SATURATION_PLOT_ONTARGET_LABEL: &str = "Targeted";
const SEQUENCING_SATURATION_PLOT_OFFTARGET_LABEL: &str = "Non-Targeted";
const LIBRARY_SEQUENCING_SATURATION_PLOT_HELP_TEXT: &str = "This plot shows the Sequencing Saturation metric as a function of downsampled sequencing depth (measured in mean reads per cell), up to the observed sequencing depth. Sequencing Saturation is a measure of the observed library complexity, and approaches 1.0 (100%) when all converted mRNA transcripts (or ligation products in the context of Fixed RNA Profiling libraries) have been sequenced. The slope of the curve near the endpoint can be interpreted as an upper bound to the benefit to be gained from increasing the sequencing depth beyond this point. The dotted line is drawn at a value reasonably approximating the saturation point.";

const MEDIAN_GENES_PLOT_X_LABEL: &str = "Mean Reads per Cell";
const MEDIAN_GENES_PLOT_Y_LABEL: &str = "Median Genes per Cell";
const MEDIAN_GENES_PLOT_Y_LABEL_TARGETED: &str = "Median Targeted Genes per Cell";
const MEDIAN_GENES_PLOT_TARGETED_SUFFIX: &str = "_Targeted";
const LIBRARY_MEDIAN_GENES_PLOT_HELP_TEXT: &str = "This plot shows the Median Genes per Cell as a function of downsampled sequencing depth in mean reads per cell, up to the observed sequencing depth. The slope of the curve near the endpoint can be interpreted as an upper bound to the benefit to be gained from increasing the sequencing depth beyond this point.";
const SAMPLE_MEDIAN_GENES_PLOT_HELP_TEXT: &str = "This plot shows the Median Genes per Cell as a function of downsampled sequencing depth in mean reads per cell, up to the observed sequencing depth. The mean reads per cell is calculated using reads from cells in this sample only. The slope of the curve near the endpoint can be interpreted as an upper bound to the benefit to be gained from increasing the sequencing depth beyond this point.";

pub fn standard_layout(x_label: &str, y_label: &str) -> Layout {
    Layout::new()
        .margin(Margin::new().left(70).right(65).top(30).bottom(70))
        .show_legend(true)
        .hover_mode(HoverMode::Closest)
        .x_axis(Axis::new().title(Title::new(x_label)))
        .y_axis(Axis::new().title(Title::new(y_label)))
        .legend(
            Legend::new()
                .y_anchor(Anchor::Bottom)
                .y(0.1)
                .x_anchor(Anchor::Right)
                .x(0.99)
                .background_color("#ffffff"),
        )
}

pub fn sequencing_saturation_layout(x_label: &str, y_label: &str, x_max: f64) -> Layout {
    standard_layout(x_label, y_label)
        .y_axis(Axis::new().title(Title::new(y_label)).range(vec![0.0, 1.0]))
        .shapes(vec![Shape::new()
            .shape_type(ShapeType::Line)
            .x0(0)
            .y0(SATURATION_LINE)
            .x1(x_max)
            .y1(SATURATION_LINE)
            .line(
                ShapeLine::new()
                    .color("#999999")
                    .width(4.0)
                    .dash(DashType::Dot),
            )])
}

pub fn format_jibes_biplots(plot: &Value) -> RawChartWithHelp {
    RawChartWithHelp {
        plot: plot.clone(),
        help: TitleWithHelp {
            help: CMO_JIBES_BIPLOT_HELP_TEXT.to_string(),
            title: CMO_JIBES_BIPLOT_TITLE.to_string(),
        },
    }
}

pub fn format_umi_on_tsne_plot(plot: &Value, library_type: &str, title: &str) -> RawChartWithHelp {
    RawChartWithHelp {
        plot: plot.clone(),
        help: TitleWithHelp {
            help: format!(
                "Shown here are the total {library_type} UMI counts for each cell-barcode. The axes correspond to the 2-dimensional embedding produced by the t-SNE algorithm over the {library_type} features. In this space, pairs of cells that are close to each other have more similar {library_type} profiles than cells that are distant from each other. The display is limited to a random subset of cells."
            ),
            title: title.to_string(),
        },
    }
}

pub fn format_tags_on_tsne_plot(plot: &Value) -> RawChartWithHelp {
    RawChartWithHelp {
        plot: plot.clone(),
        help: TitleWithHelp {
            help: CMO_TAGS_ON_TSNE_PLOT_HELP_TEXT.to_string(),
            title: CMO_TAGS_ON_TSNE_PLOT_TITLE.to_string(),
        },
    }
}

pub fn format_histogram(plot: &Value, feature: &str) -> RawChartWithHelp {
    RawChartWithHelp {
        plot: plot.clone(),
        help: TitleWithHelp {
            help: format!(
                "Histogram of {feature} counts per cell, for each {feature}.  The X-axis is the UMI counts in the log scale, while the Y-axis is the number of cells."
            ),
            title: format!("Histogram of {feature} Count"),
        },
    }
}

pub fn format_barcode_rank_plot(plot: &PlotlyChart, library: &str) -> ChartWithHelp {
    ChartWithHelp {
        plot: plot.clone(),
        help: TitleWithHelp {
            help: BARCODE_RANK_PLOT_HELP_TEXT.to_string(),
            title: format!("{library} {BARCODE_RANK_PLOT_TITLE}"),
        },
    }
}

/// for a plot drawn from metrics where the y values suddenly drop to 0
/// trim off the part that drops to 0
/// requires that vector of (x,y) used as input is sorted by x
fn trim_plot(data: Vec<(f64, f64)>) -> (Vec<f64>, Vec<f64>) {
    let mut y_prev = 0.0;
    let y_trimmed: Vec<_> = data
        .iter()
        .map(|(_, y)| *y)
        .take_while(|y| {
            if *y <= 0.0 && y_prev > 0.0 {
                return false;
            }
            y_prev = *y;
            true
        })
        .collect();
    let x_trimmed = data
        .into_iter()
        .take(y_trimmed.len())
        .map(|(x, _)| x)
        .collect();
    (x_trimmed, y_trimmed)
}

pub enum PlotType {
    LibraryPlot,
    SamplePlot,
}

// TODO: Get rid of this and create a simple deserializeable data structure representing this data inside subsample.py and pass that forward
pub fn library_sequencing_saturation_plot_from_metrics(
    metrics: &TxHashMap<String, Value>,
    is_targeted: bool,
) -> ChartWithHelp {
    let (re_untargeted, re_ontarget, re_offtarget) = (
        r"^multi_raw_rpc_(\d+)_subsampled_duplication_frac$",
        r"^multi_raw_rpc_(\d+)_subsampled_duplication_frac_ontarget$",
        r"^multi_raw_rpc_(\d+)_subsampled_duplication_frac_offtarget$",
    );
    let get_xy_data = |re: Regex| -> (Vec<f64>, Vec<f64>) {
        let mut xy_data: Vec<(f64, f64)> = vec![(0.0, 0.0)];
        xy_data.extend(
            metrics
                .iter()
                .filter_map(move |(metric_name, metric_value)| {
                    re.captures(metric_name).map(move |cap| {
                        (
                            cap.get(1).unwrap().as_str().parse::<usize>().unwrap() as f64,
                            metric_value.as_f64().unwrap(),
                        )
                    })
                }),
        );

        xy_data.sort_by(|a, b| {
            a.partial_cmp(b)
                .expect("Error sorting plot data, possibly due to NaN value")
        });
        trim_plot(xy_data)
    };

    let help_text = LIBRARY_SEQUENCING_SATURATION_PLOT_HELP_TEXT;

    if is_targeted {
        let re_ontarget = Regex::new(re_ontarget).unwrap();
        let (x_data_ontarget, y_data_ontarget) = get_xy_data(re_ontarget);
        let re_offtarget = Regex::new(re_offtarget).unwrap();
        let (x_data_offtarget, y_data_offtarget) = get_xy_data(re_offtarget);
        let layout = sequencing_saturation_layout(
            SEQUENCING_SATURATION_PLOT_X_LABEL,
            SEQUENCING_SATURATION_PLOT_Y_LABEL,
            *x_data_ontarget.last().unwrap_or(&0.),
        );
        let data = vec![
            Scatter::new(x_data_ontarget, y_data_ontarget)
                .name(SEQUENCING_SATURATION_PLOT_ONTARGET_LABEL)
                .mode(Mode::Lines)
                .fill_color(ON_TARGET_COLOR)
                .marker(Marker::new().color(ON_TARGET_COLOR))
                .line(Line::new().width(3.0)),
            Scatter::new(x_data_offtarget, y_data_offtarget)
                .name(SEQUENCING_SATURATION_PLOT_OFFTARGET_LABEL)
                .mode(Mode::Lines)
                .fill_color(OFF_TARGET_COLOR)
                .marker(Marker::new().color(OFF_TARGET_COLOR))
                .line(Line::new().width(3.0)),
        ];
        ChartWithHelp {
            plot: PlotlyChart::with_layout_and_data(layout.show_legend(true), data),
            help: TitleWithHelp {
                help: help_text.to_string(),
                title: SEQUENCING_SATURATION_PLOT_Y_LABEL.to_string(),
            },
        }
    } else {
        let re_untargeted = Regex::new(re_untargeted).unwrap();
        let (x_data_untargeted, y_data_untargeted) = get_xy_data(re_untargeted);
        let layout = sequencing_saturation_layout(
            SEQUENCING_SATURATION_PLOT_X_LABEL,
            SEQUENCING_SATURATION_PLOT_Y_LABEL,
            *x_data_untargeted.last().unwrap_or(&0.),
        );
        let data = vec![Scatter::new(x_data_untargeted, y_data_untargeted)
            .mode(Mode::Lines)
            .line(Line::new().width(3.0))];
        ChartWithHelp {
            plot: PlotlyChart::with_layout_and_data(layout.show_legend(false), data),
            help: TitleWithHelp {
                help: help_text.to_string(),
                title: SEQUENCING_SATURATION_PLOT_Y_LABEL.to_string(),
            },
        }
    }
}

// TODO: Get rid of this and create a simple deserializeable data structure representing this data inside subsample.py and pass that forward
pub fn sample_median_genes_plot_from_metrics(
    metrics: &TxHashMap<String, Value>,
    genomes: TxHashSet<String>,
    plot_type: PlotType,
    targeting_method: Option<TargetingMethod>,
) -> ChartWithHelp {
    median_genes_plot_from_metrics(metrics, genomes, plot_type, targeting_method, true)
}

pub fn library_median_genes_plot_from_metrics(
    metrics: &TxHashMap<String, Value>,
    genomes: TxHashSet<String>,
    plot_type: PlotType,
    targeting_method: Option<TargetingMethod>,
) -> ChartWithHelp {
    median_genes_plot_from_metrics(metrics, genomes, plot_type, targeting_method, false)
}

fn median_genes_plot_from_metrics(
    metrics: &TxHashMap<String, Value>,
    genomes: TxHashSet<String>,
    plot_type: PlotType,
    targeting_method: Option<TargetingMethod>,
    is_sample: bool,
) -> ChartWithHelp {
    let (re_untargeted, re_ontarget) = if is_sample {
        (
            r"^(.+)_raw_barcoded_filtered_bc_rpc_(\d+)_subsampled_filtered_bcs_median_unique_genes_detected$",
            r"^(.+)_raw_barcoded_filtered_bc_rpc_(\d+)_subsampled_filtered_bcs_median_unique_genes_detected_ontarget$",
        )
    } else {
        (
            r"^(.+)_raw_rpc_(\d+)_subsampled_filtered_bcs_median_unique_genes_detected$",
            r"^(.+)_raw_rpc_(\d+)_subsampled_filtered_bcs_median_unique_genes_detected_ontarget$",
        )
    };
    let get_xy_data = move |re: Regex| -> TxHashMap<&str, Vec<(f64, f64)>> {
        let mut xy_data: TxHashMap<&str, Vec<(f64, f64)>> = TxHashMap::default();

        for (metric_name, metric_value) in metrics {
            if let Some(cap) = re.captures(metric_name) {
                let genome = cap.get(1).unwrap().as_str();
                // don't want to plot "multi" genome
                // also skips "ANTIBODY_", "MULTIPLEXING",
                if !genomes.contains(genome) {
                    continue;
                }

                xy_data.entry(genome).or_default().push((
                    cap.get(2).unwrap().as_str().parse::<usize>().unwrap() as f64,
                    metric_value.as_f64().unwrap(),
                ));
            }
        }
        xy_data
    };

    fn make_scatter_data(
        xy_data: TxHashMap<&str, Vec<(f64, f64)>>,
        label_suffix: &str,
    ) -> Vec<Scatter<f64, f64>> {
        xy_data
            .into_iter()
            .map(|(genome, mut xy_data_genome)| {
                xy_data_genome.push((0.0, 0.0));
                xy_data_genome.sort_by(|a, b| {
                    a.partial_cmp(b)
                        .expect("Error sorting plot data, possibly due to NaN value")
                });
                let (x_data_trimmed, y_data_trimmed) = trim_plot(xy_data_genome);

                let mut label = genome.to_string();
                label.push_str(label_suffix);
                *Scatter::new(x_data_trimmed, y_data_trimmed)
                    .name(&label)
                    .mode(Mode::Lines)
                    .line(Line::new().width(3.0))
            })
            .collect()
    }

    let (data, y_label) = if let Some(targeting_method) = targeting_method {
        let xy_data_ontarget = get_xy_data(Regex::new(re_ontarget).unwrap());
        if targeting_method == TargetingMethod::TemplatedLigation {
            (
                make_scatter_data(xy_data_ontarget, ""),
                MEDIAN_GENES_PLOT_Y_LABEL,
            )
        } else {
            (
                make_scatter_data(xy_data_ontarget, MEDIAN_GENES_PLOT_TARGETED_SUFFIX),
                MEDIAN_GENES_PLOT_Y_LABEL_TARGETED,
            )
        }
    } else {
        let xy_data_untargeted = get_xy_data(Regex::new(re_untargeted).unwrap());
        (
            make_scatter_data(xy_data_untargeted, ""),
            MEDIAN_GENES_PLOT_Y_LABEL,
        )
    };

    let help_text = match plot_type {
        PlotType::LibraryPlot => LIBRARY_MEDIAN_GENES_PLOT_HELP_TEXT,
        PlotType::SamplePlot => SAMPLE_MEDIAN_GENES_PLOT_HELP_TEXT,
    };

    let layout = standard_layout(MEDIAN_GENES_PLOT_X_LABEL, y_label);

    ChartWithHelp {
        plot: PlotlyChart::with_layout_and_data(layout, data),
        help: TitleWithHelp {
            help: help_text.to_string(),
            title: y_label.to_string(),
        },
    }
}

#[derive(Copy, Debug, Deserialize, PartialEq, Clone)]
enum PythonBool {
    False,
    True,
}

impl From<PythonBool> for bool {
    fn from(item: PythonBool) -> Self {
        match item {
            PythonBool::False => false,
            PythonBool::True => true,
        }
    }
}

#[derive(Debug, Deserialize, PartialEq)]
struct TargetedPerFeatureMetricsRow {
    feature_type: String,
    feature_id: String,
    feature_name: String,
    feature_idx: usize,
    num_umis: usize,
    num_reads: usize,
    num_umis_cells: usize,
    num_reads_cells: usize,
    dup_frac: f64,
    dup_frac_cells: f64,
    is_targeted: PythonBool,
    mean_reads_per_umi: Option<f64>,
    mean_reads_per_umi_log10: Option<f64>,
    mean_reads_per_umi_cells: Option<f64>,
    mean_reads_per_umi_cells_log10: Option<f64>,
    enriched_rpu: Option<PythonBool>,
}

#[derive(Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct TargetedEnrichmentTableAndPlot {
    pub plot: PlotlyChart,
    pub gene_targeting_table: CardWithMetric,
}

// Rust implementation of the reads_per_umi_plot in analysis_tab.py
// does NOT support spatial whatsoever
pub fn targeted_enrichment_plot(
    per_feature_metrics_csv: Option<&Path>,
    log_rpu_threshold: Option<f64>, // from the metrics json
) -> Result<Option<ChartWithHelp>> {
    let mut rdr = match per_feature_metrics_csv {
        Some(csv) => Reader::from_path(csv)?,
        None => return Ok(None),
    };
    // get the rows records of per-feature-metrics CSV
    let per_feature_metrics: Vec<TargetedPerFeatureMetricsRow> = rdr.deserialize().try_collect()?;

    let (targeted_per_feature_metrics_subset, untargeted_per_feature_metrics_subset) =
        per_feature_metrics
            .iter()
            .filter(|&x| x.num_reads_cells > 0)
            .partition::<Vec<&TargetedPerFeatureMetricsRow>, _>(|x| bool::from(x.is_targeted));
    let groups = [
        (
            OFF_TARGET_LABEL.to_string(),
            untargeted_per_feature_metrics_subset,
            OFF_TARGET_COLOR.to_string(),
        ),
        (
            ON_TARGET_LABEL.to_string(),
            targeted_per_feature_metrics_subset,
            ON_TARGET_COLOR.to_string(),
        ),
    ];
    let mut data = Vec::with_capacity(groups.len() + 1);
    for (name, per_feature_metrics_subset, color) in groups {
        // remove totally overlapping points
        let per_feature_metrics_subset_deduped: Vec<&TargetedPerFeatureMetricsRow> =
            per_feature_metrics_subset
                .into_iter()
                .unique_by(|x| (x.num_umis_cells, x.num_reads_cells))
                .collect();

        let gene_name_labels: Vec<String> = per_feature_metrics_subset_deduped
            .iter()
            .map(|x| format!("Gene: {}", x.feature_name))
            .collect();
        let x: Vec<f64> = per_feature_metrics_subset_deduped
            .iter()
            .map(|x| (x.num_umis_cells as f64).log10())
            .collect();
        let y: Vec<f64> = per_feature_metrics_subset_deduped
            .iter()
            .map(|x| (x.num_reads_cells as f64).log10())
            .collect();

        data.push(
            Scatter::new(x, y)
                .name(name)
                .mode(Mode::Markers)
                .fill_color(ON_TARGET_COLOR)
                .marker(Marker::new().color(color).size(5).opacity(0.5))
                .hover_text_array(gene_name_labels),
        );
    }
    // plot separation boundary on top, if it was determined
    if let Some(b) = log_rpu_threshold {
        let a = 1.0;
        let x0 = 0.0;
        let x1 = (per_feature_metrics
            .into_iter()
            .map(|x| x.num_umis_cells)
            .max()
            .unwrap() as f64
            + 1.0)
            .log10();

        data.push(
            Scatter::new(vec![x0, x1], vec![x0 * a + b, x1 * a + b])
                .name(format!("Reads per UMI threshold {:.2}", 10.0_f64.powf(b)))
                .mode(Mode::Lines)
                .fill_color(OFF_TARGET_COLOR)
                .marker(Marker::new().color(TARGETED_ENRICHMENT_PLOT_SEPARATION_BOUNDARY_COLOR))
                .line(Line::new().width(3.0)),
        )
    }

    let layout = Layout::new()
        .x_axis(Axis::new().title(Title::new(TARGETED_ENRICHMENT_PLOT_X_LABEL)))
        .y_axis(Axis::new().title(Title::new(TARGETED_ENRICHMENT_PLOT_Y_LABEL)))
        .legend(Legend::new().orientation(Orientation::Horizontal).y(-0.2))
        .hover_mode(HoverMode::Closest);

    Ok(Some(ChartWithHelp {
        plot: PlotlyChart::with_layout_and_data(layout, data),
        help: TitleWithHelp {
            help: TARGETED_ENRICHMENT_PLOT_HELP_TEXT.to_string(),
            title: TARGETED_ENRICHMENT_PLOT_TITLE.to_string(),
        },
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_trim_plot_two_trailing() {
        let expected_result = (vec![0.0, 100.0, 200.0], vec![0.0, 100.0, 200.0]);
        let xy_data = vec![
            (0.0, 0.0),
            (100.0, 100.0),
            (200.0, 200.0),
            (500.0, 0.0),
            (600.0, 0.0),
        ];
        assert_eq!(trim_plot(xy_data), expected_result);
    }

    #[test]
    fn test_trim_plot_data_one_trailing() {
        let expected_result = (vec![0.0, 100.0, 200.0], vec![0.0, 100.0, 200.0]);
        let xy_data = vec![(0.0, 0.0), (100.0, 100.0), (200.0, 200.0), (500.0, 0.0)];
        assert_eq!(trim_plot(xy_data), expected_result);
    }

    #[test]
    fn test_trim_plot_data_no_trailing() {
        let expected_result = (vec![0.0, 100.0, 200.0], vec![0.0, 100.0, 200.0]);
        let xy_data = vec![(0.0, 0.0), (100.0, 100.0), (200.0, 200.0)];
        assert_eq!(trim_plot(xy_data), expected_result);
    }

    #[test]
    fn test_trim_plot_unsorted() {
        // trim_plot should require pre-sorted data to yield desired result
        let expected_result = (vec![0.0, 100.0, 200.0], vec![0.0, 100.0, 200.0]);
        let xy_data = vec![(100.0, 100.0), (0.0, 0.0), (200.0, 200.0)];
        assert_ne!(trim_plot(xy_data), expected_result);
    }
}
