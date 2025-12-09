#![expect(missing_docs)]
use crate::{ChartWithHelp, PlotlyChart, RawChartWithHelp, TitleWithHelp};
use cr_types::{GenomeName, LibraryType, TargetingMethod};
use metric::{AsMetricPrefix, TxHashMap, TxHashSet};
use plotly::Layout;
use plotly::common::{Anchor, Line, Mode};
use plotly::layout::{Axis, HoverMode, Legend, Margin};
use plotly::traces::Scatter;
use regex::Regex;
use serde::Deserialize;
use serde_json::Value;
use serde_json::value::RawValue;

// CONSTANTS / LABELS

const BARCODE_RANK_PLOT_TITLE: &str = "Barcode Rank Plot";
const BARCODE_RANK_PLOT_HELP_TEXT: &str = "The plot shows filtered UMI counts mapped to each GEM barcode. Barcode-cell associations can be determined by UMI count or expression profile, or removed by Protein Aggregate Detection and Filtering and/or High Occupancy GEM Filtering steps. Therefore, some regions of the graph contain both cell-associated and background-associated barcodes. When present, Gene Expression data is used to identify these barcode populations. The color of the graph is based on the local density of barcodes that are cell-associated in these regions. Hovering over the plot displays the total number and percentage of barcodes in that region called as cells along with the number of UMI counts for those barcodes and barcode rank, ordered in descending order of UMI counts.";

const BARNYARD_BIPLOT_TITLE: &str = "Barnyard Count Bi-plot";
const BARNYARD_BIPLOT_HELP_TEXT: &str = "The plot shows the UMI counts of each species for each cell. Cells are classified as singlets belonging to one species or as multiplets.";

const SEQUENCING_SATURATION_PLOT_X_LABEL: &str = "Mean Reads per Cell";
const SEQUENCING_SATURATION_PLOT_Y_LABEL: &str = "Sequencing Saturation";
const LIBRARY_SEQUENCING_SATURATION_PLOT_HELP_TEXT: &str = "This plot shows the Sequencing Saturation metric as a function of downsampled sequencing depth (measured in mean reads per cell), up to the observed sequencing depth. Sequencing Saturation is a measure of the observed library complexity, and approaches 1.0 (100%) when all converted mRNA transcripts (or ligation products in the context of Flex libraries) have been sequenced. The slope of the curve near the endpoint can be interpreted as an upper bound to the benefit to be gained from increasing the sequencing depth beyond this point.";

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
        .x_axis(Axis::new().title(x_label))
        .y_axis(Axis::new().title(y_label))
        .legend(
            Legend::new()
                .y_anchor(Anchor::Bottom)
                .y(0.1)
                .x_anchor(Anchor::Right)
                .x(0.99)
                .background_color("#ffffff"),
        )
}

pub fn sequencing_saturation_layout(x_label: &str, y_label: &str) -> Layout {
    standard_layout(x_label, y_label).y_axis(Axis::new().title(y_label).range(vec![0.0, 1.0]))
}

pub fn format_jibes_biplots(plot: Box<RawValue>, feature: &str) -> RawChartWithHelp {
    RawChartWithHelp {
        plot,
        help: TitleWithHelp {
            title: format!("Biplots of {feature} UMI Counts"),
            help: format!(
                "The plot shows the relationships between {feature} UMI counts for cells. Each point is a cell and the X and Y axes are UMI counts for a given {feature} in the log10 scale. The {feature}s on the axes can be changed with the selector. The cells which are not confidently assigned to {feature}s are indicated. The number of cells has been downsampled to a maximum count of 100,000."
            ),
        },
    }
}

pub fn format_umi_on_umap_plot(
    plot: Box<RawValue>,
    library_type: &str,
    title: &str,
) -> RawChartWithHelp {
    RawChartWithHelp {
        plot,
        help: TitleWithHelp {
            help: format!(
                "Shown here are the total {library_type} UMI counts for each cell-barcode. The axes correspond to the 2-dimensional embedding produced by the UMAP algorithm over the {library_type} features. In this space, pairs of cells that are close to each other have more similar {library_type} profiles than cells that are distant from each other. The display is limited to a random subset of cells."
            ),
            title: title.to_string(),
        },
    }
}

pub fn format_tags_on_umap_plot(plot: Box<RawValue>, library_type: &str) -> RawChartWithHelp {
    RawChartWithHelp {
        plot,
        help: TitleWithHelp {
            help: format!(
                "Shown here are the tag assignments for each cell-barcode. The axes correspond to the 2-dimensional embedding produced by the UMAP algorithm over {library_type} features. In this space, pairs of cells that are close to each other have more similar {library_type} profiles than cells that are distant from each other. The display is limited to a random subset of cells."
            ),
            title: "UMAP Projection of Cells Colored by Tag Assignment".to_string(),
        },
    }
}

pub fn format_histogram(plot: Box<RawValue>, feature: &str) -> RawChartWithHelp {
    RawChartWithHelp {
        plot,
        help: TitleWithHelp {
            help: format!(
                "Histogram of {feature} UMI counts per cell, for each {feature}.  The X-axis is the UMI counts in the log scale, while the Y-axis is the number of cells."
            ),
            title: format!("Histogram of {feature} UMI Counts"),
        },
    }
}

pub fn format_barcode_rank_plot(plot: Box<RawValue>, library: &str) -> RawChartWithHelp {
    RawChartWithHelp {
        plot,
        help: TitleWithHelp {
            help: BARCODE_RANK_PLOT_HELP_TEXT.to_string(),
            title: format!("{library} {BARCODE_RANK_PLOT_TITLE}"),
        },
    }
}

pub fn format_barnyard_biplot(plot: Box<RawValue>, library: &str) -> RawChartWithHelp {
    RawChartWithHelp {
        plot,
        help: TitleWithHelp {
            help: BARNYARD_BIPLOT_HELP_TEXT.to_string(),
            title: format!("{library} {BARNYARD_BIPLOT_TITLE}"),
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
// Note: these metrics originate in /lib/python/cellranger/subsample.py:calculate_subsampling_metrics
pub fn library_sequencing_saturation_plot_from_metrics(
    library_type: LibraryType,
    metrics: &TxHashMap<String, Value>,
) -> ChartWithHelp {
    let re = {
        let (prefix, sep) = library_type
            .as_metric_prefix()
            .map_or(("", ""), |prefix| (prefix, "_"));
        format!(r"^{prefix}{sep}multi_raw_rpc_(\d+)_subsampled_duplication_frac$")
    };
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

    let re = Regex::new(&re).unwrap();
    let (x_data, y_data) = get_xy_data(re);
    let layout = sequencing_saturation_layout(
        SEQUENCING_SATURATION_PLOT_X_LABEL,
        SEQUENCING_SATURATION_PLOT_Y_LABEL,
    );
    let data = vec![
        Scatter::new(x_data, y_data)
            .mode(Mode::Lines)
            .line(Line::new().width(3.0)),
    ];
    ChartWithHelp {
        plot: PlotlyChart::with_layout_and_data(layout.show_legend(false), data),
        help: TitleWithHelp {
            help: help_text.to_string(),
            title: SEQUENCING_SATURATION_PLOT_Y_LABEL.to_string(),
        },
    }
}

// TODO: Get rid of this and create a simple deserializeable data structure representing this data inside subsample.py and pass that forward
pub fn sample_median_genes_plot_from_metrics(
    metrics: &TxHashMap<String, Value>,
    genomes: TxHashSet<GenomeName>,
    plot_type: PlotType,
    targeting_method: Option<TargetingMethod>,
) -> ChartWithHelp {
    median_genes_plot_from_metrics(metrics, genomes, plot_type, targeting_method, true)
}

pub fn library_median_genes_plot_from_metrics(
    metrics: &TxHashMap<String, Value>,
    genomes: TxHashSet<GenomeName>,
    plot_type: PlotType,
    targeting_method: Option<TargetingMethod>,
) -> ChartWithHelp {
    median_genes_plot_from_metrics(metrics, genomes, plot_type, targeting_method, false)
}

fn median_genes_plot_from_metrics(
    metrics: &TxHashMap<String, Value>,
    genomes: TxHashSet<GenomeName>,
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
    let get_xy_data = move |re: Regex| -> TxHashMap<GenomeName, Vec<(f64, f64)>> {
        let mut xy_data: TxHashMap<GenomeName, Vec<(f64, f64)>> = TxHashMap::default();

        for (metric_name, metric_value) in metrics {
            if let Some(cap) = re.captures(metric_name) {
                let genome = GenomeName::from(cap.get(1).unwrap().as_str());
                // don't want to plot "multi" genome
                // also skips "ANTIBODY_", "MULTIPLEXING",
                if !genomes.contains(&genome) {
                    continue;
                }

                xy_data.entry(genome).or_default().push((
                    cap.get(2).unwrap().as_str().parse::<usize>().unwrap() as f64,
                    metric_value.as_f64().unwrap_or(0_f64),
                ));
            }
        }
        xy_data
    };

    fn make_scatter_data(
        xy_data: TxHashMap<GenomeName, Vec<(f64, f64)>>,
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
