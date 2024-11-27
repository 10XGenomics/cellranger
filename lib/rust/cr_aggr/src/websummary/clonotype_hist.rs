use crate::websummary::{ChartWithHelp, PlotlyChart, TitleWithHelp};
use plotly::layout::{Axis, AxisType, BarMode, HoverMode, Margin};
use plotly::{Bar, Layout};
use serde::Serialize;

#[derive(Debug, Serialize, PartialEq, Default, Clone)]
#[serde(into = "ChartWithHelp")]
pub struct ClonotypeHist {
    pub ids: Vec<usize>,
    pub proportions: Vec<(String, Vec<f64>)>, // One per origin
}

const TITLE: &str = "Top 10 Clonotype Frequencies";
const HELP: &str = "This histogram displays the fraction of cells (percentage of cells) occupied by the 10 most abundant clonotypes in the aggregated datasets. The histogram is colored on the basis of origin. The clonotype IDs on the X axis correspond to the clonotype IDs listed in the Top 10 Clonotype CDR3 Sequences below.";
const XLABEL: &str = "Clonotype ID";
const YLABEL: &str = "Fraction of Cells";

impl From<ClonotypeHist> for ChartWithHelp {
    fn from(hist: ClonotypeHist) -> ChartWithHelp {
        let layout = Layout::new()
            .show_legend(true)
            .bar_mode(BarMode::Stack)
            .hover_mode(HoverMode::Closest)
            .x_axis(Axis::new().type_(AxisType::Category).title(XLABEL))
            .y_axis(Axis::new().title(YLABEL).tick_format(",.2%"))
            .margin(Margin::new().left(60).right(40).top(25));

        let ids = hist.ids;
        let data = hist
            .proportions
            .into_iter()
            .map(|(name, y)| {
                Bar::new(ids.clone(), y)
                    .name(name)
                    .hover_template("%{y:.2%}")
            })
            .collect();

        ChartWithHelp {
            plot: PlotlyChart::with_layout_and_data(layout, data),
            help: TitleWithHelp {
                title: TITLE.into(),
                help: HELP.into(),
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::websummary::check_eq_json;

    #[test]
    fn test_clonotype_hist() {
        let hist = ClonotypeHist {
            ids: vec![1, 21],
            proportions: vec![
                ("PreVac".to_string(), vec![0.47, 0.12]),
                ("PostVac".to_string(), vec![0.0, 0.2]),
            ],
        };
        check_eq_json(
            &serde_json::to_string(&hist).unwrap(),
            r#"{
                "plot": {
                    "config": {
                        "displayModeBar": true,
                        "staticPlot": false,
                        "dragmode": "zoom",
                        "modeBarButtons": [
                            [
                                "toImage"
                            ]
                        ]
                    },
                    "data": [
                        {
                            "y": [
                                0.47,
                                0.12
                            ],
                            "x": [
                                1,
                                21
                            ],
                            "type": "bar",
                            "name": "PreVac",
                            "hovertemplate": "%{y:.2%}"
                        },
                        {
                            "y": [
                                0.0,
                                0.2
                            ],
                            "x": [
                                1,
                                21
                            ],
                            "type": "bar",
                            "name": "PostVac",
                            "hovertemplate": "%{y:.2%}"
                        }
                    ],
                    "layout": {
                        "showlegend": true,
                        "yaxis": {
                        	"tickformat": ",.2%",
                            "title": {
                                "text": "Fraction of Cells"
                            }
                        },
                        "barmode": "stack",
                        "xaxis": {
                            "type": "category",
                            "title": {
                                "text": "Clonotype ID"
                            }
                        },
                        "hovermode": "closest",
                        "margin": {
                            "r": 40,
                            "l": 60,
                            "t": 25
                        }
                    }
                },
                "help": {
                    "helpText": "This histogram displays the fraction of cells (percentage of cells) occupied by the 10 most abundant clonotypes in the aggregated datasets. The histogram is colored on the basis of origin. The clonotype IDs on the X axis correspond to the clonotype IDs listed in the Top 10 Clonotype CDR3 Sequences below.",
                    "title": "Top 10 Clonotype Frequencies"
                }
            }"#,
        );
    }
}
