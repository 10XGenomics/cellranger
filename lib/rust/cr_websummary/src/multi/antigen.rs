#![expect(missing_docs)]
use crate::{ChartWithHelp, PlotlyChart, TitleWithHelp};
use anyhow::Result;
use cr_types::clonotype::ClonotypeId;
use hclust::{ClusterDirection, DistanceMetric, HierarchicalCluster, LeafOrdering, LinkageMethod};
use itertools::Itertools;
use martian_filetypes::FileTypeRead;
use martian_filetypes::tabular_file::CsvFile;
use ndarray::Array2;
use ordered_float::OrderedFloat;
use plotly::common::{ColorBar, ColorScale, ColorScaleElement};
use plotly::layout::{Axis, AxisType};
use plotly::{HeatMap, Layout};
use serde::{Deserialize, Serialize};
use serde_json::json;
use statrs::statistics::{Data, Median};
use std::cmp::Reverse;
use std::collections::{BTreeSet, HashMap};

const TOP_N: usize = 10;
const MIN_SPECIFICITY: f64 = 25.0;

#[derive(Serialize, Deserialize)]
pub struct AntigenSpecificityRow {
    barcode: String,
    antigen: String,
    antigen_umi: u32,
    control: String,
    control_umi: u32,
    antigen_specificity_score: f64,
    mhc_allele: Option<String>,
    raw_clonotype_id: Option<String>,
    exact_subclonotype_id: Option<String>,
}

impl AntigenSpecificityRow {
    fn clonotype_id(&self) -> Option<String> {
        self.raw_clonotype_id
            .as_ref()
            .and_then(|id| if id == "None" { None } else { Some(id.clone()) })
    }
    fn sort_key(&self) -> (Option<&str>, &str, OrderedFloat<f64>) {
        (
            self.raw_clonotype_id.as_deref(),
            &self.antigen,
            self.antigen_specificity_score.into(),
        )
    }
}

fn make_index(iter: impl Iterator<Item = String>) -> (Vec<String>, HashMap<String, usize>) {
    let values: Vec<_> = iter.unique().sorted().collect();
    let index = values
        .iter()
        .cloned()
        .enumerate()
        .map(|(i, v)| (v, i))
        .collect();
    (values, index)
}

fn order_direction(array: &Array2<f64>, direction: ClusterDirection) -> Vec<usize> {
    let n = direction.n(array); // Size along the direction
    if n > 2 {
        HierarchicalCluster::new(
            array,
            DistanceMetric::Euclidean,
            LinkageMethod::Ward,
            direction,
        )
        .leaves(LeafOrdering::ModularSmallest)
    } else {
        (0..n).collect()
    }
}

fn tenx_blue_colorscale() -> ColorScale {
    // See CELLRANGER-6596
    ColorScale::Vector(
        [
            (0.0, "#0044B2"),
            (0.3385, "#0071D9"),
            (0.6927, "#B3D4F4"),
            (1.0, "#F0F2F5"),
        ]
        .into_iter()
        .map(|(v, c)| ColorScaleElement(v, c.to_string()))
        .collect(),
    )
}

struct ClonotypeSpecificity {
    clonotypes: Vec<String>, // columns
    antigens: Vec<String>,   // rows
    median_specificity: Array2<f64>,
    clonotype_sizes: Vec<usize>,
}

impl ClonotypeSpecificity {
    fn new(antigen_specificity: CsvFile<AntigenSpecificityRow>) -> Result<Self> {
        Self::from_csv_rows(antigen_specificity.read()?)
    }
    fn from_csv_rows(mut rows: Vec<AntigenSpecificityRow>) -> Result<Self> {
        rows.sort_by(|a, b| a.sort_key().cmp(&b.sort_key()));

        let (clonotypes, clonotype_index) =
            make_index(rows.iter().filter_map(AntigenSpecificityRow::clonotype_id));
        let (antigens, antigen_index) = make_index(rows.iter().map(|r| r.antigen.clone()));

        let mut clonotype_sizes = vec![0; clonotypes.len()];
        let mut specificity_mat = Array2::zeros((antigens.len(), clonotypes.len()));

        for ((clonotype_id, antigen), group_iter) in &rows
            .into_iter()
            .chunk_by(|r| (r.clonotype_id(), r.antigen.clone()))
        {
            if let Some(clonotype_id) = clonotype_id {
                let specificites: Vec<_> =
                    group_iter.map(|r| r.antigen_specificity_score).collect();
                let row = antigen_index[&antigen];
                let col = clonotype_index[&clonotype_id];
                clonotype_sizes[col] = specificites.len();
                specificity_mat[[row, col]] = Data::new(specificites).median();
            }
        }

        Ok(ClonotypeSpecificity {
            clonotypes,
            antigens,
            median_specificity: specificity_mat,
            clonotype_sizes,
        })
    }

    // Only retain upto top 10 clonotypes for each antigen. The order is defined by
    // (Specificity > 50, Clonotype size)
    fn _top_clonotypes(self) -> Self {
        let mut retain_clonotypes = BTreeSet::new();
        for row in self.median_specificity.rows() {
            let column_indices = (0..row.len())
                .filter(|&i| row[i] > MIN_SPECIFICITY)
                .sorted_by_key(|&i| Reverse(self.clonotype_sizes[i]))
                .take(TOP_N);
            retain_clonotypes.extend(column_indices);
        }

        let retain_clonotypes: Vec<_> = retain_clonotypes.into_iter().collect();

        let mut median_specificity = Array2::zeros((self.antigens.len(), retain_clonotypes.len()));
        for ((row, col), value) in median_specificity.indexed_iter_mut() {
            *value = self.median_specificity[[row, retain_clonotypes[col]]];
        }

        let (clonotypes, clonotype_sizes) = retain_clonotypes
            .iter()
            .map(|&i| (self.clonotypes[i].clone(), self.clonotype_sizes[i]))
            .unzip();

        Self {
            clonotypes,
            antigens: self.antigens,
            median_specificity,
            clonotype_sizes,
        }
    }

    fn _hierarchical_reorder(self) -> Self {
        // Reorder rows and columns based on hierarchical clustering
        let row_order = order_direction(&self.median_specificity, ClusterDirection::Rows);
        let col_order = order_direction(&self.median_specificity, ClusterDirection::Columns);

        let antigens = row_order
            .iter()
            .map(|i| self.antigens[*i].clone())
            .collect();

        let (clonotypes, clonotype_sizes) = col_order
            .iter()
            .map(|i| (self.clonotypes[*i].clone(), self.clonotype_sizes[*i]))
            .unzip();

        let mut median_specificity = Array2::zeros(self.median_specificity.raw_dim());
        for ((row, col), value) in median_specificity.indexed_iter_mut() {
            *value = self.median_specificity[[row_order[row], col_order[col]]];
        }

        Self {
            antigens,
            clonotypes,
            clonotype_sizes,
            median_specificity,
        }
    }
    fn clustermap(self) -> Option<ChartWithHelp> {
        if self.clonotype_sizes.is_empty() || self.antigens.is_empty() {
            return None;
        }
        const MAX_GAPS: usize = 70;

        let total_clonotypes = self.clonotypes.len();

        let ClonotypeSpecificity {
            clonotypes,
            antigens,
            median_specificity,
            clonotype_sizes,
        } = self._top_clonotypes()._hierarchical_reorder();

        let clonotype_sizes: Vec<_> = clonotype_sizes.into_iter().map(|i| i.to_string()).collect();

        let n_clonotypes = clonotypes.len();

        let x = clonotypes
            .into_iter()
            .map(|cl| ClonotypeId::parse(&cl).unwrap().id.to_string())
            .collect();

        // Unnecessary duplication here since the text needs to be per z-value
        // but is only a function of clonotype. We could instead just embed number of
        // cells into "x" itself and show that on hover, but we don't want it to show up in
        // the x labels.
        let (z, text): (Vec<Vec<_>>, Vec<_>) = median_specificity
            .rows()
            .into_iter()
            .map(|r| (r.iter().copied().collect(), clonotype_sizes.clone()))
            .unzip();

        let mut data = serde_json::to_value(
            HeatMap::new(x, antigens, z)
                .reverse_scale(true)
                .hover_template(
                    "Clonotype: %{x} (%{text} cells)<br>Antigen: %{y}<br>Median specificity score: %{z}",
                )
                .color_bar(ColorBar::new().title("Antigen Specificity Score").outline_width(0))
                .color_scale(tenx_blue_colorscale())
                .name(""),
        )
        .unwrap();

        // Due to plotly API limitations. Should ideally be fixed upstream
        let map = data.as_object_mut().unwrap();
        map.insert("xgap".into(), i32::from(n_clonotypes < MAX_GAPS).into());
        map.insert("ygap".into(), 1.into());
        map.insert("text".into(), text.into());
        map.insert("zmin".into(), 0.0.into());
        map.insert("zmax".into(), 100.0.into());

        let mut plot = PlotlyChart::with_layout_and_data(
            Layout::new()
                .x_axis(
                    Axis::new()
                        .fixed_range(false)
                        .title(format!("Clonotype ID ({n_clonotypes} clonotypes)"))
                        .type_(AxisType::Category),
                )
                .y_axis(Axis::new().auto_margin(true).title("Antigens")),
            vec![data],
        );

        plot.config = json!({
            "displayModeBar": true,
            "modeBarButtons": [
                [
                    "toImage", "resetScale2d"
                ]
            ],
            "scrollZoom": false,
            "showAxisDragHandles": true,
            "staticPlot": false,
        });

        Some(ChartWithHelp {
            plot,
            help: TitleWithHelp {
                help: format!(
                    "The hierarchically-clustered heatmap shows the antigen specificity for \
                {n_clonotypes} clonotypes out of {total_clonotypes} clonotypes in this sample. \
                For each antigen with >0 cellular UMI counts, we select upto {TOP_N} clonotypes \
                (ordered by size) that have a median antigen specificity score >{MIN_SPECIFICITY}. \
                The set of all such clonotypes across the antigens is show in the x-axis."
                ),
                title: "Antigen Specificity Clustermap".into(),
            },
        })
    }
}

pub fn clonotype_specificity_heatmap(
    antigen_specificity: CsvFile<AntigenSpecificityRow>,
) -> Result<Option<ChartWithHelp>> {
    Ok(ClonotypeSpecificity::new(antigen_specificity)?.clustermap())
}

#[cfg(test)]
mod tests {
    use super::*;
    use cr_types::clonotype::ClonotypeId;

    #[test]
    fn test_no_clonotypes() {
        assert!(
            ClonotypeSpecificity::from_csv_rows(vec![AntigenSpecificityRow {
                barcode: "TTTGTCATCACCGGGT-1".into(),
                antigen: "BEAM01".into(),
                antigen_umi: 180,
                control: "BEAM03".into(),
                control_umi: 2,
                antigen_specificity_score: 99.852,
                mhc_allele: None,
                raw_clonotype_id: None,
                exact_subclonotype_id: None,
            }])
            .unwrap()
            .clustermap()
            .is_none()
        );
    }

    #[test]
    fn test_one_clonotype_antigen() {
        assert!(
            ClonotypeSpecificity::from_csv_rows(vec![AntigenSpecificityRow {
                barcode: "TTTGTCATCACCGGGT-1".into(),
                antigen: "BEAM01".into(),
                antigen_umi: 180,
                control: "BEAM03".into(),
                control_umi: 2,
                antigen_specificity_score: 99.852,
                mhc_allele: None,
                raw_clonotype_id: Some(
                    ClonotypeId {
                        id: 1,
                        sample_id: None,
                    }
                    .to_string()
                ),
                exact_subclonotype_id: Some("1".into()),
            }])
            .unwrap()
            .clustermap()
            .is_some()
        );
    }
}
