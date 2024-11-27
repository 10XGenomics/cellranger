#!/usr/bin/env python
#
# Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
#
"""This file contains plots in the analysis tab that are directly related to the outputs of the SC_RNA_ANALYZER pipeline for.

gene expression and/or antibody capture, e.g. TSNE plots, clustering plots, differential expression table, etc
splitting things up in this way allows turing repo to utilize TSNE/clustering/diffexp plots without pulling in a lot of transitive dependencies
"""

from __future__ import annotations

import numpy as np

import cellranger.rna.library as rna_library
import cellranger.webshim.constants.gex as ws_gex_constants
import cellranger.websummary.plotly_tools as pltly
from cellranger.analysis.analysis_types import DifferentialExpressionWithFeatures
from cellranger.analysis.clustering import AB_PREFIX, GEX_PREFIX
from cellranger.analysis.singlegenome import (
    PROJECTION_TITLE,
    TSNE_NAME,
    UMAP_NAME,
    Projection,
    SingleGenomeAnalysis,
)
from cellranger.websummary.helpers import get_projection_key
from cellranger.websummary.numeric_converters import round_floats_in_list
from cellranger.websummary.react_components import (
    ClusteringData,
    Clusterings,
    ClusteringSelector,
    DifferentialExpressionTableValue,
    SharedCoordinatePlotCollection,
)

DIFFEXP_TABLE_HELP = {
    "helpText": "The differential expression analysis seeks to find, for each cluster, features that are more highly expressed in that cluster relative to the rest of the sample. "
    "Here a differential expression test was performed between each cluster and the rest of the sample for each feature. "
    "The Log2 fold-change (L2FC) is an estimate of the log2 ratio of expression in a cluster to that in all other cells. "
    "A value of 1.0 indicates 2-fold greater expression in the cluster of interest. "
    "The p-value is a measure of the statistical significance of the expression difference and is based on a negative binomial test. "
    "The p-value reported here has been adjusted for multiple testing via the Benjamini-Hochberg procedure. "
    "In this table you can click on a column to sort by that value. "
    "Also, in this table features were filtered by (Mean UMI counts > 1.0) and the top N features by L2FC for each cluster were retained. "
    "Features with L2FC < 0 or adjusted p-value >= 0.10 were grayed out. "
    "The number of top features shown per cluster, N, is set to limit the number of table entries shown to 10,000; N=10,000/K^2 where K is the number of clusters. "
    "N can range from 1 to 50. "
    "For the full table, please refer to the 'differential_expression.csv' files produced by the pipeline.",
    "title": "Top Features by Cluster (Log2 fold-change, p-value)",
}

SPATIAL_DIFFEXP_TABLE_HELP = {
    "helpText": "The differential expression analysis seeks to find, for each cluster, features that are more highly expressed in that cluster relative to the rest of the sample. "
    "Here a differential expression test was performed between each cluster and the rest of the sample for each feature. "
    "The Log2 fold-change (L2FC) is an estimate of the log2 ratio of expression in a cluster to that in all other spots. "
    "A value of 1.0 indicates 2-fold greater expression in the cluster of interest. "
    "The p-value is a measure of the statistical significance of the expression difference and is based on a negative binomial test. "
    "The p-value reported here has been adjusted for multiple testing via the Benjamini-Hochberg procedure. "
    "In this table you can click on a column to sort by that value. "
    "Also, in this table features were filtered by (Mean UMI counts > 1.0) and the top N features by L2FC for each cluster were retained. "
    "Features with L2FC < 0 or adjusted p-value >= 0.10 were grayed out. "
    "The number of top features shown per cluster, N, is set to limit the number of table entries shown to 10,000; N=10,000/K^2 where K is the number of clusters. "
    "N can range from 1 to 50. "
    "For the full table, please refer to the 'differential_expression.csv' files produced by the pipeline.",
    "title": "Top Features by Cluster (Log2 fold-change, p-value)",
}

# exactly the same as DIFFEXP_TABLE_HELP except UMI->object
INSITU_DIFFEXP_TABLE_HELP = {
    "helpText": "The differential expression analysis seeks to find, for each cluster, features that are more highly expressed in that cluster relative to the rest of the sample. "
    "Here a differential expression test was performed between each cluster and the rest of the sample for each feature. "
    "The Log2 fold-change (L2FC) is an estimate of the log2 ratio of expression in a cluster to that in all other cells. "
    "A value of 1.0 indicates 2-fold greater expression in the cluster of interest. "
    "The p-value is a measure of the statistical significance of the expression difference and is based on a negative binomial test. "
    "The p-value reported here has been adjusted for multiple testing via the Benjamini-Hochberg procedure. "
    "In this table you can click on a column to sort by that value. "
    "Also, in this table features were filtered by (Mean object counts > 1.0) and the top N features by L2FC for each cluster were retained. "
    "Features with L2FC < 0 or adjusted p-value >= 0.10 were grayed out. "
    "The number of top features shown per cluster, N, is set to limit the number of table entries shown to 10,000; N=10,000/K^2 where K is the number of clusters. "
    "N can range from 1 to 50. "
    "For the full table, please refer to the 'differential_expression.csv' files produced by the pipeline.",
    "title": "Top Features by Cluster (Log2 fold-change, p-value)",
}

SPATIAL_TSNE_CLUSTERING_PLOT_HELP = {
    "data": [
        [
            "",
            [
                "(left) These are the assignments of each spot-barcode to clusters by an automated clustering algorithm. "
                "The clustering groups together spots that have similar expression profiles. "
                "In this plot, spots are colored according to their cluster assignment and projected on to the tissue image. "
                "Only spots under tissue are used in the clustering algorithm.",
                "(right) Spots are colored by clustering assignment and shown in t-SNE space. "
                "The axes correspond to the 2-dimensional embedding produced by the t-SNE algorithm. "
                "In this space, pairs of spots that are close to each other have more similar gene expression profiles than spots that are distant from each other. ",
                # "The display is limited to a random subset of spots.", # TODO: We should set a max # of cells in insitu and subsample if necessary
            ],
        ]
    ],
    "title": "Clustering",
}

INSITU_TSNE_CLUSTERING_PLOT_HELP = {
    "data": [
        [
            "",
            [
                "(left) These are the assignments of each cell to clusters by an automated clustering algorithm. "
                "The clustering groups together cells that have similar expression profiles. "
                "In this plot, cells are colored according to their cluster assignment and plotted in their spatial location. "
                "Only cells with a Nucleus detected in the DAPI stain are used in the clustering algorithm.",
                "(right) Cells are colored by clustering assignment and shown in t-SNE space. "
                "The axes correspond to the 2-dimensional embedding produced by the t-SNE algorithm. "
                "In this space, pairs of cells that are close to each other have more similar gene expression profiles than cells that are distant from each other. ",
                # "The display is limited to a random subset of spots.", # TODO: We should set a max # of cells in insitu and subsample if necessary
            ],
        ]
    ],
    "title": "Clustering",
}

UMI_TSNE_PLOT = "umi_tsne_plot"


def projection_clustering_plot_help(projection: Projection):
    return {
        "data": [
            [
                "",
                [
                    "(left) Shown here are the total UMI counts for each cell-barcode. "
                    "Cells with greater UMI counts likely have higher RNA content than cells with fewer UMI counts. "
                    f"The axes correspond to the 2-dimensional embedding produced by the {PROJECTION_TITLE[projection]} algorithm. "
                    "In this space, pairs of cells that are close to each other have more similar gene expression profiles than cells that are distant from each other. "
                    "The display is limited to a random subset of cells.",
                    "(right) These are the assignments of each cell-barcode to clusters by an automated clustering algorithm. "
                    "The clustering groups together cells that have similar expression profiles. "
                    f"The axes correspond to the 2-dimensional embedding produced by the {PROJECTION_TITLE[projection]} algorithm. "
                    "In this space, pairs of cells that are close to each other have more similar gene expression profiles than cells that are distant from each other. "
                    "The display is limited to a random subset of cells. "
                    "Please use Loupe Browser to view the entire dataset.",
                ],
            ]
        ],
        "title": f"{PROJECTION_TITLE[projection]} Projection",
    }


def projection_layout_config(projection: Projection):
    return {
        "xaxis": {
            "type": "linear",
            "title": f"{PROJECTION_TITLE[projection]}1",
            "showline": False,
            "zeroline": True,
            "fixedrange": False,
        },
        "yaxis": {
            "type": "linear",
            "title": f"{PROJECTION_TITLE[projection]}2",
            "showline": False,
            "zeroline": True,
            "fixedrange": False,
        },
        "margin": {"t": 30},  # needed to keep screenshots from hitting the top
        "hovermode": "closest",
    }


def diffexp_table_from_analysis(
    diffexp,
    clustering,
    analysis: SingleGenomeAnalysis,
):
    """Show diffexp table for the given clustering."""
    if not analysis:
        raise ValueError("A analysis object must be provided.")
    if not clustering:
        raise ValueError("A clustering must be provided.")

    n_clusters: int = clustering.clusters.max()
    cluster_names = [f"Cluster {i + 1}" for i in range(n_clusters)]
    diffexp_with_features = DifferentialExpressionWithFeatures.from_cmatrix_and_diffexp(
        cmatrix=analysis.matrix, diffexp=diffexp
    )
    return diffexp_table(diffexp_with_features=diffexp_with_features, cluster_names=cluster_names)


def diffexp_table(
    diffexp_with_features: DifferentialExpressionWithFeatures, cluster_names
):  # pylint: disable=too-many-locals
    """Show diffexp table for the given clustering."""
    # Limit the number of entries in the DE table
    diffexp = diffexp_with_features.diffexp
    feature_id_to_name = dict(
        zip(diffexp_with_features.feature_ids, diffexp_with_features.feature_names)
    )
    feature_id_to_int = {x: i for i, x in enumerate(diffexp_with_features.feature_ids)}

    n_clusters = len(cluster_names)
    n_genes = int(np.floor(float(ws_gex_constants.MAX_DE_TABLE_ENTRIES) / (n_clusters**2)))
    if n_genes < 1:
        n_genes = 1
    elif n_genes > ws_gex_constants.MAX_TOP_N_GENES:
        n_genes = ws_gex_constants.MAX_TOP_N_GENES

    columns = [
        {
            "Header": "Feature",
            "columns": [
                {"Header": "ID", "accessor": "feature.id"},
                {"Header": "Name", "accessor": "feature.fn"},
            ],
        }
    ]

    # Get the union of top DE genes
    top_genes = set()
    for i, cluster_name in enumerate(cluster_names):
        means = diffexp.data[:, 0 + 3 * i]
        log2fcs = diffexp.data[:, 1 + 3 * i]

        keep_indices = np.flatnonzero(means >= ws_gex_constants.TOP_DE_GENES_MIN_MEAN)
        top_gene_indices = keep_indices[log2fcs[keep_indices].argsort()[::-1]][:n_genes]

        for j in top_gene_indices:
            top_genes.add(diffexp_with_features.feature_ids[j])

        columns.append(
            {
                "Header": cluster_name,
                "columns": [
                    {
                        "Header": "L2FC",
                        "accessor": f"c{i + 1}.l",
                        "greyedout": f"c{i + 1}.g",
                    },
                    {
                        "Header": "p-value",
                        "accessor": f"c{i + 1}.p",
                        "greyedout": f"c{i + 1}.g",
                    },
                ],
            }
        )

    table = []
    for gene_id in top_genes:
        i = feature_id_to_int[gene_id]
        gene_name = feature_id_to_name[gene_id]

        row = {"feature": {"fn": gene_name, "id": gene_id}}

        for j in range(n_clusters):
            log2fc = diffexp.data[i, 1 + (3 * j)]
            adj_p_value = diffexp.data[i, 2 + (3 * j)]

            greyed = bool(log2fc <= 0 or adj_p_value >= ws_gex_constants.PVALUE_DEEMPHASIS_CUTOFF)

            cluster_num = f"c{j + 1}"
            row[cluster_num] = DifferentialExpressionTableValue(adj_p_value, log2fc, greyed)

        table.append(row)

    # Sort by log2fc, descending, in first cluster
    if n_clusters > 0:
        table = sorted(table, key=lambda row: row["c1"].log2_fold_change, reverse=True)

    return {
        "columns": columns,
        "data": table,
    }


def sort_analysis_clusterings(clusterings):
    return sorted(clusterings.items(), key=lambda k_v: k_v[1].global_sort_key)


def _get_unit_and_plt_type(is_spatial: bool = False, is_insitu: bool = False):
    """Get the plot type and unit to use."""
    assert not (is_insitu and is_spatial)
    if is_insitu:
        unit = "Cells"
        plt_type = "scatter"
    elif is_spatial:
        unit = "Spots"
        plt_type = "scatter"
    else:
        unit = "Cells"
        plt_type = "scattergl"
    return unit, plt_type


def analysis_by_clustering_helper(
    analysis: SingleGenomeAnalysis,
    projection: Projection,
    spatial: bool = False,
    insitu: bool = False,
    library_type: str = rna_library.GENE_EXPRESSION_LIBRARY_TYPE,
    projection_marker_opacity: float = 0.9,
    projection_marker_size: float = 4,
    downsample_mask: np.ndarray | None = None,
):
    """Constructs a ClusteringSelector object from a SingleGenomeAnalyis.

    Args:
        analysis (SingleGenomeAnalysis): secondary analysis
        spatial (bool, optional): is it a spatial web summary. Defaults to False.
        insitu (bool, optional): is it an insitu web summary. Defaults to False.
        library_type (str, optional): Library type. Defaults to
            'Gene Expression'.
        projection_marker_opacity (float, optional): plotting marker opacity.
            Defaults to 0.9.
        projection_marker_size (float, optional): plotting marker size.
            Defaults to 4.
        projection (Projection): Which manifold projection to use in
            plotting: 'tsne' or 'umap'.
        downsample_mask (Optional[np.ndarray], optional): downsample points.
            Defaults to None.

    Raises:
        NotImplementedError: if projection is not umap or tsne

    Returns:
        ClusteringSelector|None: clustering selector object or None if no
            projection could be found
    """
    assert projection in (TSNE_NAME, UMAP_NAME)

    # pylint: disable=too-many-function-args,too-many-locals
    key = get_projection_key(library_type, 2)

    unit, plt_type = _get_unit_and_plt_type(spatial, insitu)

    if projection == TSNE_NAME:
        if key not in analysis.tsne:
            return None
        proj_coords = analysis.get_tsne(key=key).transformed_tsne_matrix
    elif projection == UMAP_NAME:
        if key not in analysis.umap:
            return None
        proj_coords = analysis.get_umap(key=key).transformed_umap_matrix

    layout = projection_layout_config(projection=projection).copy()
    layout["title"] = f"{PROJECTION_TITLE[projection]} Projection of {unit} by Clustering"
    layout["hover_mode"] = "closest"

    if downsample_mask is not None:
        proj_coords = proj_coords[downsample_mask, :]

    proj_data = SharedCoordinatePlotCollection(
        proj_coords[:, 0],
        proj_coords[:, 1],
        pltly.PLOT_CONFIG,
        layout,
        plt_type,
        {"opacity": projection_marker_opacity, "size": projection_marker_size},
    )

    diff_exp_tables = Clusterings()
    # order clustering by order: graph, kmeans=2, =3 etc
    clusterings = sort_analysis_clusterings(analysis.clusterings)
    for clustering_key, clustering in clusterings:
        if library_type == rna_library.ANTIBODY_LIBRARY_TYPE:
            if not clustering_key.startswith(AB_PREFIX):
                continue
        elif library_type == rna_library.GENE_EXPRESSION_LIBRARY_TYPE:
            if not clustering_key.startswith(GEX_PREFIX):
                continue
        # Add differential expression table for clustering
        diftblclust = ClusteringData(
            key=clustering_key,
            clustering=clustering,
            data=diffexp_table_from_analysis(
                analysis.differential_expression[clustering_key], clustering, analysis
            ),
        )
        diff_exp_tables.add_clustering(diftblclust)

        # Optionally downsample the clustered data before plotting
        if downsample_mask is not None:
            clustering = clustering._replace(clusters=clustering.clusters[downsample_mask])
        # Add t-sne plot for clustering
        proj_data.add_clustering(clustering_key, clustering)

    assert not (insitu and spatial)

    if insitu:
        clust_select = ClusteringSelector(
            INSITU_TSNE_CLUSTERING_PLOT_HELP, INSITU_DIFFEXP_TABLE_HELP
        )
    elif spatial:
        clust_select = ClusteringSelector(
            SPATIAL_TSNE_CLUSTERING_PLOT_HELP, SPATIAL_DIFFEXP_TABLE_HELP
        )
    else:
        clust_select = ClusteringSelector(
            projection_clustering_plot_help(projection=projection), DIFFEXP_TABLE_HELP
        )

    clust_select.right_plots = proj_data
    clust_select.tables = diff_exp_tables
    return clust_select


def analysis_by_clustering(
    sample_data,
    spatial=False,
    library_type=rna_library.GENE_EXPRESSION_LIBRARY_TYPE,
    *,
    projection: Projection,
):
    """Get the tSNE/UMAP (colored by clustering) and diffexp table for each clustering."""
    if sample_data is None:
        return None

    analysis = sample_data.get_analysis(SingleGenomeAnalysis)
    if analysis is None:
        return None
    return analysis_by_clustering_helper(
        analysis,
        projection=projection,
        spatial=spatial,
        library_type=library_type,
        projection_marker_size=(
            3 if analysis.matrix.bcs_dim > ws_gex_constants.MAX_WEBSHIM_BCS_DIM else 4
        ),
    )


def umi_on_projection_helper(
    analysis: SingleGenomeAnalysis,
    spatial: bool = False,
    library_type=None,
    tag_type=None,
    *,
    projection: Projection,
):
    """Get the t-SNE/UMAP projections for a given library and retun the left-side plot as json."""
    if library_type is None:

        if projection == TSNE_NAME:
            proj_coords = analysis.get_tsne().transformed_tsne_matrix
        else:
            proj_coords = analysis.get_umap().transformed_umap_matrix
        reads_per_bc = analysis.matrix.get_counts_per_bc()
    else:
        key = get_projection_key(library_type, 2)
        if projection == TSNE_NAME:
            if key not in analysis.tsne:
                return None
            proj_coords = analysis.get_tsne(key=key).transformed_tsne_matrix
        else:
            assert projection == UMAP_NAME
            if key not in analysis.umap:
                return None
            proj_coords = analysis.get_umap(key=key).transformed_umap_matrix
        if tag_type is None:
            matrix = analysis.matrix.select_features_by_type(library_type)
            reads_per_bc = matrix.get_counts_per_bc()
        else:
            matrix = analysis.matrix.select_features_by_type_and_tag(library_type, tag_type)
            reads_per_bc = matrix.get_counts_per_bc()

    color = get_umi_color(reads_per_bc)
    unit, plt_type = _get_unit_and_plt_type(spatial)
    data = [
        {
            "name": unit,
            "x": round_floats_in_list(proj_coords[:, 0]),
            "y": round_floats_in_list(proj_coords[:, 1]),
            "type": plt_type,
            "mode": "markers",
            "marker": {
                "opacity": 0.9,
                "size": 3 if analysis.matrix.bcs_dim > ws_gex_constants.MAX_WEBSHIM_BCS_DIM else 4,
                "color": color,
                "colorscale": "Jet",
                "colorbar": {"title": "UMI counts"},
            },
            "hovertemplate": "(%{x:.3f}, %{y:.3f}) <br> " + "UMI counts: %{text}<extra></extra>",
            "text": [f"{reads:,d}" for reads in reads_per_bc],
        }
    ]
    # Note: the help text has been included in tsne_cluster plot
    layout = projection_layout_config(projection=projection).copy()
    layout["title"] = f"{PROJECTION_TITLE[projection]} Projection of {unit} Colored by UMI Counts"
    umi_projection_plot = {
        "config": pltly.PLOT_CONFIG,
        "layout": layout,
        "data": data,
    }
    return umi_projection_plot


def get_umi_color(umi_per_bc):
    """Color barcodes by their UMI counts."""
    vmin, vmax = np.percentile(umi_per_bc, ws_gex_constants.TSNE_TOTALCOUNTS_PRCT_CLIP)
    vmin = int(vmin)
    vmax = int(vmax)
    color = [min(vmax, max(vmin, int(v))) for v in umi_per_bc]
    return color


def umi_on_projection_plot(
    sample_data,
    *,
    spatial,
    library_type,
    tag_type,
    projection: Projection,
):
    """UMI count on tSNE plot."""
    if sample_data is None:
        return None

    analysis = sample_data.get_analysis(SingleGenomeAnalysis)

    if analysis is None:
        return None

    return umi_on_projection_helper(
        analysis,
        spatial=spatial,
        library_type=library_type,
        tag_type=tag_type,
        projection=projection,
    )


########################### Only used in multi stages #######################################
def projection_diffexp_plots_from_path(analysis_dir, *, projection: Projection, library_type=None):
    """Given the base analysis directory for a single genome analysis h5 file,.

    Get the tSNE/UMAPs (left-right plots colored by clustering and umi counts)
    and diffexp table for each clustering
    """
    if analysis_dir is None:
        return None

    analysis = SingleGenomeAnalysis.load_default_format(
        analysis_dir, method="pca", projections=(projection,)
    )

    if analysis is None:
        return None

    # right (colored by clustering) tsne plot and gene expression table
    plots = analysis_by_clustering_helper(
        analysis, projection=projection, spatial=False, library_type=library_type
    )

    # left (colored by umi count) tsne plot
    if plots:
        left_plots = umi_on_projection_helper(
            analysis, spatial=False, library_type=library_type, projection=projection
        )
        if left_plots:
            plots.left_plots = left_plots

    return plots


def umi_on_projection_from_path(
    analysis_dir, projection: Projection, library_type=None, tag_type=None
):
    """Given the base analysis directory for a single genome analysis h5,.

    Get the t-SNE/UMAP colored by umi counts, optionally filtered for a specific library type
    """
    # pylint: disable=missing-function-docstring
    if analysis_dir is None:
        return None

    analysis = SingleGenomeAnalysis.load_default_format(
        analysis_dir, method="pca", projections=(projection,)
    )

    if analysis is None:
        return None

    # left (colored by umi count) tsne plot
    return umi_on_projection_helper(
        analysis, spatial=False, library_type=library_type, tag_type=tag_type, projection=projection
    )
