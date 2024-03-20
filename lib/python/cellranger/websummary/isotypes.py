#!/usr/bin/env python
#
# Copyright (c) 2023 10X Genomics, Inc. All rights reserved.
#
"""Code to produce isotype control plots for the websummary."""
from __future__ import annotations

from collections import OrderedDict
from typing import TYPE_CHECKING

import altair as alt
import numpy as np
import pandas as pd
import scipy.sparse as sp_sparse
from scipy.cluster import hierarchy

import cellranger.altair_utils as alt_utils
from cellranger.rna.library import ANTIBODY_LIBRARY_TYPE, GENE_EXPRESSION_LIBRARY_TYPE
from cellranger.websummary.numeric_converters import round_floats_in_list

if TYPE_CHECKING:
    import cellranger.matrix as cr_matrix
    from cellranger.feature_ref import FeatureDef

alt.data_transformers.disable_max_rows()

LINK_HELP_TEXT = "<a href='https://10xgen.com/spaceranger-antibody-algo' target='_blank' title='Antibody Algorithms' rel='noopener noreferrer'>isotype normalized</a>"
ALT_LINK_HELP_TEXT = "<a href='https://10xgen.com/spaceranger-antibody-algo' target='_blank' title='Antibody Algorithms' rel='noopener noreferrer'>isotype normalization</a>"
# Number of clusters used to cluster FBc features to reorder the correlation matrix
NUM_HIERARCHICAL_CLUSTERS = 10


def get_antibody_feature_dict(
    raw_matrix: cr_matrix.CountMatrix, filtered_matrix: cr_matrix.CountMatrix
) -> tuple[OrderedDict, OrderedDict, int]:
    """Returns antibody feature names present in both matricies.

    Given a raw and filtered matrix, returns an ordered dictionary mapping each
    antibody feature name present in both matrices to its index in the respectively matrix.
    The isotype controls are maintained in the end of the OrderedDict

    Args:
        raw_matrix (cr_matrix.CountMatrix): the raw matrix containing antibody feature data
        filtered_matrix (cr_matrix.CountMatrix): the filtered matrix containing antibody feature data

    Returns:
        tuple[OrderedDict, OrderedDict, int]: Two OrderedDicts representing the antibody feature
        indices for the raw and filtered matrices, respectively; and the number of control features
        in the OrderedDict. If the antibody feature type is not present in a matrix, the both OrderedDicts are empty.
    """

    def get_antibody_dict(
        feature_defs: list[FeatureDef],
        non_control_feature_ids: set[bytes],
        control_feature_ids: set[bytes],
    ) -> OrderedDict:
        """Helper function to extract the a dictionary of common antibody features from a Featurereference.

        Ordered dict is ordered by non-control features sorted by id followed by control features sorted
        by id

        Args:
            feature_defs (list[FeatureDef]): the list of feature def to extract the antibody dictionary from
            non_control_feature_ids (set[bytes]): set of ids of non control features to extract
            control_feature_ids (set[bytes]):set of ids of control features to extract

        Returns:
            OrderedDict: an OrderedDict mapping each antibody feature name to its index in the matrix
        """
        control_ab_features = sorted(
            [
                f
                for f in feature_defs
                if f.feature_type == ANTIBODY_LIBRARY_TYPE and f.id in control_feature_ids
            ],
            key=lambda x: x.id,
        )
        non_control_ab_features = sorted(
            [
                f
                for f in feature_defs
                if f.feature_type == ANTIBODY_LIBRARY_TYPE and f.id in non_control_feature_ids
            ],
            key=lambda x: x.id,
        )
        return OrderedDict(
            (f.id.decode(), f.index) for f in non_control_ab_features + control_ab_features
        )

    raw_antibody_dict = OrderedDict()
    filtered_antibody_dict = OrderedDict()
    num_control_features = 0

    if ANTIBODY_LIBRARY_TYPE in {
        f.feature_type for f in raw_matrix.feature_ref.feature_defs
    }.intersection({f.feature_type for f in filtered_matrix.feature_ref.feature_defs}):
        # Figuring out control and non-control ABs in both the filtered
        # and the raw matrix
        filtered_feature_defs = filtered_matrix.feature_ref.feature_defs
        non_control_filtered_abs_ids = {
            f.id
            for f in filtered_feature_defs
            if f.feature_type == ANTIBODY_LIBRARY_TYPE
            and f.tags.get("isotype_control", "FALSE") != "TRUE"
        }
        control_filtered_abs_ids = {
            f.id
            for f in filtered_feature_defs
            if f.feature_type == ANTIBODY_LIBRARY_TYPE
            and f.tags.get("isotype_control", "FALSE") == "TRUE"
        }

        raw_feature_defs = raw_matrix.feature_ref.feature_defs
        non_control_raw_abs_ids = {
            f.id
            for f in raw_feature_defs
            if f.feature_type == ANTIBODY_LIBRARY_TYPE
            and f.tags.get("isotype_control", "FALSE") != "TRUE"
        }
        control_raw_abs_ids = {
            f.id
            for f in raw_feature_defs
            if f.feature_type == ANTIBODY_LIBRARY_TYPE
            and f.tags.get("isotype_control", "FALSE") == "TRUE"
        }

        non_control_feature_ids = non_control_filtered_abs_ids.intersection(non_control_raw_abs_ids)
        control_feature_ids = control_filtered_abs_ids.intersection(control_raw_abs_ids)

        num_control_features = len(control_feature_ids)
        raw_antibody_dict = get_antibody_dict(
            raw_feature_defs, non_control_feature_ids, control_feature_ids
        )
        filtered_antibody_dict = get_antibody_dict(
            filtered_feature_defs, non_control_feature_ids, control_feature_ids
        )

    return raw_antibody_dict, filtered_antibody_dict, num_control_features


def get_ab_ab_correlation_matrix(
    raw_matrix: cr_matrix.CountMatrix,
    filtered_matrix: cr_matrix.CountMatrix,
) -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    """Compute the correlation matrix between all FBC features from a raw and isotype-normalized matrix.

    Args:
        raw_matrix (cr_matrix.CountMatrix): a raw matrix
        filtered_matrix (cr_matrix.CountMatrix): a filtered isotype-normalized matrix

    Returns:
        tuple[pd.DataFrame, pd.DataFrame, list[str]]: This includes the following
            - raw_cross_corr (pd.DataFrame): A Pandas DataFrame containing the correlation matrix between raw FBC features
            - filtered_cross_corr (pd.DataFrame): A Pandas DataFrame containing the correlation matrix between isotype
                normalized FBC features
            - order_list (list[str]): Order to plot FBc in
    """
    raw_antibody_dict, filtered_antibody_dict, num_control_abs = get_antibody_feature_dict(
        raw_matrix, filtered_matrix
    )
    # Create the correlation matrix

    def create_wide_corr_df(
        matrix: cr_matrix.CountMatrix, antibody_dict: OrderedDict[str, int], num_control_abs: int
    ) -> tuple[pd.DataFrame, int]:
        """Get a df of correlation of features in the dict in the matrix.

        Args:
            matrix (cr_matrix.CountMatrix): A countmatrix
            antibody_dict (OrderedDict[str, int]): Ordered dict contating AB feature name: index in matrix
            num_control_abs (int): number of control ABs seen

        Returns:
            pd.DataFrame: DataFrame of correlation in wide form
            num_non_zero_control_abs: int with the number control ABs with non-zero UMIs
        """
        ab_matrix = np.log(1 + matrix.m[list(antibody_dict.values()), :].toarray())
        row_sums = ab_matrix.sum(axis=1).flatten()
        nonzero_rows = np.flatnonzero(row_sums)
        num_non_zero_control_abs = np.count_nonzero(row_sums[-num_control_abs:])
        ab_matrix = ab_matrix[nonzero_rows, :]
        num_ab = ab_matrix.shape[0]
        corr_mat = np.corrcoef(ab_matrix, ab_matrix)[:num_ab, -num_ab:]
        antibody_names = list(antibody_dict)
        nonzero_antibody_names = [antibody_names[j] for j in nonzero_rows]
        cross_corr = pd.DataFrame(
            corr_mat,
            columns=nonzero_antibody_names,
            index=nonzero_antibody_names,
        )
        return (cross_corr, num_non_zero_control_abs)

    def get_long_cross_corr_df(x_corr: pd.DataFrame) -> pd.DataFrame:
        """Convert wide cross correlation df to a long df to plot.

        Args:
            x_corr (pd.DataFrame): wide cross correlation df

        Returns:
            pd.DataFrame: Long df to enable easy plotting
        """
        nonzero_antibody_names = list(x_corr.columns)
        cross_corr = x_corr.reset_index()
        cross_corr = pd.melt(
            cross_corr,
            id_vars=["index"],
            value_vars=nonzero_antibody_names,
            var_name="Antibody(y)",
            value_name="Corr",
        )
        cross_corr = cross_corr.rename(columns={"index": "Antibody(x)"})
        cross_corr["Corr"] = round_floats_in_list(cross_corr["Corr"])
        return cross_corr

    raw_cross_corr_wide, _ = create_wide_corr_df(raw_matrix, raw_antibody_dict, num_control_abs)
    filtered_cross_corr_wide, num_non_zero_control_abs = create_wide_corr_df(
        filtered_matrix, filtered_antibody_dict, num_control_abs
    )

    # If the correlation df is empty, skip hierarchical clustering.
    order_list = list(filtered_cross_corr_wide.columns)
    if not filtered_cross_corr_wide.iloc[
        :-num_non_zero_control_abs, :-num_non_zero_control_abs
    ].empty:
        linkage = hierarchy.linkage(
            filtered_cross_corr_wide.iloc[:-num_non_zero_control_abs, :-num_non_zero_control_abs],
            method="average",
        )
        idx_to_cluster_array = hierarchy.fcluster(
            linkage, NUM_HIERARCHICAL_CLUSTERS, criterion="maxclust"
        )
        order_list = list(
            filtered_cross_corr_wide.columns[np.argsort(idx_to_cluster_array)]
        ) + list(filtered_cross_corr_wide.columns[-num_non_zero_control_abs:])

    raw_cross_corr = get_long_cross_corr_df(raw_cross_corr_wide)
    filtered_cross_corr = get_long_cross_corr_df(filtered_cross_corr_wide)
    return raw_cross_corr, filtered_cross_corr, order_list


def gex_fbc_correlation_matrix(matrix: cr_matrix.CountMatrix):
    """Compute the correlation matrix between gene expression and FBC features from SingleGenomeAnalysis data.

    Parameters:
        matrix (CountMatrix): A isotype-normalized FBC matrix

    Returns:
        cross_corr (pd.DataFrame): A Pandas DataFrame containing the correlation matrix between gene expression and FBC features.
    """
    feature_defs = matrix.feature_ref.feature_defs
    gene_expression_features = list(
        filter(lambda fd: fd.feature_type == GENE_EXPRESSION_LIBRARY_TYPE, feature_defs)
    )

    antibody_features = list(
        filter(lambda fd: fd.feature_type == ANTIBODY_LIBRARY_TYPE, feature_defs)
    )
    # Get the intersection of the FBC and GEX dictionaries by name
    gene_expression_dict = {f.name: f.index for f in gene_expression_features}
    antibody_dict = {f.name: f.index for f in antibody_features}
    intersection = {
        name: gene_expression_dict[name] for name in gene_expression_dict if name in antibody_dict
    }
    # Create two new dictionaries based on the intersection of the two original dictionaries
    gene_expression_intersect = {
        name: index for name, index in gene_expression_dict.items() if name in intersection
    }
    antibody_intersect = {
        name: index for name, index in antibody_dict.items() if name in intersection
    }

    # Sort the dictionaries
    antibody_intersect = {k: antibody_intersect[k] for k in sorted(antibody_intersect)}
    gene_expression_intersect = {
        k: gene_expression_intersect[k] for k in sorted(gene_expression_intersect)
    }

    ab_matrix = sp_sparse.csr_matrix.log1p(matrix.m[list(antibody_intersect.values()), :])
    gex_matrix = sp_sparse.csr_matrix.log1p(matrix.m[list(gene_expression_intersect.values()), :])

    num_gex = gex_matrix.shape[0]
    corr_mat = np.corrcoef(gex_matrix.toarray(), ab_matrix.toarray())[:num_gex, -num_gex:]
    cross_corr = pd.DataFrame(
        corr_mat,
        columns=list(gene_expression_intersect.keys()),
        index=list(gene_expression_intersect.keys()),
    )

    cross_corr = cross_corr.reset_index()
    cross_corr = pd.melt(
        cross_corr,
        id_vars=["index"],
        value_vars=list(gene_expression_intersect.keys()),
        var_name="Gene",
        value_name="Corr",
    )
    cross_corr = cross_corr.rename(columns={"index": "Antibody"})
    cross_corr["Corr"] = round_floats_in_list(cross_corr["Corr"])
    return cross_corr


def fbc_isotype_correlation(
    matrix: cr_matrix.CountMatrix, filtered_barcodes: list, log_transform: bool
):
    """Calculates the correlation between counts of barcodes expressing isotype controls non-isotype controls.

    Parameters:
    -----------
    matrix: CountMatrix
        A `CountMatrix` object representing the raw matrix.
    filtered_barcodes: list
        A list of filtered barcodes to subset the raw matrix
    log_transform: bool
        log transform the data or not

    Returns:
    --------
    pandas.DataFrame:
        A DataFrame with two columns, "Isotype Antibody Counts" and "Non-Isotype Antibody Counts", containing the total
        counts of cells expressing isotype controls and non-isotype controls, respectively.
    """
    # Filter the raw matrix to only tissue covered barcodes
    matrix = matrix.select_barcodes_by_seq(list(filtered_barcodes))

    feature_defs = matrix.feature_ref.feature_defs

    antibody_features = list(
        filter(lambda fd: fd.feature_type == ANTIBODY_LIBRARY_TYPE, feature_defs)
    )
    isotype_index = [
        feat.index for feat in antibody_features if feat.tags.get("isotype_control") == "TRUE"
    ]
    antibody_index = [
        feat.index for feat in antibody_features if feat.tags.get("isotype_control") == "FALSE"
    ]
    if log_transform:
        isotype_matrix = sp_sparse.csr_matrix.log1p(matrix.m[isotype_index, :])
        antibody_matrix = sp_sparse.csr_matrix.log1p(matrix.m[antibody_index, :])
    else:
        isotype_matrix = matrix.m[isotype_index, :]
        antibody_matrix = matrix.m[antibody_index, :]

    isotype_bc_counts = round_floats_in_list(isotype_matrix.sum(axis=0).tolist()[0])
    antibody_bc_counts = round_floats_in_list(antibody_matrix.sum(axis=0).tolist()[0])

    isotype_antibody_counts_df = pd.DataFrame(
        {
            "Isotype Antibody": isotype_bc_counts,
            "Non-Isotype Antibody": antibody_bc_counts,
        }
    )

    return isotype_antibody_counts_df


def make_gex_fbc_correlation_heatmap(matrix: cr_matrix.CountMatrix):
    """Create a heatmap showing the Pearson correlation counts from gene expression features and antibody features.

    Parameters:
        matrix (cr_matrix.CountMatrix): an isotype-normalized matrix

    Returns:
        corr_plot_data (dict): A dictionary containing the plot data and accompanying help text and title.
    """
    corr_mat = gex_fbc_correlation_matrix(matrix)
    corr_plot = (
        alt.Chart(corr_mat)
        .mark_rect()
        .encode(
            x="Antibody",
            y=alt.Y("Gene", scale=alt.Scale(reverse=True)),
            color=alt.Color(
                "Corr:Q",
                scale=alt.Scale(scheme="blueorange", domain=[-1, 1]),
            ),
            tooltip=["Antibody", "Gene", "Corr"],
        )
        .properties(width=300, height=300)
        .interactive()
    )

    corr_plot = alt_utils.chart_to_json(corr_plot)
    corr_plot_help = {
        "helpText": f"Pearson correlations between log transformed gene expression counts and log transformed {LINK_HELP_TEXT} antibody counts "
        f"for all gene:antibody pairs that have matching names.",
        "title": "Gene:Antibody Correlations",
    }
    corr_plot_data = {"gex_fbc_correlation_plot": {"help": corr_plot_help, "spec": corr_plot}}
    return corr_plot_data


def make_ab_ab_correlation_heatmap(
    raw_matrix: cr_matrix.CountMatrix,
    filtered_matrix: cr_matrix.CountMatrix,
):
    """Create a heatmap showing the Pearson correlation between all antibody features for a given sample.

    Parameters:
        filtered_matrix (CountMatrix): An isotype-normalized FBC matrix
        raw_matrix (CountMatrix): A raw fbc count matrix.

    Returns:
        raw_normalized_heatmap_data: A dictionary containing the plot data and accompanying help text and title.
    """
    # make the correlation matricies
    raw_corr_mat, filtered_corr_mat, order_list = get_ab_ab_correlation_matrix(
        raw_matrix, filtered_matrix
    )
    merged_corr = pd.merge(
        raw_corr_mat,
        filtered_corr_mat,
        how="inner",
        on=["Antibody(x)", "Antibody(y)"],
        suffixes=(" Raw", " Filtered"),
    )

    raw_corr_plot = (
        alt.Chart(merged_corr, title=alt.TitleParams("Raw Antibody Counts", anchor="middle"))
        .mark_rect()
        .encode(
            x=alt.X("Antibody(x)", sort=order_list),
            y=alt.Y("Antibody(y)", scale=alt.Scale(reverse=True), sort=order_list),
            color=alt.Color(
                "Corr Raw:Q",
                title="Corr",
                scale=alt.Scale(scheme="blueorange", domain=[-1, 1]),
            ),
            tooltip=["Antibody(x)", "Antibody(y)", "Corr Raw"],
        )
        .properties(width=320, height=300)
    )

    filtered_corr_plot = (
        alt.Chart(merged_corr, title=alt.TitleParams("Normalized Antibody Counts", anchor="middle"))
        .mark_rect()
        .encode(
            x=alt.X("Antibody(x)", sort=order_list),
            y=alt.Y("Antibody(y)", scale=alt.Scale(reverse=True), sort=order_list),
            color=alt.Color(
                "Corr Filtered:Q",
                title="Corr",
                scale=alt.Scale(scheme="blueorange", domain=[-1, 1]),
            ),
            tooltip=["Antibody(x)", "Antibody(y)", "Corr Filtered"],
        )
        .properties(width=320, height=300)
    )

    raw_normalized_heatmap = alt.hconcat(raw_corr_plot, filtered_corr_plot)

    # build the plot data
    raw_normalized_heatmap = alt_utils.chart_to_json(raw_normalized_heatmap)
    raw_normalized_heatmap_help = {
        "helpText": f"Heatmap of Pearson correlation coefficients of raw (left) and isotype normalized antibody (right) counts across spots under tissue. "
        f"Antibodies are clustered using the correlation structures and then ordered according this clustering. "
        f"Compare the raw count correlation plot with the normalized count plot to see the effect of {ALT_LINK_HELP_TEXT}.",
        "title": "Antibody Correlations",
    }
    raw_normalized_heatmap_data = {
        "raw_normalized_heatmap": {
            "help": raw_normalized_heatmap_help,
            "spec": raw_normalized_heatmap,
        }
    }

    return raw_normalized_heatmap_data


# need to use the raw matrix
def make_fbc_isotype_correlation_scatter(
    raw_matrix: cr_matrix.CountMatrix, filtered_barcodes: list, log_transform: bool = False
):
    """Generates a scatter plot showing the correlation between isotype and non-isotype antibodies.

    Parameters:
        raw_matrix (cr_matrix.CountMatrix): raw unnormalized FBC matrix
        filtered_barcodes (list): a list of barcodes under tissue
        log_transform (bool): whether to log transform the data

    Returns:
        dict: A dictionary containing the plot data in JSON format, along with help text and a title.
    """
    isotype_corr = fbc_isotype_correlation(
        matrix=raw_matrix, filtered_barcodes=filtered_barcodes, log_transform=log_transform
    )

    # calculate 99th percentile for each column
    isotype_quantile = isotype_corr.quantile(0.99)
    isotype_corr["Above 99th percentile"] = (
        (isotype_corr["Isotype Antibody"] > isotype_quantile["Isotype Antibody"])
        | (isotype_corr["Non-Isotype Antibody"] > isotype_quantile["Non-Isotype Antibody"])
    ).astype(int)

    # filter out rows where either column is above the 99th percentile value
    isotype_corr_trimmed = isotype_corr[
        (isotype_corr["Isotype Antibody"] <= isotype_quantile["Isotype Antibody"])
        & (isotype_corr["Non-Isotype Antibody"] <= isotype_quantile["Non-Isotype Antibody"])
    ]

    r_squared = round(
        isotype_corr_trimmed["Isotype Antibody"].corr(isotype_corr_trimmed["Non-Isotype Antibody"])
        ** 2,
        ndigits=3,
    )

    chart_title = alt.TitleParams(
        f"r² = {r_squared}",
        anchor="middle",
    )
    isotype_corr_plot = (
        alt.Chart(isotype_corr, title=chart_title)
        .mark_point()
        .encode(
            x=alt.X("Isotype Antibody", title="Sum of Isotype Antibody Raw Counts"),
            y=alt.Y("Non-Isotype Antibody", title="Sum of Non-Isotype Antibody Raw Counts"),
            color=alt.Color("Above 99th percentile:N", legend=None),
            opacity=alt.value(0.3),
            tooltip=["Isotype Antibody", "Non-Isotype Antibody"],
        )
    )

    transform_text = "log transformed" if log_transform else "raw"

    isotype_corr_plot = alt_utils.chart_to_json(isotype_corr_plot)
    isotype_corr_plot_plot_help = {
        "helpText": f"Total {transform_text} counts from isotype antibody features and non-isotype antibody features for tissue-associated barcodes. "
        f"A high pearson correlation (r) validates the usage of isotypes for {ALT_LINK_HELP_TEXT}. "
        f"Barcodes in the top 1% of isotype antibody or non-isotype antibody counts are colored in orange and discarded before calculating r².",
        "title": "Isotype vs. Non-Isotype Antibodies",
    }
    isotype_corr_plot_plot_data = {
        "isotype_corr_plot": {"help": isotype_corr_plot_plot_help, "spec": isotype_corr_plot}
    }
    return isotype_corr_plot_plot_data, r_squared
