# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
"""Code to make the web summary for AGGR runs of (Space|Cell) Ranger."""


from __future__ import annotations

import numpy as np

import cellranger.rna.library as rna_library
import cellranger.websummary.sample_properties as wsp
from cellranger import utils as cr_util
from cellranger.analysis.singlegenome import (
    PROJECTION_TITLE,
    TSNE_NAME,
    Projection,
    SingleGenomeAnalysis,
)
from cellranger.targeted.targeted_constants import TARGETING_METHOD_TL
from cellranger.webshim.common import load_sample_data
from cellranger.webshim.constants.shared import CELLRANGER_COMMAND_NAME, PIPELINE_AGGR
from cellranger.websummary import plotly_tools as pltly
from cellranger.websummary.analysis_tab_core import projection_layout_config
from cellranger.websummary.helpers import get_projection_key
from cellranger.websummary.metrics import (
    MetricAnnotations,
    SpatialAggrMetricAnnotations,
    SpatialTargetedAggrMetricAnnotations,
    SpatialTemplateLigationAggrMetricAnnotations,
    TargetedAggrMetricAnnotations,
    TemplateLigationAggrMetricAnnotations,
)
from cellranger.websummary.numeric_converters import round_floats_in_list
from cellranger.websummary.react_components import (
    PlotWithHeader,
)
from cellranger.websummary.react_summarize import write_html_file
from cellranger.websummary.spatial_utils import SPATIAL_COMMAND_NAME
from cellranger.websummary.web_summary_builder import build_web_summary_html_sc_and_aggr


def library_on_projection_plot(sample_data, gg_id_to_name, *, projection: Projection):
    """Get the tSNE/UMAP colored by batch; only makes sense in aggregataions."""
    if sample_data is None:
        return None

    analysis = sample_data.get_analysis(SingleGenomeAnalysis)
    if analysis is None:
        return None

    # get batch names
    bcs = analysis.matrix.bcs
    library_names = cr_util.bcs_suffices_to_names(bcs, gg_id_to_name)

    library_types = analysis.matrix.get_library_types()
    if rna_library.GENE_EXPRESSION_LIBRARY_TYPE in library_types:
        key = get_projection_key(rna_library.GENE_EXPRESSION_LIBRARY_TYPE, 2)
    elif rna_library.ANTIBODY_LIBRARY_TYPE in library_types:
        key = get_projection_key(rna_library.ANTIBODY_LIBRARY_TYPE, 2)

    if projection == TSNE_NAME:
        projection_coordinates = analysis.get_tsne(key=key).transformed_tsne_matrix
    else:
        projection_coordinates = analysis.get_umap(key=key).transformed_umap_matrix
    layout = projection_layout_config(projection=projection).copy()

    data = []
    for library in gg_id_to_name.values():
        indices = np.where(np.array(library_names) == library)[0]
        if len(indices) > 0:
            # Plots with legend values >~20 characters seem to goof plotly's layout.  To avoid such
            # issues, we'll rename any id with >20 characters.
            name = library
            if len(name) > 15:
                name = name[:6] + "..." + name[-6:]
            data.append(
                {
                    "name": name + f" ({len(indices)})",
                    "x": round_floats_in_list(list(projection_coordinates[indices, 0])),
                    "y": round_floats_in_list(list(projection_coordinates[indices, 1])),
                    "indices": list(indices),
                    "type": "scattergl",
                    "mode": "markers",
                    "marker": {"opacity": 0.6, "size": 4},
                }
            )
    assert "title" not in layout
    batch_projection_plot = {
        "config": pltly.PLOT_CONFIG,
        "layout": layout,
        "data": data,
    }
    return batch_projection_plot


def build_web_summary_html_aggr(
    filename,
    sample_properties,
    gg_id_to_name_map,
    sample_data_paths,
    projection: Projection,
    sample_defs=None,
):
    """Build a web summary file for an AGGR Run.

    Args:
        filename: output file name
        sample_properties: Instance of AggrCountSampleProperties
        gg_id_to_name_map: dictionary to map gem group ids in barcode to library names
        sample_data_paths: Instance of SampleDataPaths
        projection: Projection used in websummary (TSNE_NAME or UMAP_NAME)
        sample_defs: Map of sample defs passed in from PARSE_CSV

    Returns:
        None
    """
    web_sum_data = build_web_summary_data_aggr(
        sample_properties, gg_id_to_name_map, sample_data_paths, projection, sample_defs
    )

    write_html_file(filename, web_sum_data)


def build_web_summary_data_aggr(
    sample_properties,
    gg_id_to_name_map,
    sample_data_paths,
    projection: Projection,
    sample_defs=None,
):
    """Build a web summary file for an AGGR Run.

    Args:
        sample_properties: Instance of AggrCountSampleProperties
        gg_id_to_name_map: dictionary to map gem group ids in barcode to library names
        sample_data_paths: Instance of SampleDataPaths
        projection: Tuple including all available projections
        sample_defs: Map of sample defs from PARSE_CSV

    Returns:
        web_sum_data: Web summary data
    """
    assert isinstance(sample_properties, wsp.AggrCountSampleProperties)
    # add metadata info that will be needed to draw batch tsne and determine targeting method

    sample_data = load_sample_data(
        sample_properties,
        sample_data_paths,
        (projection,),
    )
    if sample_properties.is_spatial:
        command = SPATIAL_COMMAND_NAME
        if sample_properties.is_targeted:
            if sample_data.targeting_method == TARGETING_METHOD_TL:
                metadata = SpatialTemplateLigationAggrMetricAnnotations()
            else:
                metadata = SpatialTargetedAggrMetricAnnotations()
        else:
            metadata = SpatialAggrMetricAnnotations()
    else:
        command = CELLRANGER_COMMAND_NAME
        if sample_properties.is_targeted:
            if sample_data.targeting_method == TARGETING_METHOD_TL:
                metadata = TemplateLigationAggrMetricAnnotations()
            else:
                metadata = TargetedAggrMetricAnnotations()
        else:
            metadata = MetricAnnotations(intron_mode_alerts=sample_properties.include_introns)

    web_sum_data = build_web_summary_html_sc_and_aggr(
        sample_properties,
        sample_data_paths,
        PIPELINE_AGGR,
        metadata,
        command,
        projection=projection,
        sample_defs=sample_defs,
    )

    clustering_plot = library_on_projection_plot(
        sample_data, gg_id_to_name_map, projection=projection
    )
    clustering_plot_title = f"{PROJECTION_TITLE[projection]} Projection Colored by Sample ID"

    if clustering_plot:
        web_sum_data.summary_tab["aggr_batch_projection"] = PlotWithHeader(
            clustering_plot_title,
            clustering_plot,
            help_text="Number of cells for each sample is provided in parentheses.",
        )

    return web_sum_data
