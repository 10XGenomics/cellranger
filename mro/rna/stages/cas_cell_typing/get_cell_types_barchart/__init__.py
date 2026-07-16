#
# Copyright (c) 2024 10X Genomics, Inc. All rights reserved.
#

"""Generate cell types interactive barchart."""

import json
import textwrap

import altair as alt
import martian
import polars as pl

from cellranger.altair_utils import chart_to_json
from cellranger.cell_typing.broad_tenx.cas_postprocessing import LOW_UMI_BARCODE_KEY
from cellranger.cell_typing.common_cell_typing import (
    BARCODE_KEY,
    DISPLAY_MAP_KEY,
    FINE_CELL_TYPES_KEY,
)

__MRO__ = """
stage GET_CELL_TYPES_BARCHART(
    in  csv  cell_types,
    out json cell_type_interactive_bar_chart,
    src py   "stages/cas_cell_typing/get_cell_types_barchart",
) using (
    volatile = strict,
)
"""


# Keys of things shown in the graph
NUM_BARCODES_KEY = "num_barcodes"
NUM_BC_IN_COARSE_CELL_TYPE_KEY = "num_bc_in_coarse_cell_type"
FRACTION_OF_COARSE_CELL_TYPE_KEY = "fraction_of_coarse_cell_type"

# Hyperparameters of figure
CELL_TYPE_WIDTH_IN_LEGEND = 50
COARSE_CELL_TYPE_LOWER_BOUND = 10
MAX_SUBTYPES_TO_SHOW = 5
MIN_FRACTION_TO_LIST = 0.000051
NUMBER_OF_BARCODES_LABEL = "Number of Barcodes"
CELL_TYPES_LABEL = "Cell Type"
CELL_SUB_TYPE_LABEL = "Sub-type"
FRACTION_BARCODES_LABEL = "Fraction Barcodes"
RIGHT_CHART_HEADING_FONTSIZE = 18
RIGHT_CHART_TABLE_FONTSIZE = 14
RIGHT_CHART_CELL_TYPES_WIDTH = 175
RIGHT_CHART_CELL_TYPES_X = 0
RIGHT_CHART_CELL_TYPE_STRING_LIMIT = RIGHT_CHART_CELL_TYPES_WIDTH
RIGHT_CHART_CELL_TYPES_ALIGNMENT = "left"
RIGHT_CHART_FRAC_BCS_X = 0
RIGHT_CHART_FRAC_BCS_WIDTH = 30
RIGHT_CHART_FRAC_BCS_ALIGNMENT = "right"

# Schema of postprocessed cell annotation out
cell_typing_csv_schema = {
    BARCODE_KEY: pl.String,
    DISPLAY_MAP_KEY: pl.String,
    FINE_CELL_TYPES_KEY: pl.String,
}


def main(args, outs):  # pylint: disable=too-many-locals
    if not args.cell_types:
        martian.clear(outs)
        return

    title_coarse_cell_type_key = f"title_{DISPLAY_MAP_KEY}"

    pl.set_random_seed(seed=7)
    lazy_df = pl.scan_csv(args.cell_types, schema_overrides=cell_typing_csv_schema).select(
        [BARCODE_KEY, DISPLAY_MAP_KEY, FINE_CELL_TYPES_KEY]
    )
    df = (
        lazy_df.group_by((FINE_CELL_TYPES_KEY, DISPLAY_MAP_KEY))
        .agg(
            pl.len().alias(NUM_BARCODES_KEY),
        )
        .sort(pl.col(NUM_BARCODES_KEY), descending=True)
        .collect()
    )

    df = df.lazy().with_columns(pl.col(DISPLAY_MAP_KEY).alias(title_coarse_cell_type_key)).collect()

    list_of_dicts = (
        df.lazy()
        .group_by(DISPLAY_MAP_KEY)
        .agg(pl.col(NUM_BARCODES_KEY).sum().alias(NUM_BC_IN_COARSE_CELL_TYPE_KEY))
        .collect()
        .to_dicts()
    )
    map_of_total_coarce_cell_types = {
        x[DISPLAY_MAP_KEY]: x[NUM_BC_IN_COARSE_CELL_TYPE_KEY] for x in list_of_dicts
    }
    df = df.with_columns(
        [
            pl.struct([NUM_BARCODES_KEY, DISPLAY_MAP_KEY])
            .map_elements(
                lambda x: x[NUM_BARCODES_KEY] / map_of_total_coarce_cell_types[x[DISPLAY_MAP_KEY]],
                return_dtype=pl.Float64,
            )
            .alias(FRACTION_OF_COARSE_CELL_TYPE_KEY)
        ]
    )

    df = df.with_columns(
        [
            pl.col(FINE_CELL_TYPES_KEY).map_elements(
                lambda x: "\n".join(textwrap.wrap(x, width=CELL_TYPE_WIDTH_IN_LEGEND)),
                return_dtype=pl.String,
            )
        ]
    )

    df = df.with_columns(
        [
            pl.lit(CELL_SUB_TYPE_LABEL).alias(CELL_SUB_TYPE_LABEL),
            pl.lit(FRACTION_BARCODES_LABEL).alias(FRACTION_BARCODES_LABEL),
        ]
    )

    default_value = (
        df.group_by(DISPLAY_MAP_KEY)
        .agg(pl.sum(NUM_BARCODES_KEY))
        .sort(pl.col(NUM_BARCODES_KEY), descending=True)
        .select(pl.col(DISPLAY_MAP_KEY).first())
        .to_series()
        .to_list()[0]
    )

    brush = alt.selection_point(
        fields=[DISPLAY_MAP_KEY],
        empty=False,
        value=default_value,
        clear="dblclick",
        toggle=False,
    )

    hover_brush = alt.selection_point(on="mouseover", empty=False, value=default_value)

    base_chart = alt.Chart(df)

    base_left_chart = base_chart.transform_aggregate(
        groupby=[DISPLAY_MAP_KEY], num_barcodes=f"sum({NUM_BARCODES_KEY})"
    )

    num_cell_types = (
        df.group_by(DISPLAY_MAP_KEY)
        .agg(pl.col("num_barcodes").sum())
        .filter(pl.col("num_barcodes") > COARSE_CELL_TYPE_LOWER_BOUND)
        .height
    )

    chart_height = max(15 * num_cell_types, 250)
    left_chart = (
        base_left_chart.transform_filter(f"datum.num_barcodes > {COARSE_CELL_TYPE_LOWER_BOUND}")
        .transform_calculate(
            sort_priority=f"datum.{DISPLAY_MAP_KEY} === '{LOW_UMI_BARCODE_KEY}' ? 1 : 0"
        )
        .transform_window(
            sort_index="rank()",
            sort=[
                alt.SortField("sort_priority", order="ascending"),
                alt.SortField("num_barcodes", order="descending"),
            ],
        )
        .mark_bar(
            stroke="black",
            color="#0071D9",
            cursor=alt.Cursor("pointer"),
        )
        .encode(
            y=alt.Y(f"{DISPLAY_MAP_KEY}:N")
            .sort(alt.EncodingSortField(field="sort_index", order="ascending"))
            .axis(title=CELL_TYPES_LABEL, labelFontSize=12, titleFontSize=14),
            x=alt.X("num_barcodes:Q").axis(
                title=NUMBER_OF_BARCODES_LABEL, labelFontSize=12, titleFontSize=14
            ),
            tooltip=[
                alt.Tooltip(DISPLAY_MAP_KEY, title=CELL_TYPES_LABEL),
                alt.Tooltip("num_barcodes", title=NUMBER_OF_BARCODES_LABEL),
            ],
            opacity=alt.condition(brush, alt.value(1), alt.value(0.5)),
            strokeWidth=alt.condition(hover_brush, alt.value(2), alt.value(0)),
        )
        .add_params(brush, hover_brush)
    ).properties(height=chart_height)

    transformed_chart = base_chart.transform_filter(brush)

    heading = (
        transformed_chart.transform_calculate(major_type_tooltip_text="'Major cell type'")
        .transform_aggregate(
            groupby=[title_coarse_cell_type_key, "major_type_tooltip_text"],
        )
        .mark_text(
            align=RIGHT_CHART_CELL_TYPES_ALIGNMENT,
            fontSize=RIGHT_CHART_HEADING_FONTSIZE,
            fontWeight="bold",
            x=RIGHT_CHART_CELL_TYPES_X,
        )
        .encode(
            text=f"{title_coarse_cell_type_key}:N",
            tooltip=alt.Tooltip("major_type_tooltip_text:N"),
        )
        .properties(width=RIGHT_CHART_CELL_TYPES_WIDTH)
    )

    cell_type_subheading = (
        transformed_chart.transform_calculate(
            sub_type_tooltip_text="'Cell sub-types that compose the major cell type'"
        )
        .transform_aggregate(
            groupby=[CELL_SUB_TYPE_LABEL, "sub_type_tooltip_text"],
        )
        .mark_text(
            align=RIGHT_CHART_CELL_TYPES_ALIGNMENT,
            fontSize=RIGHT_CHART_TABLE_FONTSIZE,
            fontWeight="bold",
            x=RIGHT_CHART_CELL_TYPES_X,
        )
        .encode(
            text=f"{CELL_SUB_TYPE_LABEL}:N",
            tooltip=alt.Tooltip("sub_type_tooltip_text:N"),
        )
        .properties(width=RIGHT_CHART_CELL_TYPES_WIDTH)
    )

    frac_bcs_subheading = (
        transformed_chart.transform_calculate(
            frac_bc_tooltip_text="'The fraction of barcodes in the major cell type'"
        )
        .transform_aggregate(
            groupby=[FRACTION_BARCODES_LABEL, "frac_bc_tooltip_text"],
        )
        .mark_text(
            align=RIGHT_CHART_FRAC_BCS_ALIGNMENT,
            fontSize=RIGHT_CHART_TABLE_FONTSIZE,
            fontWeight="bold",
            x=RIGHT_CHART_FRAC_BCS_X,
        )
        .encode(
            text=f"{FRACTION_BARCODES_LABEL}:N",
            tooltip=alt.Tooltip("frac_bc_tooltip_text:N"),
        )
        .properties(width=RIGHT_CHART_FRAC_BCS_WIDTH)
    )

    text_chart = (
        transformed_chart.transform_aggregate(
            groupby=[FINE_CELL_TYPES_KEY],
            frac_barcodes_in_subtype=f"sum({FRACTION_OF_COARSE_CELL_TYPE_KEY})",
        )
        .transform_window(
            row_number="row_number()",
            sort=[
                alt.SortField("frac_barcodes_in_subtype", order="descending"),  # Primary sort
                alt.SortField(FINE_CELL_TYPES_KEY, order="ascending"),  # Secondary sort
            ],
        )
        .transform_filter(f"datum.row_number < {MAX_SUBTYPES_TO_SHOW}")
        .transform_filter(f"datum.frac_barcodes_in_subtype >= {MIN_FRACTION_TO_LIST}")
        .encode(y=alt.Y("row_number:O", axis=None))
    )

    cell_type_name = cell_type_subheading & text_chart.mark_text(
        x=RIGHT_CHART_CELL_TYPES_X,
        limit=RIGHT_CHART_CELL_TYPE_STRING_LIMIT,
        align=RIGHT_CHART_CELL_TYPES_ALIGNMENT,
        fontSize=RIGHT_CHART_TABLE_FONTSIZE,
    ).encode(
        text=alt.Text(f"{FINE_CELL_TYPES_KEY}:N"),
        tooltip=[
            alt.Tooltip(f"{FINE_CELL_TYPES_KEY}:N", title=CELL_SUB_TYPE_LABEL),
            alt.Tooltip(
                "frac_barcodes_in_subtype:O",
                format="0.02%",
                title=FRACTION_BARCODES_LABEL,
            ),
        ],
    ).properties(
        width=RIGHT_CHART_CELL_TYPES_WIDTH,
    )  # Combine data tables
    num_subtype_bcs = frac_bcs_subheading & text_chart.mark_text(
        align=RIGHT_CHART_FRAC_BCS_ALIGNMENT,
        fontSize=RIGHT_CHART_TABLE_FONTSIZE,
        x=RIGHT_CHART_FRAC_BCS_X,
    ).encode(
        text=alt.Text(
            "frac_barcodes_in_subtype:O",
            format="0.02%",
        ),
    ).properties(
        width=RIGHT_CHART_FRAC_BCS_WIDTH
    )
    text = alt.vconcat(
        heading,
        alt.hconcat(cell_type_name, num_subtype_bcs, spacing=10).resolve_axis(x="independent"),
    )

    chart = (
        alt.hconcat(left_chart, text, spacing=100)
        .configure_concat(spacing=10)
        .resolve_legend(color="independent")
        .configure_view(strokeWidth=0)
    )

    with open(outs.cell_type_interactive_bar_chart, "w") as f:
        json.dump(chart_to_json(chart), f)
