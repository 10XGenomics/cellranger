#!/usr/bin/env python
#
# Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
#
# Calculate data required for rarefaction and extrapolation curves.
#

"""A stage for calculating clonotype diversity plot from clonotypes.csv.

Classes:
    Clonotypes: used to load information from clonotypes.csv for the purpose of
    plotting clonotype diversity plots.
"""


import json

import cellranger.stats as cr_stats

__MRO__ = """
stage CALCULATE_CLONOTYPE_DIVERSITY(
    in  json per_origin_clonotype_hist,
    in  json diversity_chart,
    out json plotly_diversity_chart,
    src py   "stages/vdj/clonotype_diversity",
)
"""

# import plotly.express as px
# print(px.colors.qualitative.D3)
PLOTLY_COLORS = [
    "#1f77b4",
    "#ff7f0e",
    "#2ca02c",
    "#d62728",
    "#9467bd",
    "#8c564b",
    "#e377c2",
    "#7f7f7f",
    "#bcbd22",
    "#17becf",
]

# TODO: This is the only element in the web summary created by python code because we
# don't have the rarefaction math implemented in Rust. Would be cleaner to move this
DIVERSITY_CHART_TEMPLATE = {
    "plot": {
        "config": {
            "displayModeBar": True,
            "staticPlot": False,
            "dragmode": "zoom",
            "modeBarButtons": [["toImage"]],
        },
        "data": [],
        "layout": {
            "showlegend": True,
            "margin": {"r": 40, "l": 60, "t": 25},
            "hovermode": "closest",
            "xaxis": {"title": "Number of Cells"},
            "yaxis": {"title": "Number of Unique Clonotypes"},
            "legend": dict(yanchor="top", y=0.99, xanchor="left", x=0.01, bgcolor="rgba(0,0,0,0)"),
        },
    },
    "help": {
        "helpText": "This plot displays rarefaction curves of clonotype diversity. Each curve is from 1 origin.",
        "title": "Clonotype diversity",
    },
}

# class Clonotypes(object):
#     """Clonotypes class used to load information from clonotypes.csv for the
#     purpose of plotting clonotype diversity plots.
#     Attributes:
#         clonotypes          (pandas dataframe): Pandas dataframe created
#                             by loading clonotypes.csv
#         hist                (list of int): An unsorted list of clonotype
#                             abundances
#         diversity           (Diversity): A Diversity object to calculate and
#                             save divesity curves
#     """

#     def __init__(self, clonotypes_csv):
#         self.clonotypes = pd.read_csv(clonotypes_csv)
#         # Extracting (clonotype_id, count) tuples
#         clon_unsorted = list(
#             self.clonotypes[["clonotype_id", "frequency"]].itertuples(index=False, name=None)
#         )
#         self.hist = [x[1] for x in clon_unsorted]
#         self.diversity = cr_stats.Diversity(self.hist)

#     def create_and_save_curves(self, outfile, rc_steps, ec_steps, max_n):
#         self.diversity.calc_rarefaction_curve()
#         self.diversity.calc_extrapolation_curve()
#         self.diversity.create_both_curves(outfile, rc_steps, ec_steps, max_n)

#     def can_draw_curve(self):
#         return self.diversity.is_diversity_curve_possible()


def main(args, outs):
    with open(args.per_origin_clonotype_hist) as f:
        per_origin_clonotype_hist = json.load(f)

    data = []
    for index, (origin, hist) in enumerate(per_origin_clonotype_hist):
        # In the very unlikely case of more than 10 origins, we will cycle through colors
        index = index % len(PLOTLY_COLORS)
        values = hist.values()
        data.extend(
            cr_stats.Diversity(values).calc_rarefaction_curve_plotly(
                origin, PLOTLY_COLORS[index], num_steps=50, stdev=False
            )
        )

    chart = DIVERSITY_CHART_TEMPLATE
    chart["plot"]["data"] = data

    with open(outs.plotly_diversity_chart, "w") as out_file:
        json.dump(chart, out_file, indent=4)
