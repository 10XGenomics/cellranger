# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
"""Generate aggr web summary from the individual json data."""

__MRO__ = """
stage BUILD_ANTIGEN_AGGR_WEB_SUMMARY(
    in  json antigen_histograms_json,
    in  json antigen_treemap_json,
    in  json clonotype_clustermap_json,
    out json antigen_tab,
    src py   "stages/aggregator/build_antigen_aggr_web_summary",
)
"""

import json

ANTIGEN_HISTOGRAM = "antigen_histogram_plot"
ANTIGEN_TREEMAP = "antigen_treemap_plot"
ANTIGEN_CLONOTYPE_CLUSTERMAP = "clonotype_clustermap"


def main(args, outs):
    antigen_tab = {}

    if args.antigen_histograms_json is not None:
        with open(args.antigen_histograms_json) as f:
            antigen_histograms_json = json.load(f)
        antigen_tab[ANTIGEN_HISTOGRAM] = antigen_histograms_json

    if args.antigen_treemap_json is not None:
        with open(args.antigen_treemap_json) as f:
            antigen_treemap_json = json.load(f)
        antigen_tab[ANTIGEN_TREEMAP] = antigen_treemap_json

    if args.clonotype_clustermap_json is not None:
        with open(args.clonotype_clustermap_json) as f:
            clonotype_clustermap_json = json.load(f)
        antigen_tab[ANTIGEN_CLONOTYPE_CLUSTERMAP] = clonotype_clustermap_json

    if len(antigen_tab.keys()) == 0:
        outs.antigen_tab = None
        return

    with open(outs.antigen_tab, "w") as f:
        json.dump(antigen_tab, f, indent=4)
