# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
"""Generate aggr web summary from the individual json data."""

from __future__ import annotations

import json
from typing import TYPE_CHECKING

import martian

if TYPE_CHECKING:
    import cellranger.mro_types.filetypes as mro_filetypes

__MRO__ = """
stage BUILD_ANTIGEN_AGGR_WEB_SUMMARY(
    in  json antigen_histograms_json,
    in  json antigen_treemap_json,
    in  json clonotype_clustermap_json,
    out json antigen_tab,
    src py   "stages/aggregator/build_antigen_aggr_web_summary",
)
"""


class StageInputs(martian.Record):  # pylint: disable=too-few-public-methods
    """BUILD_ANTIGEN_AGGR_WEB_SUMMARY stage inputs."""

    antigen_histograms_json: mro_filetypes.json | None
    antigen_treemap_json: mro_filetypes.json | None
    clonotype_clustermap_json: mro_filetypes.json | None


class StageOutputs(martian.Record):  # pylint: disable=too-few-public-methods
    """BUILD_ANTIGEN_AGGR_WEB_SUMMARY stage outputs."""

    antigen_tab: mro_filetypes.json | None


ChunkInputs = StageInputs
ChunkOutputs = StageOutputs

ANTIGEN_HISTOGRAM = "antigen_histogram_plot"
ANTIGEN_TREEMAP = "antigen_treemap_plot"
ANTIGEN_CLONOTYPE_CLUSTERMAP = "clonotype_clustermap"


def main(args: StageInputs, outs: StageOutputs):
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
