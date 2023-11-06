# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
"""Generate aggr web summary from the individual json data."""

__MRO__ = """
stage BUILD_MULTI_GRAPH_VIEW(
    in  json multi_graph,
    out svg  view,
    src py   "stages/multi/build_multi_graph_view",
)
"""

import json

import cellranger.multi.config as cfg


def main(args, outs):
    render = False
    if args.multi_graph is not None:
        with open(args.multi_graph) as f:
            data = json.load(f)
        graph = cfg.CrMultiGraph.from_json_val(data)
        render = graph.render_to_svg(outs.view)  # Only rendered for multiplexed runs for now

    if not render:
        outs.view = None
