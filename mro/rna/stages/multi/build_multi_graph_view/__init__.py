#!/usr/bin/env python
#
# Copyright (c) 2023 10X Genomics, Inc. All rights reserved.
#
"""Creating the multi graph SVG for websummary."""

import os
from subprocess import PIPE, CalledProcessError, Popen

from cellranger.fast_utils import MultiGraph

__MRO__ = """
stage BUILD_MULTI_GRAPH_VIEW(
    in  json multi_graph,
    out svg  view,
    src py   "stages/multi/build_multi_graph_view",
)
"""


def main(args, outs):
    if args.multi_graph is None:
        outs.view = None
        return

    create_graphviz_digraph(MultiGraph.from_path(args.multi_graph)).draw(outs.view)


class GraphVizProps:
    # pylint: disable=too-few-public-methods
    """Graphviz property constants."""

    RANKSEP = "0.5"
    NODESEP = "0.75"
    ORDERING = "out"
    FILLEDSTYLE = "filled"
    NODESHAPE = "rectangle"
    RANKSAME = "same"


class GraphVizColors:
    # pylint: disable=too-few-public-methods
    """Graphviz color constants."""

    BLACK = "black"
    WHITE = "white"
    GREY = "lightgrey"


class AGraph:
    """Minimal implementation of the pygraphviz.AGraph class.

    Only the methods we need for the MultiGraph SVG.
    """

    def __init__(self, *_, **kwargs):
        kwargs.pop("strict", None)
        kwargs.pop("directed", None)

        self.graph_attr = kwargs
        self.node_attr = {}
        self.subgraphs = {}
        self.nodes = {}
        self.edges = {}

    def add_subgraph(self, *, name="", **kwargs) -> "AGraph":
        subg = AGraph(**kwargs)
        self.subgraphs[name] = subg
        return subg

    def add_node(self, name, **kwargs):
        self.nodes[name] = kwargs

    def add_edge(self, src, dst):
        """Add an edge from src to dst in the graph."""
        if src not in self.edges:
            self.edges[src] = []
        self.edges[src].append(dst)

    def string(self) -> str:
        """Generate a dot format representation of the graph, ready for graphviz consumption."""
        g_str = ['digraph "" {']
        if attrs := self.graph_attr:
            g_str += ["graph ["]
            for prop, value in attrs.items():
                g_str += [f'{prop}="{value}",']
            g_str += ["];"]

        # global node attributes
        g_str += ["node ["]
        if label := self.node_attr.pop("label", r"\N"):
            g_str += [f'label="{label}",']

        if self.node_attr:
            for prop, value in self.node_attr.items():
                g_str += [f'{prop}="{value}",']

        g_str += ["];"]

        # add subgraphs
        for name, subgraph in self.subgraphs.items():
            g_str += [f"subgraph {name} {{"]
            if attrs := subgraph.graph_attr:
                g_str += ["graph ["]
                for prop, value in attrs.items():
                    g_str += [f'{prop}="{value}",']
                g_str += ["];"]

            for node_name, props in subgraph.nodes.items():
                flat_value = ",".join([f'{key}="{value}"' for key, value in props.items()])
                g_str += [f'"{node_name}" [{flat_value}];']

            g_str += ["}"]

        for src, dsts in self.edges.items():
            for dst in dsts:
                g_str += [f"{src}->{dst}"]

        for subgraph in self.subgraphs.values():
            if edges := subgraph.edges:
                for src, dsts in edges.items():
                    for dst in dsts:
                        g_str += [f"{src}->{dst}"]

        g_str += ["}"]

        return "\n".join(g_str)

    def draw(self, fname):
        """Run graphviz and generate the graph in the SVG format."""
        g = self.string()
        dot_cmd = os.path.abspath(
            os.path.join(os.path.dirname(__file__), "../../../../../lib/bin/dot")
        )

        cmd = [dot_cmd, "-Tsvg", "-o", fname]
        with Popen(cmd, stdin=PIPE) as dot:
            dot.communicate(input=g.encode("utf-8"))
            if dot.returncode != 0:
                raise CalledProcessError(dot.returncode, cmd)


def create_graphviz_digraph(multi_graph: MultiGraph) -> AGraph:
    """Convert MultiGraph to a graphviz `digraph`.

    This is what shows up in the experimental design tab of the multi websummary.
    """
    # pylint: disable=too-many-locals

    g = AGraph(
        strict=False,
        directed=True,
        ranksep=GraphVizProps.RANKSEP,
        nodesep=GraphVizProps.NODESEP,
        ordering=GraphVizProps.ORDERING,
    )
    g.node_attr.update(
        style=GraphVizProps.FILLEDSTYLE,
        color=GraphVizColors.BLACK,
        fillcolor=GraphVizColors.WHITE,
        shape=GraphVizProps.NODESHAPE,
    )

    ## The samples cluster
    g_samp = g.add_subgraph(
        name="cluster_samples",
        style=GraphVizProps.FILLEDSTYLE,
        color=GraphVizColors.GREY,
        rank=GraphVizProps.RANKSAME,
    )
    g_samp.add_node(
        "Samples",
        style=GraphVizProps.FILLEDSTYLE,
        color=GraphVizColors.GREY,
        fillcolor=GraphVizColors.GREY,
    )

    ## The tags cluster
    g_tags: AGraph | None = None
    if multi_graph.is_multiplexed():
        g_tags = g.add_subgraph(
            name="cluster_tags",
            style=GraphVizProps.FILLEDSTYLE,
            color=GraphVizColors.GREY,
            rank=GraphVizProps.RANKSAME,
        )

        if multi_graph.is_cmo_multiplexed():
            g_tags.add_node(
                "CMO Tags",
                style=GraphVizProps.FILLEDSTYLE,
                color=GraphVizColors.GREY,
                fillcolor=GraphVizColors.GREY,
            )
        elif multi_graph.is_hashtag_multiplexed():
            g_tags.add_node(
                "Hashtags",
                style=GraphVizProps.FILLEDSTYLE,
                color=GraphVizColors.GREY,
                fillcolor=GraphVizColors.GREY,
            )
        elif multi_graph.is_rtl_multiplexed():
            g_tags.add_node(
                "Probe Barcode IDs",
                style=GraphVizProps.FILLEDSTYLE,
                color=GraphVizColors.GREY,
                fillcolor=GraphVizColors.GREY,
            )
        elif multi_graph.is_oh_multiplexed():
            g_tags.add_node(
                "OCM Barcode IDs",
                style=GraphVizProps.FILLEDSTYLE,
                color=GraphVizColors.GREY,
                fillcolor=GraphVizColors.GREY,
            )
        else:
            raise ValueError("Expected to be multiplexed.")

    ## The gem wells cluster
    g_gw = g.add_subgraph(
        name="cluster_gem_wells",
        style=GraphVizProps.FILLEDSTYLE,
        color=GraphVizColors.GREY,
        rank=GraphVizProps.RANKSAME,
    )
    g_gw.add_node(
        "GEM wells",
        style=GraphVizProps.FILLEDSTYLE,
        color=GraphVizColors.GREY,
        fillcolor=GraphVizColors.GREY,
    )

    # TODO: Not generalized in multi gem well
    for i, (sample_id, fingerprints) in enumerate(multi_graph.sample_info_for_graphviz()):
        samp_node_id = f"sample_{i + 1}"
        g_samp.add_node(samp_node_id, id=samp_node_id, label=sample_id)
        for j, fingerprint in enumerate(fingerprints):
            gw_node_id = f"gem_well_{fingerprint.gem_well}"
            g_gw.add_node(gw_node_id, id=gw_node_id, label=fingerprint.gem_well)
            if multi_graph.is_multiplexed():
                assert g_tags is not None
                tag_node_id = f"tag_{i + 1}_{j + 1}"
                label = "+".join(fingerprint.tag_names)
                g_tags.add_node(tag_node_id, id=tag_node_id, label=label)
                g.add_edge(samp_node_id, tag_node_id)
                g.add_edge(tag_node_id, gw_node_id)
            else:
                g.add_edge(samp_node_id, gw_node_id)

    ## The libraries cluster
    g_lib = g.add_subgraph(
        name="cluster_libraries",
        style=GraphVizProps.FILLEDSTYLE,
        color=GraphVizColors.GREY,
        rank=GraphVizProps.RANKSAME,
    )
    g_lib.add_node(
        "Physical Libraries",
        style=GraphVizProps.FILLEDSTYLE,
        color=GraphVizColors.GREY,
        fillcolor=GraphVizColors.GREY,
    )

    ## The fastq_id cluster
    g_fq = g.add_subgraph(
        name="cluster_fastq",
        style=GraphVizProps.FILLEDSTYLE,
        color=GraphVizColors.GREY,
        rank=GraphVizProps.RANKSAME,
    )
    g_fq.add_node(
        "Fastq IDs",
        style=GraphVizProps.FILLEDSTYLE,
        color=GraphVizColors.GREY,
        fillcolor=GraphVizColors.GREY,
    )

    for i, (physical_library_id, gem_well, fastq_ids) in enumerate(
        multi_graph.library_info_for_graphviz()
    ):
        lib_node_id = f"library_{i + 1}"
        g_lib.add_node(lib_node_id, id=lib_node_id, label=physical_library_id)
        gw_node_id = f"gem_well_{gem_well}"
        g.add_edge(gw_node_id, lib_node_id)
        for j, fastq_id in enumerate(fastq_ids):
            fq_node_id = f"fastq_{i + 1}_{j + 1}"
            g_fq.add_node(fq_node_id, id=fq_node_id, label=fastq_id)
            g.add_edge(lib_node_id, fq_node_id)

    return g
