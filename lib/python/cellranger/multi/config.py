#!/usr/bin/env python
#
# Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
#

"""Python API for the CrMultiConfig defined in the rust crate `cr_types`.

Snapshot tests in Rust in conjunction with `test_multi_config.py` ensures that
the API is up to date
"""


from __future__ import annotations

import json
import os
from dataclasses import dataclass
from enum import Enum
from subprocess import PIPE, CalledProcessError, Popen

################################################################################
# NOTE: The schema of json (including the multi graph) that is parsed in the
# classes defined below comes from the rust code in `cr_types`. We want to
# have flexibility in modifying the rust code, hence, there is NO guarantee
# that the schema would not change and it is NOT RECOMMENDED to parse the
# json manually by using the constants below. Instead, use the classes defined
# in this module to read the json. There are tests defined in
# `test_multi_config.py` which ensures that the classes defined in this module
# are updated in sync with any changes in the rust code.

# For parsing LibraryFeatures
GENES_KEY = "GeneExpression"
VDJ_KEY = "Vdj"
FB_KEY = "FeatureBarcodes"

# For parsing Sample
SAMPLE_ID_KEY = "sample_id"
DESC_KEY = "description"
FINGERPRINTS_KEY = "fingerprints"

# For parsing CrMultiGraph
LIBRARIES_KEY = "libraries"
SAMPLES_KEY = "samples"
################################################################################


class CellMultiplexingType(Enum):
    """Cell multiplexing types supported across products."""

    CMO = "CMO"
    RTL = "RTL"
    OH = "OH"

    @staticmethod
    def all():
        """Return all FeatureTypes in the order consistent with the Rust enum definition."""
        return [CellMultiplexingType.CMO, CellMultiplexingType.RTL, CellMultiplexingType.OH]

    @staticmethod
    def from_str(val):
        """Create a CellMultiplexingType object from a string."""
        for data in CellMultiplexingType:
            if data.value == val:
                return data
        raise ValueError(f"Could not parse CellMultiplexingType from {val}")

    @staticmethod
    def from_int(val):
        """Create a FeatureType object from an int."""
        return CellMultiplexingType.all()[val]


class FeatureType(Enum):
    """The feature types which are supported in our matrix."""

    # TODO: re-enable after pylint 2.8 is released.
    # pylint: disable=invalid-name
    Antibody = "Antibody Capture"
    Antigen = "Antigen Capture"
    CRISPR = "CRISPR Guide Capture"
    Multiplexing = "Multiplexing Capture"
    Custom = "Custom"
    Gene = "Gene Expression"

    @staticmethod
    def all():
        """Return all FeatureTypes in the order consistent with the Rust enum definition."""
        return [
            FeatureType.Antibody,
            FeatureType.Antigen,
            FeatureType.CRISPR,
            FeatureType.Multiplexing,
            FeatureType.Custom,
            FeatureType.Gene,
        ]

    @staticmethod
    def from_str(val):
        """Create a FeatureType object from a string."""
        for data in FeatureType:
            if data.value == val:
                return data
        raise ValueError(f"Could not parse feature type from {val}")

    @staticmethod
    def from_int(val):
        """Create a FeatureType object from an int."""
        return FeatureType.all()[val]


def decode_feature_types(byte):
    """Feature types are bit encoded.

    This function decodes them for ease of use.
    """
    feat_types = []
    for i, feat in enumerate(FeatureType.all()):
        if byte & (1 << i):
            feat_types.append(feat)
    return feat_types


class VdjChainType(Enum):
    """VDJ chain type input by the user.

    The default is to auto detect the chain type.
    """

    # TODO: re-enable after pylint 2.8 is released.
    # pylint: disable=invalid-name
    Tcr = "T Cell Receptor"
    Ig = "B Cell Receptor"
    Auto = "Auto"

    @staticmethod
    def from_str(val):
        """Create a VdjChainType object from a string."""
        for data in VdjChainType:
            if data.value == val:
                return data
        raise ValueError(f"Could not parse vdj chain type from {val}")


class LibraryType(Enum):
    """TODO: Should targeting be a different library type?"""

    # TODO: re-enable after pylint 2.8 is released.
    # pylint: disable=invalid-name
    GeneExpression = "Gene Expression"
    ImmuneProfiling = "Immune Profiling"
    FeatureBarcoding = "Feature Barcoding"

    @staticmethod
    def from_str(val):
        """Create a LibraryType object from a string."""
        for data in LibraryType:
            if data.value == val:
                return data
        raise ValueError(f"Could not parse library type from {val}")


class LibraryFeatures:
    """Encapsulate the library type and associated feature types."""

    def __init__(self, gex=None, vdj=None, feature_types=None):
        self.gex = gex
        self.vdj = vdj
        if feature_types is not None:
            self.feature_types = set(feature_types)
        else:
            self.feature_types = None
        self.validate()

    def __eq__(self, other):
        return (
            self.gex == other.gex
            and self.vdj == other.vdj
            and self.feature_types == other.feature_types
        )

    def __hash__(self):
        return hash((self.gex, self.vdj, self.feature_types))

    def validate(self):
        """Ensure that only one of GEX/VDJ/FB is not None."""
        if self.gex is not None:
            assert self.vdj is None
            assert self.feature_types is None
        elif self.vdj is not None:
            assert self.gex is None
            assert self.feature_types is None
        elif self.feature_types is not None:
            assert self.gex is None
            assert self.vdj is None
        else:
            raise ValueError("All features are none")

    def library_type(self):
        """The LibraryType associated with this object."""
        if self.gex is not None:
            return LibraryType.GeneExpression
        elif self.vdj is not None:
            return LibraryType.ImmuneProfiling
        elif self.feature_types is not None:
            return LibraryType.FeatureBarcoding
        raise ValueError("Library type cannot be inferred in the absence of any features")

    def library_features(self):
        """Return the FeatureType if any associated with this object."""
        if self.feature_types is not None:
            return [i.value for i in self.feature_types]
        else:
            return []

    @staticmethod
    def from_json_val(val):
        """Create a LibraryFeatures object from a Json value.

        Value can be a string or a dictionary.
        """
        assert isinstance(val, dict)
        assert len(val) == 1
        if GENES_KEY in val:
            return LibraryFeatures(gex=True)
        elif VDJ_KEY in val:
            return LibraryFeatures(vdj=VdjChainType.from_str(val[VDJ_KEY]))
        elif FB_KEY in val:
            return LibraryFeatures(feature_types=decode_feature_types(val[FB_KEY]))
        raise ValueError(f"Could not parse {val} as LibraryFeatures")


@dataclass(frozen=True)
class LibraryDef:
    """Mirrors types::LibraryDef in lib/rust/cr_types.

    For now, we are not tracking the entire FastqDef, only the fastq IDs.
    """

    physical_library_id: str
    gem_well: int
    library_features: LibraryFeatures
    fastq_ids: list[str]

    @staticmethod
    def from_json_val(val):
        """Create a LibraryDef object from a dictionary."""
        assert isinstance(val, dict)
        physical_library_id = val["physical_library_id"]
        gem_well = val["gem_well"]
        library_features = LibraryFeatures.from_json_val(val["library_features"])
        fastq_ids = [fastq["id"] for fastq in val["fastqs"]]
        return LibraryDef(physical_library_id, gem_well, library_features, fastq_ids)


@dataclass(frozen=True)
class Fingerprint:
    """Refers to a set of tagged or untagged cells that went through a specific gem well."""

    gem_well: int
    tag_name: str | None = None
    translated_tag_names: list[str] | None = None
    cell_multiplexing_type: CellMultiplexingType | None = None

    def __post_init__(self):
        """Validate that we have an expected combination of fields."""
        fields_count = sum(
            [
                self.tag_name is not None,
                self.translated_tag_names is not None,
                self.cell_multiplexing_type is not None,
            ]
        )
        if fields_count not in (0, 3):
            raise ValueError(f"Unexpected combination of fields in fingerprint: {self}")

    @classmethod
    def from_json_val(cls, val: dict):
        """Create a Fingerprint object from a dictionary."""
        assert isinstance(val, dict)
        assert len(val) == 1
        tagged = val.get("Tagged", None)
        if tagged is not None:
            assert isinstance(tagged, dict)

            if len(tagged) != 4:
                raise ValueError(f"Invalid number of fields in tagged fingerprint: {tagged}")
            kwargs = tagged.copy()
            kwargs["cell_multiplexing_type"] = CellMultiplexingType.from_str(
                tagged["cell_multiplexing_type"]
            )
        else:
            untagged = val["Untagged"]
            assert len(untagged) == 1
            kwargs = untagged
        return cls(**kwargs)

    def is_multiplexed(self) -> bool:
        return self.tag_name is not None

    def get_cell_multiplexing_type(self) -> CellMultiplexingType | None:
        return self.cell_multiplexing_type


class Sample:
    """A single sample in a multi run."""

    def __init__(self, sample_id, desc, fingerprints):
        self.sample_id = sample_id
        self.description = desc
        self.fingerprints = fingerprints

    def __eq__(self, other):
        return (
            self.sample_id == other.sample_id
            and self.description == other.description
            and self.fingerprints == other.fingerprints
        )

    def __hash__(self):
        return hash((self.sample_id, self.description, self.fingerprints))

    @staticmethod
    def from_json_val(val):
        """Create a Sample object from a dictionary."""
        assert isinstance(val, dict)
        sample_id = val[SAMPLE_ID_KEY]
        desc = val[DESC_KEY]
        fingerprints = [Fingerprint.from_json_val(fp) for fp in val[FINGERPRINTS_KEY]]
        return Sample(sample_id, desc, fingerprints)

    def is_multiplexed(self):
        fp_multiplexed = [fp.is_multiplexed() for fp in self.fingerprints]
        assert all(x == fp_multiplexed[0] for x in fp_multiplexed)  # Sanity check
        return fp_multiplexed[0]

    def get_cell_multiplexing_type(self) -> CellMultiplexingType:
        """Return the CellMultiplexingType for sample.

        Returns:
            CellMultiplexingType: cell multiplexing type for sample
        """
        fp_cell_multiplexing_type = [fp.get_cell_multiplexing_type() for fp in self.fingerprints]
        assert all(
            x == fp_cell_multiplexing_type[0] for x in fp_cell_multiplexing_type
        )  # Sanity check
        return fp_cell_multiplexing_type[0]


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

    Only the methods we need for the CrMultiGraph SVG.
    """

    def __init__(self, *_, **kwargs):
        kwargs.pop("strict", None)
        kwargs.pop("directed", None)

        self.graph_attr = kwargs
        self.node_attr = {}
        self.subgraphs = {}
        self.nodes = {}
        self.edges = {}

    def add_subgraph(self, *, name="", **kwargs) -> AGraph:
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
        dot_cmd = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../bin/dot"))

        cmd = [dot_cmd, "-Tsvg", "-o", fname]
        with Popen(cmd, stdin=PIPE) as dot:
            dot.communicate(input=g.encode("utf-8"))
            if dot.returncode != 0:
                raise CalledProcessError(dot.returncode, cmd)


class CrMultiGraph:
    """LibraryDefs and Samples in this multi run."""

    def __init__(self, libraries: list[LibraryDef], samples: list[Sample]):
        self.libraries = libraries
        self.samples = samples

    def __eq__(self, other):
        return self.libraries == other.libraries and self.samples == other.samples

    def __hash__(self):
        return hash((self.libraries, self.samples))

    @staticmethod
    def from_json_val(val):
        """Create a CrMultiGraph object from a dictionary."""
        assert isinstance(val, dict)
        libraries = [LibraryDef.from_json_val(v) for v in val[LIBRARIES_KEY]]
        samples = [Sample.from_json_val(v) for v in val[SAMPLES_KEY]]
        return CrMultiGraph(libraries, samples)

    @staticmethod
    def from_json_file(multi_graph_json):
        """Create a CrMultiGraph object from a json file."""
        with open(multi_graph_json) as f:
            data = json.load(f)
        return CrMultiGraph.from_json_val(data)

    def is_multiplexed(self) -> bool:
        samp_multiplexed = [samp.is_multiplexed() for samp in self.samples]
        assert all(x == samp_multiplexed[0] for x in samp_multiplexed)  # Sanity check
        return samp_multiplexed[0]

    def get_cell_multiplexing_type(self) -> CellMultiplexingType:
        """Return the CellMultiplexingType for multigraph.

        Returns:
            CellMultiplexingType: cell multiplexing type for multigraph
        """
        samp_cell_multiplexing_type = [samp.get_cell_multiplexing_type() for samp in self.samples]
        assert all(
            x == samp_cell_multiplexing_type[0] for x in samp_cell_multiplexing_type
        )  # Sanity check
        return samp_cell_multiplexing_type[0]

    def to_graphviz_digraph(self) -> AGraph | None:
        """Convert CrMultiGraph to a graphviz `digraph`.

        This is what shows up in the experimental design tab of the multi websummary.

        TODO: This code makes simplifications which are only true in a single
        gem well case and needs to be updated when multi-gem well is supported
        """
        # pylint: disable=too-many-locals
        is_multiplexed = self.is_multiplexed()
        cell_multiplexing_type = self.get_cell_multiplexing_type()

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
        if is_multiplexed:
            g_tags = g.add_subgraph(
                name="cluster_tags",
                style=GraphVizProps.FILLEDSTYLE,
                color=GraphVizColors.GREY,
                rank=GraphVizProps.RANKSAME,
            )

            if cell_multiplexing_type == CellMultiplexingType.CMO:
                g_tags.add_node(
                    "CMO Tags",
                    style=GraphVizProps.FILLEDSTYLE,
                    color=GraphVizColors.GREY,
                    fillcolor=GraphVizColors.GREY,
                )
            elif cell_multiplexing_type == CellMultiplexingType.RTL:
                g_tags.add_node(
                    "Probe Barcode IDs",
                    style=GraphVizProps.FILLEDSTYLE,
                    color=GraphVizColors.GREY,
                    fillcolor=GraphVizColors.GREY,
                )
            elif cell_multiplexing_type == CellMultiplexingType.OH:
                g_tags.add_node(
                    "Overhang IDs",
                    style=GraphVizProps.FILLEDSTYLE,
                    color=GraphVizColors.GREY,
                    fillcolor=GraphVizColors.GREY,
                )
            else:
                raise ValueError("Unsupported cell multiplexing type: %s" % cell_multiplexing_type)

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
        for i, samp in enumerate(self.samples):
            samp_node_id = f"sample_{i + 1}"
            g_samp.add_node(samp_node_id, id=samp_node_id, label=samp.sample_id)
            for j, fingerprint in enumerate(samp.fingerprints):
                gw_node_id = f"gem_well_{fingerprint.gem_well}"
                g_gw.add_node(gw_node_id, id=gw_node_id, label=fingerprint.gem_well)
                if is_multiplexed:
                    assert g_tags is not None
                    tag_node_id = f"tag_{i + 1}_{j + 1}"
                    label = "+".join([fingerprint.tag_name] + fingerprint.translated_tag_names)
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

        for i, lib in enumerate(self.libraries):
            lib_node_id = f"library_{i + 1}"
            g_lib.add_node(lib_node_id, id=lib_node_id, label=lib.physical_library_id)
            gw_node_id = f"gem_well_{lib.gem_well}"
            g.add_edge(gw_node_id, lib_node_id)
            for j, fastq_id in enumerate(lib.fastq_ids):
                fq_node_id = f"fastq_{i + 1}_{j + 1}"
                g_fq.add_node(fq_node_id, id=fq_node_id, label=fastq_id)
                g.add_edge(lib_node_id, fq_node_id)

        return g

    def render_to_svg(self, fname):
        """Render the view svg.

        Returns:
            bool: whether it was rendered
        """
        g = self.to_graphviz_digraph()
        if g is not None:
            # print(g.string())  # Useful for debugging
            g.draw(fname)
            return True

        return False
