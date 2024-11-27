#!/usr/bin/env python
#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#
# Code for estimating tile-extents on a lane
#
from __future__ import annotations

import re
import xml.etree.ElementTree as etree
from collections.abc import Iterable
from typing import ClassVar, Literal, NamedTuple

import tenkit.stats as tk_stats


class ReadLocation(NamedTuple):
    flowcell: str
    lane: str
    surface: int
    swath: int
    tile: int
    x: int
    y: int


READ_NAME_REGEX: str = (
    r"^[a-zA-Z0-9\-]+:[a-z0-9\-]+:([a-zA-Z0-9\-]+):([0-9]):([0-9])([0-9])([0-9]+):([0-9]+):([0-9]+)$"
)

NUM_READS_TO_ESTIMATE_TILE_EXTENTS: int = 2000


def extract_read_position(read) -> ReadLocation | None:
    """Get read position info from a read's qname."""
    match = re.match(READ_NAME_REGEX, read.qname)
    if match is None:
        return None
    groups = match.groups()
    return ReadLocation(
        flowcell=groups[0],
        lane=groups[1],
        surface=int(groups[2]),
        swath=int(groups[3]),
        tile=int(groups[4]),
        x=int(groups[5]),
        y=int(groups[6]),
    )


def get_flowcell_lane_count(run_info_xml) -> int:
    tree = etree.parse(run_info_xml)
    run_node = tree.getroot().find("Run")
    assert run_node is not None
    flowcell_layout = run_node.find("FlowcellLayout")
    assert flowcell_layout is not None
    return int(flowcell_layout.attrib["LaneCount"])


def get_flowcell_layout(run_info_xml) -> DeclaredFlowcellLayout | None:
    """:param run_info_xml: Path to RunInfo.xml.

    :return: PatternedFlowcellLayout info, or none if the run isn't from a patterned flowcell
    """
    tree = etree.parse(run_info_xml)
    run_node = tree.getroot().find("Run")
    assert run_node is not None
    flowcell_layout = run_node.find("FlowcellLayout")
    if flowcell_layout is not None:
        return DeclaredFlowcellLayout.from_flowcell_layout_node(flowcell_layout)
    else:
        return None


class DeclaredFlowcellLayout:
    """Flowcell layout, as parsed from xml."""

    NAMING_CONVENTION_TILE_LENGTHS: ClassVar[dict[str, int]] = {"FiveDigit": 5, "FourDigit": 4}

    def __init__(
        self,
        *,
        lane_count: int | None = None,
        surface_count: int | None = None,
        swath_count: int | None = None,
        tile_count: int | None = None,
        section_per_lane: int | None = None,
        lane_per_section: int | None = None,
        tile_length: int | None = None,
    ) -> None:
        self.lane_count = lane_count
        self.surface_count = surface_count
        self.swath_count = swath_count
        self.tile_count = tile_count
        self.section_per_lane = section_per_lane
        self.lane_per_section = lane_per_section
        self.tile_length = tile_length

    @classmethod
    def from_flowcell_layout_node(
        cls: type[DeclaredFlowcellLayout], node: etree.Element
    ) -> DeclaredFlowcellLayout:
        """Construct a PatternedFlowcellLayout object from RunInfo.xml.

        Reads the FlowcellLayout node on RunInfo.xml.
        """

        def int_attr_or_none(key: str) -> int | None:
            val = node.attrib.get(key)
            return int(val) if val is not None else None

        tile_length = None

        tile_set = node.find("TileSet")
        if tile_set is not None:
            tile_length = cls.NAMING_CONVENTION_TILE_LENGTHS.get(
                tile_set.attrib.get("TileNamingConvention")
            )

        return cls(
            lane_count=int_attr_or_none("LaneCount"),
            surface_count=int_attr_or_none("SurfaceCount"),
            swath_count=int_attr_or_none("SwathCount"),
            tile_count=int_attr_or_none("TileCount"),
            section_per_lane=int_attr_or_none("SectionPerLane"),
            lane_per_section=int_attr_or_none("LanePerSection"),
            tile_length=tile_length,
        )


class LaneLayout:
    """Extents of a particular flowcell, lane combination."""

    def __init__(self) -> None:
        self.tile_width: int = 0
        self.tile_height: int = 0
        self.num_swaths: int = 0  # vertical stripes of tiles
        self.tiles_per_swath: int = 0  # number of tiles in a swath
        self.num_reads_observed: int = 0  # number of reads used to estimate extents

    def estimate_extents(self, read_loc: ReadLocation) -> None:
        """Cumulatively estimate extents from read observations."""
        self.tile_width = max(self.tile_width, read_loc.x)
        self.tile_height = max(self.tile_height, read_loc.y)
        self.num_swaths = max(self.num_swaths, read_loc.swath)
        self.tiles_per_swath = max(self.tiles_per_swath, read_loc.tile)
        self.num_reads_observed += 1

    def has_diffusion_duplicates(self, diffusion_radius: float) -> bool:
        # Only screen lanes where the diffusion area is <2% of the flowcell
        fc_area = self.tile_width * self.tile_height * self.num_swaths * self.tiles_per_swath
        diffusion_area = 3.141519 * diffusion_radius**2

        return tk_stats.robust_divide(diffusion_area, fc_area) < 0.025

    def to_dict(self) -> dict[str, int]:
        """Serialize to a dict."""
        return self.__dict__

    @staticmethod
    def from_dict(d: dict[str, int]):
        """Read from a dict."""
        result = LaneLayout()
        for key, value in d.items():
            setattr(result, key, int(value))
        return result


class LaneCoordinateSystem:
    """Coordinate system for a collection of lanes."""

    def __init__(self) -> None:
        self.lanes: dict[str, LaneLayout] = {}  # maps flowcell_lane to a LaneLayout

    def estimate_tile_extents(self, bam: Iterable) -> None:
        """Estimate tile extents from a sampling of within-tile x,y coords."""
        for read in bam:
            read_loc = extract_read_position(read)

            # extreme pathology case but possible
            if not read_loc:
                continue

            fc_lane = LaneCoordinateSystem.fc_lane_key(read_loc.flowcell, read_loc.lane)
            if fc_lane not in self.lanes:
                self.lanes[fc_lane] = LaneLayout()

            self.lanes[fc_lane].estimate_extents(read_loc)

            if all(
                x.num_reads_observed >= NUM_READS_TO_ESTIMATE_TILE_EXTENTS
                for x in self.lanes.values()
            ):
                break

    @staticmethod
    def fc_lane_key(flowcell: str, lane: str) -> str:
        """Build a key from a flowcell and lane."""
        return f"{flowcell}_{lane}"

    def get_layout_for_read_loc(self, read_loc: ReadLocation) -> LaneLayout:
        fc_lane = LaneCoordinateSystem.fc_lane_key(read_loc.flowcell, read_loc.lane)
        # If we encounter a lane we've never seen before, assume it's like
        # one we've already observed
        if fc_lane not in self.lanes:
            fc_lane = next(iter(self.lanes))
        layout = self.lanes[fc_lane]

        return layout

    def convert_to_lane_coords(
        self, read_loc: ReadLocation
    ) -> tuple[None, None, None] | tuple[int, int, Literal[0]]:
        """Convert within-tile coords to lane coords."""
        # If we know nothing, give up
        if not self.lanes:
            return (None, None, None)

        fc_lane = LaneCoordinateSystem.fc_lane_key(read_loc.flowcell, read_loc.lane)
        # If we encounter a lane we've never seen before, assume it's like
        # one we've already observed
        if fc_lane not in self.lanes:
            fc_lane = next(iter(self.lanes))
        layout = self.lanes[fc_lane]

        lane_x = (read_loc.swath - 1) * layout.tile_width + read_loc.x
        lane_y = (read_loc.tile - 1) * layout.tile_height + read_loc.y
        lane_z = 0  # Ignore surface for now
        return (lane_x, lane_y, lane_z)

    def to_dict(self) -> dict[str, dict[str, int]]:
        """Serialize to a dict."""
        return {key: layout.to_dict() for key, layout in self.lanes.items()}

    @staticmethod
    def from_dict(d: dict[str, dict[str, int]]):
        """Read from a dict."""
        result = LaneCoordinateSystem()
        for layout_key, layout_dict in d.items():
            result.lanes[layout_key] = LaneLayout.from_dict(layout_dict)
        return result
