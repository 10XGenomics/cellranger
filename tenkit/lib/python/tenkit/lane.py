#!/usr/bin/env python
#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#
# Code for estimating tile-extents on a lane
#
import tenkit.stats as tk_stats
from collections import namedtuple
import re
import xml.etree.ElementTree as etree

ReadLocation = namedtuple('ReadLocation', 'flowcell, lane, surface, swath, tile, x, y')
READ_NAME_REGEX="^[a-zA-Z0-9\-]+:[a-z0-9\-]+:([a-zA-Z0-9\-]+):([0-9]):([0-9])([0-9])([0-9]+):([0-9]+):([0-9]+)$"

NUM_READS_TO_ESTIMATE_TILE_EXTENTS = 2000

def extract_read_position(read):
    ''' Get read position info from a read's qname '''
    match = re.match(READ_NAME_REGEX, read.qname)
    if match is None:
        return None
    groups = match.groups()
    return ReadLocation(flowcell=groups[0], lane=groups[1],
                        surface=int(groups[2]), swath=int(groups[3]),
                        tile=int(groups[4]),
                        x=int(groups[5]), y=int(groups[6]))


def get_flowcell_lane_count(run_info_xml):
    tree = etree.parse(run_info_xml)
    run_node = tree.getroot().find("Run")
    flowcell_layout = run_node.find("FlowcellLayout")
    return int(flowcell_layout.attrib["LaneCount"])

def get_flowcell_layout(run_info_xml):
    """
    :param run_info_xml: Path to RunInfo.xml
    :return: PatternedFlowcellLayout info, or none if the run isn't from a patterned flowcell
    """
    tree = etree.parse(run_info_xml)
    run_node = tree.getroot().find("Run")
    flowcell_layout = run_node.find("FlowcellLayout")
    if flowcell_layout is not None:
        return DeclaredFlowcellLayout.from_flowcell_layout_node(flowcell_layout)
    else:
        return None


class DeclaredFlowcellLayout(object):
    NAMING_CONVENTION_TILE_LENGTHS = {
        "FiveDigit": 5,
        "FourDigit": 4
    }
    def __init__(self, *args, **kwargs):
        self.lane_count = kwargs.get('lane_count', None)
        self.surface_count = kwargs.get('surface_count', None)
        self.swath_count = kwargs.get('swath_count', None)
        self.tile_count = kwargs.get('tile_count', None)
        self.section_per_lane = kwargs.get('section_per_lane', None)
        self.lane_per_section = kwargs.get('lane_per_section', None)
        self.tile_length = kwargs.get('tile_length', None)

    @classmethod
    def from_flowcell_layout_node(cls, node):
        """
        Construct a PatternedFlowcellLayout object from the FlowcellLayout
        node on RunInfo.xml.
        """
        intAttrOrNone = lambda key: int(node.attrib.get(key)) if node.attrib.get(key) is not None else None
        tile_length = None

        tile_set = node.find("TileSet")
        if tile_set is not None:
            tile_length = cls.NAMING_CONVENTION_TILE_LENGTHS.get(tile_set.attrib.get("TileNamingConvention"))

        return cls(
            lane_count=intAttrOrNone('LaneCount'),
            surface_count=intAttrOrNone('SurfaceCount'),
            swath_count=intAttrOrNone('SwathCount'),
            tile_count=intAttrOrNone('TileCount'),
            section_per_lane=intAttrOrNone('SectionPerLane'),
            lane_per_section = intAttrOrNone('LanePerSection'),
            tile_length=tile_length
        )

class LaneLayout:
    ''' Extents of a particular flowcell, lane combination '''
    def __init__(self):
        self.tile_width = 0
        self.tile_height = 0
        self.num_swaths = 0          # vertical stripes of tiles
        self.tiles_per_swath = 0     # number of tiles in a swath
        self.num_reads_observed = 0  # number of reads used to estimate extents

    def estimate_extents(self, read_loc):
        ''' Cumulatively estimate extents from read observations '''
        self.tile_width = max(self.tile_width, read_loc.x)
        self.tile_height = max(self.tile_height, read_loc.y)
        self.num_swaths = max(self.num_swaths, read_loc.swath)
        self.tiles_per_swath = max(self.tiles_per_swath, read_loc.tile)
        self.num_reads_observed += 1

    def has_diffusion_duplicates(self, diffusion_radius):
        # Only screen lanes where the diffusion area is <2% of the flowcell
        fc_area = self.tile_width * self.tile_height * self.num_swaths * self.tiles_per_swath
        diffusion_area = 3.141519 * diffusion_radius**2

        return tk_stats.robust_divide(diffusion_area, fc_area) < 0.025

    def to_dict(self):
        ''' Serialize to a dict '''
        return self.__dict__

    @staticmethod
    def from_dict(d):
        ''' Read from a dict '''
        result = LaneLayout()
        for key, value in d.iteritems():
            setattr(result, key, int(value))
        return result

class LaneCoordinateSystem:
    def __init__(self):
        self.lanes = {}        # maps flowcell_lane to a LaneLayout

    def estimate_tile_extents(self, bam):
        ''' Estimate tile extents from a sampling of within-tile x,y coords '''
        for read in bam:
            read_loc = extract_read_position(read)

            # extreme pathology case but possible
            if not read_loc:
                continue

            fc_lane = LaneCoordinateSystem.fc_lane_key(read_loc.flowcell, read_loc.lane)
            if fc_lane not in self.lanes.keys():
                self.lanes[fc_lane] = LaneLayout()

            self.lanes[fc_lane].estimate_extents(read_loc)

            if all([x.num_reads_observed >= NUM_READS_TO_ESTIMATE_TILE_EXTENTS for x in self.lanes.values()]):
                break

    @staticmethod
    def fc_lane_key(flowcell, lane):
        ''' Build a key from a flowcell and lane '''
        return "%s_%s" % (flowcell, lane)


    def get_layout_for_read_loc(self, read_loc):
        fc_lane = LaneCoordinateSystem.fc_lane_key(read_loc.flowcell, read_loc.lane)
        # If we encounter a lane we've never seen before, assume it's like
        # one we've already observed
        if fc_lane not in self.lanes.keys():
            fc_lane = self.lanes.keys()[0]
        layout = self.lanes[fc_lane]

        return layout


    def convert_to_lane_coords(self, read_loc):
        ''' Convert within-tile coords to lane coords '''
        # If we know nothing, give up
        if not self.lanes:
            return (None, None, None)

        fc_lane = LaneCoordinateSystem.fc_lane_key(read_loc.flowcell, read_loc.lane)
        # If we encounter a lane we've never seen before, assume it's like
        # one we've already observed
        if fc_lane not in self.lanes.keys():
            fc_lane = self.lanes.keys()[0]
        layout = self.lanes[fc_lane]

        lane_x = (read_loc.swath - 1) * layout.tile_width + read_loc.x
        lane_y = (read_loc.tile - 1) * layout.tile_height + read_loc.y
        lane_z = 0  # Ignore surface for now
        return (lane_x, lane_y, lane_z)

    def to_dict(self):
        ''' Serialize to a dict '''
        result = {}
        for layout_key, layout in self.lanes.iteritems():
            result[layout_key] = layout.to_dict()
        return result

    @staticmethod
    def from_dict(d):
        ''' Read from a dict '''
        result = LaneCoordinateSystem()
        for layout_key, layout_dict in d.iteritems():
            result.lanes[layout_key] = LaneLayout.from_dict(layout_dict)
        return result

