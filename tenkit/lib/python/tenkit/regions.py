#!/usr/bin/env python
#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#
# Useful data structures
#

import bisect
import numpy as np
from collections import namedtuple

# Set up basic region tuple
Region = namedtuple('Region', 'start end')
NamedRegion = namedtuple('NamedRegion', 'start end name')

class Dirs(object):
    left = -1
    right = 1

    @staticmethod
    def from_str(s):
        if s == '+':
            return Dirs.right
        elif s == '-':
            return Dirs.left
        return None


class Regions(object):
    """Defines a set of regions and allows for determining whether a point
    is in ths sets or another region overlaps those sets.
    """
    def __init__(self, regions=None, merge = True):
        if regions is not None and not(regions == []):
            sorted_regions = sorted(regions)
            self.starts = [r[0] for r in sorted_regions]
            self.ends = [r[1] for r in sorted_regions]
            if merge:
                self.make_regions_non_overlapping()
        else:
            self.starts = []
            self.ends = []
        self.current_iter = 0

    def __iter__(self):
        return self

    def next(self):
        if self.current_iter >= len(self.starts):
            self.current_iter = 0
            raise StopIteration
        else:
            self.current_iter += 1
            return Region(self.starts[self.current_iter-1], self.ends[self.current_iter - 1])

    def get_region_list(self):
        """ Returns list of non-overlapping regions
        """

        region_list = []
        for start, end in zip(self.starts, self.ends):
            region_list.append(Region(start, end))
        return region_list

    def add_region(self, reg):
        """ Adds single region to set while removing overlap """

        # in case not given as named tuple
        region = Region(reg[0], reg[1])

        # find where this region should go
        start_index = bisect.bisect_left(self.starts, region.start)
        end_index = bisect.bisect_left(self.ends, region.end)

        # merge if required
        if start_index > 0 and self.ends[start_index-1] > region.start:
            region = Region(self.starts[start_index-1], region.end)
            start_index -= 1
        if end_index < len(self.starts) and self.starts[end_index] < region.end:
            region = Region(region.start, self.ends[end_index])
            end_index += 1

        # if no merge, insert this region in the lists
        if start_index == end_index:
            self.starts = self.starts[:start_index] + [region.start] + self.starts[start_index:]
            self.ends = self.ends[:start_index] + [region.end] + self.ends[start_index:]
        else:
            self.starts = self.starts[:start_index] + [region.start] + self.starts[end_index:]
            self.ends = self.ends[:start_index] + [region.end] + self.ends[end_index:]

    def make_regions_non_overlapping(self):
        """ Merges overlapping regions
        """
        new_starts = [self.starts[0]]
        new_ends = [self.ends[0]]
        last_index = 0
        for n in xrange(1, len(self.starts)):
            this_start = self.starts[n]
            this_end = self.ends[n]
            last_end = new_ends[last_index]

            if this_start >= last_end: # JI:  Should be >= for gap-numbered
                new_starts.append(this_start)
                new_ends.append(this_end)
                last_index += 1
            else:
                new_ends[last_index] = max(this_end, last_end)

        self.starts = new_starts
        self.ends = new_ends

    def intersect(self, regions):
        intersected_regions = Regions()
        for region in self.get_region_list():
            (start, end) = region
            overlapping_regions = regions.overlapping_regions(start, end)
            mini_intersect = [(max(start, overlapping_start), min(end, overlapping_end)) for (overlapping_start, overlapping_end) in overlapping_regions]
            for mini_intersected_region in mini_intersect:
                intersected_regions.add_region(mini_intersected_region)
        return intersected_regions

    def get_total_size(self):
        return sum([self.ends[n] - self.starts[n] for n in xrange(len(self.starts))])

    def contains_point(self, pt):
        """ Determines whether a point is contained in one of the regions
        """
        check_index = bisect.bisect(self.starts, pt) - 1
        if check_index == -1:
            return False

        if pt >= self.starts[check_index] and pt < self.ends[check_index]: # JI: Should be < for gap-numbered
            return True
        return False

    def get_region_containing_point(self, pt):
        """ Determines which (if any) region contains a point
        """
        check_index = bisect.bisect(self.starts, pt) - 1
        if check_index == -1:
            return None

        if pt >= self.starts[check_index] and pt < self.ends[check_index]: # JI: Should be < for gap-numbered
            return Region(self.starts[check_index], self.ends[check_index])
        return None

    def get_closest_region(self, pt):
        """ Returns the start and end of the closest region, and whether the region contains the point
        """
        containining_region = self.get_region_containing_point(pt)
        if not(containining_region is None):
            (start, end) = containining_region
            return (start, end, True)

        right_index = bisect.bisect(self.starts, pt)
        left_index = right_index - 1

        if right_index == 0:
            return (self.starts[right_index], self.ends[right_index], False)
        elif right_index == len(self.starts):
            return (self.starts[left_index], self.ends[left_index], False)
        else:
            left_dist = pt - self.ends[left_index]
            right_dist = self.starts[right_index] - pt

            if left_dist < right_dist:
                return (self.starts[left_index], self.ends[left_index], False)
            else:
                return (self.starts[right_index], self.ends[right_index], False)


    def get_closest_region_to_region(self, start, stop, direction = None):
        """Closest region to a given region and distance.
        direction: 1 means try to get the closest upstream region, -1 means
        closest downstream, None means in any direction. If there are
        overlapping regions, the direction is ignored.

        Return value:
        (start, stop, dist) where (start, stop) is the closest region (in the
        requested direction) and dist is the distance.
        Distance is set to 0 if the region is overlapping.
        If there are no regions (or no regions in the requested direction),
        then the return value is (None, None, np.inf).
        """
        if len(self.starts) == 0:
            return (None, None, np.inf)

        ov_regions = self.overlapping_regions(start, stop)
        if len(ov_regions) > 0:
            # Return the first overlapping region
            (start, stop) = ov_regions[0]
            return (start, stop, 0)

        left_start, left_stop, _ = self.get_closest_region(start)
        left_dist = start - left_stop
        right_start, right_stop, _ = self.get_closest_region(stop)
        right_dist = right_start - stop

        if direction == Dirs.right:
            if right_dist >= 0:
                # closest region (ignoring direction) happens to be
                # in the right direction
                return (right_start, right_stop, right_dist)
            else:
                # closest region is the one with the closest start
                idx = bisect.bisect_left(self.starts, stop)
                if idx == len(self.starts):
                    # all starts are smaller than stop
                    return (None, None, np.inf)
                return (self.starts[idx], self.ends[idx], self.starts[idx] - stop)
        elif direction == Dirs.left:
            if left_dist >= 0:
                return (left_start, left_stop, left_dist)
            else:
                # closest regions is the one with the closest end
                if start < self.starts[0]:
                    return (None, None, np.inf)
                new_regions = sorted(self.get_region_list(), key = lambda x: x.end)
                new_starts = [n.start for n in new_regions]
                new_ends = [n.end for n in new_regions]
                idx = bisect.bisect_left(new_ends, start) - 1
                return (new_starts[idx], new_ends[idx], start - new_ends[idx])
        else:
            # If this is negative, the closest region is to the right of stop
            # (so it has to be stop's closest region)
            if left_dist < 0:
                return (left_start, left_stop, left_start - stop)
            if right_dist > 0 and (left_dist > right_dist):
                return (right_start, right_stop, right_dist)
            return (left_start, left_stop, left_dist)


    def overlaps_region(self, start, end):
        """ Determines whether a region described by start and end overlaps
        any of the regions.
        """
        check_index = bisect.bisect_left(self.ends, start)
        if check_index == len(self.starts):
            return False
        if end > self.starts[check_index]:
            return True
        return False

    def overlapping_regions(self, start, end):
        """ Return regions overlapping the given interval """
        idx_left = bisect.bisect(self.starts, start)
        idx_right = bisect.bisect(self.starts, end)

        if idx_left == -1:
            return []

        if idx_left > 0 and self.ends[idx_left - 1] > start:
            idx_left -= 1

        #if idx_right > 0  and self.ends[idx_right] > start:
        #    idx_right += 1

        return [(self.starts[i], self.ends[i]) for i in range(idx_left, idx_right)]


    def merge(self, other_regions):
        """Takes the other regions and merges
        """
        for other_start, other_end in zip(other_regions.starts, other_regions.ends):
            self.add_region((other_start, other_end))


class NamedRegions(Regions):
    """Defines a set of regions and allows for determining whether a point
    is in ths sets or another region overlaps those sets.
    """
    def __init__(self, regions = None, add_names = False):
        if regions is not None and not(regions == []):
            if add_names:
                regions = [(r[0], r[1], r[2], str(i)) for i, r in enumerate(regions)]
            sorted_regions = sorted(regions)
            self.starts = [r[0] for r in sorted_regions]
            self.ends = [r[1] for r in sorted_regions]
            self.names = [r[2] for r in sorted_regions]

            end_sorted_regions = sorted(regions, key = lambda x:(x[1], x[0], x[2]))
            self.e_starts = [r[0] for r in end_sorted_regions]
            self.e_ends = [r[1] for r in end_sorted_regions]
            self.e_names = [r[2] for r in end_sorted_regions]
        else:
            self.starts = []
            self.ends = []
            self.names = []
            self.e_starts = []
            self.e_ends = []
            self.e_names = []
        self.current_iter = 0


    def overlapping_region_names(self, start, end):
        """Returns a set of names of regions overlapping the given interval."""
        if len(self.names) == 0 or start > self.e_ends[-1] or end < self.starts[0]:
            return set([])

        # e_ends[idx_left:] >= start
        idx_left = bisect.bisect_left(self.e_ends, start)
        # starts[idx_right:] > end
        idx_right = bisect.bisect(self.starts, end)

        names1 = set([self.names[i] for i in range(0, idx_right)])
        names2 = set([self.e_names[i] for i in range(idx_left, len(self.names))])

        return names1.intersection(names2)


    def __iter__(self):
        return self


    def next(self):
        if self.current_iter >= len(self.starts):
            self.current_iter = 0
            raise StopIteration
        else:
            self.current_iter += 1
            return NamedRegion(self.starts[self.current_iter - 1], self.ends[self.current_iter - 1],
                               self.names[self.current_iter - 1])
