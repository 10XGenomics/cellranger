#!/usr/bin/env python
#
# Copyright (c) 2018 10X Genomics, Inc. All rights reserved.
#

from collections import OrderedDict, defaultdict
import cellranger.library_constants as lib_constants

def get_library_type(sample_def):
    return sample_def.get('library_type') or lib_constants.DEFAULT_LIBRARY_TYPE


def get_gem_group(sample_def):
    """ Get the GEM group from a sample def.

    Defaults to 1 if the gem_group is specified as None.

    Args:
      sample_def (dict): Sample def
    Returns:
      int: GEM group
    """
    gg = sample_def['gem_group'] or 1
    return int(gg)


def get_library_id(sample_def):
    return sample_def.get('library_id')


def assign_library_ids(sample_defs):
    """ Assign library ids to a list of sample defs if not given already.

    Each (gem group, library_type) must have a distinct library id. If missing,
    the library id becomes an increasing integer.

    Args:
      sample_defs (list of dict): Sample defs, some possibly already containing a library_id key
    Returns:
      list of str: Library IDs, one per sample def.
    """
    sd_libs = []
    ids = OrderedDict()

    # Map library ids to libraries
    for sd in sample_defs:
        gem_group = get_gem_group(sd)
        lib_type = get_library_type(sd)
        lib = (gem_group, lib_type)
        sd_libs.append(lib)

        # If a library type isn't given, assign it an integer that maps
        # uniquely to (gem_group, library_type)
        default_id = str(len(ids))

        if lib not in ids:
            ids[lib] = []

        lib_id = sd.get('library_id') or default_id

        if ':' in lib_id or '\t' in lib_id:
            raise ValueError('Invalid library ID: "%s". Library IDs may not contain ":" or tab characters.' % lib_id)

        ids[lib].append(lib_id)

    # Verify that there is one library_id per library
    for lib, lib_ids in ids.iteritems():
        if len(lib_ids) > 1:
            raise ValueError('Library (gem_group=%d, type=%s) was labeled with multiple library IDs: %s .' % \
                             (lib[0], lib[1], ', '.join(lib_ids)))

    # Verify that there is one library per library_id
    id2libs = defaultdict(set)
    for lib, lib_ids in ids.iteritems():
        assert len(lib_ids) == 1
        id2libs[lib_ids[0]].add(lib)

    for lib_id, libs in id2libs.iteritems():
        if len(libs) > 1:
            lib_list = ', '.join(map(lambda x: '(gem_group=%s, library_type=%s)' % x, libs))
            raise ValueError('Library ID is not unique: "%s". It maps to the following libraries: %s' % (lib_id, lib_list))
        assert len(libs) == 1

    # Map sample defs to library ids
    return [ids[lib][0] for lib in sd_libs]
