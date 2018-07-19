#!/usr/bin/env python
#
# Copyright (c) 2018 10X Genomics, Inc. All rights reserved.
#

from collections import OrderedDict
import cellranger.library_constants as lib_constants

def get_library_type(sample_def):
    return sample_def.get('library_type')

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


def assign_library_ids(sample_defs, default_library_type=lib_constants.DEFAULT_LIBRARY_TYPE):
    """ Assign library ids to a list of sample defs if not given already.

    Each (gem group, library_type) must have a distinct library id. If missing,
    the library id becomes an increasing integer.

    Args:
      sample_defs (list of dict): Sample defs, some possibly already containing a library_id key
      default_library_type (str): Use if a library_type is None or empty
    Returns:
      list of str: Library IDs, one per sample def.
    """
    # Resulting library ID for each sample def
    sd_libs = []

    # Map tuple of (gem_group (int), library_type (str)) to library ID (str)
    library_to_id = OrderedDict()

    # Assign a library ID to each sample def (most commonly by generating a new, unique one)
    for sd in sample_defs:
        gem_group = get_gem_group(sd)
        lib_type = get_library_type(sd) or default_library_type
        lib_tuple = (gem_group, lib_type) # library
        sd_libs.append(lib_tuple) # library for each sample def

        # If a library ID wasn't given by the user, assign it an integer that maps
        # uniquely to (gem_group, library_type)
        default_id = library_to_id.get(lib_tuple, str(len(library_to_id)))

        lib_id = sd.get('library_id') or default_id

        if ':' in lib_id or '\t' in lib_id:
            raise ValueError('Invalid library ID: "%s". Library IDs may not contain ":" or tab characters.' % lib_id)

        # Check if this library ID is already assigned to a different library.
        if lib_id in set(library_to_id.values()):
            # Find the first colliding library
            other_lib = [other_lib for other_lib, other_lib_id in library_to_id.iteritems() if other_lib_id==lib_id][0]
            if other_lib != lib_tuple:
                raise ValueError('Library ID "%s" is already associated with library "%s." A library ID must identify exactly one library.' % (lib_id, other_lib))

        if lib_tuple not in library_to_id:
            library_to_id[lib_tuple] = lib_id

        elif library_to_id[lib_tuple] != lib_id:
            # Library already has an ID
            raise ValueError('Library "%s" already has ID "%s." Cannot assign it the different ID "%s."' %
                             (str(lib_tuple),
                              library_to_id[lib_tuple],
                              lib_id))

    # Map sample defs to library ids
    return [library_to_id[lib] for lib in sd_libs]
