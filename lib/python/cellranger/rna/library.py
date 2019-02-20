#!/usr/bin/env python
#
# Copyright (c) 2018 10X Genomics, Inc. All rights reserved.
#
import json
import re
import cellranger.library_constants as lib_constants
# TODO: Move lib_constants here
CRISPR_LIBRARY_TYPE = 'CRISPR Guide Capture'
ANTIBODY_LIBRARY_TYPE = 'Antibody Capture'
CUSTOM_LIBRARY_TYPE_PREFIX = "Custom"
CRISPR_METRIC_PREFIX = 'CRISPR'
ANTIBODY_METRIC_PREFIX = 'ANTIBODY'
FEATURE_LIBRARY_TYPE_MAP = {
CRISPR_METRIC_PREFIX : CRISPR_LIBRARY_TYPE,
ANTIBODY_METRIC_PREFIX : ANTIBODY_LIBRARY_TYPE,
}
RECOGNIZED_FEATURE_TYPES = [lib_constants.GENE_EXPRESSION_LIBRARY_TYPE,
                            CRISPR_LIBRARY_TYPE,
                            ANTIBODY_LIBRARY_TYPE,
]

def get_library_type_metric_prefix(lib_type):
    """Get the metric prefix for a given library type.
    """
    if lib_type == lib_constants.GENE_EXPRESSION_LIBRARY_TYPE:
        return ''
    elif lib_type == CRISPR_LIBRARY_TYPE:
        return CRISPR_METRIC_PREFIX + '_'
    elif lib_type == ANTIBODY_LIBRARY_TYPE:
        return ANTIBODY_METRIC_PREFIX + '_'
    else:
        return '%s_' % lib_type

def get_bam_library_info(bam):
    """Get the library info from a BAM's comment lines.
    Args:
      bam (pysam.AlignmentFile): BAM file
    Returns:
      list of dicts
    """
    comments = bam.header['CO']
    libraries = []
    for comment in comments:
        m = re.match(r'^library_info:(.+)$', comment)
        if m:
            libraries.append(json.loads(m.group(1)))
    return libraries

def has_genomes(library_type):
    """Do genomes make sense for a library type"""
    return library_type == lib_constants.GENE_EXPRESSION_LIBRARY_TYPE
