#!/usr/bin/env python
#
# Copyright (c) 2018 10X Genomics, Inc. All rights reserved.
#
import json
import re
import cellranger.constants as cr_constants

# TODO: Move lib_constants here
CRISPR_LIBRARY_TYPE = 'CRISPR Guide Capture'
ANTIBODY_LIBRARY_TYPE = 'Antibody Capture'
FEATURE_LIBRARY_TYPE_MAP = {
'CRISPR' : CRISPR_LIBRARY_TYPE,
'ANTIBODY' : ANTIBODY_LIBRARY_TYPE,
}

def get_library_type_metric_prefix(lib_type):
    """Get the metric prefix for a given library type.
       Some are hardcoded for historical reasons.
    """
    if lib_type == 'Gene Expression':
        return ''
    elif lib_type == 'CRISPR Guide Capture':
        return 'CRISPR_'
    elif lib_type == 'Antibody Capture':
        return 'ANTIBODY_'
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
    return library_type == cr_constants.GENE_EXPRESSION_LIBRARY_TYPE
