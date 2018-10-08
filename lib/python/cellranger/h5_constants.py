#!/usr/bin/env python
#
# Copyright (c) 2018 10X Genomics, Inc. All rights reserved.
#

MATRIX_MEM_GB_MULTIPLIER = 2  # TODO reduce this once we're confident about the actual memory bounds
NUM_MATRIX_ENTRIES_PER_MEM_GB = 50e6

# Empirical obs: with new CountMatrix setup, take ~ 50 bytes/bc
NUM_MATRIX_BARCODES_PER_MEM_GB = 2000000 
MIN_MEM_GB = 3
MIN_MEM_GB_NOWHITELIST = 64

GZIP_SUFFIX = '.gz'
LZ4_SUFFIX = '.lz4'

H5_COMPRESSION_LEVEL = 1
H5_FILETYPE_KEY = 'filetype'
H5_FEATURE_REF_ATTR = 'features'
H5_BCS_ATTR = 'barcodes'
H5_MATRIX_DATA_ATTR = 'data'
H5_MATRIX_INDICES_ATTR = 'indices'
H5_MATRIX_INDPTR_ATTR = 'indptr'
H5_MATRIX_SHAPE_ATTR = 'shape'
H5_MATRIX_ATTRS = {H5_MATRIX_DATA_ATTR: 'int32', H5_MATRIX_INDICES_ATTR: 'int64', H5_MATRIX_INDPTR_ATTR: 'int64', H5_MATRIX_SHAPE_ATTR: 'int32'}
H5_CHEMISTRY_DESC_KEY = 'chemistry_description'
H5_LIBRARY_ID_MAPPING_KEY = 'library_ids'
H5_ORIG_GEM_GROUP_MAPPING_KEY = 'original_gem_groups'
H5_METADATA_ATTRS = [H5_LIBRARY_ID_MAPPING_KEY, H5_ORIG_GEM_GROUP_MAPPING_KEY, H5_CHEMISTRY_DESC_KEY]
