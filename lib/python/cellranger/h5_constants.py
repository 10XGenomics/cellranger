#!/usr/bin/env python
#
# Copyright (c) 2018 10X Genomics, Inc. All rights reserved.
#

MATRIX_MEM_GB_MULTIPLIER = 2.6  # Increased from 2.0 to enable high-diversity samples
VIS_HD_MATRIX_MEM_GB_MULTIPLIER = 1.0
NUM_MATRIX_ENTRIES_PER_MEM_GB = 50e6

# Empirical obs: with new CountMatrix setup, take ~ 50 bytes/bc
NUM_MATRIX_BARCODES_PER_MEM_GB = 10000000
MIN_MEM_GB = 3
MIN_MEM_GB_NOWHITELIST = 64

# Empirical obs of mem usage with anndatas created for running CAS
MEM_BYTES_PER_MATRIX_NNZ_H5AD = 32  # 2^5
MEM_BYTES_PER_MATRIX_BARCODE_H5AD = 512  # 2^9
MEM_BYTES_PER_MATRIX_FEATURE_H5AD = 1024  # 2^10
MEM_BYTES_CONSTANT_H5AD = 8_388_608  # 2^23

GZIP_SUFFIX = b".gz"
LZ4_SUFFIX = b".lz4"

H5_COMPRESSION_LEVEL = 1
H5_FILETYPE_KEY = "filetype"
H5_FEATURE_REF_ATTR = "features"
H5_TARGET_SET_ATTR = b"target_sets"
H5_BCS_ATTR = "barcodes"
H5_MATRIX_DATA_ATTR = "data"
H5_MATRIX_INDICES_ATTR = "indices"
H5_MATRIX_INDPTR_ATTR = "indptr"
H5_MATRIX_SHAPE_ATTR = "shape"
H5_MATRIX_ATTRS = {
    H5_MATRIX_DATA_ATTR: "int32",
    H5_MATRIX_INDICES_ATTR: "int64",
    H5_MATRIX_INDPTR_ATTR: "int64",
    H5_MATRIX_SHAPE_ATTR: "int32",
}
H5_CHEMISTRY_DESC_KEY = "chemistry_description"
H5_LIBRARY_ID_MAPPING_KEY = "library_ids"
H5_ORIG_GEM_GROUP_MAPPING_KEY = "original_gem_groups"
H5_METADATA_ATTRS = [
    H5_LIBRARY_ID_MAPPING_KEY,
    H5_ORIG_GEM_GROUP_MAPPING_KEY,
    H5_CHEMISTRY_DESC_KEY,
]
