#!/usr/bin/env python
#
# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
#

"""Constants related to BAM file generation."""

# Bam tags
# NUM_HITS_TAG = "NH"
MULTIMAPPER_TAG = "mm"
ANTISENSE_TAG = "AN"
RAW_BARCODE_TAG = "CR"
PROCESSED_BARCODE_TAG = "CB"
RAW_BARCODE_QUAL_TAG = "CY"
RAW_UMI_TAG = "UR"
PROCESSED_UMI_TAG = "UB"
UMI_QUAL_TAG = "UY"
TRANSCRIPTS_TAG = "TX"
GENE_IDS_TAG = "GX"
# GENE_NAMES_TAG = "GN"
MAPPING_REGION_TAG = "RE"
RAW_FEATURE_BARCODE_TAG = "fr"
FEATURE_IDS_TAG = "fx"
EXTRA_FLAGS_TAG = "xf"
# MATE_RESCUE_TAG = "MR"

BAM_CHUNK_SIZE_GB = 0.5
MAX_BAM_CHUNKS = 256

# matches variants of enum cr_types::types::AlignerParam
STAR_ALIGNER = "star"
MINIMAP2_ALIGNER = "minimap2"

# matches values of cr_types::types::AlignerParam.high_conf_mapq()
STAR_DEFAULT_HIGH_CONF_MAPQ = 255
MINIMAP_DEFAULT_HIGH_CONF_MAPQ = 60
