#!/usr/bin/env python
#
# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
#

"""Commonly-used constants which have not yet been organized."""

######################################################
# DO NOT add new items to this file.
#
# - If a constant is only used from a single module, put it in that module.
# - If a constant is only used in association with a particular module, put it
#   in that module.
# - If a constant is used in a bunch of places, create a module with a more
#   specifically descriptive name to put those constants.
######################################################
# from __future__ import annotations

import os
from typing import NamedTuple

BARCODE_WHITELIST_PATH = os.path.join(os.path.dirname(__file__), "barcodes")
BARCODE_WHITELIST_TRANSLATE_PATH = os.path.join(BARCODE_WHITELIST_PATH, "translation")

# Product types defs
SPATIAL_PRODUCT_TYPE = "sp"
SINGLE_CELL_PRODUCT_TYPE = "sc"

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


# Cell Ranger read types:
#     all: all reads
#     mapped: defined per region
#     conf_mapped: defined per region
#     conf_mapped_barcoded: read must be conf_mapped and have valid BC if BCs are used
#     conf_mapped_deduped: read must be conf_mapped_barcoded, non-duplicate and have valid UMI if UMIs are used
# ALL_READ_TYPE = "all"
# MAPPED_READ_TYPE = "mapped"
CONF_MAPPED_READ_TYPE = "conf_mapped"
CONF_MAPPED_BC_READ_TYPE = "conf_mapped_barcoded"
CONF_MAPPED_DEDUPED_READ_TYPE = "conf_mapped_deduped_barcoded"

# Molecule types:
#     insert: a distinct sequenced read.
#     fragment: a called siPCR fragment molecule; a distinct (BC, UMI, gene, pos, strand) tuple
#     cDNA: a called cDNA molecule; a distinct (BC, UMI, gene) tuple
#     cDNA_candidate: a candidate cDNA molecule; a distinct (BC, UMI) tuple
INSERT_MOLECULE_TYPE = "insert"
FRAGMENT_MOLECULE_TYPE = "fragment"
CDNA_MOLECULE_TYPE = "cdna"
# CDNA_MOLECULE_CANDIDATE_TYPE = "cdna_candidate"

# Barcode types:
#    all_bcs: all barcodes
#    filtered_bcs: cell-containing barcodes
ALL_BARCODES = "all_bcs"
FILTERED_BARCODES = "filtered_bcs"

# Regions:
#     transcriptome:
#       mapped: read must be mapped to sense strand of at least one transcript
#       conf_mapped: read must be confidently mapped along the sense strand to a unique gene
#     genome:
#       mapped: read must be mapped to at least one chromosome
#       conf_mapped: read must be confidently mapped to at least one chromosome
#     exonic:
#       mapped: read must be mapped to at least one exonic region
#       conf_mapped: read must be confidently mapped to at least one exonic region
#     intergenic:
#       mapped: read must be mapped to at least one intergenic region
#       conf_mapped: read must be confidently mapped to at least one intergenic region
#     intronic:
#       mapped: read must be mapped to at least one intronic region
#       conf_mapped: read must be confidently mapped to at least one intronic region
TRANSCRIPTOME_REGION = "transcriptome"
# GENOME_REGION = "genome"
EXONIC_REGION = "exonic"
INTERGENIC_REGION = "intergenic"
INTRONIC_REGION = "intronic"

# Duplicate types:
#     cdna_pcr_uncorrected: determined by gene ID, barcode, strand, uncorrected umi
#     cdna_pcr: determined by gene ID, barcode, strand, corrected umi
#     si_pcr: determined by gene ID, pos, barcode, strand, corrected umi
# CDNA_PCR_UNCORRECTED_DUPE_TYPE = "cdna_pcr_uncorrected"
CDNA_PCR_DUPE_TYPE = "cdna_pcr"
# SI_PCR_DUPE_TYPE = "si_pcr"

ON_TARGET_SUBSAMPLE = "ontarget"
OFF_TARGET_SUBSAMPLE = "offtarget"

FORWARD_STRAND = b"+"
REVERSE_STRAND = b"-"
STRANDS = [FORWARD_STRAND, REVERSE_STRAND]

THREE_PRIME = "three_prime"
FIVE_PRIME = "five_prime"

INSERT_SIZE_CUTOFFS = list(range(0, 1550, 50))
MIN_COUNTS_PER_GENE = 1
DEFAULT_RECOVERED_CELLS_PER_GEM_GROUP = 3000

TOP_N = 5
HOMOPOLYMER_LENGTH = 15

H5_BC_SEQUENCE_COL = "bc_sequence"

CELLRANGER_VERSION_KEY = "cellranger_version"

REFERENCE_METADATA_FILE = "reference.json"
REFERENCE_STAR_PATH = "star"
REFERENCE_FASTA_PATH = "fasta/genome.fa"
REFERENCE_GENES_GTF_PATH = "genes/genes.gtf"
REFERENCE_GENOMES_KEY = "genomes"
REFERENCE_MEM_GB_KEY = "mem_gb"
# REFERENCE_NUM_THREADS_KEY = "threads"

# Ref metadata keys used by GEX and VDJ
REFERENCE_FASTA_HASH_KEY = "fasta_hash"
REFERENCE_GTF_HASH_KEY = "gtf_hash"
REFERENCE_INPUT_FASTA_KEY = "input_fasta_files"
REFERENCE_INPUT_GTF_KEY = "input_gtf_files"
REFERENCE_MKREF_VERSION_KEY = "mkref_version"
REFERENCE_VERSION_KEY = "version"
REFERENCE_TYPE_KEY = "type"
REFERENCE_TYPE = "Transcriptome"
REFERENCE_METRIC_PREFIX = "reference_"

BAM_CHUNK_SIZE_GB = 0.5
MAX_BAM_CHUNKS = 256

COUNT_GENES_MAX_MEM_GB = 64
NUM_MOLECULE_INFO_ENTRIES_PER_CHUNK = 40000000

# Ensure this value matches cr_types::rna_read::HIGH_CONF_MAPQ.
STAR_DEFAULT_HIGH_CONF_MAPQ = 255

NORM_MODE_NONE = "none"

AGG_ID_FIELD = "library_id"
AGG_H5_FIELD = "molecule_h5"
AGG_BATCH_FIELD = "batch"
AGG_DESCRIPTION_FIELD = "description"

# Cloupe field
AGG_CLOUPE_FIELD = "cloupe_file"


# Namedtuples definitions
class aggr_files(NamedTuple):
    paths: list[str]
    required: bool
    default_location: str


# Definition of the sc aggr files
SC_AGGR_FILES = {
    AGG_H5_FIELD: aggr_files(paths=["molecule_info.h5"], required=True, default_location="")
}

MAX_INSERT_SIZE = 1000

# PIPELINE NAMES
PIPELINE_VDJ = "vdj"
PIPELINE_COUNT = "count"

# Barcode for non-whitelisted barcodes in the per barcode metrics CSV
NO_BARCODE = b"NO_BARCODE"

# Spatial, tissue color
TISSUE_COLOR = "#E84B50"  # a color as understood by matplotlib, see https://matplotlib.org/2.0.2/api/colors_api.html

# alpha blending to apply to spots under tissue; used in both, the automatic and the manual alignment paths.
TISSUE_SPOTS_ALPHA = 0.50

# alpha blending to apply to tissue area outside of spots; used only in the automatic alignment.
# Has to be smaller than TISSUE_SPOTS_ALPHA.
TISSUE_NOSPOTS_ALPHA = 0.25

# color for the bounding box drawn around the tissue detection area
# this is an 8-bit RGB triple
TISSUE_BBOX_COLOR = (29, 67, 122)

# the color for the ring traced around fiducial spots on the QC image
FIDUCIAL_SPOT_COLOR = "#1D437A"
FILTER_LIST = [
    "None",
    "Non-Targeting",
    "Ignore",
]  # List of targets that are considered filter-able ie
