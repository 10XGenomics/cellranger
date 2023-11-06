# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.

"""Values used in many other modules."""

from __future__ import annotations

######################################################
# DO NOT add new items to this file.
#
# - If a constant is only used from a single module, put it in that module.
# - If a constant is only used in association with a particular module, put it
#   in that module.
# - If a constant is used in a bunch of places, create a module with a more
#   specifically descriptive name to put those constants.
######################################################

HIGH_CONF_MAPQ: int = 60
"""What is considered a high confidence mapped read pair."""

MIN_MATE_OFFSET_DUP_FILTER: int = 20
"""Distance to mate to ensure the"""

# Sequencing settings
ILLUMINA_QUAL_OFFSET: int = ord(b"!")
DEFAULT_HIGH_MAPQ: int = ord(b"<")

# Demultiplex settings
DEMULTIPLEX_ACCEPTED_SAMPLE_INDEX_LENGTH: list[int] = [8, 10]
"""New defaults for dual-index."""
DEMULTIPLEX_BARCODE_LENGTH: int = 14
DEMULTIPLEX_INVALID_SAMPLE_INDEX: str = "X"
MAX_INDICES_TO_DEMUX: int = 1500
"""Maximum number of sample indices allowed."""
SAMPLE_INDEX_MAX_HAMMING_DISTANCE: int = 1
"""Maximum separation in HD between SIs of the same kind i7/i5."""

# TAG names
PROCESSED_BARCODE_TAG: str = "BX"
SAMPLE_INDEX_TAG: str = "BC"
SAMPLE_INDEX_QUAL_TAG: str = "QT"
TRIM_TAG: str = "TR"
TRIM_QUAL_TAG: str = "TQ"

PARALLEL_LOCUS_SIZE: int = int(4e7)
"""Parallelization settings."""

#
# Settings for metrics computation
#

READ_MATE_FAR_DIST: int = 5000
"""Distance to consider reads far away for far chimeras."""

# Preflight constants
BCL_PROCESSOR_FASTQ_MODE: str = "BCL_PROCESSOR"
ILMN_BCL2FASTQ_FASTQ_MODE: str = "ILMN_BCL2FASTQ"
