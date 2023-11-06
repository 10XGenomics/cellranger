#
# Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
#
"""Constant values used by some code (e.g. websummary) that doesn't need the rest of the infrastructure imported."""
import cellranger.analysis.constants as analysis_constants

MULTIPLETS_FACTOR_NAME = analysis_constants.GEM_CLASS_MULTIPLET.decode()
UNASSIGNED_FACTOR_NAME = "Unassigned"
BLANK_FACTOR_NAME = "Blank"
ASSIGNMENT_COL_NAME = "Assignment"
ASSIGNMENT_PROB_COL_NAME = "Assignment_Probability"
JIBES_MIN_CONFIDENCE = 0.9  # the minimum posterior probability a state needs to have for the barcode in question to be assigned to that state.
# If the most likely state has lower probability than this value, the barcode will be declared UNASSIGNED_FACTOR_NAME
