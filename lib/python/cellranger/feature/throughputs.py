#
# Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
#
"""Constant values used by compute code and websummary code, put here to allow their use without lots of other.

imports
"""

MT_THROUGHPUT = "MT"
LT_THROUGHPUT = "LT"
HT_THROUGHPUT = "HT"
GEMX_THROUGHPUT = "GEMX"
THROUGHPUT_INFERRED_METRIC = "throughput_inferred"
INCONSISTENT_THROUGHPUT_METRIC = "inconsistent_throughput"  # keep in sync with name in metrics.csv
N_G = 95000  # Number of GEMs in a single NextGem gem-well
CORR_FACTOR = 1.54  # correction factor - if N_L  = cells loaded, we recover N_L/C cells
G19_N_GEMS = {
    MT_THROUGHPUT: N_G,
    LT_THROUGHPUT: 9500,
    HT_THROUGHPUT: 190000,
    GEMX_THROUGHPUT: 190000,
}
