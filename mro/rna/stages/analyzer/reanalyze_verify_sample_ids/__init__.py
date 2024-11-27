#!/usr/bin/env python
#
# Copyright (c) 2021 10x Genomics, Inc. All rights reserved.
#

"""For SC_RNA_REANALYZER pipeline runs on Aggr inputs.

Sanity check that the library_ids from the matrix match the sample_ids from the aggr sample defs.
"""

import martian
from six import ensure_binary, ensure_str

import cellranger.constants as cr_constants
import cellranger.matrix as cr_matrix

__MRO__ = """
stage REANALYZE_VERIFY_SAMPLE_IDS(
    in  h5    matrix_h5,
    in  map[] sample_defs,
    out map[] sample_defs,
    src py    "stages/analyzer/reanalyze_verify_sample_ids",
)
"""

AGG_SAMPLE_ID_FIELD = "sample_id"


def main(args, outs):
    library_map = cr_matrix.get_gem_group_index(args.matrix_h5)
    matrix_library_ids = {library_id for library_id, _ in library_map.values()}

    csv_library_ids = set()
    for row in args.sample_defs:
        if cr_constants.AGG_ID_FIELD in row:
            csv_library_ids.add(ensure_binary(row[cr_constants.AGG_ID_FIELD]))
        else:
            csv_library_ids.add(ensure_binary(row[AGG_SAMPLE_ID_FIELD]))

    if matrix_library_ids != csv_library_ids:
        str_csv_library_ids = ",".join([ensure_str(x) for x in sorted(csv_library_ids)])
        str_matrix_library_ids = ",".join([ensure_str(x) for x in sorted(matrix_library_ids)])
        this_msg = f"Sample IDs specified in CSV ({str_csv_library_ids}) do not match those contained in the input matrix ({str_matrix_library_ids})"
        martian.exit(this_msg)

    # output the sample defs so we can ensure downstream stages execute after this one
    outs.sample_defs = args.sample_defs
