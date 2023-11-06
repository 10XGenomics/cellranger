#!/usr/bin/env python
#
# Copyright (c) 2021 10x Genomics, Inc. All rights reserved.
#
"""Prepare inputs to run the PCA for batch correction."""

import martian
import numpy as np

import cellranger.h5_constants as h5_constants
import cellranger.matrix as cr_matrix
import cellranger.rna.library as rna_library
from cellranger.library_constants import ATACSEQ_LIBRARY_TYPE

__MRO__ = """
stage PCA_PREP(
    in  h5     matrix_h5,
    in  bool   is_antibody_only,
    in  bool   is_atac,
    out float  mem_gb,
    out h5     filt_matrix,
    out string library_type,
    src py     "stages/analyzer/pca_prep",
) split (
) using (
    volatile = strict,
)
"""


def split(args):
    assert (not args.is_antibody_only) or (not args.is_atac)  # cannot be atac and antibody_only
    # memory usage to load a h5 matrix
    matrix_dims = cr_matrix.CountMatrix.load_dims_from_h5(args.matrix_h5)
    (_, bcs_dim, nonzero_entries) = matrix_dims
    matrix_mem_gb = cr_matrix.CountMatrix.get_mem_gb_from_matrix_dim(bcs_dim, nonzero_entries)

    # convert from int to float,
    # and holding original matrix for conversion of PCA barcode space back costs us 3x memory
    mem_gb = max(np.ceil(3.0 * matrix_mem_gb), h5_constants.MIN_MEM_GB)
    return {"chunks": [], "join": {"__mem_gb": mem_gb}}


def join(args, outs, chunk_defs, chunk_outs):
    matrix = cr_matrix.CountMatrix.load_h5_file(args.matrix_h5)

    library_types = cr_matrix.CountMatrix.load_library_types_from_h5_file(args.matrix_h5)
    is_antibody_only = (
        rna_library.GENE_EXPRESSION_LIBRARY_TYPE not in library_types
        and rna_library.ANTIBODY_LIBRARY_TYPE in library_types
    )  # No GEX features found

    if is_antibody_only:
        matrix = matrix.select_features_by_type(rna_library.ANTIBODY_LIBRARY_TYPE)
        outs.library_type = rna_library.ANTIBODY_LIBRARY_TYPE
    elif args.is_atac:
        matrix = matrix.select_features_by_type(ATACSEQ_LIBRARY_TYPE)
        outs.library_type = ATACSEQ_LIBRARY_TYPE
    else:
        matrix = matrix.select_features_by_type(rna_library.GENE_EXPRESSION_LIBRARY_TYPE)
        outs.library_type = rna_library.GENE_EXPRESSION_LIBRARY_TYPE

    # Take features with some minimum number of total UMIs in pseudobulk across all batches
    feature_indices = np.flatnonzero(matrix.get_counts_per_feature())
    all_bcs_matrix = matrix.select_features(feature_indices)

    # filter barcodes with zero count
    bc_indices = np.flatnonzero(all_bcs_matrix.get_counts_per_bc())
    matrix = matrix.select_barcodes(bc_indices)

    matrix.save_h5_file(outs.filt_matrix)
    outs.mem_gb = martian.get_memory_allocation()
