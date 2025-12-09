#!/usr/bin/env python
#
# Copyright (c) 2021 10x Genomics, Inc. All rights reserved.
#
"""Prepare inputs to run the PCA for batch correction."""

import numpy as np

import cellranger.rna.library as rna_library
from cellranger.library_constants import ATACSEQ_LIBRARY_TYPE
from cellranger.matrix import CountMatrix

__MRO__ = """
stage PCA_PREP(
    in  h5     matrix_h5,
    in  bool   is_antibody_only,
    in  bool   is_atac,
    out h5     filt_matrix,
    out string library_type,
    src py     "stages/analyzer/pca_prep",
) split (
) using (
    volatile = strict,
)
"""


def split(args):
    assert not (args.is_antibody_only and args.is_atac)

    (num_features, num_barcodes, nnz) = CountMatrix.load_dims_from_h5(args.matrix_h5)
    matrix_mem_gib = CountMatrix.get_mem_gb_from_matrix_dim(num_barcodes, nnz, scale=1)
    mem_gib = 1 + 1.7 * matrix_mem_gib
    print(f"{num_features=},{num_barcodes=},{nnz=},{matrix_mem_gib=},{mem_gib=}")
    return {"chunks": [], "join": {"__mem_gb": mem_gib}}


def join(args, outs, _chunk_defs, _chunk_outs):
    outs.library_type = (
        rna_library.ANTIBODY_LIBRARY_TYPE
        if args.is_antibody_only
        else ATACSEQ_LIBRARY_TYPE if args.is_atac else rna_library.GENE_EXPRESSION_LIBRARY_TYPE
    )
    matrix = CountMatrix.load_data_for_library_type_from_h5(args.matrix_h5, outs.library_type)

    # Take features with some minimum number of total UMIs in pseudobulk across all batches
    feature_indices = np.flatnonzero(matrix.get_counts_per_feature())

    # filter barcodes with zero count
    bc_indices = np.flatnonzero(matrix.select_features(feature_indices).get_counts_per_bc())
    matrix.select_barcodes(bc_indices).save_h5_file(outs.filt_matrix)
