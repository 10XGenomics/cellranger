#!/usr/bin/env python
#
# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
#

import martian

import cellranger.analysis.constants as analysis_constants
import cellranger.analysis.pca as cr_pca
import cellranger.h5_constants as h5_constants
import cellranger.matrix as cr_matrix
import cellranger.rna.library as rna_library

__MRO__ = """
stage RUN_PCA(
    in  h5   matrix_h5,
    in  int  random_seed,
    in  int  num_pca_bcs,
    in  int  num_pca_genes,
    in  int  num_principal_comps,
    in  bool is_antibody_only,
    out h5   pca_h5,
    out path pca_csv,
    src py   "stages/analyzer/run_pca",
) split using (
)
"""


def split(args):
    # memory usage of just loading the full matrix

    matrix_dims = cr_matrix.CountMatrix.load_dims_from_h5(args.matrix_h5)
    (feature_dim, bcs_dim, nonzero_entries) = matrix_dims
    matrix_mem_gb = cr_matrix.CountMatrix.get_mem_gb_from_matrix_dim(bcs_dim, nonzero_entries)

    num_pcs = (
        analysis_constants.PCA_N_COMPONENTS_DEFAULT
        if args.num_principal_comps is None
        else args.num_principal_comps
    )

    irlb_mem_gb = cr_pca.get_irlb_mem_gb_from_matrix_dim(
        nonzero_entries, feature_dim, bcs_dim, num_pcs
    )
    mem_gb = max(irlb_mem_gb, matrix_mem_gb, h5_constants.MIN_MEM_GB)
    return {"chunks": [], "join": {"__mem_gb": mem_gb}}


def join(args, outs, _chunk_defs, _chunk_outs):
    matrix = cr_matrix.CountMatrix.load_h5_file(args.matrix_h5)
    if args.is_antibody_only:
        matrix = matrix.select_features_by_type(rna_library.ANTIBODY_LIBRARY_TYPE)
    else:
        matrix = matrix.select_features_by_type(rna_library.GENE_EXPRESSION_LIBRARY_TYPE)

    try:
        pca = cr_pca.run_pca(
            matrix,
            pca_features=args.num_pca_genes,
            pca_bcs=args.num_pca_bcs,
            n_pca_components=args.num_principal_comps,
            random_state=args.random_seed,
            min_count_threshold=cr_pca.DEFAULT_RUNPCA_THRESHOLD,
        )
    except (cr_matrix.NullAxisMatrixError, cr_pca.MatrixRankTooSmallException):
        martian.log_warn(
            "insufficient counts for min_count_threshold=2, downgrading to min_count_threshold=0"
        )
        try:
            pca = cr_pca.run_pca(
                matrix,
                pca_features=args.num_pca_genes,
                pca_bcs=args.num_pca_bcs,
                n_pca_components=args.num_principal_comps,
                random_state=args.random_seed,
                min_count_threshold=0,
            )
        except ValueError as e:
            martian.exit(e)

    pca_key = (
        args.num_principal_comps
        if args.num_principal_comps is not None
        else analysis_constants.PCA_N_COMPONENTS_DEFAULT
    )
    pca_map = {pca_key: pca}
    cr_pca.save_pca_h5(pca_map, outs.pca_h5)
    cr_pca.save_pca_csv(pca_map, matrix, outs.pca_csv)
