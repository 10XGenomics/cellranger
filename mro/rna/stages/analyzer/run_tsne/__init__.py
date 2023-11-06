#!/usr/bin/env python
#
# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
#

import itertools

import numpy as np

import cellranger.analysis.bhtsne as cr_tsne
import cellranger.analysis.constants as analysis_constants
import cellranger.analysis.io as analysis_io
import cellranger.analysis.pca as cr_pca
import cellranger.cr_io as cr_io
import cellranger.h5_constants as h5_constants
import cellranger.matrix as cr_matrix
import cellranger.rna.library as rna_library

__MRO__ = """
stage RUN_TSNE(
    in  h5    matrix_h5,
    in  h5    pca_h5,
    in  int   random_seed,
    in  int   perplexity,
    in  int   input_pcs,
    in  int   max_dims,
    in  int   max_iter,
    in  int   stop_lying_iter,
    in  int   mom_switch_iter,
    in  float theta,
    in  bool  is_antibody_only,
    out h5    tsne_h5,
    out path  tsne_csv,
    src py    "stages/analyzer/run_tsne",
) split using (
    in  int   tsne_dims,
    in  string feature_type,
)
"""


def split(args):
    feature_ref = cr_matrix.CountMatrix.load_feature_ref_from_h5_file(args.matrix_h5)
    feature_types = sorted(list({fd.feature_type for fd in feature_ref.feature_defs}))

    chunks = []
    matrix_mem_gb = cr_matrix.CountMatrix.get_mem_gb_from_matrix_h5(args.matrix_h5)
    min_tsne_dims = analysis_constants.TSNE_N_COMPONENTS
    max_tsne_dims = args.max_dims if args.max_dims is not None else min_tsne_dims
    for tsne_dims, feature_type in itertools.product(
        range(min_tsne_dims, max_tsne_dims + 1), feature_types
    ):
        chunks.append(
            {
                "tsne_dims": tsne_dims,
                "feature_type": feature_type,
                "__mem_gb": max(matrix_mem_gb, h5_constants.MIN_MEM_GB),
            }
        )
    return {"chunks": chunks, "join": {"__mem_gb": 1}}


def lower_no_space(s):
    return s.replace(" ", "_").lower()


def get_tsne_name(feature_type, n_components):
    return "%s_%d-d" % (lower_no_space(feature_type), n_components)


def get_tsne_key(feature_type, n_components):
    if feature_type == rna_library.DEFAULT_LIBRARY_TYPE:
        # Preserve backward HDF5 compatibility with pre-3.0 HDF5 files
        #   where the CSV directory was named "2_components" and the HDF5 dataset was named "_2"
        return str(n_components)
    else:
        return "%s_%d" % (lower_no_space(feature_type), n_components)


def main(args, outs):
    tsne_dims = args.tsne_dims

    if args.is_antibody_only:
        matrix = cr_matrix.CountMatrix.load_h5_file(args.matrix_h5)
        matrix = matrix.select_features_by_type(rna_library.ANTIBODY_LIBRARY_TYPE)
        matrix.m.data = np.log2(1 + matrix.m.data)
        tsne_input = matrix.m.transpose().toarray()

        name = get_tsne_name(rna_library.ANTIBODY_LIBRARY_TYPE, args.tsne_dims)
        key = str(args.tsne_dims)

        tsne = cr_tsne.run_tsne(
            tsne_input,
            name=name,
            key=key,
            input_pcs=args.input_pcs,
            perplexity=args.perplexity,
            theta=args.theta,
            tsne_dims=tsne_dims,
            max_iter=args.max_iter,
            stop_lying_iter=args.stop_lying_iter,
            mom_switch_iter=args.mom_switch_iter,
            random_state=args.random_seed,
        )

        cr_tsne.save_tsne_h5(tsne, outs.tsne_h5)
        cr_tsne.save_tsne_csv(tsne, matrix.bcs, outs.tsne_csv)
        return

    if args.feature_type == rna_library.GENE_EXPRESSION_LIBRARY_TYPE:
        # Use PCA for gene expression
        pca = cr_pca.load_pca_from_h5(args.pca_h5)[0]
        tsne_input = pca.transformed_pca_matrix
    else:
        # Use feature space for other feature types
        # Assumes other feature types are much lower dimension than gene expression
        matrix = cr_matrix.CountMatrix.load_h5_file(args.matrix_h5)
        matrix = matrix.select_features_by_type(args.feature_type)
        matrix.m.data = np.log2(1 + matrix.m.data)
        tsne_input = matrix.m.transpose().toarray()

    name = get_tsne_name(args.feature_type, args.tsne_dims)
    key = get_tsne_key(args.feature_type, args.tsne_dims)

    tsne = cr_tsne.run_tsne(
        tsne_input,
        name=name,
        key=key,
        input_pcs=args.input_pcs,
        perplexity=args.perplexity,
        theta=args.theta,
        tsne_dims=tsne_dims,
        max_iter=args.max_iter,
        stop_lying_iter=args.stop_lying_iter,
        mom_switch_iter=args.mom_switch_iter,
        random_state=args.random_seed,
    )

    cr_tsne.save_tsne_h5(tsne, outs.tsne_h5)

    matrix_bcs = cr_matrix.CountMatrix.load_bcs_from_h5_file(args.matrix_h5)
    cr_tsne.save_tsne_csv(tsne, matrix_bcs, outs.tsne_csv)


def join(args, outs, chunk_defs, chunk_outs):
    chunk_h5s = [chunk_out.tsne_h5 for chunk_out in chunk_outs]
    chunk_csv_dirs = [chunk_out.tsne_csv for chunk_out in chunk_outs]
    analysis_io.combine_h5_files(
        chunk_h5s, outs.tsne_h5, [analysis_constants.ANALYSIS_H5_TSNE_GROUP]
    )
    for csv_dir in chunk_csv_dirs:
        cr_io.hardlink_with_fallback(csv_dir, outs.tsne_csv)
