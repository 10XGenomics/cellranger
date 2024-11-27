#!/usr/bin/env python
#
# Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
#

import martian
import numpy as np

import cellranger.csv_io as cr_csv_io
import cellranger.h5_constants as h5_constants
import cellranger.matrix as cr_matrix
import cellranger.rna.library as rna_library
from cellranger import cell_calling_helpers
from cellranger.library_constants import ATACSEQ_LIBRARY_TYPE
from cellranger.logperf import LogPerf

__MRO__ = """
stage PREPROCESS_MATRIX(
    in  h5   matrix_h5,
    in  int  random_seed,
    in  csv  use_genes,
    in  csv  exclude_genes,
    in  csv  use_bcs,
    in  int  num_bcs,
    in  int  force_cells,
    in  bool get_peak_matrix,
    in  bool skip,
    in  bool is_visium_hd,
    in  bool is_pd,
    in  bool is_antibody_only,
    in  bool disable_run_pca,
    in  bool disable_correct_chemistry_batch,
    in  bool skip_multigenome_analysis,
    in  bool enable_tsne,
    out bool skip_antibody_analysis,
    out bool skip_antigen_analysis,
    out h5   cloupe_matrix_h5,
    out h5   preprocessed_matrix_h5,
    out bool is_multi_genome,
    out bool skip,
    out bool skip_tsne,
    out bool is_antibody_only,
    out bool disable_run_pca,
    out bool disable_correct_chemistry_batch,
    out bool skip_multigenome_analysis,
    out bool disable_hierarchical_clustering,
    src py   "stages/analyzer/preprocess_matrix",
) split (
) using (
    volatile = strict,
)
"""


def select_barcodes_and_features(
    matrix,
    num_bcs=None,
    use_bcs=None,
    use_genes=None,
    exclude_genes=None,
    force_cells=None,
):
    if force_cells is not None:
        bc_counts = matrix.get_counts_per_bc()
        bc_indices, _, _ = cell_calling_helpers.filter_cellular_barcodes_fixed_cutoff(
            bc_counts, force_cells
        )
        with LogPerf("f1"):
            matrix = matrix.select_barcodes(bc_indices)

    elif use_bcs is not None:
        bc_seqs = cr_csv_io.load_csv_rownames(use_bcs)
        missing_bcs = []
        bcs_in_matrix = set(matrix.bcs)
        for bc in bc_seqs:
            if bc not in bcs_in_matrix:
                missing_bcs.append(bc)
        if missing_bcs:  # is not empty
            martian.exit(
                "ERROR: The following barcodes are not present in the feature-barcode matrix:\n{}\n\nPlease check the input csv file which lists the barcodes to use.".format(
                    b"\n".join(missing_bcs).decode()
                )
            )
        bc_indices = matrix.bcs_to_ints(bc_seqs)
        with LogPerf("f2"):
            matrix = matrix.select_barcodes(bc_indices)

    elif num_bcs is not None and num_bcs < matrix.bcs_dim:
        bc_indices = np.sort(
            np.random.choice(np.arange(matrix.bcs_dim), size=num_bcs, replace=False)
        )
        with LogPerf("f3"):
            matrix = matrix.select_barcodes(bc_indices)

    include_indices = list(range(matrix.features_dim))
    if use_genes is not None:
        include_ids = cr_csv_io.load_csv_rownames(use_genes)
        with LogPerf("f4"):
            try:
                include_indices = matrix.feature_ids_to_ints(include_ids)
            except KeyError as e:
                err_message = str(e).strip("'")
                err_message = err_message.replace(
                    "Specified feature ID", "Feature ID specified in genes_csv"
                )
                martian.exit(err_message)

    exclude_indices = []
    if exclude_genes is not None:
        exclude_ids = cr_csv_io.load_csv_rownames(exclude_genes)
        with LogPerf("f5"):
            try:
                exclude_indices = matrix.feature_ids_to_ints(exclude_ids)
            except KeyError as e:
                err_message = str(e).strip("'")
                err_message = err_message.replace(
                    "Specified feature ID", "Feature ID specified in exclude_genes_csv"
                )
                martian.exit(err_message)

    gene_indices = np.array(sorted(list(set(include_indices) - set(exclude_indices))), dtype=int)
    with LogPerf("ff"):
        matrix = matrix.select_features(gene_indices)

    return matrix


def split(args):
    assert (not args.is_antibody_only) or (not args.get_peak_matrix)
    matrix_mem_gb = cr_matrix.CountMatrix.get_mem_gb_from_matrix_h5(args.matrix_h5)
    return {"chunks": [], "join": {"__mem_gb": max(matrix_mem_gb, h5_constants.MIN_MEM_GB)}}


def join(args, outs, chunk_defs, chunk_outs):
    outs.skip = args.skip
    outs.skip_tsne = (not args.enable_tsne) or outs.skip or args.is_visium_hd
    outs.is_antibody_only = args.is_antibody_only
    outs.disable_run_pca = args.disable_run_pca
    outs.disable_hierarchical_clustering = (
        outs.skip or not args.is_visium_hd or args.disable_run_pca or not args.is_pd
    )
    outs.disable_correct_chemistry_batch = args.disable_correct_chemistry_batch
    outs.skip_multigenome_analysis = args.skip_multigenome_analysis
    outs.skip_antibody_analysis = False
    outs.skip_antigen_analysis = False
    if args.skip:
        return

    if args.random_seed is not None:
        np.random.seed(args.random_seed)

    # detect barnyard
    genomes = cr_matrix.CountMatrix.get_genomes_from_h5(args.matrix_h5)
    if len(genomes) > 1:
        outs.is_multi_genome = True
    else:
        outs.is_multi_genome = False

    with LogPerf("load"):
        matrix = cr_matrix.CountMatrix.load_h5_file(args.matrix_h5)

    with LogPerf("select"):
        library_type = rna_library.GENE_EXPRESSION_LIBRARY_TYPE
        library_types = cr_matrix.CountMatrix.load_library_types_from_h5_file(args.matrix_h5)
        is_antibody_only = args.is_antibody_only or (
            rna_library.GENE_EXPRESSION_LIBRARY_TYPE not in library_types
            and rna_library.ANTIBODY_LIBRARY_TYPE in library_types
        )  # No GEX features found
        if is_antibody_only:
            library_type = rna_library.ANTIBODY_LIBRARY_TYPE
        if args.get_peak_matrix:
            library_type = ATACSEQ_LIBRARY_TYPE

        matrix = select_barcodes_and_features(
            matrix,
            num_bcs=args.num_bcs,
            use_bcs=args.use_bcs,
            use_genes=args.use_genes,
            exclude_genes=args.exclude_genes,
            force_cells=args.force_cells,
        )

    # Preserve original matrix attributes
    matrix_attrs = cr_matrix.load_matrix_h5_metadata(args.matrix_h5)
    # Including user-specified attrs
    matrix_attrs.update(cr_matrix.load_matrix_h5_custom_attrs(args.matrix_h5))

    # gem groups are needed for cloupe, and older versions of cellranger count
    # may not have added those to the matrix_attrs
    if cr_matrix.get_gem_group_index(args.matrix_h5) is None and (
        library_type is not ATACSEQ_LIBRARY_TYPE
    ):
        # TODO: Is this is still necessary?
        # Not sure when this broke but the sample_id is not piped to the analyzer.
        martian.throw(
            "Count matrix is missing necessary attributes, possibly because it was produced by an old version of Cellranger."
        )

        # chemistry = cr_matrix.CountMatrix.load_chemistry_from_h5(args.matrix_h5)
        # bcs = cr_matrix.CountMatrix.load_bcs_from_h5(args.matrix_h5)
        # gem_groups = list(set(cr_utils.split_barcode_seq(bc)[1] for bc in bcs))
        # count_attrs = cr_matrix.make_matrix_attrs_count(args.sample_id, gem_groups, chemistry)
        # matrix_attrs.update(count_attrs)

    # matrix h5 for cloupe (gene with zero count preserved)
    # this will only be used in reanalyzer, where user could include/exclude genes
    with LogPerf("w1"):
        matrix.save_h5_file(
            outs.cloupe_matrix_h5,
            extra_attrs=matrix_attrs,
            sw_version=martian.get_pipelines_version(),
        )

    with LogPerf("selnz"):
        matrix, _, _ = matrix.select_nonzero_axes()

    with LogPerf("w2"):
        matrix.save_h5_file(
            outs.preprocessed_matrix_h5,
            extra_attrs=matrix_attrs,
            sw_version=martian.get_pipelines_version(),
        )

    _, total_bcs = matrix.get_shape()
    # if we don't have enough barcodes after processing, skip!
    if total_bcs <= 1:
        martian.log_info(
            f"Feature-barcode matrix is tiny (num_cells = {total_bcs}) - skipping analysis."
        )
        outs.skip = True
        outs.disable_run_pca = True
        outs.disable_correct_chemistry_batch = True
        outs.skip_multigenome_analysis = True
        outs.disable_hierarchical_clustering = True
        return

    feature_cnts = matrix.get_count_of_feature_types()
    ab_cnts = feature_cnts.get(rna_library.ANTIBODY_LIBRARY_TYPE, 0)
    if ab_cnts == 0:
        outs.skip_antibody_analysis = True
    ag_cnts = feature_cnts.get(rna_library.ANTIGEN_LIBRARY_TYPE, 0)
    if ag_cnts == 0:
        outs.skip_antigen_analysis = True
    # Skip ANTIGEN_ANALYZER in SC_RNA_ANALYZER in aggr, it is called in VDJ AGGR
    if (
        library_type is not ATACSEQ_LIBRARY_TYPE
        and len(cr_matrix.get_gem_group_index(args.matrix_h5)) > 1
    ):
        outs.skip_antigen_analysis = True
    # if we don't have enough features after processing, skip!
    if all(v <= 1 for (_, v) in feature_cnts.items()):
        martian.log_info(
            "Feature-barcode matrix is tiny (fewer than 2 features per library type) - skipping analysis."
        )
        outs.skip = True
        outs.disable_run_pca = True
        outs.disable_hierarchical_clustering = True
        outs.disable_correct_chemistry_batch = True
        outs.skip_multigenome_analysis = True
