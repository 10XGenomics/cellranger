#!/usr/bin/env python3
#
# Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
#
"""Infer Gem well throughput if there's no top level argument."""


import h5py

import cellranger.constants as cr_constants
import cellranger.feature.multiplexing.infer_throughput as it
import cellranger.feature.utils as feature_utils
import cellranger.rna.library as rna_library
import cellranger.utils as cr_utils
from cellranger.chemistry import (
    CHEMISTRY_DESCRIPTION_FIELD,
    CHEMISTRY_SC3P_LT,
    SC_GEMX_CHEMISTRIES,
)
from cellranger.feature.throughputs import (
    GEMX_THROUGHPUT,
    HT_THROUGHPUT,
    LT_THROUGHPUT,
    MT_THROUGHPUT,
    THROUGHPUT_INFERRED_METRIC,
)
from cellranger.feature_ref import HASHTAG_TAG
from cellranger.matrix import CountMatrix
from cellranger.reference_paths import get_reference_genomes

__MRO__ = """
stage INFER_GEM_WELL_THROUGHPUT(
    in  map<ChemistryDef> chemistry_defs,
    in  string            throughput,
    in  h5                filtered_feature_counts_matrix,
    in  path              reference_path,
    in  h5                barcode_summary_h5,
    out string            throughput,
    out json              inferred_throughputs,
    src py                "stages/feature/infer_gem_well_throughput",
) using (
    mem_gb   = 2,
    volatile = strict,
)
"""


def main(args, outs):
    feature_ref = CountMatrix.load_feature_ref_from_h5_file(args.filtered_feature_counts_matrix)
    no_cmos = (
        feature_ref.get_count_of_feature_type(rna_library.MULTIPLEXING_LIBRARY_TYPE) == 0
        and HASHTAG_TAG not in feature_ref.all_tag_keys
    )
    if no_cmos or not feature_utils.all_files_present([args.barcode_summary_h5]):
        outs.throughput = None
        outs.inferred_throughputs = None
        return

    # get barcode counts from barcode summary
    genomes = get_reference_genomes(args.reference_path)
    barcode_summary = h5py.File(args.barcode_summary_h5, "r")
    genome = genomes[0] if len(genomes) == 1 else rna_library.MULTI_REFS_PREFIX
    lib_prefix = rna_library.get_library_type_metric_prefix(
        rna_library.GENE_EXPRESSION_LIBRARY_TYPE
    )
    key = cr_utils.format_barcode_summary_h5_key(
        lib_prefix,
        genome,
        cr_constants.TRANSCRIPTOME_REGION,
        cr_constants.CONF_MAPPED_DEDUPED_READ_TYPE,
    )
    if key not in barcode_summary:
        outs.throughput = None
        return
    counts_per_bc = barcode_summary[key][:]
    counts_per_bc[::-1].sort()  # pylint: disable=no-member

    throughput_counts = it.infer_throughput_from_background_counts(counts_per_bc)
    slope_bc_idx, throughput_gradient = it.infer_throughput_from_rankplot_gradient(counts_per_bc)

    # Final inferral is the AND operator on both algorithm outputs
    throughput_final = (
        throughput_counts
        if throughput_counts == throughput_gradient == HT_THROUGHPUT
        else MT_THROUGHPUT
    )

    gex_chemistry_def = args.chemistry_defs[rna_library.GENE_EXPRESSION_LIBRARY_TYPE]
    gex_chemistry_description = gex_chemistry_def[CHEMISTRY_DESCRIPTION_FIELD]

    if gex_chemistry_description in [
        chem[CHEMISTRY_DESCRIPTION_FIELD] for chem in SC_GEMX_CHEMISTRIES
    ]:
        # Downstream code uses "throughput" to estimate droplet counts, so GEM-X
        # data needs an annotated "throughput" to handle this, even though
        # it no longer supports actual throughput differences.
        outs.throughput = GEMX_THROUGHPUT
        return

    inferred_throughputs = {
        "throughput_specified_by_chemistry": gex_chemistry_description,
        "throughput_specified_by_user": args.throughput,
        "throughput_inferred_from_counts": throughput_counts,
        "throughput_inferred_from_gradient": throughput_gradient,
        "throughput_steepest_gradient_bc_idx": slope_bc_idx,
        THROUGHPUT_INFERRED_METRIC: throughput_final,
    }

    # We now take throughput directly from --chemistry input
    if gex_chemistry_def == CHEMISTRY_SC3P_LT:
        outs.throughput = LT_THROUGHPUT
    elif gex_chemistry_description.endswith("HT"):
        outs.throughput = HT_THROUGHPUT
    else:
        outs.throughput = throughput_final

    inferred_throughputs["throughput_final_output"] = outs.throughput
    feature_utils.write_json_from_dict(inferred_throughputs, outs.inferred_throughputs)
