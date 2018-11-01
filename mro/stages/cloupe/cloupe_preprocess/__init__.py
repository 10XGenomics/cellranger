#!/usr/bin/env python
#
# Copyright (c) 2018 10X Genomics, Inc. All rights reserved.
#

import cellranger.h5_constants as h5_constants
import cellranger.matrix as cr_matrix
import cellranger.io as cr_io
import martian
import subprocess
import os
import tenkit.log_subprocess as tk_subproc
import tenkit.safe_json as tk_json

__MRO__ = """
stage CLOUPE_PREPROCESS(
    in  string pipestance_type,
    in  string sample_id,
    in  string sample_desc,
    in  path   analysis,
    in  h5     filtered_gene_bc_matrices_h5,
    in  json   metrics_json,
    in  csv    aggregation_csv,
    in  json   gem_group_index_json,
    in  bool   no_secondary_analysis,
    out cloupe output_for_cloupe,
    out json   gem_group_index_json,
    src py     "stages/cloupe/cloupe_preprocess",
) split using (
)
"""

def get_gem_group_index_json(args, outs):
    if args.gem_group_index_json:
        cr_io.copy(args.gem_group_index_json, outs.gem_group_index_json)
        return outs.gem_group_index_json
    else:
        generated_index = cr_matrix.get_gem_group_index(args.filtered_gene_bc_matrices_h5)
        if generated_index:
            with open(outs.gem_group_index_json, 'w') as outfile:
                tk_json.dump_numpy({"gem_group_index": generated_index}, outfile)
            return outs.gem_group_index_json
        else:
            return None

def get_analysis_h5_path(args):
    return os.path.join(args.analysis, "analysis.h5")

def do_not_make_cloupe(args):
    """
    Returns True if there is a reason why this stage should not attempt to
    generate a .cloupe file
    """
    if args.no_secondary_analysis:
        martian.log_info("Skipping .cloupe generation by instruction (--no-secondary-analysis)")
        return True
    if args.analysis is None:
        martian.log_info("Skipping .cloupe generation due to missing analysis folder")
        return True
    if not os.path.exists(args.filtered_gene_bc_matrices_h5):
        martian.log_info("Skipping .cloupe generation due to missing or zero-length feature-barcode matrix")
        return True
    return False

def split(args):
    # no mem usage if skipped
    if do_not_make_cloupe(args):
        return {'chunks': [{'__mem_gb': h5_constants.MIN_MEM_GB}]}

    # Below is an old estimate (and comment) from before CELLRANGER-1667
    # CELLRANGER-762: worst case is that there's two copies of the sparse matrix,
    # one for reading, and one for writing
    old_matrix_mem_gb = 2 * cr_matrix.CountMatrix.get_mem_gb_from_matrix_h5(args.filtered_gene_bc_matrices_h5)
    old_matrix_mem_gb = max(old_matrix_mem_gb, h5_constants.MIN_MEM_GB, 6)

    # After CR-1667, we implemented a new estimate to deal with OOM errors, but to ensure backwards compatibility, we
    # select the max of the old and the new (the new underestimates for large matrices)
    new_matrix_mem_gb = cr_matrix.CountMatrix.get_mem_gb_crconverter_estimate_from_h5(args.filtered_gene_bc_matrices_h5)

    matrix_mem_gb = max(old_matrix_mem_gb, new_matrix_mem_gb)
    chunks = [{
        '__mem_gb': max(matrix_mem_gb, h5_constants.MIN_MEM_GB, 6)
    }]
    # Hack to increase min memgb - not sure why we need it
    return {'chunks': chunks, 'join': {'__mem_gb': 1}}

def join(args, outs, chunk_defs, chunk_outs):
    if chunk_outs[0].output_for_cloupe is None:
        # Set output to null if noloupe is set
        outs.output_for_cloupe = None
    else:
        cr_io.copy(chunk_outs[0].output_for_cloupe, outs.output_for_cloupe)

def main(args, outs):
    if do_not_make_cloupe(args):
        outs.output_for_cloupe = None
        return

    call = ["crconverter",
            args.sample_id,
            args.pipestance_type,
            "--matrix", args.filtered_gene_bc_matrices_h5,
            "--analysis", get_analysis_h5_path(args),
            "--output", outs.output_for_cloupe,
            "--description", args.sample_desc]

    if args.metrics_json:
        call.extend(["--metrics", args.metrics_json])
    if args.aggregation_csv:
        call.extend(["--aggregation", args.aggregation_csv])

    gem_group_index_json = get_gem_group_index_json(args, outs)
    if gem_group_index_json:
        call.extend(["--gemgroups", gem_group_index_json])
    else:
        # required argument for crconverter
        martian.exit("HDF5 matrix to be used for cloupe does not have GEM group information.")

    # the sample desc may be unicode, so send the whole
    # set of args str utf-8 to check_output
    unicode_call = [arg.encode('utf-8') for arg in call]

    # but keep the arg 'call' here because log_info inherently
    # attempts to encode the message... (TODO: should log_info
    # figure out the encoding of the input string)
    martian.log_info("Running crconverter: %s" % " ".join(call))
    try:
        results = tk_subproc.check_output(unicode_call)
        martian.log_info("crconverter output: %s" % results)
    except subprocess.CalledProcessError, e:
        outs.output_for_cloupe = None
        martian.throw("Could not generate .cloupe file: \n%s" % e.output)
