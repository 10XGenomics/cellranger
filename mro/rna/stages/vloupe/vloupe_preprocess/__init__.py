#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

import json
import os
import os.path
import subprocess
import tempfile

import martian

import tenkit.log_subprocess as tk_subproc

__MRO__ = """
stage VLOUPE_PREPROCESS(
    in  string      pipestance_type,
    in  string      sample_id,
    in  string      sample_desc,
    in  pb          enclone_output,
    in  bool        disable_vloupe,
    in  string      beam_mode,
    in  csv         feature_reference,
    in  h5          feature_barcode_matrix,
    in  csv         antigen_specificity_scores,
    in  map<string> antigen_specificity_controls,
    out vloupe      output_for_vloupe,
    src py          "stages/vloupe/vloupe_preprocess",
) using (
    mem_gb = 15,
)
"""


def main(args, outs):
    """Run the vlconverter executable with inputs that should be available in the outs.

    folder at the end of the pipeline run.  This will generate "output_for_vloupe.vloupe"
    in the stage folder.

    Memory usage not expected to be excessive with this (thus no custom split/join
    as of yet).
    """
    if args.disable_vloupe:
        outs.output_for_vloupe = None
        return

    if args.enclone_output is None or not os.path.isfile(args.enclone_output):
        martian.log_info("Proto file missing - cannot make vloupe file")
        outs.output_for_vloupe = None
        return

    call = [
        "vlconverter",
        args.sample_id,
        args.pipestance_type,
        args.enclone_output,
        "--output",
        outs.output_for_vloupe,
        "--description",
        args.sample_desc,
    ]

    if all([args.beam_mode, args.feature_reference, args.feature_barcode_matrix]):
        antigen_args = [
            "--with-antigen-counts",
            "--feature-reference",
            args.feature_reference,
            "--feature-barcode-matrix",
            args.feature_barcode_matrix,
        ]

        if args.antigen_specificity_scores:
            antigen_args.append("--feature-specificity-score")
            antigen_args.append(args.antigen_specificity_scores)

        # antigen_specificity_controls looks like:
        #     {
        #         "ALLELE_NAME1": "ANTIGEN_NAME1"
        #         "ALLELE_NAME2": "ANTIGEN_NAME2"
        #         "no_allele":    "ANTIGEN_NAME3"
        #     }
        # The ALLELE "no_allele" means there isn't an allele.
        # This datastructure is written as above as json
        if args.antigen_specificity_controls:
            prefix = "feature-specificity-controls-"
            with tempfile.NamedTemporaryFile(prefix=prefix, mode="w", delete=False) as f:
                json.dump(args.antigen_specificity_controls, f)
                antigen_args.append("--feature-specificity-controls")
                antigen_args.append(f.name)

        call.extend(antigen_args)

    # the sample desc may be unicode, so send the whole
    # set of args str utf-8 to check_output
    unicode_call = [arg.encode("utf-8") for arg in call]

    # but keep the arg 'call' here because log_info inherently
    # attempts to encode the message... (TODO: should log_info
    # figure out the encoding of the input string)
    martian.log_info("Running vlconverter: {}".format(" ".join(call)))
    try:
        results = tk_subproc.check_output(unicode_call, stderr=subprocess.PIPE)
        martian.log_info(f"vlconverter output: {results}")
    except subprocess.CalledProcessError as e:
        outs.output_for_vloupe = None
        martian.throw(
            f"Could not generate .vloupe file:\nstdout:\n{e.stdout}\n\nstderr:\n{e.stderr}"
        )
