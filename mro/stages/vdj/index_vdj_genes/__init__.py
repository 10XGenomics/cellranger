#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
import os
import subprocess
import cellranger.vdj.reference as vdj_reference
import cellranger.utils as cr_utils

__MRO__ = """
stage INDEX_VDJ_GENES(
    in  path  vdj_reference_path,
    out fasta recombinome,
    out path  recombinome_index,
    src py    "stages/vdj/index_vdj_genes",
)
"""

def main(args, outs):
    # NOOP if no vdj ref path specified
    if args.vdj_reference_path is None:
        outs.recombinome = None
        outs.recombinome_index = None
        return

    fasta_filename = vdj_reference.get_vdj_reference_fasta(args.vdj_reference_path)
    cr_utils.copy(fasta_filename, outs.recombinome)
    os.makedirs(outs.recombinome_index)

    # Build a bowtie2 index
    subprocess.check_call(['bowtie2-build',
                           outs.recombinome,
                           os.path.join(outs.recombinome_index, 'recombinome')])
