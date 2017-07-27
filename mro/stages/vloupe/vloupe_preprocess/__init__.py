#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

import os.path
import martian
import subprocess

__MRO__ = """
stage VLOUPE_PREPROCESS(
    in  string    pipestance_type,
    in  string    sample_id,
    in  string    sample_desc,
    in  bam       concat_ref_bam,
    in  bam.bai   concat_ref_bam_bai,
    in  fasta     concat_ref_fasta,
    in  fasta.fai concat_ref_fasta_fai,
    in  json      concat_ref_annotations_json,
    in  csv       clonotypes_csv,
    in  bam       consensus_bam,
    in  bam.bai   consensus_bam_bai,
    in  json      consensus_annotations_json,
    in  fasta     consensus_fasta,
    in  fasta.fai consensus_fasta_fai,
    in  string    contig_bam_relative_path,
    in  bam.bai   contig_bam_bai,
    in  json      contig_annotations_json,
    in  bed       contig_annotations_bed,
    in  fasta     contig_fasta,
    in  fasta.fai contig_fasta_fai,
    in  csv       metrics_csv,
    out vloupe    output_for_vloupe,
    src py        "stages/vloupe/vloupe_preprocess",
)
"""


def main(args, outs):
    """
    Run the vlconverter executable with inputs that should be available in the outs
    folder at the end of the pipeline run.  This will generate "output_for_vloupe.vloupe"
    in the stage folder.

    Memory usage not expected to be excessive with this (thus no custom split/join
    as of yet); it will need to load a few full files (bam.bai, fasta.fai) into memory.
    """
    if not os.path.isfile(args.concat_ref_bam) or \
       not os.path.isfile(args.consensus_bam) or \
       not os.path.isfile(args.contig_bam_bai):
        martian.log_info('One or more bam files missing - cannot make vloupe file')
        return
   
    call = ["vlconverter",
            args.sample_id,
            args.pipestance_type,
            "--output", outs.output_for_vloupe,
            "--reference-bam", args.concat_ref_bam,
            "--reference-bam-index", args.concat_ref_bam_bai,
            "--reference-fasta", args.concat_ref_fasta,
            "--reference-fasta-index", args.concat_ref_fasta_fai,
            "--reference-annotations", args.concat_ref_annotations_json,
            "--clonotypes", args.clonotypes_csv,
            "--consensus-bam", args.consensus_bam,
            "--consensus-bam-index", args.consensus_bam_bai,
            "--consensus-annotations", args.consensus_annotations_json,
            "--consensus-fasta", args.consensus_fasta,
            "--consensus-fasta-index", args.consensus_fasta_fai,
            "--contig-bam-relative-path", args.contig_bam_relative_path,
            "--contig-bam-index", args.contig_bam_bai,
            "--contig-annotations", args.contig_annotations_json,
            "--contig-bed", args.contig_annotations_bed,
            "--contig-fasta", args.contig_fasta,
            "--contig-fasta-index", args.contig_fasta_fai,
            "--description", args.sample_desc]

    martian.log_info("Running vlconverter: %s" % " ".join(call))
    try:
        results = subprocess.check_output(call)
        martian.log_info("vlconverter output: %s" % results)
    except subprocess.CalledProcessError, e:
        outs.output_for_vloupe = None
        martian.throw("Could not generate .vloupe file: \n%s" % e.output)
