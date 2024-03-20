#!/usr/bin/env python
#
# Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
#
"""Run the reference builder to generate a 10X-compatible reference."""
import martian

import cellranger.vdj.reference as cr_vdj_ref
import cellranger.vdj.reference_maker as ref_maker

__MRO__ = """
stage _MAKE_VDJ_REFERENCE(
    in  string genome_name,
    in  fasta  fasta_file,
    in  gtf[]  gtf_files,
    in  fasta  seq_file,
    in  file   remove_transcripts_file,
    in  string ref_version,
    in  string mkref_version,
    in  int    mem_gb,
    out path   reference,
    src py     "stages/make_vdj_reference",
) split (
) using (
    volatile = strict,
)
"""


def split(args):
    """Set memory resources."""
    return {"chunks": [], "join": {"__mem_gb": args.mem_gb}}


def join(args, outs, _chunk_defs, _chunk_outs):
    assert args.genome_name
    assert "/" not in args.genome_name

    try:
        if args.fasta_file:
            assert args.gtf_files
            assert not args.seq_file

            ref_maker.build_reference_fasta_from_ensembl(
                gtf_paths=args.gtf_files,
                transcripts_to_remove_path=args.remove_transcripts_file,
                genome_fasta_path=args.fasta_file,
                reference_path=outs.reference,
                reference_name=args.genome_name,
                ref_version=args.ref_version,
                mkref_version=args.mkref_version,
            )
        else:
            assert args.seq_file
            assert not args.gtf_files
            assert not args.remove_transcripts_file

            cr_vdj_ref.build_reference_fasta_from_fasta(
                fasta_path=args.seq_file,
                reference_path=outs.reference,
                reference_name=args.genome_name,
                ref_version=args.ref_version,
                mkref_version=args.mkref_version,
            )
    # Catch any errors in VDJ reference construction and print the error without
    # the stack trace.
    except cr_vdj_ref.VDJReferenceConstructionError as error:
        martian.exit(error)
