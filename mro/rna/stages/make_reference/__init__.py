#!/usr/bin/env python
#
# Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
#
"""Run the reference builder to generate a 10X-compatible reference."""
import martian

from cellranger.reference_builder import (
    GexReferenceError,
    GtfParseError,
    ReferenceBuilder,
)

__MRO__ = """
stage _MAKE_REFERENCE(
    in  string[] genome_names,
    in  fasta[]  fasta_files,
    in  gtf[]    gtf_files,
    in  string   ref_version,
    in  string   mkref_version,
    in  int      num_threads,
    in  int      mem_gb,
    out path     reference,
    src py       "stages/make_reference",
) split (
) using (
    volatile = strict,
)
"""


def split(args):
    """Set memory and CPU resources."""
    return {"chunks": [], "join": {"__mem_gb": args.mem_gb, "__threads": args.num_threads}}


def join(args, outs, _chunk_defs, _chunk_outs):
    assert args.genome_names
    assert len(args.genome_names) == len(args.fasta_files)
    assert len(args.genome_names) == len(args.gtf_files)
    assert all("/" not in genome for genome in args.genome_names)

    try:
        reference_builder = ReferenceBuilder(
            genomes=args.genome_names,
            in_fasta_fns=args.fasta_files,
            in_gtf_fns=args.gtf_files,
            out_dir=outs.reference,
            ref_version=args.ref_version,
            mkref_version=args.mkref_version,
            num_threads=args.num_threads,
            mem_gb=args.mem_gb,
        )
        reference_builder.build_gex_reference()
    except (GtfParseError, GexReferenceError) as ex:
        martian.exit(f"mkref has failed: error building reference package\n{ex}")
