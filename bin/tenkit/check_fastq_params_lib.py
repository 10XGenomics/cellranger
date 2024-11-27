#!/usr/bin/env python3
#
# Copyright (c) 2016 10x Genomics, Inc. All rights reserved.
#

from __future__ import annotations

import argparse
import sys

from six import ensure_str

from tenkit.fasta import check_fastq_types_multipath


def make_parser():
    parser = argparse.ArgumentParser(
        description="Check demux scheme and sample name, if applicable."
    )
    parser.add_argument("--fastqs", help="Input FASTQ path", required=True)
    parser.add_argument(
        "--fastqprefix", default=None, help="Prefix match for bcl2fastq/mkfastq-generated FASTQs."
    )
    return parser


def get_mro_params(fastq_paths, fastqprefixes):
    try:
        fastq_mode, sample_name = check_fastq_types_multipath(fastq_paths, fastqprefixes)
        print(f"{fastq_mode}\t{ensure_str(sample_name)}")
        sys.exit(0)
    except ValueError as ex:
        sys.stderr.write(f"{ex!s}\n")
        sys.exit(1)


def main():
    parser = make_parser()
    args = parser.parse_args()

    fastq_paths = args.fastqs.split(",")

    if args.fastqprefix:
        fastqprefixes = args.fastqprefix.split(",")
    else:
        fastqprefixes = None
    get_mro_params(fastq_paths, fastqprefixes)


if __name__ == "__main__":
    main()
