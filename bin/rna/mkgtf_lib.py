#!/usr/bin/env python3
#
# Copyright (c) 2016 10x Genomics, Inc. All rights reserved.
#

"""Genes GTF tool for 10x Genomics {product}.

Filter user-supplied GTF files for use as {product}-compatible
genes files for mkref tool.

The commands below should be preceded by '{cmd}':

Usage:
    mkgtf <input_gtf> <output_gtf> [--attribute=KEY:VALUE...]
    mkgtf -h | --help | --version

Arguments:
    input_gtf           Path to input genes GTF file.
    output_gtf          Path to filtered output genes GTF file.

Options:
    --attribute=<key:value>
                        Key-value pair in attributes field to be kept in the GTF
                            file.
    -h --help           Show this message.
    --version           Show version.
"""

from __future__ import annotations

import collections
import os
import sys

import docopt

import cellranger.cr_io as cr_io
import cellranger.reference as cr_reference
from cellranger.products import get_cmd_names


def _parse_args(product_name):
    version = "{} {} {}\n{}".format(
        product_name,
        os.getenv("TENX_SUBCMD", ""),
        os.getenv("TENX_VERSION", ""),
        os.getenv("TENX_COPYRIGHT", ""),
    )
    product, cmd = get_cmd_names(product_name)

    return docopt.docopt(__doc__.format(product=product, cmd=cmd), version=version)


def main():
    args = _parse_args(os.getenv("TENX_PRODUCT", ""))
    input_genes_file = cr_io.get_input_path(args["<input_gtf>"])
    output_genes_file = cr_io.get_output_path(args["<output_gtf>"])
    attributes_str = args["--attribute"]

    attributes = collections.defaultdict(set)
    for attribute_str in attributes_str:
        parts = attribute_str.split(":")
        if len(parts) != 2:
            sys.exit(f"Attribute option must have format <key;value>: {attribute_str}")
        key, value = parts
        attributes[key].add(value)

    gtf_builder = cr_reference.GtfBuilder(
        input_genes_file, output_genes_file, attributes=attributes
    )
    gtf_builder.build_gtf()


if __name__ == "__main__":
    main()
