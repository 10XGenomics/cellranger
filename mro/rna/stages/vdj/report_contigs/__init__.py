#!/usr/bin/env python
#
# Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

from __future__ import annotations

import json
import math
import os
import shutil
from collections import defaultdict

import pandas as pd
from six import ensure_binary, ensure_str

import cellranger.h5_constants as h5_constants
import cellranger.library_constants as lib_constants
import cellranger.utils as cr_utils
import cellranger.vdj.reference as vdj_ref
import cellranger.vdj.report as vdj_report
import cellranger.vdj.utils as vdj_utils
from cellranger.library_constants import MULTI_REFS_PREFIX

__MRO__ = """
stage REPORT_CONTIGS(
    in  path   vdj_reference_path,
    in  json   cell_barcodes,
    in  fasta  contigs,
    in  json   annotations,
    in  csv    filter_summary,
    in  tsv    umi_summary,
    in  string prefix,
    out json   summary,
    src py     "stages/vdj/report_contigs",
) split (
)
"""

MIN_CHAIN_TYPE_CONTIG_FRAC = 0.05

MEM_GB_PER_UMI_SUMMARY_GB = 7.0

LIBRARY_TYPE = lib_constants.VDJ_LIBRARY_TYPE


def split(args):
    mem_gb_annot = vdj_utils.get_mem_gb_from_annotations_json(args.annotations)

    umi_summary_bytes = os.path.getsize(args.umi_summary) if args.umi_summary else 0
    mem_gb_umi = math.ceil(MEM_GB_PER_UMI_SUMMARY_GB * float(umi_summary_bytes) / 1e9)

    mem_gb = max(mem_gb_annot, mem_gb_umi)
    mem_gb = max(mem_gb, 6)

    print("requested %d" % mem_gb)
    return {
        "chunks": [
            {
                "__mem_gb": max(h5_constants.MIN_MEM_GB, mem_gb),
            }
        ]
    }


def main(args, outs):
    prefix = ""
    if args.prefix is not None:
        prefix = ensure_str(args.prefix)
    print(prefix)
    reporter = vdj_report.VdjReporter(prefix)

    # Set a default value of 0 for number of paired cells so that it will be
    # present in the metric summary csv even when there are no paired cells
    # or in denovo mode
    reporter._get_metric_attr(
        prefix + "vdj_assembly_contig_pair_productive_full_len_bc_count", MULTI_REFS_PREFIX
    ).set_value(0)

    barcode_contigs: dict[bytes, list[tuple[bytes, bytes]]] = defaultdict(list)
    contig_annotations = {}

    # Get annotations for each contig
    for annotation in iter(json.load(open(args.annotations))):
        contig_annotations[ensure_binary(annotation["contig_name"])] = annotation

    file_size = os.path.getsize(args.umi_summary)
    if args.umi_summary and os.path.isfile(args.umi_summary) and file_size > 0:
        try:
            umi_summary = pd.read_csv(
                args.umi_summary,
                header=0,
                index_col=None,
                sep="\t",
                converters={"barcode": ensure_binary},
            )
            umi_summary["barcode"] = umi_summary["barcode"].astype("S")
            umi_summary = umi_summary.groupby("barcode")
        except KeyError as exc:
            print("Key error when reading umi_summary file.  Diagnostic information below:")
            print(f"File size is: {file_size}")
            print("Decoded columns are: " + ",".join(x for x in umi_summary.columns))
            with open(args.umi_summary) as in_file:
                data = in_file.read()
                print("Start of File:\n" + data[:50] + "\n")
                print("End of File:\n" + data[-50:] + "\n")
            print("Please contact support@10xgenomics.com to report this problem.")
            raise exc
    else:
        umi_summary = None

    if args.filter_summary:
        filter_summary = vdj_utils.load_contig_summary_table(args.filter_summary)
    else:
        filter_summary = None

    # Get contigs for each barcode
    with open(args.contigs, "rb") as contigs_file:
        for contig_hdr, contig_seq in cr_utils.get_fasta_iter(contigs_file):
            contig_name = contig_hdr.split(b" ")[0]
            if not filter_summary is None and not vdj_utils.is_contig_filtered(
                filter_summary, contig_name
            ):
                continue

            barcode = vdj_utils.get_barcode_from_contig_name(contig_name)
            barcode_contigs[barcode].append((contig_name, contig_seq))

    # Compute metrics for each barcode
    if args.cell_barcodes:
        barcodes = [
            ensure_binary(bc) for bc in vdj_utils.load_cell_barcodes_json(args.cell_barcodes)
        ]
    else:
        # Pass an empty barcode JSON for bulk
        barcodes = []

    reference = vdj_ref.VdjReference(args.vdj_reference_path)

    for barcode in barcodes:  # type: bytes
        contigs = barcode_contigs[barcode]
        annotations = [contig_annotations[contig[0]] for contig in contigs]
        reporter.vdj_barcode_contig_cb(barcode, contigs, annotations, reference, prefix)

        if not umi_summary is None and barcode in umi_summary.groups:
            bc_umi_summary = umi_summary.get_group(barcode)
        else:
            bc_umi_summary = None

        reporter.vdj_assembly_cb(bc_umi_summary, annotations, reference, prefix)

    reporter.report_summary_json(outs.summary)


def join(args, outs, chunk_defs, chunk_outs):
    shutil.copy(chunk_outs[0].summary, outs.summary)
