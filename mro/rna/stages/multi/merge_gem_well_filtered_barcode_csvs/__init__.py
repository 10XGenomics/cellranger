#!/usr/bin/env python
#
# Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
#
"""Merges the filtered barcodes csv and barcode correction CSV for multiple gem wells."""
import shutil

__MRO__ = """
stage MERGE_GEM_WELL_CSVS (
    in csv[] filtered_barcodes,
    in csv[] barcode_correction_csv,
    out csv  filtered_barcodes,
    out csv  barcode_correction_csv,
    src py   "stages/multi/merge_gem_well_filtered_barcode_csvs",
) using (
    volatile = strict,
)
"""


def main(args, outs):
    assert len(args.filtered_barcodes) == len(args.barcode_correction_csv)

    # merge the filtered barcodes CSV for multiple GEM wells
    with open(outs.filtered_barcodes, "w") as outfile:
        for unmerged_filtered_barcodes in args.filtered_barcodes:
            with open(unmerged_filtered_barcodes) as infile:
                shutil.copyfileobj(infile, outfile)

    # merge the barcode correction CSV for multiple GEM wells
    with open(outs.barcode_correction_csv, "w") as outfile:
        first = True
        for unmerged_barcode_correction_csv in args.barcode_correction_csv:
            with open(unmerged_barcode_correction_csv) as infile:
                header = infile.readline()
                assert (
                    header.strip()
                    == "library_type,barcode,reads,umis,candidate_dup_reads,umi_corrected_reads"
                )
                if first:
                    print(header, end="", file=outfile)
                    first = False

                shutil.copyfileobj(infile, outfile)
