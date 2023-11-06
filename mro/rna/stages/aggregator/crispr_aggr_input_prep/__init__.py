#!/usr/bin/env python3
#
# Copyright (c) 2022 10X Genomics, Inc. All rights reserved
#
"""Prepare input for _CRISPR_ANALYZER stage in aggr."""

import cellranger.molecule_counter as cr_mc
import cellranger.rna.library as rna_library

__MRO__ = """
stage CRISPR_AGGR_INPUT_PREP(
    in  h5   merged_molecules,
    out csv  filtered_barcodes,
    out csv  feature_reference,
    out json counter_metrics_json,
    src py   "stages/aggregator/crispr_aggr_input_prep",
) using (
    mem_gb   = 4,
    volatile = strict,
)
"""


def main(args, outs):
    with cr_mc.MoleculeCounter.open(args.merged_molecules, "r") as mc:
        # Create feature_reference csv file
        feature_ref = mc.get_feature_ref()
        with open(outs.feature_reference, "w") as file_handle:
            feature_ref.select_features_by_type(rna_library.CRISPR_LIBRARY_TYPE).to_csv(file_handle)
        # filtered_barcode
        library_info = mc.get_library_info()
        barcode_info = mc.get_barcode_info()
        barcode_seqs = mc.get_barcodes()

        filtered_bcs = mc.get_filtered_barcodes(barcode_info, library_info, barcode_seqs)
        with open(outs.filtered_barcodes, "wb") as f_bc:
            for bc in filtered_bcs:
                f_bc.write(bc + b"\n")

    # dummy empty json file
    with open(outs.counter_metrics_json, "wb") as f_metrics:
        f_metrics.write(b"{}")
