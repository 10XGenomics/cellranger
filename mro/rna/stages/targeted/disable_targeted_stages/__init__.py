#!/usr/bin/env python
#
# Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
#

"""Figures out if gDNA stages should be run."""

import martian

import cellranger.csv_io as cr_csv_io

__MRO__ = """
stage DISABLE_TARGETED_STAGES(
    in  csv  probe_set,
    in  bool is_visium_hd,
    out bool disable_targeted_gdna,
    out bool disable_sampling_stages,
    src py   "stages/targeted/disable_targeted_stages",
)
"""
# also in Rust in cr_lib/stages/collate_probe_metrics.rs
GDNA_GENE_THRESHOLD = 10


def main(args, outs):
    """Determine if the gDNA analysis is disabled.

    Disabled if there are fewer than GDNA_GENE_THRESHOLD
    genes with both spliced and unspliced probes.

    Args:
        args (An object with one input):
            args.probe_set : path to probe set reference
        outs (An object with one utput):
            outs.disable_targeted_gdna: (bool) True if gDNA analysis should be disable
    """
    outs.disable_sampling_stages = bool(args.is_visium_hd)
    if (
        (args.probe_set is None)
        or ("region" not in cr_csv_io.load_csv_columnnames(args.probe_set))
        or ("gene_id" not in cr_csv_io.load_csv_columnnames(args.probe_set))
    ):
        outs.disable_targeted_gdna = True
        martian.log_info(
            "Skipping gDNA analysis because the probe set reference does not contain region annotations"
        )
    else:
        gene_ids = cr_csv_io.read_target_features_csv(args.probe_set, header="gene_id", comment="#")
        gene_regions = cr_csv_io.read_target_features_csv(
            args.probe_set, header="region", comment="#"
        )

        unspliced_genes = {x[0] for x in zip(gene_ids, gene_regions) if x[1] == "unspliced"}
        spliced_genes = {x[0] for x in zip(gene_ids, gene_regions) if x[1] == "spliced"}

        num_genes_with_spliced_unspliced_probes = len(unspliced_genes.intersection(spliced_genes))
        if num_genes_with_spliced_unspliced_probes < GDNA_GENE_THRESHOLD:
            martian.log_info(
                "Skipping gDNA analysis because the probe set reference has "
                + f"fewer than {GDNA_GENE_THRESHOLD} genes with both spliced and"
                + "unspliced probes."
            )
        outs.disable_targeted_gdna = num_genes_with_spliced_unspliced_probes < GDNA_GENE_THRESHOLD
