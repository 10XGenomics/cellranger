# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
"""Parses the target_panel file that contains IDs of target genes.

Outputs include the
indices of target genes within the feature reference, panel metadata, and boolean
flags for disabling targeting-related stages downstream.
"""


import json
import shutil
from collections import OrderedDict

import martian

import cellranger.csv_io as cr_csv_io
import cellranger.reference as cr_reference
import cellranger.rna.library as rna_library
import cellranger.sample_def as cr_sample_def
from cellranger.targeted import simple_utils, targeted_constants

__MRO__ = """
stage PARSE_TARGET_FEATURES(
    in  map[]    sample_def,
    in  path     reference_path,
    in  json     gene_index,
    in  bool     filter_probes,
    in  bool     no_target_umi_filter,
    in  bool     no_bam,
    in  bool     is_pd,
    in  int      rps_limit,
    out fa       bait_fasta,
    out csv      target_panel,
    out csv      probe_set,
    out csv      target_panel_or_probe_set,
    out csv      target_gene_indices,
    out bool     disable_targeted,
    out bool     disable_target_umi_filter,
    out bool     no_bam,
    out int      rps_limit,
    out string   target_set_name,
    out string   targeting_method,
    out tps.json target_panel_summary,
    src py       "stages/common/parse_target_features",
) using (
    mem_gb = 4,
)
"""


def main(args, outs):
    # HACK: needed for backwards compatibility against old sample defs. ultimately this stage
    #       needs to be within count_gem_well_processor or the new (rustified) setup_chunks,
    #       which will fix this.
    for sample_def in args.sample_def:
        if cr_sample_def.get_library_type(sample_def) is None:
            sample_def[rna_library.LIBRARY_TYPE] = rna_library.GENE_EXPRESSION_LIBRARY_TYPE

    # is there any target set?
    target_sets = [
        (cr_sample_def.get_target_set_name(sd), cr_sample_def.get_target_set(sd))
        for sd in args.sample_def
        if cr_sample_def.get_target_set(sd) is not None
    ]
    target_sets = sorted(set(target_sets))

    if len(target_sets) > 1:
        martian.exit(f"Multiple target sets/target set names found:\n\t{target_sets}")
    elif len(target_sets) == 0:
        # Sample is not targeted
        ## Set relevant flags
        outs.no_bam = args.no_bam
        outs.disable_targeted = True
        outs.disable_target_umi_filter = True

        ## These will all be null if not targeted data.
        ## (Also, don't rely on deprecated martian behavior that initializes paths)
        outs.rps_limit = None
        outs.target_set_name = None
        outs.targeting_method = None
        outs.target_panel = None
        outs.probe_set = None
        outs.target_panel_or_probe_set = None
        outs.target_gene_indices = None
        outs.bait_fasta = None
        outs.target_panel_summary = None
    else:
        # Sample is targeted, gather basic info first
        target_set_name, target_set_fn = target_sets[0]

        ## load reference genes
        gene_index = cr_reference.NewGeneIndex.load_from_json(args.gene_index)
        gene_name_to_id = {gene.name: gene.id for gene in gene_index.genes}

        ## parse target set
        (target_set_metadata, target_gene_ids, bait_sequences) = simple_utils.parse_target_csv(
            target_set_fn, gene_index, gene_name_to_id, args.filter_probes
        )

        target_set_name = target_set_metadata.get("panel_name", target_set_name)
        target_set_type = target_set_metadata.get("panel_type", None)
        targeting_method = simple_utils.determine_targeting_type_from_csv(target_set_fn)

        ## get gene indices
        target_gene_indices = sorted(
            filter(
                lambda x: x is not None,
                (gene_index.gene_id_to_int(gene_id) for gene_id in target_gene_ids),
            )
        )
        assert len(target_gene_indices) > 0

        # Set up outs
        outs.no_bam = args.no_bam
        outs.disable_targeted = False
        outs.rps_limit = args.rps_limit
        outs.target_set_name = target_set_name
        outs.targeting_method = targeting_method

        if targeting_method == targeted_constants.TARGETING_METHOD_HC:
            outs.target_panel = martian.make_path("target_panel.csv")
            outs.probe_set = None
            outs.target_panel_or_probe_set = outs.target_panel
            outs.disable_target_umi_filter = args.no_target_umi_filter
            # The stages CALCULATE_TARGETED_METRICS_PD and TARGETED_WEBSUMMARY_PD
            # in the pipeline _TARGETED_ANALYZER_PD require a BAM file.
            if args.is_pd:
                outs.no_bam = False
        elif targeting_method == targeted_constants.TARGETING_METHOD_TL:
            outs.target_panel = None
            outs.probe_set = martian.make_path("probe_set.csv")
            outs.target_panel_or_probe_set = outs.probe_set
            outs.disable_target_umi_filter = True
        else:
            raise f"Unexpected targeting_method: {targeting_method}"

        ## pass through target panel file
        ## (FIXME: after moving target_set to top-level arg this isn't needed)
        shutil.copyfile(target_set_fn, outs.target_panel_or_probe_set)

        ## gene indices
        cr_csv_io.write_target_features_csv(outs.target_gene_indices, target_gene_indices)

        ## Write the bait sequence FASTA file
        simple_utils.write_bait_fasta(outs.bait_fasta, bait_sequences)

        ## Pass along some details about the target panel file in summary file
        target_panel_hash = cr_reference.compute_hash_of_file(target_set_fn)
        # TODO: This should be a type
        summary_dict = OrderedDict(
            [
                ("target_panel_hash", target_panel_hash),
                ("target_panel_name", target_set_name),
                ("target_panel_path", target_set_fn),
                ("target_panel_gene_count", len(target_gene_indices)),
                ("target_panel_type", target_set_type),
                ("targeting_method", targeting_method),
            ]
        )

        with open(outs.target_panel_summary, "w") as output:
            output.write(json.dumps(summary_dict, indent=4))
