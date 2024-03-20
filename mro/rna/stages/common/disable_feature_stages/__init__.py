# Copyright (c) 2019 10x Genomics, Inc. All rights reserved.
"""Looks at the sample def and determines which feature-counter calls can be disabled."""

import cellranger.rna.library as rna_library
from cellranger.fast_utils import MultiGraph

__MRO__ = """
stage DISABLE_FEATURE_STAGES(
    in  map[]               sample_def,
    in  bool                disable_multi,
    in  bool                disable_count,
    in  bool                is_pd,
    in  bool                in_disable_targeted,
    in  map<SampleSlfeOuts> sample_outs,
    in  json                multi_graph,
    out bool                disable_crispr,
    out bool                disable_antibody,
    out bool                disable_antigen,
    out bool                disable_multiplexing,
    out bool                disable_targeted,
    out bool                disable_legacy_stages,
    out bool                disable_library_cloupe,
    out bool                disable_gex,
    src py                  "stages/common/disable_feature_stages",
)
"""


def main(args, outs):
    # this is a workaround so that this stage doesn't have to be disabled
    # instead we just instantly return reasonable outputs
    # if this stage can be disabled it causes disabled binding problems downstream
    if args.sample_def is None:
        outs.disable_crispr = True
        outs.disable_antibody = True
        outs.disable_antigen = True
        outs.disable_multiplexing = True
        outs.disable_legacy_stages = True
        outs.disable_targeted = True
        outs.disable_library_cloupe = True
        outs.disable_gex = True
        return

    sample_def = args.sample_def
    library_types = [x.get("library_type") for x in sample_def if x.get("library_type") is not None]

    found_crispr = rna_library.CRISPR_LIBRARY_TYPE in library_types
    found_antibody = rna_library.ANTIBODY_LIBRARY_TYPE in library_types
    found_antigen = rna_library.ANTIGEN_LIBRARY_TYPE in library_types
    found_multiplexing = rna_library.MULTIPLEXING_LIBRARY_TYPE in library_types
    found_gex = rna_library.GENE_EXPRESSION_LIBRARY_TYPE in library_types

    outs.disable_crispr = not (found_crispr)
    outs.disable_antibody = not (found_antibody)
    outs.disable_antigen = not (found_antigen)
    outs.disable_multiplexing = not (found_multiplexing)
    outs.disable_gex = not (found_gex)
    outs.disable_targeted = args.in_disable_targeted

    outs.disable_legacy_stages = (not args.is_pd) and (
        (not args.disable_multi) or args.disable_count
    )

    # Read in the multi graph
    is_multiplexed = False
    if args.multi_graph:
        is_multiplexed = MultiGraph.from_path(args.multi_graph).is_multiplexed()

    # for CMO-tagged runs of multi we must run the analyzer on the whole library
    # for single sample non-muxed this is the same as the sample analyzer so we should avoid reruns
    outs.disable_library_cloupe = (
        (args.multi_graph is not None)  # is a multi run
        and (not is_multiplexed)  # is not multiplexed (includes single sample tagged)
        and (args.sample_outs is not None)
        and (len(args.sample_outs) == 1)  # is single sample
    ) or args.disable_count
