#
# Copyright (c) 2024 10X Genomics, Inc. All rights reserved.
#
"""Stage checking if cell annotation is viable but not requested in a multi run."""

import json

import martian

import cellranger.cell_typing.broad_tenx.cloud_cas_utils as cloud_utils
import cellranger.cell_typing.broad_tenx.common as ct_common
import cellranger.matrix as cr_matrix

__MRO__ = """
stage CLOUD_CELL_ANNOTATION_VIABLE_BUT_NOT_REQUESTED(
    in  h5     filtered_matrix,
    in  bool   skip_cell_annotation,
    in  string tenx_cloud_token_path,
    out bool   cloud_cell_annotation_viable_but_not_requested,
    out json   metrics,
    src py     "stages/cas_cell_typing/cloud_cell_annotation_viable_but_not_requested",
)
"""


def main(args, outs):
    if not args.filtered_matrix:
        martian.clear(outs)
        return

    ref = cr_matrix.CountMatrix.get_genomes_from_h5(args.filtered_matrix)
    # Make sure the genome is some version of GRCh38 or mm10 but there aren't more than 1 and that there is a GEX library
    valid_genome_library, _, _ = ct_common.check_valid_genome_and_library(
        matrix=args.filtered_matrix
    )
    try:
        token_info = cloud_utils.extract_token_info(args.tenx_cloud_token_path)
    except Exception as token_exc:  # pylint: disable=broad-except
        print(token_exc)
        token_info = cloud_utils.TokenInfo(successful=False)
    if (
        len(ref) == 1
        and valid_genome_library
        and not token_info.successful
        and not args.skip_cell_annotation
    ):
        outs.cloud_cell_annotation_viable_but_not_requested = True
    else:
        outs.cloud_cell_annotation_viable_but_not_requested = False

    metrics = {
        "cell_annotation_viable_but_not_requested": outs.cloud_cell_annotation_viable_but_not_requested,
    }
    with open(outs.metrics, "w") as outfile:
        json.dump(metrics, outfile)
