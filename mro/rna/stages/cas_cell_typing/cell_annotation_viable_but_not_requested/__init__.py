#
# Copyright (c) 2024 10X Genomics, Inc. All rights reserved.
#
"""Stage checking if cell annotation is viable but not requested in a multi run."""

import martian

import cellranger.cell_typing.common as ct_common
import cellranger.matrix as cr_matrix

__MRO__ = """
stage CELL_ANNOTATION_VIABLE_BUT_NOT_REQUESTED(
    in  h5   filtered_matrix,
    in  bool skip_cell_annotation,
    out bool cell_annotation_viable_but_not_requested,
    src py   "stages/cas_cell_typing/cell_annotation_viable_but_not_requested",
)
"""


def main(args, outs):
    if not args.filtered_matrix:
        martian.clear(outs)
        return

    viable_run = True
    ref = cr_matrix.CountMatrix.get_genomes_from_h5(args.filtered_matrix)
    if len(ref) > 1:
        viable_run = False

    # Make sure the genome is some version of GRCh38 or mm10 but there aren't more than 1 and that there is a GEX library
    valid_genome_library, _, _ = ct_common.check_valid_genome_and_library(
        matrix=args.filtered_matrix
    )
    if not valid_genome_library:
        viable_run = False

    if args.skip_cell_annotation:
        outs.cell_annotation_viable_but_not_requested = viable_run
        return
    else:
        outs.cell_annotation_viable_but_not_requested = False
        return
