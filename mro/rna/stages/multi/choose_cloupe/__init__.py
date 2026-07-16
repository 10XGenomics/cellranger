#!/usr/bin/env python
#
# Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
#
"""Determine if the cloupe was not generated because this is a single sample run of MULTI.

In this case we hard-link the single sample cloupe as the library-level cloupe.
"""

from __future__ import annotations

from shutil import copyfile
from typing import TYPE_CHECKING

import martian

import cellranger.cr_io as cr_io

if TYPE_CHECKING:
    import cellranger.mro_types.filetypes as mro_filetypes

__MRO__ = """
stage CHOOSE_CLOUPE(
    in  cloupe library_cloupe,
    in  map<cloupe> sample_cloupe,
    out cloupe cloupe,
    src py     "stages/multi/choose_cloupe",
) using (
    volatile = strict,
)
"""


class StageInputs(martian.Record):  # pylint: disable=too-few-public-methods
    """CHOOSE_CLOUPE stage inputs."""

    library_cloupe: None | mro_filetypes.cloupe
    sample_cloupe: None | dict[str, None | mro_filetypes.cloupe]


class StageOutputs(martian.Record):  # pylint: disable=too-few-public-methods
    """CHOOSE_CLOUPE stage outputs."""

    cloupe: None | mro_filetypes.cloupe


class ChunkInputs(StageInputs):  # pylint: disable=too-few-public-methods
    """CHOOSE_CLOUPE chunk inputs."""


class ChunkOutputs(StageOutputs):  # pylint: disable=too-few-public-methods
    """CHOOSE_CLOUPE chunk outputs."""


def main(args: ChunkInputs, outs: ChunkOutputs):
    if args.library_cloupe is None:
        if args.sample_cloupe is None:
            outs.cloupe = None
            return

        sample_cloupe = next(iter(args.sample_cloupe.values()))
        if sample_cloupe is None:
            outs.cloupe = None
            return
        # changed this from a hardlink to a copy to avoid issues with moving/deleting directories etc
        copyfile(sample_cloupe, outs.cloupe)

    else:
        cr_io.hardlink_with_fallback(args.library_cloupe, outs.cloupe)
