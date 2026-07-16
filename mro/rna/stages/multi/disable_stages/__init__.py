#!/usr/bin/env python
#
# Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
#
"""A helper stage to determine whether or not to disable the legacy bam file (holding all reads)."""

from __future__ import annotations

from typing import TYPE_CHECKING

import martian

if TYPE_CHECKING:
    import cellranger.mro_types.filetypes as mro_filetypes

__MRO__ = """
stage DISABLE_STAGES(
    in  bool  no_bam,
    in  bool  disable_multi,
    in  bool  is_visium_hd,
    in  csf[] probe_barcode_counts,
    out bool  disable_legacy_bam,
    out bool  disable_sample_bams,
    out bool  disable_subsampling,
    out bool  no_probe_barcode_counts,
    out bool  no_probe_barcode_matrix_demux,
    src py    "stages/multi/disable_stages",
) using (
    volatile = strict,
)
"""


class StageInputs(martian.Record):  # pylint: disable=too-few-public-methods
    """DISABLE_STAGES stage inputs."""

    no_bam: bool
    disable_multi: bool
    is_visium_hd: None | bool
    probe_barcode_counts: None | list[mro_filetypes.csf]


class StageOutputs(martian.Record):  # pylint: disable=too-few-public-methods
    """DISABLE_STAGES stage outputs."""

    disable_legacy_bam: bool
    disable_sample_bams: bool
    disable_subsampling: bool
    no_probe_barcode_counts: bool
    no_probe_barcode_matrix_demux: bool


class ChunkInputs(StageInputs):  # pylint: disable=too-few-public-methods
    """DISABLE_STAGES chunk inputs."""


class ChunkOutputs(StageOutputs):  # pylint: disable=too-few-public-methods
    """DISABLE_STAGES chunk outputs."""


def main(args: ChunkInputs, outs: ChunkOutputs):
    # Don't make an uber bam in multi
    outs.disable_legacy_bam = args.no_bam or not args.disable_multi

    # disable the new sample bam if no_bam was specified, or if it's not a multi run
    outs.disable_sample_bams = args.no_bam or args.disable_multi

    # Disable stages that require probe-barcode counts.
    outs.no_probe_barcode_counts = args.probe_barcode_counts is None

    # Disable demux probe-bc matrix if we are not in multi or have no probe barcode matrix
    outs.no_probe_barcode_matrix_demux = outs.no_probe_barcode_counts or args.disable_multi

    # disable subsampling iff this is a visium HD run
    outs.disable_subsampling = bool(args.is_visium_hd)
