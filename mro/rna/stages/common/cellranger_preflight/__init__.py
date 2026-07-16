#!/usr/bin/env python
#
# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
#

from __future__ import annotations

from typing import TYPE_CHECKING

import martian

import cellranger.csv_utils as cr_csv_utils
import cellranger.preflight as cr_preflight
from cellranger.feature_ref import FeatureDefException
from cellranger.sample_def import SampleDef

if TYPE_CHECKING:
    import cellranger.mro_types.filetypes as mro_filetypes
    import cellranger.mro_types.structs as mro_structs

__MRO__ = """
stage CELLRANGER_PREFLIGHT(
    in  bool             full_check,
    in  string           chemistry,
    in  SampleDef[]      sample_def,
    in  csv              target_set,
    in  path             reference_path,
    in  csv              feature_reference,
    in  CellCallingParam recovered_cells,
    in  CellCallingParam force_cells,
    in  int              r1_length,
    in  int              r2_length,
    in  string           targeting_method,
    src py               "stages/common/cellranger_preflight",
) using (
    mem_gb   = 8,
    volatile = strict,
)
"""


class StageInputs(martian.Record):  # pylint: disable=too-few-public-methods
    """CELLRANGER_PREFLIGHT stage inputs."""

    full_check: bool
    chemistry: str
    sample_def: list[mro_structs.SampleDef]
    target_set: mro_filetypes.csv
    reference_path: mro_filetypes.path
    feature_reference: mro_filetypes.csv
    recovered_cells: mro_structs.CellCallingParam
    force_cells: mro_structs.CellCallingParam
    r1_length: int
    r2_length: int
    targeting_method: str


class StageOutputs(martian.Record):  # pylint: disable=too-few-public-methods
    """CELLRANGER_PREFLIGHT stage outputs."""


class ChunkInputs(StageInputs):  # pylint: disable=too-few-public-methods
    """CELLRANGER_PREFLIGHT chunk inputs."""


class ChunkOutputs(StageOutputs):  # pylint: disable=too-few-public-methods
    """CELLRANGER_PREFLIGHT chunk outputs."""


def run_preflight_checks(args: ChunkInputs):
    cr_preflight.check_os()

    sdefs = [SampleDef.from_dict(sdef) for sdef in args.sample_def]

    cr_preflight.check_sample_info(
        sdefs, args.reference_path, args.full_check, args.feature_reference
    )

    cr_preflight.check_common_preflights(
        args.full_check,
        args.reference_path,
        args.r1_length,
        args.r2_length,
        args.recovered_cells,
        args.force_cells,
    )

    cr_preflight.check_feature_preflights(sdefs, args.feature_reference)
    cr_preflight.check_targeting_preflights(
        args.target_set,
        args.reference_path,
        args.feature_reference,
        parse_files=args.full_check,
        expected_targeting_method=args.targeting_method,
        is_spatial=False,
    )


def main(args: ChunkInputs, _outs: ChunkOutputs):
    try:
        run_preflight_checks(args)
    except (
        cr_preflight.PreflightException,
        cr_csv_utils.CSVParseException,
        FeatureDefException,
    ) as e:
        martian.exit(e.msg)

    cr_preflight.record_package_versions()
