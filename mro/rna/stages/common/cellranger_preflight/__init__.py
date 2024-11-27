#!/usr/bin/env python
#
# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
#
import martian

import cellranger.csv_utils as cr_csv_utils
import cellranger.preflight as cr_preflight
from cellranger.feature_ref import FeatureDefException

__MRO__ = """
stage CELLRANGER_PREFLIGHT(
    in  bool             full_check,
    in  string           chemistry,
    in  map[]            sample_def,
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


def run_preflight_checks(args):
    cr_preflight.check_os()

    cr_preflight.check_sample_info(
        args.sample_def, args.reference_path, args.full_check, args.feature_reference
    )

    cr_preflight.check_common_preflights(
        args.full_check,
        args.reference_path,
        args.r1_length,
        args.r2_length,
        args.recovered_cells,
        args.force_cells,
    )

    cr_preflight.check_feature_preflights(args.sample_def, args.feature_reference)
    cr_preflight.check_targeting_preflights(
        args.target_set,
        args.reference_path,
        args.feature_reference,
        parse_files=args.full_check,
        expected_targeting_method=args.targeting_method,
        is_spatial=False,
    )


def main(args, outs):
    try:
        run_preflight_checks(args)
    except (
        cr_preflight.PreflightException,
        cr_csv_utils.CSVParseException,
        FeatureDefException,
    ) as e:
        martian.exit(e.msg)

    cr_preflight.record_package_versions()
