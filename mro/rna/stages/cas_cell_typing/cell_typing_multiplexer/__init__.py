#
# Copyright (c) 2025 10X Genomics, Inc. All rights reserved.
#

"""Multiplexing cell typing methods."""

import martian

import cellranger.cr_io as cr_io
from cellranger.cell_typing.common_cell_typing import CellTypeResults

__MRO__ = """
stage CELL_TYPING_MULTIPLEXER(
    in  bool            disable_cloud_cell_types,
    in  bool            disable_azimuth_cell_types,
    in  CellTypeResults cloud_cell_type_results,
    in  bool            cloud_cas_success,
    in  bool            cloud_disable_cas_ws,
    in  bool            cloud_disable_summarize,
    in  json            cloud_summary,
    in  CellTypeResults azimuth_cell_type_results,
    in  bool            azimuth_cas_success,
    in  bool            azimuth_disable_cas_ws,
    in  bool            azimuth_disable_summarize,
    in  json            azimuth_summary,
    out CellTypeResults cell_type_results,
    out bool            cas_success,
    out bool            disable_cas_ws,
    out bool            disable_summarize,
    out json            summary,
    src py              "stages/cas_cell_typing/cell_typing_multiplexer",
)
"""


def main(args, outs):
    if not args.disable_cloud_cell_types:
        outs.cas_success = args.cloud_cas_success
        outs.disable_cas_ws = args.cloud_disable_cas_ws
        outs.disable_summarize = args.cloud_disable_summarize
        if args.cloud_summary:
            cr_io.hardlink_with_fallback(args.cloud_summary, outs.summary)
        else:
            outs.summary = None

        if args.cloud_cell_type_results:
            cloud_cell_type_results = CellTypeResults(**args.cloud_cell_type_results)

            cell_type_results = propagate_cell_type_results(cloud_cell_type_results)
            outs.cell_type_results = cell_type_results.to_outs()
        else:
            outs.cell_type_results = None

    elif not args.disable_azimuth_cell_types:
        outs.cas_success = args.azimuth_cas_success
        outs.disable_cas_ws = args.azimuth_disable_cas_ws
        outs.disable_summarize = args.azimuth_disable_summarize
        if args.azimuth_summary:
            cr_io.hardlink_with_fallback(args.azimuth_summary, outs.summary)
        else:
            outs.summary = None

        if args.azimuth_cell_type_results:
            azimuth_cell_type_results = CellTypeResults(**args.azimuth_cell_type_results)

            cell_type_results = propagate_cell_type_results(azimuth_cell_type_results)
            outs.cell_type_results = cell_type_results.to_outs()
        else:
            outs.cell_type_results = None

    else:
        martian.clear(outs=outs)


def propagate_cell_type_results(cloud_cell_type_results) -> CellTypeResults:
    """Propagate cell type results to a new struct."""
    if cloud_cell_type_results.cell_types:
        cell_types = martian.make_path("cell_types.csv").decode()
        cr_io.hardlink_with_fallback(cloud_cell_type_results.cell_types, cell_types)
    else:
        cell_types = None
    if cloud_cell_type_results.results:
        results = martian.make_path("results.json.gz").decode()
        cr_io.hardlink_with_fallback(cloud_cell_type_results.results, results)
    else:
        results = None
    if cloud_cell_type_results.metadata:
        metadata = martian.make_path("metadata.json").decode()
        cr_io.hardlink_with_fallback(cloud_cell_type_results.metadata, metadata)
    else:
        metadata = None
    return CellTypeResults(
        cell_types=cell_types,
        results=results,
        metadata=metadata,
        frac_returned_bcs=cloud_cell_type_results.frac_returned_bcs,
        model_used=cloud_cell_type_results.model_used,
        skip_downstream=cloud_cell_type_results.skip_downstream,
    )
