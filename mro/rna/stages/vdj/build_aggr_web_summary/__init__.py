# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
"""Generate VDJ aggr web summary from the json."""

from __future__ import annotations

import json
from typing import TYPE_CHECKING

import martian

if TYPE_CHECKING:
    import cellranger.mro_types.filetypes as mro_filetypes

__MRO__ = """
stage BUILD_AGGR_WEB_SUMMARY(
    in  json content,
    in  json diversity_chart,
    out json web_summary_data,
    src py   "stages/vdj/build_aggr_web_summary",
)
"""


class StageInputs(martian.Record):  # pylint: disable=too-few-public-methods
    """BUILD_AGGR_WEB_SUMMARY stage inputs."""

    content: mro_filetypes.json
    diversity_chart: mro_filetypes.json


class StageOutputs(martian.Record):  # pylint: disable=too-few-public-methods
    """BUILD_AGGR_WEB_SUMMARY stage outputs."""

    web_summary_data: mro_filetypes.json


class ChunkInputs(StageInputs):  # pylint: disable=too-few-public-methods
    """BUILD_AGGR_WEB_SUMMARY chunk inputs."""


class ChunkOutputs(StageOutputs):  # pylint: disable=too-few-public-methods
    """BUILD_AGGR_WEB_SUMMARY chunk outputs."""


def main(args: ChunkInputs, outs: ChunkOutputs):
    with open(args.content) as f:
        summary = json.load(f)

    # NOT very pretty. Would like to move this to Rust soon
    with open(args.diversity_chart) as f:
        diversity_chart = json.load(f)
    summary["summary_tab"]["vdj_diversity"] = diversity_chart
    websummary_data = {"summary": summary}
    with open(outs.web_summary_data, "w") as f:
        json.dump(websummary_data, f, indent=4)
