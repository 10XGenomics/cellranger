# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
"""Generate VDJ aggr web summary from the json."""

__MRO__ = """
stage BUILD_AGGR_WEB_SUMMARY(
    in  json content,
    in  json diversity_chart,
    out json web_summary_data,
    src py   "stages/vdj/build_aggr_web_summary",
)
"""

import json


def main(args, outs):
    with open(args.content) as f:
        summary = json.load(f)

    # NOT very pretty. Would like to move this to Rust soon
    with open(args.diversity_chart) as f:
        diversity_chart = json.load(f)
    summary["summary_tab"]["vdj_diversity"] = diversity_chart
    websummary_data = {"summary": summary}
    with open(outs.web_summary_data, "w") as f:
        json.dump(websummary_data, f, indent=4)
