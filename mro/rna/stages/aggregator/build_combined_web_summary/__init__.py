# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
"""Generate aggr web summary from the individual json data."""

__MRO__ = """
stage BUILD_COMBINED_WEB_SUMMARY(
    in  json vdj_t_data,
    in  json vdj_b_data,
    in  json vdj_t_gd_data,
    in  json count_data,
    in  json antigen_data,
    out html web_summary,
    out json web_summary_data,
    src py   "stages/aggregator/build_combined_web_summary",
)
"""

import json

from websummary import summarize

ALARM_KEY = "alarms"
SAMPLE_KEY = "sample"
TOP_LEVEL_KEY = "summary"
ANALYSIS_TAB_KEY = "analysis_tab"
ANTIGEN_TAB_KEY = "antigen_tab"
SUMMARY_TAB_KEY = "summary_tab"
VDJ_T_TAB_KEY = "vdj_t_tab"
VDJ_T_GD_TAB_KEY = "vdj_t_gd_tab"
VDJ_B_TAB_KEY = "vdj_b_tab"


def main(args, outs):
    inner_data = {}
    if args.count_data is not None:
        with open(args.count_data) as f:
            count_data = json.load(f)
        inner_data = count_data[TOP_LEVEL_KEY]

    for key, data_file in [
        (VDJ_T_TAB_KEY, args.vdj_t_data),
        (VDJ_B_TAB_KEY, args.vdj_b_data),
        (VDJ_T_GD_TAB_KEY, args.vdj_t_gd_data),
    ]:
        if data_file is not None:
            with open(data_file) as f:
                data = json.load(f)
            inner_data[SAMPLE_KEY] = data[TOP_LEVEL_KEY][SAMPLE_KEY]
            inner_data[key] = data[TOP_LEVEL_KEY][SUMMARY_TAB_KEY]

    if args.antigen_data is not None:
        assert args.count_data is not None
        with open(args.antigen_data) as f:
            antigen_data = json.load(f)
        inner_data[ANTIGEN_TAB_KEY] = antigen_data

    websummary_data = {
        TOP_LEVEL_KEY: inner_data,
    }

    with open(outs.web_summary_data, "w") as f:
        json.dump(websummary_data, f, indent=4)

    contents = """<div data-key="summary" data-component="AggrSummary">"""
    with open(outs.web_summary, "w") as outfile:
        summarize.generate_html_summary(
            websummary_data,
            contents,
            None,
            outfile,
        )
