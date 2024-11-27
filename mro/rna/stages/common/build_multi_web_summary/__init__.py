# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
"""Generate aggr web summary from the individual json data."""

__MRO__ = """
stage BUILD_MULTI_WEB_SUMMARY(
    in  map<json> web_summary_data,
    in  map<csv>  metrics_summary_csvs,
    out map<html> web_summaries,
    out map<csv>  metrics_summary_csvs,
    src py        "stages/common/build_multi_web_summary",
)
"""

import json
import shutil

import martian

from websummary import summarize


def main(args, outs):
    if args.web_summary_data is None:
        outs.web_summaries = None
        return

    outs.web_summaries = {}
    for sample, web_summary_data in args.web_summary_data.items():
        with open(web_summary_data) as f:
            web_summary_data = json.load(f)

        contents = """<div data-component="OptionalMultiplexedSummary">"""

        web_summary_path = martian.make_path(sample + "_web_summary.html")
        with open(web_summary_path, "w") as outfile:
            summarize.generate_html_summary(
                web_summary_data,
                contents,
                None,
                outfile,
            )

        outs.web_summaries[sample] = web_summary_path

    # need to copy metrics summaries into typed map so they're map-callable
    outs.metrics_summary_csvs = {}
    for sample, metrics in args.metrics_summary_csvs.items():
        path = martian.make_path(sample + "_metrics_summary.csv")
        shutil.copyfile(metrics, path)
        outs.metrics_summary_csvs[sample] = path
