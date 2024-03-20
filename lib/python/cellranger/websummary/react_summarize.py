#
# Copyright (c) 2023 10X Genomics, Inc. All rights reserved.
#
"""Write out the html, depends on websummary."""
from __future__ import annotations

from cellranger.websummary.react_components import ReactComponentEncoder, WebSummaryData
from websummary import summarize


def write_html_file(filename: str | bytes, websummary_data: WebSummaryData, contents=None):
    """Write out a web summary to an HTML file.

    Args:
        filename: string of output filename
        websummary_data: an instance of WebSummaryData
        contents: A string with the component to render
    """
    isinstance(websummary_data, WebSummaryData)
    if not contents:
        contents = """<div data-key="summary" data-component="CellRangerSummary">"""
    with open(filename, "w") as outfile:
        summarize.generate_html_summary(
            websummary_data,
            contents,
            None,
            outfile,
            cls=ReactComponentEncoder,
        )
