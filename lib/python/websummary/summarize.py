# Copyright (c) 2020 10x Genomics, Inc. All rights reserved.


"""Basic python interface for creating a websummary.

If run as "python summarize.py example" it
creates an example HTML websummary using the template in
example/summary.html and the data in example/data.json.
"""

from __future__ import annotations

import json
import os
import re
import sys
from typing import Optional, TextIO, Union

# Hard-coded template file to be used by `generate_html_summary`
DIST_TEMPLATE: str = "dist/template.html"

# The name of the default font to use.  This is exported for use by callers
# who for whatever reason need to put it in other places (e.g. plotly layouts).
DEFAULT_FONT: str = "Overpass"

FILES_TO_INLINE: dict[str, str] = {
    "tenx-websummary-script.min.js": "dist/tenx-websummary-script.min.js",
    "tenx-websummary-styles.min.css": "dist/tenx-websummary-styles.min.css",
}

# The regex to scan for includes into the summary template.  Syntax is
# [[ include LOCAL_FILE_PATH ]] or
# [[ include LOCAL_FILE_PATH ifset DATAKEY ]] or
# [[ include LOCAL_FILE_PATH ifset DATAKEY else OTHER_FILE_PATH ]]
# File path can contain letters, digits, underscores, forward slashes, and dashes.
# File path is local relative to the input template_dir.
# Data key can contain letters, underscores, and dashes.
TEMPLATE_REGEX = (
    r"\[\[ include (?P<filename>{filepath_regex})"
    r"( ifset (?P<key>{datakey_regex})"
    r"( else (?P<else_filename>{filepath_regex}))?)? \]\]"
).format(filepath_regex=r"[a-zA-Z.\/_\d-]+", datakey_regex=r"[a-zA-Z_-]+")


def find_local_path(relative: str) -> str:
    """Compute path name of various known resource files.

    Paths are relative to the location of this python file.
    """
    return os.path.join(os.path.dirname(__file__), relative)


_SANITIZE_TRANSLATION: dict[int, str] = str.maketrans(
    {ord(c): rf"\u{ord(c):04x}" for c in ("<", ">", "&", "\u2028", "\u2029")}
)


def _sanitize_json(obj, **kwargs) -> str:
    """Serializes an object to json that is safe to embed in an html document.

    When a string is serialized to json and then embedded in an html document,
    as we do here, certain characters can cause problems.

    The issue this is trying to prevent is ending up with something like the
    following in the resulting html document.

    .. code-block:: html

        <script type="application/javascript>
            const data = {"description":"</script><script>evil();"}
        </script>

    The specific list of characters we sanitize here is the same as what golang
    escapes by default: `<`, `>`, `&`, U+2028 and U+2029 (line and paragraph
    separator characters).
    """
    return json.dumps(obj, separators=kwargs.pop("separators", (",", ":")), **kwargs).translate(
        _SANITIZE_TRANSLATION
    )


def generate_html_summary(
    json_data: dict[str, Union[str, int, float, dict, list, None]],
    summary_contents: str,
    template_dir: str,
    output_filehandle: TextIO,
    cls: Optional[type[json.JSONEncoder]] = None,
) -> None:
    """Generate an HTML summary.

    Args:
        json_data:          a json-dumpable object containing data to be injected into the summary.
        summary_contents:   the inner HTML for the web summary
        template_dir:       a path to the folder containing the templates
        output_filehandle:  output where the final HTML will be written
        cls:                a JSONEncoder object to be passed to json.dumps(..., cls=cls)
    """
    template_file = DIST_TEMPLATE
    if template_dir and os.path.isfile(
        find_local_path(os.path.join(template_dir, "template.html"))
    ):
        template_file = find_local_path(os.path.join(template_dir, "template.html"))
    with open(find_local_path(template_file)) as template_fh:
        template = template_fh.read()

    for search_string, filename in FILES_TO_INLINE.items():
        with open(find_local_path(filename)) as inline_fh:
            file_contents = inline_fh.read()
        template = template.replace(("[[ %s ]]" % search_string), file_contents)

    # Replace inline subtemplates recursively from inner html.
    # We limit the number of loops to avoid infinite recursions.
    count = 0
    while True:
        if count > 100:
            raise ValueError("Maximum template recursion depth exceeded")
        count += 1
        match = re.search(TEMPLATE_REGEX, summary_contents)
        if match is None:
            break
        filename = match.group("filename")
        key_query = json_data.get(match.group("key"), None)
        if key_query is not None and not key_query:
            # Key was required to be set but is False/None
            if match.group("else_filename") is None:
                # No fallback filename, don't import and just erase the match
                summary_contents = summary_contents.replace(match.group(), "")
                continue
            else:
                filename = match.group("else_filename")
        with open(os.path.join(template_dir, filename)) as infile:
            file_contents = infile.read()
        summary_contents = summary_contents.replace(match.group(), file_contents)

    template = template.replace("[[ data.js ]]", _sanitize_json(json_data, cls=cls))
    template = template.replace("[[ summary.html ]]", summary_contents)

    output_filehandle.write(template)


def main(argv: list[str], outfile: TextIO) -> None:
    """Executable entry point.

    Parses the command line, loads the data and content files, and writes
    the result to the given output.
    """
    json_file = "example/data.json"
    contents_file = "example/summary.html"
    if len(argv) > 1 and argv[1] == "cr":
        json_file = "tests/cr_tests/data/count_small.json"
        contents_file = "tests/cr_tests/summary.html"
    elif len(argv) > 1 and argv[1] == "cr_aggr_1":
        json_file = "cr_example/cr_aggr_1.json"
        contents_file = "cr_example/aggr_summary.html"
    elif len(argv) > 1 and argv[1] == "cr_aggr_2":
        json_file = "cr_example/cr_aggr_2.json"
        contents_file = "cr_example/aggr_summary.html"
    elif len(argv) > 1 and argv[1] == "cr_multi":
        json_file = "cr_example/multi.json"
        contents_file = "cr_example/multi_summary.html"
    elif len(argv) > 1 and argv[1] == "cs_example":
        json_file = "cs_example/data.json"
        contents_file = "cs_example/web_summary.html"
    elif len(argv) > 1 and argv[1] == "cmo_multi":
        json_file = "cr_example/cmo.json"
        contents_file = "cr_example/multi_summary.html"
    elif len(argv) > 1 and argv[1] == "hashtag_multi":
        json_file = "cr_example/hashtag.json"
        contents_file = "cr_example/multi_summary.html"
    elif len(argv) > 1 and argv[1] == "cas_multi":
        json_file = "cr_example/multi_cell_annotation.json"
        contents_file = "cr_example/multi_summary.html"
    elif len(argv) > 1 and argv[1] == "ab_multi":
        json_file = "cr_example/ab2.json"
        contents_file = "cr_example/multi_summary.html"
    elif len(argv) > 1 and argv[1] == "sr_example":
        json_file = "sr_example/data.json"
        contents_file = "sr_example/web_summary.html"
    elif len(argv) > 1 and argv[1] == "cr_aggr_ag":
        json_file = "cr_example/cr_aggr_ag.json"
        contents_file = "cr_example/aggr_summary.html"
    with open(find_local_path(json_file)) as infile:
        data = json.load(infile)
    with open(find_local_path(contents_file)) as infile:
        contents = infile.read()
    generate_html_summary(
        data, contents, os.path.join(os.path.dirname(__file__), "example"), outfile
    )


if __name__ == "__main__":
    main(sys.argv, sys.stdout)
