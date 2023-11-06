# Copyright (c) 2023 10X Genomics, Inc. All rights reserved.
"""Tools for command line parsing."""
from __future__ import annotations

import os
import shlex


def preprocess_cmdline(cmdline: str) -> str:
    """Preprocess the command-line string to ensure correct argument formatting.

    Args:
        cmdline (str): The input command-line string.

    Returns:
        str: The preprocessed command-line string.
    """
    parts = cmdline.split(" ")
    preprocessed_parts = []

    i = 0
    while i < len(parts):
        part = parts[i]
        if part.startswith("--"):
            if i + 1 < len(parts) and not parts[i + 1].startswith("--"):
                # Check if the next part is not starting with "--"
                preprocessed_parts.append(part + "=" + parts[i + 1].rstrip("/"))
                i += 1  # Skip the next part
            else:
                preprocessed_parts.append(part)
        else:
            preprocessed_parts.append(part)
        i += 1

    return " ".join(preprocessed_parts)


def parse_cmdline_basename(cmdline: str) -> str:
    """Parse a command-line string while preserving basenames for paths with and without spaces.

    Args:
        cmdline (str): The input command-line string.

    Returns:
        str: A modified command-line string where basenames are preserved.
    """
    preprocessed_cmdline = preprocess_cmdline(cmdline)
    cmd_split = shlex.split(preprocessed_cmdline)
    result_list = []

    for item in cmd_split:
        parts = item.split("=")
        if len(parts) == 2:
            argument, path = parts
            basename = os.path.basename(path)
            result_list.append(f"{argument}={basename}")
        else:
            basename = os.path.basename(item)
            result_list.append(basename)
    result_string = " ".join(result_list)
    return result_string
