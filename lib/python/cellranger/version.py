#
# Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
#

"""Utility for getting the version of the pipeline."""


from __future__ import annotations

import os.path

import tenkit.log_subprocess as tk_subproc


def get_version():
    """Get the version for the pipeline.

    Normally it will read this from the .version file in the package root, but
    if that is not found it will attempt to query git for the version.
    """
    # NOTE: this makes assumptions about the directory structure
    repo_base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
    script_dir = os.path.join(repo_base_dir, "bin")
    version_fn = os.path.join(repo_base_dir, ".version")

    if os.path.exists(version_fn):
        with open(version_fn) as f:
            output = f.read()
    else:
        try:
            output = tk_subproc.check_output(
                ["git", "describe", "--tags", "--always", "--dirty"], cwd=script_dir
            )
        except OSError:
            output = "unknown_version"
    return output.strip()
