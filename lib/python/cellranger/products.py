#
# Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
#

"""Methods for identifying and working with 10X product releases."""

from __future__ import annotations


def get_cmd_names(product_name: str) -> tuple[str, str]:
    """For a given product name (must be either cellranger or spaceranger) returns the human.

    (Cell Ranger) and the command line (cellranger) name.

    :param product_name: a string either spaceranger|cellranger|cellranger-arc|cellranger-atac
    :return: A tuple with the Human/CLI name, e.g. ("Space Ranger", "spaceranger")
    """
    if product_name == "spaceranger":
        return "Space Ranger", "spaceranger"
    elif product_name == "cellranger":
        return "Cell Ranger", "cellranger"
    elif product_name == "cellranger-arc":
        return "Cell Ranger Multiome ATAC + Gene Expression", "cellranger-arc"
    elif product_name == "cellranger-atac":
        return "Cell Ranger ATAC", "cellranger-atac"
    else:
        # will never happen, but if it does we need to catch in testing with an exception
        raise RuntimeError(
            f"Unrecognized product name: '{product_name}'. Expected 'spaceranger' or 'cellranger' or "
            "'cellranger-arc' or 'cellranger-atac'"
        )
