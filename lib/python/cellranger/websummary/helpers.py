#
# Copyright (c) 2023 10X Genomics, Inc. All rights reserved.
#
"""Small functions that can be used without adding dependencies if they were elsewhere."""
from __future__ import annotations


def get_projection_key(feature_type, n_components):
    """Return the correct tsne/umap key for each feature type and number of components, e.g. gene_expression_2."""
    return "{}_{}".format(feature_type.replace(" ", "_").lower(), n_components).encode()
