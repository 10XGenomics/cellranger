# Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
"""Type definitions for various analysis types.

This allows others to use them without taking a dependency on
the actual analysis code.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, NamedTuple

if TYPE_CHECKING:
    import numpy as np

# pylint: disable=invalid-name


class PCA(NamedTuple):
    """Holds information related to PCA dimensionality reduction."""

    transformed_pca_matrix: np.ndarray[tuple[int, int], np.dtype[np.float64]]
    components: np.ndarray[tuple[int, int], np.dtype[np.float64]]
    variance_explained: np.ndarray[tuple[int, int], np.dtype[np.float64]]
    dispersion: np.ndarray[tuple[int, int], np.dtype[np.float64]]
    features_selected: np.ndarray[tuple[int, int], np.dtype[np.float64]]


class LSA(NamedTuple):
    """Hold information related to LSA dimensionality reduction."""

    transformed_lsa_matrix: np.ndarray[tuple[int, int], np.dtype[np.float64]]
    components: np.ndarray[tuple[int, int], np.dtype[np.float64]]
    variance_explained: np.ndarray[tuple[int, int], np.dtype[np.float64]]
    dispersion: np.ndarray[tuple[int, int], np.dtype[np.float64]]
    features_selected: np.ndarray[tuple[int, int], np.dtype[np.float64]]


class TSNE(NamedTuple):
    """Hold information related to TSNE projection."""

    transformed_tsne_matrix: np.ndarray[tuple[int, int], np.dtype[np.float64]]
    name: str  # Human readable form
    key: str  # Machine queryable form, must be unique


class UMAP(NamedTuple):
    """Hold information related to UMAP projection."""

    transformed_umap_matrix: np.ndarray[tuple[int, int], np.dtype[np.float64]]
    name: str  # Human readable form
    key: str  # Machine queryable form, must be unique


class DifferentialExpression(NamedTuple):
    """Hold information related to differential expression analysis."""

    data: np.ndarray[tuple[int, int], np.dtype[np.float64]]


class PLSA(NamedTuple):
    """Hold information related to PLSA dimensionality reduction."""

    transformed_plsa_matrix: np.ndarray[tuple[int, int], np.dtype[np.float64]]
    components: np.ndarray[tuple[int, int], np.dtype[np.float64]]
    variance_explained: np.ndarray[tuple[int, int], np.dtype[np.float64]]
    dispersion: np.ndarray[tuple[int, int], np.dtype[np.float64]]
    features_selected: np.ndarray[tuple[int, int], np.dtype[np.float64]]
