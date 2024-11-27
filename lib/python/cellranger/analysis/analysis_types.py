# Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
"""Type definitions for various analysis types.

This allows others to use them without taking a dependency on
the actual analysis code.
"""

from __future__ import annotations

import csv
import itertools
from dataclasses import dataclass
from typing import TYPE_CHECKING, NamedTuple

import numpy as np
from six import ensure_str
from typing_extensions import Self

if TYPE_CHECKING:
    import cellranger.matrix as cr_matrix

# pylint: disable=invalid-name

# Diffexp table keys
DIFFEXP_TABLE_FEATURE_ID_KEY = "Feature ID"
DIFFEXP_TABLE_FEATURE_NAME_KEY = "Feature Name"


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

    @classmethod
    def read_diffexp_from_csv(cls, csv_path: str | bytes) -> Self:
        """Read differential expression from a CSV."""
        with open(csv_path) as f:
            reader = csv.DictReader(f)
            field_names = reader.fieldnames
            if (
                not field_names
                or len(field_names) < 2
                or DIFFEXP_TABLE_FEATURE_ID_KEY != field_names[0]
                or DIFFEXP_TABLE_FEATURE_NAME_KEY != field_names[1]
            ):
                raise ValueError(
                    f"CSV did not begin with fields {DIFFEXP_TABLE_FEATURE_ID_KEY} and {DIFFEXP_TABLE_FEATURE_NAME_KEY}"
                )

        array_list_of_lists = []
        with open(csv_path) as f:
            csvreader = csv.reader(f)
            next(csvreader)
            for row in csvreader:
                array_list_of_lists.append(list(map(float, row[2:])))

        return cls(data=np.array(array_list_of_lists))


@dataclass
class DifferentialExpressionWithFeatures:
    """Class with differential expression and feature stuff."""

    diffexp: DifferentialExpression
    feature_names: list[str]
    feature_ids: list[str]

    @classmethod
    def from_cmatrix_and_diffexp(
        cls, cmatrix: cr_matrix.CountMatrix, diffexp: DifferentialExpression
    ) -> Self:
        """Get differential expression with features from matrix and diffexp."""
        feature_names = cmatrix.feature_ref.get_feature_names()
        feature_ids = [ensure_str(x) for x in cmatrix.feature_ref.get_feature_ids()]
        return cls(feature_names=feature_names, feature_ids=feature_ids, diffexp=diffexp)

    @classmethod
    def from_diffexp_csv(cls, csv_path: str | bytes) -> Self:
        """Get differential expression with features from matrix and diffexp."""
        diffexp = DifferentialExpression.read_diffexp_from_csv(csv_path=csv_path)
        with open(csv_path) as f:
            feature_ids, feature_names = itertools.tee(
                (ensure_str(row[DIFFEXP_TABLE_FEATURE_ID_KEY]), row[DIFFEXP_TABLE_FEATURE_NAME_KEY])
                for row in csv.DictReader(f)
            )
            feature_ids, feature_names = [x[0] for x in feature_ids], [x[1] for x in feature_names]
        return cls(feature_names=feature_names, feature_ids=feature_ids, diffexp=diffexp)


class PLSA(NamedTuple):
    """Hold information related to PLSA dimensionality reduction."""

    transformed_plsa_matrix: np.ndarray[tuple[int, int], np.dtype[np.float64]]
    components: np.ndarray[tuple[int, int], np.dtype[np.float64]]
    variance_explained: np.ndarray[tuple[int, int], np.dtype[np.float64]]
    dispersion: np.ndarray[tuple[int, int], np.dtype[np.float64]]
    features_selected: np.ndarray[tuple[int, int], np.dtype[np.float64]]
