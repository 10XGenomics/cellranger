from __future__ import annotations

import os

import pandas as pd
import scipy.io as sp_io
from six import ensure_binary, ensure_str

from cellranger.feature_ref import FeatureDef, FeatureReference
from cellranger.library_constants import GENE_EXPRESSION_LIBRARY_TYPE
from cellranger.matrix import FEATURES_TSV_GZ, CountMatrix

GENES_TSV = "genes.tsv"


def load_mtx(mtx_dir):
    legacy_fn = os.path.join(mtx_dir, GENES_TSV)
    v3_fn = os.path.join(mtx_dir, FEATURES_TSV_GZ)
    if os.path.exists(legacy_fn):
        return from_legacy_mtx(mtx_dir)

    if os.path.exists(v3_fn):
        return from_v3_mtx(mtx_dir)

    raise OSError(f"Not a valid path to a feature-barcode mtx directory: '{mtx_dir!s}'")


def save_dense_csv(mat: CountMatrix, filename):
    """Save this matrix to a dense CSV file."""
    dense_cm = pd.DataFrame(
        mat.m.toarray(),
        index=[ensure_str(f.id) for f in mat.feature_ref.feature_defs],
        columns=[ensure_str(bc) for bc in mat.bcs],
    )
    dense_cm.to_csv(filename, index=True, header=True)


def from_legacy_mtx(genome_dir):
    barcodes_tsv = ensure_binary(os.path.join(genome_dir, "barcodes.tsv"))
    genes_tsv = ensure_binary(os.path.join(genome_dir, GENES_TSV))
    matrix_mtx = ensure_binary(os.path.join(genome_dir, "matrix.mtx"))
    for filepath in [barcodes_tsv, genes_tsv, matrix_mtx]:
        if not os.path.exists(filepath):
            raise OSError(f"Required file not found: {filepath}")
    barcodes = pd.read_csv(
        barcodes_tsv.encode(),
        delimiter="\t",
        header=None,
        usecols=[0],
        dtype=bytes,
        converters={0: ensure_binary},
    ).values.squeeze()
    genes = pd.read_csv(
        genes_tsv.encode(),
        delimiter="\t",
        header=None,
        usecols=[0],
        dtype=bytes,
        converters={0: ensure_binary},
    ).values.squeeze()
    feature_defs = [
        FeatureDef(idx, gene_id, None, GENE_EXPRESSION_LIBRARY_TYPE, {})
        for (idx, gene_id) in enumerate(genes)
    ]
    feature_ref = FeatureReference(feature_defs, [])

    matrix = sp_io.mmread(matrix_mtx)
    mat = CountMatrix(feature_ref, barcodes, matrix)
    return mat


def from_v3_mtx(genome_dir):
    barcodes_tsv = ensure_str(os.path.join(genome_dir, "barcodes.tsv.gz"))
    features_tsv = ensure_str(os.path.join(genome_dir, FEATURES_TSV_GZ))
    matrix_mtx = ensure_str(os.path.join(genome_dir, "matrix.mtx.gz"))
    for filepath in [barcodes_tsv, features_tsv, matrix_mtx]:
        if not os.path.exists(filepath):
            raise OSError(f"Required file not found: {filepath}")
    barcodes = pd.read_csv(
        barcodes_tsv, delimiter="\t", header=None, usecols=[0], dtype=bytes
    ).values.squeeze()
    features = pd.read_csv(features_tsv, delimiter="\t", header=None)

    feature_defs = []
    for idx, (_, r) in enumerate(features.iterrows()):
        fd = FeatureDef(idx, r[0], r[1], r[2], {})
        feature_defs.append(fd)

    feature_ref = FeatureReference(feature_defs, [])

    matrix = sp_io.mmread(matrix_mtx)
    mat = CountMatrix(feature_ref, barcodes, matrix)
    return mat
