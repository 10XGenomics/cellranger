#!/usr/bin/env python
#
# Copyright (c) 2018 10x Genomics, Inc. All rights reserved.
#
"""Library functions for performing batch correction."""


from __future__ import annotations

import struct
from collections import Counter

import numpy as np
import sklearn.neighbors as sk_neighbors
from sklearn.metrics.pairwise import rbf_kernel

DEFAULT_BALLTREE_LEAFSIZE = 40


def batch_effect_score(
    dimred_matrix: np.ndarray[tuple[int, int], np.dtype[np.float64]],
    batch_ids: np.ndarray[int, np.dtype[np.bytes_]],
    knn_neighbors: int | None = None,
    knn_frac: float | None = 0.01,
    max_num_bcs: int | None = 10000,
):
    """Compute batch effect score on an aggregated dimension-reduced matrix.

    The batch effect score quantifies the degree of separation between groups of barcodes (batches).
    For each barcode, the fraction of its k-nearest neighbors with the same batch id is calculated.
    Local batch scores are computed by scaling and shifting the same-batch fractions to be strictly
    less than the number of batches and have expectation equal to 1 under identical batches.
    Finally, the average of these local batch scores is computed to obtain the overall batch effect
    score. A batch score of 1 indicates no separation between batches, and a batch score equal to
    the number of batches indicates perfect separation. A batch score less than 1 may occur in rare
    cases and is consistent with no batch separation.

    Args:
        dimred_matrix (np.ndarray[tuple[int, int], np.dtype[np.float64]]): dimension reduced matrix.
            Rows are samples and columns are features (e.g. principal components)
        batch_ids (np.ndarray[int, np.dtype[np.string_]]): 1-D array of batch IDs. Must have at
            least 2 observations per batch
        knn_frac (Optional[float], optional): Sets k in the k-nearest neighbors calculation to a
            fraction of the total barcodes (after subsampling to max_num_bcs). Defaults to 0.01.
        max_num_bcs (Optional[int], optional): Maximum number of barcodes to use. Larger data sets
            will be subsampled. Defaults to 10000.
        knn_neighbors (Optional[int], optional): Use a fixed k rather than knn_frac
    """
    if knn_neighbors is None and knn_frac is None:
        raise ValueError("One of knn_neighbors or knn_frac must be specified")

    num_bcs = dimred_matrix.shape[0]
    if num_bcs != len(batch_ids):
        raise ValueError("Length of batch_ids must equal number of rows in dimred_matrix")

    batch_counts_orig = Counter(batch_ids)

    # subsample barcodes if greater than the specified max_num_bcs
    if max_num_bcs is not None and num_bcs > max_num_bcs:
        np.random.seed(0)
        select_bc_idx = np.random.choice(num_bcs, max_num_bcs)
        select_bc_idx.sort()
        dimred_matrix = dimred_matrix[select_bc_idx]
        batch_ids = batch_ids[select_bc_idx]
        num_bcs = dimred_matrix.shape[0]

    # Return NaN if we completely dropped any batch or if not enough data
    # NOTE: this can happen due to subsampling but only for severely imbalanced cases,
    #       e.g. 10,000x more cells in batch 1 vs 2)
    batch_counts = Counter(batch_ids)
    if len(batch_counts) != len(batch_counts_orig) or min(batch_counts.values()) < 2:
        return np.nan

    if knn_neighbors is not None:
        num_neighbors = knn_neighbors
    else:
        num_neighbors = int(np.ceil(knn_frac * num_bcs))

    # For a cell in a given batch, what fraction of other cells share the same batch?
    # This is the expectation of same_batch_frac below, given the null of identical batch mixing
    batch_counts = Counter(batch_ids)
    num_batches = len(batch_counts)
    batch_to_frac = {batch: (count - 1) / (num_bcs - 1) for batch, count in batch_counts.items()}
    null_same_batch_frac = np.fromiter((batch_to_frac[i] for i in batch_ids), dtype=np.float64)
    # What is largest that same_batch_frac can be? Correction only relevant for very small batches
    # (e.g. 10k total cells, 100 or fewer cells in a single batch)
    batch_to_max_frac = {
        batch: min(count - 1, num_neighbors) / num_neighbors
        for batch, count in batch_counts.items()
    }
    max_same_batch_frac = np.fromiter((batch_to_max_frac[i] for i in batch_ids), dtype=np.float64)

    balltree = sk_neighbors.BallTree(dimred_matrix, leaf_size=DEFAULT_BALLTREE_LEAFSIZE)
    knn_idx = balltree.query(dimred_matrix, k=num_neighbors + 1, return_distance=False)

    same_batch_frac = np.mean(batch_ids[:, None] == batch_ids[knn_idx[:, 1:]], axis=1)
    # for identical batches, local_batch_score is equal to null_same_batch_frac in expectation
    # we rescale to be between 1 and num_batches
    local_batch_score = 1 + (num_batches - 1) * (same_batch_frac - null_same_batch_frac) / (
        max_same_batch_frac - null_same_batch_frac
    )

    return np.mean(local_batch_score)


def find_knn(curr_matrix, ref_matrix, knn):
    """For each row in curr_matrix, find k nearest neighbors in ref_matrix,.

    return an array of shape=[curr_matrix.shape[0] * knn, ], which stores
    the index of nearest neighbors in ref_matrix
    """
    balltree = sk_neighbors.BallTree(ref_matrix, leaf_size=DEFAULT_BALLTREE_LEAFSIZE)
    num_neighbors = min(ref_matrix.shape[0], knn)
    nn_idx = balltree.query(curr_matrix, k=num_neighbors, return_distance=False)
    return nn_idx.ravel().astype(int)


def serialize_batch_nearest_neighbor(fp, batch_nearest_neighbor):
    for (a, b), s in batch_nearest_neighbor.items():
        fp.write(struct.pack("qqQ", a, b, len(s)))
        for i, j in s:
            fp.write(struct.pack("qq", i, j))


def deserialize_batch_nearest_neighbor(fp):
    """>>> from cStringIO import StringIO.

    >>> batch1 = dict()
    >>> batch1[(0, 1)] = set([(1, 2), (3, 4), (5, 6)])
    >>> batch1[(1, 2)] = set([(7, 8), (9, 10)])
    >>> batch1[(3, 4)] = set([(11, 12)])
    >>> fp = StringIO()
    >>> serialize_batch_nearest_neighbor(fp, batch1)
    >>> fp.seek(0)
    >>> batch2 = deserialize_batch_nearest_neighbor(fp)
    >>> batch1 == batch2
    True
    """
    batch_nearest_neighbor = {}
    while True:
        fmt = "qqQ"
        sz = struct.calcsize("qqQ")
        buf = fp.read(sz)
        if len(buf) == 0:
            break
        elif len(buf) != sz:
            raise RuntimeError("corrupted batch_nearest_neighbor stream (key)")
        a, b, slen = struct.unpack(fmt, buf)
        fmt = "qq"
        sz = struct.calcsize("qq")
        s = set()
        for _ in range(slen):
            buf = fp.read(sz)
            if len(buf) != sz:
                raise RuntimeError("corrupted batch_nearest_neighbor stream (set)")
            i, j = struct.unpack(fmt, buf)
            s.add((i, j))
        batch_nearest_neighbor[(a, b)] = s
    return batch_nearest_neighbor


def correction_vector(dimred_matrix, cur_submatrix_idx, mnn_cur_idx, mnn_ref_idx, sigma):
    """Compute the batch-correction vector.

    1. For each MNN pair in current dataset and the reference, a pair-specific
    batch-correction vector is computed as the vector difference between the
    paired cells.
    2. For each barcode in cur dataset, a batch-correction vector is calculated
    as a weighted average of these pair-specific vectors, as computed with a
    Gaussian kernel.
    """
    num_pcs = dimred_matrix.shape[1]
    corr_vector = np.zeros((0, num_pcs))

    # the number of mnn and submatrix dim might be very large, process by chunk to save memory
    cur_submatrix_size = len(cur_submatrix_idx)
    mnn_size = len(mnn_cur_idx)
    # based on empirical testing
    cur_submatrix_chunk_size = int(1e6 / num_pcs)
    mnn_chunk_size = int(2e7 / num_pcs)

    for i in range(0, cur_submatrix_size, cur_submatrix_chunk_size):
        cur_submatrix_chunk = cur_submatrix_idx[i : i + cur_submatrix_chunk_size]
        cur_submatrix = dimred_matrix[cur_submatrix_chunk]

        weighted_sum, weights_sum = np.zeros(cur_submatrix.shape), np.zeros(cur_submatrix.shape)

        for j in range(0, mnn_size, mnn_chunk_size):
            mnn_cur_chunk = mnn_cur_idx[j : j + mnn_chunk_size]
            mnn_ref_chunk = mnn_ref_idx[j : j + mnn_chunk_size]

            mnn_cur = dimred_matrix[mnn_cur_chunk]
            weights = rbf_kernel(cur_submatrix, mnn_cur, gamma=0.5 * sigma)
            bias = dimred_matrix[mnn_ref_chunk] - mnn_cur
            weighted_sum += np.dot(weights, bias)
            weights_sum += np.tile(np.sum(weights, axis=1), (num_pcs, 1)).T

        corr_vector = np.vstack((corr_vector, weighted_sum / weights_sum))

    return corr_vector
