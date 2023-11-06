#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#
"""Louvain graph clustering interface."""
from __future__ import annotations

import errno
import os
import pickle as cPickle
import subprocess
import time
from collections.abc import Iterable

import h5py
import numpy as np
import scipy.sparse as sp_sparse
import sklearn.neighbors as sk_neighbors

import cellranger.analysis.clustering as cr_clustering
import cellranger.analysis.constants as analysis_constants
import cellranger.analysis.io as analysis_io
import tenkit.log_subprocess as tk_subproc
from cellranger.logperf import LogPerf
from cellranger.wrapped_tables import tables

# pylint: disable=invalid-name


LOUVAIN_CONVERT_BINPATH = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))), "bin", "convert"
)
LOUVAIN_BINPATH = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))), "bin", "louvain"
)
LOUVAIN_DEFAULT_SEED = 0x00C0FFEE


def compute_nearest_neighbors(submatrix, balltree, k: int, row_start: int):
    """Compute k nearest neighbors on a submatrix.

    Args:
        submatrix (np.ndarray): Data submatrix
        balltree: Nearest neighbor index (from sklearn)
        k: number of nearest neigbors to compute
        row_start: row offset into larger matrix

    Returns:
        a COO sparse adjacency matrix of nearest neighbor relations as (i,j,x).
    """
    nn_dist, nn_idx = balltree.query(submatrix, k=k + 1)

    # Remove the self-as-neighbors
    nn_idx = nn_idx[:, 1:]
    nn_dist = nn_dist[:, 1:]

    # Construct a COO sparse matrix of edges and distances
    i: np.ndarray[tuple[int, int], np.dtype[np.int_]] = np.repeat(
        row_start + np.arange(nn_idx.shape[0]), k
    )
    j: np.ndarray[tuple[int, int], np.dtype[np.int_]] = nn_idx.ravel().astype(int)
    return (i, j, nn_dist.ravel())


def write_nearest_neighbors(ijx, filename):
    """Write adjacency matrix to an HDF5 file.

    Args:
        ijx (tuple of i,j,x arrays of a COO sparse adjacency matrix)
    """
    i, j, x = ijx
    with h5py.File(filename, "w") as f:
        f.create_dataset("i", data=i)
        f.create_dataset("j", data=j)
        f.create_dataset("distance", data=x)


def write_adjacency_matrix_h5(matrix, h5_out_file, data_key="data"):
    """Write adjacency matrix to an HDF5 file from a sparse matrix.

    Specs::

        data    (array of matrix values, float64)
        indices (int64)
        indptr (int64)
        shape matrix.shape
    """
    with h5py.File(h5_out_file, "w") as h5_file:
        h5_file.create_dataset("indices", data=matrix.indices, dtype=np.int64)
        h5_file.create_dataset("indptr", data=matrix.indptr, dtype=np.int64)
        h5_file.create_dataset(data_key, data=matrix.data, dtype=np.float64)
        h5_file.create_dataset("shape", data=matrix.shape, dtype=np.int64)


def merge_nearest_neighbors(filenames, total_rows):
    """Merge nearest neighbor adjacency matrix HDF files.

    Returns:
        A sparse adjacency matrix
    """
    nn = sp_sparse.coo_matrix((total_rows, total_rows))
    for filename in filenames:
        with h5py.File(filename, "r") as h5:
            i: h5py.Dataset = h5["i"]
            nn += sp_sparse.coo_matrix((np.ones(len(i)), (i[:], h5["j"][:])), shape=nn.shape)
    return nn


def pipe_weighted_edgelist_to_convert(matrix, bin_filename, weight_filename):
    """Pipe a weighted edgelist (COO sparse matrix) to Louvain's convert utility."""
    raise ValueError("Unsupported method at the moment")


def run_louvain_weighted_clustering(
    bin_filename: str, weight_filename: str, louvain_out, seed=None
):
    """Run Louvain clustering on a weighted edge-list."""
    if seed is None:
        seed = LOUVAIN_DEFAULT_SEED
    with open(louvain_out, "wb") as f:
        tk_subproc.check_call(
            [
                LOUVAIN_BINPATH,
                bin_filename,
                "-w",
                weight_filename,
                "-q",
                "0",
                "-l",
                "-1",
                "-s",
                str(seed),
            ],
            stdout=f,
        )


def _write_to_proc(proc: subprocess.Popen[bytes], ijs: Iterable[tuple[int, int]]):
    assert proc.stdin is not None
    try:
        for ij in ijs:
            proc.stdin.write(b"%d\t%d\n" % ij)
    finally:
        proc.stdin.close()
        proc.wait()


def pipe_unweighted_edgelist_to_convert(matrix, bin_filename):
    """Pipe an unweighted edgelist (COO sparse matrix) to Louvain's convert utility."""
    with tk_subproc.Popen(
        [
            LOUVAIN_CONVERT_BINPATH,
            "-i",
            "-",
            "-o",
            bin_filename,
        ],
        stdin=subprocess.PIPE,
    ) as proc:
        # Check if the process terminated early
        time.sleep(3)
        retcode = proc.poll()
        if retcode is not None:
            proc.stdin.close()
            proc.wait()
            raise ChildProcessError(
                "'convert' command terminated early with exit code %d" % proc.returncode
            )

        # Stream text triplets to 'convert'
        print("Writing %d elements." % len(matrix.row))

        try:
            _write_to_proc(proc, zip(matrix.row, matrix.col))
        except OSError as e:
            if e.errno == errno.EPIPE:
                raise ChildProcessError(
                    "'convert' binary closed the pipe before we finished "
                    "writing to it. It terminated with exit code %d" % proc.returncode
                ) from e
            else:
                raise

        if proc.returncode != 0:
            raise ChildProcessError("'convert' command failed with exit code %d" % proc.returncode)

        if not os.path.exists(bin_filename):
            raise ChildProcessError(
                "'convert' failed to write the matrix file. Please see the "
                "standard error file (_stderr) to see if it emitted any errors."
            )


def run_louvain_unweighted_clustering(bin_filename, louvain_out, seed=None):
    """Run Louvain clustering on an unweighted edge-list."""
    if seed is None:
        seed = LOUVAIN_DEFAULT_SEED
    with open(louvain_out, "w") as f:
        tk_subproc.check_call(
            [
                LOUVAIN_BINPATH,
                bin_filename,
                "-q",
                "0",
                "-l",
                "-1",
                "-s",
                str(seed),
            ],
            stdout=f,
        )


def compute_snn_matrix(nn: sp_sparse.spmatrix, k_nearest: float | int) -> sp_sparse.coo_matrix:
    """Compute shared-nearest-neighbor matrix from a nearest-neighbor boolean matrix."""
    with LogPerf("tocsr"):
        nn = nn.tocsr(copy=False)

    # The SNN (shared nearest neighbor) similarity is
    #   The length of the nearest-neighbor intersection between two rows
    #   (divided by the max number of neighbors)
    # This can be computed via the dot products of rows in the boolean NN matrix
    with LogPerf("snn"):
        snn = (nn.dot(nn.T)) / float(k_nearest)

    # Use the SNN similarity in the modularity optimization algorithm
    # Louvain takes a text edge-list and converts to its own binary format
    with LogPerf("tocoo"):
        snn = snn.tocoo(copy=False)

    return snn


def load_louvain_results(
    num_barcodes: int, use_bcs: np.ndarray[int, np.dtype[np.uintp]], louvain_out: str | bytes
) -> np.ndarray[int, np.dtype[np.int64]]:
    """Load Louvain modularity results.

    Args:
        num_barcodes: total number of cell barcodes
        use_bcs: indices of bcs used for clustering
        louvain_out: path to louvain output file.

    Returns:
        Array of cluster labels.
    """
    labels: np.ndarray[int, np.dtype[np.int64]] = np.zeros(num_barcodes, dtype=np.int64)
    seen_idx = set()
    with open(louvain_out) as f:
        for line in f:
            used_bc_idx, cluster = map(int, line.strip().split(" "))

            # If we are seeing this index again, it's referring to
            # membership one or more levels up in the cluster hierarchy.
            if used_bc_idx in seen_idx:
                continue
            seen_idx.add(used_bc_idx)

            # 1-based cluster ids. Report 0 if barcode wasn't used.
            bc_idx = use_bcs[used_bc_idx]
            labels[bc_idx] = 1 + cluster
    return labels


def matrix_density(m: sp_sparse.spmatrix):
    return m.nnz / float(m.shape[0] * m.shape[1])


def build_neighbor_index(
    x: np.ndarray[tuple[int, int], np.dtype[np.number]], leaf_size: int
) -> sk_neighbors.BallTree:
    return sk_neighbors.BallTree(x, leaf_size=leaf_size)


def save_neighbor_index(index: sk_neighbors.BallTree, filename: str | bytes):
    with open(filename, "wb") as f:
        cPickle.dump(index, f, cPickle.HIGHEST_PROTOCOL)


def load_neighbor_index(filename: str | bytes) -> sk_neighbors.BallTree:
    with open(filename, "rb") as f:
        return cPickle.load(f)


def save_graphclust_h5(
    f: tables.File,
    labels: np.ndarray[int, np.dtype[np.int64]],
    clustering_type: str = cr_clustering.CLUSTER_TYPE_GRAPHCLUST,
):
    """Save graph clustering results in labels to h5py file f."""
    clustering_key = cr_clustering.format_clustering_key(clustering_type, 0)

    clustering = cr_clustering.create_clustering(
        clusters=labels,
        num_clusters=max(labels),
        cluster_score=0,
        clustering_type=clustering_type,
        global_sort_key=float("-inf"),  # always first
        description=cr_clustering.humanify_clustering_key(clustering_key),
    )

    group = f.create_group(f.root, analysis_constants.ANALYSIS_H5_CLUSTERING_GROUP)
    analysis_io.save_h5(f, group, clustering_key, clustering)


def save_ndarray_h5(data, path, dataset_name):
    """Save a numpy array to an hdf5 file."""
    with h5py.File(path, "w") as f:
        f.create_dataset(dataset_name, data=data)


def load_ndarray_h5(path, dataset_name):
    """Load an entire dataset into memory."""
    with h5py.File(path, "r") as f:
        return f[dataset_name][:]


def load_graphclust_from_h5(filename):
    """Read graph clustering data from HDF5 file in filename."""
    with tables.open_file(filename, "r") as f:
        # pylint: disable=protected-access
        group = f.root._v_groups[analysis_constants.ANALYSIS_H5_CLUSTERING_GROUP]

        # Take the first entry
        for key, clustering in analysis_io.load_h5_iter(group, cr_clustering.CLUSTERING):
            clustering_type, _ = cr_clustering.parse_clustering_key(key)
            if clustering_type == cr_clustering.CLUSTER_TYPE_GRAPHCLUST:
                return clustering
            elif clustering_type == cr_clustering.CLUSTER_TYPE_ATAC_GRAPHCLUST:
                return clustering
    raise KeyError(cr_clustering.CLUSTER_TYPE_GRAPHCLUST)
