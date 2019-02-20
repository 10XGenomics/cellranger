#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#
import collections
import cPickle
import errno
import h5py
import itertools
import numpy as np
import os
import scipy.sparse as sp_sparse
import sklearn.neighbors as sk_neighbors
import subprocess
import tables
import time
import tenkit.log_subprocess as tk_subproc
import cellranger.analysis.io as analysis_io
import cellranger.analysis.clustering as cr_clustering
import cellranger.analysis.constants as analysis_constants
from cellranger.logperf import LogPerf

GRAPHCLUST = collections.namedtuple('GRAPHCLUST', ['clusters'])

LOUVAIN_CONVERT_BINPATH = 'convert'
LOUVAIN_BINPATH = 'louvain'

def compute_nearest_neighbors(submatrix, balltree, k, row_start):
    """ Compute k nearest neighbors on a submatrix
    Args: submatrix (np.ndarray): Data submatrix
          balltree: Nearest neighbor index (from sklearn)
          k: number of nearest neigbors to compute
          row_start: row offset into larger matrix
    Returns a COO sparse adjacency matrix of nearest neighbor relations as (i,j,x)"""

    nn_dist, nn_idx = balltree.query(submatrix, k=k+1)

    # Remove the self-as-neighbors
    nn_idx = nn_idx[:,1:]
    nn_dist = nn_dist[:,1:]

    # Construct a COO sparse matrix of edges and distances
    i = np.repeat(row_start + np.arange(nn_idx.shape[0]), k)
    j = nn_idx.ravel().astype(int)
    return (i, j, nn_dist.ravel())

def write_nearest_neighbors(ijx, filename):
    """ Write adjacency matrix to an HDF5 file.
    Args: ijx (tuple of i,j,x arrays of a COO sparse adjacency matrix)
    """
    i, j, x = ijx
    with h5py.File(filename, 'w') as f:
        f.create_dataset('i', data=i)
        f.create_dataset('j', data=j)
        f.create_dataset('distance', data=x)

def merge_nearest_neighbors(filenames, total_rows):
    """ Merge nearest neighbor adjacency matrix HDF files.
    Returns: A sparse adjacency matrix """
    nn = sp_sparse.coo_matrix((total_rows, total_rows))
    for filename in filenames:
        h5 = h5py.File(filename, 'r')
        nn += sp_sparse.coo_matrix((np.ones(len(h5['i'])),
                                    (h5['i'][:], h5['j'][:])),
                                   shape=nn.shape)
        h5.close()
    return nn


def pipe_weighted_edgelist_to_convert(matrix, bin_filename, weight_filename):
    """ Pipe a weighted edgelist (COO sparse matrix) to Louvain's convert utility """
    raise ValueError('Unsupported method at the moment')

    devnull = open(os.devnull, 'w')

    proc = tk_subproc.Popen([LOUVAIN_CONVERT_BINPATH,
                           '-i', '/dev/stdin',
                           '-o', bin_filename,
                           '-w', weight_filename,
                         ], stdin=subprocess.PIPE, stdout=devnull, stderr=devnull)

    # Stream text triplets to 'convert'
    for ijx in itertools.izip(matrix.row, matrix.col, matrix.data):
        proc.stdin.write('%d\t%d\t%f\n' % ijx)

    proc.stdin.close()
    proc.wait()
    devnull.close()

def run_louvain_weighted_clustering(bin_filename, weight_filename, louvain_out):
    """ Run Louvain clustering on a weighted edge-list """
    with open(louvain_out, 'w') as f:
        tk_subproc.check_call([LOUVAIN_BINPATH,
                               bin_filename,
                               '-w', weight_filename,
                               '-q', '0',
                               '-l', '-1',
                           ], stdout=f)

def pipe_unweighted_edgelist_to_convert(matrix, bin_filename):
    """ Pipe an unweighted edgelist (COO sparse matrix) to Louvain's convert utility """

    proc = tk_subproc.Popen([LOUVAIN_CONVERT_BINPATH,
                           '-i', '-',
                           '-o', bin_filename,
                         ], stdin=subprocess.PIPE)


    # Check if the process terminated early
    time.sleep(3)
    retcode = proc.poll()
    if retcode is not None:
        proc.stdin.close()
        proc.wait()
        raise Exception("'convert' command terminated early with exit code %d" % proc.returncode)

    # Stream text triplets to 'convert'
    print 'Writing %d elements.' % len(matrix.row)

    try:
        for ij in itertools.izip(matrix.row, matrix.col):
            proc.stdin.write('%d\t%d\n' % ij)
        proc.stdin.close()
    except IOError as e:
        if e.errno == errno.EPIPE:
            proc.stdin.close()
            proc.wait()
            raise Exception("'convert' binary closed the pipe before we finished writing to it. It terminated with exit code %d" % proc.returncode)

    proc.wait()

    if proc.returncode != 0:
        raise Exception("'convert' command failed with exit code %d" % proc.returncode)

    if not os.path.exists(bin_filename):
        raise Exception("'convert' failed to write the matrix file. Please see the standard error file (_stderr) to see if it emitted any errors.")

def run_louvain_unweighted_clustering(bin_filename, louvain_out):
    """ Run Louvain clustering on an unweighted edge-list """
    with open(louvain_out, 'w') as f:
        tk_subproc.check_call([LOUVAIN_BINPATH,
                               bin_filename,
                               '-q', '0',
                               '-l', '-1',
                           ], stdout=f)

def compute_snn_matrix(nn, k_nearest):
    """ Compute shared-nearest-neighbor matrix from a nearest-neighbor boolean matrix """
    with LogPerf('tocsr'):
        nn = nn.tocsr(copy=False)

    # The SNN (shared nearest neighbor) similarity is
    #   The length of the nearest-neighbor intersection between two rows
    #   (divided by the max number of neighbors)
    # This can be computed via the dot products of rows in the boolean NN matrix
    with LogPerf('snn'):
        snn = (nn.dot(nn.T)) / float(k_nearest)

    # Use the SNN similarity in the modularity optimization algorithm
    # Louvain takes a text edge-list and converts to its own binary format
    with LogPerf('tocoo'):
        snn = snn.tocoo(copy=False)

    return snn

def load_louvain_results(num_barcodes, use_bcs, louvain_out):
    """ Load Louvain modularity results.
    Args: num_barcodes - total number of cell barcodes
          use_bcs - indices of bcs used for clustering
          louvain_out - path to louvain output file.
    Returns: Array of cluster labels  """

    labels = np.zeros(num_barcodes, dtype=int)
    with open(louvain_out) as f:
        for line in f:
            used_bc_idx, cluster = map(int, line.strip().split(' '))
            # Take max of community cluster ids reported by Louvain.
            # 1-based cluster ids. Report 0 if barcode wasn't used.

            bc_idx = use_bcs[used_bc_idx]
            labels[bc_idx] = max(labels[bc_idx], 1 + cluster)
    return labels

def matrix_density(m):
    return m.nnz / float(m.shape[0]*m.shape[1])

def build_neighbor_index(x, leaf_size):
    return sk_neighbors.BallTree(x, leaf_size=leaf_size)

def save_neighbor_index(index, filename):
    with open(filename, 'wb') as f:
        cPickle.dump(index, f, cPickle.HIGHEST_PROTOCOL)

def load_neighbor_index(filename):
    with open(filename) as f:
        return cPickle.load(f)

def save_graphclust_h5(f, labels):
    clustering_key = cr_clustering.format_clustering_key(cr_clustering.CLUSTER_TYPE_GRAPHCLUST, 0)

    clustering = cr_clustering.create_clustering(clusters=labels,
                                                 num_clusters=max(labels),
                                                 cluster_score=0,
                                                 clustering_type=cr_clustering.CLUSTER_TYPE_GRAPHCLUST,
                                                 global_sort_key=float('-inf'), # always first
                                                 description=cr_clustering.humanify_clustering_key(clustering_key))

    group = f.create_group(f.root, analysis_constants.ANALYSIS_H5_CLUSTERING_GROUP)
    analysis_io.save_h5(f, group, clustering_key, clustering)

def save_ndarray_h5(data, path, dataset_name):
    """ Save a numpy array to an hdf5 file """
    with h5py.File(path, 'w') as f:
        f.create_dataset(dataset_name, data=data)

def load_ndarray_h5(path, dataset_name):
    """ Load an entire dataset into memory """
    with h5py.File(path, 'r') as f:
        return f[dataset_name][:]

def load_graphclust_from_h5(filename):
    with tables.open_file(filename, 'r') as f:
        group = f.root._v_groups[analysis_constants.ANALYSIS_H5_CLUSTERING_GROUP]

        # Take the first entry
        for key, clustering in analysis_io.load_h5_iter(group, cr_clustering.CLUSTERING):
            clustering_type, _ = cr_clustering.parse_clustering_key(key)
            if clustering_type == cr_clustering.CLUSTER_TYPE_GRAPHCLUST:
                return clustering
