#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#
import martian
import numpy as np
import sys
import cellranger.analysis.clustering as cr_clustering
import cellranger.analysis.graphclust as cr_graphclust
import cellranger.analysis.io as cr_io
from cellranger.analysis.singlegenome import SingleGenomeAnalysis
import cellranger.constants as cr_constants
from cellranger.logperf import LogPerf
import cellranger.utils as cr_utils

__MRO__ = """
stage RUN_GRAPH_CLUSTERING(
    in  h5     matrix_h5,
    in  h5     pca_h5,
    in  int    num_neighbors       "Use this many neighbors",
    in  float  neighbor_a          "Use larger of (a+b*log10(n_cells) neighbors or num_neighbors",
    in  float  neighbor_b          "Use larger of (a+b*log10(n_cells) neighbors or num_neighbors",
    in  int    num_bcs             "Use this many cell-barcodes in clustering",
    in  int    input_pcs           "Use top N PCs",
    in  int    balltree_leaf_size,
    in  string similarity_type     "Type of similarity to use (nn or snn)",
    in  bool   skip,
    in  bool   is_multi_genome,
    out h5     chunked_neighbors,
    out h5     clusters_h5,
    out path   clusters_csv,
    src py     "stages/analyzer/run_graph_clustering",
) split using (
    in  pickle neighbor_index,
    in  h5     submatrix,
    in  int    row_start,
    in  int    total_rows,
    in  int    k_nearest,
    in  h5     use_bcs,
)
"""

# 1e6 cells => ~64 chunks
NN_QUERIES_PER_CHUNK = 15000

DEFAULT_BALLTREE_LEAFSIZE = 40

# Memory usage in join, empirically determined
NN_ENTRIES_PER_MEM_GB = 10000000

# Unweighted nearest neighbor (boolean: is-nearest-neighbor)
NN_SIMILARITY = 'nn'

# Shared nearest neighbor (fraction of neighbors shared)
SNN_SIMILARITY = 'snn'

SIMILARITY_TYPES = [NN_SIMILARITY, SNN_SIMILARITY]

# TODO: Martian needs to provide a way to give split more memory.
# Workaround is mrp --overrides
def split(args):
    np.random.seed(0)

    if args.skip or args.is_multi_genome:
        return {'chunks': [{'__mem_gb': cr_constants.MIN_MEM_GB}]}

    if args.similarity_type not in SIMILARITY_TYPES:
        martian.exit("Unsupported similarity type: %s. Must be one of: %s" % (args.similarity_type, ','.join(SIMILARITY_TYPES)))

    with LogPerf('load'):
        pca_mat = SingleGenomeAnalysis.load_pca_from_h5(args.pca_h5).transformed_pca_matrix

    # Subselect barcodes if desired
    if args.num_bcs is None:
        use_bcs = np.arange(pca_mat.shape[0])
    else:
        use_bcs = np.random.choice(pca_mat.shape[0], args.num_bcs, replace=False)
        pca_mat = pca_mat[use_bcs,:]

    # Record indices of selected barcodes
    use_bcs_path = martian.make_path('use_bcs.h5')
    cr_graphclust.save_ndarray_h5(use_bcs, use_bcs_path, 'use_bcs')

    # Subselect PCs if desired
    if args.input_pcs is not None:
        n_pcs = min(pca_mat.shape[1], args.input_pcs)
        pca_mat = pca_mat[:,np.arange(n_pcs)]

    # Build the nearest neighbor query index
    with LogPerf('nn_build'):
        balltree = cr_graphclust.build_neighbor_index(pca_mat,
                                                      args.balltree_leaf_size or DEFAULT_BALLTREE_LEAFSIZE)
        neighbor_index = martian.make_path('neighbor_index.pickle')
        cr_graphclust.save_neighbor_index(balltree, neighbor_index)

    # Compute the actual number of nearest neighbors we'll use
    given_num_neighbors = args.num_neighbors if args.num_neighbors is not None else cr_constants.GRAPHCLUST_NEIGHBORS_DEFAULT
    given_neighbor_a = args.neighbor_a if args.neighbor_a is not None else cr_constants.GRAPHCLUST_NEIGHBOR_A_DEFAULT
    given_neighbor_b = args.neighbor_b if args.neighbor_b is not None else cr_constants.GRAPHCLUST_NEIGHBOR_B_DEFAULT

    # Take max of {num_neighbors, a + b*log10(n)}
    use_neighbors = int(max(given_num_neighbors, np.round(given_neighbor_a + given_neighbor_b * np.log10(len(use_bcs)))))

    # Clamp to [1, n - 1]
    num_neighbors = max(1, min(use_neighbors, len(use_bcs)-1))
    print "Using %d neighbors" % num_neighbors

    # Divide the PCA matrix up into rows for NN queries
    with LogPerf('chunk_pca'):
        chunks = []
        for row_start in xrange(0, pca_mat.shape[0], NN_QUERIES_PER_CHUNK):
            row_end = min(row_start + NN_QUERIES_PER_CHUNK, pca_mat.shape[0])

            # Write the pca submatrix to an h5 file
            submatrix_path = martian.make_path('%d_submatrix.h5' % row_start)
            cr_graphclust.save_ndarray_h5(pca_mat[row_start:row_end, :], submatrix_path, 'submatrix')

            chunks.append({
                'neighbor_index': neighbor_index,
                'submatrix': submatrix_path,
                'row_start': row_start,
                'total_rows': pca_mat.shape[0],
                'k_nearest': num_neighbors,
                'use_bcs': use_bcs_path,
            })

    if args.similarity_type == SNN_SIMILARITY:
        join_mem_gb = 64
        join_threads = 4 # Overallocate
    else:
        # Scale memory with size of nearest-neighbor adjacency matrix
        join_mem_gb = max(cr_constants.MIN_MEM_GB, int(np.ceil((num_neighbors * len(use_bcs)) / NN_ENTRIES_PER_MEM_GB)))
        # HACK: use more threads for bigger mem requests to avoid mem oversubscription on clusters that don't enforce it
        join_threads = cr_utils.get_thread_request_from_mem_gb(join_mem_gb)

    return {
        'chunks': chunks,
        'join': {
            '__mem_gb': join_mem_gb,
            '__threads': join_threads,
        }}


def main(args, outs):
    np.random.seed(0)

    if args.skip or args.is_multi_genome:
        return

    with LogPerf('submatrix_load'):
        submatrix = cr_graphclust.load_ndarray_h5(args.submatrix, 'submatrix')

    with LogPerf('nn_idx_load'):
        balltree = cr_graphclust.load_neighbor_index(args.neighbor_index)

    with LogPerf('nn_query'):
        nn_matrix = cr_graphclust.compute_nearest_neighbors(submatrix, balltree, args.k_nearest, args.row_start)
        cr_graphclust.write_nearest_neighbors(nn_matrix, outs.chunked_neighbors)

def join(args, outs, chunk_defs, chunk_outs):
    if args.skip or args.is_multi_genome:
        return
    # Merge the neighbor matrices
    with LogPerf('merge_nn'):
        nn = cr_graphclust.merge_nearest_neighbors([chunk.chunked_neighbors for chunk in chunk_outs],
                                                   chunk_defs[0].total_rows)
    print 'nn\tnn_nodes\t%0.4f' % nn.shape[0]
    print 'nn\tnn_links\t%0.4f' % nn.nnz
    print 'nn\tnn_density\t%0.4f' % cr_graphclust.matrix_density(nn)
    sys.stdout.flush()

    matrix_bin = martian.make_path('matrix.bin')
    matrix_weights = martian.make_path('matrix.weights')
    louvain_out = martian.make_path('louvain.out')

    if args.similarity_type == 'snn':
        snn = cr_graphclust.compute_snn_matrix(nn, chunk_defs[0].k_nearest)

        print 'snn\tsnn_nodes\t%d' % snn.shape[0]
        print 'snn\tsnn_links\t%d' % (snn.nnz/2)
        print 'snn\tsnn_density\t%0.4f' % ((snn.nnz) / float(snn.shape[0]*(snn.shape[0]-1)))
        sys.stdout.flush()

        with LogPerf('convert'):
            cr_graphclust.pipe_weighted_edgelist_to_convert(snn, matrix_bin, matrix_weights)

        with LogPerf('louvain'):
            cr_graphclust.run_louvain_weighted_clustering(matrix_bin, matrix_weights, louvain_out)

    else:
        with LogPerf('tocoo'):
            nn = nn.tocoo(copy=False)

        with LogPerf('convert'):
            cr_graphclust.pipe_unweighted_edgelist_to_convert(nn, matrix_bin)

        with LogPerf('louvain'):
            cr_graphclust.run_louvain_unweighted_clustering(matrix_bin, louvain_out)

    with LogPerf('load_bcs'):
        barcodes = SingleGenomeAnalysis.load_bcs_from_matrix_h5(args.matrix_h5)

    use_bcs = cr_graphclust.load_ndarray_h5(chunk_defs[0].use_bcs, 'use_bcs')

    labels = cr_graphclust.load_louvain_results(len(barcodes), use_bcs, louvain_out)

    labels = cr_clustering.relabel_by_size(labels)

    # Save cluster results
    with cr_io.open_h5_for_writing(outs.clusters_h5) as f:
        cr_graphclust.save_graphclust_h5(f, labels)

    clustering_key = cr_clustering.format_clustering_key(cr_clustering.CLUSTER_TYPE_GRAPHCLUST, 0)

    cr_clustering.save_clustering_csv(outs.clusters_csv, clustering_key, labels, barcodes)

    outs.chunked_neighbors = None
