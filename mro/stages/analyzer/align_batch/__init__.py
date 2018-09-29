#!/usr/bin/env python
#
# Copyright (c) 2018 10x Genomics, Inc. All rights reserved.
#
import numpy as np
import cPickle
import martian

import sklearn.neighbors as sk_neighbors
from sklearn.metrics.pairwise import rbf_kernel
from collections import defaultdict, OrderedDict, Counter
import tables
import h5py

import cellranger.matrix as cr_matrix
import cellranger.utils as cr_util
import cellranger.analysis.pca as cr_pca
import cellranger.h5_constants as h5_constants
import cellranger.analysis.constants as analysis_constants

__MRO__  = """
stage ALIGN_BATCH(
    in  h5     matrix_h5,
    in  pickle dimred_matrix,
    in  map[]  library_info,
    in  int    alignment_knn,
    in  float  alignment_alpha,
    in  float  alignment_sigma,
    in  bool   alignment_realign_panorama,
    in  bool   skip,
    out float  batch_score_before_alignment,
    out float  batch_score_after_alignment,
    out h5     aligned_pca_h5,
    out path   aligned_pca_csv,
    src py     "stages/analyzer/align_batch",
) split (
    in  int    batch_id,
    in  map    batch_to_bc_indices,
    in  pickle dimred_matrix,
    in  pickle idx_to_batch_id,
    in  bool   need_reorder_barcode,
    in  pickle barcode_reorder_index,
    out pickle batch_nearest_neighbor,
)
"""

MAX_MEM_GB = 64
NUM_ENTRIES_PER_MEM_GB = 2000000

def option(arg, default):
    return arg if arg is not None else default

def batch_effect_score(dimred_matrix, batch_ids, knn=100, subsample=0.1):
    """
    For each cell, serach KNN and calculate the proportion of cells from
    the same batch. Then compute the ratio of the proportion to the
    percentage (number of cell) of this batch. The batch_effect_score is
    defined as the average of the ratio. Closer to 1 means no batch effect.
    """
    num_bcs = dimred_matrix.shape[0]
    assert num_bcs == len(batch_ids)

    # batch percentage
    counter = Counter(batch_ids)
    batch_to_percentage = {batch: count*1.0/sum(counter.values()) for batch, count in counter.items()}

    # BallTree for KNN
    balltree  = sk_neighbors.BallTree(dimred_matrix, leaf_size=knn)

    select_bc_idx = np.array([i for i in range(num_bcs) if np.random.uniform() < subsample])
    knn_idx = balltree.query(dimred_matrix[select_bc_idx], k=knn+1, return_distance=False)

    same_batch_ratio = []
    for bc, neighbors in zip(select_bc_idx, knn_idx):
        batch_id = batch_ids[bc]
        same_batch = len([i for i in neighbors[1:] if batch_ids[i] == batch_id])
        same_batch_ratio.append((same_batch*1.0/knn)/batch_to_percentage[batch_id])

    return np.mean(same_batch_ratio)

def split(args):
    if args.skip:
        return {'chunks': [{'__mem_gb': h5_constants.MIN_MEM_GB}]}

    gg_id_to_batch_id, batch_id_to_name = {}, {}

    for lib in args.library_info:
        gg_id_to_batch_id[lib['gem_group']] = lib['batch_id']
        batch_id_to_name[lib['batch_id']] = lib['batch_name']

    # load just the barcodes from a matrix h5
    with h5py.File(args.matrix_h5, 'r') as f:
        group_name = f.keys()[0]
        bcs = cr_matrix.CountMatrix.load_bcs_from_h5_group(f[group_name])

    batch_ids = np.array([gg_id_to_batch_id[cr_util.split_barcode_seq(bc)[1]] for bc in bcs])

    with open(args.dimred_matrix) as fp:
        dimred_matrix = cPickle.load(fp)

    # re-order matrix such that barcodes from same batch are grouped together
    new_bc_indices = None
    batch_to_bc_indices = OrderedDict()
    idx_to_batch_id = np.full(dimred_matrix.shape[0], 0, dtype=np.int8)

    base = 0
    for b_id in range(len(batch_id_to_name)):
        batch_bc_indices = np.where(batch_ids == b_id)[0]
        if batch_bc_indices.shape[0] == 0:
            continue

        new_bc_indices = batch_bc_indices if new_bc_indices is None else np.append(new_bc_indices, batch_bc_indices)
        batch_to_bc_indices[b_id] = (base, base + batch_bc_indices.shape[0])
        idx_to_batch_id[base : base + batch_bc_indices.shape[0]] = b_id
        base += len(batch_bc_indices)

    # 1. check if needs re-order; 2. if needs re-order, store the original order
    need_reorder_barcode = (not np.all(np.diff(new_bc_indices)>=0))
    if need_reorder_barcode:
        dimred_matrix = dimred_matrix[new_bc_indices]
        barcode_reorder_index = np.argsort(new_bc_indices)

        barcode_reorder_index_file = martian.make_path('barcode_reorder_index.pickle')
        with open(barcode_reorder_index_file, 'wb') as fp:
            cPickle.dump(barcode_reorder_index, fp, cPickle.HIGHEST_PROTOCOL)
    else:
        barcode_reorder_index_file = None

    dimred_matrix_file = martian.make_path('dimred_matrix.pickle')
    with open(dimred_matrix_file, 'wb') as fp:
        cPickle.dump(dimred_matrix, fp, cPickle.HIGHEST_PROTOCOL)

    idx_to_batch_id_file = martian.make_path('idx_to_batch_id.pickle')
    with open(idx_to_batch_id_file, 'wb') as fp:
        cPickle.dump(idx_to_batch_id, fp, cPickle.HIGHEST_PROTOCOL)

    matrix_mem_gb = float(dimred_matrix.shape[0] * dimred_matrix.shape[1]) / NUM_ENTRIES_PER_MEM_GB
    mem_gb = min(MAX_MEM_GB, max(int(matrix_mem_gb), h5_constants.MIN_MEM_GB))

    chunks = []
    for batch_id, (start_idx, end_idx) in batch_to_bc_indices.items():
        chunks.append({
            '__mem_gb': mem_gb,
            'batch_id': str(batch_id),
            'batch_to_bc_indices': batch_to_bc_indices,
            'dimred_matrix': dimred_matrix_file,
            'idx_to_batch_id': idx_to_batch_id_file,
            'need_reorder_barcode': need_reorder_barcode,
            'barcode_reorder_index': barcode_reorder_index_file,
        })

    return {
        'chunks': chunks, 
        'join': {'__mem_gb': MAX_MEM_GB}
    }

def find_knn(curr_matrix, ref_matrix, knn):
    """
    for each row in cur_matrix, find k nearest neighbors in ref_matrix,
    return an array of shape=[curr_matrix.shape[0] * knn, ], which stores
    the index of nearest neighbors in ref_matrix
    """
    balltree = sk_neighbors.BallTree(ref_matrix, leaf_size=knn)
    nn_idx = balltree.query(curr_matrix, k=knn, return_distance=False)
    return nn_idx.ravel().astype(int)

def main(args, outs):
    if args.skip:
        return

    with open(args.dimred_matrix) as fp:
        dimred_matrix = cPickle.load(fp)

    alignment_knn = option(args.alignment_knn, analysis_constants.ALIGNMENT_KNN)

    batch_to_bc_indices = args.batch_to_bc_indices
    batch_start_idx  = batch_to_bc_indices[args.batch_id][0]
    batch_end_idx = batch_to_bc_indices[args.batch_id][1]
    cur_matrix = dimred_matrix[batch_start_idx:batch_end_idx,:]

    # nearest neighbor pair: stores the nearest neighbors from match_i to match_j
    # key = (batch_i, batch_j), values = set((idx_i, idx_j), ...), the index here is the global index
    batch_nearest_neighbor = defaultdict(set)

    from_idx, to_idx = None, None
    # Batch balanced KNN
    for batch in args.batch_to_bc_indices:
        if batch == args.batch_id:
            continue

        ref_matrix = dimred_matrix[batch_to_bc_indices[batch][0]:batch_to_bc_indices[batch][1],]
        nn_idx_right = find_knn(cur_matrix, ref_matrix, alignment_knn)

        # convert index (in cur_matrix and ref_matrix) to global index (in dimred_matrix)
        nn_idx_left = np.repeat(np.arange(cur_matrix.shape[0]) + batch_start_idx, alignment_knn)
        nn_idx_right += batch_to_bc_indices[batch][0]

        from_idx = nn_idx_left if from_idx is None else np.concatenate([from_idx, nn_idx_left])
        to_idx = nn_idx_right if to_idx is None else np.concatenate([to_idx, nn_idx_right])

        for i, j in zip(from_idx, to_idx):
            batch_nearest_neighbor[(args.batch_id, batch)].add((i,j))

    outs.batch_nearest_neighbor = martian.make_path('batch_nearest_neighbor.pickle')
    with open(outs.batch_nearest_neighbor, 'wb') as fp:
        cPickle.dump(batch_nearest_neighbor, fp, cPickle.HIGHEST_PROTOCOL)

    return

def correction_vector(mnn_cur, mnn_ref, cur_submatrix, sigma):
    """
    Compute the batch-correction vector

    1. For each MNN pair in current dataset and the reference, a pair-specific
    batch-correction vector is computed as the vector difference between the
    paired cells.
    2. For each barcode in cur dataset, a batch-correction vector is calculated
    as a weighted average of these pair-specific vectors, as computed with a
    Gaussian kernel.
    """
    bias = mnn_ref - mnn_cur
    weights = rbf_kernel(cur_submatrix, mnn_cur, gamma=0.5*sigma)
    return np.dot(weights, bias) / np.tile(np.sum(weights, axis=1), (cur_submatrix.shape[1], 1)).T

def join(args, outs, chunk_defs, chunk_outs):
    if args.skip:
        return

    batch_to_bc_indices = chunk_defs[0].batch_to_bc_indices

    with open(chunk_defs[0].dimred_matrix) as fp:
        dimred_matrix = cPickle.load(fp)

    with open(chunk_defs[0].idx_to_batch_id) as fp:
        idx_to_batch_id = cPickle.load(fp)

    # batch score before alignment
    outs.batch_score_before_alignment = batch_effect_score(dimred_matrix, idx_to_batch_id)

    nn_pairs = {}
    for chunk_out in chunk_outs:
        with open(chunk_out.batch_nearest_neighbor) as fp:
            batch_nearest_neighbor = cPickle.load(fp)
        for k,v in batch_nearest_neighbor.items():
            nn_pairs[k] = v

    mutual_nn = {} # mnn between batches
    overlap_percentage = {} # percentage of matching cells (max of batch i and batch j)

    for i in batch_to_bc_indices:
        batch_i_size = batch_to_bc_indices[i][1] - batch_to_bc_indices[i][0]
        for j in batch_to_bc_indices:
            if i >= j:
                continue

            if (i,j) not in nn_pairs or (j,i) not in nn_pairs:
                continue

            nn_ij = nn_pairs[(i,j)]
            nn_ji = set([(y,x) for x,y in nn_pairs[(j,i)]])
            mutual_nn[(i,j)] = nn_ij & nn_ji

            batch_j_size = batch_to_bc_indices[j][1] - batch_to_bc_indices[j][0]
            overlap_percentage[(i,j)] = max(
                float(len(set([idx for idx, _ in mutual_nn[(i,j)]]))) / batch_i_size,
                float(len(set([idx for _, idx in mutual_nn[(i,j)]]))) / batch_j_size
            )

    alignment_alpha = option(args.alignment_alpha, analysis_constants.ALIGNMENT_ALPHA)
    alignment_realign_panorama = option(args.alignment_realign_panorama, analysis_constants.ALIGNMENT_REALIGN_PANORAMA)
    alignment_sigma = option(args.alignment_sigma, analysis_constants.ALIGNMENT_SIGMA)

    alignments = [k for k,v in sorted(overlap_percentage.items(), key=lambda x: x[1], reverse=True) if v > alignment_alpha]

    ## panorama stitch ##
    aligned_dimred_matrix = dimred_matrix
    batch_id_to_alignment_count = defaultdict(int)
    panoramas = [] # a list of stitched panoramas
    batch_id_to_alignment_count = defaultdict(int)

    for (i,j) in alignments:
        panorama_idx_i, panorama_idx_j = None, None
        for idx, panorama in enumerate(panoramas):
            if i in panorama:
                assert panorama_idx_i is None, "batch is found in multiple panorama"
                panorama_idx_i = idx

            if j in panorama:
                assert panorama_idx_j is None, "batch is found in multiple panorama"
                panorama_idx_j = idx

        if panorama_idx_i is None:
            panoramas.append(set([i]))
            panorama_idx_i = len(panoramas) - 1

        if panorama_idx_j is None:
            panoramas.append(set([j]))
            panorama_idx_j = len(panoramas) - 1

        # re-align within panorama is enabled
        if alignment_realign_panorama is True:
            batch_id_to_alignment_count[i] += 1
            batch_id_to_alignment_count[j] += 1
            if batch_id_to_alignment_count[i] > 3 and batch_id_to_alignment_count[j] > 3:
                continue
        # re-align within panorama is disable
        else:
            if panorama_idx_i == panorama_idx_j:
                continue

        # use the panorama with the larger shape as reference
        panorama_i_size = sum([batch_to_bc_indices[b][1] - batch_to_bc_indices[b][0] for b in panoramas[panorama_idx_i]])
        panorama_j_size = sum([batch_to_bc_indices[b][1] - batch_to_bc_indices[b][0] for b in panoramas[panorama_idx_j]])
        if panorama_i_size < panorama_j_size:
            panorama_idx_i, panorama_idx_j = panorama_idx_j, panorama_idx_i

        # merge panorama j (current) to panorama i (reference)
        # ordered batch id <-> ordered index range
        batches_in_panorama_j = sorted(list(panoramas[panorama_idx_j]))
        cur_submatrix_idx = np.concatenate([np.arange(batch_to_bc_indices[b][0], batch_to_bc_indices[b][1])
                                            for b in batches_in_panorama_j])
        cur_submatrix = aligned_dimred_matrix[cur_submatrix_idx]

        matches = []
        for ref in panoramas[panorama_idx_i]:
            for cur in panoramas[panorama_idx_j]:
                if ref < cur and (ref, cur) in mutual_nn:
                    matches.extend([(c,r) for r,c in mutual_nn[(ref, cur)]])
                if ref > cur and (cur, ref) in mutual_nn:
                    matches.extend(mutual_nn[(cur, ref)])

        mnn_cur = aligned_dimred_matrix[[a for a, _ in matches], :]
        mnn_ref = aligned_dimred_matrix[[b for _, b in matches], :]

        corr_vector = correction_vector(mnn_cur, mnn_ref, cur_submatrix, alignment_sigma)
        assert(corr_vector.shape == cur_submatrix.shape)
        cur_submatrix += corr_vector

        base = 0
        for b in batches_in_panorama_j:
            batch_size = batch_to_bc_indices[b][1] - batch_to_bc_indices[b][0]
            aligned_dimred_matrix[batch_to_bc_indices[b][0]:batch_to_bc_indices[b][1], :] = cur_submatrix[base:(base + batch_size), :]
            base += batch_size

        # merge panoramas and delete panorama j
        if panorama_idx_i != panorama_idx_j:
            panoramas[panorama_idx_i].update(panoramas[panorama_idx_j])
            panoramas.pop(panorama_idx_j)

    # batch score after alignment
    outs.batch_score_after_alignment = batch_effect_score(aligned_dimred_matrix, idx_to_batch_id)

    if chunk_defs[0].need_reorder_barcode:
        with open(chunk_defs[0].barcode_reorder_index) as fp:
            barcode_reorder_index = cPickle.load(fp)
            aligned_dimred_matrix = aligned_dimred_matrix[barcode_reorder_index]

    # TODO: currently save the aligned matrix (dimred) into pca h5 format, for downstream analysis
    aligned_pca = cr_pca.PCA(aligned_dimred_matrix, np.zeros(0), np.zeros(0), np.zeros(0), np.zeros(0))
    pca_map = {aligned_dimred_matrix.shape[1] : aligned_pca}
    filters = tables.Filters(complevel = h5_constants.H5_COMPRESSION_LEVEL)
    with tables.open_file(outs.aligned_pca_h5, 'w', filters = filters) as h5_file:
        cr_pca.save_pca_h5(pca_map, h5_file)

    matrix = cr_matrix.CountMatrix.load_h5_file(args.matrix_h5)
    cr_pca.save_pca_csv(pca_map, matrix, outs.aligned_pca_csv)
