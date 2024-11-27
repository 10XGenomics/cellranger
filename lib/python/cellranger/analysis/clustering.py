#!/usr/bin/env python
#
# Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
#
from __future__ import annotations

import os
import sys
from collections.abc import Iterable, Sequence
from typing import TYPE_CHECKING, ClassVar, NamedTuple, Protocol

import numpy as np

import cellranger.analysis.io as analysis_io

if TYPE_CHECKING:
    from cellranger.wrapped_tables import tables

    class NamedTupleProtocol(Protocol):
        _fields: ClassVar[Sequence[str]]


GEX_PREFIX = "gene_expression"
AB_PREFIX = "antibody_capture"
ATAC_PREFIX = "peaks"
CLUSTER_TYPE_KMEANS = "kmeans"
CLUSTER_TYPE_ANTIBODY_KMEANS = AB_PREFIX + "_" + CLUSTER_TYPE_KMEANS
CLUSTER_TYPE_ATAC = ATAC_PREFIX + "_" + CLUSTER_TYPE_KMEANS
CLUSTER_TYPE_GRAPHCLUST = GEX_PREFIX + "_" + "graphclust"
CLUSTER_TYPE_ATAC_GRAPHCLUST = ATAC_PREFIX + "_" + "graphclust"
CLUSTER_TYPE_ANTIBODY_GRAPHCLUST = f"{AB_PREFIX}_graphclust"
CLUSTER_TYPE_KMEDOIDS = "kmedoids"
CLUSTER_TYPE_CELLTYPES = "celltype"


class CLUSTERING(NamedTuple):
    clusters: np.ndarray[int, np.dtype[np.int64]]
    num_clusters: np.int64
    cluster_score: np.float64
    clustering_type: np.bytes_
    global_sort_key: np.float64
    description: np.bytes_


def format_legacy_clustering_key(cluster_type: str, cluster_param: int):
    assert cluster_type == CLUSTER_TYPE_KMEANS
    return "_%d" % cluster_param


def format_clustering_key(cluster_type: str, cluster_param: int):
    """Generate a machine-readable string that describes a particular clustering."""
    if cluster_type == CLUSTER_TYPE_KMEANS:
        return f"{GEX_PREFIX}_{CLUSTER_TYPE_KMEANS}_{cluster_param}_clusters"
    elif cluster_type == CLUSTER_TYPE_ANTIBODY_KMEANS:
        return "%s_%s_%d_clusters" % (AB_PREFIX, CLUSTER_TYPE_KMEANS, cluster_param)
    elif cluster_type == CLUSTER_TYPE_ATAC:
        return "%s_%s_%d_clusters" % (ATAC_PREFIX, CLUSTER_TYPE_KMEANS, cluster_param)
    elif cluster_type == CLUSTER_TYPE_KMEDOIDS:
        return "%s_%d_clusters" % (CLUSTER_TYPE_KMEDOIDS, cluster_param)
    elif cluster_type in (
        CLUSTER_TYPE_GRAPHCLUST,
        CLUSTER_TYPE_ATAC_GRAPHCLUST,
        CLUSTER_TYPE_CELLTYPES,
    ):
        return cluster_type
    else:
        raise ValueError(f"Unsupported cluster type: {cluster_type}")


def _parse_number_of_clusters(clustering_key: str):
    """Given clustering keys, output the integer number of clusters.

    Keys are of a form like _antibody_capture_kmeans_10_clusters or
    _gene_expression_kmeans_3_clusters.

    >>> _parse_number_of_clusters("_antibody_capture_kmeans_10_clusters")
    10
    >>> _parse_number_of_clusters("_gene_expression_kmeans_3_clusters")
    3
    >>> _parse_number_of_clusters("_peaks_kmeans_5_clusters")
    5
    """
    n_clusters = filter(str.isdigit, clustering_key)
    n_clusters = "".join([s for s in n_clusters])
    return int(n_clusters)


def parse_clustering_key(clustering_key: str):
    """Parse the output of format_clustering_key()."""
    if clustering_key in (
        CLUSTER_TYPE_GRAPHCLUST,
        CLUSTER_TYPE_ATAC_GRAPHCLUST,
        CLUSTER_TYPE_ANTIBODY_GRAPHCLUST,
    ):
        return (clustering_key, 0)
    elif clustering_key.startswith(GEX_PREFIX):
        n_clusters = _parse_number_of_clusters(clustering_key)
        return (CLUSTER_TYPE_KMEANS, n_clusters)
    elif clustering_key.startswith(AB_PREFIX):
        n_clusters = _parse_number_of_clusters(clustering_key)
        return (CLUSTER_TYPE_ANTIBODY_KMEANS, n_clusters)
    elif clustering_key.startswith(ATAC_PREFIX):
        n_clusters = _parse_number_of_clusters(clustering_key)
        return (CLUSTER_TYPE_ATAC, n_clusters)
    elif clustering_key.startswith(CLUSTER_TYPE_KMEDOIDS):
        _, n_clusters, _ = clustering_key.split("_")
        return (CLUSTER_TYPE_KMEDOIDS, int(n_clusters))
    elif clustering_key == CLUSTER_TYPE_CELLTYPES:
        return (clustering_key, 0)
    else:
        raise ValueError(f"Unsupported clustering type for clustering key: {clustering_key}")


def humanify_clustering_key(clustering_key: str):
    """Make a cluster_key string human-readable."""
    cluster_type, cluster_param = parse_clustering_key(clustering_key)
    if cluster_type == CLUSTER_TYPE_GRAPHCLUST:
        return "Gene Expression Graph-based"
    elif cluster_type == CLUSTER_TYPE_ATAC_GRAPHCLUST:
        return "Peaks Graph-based"
    elif cluster_type == CLUSTER_TYPE_KMEANS:
        return "Gene Expression K-means (K=%d)" % cluster_param
    elif cluster_type == CLUSTER_TYPE_ANTIBODY_KMEANS:
        return "Antibody Capture K-means (K=%d)" % cluster_param
    elif cluster_type == CLUSTER_TYPE_ATAC:
        return "Peaks K-means (K=%d)" % cluster_param
    elif cluster_type == CLUSTER_TYPE_KMEDOIDS:
        return "K-medoids (K=%d)" % cluster_param
    elif cluster_type == CLUSTER_TYPE_CELLTYPES:
        return "Celltypes"
    else:
        raise ValueError(
            f"Unsupported clustering type {cluster_type} for clustering key: {clustering_key}"
        )


def relabel_by_size(
    labels: np.ndarray[int, np.dtype[np._IntType]]
) -> np.ndarray[int, np.dtype[np._IntType]]:
    """Relabel clusters so they are sorted by number of members, descending.

    Args:
        labels (np.array(int)): 1-based cluster labels
    """
    order: np.ndarray[int, np.dtype[np._IntType]] = np.argsort(np.argsort(-np.bincount(labels)))
    return 1 + order[labels]


def save_clustering_csv(
    base_dir: str, clustering_key: str, labels: Iterable[int], barcodes: Iterable[bytes]
):
    out_dir = os.path.join(base_dir, clustering_key)
    os.makedirs(out_dir, exist_ok=True)

    clusters_fn = os.path.join(out_dir, "clusters.csv")

    header = ["Barcode", "Cluster"]
    analysis_io.save_matrix_csv(clusters_fn, labels, header, barcodes)


def create_legacy_kmeans_nodes(
    f: tables.File,
    new_group_name,
    legacy_group_name: str,
    namedtuple: NamedTupleProtocol,
    clustering_key,
):
    """Soft-link a legacy-structured (CR 1.2) kmeans subgroup (dest) to a new-style (CR 1.3) subgroup (src).

    The old-style was a group called 'kmeans' with subgroups named _K.
    The new-style is a group called 'clustering' with subgroups named kmeans_K_clusters, etc.
    """
    group = f.create_group(f.root, legacy_group_name)

    cluster_type, cluster_param = parse_clustering_key(clustering_key)
    if cluster_type != CLUSTER_TYPE_KMEANS:
        return

    legacy_key = format_legacy_clustering_key(cluster_type, cluster_param)
    subgroup = f.create_group(group, legacy_key)
    for field in namedtuple._fields:
        target = f"/{new_group_name}/_{clustering_key}/{field}"

        if target in f:
            f.create_soft_link(subgroup, field, target=target)
        else:
            sys.stderr.write(
                f"Skipped soft-link of legacy dataset to {target}; node doesn't exist\n"
            )


def subselect_barcodes(clustering: CLUSTERING, bc_indices: Sequence[int]):
    """Select a subset of barcodes from a clustering object.

    Args:
        clustering (CLUSTERS namedtuple)
    """
    return CLUSTERING(
        clusters=clustering.clusters[bc_indices],
        num_clusters=clustering.num_clusters,
        cluster_score=clustering.cluster_score,
        clustering_type=clustering.clustering_type,
        global_sort_key=clustering.global_sort_key,
        description=clustering.description,
    )


def get_cluster_sizes(clustering):
    """Returns a numpy array containing cell-counts for each cluster."""
    return np.bincount(clustering.clusters)[1:]


def sort_clusterings(clusterings):
    return sorted(clusterings, key=lambda x: x.global_sort_key)


def create_clustering(
    clusters: np.ndarray[int, np.dtype[np.int64]] | Sequence[int],
    num_clusters: int,
    cluster_score: float,
    clustering_type: str | bytes,
    global_sort_key: float,
    description: str | bytes,
):
    """Create a clustering namedtuple.

    Use numpy arrays/scalars to ensure h5 compatibility
    """
    return CLUSTERING(
        clusters=np.asarray(clusters, dtype=np.int64),
        num_clusters=np.int64(num_clusters),
        cluster_score=np.float64(cluster_score),
        clustering_type=np.bytes_(clustering_type),
        global_sort_key=np.float64(global_sort_key),
        description=np.bytes_(description),
    )
