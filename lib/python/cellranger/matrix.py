#!/usr/bin/env python3
#
# Copyright (c) 2018 10X Genomics, Inc. All rights reserved.
#
from __future__ import annotations

import copy
import os.path
import pathlib
import shutil
from collections import OrderedDict
from collections.abc import Callable, Collection, Container, Iterable, Mapping, Sequence
from io import BufferedIOBase
from typing import Any, TypeVar, overload

import h5py as h5
import numpy as np
import scipy.sparse as sp_sparse
from six import ensure_binary, ensure_str

import cellranger.cr_io as cr_io
import cellranger.h5_constants as h5_constants
import cellranger.hdf5 as cr_h5
import cellranger.rna.library as rna_library
import cellranger.sparse as cr_sparse
import cellranger.utils as cr_utils
import tenkit.safe_json as tk_safe_json
from cellranger.fast_utils import (  # pylint: disable=no-name-in-module, invalid-name
    MatrixBarcodeIndex,
)
from cellranger.feature_ref import GENOME_FEATURE_TAG, FeatureDef, FeatureReference
from cellranger.wrapped_tables import tables

HDF5_COMPRESSION = "gzip"
# Number of elements per chunk. Here, 1 MiB / (12 bytes)
HDF5_CHUNK_SIZE = 80000

DEFAULT_DATA_DTYPE = "int32"

FEATURES_TSV_GZ = "features.tsv.gz"

# Normalization tag in H5
NORMALIZATION_TAG_IN_H5 = "normalized"


# some helper functions from stats
def sum_sparse_matrix(matrix, axis: int = 0) -> np.ndarray:
    """Sum a sparse matrix along an axis."""
    axis_sum = np.asarray(matrix.sum(axis=axis))  # sum along given axis
    max_dim = np.prod(axis_sum.shape)  # get the max dimension
    return axis_sum.reshape((max_dim,))  # reshape accordingly


_T = TypeVar("_T", bound=np.generic)


def top_n(array: np.ndarray[Any, np.dtype[_T]], n: int) -> Iterable[tuple[int, _T]]:
    """Retrieve the N largest elements and their positions in a numpy ndarray.

    Args:
       array (numpy.ndarray): Array
       n (int): Number of elements

    Returns:
       list of tuple of (int, x): Tuples are (original index, value).
    """
    indices = np.argpartition(array, -n)[-n:]
    indices = indices[np.argsort(array[indices])]
    return zip(indices, array[indices])


MATRIX_H5_FILETYPE = "matrix"
MATRIX = "matrix"
MATRIX_H5_VERSION_KEY = "version"
MATRIX_H5_VERSION = 2
SOFTWARE_H5_VERSION_KEY = "software_version"

# used to distinguish from user-defined attrs introduced in aggr
MATRIX_H5_BUILTIN_ATTRS = [
    h5_constants.H5_FILETYPE_KEY,
    MATRIX_H5_VERSION_KEY,
] + h5_constants.H5_METADATA_ATTRS


class NullAxisMatrixError(Exception):
    pass


def _ensure_types_match(
    original: sp_sparse.spmatrix, subset: sp_sparse.spmatrix
) -> sp_sparse.spmatrix:
    """Scipy has the fantastic habit of subsetting a matrix of one dtype down to a float64 if no elements are selected.

    This function makes the new matrix have the same type as the old matrix, and can be called after subsetting.

    Args:
        original: The original matrix
        subset: The matrix after subsetting

    Returns:
        A type correct version of the matrix
    """
    if original.dtype != subset.dtype:
        subset = subset.astype(original.dtype)
    return subset


def _save_extra_attrs(
    file: h5.File, extra_attrs: None | (Mapping[str | bytes, str | bytes | Iterable[str | bytes]])
) -> None:
    """Set optional top level attributes."""
    if extra_attrs:
        for k, v in extra_attrs.items():
            cr_h5.set_hdf5_attr(file, k, v)


def _save_sw_version(file: h5.File, sw_version: None | str | bytes) -> None:
    if sw_version:
        file.attrs[SOFTWARE_H5_VERSION_KEY] = sw_version


class CountMatrixView:
    """Supports summing a sliced CountMatrix w/o copying the whole thing."""

    def __init__(
        self,
        matrix: CountMatrix,
        feature_indices: Collection[int] | None = None,
        bc_indices: Collection[int] | None = None,
    ):
        self.feature_mask = np.ones(matrix.features_dim, dtype="bool")
        self.bc_mask = np.ones(matrix.bcs_dim, dtype="bool")
        self.matrix = matrix

        if feature_indices is not None:
            self.feature_mask.fill(False)
            self.feature_mask[np.asarray(feature_indices)] = True
        if bc_indices is not None:
            self.bc_mask.fill(False)
            self.bc_mask[np.asarray(bc_indices)] = True

        self._update_feature_ref()

    @property
    def bcs_dim(self) -> int:
        return np.count_nonzero(self.bc_mask)

    @property
    def features_dim(self) -> int:
        return np.count_nonzero(self.feature_mask)

    def _copy(self) -> CountMatrixView:
        """Return a copy of this view."""
        view = CountMatrixView(self.matrix)
        view.bc_mask = np.copy(self.bc_mask)
        view.feature_mask = np.copy(self.feature_mask)
        view._update_feature_ref()
        return view

    def view(self) -> CountMatrixView:
        """Return a copy of this view."""
        return self._copy()

    def sum(self, axis: int | None = None):
        """Sum across an axis."""
        return cr_sparse.sum_masked(self.matrix.m, self.feature_mask, self.bc_mask, axis=axis)

    @overload
    def count_ge(self, axis: None, threshold: int) -> np.uint64: ...

    @overload
    def count_ge(self, axis: int, threshold: int) -> np.ndarray[int, np.dtype[np.uint64]]: ...

    def count_ge(
        self, axis: int | None, threshold: int
    ) -> np.uint64 | np.ndarray[int, np.dtype[np.uint64]]:
        """Count number of elements >= X over an axis."""
        return cr_sparse.count_ge_masked(
            self.matrix.m, self.feature_mask, self.bc_mask, threshold, axis
        )

    def select_barcodes(self, indices: list[int] | np.ndarray):
        """Select a subset of barcodes (by index in the original matrix) and return the resulting view."""
        view = self._copy()
        mask = np.bincount(indices, minlength=len(view.bc_mask)).astype(bool)
        mask[mask > 1] = 1
        view.bc_mask &= mask
        return view

    def select_barcodes_by_seq(self, barcode_seqs: list[bytes]) -> CountMatrixView:
        indices = self.matrix.bcs_to_ints(barcode_seqs)
        return self.select_barcodes(indices)

    def select_barcodes_by_gem_group(self, gem_group: int) -> CountMatrixView:
        return self.select_barcodes_by_seq(
            [bc for bc in self.matrix.bcs if gem_group == cr_utils.split_barcode_seq(bc)[1]]
        )

    def _update_feature_ref(self) -> None:
        """Make the feature reference consistent with the feature mask."""
        indices = np.flatnonzero(self.feature_mask)
        self.feature_ref = FeatureReference(
            feature_defs=[self.matrix.feature_ref.feature_defs[i] for i in indices],
            all_tag_keys=self.matrix.feature_ref.all_tag_keys,
            target_features=self.matrix.feature_ref.target_features,
        )

    def select_features(self, indices: list[int] | np.ndarray) -> CountMatrixView:
        """Select a subset of features and return the resulting view."""
        view = self._copy()
        mask = np.bincount(indices, minlength=len(view.feature_mask)).astype(bool)
        mask[mask > 1] = 1
        view.feature_mask &= mask
        view._update_feature_ref()
        return view

    def select_features_by_genome(self, genome: str):
        """Select the subset of gene-expression features for genes in a specific genome."""
        indices = []
        for feature in self.matrix.feature_ref.feature_defs:
            if feature.feature_type == rna_library.DEFAULT_LIBRARY_TYPE:
                if feature.tags["genome"] == genome:
                    indices.append(feature.index)
        return self.select_features(indices)

    def select_features_by_genome_and_types(self, genome: str, feature_types: Container[str]):
        """Subset the features by types and genome.

        Select the subset of gene-expression features for genes in a specific genome and matching
        one of the types listed in feature_types.
        """
        indices: list[int] = []
        for feature in self.matrix.feature_ref.feature_defs:
            if feature.feature_type in feature_types:
                include_feature = True
                # Genome test only applies to GEX features
                if feature.feature_type == rna_library.DEFAULT_LIBRARY_TYPE:
                    if feature.tags["genome"] != genome:
                        include_feature = False
                if include_feature:
                    indices.append(feature.index)
        return self.select_features(indices)

    def select_features_by_types(self, feature_types: Container[str]):
        """Subset the features by type.

        Select the subset of gene-expression features for genes in a specific genome and matching
        one of the types listed in feature_types.
        """
        indices: list[int] = []
        for feature in self.matrix.feature_ref.feature_defs:
            if feature.feature_type in feature_types:
                indices.append(feature.index)
        return self.select_features(indices)

    def select_features_by_type(self, feature_type):
        """Select the subset of features with a particular feature type (e.g. "Gene Expression")."""
        return self.select_features(self.get_feature_indices_by_type(feature_type))

    def select_features_by_ids(self, feature_ids: Iterable[bytes]) -> CountMatrixView:
        return self.select_features(self.matrix.feature_ids_to_ints(feature_ids))

    def get_feature_indices_by_type(self, feature_type: str) -> list[int]:
        """Return the list of indices of features corresponding a feature type."""
        return self.matrix.feature_ref.get_indices_for_type(feature_type)

    def get_genomes(self) -> list[str]:
        """Get a list of the distinct genomes represented by gene expression features."""
        return CountMatrix._get_genomes_from_feature_ref(self.feature_ref)

    def bcs_to_ints(self, bcs: list[bytes] | set[bytes]) -> list[int]:
        # Only works when we haven't masked barcodes.
        if np.count_nonzero(self.bc_mask) != self.matrix.bcs_dim:
            raise NotImplementedError(
                "Calling bcs_to_ints on a barcode-sliced matrix view is unimplemented"
            )
        return self.matrix.bcs_to_ints(bcs)

    def ints_to_bcs(self, bc_ints: list[int] | np.ndarray) -> list[bytes]:
        if bc_ints is None or len(bc_ints) == 0:
            return []
        sliced_bc_ints = np.flatnonzero(self.bc_mask)
        orig_bc_ints = sliced_bc_ints[np.asarray(bc_ints)]
        return [self.matrix.bcs[i] for i in orig_bc_ints]

    def int_to_feature_id(self, i: int) -> bytes:
        return self.feature_ref.feature_defs[i].id

    def get_shape(self) -> tuple[int, int]:
        """Return the shape of the sliced matrix."""
        return (np.count_nonzero(self.feature_mask), np.count_nonzero(self.bc_mask))

    def get_num_nonzero(self):
        """Return the number of nonzero entries in the sliced matrix."""
        return self.count_ge(axis=None, threshold=1)

    def get_counts_per_bc(self):
        return self.sum(axis=0)


class CountMatrix:
    # pylint: disable=too-many-public-methods
    def __init__(
        self,
        feature_ref: FeatureReference,
        bcs: Collection[bytes] | Collection[str],
        matrix: sp_sparse.spmatrix,
    ):
        # Features (genes, CRISPR gRNAs, antibody barcodes, etc.)
        self.feature_ref = feature_ref
        self.features_dim = len(feature_ref.feature_defs)
        self.feature_ids_map: dict[bytes, int] = {f.id: f.index for f in feature_ref.feature_defs}

        # Cell barcodes
        if isinstance(bcs, np.ndarray) and bcs.dtype.type is np.bytes_:
            bc_array = bcs.copy()
        elif len(bcs) == 0:
            bc_array = np.array([], dtype="S", copy=False)
        else:
            max_len = max(len(bc) for bc in bcs)
            bc_array = np.fromiter(bcs, count=len(bcs), dtype=np.dtype((np.bytes_, max_len)))
        del bcs
        bc_array.flags.writeable = False
        self.bcs: np.ndarray[int, np.dtype[np.bytes_]] = bc_array
        (self.bcs_dim,) = self.bcs.shape

        self.bcs_idx = MatrixBarcodeIndex.from_raw_bytes(
            self.bcs.tobytes(), self.bcs.dtype.itemsize
        )

        self.m: sp_sparse.spmatrix = matrix
        assert self.m.shape[1] == len(self.bcs), "Barcodes must be equal to cols of matrix"

    def get_shape(self):
        """Return the shape of the sliced matrix."""
        return self.m.shape

    def get_num_nonzero(self):
        """Return the number of nonzero entries in the sliced matrix."""
        return self.m.nnz

    def view(self):
        """Return a view on this matrix."""
        return CountMatrixView(self)

    @classmethod
    def empty(cls, feature_ref: FeatureReference, bcs: Collection[bytes], dtype=DEFAULT_DATA_DTYPE):
        """Create an empty matrix."""
        matrix = sp_sparse.lil_matrix((len(feature_ref.feature_defs), len(bcs)), dtype=dtype)
        return cls(feature_ref=feature_ref, bcs=bcs, matrix=matrix)

    @staticmethod
    def from_legacy_v1_h5(h5_file: h5.File) -> CountMatrix:
        """Create a CountMatrix from a legacy h5py.File (format version 1)."""
        genome_arrays = []
        gene_id_arrays = []
        gene_name_arrays = []
        bc_idx_arrays = []
        feat_idx_arrays = []
        data_arrays = []

        # Map barcode string to column index in new matrix
        barcode_map = OrderedDict()

        # Construct a genome-concatenated matrix and FeatureReference
        for genome_idx, genome in enumerate(h5_file.keys()):
            g = h5_file[genome]

            n_genes = sum(len(x) for x in gene_id_arrays)

            # Offset the row (gene) indices by the number of genes seen so far
            feat_idx_arrays.append(g["indices"][:] + n_genes)

            # Offset the col (barcode) indices by the number of nonzero elements seen so far

            # Map barcode (column) indices to a single unique barcode space
            barcodes = g["barcodes"][:]
            for bc in barcodes:
                if bc not in barcode_map:
                    barcode_map[bc] = len(barcode_map)

            remapped_col_inds = np.fromiter(
                (barcode_map[bc] for bc in barcodes), count=len(barcodes), dtype="uint64"
            )

            indptr = g["indptr"][:]
            assert len(indptr) == 1 + len(remapped_col_inds)

            if genome_idx == 0:
                # For the first set of barcodes encountered, there should
                # be no change in their new indices.
                assert np.array_equal(remapped_col_inds, np.arange(len(indptr) - 1))

            # Convert from CSC to COO by expanding the indptr array out

            nz_elems_per_bc = np.diff(indptr)
            assert len(nz_elems_per_bc) == len(g["barcodes"])

            bc_idx = np.repeat(remapped_col_inds, nz_elems_per_bc)
            assert len(bc_idx) == len(g["indices"])
            assert len(bc_idx) == len(g["data"])

            bc_idx_arrays.append(bc_idx)
            data_arrays.append(g["data"][:])

            gene_id_arrays.append(g["genes"][:])
            gene_name_arrays.append(g["gene_names"][:])
            genome_arrays.append(np.repeat(genome, len(g["genes"])))

        genomes = np.concatenate(genome_arrays)
        gene_ids = np.concatenate(gene_id_arrays)
        gene_names = np.concatenate(gene_name_arrays)

        # Construct FeatureReference
        feature_defs = []
        for gene_id, gene_name, genome in zip(gene_ids, gene_names, genomes):
            feature_defs.append(
                FeatureDef(
                    index=len(feature_defs),
                    id=gene_id,
                    name=gene_name,
                    feature_type=rna_library.GENE_EXPRESSION_LIBRARY_TYPE,
                    tags={GENOME_FEATURE_TAG: genome},
                )
            )
        feature_ref = FeatureReference(feature_defs, [GENOME_FEATURE_TAG])

        i = np.concatenate(feat_idx_arrays)
        j = np.concatenate(bc_idx_arrays)
        data = np.concatenate(data_arrays)

        assert isinstance(barcode_map, OrderedDict)

        matrix = sp_sparse.csc_matrix((data, (i, j)), shape=(len(gene_ids), len(barcode_map)))

        return CountMatrix(feature_ref, barcode_map.keys(), matrix)

    def feature_id_to_int(self, feature_id: bytes) -> int:
        if not isinstance(feature_id, bytes):
            raise KeyError(f"feature_id {feature_id} must be bytes, but was {type(feature_id)}")
        if feature_id not in self.feature_ids_map:
            raise KeyError(f"Specified feature ID not found in matrix: {feature_id.decode()}")
        return self.feature_ids_map[feature_id]

    def feature_ids_to_ints(self, feature_ids: Iterable[bytes]) -> list[int]:
        return sorted(self.feature_id_to_int(fid) for fid in feature_ids)

    def feature_id_to_name(self, feature_id: bytes) -> str:
        idx = self.feature_id_to_int(feature_id)
        return self.feature_ref.feature_defs[idx].name

    def int_to_feature_id(self, i: int) -> bytes:
        return self.feature_ref.feature_defs[i].id

    def int_to_feature_name(self, i: int) -> str:
        return self.feature_ref.feature_defs[i].name

    def bc_to_int(self, bc: bytes) -> int:
        """Get the integer index for a barcode.

        Args:
            bc (bytes): The barcode to search for.

        Raises:
            ValueError: `barcode` was not bytes.
            KeyError: `barcode` was not found in the set.

        Returns:
            int: the barcode index.
        """
        return self.bcs_idx.bc_to_int(bc)

    def bcs_to_ints(self, bcs: Sequence[bytes], return_sorted: bool = True) -> list[int]:
        if isinstance(bcs, set):
            bcs = list(bcs)
        return self.bcs_idx.bcs_to_ints(bcs, sort=return_sorted)

    def int_to_bc(self, j: int) -> bytes:
        return self.bcs[j]

    def ints_to_bcs(self, jj: Iterable[int]) -> list[bytes]:
        return [self.int_to_bc(j) for j in jj]

    def add(self, feature_id: bytes, bc: bytes, value=1) -> None:
        """Add a count."""
        i, j = self.feature_id_to_int(feature_id), self.bc_to_int(bc)
        self.m[i, j] += value

    def get(self, feature_id: bytes, bc: bytes) -> None:
        i, j = self.feature_id_to_int(feature_id), self.bc_to_int(bc)
        return self.m[i, j]

    def merge(self, other: CountMatrix) -> None:
        """Merge this matrix with another CountMatrix.

        Works by addition, dimensions must be the same.
        """
        assert self.features_dim == other.features_dim
        # Also requires nrows to be the same
        self.m += other.m

    def sort_indices(self) -> None:
        self.tocsc()
        if not self.m.has_sorted_indices:
            self.m.sort_indices()

    def save_h5_file(
        self,
        filename,
        extra_attrs: None | (Mapping[str | bytes, str | bytes | Iterable[str | bytes]]) = None,
        sw_version=None,
    ):
        """Save this matrix to an HDF5 file, optionally with SW version."""
        with h5.File(ensure_binary(filename), "w") as f:
            f.attrs[h5_constants.H5_FILETYPE_KEY] = MATRIX_H5_FILETYPE
            f.attrs[MATRIX_H5_VERSION_KEY] = MATRIX_H5_VERSION
            # set the software version key only if it is supplied
            _save_sw_version(f, sw_version)
            _save_extra_attrs(f, extra_attrs)

            group = f.create_group(MATRIX)
            self.save_h5_group(group)

    def save_h5_group(self, group: h5.Group) -> None:
        """Save this matrix to an HDF5 (h5py) group and converts the matrix to csc format if not already."""
        self.sort_indices()

        # Save the feature reference
        feature_ref_group = group.create_group(h5_constants.H5_FEATURE_REF_ATTR)
        self.feature_ref.to_hdf5(feature_ref_group)

        # Store barcode sequences as array of ASCII strings
        cr_h5.create_hdf5_string_dataset(
            group, h5_constants.H5_BCS_ATTR, self.bcs, compression=True
        )

        for attr, dtype in h5_constants.H5_MATRIX_ATTRS.items():
            arr = np.array(getattr(self.m, attr), dtype=dtype)
            group.create_dataset(
                attr,
                data=arr,
                chunks=(HDF5_CHUNK_SIZE,),
                maxshape=(None,),
                compression=HDF5_COMPRESSION,
                shuffle=True,
            )

    @staticmethod
    def load_dims(group: h5.Group) -> tuple[int, int, int]:
        """Load the matrix shape from an HDF5 group."""
        (rows, cols) = group[h5_constants.H5_MATRIX_SHAPE_ATTR][:]
        entries = len(group[h5_constants.H5_MATRIX_DATA_ATTR])
        return (int(rows), int(cols), entries)

    @staticmethod
    def load_dims_from_h5(filename) -> tuple[int, int, int]:
        """Load the matrix shape from an HDF5 file."""
        return CountMatrix.load_dims_from_h5_file_handle(ensure_binary(filename))

    @staticmethod
    def _load_dims_from_legacy_v1_h5_handle(f: h5.File) -> tuple[int, int, int]:
        # legacy format
        genomes = f.keys()
        num_nonzero_entries = 0
        num_gene_ids = 0
        barcodes = set()
        for genome in genomes:
            g = f[genome]
            num_nonzero_entries += len(g["data"])
            num_gene_ids += len(g["genes"])
            barcodes.update(g["barcodes"])
        return (num_gene_ids, len(barcodes), num_nonzero_entries)

    @staticmethod
    def _load_dims_from_legacy_v1_h5(filename) -> tuple[int, int, int]:
        """Load the matrix shape from a legacy h5py.File (format version 1)."""
        with h5.File(ensure_binary(filename), "r") as f:
            return CountMatrix._load_dims_from_legacy_v1_h5_handle(f)

    @staticmethod
    def get_anndata_mem_gb_from_matrix_h5(
        filename: str | bytes,
    ) -> float:
        """Estimate memory usage of anndata object from the matrix."""
        num_features, num_bcs, nonzero_entries = CountMatrix.load_dims_from_h5(filename)
        return (
            num_features * h5_constants.MEM_BYTES_PER_MATRIX_FEATURE_H5AD
            + num_bcs * h5_constants.MEM_BYTES_PER_MATRIX_BARCODE_H5AD
            + nonzero_entries * h5_constants.MEM_BYTES_PER_MATRIX_NNZ_H5AD
            + h5_constants.MEM_BYTES_CONSTANT_H5AD
        ) / 1024**3

    @staticmethod
    def get_mem_gb_from_matrix_dim(
        num_barcodes: int,
        nonzero_entries: int,
        scale: float = h5_constants.MATRIX_MEM_GB_MULTIPLIER,
        ceil: bool = True,
    ) -> float:
        """Estimate memory usage of loading a matrix."""
        matrix_mem_gb = float(nonzero_entries) / h5_constants.NUM_MATRIX_ENTRIES_PER_MEM_GB
        # We store a list and a dict of the whitelist. Based on empirical obs.
        matrix_mem_gb += float(num_barcodes) / h5_constants.NUM_MATRIX_BARCODES_PER_MEM_GB
        return scale * np.ceil(matrix_mem_gb) if ceil else scale * matrix_mem_gb

    @staticmethod
    def get_mem_gb_from_group(
        group: h5.Group,
        scale: float = h5_constants.MATRIX_MEM_GB_MULTIPLIER,
        ceil: bool = True,
    ) -> float:
        """Estimate memory usage from an HDF5 group."""
        _, num_bcs, nonzero_entries = CountMatrix.load_dims(group)
        return CountMatrix.get_mem_gb_from_matrix_dim(num_bcs, nonzero_entries, scale, ceil)

    @staticmethod
    def get_mem_gb_from_matrix_h5(
        filename: str | bytes,
        scale: float = h5_constants.MATRIX_MEM_GB_MULTIPLIER,
        ceil: bool = True,
    ) -> float:
        """Estimate memory usage from an HDF5 file."""
        return CountMatrix.get_mem_gb_from_matrix_h5_file_handle(
            ensure_binary(filename), scale, ceil
        )

    @staticmethod
    def get_mem_gb_from_matrix_h5_file_handle(
        file_handle: str | bytes | BufferedIOBase,
        scale: float = h5_constants.MATRIX_MEM_GB_MULTIPLIER,
        ceil: bool = True,
    ) -> float:
        """Estimate memory usage from an HDF5 file."""
        _, num_bcs, nonzero_entries = CountMatrix.load_dims_from_h5_file_handle(file_handle)
        return CountMatrix.get_mem_gb_from_matrix_dim(num_bcs, nonzero_entries, scale, ceil)

    @staticmethod
    def load_dims_from_h5_file_handle(
        file_handle: str | bytes | BufferedIOBase,
    ) -> tuple[int, int, int]:
        """Load the matrix shape from an HDF5 file."""
        with h5.File(file_handle, "r") as f:
            h5_version = CountMatrix._get_format_version_from_handle(f)
            if h5_version == 1:
                return CountMatrix._load_dims_from_legacy_v1_h5_handle(f)
            else:
                return CountMatrix.load_dims(f[MATRIX])

    @staticmethod
    def load_bcs_from_h5_file_handle(file_handle) -> list[bytes]:
        """Load just the barcode sequences from an HDF5 file hadle."""
        with h5.File(file_handle, "r") as f:
            h5_version = CountMatrix._get_format_version_from_handle(f)
            if h5_version == 1:
                return CountMatrix._load_bcs_from_legacy_v1_h5_file_handle(f)
            else:
                return CountMatrix.load_bcs_from_h5_group(f[MATRIX])

    @staticmethod
    def _load_indptr_from_matrix_group(group: h5.Group):
        return group[h5_constants.H5_MATRIX_INDPTR_ATTR][:]

    @staticmethod
    def load_indptr_from_file(filename):
        fn, version = CountMatrix._validate_h5_file(filename)
        if version < MATRIX_H5_VERSION:
            raise ValueError(
                "Matrix HDF5 file format version (%d) is an older version that is no longer supported."
                % version
            )
        with h5.File(fn, "r") as f:
            mat_group = f[MATRIX]
            return CountMatrix._load_indptr_from_matrix_group(mat_group)

    @classmethod
    def load(cls, group: h5.Group) -> CountMatrix:
        """Load from an HDF5 group."""
        feature_ref = CountMatrix.load_feature_ref_from_h5_group(group)
        bcs = cls.load_bcs_from_h5_group(group)

        shape = group[h5_constants.H5_MATRIX_SHAPE_ATTR][:]
        data = group[h5_constants.H5_MATRIX_DATA_ATTR][:]
        indices = group[h5_constants.H5_MATRIX_INDICES_ATTR][:]
        indptr = CountMatrix._load_indptr_from_matrix_group(group)

        # Check to make sure indptr increases monotonically (to catch overflow bugs)
        assert np.all(np.diff(indptr) >= 0)

        matrix = sp_sparse.csc_matrix((data, indices, indptr), shape=shape)

        return cls(feature_ref=feature_ref, bcs=bcs, matrix=matrix)

    @classmethod
    def load_columns_from_file(cls, group: h5.Group, start: int, end: int) -> CountMatrix:
        """Load from an HDF5 group only a subset of the columns.

        Very similar to CountMatrix.load().

        Args:
            group: The matrix H5 group
            start: the first column to select
            end: the column to select up to.

        Returns:
            A CountMatrix with the relevant columns
        """
        feature_ref = CountMatrix.load_feature_ref_from_h5_group(group)
        bcs = cls.load_bcs_from_h5_group(group)
        n_rows, n_cols, _ = cls.load_dims(group)
        if start < 0 or end > n_cols or start > end:
            raise ValueError(
                f"The column range you've specified is invalid.  Start Column = {start}, End = {end}, Total Columns = {n_cols}"
            )
        indptr = group[h5_constants.H5_MATRIX_INDPTR_ATTR][start : (end + 1)]
        col_start = indptr[0]
        col_end = indptr[-1]
        shape = np.array([n_rows, (end - start)])
        data = group[h5_constants.H5_MATRIX_DATA_ATTR][col_start:col_end]
        indices = group[h5_constants.H5_MATRIX_INDICES_ATTR][col_start:col_end]
        indptr = indptr - indptr[0]
        # Check to make sure indptr increases monotonically (to catch overflow bugs)
        assert np.all(np.diff(indptr) >= 0)
        matrix = sp_sparse.csc_matrix((data, indices, indptr), shape=shape)
        return cls(feature_ref=feature_ref, bcs=bcs[start:end], matrix=matrix)

    @staticmethod
    def load_bcs_from_h5_group(group: h5.Group) -> list[bytes]:
        """Load just the barcode sequences from an h5 group."""
        if group[h5_constants.H5_BCS_ATTR].shape is not None:
            return list(group[h5_constants.H5_BCS_ATTR][:])
        return []

    @staticmethod
    def load_bcs_from_h5(filename) -> list[bytes]:
        """Load just the barcode sequences from an HDF5 group."""
        filename = ensure_binary(filename)
        h5_version = CountMatrix.get_format_version_from_h5(filename)
        if h5_version == 1:
            return CountMatrix._load_bcs_from_legacy_v1_h5(filename)
        else:
            with h5.File(filename, "r") as f:
                return CountMatrix.load_bcs_from_h5_group(f[MATRIX])

    @staticmethod
    def _load_bcs_from_legacy_v1_h5(filename) -> list[bytes]:
        """Load just the barcode sequences from a legacy h5py.File (format version 1)."""
        with h5.File(ensure_binary(filename), "r") as f:
            return CountMatrix._load_bcs_from_legacy_v1_h5_file_handle(f)

    @staticmethod
    def _load_bcs_from_legacy_v1_h5_file_handle(file_handle) -> list[bytes]:
        """Load just the barcode sequences from a legacy h5py.File (format version 1)."""
        genomes = file_handle.keys()
        barcodes: set[bytes] = set()
        for genome in genomes:
            group = file_handle[genome]
            barcodes.update(group["barcodes"])
        return list(barcodes)

    @staticmethod
    def load_bcs_from_h5_file(filename) -> list[bytes]:
        with h5.File(ensure_binary(filename), "r") as f:
            if (
                h5_constants.H5_FILETYPE_KEY not in f.attrs
                or f.attrs[h5_constants.H5_FILETYPE_KEY] != MATRIX_H5_FILETYPE
            ):
                raise ValueError("HDF5 file is not a valid matrix HDF5 file.")

            if MATRIX_H5_VERSION_KEY in f.attrs:
                version = f.attrs[MATRIX_H5_VERSION_KEY]
            else:
                version = 1

            if version > MATRIX_H5_VERSION:
                raise ValueError(
                    "Matrix HDF5 file format version (%d) is a newer version that is not supported by this version of the software."
                    % version
                )
            if version < MATRIX_H5_VERSION:
                raise ValueError(
                    "Matrix HDF5 file format version (%d) is an older version that is no longer supported."
                    % version
                )

            if "matrix" not in f.keys():
                raise ValueError('Could not find the "matrix" group inside the matrix HDF5 file.')

            return CountMatrix.load_bcs_from_h5_group(f["matrix"])

    @staticmethod
    def load_library_types_from_h5_file(filename) -> set[str]:
        """Return a set of all library types defined in the Feature Reference."""
        with h5.File(ensure_binary(filename), "r") as f:
            version = CountMatrix._get_format_version_from_handle(f)
            if version < MATRIX_H5_VERSION:
                # Only GEX supported, check for any data and return GEX if any exists
                (gene_count, _, _) = CountMatrix._load_dims_from_legacy_v1_h5_handle(f)
                if gene_count > 0:
                    return {rna_library.GENE_EXPRESSION_LIBRARY_TYPE}
                else:
                    return set()
            else:
                feature_ref = CountMatrix.load_feature_ref_from_h5_group(f[MATRIX])
                return {f.feature_type for f in feature_ref.feature_defs}

    def get_library_types(self) -> set[str]:
        """Get the list of feature types."""
        return {f.feature_type for f in self.feature_ref.feature_defs}

    @staticmethod
    def load_feature_ref_from_h5_group(group: h5.Group) -> FeatureReference:
        """Load just the FeatureRef from an h5py.Group."""
        feature_group = group[h5_constants.H5_FEATURE_REF_ATTR]
        return FeatureReference.from_hdf5(feature_group)

    @staticmethod
    def load_feature_ref_from_h5_file(filename) -> FeatureReference:
        """Load just the FeatureRef from a matrix HDF5 file."""
        with h5.File(ensure_binary(filename), "r") as f:
            version = CountMatrix._get_format_version_from_handle(f)
            if version < MATRIX_H5_VERSION:
                raise OSError("Direct Feature Ref reading not supported for older H5 files.")
            else:
                return CountMatrix.load_feature_ref_from_h5_group(f[MATRIX])

    def tolil(self):
        if type(self.m) is not sp_sparse.lil_matrix:
            self.m = self.m.tolil()

    def tocoo(self):
        if type(self.m) is not sp_sparse.coo_matrix:
            self.m = self.m.tocoo()

    def tocsc(self):
        # Convert from lil to csc matrix for efficiency when analyzing data
        if type(self.m) is not sp_sparse.csc_matrix:
            self.m = self.m.tocsc()

    def tocsr(self):
        """Convert to a csr matrix if not already.

        Returns:
            None, mutates in place
        """
        # Convert from lil to csc matrix for efficiency when analyzing data
        if type(self.m) is not sp_sparse.csr_matrix:
            self.m = self.m.tocsr()

    def select_axes_above_threshold(
        self, threshold: int = 0
    ) -> tuple[CountMatrix, np.ndarray, np.ndarray]:
        """Select axes with sums greater than the threshold value.

        Returns:
            (CountMatrix, np.array of int, np.array of int):
                New count matrix, non-zero bc indices, feat indices
        """
        new_mat = copy.deepcopy(self)

        nonzero_bcs = np.flatnonzero(new_mat.get_counts_per_bc() > threshold)
        if self.bcs_dim > len(nonzero_bcs):
            new_mat = new_mat.select_barcodes(nonzero_bcs)

        nonzero_features = np.flatnonzero(new_mat.get_counts_per_feature() > threshold)
        if new_mat.features_dim > len(nonzero_features):
            new_mat = new_mat.select_features(nonzero_features)

        if len(nonzero_bcs) == 0 or len(nonzero_features) == 0:
            raise NullAxisMatrixError()

        return new_mat, nonzero_bcs, nonzero_features

    def select_axis_above_threshold(self, axis: int, threshold: int = 0) -> CountMatrix:
        """Select axis given a threshold of features per barcode or total counts per feature.

        Args:
            axis (int): 0 for rows (features) and 1 for columns (barcodes).

        Returns:
            (CountMatrix):
                New count matrix
        """
        if axis == 0:
            counts = self.get_counts_per_feature()
            indices = np.flatnonzero(counts > threshold)
            return self.select_features(indices)

        elif axis == 1:
            counts = self.get_counts_per_bc()
            indices = np.flatnonzero(counts > threshold)
            return self.select_barcodes(indices)

        else:
            raise ValueError("axis out of range")

    def select_nonzero_axes(self) -> tuple[CountMatrix, np.ndarray, np.ndarray]:
        """Select axes with nonzero sums.

        Returns:
            (CountMatrix, np.array of int, np.array of int):
                New count matrix, non-zero bc indices, feat indices
        """
        return self.select_axes_above_threshold(0)

    def select_nonzero_axis(self, axis: int) -> tuple[CountMatrix, np.ndarray]:
        """Select axis with nonzero sums.

        Args:
            axis (int): 0 for rows (features) and 1 for columns (barcodes).

        Returns:
            (CountMatrix, np.array of int, np.array of int):
                New count matrix, selected indices
        """
        if axis == 0:
            counts = self.get_counts_per_feature()
            indices = np.flatnonzero(counts > 0)
            return self.select_features(indices), indices

        elif axis == 1:
            counts = self.get_counts_per_bc()
            indices = np.flatnonzero(counts > 0)
            return self.select_barcodes(indices), indices

        else:
            raise ValueError("axis out of range")

    def select_barcodes(self, indices: Sequence[int]) -> CountMatrix:
        """Select a subset of barcodes and return the resulting CountMatrix."""
        submat = self.m[:, indices]
        submat = _ensure_types_match(self.m, submat)
        return CountMatrix(
            feature_ref=self.feature_ref,
            bcs=[self.bcs[i] for i in indices],
            matrix=submat,
        )

    def select_barcodes_by_seq(self, barcode_seqs: Sequence[bytes]) -> CountMatrix:
        indices = self.bcs_to_ints(barcode_seqs, False)
        return self.select_barcodes(indices)

    def select_barcodes_by_gem_group(self, gem_group: int) -> CountMatrix:
        return self.select_barcodes_by_seq(
            [bc for bc in self.bcs if gem_group == cr_utils.get_gem_group_from_barcode(bc)]
        )

    def select_features(self, indices: Iterable[int]) -> CountMatrix:
        """Select a subset of features and return the resulting matrix.

        We also update FeatureDefs to keep their indices consistent with their new position.
        """
        feature_ref = self.feature_ref.select_features(indices)
        submat = self.m[indices, :]
        submat = _ensure_types_match(self.m, submat)
        return CountMatrix(feature_ref=feature_ref, bcs=self.bcs, matrix=submat)

    def select_features_by_ids(self, feature_ids: Iterable[bytes]) -> CountMatrix:
        return self.select_features(self.feature_ids_to_ints(feature_ids))

    def remove_genes_not_on_list(self, gene_indices_to_keep: Iterable[int]) -> CountMatrix:
        """Removes all features that are GEX and not on this list, keeping all others.

        Used to subset the matrix down in the targeted assay

        Args:
            gene_indices_to_keep: list of indices of the GEX features to keep

        Returns:
            A copy of this matrix subset to the GEX feature indices requested and all other
            feature types.
        """
        indices = list(gene_indices_to_keep)
        for feature in self.feature_ref.feature_defs:
            if feature.feature_type != rna_library.GENE_EXPRESSION_LIBRARY_TYPE:
                indices.append(feature.index)
        return self.select_features(indices)

    def select_features_by_genome(self, genome: str) -> CountMatrix:
        """Select the subset of gene-expression features for genes in a specific genome."""
        indices = []
        for feature in self.feature_ref.feature_defs:
            if feature.feature_type == rna_library.DEFAULT_LIBRARY_TYPE:
                if feature.tags["genome"] == genome:
                    indices.append(feature.index)
        return self.select_features(indices)

    def select_features_by_types(self, feature_types: Container[str]) -> CountMatrix:
        """Select the subset of gene-expression features by genome and feature type."""
        indices = []
        for feature in self.feature_ref.feature_defs:
            if feature.feature_type in feature_types:
                indices.append(feature.index)
        return self.select_features(indices)

    def select_features_by_type(self, feature_type: str) -> CountMatrix:
        """Select the subset of features with a particular feature type (e.g. "Gene Expression")."""
        indices = []
        for feature in self.feature_ref.feature_defs:
            if feature.feature_type == feature_type:
                indices.append(feature.index)
        return self.select_features(indices)

    def select_features_by_type_and_tag(self, feature_type: str, tag_type: str) -> CountMatrix:
        """Select the subset of features with a particular feature type (e.g. "Antibody Capture") and tag (e.g. "Hashtag")."""
        indices = []
        for feature in self.feature_ref.feature_defs:
            if feature.feature_type == feature_type and tag_type in feature.tags:
                indices.append(feature.index)
        return self.select_features(indices)

    def get_feature_ids_by_type(self, feature_type: str) -> list[bytes]:
        """Return a list of feature ids of a particular feature type (e.g. "Gene Expression")."""
        return self.feature_ref.get_feature_ids_by_type(feature_type)

    def get_count_of_feature_type(self, feature_type: str) -> int:
        """Count how many features in the matrix are of a given type.

        (e.g. "Gene Expression")
        """
        return self.feature_ref.get_count_of_feature_type(feature_type)

    def get_count_of_feature_types(self) -> dict[str, int]:
        """Get count of each feature by type."""
        return {ft: self.get_count_of_feature_type(ft) for ft in self.get_library_types()}

    @staticmethod
    def _get_genomes_from_feature_ref(feature_ref: FeatureReference) -> list[str]:
        """Get a list of the distinct genomes represented by gene expression features."""
        return feature_ref.get_genomes(feature_type=rna_library.DEFAULT_LIBRARY_TYPE)

    def get_genomes(self):
        """Get a list of the distinct genomes represented by gene expression features."""
        return CountMatrix._get_genomes_from_feature_ref(self.feature_ref)

    @staticmethod
    def get_genomes_from_h5(filename: str | bytes) -> list[str]:
        """Get a list of the distinct genomes from a matrix HDF5 file."""
        filename = ensure_binary(filename)
        h5_version = CountMatrix.get_format_version_from_h5(filename)
        if h5_version == 1:
            return CountMatrix._get_genomes_from_legacy_v1_h5(filename)
        else:
            with h5.File(filename, "r") as f:
                feature_ref = CountMatrix.load_feature_ref_from_h5_group(f[MATRIX])
                return CountMatrix._get_genomes_from_feature_ref(feature_ref)

    @staticmethod
    def _get_genomes_from_legacy_v1_h5(filename):
        """Get a list of the distinct genomes from a legacy h5py.File (format version 1)."""
        with h5.File(ensure_binary(filename), "r") as f:
            return list(f.keys())

    def get_numfeatures_per_bc(self) -> np.ndarray[int, np.dtype[np.int_]]:
        self.m.eliminate_zeros()
        return self.m.getnnz(axis=0)

    def get_counts_per_bc(self) -> np.ndarray[int, np.dtype[np.int_]]:
        return sum_sparse_matrix(self.m, axis=0)

    def get_counts_per_barcode_for_genome(
        self, genome: str, feature_type: str | None = None
    ) -> np.ndarray[int, np.dtype[np.int_]]:
        """Sum the count matrix across feature rows with a given genome tag.

        The feature reference
        must contain a 'genome' tag. If feature_type is not null filter on it as well.
        """
        assert "genome" in self.feature_ref.all_tag_keys, "feature reference missing 'genome' tag"
        if feature_type:
            indices = [
                i
                for i, fdef in enumerate(self.feature_ref.feature_defs)
                if fdef.tags["genome"] == genome and fdef.feature_type == feature_type
            ]
        else:
            indices = [
                i
                for i, fdef in enumerate(self.feature_ref.feature_defs)
                if fdef.tags["genome"] == genome
            ]
        if indices:
            view = CountMatrixView(self, feature_indices=indices, bc_indices=None)
            return view.sum(axis=0)
        return np.array([], dtype=self.m.dtype)

    def get_counts_per_feature(self) -> np.ndarray[int, np.dtype[np.int_]]:
        return sum_sparse_matrix(self.m, axis=1)

    def get_frac_counts_per_feature(self) -> np.ndarray[int, np.dtype[np.int_]]:
        total_umi = sum(sum_sparse_matrix(self.m))
        return sum_sparse_matrix(self.m, axis=1) / total_umi

    def get_mean_and_var_per_feature(
        self,
    ) -> tuple[np.ndarray[int, np.dtype[np.float64]], np.ndarray[int, np.dtype[np.float64]]]:
        """Calculate the mean and variance on the sparse matrix efficiently.

        :return: a tuple with numpy arrays for mean and var
        """
        assert isinstance(self.m, sp_sparse.csc_matrix)
        mean_per_feature = self.m.mean(axis=1)
        second_moment = self.m.copy()
        second_moment = second_moment.power(2.0)
        var_per_feature = second_moment.sum(axis=1) / second_moment.shape[1] - np.power(
            mean_per_feature, 2.0
        )
        var_per_feature = np.asarray(var_per_feature)
        return (mean_per_feature, var_per_feature)

    def get_subselected_counts(
        self, list_feature_ids=None, list_barcodes=None, log_transform=False, library_type=None
    ) -> np.ndarray[int, np.dtype[np.int_]]:
        """Get counts per barcode, sliced various ways.

        - subset by list of feature IDs
        - subset by list of barcodes
        - subset by library_type
        """
        subselect_matrix = self

        if library_type is not None:
            assert (
                library_type in rna_library.RECOGNIZED_FEATURE_TYPES
            ), "library_type not recognized"
            subselect_matrix = subselect_matrix.select_features_by_type(library_type)

        if list_feature_ids is not None:
            subselect_matrix = subselect_matrix.select_features_by_ids(list_feature_ids)

        if list_barcodes is not None:
            subselect_matrix = subselect_matrix.select_barcodes_by_seq(list_barcodes)

        counts_feature = subselect_matrix.get_counts_per_bc()

        if log_transform:
            return np.log10(1.0 + counts_feature)
        return counts_feature

    def get_numbcs_per_feature(self) -> np.ndarray[int, np.dtype[np.int_]]:
        return sum_sparse_matrix(self.m > 0, axis=1)

    def get_top_bcs(self, cutoff: int) -> np.ndarray:
        reads_per_bc = self.get_counts_per_bc()
        index = max(0, min(reads_per_bc.size, cutoff) - 1)
        value = sorted(reads_per_bc, reverse=True)[index]
        return np.nonzero(reads_per_bc >= value)[0]

    def save_mex(
        self,
        base_dir: str,
        save_features_func: Callable[[FeatureReference, str, bool], None],
        metadata: dict | None = None,
        compress: bool = True,
    ):
        """Save in Matrix Market Exchange format.

        Note:
            This operation modifies the matrix by
            converting to a coordinate representation by calling scipy.sparse.csc_matrix.tocoo().

        Args:
          base_dir (str): Path to directory to write files in.
          save_features_func (func): Func that takes (FeatureReference, base_dir, compress) and writes
                                     a file describing the features.
          metadata (dict): Optional metadata to encode into the comments as JSON.
          compress: Whether to compress output files with gzip.
        """
        self.sort_indices()
        self.tocoo()

        os.makedirs(base_dir, exist_ok=True)

        out_matrix_fn = os.path.join(base_dir, "matrix.mtx")
        out_barcodes_fn = os.path.join(base_dir, "barcodes.tsv")
        if compress:
            out_matrix_fn += ".gz"
            out_barcodes_fn += ".gz"

        if self.m.dtype in ["uint32", "int32", "uint64", "int64"]:
            field = "integer"
            fmt = b"%i %i %i\n"
        elif self.m.dtype in ["float", "double"]:
            field = "real"
            fmt = b"%i %i %15g\n"
        else:
            raise ValueError(f"Unsupported data type for the matrix: {self.m.dtype}")

        assert isinstance(self.m, sp_sparse.coo.coo_matrix)

        rows, cols = self.m.shape
        # Header fields in the file
        rep = "coordinate"
        symmetry = "general"

        metadata = metadata or {}
        metadata["format_version"] = MATRIX_H5_VERSION

        comment = b"metadata_json: " + tk_safe_json.safe_jsonify(metadata).encode()

        with cr_io.open_maybe_gzip(out_matrix_fn, "wb") as stream:
            # write initial header line
            stream.write(f"%%MatrixMarket matrix {rep} {field} {symmetry}\n".encode())

            # write comments
            for line in comment.split(b"\n"):
                stream.write(b"%")
                stream.write(line)
                stream.write(b"\n")

            # write shape spec
            stream.write(b"%i %i %i\n" % (rows, cols, self.m.nnz))
            # write row, col, val in 1-based indexing
            for r, c, d in zip(self.m.row + 1, self.m.col + 1, self.m.data):
                stream.write(fmt % (r, c, d))

        # both GEX and ATAC provide an implementation of this in respective feature_ref.py
        save_features_func(self.feature_ref, base_dir, compress=compress)

        with cr_io.open_maybe_gzip(out_barcodes_fn, "wb") as f:
            for bc in self.bcs:
                f.write(bc + b"\n")

    @staticmethod
    def _get_format_version_from_handle(ofile: h5.File) -> int:
        if MATRIX_H5_VERSION_KEY in ofile.attrs:
            version = ofile.attrs[MATRIX_H5_VERSION_KEY]
        else:
            version = 1
        return version

    @staticmethod
    def get_format_version_from_h5(filename: bytes | str):
        with h5.File(ensure_binary(filename), "r") as f:
            return CountMatrix._get_format_version_from_handle(f)

    @staticmethod
    def _validate_h5_file(filename):
        if isinstance(filename, pathlib.PosixPath):
            fn = filename
        else:
            fn = ensure_binary(filename)
        with h5.File(fn, "r") as f:
            if (
                h5_constants.H5_FILETYPE_KEY not in f.attrs
                or ensure_str(f.attrs[h5_constants.H5_FILETYPE_KEY]) != MATRIX_H5_FILETYPE
            ):
                raise ValueError("HDF5 file is not a valid matrix HDF5 file.")
            version = CountMatrix._get_format_version_from_handle(f)
            if version > MATRIX_H5_VERSION:
                raise ValueError(
                    "Matrix HDF5 file format version (%d) is a newer version that is not supported by this version of the software."
                    % version
                )
            if version >= MATRIX_H5_VERSION and MATRIX not in f.keys():
                raise ValueError('Could not find the "matrix" group inside the matrix HDF5 file.')
        return fn, version

    @staticmethod
    def load_h5_file(filename, col_start: int | None = None, col_end: int | None = None):
        """Load a matrix H5 file, optionally subsetting down to a particular range of columns if requests.

        Args:
            filename: The name of the H5 file
            col_start: (Optional) The column to select
            col_end: (Optional) End of column select range

        Returns:
            Instance of a CountMatrix

        """
        fn, version = CountMatrix._validate_h5_file(filename)

        with h5.File(fn, "r") as f:
            if version < MATRIX_H5_VERSION:
                if col_start is not None or col_end is not None:
                    raise ValueError(
                        "Subsetting columns when loading legacy H5 files is not supported."
                    )
                # raise ValueError('Matrix HDF5 file format version (%d) is an older version that is no longer supported.' % version)
                return CountMatrix.from_legacy_v1_h5(f)
            if col_start is not None or col_end is not None:
                if col_start is None or col_end is None:
                    raise ValueError("Both or neither argument col_start/col_end must be provided.")
                return CountMatrix.load_columns_from_file(f[MATRIX], col_start, col_end)
            else:
                return CountMatrix.load(f[MATRIX])

    @staticmethod
    def count_cells_from_h5(filename):
        _, bcs, _ = CountMatrix.load_dims_from_h5(filename)
        return bcs

    @staticmethod
    def load_chemistry_from_h5(filename):
        with tables.open_file(ensure_str(filename), "r") as f:
            try:
                chemistry = f.get_node_attr("/", h5_constants.H5_CHEMISTRY_DESC_KEY)
            except AttributeError:
                chemistry = "Unknown"
        return chemistry

    def filter_barcodes(self, bcs_per_genome: dict[Any, Iterable[bytes]]) -> CountMatrix:
        """Return CountMatrix containing only the specified barcodes.

        Args:
            bcs_per_genome (dict of str to list): Maps genome to cell-associated barcodes.

        Returns:
            CountMatrix w/ the specified barcodes.
        """
        # Union all the cell-associated barcodes
        bcs = set()
        for x in bcs_per_genome.values():
            bcs |= x
        bcs = list(sorted(bcs))
        return self.select_barcodes_by_seq(bcs)

    @staticmethod
    def h5_path(base_path):
        return os.path.join(base_path, "hdf5", "matrices.hdf5")

    def normalise_library_type_by_control_features(
        self,
        library_type: str,
        scale_factor: int | float,
        control_feature_names: list[str] | list[bytes],
    ) -> int:
        """Normalise the features of a library type by control features in the matrix.

        scale by the scale_factor.

        The features in the library are all
        1. divided by 1 plus the sum of control features in each barcode
        2. multiplied by scale factor
        3. floored to be an integer

        If the control features provided are not in the matrix, the division is by
        the sum of all features in the library provided.
        The features that are normalized have NORMALIZATION_TAG_IN_H5 set to TRUE.
        If this tag was not in the H5, all other features have the tag set to FALSE.
        If it was present in the H5, its values for other features is unchanged.

        Args:
            library_type (str): string showing library type
            scale_factor (Union[int, float]): factor to scale by
            control_feature_names (list[Union[str, bytes]]): feature names that are used as control

        Returns:
            int:  Number of entries clipped due to overflow
        """
        number_of_overflow_entries = 0
        assert library_type in self.get_library_types(), (
            f"Library type to normalise not in matrix. Library type passed in {library_type}. "
            + f"Library types in matrix {self.get_library_types()}"
        )

        # If the control features exist, then size factor using those.
        # If not use the total feature count in library as control
        control_features_in_matrix = [x for x in control_feature_names if x in self.feature_ids_map]
        if control_features_in_matrix:
            size_factor = (
                self.view().select_features_by_ids(control_features_in_matrix).get_counts_per_bc()
                + 1
            )
        else:
            size_factor = self.view().select_features_by_type(library_type).get_counts_per_bc() + 1

        # Get row-indices of the library in the matrix
        indices_in_library = sorted(
            self.feature_ids_to_ints(self.get_feature_ids_by_type(library_type))
        )
        # Rescale the rows corresponding to the  by a size factor and upscale by scale factor
        # enables a de facto fixed point representation. Account for the potential of overflow
        # by first computing the matrix in floats and then converting to int.
        max_value = np.iinfo(self.view().matrix.m.dtype).max
        tmp_normalised_matrix = self.m[indices_in_library, :].dot(
            sp_sparse.diags(scale_factor / size_factor).astype(np.float64)
        )
        # Clip at max_value - 1 because we often do log1p where the 1 plus will lead to overflow
        number_of_overflow_entries = (
            (tmp_normalised_matrix.data < 0) | (tmp_normalised_matrix.data > max_value - 1)
        ).sum()
        tmp_normalised_matrix.data = np.clip(tmp_normalised_matrix.data, 0, max_value - 1)
        self.m[indices_in_library, :] = tmp_normalised_matrix.astype(self.view().matrix.m.dtype)

        # Update normalize tag in the feature ref of the matrix
        features_normalized_set = set(self.get_feature_ids_by_type(library_type))
        features_label_dict = {x: "TRUE" for x in features_normalized_set}
        if NORMALIZATION_TAG_IN_H5 not in self.feature_ref.all_tag_keys:
            self.feature_ref.add_tag(NORMALIZATION_TAG_IN_H5, features_label_dict, "FALSE")
        else:
            self.feature_ref.update_tag(NORMALIZATION_TAG_IN_H5, features_label_dict)

        return number_of_overflow_entries


def merge_matrices(h5_filenames: list[str]) -> CountMatrix | None:
    """Merge multiple matrices into a single matrix."""
    matrix = None
    for h5_filename in h5_filenames:
        if matrix is None:
            matrix = CountMatrix.load_h5_file(h5_filename)
        else:
            matrix.merge(CountMatrix.load_h5_file(h5_filename))
    if matrix is not None:
        matrix.tocsc()
    return matrix


def create_merged_matrix_from_col_concat(
    in_h5_filenames, out_h5_filename, extra_attrs=None, sw_version=None
):
    """Merge several h5 files into one larger matrix file.

    An efficient method for doing column concatenation of H5 files.  Assumes the input files all have the same
    features/barcodes in their matrix, and that each only contains a non-overlapping subset of the columns.  Strategy
    is to append the column data directly, rather than merging arbitrary non-distinct columns (As is done in
    merge_matrices).

    Args:
        in_h5_filenames: A list of filenames to column distinct sets of a larger matrix
        out_h5_filename: The desired output filename
        extra_attrs: Similar to the argument to save_h5_file
        sw_version: A version of the software to add into the attributes.

    Returns:
        Nothing
    """
    # Quick sanity check that the dimensions are the same (and hope that implies barcodes/features are the same)
    # Dimensions are tuple of row/column/nnz
    dimensions = [CountMatrix.load_dims_from_h5(f) for f in in_h5_filenames]
    total_nnz = sum(x[2] for x in dimensions)
    assert np.all(
        x[0] == dimensions[0][0] for x in dimensions
    ), "Not all row dimensions were the same"
    assert np.all(
        x[1] == dimensions[0][1] for x in dimensions
    ), "Not all col dimensions were the same"

    # Load the indptrs and figure out which columns each file has
    fn_start_ends = []
    for fn in in_h5_filenames:
        indptr = CountMatrix.load_indptr_from_file(fn)
        col_spans = indptr[1:] - indptr[:-1]
        active_cols = np.nonzero(col_spans)
        # TODO: Double pass a bit annoying here, get min/max at once
        start = np.min(active_cols)
        end = np.max(active_cols)
        fn_start_ends.append((fn, start, end))
    fn_start_ends.sort(key=lambda x: x[1])
    # Now check they don't overlap and that we have independent blocks of cols
    last_end = -1
    last_start = -1
    for _, start, end in fn_start_ends:
        if start < last_end or end <= last_end or start < last_start:
            raise ValueError("H5 files to be concatenated were not unique sets of columns")
        last_end = end
        last_start = start

    # Columnwise concatenate into a new file
    start_file = fn_start_ends[0][0]
    shutil.copyfile(start_file, out_h5_filename)
    to_append = (h5_constants.H5_MATRIX_DATA_ATTR, h5_constants.H5_MATRIX_INDICES_ATTR)
    with h5.File(out_h5_filename, "a") as outfile:
        if sw_version:
            _save_sw_version(outfile, sw_version)
        if extra_attrs:
            _save_extra_attrs(outfile, extra_attrs)
        matrix = outfile[MATRIX]
        ind_ptr = CountMatrix._load_indptr_from_matrix_group(matrix).astype(np.uint64)
        end = len(matrix[to_append[0]])
        # Resize everything to expected total size
        for dset_name in to_append:
            dset = matrix[dset_name]
            dset.resize((total_nnz,))
        # Now append new data to this
        for small_file, s_start, s_end in fn_start_ends[1:]:
            new_mat = CountMatrix.load_h5_file(small_file)
            start = end
            end = len(new_mat.m.data) + end
            matrix[h5_constants.H5_MATRIX_DATA_ATTR][start:end] = new_mat.m.data
            matrix[h5_constants.H5_MATRIX_INDICES_ATTR][start:end] = new_mat.m.indices
            # Update the column index as well, by adding an offset to account for the data we're appending.
            cur_indptr = new_mat.m.indptr.astype(np.uint64)
            cur_indptr += start
            num_entries_before = (  # suspect is zero most of the time as we grouped by columns
                ind_ptr[s_end + 1] - ind_ptr[s_end]
            )
            assert num_entries_before == 0, "Data was not cleanly separated into distinct columns"
            ind_ptr[s_start:] = cur_indptr[s_start:]  # +num_entries_before
            del new_mat
        matrix[h5_constants.H5_MATRIX_INDPTR_ATTR][:] = ind_ptr


def make_matrix_attrs_count(sample_id, gem_groups, chemistry):
    matrix_attrs = make_library_map_count(sample_id, gem_groups)
    matrix_attrs[h5_constants.H5_CHEMISTRY_DESC_KEY] = chemistry
    return matrix_attrs


def load_matrix_h5_metadata(filename):
    """Get matrix metadata attributes from an HDF5 file."""
    # TODO: Consider moving these to the MATRIX key instead of the root group
    filename = ensure_binary(filename)
    h5_version = CountMatrix.get_format_version_from_h5(filename)
    if h5_version == 1:
        return _load_matrix_legacy_v1_h5_metadata(filename.decode())
    else:
        attrs = {}
        with h5.File(filename, "r") as f:
            for key in h5_constants.H5_METADATA_ATTRS:
                val = f.attrs.get(key)
                if val is not None:
                    if np.isscalar(val) and hasattr(val, "item"):
                        # Coerce numpy scalars to python types.
                        # In particular, force np.unicode_ to unicode (python2)
                        #    or str (python 3)
                        # pylint: disable=no-member
                        attrs[key] = val.item()
                    else:
                        attrs[key] = val
            return attrs


# (needed for compatibility with v1 matrices)
def _load_matrix_legacy_v1_h5_metadata(filename: str):
    attrs = {}
    with tables.open_file(filename, "r") as f:
        all_attrs = f.get_node("/")._v_attrs
        for key in h5_constants.H5_METADATA_ATTRS:
            if hasattr(all_attrs, key):
                val = getattr(all_attrs, key)
                if np.isscalar(val) and hasattr(val, "item"):
                    # Coerce numpy scalars to python types.
                    # In particular, force np.unicode_ to unicode (python2)
                    #    or str (python 3)
                    attrs[key] = val.item()
                else:
                    attrs[key] = val
    return attrs


def load_matrix_h5_custom_attrs(filename):
    """Get matrix metadata attributes from an HDF5 file."""
    filename = ensure_binary(filename)
    h5_version = CountMatrix.get_format_version_from_h5(filename)
    if h5_version == 1:
        # no support for custom attrs in older versions
        return {}

    attrs = {}
    with h5.File(filename, "r") as f:
        for key, val in f.attrs.items():
            if key not in MATRIX_H5_BUILTIN_ATTRS:
                attrs[key] = val
        return attrs


def make_library_map_count(sample_id, gem_groups):
    # store gem group mapping for use by Cell Loupe
    unique_gem_groups = sorted(set(gem_groups))
    library_map = {
        h5_constants.H5_LIBRARY_ID_MAPPING_KEY: np.array(
            [sample_id] * len(unique_gem_groups), dtype="S"
        ),
        h5_constants.H5_ORIG_GEM_GROUP_MAPPING_KEY: np.array(unique_gem_groups, dtype=int),
    }
    return library_map


def make_library_map_aggr(gem_group_index) -> dict[str, np.ndarray]:
    # store gem group mapping for use by Cell Loupe
    library_ids = []
    original_gem_groups = []
    # Sort numerically by new gem group
    for _, (lid, og) in sorted(gem_group_index.items(), key=lambda pair: int(pair[0])):
        library_ids.append(lid)
        original_gem_groups.append(og)
    library_map = {
        h5_constants.H5_LIBRARY_ID_MAPPING_KEY: np.array(library_ids, dtype="S"),
        h5_constants.H5_ORIG_GEM_GROUP_MAPPING_KEY: np.array(original_gem_groups, dtype=int),
    }
    return library_map


def get_gem_group_index(matrix_h5: str):
    with tables.open_file(matrix_h5, mode="r") as f:
        try:
            library_ids = f.get_node_attr("/", h5_constants.H5_LIBRARY_ID_MAPPING_KEY)
            original_gem_groups = f.get_node_attr("/", h5_constants.H5_ORIG_GEM_GROUP_MAPPING_KEY)
        except AttributeError:
            return None
    library_map = {}
    for ng, (lid, og) in enumerate(zip(library_ids, original_gem_groups), start=1):
        library_map[ng] = (lid, og)
    return library_map


def inplace_csc_column_normalize_l2(X):
    """Perform in-place column L2-normalization of input matrix X.

    >>> import numpy as np
    >>> import scipy.sparse as sp
    >>> from sklearn.preprocessing import normalize
    >>> a = np.arange(12, dtype='float').reshape((3, 4))
    >>> b = sp.csc_matrix(a)
    >>> inplace_csc_column_normalize_l2(b)
    >>> np.all(normalize(a, axis=0) == b)
    True
    """
    assert X.getnnz() == 0 or isinstance(X.data[0], np.float32 | float)
    assert isinstance(X, sp_sparse.csc_matrix)
    for i in range(X.shape[1]):
        s = 0.0
        for j in range(X.indptr[i], X.indptr[i + 1]):
            s += X.data[j] * X.data[j]
        if s == 0.0:
            continue
        s = np.sqrt(s)
        for j in range(X.indptr[i], X.indptr[i + 1]):
            X.data[j] /= s
