#!/usr/bin/env python
#
# Copyright (c) 2019 10X Genomics, Inc. All rights reserved
#
"""Utils for summarizing and computing statistics on RNA reads and UMIs from molecule_info."""

from __future__ import annotations

from collections.abc import Iterable
from copy import deepcopy
from typing import Any

import h5py as h5
import numpy as np
import pandas as pd
from six import ensure_binary, ensure_str

import cellranger.utils as cr_utils
from cellranger import molecule_counter_extensions as cr_mce
from cellranger.feature_ref import GENOME_FEATURE_TAG, FeatureReference
from cellranger.molecule_counter import (
    BARCODE_IDX_COL_NAME,
    COUNT_COL_NAME,
    FEATURE_IDX_COL_NAME,
    GEM_GROUP_COL_NAME,
    LIBRARY_IDX_COL_NAME,
    MOLECULE_INFO_COLUMNS,
    UMI_COL_NAME,
    UMI_TYPE_COL_NAME,
    BarcodeInfo,
    MoleculeCounter,
)

# column header string constants useful for building df
FEATURE_REF_COLS = ["feature_type", "id", "name", "index"]
MOL_INFO_CELL_COL = "is_cell"
IS_CELL_FORMAT_STRING = "{}_cells"
FEATURE_DF_COUNT_COL = "num_reads"
FEATURE_DF_UMI_COL = "num_umis"
FEATURE_DF_DUP_COL = "dup_frac"
FEATURE_DF_BARCODE_SEQ_COL = "barcode"
FEATURE_DF_BARCODE_COUNT_COL = "num_barcodes"
FEATURE_DF_FEATURE_COUNT_COL = "num_features"
FEATURE_ID_COL = "feature_id"
FEATURE_NAME_COL = "feature_name"

# df columns with barcode information useful for summarizing by barcode
BARCODE_METADATA = [
    BARCODE_IDX_COL_NAME,
    FEATURE_DF_BARCODE_SEQ_COL,
    MOL_INFO_CELL_COL,
    GEM_GROUP_COL_NAME,
]

# mol_info chunk size used in collapse_feature_counts and collapse_barcode_chunks
CHUNK_SIZE = 20_000_000


# pylint: disable=invalid-name
def _downcast_large_uint_vector(x):
    """A roll our own version of `pd.to_numeric(...,downcast="unsigned")`.

    Profiling showed
    that the pandas version was incredibly memory and CPU intensive due to some inefficient code
    leading to an expensive call to `np.allclose(new_result, result, rtol=0)`
    """
    levels = [np.uint8, np.uint16, np.uint32, np.uint64]
    if x.dtype not in levels:
        type_name = str(type(x))
        if isinstance(x, np.ndarray):
            type_name = type_name + str(x.dtype)
        print(f"Downcasting function expected a uint64 or smaller, type was: {type_name}")
        return pd.to_numeric(x, downcast="unsigned")
    to_check = levels[: levels.index(x.dtype)]
    for smaller_type in to_check:
        ii = np.iinfo(smaller_type)
        max_val = ii.max
        # prevent overflow in the very common log(1 + x)
        below = np.all(x < max_val)
        if below:
            return x.astype(smaller_type)
    return x


def to_numeric_safe(
    x: float | int | tuple | np.ndarray | pd.Series, downcast: str | None = "unsigned"
):
    """Wrapper around pandas.to_numeric to avoid bugs and slow performance.

    Args:
        x: value to be converted
        downcast: one of (None, 'unsigned', 'float')

    Returns:
        (depends on input type): input value as numeric type, possibly downcast
    """
    if downcast == "unsigned":
        if len(x) == 0:
            # pd.to_numeric gives ValueError
            return x.astype("uint8")
        else:
            # pd.to_numeric is very slow here
            return _downcast_large_uint_vector(x)
    elif downcast in (None, "float"):
        # issues with 'to_numeric' do not apply
        return pd.to_numeric(x, downcast=downcast)
    else:
        # _downcast_large_int_vector could be implemented for
        # downcast options 'signed' and 'integer'
        raise NotImplementedError(f"to_numeric wrapper not implemented for downcast='{downcast}'")


def _conv_bytes(v):
    if isinstance(v, bytes):
        return v.decode()
    else:
        return v


def sanitize_dataframe(df: pd.DataFrame, inplace=False):
    """Convert any bytes in the DataFrame to str."""
    if not inplace:
        df = df.copy(deep=False)

    if df.index.dtype == bytes:
        df.rename(ensure_str, axis="index", inplace=True)
    elif df.index.dtype == object:
        df.rename(_conv_bytes, axis="index", inplace=True)

    df.rename(ensure_str, axis="columns", inplace=True)

    for col in df.columns:
        if df[col].dtype == object:
            df[col] = df[col].apply(_conv_bytes)
        elif df[col].dtype == bytes:
            df[col] = df[col].apply(ensure_str)
        elif isinstance(df[col].dtype, pd.CategoricalDtype):
            if df[col].cat.categories.dtype == object:
                df[col] = df[col].cat.rename_categories(ensure_str)
            elif df[col].cat.categories.dtype == bytes:
                df[col] = df[col].cat.rename_categories(_conv_bytes)

    return df


def _get_idx_mol(
    mc: MoleculeCounter,
    filter_library_idx: int | list[int] | None = None,
    filter_feature_idx: int | list[int] | None = None,
    filter_barcode_idx: int | list[int] | None = None,
) -> np.ndarray[int, np.dtype[np.bool_]]:
    if filter_library_idx is None:
        idx_mol = np.ones(mc.nrows(), dtype=bool)
    else:
        idx_mol = np.isin(mc.get_column(LIBRARY_IDX_COL_NAME), filter_library_idx)

    if filter_barcode_idx is not None:
        idx_mol &= np.isin(mc.get_column_lazy(BARCODE_IDX_COL_NAME), filter_barcode_idx)

    # filter features if needed
    if filter_feature_idx is not None:
        idx_mol &= cr_mce.get_indices_for_values(
            mc, [FEATURE_IDX_COL_NAME], [(x,) for x in filter_feature_idx]
        )
    return idx_mol


def _get_is_cell(
    mc: MoleculeCounter,
    idx_mol: np.ndarray[int, np.dtype[np.bool_]],
    exclude_cells: bool,
    exclude_noncells: bool,
    with_cell_call: bool,
) -> np.ndarray[int, np.dtype[np.bool_]] | None:
    """Get vector of cell calls if needed."""
    if exclude_cells and exclude_noncells:
        raise ValueError("Can't exclude both cells and non-cells")
    elif exclude_cells or exclude_noncells or with_cell_call:
        pass_filter = mc.get_barcode_info().pass_filter
        is_cell = np.zeros(np.sum(idx_mol), dtype=bool)
        for library_idx in set(pass_filter[:, BarcodeInfo.PASS_FILTER_LIBRARY_IDX]):
            # which are the cell-associated barcodes from this library?
            which_barcodes = pass_filter[
                pass_filter[:, BarcodeInfo.PASS_FILTER_LIBRARY_IDX] == library_idx, 0
            ]
            is_bc_and_lib = cr_mce.get_indices_for_values(
                mc, [LIBRARY_IDX_COL_NAME], [(library_idx,)]
            )[idx_mol]
            is_bc_and_lib &= cr_mce.get_indices_for_values(
                mc, [BARCODE_IDX_COL_NAME], which_barcodes
            )[idx_mol]
            is_cell |= is_bc_and_lib
            del is_bc_and_lib
        if exclude_cells:
            idx_mol[idx_mol == 1] &= ~is_cell
            return is_cell[~is_cell]
        elif exclude_noncells:
            idx_mol[idx_mol == 1] &= is_cell
            return is_cell[is_cell]
        return is_cell
    return None


def _mol_info_df_from_h5(  # pylint: disable=too-many-arguments
    mc: MoleculeCounter,
    exclude_noncells: bool,
    exclude_cells: bool,
    with_umi: bool,
    with_gem_group: bool,
    with_library_idx: bool,
    with_cell_call: bool,
    with_umi_type: bool,
    filter_library_idx: int | list[int] | None,
    filter_feature_idx: int | list[int] | None,
    filter_barcode_idx: int | list[int] | None,
    downsample: float | None,
):
    # select molecules from the specified sample(s)
    idx_mol = _get_idx_mol(mc, filter_library_idx, filter_feature_idx, filter_barcode_idx)

    is_cell = _get_is_cell(mc, idx_mol, exclude_cells, exclude_noncells, with_cell_call)

    dict_const = {}  # At about 2 GB now

    # if need be, downsample selected molecules at a certain rate
    if downsample is not None:
        assert 0 <= downsample < 1
        new_read_counts = mc.get_column(COUNT_COL_NAME)[idx_mol == 1]
        new_read_counts = np.random.binomial(new_read_counts, downsample)
        if is_cell is not None and with_cell_call:
            is_cell = is_cell[new_read_counts > 0]
        idx_mol[idx_mol == 1] &= new_read_counts > 0
        new_read_counts = (new_read_counts[new_read_counts > 0]).astype(np.uint32)

    if with_cell_call:
        dict_const[MOL_INFO_CELL_COL] = is_cell
    del is_cell

    mol_info_cols = deepcopy(MOLECULE_INFO_COLUMNS)
    if not with_umi:
        del mol_info_cols[UMI_COL_NAME]
    if not with_library_idx:
        del mol_info_cols[LIBRARY_IDX_COL_NAME]
    if not with_gem_group:
        del mol_info_cols[GEM_GROUP_COL_NAME]
    if not with_umi_type:
        del mol_info_cols[UMI_TYPE_COL_NAME]

    for x in mol_info_cols:
        # skip if downsample was used
        if x == COUNT_COL_NAME and downsample is not None:
            val = new_read_counts
            del new_read_counts
        else:
            val = mc.get_column(x)[idx_mol]
        val = to_numeric_safe(val)
        if x == COUNT_COL_NAME:
            x = FEATURE_DF_COUNT_COL
        dict_const[x] = val
    del idx_mol
    return pd.DataFrame(dict_const)  # very expensive in terms of memory, assume it doubles it


def mol_info_from_h5(  # pylint: disable=too-many-arguments
    mol_info,
    exclude_noncells: bool = False,
    exclude_cells: bool = False,
    with_umi: bool = True,
    with_barcode_seq: bool = False,
    with_gem_group: bool = True,
    with_library_idx: bool = True,
    with_cell_call: bool = True,
    with_feature_info: bool = False,
    with_umi_type: bool = False,
    filter_library_idx: int | list[int] | None = None,
    filter_feature_idx: int | list[int] | None = None,
    filter_barcode_idx: int | list[int] | None = None,
    downsample: float | None = None,
    genome_col: bool = False,
) -> pd.DataFrame:
    """Load molecule_info.h5 file as a pandas DataFrame.

    Returns a data frame containing the per-molecule columns from the molecule_info
    file (barcode_idx, count, feature_idx, gem_group, library_idx, umi) along with
    a column 'is_cell'. Numeric columns are automatically downcast to the
    smallest appropriate type.

    Args:
        mol_info (str or object): path to a cellranger molecule_info.h5 file or instance of
            MoleculeCounter object
        exclude_noncells (bool): exclude molecules from non-cell-associated barcodes
        exclude_cells (bool): exclude molecules from cell-associated barcodes
        with_umi (bool): include a column 'umi' containing the 2-bit encoded UMI
        with_barcode_seq (bool): include a column 'barcode' containing the full
            barcode sequences
        with_gem_group (bool): include a column gem_group
        with_library_idx (bool): include `library_idx` column
        with_cell_call (bool): include a column 'is_cell' indicating whether
            the barcode associated with each molecule was called as a cell
        with_feature_info (bool): merge the entire result with the feature reference
            (not memory-efficient, but useful for small analysis tasks)
        with_umi_type (bool): include a column `umi_type` indicating whether the UMI is
            transcriptomic or not [default: False]
        filter_library_idx: library indices (i.e., into
            the library_info dataset) to include, or None to include all of them
        filter_barcode_idx: barcode indices (i.e., into
            the barcodes dataset) to include, or None to include all of them. Note
            that exclude_noncells and exclude_cells will still be applied if True
        filter_feature_idx: feature indices (i.e., into
            the features dataset) to include, or None to include all of them
        downsample: between 0 and 1 exclusive, or None.  The rate at which to
            downsample reads in the molecule info.
        genome_col: If with_feature_info is set, will add a genome column to the dataframe.

    Returns:
        pd.DataFrame: molecule_info data (see above)
    """
    # pylint: disable=invalid-unary-operand-type
    if isinstance(mol_info, MoleculeCounter):
        obj_passed = True
        mc = mol_info
    else:
        obj_passed = False
        mc = MoleculeCounter.open(mol_info, "r")

    df = _mol_info_df_from_h5(
        mc,
        exclude_noncells,
        exclude_cells,
        with_umi,
        with_gem_group,
        with_library_idx,
        with_cell_call,
        with_umi_type,
        filter_library_idx,
        filter_feature_idx,
        filter_barcode_idx,
        downsample,
    )

    if with_barcode_seq:
        df[FEATURE_DF_BARCODE_SEQ_COL] = mc.get_barcodes()[df[BARCODE_IDX_COL_NAME]]
        df[FEATURE_DF_BARCODE_SEQ_COL] = df[FEATURE_DF_BARCODE_SEQ_COL].astype("category")

    if with_feature_info:
        feature_ref = feature_ref_from_h5(mol_info, genome_col)
        feature_ref.rename(columns={"index": "feature_idx"}, inplace=True)
        df = df.merge(feature_ref, how="left", on="feature_idx")
    if not obj_passed:
        mc.close()
    return df


def feature_ref_from_h5(fname_or_mol_info, genome_col: bool = False, filter_feature_types=None):
    """Load feature reference from h5 file as a pandas DataFrame.

    Args:
        fname_or_mol_info (str or object): path to a cellranger molecule_info.h5 file or instance of
            MoleculeCounter object, or a path to a matrix h5 file
        genome_col (bool): Do we want a column with the "Genome" in it? This is to preserve
        older behaviour used by puppy.panda_utils
        filter_feature_types (list or None): list of feature types to select from feature reference,
            all are kept if None

    Returns:
        pd.DataFrame: the feature reference
    """
    if isinstance(fname_or_mol_info, MoleculeCounter):
        feature_ref = fname_or_mol_info.get_feature_ref()
    else:
        f = h5.File(ensure_binary(fname_or_mol_info), "r")
        if "matrix" in f.keys():
            group = f["matrix"]["features"]
        else:
            group = f["features"]
        feature_ref = FeatureReference.from_hdf5(group)
        f.close()
    feature_ref = feature_ref.id_map
    columns = deepcopy(FEATURE_REF_COLS)
    vals = [[getattr(v, col) for col in FEATURE_REF_COLS] for v in feature_ref.values()]
    if genome_col:
        # Used by PD only code at the moment.
        for i, val in enumerate(feature_ref.values()):
            vals[i].append(val.tags.get(GENOME_FEATURE_TAG, None))
        columns.append(GENOME_FEATURE_TAG)

    df = pd.DataFrame(
        vals,
        columns=columns,
    )

    del vals
    del feature_ref

    if filter_feature_types is not None:
        df = df.loc[df["feature_type"].isin(filter_feature_types)]

    # reorder df by feature indices
    df.sort_values("index", inplace=True)
    df.reset_index(drop=True, inplace=True)

    # FIXME in pandas 0.24: df = df.astype('category')
    for col in df.columns:
        # don't convert 'id' to category -- it stuff "b'" into the string
        if col not in ["index", "id"]:
            df[col] = df[col].astype("category")

    df.rename(columns={"id": FEATURE_ID_COL, "name": FEATURE_NAME_COL}, inplace=True)
    return df


def summarize_by_feature(feature_ref_df: pd.DataFrame, mol_info_df: pd.DataFrame):
    """Compute a per-feature summary from a molecule_info DataFrame.

    For each feature, these metrics are computed:
        num_umis - # UMIs associated with the feature
        num_reads - # reads associated with the feature
        num_barcodes - # distinct barcodes containing at least one read/molecule
            mapped to the feature
        dup_frac - proportion of reads corresponding to already observed molecules,
            i.e., sequencing saturation. This is (1 - num_umis/num_reads),
            and 0 in case there are no reads.

    Args:
        feature_ref_df (pandas.DataFrame)
        mol_info_df (pandas.DataFrame)

    Returns:
        pandas.DataFrame: per-feature summary statistics from mol_info_df
                            with the following columns
                            * feature_type category
                            * feature_id   object
                            * feature_name category
                            * index        int64
                            * num_umis     int64
                            * num_reads    int64
                            * num_barcodes int64
                            * dup_frac     float64
                            * num_umis_cells int64
                            * num_reads_cells int64
                            * num_barcodes_cells int64
                            * dup_frac_cells float64
    """
    raw = _feature_summary(mol_info_df)
    filtered = _feature_summary(mol_info_df.query(MOL_INFO_CELL_COL))
    filtered.columns = [IS_CELL_FORMAT_STRING.format(x) for x in filtered.columns]

    combined = pd.concat([raw, filtered], axis=1)
    combined_cols = combined.columns

    result = feature_ref_df.merge(
        combined, how="left", left_on="index", right_on=FEATURE_IDX_COL_NAME
    )
    result[combined_cols] = result[combined_cols].fillna(0, downcast="infer")

    return result


def summarize_by_barcode(
    mol_info_df: pd.DataFrame, feature_idx_subsets: dict[str, list[int]] | None = None
):
    """Compute a per-barcode summary from a molecule_info DataFrame.

    For each feature, these metrics are computed::

        num_umis - # UMIs associated with the feature
        num_reads - # reads associated with the feature
        num_features - # distinct features containing at least one read/molecule
            assigned to this barcode
        dup_frac - proportion of reads corresponding to already observed molecules,
            i.e., sequencing saturation. This is (1 - num_umis/num_reads),
            and 0 in case there are no reads.

    All features are currently treated the same (regardless of feature_type or genome).
    The subset_feature_idx argument can be used to produce separate metrics on a
    defined subset of features.
    Any metadata columns will be propagated (barcode_idx, barcode sequence,
    and is_cell, as defined by BARCODE_METADATA)

    Args:
        mol_info_df (pandas.DataFrame): TODO
        feature_idx_subsets (Dict[str, List[int]]): produce
            separate metrics for certain named subsets of features, specified
            by their indices within the feature reference (feature_idx column
            in MoleculeCounter). There will be an additional output column
            called 'feature_subset'; the rows marked 'all_features' contain the
            metrics over all features.

    Returns:
        pandas.DataFrame: per-feature summary statistics from mol_info_df
    """
    result = _barcode_summary(mol_info_df)
    if feature_idx_subsets:
        result["feature_subset"] = "all_features"

        metadata_cols = [x for x in result.columns if x in BARCODE_METADATA]
        metadata = result[metadata_cols]

        for subset_name, subset in feature_idx_subsets.items():
            idx_subset = np.isin(mol_info_df.feature_idx, subset)
            summary_subset = _barcode_summary(mol_info_df[idx_subset])
            summary_subset.drop(metadata_cols, axis=1, inplace=True)
            summary_subset = metadata.merge(
                summary_subset, how="left", left_index=True, right_index=True
            )
            summary_subset.fillna(0, downcast="infer", inplace=True)
            summary_subset["feature_subset"] = subset_name
            result = pd.concat([result, summary_subset])
        # cast feature subset to category
        result.feature_subset = result.feature_subset.astype("category")
        # reorder columns
        main_cols = [x for x in result.columns if x != "feature_subset"]
        result = result[["feature_subset"] + main_cols]
    return result


def calculate_duplicate_fraction(df: pd.DataFrame, umi_col_name, read_col_name):
    """Calculate the fraction of duplication/sequencing saturation.

    Args:
        df: pandas dataframe
        umi_col_name: column with umi counts
        read_col_name: column with read counts

    Returns:
        A vector of duplicate fractions
    """
    return 1 - (df[umi_col_name].clip(lower=1) / df[read_col_name].clip(lower=1))


def _feature_summary(df: pd.DataFrame):
    """Aggregation for summarize_by_feature."""
    grouped = df.groupby(FEATURE_IDX_COL_NAME)

    result = pd.DataFrame(
        {
            FEATURE_DF_UMI_COL: to_numeric_safe(grouped[FEATURE_DF_COUNT_COL].count()),
            FEATURE_DF_COUNT_COL: to_numeric_safe(grouped[FEATURE_DF_COUNT_COL].sum()),
            FEATURE_DF_BARCODE_COUNT_COL: to_numeric_safe(grouped[BARCODE_IDX_COL_NAME].nunique()),
        }
    )

    # dup frac defined to be 0 when there are 0 reads
    result[FEATURE_DF_DUP_COL] = calculate_duplicate_fraction(
        result, FEATURE_DF_UMI_COL, FEATURE_DF_COUNT_COL
    )
    return result


def _barcode_summary(df: pd.DataFrame):
    """Aggregation for summarize_by_barcode."""
    grouped = df.groupby(BARCODE_IDX_COL_NAME)

    # keep barcode_idx even if it's the df index, useful for selection without
    # special handling
    metadata_present = [x for x in df.columns if x in BARCODE_METADATA]
    metadata = pd.DataFrame({x: grouped[x].first() for x in metadata_present})

    result = pd.DataFrame(
        {
            FEATURE_DF_UMI_COL: to_numeric_safe(grouped[FEATURE_DF_COUNT_COL].count()),
            FEATURE_DF_COUNT_COL: to_numeric_safe(grouped[FEATURE_DF_COUNT_COL].sum()),
            FEATURE_DF_FEATURE_COUNT_COL: to_numeric_safe(grouped[FEATURE_IDX_COL_NAME].nunique()),
        }
    )
    # dup frac defined to be 0 when there are 0 reads
    result[FEATURE_DF_DUP_COL] = 1 - (
        result[FEATURE_DF_UMI_COL].clip(lower=1) / result[FEATURE_DF_COUNT_COL].clip(lower=1)
    )

    return pd.concat([metadata, result], axis=1)


def collapse_feature_counts(
    mc: MoleculeCounter,
    filter_library_idx: list[int] | None = None,
    barcode_genome: str | None = None,
    tgt_chunk_len: int = CHUNK_SIZE,
    downsample: float | None = None,
):
    """Get read and UMI counts per feature by aggregating counts over mol_info.

    in chunks.

    Args:
        mc (MoleculeCounter): instance of MoleculeCounter object
        filter_library_idx (list of ints or None): specifies whether to restrict counts to only certain libraries
        barcode_genome (str | None): Use only barcodes from this genome. If None, use all barcodes.
        tgt_chunk_len (int): number of rows by which to chunk mol_info
        downsample (float): downsample fraction

    Returns:
        pandas.DataFrame containing per-feature summary statistics (UMIs, reads, dup_rate, and feature metadata)
    """
    barcode_info = mc.get_barcode_info()
    barcode_passing_filter = barcode_info.pass_filter
    if barcode_genome is not None:
        genome_idx = barcode_info.genomes.index(barcode_genome)
        barcode_passing_filter = barcode_passing_filter[
            barcode_passing_filter[:, BarcodeInfo.PASS_FILTER_GENOME_IDX] == genome_idx
        ]
    if filter_library_idx is None:
        filter_library_idx = [int(lib["library_id"]) for lib in mc.get_library_info()]

    num_features = mc.get_feature_ref().get_num_features()
    read_counts_by_index = np.zeros(num_features, dtype=np.uint64)
    umi_counts_by_index = np.zeros(num_features, dtype=np.uint64)
    read_counts_in_cells_by_index = np.zeros(num_features, dtype=np.uint64)
    umi_counts_in_cells_by_index = np.zeros(num_features, dtype=np.uint64)

    def collapse_chunk(chunk_start: int, chunk_len: int, filter_library_idx: list[int]):
        chunk_stop = chunk_start + chunk_len
        idx_mol = np.ones(chunk_len, dtype="bool")

        if filter_library_idx:
            idx_mol &= np.isin(
                mc.get_column_lazy(LIBRARY_IDX_COL_NAME)[chunk_start:chunk_stop], filter_library_idx
            )

        is_cell = np.zeros(idx_mol.shape, dtype="bool")
        for library_idx in filter_library_idx:
            cell_barcode_indices = barcode_passing_filter[
                barcode_passing_filter[:, BarcodeInfo.PASS_FILTER_LIBRARY_IDX] == library_idx, 0
            ]
            is_cell |= np.isin(
                mc.get_column_lazy(BARCODE_IDX_COL_NAME)[chunk_start:chunk_stop],
                cell_barcode_indices,
            )
        is_cell = is_cell[idx_mol]

        # increment feature counts
        feature_indices = mc.get_column_lazy(FEATURE_IDX_COL_NAME)[chunk_start:chunk_stop][idx_mol]
        read_counts = mc.get_column_lazy(COUNT_COL_NAME)[chunk_start:chunk_stop][idx_mol]
        del idx_mol

        if downsample is not None and downsample < 1.0:
            read_counts = np.random.binomial(read_counts, downsample)
            feature_indices = feature_indices[read_counts > 0]
            is_cell = is_cell[read_counts > 0]
            read_counts = read_counts[read_counts > 0]

        np.add.at(umi_counts_by_index, feature_indices, np.ones(feature_indices.shape[0]))
        np.add.at(read_counts_by_index, feature_indices, read_counts)
        # increment feature counts in cells
        feature_indices = feature_indices[is_cell]
        read_counts = read_counts[is_cell]
        np.add.at(umi_counts_in_cells_by_index, feature_indices, np.ones(feature_indices.shape[0]))
        np.add.at(read_counts_in_cells_by_index, feature_indices, read_counts)

    for chunk_start, chunk_len in mc.get_chunks(tgt_chunk_len, preserve_boundaries=True):
        collapse_chunk(chunk_start, chunk_len, filter_library_idx)

    result = pd.DataFrame(
        {
            FEATURE_DF_UMI_COL: to_numeric_safe(umi_counts_by_index),
            FEATURE_DF_COUNT_COL: to_numeric_safe(read_counts_by_index),
            f"{FEATURE_DF_UMI_COL}_cells": to_numeric_safe(umi_counts_in_cells_by_index),
            f"{FEATURE_DF_COUNT_COL}_cells": to_numeric_safe(read_counts_in_cells_by_index),
            FEATURE_IDX_COL_NAME: to_numeric_safe(np.arange(umi_counts_by_index.shape[0])),
        }
    )
    for suffix in ["", "_cells"]:
        result[FEATURE_DF_DUP_COL + suffix] = 1 - (
            result[FEATURE_DF_UMI_COL + suffix].clip(lower=1)
            / result[FEATURE_DF_COUNT_COL + suffix].clip(lower=1)
        )

    feature_ref_df = feature_ref_from_h5(mc, genome_col=True)
    feature_ref_df.rename(columns={"index": FEATURE_IDX_COL_NAME}, inplace=True)
    feature_ref_df = feature_ref_df.merge(result, how="left", on=FEATURE_IDX_COL_NAME)

    return feature_ref_df


def _barcode_is_cell(
    is_cell,
    library_info: Iterable[dict[str, Any]],
    bcs_passing_filter: np.ndarray,
    filter_library_idx: list[int] | None,
    num_barcodes,
):
    for lib in library_info:
        lib_idx = int(lib["library_id"])
        if filter_library_idx is not None and lib["library_id"] not in filter_library_idx:
            continue
        gg = lib[GEM_GROUP_COL_NAME]
        cell_indices = bcs_passing_filter[bcs_passing_filter[:, 1] == lib_idx, 0]
        # offset barcode indices according to gem group index
        cell_indices = np.add(cell_indices, (gg - 1) * num_barcodes)
        is_cell[cell_indices] = True
    return is_cell


def collapse_barcode_counts(
    mc: MoleculeCounter,
    filter_library_idx: list[int] | None = None,
    filter_feature_idx=None,
    tgt_chunk_len: int = CHUNK_SIZE,
    downsample: float | None = None,
):
    """Get read and UMI counts per feature by aggregating counts over mol_info in chunks.

    Args:
        mc (MoleculeCounter): instance of MoleculeCounter object
        filter_library_idx (list of ints or None): specifies whether to restrict counts to only certain libraries
        filter_feature_idx: TODO: document
        tgt_chunk_len (int): number of rows by which to chunk mol_info
        downsample (float): downsample fraction

    Returns:
        pandas.DataFrame containing per-barcode summary statistics (UMIs, reads, dup_rate, and barcode metadata)
    """
    gem_wells = sorted(list(set(mc.get_gem_groups())))
    num_barcodes = len(mc.get_barcodes())
    ordered_barcodes = cr_utils.format_barcode_seqs(mc.get_barcodes(), gem_wells)
    read_counts_by_index = np.zeros(len(ordered_barcodes), dtype=np.uint64)
    umi_counts_by_index = np.zeros(len(ordered_barcodes), dtype=np.uint64)

    def collapse_chunk(
        chunk_start: int,
        chunk_len: int,
        umi_counts_by_index: np.ndarray[int, np.dtype[np.uint64]],
        read_counts_by_index: np.ndarray[int, np.dtype[np.uint64]],
    ):
        chunk_stop = chunk_start + chunk_len
        idx_mol = np.ones(chunk_len, dtype=bool)

        if filter_library_idx:
            idx_mol &= np.isin(
                mc.get_column_lazy(LIBRARY_IDX_COL_NAME)[chunk_start:chunk_stop], filter_library_idx
            )

        if filter_feature_idx:
            idx_mol &= np.isin(
                mc.get_column_lazy(FEATURE_IDX_COL_NAME)[chunk_start:chunk_stop], filter_feature_idx
            )

        # increment feature counts
        barcode_indices = mc.get_column_lazy(BARCODE_IDX_COL_NAME)[chunk_start:chunk_stop][idx_mol]
        gg_indices = mc.get_column_lazy(GEM_GROUP_COL_NAME)[chunk_start:chunk_stop][idx_mol]
        read_counts = mc.get_column_lazy(COUNT_COL_NAME)[chunk_start:chunk_stop][idx_mol]
        if downsample is not None and downsample < 1.0:
            read_counts = np.random.binomial(read_counts, downsample)
            barcode_indices = barcode_indices[read_counts > 0]
            gg_indices = gg_indices[read_counts > 0]
            read_counts = read_counts[read_counts > 0]

        del idx_mol

        # barcodes must be offset by num_bcs * gem_well_index in ordered_barcodes
        bc_idx_offset = np.multiply(np.subtract(gg_indices, 1), num_barcodes)
        np.add.at(
            umi_counts_by_index,
            np.add(barcode_indices, bc_idx_offset),
            np.ones(barcode_indices.shape[0]),
        )
        np.add.at(
            read_counts_by_index,
            np.add(barcode_indices, bc_idx_offset),
            read_counts,
        )

    for chunk_start, chunk_len in mc.get_chunks(tgt_chunk_len, preserve_boundaries=True):
        collapse_chunk(chunk_start, chunk_len, umi_counts_by_index, read_counts_by_index)

    result = pd.DataFrame(
        {
            FEATURE_DF_UMI_COL: to_numeric_safe(umi_counts_by_index),
            FEATURE_DF_COUNT_COL: to_numeric_safe(read_counts_by_index),
        }
    )
    del read_counts_by_index, umi_counts_by_index

    result[FEATURE_DF_BARCODE_SEQ_COL] = ordered_barcodes
    del ordered_barcodes

    # add col for whether or not it's a cell
    result[MOL_INFO_CELL_COL] = _barcode_is_cell(
        np.full(shape=result.shape[0], fill_value=False, dtype=bool),
        mc.get_library_info(),
        mc.get_barcode_info().pass_filter,
        filter_library_idx,
        num_barcodes,
    )

    # to save space, keep only rows that are cells or have non-zero counts
    result = result[(result[FEATURE_DF_UMI_COL] > 0) | (result[MOL_INFO_CELL_COL])]

    # dup frac defined to be 0 when there are 0 reads
    result[FEATURE_DF_DUP_COL] = 1 - (
        result[FEATURE_DF_UMI_COL].clip(lower=1) / result[FEATURE_DF_COUNT_COL].clip(lower=1)
    )

    return result
