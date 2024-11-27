#!/usr/bin/env python
#
# Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
#
#

"""Feature Assignments object."""
# pylint: disable=invalid-name, inconsistent-return-statements
from __future__ import annotations

import itertools
import json
from collections.abc import MutableMapping
from doctest import testmod
from functools import cmp_to_key
from typing import TYPE_CHECKING, NamedTuple

import numpy as np
import pandas as pd
from six import ensure_binary, ensure_str

import tenkit.stats as tk_stats
from cellranger import utils as cr_utils
from tenkit.safe_json import safe_jsonify

if TYPE_CHECKING:
    from cellranger.matrix import CountMatrix

NUM_FEATURES = "num_features"
FEATURE_CALL = "feature_call"
NUM_UMIS = "num_umis"
CELL_BARCODE = "cell_barcode"


class CellsPerFeature(MutableMapping):
    """A class to hold the cells per feature.

    Basically a dictionary like:

            "Tag1": ["BC1", "BC2", "BC3"],
            "Tag2": ["BC3", "BC4"],
            "Tag3": ["BC1", "BC3"],
    """

    def __init__(self):
        """Initialize like a dictionary."""
        self._data: dict[bytes, list[bytes]] = {}

    def __getitem__(self, item: bytes):
        """Returns a list for a given assignment."""
        assert isinstance(item, bytes)
        return self._data[item]

    def __setitem__(self, key: bytes, value: list[bytes]):
        assert isinstance(value, list)
        self._data[key] = value

    def __delitem__(self, key: bytes):
        assert isinstance(key, bytes)
        del self._data[key]

    def __iter__(self):
        return iter(self._data)

    def __len__(self):
        return len(self._data)

    def save_to_file(self, filename):
        """Encodes the file.

        Args:
            filename: a filename with a json suffix

        Returns:
            None
        """
        with open(filename, "w") as f:
            f.write(safe_jsonify(self._data, pretty=True))

    @staticmethod
    def load_from_file(filename: str | bytes) -> CellsPerFeature:
        """Load the data from a file.

        Args:
            filename: the file to load from

        Returns:
            a new CallsPerFeature object
        """
        if filename is None:
            return CellsPerFeature.from_dictionary({})
        with open(filename) as f:
            data = json.load(f)
        return CellsPerFeature.from_dictionary(data)

    @staticmethod
    def from_dictionary(data: dict[str | bytes, list[str | bytes]]) -> CellsPerFeature:
        """Loads from a dictionary.

        Args:
            data (Dict[str|bytes,List[str|bytes]]): [TODO]

        Returns:
            CallsPerFeature: object with str->bytes
        """
        return_value = CellsPerFeature()

        for k, v in data.items():
            k = ensure_binary(k)
            v = [ensure_binary(x) for x in v]
            return_value[k] = v
        return return_value

    def get_multiplets(self):
        """Returns barcodes that occur under multiple keys in the dictionary.

        Returns:
            multiplets: set of bytes representing multiplet barcodes
        """
        multiplets = set()  # set(bytes)
        seen_bcs = set()  # set(bytes)
        for tag in self:
            for bc in self[tag]:
                if bc in seen_bcs:
                    # barcode was seen in another sample or fingerprint
                    multiplets.add(bc)
                seen_bcs.add(bc)

        return multiplets

    def get_set_of_unique_bcs(self):
        """Get a list of all the unique barcodes from all assigned features.

        Returns:
            A set of all assigned barcodes
        """
        return set(itertools.chain.from_iterable(self._data.values()))


class SampleBarcodes(CellsPerFeature):  # pylint: disable=too-many-ancestors
    """A class to hold the barcodes per sample.

    Acts as a dictionary like:

            "Sample1": ["BC1", "BC2", "BC3"],
            "Sample2": ["BC3", "BC4"],
            "Sample3": ["BC1", "BC3"],

    We want the contents to be bytes, following convention.
    We are just renaming CellsPerFeature because it does exactly what we need.
    """

    def get_multiplets(self):
        """Override the parent class' get_multiplets.

        Tag data are expected to overlap (multiplets)
        but SampleBarcodes should have non-overlapping lists of barcodes per sample.

        Returns:
            multiplets: set of bytes representing multiplet barcodes
        """
        return set()


class CategoryNames(NamedTuple):
    """To name outputs in FeatureAssignmentsMatrix.get_features_per_cell_table."""

    background: str
    ambiguous: str
    singles: str
    multiples: str


CMO_CATEGORY_NAMES = CategoryNames(
    "No tag molecules", "No tag assigned", "1 tag assigned", "More than 1 tag assigned"
)


class FeatureAssignmentsMatrix:
    """This is a sparse Pandas dataframe where rows are features and columns are barcodes.

    A (feature, barcode) pair has value 1 if that barcode is assigned that feature, else 0.
    """

    # Dataype used for the dataframe
    FEATURE_ASSIGNMENT_DTYPE = "uint8"

    def __init__(self, df: pd.DataFrame, matrix: CountMatrix):
        """Constructor for children classes."""
        self.df = df
        self.matrix = matrix
        self.features_per_cell_table = None

    def get_shape(self):
        """Return the shape of the matrix."""
        return self.df.shape

    def get_num_features_per_barcode(self):
        """Return number of assigned features per barcode."""
        num_features_per_barcode = self.df.sum(axis=0)
        return num_features_per_barcode

    def get_barcodes_per_feature(self, feature):
        """Return barcodes assigned to the given feature."""
        bcs = self.df.loc[feature]
        return bcs[bcs > 0].index.values

    def get_nlets(self, n):
        """Return barcodes assigned to n features.

        Unassigned barcodes would be n=0, doublets would be n=2, etc.
        """
        num_features_per_barcode = self.get_num_features_per_barcode()
        nlets = num_features_per_barcode[num_features_per_barcode == n].index.values
        return nlets

    def get_multiplets(self):
        """Return multiplet barcodes."""
        num_features_per_barcode = self.get_num_features_per_barcode()
        multiplets = num_features_per_barcode[num_features_per_barcode > 1].index.values
        return multiplets

    def filter_multiplets(self):
        """Return a new dataframe with multiplet barcodes filtered out."""
        multiplets = self.get_multiplets()
        singlets_df = self.df.drop(columns=multiplets)
        return singlets_df

    def get_cells_per_feature(self):
        """Convenience function that returns a CellsPerTag object.

        Basically a dict of (feature_id: [cell]) pairs.
        """
        cells_per_feature = CellsPerFeature()
        for feature in self.df.index:
            barcodes_per_feature = self.get_barcodes_per_feature(feature)
            cells_per_feature[feature] = barcodes_per_feature.tolist()
        return cells_per_feature

    def get_features_per_cell_table(self, sep=b"|"):
        """Returns a list of features assigned, with metadata.

        Returns a Pandas dataframe, indexed by bc, that provides a list of features assigned and
        metadata. Columns are ['num_features', 'feature_call', 'num_umis']
        Creates the "tag_calls_per_cell" output for Multiplexing Capture library
        """
        columns = [NUM_FEATURES, FEATURE_CALL, NUM_UMIS]
        features_per_cell_table = pd.DataFrame(columns=columns)
        self.matrix.tocsc()
        for cell in self.df.columns.values:
            calls = self.df[cell].to_numpy().nonzero()
            calls = self.df.index[calls].values

            num_features = len(calls)
            # TODO: Make a separate PR to include these cells
            if num_features == 0:
                continue

            # These three lines represent a fast way of populating the UMIS
            features = sep.join(calls)
            cell_index = self.matrix.bc_to_int(cell)
            bc_data = np.ravel(
                self.matrix.m.getcol(cell_index).toarray()
            )  # get data for a single bc, densify it, then flatten it
            umis = sep.join(
                [str(bc_data[self.matrix.feature_id_to_int(f_id)]).encode() for f_id in calls]
            )

            features_per_cell_table.loc[cell] = (num_features, features, umis)

        features_per_cell_table.index.name = CELL_BARCODE
        features_per_cell_table.sort_values(by=[FEATURE_CALL], inplace=True, kind="mergesort")
        return features_per_cell_table

    def _get_percentage_cells(self, numerator: int):
        """Generate values for pct_cells of the summary table.

        This is done by diving the number of cells in whichever
        category by the total number of cells, and round up to 1 decimal point
        of precision.
        """
        num_cells = len(self.matrix.bcs)
        percentage = 100 * tk_stats.robust_divide(numerator, num_cells)
        return np.round(percentage, 1)

    def get_feature_calls_summary(
        self,
        sep: bytes = b"|",
        cat_names: CategoryNames = CategoryNames(
            "No feature molecules",
            "No confident call",
            "1 feature assigned",
            "More than 1 feature assigned",
        ),
    ) -> pd.DataFrame:
        """Feature calls summary table."""
        if self.features_per_cell_table is None:
            self.features_per_cell_table = self.get_features_per_cell_table()
        num_cells = "num_cells"
        column_titles = [num_cells, "pct_cells", "median_umis", "stddev_umis"]
        feature_calls_summary = pd.DataFrame(columns=column_titles)

        num_cells_without_feature_umis = self.get_num_cells_without_molecules()
        feature_calls_summary.loc[cat_names.background] = (
            num_cells_without_feature_umis,
            self._get_percentage_cells(num_cells_without_feature_umis),
            "None",
            "None",
        )

        num_cells_without_features = self.get_nlets(0)
        feature_calls_summary.loc[cat_names.ambiguous] = (
            len(num_cells_without_features),
            self._get_percentage_cells(len(num_cells_without_features)),
            "None",
            "None",
        )

        num_singlets = self.get_nlets(1)
        feature_calls_summary.loc[cat_names.singles] = (
            len(num_singlets),
            self._get_percentage_cells(len(num_singlets)),
            "None",
            "None",
        )

        num_multiplets = self.get_multiplets()
        feature_calls_summary.loc[cat_names.multiples] = (
            len(num_multiplets),
            self._get_percentage_cells(len(num_multiplets)),
            "None",
            "None",
        )

        for f_call, table_iter in itertools.groupby(
            self.features_per_cell_table.itertuples(), key=lambda row: row.feature_call
        ):
            if isinstance(f_call, bytes):
                f_call = f_call.decode()
            assert isinstance(f_call, str)
            if f_call == "None":
                continue
            this_cells = 0
            umis_per_cell = []

            for row in table_iter:
                this_cells += 1
                if row.num_features > 1:
                    umis_per_cell.append(sum(float(x) for x in row.num_umis.split(sep)))
                else:
                    umis_per_cell.append(float(row.num_umis))

            feature_calls_summary.loc[f_call] = (
                this_cells,
                self._get_percentage_cells(this_cells),
                np.median(umis_per_cell),
                np.round(np.std(umis_per_cell), 1),
            )
        feature_calls_summary.index.name = "Category"
        feature_calls_summary[num_cells] = feature_calls_summary[num_cells].astype(dtype=np.int32)

        ## Sort groups with singlets first
        def _cmp(a, b):
            """Python 2 cmp."""
            return (a > b) - (a < b)

        sep_str = ensure_str(sep)

        def _cmp_mult(x, y):
            """Compare by numbers of sep then by string."""
            c1 = _cmp(x.count(sep_str), y.count(sep_str))
            if c1 == 0:
                return _cmp(x, y)
            else:
                return c1

        ordered_cats = sorted((x for x in feature_calls_summary.index), key=cmp_to_key(_cmp_mult))
        # Remove the category names from the sorted list and add back in later
        for cat in cat_names:
            ordered_cats.remove(cat)
        ordered_cats = list(cat_names) + ordered_cats
        feature_calls_summary = feature_calls_summary.reindex(ordered_cats)

        return feature_calls_summary

    def get_num_cells_without_molecules(self) -> int:
        """Returns the number of cells without molecules.

        Returns:
            int32: number of barcodes that do not have any feature umis
        """
        counts = self.matrix.get_subselected_counts(log_transform=False)
        return np.sum(counts == 0)


class BarcodeAssignmentException(Exception):
    """Exception when barcode overrides don't match cell calling."""


def validate_force_sample_barcodes(
    filtered_barcodes_csv, cells_per_tag_json, non_singlet_barcodes_json
):
    """Validate force_sample_barcodes.

    If the user overrides the sample deplex using a CSV file, we need to verify the barcodes in the
    CSV match the results from cell calling, or the behavior is undefined (and as of 3/14/2022
    leads to a downstream error).

    Returns:
        None or raises Exception
    """
    # Loads cell calling barcodes.
    filtered_barcodes = cr_utils.get_cell_associated_barcode_set(filtered_barcodes_csv)
    # Loads barcodes assigned to CMOs
    cmo_barcodes = SampleBarcodes.load_from_file(cells_per_tag_json).get_set_of_unique_bcs()
    # Loads Multiplet/Unassigned/Blank
    non_singlet_barcodes = SampleBarcodes.load_from_file(
        non_singlet_barcodes_json
    ).get_set_of_unique_bcs()
    all_forced_barcodes = cmo_barcodes.union(non_singlet_barcodes)
    # Forced but not cell called
    force_no_cell_call = all_forced_barcodes.difference(filtered_barcodes)
    cell_call_no_force = filtered_barcodes.difference(all_forced_barcodes)
    if force_no_cell_call or cell_call_no_force:
        msg = (
            "Mismatch between cell calling and Barcode Assignment CSV.  The Barcode Assignment CSV must have an entry "
            "for each barcode called by cell calling, and cannot have other barcodes within it.\n"
        )
        if force_no_cell_call:
            msg += "\nThe following barcodes were in the Barcode Assignment CSV but were not called as cells:\n"
            msg += "\n".join(x.decode() for x in force_no_cell_call)
        if cell_call_no_force:
            msg += "\nThe following barcodes were called as cells but were not in the Barcode Assignment CSV:\n"
            msg += "\n".join(x.decode() for x in cell_call_no_force)
        raise BarcodeAssignmentException(msg)


def recursive_logical_or(arrays):
    """Recursively run OR operation on array of arrays.

    >>> recursive_logical_or([[1]])
    [1]
    >>> recursive_logical_or([[0, 0, 1], [0, 1, 1]]).tolist()
    [0, 1, 1]
    >>> recursive_logical_or([[0, 0, 1], [0, 1, 1], [1, 0, 1]]).tolist()
    [1, 1, 1]
    """
    assert len(arrays) > 0
    if len(arrays) == 1:
        return arrays[0]
    elif len(arrays) == 2:
        return np.logical_or(arrays[0], arrays[1]).astype(int)
    elif len(arrays) > 2:
        return np.logical_or(arrays[0], recursive_logical_or(arrays[1:])).astype(int)


if __name__ == "__main__":
    testmod(name="recursive_logical_or", verbose=False)
