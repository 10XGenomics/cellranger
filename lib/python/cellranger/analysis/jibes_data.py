# Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
"""The JibesData class used by both the Rust and Python implementations."""
from __future__ import annotations

import numpy as np
import pandas as pd

BARCODE_COL = "Barcode"


class JibesData:
    """The raw data for a JIBES model."""

    # pylint: disable=too-few-public-methods
    def __init__(self, data, column_names, barcodes):
        """Initialize the data with a numpy matrix, with columns given the column names and.

        barcodes the row names.

        Basically a lightweight version of a pandas data.frame

        :param data:
        :param column_names:
        :param barcodes:
        """
        self.counts = data
        self.column_names = column_names
        self.barcodes = barcodes

    def subset_to_barcodes(self, valid_bcs):
        """Take only a subset of the barcodes to use.

        :param valid_bcs: an iterator of valid_bcs
        :return: a JibesData instance with a subset of all barcodes
        """
        okay = set(valid_bcs)
        to_use = [bc in okay for bc in self.barcodes]
        barcodes = self.barcodes[to_use]
        counts = self.counts[to_use, :]
        return JibesData(counts, self.column_names, barcodes)

    def get_df_from_data(self) -> pd.DataFrame:
        """Get a pandas data from of the counts and barcodes.

        Return:
          A data frame with counts.
        """
        data_dict = {
            x: np.asarray(self.counts[:, i]).ravel() for i, x in enumerate(self.column_names)
        }
        data_dict[BARCODE_COL] = self.barcodes
        df = pd.DataFrame.from_dict(data_dict)
        return df
