# Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

import cellranger.rna.library as rna_library
from cellranger.feature.utils import get_feature_counts_as_df

if TYPE_CHECKING:
    import pandas as pd

    from cellranger.matrix import CountMatrix

# fraction of reads necessary for a bc to be deemed highly corrected
HIGH_UMI_CORRECTION_THRESHOLD = 0.5
NUM_READS_THRESHOLD = 10000  # number of reads necessary for a bc to be deemed highly corrected
BACKGROUND_ANTIBODY_UMI_THRESHOLD = 1000  # if less UMIs, antibody not used for aggregate detection
TOP_UMI_BCS = 25  # top barcodes by antibody counts to look for aggregates

FRACTION_CORRECTED_READS = "frac_corrected_reads"
FRACTION_TOTAL_READS = "frac_total_reads"
FRACTION_SAMPLE_READS = "frac_sample_reads"


def augment_correction_table_with_corrected_reads_fraction(correction_data: pd.DataFrame):
    """Add a fractional corrected barcodes column to the barcode correction table.

    The dataframe is modified in-place.
    """
    correction_data[FRACTION_CORRECTED_READS] = correction_data["umi_corrected_reads"].divide(
        correction_data["reads"]
    )


def augment_correction_table_with_read_fraction(
    correction_data: pd.DataFrame, col_name: str = FRACTION_TOTAL_READS
):
    """Add a column to the correction data with fractional read data by library type.

    The dataframe is modified in-place.

    Adds fractional reads for each barcode, computed per library type, such that
    for each barcode this value is the fraction of total reads in the associated
    library type.
    """
    by_library_type = correction_data.groupby(rna_library.LIBRARY_TYPE)
    reads_by_library = by_library_type["reads"].sum().to_dict()

    correction_data["total_reads"] = correction_data["library_type"].map(reads_by_library)
    correction_data[col_name] = correction_data["reads"].divide(correction_data["total_reads"])
    del correction_data["total_reads"]


def filter_correction_table(correction_data: pd.DataFrame, library_type: str) -> pd.DataFrame:
    """Filter the correction table by library type.

    Also indexes the resulting dataframe by barcode.
    """
    filtered = correction_data[correction_data[rna_library.LIBRARY_TYPE] == library_type]
    filtered.set_index("barcode", inplace=True)
    return filtered


def subselect_augmented_table(barcode_list, augmented_table: pd.DataFrame):
    """Make a slice of the barcode summary csv with the given list of barcodes."""
    if isinstance(barcode_list, set):
        barcode_list = list(barcode_list)
    subselected_table = augmented_table.loc[barcode_list].drop(
        ["reads", "candidate_dup_reads"], axis=1
    )
    return subselected_table


def detect_outlier_umis_bcs(matrix, multiplier=3, lib_type=rna_library.ANTIGEN_LIBRARY_TYPE):
    """Find barcodes with UMI counts in excess of IQR * multiplier defined on the top 100 barcodes."""
    counts = matrix.select_features_by_type(lib_type).get_counts_per_bc()
    top100_idx = np.argsort(-counts)[:100]
    q3 = np.quantile(counts[top100_idx], 0.75)
    q1 = np.quantile(counts[top100_idx], 0.25)
    threshold = q3 + (q3 - q1) * multiplier
    # min cutoff=1000 umis to be labeled as aggregate
    if threshold < 1000:
        return []
    outlier_idx = top100_idx[counts[top100_idx] >= threshold]
    return matrix.select_barcodes(outlier_idx).bcs


def detect_highly_corrected_bcs(augmented_table: pd.DataFrame):
    """Given the augmented correction table, find barcodes that exceed.

    pre-defined thresholds of fraction corrected reads and fraction total reads
    """
    high = augmented_table[FRACTION_CORRECTED_READS].values > HIGH_UMI_CORRECTION_THRESHOLD
    high &= augmented_table["reads"] > NUM_READS_THRESHOLD
    highly_corrected_bcs = set(augmented_table.index[high])
    return highly_corrected_bcs


def _calculate_fraction_to_use(num_total_signal_antibodies):
    """Given the total number of antibodies used in this experiments, decide how many to use for aggregate detection.

    Use a linear model where coefficients satisfy the following constraints: 100% of the 5-antibody panel and
    60% of the 25 antibody panel need to be selected, and the fraction needs to decrease monotonically in between.
    m * 5 + b = 1.0
    m * 25 + b = 0.6

    Which amounts to the following conversion table:
    Total antibodies: 5  6  7  8  9  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25
    How many to keep: 5  6  7  8  8   9  10  11  11  12  12  13  13  14  14  14  15  15  15  15  15

    >>> _calculate_fraction_to_use(5)
    1.0
    >>> _calculate_fraction_to_use(17)
    0.76
    >>> _calculate_fraction_to_use(25)
    0.6
    >>> _calculate_fraction_to_use(30)
    0.6
    """
    assert num_total_signal_antibodies >= 5
    if num_total_signal_antibodies > 26:
        fraction_to_use = 0.6
    else:
        m = -0.02
        b = 1.1
        fraction_to_use = m * num_total_signal_antibodies + b
    return fraction_to_use


def detect_aggregate_barcodes(matrix: CountMatrix, *, num_probe_barcodes: int | None = None):
    """Given the antibody UMI counts, find barcodes which are enriched in a certain number of antibodies,.

    which is a dynamically adjusted fraction of the number of signal antibodies from the panel
    """
    ab_matrix = matrix.select_features_by_type(rna_library.ANTIBODY_LIBRARY_TYPE)
    raw_ab_counts = get_feature_counts_as_df(ab_matrix, get_transpose=True)
    del ab_matrix  # No longer used, save the memory
    # Not all antibodies in the panel get used, some do not get stained well, so drop such antibodies
    # Note: in code below we do df.values.sum(axis=0) instead of df.sum() to create making a copy of the
    # data that hogs memory.
    # See: https://github.com/pandas-dev/pandas/issues/16788
    background_columns = raw_ab_counts.columns[
        np.where(raw_ab_counts.values.sum(axis=0) < BACKGROUND_ANTIBODY_UMI_THRESHOLD)
    ]
    raw_ab_counts = raw_ab_counts.drop(labels=background_columns, axis=1)
    num_signal_antibodies = len(raw_ab_counts.columns)
    # cannot reliably differentiate low-order multiplets from aggregates, so do not run the algo
    if num_signal_antibodies < 5:
        return []

    # select a list of barcodes based on their total antibody UMI counts
    top_umi_barcodes = (num_probe_barcodes or 1) * TOP_UMI_BCS
    total_sum = raw_ab_counts.values.sum(axis=1)
    top_indices = np.argsort(total_sum)[-top_umi_barcodes:]
    del total_sum
    candidate_barcodes = raw_ab_counts.index[top_indices]
    del top_indices
    # loop over each antibody, and see if any of the suspicious barcodes are enriched in it
    # creates the following dictionary, where 1 means that antibody was found in high quantities
    # {TACGGGCAGGTATAGT-1: [1, 1, 1, 0, 1]}
    #                     CD3  CD4  CD8a CD19  PD-1
    # TACGGGCAGGTATAGT-1   1    1     1    0    1
    found_in_top_ab = {key: [] for key in candidate_barcodes}
    for col in raw_ab_counts.columns:
        # select top barcodes enriched in this antibody
        barcodes_high_in_this_ab = raw_ab_counts[col].sort_values()[-top_umi_barcodes:]
        for bc in candidate_barcodes:
            if bc in barcodes_high_in_this_ab:
                # if this candidate barcode is among those top enriched barcodes, add 1
                found_in_top_ab[bc].append(1)
            else:
                found_in_top_ab[bc].append(0)

    # now we know which antibodies every GEM is enriched in
    # the ones that have most such antibodies are aggregates
    fraction_to_use = _calculate_fraction_to_use(num_signal_antibodies)
    aggregate_barcodes = set()
    for bc in candidate_barcodes:
        if sum(found_in_top_ab[bc]) >= int(np.round(num_signal_antibodies * fraction_to_use)):
            aggregate_barcodes.add(bc)

    return aggregate_barcodes
