#!/usr/bin/env python
#
# Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
#
#
"""Generalization of the `FeatuerAssigner` class.

Generalizes the code for assigning features to cells.
"""
from __future__ import annotations

import itertools
from collections import OrderedDict, defaultdict

import numpy as np
import pandas as pd
import scipy.stats as sp_stats
from pandas.arrays import SparseArray
from scipy import optimize
from scipy.special import comb
from sklearn import mixture

import cellranger.matrix as cr_matrix
import cellranger.rna.library as rna_library
import tenkit.stats as tk_stats
from cellranger.analysis.combinatorics import generate_solutions, multinomial_comb
from cellranger.feature.feature_assignments import CellsPerFeature, FeatureAssignmentsMatrix
from cellranger.feature.throughputs import CORR_FACTOR, N_G
from cellranger.pandas_utils import sanitize_dataframe

SUPPORTED_FEATURE_TYPES = [
    rna_library.ANTIBODY_LIBRARY_TYPE,
    rna_library.ANTIGEN_LIBRARY_TYPE,
    rna_library.CRISPR_LIBRARY_TYPE,
    rna_library.MULTIPLEXING_LIBRARY_TYPE,
]

SUPPORTED_METHODS = ["GMM"]
FEATURE_SEPARATOR = b"|"

UMI_NUM_TRIES = 10  # Number of initial points to try for GMM-fitting
UMI_MIX_INIT_SD = 0.25  # Initial standard deviation for GMM components
MIN_COUNTS_PER_ANTIBODY = 1000  # Filter out background antibodies

# Filtering feature assignments based on UMI thresholds and correlation with other tags
COUNTS_DYNAMIC_RANGE = 50.0
MIN_COUNTS_PER_TAG = 1000
MAX_ALLOWED_PEARSON_CORRELATION = 0.5

_IS_CONTAMINANT_LC = "is_contaminant"
FEATURE_METRICS = {
    "UMI_COUNTS": "total_umi_counts",
    "UMI_THRESHOLD": "umi_threshold",
    "SNR": "log10_snr",
    "IS_CONTAMINANT": _IS_CONTAMINANT_LC,
    "CORRELATED": "source_tags",
}

OBSERVED_COL = "observed"
MEASURABLE_COL = "measurable"
POISSON_COL = "poisson_loading"

MINIMUM_DEPTH_TO_FILTER = 1000000
# if tag reads usable per cell is below this, do NOT do UMI filtering
# Essentially turns off UMI filtering as we won't get 1M usable reads per cell

NUM_TOTAL_TAGS = 14

# set measurable probabilities to 0 for poisson probabilities < this value
# else takes forever to compute
MIN_POISSON_PROB = 1e-5
MULTIPLET_ROW_INDEX = "all_multiplets"

FEATURE_CALL_COL = "feature_call"
CELL_BARCODE_COL = "cell_barcode"


def calculate_fat_tail_frac(
    freqs_df,
    start_at=4,
    report_prefix="",
    depth_suffix="",
):
    """Takes the marginal_tag_frequencies df and computes the fat-tail metric."""
    freqs_df.drop(labels=[MULTIPLET_ROW_INDEX], inplace=True)

    expected_freqs = freqs_df[MEASURABLE_COL].values
    observed_freqs = freqs_df[OBSERVED_COL].values
    num_lets = len(observed_freqs)
    num_recovered_cells = observed_freqs[0] + np.sum(
        observed_freqs[1:num_lets] * range(1, num_lets)
    )
    sub_lets = range(start_at, num_lets)
    diff_freqs = observed_freqs[sub_lets] - expected_freqs[sub_lets]
    diff_freqs[diff_freqs < 0] = 0
    mass_fat_tail = np.sum(diff_freqs * sub_lets)
    frac_fat_tail = float(mass_fat_tail) / num_recovered_cells

    return {f"{report_prefix}frac_fat_tail{depth_suffix}": frac_fat_tail}


def get_multiplet_counts_unrounded(obs_cells: float, n_gems: int = N_G) -> np.ndarray:
    """Returns unrounded counts for multiplets 1 and above."""
    num_loaded_cells = CORR_FACTOR * obs_cells
    # estimate the numnber of cells loaded from the number of recovered cells
    loading_rate = float(num_loaded_cells) / n_gems  # Poisson loading rate
    # expected counts of i-lets just due to Poisson loading
    poisson_fracs = TagAssigner.get_poisson_multiplets_fractions(loading_rate)
    poisson_counts = poisson_fracs * n_gems / CORR_FACTOR
    return poisson_counts


def get_multiplet_counts(obs_cells: float, n_gems: int = N_G) -> np.ndarray:
    return np.round(get_multiplet_counts_unrounded(obs_cells, n_gems))


class CellEstimateCantConvergeException(Exception):
    """Used to signal errors due to estimating the expected number of cell."""


def calculate_expected_total_cells(obs_barcodes: int, n_gems: int = N_G) -> float:
    """Solve for the observed cell counts.

    Given an observed number of cell barcodes, try to solve for
    the number of likely observed cells in the original mixture by estimating the number
    of different k-lets implied by that total number of barcodes.
    """

    def to_minimize(x):
        return np.power(obs_barcodes - np.sum(get_multiplet_counts_unrounded(x, n_gems)), 2.0)

    z = optimize.minimize(to_minimize, x0=obs_barcodes * 1.1)
    if not z.success:
        raise CellEstimateCantConvergeException(
            "Could not estimate number of cells based on number of barcodes"
        )
    if z.fun > 2:
        msg = "Could not converge on an accurate cell estimate."
        if obs_barcodes > 62500:
            msg += (
                " There is an unexpectedly high number of cells detected. "
                "Consider using or changing the parameter --force-cells or --expect-cells and rerunning the pipeline."
            )
        raise CellEstimateCantConvergeException(msg)
    return z.x[0]


class FeatureAssignments:  # pylint: disable=too-few-public-methods
    """A class for values returned by `FeatureAssigner.get_feature_assignments`.

    Replaces a dictionary with the following properties:

        {'bc_indices': np.array,
        'umi_counts': int,
        'is_contaminant': bool,
        'correlated_features': list}},

    where bc_indices is a np array that indexes into matrix.bcs,
    correlated_features is a list of feature_ids
    """

    def __init__(
        self,
        bc_indices: np.ndarray,
        umi_counts: int,
        is_contaminant: bool,
        correlated_features=None,
    ):
        self.bc_indices = bc_indices
        self.umi_counts = umi_counts
        self.is_contaminant = is_contaminant
        self.correlated_features = correlated_features

    def __getitem__(self, item):
        """This class replaces a bunch of untyped dictionaries.

        For backwards compatability we enable attr access this way.
        """
        return getattr(self, item)


def _calculate_n_let_diversity_probs(n_let: int, n_tags: int) -> np.ndarray:
    """Returns the probabilities that an nlet contains tag types.

    Given n_tags = len(tag_fractions) and a n_let,
    return a vector of probabilities that an n-let contains
    1, 2, 3, 4, etc. tag types

    Args:
        n_let: the integer number of tags in this multiplet
        n_tags: The number of tags

    Returns:
        A vector of probabilities
    """
    probs = np.zeros(n_let)
    if n_tags == 0:
        return probs
    base_prob = np.power(1.0 / float(n_tags), n_let)
    for x in range(1, n_let + 1):
        ways_to_get_k_let_with_x_tag_types: np.ndarray = comb(n_tags, x)
        total_combos = 0.0
        for combo in generate_solutions(x, n_let, True):
            total_combos += multinomial_comb(combo)
        probs[x - 1] = ways_to_get_k_let_with_x_tag_types * total_combos * base_prob
    return probs


def call_presence_with_gmm_ab(umi_counts: np.ndarray, *, min_umi_threshold: int = 1) -> np.ndarray:
    """Given the UMI counts for a specific antibody, separate signal from background.

    A cell must have at least `min_umi_threshold` UMIs for this feature to be considered positive.
    """
    if np.max(umi_counts) == 0 or max(umi_counts.shape) < 2:
        # there are no UMIs, or only one UMI, each barcode has 0 count
        return np.repeat(False, len(umi_counts))

    log_umi_counts = np.log10(1 + umi_counts).reshape(-1, 1)

    # Initialize and fit the gaussian Mixture Model
    gmm = mixture.GaussianMixture(n_components=2, n_init=10, covariance_type="tied", random_state=0)
    gmm.fit(log_umi_counts)
    positive_component = np.argmax(gmm.means_)

    # Classify each cell
    return (umi_counts >= min_umi_threshold) & (gmm.predict(log_umi_counts) == positive_component)


# This cannot use a namedtuple because those are immutable.
class AssignmentMetadata:
    """Metadata for feature assignments."""

    # pylint: disable=too-many-instance-attributes,too-few-public-methods

    def __init__(
        self,
        num_cells,
        num_cells_with_features,
        num_cells_without_feature_umis,
        num_singlets,
        num_multiplets,
    ):
        super().__init__()
        self.num_cells = num_cells
        self.num_cells_without_features = num_cells - num_cells_with_features
        assert self.num_cells_without_features >= 0
        self.frac_cells_without_features = tk_stats.robust_divide(
            self.num_cells_without_features, self.num_cells
        )
        self.num_cells_blanks = None
        self.frac_cells_blanks = None
        self.num_cells_unassigned = None
        self.frac_cells_unassigned = None
        self.num_cells_without_feature_umis = num_cells_without_feature_umis
        self.frac_cells_without_feature_umis = tk_stats.robust_divide(
            self.num_cells_without_feature_umis, self.num_cells
        )
        self.num_singlets = num_singlets
        self.frac_singlets = tk_stats.robust_divide(self.num_singlets, self.num_cells)
        self.num_multiplets = num_multiplets
        self.frac_multiplets = tk_stats.robust_divide(self.num_multiplets, self.num_cells)
        self.umi_thresholds = None
        self.threshold_cutoff = None
        self.freq_counts = None


class FeatureAssigner:
    """Assigns features to barcodes.

    Determines "true presence" of a feature, above background.
    """

    def __init__(
        self,
        sub_matrix: cr_matrix.CountMatrix,
        assignments: OrderedDict[bytes, FeatureAssignments] | None = None,
        features_per_cell: dict[bytes, list[bytes]] | None = None,
        features_per_cell_table: pd.DataFrame | None = None,
        assignment_metadata: AssignmentMetadata | None = None,
    ):
        """Constructs a `FeatureAssigner` object.

        A `FeatureAssigner` object identifies which barcodes "truly" contain
        the presence of features, above background. Serves as a parent class for specialized
        sub-classes such as `GuideAssigner`, `AntibodyOrAntigenAssigner` and `TagAssigner`.

        Args:
            sub_matrix (CountMatrix): filtered count matrix subset to features of interest.
            assignments (dict): initialized in child classes
            features_per_cell (dict): contains (barcode: [feature_assigned]) pairs
            features_per_cell_table (pd.DataFrame): indexed by bc, provides
                a list of features assigned and metadata (umi counts, num_features, etc.).
                Columns: ['num_features', 'feature_call', 'num_umis']
            assignment_metadata (AssignmentMetadata): contains metadata about
                feature assignments, such as the fraction with 1 feature,
                multiple features, etc.
        """
        assert isinstance(
            sub_matrix, cr_matrix.CountMatrix
        ), "matrix should be a CountMatrix object"
        self.sub_matrix = sub_matrix
        feature_types = self.sub_matrix.get_library_types()
        assert len(feature_types) == 1
        feature_type = next(iter(feature_types))
        assert feature_type in SUPPORTED_FEATURE_TYPES, f"feature_type {feature_type} not supported"

        self.assignments = assignments
        self.features_per_cell = features_per_cell
        self.features_per_cell_table = features_per_cell_table
        self.assignment_metadata = assignment_metadata

    def get_feature_assignments(self) -> dict[bytes, FeatureAssignments]:
        """Determines which bcs truly contain each feature above background.

        Returns:
            assignments (dict): {feature_id : FeatureAssignments object,
        """
        raise NotImplementedError("this method must be overridden")

    # TODO: this is now just a column in TagAssignmentsMatrix
    def get_features_per_cell(self) -> defaultdict[bytes, list[bytes]]:
        """Returns a dict of (barcode: [feature_assigned]) pairs."""
        if self.assignments is None:
            self.assignments = self.get_feature_assignments()

        features_per_cell: defaultdict[bytes, list[bytes]] = defaultdict(list)
        for feature_id, asst in self.assignments.items():
            assert isinstance(feature_id, bytes)
            # check for contaminat tags in tag librares
            if asst.is_contaminant:
                continue
            for i in asst.bc_indices:
                features_per_cell[self.sub_matrix.bcs[i]].append(feature_id)

        return features_per_cell

    # TODO: this can now be replaced with the same function from TagAssignmentsMatrix
    # TODO: do NOT modify this function further. it will be eventually removed.
    def get_cells_per_feature(self):
        """Convenience function that returns a dict of (feature_id: [cell]) pairs."""
        if self.assignments is None:
            self.assignments = self.get_feature_assignments()

        cells_per_feature = CellsPerFeature()
        for f_id, asst in self.assignments.items():
            # check for contaminat tags in tag librares
            if asst.is_contaminant:
                continue
            cells_per_feature[f_id] = [self.sub_matrix.bcs[i] for i in asst.bc_indices]
        return cells_per_feature

    # TODO: this can now be replaced with the same function from TagAssignmentsMatrix
    # TODO: do NOT modify this function further. it will be eventually removed.
    def get_features_per_cell_table(self, sep: bytes = FEATURE_SEPARATOR) -> pd.DataFrame:
        """Returns a Pandas dataframe of assigned features.

        The dataframe is indexed by bc, and provides a list of features assigned and
        metadata.

        Columns are ['num_features', 'feature_call', 'num_umis']
        """
        if self.features_per_cell is None:
            self.features_per_cell = self.get_features_per_cell()

        columns = ["num_features", FEATURE_CALL_COL, "num_umis"]
        features_per_cell_table = pd.DataFrame(columns=columns)
        for cell in self.features_per_cell:
            calls = self.features_per_cell.get(cell)
            num_features = len(calls)
            features = sep.join(calls)

            # These three lines represent a fast way of populating the UMIS
            cell_bc_index = self.sub_matrix.bc_to_int(cell)
            bc_data = np.ravel(
                self.sub_matrix.m.getcol(cell_bc_index).toarray()
            )  # get data for a single bc, densify it, then flatten it
            umis = sep.join(
                str(bc_data[self.sub_matrix.feature_id_to_int(f_id)]).encode() for f_id in calls
            )

            features_per_cell_table.loc[cell] = (num_features, features, umis)

        features_per_cell_table.index.name = CELL_BARCODE_COL
        features_per_cell_table.sort_values(by=[FEATURE_CALL_COL], inplace=True, kind="mergesort")
        return features_per_cell_table

    def compute_generic_assignment_metadata(self):
        """Returns metadata for feature assignments."""
        if self.features_per_cell_table is None:
            self.features_per_cell_table = self.get_features_per_cell_table()

        return AssignmentMetadata(
            num_cells=self.sub_matrix.bcs_dim,
            num_cells_with_features=self.features_per_cell_table.shape[0],
            num_cells_without_feature_umis=self.get_num_cells_without_molecules(),
            num_singlets=self.features_per_cell_table[
                self.features_per_cell_table["num_features"] == 1
            ].shape[0],
            num_multiplets=self.features_per_cell_table[
                self.features_per_cell_table["num_features"] > 1
            ].shape[0],
        )

    @staticmethod
    def _compute_umi_thresholds(
        matrix: cr_matrix.CountMatrix,
        assts: dict[bytes, FeatureAssignments],
        multiplexing: bool = False,
        sep: bytes = FEATURE_SEPARATOR,
    ):
        """Computes, for each feature, the lowest count of that feature.

        Taken from amongst cells assigned that feature.

        Args:
            matrix (CountMatrix): object
            assignments (dict): {feature_id : {'bc_indices': np.array,
                                               'umi_counts': int,
                                               'is_contaminant': bool,
                                               'correlated_features': list}},
                where bc_indices is a np array that indexes into matrix.bcs,
                correlated_features is a list of feature_ids

        Returns:
            pd.DataFrame: A summary of feature assignment thresholds
                 for both assigned features and features that are considered as
                 contaminants. The values listed in the "UMI threshold" column
                 are minimal UMI count for a feature amongst cells assigned to
                 that feature
        """
        # CAUTION! The dataframe is parsed on the rust side, any changes
        # to the format should lead to updates in write_websummary_json.rs

        if multiplexing:
            columns = [
                FEATURE_METRICS["UMI_COUNTS"],
                FEATURE_METRICS["UMI_THRESHOLD"],
                FEATURE_METRICS["SNR"],
                FEATURE_METRICS["IS_CONTAMINANT"],
                FEATURE_METRICS["CORRELATED"],
            ]
        else:
            columns = [FEATURE_METRICS["UMI_COUNTS"], FEATURE_METRICS["UMI_THRESHOLD"]]

        f_ids = sorted(assts.keys())
        umi_thresholds = sanitize_dataframe(pd.DataFrame(columns=columns))
        for f_id in f_ids:
            if isinstance(f_id, str):
                str_f_id = f_id
            else:
                assert isinstance(f_id, bytes)
                str_f_id = f_id.decode()
            if len(assts[f_id].bc_indices) == 0:
                continue

            umi_counts = matrix.get_subselected_counts(list_feature_ids=[f_id], log_transform=False)
            thres = int(np.amin(umi_counts[assts[f_id].bc_indices]))
            umi_thresholds.loc[str_f_id, FEATURE_METRICS["UMI_THRESHOLD"]] = thres
            umi_thresholds.loc[str_f_id, FEATURE_METRICS["UMI_COUNTS"]] = assts[f_id].umi_counts
            if multiplexing:
                peak_interval = FeatureAssigner._get_umi_interval(umi_counts, thres)
                umi_thresholds.loc[str_f_id, FEATURE_METRICS["SNR"]] = np.round(
                    np.log10(peak_interval[1] + 1) - np.log10(peak_interval[0] + 1), 3
                )
                umi_thresholds.loc[str_f_id, FEATURE_METRICS["IS_CONTAMINANT"]] = assts[
                    f_id
                ].is_contaminant
                if not assts[f_id].is_contaminant:
                    umi_thresholds.loc[str_f_id, FEATURE_METRICS["CORRELATED"]] = "NA"
                elif len(assts[f_id].correlated_features) > 0:
                    umi_thresholds.loc[str_f_id, FEATURE_METRICS["CORRELATED"]] = sep.join(
                        assts[f_id].correlated_features
                    )
                else:
                    umi_thresholds.loc[str_f_id, FEATURE_METRICS["CORRELATED"]] = "unknown"

        umi_thresholds = umi_thresholds.astype(
            {FEATURE_METRICS["UMI_THRESHOLD"]: int, FEATURE_METRICS["UMI_COUNTS"]: int}
        )
        return umi_thresholds

    @staticmethod
    def _get_umi_interval(counts, threshold, pct_width=25) -> tuple[float, float]:
        """Returns the Q3 and Q1 percentiles.

        Returns:
            float: the Q3 of counts < threshold
            float: the Q1 of the counts > threshold
        """
        noise = [x for x in counts if x < threshold]
        signal = [x for x in counts if x >= threshold]
        # set value to 0 if there is no noise or no signal (SNR would be negative in this case)
        q3_noise = np.percentile(noise, 100 - pct_width) if len(noise) > 0 else 0
        q1_signal = np.percentile(signal, pct_width) if len(signal) > 0 else 0
        return q3_noise, q1_signal

    @staticmethod
    def get_generic_feature_assignment_metrics(
        assignment_metadata: AssignmentMetadata,
        feature_bc_name: str,
        report_prefix: str,
        depth_suffix: str = "",
    ) -> dict[str, float]:
        """Compute summary metrics on feature assignments.

        Args:
            assignment_metadata: source metadata
            feature_bc_name: TODO
            report_prefix: TODO
            depth_suffix (str): the depth, if any, subsampled to in
                units of 1000 reads per cell. E.g. "1.0k". "" if full-depth

        Returns:
            dict: json_key of a metric and the corrosponding values.
        """
        name_fr_ce_w_features = report_prefix + f"frac_cells_with_{feature_bc_name}{depth_suffix}"

        frac_cells_with_features = (
            assignment_metadata.frac_singlets + assignment_metadata.frac_multiplets
        )

        name_fr_ce_w_sg_features = (
            report_prefix + f"frac_cells_with_single_{feature_bc_name}{depth_suffix}"
        )
        frac_cells_with_single_features = assignment_metadata.frac_singlets

        name_fr_ce_w_mult_features = (
            report_prefix + f"frac_cells_with_multiple_{feature_bc_name}{depth_suffix}"
        )
        frac_cells_multiple_features = assignment_metadata.frac_multiplets

        res = {
            name_fr_ce_w_features: frac_cells_with_features,
            name_fr_ce_w_sg_features: frac_cells_with_single_features,
            name_fr_ce_w_mult_features: frac_cells_multiple_features,
        }

        return res

    # TODO: do NOT modify this function further. it will be eventually removed.
    def get_generic_feature_calls_summary(
        self, assignment_metadata, feature_mol_name, feature_bc_name, sep=FEATURE_SEPARATOR
    ):
        """Returns a Pandas dataframe that summarizes feature assignment metrics.

        Creates the "protospacer_calls_summary" output for CRISPR
        """
        if self.features_per_cell_table is None:
            self.features_per_cell_table = self.get_features_per_cell_table()

        column_titles = ["num_cells", "pct_cells", "median_umis", "stddev_umis"]
        feature_calls_summary = pd.DataFrame(columns=column_titles)

        feature_calls_summary.loc[f"No {feature_mol_name} molecules"] = (
            assignment_metadata.num_cells_without_feature_umis,
            100 * assignment_metadata.frac_cells_without_feature_umis,
            "None",
            "None",
        )
        feature_calls_summary.loc["No confident call"] = (
            assignment_metadata.num_cells_without_features,
            100 * assignment_metadata.frac_cells_without_features,
            "None",
            "None",
        )

        feature_calls_summary.loc[f"1 {feature_bc_name} assigned"] = (
            assignment_metadata.num_singlets,
            100 * assignment_metadata.frac_singlets,
            "None",
            "None",
        )

        feature_calls_summary.loc[f"More than 1 {feature_bc_name} assigned"] = (
            assignment_metadata.num_multiplets,
            100 * assignment_metadata.frac_multiplets,
            "None",
            "None",
        )

        for f_call, table_iter in itertools.groupby(
            self.features_per_cell_table.itertuples(), key=FeatureAssigner.sort_by_feature_call
        ):
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
                100 * tk_stats.robust_divide(this_cells, self.assignment_metadata.num_cells),
                np.median(umis_per_cell),
                np.std(umis_per_cell),
            )
        feature_calls_summary.index.name = FEATURE_CALL_COL
        return feature_calls_summary

    def get_num_cells_without_molecules(self) -> int:
        """Returns the number of cells without molecules.

        Returns:
            int32: number of barcodes that do not have any feature umis
        """
        counts = self.sub_matrix.get_subselected_counts(log_transform=False, library_type=None)
        return np.sum(counts == 0)

    @staticmethod
    def sort_by_feature_call(row):
        """Returns feature_call field of a row from a Pandas df.

        Used in sorting.
        """
        return row.feature_call

    def create_feature_assignments_matrix(self) -> FeatureAssignmentsMatrix:
        """Create FeatureAssignments object.

        This is a data structure where rows are features and
        columns are barcodes.  In each cell, the value is 1 if barcode was
        assigned this feature, otherwise 0.
        """
        if self.assignments is None:
            self.assignments = self.get_feature_assignments()

        feature_assignments_dict = {}
        for f_id, asst in self.assignments.items():
            if not isinstance(f_id, bytes):
                raise ValueError(f"Feature ID {f_id} must be bytes, but was {type(f_id)}")
            if asst.is_contaminant:  # only relevant for multiplexing library
                continue

            mask = np.zeros(len(self.sub_matrix.bcs))
            if len(asst.bc_indices) > 0:
                # Put 1 for barcodes that were assigned to this feature. Basically, if
                # bc_indices = [1, 4], go from  [0, 0, 0, 0, 0] -> [0, 1, 0, 0, 1]
                np.put(mask, asst.bc_indices, 1)
            # instead of creating a dense df and then putting it in a sparse matrix,
            # better to create sparse arrays and put those into sparse df
            sparse_mask = SparseArray(mask, dtype=FeatureAssignmentsMatrix.FEATURE_ASSIGNMENT_DTYPE)
            feature_assignments_dict[f_id] = sparse_mask
            del mask, sparse_mask
        feature_assignments_df = pd.DataFrame(feature_assignments_dict)
        del feature_assignments_dict
        # when no assignments for whatever reason, simply return the object with empty df
        if len(feature_assignments_df) == 0:
            return FeatureAssignmentsMatrix(feature_assignments_df, self.sub_matrix)

        feature_assignments_df.index = self.sub_matrix.bcs
        feature_assignments_df = feature_assignments_df.T
        feature_assignments_matrix = FeatureAssignmentsMatrix(
            feature_assignments_df, self.sub_matrix
        )
        return feature_assignments_matrix


class GuideAssigner(FeatureAssigner):
    """Sub-class of FeatureAssigner specific to CRISPR Library features."""

    def __init__(
        self,
        matrix: cr_matrix.CountMatrix,
        *,
        min_crispr_umi_threshold: int,
    ):
        super().__init__(
            matrix,
            assignments=None,
            features_per_cell=None,
            features_per_cell_table=None,
            assignment_metadata=None,
        )
        self.report_prefix = rna_library.get_library_type_metric_prefix(
            rna_library.CRISPR_LIBRARY_TYPE
        )
        self.feature_mol_name = "guide"
        self.feature_bc_name = "protospacer"
        self.min_crispr_umi_threshold = min_crispr_umi_threshold

        self.method = "GMM"
        assert self.method in SUPPORTED_METHODS, f"Method {self.method} not supported"

    def get_feature_assignments(self) -> dict[bytes, FeatureAssignments]:
        return self.get_guide_assignments()

    def get_guide_assignments(self) -> dict[bytes, FeatureAssignments]:
        """Get feature assignments for CRISPR library."""
        assignments: dict[bytes, FeatureAssignments] = OrderedDict()
        feature_ids = [f.id for f in self.sub_matrix.feature_ref.feature_defs]
        self.sub_matrix.tocsr()
        for feature_id in feature_ids:
            if not isinstance(feature_id, bytes):
                raise ValueError(
                    f"Feature ID {feature_id} must be bytes but was {type(feature_id)}"
                )
            umi_counts = self.sub_matrix.get_subselected_counts(
                log_transform=False, list_feature_ids=[feature_id]
            )

            in_high_umi_component = GuideAssigner._call_presence(
                umi_counts, self.method, min_crispr_umi_threshold=self.min_crispr_umi_threshold
            )
            assignments[feature_id] = FeatureAssignments(
                np.flatnonzero(np.array(in_high_umi_component)), sum(umi_counts), False, None
            )
        return assignments

    @staticmethod
    def _call_presence(
        counts: np.ndarray,
        method: str = "GMM",
        *,
        min_crispr_umi_threshold: int,
    ) -> np.ndarray:
        """Classify each cell as positive/negative for a CRISPR feature using a GMM.

        A cell must have at least `min_crispr_umi_threshold` UMIs for this CRISPR feature to be
        considered positive. This threshold is used to exclude CRISPR features with only background
        signal and no foreground signal. Without this filter, a CRISPR feature with 0 or 1 UMIs
        in each cell would call all the cells with one UMI as positive, which renders meaningless
        the metric `Cells with one or more protospacers detected`. The threshold value was chosen
        by being the smallest sufficient value on experimental data.

        Args:
            counts: feature counts
            method: the method to use. Only supports GMMs right now

        Returns:
            Booleans indicating whether feature is present above background
        """
        if method == "GMM":
            return call_presence_with_gmm_ab(counts, min_umi_threshold=min_crispr_umi_threshold)
        raise ValueError(f"Method {method} is not supported")

    def create_guide_assignments_matrix(self) -> FeatureAssignmentsMatrix:
        return self.create_feature_assignments_matrix()

    def compute_assignment_metadata(self) -> AssignmentMetadata:
        """Compute assignment metadata for CRISPR library."""
        assignment_metadata = self.compute_generic_assignment_metadata()

        assert self.assignments is not None
        umi_thresholds = self._compute_umi_thresholds(self.sub_matrix, self.assignments)
        # keep the crispr output unchanged
        assignment_metadata.umi_thresholds = umi_thresholds["umi_threshold"].to_dict()
        return assignment_metadata

    def get_feature_assignment_metrics(self, depth_suffix: str = ""):
        """Feature assignment metrics for CRISPR library."""
        if self.assignment_metadata is None:
            self.assignment_metadata = self.compute_assignment_metadata()
        return self.get_generic_feature_assignment_metrics(
            self.assignment_metadata, self.feature_bc_name, self.report_prefix, depth_suffix
        )

    def get_feature_calls_summary(self) -> pd.DataFrame:
        """Feature calls summary table for CRISPR library."""
        if self.assignment_metadata is None:
            self.assignment_metadata = self.compute_assignment_metadata()
        return FeatureAssigner.get_generic_feature_calls_summary(
            self, self.assignment_metadata, self.feature_mol_name, self.feature_bc_name
        )


class AntibodyOrAntigenAssigner(FeatureAssigner):
    """Sub-class of FeatureAssigner specific to Antibody or Antigen Library features."""

    def __init__(self, matrix: cr_matrix.CountMatrix, feature_type: str):

        super().__init__(
            matrix.select_features_by_type(feature_type),
            assignments=None,
            features_per_cell=None,
            features_per_cell_table=None,
            assignment_metadata=None,
        )
        self.report_prefix = rna_library.get_library_type_metric_prefix(feature_type)

        if feature_type == rna_library.ANTIGEN_LIBRARY_TYPE:
            self.feature_mol_name = self.feature_bc_name = "antigen"
        else:
            self.feature_mol_name = self.feature_bc_name = "antibody"
        self.method = "GMM"
        assert self.method in SUPPORTED_METHODS, f"Method {self.method} not supported"
        self.assignment_list = [b"CD3", b"CD19", b"CD14"]

    def get_feature_assignments(self):
        return self.get_antibody_assignments()

    def get_antibody_assignments(self) -> dict[bytes, FeatureAssignments]:
        """Get feature assignments for Antibody library."""
        assignments: OrderedDict[bytes, FeatureAssignments] = OrderedDict()
        feature_ids = [f.id for f in self.sub_matrix.feature_ref.feature_defs]

        for feature_id in feature_ids:
            assert isinstance(feature_id, bytes)
            if feature_id not in self.assignment_list:
                continue
            umi_counts = self.sub_matrix.get_subselected_counts(
                log_transform=False, list_feature_ids=[feature_id]
            )
            total_count = umi_counts.sum()
            if total_count < MIN_COUNTS_PER_ANTIBODY:
                continue

            in_high_umi_component = AntibodyOrAntigenAssigner._call_presence(
                umi_counts, self.method
            )
            assignments[feature_id] = FeatureAssignments(
                np.flatnonzero(np.array(in_high_umi_component)), total_count, False, None
            )

        return assignments

    @staticmethod
    def _call_presence(counts, method: str = "GMM") -> np.ndarray:
        """Calls the appropriate method on counts.

        Args:
            counts (np.array): feature counts (int32)
            method (str): the method to use. Only supports GMMs right now

        Returns:
            np.array: Booleans indicating whether feature is "truly present", above background
        """
        if method == "GMM":
            return call_presence_with_gmm_ab(counts)
        raise ValueError(f"Method {method} is not supported")

    def create_ab_assignments_matrix(self):
        return self.create_feature_assignments_matrix()

    def compute_assignment_metadata(self):
        """Compute assignment metadata for Antibody library."""
        assignment_metadata = self.compute_generic_assignment_metadata()

        umi_thresholds = self._compute_umi_thresholds(self.sub_matrix, self.assignments)
        # keep the antibody output unchanged
        assignment_metadata.umi_thresholds = umi_thresholds["umi_threshold"].to_dict()
        return assignment_metadata

    def get_feature_assignment_metrics(self, depth_suffix: str = "") -> dict[str, float]:
        """Feature assignment metrics for Antibody library."""
        if self.assignment_metadata is None:
            self.assignment_metadata = self.compute_assignment_metadata()
        return self.get_generic_feature_assignment_metrics(
            self.assignment_metadata, self.feature_bc_name, self.report_prefix, depth_suffix
        )


class TagAssigner(FeatureAssigner):
    """Sub-class of FeatureAssigner specific to Multiplexing Capture Library features.

    Args:
        assignments:  {feature_id : {'bc_indices': np.array,
                                 'umi_counts': int,
                                 'is_contaminant': bool,
                                 'correlated_features': list}},
                where `bc_indices` is a np array that indexes into `matrix.bcs`,
                `correlated_features` is a list of feature_ids
    """

    def __init__(self, matrix: cr_matrix.CountMatrix, feature_type: str, n_gems: int = N_G):
        FeatureAssigner.__init__(
            self,
            matrix,
            assignments=None,
            features_per_cell=None,
            features_per_cell_table=None,
            assignment_metadata=None,
        )
        self.report_prefix = rna_library.get_library_type_metric_prefix(feature_type)
        self.feature_mol_name = self.feature_bc_name = "tag"
        self.method = "GMM"
        self.n_gems = n_gems
        assert self.method in SUPPORTED_METHODS, f"Method {self.method} not supported"

    def get_feature_assignments(self) -> dict[bytes, FeatureAssignments]:
        return self.get_tag_assignments()

    def get_tag_assignments(self) -> dict[bytes, FeatureAssignments]:
        """Get feature assignments for Tag library."""
        assignments: OrderedDict[bytes, FeatureAssignments] = OrderedDict()
        feature_ids = [f.id for f in self.sub_matrix.feature_ref.feature_defs]

        for feature_id in feature_ids:
            assert isinstance(feature_id, bytes)
            umi_counts = self.sub_matrix.get_subselected_counts(
                log_transform=False, list_feature_ids=[feature_id]
            )

            in_high_umi_component = TagAssigner._call_presence(umi_counts, self.method)
            assignments[feature_id] = FeatureAssignments(
                np.flatnonzero(np.array(in_high_umi_component)), sum(umi_counts), False, []
            )

        assignments = self.identify_contaminant_tags(assignments)
        return assignments

    @staticmethod
    def _call_presence(counts, method: str = "GMM") -> np.ndarray:
        """Calls the appropriate method on counts.

        Args:
            counts (np.array): feature counts (int32)
            method (str): the method to use. Only supports GMMs right now

        Returns:
            np.array: Booleans indicating whether feature is "truly present", above background
        """
        if method == "GMM":
            return call_presence_with_gmm_ab(counts)
        raise ValueError(f"Method {method} is not supported")

    def create_tag_assignments_matrix(self) -> FeatureAssignmentsMatrix:
        return self.create_feature_assignments_matrix()

    def compute_assignment_metadata(self) -> AssignmentMetadata:
        """Compute assignment metadata for Tag library."""
        assignment_metadata = self.compute_generic_assignment_metadata()

        freqs_df = self.compute_assignment_freqs(assignment_metadata.num_cells_without_features)
        assignment_metadata.freq_counts = freqs_df

        umi_thresholds = self._compute_umi_thresholds(
            self.sub_matrix, self.assignments, multiplexing=True
        )
        assignment_metadata.umi_thresholds = umi_thresholds

        return assignment_metadata

    def get_feature_assignment_metrics(
        self,
        depth_suffix: str = "",
    ) -> dict[str, float]:
        """Feature assignment metrics for Tag library."""
        if self.assignment_metadata is None:
            self.assignment_metadata = self.compute_assignment_metadata()

        res = self.get_generic_feature_assignment_metrics(
            self.assignment_metadata, self.feature_bc_name, self.report_prefix, depth_suffix
        )

        freqs_df = self.assignment_metadata.freq_counts.drop(
            labels=[MULTIPLET_ROW_INDEX], inplace=False
        )
        res.update(
            TagAssigner._compute_efficiency_metrics(
                freqs_df[MEASURABLE_COL].values,
                freqs_df[OBSERVED_COL].values,
                self.report_prefix,
                depth_suffix,
            )
        )
        res.update(self._compute_no_tag_metric(self.report_prefix, depth_suffix))
        return res

    def _compute_no_tag_metric(self, report_prefix: str = "", depth_suffix: str = ""):
        return {
            f"{report_prefix}frac_no_tag_assigned{depth_suffix}": self.assignment_metadata.frac_cells_without_features
        }

    def _compute_no_tag_molecules_metric(self, report_prefix=""):
        return {
            f"{report_prefix}frac_no_tag_molecules": self.assignment_metadata.frac_cells_without_feature_umis
        }

    @staticmethod
    def _compute_efficiency_metrics(
        expected_freqs: np.ndarray,
        observed_freqs: np.ndarray,
        report_prefix: str = "",
        depth_suffix: str = "",
    ):
        """Computes sc recovery efficiency.

        Given expected and observed freqs.

        Assumes inputs are arrays that list expected and observed frequencies
        starting from 0 tags.
        """
        singlets_expected = expected_freqs[1]
        singlets_observed = observed_freqs[1]

        multiplets_expected = sum(expected_freqs[2:])
        multiplets_observed = sum(observed_freqs[2:])

        # expect < 1 because of loss due to fat tail and some singlets being assigned 0 tags
        # if > 1 we are not calling multiplets sensitively enough
        sc_rec_efficiency = tk_stats.robust_divide(singlets_observed, singlets_expected)
        mtp_rec_efficiency = tk_stats.robust_divide(multiplets_observed, multiplets_expected)

        res = {
            f"{report_prefix}expected_num_singlets{depth_suffix}": singlets_expected,
            f"{report_prefix}observed_num_singlets{depth_suffix}": singlets_observed,
            f"{report_prefix}expected_num_multiplets{depth_suffix}": multiplets_expected,
            f"{report_prefix}observed_num_multiplets{depth_suffix}": multiplets_observed,
            f"{report_prefix}sc_rec_efficiency{depth_suffix}": sc_rec_efficiency,
            f"{report_prefix}mtp_rec_efficiency{depth_suffix}": mtp_rec_efficiency,
        }

        return res

    def get_feature_calls_summary(self):
        """Feature calls summary table for Tag library."""
        if self.assignment_metadata is None:
            self.assignment_metadata = self.compute_assignment_metadata()
        return FeatureAssigner.get_generic_feature_calls_summary(
            self, self.assignment_metadata, self.feature_mol_name, self.feature_bc_name
        )

    def compute_assignment_freqs(
        self, num_zero_features: int, num_total_tags: int = NUM_TOTAL_TAGS
    ) -> pd.DataFrame:
        """Compute assignment frequencies.

        Args:
            num_zero_features: the pre-computed number of cells assigned exactly 0 tags.
            num_total_tags: The number of rows to produce.

        Returns:
            pd.DataFrame: with the following columns:
                n: ranges from 0 to num_total_tags and
                poisson_loading: number of cells being n-lets in theory
                measurable: number of cells that should be measured to have n tags
                observed: number of cells assigned with n tags
        """
        freq_df = pd.DataFrame(
            columns=[POISSON_COL, MEASURABLE_COL, OBSERVED_COL], index=range(num_total_tags + 1)
        )

        freq_df.loc[0, :] = [0, 0, num_zero_features]

        # Calculate observed frequencies (note this does not include an entry for i = 0)
        obs_freqs = np.array([self._get_freq_num_features(i) for i in range(1, num_total_tags + 1)])

        num_cells_called = len(self.sub_matrix.bcs)
        estimated_cells = calculate_expected_total_cells(num_cells_called, self.n_gems)
        # estimate of total cells actually present, while accounting for the fact that
        # a k-let probably contains k-cells in a single partition

        poisson_counts = get_multiplet_counts(
            estimated_cells, self.n_gems
        )  # expected counts, according to Poisson statistics

        num_loaded_cells = CORR_FACTOR * estimated_cells
        # estimate the numnber of cells loaded from the number of recovered cells
        loading_rate = float(num_loaded_cells) / self.n_gems  # Poisson loading rate

        measurable_fracs = self.get_measurable_multiplets_fractions(loading_rate)
        measurable_counts = np.round(
            TagAssigner.get_partitions_from_fractions(measurable_fracs, num_gems=self.n_gems)
        )

        for i in range(1, num_total_tags + 1):
            freq_df.loc[i, POISSON_COL] = poisson_counts[i - 1]
            freq_df.loc[i, MEASURABLE_COL] = measurable_counts[i - 1]
            freq_df.loc[i, OBSERVED_COL] = obs_freqs[i - 1]

        freq_df.loc[MULTIPLET_ROW_INDEX] = [
            np.sum(freq_df[POISSON_COL]) - freq_df.loc[1, POISSON_COL],
            np.sum(freq_df[MEASURABLE_COL]) - freq_df.loc[1, MEASURABLE_COL],
            np.sum(freq_df[OBSERVED_COL]) - freq_df.loc[1, OBSERVED_COL] - num_zero_features,
        ]
        return freq_df

    def _get_freq_num_features(self, num_features: int) -> int:
        """Computes the number of cells assigned a specified number of features.

        e.g. if num_features = 3, returns the number of cells assigned exactly 3 features

        Note: does not work correctly for `num_features = 0` because
              `features_per_cell_table` only has entries for cells assigned
              at least 1 feature
        """
        assert isinstance(num_features, int), "num_features should be of type int"
        assert num_features > 0, "method only works for num_features > 0"

        if self.features_per_cell_table is None:
            self.features_per_cell_table = self.get_features_per_cell_table()

        return self.features_per_cell_table[
            self.features_per_cell_table["num_features"] == num_features
        ].shape[0]

    def identify_contaminant_tags(
        self,
        assignments: dict[bytes, FeatureAssignments],
        min_counts: int = MIN_COUNTS_PER_TAG,
        dynamic_range: float = COUNTS_DYNAMIC_RANGE,
        max_cor: float = MAX_ALLOWED_PEARSON_CORRELATION,
    ) -> dict[bytes, FeatureAssignments]:
        """Filter out tags that are likely contaminants based on total UMI counts.

        Fileters across the
        cell barcodes and the correlation of counts per cell with tags with more UMI counts.

        Args:
            assignments (Dict[bytes, FeatureAssignments]): with keys being tag IDs
            min_counts (int, optional): tags with total counts less than minimal will be
            contaminants. Defaults to MIN_COUNTS_PER_TAG.
            dynamic_range (float, optional): tags with total counts below this range from the
            top will be contaminants. Defaults to COUNTS_DYNAMIC_RANGE.
            max_cor (float, optional): tags with a Pearson correlation with a more abundant
            tag larger than it will be contaminants. Defaults to MAX_ALLOWED_PEARSON_CORRELATION.

        Returns:
            Dict[bytes, FeatureAssignments]
        """
        max_count = max(assignments[x].umi_counts for x in assignments)
        contaminant_tags = set()

        all_features = sorted(assignments.keys())
        for feature_id in all_features:
            count = assignments[feature_id].umi_counts
            if count == 0:
                assignments.pop(feature_id, None)
                contaminant_tags.add(feature_id)
            elif count < min_counts or count * dynamic_range < max_count:
                assignments[feature_id].is_contaminant = True
                contaminant_tags.add(feature_id)

        all_features = sorted(assignments.keys())

        # correlation check
        df = pd.DataFrame(
            self.sub_matrix.m.toarray(),
            index=[f.id for f in self.sub_matrix.feature_ref.feature_defs],
            columns=self.sub_matrix.bcs,
        ).T
        df = df.loc[:, all_features]
        cor_df = df.corr()
        for i in range(len(all_features) - 1):
            for j in range(i + 1, len(all_features)):
                this_cor = cor_df.loc[all_features[i], all_features[j]]
                if np.isnan(this_cor) or this_cor < max_cor:
                    continue
                if assignments[all_features[i]].umi_counts < (
                    assignments[all_features[j]].umi_counts
                ):
                    assignments[all_features[i]].is_contaminant = True
                    assignments[all_features[i]].correlated_features.append(all_features[j])
                else:
                    assignments[all_features[j]].is_contaminant = True
                    assignments[all_features[j]].correlated_features.append(all_features[i])
        return assignments

    def _infer_cell_fractions_by_tags(self) -> list[float]:
        """Get relative fraction of cells being stained by each tag."""
        if self.assignments is None:
            self.assignments = self.get_tag_assignments()
        cell_counts = []
        for _, asst in self.assignments.items():
            if asst.is_contaminant:
                continue
            cell_counts.append(len(asst.bc_indices))
        total_cell_counts = float(sum(cell_counts))  # doublets get count twice and so forth
        cell_fractions = [tk_stats.robust_divide(x, total_cell_counts) for x in cell_counts]
        return cell_fractions

    @staticmethod
    def get_poisson_multiplets_fractions(loading_rate, max_klet: int = NUM_TOTAL_TAGS):
        """Returns the fractions of GEMs for each n-let.

        Given a loading_rate and the total number of tags,
        returns fractions of GEMs that would be 1-let, 2-let ..., n-let
        """
        rate_v = sp_stats.poisson(loading_rate)
        k_lets = range(1, max_klet + 1)
        return rate_v.pmf(k_lets)

    @staticmethod
    def get_partitions_from_fractions(
        fractions: np.ndarray, num_gems: int = N_G, corr_factor: float = CORR_FACTOR
    ) -> np.ndarray:
        """Return the scaled fractions.

        Args:
            fractions (np.ndarray): TODO
            num_gems (int): TODO
            corr_factor (float): TODO

        Returns:
            np.ndarray:
        """
        # we only recover 1/corr_factor of those GEMs from the emulsion
        return fractions * (num_gems / corr_factor)

    def _helper_get_measurable_multiplets(self, n_let: int) -> np.ndarray:
        """A helper function to get a list of probability to detect a GEM with n_let cells.

        :return: np.array of float with a length of n_let
        """
        if self.assignments is None:
            self.assignments = self.get_tag_assignments()
        n_tags = sum(not asst.is_contaminant for asst in self.assignments.values())
        return _calculate_n_let_diversity_probs(n_let, n_tags)

    def get_measurable_multiplets_fractions(self, loading_rate) -> np.ndarray:
        """Predict the fractions of GEMs containing different number of cells.

        Uses given number of cells and tags.

        :return: np.array of non-zero float (variable length)
        """
        raw_poisson_fracs = TagAssigner.get_poisson_multiplets_fractions(loading_rate)

        measurable_fractions = np.array([0.0] * len(raw_poisson_fracs))
        for i, tfrac in enumerate(raw_poisson_fracs):
            if tfrac < MIN_POISSON_PROB:
                prob_measuring_this = 0.0
            else:
                prob_measuring_this = self._helper_get_measurable_multiplets(i + 1)
            measurable_fractions[: i + 1] += tfrac * prob_measuring_this
        return measurable_fractions
