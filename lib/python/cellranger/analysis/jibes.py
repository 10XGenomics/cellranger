# Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
"""Code for simulating and fitting data from a Joint Inference By Exploiting Stoichiometry model.

The model assumes that for a given number of labeled cells going through a GEM Well, we have a
poisson based expectation around the number of k-lets (1,2,3,etc.) we might expect, and that in
log space the UMI (or read) counts for multiplets are a linear combination of the singlet states.

An EM algorithm is implemented to jointly infer the background, foreground, and variance of
counts across the entire dataset that is given.
"""
# allow zi, ps, X_full, one_EM_step and perform_EM names
# pylint: disable=invalid-name
from __future__ import annotations

import numpy as np
import pandas as pd
from six import ensure_binary

import cellranger.analysis.multigenome as cr_mga
import cellranger.feature.feature_assigner as cr_fa
import cellranger.feature.throughputs
import cellranger.pandas_utils as pu
import cellranger.rna.library as rna_library
from cellranger.analysis.jibes_constants import (
    ASSIGNMENT_COL_NAME,
    ASSIGNMENT_PROB_COL_NAME,
    BLANK_FACTOR_NAME,
    JIBES_MIN_CONFIDENCE,
    MULTIPLETS_FACTOR_NAME,
    UNASSIGNED_FACTOR_NAME,
)
from cellranger.analysis.jibes_data import JibesData
from cellranger.analysis.jibes_py import JibesModelPy
from cellranger.matrix import CountMatrix
from cellranger.molecule_counter import FEATURE_IDX_COL_NAME, MoleculeCounter

USE_PYO3 = True
if USE_PYO3:
    # pylint: disable=no-name-in-module,import-error
    from cellranger.analysis.jibes_o3 import JibesEMO3 as JibesEM
    from cellranger.analysis.jibes_o3 import JibesModelO3 as JibesModel
else:
    JibesModel = JibesModelPy
    from cellranger.analysis.jibes_py import JibesEMPy as JibesEM


_BARCODE = pu.FEATURE_DF_BARCODE_SEQ_COL
_COUNT = "count"


def get_valid_tags(tag_calls_per_cell_fn):
    """Loads the list of.

    :param tag_calls_per_cell_fn: Filename of tag_calls_per_cell file
    :return: set of valid CMO tags
    """
    cur_calls = pd.read_csv(tag_calls_per_cell_fn, usecols=[cr_fa.FEATURE_CALL_COL])
    valid = cur_calls[cr_fa.FEATURE_CALL_COL].astype("category")
    okay = valid.cat.categories
    all_tags = []
    for assignment in okay:
        assignment = ensure_binary(assignment)
        all_tags += assignment.split(cr_fa.FEATURE_SEPARATOR)
    return set(all_tags)


def load_tag_counts_from_molecule_info(mol_info_file):
    """Load all the read counts for the multiplexing library type from a molecule info file.

    Each cell-associated barcode and return as a data frame following log10 transformation
    along with a description of the features loaded.
    :param mol_info_file: A path to a molecule info file.
    :return: A tuple with a Pandas DataFrame and a Dictionary with Feature info.
    """
    with MoleculeCounter.open(mol_info_file, "r") as mc:
        feature_indices = mc.feature_reference.get_indices_for_type(
            rna_library.MULTIPLEXING_LIBRARY_TYPE
        )
        df = pu.mol_info_from_h5(
            mc,
            exclude_noncells=True,
            with_umi=False,
            with_barcode_seq=True,
            with_gem_group=False,
            with_library_idx=False,
            filter_feature_idx=feature_indices,
        )
    feature_map = {}
    for i in feature_indices:
        feature_map[i] = mc.feature_reference.feature_defs[i]
    grouped = df.groupby([_BARCODE, FEATURE_IDX_COL_NAME])
    dfs = pd.DataFrame({_COUNT: grouped[pu.FEATURE_DF_COUNT_COL].sum()})
    dfs = dfs.reset_index()
    df = dfs.pivot(index=_BARCODE, columns=FEATURE_IDX_COL_NAME, values=_COUNT)
    df = df.reset_index()
    cnames = list(df.columns)[1:]
    df.fillna({x: 0 for x in cnames}, inplace=True)
    for col in cnames:
        df[col] = np.log10(df[col] + 1)
    return (df, feature_map)


def create_initial_parameters(
    data: JibesData, assignments, n_gems: int = cellranger.feature.throughputs.N_G
):
    """Create initial conditions based on initial assignments.

    (or irrelevant) create initial conditions for the jibes model to fit

    Args:
        data:  a JibesData instance
        assignments: an ndarray ndarray where each element is a string matching the
                     column name of data
        n_gems: [TODO]

    Returns:
        JibesModel: initial conditions
    """
    # pylint: disable=too-many-locals
    k = len(data.column_names)
    foreground_means = np.zeros(k)
    std_devs = np.zeros(k)
    backgrounds = np.zeros(k)
    bad_indices = set()
    singletons = np.isin(assignments, data.column_names)
    for i, name in enumerate(data.column_names):
        vals = assignments == name
        if np.sum(vals) < 2:
            foreground_means[i] = np.nan
            backgrounds[i] = np.nan
            std_devs[i] = np.nan
            bad_indices.add(i)
        else:
            # TODO: Assuming only one library in the data
            tmp = assignments != name
            other_singletons = tmp & singletons
            if np.sum(other_singletons) > 0:
                backgrounds[i] = data.counts[other_singletons, i].mean()
            else:
                backgrounds[i] = data.counts[:, i].mean()
            foreground_counts = data.counts[vals, i].squeeze()
            bg = backgrounds[i]
            foreground_means[i] = np.maximum(0.6 + bg, np.mean(foreground_counts)) - bg
            std_devs[i] = foreground_counts.std()
    # If the marginal caller goes goofy, sometimes we have no singlet calls for a tag,
    # and we'll fill in with the average of other values
    if len(bad_indices) > 0:
        print(
            f"A total of {len(bad_indices)} tags did not have enough assigned samples during initialization."
        )
        bad_indices = list(bad_indices)
        if len(bad_indices) == len(data.column_names):
            # TODO: Need Better Initial Conditions
            foreground_means[bad_indices] = 1.0
            backgrounds[bad_indices] = 0.5
            std_devs[bad_indices] = 0.3
            # raise RuntimeError("None of the Tags had singleton calls in the earlier data.")
        else:
            good = list(x for x in range(len(data.column_names)) if x not in bad_indices)
            fm = np.mean(foreground_means[good])
            bm = np.mean(backgrounds[good])
            vm = np.mean(std_devs[good])
            foreground_means[bad_indices] = fm
            backgrounds[bad_indices] = bm
            std_devs[bad_indices] = vm
    # TODO: Improve this edge case handling
    zero_std_devs = np.where(std_devs == 0.0)[0]
    if len(zero_std_devs) > 0:
        print("Edge case with no variance in initial calls detected!!")
        std_devs[zero_std_devs] = 0.2
    model = JibesModel(backgrounds, foreground_means, std_devs, n_gems=n_gems)
    return model


def _order_initial_parameters(assignments, data):
    """Helper function to make sure the barcodes in two dataframes used for initialization.

    are in the same order.

    :param assignments:
    :param data:
    :return:
    """
    bc_to_index = {x: i for i, x in enumerate(assignments[cr_fa.CELL_BARCODE_COL])}
    indices = np.fromiter(
        (bc_to_index[x] for x in data.barcodes), dtype=int, count=len(data.barcodes)
    )
    return assignments.take(indices).reset_index(drop=True)


def initial_params_from_tag_calls_df(assignments, data, n_gems=cellranger.feature.throughputs.N_G):
    """Make initial parameter guesses to fit a JIBES model.

    Given a DataFrame loaded from tag_calls_per_cell_fn created by CALL_TAGS_MARGINAL,
    make initial parameter guesses to fit a JIBES model with by utilizing the
    earlier cell calls.

    Args:
        assignments (pd.DataFrame): produced from loading a tag_calls_per_cell file
        data (JibesData): A JibesData object
        n_gems (int): [TODO]

    Returns:
        JibesModel: initial conditions
    """
    # pylint: disable=too-many-locals
    assert isinstance(data, JibesData)
    assert assignments.shape[0] == len(data.barcodes)
    assignments = _order_initial_parameters(assignments, data)
    return create_initial_parameters(
        data, assignments[cr_fa.FEATURE_CALL_COL].values, n_gems=n_gems
    )


def get_cols_associated_with_assignments(fit):
    """From a fit model, determine which cols are associated with each state.

    Blank, singlet (or doublets of same type) and multiplet latent states.
    :param fit:
    :return: a dictionary that goes from assignment (numeric for singlets) ->row_indices
    """
    isinstance(fit, JibesEM)
    X = fit.X[:, 1:]
    row_sums = np.sum(X, axis=1)
    all_indices = set(range(1, X.shape[0]))
    matches = {}
    for k in range(fit.k):
        all_tag_rows = [i for i in range(1, X.shape[0]) if X[i, k] == row_sums[i]]
        matches[k] = all_tag_rows
        for x in all_tag_rows:
            all_indices.discard(x)
    # Now the zero category, silly cut/paste, it's the first row
    all_tag_rows = [i for i in range(X.shape[0]) if row_sums[i] == 0]
    matches[BLANK_FACTOR_NAME] = all_tag_rows
    for x in all_tag_rows:
        all_indices.discard(x)
    matches[MULTIPLETS_FACTOR_NAME] = list(all_indices)
    return matches


def format_count_column_name(name: bytes | str) -> str:
    """Used for accessing the counts in an assignment dataframe.

    Args:
        name:

    Returns:
        f"{name}_cnts
    """
    if isinstance(name, bytes):
        name = name.decode()
    return f"{name}_cnts"


def get_assignment_df(
    fit: JibesEM,
    is_jibes_multi_genome: bool = False,
    confidence: float = JIBES_MIN_CONFIDENCE,
    exclude_blanks_from_posterior: bool = False,
    blank_data_to_append: JibesData = None,
) -> pd.DataFrame:
    """Get the assignments as a dataframe.

    Returns a Pandas dataframe where each row is a barcode, and each column is either a feature ID, one of
    [Barcode, Multiplet, Blank, Assignment, Assignment_Probability].
    Each value is the posterior probability for that barcode to be in that state.
    Note that the df is not indexed by barcodes, that info is in the barcodes column.
    Also, Assignment is either a feature ID (singlet) or Blank (no CMO state) or
    Unassigned (no single state is likely enough to exceed the minimum confidence required) or
    Multiplet (exceeds posterior probability for >1 CMO state).
    Assignment_Probability contains the probability of the most likely state,
    regardless of whether or not it meets the minimum threshold).

    Returns:
      pd.DataFrame: each row is a barcode, and each column is either a feature ID, one of
            [Barcode, Multiplet, Blank, Assignment, Assignment_Probability].
    """
    assert isinstance(fit, JibesEM)
    df = fit.data.get_df_from_data()
    matches = get_cols_associated_with_assignments(fit)
    all_names = []
    posterior = fit.posterior  # only get this once as it's a copy operation from the Rust side
    cnts = pd.DataFrame()
    for i, col in enumerate(fit.data.column_names):
        indices = matches[i]
        name = col
        all_names.append(name)
        probs = posterior[:, indices].sum(axis=1)
        cnts[format_count_column_name(name)] = df[name]
        df[name] = probs

    indices = matches[MULTIPLETS_FACTOR_NAME]
    all_names.append(MULTIPLETS_FACTOR_NAME)
    probs = posterior[:, indices].sum(axis=1)
    df[MULTIPLETS_FACTOR_NAME] = probs

    indices = matches[BLANK_FACTOR_NAME]
    all_names.append(BLANK_FACTOR_NAME)
    probs = posterior[:, indices].sum(axis=1)
    df[BLANK_FACTOR_NAME] = probs

    df[ASSIGNMENT_COL_NAME] = df[all_names].idxmax(axis=1)
    df[ASSIGNMENT_PROB_COL_NAME] = df[all_names].max(axis=1)
    df = pd.concat([df, cnts], axis=1)
    if is_jibes_multi_genome:
        return df

    if blank_data_to_append is not None and len(blank_data_to_append.barcodes) > 0:
        blank_cnts = pd.DataFrame()
        blank_df = blank_data_to_append.get_df_from_data()
        for name in fit.data.column_names:
            blank_cnts[format_count_column_name(name)] = blank_df[name]
            blank_df[name] = np.repeat(0.0, blank_df.shape[0])
        blank_df[MULTIPLETS_FACTOR_NAME] = np.repeat(0.0, blank_df.shape[0])
        blank_df[BLANK_FACTOR_NAME] = np.repeat(1.0, blank_df.shape[0])

        blank_df[ASSIGNMENT_COL_NAME] = np.repeat(BLANK_FACTOR_NAME, blank_df.shape[0])
        blank_df[ASSIGNMENT_PROB_COL_NAME] = np.repeat(1.0, blank_df.shape[0])
        blank_df = pd.concat([blank_df, blank_cnts], axis=1)
        df = pd.concat([df, blank_df])
        df = df.reset_index(drop=True)

    df[ASSIGNMENT_COL_NAME] = df.apply(
        lambda row: _enforce_min_confidence(
            row[ASSIGNMENT_COL_NAME],
            row[ASSIGNMENT_PROB_COL_NAME],
            row[BLANK_FACTOR_NAME],
            exclude_blanks_from_posterior,
            confidence=confidence,
        ),
        axis=1,
    )
    df.index.name = "Row"
    return df


def _enforce_min_confidence(
    candidate_asst,
    asst_probability,
    blank_probability,
    exclude_blanks_from_posterior,
    confidence=JIBES_MIN_CONFIDENCE,
    unassigned_name=UNASSIGNED_FACTOR_NAME,
):
    if candidate_asst != BLANK_FACTOR_NAME and exclude_blanks_from_posterior:
        asst_probability = asst_probability / (1.0 - blank_probability)
    if asst_probability < confidence:
        return unassigned_name
    return candidate_asst


class JIBESMultiGenome(cr_mga.MultiGenomeAnalysis):
    """A class to enable the JIBES model to classify bar."""

    def __init__(
        self,
        filtered_matrix: CountMatrix | None = None,
        n_gems: int = cellranger.feature.throughputs.N_G,
    ):
        """Initialize the model.

        Args:
            filtered_matrix:
        """
        super().__init__(filtered_matrix)
        if self.filtered_matrix is not None:
            assert isinstance(filtered_matrix, CountMatrix)
        self.suffix = "_jibes"
        self.n_gems = n_gems

    def _infer_multiplets(self, counts0, counts1, bootstraps=None, **kwargs):
        """See base class for doc_string.

        Args:
            counts0:
            counts1:
            bootstraps: Ignored

        Returns:
            int: observed multiplets
            np.ndarray:
            np.ndarray:
        """
        assert len(counts0) == len(counts1)
        n_gems = kwargs.get("n_gems", cellranger.feature.throughputs.N_G)
        barcodes = np.arange(len(counts0))
        cnt0 = np.log10(counts0 + 1)
        cnt1 = np.log10(counts1 + 1)
        data = JibesData(
            np.vstack((cnt0, cnt1)).T,
            [
                cr_mga.analysis_constants.GEM_CLASS_GENOME0,
                cr_mga.analysis_constants.GEM_CLASS_GENOME1,
            ],
            barcodes,
        )
        # Initialize the model
        assignments = np.array(cr_mga.classify_gems(counts0, counts1))
        # return as observed when only one or none genome was observed
        if len(set(assignments)) <= 1:
            return 0, np.array(0), assignments
        init = create_initial_parameters(data, assignments, n_gems=self.n_gems)
        fitter = JibesEM(data, init, n_gems=n_gems)
        fitter.perform_EM()
        df = get_assignment_df(fitter, is_jibes_multi_genome=True)
        obs_mults = sum(df[ASSIGNMENT_COL_NAME] == MULTIPLETS_FACTOR_NAME)
        c1 = sum(df[ASSIGNMENT_COL_NAME] == cr_mga.analysis_constants.GEM_CLASS_GENOME0)
        c2 = sum(df[ASSIGNMENT_COL_NAME] == cr_mga.analysis_constants.GEM_CLASS_GENOME1)
        # TODO: The multigenome code currently does not have a "Blank" state, though perhaps
        # it should be included, not doing so now to stay consistent though
        # n = float(sum(df['Assignment'] != "Blanks"))

        return (
            obs_mults,
            np.array(cr_mga.infer_multiplets_from_observed(obs_mults, c1, c2)),
            df[ASSIGNMENT_COL_NAME].values,
        )


# pylint: disable=too-few-public-methods
class FitterPickle:
    """Holds a reduced set of JibesEM object data for downstream usage."""

    def __init__(self, current_instance):
        """Create a new instance.

        :param current_instance: Either a JibesEMO3 or JibesEMPy object
        """
        self.data = current_instance.data
        # Since the counts in the model no longer include the blanks, delete this so they can't be used again
        # inappropriately
        self.data.counts = None
        self.LL = current_instance.LL
        self.n = current_instance.n
        self.z = current_instance.z
        self.iterations = current_instance.iterations
        cur_model = current_instance.model
        if USE_PYO3:
            cur_model = JibesModelPy(
                cur_model.background,
                cur_model.foreground,
                cur_model.std_devs,
                cur_model.frequencies,
                cur_model.blank_freq,
                cur_model.n_gems,
            )
        self.model = cur_model
