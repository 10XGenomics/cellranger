"""A TagAssigner class based on the JIBES model."""

# Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

from __future__ import annotations

import json
import resource
import time
from enum import Enum

import numpy as np
import pandas as pd
from six import ensure_binary, ensure_str

import cellranger.analysis.jibes as jibes
import cellranger.analysis.jibes_data as jibes_data
import cellranger.feature.feature_assigner as cr_fa
import cellranger.rna.library as rna_library
import tenkit.stats as tk_stats
from cellranger.analysis.jibes_constants import (
    ASSIGNMENT_COL_NAME,
    ASSIGNMENT_PROB_COL_NAME,
    BLANK_FACTOR_NAME,
    JIBES_MIN_CONFIDENCE,
    MULTIPLETS_FACTOR_NAME,
    UNASSIGNED_FACTOR_NAME,
)
from cellranger.analysis.jibes_data import BARCODE_COL
from cellranger.feature.throughputs import G19_N_GEMS, N_G

MODEL_RESULTS_HEADER = [
    "ID",
    "iterations",
    "LL",
    "num_obs",
    "est_cells",
    "singlets",
    "Over0.9",
    "Over0.95",
    "Over0.99",
    "latent_states",
    "elapsed",
    "max_rss",
    "singlet_capture_ratio",
]


class JibesError(Exception):
    """Exception class to indicate an error in the algorithm.

    that should lead to an exit instead of a stack trace.
    """


def make_model_results_list(name, samp_id, assigner):
    """Make a list of list with results in it.

    :param name: Sample Name
    :param samp_id: Sample Id
    :param assigner: instance of JibesTagAssigner
    :return:
    """
    isinstance(assigner, JibesTagAssigner)
    results = [MODEL_RESULTS_HEADER]
    res = assigner.jibes_assignments
    singlets = (res[ASSIGNMENT_COL_NAME] != MULTIPLETS_FACTOR_NAME) & (
        res[ASSIGNMENT_COL_NAME] != BLANK_FACTOR_NAME
    )
    total_singlets = np.sum(singlets)
    o90 = np.sum(res[ASSIGNMENT_PROB_COL_NAME][singlets] > 0.90)
    o95 = np.sum(res[ASSIGNMENT_PROB_COL_NAME][singlets] > 0.95)
    o99 = np.sum(res[ASSIGNMENT_PROB_COL_NAME][singlets] > 0.99)
    max_rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    fitter = assigner.fitter
    metrics = assigner.get_feature_assignment_metrics()
    report_prefix = assigner.report_prefix
    results += [
        [
            name,
            samp_id,
            fitter.iterations,
            fitter.LL,
            fitter.data.counts.shape[0],
            fitter.estimated_cells,
            total_singlets,
            o90,
            o95,
            o99,
            fitter.z,
            assigner.fit_time_s,
            max_rss,
            metrics[f"{report_prefix}sc_rec_efficiency"],
        ]
    ]
    return results


PARAM_RESULTS_HEADER = ["samp_name", "ID", "feature", "background", "foreground", "variance"]


def make_parameter_results_list(jibes_fit, sample_id, samp_name):
    """Returns a list of list giving the estimated parameters for each tag.

    :param jibes_fit: instance of JibesEM
    :param sample_id:
    :param samp_name:
    :return: A list of lists with
    """
    assert isinstance(jibes_fit, jibes.JibesEM)
    cols = jibes_fit.data.column_names
    assert len(cols) == len(jibes_fit.model.background)
    return [PARAM_RESULTS_HEADER] + [
        [samp_name, sample_id, name, b, f, v]
        for (name, b, f, v) in zip(
            cols, jibes_fit.model.background, jibes_fit.model.foreground, jibes_fit.model.std_devs
        )
    ]


def make_parameter_table_rows(jibes_fit: jibes.JibesEM, output_json) -> None:
    """Get a jibes parameter table.

    Generates the table in a form that Rust can read as rows of the table.

    Args:
        jibes_fit: instance of JibesEM
        output_json: Filename for json with a list of lists with results.
    """
    assert isinstance(jibes_fit, jibes.JibesEM)
    cols = jibes_fit.data.column_names
    assert len(cols) == len(jibes_fit.model.background)
    res = []
    for i, name in enumerate(cols):
        # pylint: disable=invalid-name
        b = jibes_fit.model.background[i]
        f = jibes_fit.model.foreground[i]
        v = jibes_fit.model.std_devs[i]
        res.append(
            {
                "cmo_name": ensure_str(name),
                "background": b,
                "foreground": f,
                "std_dev": v,
            }
        )

    with open(output_json, "w") as outf:
        json.dump(res, outf)


class FitStatus(Enum):
    """Specifies the status of a fit."""

    # pylint:disable=invalid-name
    NotStarted = 0
    Finished = 1
    FailedToConverge = 2


class JibesTagAssigner(cr_fa.TagAssigner):
    """A class to do Tag assignments using the JIBES algorithm instead of a GMM.

    Currently inherits from TagAssigner to enable better code sharing, although the other
    FeatureAssigner classes are all marginal callers designed to represent values above
    background, not directly infer the multiplet and blanks state as the JIBES caller does.
    """

    def __init__(
        self,
        matrix,
        tag_calls_fn,
        library_type=rna_library.MULTIPLEXING_LIBRARY_TYPE,
        n_gems=N_G,
        confidence=JIBES_MIN_CONFIDENCE,
        exclude_blanks_from_posterior=True,
    ):
        """A Tag Assigner class based on the JIBES model.

        :param matrix: a CountMatrix object
        :param tag_calls_fn: File name to a tag calls file produced by the marginal tag caller
        :param confidence: Minimum confidence required to assign a cell
        :param exclude_blanks_from_posterior: When true, don't consider the blank state when determining singlet/multiplet
        confidence levels for thresholding.
        """
        super().__init__(matrix, library_type, n_gems)
        self.exclude_blanks_from_posterior = exclude_blanks_from_posterior

        self.valid_tags: set[bytes] = jibes.get_valid_tags(tag_calls_fn)
        if len(self.valid_tags) == 0:
            err_msg = """No cell multiplexing tag sequences were detected in the
Multiplexing Capture library. Common causes include:
1. Wrong pattern or sequences provided in the feature reference (CMO reference) csv file.
2. Corrupt or low quality reads.
3. Incorrect input fastq files for the Multiplexing Capture library.
Contact support@10xgenomics.com for additional help with this error."""
            raise JibesError(err_msg)
        assert isinstance(next(iter(self.valid_tags)), bytes)

        # TODO: Remove dependence on upstream tag_calls file
        dense_matrix = self.sub_matrix.m.toarray().astype(np.float64).T
        # With customer data we've occasionally seen barcodes with near 0 counts,
        # we plan to exclude these from fitting and just call them blank with prob = 1
        row_sums = np.sum(dense_matrix, axis=1)
        row_sums = np.ravel(row_sums)
        near_zero_indices = np.array(row_sums <= 12)
        for i in range(dense_matrix.shape[1]):
            # Note: The reverse transform is applied inside
            dense_matrix[:, i] = np.log10(dense_matrix[:, i] + 1.0)

        self.all_tags: list[bytes] = [x.id for x in self.sub_matrix.feature_ref.feature_defs]
        assert isinstance(self.all_tags[0], bytes)
        self.dropped_tags = set(self.all_tags) - self.valid_tags
        good_indices = [(i, x) for i, x in enumerate(self.all_tags) if x not in self.dropped_tags]
        assert good_indices, (self.all_tags, self.valid_tags)
        dense_matrix = dense_matrix[:, [x[0] for x in good_indices]]
        col_names = [x[1] for x in good_indices]

        data = jibes_data.JibesData(
            dense_matrix[~near_zero_indices, :], col_names, self.sub_matrix.bcs[~near_zero_indices]
        )
        self.near_zero_data = jibes_data.JibesData(
            dense_matrix[near_zero_indices, :], col_names, self.sub_matrix.bcs[near_zero_indices]
        )

        self.fitter = self._load_preliminary_assignments(tag_calls_fn, data, n_gems)

        self.fit_status = FitStatus.NotStarted
        self.fit_time_s = np.nan
        self.jibes_assignments = None
        self.non_singlet_barcodes = None
        self.confidence = confidence

    def _load_preliminary_assignments(self, tag_calls_fn: str, data, n_gems):
        prelim_assignments = pd.read_csv(
            tag_calls_fn,
            converters={
                cr_fa.FEATURE_CALL_COL: ensure_binary,
                cr_fa.CELL_BARCODE_COL: ensure_binary,
            },
            usecols=[cr_fa.FEATURE_CALL_COL, cr_fa.CELL_BARCODE_COL],
        )
        # Now remove near-zero barcodes from initial calls made upstream
        # keeping these two stages in sync is a bit of a pain.
        all_bcs = prelim_assignments[cr_fa.CELL_BARCODE_COL].astype("S")
        mask = np.fromiter((not (bc in self.near_zero_data.barcodes) for bc in all_bcs), dtype=bool)
        prelim_assignments = prelim_assignments[mask]
        init_data = data.subset_to_barcodes(all_bcs[mask])
        del mask
        if len(self.valid_tags) > 0 or data.counts.shape[0] == 0:
            starting_model = jibes.initial_params_from_tag_calls_df(
                prelim_assignments, init_data, n_gems=n_gems
            )
            return jibes.JibesEM(data, starting_model, n_gems=n_gems)
        else:
            return None

    # pylint: disable=too-many-locals
    def get_tag_assignments(self):
        """Get Tag assignments for this dataset after running fitting.

        There is some nuance in how these values are computed for the JIBES caller as oppposed
        to the older marginal ones.  In particular, this function normally lists all the barcodes
        that are "above" background for a given tag.  In the JIBES caller, instead of above/below background for
        each tag, we estimate which latent state (either a singlet or multiplet) a barcode might
        be in and assign a probability value to every such state.

        To transform the joint inference back into the marginal inference, we do the following
        operations:

        - If a barcode is "Blanks" it is not assigned to any FeatureAssignment
        - If a barcode is "Unassigned", it is not assigned to any FeatureAssignment due to lack of confidence
        - If a barcode is a Singlet, it is assigned to to the FeatureAssignment for that tag.
        - If a barcode is a multiplet, we find the most likely multiplet state and assign the barcode to all tags > 0 for that state.

        Note that in JIBES a barcode is assigned as a multiplet if the summed probability over
        all multiplet latent states is greater than any singlet state, so the `most likely multiplet`
        state, might not be the most likely state.

        :return: A dictionary of {tag_name:[FeatureAssignmentObjects]
        """
        if self.assignments is not None or self.fitter is None:
            return self.assignments
        start = time.time()
        self.fitter.perform_EM(abs_tol=0.1)
        self.fit_time_s = time.time() - start
        if not self.fitter.converged:
            self.fit_status = FitStatus.FailedToConverge
        self.jibes_assignments = jibes.get_assignment_df(
            self.fitter,
            confidence=self.confidence,
            exclude_blanks_from_posterior=self.exclude_blanks_from_posterior,
            blank_data_to_append=self.near_zero_data,
        )

        col_names = self.fitter.data.column_names
        name_to_index = {name: i for i, name in enumerate(col_names)}
        indices = [[] for _ in range(len(col_names))]
        multiplet_latent_state_cols = jibes.get_cols_associated_with_assignments(self.fitter)[
            MULTIPLETS_FACTOR_NAME
        ]
        # Posterior is dimenstions n_bcs x n_latent_states
        multiplet_posterior = self.fitter.posterior[:, multiplet_latent_state_cols]
        # X is dimensions n_latent_states x n_tags
        mutiplet_latent_states = self.fitter.X[multiplet_latent_state_cols, 1:]

        def get_most_likely_multiplet_state(row_index):
            max_prob_index = np.argmax(multiplet_posterior[row_index, :])
            latent_states = mutiplet_latent_states[max_prob_index, :]
            return latent_states

        non_singlet_barcodes = {
            BLANK_FACTOR_NAME: [],
            UNASSIGNED_FACTOR_NAME: [],
            MULTIPLETS_FACTOR_NAME: [],
        }

        for i in range(self.fitter.n):
            asn = self.jibes_assignments[ASSIGNMENT_COL_NAME][i]
            # the jibes_assignment dataframe is not guaranteed to have the same order
            # as when the data was loaded, so we want to check that it matches still
            bc = self.jibes_assignments[BARCODE_COL][i]
            assert bc == self.fitter.data.barcodes[i], "Dataframe is out of order!"
            bc_index = self.sub_matrix.bc_to_int(self.jibes_assignments[BARCODE_COL][i])
            if asn in non_singlet_barcodes:
                non_singlet_barcodes[asn].append(bc_index)
                if asn == MULTIPLETS_FACTOR_NAME:
                    most_probable_multiple_state = get_most_likely_multiplet_state(i)
                    for k, count in enumerate(most_probable_multiple_state):
                        if count > 0.0:
                            indices[k].append(bc_index)
            else:  # Singlet
                col = name_to_index[asn]
                indices[col].append(bc_index)

        blanks = non_singlet_barcodes[BLANK_FACTOR_NAME]
        for bc in self.near_zero_data.barcodes:
            bc_index = self.sub_matrix.bc_to_int(bc)
            blanks.append(bc_index)

        self.non_singlet_barcodes = non_singlet_barcodes

        assignments = {}
        for j, c_name in enumerate(col_names):
            umi_counts = self.sub_matrix.get_subselected_counts(
                log_transform=False, list_feature_ids=[c_name]
            )
            assignments[c_name] = cr_fa.FeatureAssignments(
                indices[j],
                umi_counts.sum(),
                False,
                [],
            )

        self.assignments = assignments
        return assignments

    def compute_assignment_metadata(self):
        """Compute assignment metadata for Tag library."""
        assignment_metadata = self.compute_generic_assignment_metadata()
        freqs_df = self.compute_assignment_freqs(assignment_metadata.num_cells_without_features)
        assignment_metadata.freq_counts = freqs_df
        umi_thresholds = self._compute_umi_thresholds(
            self.sub_matrix, self.assignments, multiplexing=True
        )
        assignment_metadata.umi_thresholds = umi_thresholds
        assignment_metadata.num_cells_blanks = len(self.non_singlet_barcodes[BLANK_FACTOR_NAME])
        assignment_metadata.frac_cells_blanks = tk_stats.robust_divide(
            assignment_metadata.num_cells_blanks, assignment_metadata.num_cells
        )
        assignment_metadata.num_cells_unassigned = len(
            self.non_singlet_barcodes[UNASSIGNED_FACTOR_NAME]
        )
        assignment_metadata.frac_cells_unassigned = tk_stats.robust_divide(
            assignment_metadata.num_cells_unassigned, assignment_metadata.num_cells
        )

        return assignment_metadata

    def _compute_no_tag_metric(self, report_prefix="", depth_suffix=""):
        return {
            f"{report_prefix}frac_no_tag_assigned{depth_suffix}": self.assignment_metadata.frac_cells_without_features,
            f"{report_prefix}observed_num_no_tag_assigned{depth_suffix}": self.assignment_metadata.num_cells_without_features,
            f"{report_prefix}frac_blanks{depth_suffix}": self.assignment_metadata.frac_cells_blanks,
            f"{report_prefix}observed_num_blanks{depth_suffix}": self.assignment_metadata.num_cells_blanks,
            f"{report_prefix}frac_unassigned{depth_suffix}": self.assignment_metadata.frac_cells_unassigned,
            f"{report_prefix}observed_num_unassigned{depth_suffix}": self.assignment_metadata.num_cells_unassigned,
        }

    def prepare_for_pickle(self):
        """Until pickle is no longer used, will convert fitter object.

        Fitter object is either a lightweight rust class or a heavy weight
        python class.  Either will be converted into a smaller data structure
        to pass along.

        :return:
        """
        self.fitter = jibes.FitterPickle(self.fitter)


def run_assignment_stage(
    filtered_matrix,
    marginal_tag_calls_per_cell,
    library_type=rna_library.MULTIPLEXING_LIBRARY_TYPE,
    throughput="MT",
    confidence=None,
    exclude_blanks_from_posterior=False,
):
    """Runs the assignment stage.

    :return:
    """
    if confidence is None:
        confidence = JIBES_MIN_CONFIDENCE
    n_gems = G19_N_GEMS.get(throughput, N_G)
    assigner = JibesTagAssigner(
        filtered_matrix,
        marginal_tag_calls_per_cell,
        library_type,
        n_gems,
        confidence,
        exclude_blanks_from_posterior=exclude_blanks_from_posterior,
    )
    assigner.get_tag_assignments()
    return assigner
