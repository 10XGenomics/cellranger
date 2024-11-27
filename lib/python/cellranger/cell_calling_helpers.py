# Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
from __future__ import annotations

import enum
import itertools
from collections import Counter, OrderedDict
from collections.abc import Collection, Container, Iterable, Mapping
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, NamedTuple

import martian
import numpy as np
import pandas as pd
import scipy.stats as sp_stats
from scipy import interpolate

import cellranger.cell_calling as cr_cell
import cellranger.feature.antibody.analysis as ab_utils
import cellranger.rna.library as rna_library
from cellranger.cell_barcodes import CellBarcodes, InvalidCellBarcode
from cellranger.cell_calling import get_empty_drops_range
from cellranger.chemistry import CHEMISTRY_DESCRIPTION_FIELD, CHEMISTRY_SC3P_LT
from cellranger.constants import DEFAULT_RECOVERED_CELLS_PER_GEM_GROUP
from cellranger.metrics import BarcodeFilterResults
from cellranger.pandas_utils import sanitize_dataframe
from cellranger.rna.report_matrix import load_per_barcode_metrics
from tenkit.stats import robust_divide

if TYPE_CHECKING:
    from cellranger import matrix as cr_matrix

## Cell calling constants
ORDMAG_NUM_BOOTSTRAP_SAMPLES = 100
ORDMAG_RECOVERED_CELLS_QUANTILE = 0.99
NP_SORT_KIND = "stable"
N_CANDIDATE_BARCODES_GRADIENT = 20000
N_CANDIDATE_BARCODES_GRADIENT_LT = 2000
MIN_RECOVERED_CELLS_PER_GEM_GROUP = 50
MAX_RECOVERED_CELLS_PER_GEM_GROUP = 1 << 18


## High occupancy GEM filtering constants (multiplex FRP only)
RECOVERY_FACTOR = 1 / 1.65  # suggested by assaydev
TOTAL_INSTRUMENT_PARTITIONS = 115000  # assumes typical standard kit run


class MitoStats:  # pylint: disable=too-few-public-methods
    """Holds mitochondrial stats for diagnotic plots to be used in the websummary."""

    def __init__(self, matrix: cr_matrix.CountMatrix):
        """Class initialized with various stats.

        Args:
            matrix (cr_matrix.CountMatrix): Count matrix.
            mt_index (list): List of indices of mitochondrial genes in the count matrix.
        """
        self._get_mito_loc(matrix)
        self.total_umi = matrix.get_counts_per_bc()
        self.mt_umi = matrix.select_features(self.mt_index).get_counts_per_bc()
        self.mt_pct = (
            100.0 * self.mt_umi / self.total_umi
        )  # Percentage of UMIs per cell mapping to MT genes
        self.total_umi_log = np.log1p(self.total_umi) / np.log(10)
        # TO DO: add genome information for BY samples.

    def _get_mito_loc(self, matrix: cr_matrix.CountMatrix):
        """Function to get indices of mitochondrial genes."""
        mt_ensembl_ids = [
            # human
            "ENSG00000198888",
            "ENSG00000198763",
            "ENSG00000198804",
            "ENSG00000198712",
            "ENSG00000228253",
            "ENSG00000198899",
            "ENSG00000198938",
            "ENSG00000198840",
            "ENSG00000212907",
            "ENSG00000198886",
            "ENSG00000198786",
            "ENSG00000198695",
            "ENSG00000198727",
            # mouse
            "ENSMUSG00000064341",
            "ENSMUSG00000064345",
            "ENSMUSG00000064351",
            "ENSMUSG00000064354",
            "ENSMUSG00000064356",
            "ENSMUSG00000064357",
            "ENSMUSG00000064358",
            "ENSMUSG00000064360",
            "ENSMUSG00000065947",
            "ENSMUSG00000064363",
            "ENSMUSG00000064367",
            "ENSMUSG00000064368",
            "ENSMUSG00000064370",
        ]
        # The split on the underscore is to accommodate barnyard samples with prefixes on feature IDs
        ensembl_ids = [x.id.decode("UTF-8").split("_")[-1] for x in matrix.feature_ref.feature_defs]
        self.mt_index = [idx for idx, eid in enumerate(ensembl_ids) if eid in mt_ensembl_ids]


@dataclass
class HighOccupancyGemSummary:
    rtl_multiplexing_cells_per_gem_histogram: None | dict[str, int] = None
    rtl_multiplexing_estimated_lambda: float | None = 0.0
    rtl_multiplexing_total_probe_barcodes: int | None = 0
    rtl_multiplexing_high_occupancy_probe_barcode_count_threshold: int | None = 0
    rtl_multiplexing_fraction_cell_gems_high_occupancy: float | None = 0.0
    rtl_multiplexing_total_cells_in_high_occupancy_gems: int | None = 0
    rtl_multiplexing_fraction_cells_in_high_occupancy_gems: float | None = 0.0
    rtl_multiplexing_fraction_cell_reads_high_occupancy_gems: float | None = 0.0
    rtl_multiplexing_fraction_reads_high_occupancy_gems: float | None = 0.0
    rtl_multiplexing_fraction_cell_reads_in_cell_gems: float | None = 0.0


class FilterMethod(enum.Enum):
    # Caller-provided list of barcodes
    MANUAL = 0
    # Take the top N barcodes by count
    TOP_N_BARCODES = 1
    # Take barcodes within an order of magnitude of the max by count
    ORDMAG = 2
    # The above (ORDMAG), then find barcodes that differ from the ambient profile
    ORDMAG_NONAMBIENT = 3
    # The above (ORDMAG), then find barcodes above the steepest gradient in the log-log rank plot
    GRADIENT = 4
    # Apply (GRADIENT) to total counts, then keep barcodes with target gene counts above minimum
    TARGETED = 5


def get_filter_method_name(fm: FilterMethod) -> str:
    if fm == FilterMethod.MANUAL:
        return "manual"
    elif fm == FilterMethod.TOP_N_BARCODES:
        return "topn"
    elif fm == FilterMethod.ORDMAG:
        return "ordmag"
    elif fm == FilterMethod.ORDMAG_NONAMBIENT:
        return "ordmag_nonambient"
    elif fm == FilterMethod.GRADIENT:
        return "gradient"
    elif fm == FilterMethod.TARGETED:
        return "targeted"
    else:
        raise ValueError("Unsupported filter method value %d" % fm)


def get_filter_method_from_string(name: str) -> FilterMethod:
    if name == "manual":
        return FilterMethod.MANUAL
    elif name == "topn":
        return FilterMethod.TOP_N_BARCODES
    elif name == "ordmag":
        return FilterMethod.ORDMAG
    elif name == "ordmag_nonambient":
        return FilterMethod.ORDMAG_NONAMBIENT
    elif name == "gradient":
        return FilterMethod.GRADIENT
    elif name == "targeted":
        return FilterMethod.TARGETED
    else:
        raise ValueError("Unknown filter method value %d" % name)


def validate_cell_calling_args(
    recovered_cells, force_cells, cell_barcodes, method_name: str, feature_types
):
    def throw_err(msg):
        raise ValueError(msg)

    # this will throw a ValueError if it's not a recognized error
    method = get_filter_method_from_string(method_name)

    if method in (FilterMethod.ORDMAG, FilterMethod.ORDMAG_NONAMBIENT):
        pass

    elif method == FilterMethod.MANUAL:
        if cell_barcodes is None:
            throw_err(f"'cell_barcodes' must be specified when method is '{method_name}'")

    elif method == FilterMethod.TOP_N_BARCODES:
        if force_cells is None:
            throw_err(f"'force_cells' must be specified when method is '{method_name}'")


###############################################################################
def detect_antibody_antigen_aggreagtes(
    correction_data: pd.DataFrame, matrix, lib_type, metrics, *, num_probe_barcodes: int | None
):
    augmented_table = ab_utils.filter_correction_table(correction_data, lib_type)
    if lib_type == rna_library.ANTIBODY_LIBRARY_TYPE:
        # First, detect highly corrected barcodes
        highly_corrected_bcs = ab_utils.detect_highly_corrected_bcs(augmented_table)
        # Next, detect barcodes with high antibody UMI counts
        aggregate_bcs = ab_utils.detect_aggregate_barcodes(
            matrix, num_probe_barcodes=num_probe_barcodes
        )
        # Combine both sets and report a df of aggregate barcodes
        bcs_to_remove = list(set.union(highly_corrected_bcs, aggregate_bcs))
        assert len(bcs_to_remove) == len(set(bcs_to_remove))
    else:
        bcs_to_remove = ab_utils.detect_outlier_umis_bcs(matrix)
    removed_bcs_df = ab_utils.subselect_augmented_table(bcs_to_remove, augmented_table)

    ### report how many aggregates were found, and the fraction of reads those accounted for
    report_prefix = rna_library.get_library_type_metric_prefix(lib_type)
    metrics[report_prefix + "number_aggregate_GEMs"] = len(bcs_to_remove)
    frac_reads_removed = removed_bcs_df[ab_utils.FRACTION_TOTAL_READS].sum()
    metrics[report_prefix + "reads_lost_to_aggregate_GEMs"] = frac_reads_removed
    return bcs_to_remove, removed_bcs_df


def remove_antibody_antigen_aggregates(
    correction_data: pd.DataFrame,
    matrix,
    lib_types,
    disable_aggreagtion: bool,
    *,
    num_probe_barcodes: int | None,
):
    """Given a CountMatrix and and csv file containing information about umi corrected reads, detect.

    antibody protein aggreagtes:
        1) all barcodes with unusually high fraction of corrected reads, and
        2) all barcodes with unusually high antibody UMI counts
    antigen protein aggregates:
        1) all barcodes with unusually high antigen UMI counts
    remove both from the CoutMatrix
    """
    metrics_to_report = {}
    bcs_to_remove = set()
    removed_bcs_df = pd.DataFrame()

    if rna_library.ANTIBODY_LIBRARY_TYPE in lib_types:
        ab_to_remove, removed_ab_df = detect_antibody_antigen_aggreagtes(
            correction_data,
            matrix,
            rna_library.ANTIBODY_LIBRARY_TYPE,
            metrics_to_report,
            num_probe_barcodes=num_probe_barcodes,
        )
        bcs_to_remove = set.union(bcs_to_remove, ab_to_remove)
        removed_bcs_df = pd.concat([removed_bcs_df, removed_ab_df])
    if rna_library.ANTIGEN_LIBRARY_TYPE in lib_types:
        ag_to_remove, removed_ag_df = detect_antibody_antigen_aggreagtes(
            correction_data,
            matrix,
            rna_library.ANTIGEN_LIBRARY_TYPE,
            metrics_to_report,
            num_probe_barcodes=num_probe_barcodes,
        )
        bcs_to_remove = set.union(bcs_to_remove, ag_to_remove)
        removed_bcs_df = pd.concat([removed_bcs_df, removed_ag_df])
    removed_bcs_df.drop_duplicates(inplace=True)

    bcs_to_remove = {matrix.bc_to_int(bc) for bc in bcs_to_remove}
    # make sure filtered_bcs is in deterministic order or any later bootstrap sampling will not be deterministic
    filtered_bcs = [i for i in range(matrix.bcs_dim) if i not in bcs_to_remove]
    cleaned_matrix = matrix.select_barcodes(filtered_bcs)

    # Sort and sanitize
    removed_bcs_df.sort_index(inplace=True)
    sanitize_dataframe(removed_bcs_df, inplace=True)

    if disable_aggreagtion:
        # Run all the code above and generate metrics,
        # but do not actually remove detected aggregate barcodes
        return matrix, metrics_to_report, removed_bcs_df
    return cleaned_matrix, metrics_to_report, removed_bcs_df


def _get_high_occupancy_gem_threshold(
    estimated_lambda: float,
    probe_barcodes_observed: list[str],
    total_simulated_gems: int = 1000000,
) -> int:
    """Determine threshold for classifying GEMs as having more probe barcodes than expected.

    Args:
        estimated_lambda (float): estimated lambda for the run
        probe_barcodes_observed (Iterable[str]): full list of probe barocodes from each filtered cell barcode (one entry per filtered barcode).
            Used to get the total number of probe barcodes and account for the frequency of each probe barcode within filtered barcodes.
        total_simulated_gems (int): total number of simulated GEMs to generate at specified lambda
            note that can be arbitrary value, does not have to match GEMs/run.

    Returns:
        int: the number of probe barcodes per GEM over which a GEMs will be considered high occupancy
    """
    if estimated_lambda == 0:
        return 0

    np.random.seed(0)
    random_gems = np.random.poisson(estimated_lambda, total_simulated_gems)
    random_gems = list(random_gems[random_gems > 0])
    max_cells_in_gem = np.max(random_gems)

    barcode_probabilities = np.array(list(Counter(probe_barcodes_observed).values())) / len(
        probe_barcodes_observed
    )
    barcodes = list(range(len(barcode_probabilities)))

    # Much faster to sim random max number of cells for all and then subset
    simulated_barcodes = np.random.choice(
        barcodes, size=(len(random_gems), max_cells_in_gem), p=barcode_probabilities
    )
    simulated_barcode_counts = [
        len(set(simulated_barcodes[i, :barcode_count]))
        for i, barcode_count in enumerate(random_gems)
    ]

    return int(np.ceil(np.quantile(simulated_barcode_counts, 0.999)))


def remove_bcs_from_high_occupancy_gems(
    filtered_bcs: Collection[str],
    genome_filtered_bcs: dict[str, Iterable[str]],
    per_barcode_metrics: str,
    probe_bc_offset: int,
    total_instrument_partitions: int = TOTAL_INSTRUMENT_PARTITIONS,
    recovery_factor: float = RECOVERY_FACTOR,
) -> tuple[Iterable[str], dict[str, Iterable[str]], HighOccupancyGemSummary]:
    """Remove barcodes in GEMs that have more probe barcodes observed than would expect given loading.

    Args:
        filtered_bcs (Iterable[str]): cell calls prior to filtering
        genome_filtered_bcs (Dict[str, Iterable[str]]): genome_filtered_barcodes prior to filtering
        per_barcode_metrics (str): path to per_barcode_metrics file
        probe_bc_offset (int): starting position of probe barcode in cell barcode
        total_instrument_partitions (int, optional): Approximate number of partitions genrated by instrument. Defaults to total_instrument_partitions.
        recovery_factor (float, optional): Approximate recovery factor of generated partitions. Defaults to RECOVERY_FACTOR.

    Returns:
        Tuple[Iterable[str], Dict[str, Iterable[str]], HighOccupancyGemSummary]: tuple with
            input filtered_bcs with high occupancy GEMs removed),
            input genome_filtered_bcs with high occupancy GEMs removed),
            HighOccupancyGemSummary with filtering summary metrics
    """
    probe_bcs = [bc[probe_bc_offset:] for bc in filtered_bcs]
    total_probe_barcodes = len(set(probe_bcs))

    # Histogram of total filtered barcodes per GEM
    gem_bcs = [bc[:probe_bc_offset] for bc in filtered_bcs]
    gems_with_filtered_barcodes = len(set(gem_bcs))
    barcodes_per_gem = Counter(gem_bcs)
    gem_bc_count_distribution = Counter(barcodes_per_gem.values())

    gem_bc_count_distribution[0] = max(
        0, int((total_instrument_partitions * recovery_factor) - gems_with_filtered_barcodes)
    )  # don't let this go negative

    # Estimate lambda
    estimated_lambda = float(
        np.average(
            list(gem_bc_count_distribution.keys()), weights=list(gem_bc_count_distribution.values())
        )
    )

    # Compute high occupancy GEM threshold
    high_occupancy_gem_threshold = _get_high_occupancy_gem_threshold(estimated_lambda, probe_bcs)

    # Stats re: high occupancy gems
    barcodes_in_high_occupancy_gem = [
        barcode
        for barcode in filtered_bcs
        if barcodes_per_gem[barcode[0:probe_bc_offset]] > high_occupancy_gem_threshold
    ]
    total_high_occupancy_gems = sum(
        count > high_occupancy_gem_threshold for _, count in barcodes_per_gem.items()
    )
    total_barcodes_in_high_occupancy_gems = len(barcodes_in_high_occupancy_gem)

    barcodes_in_high_occupancy_gem_set = set(barcodes_in_high_occupancy_gem)

    # Final summary
    summary = HighOccupancyGemSummary()
    summary.rtl_multiplexing_cells_per_gem_histogram = dict(gem_bc_count_distribution)
    summary.rtl_multiplexing_estimated_lambda = estimated_lambda
    summary.rtl_multiplexing_total_probe_barcodes = total_probe_barcodes
    summary.rtl_multiplexing_high_occupancy_probe_barcode_count_threshold = (
        high_occupancy_gem_threshold
    )
    summary.rtl_multiplexing_fraction_cell_gems_high_occupancy = robust_divide(
        total_high_occupancy_gems, gems_with_filtered_barcodes
    )
    summary.rtl_multiplexing_total_cells_in_high_occupancy_gems = (
        total_barcodes_in_high_occupancy_gems
    )
    summary.rtl_multiplexing_fraction_cells_in_high_occupancy_gems = robust_divide(
        total_barcodes_in_high_occupancy_gems, len(filtered_bcs)
    )

    # Read metrics (use of original cell calls vs. version with high occupancy GEMs removed is intentional)
    bc_to_reads = load_per_barcode_metrics(per_barcode_metrics)
    if bc_to_reads:
        # Only compute read-wise metrics if we have GEX data.
        # These are only PD-facing; we can figure out how to add them if someone
        # actually requests them.
        total_reads = sum(bc_to_reads.values())  # includes NO_BARCODE entry already
        total_cell_reads = sum(bc_to_reads[bc] for bc in filtered_bcs)
        total_reads_cell_gems = sum(
            bc_to_reads[bc] for bc in bc_to_reads if bc[0:probe_bc_offset] in barcodes_per_gem
        )
        total_cell_reads_high_occupancy_gems = sum(
            bc_to_reads[bc] for bc in filtered_bcs if bc in barcodes_in_high_occupancy_gem_set
        )
        summary.rtl_multiplexing_fraction_cell_reads_high_occupancy_gems = robust_divide(
            total_cell_reads_high_occupancy_gems, total_cell_reads
        )
        summary.rtl_multiplexing_fraction_reads_high_occupancy_gems = robust_divide(
            total_cell_reads_high_occupancy_gems, total_reads
        )
        summary.rtl_multiplexing_fraction_cell_reads_in_cell_gems = robust_divide(
            total_cell_reads, total_reads_cell_gems
        )

    # Revised filtered barcodes
    filtered_bcs = [bc for bc in filtered_bcs if bc not in barcodes_in_high_occupancy_gem_set]
    for k in genome_filtered_bcs:
        genome_filtered_bcs[k] = [
            bc for bc in genome_filtered_bcs[k] if bc not in barcodes_in_high_occupancy_gem_set
        ]

    return filtered_bcs, genome_filtered_bcs, summary


################################################################################
class MetricGroups(NamedTuple):
    gem_group: int
    genome: str
    sample: str
    method: str


def call_initial_cells(
    matrix,
    genomes: Iterable[str],
    sample: str,
    unique_gem_groups: Collection[int],
    method: FilterMethod,
    recovered_cells: int | None,
    cell_barcodes: str | bytes,
    force_cells: int,
    feature_types: Container,
    chemistry_description: str | None,
    target_features: Iterable[int] | None = None,
    has_cmo_data: bool | None = False,
    *,
    num_probe_barcodes: int | None = None,
):
    """Call initial cells using the ordmag algorithm.

    Args:
        chemistry_description: Used to update the computation for LT Targeted
        num_probe_barcodes: Used to update ambient RNA estimation for RTL-multiplexing
    """
    # make sure we have a CountMatrixView
    matrix = matrix.view()

    # (gem_group, genome) => dict
    filtered_metrics_groups: dict[MetricGroups, BarcodeFilterResults] = {}
    # (gem_group, genome) => list of barcode strings
    filtered_bcs_groups = OrderedDict()

    if recovered_cells:
        gg_recovered_cells = recovered_cells // len(unique_gem_groups)
    else:
        gg_recovered_cells = None

    if force_cells:
        gg_force_cells = force_cells // len(unique_gem_groups)
    else:
        gg_force_cells = None
    for genome in genomes:
        # All sub-selection of the matrix should happen here & be driven by genome & feature_types
        genome_matrix = matrix.select_features_by_genome_and_types(genome, feature_types)

        # Make initial cell calls for each gem group individually
        for gem_group in unique_gem_groups:
            gg_matrix = genome_matrix.select_barcodes_by_gem_group(gem_group)

            gg_filtered_metrics, gg_filtered_bcs = _call_cells_by_gem_group(
                gg_matrix,
                method,
                gg_recovered_cells,
                gg_force_cells,
                cell_barcodes,
                chemistry_description,
                target_features,
                has_cmo_data,
                num_probe_barcodes,
            )
            filtered_metrics_groups[
                MetricGroups(gem_group, genome, sample, get_filter_method_name(method))
            ] = gg_filtered_metrics
            filtered_bcs_groups[(gem_group, genome)] = gg_filtered_bcs

    return filtered_metrics_groups, filtered_bcs_groups


def _call_cells_by_gem_group(
    gg_matrix,
    method: FilterMethod,
    gg_recovered_cells: int | None,
    gg_force_cells: int | None,
    cell_barcodes: str | bytes,
    chemistry_description,
    target_features: Iterable[int] | None,
    has_cmo_data: bool | None = False,
    num_probe_barcodes: int | None = None,
) -> tuple[BarcodeFilterResults, list[bytes]]:
    # counts per barcode
    gg_bc_counts = gg_matrix.get_counts_per_bc()

    if method in (FilterMethod.ORDMAG, FilterMethod.ORDMAG_NONAMBIENT):
        gg_filtered_indices, gg_filtered_metrics, msg = filter_cellular_barcodes_ordmag(
            gg_bc_counts,
            gg_recovered_cells,
            chemistry_description,
            num_probe_barcodes,
        )
        gg_filtered_bcs = gg_matrix.ints_to_bcs(gg_filtered_indices)

    elif method == FilterMethod.MANUAL:
        bcs_in_matrix: set[bytes] = set(gg_matrix.matrix.bcs)
        cell_barcodes = CellBarcodes(cell_barcodes)
        if has_cmo_data:
            cell_barcodes.verify_bcs_match_whitelist(chemistry_description)
        gg_filtered_bcs = []
        for bc in cell_barcodes:
            if bc in bcs_in_matrix:
                gg_filtered_bcs.append(bc)
            elif has_cmo_data:
                # If a Barcode Assignment CSV is used, all barcodes should appear in the
                # raw matrix or we consider this an error, as it will mess up downstream plotting.
                # I believe we only allow this to be possible for other data types to enable spatial
                # tissue detection to expect a barcode but not get data for it.
                raise InvalidCellBarcode(f"Barcode: {bc.decode()} had no observed reads.")
        gg_filtered_metrics = BarcodeFilterResults.init_with_constant_call(len(gg_filtered_bcs))
        msg = None

    elif method == FilterMethod.TOP_N_BARCODES:
        gg_filtered_indices, gg_filtered_metrics, msg = filter_cellular_barcodes_fixed_cutoff(
            gg_bc_counts, gg_force_cells
        )
        gg_filtered_bcs = gg_matrix.ints_to_bcs(gg_filtered_indices)

    elif method == FilterMethod.GRADIENT:
        gg_filtered_indices, gg_filtered_metrics, msg = filter_cellular_barcodes_gradient(
            gg_bc_counts, gg_recovered_cells
        )
        gg_filtered_bcs = gg_matrix.ints_to_bcs(gg_filtered_indices)

    elif method == FilterMethod.TARGETED:
        gg_bc_counts_targeted = gg_matrix.select_features(target_features).get_counts_per_bc()
        if chemistry_description == CHEMISTRY_SC3P_LT[CHEMISTRY_DESCRIPTION_FIELD]:
            max_num_additional_cells = N_CANDIDATE_BARCODES_GRADIENT_LT
        else:
            max_num_additional_cells = N_CANDIDATE_BARCODES_GRADIENT
        gg_filtered_indices, gg_filtered_metrics, msg = filter_cellular_barcodes_targeted(
            gg_bc_counts, gg_recovered_cells, gg_bc_counts_targeted, max_num_additional_cells
        )
        gg_filtered_bcs = gg_matrix.ints_to_bcs(gg_filtered_indices)

    else:
        martian.exit(f"Unsupported BC filtering method: {method}")
        raise SystemExit()

    if msg is not None:
        martian.log_info(msg)

    return gg_filtered_metrics, gg_filtered_bcs


def call_additional_cells(
    matrix,
    unique_gem_groups: list[int],
    genomes: list[str],
    filtered_bcs_groups,
    feature_types,
    chemistry_description,
    probe_barcode_sample_id,
    *,
    num_probe_barcodes,
    emptydrops_minimum_umis: dict[tuple[str, int], int] | int,
    num_sims: int = cr_cell.NUM_SIMS,
):
    """Call additional cells using the emptydrops algorithm.

    Args:
        chemistry_description: Used to update ambient RNA estimation for LT GEX/GEX+Ab
        num_sims: How many multinomial log likelihoods to generate
        num_probe_barcodes: Used to update ambient RNA estimation for RTL-multiplexing

    """
    # Track these for recordkeeping
    eval_bcs_arrays = []
    umis_per_bc_arrays = []
    loglk_arrays = []
    pvalue_arrays = []
    pvalue_adj_arrays = []
    nonambient_arrays = []
    genome_call_arrays = []
    sample_id_arrays = []

    emptydrops_threshold = dict()
    matrix = matrix.select_features_by_types(feature_types)

    for gg in unique_gem_groups:
        # Do it by gem group and genome
        gg_matrix = matrix.select_barcodes_by_gem_group(gg)
        # Get all initial cell calls for this gem group.
        gg_bcs = set()
        for group, bcs in filtered_bcs_groups.items():
            if group[0] == gg:
                gg_bcs.update(bcs)
        gg_bcs = list(sorted(gg_bcs))

        # All sub-selection of the matrix should happen here & be driven by genome & feature_types
        for genome in genomes:
            genome_matrix = gg_matrix.select_features_by_genome(genome)

            if isinstance(emptydrops_minimum_umis, int):
                ed_minimum_umis = emptydrops_minimum_umis
            else:
                ed_minimum_umis = emptydrops_minimum_umis[(genome, gg)]

            result = cr_cell.find_nonambient_barcodes(
                genome_matrix,
                gg_bcs,
                chemistry_description,
                num_probe_barcodes,
                emptydrops_minimum_umis=ed_minimum_umis,
                num_sims=num_sims,
            )
            if result is None:
                print(f"Failed at attempt to call non-ambient barcodes in GEM well {gg}")
                continue

            umis_per_bc = gg_matrix.get_counts_per_bc()

            emptydrops_threshold[(genome, gg)] = result.emptydrops_minimum_umis
            eval_bcs = np.array([x.decode() for x in gg_matrix.bcs[result.eval_bcs]])
            eval_bcs_arrays.append(eval_bcs)
            umis_per_bc_arrays.append(umis_per_bc[result.eval_bcs])
            loglk_arrays.append(result.log_likelihood)
            pvalue_arrays.append(result.pvalues)
            pvalue_adj_arrays.append(result.pvalues_adj)
            nonambient_arrays.append(result.is_nonambient)
            genome_call_arrays.append([genome] * len(result.eval_bcs))  # allow doublet calling
            sample_id_arrays.append([probe_barcode_sample_id] * len(result.eval_bcs))

            # Update the lists of cell-associated barcodes
            eval_bc_strs = np.array(gg_matrix.bcs)[result.eval_bcs]
            filtered_bcs_groups[(gg, genome)].extend(eval_bc_strs[result.is_nonambient])

    if len(eval_bcs_arrays) > 0:
        nonambient_summary = pd.DataFrame(
            OrderedDict(
                [
                    ("barcode", np.concatenate(eval_bcs_arrays)),
                    ("umis", np.concatenate(umis_per_bc_arrays)),
                    ("ambient_loglk", np.concatenate(loglk_arrays)),
                    ("pvalue", np.concatenate(pvalue_arrays)),
                    ("pvalue_adj", np.concatenate(pvalue_adj_arrays)),
                    ("nonambient", np.concatenate(nonambient_arrays)),
                    ("genome", np.concatenate(genome_call_arrays)),
                    ("sample", np.concatenate(sample_id_arrays)),
                ]
            )
        )
    else:
        nonambient_summary = pd.DataFrame()

    return filtered_bcs_groups, nonambient_summary, emptydrops_threshold


def apply_mitochondrial_threshold(
    matrix: cr_matrix.CountMatrix,
    genomes: list[str],
    unique_gem_groups: list[int],
    filtered_bcs_groups: OrderedDict,
    max_mito_percent: float,
) -> tuple[OrderedDict, pd.DataFrame]:
    """Subsets to barcodes that do not exceed set threshold of mitochondrial read percentage.

    Args:
        matrix (cr_matrix.CountMatrix): feature barcode matrix
        genomes (list[str]): list of genomes
        unique_gem_groups (list[int]): list of gem groups
        filtered_bcs_groups (OrderedDict): list of filtered barcodes
        max_mito_percent (float): mitochondrial read percentage threshold

    Returns:
        tuple[OrderedDict, pd.DataFrame]: updated list of filtered barcodes, summary information dataframe
    """
    matrix = matrix.select_features_by_type(rna_library.GENE_EXPRESSION_LIBRARY_TYPE)

    bcs_removed = []
    total_umis = []
    mt_pct = []
    genome_list = []
    gem_group_list = []

    # Do it by gem group and genome
    for genome in genomes:
        genome_matrix = matrix.select_features_by_genome(genome)
        for gem_group in unique_gem_groups:
            gem_group_matrix = genome_matrix.select_barcodes_by_gem_group(gem_group)
            # Subset matrix to barcodes that have already been called as cells
            gem_group_matrix = gem_group_matrix.select_barcodes_by_seq(
                filtered_bcs_groups[(gem_group, genome)]
            )
            mito_stats = MitoStats(gem_group_matrix)
            filtered_barcode_index = [
                idx for idx, mt in enumerate(mito_stats.mt_pct) if (mt > max_mito_percent)
            ]

            n_filtered_barcodes = len(filtered_barcode_index)

            if n_filtered_barcodes > 0:
                # Update the lists of cell-associated barcodes
                filtered_bcs_groups[(gem_group, genome)] = np.delete(
                    arr=np.array(gem_group_matrix.bcs), obj=np.array(filtered_barcode_index)
                )
                bcs_removed.append(np.array(gem_group_matrix.bcs)[np.array(filtered_barcode_index)])
                total_umis.append(np.array(mito_stats.total_umi)[np.array(filtered_barcode_index)])
                mt_pct.append(np.array(mito_stats.mt_pct)[np.array(filtered_barcode_index)])
                genome_list.append(np.repeat(genome, n_filtered_barcodes))
                gem_group_list.append(np.repeat(gem_group, n_filtered_barcodes))

    all_bcs_removed = list(itertools.chain.from_iterable(bcs_removed))
    if len(all_bcs_removed) > 0:
        mitochondrial_summary = pd.DataFrame(
            {
                "Barcodes removed": all_bcs_removed,
                "Total UMIs": list(itertools.chain.from_iterable(total_umis)),
                "MT percent": list(itertools.chain.from_iterable(mt_pct)),
                "threshold": max_mito_percent,
                "genome": list(itertools.chain.from_iterable(genome_list)),
                "gem_group": list(itertools.chain.from_iterable(gem_group_list)),
            }
        )
    else:
        mitochondrial_summary = pd.DataFrame(
            {
                "Barcodes removed": [-1],
                "Total UMIs": [-1],
                "MT percent": [-1],
                "threshold": max_mito_percent,
            }
        )
    return filtered_bcs_groups, mitochondrial_summary


def apply_global_minimum_umis_threshold(
    matrix,
    unique_gem_groups: list[int],
    genomes: list[str],
    filtered_bcs_groups,
    feature_types,
    global_minimum_umis: dict[tuple[str, int], int] | int,
):
    """Apply a global minimum UMI threshold on cell calls."""
    matrix = matrix.select_features_by_types(feature_types)

    # Do it by gem group and genome
    for genome in genomes:
        genome_matrix = matrix.select_features_by_genome(genome)
        for gg in unique_gem_groups:
            gg_matrix = genome_matrix.select_barcodes_by_gem_group(gg)
            orig_cell_bcs = set(filtered_bcs_groups[(gg, genome)])
            orig_cell_idx = np.flatnonzero(
                np.fromiter(
                    (bc in orig_cell_bcs for bc in gg_matrix.bcs),
                    count=len(gg_matrix.bcs),
                    dtype=bool,
                )
            )

            if isinstance(global_minimum_umis, int):
                minimum_umis = global_minimum_umis
            else:
                minimum_umis = global_minimum_umis[(genome, gg)]

            umis_per_cell_bc = gg_matrix.get_counts_per_bc()[orig_cell_idx]
            final_cell_idx = orig_cell_idx[umis_per_cell_bc >= minimum_umis]

            # Update the lists of cell-associated barcodes
            filtered_bcs_groups[(gg, genome)] = np.array(gg_matrix.bcs)[final_cell_idx]

    return filtered_bcs_groups


################################################################################
def merge_filtered_metrics(
    filtered_metrics: Mapping[MetricGroups, BarcodeFilterResults]
) -> dict[str, float | int]:
    """Merge all the barcode filter results and return them as a dictionary."""
    total_filtered_bcs = 0
    total_filtered_bcs_var = 0
    dresult: dict[str, float | int] = {}
    for key_tuple, fm in filtered_metrics.items():
        dresult.update(
            fm.to_dict_with_prefix(key_tuple.gem_group, key_tuple.sample, key_tuple.method)
        )
        # Compute metrics over all gem groups
        total_filtered_bcs += fm.filtered_bcs
        total_filtered_bcs_var += fm.filtered_bcs_var

    # Estimate CV based on sum of variances and means
    # This metric only applies to the initial cell calls
    filtered_bcs_cv = robust_divide(np.sqrt(total_filtered_bcs_var), total_filtered_bcs)

    dresult.update({"filtered_bcs_cv": filtered_bcs_cv})
    return dresult


def combine_initial_metrics(
    genomes: Iterable[str],
    filtered_metrics_groups: Mapping[MetricGroups, BarcodeFilterResults],
    genome_filtered_bcs,
    summary,
):
    # Combine initial-cell-calling metrics
    for genome in genomes:
        # Merge metrics over all gem groups and samples for this genome
        txome_metrics = {k: v for k, v in filtered_metrics_groups.items() if k.genome == genome}
        txome_summary = merge_filtered_metrics(txome_metrics)

        prefix = genome + "_" if genome else ""
        summary.update({(f"{prefix}{key}"): summary for key, summary in txome_summary.items()})

        summary[f"{prefix}filtered_bcs"] = len(genome_filtered_bcs.get(genome, {}))

    return summary


def summarize_bootstrapped_top_n(top_n_boot, nonzero_counts):
    top_n_bcs_mean = np.mean(top_n_boot)
    top_n_bcs_var = np.var(top_n_boot)
    top_n_bcs_sd = np.sqrt(top_n_bcs_var)
    result = BarcodeFilterResults()
    result.filtered_bcs_var = top_n_bcs_var
    result.filtered_bcs_cv = robust_divide(top_n_bcs_sd, top_n_bcs_mean)
    result.filtered_bcs_lb = np.round(sp_stats.norm.ppf(0.025, top_n_bcs_mean, top_n_bcs_sd), 0)
    result.filtered_bcs_ub = np.round(sp_stats.norm.ppf(0.975, top_n_bcs_mean, top_n_bcs_sd), 0)

    nbcs = int(np.round(top_n_bcs_mean))
    result.filtered_bcs = nbcs

    # make sure that if a barcode with count x is selected, we select all barcodes with count >= x
    # this is true for each bootstrap sample, but is not true when we take the mean

    if nbcs > 0:
        order = np.argsort(nonzero_counts, kind=NP_SORT_KIND)[::-1]
        sorted_counts = nonzero_counts[order]

        cutoff = sorted_counts[nbcs - 1]
        index = nbcs - 1
        while (index + 1) < len(sorted_counts) and sorted_counts[index] == cutoff:
            index += 1
            # if we end up grabbing too many barcodes, revert to initial estimate
            if (index + 1 - nbcs) > 0.20 * nbcs:
                return result
            result.filtered_bcs = index + 1
            result.filtered_bcs_cutoff = cutoff
    return result


def find_within_ordmag(x, baseline_idx):
    x_ascending = np.sort(x)
    # Add +1 as we're getting from the other side
    baseline = x_ascending[-(baseline_idx + 1)]
    cutoff = np.maximum(1, np.round(0.1 * baseline)).astype(int)
    # Return the index corresponding to the cutoff in descending order
    return len(x) - np.searchsorted(x_ascending, cutoff)


def estimate_recovered_cells_ordmag(nonzero_bc_counts, max_expected_cells: int):
    """Estimate the number of recovered cells by trying to find ordmag(recovered) =~ filtered."""
    # - Search for a result such that some loss(recovered_cells, filtered_cells) is minimized.
    # - Here I'm using (obs - exp)**2 / exp, which approximates a proportion for small differences
    #   but blows up for large differences.
    # - Test over a log2-spaced range of values from 1..262_144
    recovered_cells = np.linspace(1, np.log2(max_expected_cells), 2000)
    recovered_cells = np.unique(np.round(np.power(2, recovered_cells)).astype(int))
    baseline_bc_idx = np.round(recovered_cells * (1 - ORDMAG_RECOVERED_CELLS_QUANTILE))
    baseline_bc_idx = np.minimum(baseline_bc_idx.astype(int), len(nonzero_bc_counts) - 1)
    filtered_cells = find_within_ordmag(nonzero_bc_counts, baseline_bc_idx)
    loss = np.power(filtered_cells - recovered_cells, 2) / recovered_cells
    idx = np.argmin(loss)
    return recovered_cells[idx], loss[idx]


def filter_cellular_barcodes_ordmag(
    bc_counts: np.ndarray[int, np.dtype[np.int_]],
    recovered_cells: int | None,
    chemistry_description: str | None = None,
    num_probe_barcodes: int | None = None,
):
    """All barcodes that are close to within an order of magnitude of a top barcode.

    Takes all barcodes that are close to within an order of magnitude of a
    top barcode that likely represents a cell.
    """
    rs = np.random.RandomState(0)

    metrics = BarcodeFilterResults(0)

    nonzero_bc_counts = bc_counts[bc_counts > 0]
    if len(nonzero_bc_counts) == 0:
        msg = "WARNING: All barcodes do not have enough reads for ordmag, allowing no bcs through"
        return [], metrics, msg

    if recovered_cells is None:
        # Set the most cells to examine based on the empty drops range for this chemistry
        max_expected_cells = (
            min(
                get_empty_drops_range(chemistry_description, num_probe_barcodes)[0],
                MAX_RECOVERED_CELLS_PER_GEM_GROUP,
            )
            if (chemistry_description is not None)
            else MAX_RECOVERED_CELLS_PER_GEM_GROUP
        )
        recovered_cells, loss = np.mean(
            np.array(
                [
                    estimate_recovered_cells_ordmag(
                        rs.choice(nonzero_bc_counts, len(nonzero_bc_counts)), max_expected_cells
                    )
                    for _ in range(ORDMAG_NUM_BOOTSTRAP_SAMPLES)
                ]
            ),
            axis=0,
        )
        recovered_cells = max(int(np.round(recovered_cells)), MIN_RECOVERED_CELLS_PER_GEM_GROUP)
        print(f"Found recovered_cells = {recovered_cells} with loss = {loss}")
    else:
        recovered_cells = max(recovered_cells, MIN_RECOVERED_CELLS_PER_GEM_GROUP)
        print(f"Using provided recovered_cells = {recovered_cells}")

    baseline_bc_idx = int(np.round(float(recovered_cells) * (1 - ORDMAG_RECOVERED_CELLS_QUANTILE)))
    baseline_bc_idx = min(baseline_bc_idx, len(nonzero_bc_counts) - 1)

    # Bootstrap sampling; run algo with many random samples of the data
    top_n_boot = np.array(
        [
            find_within_ordmag(
                rs.choice(nonzero_bc_counts, len(nonzero_bc_counts)), baseline_bc_idx
            )
            for _ in range(ORDMAG_NUM_BOOTSTRAP_SAMPLES)
        ]
    )

    metrics.update(summarize_bootstrapped_top_n(top_n_boot, nonzero_bc_counts))

    # Get the filtered barcodes
    top_n = metrics.filtered_bcs
    top_bc_idx = np.sort(np.argsort(bc_counts, kind=NP_SORT_KIND)[::-1][0:top_n])
    assert top_n <= len(nonzero_bc_counts), "Invalid selection of 0-count barcodes!"
    return top_bc_idx, metrics, None


def filter_cellular_barcodes_fixed_cutoff(bc_counts, cutoff: int):
    nonzero_bcs = len(bc_counts[bc_counts > 0])
    top_n = min(cutoff, nonzero_bcs)
    top_bc_idx = np.sort(np.argsort(bc_counts, kind=NP_SORT_KIND)[::-1][0:top_n])
    metrics = BarcodeFilterResults.init_with_constant_call(top_n)
    metrics.filtered_bcs_cutoff = np.sort(bc_counts)[::-1][top_n - 1]
    return top_bc_idx, metrics, None


def filter_cellular_barcodes_targeted(
    bc_counts,
    recovered_cells,
    bc_counts_targeted,
    max_num_additional_cells=N_CANDIDATE_BARCODES_GRADIENT,
    min_umis_additional_cells=cr_cell.TARGETED_CC_MIN_UMIS_ADDITIONAL_CELLS,
    min_targeted_umis=cr_cell.TARGETED_CC_MIN_UMIS_FROM_TARGET_GENES,
):
    """First, call cells using the gradient method (filter_cellular_barcodes_gradient) on the total bc_counts from all genes.

    Then retain all cell-associated barcodes which pass the min_targeted_umis threshold of UMI counts from target genes.
    """
    grad_cells_idx, _, msg = filter_cellular_barcodes_gradient(
        bc_counts, recovered_cells, max_num_additional_cells, min_umis_additional_cells
    )
    if msg is not None:
        martian.log_info(msg)

    targ_cells_idx = [idx for idx in grad_cells_idx if bc_counts_targeted[idx] >= min_targeted_umis]
    num_targ_cells = len(targ_cells_idx)
    metrics = BarcodeFilterResults.init_with_constant_call(num_targ_cells)
    metrics.filtered_bcs_cutoff = min_targeted_umis
    return targ_cells_idx, metrics, None


def filter_cellular_barcodes_gradient(
    bc_counts,
    recovered_cells,
    max_num_additional_cells=N_CANDIDATE_BARCODES_GRADIENT,
    min_umis_additional_cells=cr_cell.TARGETED_CC_MIN_UMIS_ADDITIONAL_CELLS,
    infer_throughput=False,
):
    """Take all barcodes with counts above the value of the steepest gradient.

    From the log-transformed barcode rank plot [log(counts) vs log(rank) sorted by descending counts].
    This minimum gradient is computed over the allowed x-range based on ordmag cutoff and
    the allowable number of additional cells to consider. Prior to computing this value, the
    barcode rank plot is fit to a smoothing spline curve, and the interpolated first derivative of this
    curve is used to identify the point associated with maximum descent.
    """
    if recovered_cells is None:
        recovered_cells = DEFAULT_RECOVERED_CELLS_PER_GEM_GROUP
    else:
        recovered_cells = max(recovered_cells, MIN_RECOVERED_CELLS_PER_GEM_GROUP)

    metrics = BarcodeFilterResults(0)

    nonzero_bc_counts = bc_counts[bc_counts > 0]
    nonzero_bc_counts = np.array(sorted(nonzero_bc_counts)[::-1])  # sort in descending order
    if len(nonzero_bc_counts) == 0:
        msg = "WARNING: All barcodes do not have enough reads for gradient, allowing no bcs through"
        return np.array([]), metrics, msg

    baseline_bc_idx = int(np.round(float(recovered_cells) * (1 - ORDMAG_RECOVERED_CELLS_QUANTILE)))
    baseline_bc_idx = min(baseline_bc_idx, len(nonzero_bc_counts) - 1)
    baseline_count_threshold = nonzero_bc_counts[baseline_bc_idx]

    if infer_throughput:
        lower_bc_idx = 0
        max_num_additional_cells = 150000  # need to search bc idx 70-220k
        min_umis_additional_cells = 3
    else:
        lower_bc_idx = np.sum(nonzero_bc_counts >= baseline_count_threshold / 10.0) - 1
        lower_bc_idx = min(lower_bc_idx, len(nonzero_bc_counts) - 1)

    upper_bc_idx = min(
        lower_bc_idx + max_num_additional_cells,
        np.sum(nonzero_bc_counts >= min_umis_additional_cells),
    )
    upper_bc_idx = max(upper_bc_idx, lower_bc_idx)
    upper_bc_idx = min(upper_bc_idx, len(nonzero_bc_counts) - 1)

    uniq_counts = sorted(set(nonzero_bc_counts))[
        ::-1
    ]  # collapse unique values and sort in descending order
    log_y_values = [np.log10(a) for a in uniq_counts]
    x_values = [np.sum(nonzero_bc_counts >= a) for a in uniq_counts]
    log_x_values = [np.log10(x) for x in x_values]
    log_x_values.append(
        np.log10(1 + sum(nonzero_bc_counts))
    )  # append end-point values to handle small input arrays
    log_y_values.append(0.0)

    # Fit spline to barcode rank plot curve, then interpolate its first derivative over log_x_values
    spline_degree = min(3, len(log_y_values) - 1)  # adjust degree parameter for small input arrays
    fit_spline = interpolate.UnivariateSpline(
        x=log_x_values, y=log_y_values, k=spline_degree, s=0, check_finite=True
    )
    # For large input arrays, reduce the number of knots for progressive smoothing
    if len(log_x_values) > 50:
        num_spline_knots = get_spline_num_knots(len(log_x_values))
        orig_knots = fit_spline.get_knots()
        if num_spline_knots < len(orig_knots):
            knot_values = [
                orig_knots[i]
                for i in np.linspace(1, len(orig_knots) - 2, num_spline_knots - 2, dtype=int)
            ]
            fit_spline = interpolate.LSQUnivariateSpline(
                x=log_x_values, y=log_y_values, t=knot_values, k=spline_degree, check_finite=True
            )

    # Find minimum gradient value for barcodes within the allowed x-range
    gradients = fit_spline(log_x_values[0:-1], 1)
    gradients = np.where(
        [x >= lower_bc_idx and x <= upper_bc_idx for x in x_values],
        gradients,
        np.repeat(0, len(gradients)),
    )
    gradient_count_cutoff = np.round(10 ** log_y_values[np.argmin(gradients)], 0)
    gradient_num_cells = max(np.sum(nonzero_bc_counts > gradient_count_cutoff), lower_bc_idx + 1)

    # Get filtered barcodes corresponding to gradient-based cutoff
    top_n = min(gradient_num_cells, len(nonzero_bc_counts))
    top_bc_idx = np.sort(np.argsort(bc_counts, kind=NP_SORT_KIND)[::-1][0:top_n])
    metrics = BarcodeFilterResults.init_with_constant_call(top_n)
    return top_bc_idx, metrics, None


def get_spline_num_knots(n):
    """Heuristic function to estimate the number of knots to be used for spline interpolation.

    as a function of the number of unique input data points n.
    """
    if n < 50:
        return int(n)
    a1 = np.log2(50)
    a2 = np.log2(100)
    a3 = np.log2(140)
    a4 = np.log2(200)
    if n < 200:
        return int(2 ** (a1 + (a2 - a1) * (n - 50) / 150))
    if n < 800:
        return int(2 ** (a2 + (a3 - a2) * (n - 200) / 600))
    if n < 3200:
        return int(2 ** (a3 + (a4 - a3) * (n - 800) / 2400))
    return int(200 + (n - 3200) ** (0.2))


################################################################################


class CellCallingParam:
    """Python equivalent of CellCallingParam struct."""

    per_gem_well: float | None = None
    per_sample: dict[str, float] | None = None

    def __init__(self, data):
        """Initialize with data provided by martian."""
        if data is not None:
            self.per_gem_well = data["per_gem_well"]
            self.per_sample = data["per_sample"]

        # We check that only one of these is set
        gw_none = self.per_gem_well is None
        ps_none = self.per_sample is None or all(v is None for v in self.per_sample.values())
        assert gw_none or ps_none


class CellCalling:
    """Python equivalent of cr_lib::parse_multi_config::CellCalling struct."""

    recovered_cells: CellCallingParam | None = CellCallingParam(None)
    force_cells: CellCallingParam | None = CellCallingParam(None)
    emptydrops_minimum_umis: CellCallingParam | None = CellCallingParam(None)
    global_minimum_umis: CellCallingParam | None = CellCallingParam(None)
    max_mito_percent: CellCallingParam | None = CellCallingParam(None)
    cell_barcodes: str | None = None  # A JSON file path
    override_mode: str | None = None
    override_library_types: list[str] | None = None
    disable_ab_aggregate_detection: bool = False
    disable_high_occupancy_gem_detection: bool = False

    def __init__(self, args: dict[str, Any] | None):
        """Converts from an args dictionary.

        Takes the untyped arguments passed by martian into the stage code and
        convert them into a typed object.
        """
        if args is not None:
            self.cell_barcodes = args["cell_barcodes"]
            self.disable_ab_aggregate_detection = args["disable_ab_aggregate_detection"]
            self.disable_high_occupancy_gem_detection = args["disable_high_occupancy_gem_detection"]
            self.recovered_cells = CellCallingParam(args["recovered_cells"])
            self.force_cells = CellCallingParam(args["force_cells"])
            self.emptydrops_minimum_umis = CellCallingParam(args["emptydrops_minimum_umis"])
            self.max_mito_percent = CellCallingParam(args["max_mito_percent"])
            self.global_minimum_umis = CellCallingParam(args["global_minimum_umis"])
            self.override_mode = args["override_mode"]
            self.override_library_types = args["override_library_types"]


def get_recovered_cells(recovered_cells: CellCallingParam, sample=None) -> int | None:
    """Extracts a value for recovered_cells from the CellCallingParam struct.

    a. A single library-level value
    b. An aggregate library-level value for multiplex FRP samples
    c. The value for a specific multiplex FRP sample
    """
    if recovered_cells.per_gem_well is not None:
        return int(recovered_cells.per_gem_well)
    elif recovered_cells.per_sample is not None and sample is None:
        acc = 0
        for val in recovered_cells.per_sample.values():
            if val is None:
                return None
            acc += val
        return int(acc)
    elif recovered_cells.per_sample is not None and sample is not None:
        return (
            int(recovered_cells.per_sample.get(sample))
            if recovered_cells.per_sample.get(sample)
            else None
        )
    else:
        return None


def get_force_cells(force_cells: CellCallingParam, sample=None) -> int | None:
    """Extracts a value for force_cells from the CellCallingParam struct."""
    if force_cells.per_gem_well is not None:
        return int(force_cells.per_gem_well)
    elif force_cells.per_sample is not None and sample is not None:
        return (
            int(force_cells.per_sample.get(sample)) if force_cells.per_sample.get(sample) else None
        )
    else:
        return None


def get_global_minimum_umis(global_minimum_umis: CellCallingParam, sample=None) -> int:
    """Extracts a value for global_minimum_umis from the CellCallingParam struct."""
    if global_minimum_umis.per_gem_well is not None:
        return int(global_minimum_umis.per_gem_well)
    elif global_minimum_umis.per_sample is not None and sample is not None:
        return int(global_minimum_umis.per_sample.get(sample) or cr_cell.MIN_GLOBAL_UMIS)
    else:
        return cr_cell.MIN_GLOBAL_UMIS


def get_max_mito_percent(
    max_mito_percent: CellCallingParam | None, sample: str | None = None
) -> float:
    """Extracts a value for maximum mitochondrial percentage from the CellCallingParam struct."""
    if max_mito_percent is not None:
        if max_mito_percent.per_gem_well is not None:
            return max_mito_percent.per_gem_well
        elif max_mito_percent.per_sample is not None and sample is not None:
            return max_mito_percent.per_sample.get(sample) or cr_cell.MAX_MITO_PCT
        else:
            return cr_cell.MAX_MITO_PCT
    else:
        return cr_cell.MAX_MITO_PCT


def get_emptydrops_minimum_umis(emptydrops_minimum_umis: CellCallingParam, sample=None) -> int:
    """Extracts a value for emptydrops_minimum_umis from the CellCallingParam struct."""
    if emptydrops_minimum_umis.per_gem_well is not None:
        return int(emptydrops_minimum_umis.per_gem_well)
    elif emptydrops_minimum_umis.per_sample is not None and sample is not None:
        return int(emptydrops_minimum_umis.per_sample.get(sample) or cr_cell.MIN_UMIS)
    else:
        return cr_cell.MIN_UMIS


def remove_cells_with_zero_targeted_counts(matrix, filtered_bcs, genome_filtered_bcs, summary):
    target_features = matrix.feature_ref.get_target_feature_indices()
    bc_counts_targeted = matrix.select_features(target_features).get_counts_per_bc()
    bcs_zero_idx = np.argwhere(bc_counts_targeted == 0).flatten()
    bcs_zero = set(matrix.ints_to_bcs(bcs_zero_idx))

    # Revised filtered barcodes
    revised_filtered_bcs = [bc for bc in filtered_bcs if bc not in bcs_zero]
    for k in genome_filtered_bcs:
        genome_filtered_bcs[k] = [bc for bc in genome_filtered_bcs[k] if bc not in bcs_zero]

    summary.update(
        {"cell_bcs_removed_with_zero_targeted_umis": len(filtered_bcs) - len(revised_filtered_bcs)}
    )
    return revised_filtered_bcs, genome_filtered_bcs
