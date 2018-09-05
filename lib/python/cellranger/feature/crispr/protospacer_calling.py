import itertools
import numpy as np
import pandas as pd
pd.set_option("compute.use_numexpr", False)
import tenkit.stats as tk_stats
import sys
import cellranger.stats as cr_stats

UMI_NUM_TRIES = 10 # Number of initial points to try for GMM-fitting
UMI_MIX_INIT_SD = 0.25 # Intial standard deviation for GMM components

ESSENTIAL_CRISPR_FILES = ["cells_per_protospacer.json", "protospacer_calls_summary.csv", "protospacer_umi_thresholds_json.json"]
CRISPR_ANALYSIS_FILE_NAMES =  ["protospacer_calls_summary.csv", "protospacer_calls_per_cell.csv", "cells_per_protospacer.json",
                                    "protospacer_umi_thresholds_csv.csv", "protospacer_umi_thresholds_json.json",
                                    "perturbation_efficiencies_by_feature.csv", "perturbation_efficiencies_by_target.csv"]

def get_protospacer_call_metrics(ps_calls_summary, num_gex_cbs, report_prefix):
    metrics_dict = {}
    num_cells_with_multiple_protospacers = ps_calls_summary.loc['> 1 protospacer expressed', '# Cells']
    num_cells_with_protospacer =  (ps_calls_summary.loc['1 protospacer expressed', '# Cells'] +
                                        num_cells_with_multiple_protospacers)

    frac_cells_with_protospacer = tk_stats.robust_divide(num_cells_with_protospacer, num_gex_cbs)
    frac_cells_with_multiple_protospacer = tk_stats.robust_divide(num_cells_with_multiple_protospacers, num_gex_cbs)

    metrics_dict.update({
                report_prefix + 'frac_cells_with_protospacer': frac_cells_with_protospacer,
                report_prefix + 'frac_cells_with_multiple_protospacer': frac_cells_with_multiple_protospacer,
                })

    return metrics_dict

def _call_presence_with_gmm(umi_counts):
    log_umi = np.log(1. + umi_counts)

    umi_try_means = [[0, x] for x in np.linspace(0, np.max(log_umi), 1+UMI_NUM_TRIES)[1:]]
    umi_gmm = cr_stats.multistart_gmm(data = cr_stats.to_col_vec(log_umi),
                                                     weights=[0.5, 0.5],
                                                     means_list=umi_try_means,
                                                     sd=UMI_MIX_INIT_SD)

    if not umi_gmm.converged_:
        sys.stdout.write("Warning: EM did not converge for UMIs!")

    umi_posterior = umi_gmm.predict_proba(cr_stats.to_col_vec(log_umi))
    high_umi_component = np.argmax(umi_gmm.means_)
    in_high_umi_component = umi_posterior[:,high_umi_component] > 0.5

    return in_high_umi_component

def get_ps_calls_and_summary(filtered_guide_counts_matrix, f_map):
    """Calculates protospacer calls per cell and summarizes them
    Args:
        filtered_guide_counts_matrix: CountMatrix - obtained by selecting features by CRISPR library type on the feature counts matrix
        f_map: dict - map of feature ID:feature sequence pairs

    Returns:
        First 3 outputs as specified in docstring for get_perturbation_calls
        ps_calls_summary is a Pandas dataframe summarizing descriptive statistics for each perturbation_call (unique combination of protospacers) found in
        the dataset, along with some overall summary statistics about the mulitplicty of infection

    """
    if filtered_guide_counts_matrix is None:
        return (None, None, None, None)
    (ps_calls_table, presence_calls, cells_with_ps, umi_thresholds)  = get_perturbation_calls(filtered_guide_counts_matrix,
                                                                                f_map,)

    ps_calls_table.sort_values(by=['feature_call'], inplace=True, kind='mergesort')
    ps_calls_summary = get_ps_calls_summary(ps_calls_table, filtered_guide_counts_matrix)

    return (ps_calls_table, presence_calls, cells_with_ps, ps_calls_summary, umi_thresholds)

def get_perturbation_calls(filtered_guide_counts_matrix, feature_map):
    """For each barcode in the gex cell list, calculates the protospacers present in it.
    Args:
        filtered_guide_counts_matrix: CountMatrix - obtained by selecting features by CRISPR library type on the feature counts matrix
        feature_map: dict - (feature_name:feature_barcode) pairs

    Returns:
        (calls_df, presence_calls, cells_with_ps)
            calls_df: pandas DataFrame indexed by cell-barcode with 3 columns:(num_features, feature_call, num_umis)
                        num_features: int - number of features in that cell
                        feature_call = "|" de-limited string of feature_names
                        num_umis "|" de-limited string of umi_counts for each feature
            presence_calls: dict of (feature_name: list_of_presence_calls) where list_of_presence_calls[i] is True if cell_list[i] contains said feature
            cells_with_ps: dict of (feature_name: list_of_cells_with_that_feature) pairs
    """
    presence_calls = {}
    cells_with_ps = {}
    if filtered_guide_counts_matrix is None:
        return (None, None, None)

    filtered_bcs = filtered_guide_counts_matrix.bcs
    umi_thresholds = {}

    for feature_id in sorted(feature_map.keys()):
        feature_id_int = filtered_guide_counts_matrix.feature_id_to_int(feature_id)
        this_feature_view = filtered_guide_counts_matrix.select_features([feature_id_int])
        umi_counts = this_feature_view.get_counts_per_bc()

        in_high_umi_component = _call_presence_with_gmm(umi_counts)
        presence_calls[feature_id] = in_high_umi_component
        cells_with_ps[feature_id] = [filtered_bcs[i] for i in np.flatnonzero(np.array(in_high_umi_component))]

        if np.any(in_high_umi_component):
            umi_thresholds[feature_id] = np.amin(umi_counts[in_high_umi_component])

    calls_per_cell = _get_calls_per_cell(cells_with_ps)
    calls_df = _get_cell_calls_df(calls_per_cell, filtered_guide_counts_matrix)
    return (calls_df, presence_calls, cells_with_ps, umi_thresholds)

def _get_num_cells_without_guide_umis(guide_counts_matrix):
    counts_per_bc = guide_counts_matrix.get_counts_per_bc()
    return np.sum(counts_per_bc==0)

def get_ps_calls_summary(ps_calls_table, guide_counts_matrix):
    # assumes ps_calls_table is sorted by feature_call
    num_gex_cbs = len(guide_counts_matrix.bcs)
    num_cells_with_ps_calls = len(set(ps_calls_table.index))
    num_cells_without_ps_calls = num_gex_cbs - num_cells_with_ps_calls
    assert num_cells_without_ps_calls >= 0
    frac_cells_without_ps_calls = tk_stats.robust_divide(num_cells_without_ps_calls, num_gex_cbs)

    num_cells_without_guide_umis = _get_num_cells_without_guide_umis(guide_counts_matrix)
    frac_cells_without_guide_umis = tk_stats.robust_divide(num_cells_without_guide_umis, num_gex_cbs)

    column_titles = ['# Cells', '% Cells', 'Median UMIs', 'Std. Dev UMIs']
    ps_summary = pd.DataFrame(columns = column_titles)

    ps_summary.loc['No guide molecules'] = (num_cells_without_guide_umis,
                                                    100*frac_cells_without_guide_umis,
                                                    'N/A',
                                                    'N/A')
    ps_summary.loc['No confident call'] = (num_cells_without_ps_calls,
                                                    100*frac_cells_without_ps_calls,
                                                    'N/A',
                                                    'N/A')

    singlets_table = ps_calls_table[ps_calls_table['num_features']==1]
    num_cells_with_singlets = len(set(singlets_table.index))
    frac_cells_with_singlets = tk_stats.robust_divide(num_cells_with_singlets, num_gex_cbs)

    ps_summary.loc['1 protospacer expressed'] = (num_cells_with_singlets,
                                                        100*frac_cells_with_singlets,
                                                        'N/A',
                                                        'N/A')

    multiplets_table = ps_calls_table[ps_calls_table['num_features']>1]
    num_cells_with_multiplets = len(set(multiplets_table.index))
    frac_cells_with_multiplets = tk_stats.robust_divide(num_cells_with_multiplets, num_gex_cbs)

    ps_summary.loc['> 1 protospacer expressed'] = (num_cells_with_multiplets,
                                                        100*frac_cells_with_multiplets,
                                                        'N/A',
                                                        'N/A')

    for f_call, table_iter in itertools.groupby(ps_calls_table.itertuples(), key = sort_by_feature_call):
        if f_call=='None':
            continue
        num_cells = 0
        umis_per_cell = []

        for row in table_iter:
            num_cells += 1
            if row.num_features > 1:
                umis_per_cell.append(sum([ float(x) for x in (row.num_umis).split('|') ] ))
                this_feature_call = " - ".join(f_call.split('|'))
            else:
                umis_per_cell.append(float(row.num_umis))
                this_feature_call = f_call

        ps_summary.loc[this_feature_call] = (num_cells,
                                                100*tk_stats.robust_divide(num_cells, num_gex_cbs),
                                                np.median(umis_per_cell),
                                                np.std(umis_per_cell)
                                            )
    return ps_summary

def sort_by_feature_call(row):
    return row.feature_call

def _get_calls_per_cell(cells_with_ps):
    calls_per_cell = {}

    for ps in cells_with_ps:
        bc_list = cells_with_ps[ps]

        for bc in bc_list:
            this_bc_calls = calls_per_cell.get(bc, [])
            this_bc_calls.append(ps)
            calls_per_cell[bc] = this_bc_calls

    return calls_per_cell

def _get_cell_calls_df(calls_per_cell, filtered_guide_counts_matrix):
    columns = ['num_features', 'feature_call', 'num_umis']
    calls_per_cell_df = pd.DataFrame(columns = columns)

    for cell in calls_per_cell:
        calls = calls_per_cell.get(cell)
        num_features = len(calls)

        if num_features > 1:
            features = "|".join(calls)
            umis = "|".join([str(filtered_guide_counts_matrix.get(ps_id, cell))
                                    for ps_id in calls])
            calls_per_cell_df.loc[cell] = (num_features, features, umis)
        else:
            calls_per_cell_df.loc[cell] = (num_features,
                                           calls[0],
                                          str(filtered_guide_counts_matrix.get(calls[0], cell)))

    return calls_per_cell_df
