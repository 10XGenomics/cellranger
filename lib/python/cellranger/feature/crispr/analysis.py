import itertools
import numpy as np
import pandas as pd
pd.set_option("compute.use_numexpr", False)
import tenkit.stats as tk_stats
import collections
import sys
import cellranger.analysis.diffexp as cr_diffexp
import cellranger.stats as cr_stats
import random
import cellranger.library_constants as lib_constants

NUM_BOOTSTRAPS = 500 # number of bootstrap draws to do for calculating
                     # empirical confidence intervals for perturbation efficiencies
CI_LOWER_BOUND = 5.0   # CI lower bound (ie percentile value) for perturbation efficiencies
CI_UPPER_BOUND = 95.0  # CI upper bound (ie percentile value) for perturbation efficiencies
FILTER_LIST = ['None', 'Non-Targeting', 'Ignore'] # List of targets that are considered filter-able ie
                                                    # for which we can't or won't compute perturbation efficiencies
CONTROL_LIST = ['Non-Targeting'] # Target IDs used for specifying control perturbations
MIN_NUMBER_CELLS_PER_PERTURBATION = 10 # Minimum number of cells a perturbation has to have before we compute differential expression for it

UMI_NUM_TRIES = 10 # Number of initial points to try for GMM-fitting
UMI_MIX_INIT_SD = 0.25 # Intial standard deviation for GMM components

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

def _add_bcs_without_ps_calls(bc_targets_dict, bcs):
    bcs_without_ps_calls = list(set(bcs).difference(set(bc_targets_dict.keys())))
    for bc in bcs_without_ps_calls:
        bc_targets_dict[bc] = {'target_name': 'None', 'num_features': -1, 'target_id': 'None', 'feature_call': 'None'}

    return bc_targets_dict

def _get_bc_targets_dict(target_info, protospacers_per_cell, ignore_multiples=True,
                                 control_list = ['Non-Targeting'], filter_list = FILTER_LIST):
    bc_targets_dict = {}

    for this_row in protospacers_per_cell.itertuples():
        feature_call = this_row.feature_call

        if feature_call in target_info:
            this_target_id = target_info.get(feature_call).get('target_id')
            this_target_name = target_info.get(feature_call).get('target_name')
        else:
            this_target_id = this_target_name = feature_call

        if this_row.num_features > 1:
            if ignore_multiples:
                this_target_id = this_target_name = 'Ignore'
            else:
                this_features = feature_call.split('|')
                this_target_ids = [target_info.get(x).get('target_id') for x in this_features]
                this_target_names = [target_info.get(x).get('target_name') for x in this_features]

                if set(this_target_ids)==set(control_list):
                    this_target_id = this_target_name = 'Non-Targeting'
                else:
                    this_target_ids = list(set(this_target_ids).difference(set(filter_list)))
                    this_target_names = list(set(this_target_names).difference(set(filter_list)))
                    if this_target_ids:
                        this_target_id = "|".join(this_target_ids)
                        this_target_name = "|".join(this_target_names)
                    else:
                        this_target_id = this_target_name = 'Ignore'


        bc_targets_dict[this_row.Index] = {'num_features' : this_row.num_features,
                                           'feature_call' : feature_call,
                                            'target_id': this_target_id,
                                            'target_name': this_target_name,
                                           }
    return bc_targets_dict

def _get_target_name_from_id(target_id_name_map, target_id, sep="|"):
    if sep not in target_id:
        return target_id_name_map.get(target_id, target_id)
    return "|".join([target_id_name_map.get(x, x) for x in target_id.split(sep) ])

def _should_filter(perturbation_name,
                    feature_ref_table,
                    target_info,
                    filter_list = FILTER_LIST,
                    by_feature=False):
    target_id_name_map = _get_target_id_name_map(feature_ref_table, by_feature)
    target_tuple = _get_target_id_from_name(
                                            perturbation_name,
                                            target_id_name_map,
                                            target_info,
                                            by_feature,)
    return all([x in filter_list for x in target_tuple[1]])


def _get_target_id_from_name(this_perturbation_name, target_id_name_map,
                               target_info, by_feature):
    if "|" not in this_perturbation_name:
        return ([this_perturbation_name], [target_id_name_map[this_perturbation_name]])

    p_names = this_perturbation_name.split('|')
    if not(by_feature):
        return (p_names, [target_id_name_map[p_name] for p_name in p_names])

    return (p_names, [target_id_name_map[p_name] for p_name in p_names])

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
    (ps_calls_table, presence_calls, cells_with_ps)  = get_perturbation_calls(filtered_guide_counts_matrix,
                                                                                f_map,)

    ps_calls_table.sort_values(by=['feature_call'], inplace=True, kind='mergesort')
    ps_calls_summary = get_ps_calls_summary(ps_calls_table, len(filtered_guide_counts_matrix.bcs))

    return (ps_calls_table, presence_calls, cells_with_ps, ps_calls_summary)

def get_perturbation_efficiency(  feature_ref_table,
                                    protospacers_per_cell,
                                    feature_count_matrix,
                                    by_feature=False,
                                    ignore_multiples=False,
                                    filter_list = FILTER_LIST,
                                    num_bootstraps = NUM_BOOTSTRAPS,
                                    ci_lower = CI_LOWER_BOUND,
                                    ci_upper = CI_UPPER_BOUND,
                                    sseq_params=None,
                                 ):

    """ Calculates log2 fold change and empirical confidence intervals for log2 fold change for target genes.
    Args:
        feature_ref_table: pandas DataFrame - obtained by reading the feature reference file
        protospacers_per_cell: pandas DataFrame - protospacer calls per cell. Output of protospacer calling stage
        feature_count_matrix: CountMatrix - Feature Barcode Matrix
        by_feature: bool - if True, cells are grouped by the combination of protospacers present in them, rather than the gene targets of those protospacers.
                            Right now, the pipeline measures perturbations both by feature (by_feature=True) or by target (by_feature=False) by calling
                            the MEASURE_PERTURBATIONS_PD stage twice using the "as" keyword in Martian 3.0
        ignore_multiples: bool - if True, cells with multiple protospacers are ignored from the analysis. Default is False
        filter_list: list - list of target id values to be filtered out of differential expression calculations
        num_bootstraps: int - number of bootstrap draws for computing empirical confideence interval for log2 fold change
        ci_lower: float -  percentile value for calculating the lower bound of the emprirical confidence interval for log2 fold change
        ci_upper: float - percentile value for calculating the upper bound of the emprirical confidence interval for log2 fold change
        sseq_params: dict - global parameters for differential expression calculations. If None (default),
                    computed by cr_diffexp.compute_sseq_params

    Returns:
        (log2_fold_change, fold_change_CI): tuple of two dictionaries -
            log2_fold_change: keys are perturbations, grouped by features or by targets.
                              Values are a dictionary where the keys are individual features or targets,
                                and values are a tuple of floats, each with 4 entries.
                              First entry is the log2_fold_change, second entry is the p-value.=,
                              third entry is the sum of UMI counts for that target gene across all cells in the condition,
                              fourth entry is the sum of UMI counts for that target gene across all cells in the control.
                              Eg. {'COA3': {'COA3': (-1.4, 0.02, 5, 20)},
                                   'COX18|COX6C': {'COX18': (-0.87, 1.0, 3, 40), ...} ...}

            fold_change_CI: similar nested dict, except the tuple indicates the confidence interval.
                            First entry is the lower percentile bound and second is the
                            upper percentile bound.
                            Eg. {'COA3': {'COA3': (-2.12, -0.93)},
                                 'COX18|COX6C': {'COX18': (-1.16, -0.57),...}... }
    """

    matrix = feature_count_matrix.select_features_by_type(lib_constants.GENE_EXPRESSION_LIBRARY_TYPE)
    target_id_name_map = _get_target_id_name_map(feature_ref_table, by_feature)
    target_info = _get_target_info(feature_ref_table)
    (target_calls, perturbation_keys) = _get_ps_clusters(
                                                feature_ref_table,
                                                protospacers_per_cell,
                                                matrix.bcs,
                                                by_feature,
                                                ignore_multiples)

    num_cells_per_perturbation = {name : sum(target_calls==number)  for (number, name) in perturbation_keys.iteritems() }

    diff_exp_results =  _analyze_transcriptome(matrix,
                                                      target_calls,
                                                      perturbation_keys,
                                                      feature_ref_table,
                                                      by_feature,
                                                      target_info,
                                                      sseq_params,
                                                      filter_list)
    if diff_exp_results is None:
        return (None, None, None)

    results_per_perturbation = diff_exp_results['results_per_perturbation']
    log2_fold_change = _get_log2_fold_change(
                                                                results_per_perturbation,
                                                                target_id_name_map,
                                                                target_info,
                                                                by_feature)

    computed_params = diff_exp_results['sseq_params']



    fold_change_CI = _get_confidence_intervals(matrix, target_info,
                                    results_per_perturbation, computed_params,
                                                            perturbation_keys, target_calls,
                                                            target_id_name_map, filter_list,
                                                            num_bootstraps, ci_lower, ci_upper,
                                                            by_feature)


    return (log2_fold_change, fold_change_CI, num_cells_per_perturbation)

def construct_df(f_change, f_change_ci, num_cells_per_perturbation, by_feature):
    if (f_change is None) or (f_change_ci is None):
        return None

    if by_feature:
        target_string = 'Target Guide'
    else:
        target_string = 'Target Gene'

    column_list = ['Perturbation',
                    target_string,
                    'Log2 Fold Change',
                    'p Value',
                    'Log2 Fold Change Lower Bound',
                    'Log2 Fold Change Upper Bound',
                    'Cells with Perturbation',
                    'Mean UMI Count Among Cells with Perturbation',
                    'Cells with Non-Targeting Guides',
                    'Mean UMI Count Among Cells with Non-Targeting Guides',
                  ]
    this_df = pd.DataFrame(columns = column_list)
    counter = 0
    control_num_cells = num_cells_per_perturbation['Non-Targeting']
    for key in sorted(f_change.keys()):
        this_key_results = f_change.get(key)
        this_key_ci = f_change_ci.get(key)
        if this_key_results is None:
            continue

        this_num_cells = num_cells_per_perturbation[key]

        for (ps, results) in this_key_results.iteritems():
            lower_bound = this_key_ci.get(ps)[0]
            upper_bound = this_key_ci.get(ps)[1]
            this_df.loc[counter] = (key,
                                    ps,
                                    results[0],
                                    results[1],
                                    lower_bound,
                                    upper_bound,
                                    this_num_cells,
                                    tk_stats.robust_divide(results[2], this_num_cells),
                                    control_num_cells,
                                    tk_stats.robust_divide(results[3], control_num_cells)
                                    )
            counter += 1
    this_df.sort_values(by=['Log2 Fold Change'], ascending = True, inplace = True)
    return this_df

def _analyze_transcriptome(matrix,
                                    target_calls,
                                    perturbation_keys,
                                    feature_ref_table,
                                    by_feature,
                                    target_info,
                                    sseq_params=None,
                                    filter_list = FILTER_LIST):

    """ Compute differential expression for each perturbation vs non-targeting control
        Args: matrix : CountMatrix (gene expression data, sub-selected from all feature expression data in the CountMatrix)
              target_calls: np.array(int) - integer perturbation target labels.
                                            Numbering starts at 1
              perturbation_keys: dict - (cluster_number:perturbation_name) pairs
              feature_ref_table: pandas DataFrame - obtained by reading the feature reference file
              by_feature: bool - if True, cells are grouped by features (rather than by target)
              target_info: dict - Nested dict: {feature1: {'target_id': value1, 'target_name': value2}, ...}
              sseq_params: dict - params from compute_sseq_params
              filter_list: list - list of target ids to be filtered out of the analysis

        Outs:
            dict -
                    {'results_all_perturbations': named tuple containing a dataframe of differential expression outputs,
                    'results_per_perturbation': dict of perturbation_name:diffexp_dataframe pairs,
                    'sseq_params': dict of global SSEQ diffexp params}
    """

    DIFFERENTIAL_EXPRESSION = collections.namedtuple('DIFFERENTIAL_EXPRESSION', ['data'])
    results_per_perturbation = {}
    n_clusters = len(perturbation_keys.keys())
    filter_cluster_indices = [x for x in perturbation_keys if perturbation_keys[x] in filter_list]
    n_effective_clusters = n_clusters - len(filter_cluster_indices)
    if sseq_params is None:
        print "Computing params..."
        sys.stdout.flush()
        sseq_params = cr_diffexp.compute_sseq_params(matrix.m)

    # Create a numpy array with 3*k columns, where k is the number of perturbations
    # k = n_clusters - 2 since we don't need to report results for ignore and non-targeting
    # each group of 3 columns is mean, log2, pvalue for cluster i
    all_de_results = np.zeros((matrix.features_dim, 3*n_effective_clusters))

    nt_indices = [x for x in perturbation_keys if perturbation_keys[x]=='Non-Targeting']
    if (nt_indices is None) or (nt_indices==[]):
        return
    nt_index = nt_indices[0]

    in_control_cluster = target_calls == nt_index
    column_counter = 0
    feature_defs = matrix.feature_ref.feature_defs
    gene_ids = [feature_def.id for feature_def in feature_defs]
    gene_names = [feature_def.name for feature_def in feature_defs]

    for cluster in perturbation_keys:
        perturbation_name = perturbation_keys.get(cluster)
        if (cluster in filter_cluster_indices) or _should_filter(perturbation_name,
                                                                    feature_ref_table,
                                                                    target_info,
                                                                    filter_list,
                                                                    by_feature):
            continue

        in_cluster = target_calls == cluster
        group_a = np.flatnonzero(in_cluster)

        if len(group_a) < MIN_NUMBER_CELLS_PER_PERTURBATION:
            continue

        group_b = np.flatnonzero(in_control_cluster)
        print 'Computing DE for perturbation %s...' % perturbation_name
        sys.stdout.flush()

        de_result = cr_diffexp.sseq_differential_expression(
                                matrix.m, group_a, group_b, sseq_params)
        de_result['Gene ID'] = gene_ids
        de_result['Gene Name'] = gene_names
        results_per_perturbation[perturbation_name] = de_result
        all_de_results[:, 0+3*(column_counter)] = de_result['norm_mean_a']
        all_de_results[:, 1+3*(column_counter)] = de_result['log2_fold_change']
        all_de_results[:, 2+3*(column_counter)] = de_result['adjusted_p_value']
        column_counter += 1

    return {'results_all_perturbations':DIFFERENTIAL_EXPRESSION(all_de_results),
            'results_per_perturbation':results_per_perturbation,
            'sseq_params':sseq_params,
           }

def _get_confidence_intervals(matrix, target_info,
                                        results_per_perturbation,
                                        computed_params,
                                        perturbation_keys,
                                        target_calls,
                                        target_id_name_map,
                                        filter_list = FILTER_LIST,
                                        num_bootstraps=100,
                                        ci_lower = CI_LOWER_BOUND,
                                        ci_upper = CI_UPPER_BOUND,
                                        by_feature=False):

    log2_fold_change_CIs = {}
    nt_index = [x for x in perturbation_keys if perturbation_keys[x]=='Non-Targeting'][0]
    in_control_cluster = target_calls == nt_index
    cluster_from_name = {b:a for (a,b) in perturbation_keys.iteritems()}

    for perturbation in results_per_perturbation:
        (this_names, this_ids) = _get_target_id_from_name(perturbation,
                                                             target_id_name_map,
                                                            target_info,
                                                          by_feature,)
        this_cluster = cluster_from_name[perturbation]
        in_cluster = target_calls == this_cluster
        group_a = np.flatnonzero(in_cluster)
        group_b = np.flatnonzero(in_control_cluster)
        log2_fold_change_CIs[perturbation] = {name:
                                                    _get_fold_change_cis(
                                                              matrix,
                                                              target,
                                                              perturbation,
                                                              group_a,
                                                              group_b,
                                                              computed_params,
                                                              num_bootstraps,
                                                              ci_lower,
                                                              ci_upper)
                                                    for (name, target) in zip(this_names, this_ids)
                                                    if target not in filter_list
                                                 }


    return log2_fold_change_CIs

def _get_matrix_target_genes(matrix,
                                    perturbation_name,
                                    target_id_name_map,
                                    target_info,
                                    filter_list = FILTER_LIST,
                                    by_feature=False):
    target_ids = _get_target_id_from_name(perturbation_name,
                                         target_id_name_map,
                                         target_info,
                                         by_feature,
                                        )
    filtered_ids = [x for x in target_ids[1] if x not in filter_list]
    target_gene_ints = matrix.gene_ids_to_ints(filtered_ids)

    return matrix.select_genes(target_gene_ints)



def _get_log2_fold_change(results_per_perturbation,
                                     target_id_name_map,
                                     target_info,
                                     by_feature,
                                     filter_list = FILTER_LIST):

    log2_fold_change = {}
    for perturbation in results_per_perturbation:
        this_results = results_per_perturbation[perturbation]
        (this_names, this_ids) = _get_target_id_from_name(perturbation,
                                                         target_id_name_map,
                                                        target_info,
                                                       by_feature,)
        log2_fold_change[perturbation] = {name: _get_ko_per_target(this_results, target)
                                                        for (name, target) in zip(this_names, this_ids)
                                                       if target not in filter_list}

    return log2_fold_change

def _get_fold_change_cis( matrix,
                          target,
                          perturbation,
                          group_a,
                          group_b,
                          computed_params,
                          num_bootstraps,
                          ci_lower,
                          ci_upper):

    target_gene_ints = matrix.feature_ids_to_ints([target])

    this_matrix = matrix.select_features(target_gene_ints)

    return _compute_fold_change_CIs(this_matrix.m, group_a, group_b,
                                        computed_params, num_bootstraps,
                                        ci_lower, ci_upper,)

def _compute_fold_change_CIs(x, cond_a, cond_b, computed_params,
                               num_bootstraps = NUM_BOOTSTRAPS, ci_lower = CI_LOWER_BOUND, ci_upper = CI_UPPER_BOUND):
    log2_fold_change_vals = range(num_bootstraps)
    for i in xrange(num_bootstraps):
        this_cond_a = _sample_with_replacement(cond_a, len(cond_a))
        this_cond_b = _sample_with_replacement(cond_b, len(cond_b))
        x_a = x[:,this_cond_a]
        x_b = x[:,this_cond_b]

        # Size factors
        size_factor_a = np.sum(computed_params['size_factors'][cond_a])
        size_factor_b = np.sum(computed_params['size_factors'][cond_b])

        gene_sums_a = np.squeeze(np.asarray(x_a.sum(axis=1)))
        gene_sums_b = np.squeeze(np.asarray(x_b.sum(axis=1)))

        log2_fold_change_vals[i] = (np.log2((1+gene_sums_a)/(1+size_factor_a)) - \
                                        np.log2((1+gene_sums_b)/(1+size_factor_b)))


    return (np.percentile(log2_fold_change_vals, ci_lower),
            np.percentile(log2_fold_change_vals, ci_upper))

def _sample_with_replacement(list_input, num_samples):
    return [random.choice(list_input) for i in xrange(num_samples)]

def _get_ko_per_target(results, target):
    return (results.loc[results['Gene ID']==target]['log2_fold_change'].values[0],
            results.loc[results['Gene ID']==target]['p_value'].values[0],
            results.loc[results['Gene ID']==target]['sum_a'].values[0],
            results.loc[results['Gene ID']==target]['sum_b'].values[0],
           )

def _get_target_id_name_map(feature_ref_table, by_feature=False):
    """returns a dict of target_gene_name:target_gene_id pairs if by_feature = False
        and a dict of protospacer:target_gene_id pairs if by_feature = True"""
    target_info = _get_target_info(feature_ref_table)
    if not(by_feature):
        return {target_info[x]['target_name'] : target_info[x]['target_id'] for x in target_info}
    return {x : target_info[x]['target_id'] for x in target_info}

def _get_target_info(feature_ref_table):
    """Returns a nested dict. {feature1: {'target_id': value1, 'target_name': value2}, ...}
    """
    target_info = {}
    for this_row in feature_ref_table.itertuples():
        this_target_id = getattr(this_row, 'target_gene_id')
        this_target_name = getattr(this_row, 'target_gene_name', this_target_id)
        target_info[this_row.id] = {'target_id': this_target_id,
                                    'target_name': this_target_name,
                                    }

    return target_info

def _get_ps_clusters(feature_ref_table, protospacers_per_cell, bcs, by_feature=True, ignore_multiples = True):
    """Returns a tuple (target_calls, perturbation_keys).
        target_calls: np.array(int) - identifies the perturbation assigned to each cell in the gene-barcode matrix
        perturbation_keys: dict - (cluster_number:perturbation_name) pairs
    """
    if not(by_feature):
        return _get_ps_clusters_by_target(feature_ref_table, protospacers_per_cell, bcs, ignore_multiples)
    return _get_ps_clusters_by_feature(feature_ref_table, protospacers_per_cell, bcs, ignore_multiples)

def _get_ps_clusters_by_target(feature_ref_table, protospacers_per_cell, bcs, ignore_multiples=True):
    target_info = _get_target_info(feature_ref_table)
    all_target_names = [target_info.get(x).get('target_id') for x in target_info]
    assert "Non-Targeting" in all_target_names, "Perturbation calculations require a Non-Targeting control"

    target_id_name_map = {v.get('target_id'):v.get('target_name') for v in target_info.itervalues()}

    bc_targets_dict = _get_bc_targets_dict(target_info, protospacers_per_cell, ignore_multiples,
                                 control_list = CONTROL_LIST, filter_list = FILTER_LIST)
    bc_targets_dict = _add_bcs_without_ps_calls(bc_targets_dict, bcs)

    unique_targets_list = list(set(
                            [bc_targets_dict.get(bc).get('target_id') for bc in bc_targets_dict]
                                ))
    unique_targets_list.sort()

    target_to_int = {a:b for (a,b) in
                            zip(unique_targets_list, range(1, len(unique_targets_list)+1 ))
                         }
    perturbation_from_cluster = {b: _get_target_name_from_id(target_id_name_map, a, sep="|")
                                     for (a,b) in target_to_int.iteritems()}

    calls_vector = [target_to_int.get(bc_targets_dict.get(bc).get('target_id'))
                                for bc in bcs]

    return (np.asarray(calls_vector),
                perturbation_from_cluster
           )

def _get_ps_clusters_by_feature(feature_ref_table, protospacers_per_cell, bcs, ignore_multiples=True):
    target_info = _get_target_info(feature_ref_table)
    all_target_names = [target_info.get(x).get('target_id') for x in target_info]
    assert "Non-Targeting" in all_target_names, "Perturbation calculations require a Non-Targeting control"

    bc_targets_dict = _get_bc_targets_dict(target_info, protospacers_per_cell, ignore_multiples)
    bc_targets_dict = _add_bcs_without_ps_calls(bc_targets_dict, bcs)

    unique_feature_calls = list(set(
                            [_get_feature_from_bc(bc_targets_dict, bc) for bc in bc_targets_dict]
                                ))
    unique_feature_calls.sort()

    feature_to_int = {a:b for (a,b) in
                            zip(unique_feature_calls, range(1, len(unique_feature_calls)+1 ))
                         }


    perturbation_from_cluster = {b:a for (a,b) in feature_to_int.iteritems()}

    calls_vector = [feature_to_int.get(_get_feature_from_bc(bc_targets_dict, bc))
                                for bc in bcs]

    return (np.asarray(calls_vector),
            perturbation_from_cluster
           )

def _get_feature_from_bc(bc_targets_dict, bc, filter_list = ['Non-Targeting', 'Ignore', 'None']):
    if bc_targets_dict.get(bc).get('target_id') not in filter_list:
        return bc_targets_dict.get(bc).get('feature_call')

    if bc_targets_dict.get(bc).get('target_id') in ['None']:
        return 'Ignore'

    return bc_targets_dict.get(bc).get('target_id')

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

    for feature_id in sorted(feature_map.keys()):
        feature_id_int = filtered_guide_counts_matrix.feature_id_to_int(feature_id)
        this_feature_view = filtered_guide_counts_matrix.select_features([feature_id_int])
        umi_counts = this_feature_view.get_counts_per_bc()

        in_high_umi_component = _call_presence_with_gmm(umi_counts)
        presence_calls[feature_id] = in_high_umi_component
        cells_with_ps[feature_id] = [filtered_bcs[i] for i in np.flatnonzero(np.array(in_high_umi_component))]

    calls_per_cell = _get_calls_per_cell(cells_with_ps)
    calls_df = _get_cell_calls_df(calls_per_cell, filtered_guide_counts_matrix)
    return (calls_df, presence_calls, cells_with_ps)

def get_ps_calls_summary(ps_calls_table, num_gex_cbs):
    # assumes input is sorted by feature_call
    num_cells_with_ps_calls = len(set(ps_calls_table.index))
    num_cells_without_ps_calls = num_gex_cbs - num_cells_with_ps_calls
    frac_cells_without_ps_calls = tk_stats.robust_divide(num_cells_without_ps_calls, num_gex_cbs)
    assert num_cells_without_ps_calls >= 0

    column_titles = ['# Cells', '% Cells', 'Median UMIs', 'Std. Dev UMIs']
    ps_summary = pd.DataFrame(columns = column_titles)

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
