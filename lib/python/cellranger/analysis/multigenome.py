#!/usr/bin/env python
#
# Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
#
import cellranger.analysis.constants as analysis_constants
import cellranger.library_constants as lib_constants
import cellranger.io as cr_io
import csv
import json
import numpy as np
import os
import scipy.stats
import sys
import tenkit.safe_json as tk_safe_json
import tenkit.stats as tk_stats

class MultiGenomeAnalysis:
    """ Analysis of matrices when >1 genome is present """
    def __init__(self, raw_matrix=None, filtered_matrix=None):
        self.raw_matrix = raw_matrix
        self.filtered_matrix = filtered_matrix

        self.summary = {}
        self.result = {}
        self.top_two_txomes = None

        # Fix random seed
        np.random.seed(0)

    def is_zero_matrix(self):
        if self.filtered_matrix is None:
            return True
        if self.filtered_matrix.bcs_dim == 0 or self.filtered_matrix.features_dim == 0:
            return True
        return False

    @staticmethod
    def _bootstrap_inferred_multiplets(counts0, counts1):
        n_multiplet_obs = sum(MultiGenomeAnalysis._classify_gems(counts0, counts1) == analysis_constants.GEM_CLASS_MULTIPLET)

        n_multiplet_boot = np.zeros(analysis_constants.NUM_MULTIPLET_BOOTSTRAP_SAMPLES)
        for i in xrange(0, analysis_constants.NUM_MULTIPLET_BOOTSTRAP_SAMPLES):
            boot_idx = np.random.choice(len(counts0), len(counts0))
            counts0_boot = counts0[boot_idx]
            counts1_boot = counts1[boot_idx]
            gem_cls_boot = MultiGenomeAnalysis._classify_gems(counts0_boot, counts1_boot)
            n_obs_multiplet_boot = sum(gem_cls_boot == analysis_constants.GEM_CLASS_MULTIPLET)
            n_cells0_boot = sum(gem_cls_boot == analysis_constants.GEM_CLASS_GENOME0)
            n_cells1_boot = sum(gem_cls_boot == analysis_constants.GEM_CLASS_GENOME1)
            n_multiplet_boot[i] = MultiGenomeAnalysis._infer_multiplets_from_observed(n_obs_multiplet_boot, n_cells0_boot, n_cells1_boot)
        return n_multiplet_obs, n_multiplet_boot

    @staticmethod
    def _classify_gems(counts0, counts1):
        """ Infer number of distinct transcriptomes present in each GEM (1 or 2) and
            report analysis_constants.GEM_CLASS_GENOME0 for a single cell w/ transcriptome 0,
            report analysis_constants.GEM_CLASS_GENOME1 for a single cell w/ transcriptome 1,
            report analysis_constants.GEM_CLASS_MULTIPLET for multiple transcriptomes """
        # Assumes that most of the GEMs are single-cell; model counts independently
        thresh0, thresh1 = [analysis_constants.DEFAULT_MULTIPLET_THRESHOLD] * 2
        if sum(counts0 > counts1) >= 1 and sum(counts1 > counts0) >= 1:
            thresh0 = np.percentile(counts0[counts0 > counts1], analysis_constants.MULTIPLET_PROB_THRESHOLD*100.0)
            thresh1 = np.percentile(counts1[counts1 > counts0], analysis_constants.MULTIPLET_PROB_THRESHOLD*100.0)

        doublet = np.logical_and(counts0 >= thresh0, counts1 >= thresh1)
        dtype = np.dtype('|S%d' % max(len(cls) for cls in analysis_constants.GEM_CLASSES))
        result = np.where(doublet, analysis_constants.GEM_CLASS_MULTIPLET, analysis_constants.GEM_CLASS_GENOME0).astype(dtype)
        result[np.logical_and(np.logical_not(result == analysis_constants.GEM_CLASS_MULTIPLET), counts1 > counts0)] = analysis_constants.GEM_CLASS_GENOME1

        return result

    @staticmethod
    def _infer_multiplets_from_observed(n_obs_multiplets, n_cells0, n_cells1):
        """ Given a number of observed multiplets and cell counts for two transcriptomes,
        infer the total number of multiplets (observed + unobserved) """

        if n_cells0 == 0 or n_cells1 == 0:
            return 0

        # Prior probability of a doublet given counts for each cell type (ignore N_cells > 2)
        p_obs_multiplet = 2 * (float(n_cells0) / float(n_cells0 + n_cells1)) * (float(n_cells1) / float(n_cells0 + n_cells1))

        # Analytical MOM/MLE of binomial N given p, k
        mle = float(n_obs_multiplets) / p_obs_multiplet
        # In some (artificial) datasets, the mle can be higher than the total number of cells
        # observed.  This occurs when n_obs_multiplets > n_cells0|1. The right way to fix that would be to
        # do inference in a full model that didn't fix some parameters.  In practice, multigenomes
        # are a rare analysis and most data isn't artificial, so we are implementing
        # a small hack instead.
        return min(mle, float(n_obs_multiplets + n_cells0 + n_cells1))

    @staticmethod
    def _compute_count_purity(counts0, counts1):
        """ Compute fraction of counts in putative single-cell GEMs
        originating from the non-cell transcriptome """
        gem_occupancy = MultiGenomeAnalysis._classify_gems(counts0, counts1)
        frac0 = counts0.astype(float) / (counts0 + counts1).astype(float)
        purity0 = frac0[gem_occupancy == analysis_constants.GEM_CLASS_GENOME0]
        purity1 = 1 - frac0[gem_occupancy == analysis_constants.GEM_CLASS_GENOME1]

        # Compute number of purity outliers.
        # Note: This is not used for any important metrics, just for diagnostics.
        threshold0, threshold1 = 1.0, 1.0
        fit_purity0 = purity0[np.logical_and(purity0 > 0, purity0 < 1)]
        fit_purity1 = purity1[np.logical_and(purity1 > 0, purity1 < 1)]
        if len(fit_purity0) > 1 and len(fit_purity1) > 1:
            try:
                alpha0, beta0, _, _ = scipy.stats.beta.fit(fit_purity0, floc=0, fscale=1)
                alpha1, beta1, _, _ = scipy.stats.beta.fit(fit_purity1, floc=0, fscale=1)
                threshold0 = scipy.stats.beta.ppf(analysis_constants.COUNT_PURITY_OUTLIER_PROB_THRESHOLD, alpha0, beta0)
                threshold1 = scipy.stats.beta.ppf(analysis_constants.COUNT_PURITY_OUTLIER_PROB_THRESHOLD, alpha1, beta1)
            except scipy.stats._continuous_distns.FitSolverError as e:
                print >> sys.stderr, e
                threshold0, threshold1 = 1.0, 1.0
            except scipy.stats._continuous_distns.FitDataError as e:
                print >> sys.stderr, e
                threshold0, threshold1 = 1.0, 1.0

        outlier0 = np.logical_and(gem_occupancy == analysis_constants.GEM_CLASS_GENOME0,
                                  frac0 < threshold0)
        outlier1 = np.logical_and(gem_occupancy == analysis_constants.GEM_CLASS_GENOME1,
                                  (1 - frac0) < threshold1)
        n_outlier0 = sum(outlier0)
        n_outlier1 = sum(outlier1)
        frac_outlier0 = tk_stats.robust_divide(n_outlier0, len(purity0))
        frac_outlier1 = tk_stats.robust_divide(n_outlier1, len(purity1))
        is_outlier = np.logical_or(outlier0, outlier1).astype(int)

        # Let the UMI count purity be the total fraction of counts that don't belong,
        #   for all barcodes classified as non-multiplets.
        gems0 = gem_occupancy == analysis_constants.GEM_CLASS_GENOME0
        mean_purity0 = tk_stats.robust_divide(counts0[gems0].sum(),
                                              (counts0[gems0] + counts1[gems0]).sum())

        gems1 = gem_occupancy == analysis_constants.GEM_CLASS_GENOME1
        mean_purity1 = tk_stats.robust_divide(counts1[gems1].sum(),
                                              (counts0[gems1] + counts1[gems1]).sum())

        single_cell = (gem_occupancy == analysis_constants.GEM_CLASS_GENOME0) | \
                      (gem_occupancy == analysis_constants.GEM_CLASS_GENOME1)
        mean_overall_purity = tk_stats.robust_divide(np.maximum(counts0[single_cell],
                                                                counts1[single_cell]).sum(),
                                                     (counts0+counts1)[single_cell].sum())

        return (mean_purity0, mean_purity1, mean_overall_purity,
                n_outlier0, n_outlier1, frac_outlier0, frac_outlier1,
                is_outlier)

    def run_all(self):
        d = {}

        # Compute N_cells > 1 between top two genomes
        genomes = self.raw_matrix.get_genomes()
        genome_raw_mats = [self.raw_matrix.select_features_by_genome(g) for g in genomes]
        genome_filt_mats = [self.filtered_matrix.select_features_by_genome(g) for g in genomes]

        txome_counts = [mat.m.sum() for mat in genome_raw_mats]
        top_txome_idx = sorted(np.argsort(txome_counts)[::-1][0:2])

        top_two_txomes = [genomes[i] for i in top_txome_idx]
        self.top_two_txomes = top_two_txomes

        top_txome_cell_bc_seqs = [genome_filt_mats[i].bcs for i in top_txome_idx]
        use_barcodes = sorted(list(reduce(lambda a, x: a | set(x), top_txome_cell_bc_seqs, set())))

        # Don't compute multiplet / purity metrics if no cells detected
        if len(use_barcodes) == 0:
            return

        # Sum the raw matrices by barcode and stack the two genomes together
        top_two_raw_mats = [genome_raw_mats[i] for i in top_txome_idx]
        top_txome_reads_per_bc = np.vstack(tuple([m.select_barcodes_by_seq(use_barcodes).get_counts_per_bc() for m in top_two_raw_mats]))

        n_multiplet_obs, n_multiplet_boot = MultiGenomeAnalysis._bootstrap_inferred_multiplets(
            top_txome_reads_per_bc[0, :],
            top_txome_reads_per_bc[1, :])
        d['filtered_bcs_observed_multiplets'] = int(round(n_multiplet_obs))
        d['filtered_bcs_inferred_multiplets'] = int(round(n_multiplet_boot.mean()))
        multiplet_rate = tk_stats.robust_divide(n_multiplet_boot.mean(), len(use_barcodes))
        d['filtered_bcs_inferred_multiplet_rate'] = multiplet_rate
        d['filtered_bcs_inferred_normalized_multiplet_rate'] = 1000 * tk_stats.robust_divide(multiplet_rate, len(use_barcodes))
        d['filtered_bcs_inferred_multiplet_rate_lb'] = tk_stats.robust_divide(np.percentile(n_multiplet_boot, 2.5), len(use_barcodes))
        d['filtered_bcs_inferred_multiplet_rate_ub'] = tk_stats.robust_divide(np.percentile(n_multiplet_boot, 97.5), len(use_barcodes))

        (purity0, purity1, overall_purity,
         n_purity_outlier0, n_purity_outlier1,
         frac_purity_outlier0, frac_purity_outlier1,
         is_outlier) = MultiGenomeAnalysis._compute_count_purity(top_txome_reads_per_bc[0, :],
                                                                 top_txome_reads_per_bc[1, :])
        d['%s_filtered_bcs_mean_count_purity' % top_two_txomes[0]] = purity0
        d['%s_filtered_bcs_mean_count_purity' % top_two_txomes[1]] = purity1
        d['%s_filtered_bcs_mean_count_purity' % lib_constants.MULTI_REFS_PREFIX] = overall_purity
        d['%s_filtered_bcs_purity_outliers' % top_two_txomes[0]] = n_purity_outlier0
        d['%s_filtered_bcs_purity_outliers' % top_two_txomes[1]] = n_purity_outlier1
        d['%s_filtered_bcs_frac_purity_outlier' % top_two_txomes[0]] = frac_purity_outlier0
        d['%s_filtered_bcs_frac_purity_outlier' % top_two_txomes[1]] = frac_purity_outlier1
        d['%s_filtered_bcs_frac_purity_outlier' % lib_constants.MULTI_REFS_PREFIX] = frac_purity_outlier0 + frac_purity_outlier1

        # Report N_cell classifications on observed data
        gem_class_call = MultiGenomeAnalysis._classify_gems(
            top_txome_reads_per_bc[0, :],
            top_txome_reads_per_bc[1, :])
        self.result = {
            'barcode': use_barcodes,
            'call': gem_class_call.tolist(),
            'count0': top_txome_reads_per_bc[0, :].tolist(),
            'count1': top_txome_reads_per_bc[1, :].tolist(),
            'genome0': top_two_txomes[0],
            'genome1': top_two_txomes[1],
            'purity_outlier': is_outlier.tolist(),
        }
        self.summary = d

    def save_gem_class_csv(self, base_dir):
        csv_file_path = os.path.join(base_dir, 'gem_classification.csv')
        cr_io.makedirs(os.path.dirname(csv_file_path), allow_existing=True)
        with open(csv_file_path, 'wb') as f:
            writer = csv.writer(f, lineterminator=os.linesep)
            writer.writerow(['barcode',
                             self.result['genome0'],
                             self.result['genome1'],
                             'call'
                             ])
            for i in xrange(len(self.result['barcode'])):
                call = self.result['call'][i]
                call = call.replace(analysis_constants.GEM_CLASS_GENOME0, self.result['genome0'])
                call = call.replace(analysis_constants.GEM_CLASS_GENOME1, self.result['genome1'])
                writer.writerow([
                    self.result['barcode'][i],
                    self.result['count0'][i],
                    self.result['count1'][i],
                    call,
                ])

    def save_gem_class_json(self, base_dir):
        json_file_path = MultiGenomeAnalysis.json_path(base_dir)
        cr_io.makedirs(os.path.dirname(json_file_path), allow_existing=True)
        with open(json_file_path, 'w') as f:
            json.dump(tk_safe_json.json_sanitize(self.result), f, indent=4, sort_keys=True)

    def save_summary_json(self, filename):
        with open(filename, 'w') as f:
            json.dump(tk_safe_json.json_sanitize(self.summary), f, indent=4, sort_keys=True)

    @staticmethod
    def load_default_format(base_dir):
        json_file_path = MultiGenomeAnalysis.json_path(base_dir)
        if os.path.exists(json_file_path):
            return MultiGenomeAnalysis.load_json(json_file_path)
        else:
            return None

    @staticmethod
    def load_json(filename):
        analysis = MultiGenomeAnalysis()
        with open(filename, 'r') as f:
            analysis.result = json.load(f)
        return analysis

    @staticmethod
    def json_path(base_path):
        return os.path.join(base_path, "analysis.json")
