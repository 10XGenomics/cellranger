#!/usr/bin/env python
#
# Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
#
import cellranger.constants as cr_constants
import cellranger.utils as cr_utils
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
    def __init__(self, raw_matrices=None, filtered_matrices=None):
        self.raw_matrices = raw_matrices
        self.filtered_matrices = filtered_matrices

        self.summary = {}
        self.result = {}

        # Fix random seed
        np.random.seed(0)

    def is_zero_matrix(self):
        if self.filtered_matrices is None:
            return True
        if all(m.bcs_dim == 0 or m.genes_dim == 0 for m in self.filtered_matrices.matrices.itervalues()):
            return True
        return False

    @staticmethod
    def _bootstrap_inferred_multiplets(counts0, counts1):
        n_multiplet_obs = sum(MultiGenomeAnalysis._classify_gems(counts0, counts1) == cr_constants.GEM_CLASS_MULTIPLET)

        n_multiplet_boot = np.zeros(cr_constants.NUM_MULTIPLET_BOOTSTRAP_SAMPLES)
        for i in xrange(0, cr_constants.NUM_MULTIPLET_BOOTSTRAP_SAMPLES):
            boot_idx = np.random.choice(len(counts0), len(counts0))
            counts0_boot = counts0[boot_idx]
            counts1_boot = counts1[boot_idx]
            gem_cls_boot = MultiGenomeAnalysis._classify_gems(counts0_boot, counts1_boot)
            n_obs_multiplet_boot = sum(gem_cls_boot == cr_constants.GEM_CLASS_MULTIPLET)
            n_cells0_boot = sum(gem_cls_boot == cr_constants.GEM_CLASS_GENOME0)
            n_cells1_boot = sum(gem_cls_boot == cr_constants.GEM_CLASS_GENOME1)
            n_multiplet_boot[i] = MultiGenomeAnalysis._infer_multiplets_from_observed(n_obs_multiplet_boot, n_cells0_boot, n_cells1_boot)
        return n_multiplet_obs, n_multiplet_boot

    @staticmethod
    def _classify_gems(counts0, counts1):
        """ Infer number of distinct transcriptomes present in each GEM (1 or 2) and
            report cr_constants.GEM_CLASS_GENOME0 for a single cell w/ transcriptome 0,
            report cr_constants.GEM_CLASS_GENOME1 for a single cell w/ transcriptome 1,
            report cr_constants.GEM_CLASS_MULTIPLET for multiple transcriptomes """
        # Assumes that most of the GEMs are single-cell; model counts independently
        thresh0, thresh1 = [cr_constants.DEFAULT_MULTIPLET_THRESHOLD] * 2
        if sum(counts0 > counts1) >= 1 and sum(counts1 > counts0) >= 1:
            thresh0 = np.percentile(counts0[counts0 > counts1], cr_constants.MULTIPLET_PROB_THRESHOLD)
            thresh1 = np.percentile(counts1[counts1 > counts0], cr_constants.MULTIPLET_PROB_THRESHOLD)

        doublet = np.logical_and(counts0 >= thresh0, counts1 >= thresh1)
        dtype = np.dtype('|S%d' % max(len(cls) for cls in cr_constants.GEM_CLASSES))
        result = np.where(doublet, cr_constants.GEM_CLASS_MULTIPLET, cr_constants.GEM_CLASS_GENOME0).astype(dtype)
        result[np.logical_and(np.logical_not(result == cr_constants.GEM_CLASS_MULTIPLET), counts1 > counts0)] = cr_constants.GEM_CLASS_GENOME1

        return result

    @staticmethod
    def _infer_multiplets_from_observed(n_obs_multiplets, n_cells0, n_cells1):
        """ Given a number of observed multiplets and cell counts for two transcriptomes,
        infer the total number of multiplets (observed + unobserved) """

        if n_cells0 == 0 or n_cells1 == 0:
            return 0

        # Prior probability of a doublet given counts for each cell type (ignore N_cells > 2)
        p_obs_multiplet = 2*(float(n_cells0)/float(n_cells0+n_cells1))*(float(n_cells1)/float(n_cells0+n_cells1))

        # Brute force MLE of binomial n
        n_mle = 0
        if n_obs_multiplets > 0:
            likelihood = scipy.stats.binom.pmf(n_obs_multiplets, xrange(0, n_cells0 + n_cells1), p_obs_multiplet)
            n_mle = np.argmax(likelihood)
        return n_mle

    @staticmethod
    def _compute_count_purity(counts0, counts1):
        """ Compute fraction of counts in putative single-cell GEMs
        originating from the non-cell transcriptome """
        gem_occupancy = MultiGenomeAnalysis._classify_gems(counts0, counts1)
        frac0 = counts0.astype(float) / (counts0 + counts1).astype(float)
        purity0 = frac0[gem_occupancy == cr_constants.GEM_CLASS_GENOME0]
        purity1 = 1 - frac0[gem_occupancy == cr_constants.GEM_CLASS_GENOME1]
        overall_purity = np.concatenate([purity0, purity1])

        # Compute number of purity outliers
        threshold0, threshold1 = 1.0, 1.0
        fit_purity0 = purity0[np.logical_and(purity0 > 0, purity0 < 1)]
        fit_purity1 = purity1[np.logical_and(purity1 > 0, purity1 < 1)]
        if len(fit_purity0) > 1 and len(fit_purity1) > 1:
            try:
                alpha0, beta0, _, _ = scipy.stats.beta.fit(fit_purity0, floc=0, fscale=1)
                alpha1, beta1, _, _ = scipy.stats.beta.fit(fit_purity1, floc=0, fscale=1)
                threshold0 = scipy.stats.beta.ppf(cr_constants.COUNT_PURITY_OUTLIER_PROB_THRESHOLD, alpha0, beta0)
                threshold1 = scipy.stats.beta.ppf(cr_constants.COUNT_PURITY_OUTLIER_PROB_THRESHOLD, alpha1, beta1)
            except scipy.stats._continuous_distns.FitSolverError as e:
                print >> sys.stderr, e
                threshold0, threshold1 = 1.0, 1.0
            except scipy.stats._continuous_distns.FitDataError as e:
                print >> sys.stderr, e
                threshold0, threshold1 = 1.0, 1.0

        outlier0 = np.logical_and(gem_occupancy == cr_constants.GEM_CLASS_GENOME0,
                                  frac0 < threshold0)
        outlier1 = np.logical_and(gem_occupancy == cr_constants.GEM_CLASS_GENOME1,
                                  (1-frac0) < threshold1)
        n_outlier0 = sum(outlier0)
        n_outlier1 = sum(outlier1)
        frac_outlier0 = tk_stats.robust_divide(n_outlier0, len(purity0))
        frac_outlier1 = tk_stats.robust_divide(n_outlier1, len(purity1))
        is_outlier = np.logical_or(outlier0, outlier1).astype(int)

        return (purity0.mean(), purity1.mean(), overall_purity.mean(),
                n_outlier0, n_outlier1, frac_outlier0, frac_outlier1,
                is_outlier)

    def run_all(self):
        d = {}

        # Compute N_cells > 1 between top two genomes
        txome_counts = [mat.m.sum() for mat in self.raw_matrices.matrices.values()]
        top_txome_idx = sorted(np.argsort(txome_counts)[::-1][0:2])
        top_two_txomes = [self.raw_matrices.matrices.keys()[i] for i in top_txome_idx]
        top_txome_cell_bc_seqs = [self.filtered_matrices.matrices.values()[i].bcs for i in top_txome_idx]
        use_barcodes = self.raw_matrices.union_barcodes(top_txome_cell_bc_seqs, top_two_txomes)

        # Don't compute multiplet / purity metrics if no cells detected
        if len(use_barcodes) == 0:
            return

        top_txome_reads_per_bc = self.raw_matrices._get_stacked_reads_per_bc(top_two_txomes, use_barcodes)

        n_multiplet_obs, n_multiplet_boot = MultiGenomeAnalysis._bootstrap_inferred_multiplets(
            top_txome_reads_per_bc[0,:],
            top_txome_reads_per_bc[1,:])
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
         is_outlier) = MultiGenomeAnalysis._compute_count_purity(top_txome_reads_per_bc[0,:],
                                                                 top_txome_reads_per_bc[1,:])
        d['%s_filtered_bcs_mean_count_purity' % top_two_txomes[0]] = purity0
        d['%s_filtered_bcs_mean_count_purity' % top_two_txomes[1]] = purity1
        d['%s_filtered_bcs_mean_count_purity' % cr_constants.MULTI_REFS_PREFIX] = overall_purity
        d['%s_filtered_bcs_purity_outliers' % top_two_txomes[0]] = n_purity_outlier0
        d['%s_filtered_bcs_purity_outliers' % top_two_txomes[1]] = n_purity_outlier1
        d['%s_filtered_bcs_frac_purity_outlier' % top_two_txomes[0]] = frac_purity_outlier0
        d['%s_filtered_bcs_frac_purity_outlier' % top_two_txomes[1]] = frac_purity_outlier1
        d['%s_filtered_bcs_frac_purity_outlier' % cr_constants.MULTI_REFS_PREFIX] = frac_purity_outlier0 + frac_purity_outlier1

        # Report N_cell classifications on observed data
        gem_class_call = MultiGenomeAnalysis._classify_gems(
            top_txome_reads_per_bc[0,:],
            top_txome_reads_per_bc[1,:])
        self.result = {
            'barcode': use_barcodes,
            'call': gem_class_call.tolist(),
            'count0': top_txome_reads_per_bc[0,:].tolist(),
            'count1': top_txome_reads_per_bc[1,:].tolist(),
            'genome0': top_two_txomes[0],
            'genome1': top_two_txomes[1],
            'purity_outlier': is_outlier.tolist(),
        }
        self.summary = d

    def save_gem_class_csv(self, base_dir):
        csv_file_path = os.path.join(base_dir, 'gem_classification.csv')
        cr_utils.makedirs(os.path.dirname(csv_file_path), allow_existing=True)
        with open(csv_file_path, 'wb') as f:
            writer = csv.writer(f)
            writer.writerow(['barcode',
                             self.result['genome0'],
                             self.result['genome1'],
                             'call'
                         ])
            for i in xrange(len(self.result['barcode'])):
                call = self.result['call'][i]
                call = call.replace(cr_constants.GEM_CLASS_GENOME0, self.result['genome0'])
                call = call.replace(cr_constants.GEM_CLASS_GENOME1, self.result['genome1'])
                writer.writerow([
                    self.result['barcode'][i],
                    self.result['count0'][i],
                    self.result['count1'][i],
                    call,
                ])

    def save_gem_class_json(self, base_dir):
        json_file_path = MultiGenomeAnalysis.json_path(base_dir)
        cr_utils.makedirs(os.path.dirname(json_file_path), allow_existing=True)
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
