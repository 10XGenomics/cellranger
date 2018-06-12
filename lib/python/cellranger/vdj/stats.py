#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#
import itertools
import numpy as np
import tenkit.stats as tk_stats
import cellranger.stats as cr_stats

def call_vdj_cells(umi_barcode_idx,
                   umi_read_pairs,
                   barcodes,
                   rpu_mix_init_sd, umi_mix_init_sd, verbosity=0):
    """ 1) Take all UMIs with >1 readpair (including single-readpair umis tends to create a degenerate 0-variance component)
        2) Fit a 2-component Gaussian mixture model to the log(N50 [readpairs per UMI] per barcode)
        3) Let the smallest RPU in the higher-mean component be the RPU threshold
        4) Take all UMIs that exceed the RPU threshold
        5) Fit a 2-component Gaussian mixture model to the log(Filtered UMIs per barcode)
        6) Call as cells the BCs with posterior probability > 0.5 of belonging to the higher-mean component

        Args: umi_barcode_idx (np.array(int)) - barcode index for each UMI
              umi_read_pairs (np.array(int)) - readpair count for each UMI
              barcodes (np.array(str)) - barcode strings
              rpu_mix_init_sd - initialize the RPU GMM standard deviation to this value
              umi_mix_init_sd - initialize the UMI GMM standard deviation to this value

        Returns tuple of:
              (rpu_threshold (int) - read pairs per UMI threshold used,
               umi_threshold (int) - UMIs per barcode threshold to call cells
               barcode_support (dict of string:int) - filtered UMI count per barcode
              )
"""
    assert len(umi_read_pairs) == len(umi_barcode_idx)
    assert rpu_mix_init_sd > 0
    assert umi_mix_init_sd > 0

    # Take UMIs with readpairs > 1
    # Col 0: Barcode idx; Col 1: Read pairs
    use_umis = umi_read_pairs > 1
    umi_info = np.hstack((cr_stats.to_col_vec(umi_barcode_idx[use_umis]),
                          cr_stats.to_col_vec(umi_read_pairs[use_umis])))

    # Check if enough UMIs to call anything
    if umi_info.shape[0] < 1:
        print "Warning: Not enough UMIs to call cells"
        return (1, 1, {}, 0.0)

    # Sort by barcode (in-place)
    view = umi_info.view(dtype=[('bc', umi_info.dtype),
                                ('reads', umi_info.dtype)])
    view.sort(order=['bc'], axis=0)

    # Summarize barcodes
    bc_rpu = []
    bc_reads = []

    for _, group_iter in itertools.groupby(enumerate(umi_info[:,0]), key=lambda idx_bc: idx_bc[1]):
        idx, _ = next(group_iter)
        group_len = 1 + sum(1 for _ in group_iter)
        read_pairs = umi_info[idx:(idx + group_len), 1]

        # Require the BC to have at least 2 raw UMIs
        # to reduce the influence of noise BCs.
        if len(read_pairs) < 2:
            continue

        # Estimate local RPU using N50 to reduce the influence
        # of noise UMIs.
        bc_rpu.append(tk_stats.NX(read_pairs, 0.5))
        bc_reads.append(np.sum(read_pairs))

    if len(bc_rpu) < 2:
        print "Warning: Not enough barcodes with enough UMIs to call cells"
        return (1, 1, {}, 0.0)

    bc_rpu = np.array(bc_rpu, dtype=float)
    bc_reads = np.array(bc_reads, dtype=int)

    ## Fit a GMM to log(RPU) per BC
    log_rpu = np.log(bc_rpu)

    # Trim RPU above the 99th percentile RPU of BCs w/ >= N50 reads
    RPU_TRIM_QUANTILE = 0.99
    n50_reads = tk_stats.NX(bc_reads, 0.5)
    p99_rpu = np.percentile(bc_rpu[bc_reads >= n50_reads], 100.0*RPU_TRIM_QUANTILE)
    log_rpu_trimmed = log_rpu[bc_rpu <= p99_rpu]

    if verbosity > 0:
        print 'n50=%0.4f, p99=%0.4f, mean_rpu=%0.4f, raw_mean=%0.4f' % (n50_reads, p99_rpu, np.mean(np.exp(log_rpu_trimmed)), np.mean(bc_rpu))

    # Try various starting means
    RPU_NUM_TRIES = 10
    rpu_try_means = [[0, x] for x in np.linspace(0, np.max(log_rpu_trimmed), 1+RPU_NUM_TRIES)[1:]]
    rpu_gmm = cr_stats.multistart_gmm(data=cr_stats.to_col_vec(log_rpu_trimmed),
                             weights=[0.5, 0.5],
                             means_list=rpu_try_means,
                             sd=rpu_mix_init_sd)

    if not rpu_gmm.converged_:
        print "Warning: EM did not converge for RPU!"

    rpu_posterior = rpu_gmm.predict_proba(cr_stats.to_col_vec(log_rpu))
    high_rpu_component = np.argmax(rpu_gmm.means_)
    in_high_rpu_component = rpu_posterior[:,high_rpu_component] > 0.5

    if np.sum(in_high_rpu_component) > 0:
        rpu_threshold = max(1, int(np.round(np.exp(np.min(log_rpu[in_high_rpu_component])))))
    else:
        rpu_threshold = 1

    if verbosity > 0:
        print "RPPU threshold: %d" % rpu_threshold
        top_reads = np.sum(bc_reads[log_rpu >= np.log(rpu_threshold)])
        frac_top = tk_stats.robust_divide(top_reads, np.sum(bc_reads))
        print "Frac reads in top RPU component: %0.4f" % frac_top


    # Filter UMIs by the computed threshold
    filt_umi_info = umi_info[umi_info[:,1] >= rpu_threshold, :]

    # Count the number of filtered UMIs for each BC
    bc_filt_umis = np.bincount(filt_umi_info[:,0], minlength=len(barcodes))
    bc_filt_umis_nz = bc_filt_umis[np.flatnonzero(bc_filt_umis)]


    ## Fit a GMM to log(UMIs per BC)

    # Add some noise

    N_PARTITIONS = 90000
    ADD_NOISE_RATE = 10000/float(N_PARTITIONS)
    umis = bc_filt_umis_nz

    n_empties = max(0, N_PARTITIONS - len(bc_filt_umis_nz))
    bc_filt_umis_nz = np.concatenate((umis, np.zeros(n_empties, dtype=bc_filt_umis_nz.dtype)))
    bc_filt_umis_nz += np.random.poisson(ADD_NOISE_RATE, len(bc_filt_umis_nz))
    bc_filt_umis_nz = bc_filt_umis_nz[np.flatnonzero(bc_filt_umis_nz)]

    log_umi = np.log(bc_filt_umis_nz)

    # Try various starting means
    UMI_NUM_TRIES = 10
    umi_try_means = [[0, x] for x in np.linspace(0, np.max(log_umi), 1+UMI_NUM_TRIES)[1:]]
    umi_gmm = cr_stats.multistart_gmm(data=cr_stats.to_col_vec(log_umi),
                             weights=[0.5, 0.5],
                             means_list=umi_try_means,
                             sd=umi_mix_init_sd)

    if not umi_gmm.converged_:
        print "Warning: EM did not converge for UMIs!"

    umi_posterior = umi_gmm.predict_proba(cr_stats.to_col_vec(log_umi))
    high_umi_component = np.argmax(umi_gmm.means_)
    in_high_umi_component = umi_posterior[:,high_umi_component] > 0.5

    # Require at least 2 UMIs per barcode
    if np.sum(in_high_umi_component) > 0:
        umi_threshold = max(2, int(np.round(np.exp(np.min(log_umi[in_high_umi_component])))))
    else:
        umi_threshold = 2

    # Fraction of cell calls with high posterior probability
    cell_posterior = umi_posterior[log_umi >= np.log(umi_threshold), high_umi_component]
    confidence = np.mean(cell_posterior >= 0.99)

    # Record filtered UMI counts for each barcode
    bc_support = {barcodes[i]:bc_filt_umis[i] for i in xrange(len(bc_filt_umis))}

    if verbosity > 0:
        print "UMI threshold: %d" % umi_threshold
        good_bc = bc_filt_umis >= umi_threshold
        all_bc_reads = np.bincount(umi_info[:,0], weights=umi_info[:,1], minlength=len(barcodes))

        good_bc_reads = np.sum(all_bc_reads[good_bc])
        frac_top = tk_stats.robust_divide(good_bc_reads, np.sum(all_bc_reads))
        print "Frac reads in top UMI component: %0.4f" % frac_top

        good_bc_umis = np.sum(bc_filt_umis[good_bc])
        frac_top = tk_stats.robust_divide(good_bc_umis, np.sum(bc_filt_umis))
        print "Frac good UMIs in top UMI component: %0.4f" % frac_top

        print "Confidence: %0.2f" % confidence

    return (rpu_threshold, umi_threshold, bc_support, confidence)
