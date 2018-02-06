#!/usr/bin/env python
#
# Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
#
# Compute general quality metrics from sequencing runs
#
import illuminate
import numpy as np
import logging

QC_KEY_READ1='read1'
QC_KEY_READ2='read2'
QC_KEY_READ3='read3'
QC_KEY_READ4='read4'
QC_KEY_INDEX1='index1'
QC_KEY_INDEX2='index2'
QC_KEY_I7='index1'
QC_KEY_I5='index2'

QC_SHORTCUT_MAP = {
    'R1': QC_KEY_READ1,
    'R2': QC_KEY_READ2,
    'R3': QC_KEY_READ3,
    'R4': QC_KEY_READ4,
    'I1': QC_KEY_INDEX1,
    'I2': QC_KEY_INDEX2
}


def get_illumina_sequencing_metrics(flowcell_path):
    """
    Return a dictionary of sequencing metrics from the flowcell located at run_path.

    Generates the following keys (| indicates keys with combinations)

    PhiX_error_worst_tile
    experiment_name
    fwhm_A|C|G|T
    index1|2_PhiX_error_by_cycle
    index1|2_mean_phasing
    index1|2_mean_prephasing
    index1|2_q20_fraction
    index1|2_q20_fraction_by_cycle
    index1|2_q30_fraction
    index1|2_q30_fraction_by_cycle
    intensity_A|C|G|T
    lanecount
    mean_cluster_density
    mean_cluster_density_pf
    percent_pf_clusters
    read1|2_PhiX_error_by_cycle
    read1|2_mean_phasing
    read1|2_mean_prephasing
    read1|2_q20_fraction
    read1|2_q20_fraction_by_cycle
    read1|2_q30_fraction
    read1|2_q30_fraction_by_cycle
    rta_version
    run_id
    signoise_ratio
    start_datetime
    surfacecount
    swathcount
    tilecount
    total_cluster_density
    total_cluster_density_pf

    :param flowcell_path:  The path to the flowcell directory.
    :rtype: dict[string->T]
    """
    results = {}

    ds = illuminate.InteropDataset(flowcell_path)

    # Copy over flowcell layout params
    for (k, v) in ds.meta.flowcell_layout.items():
        results[k] = v

    try:
        results['rta_version'] = ds.meta.rta_version
    except AttributeError:
        results['rta_version'] = "unknown"

    try:
        results['run_id'] = ds.meta.runID
    except AttributeError:
        results['run_id'] = "unknown"

    try:
        results['experiment_name'] = ds.meta.experiment_name
    except AttributeError:
        results['experiment_name'] = "unknown"

    try:
        results['start_datetime'] = str(ds.meta.start_datetime)
    except AttributeError:
        results['start_datetime'] = "unknown"

    # Determine the indices we have for our 'semantic' reads
    read_config = ds.meta.read_config
    n_reads = len(read_config)

    read_length = np.array([read_config[x]['cycles'] for x in range(0, n_reads)])
    end_cycle = list(read_length.cumsum())
    start_cycle = [0] + end_cycle[:-1]

    main_reads = [x['read_num'] for x in read_config if x['is_index'] == False]
    r1_id = main_reads[0]
    if len(main_reads) > 1:
        r2_id = main_reads[1]
    else:
        r2_id = None

    idx_reads = [x['read_num'] for x in read_config if x['is_index'] == True]
    if len(idx_reads) > 0:
        i1_id = idx_reads[0]
    else:
        i1_id = None

    if len(idx_reads) > 1:
        i2_id = idx_reads[1]
    else:
        i2_id = None

    read_names = {'read1':r1_id, 'read2':r2_id, 'index1':i1_id, 'index2':i2_id}

    # Qscore fractions by read
    try:
        qm = ds.QualityMetrics()
        q_levels = [ 20, 30 ]

        for q in q_levels:

            good_q_cols = [x for x in qm.df.keys() if x[0]=='q' and int(x[1:]) >= q]
            all_q_cols = [x for x in qm.df.keys() if x[0]=='q']

            # Group everything by cycle, and sum out tile, etc.
            grp = qm.df.groupby('cycle')
            grs = grp.sum()

            # Bases per cycle, all qualities
            qall = grs[all_q_cols].sum(axis=1)

            # Bases per cycle at requested quality
            qgood = grs[good_q_cols].sum(axis=1)

            fraction = qgood.values.astype('float') / qall

            for (read_name, rid) in read_names.items():
                tag = read_name + "_q" + str(q) + "_fraction"
                tag_by_cycle = tag + "_by_cycle"

                if rid is not None:
                    qfrac_by_cycle = fraction[start_cycle[rid-1]:end_cycle[rid-1]].values
                    results[tag] = qfrac_by_cycle.mean()
                    results[tag_by_cycle] = list(qfrac_by_cycle)
                else:
                    results[tag] = None
                    results[tag_by_cycle] = None
    except BaseException, e:
        print "Cannot read InterOp Quality Metrics: %s" % e.message
        # error case-- null out all the results
        q_levels = [ 20, 30 ]
        for q in q_levels:
            for (read_name, rid) in read_names.items():
                tag = read_name + "_q" + str(q) + "_fraction"
                tag_by_cycle = tag + "_by_cycle"
                results[tag] = None
                results[tag_by_cycle] = None

    try:
        em = ds.ErrorMetrics()
        err_by_cycle = em.df.groupby(em.df.cycle).mean()
        err_by_cycle.reset_index(inplace=True)

        # Get the PhiX error rate by cycle -
        for (read_name, rid) in read_names.items():
            tag = read_name + "_PhiX_error_by_cycle"

            if rid is None or read_config[rid-1]['is_index']:
                results[tag] = None
            else:
                cycle = err_by_cycle.cycle
                cycle_slice = np.logical_and(cycle >= start_cycle[rid-1], cycle < end_cycle[rid-1])
                results[tag] = list(err_by_cycle[cycle_slice].rate.values / 100)

        # PhiX error rate of worst tile
        worst_tile_phix = em.df.groupby(em.df.tile).mean().rate.max()
        results["PhiX_error_worst_tile"] = worst_tile_phix
    except BaseException, e:
        print "Cannot read InterOp Error Metrics: %s" % e.message
        # PhiX was not run -- fill in None
        for (read_name, rid) in read_names.items():
            tag = read_name + "_PhiX_error_by_cycle"
            results[tag] = None

        results["PhiX_error_worst_tile"] = None


    # Cluster summary stats
    try:
        tm = ds.TileMetrics()

        results['mean_cluster_density'] = tm.mean_cluster_density
        results['mean_cluster_density_pf'] = tm.mean_cluster_density_pf
        results['total_cluster_density'] = tm.total_cluster_density
        results['total_cluster_density_pf'] = tm.total_cluster_density_pf
        results['num_clusters'] = tm.num_clusters
        results['num_clusters_pf'] = tm.num_clusters_pf
        results['percent_pf_clusters'] = tm.percent_pf_clusters
        results['PhiX_aligned'] = tm.aligned

        # Split out (pre)phasing by read type
        for (read_name, rid) in read_names.items():
            tag = read_name + "_mean_phasing"
            pretag = read_name + "_mean_prephasing"

            # TileMetrics v3 lookahead: this is not available
            if rid is None or not tm.mean_phasing:
                results[tag] = None
                results[pretag] = None
            else:
                results[tag] = tm.mean_phasing[rid-1]
                results[pretag] = tm.mean_prephasing[rid-1]
    except BaseException, e:
        print "Cannot read InterOp Error Metrics: %s" % e.message
        # for some reason we don't have tile metrics -- blank out everything
        results['mean_cluster_density'] = None
        results['mean_cluster_density_pf'] = None
        results['total_cluster_density'] = None
        results['total_cluster_density_pf'] = None
        results['num_clusters'] = None
        results['percent_pf_clusters'] = None
        results['PhiX_aligned'] = None

        for (read_name, rid) in read_names.items():
            tag = read_name + "_mean_phasing"
            pretag = read_name + "_mean_prephasing"
            results[tag] = None
            results[pretag] = None

    # Corrected Intensity metrics
    # note: deprecated in HiSeq X / 4k
    try:
        cim = ds.CorrectedIntensityMetrics()
        if 'signoise_ratio' in cim.df:
            results['signoise_ratio'] = cim.df['signoise_ratio'].mean()
        else:
            results['signoise_ratio'] = None
    except BaseException, e:
        print "Cannot read InterOp CorrectedIntensity Metrics: %s" % e.message
        results['signoise_ratio'] = None

    try:
        exm = ds.ExtractionMetrics()
        ex_df = exm.df
        # looking ahead to ExtractionMetrics v4
        if hasattr(exm, 'channelcount'):
            channelcount = exm.channelcount
            for idx in range(channelcount):
                results['fwhm_%d' % idx] = ex_df['fwhm_%d' % idx].mean()
                results['intensity_%d' % idx] = ex_df['intensity_%d' % idx].mean()
        else:
            results['fwhm_A'] = ex_df['fwhm_A'].mean()
            results['fwhm_C'] = ex_df['fwhm_C'].mean()
            results['fwhm_G'] = ex_df['fwhm_G'].mean()
            results['fwhm_T'] = ex_df['fwhm_T'].mean()

            results['intensity_A'] = ex_df['intensity_A'].mean()
            results['intensity_C'] = ex_df['intensity_C'].mean()
            results['intensity_G'] = ex_df['intensity_G'].mean()
            results['intensity_T'] = ex_df['intensity_T'].mean()
    except BaseException, e:
        print "Cannot read InterOp Extraction Metrics: %s" % e.message
        results['fwhm_A'] = None
        results['fwhm_C'] = None
        results['fwhm_G'] = None
        results['fwhm_T'] = None
        results['intensity_A'] = None
        results['intensity_C'] = None
        results['intensity_G'] = None
        results['intensity_T'] = None

    results['yield'] = _compute_yield_raw(results)
    results['yield_pf'] = _compute_yield_pf(results)
    results['yield_pf_q30'] = _compute_yield_pf_q30(results)

    return results


def _compute_yield_raw(metrics):
    """
    Compute raw yield (number of bases total) from a populated
    metrics object generated by get_illumina_sequencing_metrics()
    """
    if not metrics.get("num_clusters"):
        return None

    total_length = 0
    total_length += len(metrics.get('read1_q30_fraction_by_cycle') or [])
    total_length += len(metrics.get('read2_q30_fraction_by_cycle') or [])
    total_length += len(metrics.get('index1_q30_fraction_by_cycle') or [])
    total_length += len(metrics.get('index2_q30_fraction_by_cycle') or [])

    return int(metrics.get('num_clusters')*total_length)


def _compute_yield_pf(metrics):
    """
    Compute the yield (number of bases passing filter) from a populated
    metrics object generated by get_illumina_sequencing_metrics()
    """
    if not metrics.get('num_clusters_pf'):
        return None

    total_length = 0
    total_length += len(metrics.get('read1_q30_fraction_by_cycle') or [])
    total_length += len(metrics.get('read2_q30_fraction_by_cycle') or [])
    total_length += len(metrics.get('index1_q30_fraction_by_cycle') or [])
    total_length += len(metrics.get('index2_q30_fraction_by_cycle') or [])

    return int(metrics.get('num_clusters_pf')*total_length)


def _compute_yield_pf_q30(metrics):
    """
    Compute the number of bases passing filter + Q30 from a populated
    metrics dictionary generated by get_illumina_sequencing_metrics()
    """
    if not metrics.get('num_clusters_pf'):
        return None

    num_q30_bases = 0
    num_clusters_pf = metrics.get('num_clusters_pf')
    for read_type in ('read1', 'read2', 'index1', 'index2'):
        if metrics.get('%s_q30_fraction' % read_type):
            num_q30_bases += metrics.get('%s_q30_fraction' % read_type)\
                             *len(metrics['%s_q30_fraction_by_cycle' % read_type])\
                             *num_clusters_pf
    return int(round(num_q30_bases))


def infer_symbolic_read_metrics(metrics, key_prefix, read_type, start_idx=0, length=None):
    """
    Take an existing metrics object and return a copy of the dictionary with metrics
    computed from the underlying raw read metrics.  Copies the per-cycle metrics
    for (length) bases starting from start_idx, and where applicable, computes
    aggregate statistics over those bases.

    Returns a new dictionary object with the additional metrics added.

    Phasing and prephasing metrics will not be copied if only a subset of the read length
    is selected.

    Example usage:

    Chromium Genome barcode metrics computation:
    metrics = infer_symbolic_read_metrics(metrics, 'barcode', QC_KEY_READ1, 0, 16)

    :param metrics: A metrics dictionary emitted by get_illumina_sequencing_metrics()
    :param key_prefix: The key prefix of the metrics to be generated
    :param read_type: The read type to draw from (must be a QC_KEY read value)
    :param start_idx: The index of the read (0-based) to start reading from.  Defaults
                      to the first read.
    :param end_idx: The number of bases to consider.  If None, consider all bases to the end.
    :return: The union of old and new metrics.
    """
    # if read_type is shorthand/file type, expand to full name from Illuminate
    if read_type in QC_SHORTCUT_MAP:
        read_type = QC_SHORTCUT_MAP[read_type]

    if "%s_q30_fraction" % read_type not in metrics:
        logging.warning("Read type Q30 not found in metrics, ignoring: %s" % read_type)
        return metrics
    elif not metrics.get("%s_q30_fraction_by_cycle" % read_type, None):
        logging.warning("Read type not run, ignoring: %s" % read_type)
        return metrics

    mean_phasing = metrics.get('%s_mean_phasing' % read_type, 0)
    mean_prephasing = metrics.get('%s_mean_prephasing' % read_type, 0)
    q30_fraction_by_cycle = metrics.get('%s_q30_fraction_by_cycle' % read_type, None)

    # if only a subset of bases are used, do not propagate phasing/prephasing averages
    copy_phasing_metrics = True
    if start_idx > 0 or q30_fraction_by_cycle is None or (length is not None and length < len(q30_fraction_by_cycle)):
        copy_phasing_metrics = False

    # default the last idx to use
    if length is None:
        end_idx = len(q30_fraction_by_cycle)
    else:
        end_idx = start_idx+length

    new_metrics = {}
    if copy_phasing_metrics:
        new_metrics['%s_mean_phasing' % key_prefix] = mean_phasing
        new_metrics['%s_mean_prephasing' % key_prefix] = mean_prephasing

    # copy cycles by interval
    for cycle_metric in ('PhiX_error_by_cycle', 'q20_fraction_by_cycle', 'q30_fraction_by_cycle'):
        orig_cycle_metrics = metrics.get('%s_%s' % (read_type, cycle_metric))
        new_metrics_key = '%s_%s' % (key_prefix, cycle_metric)
        if orig_cycle_metrics:
            new_metrics[new_metrics_key] = orig_cycle_metrics[start_idx:end_idx]
        else:
            new_metrics[new_metrics_key] = None

    # aggregate cycles where appropriate
    for mean_metric in ('q20_fraction', 'q30_fraction'):
        cycle_metrics = new_metrics['%s_%s_by_cycle' % (key_prefix, mean_metric)]
        if cycle_metrics:
            new_metrics["%s_%s" % (key_prefix, mean_metric)] = np.mean(cycle_metrics)
        else:
            new_metrics["%s_%s" % (key_prefix, mean_metric)] = None

    metrics_copy = dict(metrics)
    metrics_copy.update(new_metrics)
    return metrics_copy



