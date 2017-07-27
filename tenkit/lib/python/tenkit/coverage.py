#!/usr/bin/env python
#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#
# Measure coverage in give regions
#
import tenkit.pandas as pd
import numpy as np
import tenkit.bio_io as tk_io
import tenkit.bam as tk_bam
import tenkit.stats as tk_stats
import martian
import tenkit.hdf5 as tk_hdf5

def mean_coverage_region(bam, region, read_filter=lambda x: True):
    ''' Measure the coverage mean in a region '''
    (chrom, start, end) = region
    reads_iter = tk_bam.filter_bam(bam.fetch(chrom, start, end), remove_unmapped=True, read_filter=read_filter)

    depth_df = get_depth_info(reads_iter, chrom, start, end)
    return depth_df.coverage.mean()


def get_depth_info(read_iter, chrom, cstart, cend):

    depths = np.zeros(cend-cstart, np.int32)

    for read in read_iter:
        pos = read.pos
        rstart = max(pos, cstart)

        # Increment to the end of the window or the end of the
        # alignment, whichever comes first
        rend = min(read.aend, cend)
        depths[(rstart-cstart):(rend-cstart)] += 1

    positions = np.arange(cstart, cend, dtype=np.int32)

    depth_df = pd.DataFrame({"chrom": chrom, "pos": positions, "coverage": depths})
    return depth_df


def get_depth_info_json(info):
    fixed_info = {int(x): y for (x, y) in info.iteritems()}

    total_depth_counts = sum(fixed_info.values())
    median_depth = None
    sorted_depths = sorted(fixed_info.keys())
    seen_depth_count = 0
    mean_depth = 0.0
    for depth in sorted_depths:
        seen_depth_count += fixed_info[depth]
        mean_depth += float(depth*fixed_info[depth])/float(total_depth_counts)
        if seen_depth_count > total_depth_counts/2 and median_depth is None:
            median_depth = depth
    zero_cov_fract = tk_stats.robust_divide(float(fixed_info.get(0, 0.0)), float(total_depth_counts))

    return (mean_depth, median_depth, zero_cov_fract)


EXONS_SAMPLE_COVERAGE = 3000
WGS_WINDOWS_SAMPLE_COVERAGE = 3000
WGS_WINDOW_SIZE = 10000
WGS_WINDOWS_SMALL_GENOME = 5

def allow_all(x):
    return True

def estimate_mean_coverage(targets_file, bam_in, read_filter=lambda x: True):

    if targets_file is not None:
        target_regions_dict = tk_io.get_target_regions_dict(open(targets_file))

        # Pick a random sample of target regions to estimate overall depth on
        targets = [(chrom, start, end) for (chrom, regions) in target_regions_dict.items() for (start, end) in regions if end-start > 0]
        if len(targets) == 0:
            martian.log_info("No non-empty target regions")
            return 1.0

        np.random.seed(0)
        regions_to_sample = min(EXONS_SAMPLE_COVERAGE, len(targets))
        region_indices = np.random.choice(len(targets), regions_to_sample, replace=False)
        sample_targets = [targets[idx] for idx in region_indices]

    else:
        # Pick a series of random intervals on the genome to measure coverage
        np.random.seed(0)

        if sum(bam_in.lengths) < 1e6:
            num_windows = WGS_WINDOWS_SMALL_GENOME
        else:
            num_windows = WGS_WINDOWS_SAMPLE_COVERAGE

        chrom_probs = np.array(bam_in.lengths, dtype=np.float) / sum(bam_in.lengths)
        rand_chroms = np.random.choice(len(bam_in.lengths), num_windows, replace=True, p=chrom_probs)

        starts = [np.random.randint(max(bam_in.lengths[chrom]-WGS_WINDOW_SIZE, 1)) for chrom in rand_chroms]
        sample_targets = [(bam_in.references[chrom], start, min(start+WGS_WINDOW_SIZE, bam_in.lengths[chrom])) for (chrom, start) in zip(rand_chroms, starts)]

    mean_depth = float(np.mean([mean_coverage_region(bam_in, region, read_filter) for region in sample_targets]))
    return mean_depth


def get_hap_coverage(in_bam, ps_h5, chrom, start, stop, cov_quals):
    """Return a dataframe with coverage per haplotype.

    Args:
    - in_bam: reader for a position sorted bam
    - ps_h5: HDF5 with phase set coordinates
    - chrom, start, stop: region to get coverage
    - cov_quals: Array of MAPQ cutoffs.

    Return value:
    A dataframe with columns:
    - chrom
    - pos
    - cov_q<M>_hap<H> for all M in cov_quals and for H in [0, 1, 2]: This is the
    coverage on haplotype H using reads of MAPQ >= M. Haplotype 2 corresponds to
    unphased.
    - phase_set: null if ps_h5 is missing.
    """
    coverages = [np.zeros((stop - start, 3)) for _ in cov_quals]

    for _, read in enumerate(in_bam.fetch(str(chrom), int(start), int(stop))):
        if not read.is_unmapped and not read.aend is None and not read.is_secondary and not read.is_duplicate:
            hap = tk_io.get_read_haplotype(read)
            hap_idx = 2 if hap is None else hap - 1
            range_start = max(0, read.pos - start)
            range_stop = min(stop, read.aend) - start
            for qi, q in enumerate(cov_quals):
                if read.mapq >= q:
                    coverages[qi][range_start:range_stop + 1, hap_idx] += 1

    base_df = pd.DataFrame({'chrom':chrom, 'pos':np.arange(start, stop)})
    dfs = map(lambda x: pd.DataFrame(x[0], columns=['cov_q' + str(x[1]) + '_hap' + str(i) for i in range(3)]),
              zip(coverages, cov_quals))
    df = pd.concat([base_df, pd.concat(dfs, axis=1)], axis=1)

    phase_sets = -np.ones((stop - start, ), dtype=np.int)

    # This can be None if for example the input is unbarcoded.
    if not ps_h5 is None:
        ps_df = tk_hdf5.read_data_frame(ps_h5)
        ps_df = ps_df[np.logical_and(ps_df.chrom == chrom, np.logical_and(ps_df.end >= start, ps_df.start < stop))]

        for _, row in ps_df.iterrows():
            range_start = max(0, row.start - start)
            range_stop = min(stop, row.end) - start
            phase_sets[range_start:range_stop + 1] = row.phase_set

    df['phase_set'] = phase_sets
    return df
