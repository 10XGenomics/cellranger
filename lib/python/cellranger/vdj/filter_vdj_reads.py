#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#
# Utility functions for vdj read filtering

import itertools
import os
import sys
import tenkit.fasta as tk_fasta
import tenkit.log_subprocess as tk_subproc
import tenkit.seq as tk_seq

def write_bam_read_fastq(out, read):
    if read.is_reverse:
        seq, qual = tk_seq.get_rev_comp(read.seq), read.qual[::-1]
    else:
        seq, qual = read.seq, read.qual
    tk_fasta.write_read_fastq(out, read.qname, seq, qual)

def run_read_match(read1_path, read2_path, fasta_path, out_bam_filename, strand, sw_params):
    assert strand in ('+', '-')
    cmd = ['vdj_asm', 'read-match',
           '--ref', fasta_path,
           '--r1', read1_path,
           '--outbam', out_bam_filename,
           '--seed=' + str(sw_params['seed']),
           '--min-sw-score=' + str(sw_params['min_sw_score'])]

    if strand == '-':
        cmd.append('--rev-strand')

    if read2_path:
        cmd.extend(['--r2', read2_path])

    print >> sys.stderr, 'Running', ' '.join(cmd)
    tk_subproc.check_call(cmd, cwd=os.getcwd())

def get_pair_iter(iterable):
    """ Return (x_i, x_(i+1)) for i in {0,2,4,...} """
    return itertools.izip(iterable, iterable)
