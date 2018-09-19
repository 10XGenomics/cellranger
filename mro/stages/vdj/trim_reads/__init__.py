#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# Trim adapter sequences
#
import itertools
import re
import os.path
import sys
import martian
import tenkit.log_subprocess as tk_subproc
import tenkit.seq as tk_seq
import cellranger.h5_constants as h5_constants
import cellranger.chemistry as cr_chem
import cellranger.report as cr_report
import cellranger.utils as cr_utils
import cellranger.io as cr_io
import cellranger.vdj.report as vdj_report

__MRO__ = """
stage TRIM_READS(
    in  fastq[] read1s,
    in  fastq[] read2s,
    in  map[]   primers,
    in  map     chemistry_def,
    out fastq[] read1s,
    out fastq[] read2s,
    out json    summary,
    out pickle  chunked_reporter,
    src py      "stages/vdj/trim_reads",
) split using (
    in  fastq   read1s_chunk,
    in  fastq   read2s_chunk,
)
"""

# Names of "primer" sequences expected to be anchored at the 5' end of R1
R1_ANCHORED_FIVE_PRIME_SEQS = ['spacer']

# Names of "primer" sequences expected at the 3' end of R1 and reverse-complemented
R1_THREE_PRIME_REV_COMP_SEQS = ['R2', 'P7']

# Names of "primer" sequences expected at the 3' end of R1
R1_THREE_PRIME_SEQS = ['polyA', 'rt_primer']

# Names of "primer" sequences expected at the 3' end of R2 and reverse-complemented
R2_THREE_PRIME_REV_COMP_SEQS = ['spacer', 'R1', 'P5']

# Names of "primer" sequences expected at the 3' end of R2
R2_THREE_PRIME_SEQS = ['polyA', 'rt_primer']


def get_vdj_trim_metrics(reporter, stdout_file, paired_end):
    with open(stdout_file, 'r') as f:
        for line in f:
            if not line.startswith('==='):
                continue
            if line.startswith('=== Summary'):
                f.next()

                if paired_end:
                    npairs = int(re.sub(',', '', f.next().split()[4]))
                else:
                    npairs = int(re.sub(',', '', f.next().split()[3]))
            else:
                if paired_end:
                    match = re.match(r'=== (\S+) read: Adapter (\S+) ===', line)
                    read_name, seq_name = match.group(1, 2)
                else:
                    match = re.match(r'=== Adapter (\S+) ===', line)
                    read_name, seq_name = (None, match.group(1))

                if read_name in ('First', 'Second'):
                    read_name = 'read1' if read_name == 'First' else 'read2'
                elif not paired_end:
                    read_name = 'read1'
                else:
                    raise ValueError('Unrecognized read-type in cutadapt output: %s' % read_name)

                f.next()

                seq_info = f.next()
                ntrimmed = int(seq_info.split()[-2])
                metric_name = 'vdj_trim_{}_frac'.format(read_name)

                metric = reporter._get_metric_attr(metric_name, seq_name)
                metric.add(ntrimmed, filter=True)
                metric.add(npairs - ntrimmed, filter=False)


def run_cutadapt(args, out_read1s, out_read2s, chemistry_def, stdout=sys.stdout):
    paired_end = cr_chem.is_paired_end(chemistry_def)

    # If single end, determine which read the single read is (R1 or R2)
    if paired_end:
        single_read = None
    else:
        single_read = cr_chem.get_rna_read_def(chemistry_def).read_type
        assert single_read in ('R1', 'R2')

    out_r1_file = cr_io.open_maybe_gzip_subprocess(out_read1s, 'w')

    # Note: The complexity of forcing cutadapt to output a compressed file
    #       means we'll have to give up on that for now.
    cmd = ['cutadapt',
           '-e', '0.12', '--times', '3', '--overlap', '5',
           '-f', 'fastq',
           '-o', '/proc/%d/fd/%d' % (os.getpid(), out_r1_file.fileno())]

    out_r2_file = None
    if paired_end:
        out_r2_file = cr_io.open_maybe_gzip_subprocess(out_read2s, 'w')
        cmd.extend(['-p', '/proc/%d/fd/%d' % (os.getpid(), out_r2_file.fileno())])

    primers = {anno['name']:anno['seq'] for anno in args.primers}

    if paired_end or single_read == 'R1':
        # R1 adapters
        for name in R1_ANCHORED_FIVE_PRIME_SEQS:
            if name in primers:
                cmd.extend(['-g', '%s=^%s' % (name, primers[name])])

        for name in R1_THREE_PRIME_REV_COMP_SEQS:
            if name in primers:
                cmd.extend(['-a', '%s_rc=%s' % (name, tk_seq.get_rev_comp(primers[name]))])

        for name in R1_THREE_PRIME_SEQS:
            if name in primers:
                cmd.extend(['-a', '%s=%s' % (name, primers[name])])

    if paired_end or single_read == 'R2':
        for name in R2_THREE_PRIME_REV_COMP_SEQS:
            if name in primers:
                flag = '-A' if paired_end else '-a'
                cmd.extend([flag, '%s_rc=%s' % (name, tk_seq.get_rev_comp(primers[name]))])

        for name in R2_THREE_PRIME_SEQS:
            if name in primers:
                flag = '-A' if paired_end else '-a'
                cmd.extend([flag, '%s=%s' % (name, primers[name])])


    read1_file = cr_io.open_maybe_gzip_subprocess(args.read1s_chunk)
    cmd.extend(['/proc/%d/fd/%d' % (os.getpid(), read1_file.fileno())])

    read2_file = None
    if paired_end:
        read2_file = cr_io.open_maybe_gzip_subprocess(args.read2s_chunk)
        cmd.extend(['/proc/%d/fd/%d' % (os.getpid(), read2_file.fileno())])

    print cmd

    status = tk_subproc.check_call(cmd, stdout=stdout)

    # closing these files is important both because we need to wait on the
    # subprocess, if any, or else its rusage isn't accounted for for this
    # process, and because if we don't have a reference to the objects down
    # here, then python's garbage collector is free to finalize the objects
    # before cmd runs, which would result in a failure.
    out_r1_file.close()
    if out_r2_file:
        out_r2_file.close()
    read1_file.close()
    if read2_file:
        read2_file.close()

    return status


def split(args):
    chunks = []

    if cr_chem.is_paired_end(args.chemistry_def):
        assert len(args.read1s) == len(args.read2s)

    for fastq1, fastq2 in itertools.izip_longest(args.read1s, args.read2s):
        chunks.append({
            'read1s_chunk': fastq1,
            'read2s_chunk': fastq2,
        })

    return {'chunks': chunks}


def main(args, outs):
    paired_end = cr_chem.is_paired_end(args.chemistry_def)

    # Write compressed files
    outs.read1s += h5_constants.LZ4_SUFFIX
    outs.read2s += h5_constants.LZ4_SUFFIX

    cutadapt_out = os.path.join(os.path.dirname(outs.chunked_reporter), 'cutadapt_stdout')
    with open(cutadapt_out, 'w') as cut_stdout:
        status = run_cutadapt(args,
                              outs.read1s, outs.read2s,
                              args.chemistry_def,
                              cut_stdout)

    if args.read2s_chunk == None:
        outs.read2s = None

    if status != 0:
        martian.log_info('Error while running cutadapt')
    else:
        reporter = vdj_report.VdjReporter(primers=cr_utils.get_primers_from_dicts(args.primers))
        get_vdj_trim_metrics(reporter,
                             cutadapt_out,
                             paired_end)
        reporter.save(outs.chunked_reporter)


def join(args, outs, chunk_defs, chunk_outs):
    outs.read1s = [chunk_out.read1s for chunk_out in chunk_outs]

    if all(chunk_out.read2s == None for chunk_out in chunk_outs):
        outs.read2s = []
    else:
        outs.read2s = [chunk_out.read2s for chunk_out in chunk_outs]

    reporters = [chunk_out.chunked_reporter for chunk_out in chunk_outs]
    if all([reporter is None for reporter in reporters]):
        cr_io.write_empty_json(outs.summary)
    else:
        final_report = cr_report.merge_reporters(reporters)
        final_report.report_summary_json(outs.summary)
