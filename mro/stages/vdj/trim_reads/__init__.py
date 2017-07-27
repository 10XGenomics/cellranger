#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# Trim adapter sequences
#
import itertools
import re
import os.path
import martian
import subprocess
import tenkit.seq as tk_seq
import cellranger.chemistry as cr_chem
import cellranger.report as cr_report
import cellranger.utils as cr_utils
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


def get_vdj_trim_metrics(reporter, stdout_file):
    with open(stdout_file, 'r') as f:
        for line in f:
            if not line.startswith('==='):
                continue
            if line.startswith('=== Summary'):
                f.next()
                npairs = int(re.sub(',', '', f.next().split()[4]))
            else:
                match = re.match(r'=== (\S+) read: Adapter (\S+) ===', line)
                read_name, seq_name = match.group(1, 2)

                if read_name in ('First', 'Second'):
                    read_name = 'read1' if read_name == 'First' else 'read2'
                else:
                    raise ValueError('Unrecognized read-type in cutadapt output: %s' % read_name)

                f.next()

                seq_info = f.next()
                ntrimmed = int(seq_info.split()[-2])
                metric_name = 'vdj_trim_{}_frac'.format(read_name)

                metric = reporter._get_metric_attr(metric_name, seq_name)
                metric.add(ntrimmed, filter=True)
                metric.add(npairs - ntrimmed, filter=False)


def run_cutadapt(args, out_read1s, out_read2s):
    paired_end = cr_chem.is_paired_end(args.chemistry_def)

    cmd = ['cutadapt',
           '--minimum-length', '50', '--times', '3', '--overlap', '5',
           '-o', out_read1s]

    if paired_end:
        cmd.extend(['-p', out_read2s])

    primers = {anno['name']:anno['seq'] for anno in args.primers}

    for name in R1_ANCHORED_FIVE_PRIME_SEQS:
        if name in primers:
            cmd.extend(['-g', '%s=^%s' % (name, primers[name])])

    for name in R1_THREE_PRIME_REV_COMP_SEQS:
        if name in primers:
            cmd.extend(['-a', '%s_rc=%s' % (name, tk_seq.get_rev_comp(primers[name]))])

    for name in R1_THREE_PRIME_SEQS:
        if name in primers:
            cmd.extend(['-a', '%s=%s' % (name, primers[name])])

    if paired_end:
        for name in R2_THREE_PRIME_REV_COMP_SEQS:
            if name in primers:
                cmd.extend(['-A', '%s_rc=%s' % (name, tk_seq.get_rev_comp(primers[name]))])

        for name in R2_THREE_PRIME_SEQS:
            if name in primers:
                cmd.extend(['-A', '%s=%s' % (name, primers[name])])

    cmd.extend([args.read1s_chunk])

    if paired_end:
        cmd.extend([args.read2s_chunk])

    status = subprocess.check_call(cmd)
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
    status = run_cutadapt(args, outs.read1s, outs.read2s)

    if args.read2s_chunk == None:
        outs.read2s = None

    if status != 0:
        martian.log_info('Error while running cutadapt')
    else:
        reporter = vdj_report.VdjReporter(primers=cr_utils.get_primers_from_dicts(args.primers))
        get_vdj_trim_metrics(reporter,
                             os.path.join(os.path.dirname(outs.chunked_reporter), '..', '_stdout'))
        reporter.save(outs.chunked_reporter)


def join(args, outs, chunk_defs, chunk_outs):
    outs.read1s = [chunk_out.read1s for chunk_out in chunk_outs]

    if all(chunk_out.read2s == None for chunk_out in chunk_outs):
        outs.read2s = []
    else:
        outs.read2s = [chunk_out.read2s for chunk_out in chunk_outs]

    reporters = [chunk_out.chunked_reporter for chunk_out in chunk_outs]
    if all([reporter is None for reporter in reporters]):
        cr_utils.write_empty_json(outs.summary)
    else:
        final_report = cr_report.merge_reporters(reporters)
        final_report.report_summary_json(outs.summary)
