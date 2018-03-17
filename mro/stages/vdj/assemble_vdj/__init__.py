#!/usr/bin/env python
#
# Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
#

import martian
import numpy as np
import itertools
import json
import os
import os.path
import sys
import tenkit.bam as tk_bam
import tenkit.log_subprocess as tk_subproc
import tenkit.safe_json as tk_safe_json
import cellranger.chemistry as cr_chem
import cellranger.utils as cr_utils

__MRO__ = """
stage ASSEMBLE_VDJ(
    in  map       chemistry_def,
    in  bam[]     barcode_chunked_bams   "reads split by barcode for chunking",
    in  int[]     chunk_gem_groups,
    in  tsv[]     reads_per_bc,
    in  int       mem_gb,
    in  int       min_kmer_count         "min number of UMIs in which a kmer occurs to consider it",
    in  int       min_contig_len         "min length of output sequences",
    in  int       min_qual               "min extension quality to consider a branch",
    in  float     score_factor           "assign UMIs to sequences with score_factor * best_umi_path_score",
    in  float     qual_factor            "consider branches with qual_factor * best extension quality",
    in  float     min_sw_score           "min SW score for read-contig alignments",
    in  int       min_readpairs_per_umi  "min readpairs/UMI (absolute count)",
    in  float     rt_error               "RT error for base qual computation",
    in  bool      use_unmapped           "Use unmapped reads in input",
    out fasta     contig_fasta,
    out fasta.fai contig_fasta_fai,
    out fastq     contig_fastq,
    out bam       contig_bam,
    out bam.bai   contig_bam_bai,
    out tsv       summary_tsv,
    out tsv       umi_summary_tsv,
    out json      metrics_summary_json,
    src py        "stages/vdj/assemble_vdj",
) split using (
    in  bam       chunked_bam,
    in  int       gem_group,
)
"""

# Empirical maxrss seems to scale sublinearly
# None of this is used if mem_gb is manually specified
MEM_BYTES_PER_READ = 15000
MAX_READS_PER_BC = 200000

# Merge output bams in chunks of N
MERGE_BAMS_N = 64

def run_assembly(fastq_pref, fasta_pref, args):
    cmd = [
        'vdj_asm', 'asm',
        fastq_pref,
        fasta_pref,
        '--kmers=' + str(args.min_kmer_count),
        '--min-contig=' + str(args.min_contig_len),
        '--min-qual=' + str(args.min_qual),
        '--score-factor=' + str(args.score_factor),
        '--qual-factor=' + str(args.qual_factor),
        '--min-sw-score=' + str(args.min_sw_score),
        '--rt-error=' + str(args.rt_error)]

    if not cr_chem.has_umis(args.chemistry_def):
        martian.log_info('Assembly without UMIs is not fully supported.')

    if cr_chem.is_paired_end(args.chemistry_def):
        cmd.append('--min-umi-reads=' + str(2 * args.min_readpairs_per_umi))
    else:
        cmd.append('--min-umi-reads=' + str(args.min_readpairs_per_umi))
        cmd.append('--single-end')

    if args.use_unmapped:
        cmd.append('--use-unmapped')

    #cmd.append('--mixture-filter')

    print >> sys.stderr, 'Running', ' '.join(cmd)
    tk_subproc.check_call(cmd, cwd=os.getcwd())


def split(args):
    chunks = []

    for reads_per_bc_file, bam, gem_group in itertools.izip(args.reads_per_bc,
                                                            args.barcode_chunked_bams,
                                                            args.chunk_gem_groups):
        subsample_rate = 1.0

        with open(reads_per_bc_file) as f:
            reads_per_bc = []
            for line in f:
                _, reads = line.strip().split()
                reads_per_bc.append(float(reads) * subsample_rate)

        max_reads = np.max(reads_per_bc + [0.0])

        # vdj_asm is hard-coded to use a maximum of 200k reads / BC.
        max_reads = min(MAX_READS_PER_BC, max_reads)

        # The assembly step takes roughly num_reads * MEM_BYTES_PER_READ bytes of memory to complete each BC.
        if args.mem_gb is None:
            mem_gb = max(2, int(np.ceil(MEM_BYTES_PER_READ * max_reads / 1e9)))
        else:
            mem_gb = args.mem_gb

        chunks.append({
            'chunked_bam': bam,
            'gem_group': gem_group,
            '__mem_gb': mem_gb,
        })

    if not chunks:
        # No input reads
        return {'chunks': []}
    return {'chunks': chunks, 'join': {'__threads': 4, '__mem_gb': 16}}


def main(args, outs):
    run_assembly(args.chunked_bam, martian.make_path(''), args)

    out_pref = os.path.splitext(os.path.basename(args.chunked_bam))[0]
    out_pref = martian.make_path(out_pref)
    cr_utils.move(out_pref + '.fasta', outs.contig_fasta)
    cr_utils.move(out_pref + '.fastq', outs.contig_fastq)
    cr_utils.move(out_pref + '_summary.tsv', outs.summary_tsv)
    cr_utils.move(out_pref + '_umi_summary.tsv', outs.umi_summary_tsv)
    cr_utils.move(out_pref + '_sorted.bam', outs.contig_bam)
    cr_utils.move(out_pref + '_sorted.bam.bai', outs.contig_bam_bai)
    cr_utils.move(out_pref + '_metrics_summary.json', outs.metrics_summary_json)


def join(args, outs, chunk_defs, chunk_outs):
    contigs = []
    contig_fastqs = []
    contig_bams = []

    if len(chunk_outs) == 0:
        # No input reads
        # Create empty BAM file
        with open(outs.contig_bam, 'w') as f:
            pass
        outs.contig_bam_bai = None
        # Create empty contig FASTA
        with open(outs.contig_fasta, 'w') as f:
            pass
        outs.contig_fasta_fai = None
        # Create empty contig FASTQ
        with open(outs.contig_fastq, 'w') as f:
            pass
        outs.metrics_summary_json = None
        outs.summary_tsv = None
        outs.umi_summary_tsv = None
        return

    summary_tsvs = []
    umi_summary_tsvs = []

    for chunk_out in chunk_outs:
        if not os.path.isfile(chunk_out.contig_fasta):
            continue
        contigs.append(chunk_out.contig_fasta)

        contig_fastqs.append(chunk_out.contig_fastq)
        contig_bams.append(chunk_out.contig_bam)

        summary_tsvs.append(chunk_out.summary_tsv)
        umi_summary_tsvs.append(chunk_out.umi_summary_tsv)

    cr_utils.concatenate_files(outs.contig_fasta, contigs)

    if os.path.getsize(outs.contig_fasta) > 0:
        tk_subproc.check_call('samtools faidx %s' % outs.contig_fasta, shell=True)
        outs.contig_fasta_fai = outs.contig_fasta + '.fai'

    cr_utils.concatenate_files(outs.contig_fastq, contig_fastqs)

    if len(summary_tsvs) > 0:
        cr_utils.concatenate_headered_files(outs.summary_tsv, summary_tsvs)
    if len(umi_summary_tsvs) > 0:
        cr_utils.concatenate_headered_files(outs.umi_summary_tsv, umi_summary_tsvs)

    if contig_bams:
        # Merge every N BAMs. Trying to merge them all at once
        #  risks hitting the filehandle limit.
        n_merged = 0

        while len(contig_bams) > 1:
            to_merge = contig_bams[0:MERGE_BAMS_N]

            tmp_bam = martian.make_path('merged-%04d.bam' % n_merged)
            n_merged += 1

            print "Merging %d BAMs into %s ..." % (len(to_merge), tmp_bam)
            tk_bam.merge(tmp_bam, to_merge, threads=args.__threads)

            # Delete any temporary bams that have been merged
            for in_bam in to_merge:
                if os.path.basename(in_bam).startswith('merged-'):
                    cr_utils.remove(in_bam)

            # Pop the input bams and push the merged bam
            contig_bams = contig_bams[len(to_merge):] + [tmp_bam]

        if os.path.basename(contig_bams[0]).startswith('merged-'):
            # We merged at least two chunks together.
            # Rename it to the output bam.
            cr_utils.move(contig_bams[0], outs.contig_bam)
        else:
            # There was only a single chunk, so copy it from the input
            cr_utils.copy(contig_bams[0], outs.contig_bam)

        tk_bam.index(outs.contig_bam)

        # Make sure the Martian out matches the actual index filename
        outs.contig_bam_bai = outs.contig_bam + '.bai'

    # Merge the assembler summary jsons
    merged_summary = cr_utils.merge_jsons_single_level([out.metrics_summary_json for out in chunk_outs])

    with open(outs.metrics_summary_json, 'w') as f:
        json.dump(tk_safe_json.json_sanitize(merged_summary), f, indent=4, sort_keys=True)
