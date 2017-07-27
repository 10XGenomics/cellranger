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
import pandas as pd
import subprocess
import sys
import tenkit.bam as tk_bam
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
    in  float     nx                     "min number of reads/UMI (90.0 means take N90)",
    in  int       npaths                 "num of paths per graph component to consider",
    in  float     score_factor           "assign UMIs to sequences with score_factor * best_umi_path_score",
    in  float     qual_factor            "consider branches with qual_factor * best extension quality",
    in  bool      use_sw                 "align reads on contigs using SW or fast k-mer extension",
    in  float     min_sw_score           "min SW score for read-contig alignments",
    in  map       min_readpairs_per_umi  "per-gem-group min readpairs/UMI (absolute count, unlike nx above)",
    in  map       subsample_rate         "per-gem-group read subsampling rate",
    in  float     rt_error               "RT error for base qual computation",
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
MEM_BYTES_PER_READ = 15000
MAX_READS_PER_BC = 200000

def run_assembly(fastq_pref, fasta_pref, args):
    cmd = [
        'vdj_asm', 'asm',
        fastq_pref,
        fasta_pref,
        '--kmers=' + str(args.min_kmer_count),
        '--min-contig=' + str(args.min_contig_len),
        '--npaths=' + str(args.npaths),
        '--nx=' + str(args.nx),
        '--min-qual=' + str(args.min_qual),
        '--score-factor=' + str(args.score_factor),
        '--qual-factor=' + str(args.qual_factor),
        '--min-sw-score=' + str(args.min_sw_score),
        '--rt-error=' + str(args.rt_error),
        '--subsample-rate=' + str(args.subsample_rate[str(args.gem_group)]),
    ]
    if not cr_chem.has_umis(args.chemistry_def):
        martian.log_info('Assembly without UMIs is not fully supported.')
    if not args.use_sw:
        cmd.append('--fast-align')
    if not args.min_readpairs_per_umi is None:
        # If only assembling with read2, adjust this cutoff
        # NOTE: Martian stores the gem_group dict keys as strings
        cutoff = args.min_readpairs_per_umi[str(args.gem_group)]
        cmd.append('--min-umi-reads=' + str(cutoff))
    print >> sys.stderr, 'Running', ' '.join(cmd)
    subprocess.check_call(cmd, cwd=os.getcwd())


def split(args):
    chunks = []

    for reads_per_bc_file, bam, gem_group in itertools.izip(args.reads_per_bc,
                                                            args.barcode_chunked_bams,
                                                            args.chunk_gem_groups):
        subsample_rate = args.subsample_rate[str(gem_group)]

        with open(reads_per_bc_file) as f:
            reads_per_bc = []
            for line in f:
                _, reads = line.strip().split()
                reads_per_bc.append(float(reads) * subsample_rate)

        max_reads = np.max(reads_per_bc + [0.0])

        # vdj_asm is hard-coded to use a maximum of 200k reads / BC.
        max_reads = min(MAX_READS_PER_BC, max_reads)

        # The assembly step takes roughly num_reads * MEM_BYTES_PER_READ bytes of memory to complete each BC.
        mem_gb = max(2.0, int(np.ceil(MEM_BYTES_PER_READ * max_reads / 1e9)))

        chunks.append({
            'chunked_bam': bam,
            'gem_group': gem_group,
            '__mem_gb': mem_gb,
        })

    # If there were no input reads, create a dummy chunk
    if not chunks:
        chunks.append({'chunked_bam': None})
    return {'chunks': chunks, 'join': {'__threads': 4}}


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

    summary_df_parts = []
    umi_summary_df_parts = []

    for chunk_out in chunk_outs:
        if not os.path.isfile(chunk_out.contig_fasta):
            continue
        contigs.append(chunk_out.contig_fasta)

        contig_fastqs.append(chunk_out.contig_fastq)
        contig_bams.append(chunk_out.contig_bam)
        summary_df_parts.append(pd.read_csv(chunk_out.summary_tsv,
                                     header=0, index_col=None, sep='\t',
                                     dtype={'component': int, 'num_reads': int,
                                            'num_pairs': int, 'num_umis': int}))

        umi_summary_df_parts.append(pd.read_csv(chunk_out.umi_summary_tsv,
                                     header=0, index_col=None, sep='\t',
                                     dtype={'umi_id': int, 'reads': int,
                                            'min_umi_reads': int, 'contigs': str}))

    summary_df = pd.concat(summary_df_parts, ignore_index=True)
    umi_summary_df = pd.concat(umi_summary_df_parts, ignore_index=True)

    cr_utils.concatenate_files(outs.contig_fasta, contigs)

    if os.path.getsize(outs.contig_fasta) > 0:
        subprocess.check_call('samtools faidx %s' % outs.contig_fasta, shell=True)
        outs.contig_fasta_fai = outs.contig_fasta + '.fai'

    cr_utils.concatenate_files(outs.contig_fastq, contig_fastqs)

    if summary_df is not None:
        summary_df.to_csv(outs.summary_tsv, header=True, index=False, sep='\t')
    if umi_summary_df is not None:
        umi_summary_df.to_csv(outs.umi_summary_tsv, header=True, index=False, sep='\t')

    if contig_bams:
        tk_bam.merge(outs.contig_bam, contig_bams, threads=args.__threads)
        tk_bam.index(outs.contig_bam)

        # Make sure the Martian out matches the actual index filename
        outs.contig_bam_bai = outs.contig_bam + '.bai'

    # Merge the assembler summary jsons
    merged_summary = cr_utils.merge_jsons_single_level([out.metrics_summary_json for out in chunk_outs])

    with open(outs.metrics_summary_json, 'w') as f:
        json.dump(tk_safe_json.json_sanitize(merged_summary), f, indent=4, sort_keys=True)
