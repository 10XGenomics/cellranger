#!/usr/bin/env python
#
# Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
#
# Builds consensus sequences for the detected clonotypes based on the sequences
# of the member cells and computes the base qualities of the consensus sequences.

from collections import defaultdict, Counter
import itertools
import json
import martian
import numpy as np
from numpy.random import choice, seed
import os
import pysam
import re
import shutil
import sys
import tenkit.bam as tk_bam
import tenkit.fasta as tk_fasta
import tenkit.log_subprocess as tk_subproc
import tenkit.seq as tk_seq
from cellranger.constants import PROCESSED_UMI_TAG, PROCESSED_BARCODE_TAG, MIN_MEM_GB
import cellranger.fastq as cr_fastq
import cellranger.report as cr_report
import cellranger.utils as cr_utils
import cellranger.vdj.annotations as vdj_annot
from cellranger.vdj.constants import VDJ_V_FEATURE_TYPES, VDJ_J_FEATURE_TYPES, VDJ_ANNOTATION_MIN_SCORE_RATIO
import cellranger.vdj.report as vdj_report
import cellranger.vdj.utils as vdj_utils

__MRO__ = """
stage ASSEMBLE_CONSENSUS(
    in  path       vdj_reference_path,
    in  map[]      primers,
    in  string     metric_prefix                "to prefix summary metrics (eg. if you're running for raw and inferred clonotypes)",
    in  json       annotations_json,
    in  json       clonotype_assignments,
    in  tsv        umi_summary_tsv,
    in  fastq      contigs_fastq,
    in  bam        contig_bam,
    in  map        min_score_ratios             "dict of min score ratios by feature type to use for filtering annotations",
    in  map        min_word_sizes               "dict of min word sizes by feature type to use for filtering annotations",
    out pickle     chunked_reporter,
    out bam[]      chunked_consensus_bams,
    out bam[]      chunked_concat_ref_bams,
    out bam        consensus_bam,
    out bam.bai    consensus_bam_bai,
    out bam        concat_ref_bam,
    out bam.bai    concat_ref_bam_bai,
    out json       consensus_annotations_json,
    out json       concat_ref_annotations_json,
    out fasta      consensus_fasta,
    out fastq      consensus_fastq,
    out fasta.fai  consensus_fasta_fai,
    out fasta      concat_ref_fasta,
    out fasta.fai  concat_ref_fasta_fai,
    out csv        consensus_annotations_csv,
    out json       summary,
    out csv        clonotypes                   "info about clonotypes",
    src py         "stages/vdj/assemble_consensus",
) split using (
    in  string[]   chunk_clonotypes,
    in  json       chunk_annotations,
)
"""


MAPPED_UNPAIRED_FLAG = int('0b1000000', 2)

MAX_CHUNKS = 200
MIN_CLONOTYPES_PER_CHUNK = 1
# Maximum number of cells to use for computing consensus base qualities.
# Using more than a few generally doesn't increase our confidence.
MAX_CELLS_FOR_BASE_QUALS = 10

# Cap the number of reads per UMI and per contig to use for base qualities
MAX_READS_PER_UMI = 5000
MAX_READS_PER_CONTIG = 20000

# If fewer contigs than this, we'll just pick one as the consensus.
MIN_CONTIGS_FOR_CONSENSUS = 2

# Merge consensus BAMs every N files
MERGE_BAMS_EVERY = 10

def rm_files(filenames):
    for filename in filenames:
        cr_utils.remove(filename, allow_nonexisting=True)

def split(args):
    seed(1)

    if not args.clonotype_assignments:
        return {'chunks':[{'chunk_clonotypes':[]}]}

    # Distribute clonotypes to chunks. The input file might be sorted by frequency
    # so make sure you shuffle the clonotypes in order to better distribute load.
    with open(args.clonotype_assignments) as f:
        clonotypes = json.load(f)
    num_clonotypes = len(clonotypes)

    clonotypes_per_chunk = max(MIN_CLONOTYPES_PER_CHUNK, np.ceil(float(num_clonotypes) / MAX_CHUNKS))
    num_chunks = int(np.ceil(float(num_clonotypes) / clonotypes_per_chunk))

    chunk_clonotypes = [[] for _ in xrange(num_chunks)]

    chunks = []

    full_annot_mem = vdj_utils.get_mem_gb_from_annotations_json(args.annotations_json)
    join_mem_gb = max(16, full_annot_mem)

    if num_clonotypes > 0:
        # Pick a chunk 0..num_chunks num_clonotypes times.
        chunk_assignments = choice(np.arange(num_chunks), num_clonotypes, replace=True)
        for idx, clonotype_id in enumerate(clonotypes.iterkeys()):
            chunk_clonotypes[chunk_assignments[idx]].append(clonotype_id)

        # Remove empty chunks
        chunk_clonotypes = [x for x in chunk_clonotypes if len(x) > 0]

        # Map clonotype to chunk
        clo2chunk = {}
        for chunk_idx, clo_ids in enumerate(chunk_clonotypes):
            for clo_id in clo_ids:
                clo2chunk[clo_id] = chunk_idx

        ### Split the contig annotation json into chunks
        chunk_annot_fns = [martian.make_path('annot-%04d.json' % i) for i,_ in enumerate(chunk_clonotypes)]

        clo_key = '%s_clonotype_id' % args.metric_prefix
        cons_key = '%s_consensus_id' % args.metric_prefix
        with open(args.annotations_json) as input_fh, \
             vdj_utils.CachedJsonDictListWriters(chunk_annot_fns) as writers:
            for contig_dict in vdj_utils.get_json_obj_iter(input_fh):
                clo_id = contig_dict.get('info', {}).get(clo_key)
                cons_id = contig_dict.get('info', {}).get(cons_key)
                if clo_id is not None and cons_id is not None:
                    writers.write(contig_dict, clo2chunk[clo_id])


        for clo_list, annot_fn in itertools.izip(chunk_clonotypes, chunk_annot_fns):
            chunk_mem_gb = max(MIN_MEM_GB, vdj_utils.get_mem_gb_from_annotations_json(annot_fn))

            chunks.append({
                'chunk_clonotypes': clo_list,
                'chunk_annotations': annot_fn,
                '__mem_gb': chunk_mem_gb,
            })

    if len(chunks) == 0:
        return {'chunks': []}

    return {'chunks': chunks, 'join': {'__mem_gb': join_mem_gb}}


def main(args, outs):
    outs.chunked_consensus_bams = []
    outs.chunked_concat_ref_bams = []

    chunk_clonotypes = set(args.chunk_clonotypes)

    reporter = vdj_report.VdjReporter()
    if not args.clonotype_assignments or not vdj_utils.bam_has_seqs(args.contig_bam):
        # always produce an empty summary
        reporter.save(outs.chunked_reporter)
        return

    # Get the clonotype-barcode assignments
    with open(args.clonotype_assignments) as f:
        clonotypes = json.load(f)

    # Partition contig annotations by consensus id
    consensus_to_contigs = defaultdict(list)
    relevant_contig_ids = set()

    with open(args.chunk_annotations) as f:
        contigs = vdj_annot.load_contig_list_from_json(f, args.vdj_reference_path)

    clo_key = '%s_clonotype_id' % args.metric_prefix
    cons_key = '%s_consensus_id' % args.metric_prefix

    for contig in contigs:
        clo_id = contig.info_dict.get(clo_key)
        cons_id = contig.info_dict.get(cons_key)
        assert clo_id in chunk_clonotypes and cons_id is not None

        consensus_to_contigs[cons_id].append(contig)
        relevant_contig_ids.add(contig.contig_name)

    assert len(consensus_to_contigs) > 0

    in_bam = tk_bam.create_bam_infile(args.contig_bam)

    n_merged_bams = 0

    # For all contigs relevant to this chunk,
    #   get the assembler umi data required for base qual recalculation.
    # Do not attempt to read into a pandas object because it can be huge.
    contig_umis = defaultdict(set)
    with open(args.umi_summary_tsv, 'r') as umi_file:
        for line in umi_file:
            fields = line.strip().split('\t')
            umi = fields[2]
            if umi == 'umi' or len(fields) < 7:
                continue
            good_umi = fields[5].lower() == 'true'
            contig_ids = set(fields[6].split(','))
            if good_umi and len(contig_ids & relevant_contig_ids) > 0:
                for c in contig_ids:
                    contig_umis[c].add(umi)

    consensus_fastq = open(outs.consensus_fastq, 'w')
    consensus_fasta = open(outs.consensus_fasta, 'w')
    ref_fasta = open(outs.concat_ref_fasta, 'w')

    consensus_contigs = []
    ref_contigs = []

    assert(args.metric_prefix in reporter.vdj_clonotype_types)

    # Iterate over clonotype assignments
    for clonotype_id, clonotype in clonotypes.iteritems():
        if not clonotype_id in chunk_clonotypes:
            continue

        for consensus_id, consensus in clonotype['consensuses'].iteritems():
            cdr = consensus['cdr3_seq']

            # Verify that the contig annotation data are consistent with the clonotype assignment data
            assert set(consensus['cell_contigs']) == \
                set(c.contig_name for c in consensus_to_contigs[consensus_id])
            sel_contigs = consensus_to_contigs[consensus_id]
            sel_contig_ids = [c.contig_name for c in sel_contigs]

            # Keep track of the "best" contig. This will be used in case the
            # merging fails.
            best_contig = None

            # Keep track of the set of distinct annotations of the contigs to merge.
            # Will use to report rate of discrepancies.
            feature_annotations = defaultdict(set)

            for contig in sel_contigs:
                for anno in contig.annotations:
                    feature_annotations[anno.feature.region_type].add(anno.feature.gene_name)

                # Always choose a productive over a non-productive. Between
                # contigs with the same productivity, choose the one that had more UMIs.
                if best_contig is None or (not best_contig.productive and contig.productive) or \
                   (best_contig.productive == contig.productive and \
                    best_contig.umi_count < contig.umi_count):

                    best_contig = contig

            assert best_contig is not None

            anno_count = np.max([len(feature_annotations[v]) for v in VDJ_V_FEATURE_TYPES])
            metric = reporter._get_metric_attr('vdj_clonotype_gt1_v_annotations_contig_frac', args.metric_prefix)
            metric.add(1, filter=anno_count > 1)

            anno_count = np.max([len(feature_annotations[v]) for v in VDJ_J_FEATURE_TYPES])
            metric = reporter._get_metric_attr('vdj_clonotype_gt1_j_annotations_contig_frac', args.metric_prefix)
            metric.add(1, filter=anno_count > 1)

            wrong_cdr_metric = reporter._get_metric_attr('vdj_clonotype_consensus_wrong_cdr_contig_frac', args.metric_prefix)

            tmp_dir = martian.make_path(consensus_id + '_outs')
            cr_utils.mkdir(tmp_dir, allow_existing=True)

            res = get_consensus_seq(consensus_id, sel_contig_ids, best_contig.contig_name, tmp_dir, args)
            (best_seq, best_quals, consensus_seq, contig_to_cons_bam, contig_fastq, contig_fasta) = res

            outs.chunked_consensus_bams.append(contig_to_cons_bam)

            # make sure the bam file has the right header (single sequence with this consensus name)
            tmp_bam = tk_bam.create_bam_infile(contig_to_cons_bam)
            if list(tmp_bam.references) != [consensus_id]:
                # Print some info to help us debug
                print tmp_bam.references, consensus_id
                assert(list(tmp_bam.references) == [consensus_id])
            tmp_bam.close()

            if consensus_seq:
                # If this is not None, we actually built a consensus, so we have to compute the quals from scratch.
                # Use a subset of the contigs for computing quals.
                contig_ids = map(lambda c: c.contig_name, sorted(sel_contigs, key=lambda c:c.umi_count, reverse=True))
                contig_ids = contig_ids[0:MAX_CELLS_FOR_BASE_QUALS]

                consensus_quals = get_consensus_quals(in_bam, consensus_id, contig_fasta,
                                                      contig_ids, contig_umis, tmp_dir)
            else:
                consensus_seq = best_seq
                consensus_quals = best_quals

            assert(len(consensus_seq) == len(consensus_quals))

            total_read_count = sum([c.read_count for c in sel_contigs])
            total_umi_count = sum([c.umi_count for c in sel_contigs])

            contig_info_dict = {
                'cells': clonotype['barcodes'],
                'cell_contigs': sel_contig_ids,
                'clonotype_freq': clonotype['freq'],
                'clonotype_prop': clonotype['prop'],
            }

            contig = annotate_consensus_contig(args.vdj_reference_path,
                                               args.min_score_ratios,
                                               args.min_word_sizes,
                                               consensus_id, clonotype_id,
                                               consensus_seq, consensus_quals,
                                               read_count=total_read_count,
                                               umi_count=total_umi_count,
                                               info_dict=contig_info_dict,
                                               primers=args.primers)

            wrong_cdr_metric.add(1, filter=contig.cdr3_seq is None or contig.cdr3_seq != cdr)

            if contig.cdr3_seq is None or contig.cdr3_seq != cdr:
                # Something went wrong. Use "best" contig as the consensus.
                consensus_seq = best_seq
                consensus_quals = best_quals
                contig = annotate_consensus_contig(args.vdj_reference_path,
                                                   args.min_score_ratios,
                                                   args.min_word_sizes,
                                                   consensus_id, clonotype_id,
                                                   consensus_seq, consensus_quals,
                                                   read_count=total_read_count,
                                                   umi_count=total_umi_count,
                                                   info_dict=contig_info_dict,
                                                   primers=args.primers)

            assert(not contig.cdr3_seq is None and contig.cdr3_seq == cdr)

            consensus_contigs.append(contig)

            tk_fasta.write_read_fasta(consensus_fasta, consensus_id, consensus_seq)
            tk_fasta.write_read_fastq(consensus_fastq, consensus_id,
                                      consensus_seq, consensus_quals)
            assert(len(consensus_seq) == len(consensus_quals))

            ref_seq_parts, ref_annos = contig.get_concat_reference_sequence()

            # Align the contigs and consensus to a synthetic concatenated reference
            if ref_seq_parts is not None:
                # Trim the last segment down to the annotated length
                #   to avoid including the entire (500nt) C-region
                ref_seq_parts[-1] = ref_seq_parts[-1][0:ref_annos[-1].annotation_match_end]

                # Concatenate the reference VDJC segments
                ref_seq = reduce(lambda x, y: x + y, ref_seq_parts)
                ref_name = re.sub('consensus', 'concat_ref', consensus_id)

                # Reannotate the reference sequence.
                # Restrict the annotation to the already-called segments to
                #   reduce the risk of discordance between the consensus and
                #   concat_ref annotations.
                ref_contig = annotate_consensus_contig(args.vdj_reference_path,
                                                       args.min_score_ratios,
                                                       args.min_word_sizes,
                                                       ref_name, clonotype_id,
                                                       ref_seq, 'I'*len(ref_seq),
                                                       use_features=set([a.feature.feature_id for a in ref_annos]),
                )
                ref_contigs.append(ref_contig)

                # Add the consensus sequence to the input FASTQ (next to the contigs)
                with open(contig_fastq, 'a') as contig_fq:
                    # Create a fake UMI and barcode
                    header = cr_fastq.AugmentedFastqHeader(consensus_id)
                    header.set_tag(PROCESSED_UMI_TAG, consensus_id)
                    header.set_tag(PROCESSED_BARCODE_TAG, consensus_id)
                    tk_fasta.write_read_fastq(contig_fq, header.to_string(),
                                              consensus_seq, consensus_quals)

                # Reuse this file (this had the assembly output but we don't need it anymore)
                ref_fasta_name = martian.make_path(consensus_id + '_contigs.fasta')
                with open(ref_fasta_name, 'w') as f:
                    tk_fasta.write_read_fasta(f, ref_name, ref_seq)

                # Also append to the final output
                tk_fasta.write_read_fasta(ref_fasta, ref_name, ref_seq)

                cmd = [
                    'vdj_asm', 'base-quals',
                    martian.make_path(consensus_id + '_contigs'),
                    tmp_dir,
                    '--single-end'
                ]
                sys.stderr.write('Running ' + ' '.join(cmd) + '\n')

                tk_subproc.check_call(cmd, cwd=os.getcwd())

                # Move out of tmp dir
                rec_bam = martian.make_path(consensus_id + '_reference.bam')
                cr_utils.move(os.path.join(tmp_dir, consensus_id + '_contigs.bam'), rec_bam)
                outs.chunked_concat_ref_bams.append(rec_bam)

            if os.path.isdir(tmp_dir):
                shutil.rmtree(tmp_dir)

            # Clean up unneeded files ASAP
            rm_files([consensus_id + '_contigs.fasta',
                      consensus_id + '_contigs.fastq'])

            # Merge N most recent BAM files to avoid filesystem overload
            if len(outs.chunked_consensus_bams) >= MERGE_BAMS_EVERY:
                assert len(outs.chunked_consensus_bams) == len(outs.chunked_concat_ref_bams)

                new_cons_bam = martian.make_path('merged-consensus-%03d.bam' % n_merged_bams)
                concatenate_bams(new_cons_bam, outs.chunked_consensus_bams)
                rm_files(outs.chunked_consensus_bams)
                outs.chunked_consensus_bams = [new_cons_bam]

                new_ref_bam = martian.make_path('merged-ref-%03d.bam' % n_merged_bams)
                concatenate_bams(new_ref_bam, outs.chunked_concat_ref_bams)
                rm_files(outs.chunked_concat_ref_bams)
                outs.chunked_concat_ref_bams = [new_ref_bam]

                n_merged_bams += 1

    in_bam.close()

    consensus_fastq.close()
    consensus_fasta.close()
    ref_fasta.close()

    reporter.save(outs.chunked_reporter)

    with open(outs.consensus_annotations_json, 'w') as out_file:
        vdj_annot.save_annotation_list_json(out_file, consensus_contigs)

    with open(outs.concat_ref_annotations_json, 'w') as out_file:
        vdj_annot.save_annotation_list_json(out_file, ref_contigs)


def get_consensus_seq(clonotype_name, sel_contigs, best_contig, out_dir, args):
    """Build a consensus sequence from a set of contigs.

    Args:
    - clonotype_name: Used to prefix output files.
    - sel_contigs: Names of contigs to use for consensus building.
    - best_contig: Name of "best" contig. Will search for this contig's sequence
        and base qualities.
    - out_dir: dir used for temporary results
    - args: stage args.

    - Return value:
    A tuple (best_contig_seq, best_contig_quals, consensus_seq, out_bam_name, out_fastq_name, out_fasta_name).
    - best_contig_seq/best_contig_quals: the sequence and quals of the best contig
    - consensus_seq: the consensus sequence or None if no consensus could be built.
    - out_bam_name: Path of BAM with alignments of contigs to consensus seq.
    - out_fastq_name: FASTQ with contig sequences.
    - out_fasta_name: FASTA with consensus sequence.
    enough reads for consensus.
    """

    best_contig_seq = None
    best_contig_quals = None

    # Input to base quality computation - we don't really need the
    # base qualities because we will replace them by read-based qualities
    # But we need to do this to get proper alignments of contigs against
    # the consensus.
    out_fastq_name = martian.make_path(clonotype_name + '_contigs.fastq')

    # Input to assembly
    out_bam_name = martian.make_path(clonotype_name + '_contigs.bam')

    # The reference in the output bam doesn't really matter.
    out_bam, _ = tk_bam.create_bam_outfile(out_bam_name, ['chr1'], [1])

    # Read the entire fastq (all contigs) and write the selected contigs to
    # a bam for the assembler and a fastq for the aligner.
    with open(args.contigs_fastq, 'r') as f, open(out_fastq_name, 'w') as out_fq:
        fq_iter = tk_fasta.read_generator_fastq(f)
        for (name, seq, quals) in fq_iter:
            if name in sel_contigs:
                if name == best_contig:
                    best_contig_seq = seq
                    best_contig_quals = quals

                header = cr_fastq.AugmentedFastqHeader(name)
                # Create a pseudo-UMI for each input contig
                header.set_tag(PROCESSED_UMI_TAG, name)
                # Put all reads on the same "barcode". This is important, so
                # the assembler assembles all of them together.
                header.set_tag(PROCESSED_BARCODE_TAG, clonotype_name)

                record = pysam.AlignedRead()

                record.reference_start = 0
                record.reference_id = 0
                # Wrap with str() or pysam will crash when given unicode
                record.qname = str(header.to_string())
                record.seq = seq
                record.qual = quals
                record.flag = MAPPED_UNPAIRED_FLAG

                out_bam.write(record)

                # Now change the tags. The final bam concatenation code will pull
                # the tags out of the header, so we want these to be meaningful.
                # Put the real barcode in the barcode tag. The alignment-base-qual
                # code will ignore it anyway.
                header.set_tag(PROCESSED_BARCODE_TAG, name.split('_')[0])
                tk_fasta.write_read_fastq(out_fq, header.to_string(), seq, quals)

    out_bam.close()
    assert(not best_contig_seq is None)

    out_fasta_name = martian.make_path(clonotype_name + '_contigs.fasta')

    # Run the assembler to produce a consensus sequence. Read contig-reads from out_bam_name.
    # The resulting sequences will be in out_dir/<clonotype_name>_contigs.fasta. This is the
    # only output of the assembler we care about.
    if len(sel_contigs) >= MIN_CONTIGS_FOR_CONSENSUS:
        cmd = [
            'vdj_asm', 'asm',
            out_bam_name,
            out_dir,
            '--single-end',
            '--cons', # required so we produce a single output sequence
            '--kmers=0',
            '--min-qual=0',
            '--score-factor=0.0'
        ]
        sys.stderr.write('Running ' + ' '.join(cmd) + '\n')

        tk_subproc.check_call(cmd, cwd=os.getcwd())

        with open(os.path.join(out_dir, clonotype_name + '_contigs.fasta'), 'r') as contig_f:
            lines = contig_f.readlines()
            if lines:
                out_seq = lines[1].strip()
            else:
                # In some rare cases (eg. input contigs have 0 quality), assembly might fail.
                out_seq = None
    else:
        out_seq = None

    # Write the best contig sequence on a new fasta. We need to make sure this has the
    # right contig name because this will be the name written in the bam alignments
    # of the contigs against the consensus
    with open(out_fasta_name, 'w') as f:
        tk_fasta.write_read_fasta(f, clonotype_name, out_seq if out_seq else best_contig_seq)

    # Now align the same reads that were used in vdj_asm against the consensus that you just got.
    # The output will be in out_dir/<clonotype_name> + '_contigs.bam'
    cmd = [
        'vdj_asm', 'base-quals',
        martian.make_path(clonotype_name + '_contigs'),
        out_dir,
        '--single-end'
    ]
    sys.stderr.write('Running ' + ' '.join(cmd) + '\n')

    tk_subproc.check_call(cmd, cwd=os.getcwd())

    # Move the BAM of the contigs aligned against the consensus out of the outs
    # (Will overwrite this bam which was already used as input to assembly).
    cr_utils.move(os.path.join(out_dir, clonotype_name + '_contigs.bam'), out_bam_name)

    return (best_contig_seq, best_contig_quals, out_seq, out_bam_name, out_fastq_name, out_fasta_name)


def get_consensus_quals(in_bam, clonotype_name, in_fasta, sel_contigs, contig_umis, out_dir):
    """Compute base quality scores of a sequence.

    Args:
    - in_bam: bam file to get the list of reads assigned to UMIs on the selected contigs
    - clonotype_name: Used for naming output files.
    - sel_contigs: Contigs that led to the consensus sequence above
    - contig_umis: from contig name to list of umis assigned to that contig

    Return value:
    String with base qualities (in FASTQ format).
    """

    pref = re.sub('.fasta', '', os.path.basename(in_fasta))
    fastq1 = re.sub('.fasta', '_1.fastq', in_fasta)
    fastq2 = re.sub('.fasta', '_2.fastq', in_fasta)

    sel_reads = {}

    for contig in sel_contigs:
        umi_read_count = Counter()
        barcode = contig.split('_')[0]
        contig_read_count = 0

        # Wrap contig w/ str() because pysam crashes on unicode input
        for read in in_bam.fetch(str(contig)):
            # NOTE: Assembler assumes that any tags are part of the read name
            # BUT the bam that we feed to this stage has the tags stripped out
            # of the name.
            umi = read.get_tag(PROCESSED_UMI_TAG)
            if umi in contig_umis[contig] and not read.is_secondary:
                umi_read_count[umi] += 1
                if umi_read_count[umi] >= MAX_READS_PER_UMI:
                    continue

                contig_read_count += 1
                if contig_read_count >= MAX_READS_PER_CONTIG:
                    continue

                if not read.qname in sel_reads:
                    sel_reads[read.qname] = [None, None]
                sel_reads[read.qname][read.is_read2] = read

    with open(fastq1, 'w') as f1, open(fastq2, 'w') as f2:
        for read_name, pair in sel_reads.iteritems():
            read1, read2 = pair[0], pair[1]

            if read1 is None:
                # Replace the UMI with <BC>_<UMI>.
                umi = read2.get_tag(PROCESSED_UMI_TAG)
            else:
                umi = read1.get_tag(PROCESSED_UMI_TAG)

            header = cr_fastq.AugmentedFastqHeader(read_name)
            header.set_tag(PROCESSED_UMI_TAG, barcode + '_' + umi)
            header.set_tag(PROCESSED_BARCODE_TAG, barcode)

            if read1 is None:
                out_seq1 = ""
                out_quals1 = ""
            else:
                out_seq1 = tk_seq.get_rev_comp(read1.seq) if read1.is_reverse else read1.seq
                out_quals1 = read1.qual[::-1] if read1.is_reverse else read1.qual
            tk_fasta.write_read_fastq(f1, header.to_string(), out_seq1, out_quals1)

            if read2 is None:
                out_seq2 = ""
                out_quals2 = ""
            else:
                out_seq2 = tk_seq.get_rev_comp(read2.seq) if read2.is_reverse else read2.seq
                out_quals2 = read2.qual[::-1] if read2.is_reverse else read2.qual
            tk_fasta.write_read_fastq(f2, header.to_string(), out_seq2, out_quals2)

    assert(len(sel_reads) > 0)

    cmd = [
        'vdj_asm', 'base-quals',
        re.sub('.fasta', '', in_fasta),
        out_dir
    ]
    sys.stderr.write('Running ' + ' '.join(cmd) + '\n')

    tk_subproc.check_call(cmd, cwd=os.getcwd())

    with open(os.path.join(out_dir, pref + '.fastq'), 'r') as f:
        lines = f.readlines()

    return lines[3].strip()


def annotate_consensus_contig(reference_path, min_score_ratios, min_word_sizes,
                              contig_name, clonotype_name,
                              seq, quals,
                              read_count=None, umi_count=None,
                              info_dict=None,
                              primers=None,
                              use_features=None):
    """ Given a sequence and some auxiliary info, return a populated AnnotatedContig """

    contig = vdj_annot.AnnotatedContig(contig_name,
                                       seq,
                                       quals=quals,
                                       clonotype=clonotype_name,
                                       read_count=read_count,
                                       umi_count=umi_count,
                                       info_dict=info_dict,
                                       filtered=True,
                                       high_confidence=True)

    res = vdj_annot.setup_feature_aligners(reference_path,
                                           min_score_ratios,
                                           min_word_sizes,
                                           use_features=use_features)
    feature_types, feature_aligners, feature_filters = res

    contig.annotations = contig.annotate_features(feature_types,
                                                  feature_aligners,
                                                  feature_filters)
    if primers:
        primer_aligner, primer_filter = vdj_annot.setup_primer_aligner(primers,
                                                                       VDJ_ANNOTATION_MIN_SCORE_RATIO)
        contig.primer_annotations = contig.annotate_features_by_group(primer_aligner,
                                                                      alignment_filter=primer_filter)

    contig.unannotated_intervals = contig.get_unannotated_intervals()
    contig.annotate_cdr3()

    return contig


def join(args, outs, chunk_defs, chunk_outs):
    if len(chunk_outs) == 0:
        # Set all outputs to null
        for slot in outs.slots:
            setattr(outs, slot, None)
        return

    reporters = [chunk_out.chunked_reporter for chunk_out in chunk_outs]
    final_report = cr_report.merge_reporters(reporters)
    final_report.report_summary_json(outs.summary)

    consensus_contigs = []
    ref_contigs = []
    all_bams = []
    all_ref_bams = []

    for chunk in chunk_outs:
        if chunk.consensus_annotations_json and os.path.isfile(chunk.consensus_annotations_json):
            # Collect consensus annotations
            new_contigs = vdj_annot.load_cell_contigs_from_json(chunk.consensus_annotations_json,
                                                                args.vdj_reference_path,
                                                                group_key='clonotype')
            for cl in new_contigs:
                consensus_contigs.extend(cl.chains)

            # Collect concat_ref annotations
            new_ref_contigs = vdj_annot.load_cell_contigs_from_json(chunk.concat_ref_annotations_json,
                                                                    args.vdj_reference_path,
                                                                    group_key='clonotype')
            for cl in new_ref_contigs:
                ref_contigs.extend(cl.chains)

            all_bams.extend(chunk.chunked_consensus_bams)
            all_ref_bams.extend(chunk.chunked_concat_ref_bams)

    if consensus_contigs:
        all_fastqs = [chunk_out.consensus_fastq for chunk_out in chunk_outs]
        cr_utils.concatenate_files(outs.consensus_fastq, all_fastqs)

        all_fastas = [chunk_out.consensus_fasta for chunk_out in chunk_outs]
        concatenate_and_index_fastas(outs.consensus_fasta, all_fastas)
        outs.consensus_fasta_fai = outs.consensus_fasta + '.fai'

        all_fastas = [chunk_out.concat_ref_fasta for chunk_out in chunk_outs]
        concatenate_and_index_fastas(outs.concat_ref_fasta, all_fastas)
        outs.concat_ref_fasta_fai = outs.concat_ref_fasta + '.fai'

        concatenate_sort_and_index_bams(outs.consensus_bam, all_bams)
        outs.consensus_bam_bai = outs.consensus_bam + '.bai'
        concatenate_sort_and_index_bams(outs.concat_ref_bam, all_ref_bams)
        outs.concat_ref_bam_bai = outs.concat_ref_bam + '.bai'

        # Sort contigs (and clonotypes) by frequency.
        with open(args.clonotype_assignments) as f:
            clonotypes = json.load(f)
        clonotype_freqs = {cid:c['freq'] for cid, c in clonotypes.iteritems()}

    consensus_contigs.sort(key=lambda x:clonotype_freqs[x.clonotype], reverse=True)
    ref_contigs.sort(key=lambda x:clonotype_freqs[x.clonotype], reverse=True)

    with open(outs.consensus_annotations_json, 'w') as out_file:
        vdj_annot.save_annotation_list_json(out_file, consensus_contigs)

    with open(outs.concat_ref_annotations_json, 'w') as out_file:
        vdj_annot.save_annotation_list_json(out_file, ref_contigs)

    with open(outs.consensus_annotations_csv, 'w') as out_file:
        vdj_annot.save_consensus_list_csv(out_file, consensus_contigs)

    with open(outs.clonotypes, 'w') as f:
        vdj_annot.save_clonotype_info_csv(f, consensus_contigs)

    outs.chunked_consensus_bams = []
    outs.chunked_concat_ref_bams = []


def concatenate_and_index_fastas(out_fasta, fastas):
    cr_utils.concatenate_files(out_fasta, fastas)
    tk_subproc.check_call(['samtools', 'faidx', out_fasta], cwd=os.getcwd())

def concatenate_bams(out_bam, bams):
    # Drop the UMI tag since it's useless
    vdj_utils.concatenate_and_fix_bams(out_bam, bams, drop_tags=[PROCESSED_UMI_TAG])

def concatenate_sort_and_index_bams(out_bam, bams):
    tmp_bam = out_bam.replace('.bam', '_tmp.bam')
    # Drop the UMI tag since it's useless
    vdj_utils.concatenate_and_fix_bams(tmp_bam, bams, drop_tags=[PROCESSED_UMI_TAG])
    # Wrap filenames in str() to prevent pysam crash on unicode input
    tk_bam.sort_and_index(str(tmp_bam), str(out_bam.replace('.bam', '')))
    cr_utils.remove(tmp_bam, allow_nonexisting=True)
