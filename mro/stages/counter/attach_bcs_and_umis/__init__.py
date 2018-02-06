#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
import itertools
import json
import martian
import os
import sys
import tenkit.log_subprocess as tk_subproc
import cellranger.chemistry as cr_chem
import cellranger.constants as cr_constants
import cellranger.reference as cr_reference
import cellranger.utils as cr_utils

__MRO__ = '''
stage ATTACH_BCS_AND_UMIS(
    in  bam[]    genome_inputs,
    in  path     reference_path,
    in  int[]    gem_groups,
    in  map      chemistry_def,
    in  map      annotation_params,
    in  string   barcode_whitelist,
    in  json     barcode_counts,
    in  float    barcode_confidence_threshold,
    in  int      umi_min_qual_threshold,
    in  string[] bam_comments,
    in  bool     rescue_multimappers,
    in  bool     correct_barcodes,
    in  bool     skip_metrics,
    in  bool     paired_end,
    out bam[]    output,
    out int[]    num_alignments,
    out bincode  chunked_reporter,
    out json     summary,
    out csv      barcodes_detected,
    out path     gene_index_tab,
    out path     metric_chunk_list,
    out json     chunk_metadata,
    src py       "stages/counter/attach_bcs_and_umis",
) split using (
    in  bam      chunk_genome_input,
    in  bam      chunk_trimmed_input,
    in  int      gem_group,
    in  json     bam_comments_json,
)
'''

def split(args):
    # Write BAM comments to json file
    bam_comment_fn = martian.make_path('bam_comments.json')
    with open(bam_comment_fn, 'w') as f:
        json.dump(args.bam_comments, f)

    chunk_mem_gb = cr_utils.get_mem_gb_request_from_barcode_whitelist(args.barcode_whitelist)
    join_mem_gb = cr_utils.get_mem_gb_request_from_barcode_whitelist(args.barcode_whitelist, args.gem_groups)

    chunks = []
    for chunk_genome_input, gem_group in itertools.izip_longest(
            args.genome_inputs, args.gem_groups):
        chunks.append({
            'chunk_genome_input': chunk_genome_input,
            'gem_group': gem_group,
            'bam_comments_json': bam_comment_fn,
            '__mem_gb': chunk_mem_gb,
        })
    join = {
        '__mem_gb': join_mem_gb,
    }
    return {'chunks': chunks, 'join': join}

def main(args, outs):
    convert_pickle_to_rust_index(cr_utils.get_reference_genes_index(args.reference_path), outs.gene_index_tab)

    if args.barcode_whitelist is None:
        barcode_whitelist = 'null'
    elif not os.path.exists(args.barcode_whitelist):
        barcode_whitelist = cr_utils.get_barcode_whitelist_path(args.barcode_whitelist)
    else:
        barcode_whitelist = args.barcode_whitelist

    cmd = [
        'annotate_reads', 'main',
        args.chunk_genome_input,
        outs.output,
        outs.chunked_reporter,
        args.reference_path,
        outs.gene_index_tab,
        args.barcode_counts,
        barcode_whitelist,
        str(args.gem_group),
        outs.chunk_metadata,
        cr_chem.get_strandedness(args.chemistry_def),
        '--bam-comments', args.bam_comments_json,
    ]

    if cr_chem.get_endedness(args.chemistry_def) == cr_constants.FIVE_PRIME:
        cmd.append('--fiveprime')

    print >> sys.stderr, 'Running', ' '.join(cmd)
    tk_subproc.check_call(cmd, cwd=os.getcwd())
    with open(outs.chunk_metadata) as f:
        metadata = json.load(f)
    outs.num_alignments = metadata['num_alignments']

def join(args, outs, chunk_defs, chunk_outs):
    outs.output = [str(chunk_out.output) for chunk_out in chunk_outs]
    outs.chunked_reporter = None
    outs.coerce_strings()

    with open(outs.metric_chunk_list, 'w') as f:
        for chunk_out in chunk_outs:
            f.write(chunk_out.chunked_reporter + '\n')
    cmd = [
        'annotate_reads', 'join',
        outs.metric_chunk_list,
        outs.summary,
        outs.barcodes_detected,
    ]
    print >> sys.stderr, 'Running', ' '.join(cmd)
    tk_subproc.check_call(cmd, cwd=os.getcwd())
    outs.num_alignments = [chunk_out.num_alignments for chunk_out in chunk_outs]

def convert_pickle_to_rust_index(pickle_path, out_path):
    # TODO: we could possibly avoid this by using serde-pickle
    gene_index = cr_reference.GeneIndex.load_pickle(pickle_path)
    with open(out_path, 'w') as writer:
        writer.write('\t'.join(["transcript_id", "gene_id", "gene_name", "transcript_length"]))
        for (transcript_id, transcript) in gene_index.transcripts.iteritems():
            writer.write('\n' + '\t'.join([transcript_id, transcript.gene.id, transcript.gene.name, str(transcript.length)]))
