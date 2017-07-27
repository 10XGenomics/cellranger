#!/usr/bin/env python
#
# Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
#
# Align VDJ features and primers to assembled contigs.
# Also report variants relative to the reference feaures.

import cPickle
import os.path
import pandas as pd
from pandas.io.common import EmptyDataError
import cellranger.utils as cr_utils
import cellranger.vdj.annotations as vdj_annot
import cellranger.vdj.constants as vdj_constants
import cellranger.vdj.utils as vdj_utils
import tenkit.fasta as tk_fasta


__MRO__ = """
stage ANNOTATE_CONTIGS(
    in  path       vdj_reference_path,
    in  fasta      contigs,
    in  fastq      contigs_fastq,
    in  string[][] barcodes_in_chunks,
    in  map[]      primers,
    in  csv        filter_summary,
    in  tsv        contig_summary,
    in  map        min_score_ratios     "dict of min score ratios by feature type to use for filtering annotations",
    in  map        min_word_sizes       "dict of min word sizes by feature type to use for filtering annotations",
    in  json       cell_barcodes,
    out pickle     chunked_annotations,
    out json       annotations          "Annotations for filtered contigs",
    out json       raw_annotations      "Annotations for all contigs",
    out bed        annotations_bed      "BED for IGV",
    out csv        annotations_csv,
    src py         "stages/vdj/annotate_contigs",
) split using (
    in  string[]   barcodes,
)
"""

VDJ_ANNOTATE_FEATURE_TYPES = [
    "5'UTR",
    'L-REGION',
    'C-REGION',
    "3'UTR",
] + \
    vdj_constants.VDJ_V_FEATURE_TYPES + \
    vdj_constants.VDJ_J_FEATURE_TYPES

def annotations_json_to_bed(annotations_json):
    """
    Convert JSON formatted annotations to a flattened BED format (list of lists).

    Args:
        annotations_json (list of dict): JSON formatted annotations

    Returns:
        list of lists: list of entries corresponding to BED file entries with all entries as strings.

    """

    bed_entries = []

    for annotation in annotations_json:
        annotations = annotation.get('annotations')
        primers = annotation.get('primer_annotations')

        contig_name = annotation.get('contig_name')

        for feature in annotations + primers:
            start = feature.get('contig_match_start')
            end = feature.get('contig_match_end')
            type = '%s_%s' % (feature.get('gene_name'), feature.get('region_type'))

            bed_entries.append([contig_name, str(start), str(end), type])

    return bed_entries


def split(args):
    chunks = []

    # Reuse the chunks from assemble VDJ to chunk
    if args.barcodes_in_chunks is None:
        # Assume that data are bulk
        barcodes_in_chunks = [['']]  # the get barcode name function returns '' for all bulk contig names
    else:
        barcodes_in_chunks = args.barcodes_in_chunks

    for barcodes in barcodes_in_chunks:
        chunks.append({
            'barcodes': barcodes,
        })

    join_mem = max(6.0, 12.0 * float(len(chunks)) / 1000.0)
    return {'chunks': chunks, 'join': {'__mem_gb': join_mem }}


def main(args, outs):
    if args.vdj_reference_path is None:
        outs.chunked_annotations = None
        return
    chunk_contigs = []
    barcodes_in_chunk = set(args.barcodes)

    # Set of barcodes that were called as cells
    if args.cell_barcodes:
        cell_barcodes_set = set(vdj_utils.load_cell_barcodes_json(args.cell_barcodes))
    else:
        cell_barcodes_set = set()

    # Setup feature reference sequences
    res = vdj_annot.setup_feature_aligners(args.vdj_reference_path,
                                           args.min_score_ratios,
                                           args.min_word_sizes)
    feature_types, feature_aligners, feature_filters = res

    # Setup primer reference sequnces
    if args.primers:
        primer_aligner, primer_filter = vdj_annot.setup_primer_aligner(args.primers,
                                                                       vdj_constants.VDJ_ANNOTATION_MIN_SCORE_RATIO)

    read_counts = {}
    umi_counts = {}
    if args.contig_summary and os.path.isfile(args.contig_summary):
        contig_summary = pd.read_csv(args.contig_summary, header=0, index_col=None, sep='\t')
        for _, row in contig_summary.iterrows():
            read_counts[row.contig_name] = int(row.num_reads)
            umi_counts[row.contig_name] = int(row.num_umis)

    if args.filter_summary:
        try:
            filter_summary = vdj_utils.load_contig_summary_table(open(args.filter_summary))
        except EmptyDataError:
            filter_summary = None
    else:
        filter_summary = None

    if not args.contigs_fastq is None:
        fq_iter = tk_fasta.read_generator_fastq(open(args.contigs_fastq), paired_end=False)

    for header, contig_sequence in cr_utils.get_fasta_iter(open(args.contigs)):
        if args.contigs_fastq is None:
            contig_quals = None
        else:
            header_fq, contig_sequence_fq, contig_quals = fq_iter.next()
            assert(contig_sequence_fq == contig_sequence)
            assert(header_fq == header)

        barcode = vdj_utils.get_barcode_from_contig_name(header)
        contig_name = header.split(' ')[0]

        # Only annotate barcodes assigned to this chunk and contigs with enough read support
        if barcode in barcodes_in_chunk:
            if filter_summary is not None:
                filtered = vdj_utils.is_contig_filtered(filter_summary, contig_name)
            else:
                filtered = True

            contig = vdj_annot.AnnotatedContig(contig_name,
                                               contig_sequence,
                                               quals=contig_quals,
                                               barcode=barcode,
                                               is_cell=barcode in cell_barcodes_set,
                                               filtered=filtered,
                                               read_count=read_counts.get(contig_name),
                                               umi_count=umi_counts.get(contig_name),
                                               )

            contig.annotations = contig.annotate_features(feature_types,
                                                          feature_aligners,
                                                          feature_filters)

            if args.primers:
                contig.primer_annotations = contig.annotate_features_by_group(primer_aligner,
                                                                              alignment_filter=primer_filter)

            contig.annotate_cdr3()

            chunk_contigs.append(contig)

    cPickle.dump(chunk_contigs, open(outs.chunked_annotations, 'wb'),
                 protocol=cPickle.HIGHEST_PROTOCOL)

def join(args, outs, chunk_defs, chunk_outs):
    # Merge pickles and save JSON
    all_contigs = []

    for chunk in chunk_outs:
        if chunk.chunked_annotations is not None:
            all_contigs.extend(cPickle.load(open(chunk.chunked_annotations, 'rb')))

    # Clear this temporary, chunk-specific out
    outs.chunked_annotations = None

    # Write all contigs
    with open(outs.raw_annotations, 'w') as out_file:
        vdj_annot.save_annotation_list_json(out_file, all_contigs)

    # Write filtered contigs
    with open(outs.annotations, 'w') as out_file:
        vdj_annot.save_annotation_list_json(out_file, [c for c in all_contigs if c.filtered])

    # Save a BED formatted file of a subset of annotations
    with open(outs.annotations_bed, 'w') as output_file:
        bed_lines = cr_utils.flatten_list([c.get_annotations_bed() for c in all_contigs if c.filtered])

        for bed_line in bed_lines:
            output_file.write(bed_line + '\n')

    # Write annotations CSV
    with open(outs.annotations_csv, 'w') as csv:
        vdj_annot.save_contig_list_csv(csv, all_contigs)
