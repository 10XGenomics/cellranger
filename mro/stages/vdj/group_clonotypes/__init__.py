#!/usr/bin/env python
#
# Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
#
# Assigns barcodes to clonotypes based on CDR3 identity.

import os
import cPickle
import cellranger.constants as cr_constants
import cellranger.utils as cr_utils
import cellranger.vdj.annotations as vdj_annot
import cellranger.vdj.report as vdj_report
import cellranger.vdj.utils as vdj_utils
import tenkit.safe_json as tk_safe_json

__MRO__ = """
stage GROUP_CLONOTYPES(
    in  path   vdj_reference_path,
    in  json   cell_barcodes,
    in  json   annotations,
    in  bool   use_non_productive               "use non-productive CDR sequences in clonotype definitions",
    in  bool   use_non_full_len                 "use CDR3s from non-full-length contigs",
    out json   summary,
    out json   clonotype_assignments            "info about clonotypes",
    out json   contig_annotations               "input annotations updated with clonotype assignment info",
    out csv    contig_annotations_csv           "input annotations updated with clonotype assignment info (CSV)",
    out pickle contig_annotations_pickle        "input annotations updated with clonotype assignment info (pickle)",
    out csv    filtered_contig_annotations_csv  "high-confidence contig annotations w/ clonotype info",
    src py     "stages/vdj/group_clonotypes",
) split using (
)
"""

def split(args):
    mem_gb = max(cr_constants.MIN_MEM_GB,
                 vdj_utils.get_mem_gb_from_annotations_json(args.annotations))
    return {
        'chunks': [{'__mem_gb': mem_gb}],
        'join': {'__mem_gb': mem_gb},
    }

def main(args, outs):
    reporter = vdj_report.VdjReporter()

    cell_barcodes = set(vdj_utils.load_cell_barcodes_json(args.cell_barcodes))

    barcode_contigs = vdj_annot.load_cell_contigs_from_json(args.annotations, args.vdj_reference_path,
                                                            group_key='barcode')

    # From CDR sequence to sequence id
    sequences = {}
    # From clonotype (tuple of CDR ids) to clonotype id
    clonotypes = {}

    # From barcode to clonotype id
    bc_clonotype_assignments = {}

    # First pass: Just keep track of observed CDR3s
    for contig_list in barcode_contigs:

        # This will be a tuple of sequences like "TRA_<cdr seq>"
        barcode_clonotype_tuple = contig_list.clonotype_tuple(require_productive=not args.use_non_productive,
                                                              require_full_len=True,
                                                              require_high_conf=True)

        # Give unique numerical ids to the CDR3 sequences
        if barcode_clonotype_tuple:
            for cdr_seq in barcode_clonotype_tuple:
                sequences.setdefault(cdr_seq, len(sequences))

    # From sequence id to CDR sequence
    sequence_ids = {seq_id:seq for seq, seq_id in sequences.iteritems()}

    # Do a second pass to potentially use non-full length contigs with a valid CDR3.
    for contig_list in barcode_contigs:
        if args.use_non_full_len:
            barcode_clonotype_tuple = []

            for c in contig_list.contigs():
                (_, cl_seq) = c.clonotype_seq()
                # If this contig has a CDR3 and we can infer the gene type of
                # that CDR3 (either based on the contig itself or based on
                # other full-length contigs that had this CDR3, then add this
                # to the clonotype tuple).
                if cl_seq in sequences:
                    # this will rescue contigs that have a chain and CDR3 assigned
                    # but aren't full length
                    barcode_clonotype_tuple.append(cl_seq)
        else:
            barcode_clonotype_tuple = contig_list.clonotype_tuple(require_productive=(not args.use_non_productive),
                                                                  require_full_len=True,
                                                                  require_high_conf=True)
        barcode_clonotype = tuple(sorted(list(set([sequences[s] for s in barcode_clonotype_tuple]))))

        if barcode_clonotype:
            clonotype_id = clonotypes.setdefault(barcode_clonotype, len(clonotypes))
            bc_clonotype_assignments[contig_list.name] = clonotype_id

    # From clonotype id to tuple of CDRs
    clonotype_ids = {clonotype_id:clonotype_tuple for clonotype_tuple, clonotype_id in clonotypes.iteritems()}

    out_clonotypes = vdj_annot.report_clonotypes(reporter, 'raw', cell_barcodes,
                                                 clonotype_ids, sequence_ids,
                                                 barcode_contigs, bc_clonotype_assignments)

    with open(outs.clonotype_assignments, 'w') as out_file:
        tk_safe_json.dump_numpy(tk_safe_json.json_sanitize(out_clonotypes),
                                out_file, pretty=True)

    # Add clonotype assignments to contig annotations
    del barcode_contigs
    with open(args.annotations) as f:
        all_contigs = vdj_annot.load_contig_list_from_json(f, args.vdj_reference_path)

    vdj_annot.label_contigs_with_consensus(out_clonotypes, all_contigs, 'raw')

    # Write augmented contig annotations
    with open(outs.contig_annotations, 'w') as out_file:
        vdj_annot.save_annotation_list_json(out_file, all_contigs)

    with open(outs.contig_annotations_csv, 'w') as out_file:
        vdj_annot.save_contig_list_csv(out_file, all_contigs, write_inferred=False)

    with open(outs.contig_annotations_pickle, 'w') as out_file:
        cPickle.dump(all_contigs, out_file, protocol=cPickle.HIGHEST_PROTOCOL)

    # Write filtered contig annotations
    with open(outs.filtered_contig_annotations_csv, 'w') as out_file:
        filtered_contigs = filter(lambda x: x.high_confidence and x.is_cell, all_contigs)
        vdj_annot.save_contig_list_csv(out_file, filtered_contigs, write_inferred=False)

    reporter.report_summary_json(outs.summary)


def join(args, outs, chunk_defs, chunk_outs):
    # Copy files from single chunk to join
    for out_name in ['summary',
                     'clonotype_assignments',
                     'contig_annotations',
                     'contig_annotations_csv',
                     'filtered_contig_annotations_csv',
                     'contig_annotations_pickle',
                 ]:

        src = getattr(chunk_outs[0], out_name)
        dest = getattr(outs, out_name)
        if os.path.isfile(src):
            cr_utils.copy(src, dest)
        else:
            setattr(outs, out_name, None)
