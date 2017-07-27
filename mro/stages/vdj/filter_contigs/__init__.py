#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#
# Mark contigs as high confidence

import itertools
import os
import tenkit.fasta as tk_fasta
import cellranger.constants as cr_constants
import cellranger.utils as cr_utils
import cellranger.vdj.annotations as vdj_annot
import cellranger.vdj.constants as vdj_constants
import cellranger.vdj.report as vdj_report
import cellranger.vdj.utils as vdj_utils

__MRO__ = """
stage FILTER_CONTIGS(
    in  path  vdj_reference_path,
    in  fasta contig_fasta,
    in  fastq contig_fastq,
    in  json  contig_annotations,
    out fasta filtered_contig_fasta,
    out fastq filtered_contig_fastq,
    out json  contig_annotations,
    out json  summary,
    src py    "stages/vdj/filter_contigs",
) split using (
)
"""

def split(args):
    mem_gb = max(cr_constants.MIN_MEM_GB,
                 vdj_utils.get_mem_gb_from_annotations_json(args.contig_annotations))
    return {
        'chunks': [{'__mem_gb': mem_gb}],
        'join': {'__mem_gb': mem_gb},
    }

def main(args, outs):
    reporter = vdj_report.VdjReporter()

    with open(args.contig_annotations) as f:
        contigs = vdj_annot.load_contig_list_from_json(f, args.vdj_reference_path)

    contigs.sort(key=lambda c: (c.barcode, c.get_single_chain(), not c.productive, -c.umi_count, -c.read_count, -len(c)))

    low_confidence_contigs = set()
    cell_contigs = set()

    for (bc, chain), group in itertools.groupby(contigs, key=lambda c: (c.barcode, c.get_single_chain())):
        first_cdr3 = None
        seen_cdr3s = set()

        for contig in group:
            contig.high_confidence = True

            if contig.is_cell:
                cell_contigs.add(contig.contig_name)

            if first_cdr3 is None:
                first_cdr3 = contig.cdr3_seq

            # Mark as low confidence:
            # 1) Any additional CDR3s beyond the highest-(productive,UMI,read,length) contig's CDR3
            #    with a single UMI, or
            extraneous_cdr3 = first_cdr3 is not None \
               and contig.cdr3_seq != first_cdr3 \
               and contig.umi_count == 1

            # 2) Any contigs with a repeated CDR3.
            repeat_cdr3 = contig.cdr3_seq in seen_cdr3s

            if extraneous_cdr3 or repeat_cdr3:
                contig.high_confidence = False
                low_confidence_contigs.add(contig.contig_name)

            seen_cdr3s.add(contig.cdr3_seq)

            if chain in vdj_constants.VDJ_GENES:
                reporter._get_metric_attr('vdj_high_conf_prod_contig_frac', chain).add(1, filter=contig.high_confidence)
            reporter._get_metric_attr('vdj_high_conf_prod_contig_frac', cr_constants.MULTI_REFS_PREFIX).add(1, filter=contig.high_confidence)

    # Write augmented contig annotations
    with open(outs.contig_annotations, 'w') as f:
        vdj_annot.save_annotation_list_json(f, contigs)

    # Write filtered fasta
    with open(args.contig_fasta) as in_file, \
         open(outs.filtered_contig_fasta, 'w') as out_file:
        for hdr, seq in cr_utils.get_fasta_iter(in_file):
            # Keep contigs that are high confidence & in cells
            if hdr not in low_confidence_contigs and hdr in cell_contigs:
                tk_fasta.write_read_fasta(out_file, hdr, seq)

    # Write filtered fastq
    with open(args.contig_fastq) as in_file, \
         open(outs.filtered_contig_fastq, 'w') as out_file:
        for name, seq, qual in tk_fasta.read_generator_fastq(in_file):
            if name not in low_confidence_contigs and name in cell_contigs:
                tk_fasta.write_read_fastq(out_file, name, seq, qual)

    reporter.report_summary_json(outs.summary)


def join(args, outs, chunk_defs, chunk_outs):
    # Copy files from single chunk to join
    for out_name in ['summary',
                     'contig_annotations',
                     'filtered_contig_fasta',
                     'filtered_contig_fastq',
    ]:
        src = getattr(chunk_outs[0], out_name)
        dest = getattr(outs, out_name)
        if os.path.isfile(src):
            cr_utils.copy(src, dest)
        else:
            setattr(outs, out_name, None)
