#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
import array
import collections
import itertools
import tenkit.bam as tk_bam
import tenkit.seq as tk_seq
import cellranger.constants as cr_constants
import cellranger.report as cr_report
import cellranger.utils as cr_utils

__MRO__ = '''
stage MARK_DUPLICATES(
    in  bam    input,
    in  path   reference_path,
    in  map    align,
    in  int    mem_gb,
    out bam[]  output,
    out pickle chunked_reporter,
    out json   summary,
    src py     "stages/counter/mark_duplicates",
) split using (
    in  string chunk_start,
    in  string chunk_end,
)
'''

def split(args):
    with tk_bam.create_bam_infile(args.input) as in_bam:
        chunks = tk_bam.chunk_bam_records(in_bam, chunk_bound_key=cr_utils.barcode_sort_key,
                                          chunk_size_gb=cr_constants.BAM_CHUNK_SIZE_GB,
                                          max_chunks=cr_constants.MAX_BAM_CHUNKS)
    if args.mem_gb is not None and args.mem_gb > cr_constants.MIN_MEM_GB:
        for chunk in chunks:
            chunk['__mem_gb'] = args.mem_gb
    return {'chunks': chunks}

def correct_umi(seq, counts):
    corrected_seq = seq
    count = counts.get(seq, 0)

    a = array.array('c', seq)
    for pos in xrange(len(a)):
        existing = a[pos]
        for c in tk_seq.NUCS:
            if c == existing:
                continue
            a[pos] = c
            test_str = a.tostring()

            value = counts.get(test_str, 0)
            if value > count or (value == count and corrected_seq < test_str):
                corrected_seq = test_str
                count = value

        a[pos] = existing
    return corrected_seq

def correct_umis(dupe_keys):
    corrected_dupe_keys = collections.defaultdict(dict)
    for dupe_key, umis in dupe_keys.iteritems():
        for umi in umis:
            new_umi = correct_umi(umi, umis)
            if not (new_umi == umi):
                corrected_dupe_keys[dupe_key][umi] = new_umi

    return corrected_dupe_keys

def mark_dupes(bc, gene_id, reads, args, dupe_type, dupe_func, reporter,
               corrected_dupe_keys=None, out_bam=None):
    mark = dupe_type == cr_constants.CDNA_PCR_DUPE_TYPE
    assert not mark or (corrected_dupe_keys is not None and out_bam is not None)

    dupe_key_umi_counts = collections.defaultdict(dict)

    for read in reads:
        read.is_duplicate = False
        if cr_utils.is_read_dupe_candidate(read, cr_utils.get_high_conf_mapq(args.align)):
            dupe_key = dupe_func(read)
            umi_counts = dupe_key_umi_counts[dupe_key]
            umi = cr_utils.get_read_umi(read)

            if corrected_dupe_keys and dupe_key in corrected_dupe_keys and umi in corrected_dupe_keys[dupe_key]:
                corrected_umi = corrected_dupe_keys[dupe_key][umi]
                cr_utils.set_tag(read, cr_constants.PROCESSED_UMI_TAG, umi, corrected_umi)
                umi = corrected_umi

            if umi in umi_counts:
                read.is_duplicate = True
                umi_counts[umi] += 1
            else:
                umi_counts[umi] = 1

        if mark:
            reporter.mark_dupes_corrected_cb(read)
            out_bam.write(read)

        reporter.mark_dupes_bam_cb(read, dupe_type)

    for _, umis in dupe_key_umi_counts.iteritems():
        reporter.mark_dupes_group_cb(gene_id, umis, dupe_type)

    return dupe_key_umi_counts

def main(args, outs):
    outs.coerce_strings()

    in_bam = tk_bam.create_bam_infile(args.input)
    in_bam_chunk = tk_bam.read_bam_chunk(in_bam, (args.chunk_start, args.chunk_end))
    out_bam, _ = tk_bam.create_bam_outfile(outs.output, None, None, template=in_bam)

    chroms = in_bam.references
    reporter = cr_report.Reporter(reference_path=args.reference_path,
                                  high_conf_mapq=cr_utils.get_high_conf_mapq(args.align),
                                  chroms=chroms)

    for (gg, bc, gene_ids), reads_iter in itertools.groupby(in_bam_chunk, key=cr_utils.barcode_sort_key):
        # Ignore reads w/o a valid barcode, unmapped reads and reads that map to more than 1 gene
        if bc is None or gg is None or gene_ids is None or len(gene_ids) != 1:
            for read in reads_iter:
                reporter.mark_dupes_corrected_cb(read)
                out_bam.write(read)
            continue

        reads = list(reads_iter)
        gene_id = gene_ids[0]

        # Count cDNA PCR duplicates with uncorrected UMIs
        dupe_key_umi_counts = mark_dupes(bc, gene_id, reads, args,
                                         cr_constants.CDNA_PCR_UNCORRECTED_DUPE_TYPE,
                                         cr_utils.cdna_pcr_dupe_func,
                                         reporter)

        # Record UMI corrections
        umi_corrections = correct_umis(dupe_key_umi_counts)

        # Mark duplicates for cDNA PCR duplicates with corrected UMIs
        mark_dupes(bc, gene_id, reads, args,
                   cr_constants.CDNA_PCR_DUPE_TYPE,
                   cr_utils.cdna_pcr_dupe_func,
                   reporter,
                   corrected_dupe_keys=umi_corrections,
                   out_bam=out_bam)

        # Count duplicates for SI PCR duplicates with uncorrected UMIs
        mark_dupes(bc, gene_id, reads, args,
                   cr_constants.SI_PCR_DUPE_TYPE,
                   cr_utils.si_pcr_dupe_func,
                   reporter)

    in_bam.close()
    out_bam.close()
    reporter.save(outs.chunked_reporter)

def join(args, outs, chunk_defs, chunk_outs):
    outs.chunked_reporter = None
    outs.output = [chunk_out.output for chunk_out in chunk_outs]

    reporter = cr_report.merge_reporters([chunk_out.chunked_reporter for chunk_out in chunk_outs])
    reporter.report_summary_json(outs.summary)
