#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#
# - Filter input reads based on similarity to reference VDJ segments.
#
# - Correct UMI sequences based on UMI frequencies (per barcode) post filtering.
#

from collections import defaultdict
import itertools
import martian
import numpy as np
import pysam
import cPickle
import tables
import tenkit.bam as tk_bam
import tenkit.seq as tk_seq
import cellranger.chemistry as cr_chem
import cellranger.constants as cr_constants
import cellranger.fastq as cr_fastq
import cellranger.report as cr_report
import cellranger.utils as cr_utils
import cellranger.vdj.constants as vdj_constants
import cellranger.vdj.filter_vdj_reads as vdj_filt
import cellranger.vdj.reference as vdj_reference
import cellranger.vdj.report as vdj_report
import cellranger.vdj.umi_info as vdj_umi_info

__MRO__ = """
stage FILTER_VDJ_READS(
    in  path    vdj_reference_path,
    in  fastq[] read1s                   "Per-chunk fastqs",
    in  fastq[] read2s                   "Per-chunk fastqs",
    in  json[]  chunk_barcodes           "Per-chunk barcodes",
    in  map     chemistry_def,
    in  map     sw_params                "Params for SW alignment (seed, min_sw_score)",
    in  int     min_readpairs_per_umi,
    out pickle  chunked_reporter,
    out pickle  chunked_gene_umi_counts,
    out bam[]   barcode_chunked_bams,
    out json[]  barcodes_in_chunks,
    out tsv[]   reads_per_bc,
    out json    summary,
    out h5      umi_info,
    src py      "stages/vdj/filter_vdj_reads",
) split using (
    in  fastq   read1_chunk,
    in  fastq   read2_chunk,
    in  json    barcodes_chunk,
)
"""

def get_dummy_chunk():
    read1_out_filename = martian.make_path('chunk0_1.fastq')
    read2_out_filename = martian.make_path('chunk0_2.fastq')
    with open(read1_out_filename, 'w'), open(read2_out_filename, 'w'):
        pass
    chunks = [{
        'read1_chunk': read1_out_filename,
        'read2_chunk': read2_out_filename,
        'barcodes_chunk': None,
    }]
    return {'chunks': chunks}


def split(args):
    assert args.read1s is not None and args.read2s is not None

    chunks = []

    # Ensure that data are barcoded
    assert cr_chem.get_barcode_whitelist(args.chemistry_def) is not None

    for read1_fq, read2_fq, barcodes_json in zip(args.read1s, args.read2s,
                                                 args.chunk_barcodes):
        chunks.append({
            'read1_chunk': read1_fq,
            'read2_chunk': read2_fq,
            'barcodes_chunk': barcodes_json,
            '__mem_gb': 3,
        })

    # Martian doesn't like empty chunk lists so create a chunk w/ empty data
    if len(chunks) == 0:
        return get_dummy_chunk()

    return {'chunks': chunks}

def default_dict_int():
    return defaultdict(int)

def get_bc_grouped_pair_iter(bam, paired_end):
    """ Yields (bc, pair_iter)
        where pair_iter yields (AugmentedFastqHeader, (read1, read2|None)) for the barcode """
    wrap_header = lambda pair: (cr_fastq.AugmentedFastqHeader(pair[0].qname), pair)
    get_barcode = lambda hdr_pair: hdr_pair[0].get_tag(cr_constants.PROCESSED_BARCODE_TAG)

    if paired_end:
        bam_iter = vdj_filt.get_pair_iter(bam)
    else:
        bam_iter = itertools.imap(lambda r1: (r1, None), bam)

    return itertools.groupby(itertools.imap(wrap_header, bam_iter), key=get_barcode)

def process_bam_barcode(bam, pair_iter, bc, corrected_umis, reporter,
                        gene_umi_counts_per_bc,
                        strand, out_bam,
                        asm_min_readpairs_per_umi,
                        paired_end):
    """ Process all readpairs from pair_iter, all having the same bc """

    # Note: "gene" in this function is actually "chain"
    # Readpair counts per UMI (per-gene); {gene: {UMI: count}}

    # Note: Using a lambda here breaks cPickle for some reason
    gene_umi_counts = defaultdict(default_dict_int)

    read_pairs_written = 0

    for header, (read1, read2) in pair_iter:
        (gene1, gene2) = reporter.vdj_recombinome_bam_cb(read1, read2, bam, strand)

        umi = header.get_tag(cr_constants.RAW_UMI_TAG)
        corrected_umi = corrected_umis[umi]

        # Count readpairs per UMI
        if gene1 is not None:
            gene_umi_counts[gene1][corrected_umi] += 1
        if gene2 is not None:
            gene_umi_counts[gene2][corrected_umi] += 1
        if gene1 is None and gene2 is None:
            # Allow unmapped UMIs
            gene_umi_counts["None"][corrected_umi] += 1
        gene_umi_counts[cr_constants.MULTI_REFS_PREFIX][corrected_umi] += 1

        header.set_tag(cr_constants.PROCESSED_UMI_TAG, corrected_umi)
        read1.qname = header.to_string()

        if read2 is not None:
            header2 = cr_fastq.AugmentedFastqHeader(read2.qname)
            assert(header2.get_tag(cr_constants.RAW_UMI_TAG) == umi)
            header2.set_tag(cr_constants.PROCESSED_UMI_TAG, corrected_umi)
            read2.qname = header2.to_string()

        reporter._get_metric_attr('vdj_corrected_umi_frac').add(1, filter=corrected_umis[umi] != umi)
        read_pairs_written += 1

        # Write whether this pair was filtered or not.
        out_bam.write(read1)
        if read2 is not None:
            out_bam.write(read2)

    # Report read-pairs/umi
    for gene in reporter.vdj_genes:
        tot_readpairs = 0
        asm_bad_readpairs = 0

        for reads_per_umi in gene_umi_counts[gene].itervalues():
            reporter._get_metric_attr('vdj_recombinome_readpairs_per_umi_distribution', gene).add(reads_per_umi)

            if reads_per_umi < asm_min_readpairs_per_umi:
                asm_bad_readpairs += reads_per_umi
            tot_readpairs += reads_per_umi

        reporter._get_metric_attr('vdj_recombinome_low_support_reads_frac', gene).set_value(asm_bad_readpairs, tot_readpairs)

    gene_umi_counts_per_bc[bc] = gene_umi_counts


def main(args, outs):
    reporter = vdj_report.VdjReporter(vdj_reference_path=args.vdj_reference_path)
    gene_umi_counts_per_bc = {}

    strand = cr_chem.get_strandedness(args.chemistry_def)

    paired_end = cr_chem.is_paired_end(args.chemistry_def)
    assert paired_end != (args.read2_chunk is None)

    # For the entire chunk, match reads against the V(D)J reference
    ref_fasta = vdj_reference.get_vdj_reference_fasta(args.vdj_reference_path)

    # The filtering code will write this bam. Then we'll read it, correct the UMIs
    # and write outs.chunked_bams.
    filter_bam = martian.make_path('tmp.bam')

    vdj_filt.run_read_match(args.read1_chunk, args.read2_chunk,
                            ref_fasta, filter_bam, strand, args.sw_params)

    # Make two passes over the BAM file, processing one barcode at a time
    bam1 = pysam.AlignmentFile(filter_bam, check_sq=False)
    bam2 = pysam.AlignmentFile(filter_bam, check_sq=False)
    bc_iter1 = get_bc_grouped_pair_iter(bam1, paired_end)
    bc_iter2 = get_bc_grouped_pair_iter(bam2, paired_end)

    reads_per_bc = open(outs.reads_per_bc, 'w')
    out_bam, _ = tk_bam.create_bam_outfile(outs.barcode_chunked_bams, None, None, template=bam1)

    for (bc, pair_iter1), (_, pair_iter2) in itertools.izip(bc_iter1, bc_iter2):
        nreads = 0

        # Pass 1: UMI correction
        umi_counts = defaultdict(int)
        for header, (read1, read2) in pair_iter1:
            nreads += 2
            umi_counts[header.get_tag(cr_constants.RAW_UMI_TAG)] += 1

        corrected_umis = correct_umis(umi_counts)

        # Pass 2: Write the UMI-corrected records
        process_bam_barcode(bam1, pair_iter2, bc, corrected_umis,
                            reporter, gene_umi_counts_per_bc, strand,
                            out_bam,
                            args.min_readpairs_per_umi,
                            paired_end)

        reads_per_bc.write('{}\t{}\n'.format(bc, nreads))

    bam1.close()
    bam2.close()
    out_bam.close()

    # Write bc-gene-umi counts
    cPickle.dump(gene_umi_counts_per_bc, open(outs.chunked_gene_umi_counts, 'w'))

    # Copy the input barcodes
    if args.barcodes_chunk is not None:
        cr_utils.copy(args.barcodes_chunk, outs.barcodes_in_chunks)
    else:
        outs.barcodes_in_chunks = None

    reporter.save(outs.chunked_reporter)


def write_umi_info(pickles, filename):
    """" Write an H5 with (bc, chain, read_count) tuples """
    filters = tables.Filters(complevel = cr_constants.H5_COMPRESSION_LEVEL)

    with tables.open_file(filename, 'w', filters=filters) as h5:
        umi_info = vdj_umi_info.create_arrays(h5)

        bc_to_int = {}
        chain_to_int = {}

        for pickle in pickles:
            bc_chain_umi_counts = cPickle.load(open(pickle))

            for bc, chain_umis in bc_chain_umi_counts.iteritems():
                for chain, umi_counts in chain_umis.iteritems():
                    n_umis = len(umi_counts)

                    if chain != cr_constants.MULTI_REFS_PREFIX and n_umis > 0:
                        if bc not in bc_to_int:
                            bc_to_int[bc] = len(bc_to_int)
                        if chain not in chain_to_int:
                            chain_to_int[chain] = len(chain_to_int)

                        umi_info['barcode_idx'].append(np.full(n_umis, bc_to_int[bc],
                                                               dtype=vdj_umi_info.get_dtype('barcode_idx')))
                        umi_info['chain_idx'].append(np.full(n_umis, chain_to_int[chain],
                                                             dtype=vdj_umi_info.get_dtype('chain_idx')))
                        umi_info['reads'].append(np.fromiter(umi_counts.itervalues(),
                                                             vdj_umi_info.get_dtype('reads'), count=n_umis))

        vdj_umi_info.set_ref_column(h5, 'barcodes', np.array(sorted(bc_to_int.keys(), key=bc_to_int.get)))
        vdj_umi_info.set_ref_column(h5, 'chains', np.array(sorted(chain_to_int.keys(), key=chain_to_int.get)))


def join(args, outs, chunk_defs, chunk_outs):
    outs.chunked_reporter = None
    reporter = cr_report.merge_reporters([chunk_out.chunked_reporter for chunk_out in chunk_outs])

    outs.reads_per_bc = [chunk_out.reads_per_bc for chunk_out in chunk_outs]

    outs.barcode_chunked_bams = [chunk_out.barcode_chunked_bams for chunk_out in chunk_outs]

    # Output barcodes in each chunk
    outs.barcodes_in_chunks = [co.barcodes_in_chunks for co in chunk_outs]

    # Write UMI info file
    write_umi_info([c.chunked_gene_umi_counts for c in chunk_outs], outs.umi_info)

    # Record reference info
    if args.vdj_reference_path is not None:
        reporter.store_reference_metadata(args.vdj_reference_path,
                                          vdj_constants.REFERENCE_TYPE,
                                          vdj_constants.REFERENCE_METRIC_PREFIX)

    # Write output json
    reporter.report_summary_json(outs.summary)


def correct_umis(umi_counts):
    corrected_umis = {}

    for umi, umi_count in umi_counts.iteritems():
        corrected_umi = umi
        count = umi_count

        # Try changing each position of the sequence
        for pos in xrange(len(umi)):
            test_umi_chars = [c for c in umi]
            existing = umi[pos]

            # Try every other char
            for c in tk_seq.NUCS:
                if c == existing:
                    continue
                test_umi_chars[pos] = c
                test_umi = ''.join(test_umi_chars)

                test_count = umi_counts.get(test_umi, 0)
                if test_count > count or (test_count == count and corrected_umi < test_umi):
                    corrected_umi = test_umi
                    count = test_count

        corrected_umis[umi] = corrected_umi
    return corrected_umis
