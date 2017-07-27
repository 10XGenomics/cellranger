#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#
# Filter input reads based on similarity to reference VDJ segments.
# Correct UMI sequences based on UMI frequencies (per barcode) post filtering.
#
# Input reads can be provided in two different ways (which are mutually exclusive):
#
#   1. If interleaved_reads is provided, we assume the data are barcoded. This stage
#      will split the data into chunks, with (roughly) reads_per_chunk reads per chunk.
#      Splitting is done in such a way that barcodes are never split across chunks.
#
#   2. If read1s and read2s are provided, we assume the data are non-barcoded.
#      In this case we'll just make a single chunk and copy all reads to that chunk.
#      This is because subsequent stages assume that barcodes do not span chunks.
#
# If the output_fastqs flag is set, then each chunk will write the reads passing
# filtering to a pair of fastqs. Otherwise all reads (filtered or not) will be
# written to a single BAM per chunk. Filtered readpairs have both mates unmapped.
#
from collections import defaultdict
import itertools
import json
import os
import re
import sys
import martian
import numpy as np
import cPickle
import subprocess
import tables
import tenkit.bam as tk_bam
import tenkit.fasta as tk_fasta
import tenkit.seq as tk_seq
import cellranger.chemistry as cr_chem
import cellranger.constants as cr_constants
import cellranger.fastq as cr_fastq
import cellranger.report as cr_report
import cellranger.vdj.constants as vdj_constants
import cellranger.vdj.reference as vdj_reference
import cellranger.vdj.report as vdj_report
import cellranger.vdj.umi_info as vdj_umi_info

__MRO__ = """
stage FILTER_VDJ_READS(
    in  path       vdj_reference_path,
    in  fastq[]    read1s                   "Per-chunk fastqs",
    in  fastq[]    read2s                   "Per-chunk fastqs",
    in  json[]     chunk_barcodes           "Per-chunk barcodes",
    in  map        chemistry_def,
    in  json       extract_reads_summary,
    in  map        sw_params                "Params for SW alignment (seed, min_sw_score)",
    in  bool       output_fastqs            "Output a pair of fastqs instead of a single bam",
    out pickle     chunked_reporter,
    out pickle     chunked_gene_umi_counts,
    out bam[]      barcode_chunked_bams     "Either this or the next two will be populated",
    out fastq[]    barcode_chunked_read1,
    out fastq[]    barcode_chunked_read2,
    out string[][] barcodes_in_chunks,
    out tsv[]      reads_per_bc,
    out json       summary,
    out h5         umi_info,
    src py         "stages/vdj/filter_vdj_reads",
) split using (
    in  fastq      read1_chunk,
    in  fastq      read2_chunk,
    in  string[]   barcodes_chunk,
)
"""

def write_bam_read_fastq(out, read):
    if read.is_reverse:
        seq, qual = tk_seq.get_rev_comp(read.seq), read.qual[::-1]
    else:
        seq, qual = read.seq, read.qual
    tk_fasta.write_read_fastq(out, read.qname, seq, qual)

def run_read_match(fq_pref, fasta_path, out_bam_filename, chemistry_def, sw_params):

    cmd = ['vdj_asm', 'read-match', fasta_path, fq_pref, out_bam_filename,
           '--seed=' + str(sw_params['seed']),
           '--min-sw-score=' + str(sw_params['min_sw_score'])]
    if cr_chem.get_strandedness(chemistry_def) == '-':
        cmd.append('--rev-strand')

    print >> sys.stderr, 'Running', ' '.join(cmd)
    subprocess.check_call(cmd, cwd=os.getcwd())


def get_dummy_chunk():
    read1_out_filename = martian.make_path('chunk0_1.fastq')
    read2_out_filename = martian.make_path('chunk0_2.fastq')
    with open(read1_out_filename, 'w'), open(read2_out_filename, 'w'):
        pass
    chunks = [{
        'read1_chunk': read1_out_filename,
        'read2_chunk': read2_out_filename,
        'barcodes_chunk': [],
    }]
    return {'chunks': chunks}


def split(args):
    assert args.read1s is not None and args.read2s is not None

    chunks = []

    if cr_chem.get_barcode_whitelist(args.chemistry_def) is not None:

        # Data are barcoded
        for read1_fq, read2_fq, barcodes_json in zip(args.read1s, args.read2s,
                                                     args.chunk_barcodes):
            with open(barcodes_json) as f:
                chunk_barcodes = json.load(f)

            chunks.append({
                'read1_chunk': read1_fq,
                'read2_chunk': read2_fq,
                'barcodes_chunk': chunk_barcodes,
                '__mem_gb': 3.0,
            })

    else:
        # Most stages assume that each chunk has a single barcode.
        # So unfortunately we have to put all reads in the same chunk, otherwise
        # metric computation will break.
        read1_out_filename = martian.make_path('chunk0_1.fastq')
        read2_out_filename = martian.make_path('chunk0_2.fastq')
        with open(read1_out_filename, 'w') as read1_out, open(read2_out_filename, 'w') as read2_out:
            for read1_file, read2_file in zip(args.read1s, args.read2s):
                with open(read1_file) as in1, open(read2_file) as in2:
                    fastq1_iter = tk_fasta.read_generator_fastq(in1, paired_end=False)
                    fastq2_iter = tk_fasta.read_generator_fastq(in2, paired_end=False)

                    for read1_tuple in fastq1_iter:
                        read2_tuple = fastq2_iter.next()
                        tk_fasta.write_read_fastq(read1_out, *read1_tuple)
                        tk_fasta.write_read_fastq(read2_out, *read2_tuple)

        chunks.append({
            'read1_chunk': read1_out_filename,
            'read2_chunk': read2_out_filename,
            'barcodes_chunk': [""],
        })

    # Martian doesn't like empty chunk lists so create a chunk w/ empty data
    if len(chunks) == 0:
        return get_dummy_chunk()

    return {'chunks': chunks}

def default_dict_int():
    return defaultdict(int)

def get_pair_iter(iterable):
    """ Return (x_i, x_(i+1)) for i in {0,2,4,...} """
    return itertools.izip(iterable, iterable)

def is_mapped(read1, read2):
    if read2 is not None:
        return not (read1.is_unmapped and read2.is_unmapped)
    else:
        return not read1.is_unmapped

def get_bc_grouped_pair_iter(bam):
    """ Yields (bc, pair_iter)
        where pair_iter yields (AugmentedFastqHeader, (read1, read2)) for the barcode """
    wrap_header = lambda pair: (cr_fastq.AugmentedFastqHeader(pair[0].qname), pair)
    get_barcode = lambda hdr_pair: hdr_pair[0].get_tag(cr_constants.PROCESSED_BARCODE_TAG)

    return itertools.groupby(
        itertools.imap(wrap_header, get_pair_iter(bam)),
        key=get_barcode)


def write_barcode_fastq(bam, pair_iter, bc, corrected_umis, reporter,
                        gene_umi_counts_per_bc,
                        strand, out_bam, out_fastq1, out_fastq2):
    """ Process all readpairs from pair_iter, all having the same bc """

    # Note: "gene" in this function is actually "chain"
    # Readpair counts per UMI (per-gene); {gene: {UMI: count}}
    gene_umi_counts = defaultdict(default_dict_int)

    read_pairs_written = 0

    for header, (read1, read2) in pair_iter:
        (gene1, gene2) = reporter.vdj_recombinome_bam_cb(read1, read2, bam, strand)

        if is_mapped(read1, read2):
            umi = header.get_tag(cr_constants.RAW_UMI_TAG)
            corrected_umi = corrected_umis[umi]

            # Count readpairs per UMI
            if gene1 is not None or gene2 is not None:
                for gene in set(filter(lambda x: x is not None, [gene1, gene2, cr_constants.MULTI_REFS_PREFIX])):
                    gene_umi_counts[gene][corrected_umi] += 1

            header.set_tag(cr_constants.PROCESSED_UMI_TAG, corrected_umi)
            read1.qname = header.to_string()

            header2 = cr_fastq.AugmentedFastqHeader(read2.qname)
            assert(header2.get_tag(cr_constants.RAW_UMI_TAG) == umi)
            header2.set_tag(cr_constants.PROCESSED_UMI_TAG, corrected_umi)
            read2.qname = header2.to_string()

            reporter._get_metric_attr('vdj_corrected_umi_frac').add(1, filter=corrected_umis[umi] != umi)
            read_pairs_written += 1

        if not out_bam is None:
            # Write whether this pair was filtered or not.
            out_bam.write(read1)
            out_bam.write(read2)
        elif is_mapped(read1, read2):
            write_bam_read_fastq(out_fastq1, read1)
            write_bam_read_fastq(out_fastq2, read2)

    # Report read-pairs/umi
    for gene in reporter.vdj_genes:
        for reads_per_umi in gene_umi_counts[gene].itervalues():
            reporter._get_metric_attr('vdj_recombinome_readpairs_per_umi_distribution', gene).add(reads_per_umi)

    gene_umi_counts_per_bc[bc] = gene_umi_counts


def main(args, outs):
    reporter = vdj_report.VdjReporter(vdj_reference_path=args.vdj_reference_path)
    gene_umi_counts_per_bc = {}

    strand = cr_chem.get_strandedness(args.chemistry_def)

    # For the entire chunk, match reads against the V(D)J reference
    ref_fasta = vdj_reference.get_vdj_reference_fasta(args.vdj_reference_path)
    fq_prefix = re.sub('_1.fastq', '', args.read1_chunk)
    # The filtering code will write this bam. Then we'll read it, correct the UMIs
    # and write outs.chunked_bams.
    filter_bam = martian.make_path('tmp.bam')

    run_read_match(fq_prefix, ref_fasta, filter_bam, args.chemistry_def, args.sw_params)

    # Make two passes over the BAM file, processing one barcode at a time
    bam1 = tk_bam.create_bam_infile(filter_bam)
    bam2 = tk_bam.create_bam_infile(filter_bam)
    bc_iter1 = get_bc_grouped_pair_iter(bam1)
    bc_iter2 = get_bc_grouped_pair_iter(bam2)

    reads_per_bc = open(outs.reads_per_bc, 'w')
    if args.output_fastqs:
        out_fastq1 = open(outs.barcode_chunked_read1, 'w')
        out_fastq2 = open(outs.barcode_chunked_read2, 'w')
        out_bam = None
    else:
        out_bam, _ = tk_bam.create_bam_outfile(outs.barcode_chunked_bams, None, None, template=bam1)
        out_fastq1 = None
        out_fastq2 = None

    for (bc, pair_iter1), (_, pair_iter2) in itertools.izip(bc_iter1, bc_iter2):
        nreads = 0

        # Pass 1: UMI correction
        umi_counts = defaultdict(int)
        for header, (read1, read2) in pair_iter1:
            nreads += 2
            if is_mapped(read1, read2):
                umi_counts[header.get_tag(cr_constants.RAW_UMI_TAG)] += 1

        corrected_umis = correct_umis(umi_counts)

        # Pass 2: Write the UMI-corrected records
        write_barcode_fastq(bam1, pair_iter2, bc, corrected_umis,
                            reporter, gene_umi_counts_per_bc, strand,
                            out_bam, out_fastq1, out_fastq2)

        reads_per_bc.write('{}\t{}\n'.format(bc, nreads))

    bam1.close()
    bam2.close()
    if args.output_fastqs:
        out_fastq1.close()
        out_fastq2.close()
    else:
        out_bam.close()

    # Write bc-gene-umi counts
    cPickle.dump(gene_umi_counts_per_bc, open(outs.chunked_gene_umi_counts, 'w'))

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
    if args.output_fastqs:
        outs.barcode_chunked_read1 = [chunk_out.barcode_chunked_read1 for chunk_out in chunk_outs]
        outs.barcode_chunked_read2 = [chunk_out.barcode_chunked_read2 for chunk_out in chunk_outs]
        outs.barcode_chunked_bams = []
    else:
        outs.barcode_chunked_read1 = []
        outs.barcode_chunked_read2 = []
        outs.barcode_chunked_bams = [chunk_out.barcode_chunked_bams for chunk_out in chunk_outs]

    # Output barcodes in each chunk
    outs.barcodes_in_chunks = [chunk_def.barcodes_chunk for chunk_def in chunk_defs]

    # If a single chunk w/ no barcodes, return null for chunk info
    if len(outs.barcodes_in_chunks) == 1 and outs.barcodes_in_chunks[0][0] == '':
        outs.barcodes_in_chunks = None

    # Write UMI info (only for barcoded data)
    if cr_chem.get_barcode_whitelist(args.chemistry_def) is not None:
        write_umi_info([c.chunked_gene_umi_counts for c in chunk_outs], outs.umi_info)

    reporter.store_reference_metadata(args.vdj_reference_path, vdj_constants.REFERENCE_TYPE, vdj_constants.REFERENCE_METRIC_PREFIX)

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
