#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
import itertools
import martian
import numpy as np
import pysam
import random
import tenkit.constants as tk_constants
import tenkit.seq as tk_seq
import cellranger.chemistry as cr_chem
import cellranger.constants as cr_constants
import cellranger.report as cr_report
import cellranger.utils as cr_utils
from cellranger.fastq import BarcodeCounter, FastqReader, ChunkedFastqWriter, ChunkedBamWriter, AugmentedFastqHeader, compute_trim_defs
from cellranger.fastq import infer_barcode_reverse_complement

__MRO__ = """
stage EXTRACT_READS(
    in  map[]    chunks,
    in  map      chemistry_def,
    in  string   barcode_whitelist,
    in  int      reads_per_file,
    in  float    subsample_rate,
    in  int      initial_reads,
    in  map[]    primers,
    in  map      align,
    in  int      rna_read_length,
    in  bool     skip_metrics,
    out pickle   chunked_reporter,
    out json     summary,
    out json     barcode_counts,
    out fastq[]  reads,
    out fastq[]  read2s,
    out bam[]    trimmed_seqs,
    out int[]    gem_groups,
    out string[] read_groups,
    out map      align,
    out string[] bam_comments,
    src py       "stages/common/extract_reads",
) split using (
    in  map      read_chunks,
    in  int      gem_group,
    in  bool     reads_interleaved,
    in  bool     barcode_rc,
)
"""

def split(args):
    gem_groups = [chunk['gem_group'] for chunk in args.chunks]
    chunk_mem_gb = cr_utils.get_mem_gb_request_from_barcode_whitelist(args.barcode_whitelist)
    join_mem_gb = cr_utils.get_mem_gb_request_from_barcode_whitelist(args.barcode_whitelist, gem_groups)

    chunks = []
    for chunk in args.chunks:
        chunk['__mem_gb'] = chunk_mem_gb

        if args.initial_reads is None:
            chunk['initial_reads'] = None
        else:
            chunk['initial_reads'] = args.initial_reads / len(args.chunks)

        chunks.append(chunk)
    join = {
        '__mem_gb': join_mem_gb,
    }
    return {'chunks': chunks, 'join': join}

def join(args, outs, chunk_defs, chunk_outs):
    outs.reads, outs.read2s, outs.trimmed_seqs, outs.gem_groups, outs.read_groups = [], [], [], [], []

    for chunk_out in chunk_outs:
        outs.reads += [read for read in chunk_out.reads]
        outs.read2s += [read2 for read2 in chunk_out.read2s]
        outs.trimmed_seqs += [trimmed_seq for trimmed_seq in chunk_out.trimmed_seqs]
        outs.gem_groups += [gem_group for gem_group in chunk_out.gem_groups]
        outs.read_groups += [read_group for read_group in chunk_out.read_groups]

    # Ensure consistency of BAM comments
    assert all(chunk_out.bam_comments == chunk_outs[0].bam_comments for chunk_out in chunk_outs)
    outs.bam_comments = chunk_outs[0].bam_comments

    bc_counter = BarcodeCounter(args.barcode_whitelist, outs.barcode_counts, gem_groups=outs.gem_groups)
    for chunk_def, chunk_out in zip(chunk_defs, chunk_outs):
        bc_counter.merge(chunk_def.gem_group, chunk_out.barcode_counts)
    bc_counter.close()

    outs.align = cr_utils.select_alignment_params(args.align)

    outs.chunked_reporter = None
    reporter = cr_report.merge_reporters([chunk_out.chunked_reporter for chunk_out in chunk_outs])

    reporter.store_chemistry_metadata(args.chemistry_def)

    reporter.report_summary_json(outs.summary)

def main(args, outs):
    random.seed(0)

    paired_end = cr_chem.is_paired_end(args.chemistry_def)

    # Use the chemistry to get the locations of various sequences
    rna_read_def = cr_chem.get_rna_read_def(args.chemistry_def)
    rna_read2_def = cr_chem.get_rna_read2_def(args.chemistry_def)
    bc_read_def = cr_chem.get_barcode_read_def(args.chemistry_def)
    si_read_def = cr_chem.get_si_read_def(args.chemistry_def)
    umi_read_def = cr_chem.get_umi_read_def(args.chemistry_def)

    read_defs = [rna_read_def, rna_read2_def,
                 bc_read_def, si_read_def, umi_read_def]
    read_tags = [None, None,
                 (cr_constants.RAW_BARCODE_TAG, cr_constants.RAW_BARCODE_QUAL_TAG),
                 (tk_constants.SAMPLE_INDEX_TAG, tk_constants.SAMPLE_INDEX_QUAL_TAG),
                 (cr_constants.RAW_UMI_TAG, cr_constants.UMI_QUAL_TAG),
             ]

    # Determine which trimmed sequences need to be retained
    trim_defs = compute_trim_defs(read_defs, read_tags, args.chemistry_def.get('retain_trimmed_suffix_read'))

    outs.bam_comments = sorted(set([td.bam_to_fastq for td in trim_defs.itervalues()]))

    gem_groups = [chunk['gem_group'] for chunk in args.chunks]
    reporter = cr_report.Reporter(umi_length=cr_chem.get_umi_length(args.chemistry_def),
                                  primers=cr_utils.get_primers_from_dicts(args.primers),
                                  gem_groups=gem_groups)

    # Determine if barcode sequences need to be reverse complemented.
    bc_check_rc = FastqReader(args.read_chunks, bc_read_def, args.reads_interleaved, None)
    barcode_whitelist = cr_utils.load_barcode_whitelist(args.barcode_whitelist)
    barcode_rc = infer_barcode_reverse_complement(barcode_whitelist, bc_check_rc.in_iter)
    bc_check_rc.close()


    # Determine which read_iters need to retain trimmed sequence
    # (only one per read-type e.g., one per R1, one per R2, etc.)
    read_types_with_trim_def = set()
    rna_read_trim_defs = None
    rna_read2_trim_defs = None
    bc_read_trim_defs = None
    si_read_trim_defs = None
    umi_read_trim_defs = None

    if rna_read_def.read_type not in read_types_with_trim_def:
        rna_read_trim_defs = trim_defs
        read_types_with_trim_def.add(rna_read_def.read_type)
    if rna_read2_def.read_type not in read_types_with_trim_def:
        rna_read2_trim_defs = trim_defs
        read_types_with_trim_def.add(rna_read2_def.read_type)
    if bc_read_def.read_type not in read_types_with_trim_def:
        bc_read_trim_defs = trim_defs
        read_types_with_trim_def.add(bc_read_def.read_type)
    if si_read_def.read_type not in read_types_with_trim_def:
        si_read_trim_defs = trim_defs
        read_types_with_trim_def.add(si_read_def.read_type)
    if umi_read_def.read_type not in read_types_with_trim_def:
        umi_read_trim_defs = trim_defs
        read_types_with_trim_def.add(umi_read_def.read_type)

    # Setup read iterators.
    rna_reads = FastqReader(args.read_chunks, rna_read_def, args.reads_interleaved, rna_read_trim_defs)
    rna_read2s = FastqReader(args.read_chunks, rna_read2_def, args.reads_interleaved, rna_read2_trim_defs)
    bc_reads = FastqReader(args.read_chunks, bc_read_def, args.reads_interleaved, bc_read_trim_defs)
    si_reads = FastqReader(args.read_chunks, si_read_def, args.reads_interleaved, si_read_trim_defs)

    if cr_chem.has_umis(args.chemistry_def):
        umi_reads = FastqReader(args.read_chunks, umi_read_def, args.reads_interleaved, umi_read_trim_defs)
    else:
        umi_reads = FastqReader(None, None, False, None)

    fastq_readers = (rna_reads, rna_read2s, bc_reads, si_reads, umi_reads)

    # Compute trim order of the readers; this is to ensure stability in the ordering
    # in which trimmed sequence is added to the TRIMMED_SEQ tags
    trim_order = list(np.argsort([reader.read_def.read_type for reader in fastq_readers if reader.read_def is not None]))

    read1_writer = ChunkedFastqWriter(outs.reads, args.reads_per_file)
    if paired_end:
        read2_writer = ChunkedFastqWriter(outs.read2s, args.reads_per_file)

    bc_counter = BarcodeCounter(args.barcode_whitelist, outs.barcode_counts)

    all_read_iter = itertools.izip_longest(*[reader.in_iter for reader in fastq_readers])

    # Bam file to write auxiliary data to (that won't fit in a fastq hdr / QNAME)
    trimmed_seq_writer = ChunkedBamWriter(outs.trimmed_seqs, args.reads_per_file)

    EMPTY_READ = (None, '', '')

    reporter.extract_reads_init()

    for extractions in itertools.islice(all_read_iter, args.initial_reads):
        # Downsample
        if random.random() > args.subsample_rate:
            continue

        rna_extraction, rna2_extraction, bc_extraction, si_extraction, umi_extraction = extractions

        rna_read = rna_extraction.read if rna_extraction is not None else EMPTY_READ
        rna_read2 = rna2_extraction.read if rna2_extraction is not None else EMPTY_READ
        bc_read = bc_extraction.read if bc_extraction is not None else EMPTY_READ
        si_read = si_extraction.read if si_extraction is not None else EMPTY_READ
        umi_read = umi_extraction.read if umi_extraction is not None else EMPTY_READ

        # Extra trimming for internal purposes
        if args.rna_read_length is not None:
            rna_read = (rna_read[0],
                        rna_read[1][0:args.rna_read_length],
                        rna_read[2][0:args.rna_read_length])

        # Accumulate trimmed sequence; ordering is by read-type (I1,I2,R1,R2)
        # to ensure stability
        trimmed_seq = ''
        trimmed_qual = ''
        for i in trim_order:
            if extractions[i] is None:
                continue
            trimmed_seq += extractions[i].trimmed_seq
            trimmed_qual += extractions[i].trimmed_qual

        if bc_read != EMPTY_READ:
            # Reverse complement the barcode if necessary
            if barcode_rc:
                bc_read = (bc_read[0], tk_seq.get_rev_comp(bc_read[1]), bc_read[2][::-1])
            # Track the barcode count distribution
            bc_counter.count(*bc_read)

        # Calculate metrics on raw sequences
        reporter.raw_fastq_cb(rna_read, rna_read2, bc_read, si_read, umi_read, args.gem_group, skip_metrics=args.skip_metrics)

        # Construct new fastq headers
        fastq_header1 = AugmentedFastqHeader(rna_read[0])
        fastq_header1.set_tag(tk_constants.SAMPLE_INDEX_TAG, si_read[1])
        fastq_header1.set_tag(tk_constants.SAMPLE_INDEX_QUAL_TAG, si_read[2])
        fastq_header1.set_tag(cr_constants.RAW_BARCODE_TAG, bc_read[1])
        fastq_header1.set_tag(cr_constants.RAW_BARCODE_QUAL_TAG, bc_read[2])
        fastq_header1.set_tag(cr_constants.RAW_UMI_TAG, umi_read[1])
        fastq_header1.set_tag(cr_constants.UMI_QUAL_TAG, umi_read[2])

        fastq_header_str1 = fastq_header1.to_string()

        read1_writer.write((fastq_header_str1, rna_read[1], rna_read[2]))

        # Write trimmed sequence data to a separate, unaligned BAM file
        # Note: We assume that there is only one trimmed sequence per read-pair
        trimmed_seq_data = pysam.AlignedSegment()
        trimmed_seq_data.query_name = fastq_header_str1.split(AugmentedFastqHeader.WORD_SEP)[0]
        trimmed_seq_data.flag = 4
        trimmed_seq_data.seq = trimmed_seq
        trimmed_seq_data.qual = trimmed_qual
        trimmed_seq_writer.write(trimmed_seq_data)

        if paired_end:
            fastq_header2 = AugmentedFastqHeader(rna_read2[0])
            fastq_header2.set_tag(tk_constants.SAMPLE_INDEX_TAG, si_read[1])
            fastq_header2.set_tag(tk_constants.SAMPLE_INDEX_QUAL_TAG, si_read[2])
            fastq_header2.set_tag(cr_constants.RAW_BARCODE_TAG, bc_read[1])
            fastq_header2.set_tag(cr_constants.RAW_BARCODE_QUAL_TAG, bc_read[2])
            fastq_header2.set_tag(cr_constants.RAW_UMI_TAG, umi_read[1])
            fastq_header2.set_tag(cr_constants.UMI_QUAL_TAG, umi_read[2])

            read2_writer.write((fastq_header2.to_string(), rna_read2[1], rna_read2[2]))

    reporter.extract_reads_finalize()

    # Close input and output files.
    rna_reads.close()
    if paired_end:
        rna_read2s.close()
    bc_reads.close()
    si_reads.close()
    umi_reads.close()

    read1_writer.close()
    if paired_end:
        read2_writer.close()
    bc_counter.close()

    trimmed_seq_writer.close()

    # Set stage output parameters.
    if len(read1_writer.file_paths) > 0:
        outs.reads = read1_writer.get_out_paths()
        if paired_end:
            outs.read2s = read2_writer.get_out_paths(len(outs.reads))
        else:
            outs.read2s = []
        outs.gem_groups = [args.gem_group] * len(outs.reads)
        outs.read_groups = [args.read_group] * len(outs.reads)
        outs.trimmed_seqs = trimmed_seq_writer.get_out_paths()
    else:
        outs.reads = []
        outs.read2s = []
        outs.gem_groups = []
        outs.read_groups = []
        outs.trimmed_seqs = []

    assert len(outs.gem_groups) == len(outs.reads)
    if paired_end:
        assert len(outs.reads) == len(outs.read2s)
    assert len(outs.trimmed_seqs) == len(outs.reads)

    # this is the first reporter stage, so store the pipeline metadata
    reporter.store_pipeline_metadata(martian.get_pipelines_version())

    reporter.save(outs.chunked_reporter)
