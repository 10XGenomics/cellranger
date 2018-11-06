#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
from collections import defaultdict
import itertools
import json
import martian
import numpy as np
import random
import tenkit.constants as tk_constants
import tenkit.safe_json as tk_safe_json
import tenkit.seq as tk_seq
import cellranger.chemistry as cr_chem
import cellranger.constants as cr_constants
import cellranger.rna.feature_ref as rna_feature_ref
import cellranger.rna.library as rna_library
import cellranger.report as cr_report
import cellranger.utils as cr_utils
from cellranger.fastq import BarcodeCounter, FastqReader, FastqFeatureReader, \
    ChunkedFastqWriter, AugmentedFastqHeader, get_bamtofastq_defs
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
    in  int      r1_length,
    in  int      r2_length,
    in  bool     skip_metrics,
    in  path     reference_path,
    in  csv      feature_reference,
    in  bool     augment_fastq,
    in  map[]    library_info,
    out pickle   chunked_reporter,
    out json     summary,
    out json     barcode_counts,
    out fastq[]  reads,
    out fastq[]  read2s,
    out int[]    gem_groups,
    out string[] library_types,
    out string[] read_groups,
    out map      align,
    out string[] bam_comments,
    src py       "stages/common/extract_reads",
) split using (
    in  map      read_chunks,
    in  bool     reads_interleaved,
    in  bool     barcode_rc,
)
"""

COMPRESSION = 'lz4'

def split(args):
    chunks = []
    for chunk in args.chunks:
        chunk['__threads'] = 1
        chunk['__mem_gb'] = 6

        if args.initial_reads is None:
            chunk['chunk_initial_reads'] = None
        else:
            chunk['chunk_initial_reads'] = args.initial_reads / len(args.chunks)

        # Combine the per-chunk subsample rates with the toplevel rate
        chunk['chunk_subsample_rate'] = chunk.pop('subsample_rate') * args.subsample_rate

        # NOTE: Here we implicitly convert SETUP_CHUNKS' output to Martian "chunk" definitions.
        #   Keys coming from SETUP_CHUNKS can override the stage-level args and are not
        #   explicitly typed by Martian. This is bad.

        chunks.append(chunk)
    join = {
        '__mem_gb': 8,
    }
    return {'chunks': chunks, 'join': join}

def join(args, outs, chunk_defs, chunk_outs):
    outs.reads, outs.read2s, outs.tags = [], [], []
    outs.gem_groups, outs.library_types, outs.library_ids, outs.read_groups = [], [], [], []

    for chunk_out in chunk_outs:
        outs.reads += [read for read in chunk_out.reads]
        outs.read2s += [read2 for read2 in chunk_out.read2s]
        outs.tags += [tags for tags in chunk_out.tags]
        outs.gem_groups += [gem_group for gem_group in chunk_out.gem_groups]
        outs.library_types += [lt for lt in chunk_out.library_types]
        outs.library_ids += [li for li in chunk_out.library_ids]
        outs.read_groups += [read_group for read_group in chunk_out.read_groups]

    # Ensure that we have non-zero reads
    if not outs.reads:
        martian.exit("No reads found. Check the input fastqs and/or the chemistry definition")
    # Ensure consistency of BAM comments
    assert all(chunk_out.bam_comments == chunk_outs[0].bam_comments for chunk_out in chunk_outs)
    outs.bam_comments = chunk_outs[0].bam_comments

    # Write barcode counts (merged by library_type)
    bc_counters = BarcodeCounter.merge_by([co.barcode_counts for co in chunk_outs],
                                          [cd.library_type for cd in chunk_defs],
                                          args.barcode_whitelist,
                                          outs.gem_groups)
    with open(outs.barcode_counts, 'w') as f:
        tk_safe_json.dump_numpy(bc_counters, f)

    # Write feature counts
    feature_counts = None
    for chunk_def, chunk_out in itertools.izip(chunk_defs, chunk_outs):
        with open(chunk_out.feature_counts) as f:
            chunk_counts = np.asarray(json.load(f), dtype=int)
            if feature_counts is None:
                feature_counts = chunk_counts
            else:
                feature_counts += chunk_counts

    with open(outs.feature_counts, 'w') as f:
        json.dump(tk_safe_json.json_sanitize(list(feature_counts)), f)


    outs.align = cr_utils.select_alignment_params(args.align)

    # Group reporters by library type
    outs.chunked_reporter = None
    reporter_groups = defaultdict(list)
    for chunk_def, chunk_out in zip(chunk_defs, chunk_outs):
        if not chunk_out.reads:
            continue
        chunk_lib_types = set(lt for lt in chunk_out.library_types)
        assert len(chunk_lib_types) == 1
        lib_type = list(chunk_lib_types)[0]
        reporter_groups[lib_type].append(chunk_out.chunked_reporter)

    # Merge reporters and prefix JSON keys by library type
    summary = {}
    for lib_type, reporters in reporter_groups.iteritems():
        j = cr_report.merge_reporters(reporters).to_json()

        prefix = rna_library.get_library_type_metric_prefix(lib_type)
        j_prefixed = dict((prefix + k, v) for k, v in j.iteritems())

        summary.update(j_prefixed)

    # Use a temporary reporter to generate the metadata (w/o a prefix)
    tmp_reporter = cr_report.Reporter()
    tmp_reporter.store_chemistry_metadata(args.chemistry_def)
    summary.update(tmp_reporter.to_json())

    # Write summary JSON
    with open(outs.summary, 'w') as f:
        tk_safe_json.dump_numpy(summary, f, pretty=True)

def main(args, outs):
    random.seed(0)

    paired_end = cr_chem.is_paired_end(args.chemistry_def)

    # Build the feature reference
    if args.reference_path:
        feature_ref = rna_feature_ref.from_transcriptome_and_csv(args.reference_path,
                                                                 args.feature_reference)
    else:
        feature_ref = rna_feature_ref.FeatureReference.empty()

    # Setup feature barcode extraction
    feature_extractor = rna_feature_ref.FeatureExtractor(feature_ref,
                                                         use_feature_types=[args.library_type])

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

    # Determine which trimmed sequences need to be retained for bamtofastq
    trim_defs = get_bamtofastq_defs(read_defs, read_tags)
    outs.bam_comments = sorted(set(trim_defs.itervalues()))

    num_libraries = len(args.library_info)
    reporter = cr_report.Reporter(umi_length=cr_chem.get_umi_length(args.chemistry_def),
                                  primers=cr_utils.get_primers_from_dicts(args.primers),
                                  num_libraries=num_libraries)

    # Determine if barcode sequences need to be reverse complemented.
    with FastqReader(args.read_chunks, bc_read_def, args.reads_interleaved, None, None) as bc_check_rc:
        barcode_whitelist = cr_utils.load_barcode_whitelist(args.barcode_whitelist, True)
        barcode_rc = infer_barcode_reverse_complement(barcode_whitelist, bc_check_rc.in_iter)

    # Log the untrimmed read lengths to stdout
    r1_read_def = cr_constants.ReadDef(rna_read_def.read_type, 0, None)
    r1_reader = FastqReader(args.read_chunks, r1_read_def, args.reads_interleaved, None, None)

    r1_untrimmed_len = 0
    for read in itertools.islice(r1_reader.in_iter, cr_constants.DETECT_CHEMISTRY_INITIAL_READS):
        r1_untrimmed_len = max(r1_untrimmed_len, len(read[1]))
    print "Read 1 untrimmed length = ", r1_untrimmed_len
    print "Input arg r1_length = ", args.r1_length
    r1_reader.close()

    if paired_end:
        r2_read_def = cr_constants.ReadDef(rna_read2_def.read_type, 0, None)
        r2_reader = FastqReader(args.read_chunks, r2_read_def, args.reads_interleaved, None, None)

        r2_untrimmed_len = 0
        for read in itertools.islice(r2_reader.in_iter, cr_constants.DETECT_CHEMISTRY_INITIAL_READS):
            r2_untrimmed_len = max(r2_untrimmed_len, len(read[1]))
        print "Read 2 untrimmed length = ", r2_untrimmed_len
        print "Input arg r2_length = ", args.r2_length
        r2_reader.close()

    # Setup read iterators.
    r1_length = args.r1_length
    r2_length = args.r2_length

    rna_reads = FastqReader(args.read_chunks, rna_read_def, args.reads_interleaved, r1_length, r2_length)
    rna_read2s = FastqReader(args.read_chunks, rna_read2_def, args.reads_interleaved, r1_length, r2_length)
    bc_reads = FastqReader(args.read_chunks, bc_read_def, args.reads_interleaved, r1_length, r2_length)
    si_reads = FastqReader(args.read_chunks, si_read_def, args.reads_interleaved, r1_length, r2_length)

    if cr_chem.has_umis(args.chemistry_def):
        umi_reads = FastqReader(args.read_chunks, umi_read_def, args.reads_interleaved, r1_length, r2_length)
    else:
        umi_reads = FastqReader(None, None, False, r1_length, r2_length)

    # Record feature counts:
    feature_counts = np.zeros(feature_ref.get_num_features(), dtype=int)

    # If this library type has no feature barcodes, make the reader a NOOP
    if feature_extractor.has_features_to_extract():
        feature_reads = FastqFeatureReader(args.read_chunks, feature_extractor,
                                           args.reads_interleaved,
                                           r1_length, r2_length)
    else:
        feature_reads = FastqReader(None, None, None, r1_length, r2_length)


    fastq_readers = (rna_reads, rna_read2s, bc_reads, si_reads, umi_reads, feature_reads)

    read1_writer = ChunkedFastqWriter(outs.reads, args.reads_per_file, compression=COMPRESSION)
    if paired_end:
        read2_writer = ChunkedFastqWriter(outs.read2s, args.reads_per_file, compression=COMPRESSION)

    tag_writer = None
    if not args.augment_fastq:
        tag_writer = ChunkedFastqWriter(outs.tags, args.reads_per_file, compression=COMPRESSION)

    bc_counter = BarcodeCounter(args.barcode_whitelist, outs.barcode_counts)

    all_read_iter = itertools.izip_longest(*[reader.in_iter for reader in fastq_readers])

    EMPTY_READ = (None, '', '')

    reporter.extract_reads_init()

    for extractions in itertools.islice(all_read_iter, args.chunk_initial_reads):
        # Downsample
        if random.random() > args.chunk_subsample_rate:
            continue

        rna_extraction, rna2_extraction, bc_extraction, si_extraction, umi_extraction, feature_extraction = extractions

        rna_read = rna_extraction if rna_extraction is not None else EMPTY_READ
        rna_read2 = rna2_extraction if rna2_extraction is not None else EMPTY_READ
        bc_read = bc_extraction if bc_extraction is not None else EMPTY_READ
        si_read = si_extraction if si_extraction is not None else EMPTY_READ
        umi_read = umi_extraction if umi_extraction is not None else EMPTY_READ

        if (not rna_read[1]) or (paired_end and (not rna_read2[1])):
            # Read 1 is empty or read 2 is empty (if paired_end)
            # Empty reads causes issue with STAR aligner, so eliminate
            # them here
            continue

        if bc_read != EMPTY_READ:
            # Reverse complement the barcode if necessary
            if barcode_rc:
                bc_read = (bc_read[0], tk_seq.get_rev_comp(bc_read[1]), bc_read[2][::-1])
            # Track the barcode count distribution
            bc_counter.count(*bc_read)

        # Calculate metrics on raw sequences
        lib_idx = [i for i,x in enumerate(args.library_info) if x['library_id'] == args.library_id][0]
        reporter.raw_fastq_cb(rna_read, rna_read2, bc_read, si_read, umi_read, lib_idx,
                              skip_metrics=args.skip_metrics)

        # Construct new fastq headers
        fastq_header1 = AugmentedFastqHeader(rna_read[0])
        fastq_header1.set_tag(tk_constants.SAMPLE_INDEX_TAG, si_read[1])
        fastq_header1.set_tag(tk_constants.SAMPLE_INDEX_QUAL_TAG, si_read[2])
        fastq_header1.set_tag(cr_constants.RAW_BARCODE_TAG, bc_read[1])
        fastq_header1.set_tag(cr_constants.RAW_BARCODE_QUAL_TAG, bc_read[2])
        fastq_header1.set_tag(cr_constants.RAW_UMI_TAG, umi_read[1])
        fastq_header1.set_tag(cr_constants.UMI_QUAL_TAG, umi_read[2])

        feat_raw_bc = None
        feat_proc_bc = None
        feat_qual = None
        feat_ids = None

        if feature_extraction:
            if feature_extraction.barcode:
                feat_raw_bc = feature_extraction.barcode
                feat_qual = feature_extraction.qual

            if len(feature_extraction.ids) > 0:
                feat_proc_bc = feature_extraction.barcode
                feat_ids = ';'.join(feature_extraction.ids)

                # If hit a single feature ID, count its frequency
                if len(feature_extraction.ids) == 1:
                    feature_counts[feature_extraction.indices[0]] += 1

        if feat_raw_bc:
            fastq_header1.set_tag(cr_constants.RAW_FEATURE_BARCODE_TAG, feat_raw_bc)
            fastq_header1.set_tag(cr_constants.FEATURE_BARCODE_QUAL_TAG, feat_qual)
        if feat_ids:
            fastq_header1.set_tag(cr_constants.PROCESSED_FEATURE_BARCODE_TAG, feat_proc_bc)
            fastq_header1.set_tag(cr_constants.FEATURE_IDS_TAG, feat_ids)

        if args.augment_fastq:
            read1_writer.write((fastq_header1.to_string(), rna_read[1], rna_read[2]))
        else:
            read1_writer.write((rna_read[0], rna_read[1], rna_read[2]))
            tag_writer.write((fastq_header1.to_string(), '', ''))

        if paired_end:
            fastq_header2 = AugmentedFastqHeader(rna_read2[0])
            fastq_header2.set_tag(tk_constants.SAMPLE_INDEX_TAG, si_read[1])
            fastq_header2.set_tag(tk_constants.SAMPLE_INDEX_QUAL_TAG, si_read[2])
            fastq_header2.set_tag(cr_constants.RAW_BARCODE_TAG, bc_read[1])
            fastq_header2.set_tag(cr_constants.RAW_BARCODE_QUAL_TAG, bc_read[2])
            fastq_header2.set_tag(cr_constants.RAW_UMI_TAG, umi_read[1])
            fastq_header2.set_tag(cr_constants.UMI_QUAL_TAG, umi_read[2])

            if feat_raw_bc:
                fastq_header2.set_tag(cr_constants.RAW_FEATURE_BARCODE_TAG, feat_raw_bc)
                fastq_header2.set_tag(cr_constants.FEATURE_BARCODE_QUAL_TAG, feat_qual)
            if feat_ids:
                fastq_header2.set_tag(cr_constants.PROCESSED_FEATURE_BARCODE_TAG, feat_proc_bc)
                fastq_header2.set_tag(cr_constants.FEATURE_IDS_TAG, feat_ids)

            if args.augment_fastq:
                read2_writer.write((fastq_header2.to_string(), rna_read2[1], rna_read2[2]))
            else:
                read2_writer.write((rna_read2[0], rna_read2[1], rna_read2[2]))

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
    if not args.augment_fastq:
        tag_writer.close()
    bc_counter.close()

    # Write feature BC read counts
    with open(outs.feature_counts, 'w') as f:
        json.dump(tk_safe_json.json_sanitize(list(feature_counts)), f)

    # Set stage output parameters.
    if len(read1_writer.file_paths) > 0:
        outs.reads = read1_writer.get_out_paths()

        if paired_end:
            outs.read2s = read2_writer.get_out_paths(len(outs.reads))
        else:
            outs.read2s = []

        if args.augment_fastq:
            outs.tags = []
        else:
            outs.tags = tag_writer.get_out_paths(len(outs.tags))

        libraries = args.library_info
        library = [li for li in libraries if li['library_id'] == args.library_id][0]

        outs.gem_groups = [library['gem_group']] * len(outs.reads)
        outs.library_types = [library['library_type']] * len(outs.reads)
        outs.library_ids = [library['library_id']] * len(outs.reads)
        outs.read_groups = [args.read_group] * len(outs.reads)
    else:
        outs.reads = []
        outs.read2s = []
        outs.tags = []
        outs.gem_groups = []
        outs.library_types = []
        outs.library_ids = []
        outs.read_groups = []

    assert len(outs.gem_groups) == len(outs.reads)
    assert args.augment_fastq or len(outs.tags) == len(outs.reads)

    if paired_end:
        assert len(outs.reads) == len(outs.read2s)

    # this is the first reporter stage, so store the pipeline metadata
    reporter.store_pipeline_metadata(martian.get_pipelines_version())

    reporter.save(outs.chunked_reporter)
