#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
import itertools
import time
import tenkit.bam as tk_bam
import tenkit.constants as tk_constants
import cellranger.chemistry as cr_chem
import cellranger.constants as cr_constants
from cellranger.fastq import AugmentedFastqHeader
import cellranger.reference as cr_reference
import cellranger.report as cr_report
import cellranger.stats as cr_stats
import cellranger.transcriptome as cr_transcriptome
import cellranger.utils as cr_utils

__MRO__ = '''
stage ATTACH_BCS_AND_UMIS(
    in  bam[]    genome_inputs,
    in  bam[]    trimmed_inputs,
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
    out pickle   chunked_reporter,
    out json     summary,
    out h5       barcode_summary,
    src py       "stages/counter/attach_bcs_and_umis",
) split using (
    in  bam      chunk_genome_input,
    in  bam      chunk_trimmed_input,
    in  int      gem_group,
)
'''

def split(args):
    chunk_mem_gb = cr_utils.get_mem_gb_request_from_barcode_whitelist(args.barcode_whitelist)
    join_mem_gb = cr_utils.get_mem_gb_request_from_barcode_whitelist(args.barcode_whitelist, args.gem_groups)

    chunks = []
    for chunk_genome_input, chunk_trimmed_input, gem_group in itertools.izip_longest(
            args.genome_inputs, args.trimmed_inputs or [], args.gem_groups):
        chunks.append({
            'chunk_genome_input': chunk_genome_input,
            'chunk_trimmed_input': chunk_trimmed_input,
            'gem_group': gem_group,
            '__mem_gb': chunk_mem_gb,
        })
    join = {
        '__mem_gb': join_mem_gb,
    }
    return {'chunks': chunks, 'join': join}

def main(args, outs):
    reference_star_path = cr_utils.get_reference_star_path(args.reference_path)
    star_index = cr_transcriptome.build_star_index(reference_star_path)
    chroms = star_index[0][0]
    gene_index = cr_reference.GeneIndex.load_pickle(cr_utils.get_reference_genes_index(args.reference_path))
    barcode_whitelist = cr_utils.load_barcode_whitelist(args.barcode_whitelist)
    barcode_dist = cr_utils.load_barcode_dist(args.barcode_counts, barcode_whitelist, args.gem_group)
    reporter = cr_report.Reporter(reference_path=args.reference_path,
                                  high_conf_mapq=cr_constants.STAR_DEFAULT_HIGH_CONF_MAPQ,
                                  gene_index=gene_index,
                                  chroms=chroms,
                                  barcode_whitelist=barcode_whitelist,
                                  barcode_dist=barcode_dist,
                                  gem_groups=args.gem_groups,
                                  umi_length=cr_chem.get_umi_length(args.chemistry_def),
                                  umi_min_qual_threshold=args.umi_min_qual_threshold)

    reporter.attach_bcs_init()
    outs.num_alignments = process_alignments(args.chunk_genome_input, args.chunk_trimmed_input, outs.output, args.bam_comments, reporter, gene_index, star_index, args)
    reporter.attach_bcs_finalize()
    reporter.save(outs.chunked_reporter)

def join(args, outs, chunk_defs, chunk_outs):
    outs.output = [str(chunk_out.output) for chunk_out in chunk_outs]
    outs.num_alignments = [chunk_out.num_alignments for chunk_out in chunk_outs]
    outs.chunked_reporter = None
    outs.coerce_strings()

    reporter = cr_report.merge_reporters([chunk_out.chunked_reporter for chunk_out in chunk_outs])
    reporter.store_reference_metadata(args.reference_path, cr_constants.REFERENCE_TYPE, cr_constants.REFERENCE_METRIC_PREFIX)
    reporter.report_summary_json(outs.summary)
    reporter.report_barcodes_h5(outs.barcode_summary)

# HELPER METHODS

def process_alignments(genome_bam_file, trimmed_bam_file, out_bam_file, bam_comments, reporter, gene_index, star_index, args):
    in_genome_bam = tk_bam.create_bam_infile(genome_bam_file)
    bam_open_time = time.time()
    in_trimmed_bam = tk_bam.create_bam_infile(trimmed_bam_file) if trimmed_bam_file else None
    out_bam, _ = tk_bam.create_bam_outfile(out_bam_file, None, None, template=in_genome_bam, cos=bam_comments)
    bam_iter = cr_utils.iter_by_qname(in_genome_bam, in_trimmed_bam)
    num_alignments = 0
    read_consume_time = None
    for qname, reads_iter, trimmed_iter in bam_iter:
        reads = list(reads_iter)
        if read_consume_time is None:
            read_consume_time = time.time()
            # if streaming, verify we're actually streaming
            print "Time to first read: %f seconds" % (read_consume_time - bam_open_time)
        num_alignments += len(reads)
        trimmed = list(trimmed_iter)
        trimmed_read = trimmed[0] if len(trimmed) > 0 else None
        for read in process_qname(qname, reads, trimmed_read, reporter, gene_index, star_index, args):
            out_bam.write(read)

    in_genome_bam.close()
    if in_trimmed_bam: in_trimmed_bam.close()
    out_bam.close()

    return num_alignments

def process_qname(qname, reads, trimmed_read, reporter, gene_index, star_index, args):
    stripped_qname, bc_tags, barcode_info, umi_info = make_barcode_tags(qname, reporter, args)
    trim_tags = make_trimmed_tags(trimmed_read) if trimmed_read else []

    all_genes = set()
    saw_primary = False
    primary_candidate_reads = None
    secondary_reads = []
    any_antisense = False
    any_discordant_pairs = False
    transcript_insert_sizes = []
    regions = []

    annotation_params = args.annotation_params or {}
    for read1, read2 in cr_utils.iter_read_pairs(reads, args.paired_end):
        # restore original qname
        read1.qname = stripped_qname
        if read2: read2.qname = stripped_qname

        # NOTE: STAR always marks read2 as secondary, even if read1 is primary. fix that here.
        if read2: read2.is_secondary = read1.is_secondary

        r1_ann_tags, r2_ann_tags, genes, insert_sizes, antisense, discordant_pair, r1_region, r2_region = make_annotation_tags(
                read1, read2, gene_index, star_index, reporter, annotation_params)
        all_genes.update(genes)
        transcript_insert_sizes += insert_sizes
        any_antisense |= antisense
        any_discordant_pairs |= discordant_pair
        regions.append(r1_region)
        if read2: regions.append(r2_region)

        if read1.mapq == cr_constants.STAR_DEFAULT_HIGH_CONF_MAPQ:
            saw_primary = True
        elif len(genes) > 0 and primary_candidate_reads is None:
            primary_candidate_reads = [read1]
            if read2: primary_candidate_reads.append(read2)
        else:
            secondary_reads.append(read1)
            if read2: secondary_reads.append(read2)

        read1.tags += r1_ann_tags + bc_tags + trim_tags
        if read2: read2.tags += r2_ann_tags + bc_tags + trim_tags

    rescuable = len(all_genes) == 1 and not saw_primary and primary_candidate_reads is not None
    if args.rescue_multimappers and rescuable:
        for primary_candidate in primary_candidate_reads:
            primary_candidate.mapq = cr_constants.STAR_DEFAULT_HIGH_CONF_MAPQ
            primary_candidate.is_secondary = False
            primary_candidate.tags += [(cr_constants.MULTIMAPPER_TAG, 1)]
        for secondary_read in secondary_reads:
            secondary_read.mapq = 0
            secondary_read.is_secondary = True
            secondary_read.tags += [(cr_constants.MULTIMAPPER_TAG, 0)]

    reporter.aligned_bam_cb(reads, transcript_insert_sizes, all_genes, regions, barcode_info, umi_info, any_antisense, any_discordant_pairs)

    return reads

def make_annotation_tags(read1, read2, gene_index, star_index, reporter, annotation_params):
    chrom_info, transcript_info, exon_info = star_index
    r1_tags = []
    r2_tags = []
    all_genes = set()
    transcript_insert_sizes = []
    any_antisense = False
    discordant_pair = False
    r1_region = None
    r2_region = None

    if read2 is not None: # paired end
        # align to transcriptome
        r1_hits, r1_antisense, r1_region = cr_transcriptome.align_to_transcriptome(read1, chrom_info, transcript_info, exon_info, annotation_params)
        r2_hits, r2_antisense, r2_region = cr_transcriptome.align_to_transcriptome(read2, chrom_info, transcript_info, exon_info, annotation_params)
        any_antisense = r1_antisense or r2_antisense

        # annotate both ends with the union of genes / transcripts
        r1_hits_union = []
        r2_hits_union = []
        r1_transcript_ids = set(r1_hits)
        r2_transcript_ids = set(r2_hits)
        tx_union = r1_transcript_ids | r2_transcript_ids
        tx_interection = r1_transcript_ids & r2_transcript_ids
        discordant_pair = len(tx_union) > 0 and len(tx_interection) == 0
        for transcript_id in sorted(tx_union):
            r1_hit = r1_hits.get(transcript_id, None)
            r2_hit = r2_hits.get(transcript_id, None)
            # compute insert size if possible
            if r1_hit is not None and r2_hit is not None:
                insert_size_for_transcript = compute_paired_insert_size(r1_hit, r2_hit)
                if insert_size_for_transcript <= cr_constants.MAX_INSERT_SIZE:
                    transcript_insert_sizes.append(insert_size_for_transcript)
            # add dummy alignment if necessary
            r1_hits_union.append(r1_hit or cr_transcriptome.make_dummy_alignment(transcript_id))
            r2_hits_union.append(r2_hit or cr_transcriptome.make_dummy_alignment(transcript_id))

        if len(r1_hits_union) > 0:
            gene_hits = cr_utils.get_genes_from_transcripts(r1_hits_union, gene_index)
            all_genes.update(gene_hits)
            r1_tx = cr_transcriptome.make_tx_tag(r1_hits_union)
            r2_tx = cr_transcriptome.make_tx_tag(r2_hits_union)
            gx, gn = cr_transcriptome.make_gx_tag(gene_hits)
            r1_tags += cr_utils.make_annotation_tags(r1_tx, gx, gn)
            r2_tags += cr_utils.make_annotation_tags(r2_tx, gx, gn)

    else: # single end
        r1_hits, any_antisense, r1_region = cr_transcriptome.align_to_transcriptome(read1, chrom_info, transcript_info, exon_info, annotation_params)
        r1_hits = r1_hits.values()
        for hit in r1_hits:
            if hit.strand == cr_constants.FORWARD_STRAND:
                transcript_length = gene_index.get_transcript_length(hit.transcript)
                insert_size_for_transcript = transcript_length - hit.pos
                if insert_size_for_transcript <= cr_constants.MAX_INSERT_SIZE:
                    transcript_insert_sizes.append(insert_size_for_transcript)

        if len(r1_hits) > 0:
            gene_hits = cr_utils.get_genes_from_transcripts(r1_hits, gene_index)
            all_genes.update(gene_hits)
            r1_tx = cr_transcriptome.make_tx_tag(r1_hits)
            r1_gx, r1_gn = cr_transcriptome.make_gx_tag(gene_hits)
            r1_tags += cr_utils.make_annotation_tags(r1_tx, r1_gx, r1_gn)

    if r1_region: r1_tags.append((cr_constants.MAPPING_REGION_TAG, cr_transcriptome.make_region_tag(r1_region)))
    if r2_region: r2_tags.append((cr_constants.MAPPING_REGION_TAG, cr_transcriptome.make_region_tag(r2_region)))

    return r1_tags, r2_tags, all_genes, transcript_insert_sizes, any_antisense, discordant_pair, r1_region, r2_region

def compute_paired_insert_size(r1_hit, r2_hit):
    assert r1_hit.transcript == r2_hit.transcript
    start = min(r1_hit.pos, r2_hit.pos)
    end = max(r1_hit.pos + r1_hit.alen, r2_hit.pos + r2_hit.alen)
    insert_size_for_transcript = end - start
    return insert_size_for_transcript

def make_barcode_tags(qname, reporter, args):
    gem_group = args.gem_group
    correct_barcodes = args.correct_barcodes
    barcode_confidence_threshold = args.barcode_confidence_threshold
    barcode_whitelist = reporter.barcode_whitelist
    barcode_dist = reporter.barcode_dist

    tags = []
    fastq_header = AugmentedFastqHeader(qname)

    # Barcode tags
    raw_bc_seq = fastq_header.get_tag(cr_constants.RAW_BARCODE_TAG)
    bc_qual = fastq_header.get_tag(cr_constants.RAW_BARCODE_QUAL_TAG)
    barcode_info = None

    if len(raw_bc_seq) > 0:
        processed_bc_seq = reporter.raw_barcode_cb(raw_bc_seq, bc_qual)

        # Add the gem group
        if processed_bc_seq is not None:
            processed_bc_seq = cr_utils.format_barcode_seq(processed_bc_seq, gem_group=gem_group)

        if (processed_bc_seq is None) and (barcode_whitelist is not None):
            if correct_barcodes:
                # Try to correct the barcode
                processed_bc_seq = cr_stats.correct_bc_error(barcode_confidence_threshold, raw_bc_seq, bc_qual, barcode_dist)

                # Add the gem group
                if processed_bc_seq is not None:
                    processed_bc_seq = cr_utils.format_barcode_seq(processed_bc_seq, gem_group=gem_group)

            else:
                # If the barcode was already corrected, take that (gem group is included)
                processed_bc_seq = fastq_header.get_tag(cr_constants.PROCESSED_BARCODE_TAG)
        tags.append((cr_constants.RAW_BARCODE_TAG, raw_bc_seq))
        tags.append((cr_constants.RAW_BARCODE_QUAL_TAG, bc_qual))
        if processed_bc_seq is not None:
            tags.append((cr_constants.PROCESSED_BARCODE_TAG, processed_bc_seq))
        barcode_info = cr_constants.ProcessedRead(raw_bc_seq, processed_bc_seq, bc_qual)

    # UMI tags
    raw_umi_seq = fastq_header.get_tag(cr_constants.RAW_UMI_TAG)
    umi_qual = fastq_header.get_tag(cr_constants.UMI_QUAL_TAG)
    umi_info = None

    if len(raw_umi_seq) > 0:
        processed_umi_seq = reporter.raw_umi_cb(raw_umi_seq, umi_qual)
        tags.append((cr_constants.RAW_UMI_TAG, raw_umi_seq))
        tags.append((cr_constants.UMI_QUAL_TAG, umi_qual))
        if processed_umi_seq is not None:
            tags.append((cr_constants.PROCESSED_UMI_TAG, processed_umi_seq))
        umi_info = cr_constants.ProcessedRead(raw_umi_seq, processed_umi_seq, umi_qual)

    # Sample index tags
    si_seq = fastq_header.get_tag(tk_constants.SAMPLE_INDEX_TAG)
    si_qual = fastq_header.get_tag(tk_constants.SAMPLE_INDEX_QUAL_TAG)

    if len(si_seq) > 0:
        tags.append((tk_constants.SAMPLE_INDEX_TAG, si_seq))
        tags.append((tk_constants.SAMPLE_INDEX_QUAL_TAG, si_qual))

    stripped_qname = fastq_header.fastq_header

    return stripped_qname, tags, barcode_info, umi_info

def make_trimmed_tags(trimmed_read):
    tags = []
    trimmed_seq = trimmed_read.seq
    trimmed_qual = trimmed_read.qual
    if trimmed_seq is not None and len(trimmed_seq) > 0:
        tags.append((cr_constants.TRIMMED_SEQ_TAG, trimmed_seq))
        tags.append((cr_constants.TRIMMED_QUAL_TAG, trimmed_qual))
    return tags
