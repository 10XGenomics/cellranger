#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
import tenkit.bam as tk_bam
import cellranger.chemistry as cr_chem
import cellranger.constants as cr_constants
import cellranger.matrix as cr_matrix
import cellranger.reference as cr_reference
import cellranger.report as cr_report
import cellranger.utils as cr_utils

__MRO__ = """
stage COUNT_GENES(
    in  string sample_id,
    in  bam[]  inputs,
    in  path   reference_path,
    in  map    chemistry_def,
    in  string barcode_whitelist,
    in  h5     barcode_summary,
    in  int[]  gem_groups,
    in  map    align,
    out h5     matrices_h5,
    out path   matrices_mex,
    out pickle chunked_reporter,
    out json   reporter_summary,
    out h5     barcode_summary,
    src py     "stages/counter/count_genes",
) split using (
    in  bam    chunk_input,
)
"""

def split(args):
    chunk_mem_gb = cr_utils.get_mem_gb_request_from_barcode_whitelist(args.barcode_whitelist)

    chunks = []
    for chunk_input in args.inputs:
        chunks.append({
            'chunk_input': chunk_input,
            '__mem_gb': chunk_mem_gb,
        })

    join_mem_gb = cr_utils.get_mem_gb_request_from_barcode_whitelist(args.barcode_whitelist, args.gem_groups, use_min=False)

    # Account for memory used by reporters (particularly the bc and umi diversity dicts)
    genomes = cr_utils.get_reference_genomes(args.reference_path)

    barcode_whitelist = cr_utils.load_barcode_whitelist(args.barcode_whitelist)
    if barcode_whitelist is not None:
        num_barcodes = len(barcode_whitelist) * max(args.gem_groups)
    else:
        num_barcodes = cr_utils.get_num_barcodes_from_barcode_summary(args.barcode_summary)

    max_bc_diversity_entries = num_barcodes
    max_umi_diversity_entries = 4 ** cr_chem.get_umi_length(args.chemistry_def)

    # Multiply by 2 to hold the current reporter + accumulating reporter in the merge
    bc_diversity_mem_gb = (2 * max_bc_diversity_entries * cr_constants.BYTES_PER_STR_INT_DICT_ENTRY * (len(genomes) + 1) * len(cr_constants.READ_TYPES))/1e9
    umi_diversity_mem_gb = (2 * max_umi_diversity_entries * cr_constants.BYTES_PER_STR_INT_DICT_ENTRY * (len(genomes) + 1) * len(cr_constants.READ_TYPES))/1e9
    join_mem_gb = min(cr_constants.COUNT_GENES_MAX_MEM_GB, max(cr_constants.MIN_MEM_GB, int(join_mem_gb + bc_diversity_mem_gb + umi_diversity_mem_gb)))
    join = {
        '__mem_gb': join_mem_gb,
    }
    return {'chunks': chunks, 'join': join}

def join_matrices(args, outs, chunk_defs, chunk_outs):
    chunk_h5s = [chunk_out.matrices_h5 for chunk_out in chunk_outs]
    matrices = cr_matrix.merge_matrices(chunk_h5s)
    matrix_attrs = cr_matrix.make_matrix_attrs_count(args.sample_id, args.gem_groups, cr_chem.get_description(args.chemistry_def))
    matrices.save_h5(outs.matrices_h5, extra_attrs=matrix_attrs)
    matrices.save_mex(outs.matrices_mex)

def join_reporter(args, outs, chunk_defs, chunk_outs):
    outs.chunked_reporter = None
    reporter = cr_report.merge_reporters([chunk_out.chunked_reporter for chunk_out in chunk_outs])

    reporter.report_summary_json(outs.reporter_summary)
    reporter.report_barcodes_h5(outs.barcode_summary)

def join(args, outs, chunk_defs, chunk_outs):
    join_reporter(args, outs, chunk_defs, chunk_outs)
    join_matrices(args, outs, chunk_defs, chunk_outs)

def main(args, outs):
    in_bam = tk_bam.create_bam_infile(args.chunk_input)

    chroms = in_bam.references

    barcode_whitelist = cr_utils.load_barcode_whitelist(args.barcode_whitelist)
    barcode_summary = cr_utils.load_barcode_summary(args.barcode_summary) if not barcode_whitelist else None

    gene_index = cr_reference.GeneIndex.load_pickle(cr_utils.get_reference_genes_index(args.reference_path))
    reporter = cr_report.Reporter(reference_path=args.reference_path,
                                  high_conf_mapq=cr_utils.get_high_conf_mapq(args.align),
                                  gene_index=gene_index,
                                  chroms=chroms,
                                  barcode_whitelist=barcode_whitelist,
                                  barcode_summary=barcode_summary,
                                  gem_groups=args.gem_groups)

    if barcode_whitelist:
        barcode_seqs = cr_utils.format_barcode_seqs(barcode_whitelist, args.gem_groups)
    else:
        barcode_seqs = barcode_summary

    genomes = cr_utils.get_reference_genomes(args.reference_path)
    genes = cr_utils.split_genes_by_genomes(gene_index.get_genes(), genomes)
    matrices = cr_matrix.GeneBCMatrices(genomes, genes, barcode_seqs)

    for read in in_bam:
        is_conf_mapped_deduped, genome, gene_id, bc = reporter.count_genes_bam_cb(read,
                                                                                  use_umis=cr_chem.has_umis(args.chemistry_def))
        if is_conf_mapped_deduped:
            matrices.add(genome, gene_id, bc)

    in_bam.close()

    matrices.save_h5(outs.matrices_h5)
    reporter.save(outs.chunked_reporter)
