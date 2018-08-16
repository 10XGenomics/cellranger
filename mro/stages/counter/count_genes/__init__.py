#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
import martian
import tenkit.bam as tk_bam
import cellranger.chemistry as cr_chem
import cellranger.constants as cr_constants
import cellranger.rna.feature_ref as rna_feature_ref
import cellranger.matrix as cr_matrix
import cellranger.reference as cr_reference
import cellranger.rna.library as rna_library
import cellranger.rna.matrix as rna_matrix
import cellranger.report as cr_report
import cellranger.utils as cr_utils

__MRO__ = """
stage COUNT_GENES(
    in  string sample_id,
    in  bam[]  inputs,
    in  path   reference_path,
    in  csv    feature_reference,
    in  map    chemistry_def,
    in  string barcode_whitelist,
    in  csv    barcodes_detected,
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
    chunks = []
    for chunk_input in args.inputs:
        chunks.append({
            'chunk_input': chunk_input,
            '__mem_gb': 8,
        })

    # FIXME: Need two copies of the matrix, but we don't have matrices yet
    #        Consider using info from MARK_DUPLICATES and the whitelists
    #          to estimate the matrix memory size.
    join = {
        '__mem_gb': 12,
    }
    return {'chunks': chunks, 'join': join}

def join_matrices(args, outs, chunk_defs, chunk_outs):
    chunk_h5s = [chunk_out.matrices_h5 for chunk_out in chunk_outs]
    matrix = cr_matrix.merge_matrices(chunk_h5s)
    matrix_attrs = cr_matrix.make_matrix_attrs_count(args.sample_id, args.gem_groups, cr_chem.get_description(args.chemistry_def))
    matrix.save_h5_file(outs.matrices_h5, extra_attrs=matrix_attrs)

    rna_matrix.save_mex(matrix,
                        outs.matrices_mex,
                        martian.get_pipelines_version())

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

    libraries = rna_library.get_bam_library_info(in_bam)
    distinct_library_types = sorted(list(set([x['library_type'] for x in libraries])))
    library_prefixes = map(lambda lib: rna_library.get_library_type_metric_prefix(lib['library_type']),
                           libraries)

    chroms = in_bam.references

    barcode_whitelist = cr_utils.load_barcode_whitelist(args.barcode_whitelist)
    barcode_summary = cr_utils.load_barcode_tsv(args.barcodes_detected) if not barcode_whitelist else None

    # TODO: this is redundant
    gene_index = cr_reference.GeneIndex.load_pickle(cr_utils.get_reference_genes_index(args.reference_path))
    reporter = cr_report.Reporter(reference_path=args.reference_path,
                                  high_conf_mapq=cr_utils.get_high_conf_mapq(args.align),
                                  gene_index=gene_index,
                                  chroms=chroms,
                                  barcode_whitelist=barcode_whitelist,
                                  barcode_summary=barcode_summary,
                                  gem_groups=args.gem_groups,
                                  library_types=distinct_library_types)

    feature_ref = rna_feature_ref.from_transcriptome_and_csv(args.reference_path,
                                                             args.feature_reference)

    if barcode_whitelist:
        barcode_seqs = cr_utils.format_barcode_seqs(barcode_whitelist, args.gem_groups)
    else:
        barcode_seqs = barcode_summary

    matrix = cr_matrix.CountMatrix.empty(feature_ref, barcode_seqs, dtype='int32')

    for qname, reads_iter, _ in cr_utils.iter_by_qname(in_bam, None):
        is_conf_mapped_deduped, genome, feature_id, bc = reporter.count_genes_bam_cb(reads_iter,
                                                                                     libraries,
                                                                                     library_prefixes,
                                                                                     use_umis=cr_chem.has_umis(args.chemistry_def))
        if is_conf_mapped_deduped:
            matrix.add(feature_id, bc)

    in_bam.close()

    reporter.store_reference_metadata(args.reference_path, cr_constants.REFERENCE_TYPE, cr_constants.REFERENCE_METRIC_PREFIX)

    matrix.save_h5_file(outs.matrices_h5)
    reporter.save(outs.chunked_reporter)
