#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#
import json
import martian
import tables
import cellranger.analysis.io as cr_io
import cellranger.constants as cr_constants
import cellranger.matrix as cr_matrix
import cellranger.utils as cr_utils
import cellranger.webshim.common as cr_webshim
import cellranger.webshim.data as cr_webshim_data
from cellranger.webshim.constants.shared import PIPELINE_REANALYZE

__MRO__ = """
stage SUMMARIZE_REANALYSIS(
    in  string analysis_id,
    in  string analysis_desc,
    in  h5     filtered_matrices,
    in  path   analysis,
    out html   web_summary,
    out json   summary,
    src py     "stages/analyzer/summarize_reanalysis",
) split using (
)
"""
def split(args):
    # Estimate memory usage from the matrix stored in the analysis h5
    if args.analysis:
        with tables.open_file(cr_io.h5_path(args.analysis), 'r') as f:
            matrix = getattr(f.root, cr_constants.ANALYSIS_H5_MATRIX_GROUP)
            matrix_mem_gb = cr_matrix.GeneBCMatrix.get_mem_gb_from_group(matrix)
    else:
        matrix_mem_gb = cr_constants.MIN_MEM_GB

    chunks = [{
        '__mem_gb': matrix_mem_gb,
    }]
    return {'chunks': chunks}

def main(args, outs):
    genomes = cr_matrix.GeneBCMatrices.load_genomes_from_h5(args.filtered_matrices)
    chemistry = cr_matrix.GeneBCMatrices.load_chemistry_from_h5(args.filtered_matrices)
    total_cells = cr_matrix.GeneBCMatrices.count_cells_from_h5(args.filtered_matrices)
    summary = {'chemistry_description': chemistry, 'filtered_bcs_transcriptome_union': total_cells}
    with open(outs.summary, 'w') as f:
        json.dump(summary, f, indent=4, sort_keys=True)

    sample_properties = cr_webshim.get_sample_properties(args.analysis_id, args.analysis_desc, genomes, version=martian.get_pipelines_version())

    sample_data_paths = cr_webshim_data.SampleDataPaths(
        summary_path=outs.summary,
        analysis_path=args.analysis,
    )

    sample_data = cr_webshim.load_sample_data(sample_properties, sample_data_paths)
    cr_webshim.build_web_summary_html(outs.web_summary, sample_properties, sample_data, PIPELINE_REANALYZE)

def join(args, outs, chunk_defs, chunk_outs):
    chunk_out = chunk_outs[0]

    cr_utils.copy(chunk_out.web_summary, outs.web_summary)
    cr_utils.copy(chunk_out.summary, outs.summary)
