#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

from __future__ import annotations

import json
import shutil

import h5py as h5
import martian
from six import ensure_binary

import cellranger.analysis.io as analysis_io
import cellranger.h5_constants as h5_constants
import cellranger.matrix as cr_matrix
import cellranger.rna.matrix as rna_matrix
import cellranger.websummary.sample_properties as wsp
import tenkit.safe_json as tk_json
from cellranger.webshim.constants.shared import PIPELINE_REANALYZE
from cellranger.websummary import web_summary_builder as wsb

__MRO__ = """
stage SUMMARIZE_REANALYSIS(
    in  string sample_id,
    in  string sample_desc,
    in  h5     filtered_matrices,
    in  path   analysis,
    in  json   analyze_matrices_summary,
    in  json   antibody_histograms,
    in  json   antibody_treemap,
    out html   web_summary,
    out json   summary,
    out path   feature_bc_matrix_mex,
    src py     "stages/analyzer/summarize_reanalysis",
) split using (
)
"""

SUMMARIZE_REANALYSIS_MIN_MEM_GB = 3.0


def split(args):
    if args.analysis:
        # Estimate memory usage from the matrix stored in the analysis h5
        h5_path = analysis_io.h5_path(args.analysis)
        with h5.File(ensure_binary(h5_path), "r") as f:
            matrix_mem_gb = cr_matrix.CountMatrix.get_mem_gb_from_group(f["matrix"])
    else:
        matrix_mem_gb = h5_constants.MIN_MEM_GB

    # The join stage memory was previously hard coded to 3 GB , the stage was modified to output the
    # MEX file of the matrix to, which is a task also performed by the FILTER_BARCODES stage
    # I duplicated the logic for that memory requirement here, so we could evolve it later
    # without affecting the original
    join_mem_gb = 2 * cr_matrix.CountMatrix.get_mem_gb_from_matrix_h5(args.filtered_matrices)
    join_mem_gb = max(join_mem_gb, SUMMARIZE_REANALYSIS_MIN_MEM_GB)
    chunks = [
        {
            "__mem_gb": matrix_mem_gb,
        }
    ]
    return {"chunks": chunks, "join": {"__mem_gb": join_mem_gb}}


def main(args, outs):
    genomes = cr_matrix.CountMatrix.get_genomes_from_h5(args.filtered_matrices)
    chemistry = cr_matrix.CountMatrix.load_chemistry_from_h5(args.filtered_matrices)
    total_cells = cr_matrix.CountMatrix.count_cells_from_h5(args.filtered_matrices)

    summary = {"chemistry_description": chemistry, "filtered_bcs_transcriptome_union": total_cells}
    if args.analyze_matrices_summary:
        with open(args.analyze_matrices_summary) as reader:
            analysis_summary = json.load(reader)
        summary.update(analysis_summary)

    with open(outs.summary, "w") as f:
        tk_json.dump_numpy(summary, f, indent=4, sort_keys=True)

    sample_properties = wsp.CountSampleProperties(
        sample_id=args.sample_id, sample_desc=args.sample_desc, genomes=genomes
    )

    sample_data_paths = wsp.SampleDataPaths(
        summary_path=outs.summary,
        analysis_path=args.analysis,
        antibody_histograms_path=args.antibody_histograms,
        antibody_treemap_path=args.antibody_treemap,
    )

    _ = wsb.build_web_summary_html_sc(
        outs.web_summary, sample_properties, sample_data_paths, PIPELINE_REANALYZE
    )


def join(args, outs, chunk_defs, chunk_outs):
    chunk_out = chunk_outs[0]

    shutil.copy(chunk_out.web_summary, outs.web_summary)
    shutil.copy(chunk_out.summary, outs.summary)

    matrix = cr_matrix.CountMatrix.load_h5_file(args.filtered_matrices)
    rna_matrix.save_mex(matrix, outs.feature_bc_matrix_mex, martian.get_pipelines_version())
