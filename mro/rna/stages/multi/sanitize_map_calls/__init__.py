#!/usr/bin/env python
#
# Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
#
"""A helper stage to pipe through map called stages to sidestep martian issue."""


import os

import martian

import cellranger.cr_io as cr_io

__MRO__ = """
stage SANITIZE_MAP_CALLS(
    in  map<path>   in_crispr_analysis,
    in  map<path>   in_rna_analysis,
    in  map<cloupe> in_cloupe_file,
    in  map<json>   in_metrics_summary,
    in  map<json>   in_sample_tsne_plots,
    in  map<json>   in_sample_barcode_rank_plots,
    in  map<json>   in_sample_treemap_plots,
    out map<path>   crispr_analysis,
    out map<path>   rna_analysis,
    out map<cloupe> cloupe_file,
    out map<json>   metrics_summary,
    out map<json>   sample_tsne_plots,
    out map<json>   sample_barcode_rank_plots,
    out map<json>   sample_treemap_plots,
    src py          "stages/multi/sanitize_map_calls",
) using (
    volatile = strict,
)
"""


def recursive_hard_link_dict(in_files, prefixes=None):
    """Hard link files into this stage directory.

    For a dict with type [String,path_or_file_or_dict], where the keys represent sample IDs
    or other levels of nesting, create a dict with the same keys where all the path_or_file
    values are hardlinked into this stage directory from their old paths.
    If path_or_file_or_dict is a directory, recursively hardlink the contents of the directory.
    When nesting the key used to access each upper level will be used as a prefix
    when constructing the final file name.

    For example,

      {
        "sample1": {
          "gene_expression": "/mydir/gex.txt",
          "antibody": "/mydir/ab.txt"
        },
        "sample2": {
          "gene_expression": "/mydir/gex.txt",
          "antibody": "/mydir/ab.txt"
        },
      }

    would become

      {
        "sample1": {
          "gene_expression": "sample1_gene_expression_gex.txt",
          "antibody": "sample1_antibody.ab.txt"
         },
         "sample2": {
           "gene_expression": "sample2_gene_expression_gex.txt",
           "antibody": "sample2_antibody_ab.txt"
         },
      }
    """
    if in_files is None:
        return None

    if prefixes is None:
        prefixes = []

    out_files = {}
    for k, path_or_dict in in_files.items():
        if path_or_dict is None:
            out_files[k] = None
        else:
            if isinstance(path_or_dict, dict):
                out_files[k] = recursive_hard_link_dict(in_files[k], prefixes + [k])
            elif isinstance(path_or_dict, (str)):
                final_prefixes = prefixes + [k]
                old_path = path_or_dict
                new_path = martian.make_path(
                    "_".join(final_prefixes) + "_" + os.path.basename(old_path)
                ).decode("utf8")
                cr_io.hardlink_with_fallback(old_path, new_path)
                out_files[k] = new_path
            else:
                raise ValueError(
                    "Input dictionary may not contain any elements other than dict and string: %s"
                    % path_or_dict
                )
    return out_files


def main(args, outs):
    outs.rna_analysis = recursive_hard_link_dict(args.in_rna_analysis)
    outs.crispr_analysis = recursive_hard_link_dict(args.in_crispr_analysis)
    outs.cloupe_file = recursive_hard_link_dict(args.in_cloupe_file)
    outs.metrics_summary = recursive_hard_link_dict(args.in_metrics_summary)
    outs.sample_tsne_plots = recursive_hard_link_dict(args.in_sample_tsne_plots)
    outs.sample_barcode_rank_plots = recursive_hard_link_dict(args.in_sample_barcode_rank_plots)
    outs.sample_treemap_plots = recursive_hard_link_dict(args.in_sample_treemap_plots)
