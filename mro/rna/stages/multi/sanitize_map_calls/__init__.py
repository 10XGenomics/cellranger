#!/usr/bin/env python
#
# Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
#
"""A helper stage to pipe through map called stages to sidestep martian issue."""

import cellranger.cr_io as cr_io
from cellranger.fast_utils import MultiGraph

__MRO__ = """
stage SANITIZE_MAP_CALLS(
    in  json              multi_graph,
    in  map<path>         in_crispr_analysis,
    in  map<path>         in_rna_analysis,
    in  map<cloupe>       in_cloupe_file,
    in  map<json>         in_metrics_summary,
    in  map<json>         in_sample_projection_plots,
    in  map<json>         in_sample_barcode_rank_plots,
    in  map<json>         in_sample_treemap_plots,
    in  map<VDJ_ANALYZER> in_vdj_t_analyzer,
    in  map<VDJ_ANALYZER> in_vdj_t_gd_analyzer,
    in  map<VDJ_ANALYZER> in_vdj_b_analyzer,
    out map<path>         crispr_analysis,
    out map<path>         rna_analysis,
    out map<cloupe>       cloupe_file,
    out map<json>         metrics_summary,
    out map<json>         sample_projection_plots,
    out map<json>         sample_barcode_rank_plots,
    out map<json>         sample_treemap_plots,
    out map<VDJ_ANALYZER> vdj_t_analyzer,
    out map<VDJ_ANALYZER> vdj_t_gd_analyzer,
    out map<VDJ_ANALYZER> vdj_b_analyzer,
    src py                "stages/multi/sanitize_map_calls",
) using (
    volatile = false,
) retain (
    metrics_summary,
"""


# pylint: disable=too-few-public-methods,too-many-instance-attributes
class VdjAnalyzerClonotypeOuts:
    """Python equivalent of martian BeamAnalyzerOutputs struct."""

    def __init__(self):
        """Initiate an empty martian struct."""
        self.airr_rearrangement = None
        self.all_contig_annotations_csv = None
        self.all_contig_annotations_json = None
        self.clonotypes_csv = None
        self.concat_ref_bam = None
        self.concat_ref_bam_bai = None
        self.concat_ref_fasta = None
        self.concat_ref_fasta_fai = None
        self.consensus_annotations_csv = None
        self.consensus_bam = None
        self.consensus_bam_bai = None
        self.consensus_fasta = None
        self.consensus_fasta_fai = None
        self.donor_ref_fa = None
        self.enclone_output = None
        self.filtered_contig_annotations_csv = None


# pylint: disable=too-few-public-methods
class BeamAnalyzerOutputs:
    """Python equivalent of martian BeamAnalyzerOutputs struct."""

    def __init__(self):
        """Initiate an empty martian struct."""
        self.antigen_specificity_scores = None
        self.antigen_assignment = None
        self.clonotype_concordance = None
        self.exact_subclonotype_concordance = None
        self.specificity_summary = None
        self.antigen_vdj_metrics_json = None
        self.antigen_vdj_metrics_bin = None
        self.per_barcode = None


# pylint: disable=too-few-public-methods,too-many-instance-attributes
class VdjReport:
    """Python equivalent of martian VdjReport struct."""

    def __init__(self):
        """Initiate an empty martian struct."""
        self.vdj_contig_info = None
        self.vloupe = None
        self.metrics_summary_json = None
        self.metrics_summary_csv = None
        self.web_summary = None
        self.web_summary_data = None
        self.contig_fastq = None
        self.filtered_contig_fastq = None
        self.contig_fasta = None
        self.contig_fasta_fai = None
        self.filtered_contig_fasta = None
        self.annotations_bed = None
        self.cell_barcodes = None
        self.cdr3_barcodes = None
        self.all_contig_barcodes = None
        self.productive_barcodes = None
        self.productive_cell_barcodes = None
        self.filter_summary = None
        self.filter_metrics = None
        self.per_bc_filters = None
        self.umi_summary = None
        self.barcode_brief = None
        self.report = None


# pylint: disable=invalid-name
class VDJ_ANALYZER:
    """Python equivalent of martian VDJ_ANALYZER struct."""

    def __init__(self):
        """Initiate an empty martian struct."""
        self.clonotype = VdjAnalyzerClonotypeOuts().__dict__
        self.beam_analyzer = BeamAnalyzerOutputs().__dict__
        self.report = VdjReport().__dict__


def main(args, outs):
    outs.rna_analysis = cr_io.recursive_hard_link_dict(args.in_rna_analysis)
    outs.crispr_analysis = cr_io.recursive_hard_link_dict(args.in_crispr_analysis)
    outs.cloupe_file = cr_io.recursive_hard_link_dict(args.in_cloupe_file)
    outs.metrics_summary = cr_io.recursive_hard_link_dict(args.in_metrics_summary)
    outs.sample_projection_plots = cr_io.recursive_hard_link_dict(args.in_sample_projection_plots)
    outs.sample_barcode_rank_plots = cr_io.recursive_hard_link_dict(
        args.in_sample_barcode_rank_plots
    )
    outs.sample_treemap_plots = cr_io.recursive_hard_link_dict(args.in_sample_treemap_plots)

    if args.multi_graph:
        config = MultiGraph.from_path(args.multi_graph)
        samples = config.sample_ids()

        if args.in_vdj_t_analyzer is not None:
            outs.vdj_t_analyzer = cr_io.recursive_hard_link_dict(
                args.in_vdj_t_analyzer, prefixes=["vdj_t"]
            )
        else:
            outs.vdj_t_analyzer = {k: VDJ_ANALYZER().__dict__ for k in samples}
        if args.in_vdj_t_gd_analyzer is not None:
            outs.vdj_t_gd_analyzer = cr_io.recursive_hard_link_dict(
                args.in_vdj_t_gd_analyzer, prefixes=["vdj_t_gd"]
            )
        else:
            outs.vdj_t_gd_analyzer = {k: VDJ_ANALYZER().__dict__ for k in samples}
        if args.in_vdj_b_analyzer is not None:
            outs.vdj_b_analyzer = cr_io.recursive_hard_link_dict(
                args.in_vdj_b_analyzer, prefixes=["vdj_b"]
            )
        else:
            outs.vdj_b_analyzer = {k: VDJ_ANALYZER().__dict__ for k in samples}
