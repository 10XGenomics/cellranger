#
# Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
#
"""Functions for the stage code BUILD_SAMPLE_OUTS and BUILD _SAMPLE_OUTS_PD."""
from cellranger.cr_io import hard_link

SAMPLE_OUTS_COUNT = "count"
ANTIGEN_ANALYSIS = "antigen_analysis"


def build_sample_outs(args, outs, is_pd: bool):
    """Builds the sample outs for code that is shared between PD and CS."""
    count = {}
    analysis_key = "analysis_csv" if is_pd else "analysis"

    count[analysis_key] = hard_link(args.rna_analysis)
    count["sample_cloupe"] = hard_link(args.cloupe)
    count["crispr_analysis"] = hard_link(args.crispr_analysis)
    if hasattr(args, "analysis"):
        count["analysis"] = hard_link(args.analysis)

    # Add cell_types if available
    if hasattr(args, "cell_types") and args.cell_types is not None:
        count["cell_types"] = {}
        for key, in_path in args.cell_types.items():
            count["cell_types"][key] = hard_link(in_path)
    else:
        count["cell_types"] = None

    def hard_link_sample_slfe_out(in_key: str, out_key: str, include_if: bool = True):
        in_path = (args.sample_slfe_outs or {}).get(in_key, None)
        if in_path is None or not include_if:
            count[out_key] = None
            return
        count[out_key] = hard_link(in_path)

    hard_link_sample_slfe_out("feature_reference", "feature_reference_csv")
    hard_link_sample_slfe_out("bam_file", "sample_alignments")
    hard_link_sample_slfe_out("bai_index_file", "sample_alignments_index_bai")
    hard_link_sample_slfe_out("csi_index_file", "sample_alignments_index_csi")
    hard_link_sample_slfe_out("filtered_barcodes", "sample_filtered_barcodes_csv")
    hard_link_sample_slfe_out("aggregate_barcodes", "aggregate_barcodes")
    hard_link_sample_slfe_out("filtered_matrix_h5", "sample_filtered_feature_bc_matrix")
    hard_link_sample_slfe_out("filtered_matrix_mex", "sample_filtered_feature_bc_matrix_mex")
    hard_link_sample_slfe_out(
        "raw_matrix_h5",
        "sample_raw_feature_bc_matrix",
        include_if=args.output_per_sample_raw_matrix,
    )
    hard_link_sample_slfe_out(
        "raw_matrix_mex",
        "sample_raw_feature_bc_matrix_mex",
        include_if=args.output_per_sample_raw_matrix,
    )
    hard_link_sample_slfe_out("molecule_info", "sample_molecule_info")
    hard_link_sample_slfe_out("raw_probe_bc_matrix", "sample_raw_probe_bc_matrix")
    hard_link_sample_slfe_out("target_panel", "target_panel")
    if is_pd:
        hard_link_sample_slfe_out("per_probe_metrics", "sample_per_probe_metrics")
    else:
        hard_link_sample_slfe_out("probe_set", "probe_set")

    sample_outs = {}
    sample_outs[SAMPLE_OUTS_COUNT] = count
    sample_outs["metrics_summary"] = hard_link(args.metrics_summary_csv)
    # Note that although this assignment is the same, the values being assigned differ between PD and CS code.
    sample_outs["vdj_b"] = args.vdj_b_outs
    sample_outs["vdj_t"] = args.vdj_t_outs
    sample_outs["vdj_t_gd"] = args.vdj_t_gd_outs
    sample_outs["web_summary"] = hard_link(args.web_summary)

    if args.beam_analyzer is not None:
        sample_outs["antigen_analysis"] = {}
        for key in ["antigen_specificity_scores", "per_barcode"]:
            sample_outs["antigen_analysis"][key] = hard_link(args.beam_analyzer.get(key))
    else:
        sample_outs["antigen_analysis"] = None

    outs.sample_outs = sample_outs
