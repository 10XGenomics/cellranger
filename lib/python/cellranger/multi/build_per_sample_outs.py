#
# Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
#
"""Functions for the stage code BUILD_SAMPLE_OUTS and BUILD _SAMPLE_OUTS_PD."""
from cellranger.cr_io import hard_link, recursive_hard_link_dict


def build_sample_outs(args, outs):
    """Builds the sample outs for code that is shared between PD and CS."""
    # TODO: type this method by defining a common input and output type for BUILD_SAMPLE_OUTS and BUILD_SAMPLE_OUTS_PD
    sample_outs = {}
    sample_outs["analysis"] = hard_link(args.rna_analysis)
    sample_outs["crispr_analysis"] = hard_link(args.crispr_analysis)
    sample_outs["metrics_summary"] = hard_link(args.metrics_summary_csv)
    sample_outs["sample_cloupe"] = hard_link(args.cloupe)
    sample_outs["vdj_b"] = args.vdj_b_outs
    sample_outs["vdj_t_gd"] = args.vdj_t_gd_outs
    sample_outs["vdj_t"] = args.vdj_t_outs
    sample_outs["web_summary"] = hard_link(args.web_summary)

    sample_outs["cell_types"] = recursive_hard_link_dict(args.cell_types)
    if sample_outs["cell_types"] and isinstance(sample_outs["cell_types"].get("cell_types"), dict):
        sample_outs["cell_types"] = sample_outs["cell_types"].get("cell_types")

    sample_outs.update(
        (
            dest,
            hard_link(args.sample_slfe_outs[source]) if args.sample_slfe_outs is not None else None,
        )
        for dest, source in (
            ("aggregate_barcodes", "aggregate_barcodes"),
            ("sample_alignments_index_bai", "bai_index_file"),
            ("sample_alignments_index_csi", "csi_index_file"),
            ("sample_alignments", "bam_file"),
            ("sample_filtered_barcodes_csv", "filtered_barcodes"),
            ("sample_filtered_feature_bc_matrix_mex", "filtered_matrix_mex"),
            ("sample_filtered_feature_bc_matrix", "filtered_matrix_h5"),
            ("sample_molecule_info", "molecule_info"),
            ("sample_raw_probe_bc_matrix", "raw_probe_bc_matrix"),
        )
    )

    sample_outs.update(
        (
            dest,
            (
                hard_link(args.sample_slfe_outs[source])
                if args.output_per_sample_raw_matrix is not None
                else None
            ),
        )
        for dest, source in (
            ("sample_raw_feature_bc_matrix", "raw_matrix_h5"),
            ("sample_raw_feature_bc_matrix_mex", "raw_matrix_mex"),
        )
    )

    sample_outs["antigen_analysis"] = (
        {
            key: hard_link(args.beam_analyzer[key])
            for key in ["antigen_specificity_scores", "per_barcode"]
        }
        if args.beam_analyzer is not None
        else None
    )

    outs.sample_outs = sample_outs
