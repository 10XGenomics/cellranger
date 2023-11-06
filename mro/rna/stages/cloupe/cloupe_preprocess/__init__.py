#!/usr/bin/env python
#
# Copyright (c) 2018 10X Genomics, Inc. All rights reserved.
#


import math
import os
import shutil
import subprocess

import martian

import cellranger.barcodes.utils as bc_utils
import cellranger.constants as cr_constants
import cellranger.h5_constants as h5_constants
import cellranger.matrix as cr_matrix
import tenkit.log_subprocess as tk_subproc
import tenkit.safe_json as tk_json
from cellranger.spatial.data_utils import (
    DARK_IMAGES_CHANNELS,
    DARK_IMAGES_COLORIZED,
    tissue_position_float_to_int,
)

__MRO__ = """

stage CLOUPE_PREPROCESS(
    in  string   pipestance_type,
    in  string   sample_id,
    in  string   sample_desc,
    in  path     analysis,
    in  h5       filtered_gene_bc_matrices_h5,
    in  json     metrics_json,
    in  csv      aggregation_csv,
    in  json     gem_group_index_json,
    in  string[] image_page_names,
    in  file[]   tissue_image_paths,
    in  int      dark_images,
    in  csv      tissue_positions,
    in  txt      fiducial_positions_list,
    in  json     dzi_info,
    in  path[]   dzi_tiles_paths,
    in  json     scale_factors_json,
    in  bool     no_secondary_analysis,
    in  string   barcode_whitelist,
    in  string   hd_slide_name,
    in  json     loupe_map,
    in  string   product_type,
    in  json     cells_per_sample,
    in  json     cells_per_tag,
    in  json     cells_per_protospacer,
    in  csv      spatial_enrichment,
    in  path     spatial_deconvolution_path,
    out cloupe   output_for_cloupe,
    out json     gem_group_index_json,
    src py       "stages/cloupe/cloupe_preprocess",
) split (
) using (
    volatile = strict,
)
"""


def get_gem_group_index_json(args, outs):
    if args.gem_group_index_json:
        shutil.copy(args.gem_group_index_json, outs.gem_group_index_json)
        return outs.gem_group_index_json
    else:
        generated_index = cr_matrix.get_gem_group_index(args.filtered_gene_bc_matrices_h5)
        if generated_index:
            with open(outs.gem_group_index_json, "w") as outfile:
                tk_json.dump_numpy({"gem_group_index": generated_index}, outfile)
            return outs.gem_group_index_json
        else:
            return None


def get_analysis_h5_path(args):
    return os.path.join(args.analysis, "analysis.h5")


def do_not_make_cloupe(args):
    """Returns True if there is a reason why this stage should not attempt to.

    generate a .cloupe file.
    #General:
        Case1: no_secondary_analysis is True
        Case2: No analysis folder, can't make a cloupe file
        Case3: No filtered matrix. Is there a case where an analysis folder exists but no filtered gene matrix?
    #Spatial specific
        Single sample cases:
            Case4: If Tissue positions list present but missing cloupe implementation of barcodes
            Case5: If not Tissue positions list but SPATIAL Pipestance. Image crashed but the rest worked
    """
    # General cases
    # Case1
    if args.no_secondary_analysis:
        martian.log_info("Skipping .cloupe generation by instruction (--no-secondary-analysis)")
        return True
    # Case2
    if args.analysis is None:
        martian.log_info("Skipping .cloupe generation due to missing analysis folder")
        return True
    # Case3
    if not os.path.exists(args.filtered_gene_bc_matrices_h5):
        martian.log_info(
            "Skipping .cloupe generation due to missing or zero-length feature-barcode matrix"
        )
        return True
    # Spatial specific
    if args.product_type == cr_constants.SPATIAL_PRODUCT_TYPE:
        # Single sample cases
        if args.aggregation_csv is None:
            # Case4
            if not (
                bc_utils.is_whitelist_spatial(args.barcode_whitelist)
                or args.hd_slide_name is not None
            ):
                martian.log_info("Skipping .cloupe generation due to unsupported barcode whitelist")
                return True

            # Case5
            if not args.tissue_positions:
                martian.log_info(
                    "Skipping .cloupe generation due to spatial pipeline with no image"
                )
                return True

    return False


def split(args):
    # no mem usage if skipped
    if do_not_make_cloupe(args):
        return {"chunks": [{"__mem_gb": h5_constants.MIN_MEM_GB}]}

    # See PR #1196 for a discussion of
    # Where this memory estimate is derived from.  Note that analysis did not account for imaging data.

    mem_gb_matrix = cr_matrix.CountMatrix.get_mem_gb_crconverter_estimate_from_h5(
        args.filtered_gene_bc_matrices_h5
    )

    chunks = [
        {
            "__mem_gb": math.ceil(max(mem_gb_matrix, h5_constants.MIN_MEM_GB)),
            "__threads": 2,
        }
    ]
    # Hack to increase min memgb - not sure why we need it
    return {"chunks": chunks, "join": {"__mem_gb": 1}}


def join(args, outs, chunk_defs, chunk_outs):
    if chunk_outs[0].output_for_cloupe is None:
        # Set output to null if noloupe is set
        outs.output_for_cloupe = None
    else:
        shutil.copy(chunk_outs[0].output_for_cloupe, outs.output_for_cloupe)


def main(args, outs):
    if do_not_make_cloupe(args):
        outs.output_for_cloupe = None
        return

    call = [
        "crconverter",
        args.sample_id,
        args.pipestance_type,
        "--matrix",
        args.filtered_gene_bc_matrices_h5,
        "--analysis",
        get_analysis_h5_path(args),
        "--output",
        outs.output_for_cloupe,
        "--description",
        args.sample_desc,
    ]

    if args.metrics_json:
        call.extend(["--metrics", args.metrics_json])
    if args.aggregation_csv:
        call.extend(["--aggregation", args.aggregation_csv])
        # Only give loupemap if we are aggregating and the sample is spatial
        if args.product_type == cr_constants.SPATIAL_PRODUCT_TYPE:
            call.extend(["--loupemap", args.loupe_map])
    # assume whole thing if tissue positions present
    if args.tissue_positions:
        # Revert after CLOUPE-3789 is merged
        tissue_pos_file = martian.make_path("tissue_positions.csv").decode()
        tissue_position_float_to_int(args.tissue_positions, tissue_pos_file)

        if args.dark_images == DARK_IMAGES_CHANNELS:
            spatial_image_type = "fluorescent"
        elif args.dark_images == DARK_IMAGES_COLORIZED:
            spatial_image_type = "merged"
        else:
            spatial_image_type = "brightfield"
        call.extend(
            [
                "--spatial-image-path",
                args.tissue_image_paths[0],
                "--spatial-tissue-path",
                tissue_pos_file,  # Revert after CLOUPE-3789 is merged
                "--spatial-dzi-path",
                args.dzi_info,
                "--spatial-tiles-paths",
                ",".join(args.dzi_tiles_paths),
                "--spatial-fiducials-path",
                args.fiducial_positions_list,
                "--spatial-scalefactors-path",
                args.scale_factors_json,
                "--spatial-image-type",
                spatial_image_type,
            ]
        )

    if args.image_page_names:
        call.extend(["--spatial-image-page-names", (",").join(args.image_page_names)])

    if args.cells_per_sample:
        call.extend(["--cells-per-sample", args.cells_per_sample])
    if args.cells_per_tag:
        call.extend(["--cells-per-tag", args.cells_per_tag])
    if args.cells_per_protospacer:
        call.extend(["--cells-per-protospacer", args.cells_per_protospacer])
    if args.spatial_enrichment:
        call.extend(["--spatial-enrichment", args.spatial_enrichment])
    if args.spatial_deconvolution_path:
        call.extend(["--spatial-deconvolution-path", args.spatial_deconvolution_path])

    gem_group_index_json = get_gem_group_index_json(args, outs)
    if gem_group_index_json:
        call.extend(["--gemgroups", gem_group_index_json])
    else:
        # required argument for crconverter
        martian.exit("HDF5 matrix to be used for cloupe does not have GEM group information.")

    # the sample desc may be unicode, so send the whole
    # set of args str utf-8 to check_output
    print(call)
    unicode_call = [arg.encode("utf-8") for arg in call]

    # but keep the arg 'call' here because log_info inherently
    # attempts to encode the message... (TODO: should log_info
    # figure out the encoding of the input string)
    martian.log_info("Running crconverter: %s" % " ".join(call))
    try:
        results = tk_subproc.check_output(unicode_call, stderr=subprocess.STDOUT)
        martian.log_info("crconverter output: %s" % results)
    except subprocess.CalledProcessError as e:
        outs.output_for_cloupe = None
        martian.throw("Could not generate .cloupe file: \n%s" % e.output)
