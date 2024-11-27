#!/usr/bin/env python
#
# Copyright (c) 2024 10X Genomics, Inc. All rights reserved.
#
"""Determine whether or not the cell annotation stage should run."""


__MRO__ = """
stage CELL_ANNOTATION_PREFLIGHT(
    in  bool   is_cr_annotate,
    in  bool   is_multi,
    in  string tenx_cloud_token_path,
    in  string cell_annotation_model,
    in  path   reference_path,
    in  string cas_track_name,
    in  bool   skip_cell_annotation,
    in  h5     cr_annotate_filtered_matrix,
    in  cloupe cr_annotate_sample_cloupe,
    src py     "stages/cas_cell_typing/cell_annotation_preflight",
) using (
    mem_gb   = 8,
    volatile = strict,
)
"""

import os
import re
import time

import martian

import cellranger.cell_typing.cloud_cas_client as cloud_client_lib
import cellranger.cell_typing.cloud_cas_utils as cloud_cas_utils
import cellranger.matrix as cr_matrix


def _check_connection(args, cloud_client: cloud_client_lib.CloudAPIClient):
    """Verify that we can connect to the cloud server with the current client.

    Exit if there is something wrong with the setup.  This will try 3 times in case of the
    server being unreachable, but will exit immediately if there is something wrong with the token.
    """
    num_tries = 0
    max_retries = 3
    retry_interval = 60
    preflight_succeeded = False
    while num_tries < max_retries:
        num_tries += 1
        try:
            connected, valid_user = cloud_cas_utils.test_cloud_connection(cloud_client)
        except ValueError as err:
            martian.exit(
                "Error while trying to test 10x Cloud connection. "
                f"Please check the contents of {args.tenx_cloud_token_path}, or run 'cellranger cloud auth verify' to verify that your token is valid. \n"
                f"Error seen: {err}"
            )
        if connected and valid_user:
            preflight_succeeded = True
            break
        elif connected and not valid_user:
            martian.exit(
                f"Invalid cloud access token found at {args.tenx_cloud_token_path}.  "
                f"Please check the contents of this file, or run 'cellranger cloud auth setup' to re-authenticate."
            )
        # wait momentarily for server to come back
        martian.log_info(f"Waiting {retry_interval} seconds for server response")
        time.sleep(retry_interval)

    if not preflight_succeeded:
        martian.exit(f"Could not connect to cloud server after {max_retries} attempts")


def main(args, _outs):
    if args.skip_cell_annotation:
        return

    # CHECKS COMMON TO ANNOTATE, COUNT, MULTI
    # if cloud token path is None, check local token
    if args.tenx_cloud_token_path is None:
        if not cloud_cas_utils.LOCAL_TOKEN_PATH.exists():
            martian.log_warn(
                "No 10x cloud token defined. To enable cell annotation via the 10x Cloud, please run `cellranger cloud auth setup`."
            )
    elif args.tenx_cloud_token_path and not os.path.exists(args.tenx_cloud_token_path):
        martian.exit(
            f"No 10x cloud token found at {args.tenx_cloud_token_path}. To create a token, please run `cellranger cloud auth setup`. "
            "Please check that the path and is valid and ensure the permissions are correct."
        )
    # check if the token is valid
    token_info = None
    try:
        token_info = cloud_cas_utils.extract_token_info(args.tenx_cloud_token_path)
    except Exception as token_exc:  # pylint: disable=broad-except
        if args.tenx_cloud_token_path is not None:
            martian.exit(
                f"Could not read token information from {args.tenx_cloud_token_path}: {token_exc}. "
                "Please check that the path and is valid and ensure the permissions are correct."
            )
        elif args.is_multi:
            martian.exit(
                "Failed to authenticate for Cell Ranger Cloud. "
                "Please run `cellranger cloud auth setup` to set up cloud authentication or provide a valid token via the tenx-cloud-token-path parameter in the config CSV."
            )
        else:
            martian.exit(
                "Failed to authenticate for Cell Ranger Cloud. "
                "Please run `cellranger cloud auth setup` to set up cloud authentication or provide a valid token via the --tenx-cloud-token-path parameter."
            )

    # check if the model provided is valid (and that we can connect to the cloud)
    if token_info:
        with cloud_client_lib.CloudAPIClient(
            token_info.server_url,
            token_info.cli_token,
            cfa_tokens=token_info.cfa_tokens,
        ) as cloud_client:
            _check_connection(args, cloud_client)

            if args.cell_annotation_model:
                # this could also be in a retry loop but it's unlikely to fail
                models = cloud_client.get_available_models()
                matched_models = [
                    model for model in models if model.name == f"{args.cell_annotation_model}"
                ]
                if len(matched_models) != 1:
                    martian.exit(
                        f"{args.cell_annotation_model} is not a valid cell annotation model. "
                        "Please run `cellranger cloud annotation models` to view available models."
                    )

                genomes = set(matched_models[0].genomes)
                if args.reference_path is not None:
                    ref = [os.path.basename(args.reference_path)]

                    # Check if annotation ref and genome are compatible
                    # Allows for matching GRCh38 to custom genomes like GRCh38-2024-A_and_XYZ.
                    if not any(any(genome in r for r in ref) for genome in genomes):
                        valid_models = cloud_cas_utils.get_valid_models(models, ref)
                        martian.exit(
                            f"The {args.cell_annotation_model} model cannot run data from the {ref} genome. "
                            f"Valid models for this reference are {valid_models}"
                        )

    # CHECKS SPECIFIC TO ANNOTATE
    if args.is_cr_annotate:
        # Check if cas_track_name_from_user contains invalid characters
        # that can't be displayed in loupe browser
        if args.cas_track_name:
            re_check = re.compile(r"^[a-zA-Z0-9\s_-]*$")
            if not re_check.match(args.cas_track_name):
                martian.exit(
                    f"Invalid characters in cloupe group name: {args.cas_track_name}. "
                    "Only alphanumeric characters, spaces, underscores, and dashes are allowed."
                )

        # Barnyard check
        # Checks if the matrix file exists. Guaranteed in count or multi but not annotate
        if not os.path.exists(args.cr_annotate_filtered_matrix):
            martian.exit(
                f"Filtered matrix file not found at {args.cr_annotate_filtered_matrix}. "
                "Please check that the path and is valid and ensure the permissions are correct."
            )
        try:
            ref = cr_matrix.CountMatrix.get_genomes_from_h5(args.cr_annotate_filtered_matrix)
            _, num_bcs, _ = cr_matrix.CountMatrix.load_dims_from_h5(
                args.cr_annotate_filtered_matrix
            )
        except Exception as matrix_exc:  # pylint: disable=broad-except
            martian.exit(
                f"Could not validate the matrix.h5 {args.cr_annotate_filtered_matrix}: {matrix_exc}. "
                "Please check that the path and is valid and ensure the permissions are correct "
                "and the the matrix file is a valid 10x Genomics matrix.h5 file."
            )
        if len(ref) > 1:
            martian.exit(
                f"Multi-genome references, {ref}, are not supported by cell annotation at this time."
            )
        if num_bcs > 800000:
            martian.exit(
                f"Filtered matrix file {args.cr_annotate_filtered_matrix} contains more barcodes than the maximum allowed threshold. "
                "Cell Ranger Cell Annotation has a maximum limit of 800,000 barcodes."
            )

        # check if the cloupe file is valid
        ## Check the path is valid
        if args.cr_annotate_sample_cloupe:
            if not os.path.exists(args.cr_annotate_sample_cloupe):
                martian.exit(
                    f"Cloupe file not found at {args.cr_annotate_sample_cloupe}. "
                    "Please check that the path and is valid and ensure the permissions are correct."
                )
        ## Check if the extension is valid
        if args.cr_annotate_sample_cloupe and not args.cr_annotate_sample_cloupe.endswith(
            ".cloupe"
        ):
            martian.exit(f"{args.cr_annotate_sample_cloupe} does not end with a .cloupe extension.")
