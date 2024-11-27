#!/usr/bin/env python
#
# Copyright (c) 2023 10X Genomics, Inc. All rights reserved.
#

"""Call cell types based on Broad Cell Annotation Service."""

__MRO__ = """
stage CALL_CLOUD_CELL_TYPES(
    in  string          sample_id,
    in  string          sample_desc,
    in  h5              filtered_matrix,
    in  string          cell_annotation_model,
    in  file            tenx_cloud_token_path,
    in  string          pipestance_type,
    in  bool            override_num_bc_limit,
    in  cloupe          sample_cloupe,
    out CellTypeResults cell_type_results,
    out bool            cas_success,
    out bool            disable_cas_ws,
    out bool            disable_summarize,
    out json            summary,
    src py              "stages/cas_cell_typing/call_cloud_cell_types",
) split (
) using (
)
"""

import csv
import math
from dataclasses import dataclass

import martian

import cellranger.cell_typing.cloud_cas_client as cloud_client_lib
import cellranger.cell_typing.cloud_cas_utils as cloud_utils
import cellranger.cell_typing.common as ct_common
import cellranger.matrix as cr_matrix
from tenkit import safe_json


@dataclass
class RunInfo:
    """Class storing run info."""

    num_features: int
    num_bcs: int
    num_entries: int
    is_cra: bool


def write_analysis_files(
    filtered_matrix: str | bytes,
    cell_annotation_results: cloud_utils.CasResults,
    model_used: str | None,
) -> ct_common.CellTypeResults:
    """Return cell_type_results with necessary outs and write CSV/JSON outputs."""
    results_json_gz_path = martian.make_path("results.json.gz").decode()
    cell_types_csv_path = martian.make_path("cell_types.csv").decode()
    cell_typing_metadata_path = martian.make_path("metadata.json").decode()

    cell_type_results = ct_common.CellTypeResults(model_used=model_used)
    if cell_annotation_results.cell_annotation_results_json_gz:
        cloud_utils.write_results(
            cell_annotation_results.cell_annotation_results_json_gz,
            results_json_gz_path,
        )
        cell_type_results.results = results_json_gz_path
    else:

        cell_type_results.results = None

    if cell_annotation_results.cell_types:
        with open(cell_types_csv_path, "w") as f:
            f.write(cell_annotation_results.cell_types.decode())
        cell_type_results.cell_types = cell_types_csv_path

        # compute bc fraction by comparing the return length vs input h5 bc length
        _, num_bcs, _ = cr_matrix.CountMatrix.load_dims_from_h5(filtered_matrix)
        with open(cell_types_csv_path) as f:
            reader = csv.DictReader(f)
            cell_type_results.frac_returned_bcs = round(sum(1 for _ in reader) / num_bcs, 3)

    else:
        cell_type_results.cell_types = None
        cell_type_results.frac_returned_bcs = 0.0

    if cell_annotation_results.cas_metadata:
        with open(cell_typing_metadata_path, "w") as f:
            f.write(cell_annotation_results.cas_metadata.decode())
        cell_type_results.metadata = cell_typing_metadata_path
    else:
        cell_type_results.metadata = None

    return cell_type_results


def alarm_vs_exit(message, use_exit):
    """Use martian alarm vs martian exit given some condition.

    Parameters:
    - message (str): The message to be logged with the alarm or exit.
    - use_exit (bool): If True, the function will exit. If False, it will trigger an alarm.

    Returns:
    None
    """
    if use_exit:
        martian.exit(message)
    else:
        martian.alarm(message)


def get_token_info(
    args,
    outs,
    run_info,
    error_cell_type_results,
):
    """Get token info."""
    try:
        token_info = cloud_utils.extract_token_info(args.tenx_cloud_token_path)
    except Exception as token_exc:  # pylint: disable=broad-except
        martian.clear(outs)
        outs.cell_type_results = error_cell_type_results
        alarm_vs_exit(
            message=f"Could not read token information from {args.tenx_cloud_token_path}: {token_exc}; skipping annotation",
            use_exit=run_info.is_cra,
        )
        outs.cas_success = False
        token_info = cloud_utils.TokenInfo(successful=False)
    return token_info


def get_annotation_params(
    args, outs, model_name, genome_name, token_info, run_info
) -> cloud_client_lib.AnnotationParams | None:
    """Attempt to get the parameters required to invoke against the cloud.

    This will throw an alarm/exception
    """
    error_cell_type_results = ct_common.CellTypeResults(
        model_used=model_name,
        skip_downstream=True,
    ).to_outs()

    try:
        martian.log_info(
            f"Trying to get annotation params for pipestance {args.pipestance_type}, genome {genome_name}, model {model_name}"
        )
        annotation_params = cloud_utils.run_get_annotation_defaults(
            server_url=token_info.server_url,
            cli_token=token_info.cli_token,
            cfa_tokens=token_info.cfa_tokens,
            pipestance_type=args.pipestance_type,
            genome_name=genome_name,
            pipeline_requirements=[ct_common.CELL_ANNOTATION_REQUIRES_JSON_RESULTS],
            cell_annotation_model=model_name,
        )
        martian.log_info(f"resolved annotation params: {annotation_params}")

        # no valid cloud pipeline available
        if annotation_params is None:
            outs.cas_success = False
            outs.disable_cas_ws = True
            outs.cell_type_results = error_cell_type_results
            outs.disable_summarize = True
            alarm_vs_exit(
                message=f"no model-pipeline combination found for pipeline inputs: pipestance type {args.pipestance_type}, genome {genome_name}, model {model_name}",
                use_exit=run_info.is_cra,
            )

        return annotation_params
    except cloud_utils.CloudCellAnnotationException as cloud_cas_exception:
        session_stats = cloud_cas_exception.stats
        with open(outs.summary, "w") as outfile:
            outfile.write(safe_json.safe_jsonify(session_stats.stats, pretty=True))

        outs.cas_success = False
        outs.disable_cas_ws = False
        outs.cell_type_results = error_cell_type_results
        outs.disable_summarize = True

        alarm_vs_exit(
            message=f"cell annotation param lookup failed on barcodes in sample ID: {args.sample_id}, sample description: {args.sample_desc}, cell annotation model: {model_name}. "
            f"\nCell annotation lookup failed: {cloud_cas_exception}",
            use_exit=run_info.is_cra,
        )
        return None


def run_cell_annotation(
    args,
    outs,
    pipeline_params: cloud_client_lib.AnnotationParams,
    genome_name,
    token_info,
    run_info,
):
    """Ping the end-point and obtain cell annotations."""
    error_cell_type_results = ct_common.CellTypeResults(
        model_used=pipeline_params.model,
        skip_downstream=True,
    ).to_outs()
    try:
        if run_info.is_cra:
            max_timeout = cloud_utils.MAX_CLOUD_TIMEOUT_SECONDS_ANNOTATION_ONLY
        else:
            max_timeout = cloud_utils.MAX_CLOUD_TIMEOUT_SECONDS
        cell_annotation_results, session_stats = cloud_utils.run_cloud_query(
            server_url=token_info.server_url,
            cli_token=token_info.cli_token,
            cfa_tokens=token_info.cfa_tokens,
            pipestance_type=args.pipestance_type,
            genome_name=genome_name,
            pipeline_params=pipeline_params,
            filtered_matrix_path=args.filtered_matrix,
            max_timeout=max_timeout,
        )

    except cloud_utils.CloudCellAnnotationException as cloud_cas_exception:
        session_stats = cloud_cas_exception.stats
        with open(outs.summary, "w") as outfile:
            outfile.write(safe_json.safe_jsonify(session_stats.stats, pretty=True))

        outs.cas_success = False
        outs.disable_cas_ws = False
        outs.cell_type_results = error_cell_type_results
        outs.disable_summarize = True
        alarm_vs_exit(
            message=f"cell annotation failed on barcodes in sample ID: {args.sample_id}, sample description: {args.sample_desc}, cell annotation model: {pipeline_params.model}. "
            f"num features: {run_info.num_features}, num barcodes: {run_info.num_bcs}, num non zero entries: {run_info.num_entries}"
            f"\nCell annotation failed: {cloud_cas_exception}",
            use_exit=run_info.is_cra,
        )
        return
    except ValueError as api_exc:
        martian.clear(outs)
        outs.cell_type_results = error_cell_type_results
        alarm_vs_exit(
            message=f"cell annotation failed on barcodes in sample ID: {args.sample_id}, sample description: {args.sample_desc}, cell annotation model: {pipeline_params.model}. "
            f"num features: {run_info.num_features}, num barcodes: {run_info.num_bcs}, num non zero entries: {run_info.num_entries}"
            f"\nInvalid parameters passed to cloud API client: {api_exc}",
            use_exit=run_info.is_cra,
        )
        outs.cas_success = False
        outs.disable_summarize = True
        return
    except Exception as fail_message:  # pylint: disable=broad-except
        martian.clear(outs)
        outs.cell_type_results = error_cell_type_results
        alarm_vs_exit(
            message=f"cell annotation failed on barcodes in sample ID: {args.sample_id}, sample description: {args.sample_desc}, cell annotation model: {pipeline_params.model}. "
            f"num features: {run_info.num_features}, num barcodes: {run_info.num_bcs}, num non zero entries: {run_info.num_entries}"
            f"\nUncaught message in cloud annotation: {fail_message}",
            use_exit=run_info.is_cra,
        )
        outs.cas_success = False
        outs.disable_summarize = True
        return
    outs.disable_summarize = False

    with open(outs.summary, "w") as outfile:
        outfile.write(safe_json.safe_jsonify(session_stats.stats, pretty=True))

    if not cell_annotation_results.is_complete():
        martian.clear(outs)
        outs.cell_type_results = error_cell_type_results
        alarm_vs_exit(
            message=f"cell annotation failed on barcodes in sample ID: {args.sample_id}, sample description: {args.sample_desc}, cell annotation model: {pipeline_params.model}. "
            f"num features: {run_info.num_features}, num barcodes: {run_info.num_bcs}, num non zero entries: {run_info.num_entries}"
            "\ncell annotation returned empty results.",
            use_exit=run_info.is_cra,
        )
        outs.cas_success = False
        outs.disable_summarize = True
        return

    cell_type_results = write_analysis_files(
        args.filtered_matrix,
        cell_annotation_results,
        model_used=pipeline_params.model,
    )
    outs.cell_type_results = cell_type_results.to_outs()
    outs.cas_success = True
    outs.disable_summarize = False


def split(args):
    if not args.filtered_matrix:
        return {
            "chunks:": [],
            "join": {"__mem_gb": 1},
        }

    _, num_bcs, _ = cr_matrix.CountMatrix.load_dims_from_h5(args.filtered_matrix)
    # rough estimate -- JSON representation of output will be 20k per barcode (it's verbose)
    join_mem_gb = int(math.ceil((20_000 * num_bcs) / (1024 * 1024 * 1024)))

    return {
        "chunks": [],
        "join": {
            "__mem_gb": join_mem_gb,
        },
    }


def join(
    args, outs, _chunk_defs, _chunk_outs
):  # pylint: disable=too-many-locals,too-many-statements
    is_cra = args.pipestance_type in ["CELLRANGER_ANNOTATE_PD", "CELLRANGER_ANNOTATE_CS"]

    if not args.pipestance_type:
        args.pipestance_type = "UNKNOWN_PIPESTANCE"

    # For backward compatibility run all cell annotation runs through cloud
    if not args.cell_annotation_model or len(args.cell_annotation_model.split(":", 1)) < 2:
        inferred_model_name = args.cell_annotation_model
    else:
        model_tokens = args.cell_annotation_model.split(":", 1)
        inferred_model_name = model_tokens[1] if model_tokens[1] else None

    # set model name to explicit None if blank / unset
    if not inferred_model_name:
        inferred_model_name = None

    error_cell_type_results = ct_common.CellTypeResults(
        model_used=inferred_model_name,
        skip_downstream=True,
    ).to_outs()

    # SAP-38; after SAP-39 resolution, drop THRESHOLD to 1
    num_features, num_bcs, num_entries = cr_matrix.CountMatrix.load_dims_from_h5(
        args.filtered_matrix
    )
    run_info = RunInfo(
        num_features=num_features, num_bcs=num_bcs, num_entries=num_entries, is_cra=is_cra
    )

    if run_info.num_bcs < cloud_utils.MIN_BARCODE_THRESHOLD:
        martian.clear(outs)
        outs.cell_type_results = error_cell_type_results
        alarm_vs_exit(
            message=f"Under minimum barcode count required to compute cell annotation: {run_info.num_bcs}",
            use_exit=run_info.is_cra,
        )
        outs.cas_success = False
        outs.disable_summarize = True
        return

    if not args.override_num_bc_limit and run_info.num_bcs > cloud_utils.MAX_BARCODE_THRESHOLD:
        martian.clear(outs)
        outs.cell_type_results = error_cell_type_results
        alarm_vs_exit(
            message=f"Over the maximum barcode count to compute cell annotation: {run_info.num_bcs}",
            use_exit=run_info.is_cra,
        )
        outs.cas_success = False
        outs.disable_summarize = True
        return

    if not args.filtered_matrix:
        martian.clear(outs)
        outs.cell_type_results = error_cell_type_results
        alarm_vs_exit(
            message="No filtered matrix supplied; skipping annotation", use_exit=run_info.is_cra
        )
        outs.cas_success = False
        outs.disable_cas_ws = True
        outs.disable_summarize = True
        return

    # skip annotation barnyard sample
    ref = cr_matrix.CountMatrix.get_genomes_from_h5(args.filtered_matrix)
    if len(ref) > 1:
        martian.log_warn(
            f"Multi-genome references, {ref}, are not supported by cell annotation at this time."
        )
        outs.disable_summarize = True
        return

    # Make sure the genome is some version of GRCh38 or mm10 but there aren't more than 1 and that there is a GEX library
    valid_genome_library, _, genome_name = ct_common.check_valid_genome_and_library(
        matrix=args.filtered_matrix
    )
    if not valid_genome_library:
        martian.clear(outs)
        outs.cell_type_results = error_cell_type_results
        alarm_vs_exit(
            message="Cell annotation was skipped because the reference used is not yet supported",
            use_exit=run_info.is_cra,
        )
        outs.cas_success = False
        outs.disable_cas_ws = True
        outs.disable_summarize = True
        return

    # lookup cloud access token information
    token_info = get_token_info(
        args=args,
        outs=outs,
        run_info=run_info,
        error_cell_type_results=error_cell_type_results,
    )
    if not token_info.successful:
        outs.disable_summarize = True
        return

    # look up annotation parameters (pipeline version, model name) from cloud first
    # (this sets outs if it fails and returns None-- not sure if I like the side effects)
    params = get_annotation_params(
        args=args,
        outs=outs,
        model_name=inferred_model_name,
        genome_name=genome_name,
        token_info=token_info,
        run_info=run_info,
    )

    # only run cell annotation if valid params returned
    if params is not None:
        martian.log_info(
            f"Trying to run cell annotation on barcodes. Sample ID: {args.sample_id}, Sample description: {args.sample_desc}, cell annotation model: {params.model}, "
            f"num features: {run_info.num_features}, num barcodes: {run_info.num_bcs}, num non zero entries: {run_info.num_entries}"
        )
        run_cell_annotation(
            args=args,
            outs=outs,
            pipeline_params=params,
            genome_name=genome_name,
            token_info=token_info,
            run_info=run_info,
        )
