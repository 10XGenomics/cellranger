#
# Copyright (c) 2024 10X Genomics, Inc. All rights reserved.
#

"""Utilities to run cell annotation on the cloud."""

import os
import pathlib
import time
from dataclasses import dataclass
from urllib.parse import urlparse

import martian
import requests

import cellranger.cell_typing.cloud_cas_client as cloud_client_lib
import cellranger.cell_typing.common as ct_common
from cellranger.cell_typing.cloud_cas_client import CloudSessionStats

# wait max 3 hours for completion
MAX_CLOUD_TIMEOUT_SECONDS = 60 * 60 * 3

# wait max 12 hours if annotation-only pipeline
MAX_CLOUD_TIMEOUT_SECONDS_ANNOTATION_ONLY = 60 * 60 * 12

# file names of cloud results
CLOUD_RESULTS_FILENAME = "cas_results.json.gz"
CLOUD_CELL_TYPES_FILENAME = "cas_cell_types.csv"
CLOUD_METADATA_FILENAME = "metadata.json"

# path to the token in the tarball
LOCAL_TOKEN_PATH = pathlib.Path(__file__).parent / "tokens" / "10xcloud_token.json"

# try 3 times max
MAX_RETRIES = 3
RETRY_INTERVAL = 180

# SAP-38/39
MIN_BARCODE_THRESHOLD = 100

MAX_BARCODE_THRESHOLD = 800_000

DEFAULT_HUMAN_MODEL = "human_pca_v1_beta"
DEFAULT_MOUSE_MODEL = "mouse-pca-512-log1p-zscore-v1"


@dataclass
class CasResults:
    """Stores results from CAS."""

    cell_annotation_results_json_gz: bytes | None = None
    cell_types: bytes | None = None
    cas_metadata: bytes | None = None

    def is_complete(self) -> bool:
        return all(
            (
                self.cell_annotation_results_json_gz is not None,
                self.cell_types is not None,
                self.cas_metadata is not None,
            )
        )

    def outputs_not_found(self) -> list:
        """Returns a list of outputs not found."""
        outputs_not_found = []
        if self.cell_annotation_results_json_gz is None:
            outputs_not_found.append(CLOUD_RESULTS_FILENAME)
        if self.cell_types is None:
            outputs_not_found.append(CLOUD_CELL_TYPES_FILENAME)
        if self.cas_metadata is None:
            outputs_not_found.append(CLOUD_METADATA_FILENAME)
        return outputs_not_found


@dataclass
class TokenInfo:
    """Stores token info."""

    successful: bool
    server_url: str | None = None
    cli_token: str | None = None
    cfa_tokens: dict | None = None


class CloudCellAnnotationException(Exception):
    """Non-fatal exceptions for cell annotation."""

    def __init__(self, session_stats: CloudSessionStats, *args):
        self.stats = session_stats
        super().__init__(*args)
        self.stats.mark_failure(str(self))


def extract_token_info(token_path_in) -> TokenInfo:
    """Extract CLI and (if applicable) CFA tokens from the supplied args.

    Returns (server_url, cli_token, cfa_token)
    """
    if not token_path_in:
        # look for a token in the tokens directory, use if present
        if not LOCAL_TOKEN_PATH.exists():
            raise ValueError("No token supplied")
        tenx_cloud_token_path = LOCAL_TOKEN_PATH.absolute()
    else:
        # presence of the token should also have been validated by a preflight
        if not os.path.exists(token_path_in):
            raise ValueError(f"no token exists at {token_path_in}")
        tenx_cloud_token_path = token_path_in

    # check for explicit use local token arg for additional header content
    use_local_token = os.environ.get("TENX_CLOUD_USE_LOCAL_TOKEN")
    args_cli_token: str | None = None
    if use_local_token is not None:
        # if env variable set, use CLI token from passed in path if
        # available; get other settings from embedded local token
        _, args_cli_token = ct_common.read_cli_token(tenx_cloud_token_path)
        tenx_cloud_token_path = LOCAL_TOKEN_PATH.absolute()

    server_url, token_info = ct_common.read_full_token(tenx_cloud_token_path)
    if server_url is not None and token_info is not None:
        # if the full token is valid JSON, parse it
        parse_result = urlparse(server_url)
        if not parse_result.scheme:
            server_url = f"https://{server_url}"
        cf_dict = token_info.get("cloudflare_token", None)
        if cf_dict:
            cfa_tokens = cf_dict
        else:
            cfa_tokens = None
        cli_token = token_info["cli_token"]
    else:
        # if the path contains a token only, return it
        server_url, cli_token = ct_common.read_cli_token(tenx_cloud_token_path)
        cfa_tokens = None

    # CLI token from passed in arg takes precedence over env token
    if args_cli_token is not None:
        cli_token = args_cli_token
    if server_url is None or cli_token is None:
        raise ValueError(f"Unable to read a valid token from {token_path_in}")
    return TokenInfo(
        server_url=server_url, cli_token=cli_token, cfa_tokens=cfa_tokens, successful=True
    )


def get_valid_models(models: list, ref: str) -> list:
    """Returns valid cell annotation models.

    Given a reference string and a list of class AnnotationModelDescription return the model names that can be used with the given reference.

    Args:
        models (list): List of class AnnotationModelDescription. Typically retrieved using `get_available_models`
        ref (string): Reference name as a string

    Returns:
        valid_models (list): List of valid models that can be used with the reference
    """
    valid_models = [
        model.name for model in models if any(genome in ref for genome in model.genomes)
    ]
    return valid_models


def test_cloud_connection(cloud_client: cloud_client_lib.CloudAPIClient) -> tuple[bool, bool]:
    """Test cloud connections.

    Returns (connected, valid_user).  Connected will be false if a connection could not be made;
    valid_user will be false if the user token is invalid.
    """
    cloud_client.session_stats.start(CloudSessionStats.OP_CONNECT)
    try:
        cloud_client.preflight_check()
        cloud_client.session_stats.passed(CloudSessionStats.OP_CONNECT)
        martian.log_info("cloud preflight check succeeded")
        return True, True
    except cloud_client_lib.CloudAPIError as exc:
        cloud_client.session_stats.failed(CloudSessionStats.OP_CONNECT, reason=str(exc))
        if exc.status_code == 403:
            martian.alarm(f"cloud preflight check had invalid user credentials: {exc}")
            return True, False
        else:
            martian.alarm(f"cloud preflight check failed: {exc}")
            return False, False


def run_analysis(
    filtered_matrix_path: str,
    cloud_client: cloud_client_lib.CloudAPIClient,
    pipestance_type: str,
    genome_name: str,
    pipeline_params: cloud_client_lib.AnnotationParams,
    max_timeout: float = MAX_CLOUD_TIMEOUT_SECONDS,
) -> tuple[bool, list[cloud_client_lib.AnnotationOutput] | None]:
    """Invoke the celltyper analysis and wait until a result.

    Args:
        filtered_matrix_path (str): Filtered matrix path.
        cloud_client (cloud_client_lib.CloudAPIClient): The client to use to communicate with the cloud.
        pipestance_type (str | None): The pipestance type.
        genome_name: The name of the genome to pass to the API (should agree with the H5)
        pipeline_params: The invocation parameters (model/pipeline version)
        max_timeout (float): the maximum timeout.

    Raises:
        AnnotationError: if any unexpected failure occurs.

    Returns:
        bool: If the upload and pipeline execution completed successfully.
        list[AnnotationOutputs]: The list of files to download from the completed analysis.
    """
    # upload H5 to cloud and kick off analysis
    # TODO: pass distance, dataset info vars in a "dev mode"
    martian.log_info("Uploading and invoking analysis")
    try:
        pending_analysis_id = cloud_client.annotate_file(
            filtered_matrix_path,
            pipestance_type=pipestance_type,
            genome_name=genome_name,
            pipeline_params=pipeline_params,
        )  # throws AnnotationError
    except cloud_client_lib.AnnotationError as error:
        cloud_client.session_stats.failed(CloudSessionStats.OP_UPLOAD, reason=str(error))
        return False, None

    # wait until max timeout for the analysis to complete
    # if it does not complete within the timeout time allotted, raise exception
    # and do not permit restart
    martian.log_info(f"Polling for analysis {pending_analysis_id}")
    start_time = time.time()

    analysis_complete = False
    while time.time() - start_time < max_timeout:
        annotation_completed = False
        try:
            annotation_completed = cloud_client.is_annotation_completed(pending_analysis_id)
        except cloud_client_lib.CloudAPIError:
            # API error might be transient, maintain previous sleep schedule
            time.sleep(10)

        if annotation_completed:
            analysis_complete = True
            time_taken_for_analysis = time.time() - start_time
            martian.log_info(f"Annotation done. Took {time_taken_for_analysis} sec.")
            try:
                martian.log_info("Checking if annotating was successful.")
                status = cloud_client.get_annotate_status(pending_analysis_id)
                if status.status == "COMPLETED":
                    # happy path
                    martian.log_info("Annotation successful.")
                    return True, status.outputs
                else:
                    # failed path -- can attempt to retry (could be transient)
                    martian.log_info("Annotation failed.")
                    return False, None
            except cloud_client_lib.AnnotationError:
                # this happens if a pipeline executes with no outs, or there
                # has been no runs associated with the analysis
                # assume pipeline will rerun in this manner when re-invoked;
                # do not attempt to repeat.
                return True, None
        else:
            # try every 10 sec
            time.sleep(10)

    # CELLRANGER-8729: if the operation times out, fail, do not permit restart
    if not analysis_complete:
        raise cloud_client_lib.AnnotationError("analysis timeout exceeded")


def download_analysis_results(
    cloud_client: cloud_client_lib.CloudAPIClient,
    outputs_to_download: list[cloud_client_lib.AnnotationOutput],
) -> CasResults:
    """Download the cell_annotation_results.json.gz file from the celltyper outs.

    Returns:
        None if the download was unsuccessful, a trio of byte arrays of the content if so.
    """
    download_url_pairs = [
        (output.name, cloud_client.get_output_download_url(output))
        for output in outputs_to_download
    ]
    if not outputs_to_download:
        martian.log_info("No files received to download")
    else:
        martian.log_info("Download requests for the following files received:")
        for output in outputs_to_download:
            martian.log_info(f"{output.name} with file-id {output.file_id}")

    cell_annotation_results = CasResults()
    for filename, download_url in download_url_pairs:
        # exclusively pull cell_annotation_results.json.gz
        if filename == CLOUD_RESULTS_FILENAME:
            martian.log_info(f"Downloading results file {filename} from url {download_url}")
            r = requests.get(download_url, allow_redirects=True, timeout=60)
            if r.status_code == 200:
                # the reason for the memory estimation in split
                cell_annotation_results.cell_annotation_results_json_gz = r.content

        if filename == CLOUD_CELL_TYPES_FILENAME:
            martian.log_info(f"Downloading cell types file {filename} from url {download_url}")
            r = requests.get(download_url, allow_redirects=True, timeout=60)
            if r.status_code == 200:
                cell_annotation_results.cell_types = r.content

        if filename == CLOUD_METADATA_FILENAME:
            martian.log_info(f"Downloading metadata file {filename} from url {download_url}")
            r = requests.get(download_url, allow_redirects=True, timeout=60)
            if r.status_code == 200:
                cell_annotation_results.cas_metadata = r.content

    if not cell_annotation_results.is_complete():
        outputs_not_found = cell_annotation_results.outputs_not_found()
        martian.alarm(f"No outputs with name {','.join(outputs_not_found)} found")
    return cell_annotation_results


def write_results(results_list_gz_bytes, results_json_gz_path) -> None:
    """Write results files."""
    # Write out the results with the barcode info, formatted as desired
    # (this may not be strictly necessary as celltyper outputs in this format)
    with open(results_json_gz_path, "wb") as f:
        f.write(results_list_gz_bytes)


def run_get_annotation_defaults(  # pylint: disable=too-many-arguments
    server_url: str,
    cli_token: str,
    cfa_tokens: dict | None,
    pipestance_type: str,
    genome_name: str,
    pipeline_requirements: list[str] | None,
    cell_annotation_model: str | None,
    max_retries: int = MAX_RETRIES,
    retry_interval: int = RETRY_INTERVAL,
) -> cloud_client_lib.AnnotationParams | None:
    """Get the appropriate defaults to invoke the cell annotation service."""
    with cloud_client_lib.CloudAPIClient(
        server_url,
        cli_token,
        cfa_tokens=cfa_tokens,
    ) as cloud_client:
        num_tries = 0

        # try connecting to server
        preflight_succeeded = False
        while num_tries < max_retries:
            num_tries += 1
            connected, valid_user = test_cloud_connection(cloud_client)
            if connected and valid_user:
                preflight_succeeded = True
                break
            elif connected and not valid_user:
                raise CloudCellAnnotationException(
                    cloud_client.session_stats, "Invalid cloud access token"
                )
            # if not connected, try again
            martian.log_info(f"Waiting {retry_interval} seconds to restart")
            time.sleep(retry_interval)

        if not preflight_succeeded:
            raise CloudCellAnnotationException(
                cloud_client.session_stats, "Could not connect to cloud server; skipping annotation"
            )

        num_tries = 0
        while num_tries < max_retries:
            num_tries += 1
            martian.log_info(f"Getting default model/pipeline params: attempt {num_tries}...")
            try:
                valid, params = cloud_client.get_annotation_params(
                    pipestance_type=pipestance_type,
                    genome_name=genome_name,
                    requires=pipeline_requirements,
                    cell_annotation_model=cell_annotation_model,
                )
                if valid:
                    return params
                else:
                    # no model-pipeline combination that meets requirements supplied;
                    # return None (calling stage should handle error behavior)
                    return None

            # if http 500, attempt to retry later
            except cloud_client_lib.AnnotationError:
                time.sleep(retry_interval)

    martian.log_info("Annotation parameter lookup failed.")
    raise CloudCellAnnotationException(
        cloud_client.session_stats,
        "Failed annotation parameter lookup after max retries",
    )


def run_cloud_query(  # pylint: disable=too-many-arguments
    server_url: str,
    cli_token: str,
    cfa_tokens: dict | None,
    filtered_matrix_path: str,
    pipestance_type: str,
    genome_name: str,
    pipeline_params: cloud_client_lib.AnnotationParams,
    max_retries=MAX_RETRIES,
    retry_interval=RETRY_INTERVAL,
    max_timeout=MAX_CLOUD_TIMEOUT_SECONDS,
) -> tuple[CasResults, CloudSessionStats]:
    """Run the cell annotation service and get results from the proxy server."""
    with cloud_client_lib.CloudAPIClient(
        server_url,
        cli_token,
        cfa_tokens=cfa_tokens,
    ) as cloud_client:
        # do a preflight check to fail fast (while allowing for transient cloud down-ness)
        num_tries = 0
        preflight_succeeded = False

        cloud_client.session_stats.cell_annotation_model = pipeline_params.model
        while num_tries < max_retries:
            num_tries += 1
            connected, valid_user = test_cloud_connection(cloud_client)
            if connected and valid_user:
                preflight_succeeded = True
                break
            elif connected and not valid_user:
                raise CloudCellAnnotationException(
                    cloud_client.session_stats, "Invalid cloud access token"
                )
            # if not connected, try again
            martian.log_info(f"Waiting {retry_interval} seconds to restart")
            time.sleep(retry_interval)

        if not preflight_succeeded:
            raise CloudCellAnnotationException(
                cloud_client.session_stats, "Could not connect to cloud server; skipping annotation"
            )

        # try upload/analysis up to MAX_RETRIES
        num_tries = 0
        annotation_succeeded = False
        while num_tries < max_retries:
            num_tries += 1
            try:
                annotation_succeeded, outputs_to_download = run_analysis(
                    filtered_matrix_path=filtered_matrix_path,
                    cloud_client=cloud_client,
                    pipestance_type=pipestance_type,
                    genome_name=genome_name,
                    pipeline_params=pipeline_params,
                    max_timeout=max_timeout,
                )
                if annotation_succeeded:
                    break
            except cloud_client_lib.AnnotationError as error:
                # bail on unexpected analysis condition
                cloud_client.session_stats.failed(CloudSessionStats.OP_ANALYZE, str(error))
                raise CloudCellAnnotationException(
                    cloud_client.session_stats, "Unexpected server-side analysis failure"
                ) from error
            # don't immediately retry -- give it 3 minutes to change transient condition
            martian.log_info("Waiting 3 minutes to retry analysis invocation")
            time.sleep(retry_interval)

        # if the pipeline did not succeed,
        if not annotation_succeeded:
            martian.log_info("Annotation failed.")
            raise CloudCellAnnotationException(
                cloud_client.session_stats,
                "Failed annotation after max retries",
            )

        # if the pipeline completed but there's no information, skip
        if not outputs_to_download:
            martian.log_info("The pipeline completed but there is no output to download.")
            raise CloudCellAnnotationException(
                cloud_client.session_stats,
                "Annotation run yielded no outputs",
            )

        num_tries = 0
        martian.log_info("Preparing to download results")
        while num_tries < max_retries:
            num_tries += 1
            cloud_client.session_stats.start(CloudSessionStats.OP_DOWNLOAD)
            try:
                martian.log_info(f"Attempting to download results the {num_tries}-th time.")
                results_content = download_analysis_results(cloud_client, outputs_to_download)
                if results_content.is_complete():
                    martian.log_info(f"Download of results successful on the {num_tries}-th try.")
                    cloud_client.session_stats.passed(CloudSessionStats.OP_DOWNLOAD)
                    break
            except requests.HTTPError as exc:
                martian.log_info(f"Download of results failed on the {num_tries}-th try.")
                cloud_client.session_stats.failed(CloudSessionStats.OP_DOWNLOAD, reason=str(exc))

            # don't immediately retry -- give it 3 minute to change download condition
            # (and allow for a cloud restart)
            martian.log_info("Waiting 3 minutes to retry failed download")
            time.sleep(retry_interval)

        if not results_content.is_complete():
            cloud_client.session_stats.failed(CloudSessionStats.OP_DOWNLOAD)
            raise CloudCellAnnotationException(
                cloud_client.session_stats, "Could not download results from cloud endpoint"
            )

        cloud_client.session_stats.mark_successful()
        return results_content, cloud_client.session_stats
