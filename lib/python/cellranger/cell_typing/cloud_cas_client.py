#!/usr/bin/python3
#
# Copyright (c) 2024 10X Genomics, Inc. All rights reserved.
#

"""Python Client to run cloud cell annotation."""

import base64
import functools
import hashlib
import logging
import os
import random
import string
import textwrap
import time
from collections import defaultdict
from dataclasses import dataclass
from urllib.parse import urlencode, urljoin

import requests
from requests.adapters import HTTPAdapter, Retry

import cellranger.matrix as cr_matrix

TEST_ENV = "test"
PROD_ENV = "prod"

cloud_analysis_endpoints = {
    "annotation-api": "api/cloud-analysis/annotation-api/v1",
    "files-api": "api/cloud-analysis/files-api/v1",
    "pipelines-api": "api/cloud-analysis/pipelines-api/v1",
    "projects-api": "api/cloud-analysis/projects-api/v1",
}


@dataclass
class AnnotationOutput:
    """Annotation output.

    file_id: file id
    name: file name.
    """

    file_id: str
    name: str


@dataclass
class AnnotationStatus:
    """Annotation status.

    id: the annotation request id
    status: the annotation status
    outputs: the annotation outputs on successful completion.
    """

    annotation_id: str
    status: str
    outputs: list[AnnotationOutput]


class AnnotationError(Exception):
    """Annotation Error."""


class CloudAPIError(Exception):
    """API Error."""

    def __init__(self, status_code: int, *args: object) -> None:
        """Create an exception and pass the HTTP error code as a param."""
        super().__init__(*args)
        self.status_code = status_code


@dataclass
class AnnotationModelDescription:
    """Annotation model description."""

    name: str
    displayName: str  # pylint: disable=invalid-name
    description: str
    date: str
    genomes: list[str]


@dataclass
class AnnotationParams:
    """Parameters for invoking the annotation analysis."""

    pipelineName: str  # pylint: disable=invalid-name
    pipelineVersion: str  # pylint: disable=invalid-name
    model: str


def _format_headers(headers):
    return "\n".join(f"{k}: {v}" for k, v in headers.items())


def print_roundtrip(response, *args, **kwargs):
    """Print request and response."""
    reqbody = response.request.body
    if reqbody and len(reqbody) > 100:
        reqbody = "Request body omitted"

    print(f" kwargs: {kwargs}")
    print(
        textwrap.dedent(
            """
        ---------------- request ----------------
        {req.method} {req.url}
        {reqhdrs}

        {req_body}
        ---------------- response ----------------
        {res.status_code} {res.reason} {res.url}
        {reshdrs}

        {res.text}
        """
        ).format(
            req=response.request,
            req_body=reqbody,
            res=response,
            reqhdrs=_format_headers(response.request.headers),
            reshdrs=_format_headers(response.headers),
        )
    )


def md5sum_file(file_path: str, chunk_size: int = 8192) -> str:
    """Get file md5 checksum."""
    md5 = hashlib.md5()
    with open(file_path, "rb") as f:
        while chunk := f.read(chunk_size):
            md5.update(chunk)

    return md5.hexdigest()


def check_api_errors(r):
    """Check for http response error."""
    http_error_msg = ""
    reason = r.reason

    if 400 <= r.status_code < 500:
        try:
            reason = r.json()
        except requests.exceptions.JSONDecodeError:
            # use default reason if failed to parse error response
            pass

        http_error_msg = f"Client Error {r.status_code}: {reason} for url: {r.url}"
    elif 500 <= r.status_code < 600:
        http_error_msg = f"Server Error {r.status_code}: {reason} for url: {r.url}"

    if http_error_msg:
        raise CloudAPIError(r.status_code, http_error_msg)


def http_json_response(func):
    """Wrapper to return json response."""

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        try:
            resp = func(*args, **kwargs)
            check_api_errors(resp)
            return resp.json()
        except requests.exceptions.JSONDecodeError as error:
            # could use something like HTTP 417 here but for now zero out the error code
            raise CloudAPIError(0, f"Error parsing {func.__name__} response json") from error

    return wrapper


class CloudSessionStats:
    """Cloud session statistics."""

    OP_CONNECT = "connect"
    OP_UPLOAD = "upload"
    OP_ANALYZE = "analyze"
    OP_DOWNLOAD = "download"

    operations = (OP_CONNECT, OP_UPLOAD, OP_ANALYZE, OP_DOWNLOAD)

    def __init__(self):
        """Initialize operation counters and other base statistics."""
        self._counts = defaultdict(lambda: 0)
        self._successes = defaultdict(lambda: False)
        self._starts = defaultdict(time.time)
        self._times = defaultdict(lambda: 0.0)
        self.file_size: int = 0
        self.num_barcodes: int = 0
        self.cell_annotation_model: str = ""
        self.server_url: str = ""
        self.analysis_id: str | None = None
        self._complete_success = False
        self._start_time = time.time()
        self._total_time = 0
        self._failure_mode = ""

    def __incr(self, operation):
        """Indicate that an operation was tried."""
        self._counts[operation] += 1

    def __succeed(self, operation):
        """Indicate that an operation succeeded."""
        self._successes[operation] = True

    def __fail(self, operation):
        """Indicate that an operation failed."""
        self._successes[operation] = False

    def __start(self, operation):
        """Indicate that an operation started."""
        self._starts[operation] = time.time()

    def __end(self, operation):
        """Indicate that an operation ended."""
        self._times[operation] = time.time() - self._starts[operation]

    def start(self, operation):
        """Start and count an operation."""
        self.__incr(operation)
        self.__start(operation)

    def passed(self, operation):
        """End and mark an operation as successful."""
        self.__end(operation)
        self.__succeed(operation)

    def failed(self, operation, reason: str | None = None):
        """End and mark an operation as failed.  Failure mode will be set."""
        self.__end(operation)
        self.__fail(operation)
        if reason is None:
            self._failure_mode = f"Failed {operation}"
        else:
            self._failure_mode = reason

    def __stat(self, operation):
        return {
            "tries": self._counts[operation],
            "success": self._successes[operation],
            "duration": self._times[operation],
        }

    def mark_successful(self):
        self._complete_success = True
        self._total_time = time.time() - self._start_time

    def mark_failure(self, reason: str):
        """Set the recorded failure mode to the specified explicit reason.

        This may not pertain to a particular operation.  This will also
        establish the duration of the stats lifecycle.
        """
        self._complete_success = False
        self._failure_mode = reason
        self._total_time = time.time() - self._start_time

    @property
    def stats(self) -> dict:
        """Return statistics about cloud operations as a dictionary."""
        base_stats = dict()
        # need to flatten the operation-level statistics for metrics consumption
        for operation in self.__class__.operations:
            op_stat = self.__stat(operation)
            base_stats[f"{operation}_tries"] = op_stat["tries"]
            base_stats[f"{operation}_duration"] = op_stat["duration"]
            base_stats[f"{operation}_success"] = op_stat["success"]
        base_stats["barcodes"] = self.num_barcodes
        base_stats["file_size"] = self.file_size
        base_stats["cell_annotation_model"] = self.cell_annotation_model
        base_stats["server_url"] = self.server_url
        base_stats["success"] = self._complete_success
        base_stats["failure_mode"] = self._failure_mode
        base_stats["time_elapsed"] = self._total_time
        base_stats["analysis_id"] = self.analysis_id
        return base_stats


class CloudAPIClient:
    """API Client to access 10X Cloud."""

    def __init__(
        self,
        server_url: str,
        api_token: str,
        cfa_tokens: dict | None = None,
        headers=None,
        verbosity: int = 0,
    ):
        self._base_url = server_url
        self.headers = headers
        self.verbosity = verbosity
        self.endpoints = cloud_analysis_endpoints
        self.session_stats = None

        self.urls = {k: urljoin(self._base_url, v) for (k, v) in cloud_analysis_endpoints.items()}

        if not api_token:
            raise ValueError("api_token not provided")

        self.cookies = {"session_id": api_token}

        if cfa_tokens:
            self.headers = cfa_tokens

        self._session = requests.Session()
        self._session.headers = self.headers
        self._session.cookies.update(self.cookies)

        retries = Retry(
            total=5,
            backoff_factor=0.1,
            status_forcelist=[502, 504, 524],
            allowed_methods=["GET", "POST", "PUT"],
        )
        self._session.mount("http://", HTTPAdapter(max_retries=retries))
        self._session.mount("https://", HTTPAdapter(max_retries=retries))
        if verbosity > 2:
            self._session.hooks["response"].append(print_roundtrip)

    def __enter__(self):
        self.session_stats = CloudSessionStats()
        self.session_stats.server_url = self._base_url
        return self

    def __exit__(self, *args):
        self.session_stats = None
        self._session.close()

    @property
    def stats(self) -> CloudSessionStats | None:
        return self.session_stats

    @http_json_response
    def _get_current_user(self) -> requests.Response:
        users_url = urljoin(self._base_url, "api/core/users-api/v1/user")
        return self._session.get(users_url)

    @http_json_response
    def _create_file_upload(
        self, input_file: str, file_name: str, file_size: int, md5: str, metadata=None
    ) -> requests.Response:
        """Create a file upload."""
        body = {
            "filename": file_name,
            "md5": md5,
            "size": file_size,
            "path": input_file,
            "fileType": "INPUT",
            "metadata": metadata,
        }
        files_url = self.urls["files-api"]
        return self._session.post(f"{files_url}/uploads", json=body)

    @http_json_response
    def _start_file_multiupload(self, blob_id: str) -> requests.Response:
        """Start File Multipart Upload API."""
        params = urlencode({"blobId": blob_id})
        files_url = self.urls["files-api"]
        return self._session.get(f"{files_url}/multipart/start-upload?{params}")

    @http_json_response
    def _get_part_upload_url(
        self,
        blob_id: str,
        upload_id: str,
        part_num: int,
        md5base64: str,
        expire_seconds: int,
    ) -> requests.Response:
        params = {
            "blobId": blob_id,
            "partNumber": part_num,
            "hash": md5base64,
            "expirationSeconds": expire_seconds,
            "uploadId": upload_id,
        }
        encode_params = urlencode(params)
        files_url = self.urls["files-api"]
        return self._session.get(
            f"{files_url}/multipart/get-upload-url?{encode_params}",
        )

    @http_json_response
    def _complete_file_multiupload(
        self, upload_id: str, blob_id: str, parts: dict
    ) -> requests.Response:
        files_url = self.urls["files-api"]
        return self._session.post(
            f"{files_url}/multipart/complete-upload",
            json={"parts": parts, "blobId": blob_id, "uploadId": upload_id},
        )

    @http_json_response
    def _complete_file_upload(self, file_id: str) -> requests.Response:
        files_url = self.urls["files-api"]
        return self._session.put(
            f"{files_url}/uploads/{file_id}",
            json={"status": "DONE"},
        )

    @http_json_response
    def _create_adhoc_analysis(self, req: dict) -> requests.Response:
        pipelines_url = self.urls["pipelines-api"]
        return self._session.post(
            f"{pipelines_url}/adhoc-analyses",
            json=req,
        )

    @http_json_response
    def _get_analysis(self, analysis_id: str) -> requests.Response:
        projects_url = self.urls["projects-api"]
        return self._session.get(
            f"{projects_url}/analyses/{analysis_id}",
        )

    @http_json_response
    def _get_analysis_output(self, project_id: str, analysis_id: str) -> requests.Response:
        projects_url = self.urls["files-api"]
        return self._session.get(
            f"{projects_url}/projects/{project_id}/analyses/{analysis_id}",
        )

    @http_json_response
    def _create_file_download(self, file_id: str) -> requests.Response:
        files_url = self.urls["files-api"]
        return self._session.post(
            f"{files_url}/files/{file_id}/download",
        )

    @http_json_response
    def _get_model_list(self) -> requests.Response:
        annotation_url = self.urls["annotation-api"]
        return self._session.get(
            f"{annotation_url}/models",
        )

    @http_json_response
    def _get_annotation_params(
        self,
        pipestance_type: str,
        genome_name: str,
        requires: list[str] | None = None,
        cell_annotation_model: str | None = None,
    ) -> requests.Response:
        annotation_url = self.urls["annotation-api"]
        args = {
            "pipelineName": pipestance_type,
            "genomes": [genome_name],
        }
        if requires is None:
            args["requires"] = []
        else:
            args["requires"] = requires

        if cell_annotation_model is not None and cell_annotation_model != "auto":
            args["modelName"] = cell_annotation_model
        return self._session.post(f"{annotation_url}/defaults", json=args)

    def _part_upload(self, upload_id: str, file_id: str, part: str, chunk_data: bytes):
        md5 = hashlib.md5(chunk_data)
        md5base64 = base64.b64encode(md5.digest()).decode("utf-8")
        if self.verbosity:
            logging.debug(
                "md5base64: %s, hex: %s, size: %d",
                md5base64,
                md5.hexdigest(),
                len(chunk_data),
            )

        part_response = self._get_part_upload_url(
            file_id, upload_id, part_num=part, md5base64=md5base64, expire_seconds=600
        )
        presigned_url = part_response["presignedUrl"]

        headers = {
            "Content-MD5": md5base64,
            "Content-Type": "application/octet-stream",
            "Accept": "application/json",
        }
        response = requests.put(presigned_url, data=chunk_data, headers=headers, timeout=15)
        if self.verbosity > 2:
            print_roundtrip(response)

        response.raise_for_status()

        if self.verbosity > 1:
            logging.debug("part upload: partNumber %d, Etag: %s", part, response.headers["ETag"])
        return response.headers["ETag"]

    def get_annotation_params(
        self,
        pipestance_type: str,
        genome_name: str,
        requires: list[str] | None = None,
        cell_annotation_model: str | None = None,
    ) -> tuple[bool, AnnotationParams | None]:
        """Get the default parameters (pipeline name, pipeline version, model name) for the supplied pipestance type, genome, and required pipeline capabilities.

        If no pipeline matches the supplied parameters, the return value will be
        false, with a null AnnotationParams object.
        """
        try:
            resp = self._get_annotation_params(
                pipestance_type,
                genome_name,
                requires=requires,
                cell_annotation_model=cell_annotation_model,
            )
            return True, AnnotationParams(**resp)
        except CloudAPIError as cloud_err:
            # if the server returns HTTP 400, that means there is no model/pipeline
            # available for the supplied inputs.
            if cloud_err.status_code == 400:
                return False, None
            else:
                raise AnnotationError(
                    "API error while looking up annotation pipeline defaults"
                ) from cloud_err
        # probably won't hit this error case, but to cover the bases...
        except requests.HTTPError as error:
            if error.response.status_code == 400:
                return False, None
            else:
                raise AnnotationError(
                    "HTTP error while looking up annotation pipeline defaults"
                ) from error

    def upload_input_file(
        self,
        input_file: str,
        file_name: str,
        metadata: dict,
        chunk_size: int = 5 * 1024 * 1024,
    ) -> str:
        """Upload Input File to 10x Cloud using aws presigned multipart upload.

        Args:
            input_file: path to input file
            file_name: name of the input file
            metadata: metadata to be stored with file
            chunk_size: chunk size for multipart upload
        Returns:
            file id.
        """
        md5 = md5sum_file(input_file)
        file_size = os.path.getsize(input_file)
        self.session_stats.file_size = file_size
        response = self._create_file_upload(input_file, file_name, file_size, md5, metadata)
        file_id, blob_id, status = (
            response["id"],
            response["blobId"],
            response["status"],
        )
        if status == "SKIP":
            return file_id
        if status != "READY":
            raise AnnotationError(f"input file is in unexpected cloud state: {status}")

        upload = self._start_file_multiupload(
            blob_id=blob_id,
        )

        upload_id = upload["uploadId"]
        parts = []
        with open(input_file, "rb") as f:
            part_num = 1
            while chunk_data := f.read(chunk_size):
                etag = self._part_upload(upload_id, blob_id, part_num, chunk_data)
                parts.append({"partNumber": part_num, "etag": etag})
                part_num = part_num + 1

        self._complete_file_multiupload(upload_id, blob_id, parts)
        self._complete_file_upload(file_id)
        return file_id

    def is_annotation_completed(self, annotation_id: str) -> bool:
        """Check if annotation is finished.

        Args:
            annotation_id: the annotation request id
        returns:
            boolean to indicate whether the annotation is finished or not.
        """
        resp = self._get_analysis(annotation_id)
        runs = resp["analysisRuns"]
        if len(runs) == 0:
            raise AnnotationError(f"analysis {annotation_id} has no runs")

        run = runs[0]
        state = run["status"]
        return not state in ("INPROGRESS", "NEW")

    def get_available_models(self) -> list[AnnotationModelDescription]:
        """Return the list of available models.

        Returns:
            An array of model description dictionaries.
        """
        resp = self._get_model_list()
        models: list[AnnotationModelDescription] = [
            AnnotationModelDescription(**model) for model in resp
        ]
        return models

    def annotate_file(
        self,
        h5_filepath: str,
        pipestance_type: str,
        genome_name: str,
        pipeline_params: AnnotationParams,
        with_distance_info: bool = False,
        with_dataset_info: bool = True,
        upload_chunk_size: int = 5 * 1024 * 1024,
    ) -> str:
        """Annotate h5 file.

        Args:
            h5_filepath: file path to h5 file
            pipestance_type: pipestance type
            genome_name: the name of the genome from the source h5
            pipeline_params: The invocation parameters returned from the server.
            with_distance_info: whether to include neighbor distances in the result output
            with_dataset_info: whether to include dataset source info in the result output
            upload_chunk_size: the chunk size used for multipart upload

        Raises:
            May raise an AnnotationError if any upload/initialize operation fails.

        Returns:
            str: unique annotation request id.
        """
        if not os.path.exists(h5_filepath):
            raise AnnotationError(f"input file {h5_filepath} does not exist")

        _, num_bcs, _ = cr_matrix.CountMatrix.load_dims_from_h5(h5_filepath)

        try:
            self.session_stats.start(CloudSessionStats.OP_UPLOAD)
            self.session_stats.num_barcodes = num_bcs
            self.session_stats.cell_annotation_model = pipeline_params.model
            file_id = self.upload_input_file(
                input_file=h5_filepath,
                file_name=os.path.basename(h5_filepath),
                metadata={"source": "cellranger"},
                chunk_size=upload_chunk_size,
            )
            self.session_stats.passed(CloudSessionStats.OP_UPLOAD)
        except requests.HTTPError as error:
            err_msg = "HTTP error while uploading input file"
            self.session_stats.failed(CloudSessionStats.OP_UPLOAD, reason=err_msg)
            raise AnnotationError(err_msg) from error
        except OSError as error:
            err_msg = f"OSError uploading input file: {error}"
            self.session_stats.failed(CloudSessionStats.OP_UPLOAD, reason=err_msg)
            raise AnnotationError(err_msg) from error

        alphabet = string.ascii_lowercase + string.digits
        analysis_name = "".join(random.choices(alphabet, k=8))

        if os.environ.get("TENX_CLOUD_ENV", "") == "true":
            analysis_source = "cloud"
        else:
            analysis_source = "local"

        req = {
            "analysisName": analysis_name,
            "description": f"cell annotation ({analysis_source})",
            "pipelineName": pipeline_params.pipelineName,
            "productVersion": pipeline_params.pipelineVersion,
            "productVariant": "default",
            "params": {
                "pipelineName": pipestance_type,
                "id": analysis_name,
                "modelName": pipeline_params.model,
                "filteredMatricesH5": file_id,
                "genomeName": genome_name,
                "withDistances": with_distance_info,
                "withDatasetInfo": with_dataset_info,
            },
        }

        analysis = self._create_adhoc_analysis(req)
        analysis_id = analysis["analysisID"]

        self.session_stats.analysis_id = analysis_id
        self.session_stats.start(CloudSessionStats.OP_ANALYZE)

        return analysis_id

    def get_annotate_status(self, annotation_id: str) -> AnnotationStatus:
        """Get Annotation Status.

        Args:
            annotation_id: the annotation id returned from annotation_file call
        Returns:
            AnnotationStatus: the status of annotation.
        """
        resp = self._get_analysis(annotation_id)
        runs = resp["analysisRuns"]
        if len(runs) > 0:
            run = runs[0]
            state = run["status"]
            if state == "COMPLETED":
                self.session_stats.passed(CloudSessionStats.OP_ANALYZE)
                project_id = resp["projectId"]
                res = self._get_analysis_output(project_id=project_id, analysis_id=annotation_id)
                outs = [AnnotationOutput(k["id"], k["filename"]) for k in res["files"]]
                if len(outs) == 0:
                    err_msg = "completed analysis has no output"
                    self.session_stats.failed(CloudSessionStats.OP_ANALYZE, reason=err_msg)
                    raise AnnotationError(err_msg)
                return AnnotationStatus(annotation_id, state, outs)
            elif state == "FAILED":
                self.session_stats.failed(
                    CloudSessionStats.OP_ANALYZE, reason="analysis pipeline failure"
                )
                return AnnotationStatus(annotation_id, state, None)
            else:
                return AnnotationStatus(annotation_id, state, None)

        else:
            raise AnnotationError("analysis has no runs")

    def get_output_download_url(self, file: AnnotationOutput) -> str:
        """Retrieve Presigned URL to download file.

        Args:
            file (AnnotationOutput): the output object returned from AnnotationStatus
        Returns:
            str: presigned download URL.
        """
        download = self._create_file_download(file.file_id)
        return download["downloadURL"]

    def preflight_check(self):
        """Preflight check."""
        # TODO to update
        self._get_current_user()
