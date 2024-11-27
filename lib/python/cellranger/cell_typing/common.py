#!/usr/bin/env python3
#
# Copyright (c) 2023 10X Genomics, Inc. All rights reserved.
#
"""Cell Typing utility functions."""
import json
import os
import pathlib
from dataclasses import asdict, dataclass

import cellranger.matrix as cr_matrix
from cellranger.rna.library import GENE_EXPRESSION_LIBRARY_TYPE

ALLOWED_GENOMES = ["GRCh38", "hg19", "mm10", "GRCm39"]
SCORE_THRESHOLD = 0.8

CAS_GENOME_REF_BASE_PATHS = pathlib.Path(__file__).parent / "cas_genome_references"
CAS_CSV_PATHS = {
    "GRCh38": CAS_GENOME_REF_BASE_PATHS / "cellarium_cas_tx_pca_002_grch38_2020_a.csv",
    "mm10": CAS_GENOME_REF_BASE_PATHS / "refdata-gex-mm10-2020-A.csv",
}

# pipeline requirements
CELL_ANNOTATION_REQUIRES_JSON_RESULTS = "results:json"


@dataclass
class CellTypeResults:
    """Common structure for results from remote cell annotation stages."""

    cell_types: str | None = None
    results: str | None = None
    metadata: str | None = None
    frac_returned_bcs: float = 0.0
    model_used: str | None = None
    skip_downstream: bool = False

    def to_outs(self):
        """Return the value this form needs to take to be in stage outs."""
        return asdict(self)


def check_h5_type(filtered_matrix_path: str) -> str:
    """Checks the matrix filetype by investigating the suffix.

    Args:
        filtered_matrix_path (str): The path to the matrix file.

    Returns:
        str: The type of the matrix file ('h5' or 'h5ad').
    """
    matrix_type = pathlib.Path(filtered_matrix_path).suffix
    if matrix_type not in [".h5", ".h5ad"]:
        raise ValueError(f"Unsupported matrix type: {matrix_type}")
    return matrix_type.replace(".", "")


def check_valid_genome_and_library(matrix: str) -> tuple[bool, bool, str | None]:
    """Checks if the given matrix contains valid genome/species and library type.

    Args:
        matrix (str): Path to the matrix to be checked.

    Returns:
        tuple[bool, bool, str]: True if the matrix meets the criteria for valid genome and library information.
        True if the genome is determined to be human. Genome name if there is a single genome.
    """
    genomes = cr_matrix.CountMatrix.get_genomes_from_h5(matrix)
    has_grch38_mm10 = any(
        any(allowed_genome in genome for genome in genomes) for allowed_genome in ALLOWED_GENOMES
    )
    library_types = cr_matrix.CountMatrix.load_library_types_from_h5_file(matrix)
    is_human = any("GRCh38" in item or "hg19" in item for item in genomes)
    return (
        (len(genomes) == 1 and has_grch38_mm10 and GENE_EXPRESSION_LIBRARY_TYPE in library_types),
        is_human,
        genomes[0] if len(genomes) == 1 else None,
    )


def read_cli_token(path) -> tuple[str | None, str | None]:
    """Returns the server URL and CLI token at the specified path."""
    url = os.environ.get("TENX_CLOUD_URL", "https://cloud.10xgenomics.com")
    try:
        with open(path) as infile:
            cli_token = infile.read()
            cli_token = cli_token.strip()
        return url, cli_token
    except Exception as err:  # pylint: disable=broad-except
        # could not read the token
        print(f"Error while reading the token at {path}.\nError {err}")
        return None, None


def read_full_token(path) -> tuple[str | None, dict | None]:
    """Returns the server URL and token info contained in the full token JSON at the specified path.

    If the token at the specified path isn't JSON, return None.
    """
    try:
        with open(path) as infile:
            token_info: dict = json.load(infile)
            # should only be one key -- the server url
            url = next(iter(token_info.keys()))
            return url, token_info[url]
    except json.JSONDecodeError as err:
        print(f"Error decoding token JSON at {path}.\n error seen {err}")
        return None, None
