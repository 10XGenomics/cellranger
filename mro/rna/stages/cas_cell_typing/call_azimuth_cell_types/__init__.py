#
# Copyright (c) 2025 10X Genomics, Inc. All rights reserved.
#

"""Call cell types based on Satija-Lab's Azimuth."""

import csv
import itertools
import json
import math
from dataclasses import asdict, dataclass

import martian
import numpy as np

import cellranger.matrix as cr_matrix
from cellranger.cell_typing.azimuth.annotate import AzimuthNN
from cellranger.cell_typing.azimuth.azimuth_utils import get_total_number_of_classes
from cellranger.cell_typing.common_cell_typing import CellAnnotationMetadata, CellTypeResults
from cellranger.fast_utils import total_umis_in_library
from cellranger.rna.library import GENE_EXPRESSION_LIBRARY_TYPE

__MRO__ = """
stage CALL_AZIMUTH_CELL_TYPES(
    in  string          sample_id,
    in  string          sample_desc,
    in  h5              filtered_matrix,
    in  string          cell_annotation_model,
    in  file            tenx_cloud_token_path,
    in  string          pipestance_type,
    in  bool            override_num_bc_limit,
    in  bool            is_pd,
    out CellTypeResults cell_type_results,
    out bool            cas_success,
    out bool            disable_cas_ws,
    out bool            disable_summarize,
    out json            summary,
    src py              "stages/cas_cell_typing/call_azimuth_cell_types",
) split (
    in  int             chunk_ind_start,
    in  int             chunk_ind_end,
    out csv             cell_types_csv,
) using (
    mem_gb = 4,
)
"""

PROCESSING_CHUNK_SIZE = 8_192
PROCESSING_MODEL_BATCH_SIZE = 256


@dataclass
class RunInfo:
    """Class storing run info."""

    num_features: int
    num_bcs: int
    num_entries: int
    is_cra: bool


def split(args):
    if not args.filtered_matrix:
        return {
            "chunks": [],
            "join": {"__mem_gb": 1},
        }

    genomes_in_matrix = cr_matrix.CountMatrix.get_genomes_from_h5(args.filtered_matrix)
    library_types = cr_matrix.CountMatrix.load_library_types_from_h5_file(args.filtered_matrix)
    is_human = any("GRCh38" in item or "hg19" in item for item in genomes_in_matrix)
    if (
        len(genomes_in_matrix) > 1
        or not is_human
        or GENE_EXPRESSION_LIBRARY_TYPE not in library_types
    ):
        martian.log_info(
            f"Exiting because genome: {genomes_in_matrix}, library types: {library_types} and is_human: {is_human}"
        )
        return {
            "chunks": [],
            "join": {"__mem_gb": 1},
        }

    non_zero_entries = total_umis_in_library(
        matrix_path=args.filtered_matrix,
        library=GENE_EXPRESSION_LIBRARY_TYPE,
    )
    martian.log_info(f"non-zero GEX entries = {non_zero_entries}")
    if non_zero_entries == 0:
        martian.log_info("Disabling because there are zero non-zero GEX entries!")
        return {
            "chunks": [],
            "join": {"__mem_gb": 1},
        }

    num_features, num_bcs, _nnz = cr_matrix.CountMatrix.load_dims_from_h5(args.filtered_matrix)
    matrix_mem_gb = cr_matrix.CountMatrix.get_mem_gb_from_matrix_h5(args.filtered_matrix)

    processing_boundaries = list(np.arange(0, num_bcs, PROCESSING_CHUNK_SIZE)) + [num_bcs]

    total_number_of_classes = get_total_number_of_classes()
    print(f"{total_number_of_classes=}")

    # 1GB=134_217_728 integers of 64 bits
    dense_matrix_mem_gb = math.ceil(PROCESSING_MODEL_BATCH_SIZE * num_features / 134_000_000)
    softmax_mem_gb = math.ceil(2 * PROCESSING_CHUNK_SIZE * total_number_of_classes / 134_000_000)

    chunks = [
        {
            "chunk_ind_start": int(start_ind),
            "chunk_ind_end": int(end_ind),
            "__mem_gb": matrix_mem_gb + dense_matrix_mem_gb + softmax_mem_gb,
            "__vmem_gb": 128,
        }
        for (start_ind, end_ind) in itertools.pairwise(processing_boundaries)
    ]

    return {
        "chunks": chunks,
        "join": {
            "__mem_gb": 1,
        },
    }


def main(args, outs):
    feature_bc_matrix = (
        cr_matrix.CountMatrix.load_h5_file(args.filtered_matrix)
        .select_features_by_types([GENE_EXPRESSION_LIBRARY_TYPE])
        .select_barcodes(np.arange(args.chunk_ind_start, args.chunk_ind_end))
    )
    chunk_fbm_h5_path = martian.make_path("chunk_fbm.h5").decode()
    feature_bc_matrix.save_h5_file(chunk_fbm_h5_path)

    azimuth = AzimuthNN(eval_batch_size=PROCESSING_MODEL_BATCH_SIZE)
    azimuth.query_feature_barcode_matrix_h5(matrix_path=chunk_fbm_h5_path)
    azimuth.process_query()  # get the inputs in format that the model can be run
    _ = azimuth.run_inference_model()  # run the model
    _ = azimuth.calibrate_predictions()
    _ = azimuth.process_outputs()  # postprocess model outputs
    azimuth.azimuth_refine()

    azimuth.write_results_to_csv(outs.cell_types_csv, fill_unannotated=True)
    num_bcs = len(feature_bc_matrix.bcs)
    num_bcs_annotated = len(azimuth.process_outputs()["full_hierarchical_labels"])
    frac_returned_bcs = float(num_bcs_annotated / num_bcs)
    print(f"frac_returned_bcs={frac_returned_bcs}")


def join(args, outs, _chunk_defs, chunk_outs):
    error_cell_type_results = CellTypeResults(
        model_used="Azimuth",
        skip_downstream=True,
    ).to_outs()
    if not args.filtered_matrix:
        martian.log_info("No matrix specified. Exiting!")
        martian.clear(outs)
        outs.cell_type_results = error_cell_type_results
        outs.cas_success = False
        outs.disable_cas_ws = True
        outs.disable_summarize = True
        return

    # Need to subset matrix that have multiple library types into gene expression only
    if not chunk_outs:
        martian.log_info("No valid chunks to process. Exiting!")
        martian.clear(outs)
        outs.cell_type_results = error_cell_type_results
        outs.cas_success = False
        outs.disable_cas_ws = True
        outs.disable_summarize = True
        return

    _num_features, num_bcs, _nnz = cr_matrix.CountMatrix.load_dims_from_h5(args.filtered_matrix)
    first_viable_csg = next(chnk.cell_types_csv for chnk in chunk_outs if chnk.cell_types_csv)

    with open(first_viable_csg) as f:
        csv_reader = csv.DictReader(f)
        headers = csv_reader.fieldnames

    cell_types_csv_path = martian.make_path("cell_types.csv").decode()
    num_bcs_returned = 0
    with open(cell_types_csv_path, "w", newline="") as out_csv:
        writer = csv.DictWriter(out_csv, fieldnames=headers)
        writer.writeheader()
        for chnk in chunk_outs:
            if not chnk.cell_types_csv:
                continue
            with open(chnk.cell_types_csv) as in_csv:
                reader = csv.DictReader(in_csv)
                for row in reader:
                    writer.writerow(row)
                    num_bcs_returned += 1

    metadata = CellAnnotationMetadata(
        developer="Satija Lab", tree_version_used="PanHuman", is_beta_model=True
    )
    metadata_json_path = martian.make_path("metadata.json").decode()
    with open(metadata_json_path, "w") as f:
        json.dump(asdict(metadata), f, indent=4)

    cell_type_results = CellTypeResults(
        cell_types=cell_types_csv_path,
        frac_returned_bcs=num_bcs_returned / num_bcs,
        model_used=f"Azimuth {AzimuthNN.get_azimuth_version()}",
        metadata=metadata_json_path,
    )
    print(f"cell_type_results={cell_type_results}")

    outs.cell_type_results = cell_type_results.to_outs()
    outs.cas_success = True
    outs.disable_cas_ws = False
    outs.disable_summarize = False
    outs.summary = None
