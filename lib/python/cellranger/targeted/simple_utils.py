#
# Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
#
"""Simple methods that don't require a lot of dependencies.

For things which might be needed in things like preflights.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import TYPE_CHECKING

from six import ensure_binary, ensure_str

import cellranger.csv_utils as cr_csv_utils
from cellranger.targeted.targeted_constants import (
    OLIGO_NAME,
    TARGETING_METHOD_FILE_FORMAT_REGEX,
    TARGETING_METHOD_FILE_NAMES,
    TARGETING_METHOD_HC,
    TARGETING_METHOD_HC_FILE_FORMAT,
    TARGETING_METHOD_TL,
    TARGETING_METHOD_TL_ALLOWED_REGIONS,
    TARGETING_METHOD_TL_FILE_FORMAT,
    TargetingMethod,
)

if TYPE_CHECKING:
    from cellranger.reference import GeneIndex, NewGeneIndex


def parse_bool(astring: str):
    """Return True when string is "true" (case insensitive) and False when "false" and None otherwise."""
    abool = astring.casefold()
    return abool == "true" or (False if abool == "false" else None)


def parse_target_csv(
    filename: Path | str | bytes,
    ref_gene_index: GeneIndex | NewGeneIndex | None = None,
    gene_name_to_id: dict[bytes, bytes] | None = None,
    filter_probes: bool = True,
    expected_targeting_method: str | None = None,
) -> tuple[dict[str, str], list[bytes], list[tuple[bytes, bytes]]]:
    """Return a list of gene_ids and a list (gene_id, bait_seq) tuples given a CSV input file.

    Args:
        filename: name of target panel CSV file
        ref_gene_index (GeneIndex): optional GeneIndex object -- for example,
            returned by NewGeneIndex.load_from_reference(my_reference_path).
            If using, then gene_name_to_id is also required.
        gene_name_to_id (dict): optional dict with gene names as keys and gene
            ids as values also for gene name conversion support. If using then
            ref_gene_index is also required.
        filter_probes (bool): flag indicating whether probe filtering was enabled
        expected_targeting_method (str): name of targeting method expected based on cmd line arg

    Returns:
        dict[str,str]: header metadata
        list[bytes]: gene_ids
        list[(bytes,bytes)]: gene_id, bait_seq
    """
    expected_descriptive_name = TargetingMethod.get_file_name(expected_targeting_method)

    method_info = determine_targeting_method_info_from_csv(filename, expected_targeting_method)
    descriptive_name = method_info.descriptive_name
    targeting_method = method_info.method
    if expected_targeting_method is not None and targeting_method != expected_targeting_method:
        raise cr_csv_utils.CSVParseException(
            f"The supplied file {filename} is a {descriptive_name}, but we were expecting a {expected_descriptive_name}."
        )

    (
        valid_cols,
        required_cols,
        required_metadata,
        conflicting_metadata,
    ) = get_target_panel_or_probe_set_file_format_spec(method_info.method, method_info.file_version)

    # Start reader and do basic checks on file format
    reader = cr_csv_utils.load_csv_filter_comments(
        filename,
        descriptive_name=descriptive_name,
        required_cols=required_cols,
        valid_cols=valid_cols,
    )

    check_target_csv_metadata(filename, descriptive_name, required_metadata, conflicting_metadata)

    valid_bases = b"ATGCN-"  # "-" used as separator between probe halves (and if present gap)
    min_bait_length = 1

    # Track set of unique gene IDs and gene_id, bait_sequence tuples
    gene_included: dict[bytes, bool] = {}
    bait_sequences: list[tuple[bytes, bytes]] = []
    bait_ids: set[bytes] = set()

    if (ref_gene_index is None and gene_name_to_id is not None) or (
        gene_name_to_id is None and ref_gene_index is not None
    ):
        raise ValueError(
            "ref_gene_index and gene_name_to_id must either both be unspecified "
            "or both be specified, not just one of the two."
        )

    ## Determine offset for reporting line numbers
    metadata_fields = load_target_csv_metadata(filename, descriptive_name)
    offset = 2 + len(metadata_fields)

    # More detailed per-row checks
    last_seq_length = None
    for row_num, entry in enumerate(reader):
        # Check that there aren't extra or missing columns in the CSV
        for key in entry:
            if key is None:
                msg = (
                    "Your %s file contains more columns than the header on row %d. Please use a csv file with a header for each column.\nYou might have an comma character in a field. Commas are permitted in some fields, but fields containing commas must be enclosed in quotes."
                    % (descriptive_name, (row_num + offset))
                )
                raise cr_csv_utils.CSVParseException(msg)

            if entry[key] is None or entry[key] == "":
                msg = (
                    "Your %s file contains an empty column or fewer columns than the header on row %d. You might have a missing comma. Please use a csv file with a header for each column and a value for each column in each row."
                    % (descriptive_name, (row_num + offset))
                )
                raise cr_csv_utils.CSVParseException(msg)

        gene_id: bytes = ensure_binary(entry["gene_id"])
        # Track gene IDs and bait sequences if provided
        if ref_gene_index is not None:
            # both not None
            assert gene_name_to_id is not None
            gene_id = sanitize_gene_entry(gene_id, ref_gene_index, gene_name_to_id)
            if (
                TARGETING_METHOD_TL_FILE_FORMAT not in metadata_fields
                and ref_gene_index.gene_id_to_int(gene_id) is None
            ):
                raise ValueError(f"Gene {ensure_str(gene_id)} not seen in reference")
        elif not ref_gene_index is None or not gene_name_to_id is None:
            raise ValueError(
                "ref_gene_index and gene_name_to_id must either both be unspecified or both be specified, not just one of the two."
            )

        # Check for valid region if region is part of valid_cols

        probe_region = entry.get("region")
        if probe_region is not None and probe_region not in TARGETING_METHOD_TL_ALLOWED_REGIONS:
            raise ValueError(
                f"Line {row_num + offset} has a region that is not accepted. Allowed regions are: {','.join(TARGETING_METHOD_TL_ALLOWED_REGIONS)}"
            )

        included_string = entry.get("included", "true")
        included = parse_bool(included_string)
        if included is None:
            raise ValueError(
                f'The column "included" must be "true" or "false" and saw "{ensure_str(included_string)}"'
            )
        gene_included[gene_id] = gene_included.get(gene_id, True) and included

        bait_id = entry.get("bait_id", None) or entry.get("probe_id", None)
        if not bait_id is None:
            bait_id: bytes = ensure_binary(bait_id)
            if bait_id not in bait_ids:
                bait_ids.add(bait_id)
            else:
                msg = "Your {} file contains a duplicate {} ID on row {}: '{}'. All entries in the {}_id column must be unique to their row.".format(
                    ensure_str(descriptive_name),
                    OLIGO_NAME.get(targeting_method, "probe or bait"),
                    (row_num + offset),
                    ensure_str(bait_id),
                    OLIGO_NAME.get(targeting_method, "probe_id or bait"),
                )
                raise cr_csv_utils.CSVParseException(msg)

        bait_probe_seq = entry.get("bait_seq", None) or entry.get("probe_seq", None)
        if (
            bait_probe_seq is not None
            and targeting_method == TARGETING_METHOD_TL
            and float(method_info.file_version) < 3.0
        ):
            if last_seq_length is None:
                last_seq_length = len(bait_probe_seq)
            elif last_seq_length != len(bait_probe_seq):
                msg = "The {} file contains {} sequences of different lengths starting on on row {}. Please use a csv file where all {} sequences are of the same length.".format(
                    descriptive_name,
                    OLIGO_NAME.get(targeting_method, "oligo"),
                    (row_num + offset),
                    OLIGO_NAME.get(targeting_method, "oligo"),
                )
                raise cr_csv_utils.CSVParseException(msg)
            else:
                last_seq_length = len(bait_probe_seq)

        # Check valid DNA sequences for baits
        bait_seq: bytes = ensure_binary(
            (entry.get("bait_seq") or entry.get("probe_seq", DUMMY_FASTA_SEQUENCE)).upper()
        )
        if not all(bait_seq[i : i + 1] in valid_bases for i in range(len(bait_seq))):
            raise cr_csv_utils.CSVParseException(
                f'Invalid {OLIGO_NAME.get(targeting_method, "oligo")} sequence '
                f"found in {ensure_str(descriptive_name)} file: "
                f'"{bait_seq.decode()}" on row {row_num + offset}. '
                f'May only contain "{", ".join(chr(x) for x in sorted(valid_bases))}".'
            )

        if len(bait_seq) < min_bait_length:
            raise cr_csv_utils.CSVParseException(
                'Invalid {} sequence found in {} file: "{}" on row {}. {}bp length is shorter than allowed minimum of {}bp.'.format(
                    OLIGO_NAME.get(targeting_method, "oligo"),
                    descriptive_name,
                    bait_seq,
                    row_num + offset,
                    len(bait_seq),
                    min_bait_length,
                )
            )
        if included:
            bait_sequences.append((gene_id, bait_seq))

    if filter_probes:
        gene_ids = {gene_id for gene_id, included in gene_included.items() if included}
        if len(gene_included) and not len(gene_ids):
            raise cr_csv_utils.CSVParseException(
                'All probes in the probe set file are set to false for "included" field. At least 10 probes need to have true values for included.'
            )
    else:
        gene_ids = set(gene_included.keys())

    return (
        metadata_fields,
        sorted(gene_ids),
        sorted(bait_sequences, key=lambda bait_tuple: bait_tuple[0]),
    )


def get_target_panel_or_probe_set_file_format_spec(targeting_method, file_version):
    """Returns valid and required columns and metadata (header) fields for probe set or target panel files.

    Args:
        targeting_method (str): The targeting method used.
        file_version (str): X.Y file version string matching regex TARGETING_METHOD_FILE_FORMAT_REGEX

    Returns:
        list : valid columns
        list : required columns
        set  : required metadata fields
    """
    required_cols = ["gene_id"]
    required_metadata = {
        "reference_genome",
        "reference_version",
        "panel_type",
        "panel_name",
    }
    conflicting_metadata = {}
    valid_cols = ["gene_id"]
    if targeting_method == TARGETING_METHOD_HC:
        required_metadata.update({TARGETING_METHOD_HC_FILE_FORMAT})
        valid_cols.extend(["bait_seq", "bait_id"])
        conflicting_metadata.update(
            {TARGETING_METHOD_TL_FILE_FORMAT: TARGETING_METHOD_HC_FILE_FORMAT}
        )
    elif targeting_method == TARGETING_METHOD_TL:
        required_metadata.update({TARGETING_METHOD_TL_FILE_FORMAT})
        valid_cols.extend(["probe_seq", "probe_id", "included"])
        if float(file_version) >= 2.0:  # region introduced in 2.0
            valid_cols.append("region")
        if float(file_version) >= 3.0:  # gene_name introduced in 3.0
            valid_cols.append("gene_name")
        required_cols.extend(["probe_seq", "probe_id"])
        conflicting_metadata.update(
            {TARGETING_METHOD_HC_FILE_FORMAT: TARGETING_METHOD_TL_FILE_FORMAT}
        )

    return valid_cols, required_cols, required_metadata, conflicting_metadata


def load_target_csv_metadata(filename: str, descriptive_name: str) -> dict[str, str]:
    """Return a dictionary of metadata fields contained in header of target panel CSV file.

    Args:
        filename (str): name of target panel CSV file
        descriptive_name (str): descriptive name for file for use in error messages

    Returns:
        dict (str, str): dictionary mapping metadata keys to values
    """
    metadata = {}
    for row in open(filename):  # TODO need to make this parsing more robust...
        row = row.strip()
        if row.startswith("#"):
            k_v_pair = row.strip("#").split("=")
            if len(k_v_pair) != 2:
                raise cr_csv_utils.CSVParseException(
                    f'Invalid metadata format detected in {descriptive_name} file: "{row}". Must follow "#field=value" format.'
                )
            k, v = k_v_pair
            k, v = k.strip(), v.strip()
            metadata[k] = v
        else:
            break

    return metadata


def check_target_csv_metadata(
    filename, descriptive_name, required_metadata, conflicting_metadata
) -> float:
    """[summary].

    Args:
        filename (str): Path to target panel or probe set file
        descriptive_name (str): descriptive name (probe set or target panel)
        required_metadata (set): the required metadata (header) fields
        conflicting_metadata (dict): any conflicting metadata (header) fields

    Raises:
        cr_csv_utils.CSVParseException: if metadata fields dont comply with file format spec
    """
    # Load metadata fields
    metadata_fields = load_target_csv_metadata(filename, descriptive_name)
    missing_required_fields = required_metadata.difference(metadata_fields)
    if missing_required_fields:
        raise cr_csv_utils.CSVParseException(
            f"The following metadata fields are required in the {descriptive_name} header: "
            '"{}", but were not found. Please include these fields in #field=value format at the top of the {} file.'.format(
                ", ".join(sorted(missing_required_fields)),
                descriptive_name,
            )
        )
    if "panel_name" in required_metadata:
        panel_name = metadata_fields["panel_name"]
        if "/" in panel_name:
            raise cr_csv_utils.CSVParseException(
                f'The character "/" cannot appear in the target panel or probe set CSV panel name: {panel_name}'
            )

    method_info = determine_targeting_method_info_from_metadata(metadata_fields)

    for field1, field2 in conflicting_metadata.items():
        if field1 in metadata_fields and field2 in metadata_fields:
            raise cr_csv_utils.CSVParseException(
                f'The "{field1}" metadata field in the {descriptive_name} CSV file header conflicts with the "{field2}" metadata field. Please check the {descriptive_name} csv file.'
            )

    if float(method_info.file_version) > 3.0:
        raise cr_csv_utils.CSVParseException(
            f'The {descriptive_name} file {filename} contains an unknown {method_info.file_format_tag}: "{method_info.file_version}". Must be 3.0 or less.'
        )
    return float(method_info.file_version)


def determine_targeting_method_info_from_metadata(
    target_set_metadata, expected_targeting_method=None
):
    """Determines the targeting type (either TL or HC) and the version from metadata read from load_target_csv_metadata()."""
    method_info = TargetingMethod.get_targeting_method_from_metadata(target_set_metadata)

    if method_info is None:
        expected_file_name = TargetingMethod.get_file_name(expected_targeting_method)
        expected_file_format = TargetingMethod.get_file_format(expected_targeting_method)

        raise cr_csv_utils.CSVParseException(
            f"The provided {expected_file_name} does not include the required {expected_file_format} metadata."
        )

    if re.match(TARGETING_METHOD_FILE_FORMAT_REGEX, method_info.file_version) is None:
        raise cr_csv_utils.CSVParseException(
            f'The {method_info.descriptive_name} file contains an invalid value for {method_info.file_format_tag}: "{method_info.file_version}". {method_info.file_format_tag} must conform to X.Y formatting (e.g. 1.0)'
        )

    return method_info


def determine_targeting_method_info_from_csv(filename, expected_targeting_method):
    """Reads the target panel csv metadata to determine the targeting type and returns the format version."""
    expected_descriptive_name = TARGETING_METHOD_FILE_NAMES.get(
        expected_targeting_method, "target panel or probe set"
    )
    target_set_metadata = load_target_csv_metadata(filename, expected_descriptive_name)

    return determine_targeting_method_info_from_metadata(
        target_set_metadata, expected_targeting_method
    )


def determine_targeting_type_from_csv(filename):
    """Reads the target panel csv metadata to determine the targeting type."""
    method_info = determine_targeting_method_info_from_csv(filename, "target panel or probe set")

    return method_info.method


DUMMY_FASTA_SEQUENCE = b"N" * 120


def sanitize_gene_entry(
    raw_entry: bytes, ref_gene_index: GeneIndex | NewGeneIndex, gene_name_to_id: dict[bytes, bytes]
) -> bytes:
    """Returns gene_id that corresponds to a gene_id in the specified reference or raw_entry otherwise.

    Args:
        raw_entry (str): the raw entry user provided for the gene
        ref_gene_index (GeneIndex): GeneIndex object -- for example,
            returned by NewGeneIndex.load_from_reference(my_reference_path).
        gene_name_to_id (dict): dict with gene names as keys and gene ids as
            values also for gene name conversion support
    Returns:
        str: gene_id that occurs in the reference or raw_entry otherwise
    """
    # Could be transcript ID, convert to gene if so
    if raw_entry.split(b".")[0] in ref_gene_index.transcripts:
        gene_id = ref_gene_index.get_gene_from_transcript(raw_entry.split(b".")[0].decode()).id
        return gene_id
    else:
        # we also support gene_id, gene_id without suffix, and gene_name
        possible_gene_ids: list[bytes | None] = [
            raw_entry,
            raw_entry.split(b".")[0],
            gene_name_to_id.get(raw_entry, None),
        ]

        for gene_id in possible_gene_ids:
            if gene_id is not None and ref_gene_index.gene_id_to_int(gene_id) is not None:
                return gene_id

        return raw_entry


######################################################################
# Utility functions related to target panel file
######################################################################


def write_bait_fasta(filename, bait_sequences):
    """Write a fasta with bait sequences from (gene_id, bait_seq) tuples.

    Args:
        filename (str): file name for desired fasta file
        bait_sequences (list[tuple[str,str]]): list of tuples containing (gene_id, bait_sequence)
    """
    with open(filename, "w") as output:
        for i, (gene_id, bait_seq) in enumerate(bait_sequences):
            output.write(f">{gene_id}|{i}\n{bait_seq}\n")
