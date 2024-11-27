# Copyright (c) 2019 10x Genomics, Inc. All rights reserved.

from __future__ import annotations

import csv
import os
from collections import Counter

import numpy as np
from six import ensure_binary, ensure_str


class CSVParseException(Exception):
    def __init__(self, msg):
        super().__init__(msg)
        self.msg = msg

    def __str__(self):
        return self.msg


class CSVEmptyException(CSVParseException):
    """The CSV file has no data rows."""


def _iter_file_rows(filename_bytes, filename, descriptive_name):
    """Iterates through non-comment lines in the file, returning them as unicode strings.

    Also skips the BOM if present.

    Raises:
        CSVParseException: The file has no non-comment lines or some lines can't be decoded as ascii.
    """
    with open(filename_bytes, encoding="utf-8-sig", newline=None) as f:
        # encoding utf-8-sig will skip BOM if present.
        has_data = False
        # First non-comment row is treated specially, because it is the header.
        first = True
        try:
            for row in f:
                row = row.strip()
                if not row:
                    # Skip blank rows.  The csv reader does this anyway, but doing
                    # it here saves us the trouble of decoding from ascii.
                    continue
                if not row.startswith("#"):
                    has_data = True
                    if first:
                        # Remove whitespace around column headers
                        row = ",".join(x.strip() for x in row.split(","))
                        first = False
                    # Verify is ASCII Compatible
                    if not row.isascii():
                        raise CSVParseException(
                            f"The {descriptive_name} csv file {filename} contains non-ascii characters.\n\nRow:\n{row}"
                        )
                    yield row
        # We could encounter invalid unicode when decoding any row, which has
        # some behind-the-scenes buffering. Provide an error message that
        # attempts to provide better context, showing the surrounding ~60 characters.
        except UnicodeDecodeError as err:
            context_start = max(err.start - 30, 0)
            badchar_offset = err.start - context_start
            snippet = (
                err.object[context_start : err.start + 30]
                .decode("utf-8", errors="replace")
                .replace("\r", " ")
                .replace("\n", " ")
            )
            caret = (" " * badchar_offset) + "^"
            raise CSVParseException(
                f"The {descriptive_name} csv file {filename} has one or more invalid utf-8 characters: {err.reason}. "
                f"The first bad character is at absolute position {err.start}.\n"
                "This snippet shows the offending character in context:\n"
                f"{snippet}\n{caret}"
            ) from err

    if not has_data:
        raise CSVEmptyException(f"The {descriptive_name} csv file {filename} has no data.")


def _validate_columns(reader, descriptive_name, required_cols, valid_cols):
    col_names = reader.fieldnames
    col_set = frozenset(col_names)
    if not col_set.issuperset(set(required_cols)):
        raise CSVParseException(
            """The {} file header does not contain one or more required comma-separated fields: "{}".
The following fields were found: "{}".
Please check that your file is in CSV format and has the required field names.""".format(
                descriptive_name, ", ".join(required_cols), ", ".join(col_names)
            )
        )

    if valid_cols:
        valid_set = frozenset(valid_cols)
        if not col_set.issubset(valid_set):
            extra_columns = sorted(col_set - valid_set)

            msg = """The {} file header contains invalid columns: "{}".
Only the following comma-delimited fields may be used in the file header: "{}".""".format(
                descriptive_name, ", ".join(extra_columns), ", ".join(valid_cols)
            )

            if required_cols:
                msg += """
The following columns are required: "{}".""".format(
                    ", ".join(required_cols)
                )

            msg += f"""
Please check that you have formatted the {descriptive_name} file correctly and included the appropriate column headers."""

            raise CSVParseException(msg)

    if len(col_set) < len(col_names):
        msg = "{} csv has a duplicated column: {}".format(
            descriptive_name,
            ", ".join(sorted(name for name, count in Counter(col_names).items() if count > 1)),
        )
        raise CSVParseException(msg)


def load_csv_filter_comments(
    filename: os.PathLike, descriptive_name, required_cols, valid_cols=None
):
    """Returns non-comment lines from the csv files.

    Verifies ASCII encoding and no duplicate columns.

    Args:
        filename: The csv file to open.
        descriptive_name (str): The description of the file, to use in error messages.
        required_cols (sequence of str): columns which must be present.
        valid_cols (sequence of str): columns which may be present.  If provided,
            then it is an error for other columns to be present.

    Returns:
        csv.DictReader: The non-comment rows.

    Raises:
        CSVParseException
    """
    filename_bytes = ensure_binary(os.fspath(filename))
    if not os.path.isfile(filename_bytes):
        raise CSVParseException(f"Could not find the {descriptive_name} csv file {filename}")

    if not os.access(filename_bytes, os.R_OK):
        msg = (
            f"The {descriptive_name} csv is not readable, please check file permissions: {filename}"
        )
        raise CSVParseException(msg)

    reader = csv.DictReader(_iter_file_rows(filename_bytes, filename, descriptive_name))
    _validate_columns(reader, descriptive_name, required_cols, valid_cols)
    return reader


def write_filtered_barcodes(out_csv: os.PathLike, bcs_per_genome: dict[str, list[str]]):
    """Write the barcodes to a CSV.

    Args:
        bcs_per_genome (dict of str to list): Map each genome to its cell-associated barcodes
    """
    with open(out_csv, "w") as f:
        writer = csv.writer(f, lineterminator="\n")
        for genome, bcs in bcs_per_genome.items():
            for bc in bcs:
                writer.writerow([ensure_str(genome), ensure_str(bc)])


def write_isotype_normalization_csv(
    out_csv: str | bytes,
    barcode_list: np.ndarray,
    in_tissue: list[bool],
    normalization_array: np.ndarray,
):
    """Write the normalization constant for each barcode to a CSV.

    Args:
        out_csv (Union[str, bytes]): Path of the CSV
        barcode_list (np.ndarray): list of barcodes
        in_tissue (list[bool]): boolean list of in tissue
        normalization_array (np.ndarray): normalization for each barcode - in the same order
    """
    assert barcode_list.shape == normalization_array.shape, (
        "Lengths of barcodes and normalization factors"
        + "do not match. "
        + f"barcode list length: {barcode_list.shape}. "
        + f"Normalization array shape: {normalization_array.shape}."
    )
    assert barcode_list.shape[0] == len(in_tissue), (
        "Lengths of barcodes and in tissue vector"
        + "do not match. "
        + f"barcode list length: {barcode_list.shape}. "
        + f"In tissue length: {len(in_tissue)}."
    )
    with open(out_csv, "w") as f:
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(["barcode", "in_tissue", "normalization_factor"])
        for barcode, bc_in_tissue, normalized_value in zip(
            barcode_list, in_tissue, normalization_array
        ):
            writer.writerow(
                [
                    ensure_str(barcode),
                    str(int(bc_in_tissue)),
                    np.format_float_positional(normalized_value, precision=6),
                ]
            )


def combine_csv(input_csvs, output_csv, header_lines=1):
    """Combine a list of CSV files.

    Files specified in input_csvs are combined into
    a single output csv in output_csv. It is assumed that all
    CSVs have the same header structure. The number of header
    lines is specified in header_lines.
    """
    if not input_csvs:
        output_csv = None
        return
    with open(output_csv, "w") as out:
        for i, icsv in enumerate(input_csvs):
            with open(icsv) as infile:
                header = []
                for _ in range(header_lines):
                    header.append(next(infile))
                if i == 0:
                    for hline in header:
                        out.write(hline)
                for line in infile:
                    out.write(line)
