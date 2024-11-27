#!/usr/bin/env python3
#
# Copyright (c) 2017 10x Genomics, Inc. All rights reserved.
#

# TODO(Spatial team): This stage is no longer used by Count AGGR/Reanalyze at all, and
# parse_aggr_csv.rs is used instead.
# It's currently only used in Spatial Aggr, so we should extend parse_aggr_csv.rs
# to handle Spatial Aggr CSV parsing and get rid of this stage.

from __future__ import annotations

import codecs
import csv
import os
from pathlib import Path
from typing import NamedTuple

import martian

import cellranger.constants as cr_constants
import cellranger.hdf5 as cr_h5
import cellranger.matrix as cr_matrix
import cellranger.molecule_counter as cr_mc
import cellranger.spatial.spatial_aggr_files as sa_files

__MRO__ = """
stage PARSE_CSV(
    in  path   pipestance_root,
    in  csv    aggregation_csv,
    in  bool   reanalyze,
    in  h5     matrix_h5,
    in  string product_type,
    out csv    aggregation_csv,
    out map[]  sample_defs,
    src py     "stages/aggregator/parse_csv",
)
"""


class agg_fields(NamedTuple):
    type: str
    required: bool


class agg_folder_fields(NamedTuple):
    """Aggr folder fields.

    Structured namedtuple for aggr fields
    """

    required: bool
    type: str
    files: list[str]
    default: str


AGG_SAMPLE_ID_FIELD = "sample_id"
META_PARSE_CSV_FIELDS = {
    cr_constants.SPATIAL_PRODUCT_TYPE: {
        cr_constants.AGG_ID_FIELD: agg_fields(type="string", required=True),
        cr_constants.AGG_H5_FIELD: agg_fields(type="file_path", required=True),
        cr_constants.AGG_BATCH_FIELD: agg_fields(type="string", required=False),
        cr_constants.AGG_CLOUPE_FIELD: agg_fields(type="file_path", required=True),
        sa_files.AGG_SPATIAL_FIELD: agg_folder_fields(
            type="folder_path",
            required=False,
            files=[
                sa_files.AGG_TISSUE_POSITION_FIELD,
                sa_files.AGG_SCALE_FACTORS_FIELD,
                sa_files.AGG_HIRES_IMAGES_FIELD,
                sa_files.AGG_LOWRES_IMAGES_FIELD,
            ],
            default="spatial",
        ),
    },
    cr_constants.SINGLE_CELL_PRODUCT_TYPE: {
        AGG_SAMPLE_ID_FIELD: agg_fields(type="string", required=True),
        cr_constants.AGG_H5_FIELD: agg_fields(type="file_path", required=True),
        cr_constants.AGG_BATCH_FIELD: agg_fields(type="string", required=False),
    },
}

META_PARSE_CSV_FILES = {
    cr_constants.SPATIAL_PRODUCT_TYPE: sa_files.SPATIAL_AGGR_FILES,
    cr_constants.SINGLE_CELL_PRODUCT_TYPE: cr_constants.SC_AGGR_FILES,
}


def main(args, outs):
    """PARSE_CSV checks the header of the csv file, existance of provided paths/files for aggregation."""
    if args.reanalyze and not args.aggregation_csv:
        return  # CSV not required for reanalyze

    # Get the fields to check for in the csv
    if args.product_type not in META_PARSE_CSV_FIELDS:
        martian.exit(
            f"Product type {args.product_type} for aggr not recognized. Needs to be defined in META_PARSE_CSV_FIELDS"
        )
    else:
        field_validation = META_PARSE_CSV_FIELDS[args.product_type]

    # Get file paths
    aggr_files = META_PARSE_CSV_FILES[args.product_type]
    outs.sample_defs = parse_sample_sheet(
        args.pipestance_root, args.aggregation_csv, field_validation, aggr_files
    )

    if args.reanalyze and args.matrix_h5:
        library_map = cr_matrix.get_gem_group_index(args.matrix_h5)
        matrix_library_ids = {library_id for library_id, _ in library_map.values()}
        csv_library_ids = {row[cr_constants.AGG_ID_FIELD].encode() for row in outs.sample_defs}
        if matrix_library_ids != csv_library_ids:
            str_csv_library_ids = ",".join([x.decode() for x in sorted(csv_library_ids)])
            str_matrix_library_ids = ",".join([x.decode() for x in sorted(matrix_library_ids)])
            this_msg = f"Library IDs specified in CSV ({str_csv_library_ids}) do not match those contained in the input matrix ({str_matrix_library_ids})"
            martian.exit(this_msg)
    copy_csv(args.aggregation_csv, outs.aggregation_csv)


def _replace_with_underscore(ex: UnicodeDecodeError):
    return ("_", ex.end)


codecs.register_error("underscore_replace", _replace_with_underscore)


_ILLEGAL_CHARACTER_TRANSLATION = bytes(
    i if i < 128 and i not in b"\n\r,\\/" else ord("_") for i in range(256)
)


def _text_cleanup(text: str) -> str:
    """Clean non-ascii characters with valid ones.

    Args:
        text (str): text to be corrected.

    Returns:
        str:
    """
    # replacing all non-ascii characters
    return (
        text.encode("ascii", "underscore_replace")
        .translate(_ILLEGAL_CHARACTER_TRANSLATION)
        .decode("ascii")
    )


def copy_csv(in_csv, out_csv):
    """Copy CSV, fixing formatting as needed."""
    with open(in_csv) as reader, open(out_csv, "w") as writer:
        writer.write(reader.read().replace("\r", "\n"))


def check_multiple_paths(
    paths: list[str],
    piperoot: str,
    aggr_csv_path: str,
    field: str,
    subfolder: str = "",
) -> str:
    """Find and return valid paths.

    Customers provide paths in a csv.
    The path provided might be an absolute or a relative path.

    First with absolute path, then piperoot and finally aggregation csv file relative path.
    Exists if path not found.

    Args:
        paths (list[str]): List of paths
        subfolder (str): main subfolder
        piperoot (str): Path to the pipestance root
        aggr_csv_path (str): name of the file
        field (str): which field of the aggr.csv it is

    Returns:
        str: path that has been validated
    """
    # Check if any of the existing path passes tests.
    # Files can have multiple defaults paths to be backwards compatible
    # First we check all of them

    acceptable_paths = []
    if not isinstance(paths, list):
        raise TypeError(f"Paths input must be a list of strings, not {type(paths)}")
    for file_path in paths:
        test_path = Path(subfolder, file_path)
        # If absolute path, we assume that the user has provided the best possible path.
        if not test_path.is_absolute():
            # Assume it's relative to where the aggr csv is.
            aggr_relpath = Path(aggr_csv_path) / test_path
            if os.path.exists(aggr_relpath):
                acceptable_paths.append(aggr_relpath.resolve())
            elif piperoot is not None:
                # assume path is relative to pipestance root
                pipe_relpath = Path(piperoot) / test_path
                if pipe_relpath.exists():
                    acceptable_paths.append(pipe_relpath.resolve())
        elif test_path.exists():
            acceptable_paths.append(test_path)
    # Now we check if we have exactly one path
    if len(acceptable_paths) == 1:
        final_path = acceptable_paths.pop()
        if not os.access(final_path, os.R_OK):
            martian.exit(f"Input is not readable, please check file permissions: {final_path}")
        return str(final_path)
    elif len(acceptable_paths) == 0:
        joined_paths = " or ".join(paths)
        martian.exit(
            f"Specified {field} file cannot be found at {joined_paths}.\n"
            f"Please check that the file is present and that the given path is correct."
        )
    else:
        martian.exit(
            f"Found multiple valid inputs for {field}. Please keep only one valid input file."
        )


def expand(path):
    """Returns the full path."""
    return os.path.expandvars(os.path.expanduser(path))


def _check_rows(reader_obj: csv.DictReader[str]):
    """Iterates through rows, exiting if any fields are blank.

    Args:
        reader_obj (csv.DictReader): The lines present in the input aggr CSV

    Returns:
        sequence of rows.
    """
    for i, row in enumerate(reader_obj, start=1):
        for header, value in row.items():
            if value is None or value.strip() == "":
                martian.exit(f"One of the lines in your CSV has a missing {header} field")
        yield i, row


def parse_sample_sheet(
    piperoot, aggr_csv_path, field_validation, aggr_files
) -> list[dict[str, str]]:
    """Checks that the csv carries the needed fields and that the file need files are present.

    Builds a sample_def list. Each row defines an input samples with validated inputs.
    """
    if not os.path.exists(aggr_csv_path):
        martian.exit(f"Sample sheet does not exist: {aggr_csv_path}")

    if not os.access(aggr_csv_path, os.R_OK):
        martian.exit(
            f"Sample sheet is not readable, please check file permissions: {aggr_csv_path}"
        )
    sample_defs: list[dict[str, str]] = []
    with open(aggr_csv_path) as f:
        # skip comment lines
        reader = csv.DictReader(row for row in f if not row.startswith("#"))

        # if the header row has a path in it, they forgot the header line
        fields = [field.strip() for field in reader.fieldnames]
        for field in fields:
            if os.path.exists(field):
                martian.exit(
                    f"Your CSV header contains a path: {field}\n  Did you forget to include a header line?"
                )
            if not field or field.isspace():
                martian.exit("Your CSV header has an unnamed column, please name all columns.")

        required_fields = [
            field for field in field_validation.keys() if field_validation[field].required
        ]
        optional_fields = [
            field for field in field_validation.keys() if not field_validation[field].required
        ]

        found_fields = [field.lower() for field in fields]
        for field in required_fields:
            if field not in found_fields:
                martian.exit(f"Your header row is missing a required field: {field}.")

        # If we find optional fields, we now consider those as required
        # NOTE: we rely below on an optional field (the 'spatial' folder with tissue positions)
        #       coming after the required fields (specifically the aggr h5)
        for field in optional_fields:
            if field in found_fields:
                required_fields.append(field)
        check_unique = {}
        for field in required_fields:
            check_unique[field] = set()

        # now to check each row
        for i, raw_row in _check_rows(reader):
            # convert field names to lowercase and strip whitespace
            row = {k.strip().lower(): v.strip() for (k, v) in raw_row.items()}

            # check required fields again (it could be missing for just this row)
            for field in required_fields:
                if field not in row:
                    martian.exit("Row %d is missing a required field: %s." % (i, field))
                if len(row[field]) == 0:
                    martian.exit("Row %d has an empty value for field: %s." % (i, field))

                if field_validation[field].type == "string":
                    row[field] = _text_cleanup(row[field])
                # Check path and file existence.
                elif field_validation[field].type == "file_path":
                    # Sending in a list instead of the path directly.
                    file_path = check_multiple_paths(
                        [expand(row[field])],
                        piperoot=piperoot,
                        aggr_csv_path=aggr_csv_path,
                        field=field,
                    )
                    row[field] = file_path
                # Need to add new fields that are not given by user for spatial aggr
                elif field_validation[field].type == "folder_path":
                    # Check that path exists
                    # Sending in a list instead of the path directly.
                    folder_path = check_multiple_paths(
                        [row[field]], piperoot=piperoot, aggr_csv_path=aggr_csv_path, field=field
                    )
                    if folder_path:
                        for field_name in field_validation[field].files:
                            # As we provide default locations for files to find in a folder
                            # we use check_paths instead of check_path
                            row[field_name] = check_multiple_paths(
                                paths=aggr_files[field_name].paths,
                                subfolder=folder_path,
                                piperoot=piperoot,
                                aggr_csv_path=aggr_csv_path,
                                field=field_name,
                            )
                else:
                    martian.exit(f"Field type {field} not recognized.")
                # Add the value for the field in the unique checker
                check_unique[field].add(row[field])
                # If this is the h5 file, check whether this is a Visium HD run.
                if field == cr_constants.AGG_H5_FIELD:
                    mol_info_fn = row[field]
                    if not cr_h5.is_hdf5(mol_info_fn):
                        martian.exit(f"Input molecule file is not a valid HDF5 file: {mol_info_fn}")
                    with cr_mc.MoleculeCounter.open(mol_info_fn, "r") as mol_info:
                        is_visium_hd = mol_info.get_visium_hd_slide_name() is not None
                    if is_visium_hd:
                        martian.exit(
                            f"Visium HD runs, such as in row {i+1} of your CSV, are not currently supported by spaceranger aggr."
                        )
            if AGG_SAMPLE_ID_FIELD in row:
                row[cr_constants.AGG_ID_FIELD] = row.pop(AGG_SAMPLE_ID_FIELD)

            sample_defs.append(row)
        # file has no sample rows
        if len(sample_defs) == 0:
            martian.exit("CSV file contains no sample data.")
        # Finally check that there are no duplicated entries in the csv
        n_samples = len(sample_defs)
        for field_name, field_set in check_unique.items():
            if len(field_set) != n_samples:
                martian.exit(
                    f"There is a duplicated entry in the column '{field_name}' of {aggr_csv_path}"
                )

    return sample_defs
