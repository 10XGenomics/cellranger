#!/usr/bin/env python
#
# Copyright (c) 2016 10x Genomics, Inc. All rights reserved.
#
# Utilities for reading and generating Illumina Experiment Manager sample sheets
#


import csv
import logging
import re
from collections import namedtuple

from six import ensure_text

import tenkit.pandas as pandas
import tenkit.seq as tk_seq
from tenkit.sample_index import SAMPLE_DUAL_INDEX_MAP, SAMPLE_SINGLE_INDEX_MAP

NEW_SECTION_RE = re.compile(r"^\[(?P<section>\w+)\]")
CELL_RANGE_RE = re.compile(r"^(\d+)\-(\d+)$")
OLIGO_RE = re.compile(r"^[ACGNT]+$")
SECTION_NAME_HEADER = "Header"
SECTION_NAME_READS = "Reads"
SECTION_NAME_SETTINGS = "Settings"
SECTION_NAME_DATA = "Data"

I7_INDEX_COL = "index"
I5_INDEX_COL = "index2"


class IndexAmbiguityException(Exception):
    """Special exception for when a samplesheet is entered where.

    sample indices (or sets of sample indices) may collide if
    the barcode mismatch parameter is too lenient.
    """

    def __init__(self, message):
        super().__init__()
        self.message = message

    def __str__(self):
        return repr(self.message)


class SingleIndexFlowcellException(Exception):
    """Special exception for when a 10x dual-index oligo is.

    specified on a run that appears to not have been done
    with dual-indexing.
    """

    def __init__(self, message):
        super().__init__()
        self.message = message

    def __str__(self):
        return repr(self.message)


class DualIndexFlowcellException(Exception):
    """Special exception for when a 10x single-index oligo is.

    specified without an i5 on a dual-index flowcell.
    """

    def __init__(self, message):
        super().__init__()
        self.message = message

    def __str__(self):
        return repr(self.message)


class NoFastqDataException(Exception):
    """Exception thrown when there is no FASTQ data for the specified params."""


def index_hamming(x, y):
    """Determine hamming distance between two sample indices."""
    if len(x) != len(y):
        # if SI lengths are different, for some reason, return min(len(x,y))
        return min(len(x), len(y))
    # if Ns in sample index, count as distance 0
    return sum(xi not in (yi, "N") for xi, yi in zip(x, y))


class SampleSheetSection:
    """A named section with rows of text."""

    def __init__(self, name, rows):
        self.name = name
        self.rows = [[ensure_text(c) for c in row] for row in rows] if rows else []

    def to_row_array(self):
        if self.name:
            output_rows = [[f"[{self.name}]"]]
        else:
            output_rows = []
        output_rows.extend(self.rows)
        return output_rows

    @property
    def ends_in_blank(self):
        return self.rows and not any(val for val in self.rows[-1])


def read_csv_rows(path):
    """Extract the rows from the CSV at the specified path.

    Will throw an error if the file doesn't exist.

    :type path: string
    :rtype: list[list[string]]
    """
    with open(path, encoding="utf-8-sig") as infile:
        reader = csv.reader(infile, delimiter=",")
        rows = [[col.strip() for col in row] for row in reader]
        # eliminate trailing cols that have no entries (CSI-215)
        for idx, row in enumerate(rows):
            clip_index = 0
            for col in row[::-1]:
                if not col:
                    clip_index -= 1
                else:
                    break
            if clip_index < 0:
                rows[idx] = row[:clip_index]
    return rows


def write_csv_rows(rows, path):
    """Write CSV rows in a standard format.

    :type rows: list[list[string]]
    :type path: string
    """
    with open(path, "w") as outfile:
        writer = csv.writer(outfile, delimiter=",", quoting=csv.QUOTE_MINIMAL)
        writer.writerows(rows)


def rows_are_valid_csv(rows):
    """Determine whether the rows comprise a readable simple CSV,.

    with a lane number, sample and index (in that order)
    :type rows: list[list[string]]
    :rtype: bool
    """
    if not rows:
        return False

    # TENKIT-109 change: having incorrect headers
    # is yielding weird error messages; just return
    # false if the first line isn't a header
    if row_is_simple_header(rows[0]):
        data_idx = 1
    else:
        logging.warning(
            "The first row in the CSV file doesn't contain the required fields.\nPlease check for extra spaces or typos."
        )
        return False

    pop_rows = [row for row in rows[data_idx:] if row]
    tuples = [row_is_simple_data(row) for row in pop_rows]
    for tup in tuples:
        if tup[1]:
            logging.warning(tup[1])

    return all(tup[0] for tup in tuples)


def row_get_col_index(row, col):
    """Return the case-insensitive index of the column.

    in the row, or -1 if the column is not present.
    """
    cols = [c.lower() for c in row]
    if col.lower() in cols:
        return cols.index(col)
    return -1


def row_is_simple_header(row):
    """Determine whether the row is a header row.

     The
    row must start with "lane", "sample" and "index".

    :type row: list[string]
    :rtype: bool
    """
    return (
        len(row) >= 3
        and row_get_col_index(row, "lane") == 0
        and row_get_col_index(row, "sample") == 1
        and row_get_col_index(row, "index") == 2
    )


def row_is_simple_dual_index(row):
    """Akin to row_is_simple_header, but check for dual-index."""
    return row_is_simple_header(row) and row_get_col_index(row, "index2") == 3


def cell_is_valid_lane(cell):
    if cell.isdigit() and int(cell) > 0:
        return True
    if cell in ("all", "*"):
        return True
    if CELL_RANGE_RE.match(cell):
        match = CELL_RANGE_RE.match(cell)
        first = int(match.group(1))
        second = int(match.group(2))
        return first < second
    return False


def cell_expand_lanes(cell, fc_lane_count):
    """Given the value of the lanes in the cell and a lane.

    count from the flowcell, generate the set of lanes
    represented by the cell value.  For example, '1' will
    return [1], but 'all' will return [1,2,3,4] if the
    fc_lane_count is 4.

    Do not return any illegal lane combinations; that is,
    if the flowcell has two lanes and the cell is '3', return
    a blank array.

    :rtype: list[int]
    """
    lanes = []
    if not cell_is_valid_lane(cell):
        return []
    elif cell.isdigit():
        cellint = int(cell)
        lanes = [cellint] if cellint <= fc_lane_count else []
    elif cell in ("all", "*", ""):
        lanes = range(1, 1 + fc_lane_count)
    elif CELL_RANGE_RE.match(cell):
        match = CELL_RANGE_RE.match(cell)
        first = int(match.group(1))
        second = int(match.group(2))
        lanes = [lane for lane in range(first, 1 + second) if lane <= fc_lane_count]
    return lanes


def row_is_simple_data(row):
    """Return whether row appears to match lane-sample-index criteria,.

    and why not if there is not a match.
    :type row: list[string]
    :rtype: tuple[bool, string]
    """
    if not len(row) >= 3:
        return False, "Row has less than three columns"
    if not cell_is_valid_lane(row[0]):
        return False, f"First column not a valid lane: {row[0]}"
    if not row[1]:
        return False, "Sample name blank"
    if not (
        row[2] in SAMPLE_SINGLE_INDEX_MAP
        or row[2] in SAMPLE_DUAL_INDEX_MAP
        or OLIGO_RE.match(row[2])
    ):
        return False, f"Unrecognized sample index: {row[2]}"
    if len(row) > 3 and not (row[3] in SAMPLE_DUAL_INDEX_MAP or OLIGO_RE.match(row[3])):
        return False, f"Unrecognized dual index: {row[3]}"
    return True, None


def row_is_dual_index_data(row):
    """Row contains explicit dual index data, be it from 10x or.

    another i5 oligo.
    """
    _, err = row_is_simple_data(row)
    if err is not None:
        return False, err
    return len(row) > 3, None


def row_is_10x_dual_index(row):
    """Row refers to a paired i7/i5 dual index manufactured by.

    10x.  In this case, only one index is necessary, as a
    symbolic index refers to both i7 and i5 indices.
    """
    _, err = row_is_simple_data(row)
    if err is not None:
        return False, err
    return row[2] in SAMPLE_DUAL_INDEX_MAP, None


def file_is_iem_samplesheet(path):
    """Determine whether the specified input file is an Illumina Experiment Manager (IEM).

    sample sheet.

    :type path: string
    :rtype: bool
    """
    return rows_are_iem_samplesheet(read_csv_rows(path))


def file_is_simple_samplesheet(path):
    """Determine whether the specified input file is a simple CSV sample sheet.

    with lanes/samples/indices.

    :type path: string
    :rtype: bool
    """
    return rows_are_valid_csv(read_csv_rows(path))


def file_is_simple_dual_index_samplesheet(path):
    """Determine whether the specified input file contains dual-indexed.

    data, either by virtue of having an index2 column, or by having
    rows that contain 10x paired sample indices.
    """
    if not file_is_simple_samplesheet(path):
        return False

    rows = read_csv_rows(path)
    if row_is_simple_dual_index(rows[0]):
        return True
    for row in rows[1:]:
        if row_is_10x_dual_index(row):
            return True
    return False


def file_get_iem_data_frame(path):
    """Return the IEM samplesheet data as a Pandas DataFrame,.

    to perform better slicing operations.
    """
    rows = read_csv_rows(path)
    if not rows_are_iem_samplesheet(rows):
        raise ValueError(f"Invalid IEM samplesheet format: {path}")
    section_gen = rows_iem_section_generator(rows)
    for section in section_gen:
        if section_is_valid_data(section):
            # TODO this appears to be a problem if you have data columns
            # with trailing all-blank entries (see CSI-215 fix)
            df = pandas.DataFrame(data=section.rows[1:], columns=section.rows[0])
            # skip tailing rows
            return df[df["Sample_ID"].notnull()]
    raise ValueError(f"Invalid IEM samplesheet format, no data found: {path}")


def rows_are_iem_samplesheet(rows):
    """Determine whether the rows comprise an Illumina Experiment Manager (IEM).

    sample sheet by checking for the presence of a [Data] section with
    sample header.

    :type rows: list[list[string]]
    :rtype: bool
    """
    # criteria: has to have [Data] section with recognized sample index.
    section_gen = rows_iem_section_generator(rows)
    for section in section_gen:
        if section_is_valid_data(section):
            if not iem_rows_all_have_sample_id(section.rows):
                logging.warning("Blank Sample_ID entries detected in data section")
                return False
            else:
                return True
    return False


def iem_rows_all_have_sample_id(rows):
    """Return whether all IEM rows have a filled-in Sample_ID field."""
    if not row_is_data_header(rows[0]):
        return False
    sample_id_idx = rows[0].index("Sample_ID")
    for row in rows[1:]:
        if len(row) >= sample_id_idx + 1:
            if not row[sample_id_idx]:
                return False
    return True


def iem_has_dual_index(path):
    rows = read_csv_rows(path)
    section_gen = rows_iem_section_generator(rows)
    for section in section_gen:
        if section_is_valid_data(section):
            header = section.rows[0]
            if "index2" not in header:
                # we don't have an index2
                return False
            col_idx = header.index("index2")
            for row in section.rows[1:]:
                if row and not row[col_idx]:
                    # one bad apple (blank index2 col) spoils the bunch
                    return False
            return True
    return False


def section_is_valid_data(section):
    """Return whether the specified section contains sufficient information.

    to populate bcl2fastq.
    :type rows: SampleSheetSection
    :rtype: bool
    """
    return (
        section.name == SECTION_NAME_DATA
        and len(section.rows) > 1
        and row_is_data_header(section.rows[0])
    )


def section_get_default_header():
    return SampleSheetSection(SECTION_NAME_HEADER, [["EMFileVersion", "4"]])


def row_is_section_header(row):
    """Return whether or not the specified row marks a new section in the.

    sample sheet (e.g., [Header] in first cell)

    :type row: list[string]
    :rtype: bool
    """
    return len(row) > 0 and NEW_SECTION_RE.match(row[0])


def row_get_section_name(row):
    """Return the name of the section contained in the row, if the.

    row is a section header.

    :type row: list[string]
    :rtype: string
    """
    if not row_is_section_header(row):
        return None
    return NEW_SECTION_RE.match(row[0]).group("section")


def row_is_data_header(row):
    """Returns whether or not the row of strings is an Illumina data header line.

    :type row: list[string]
    :rtype: bool
    """
    # http://support.illumina.com/content/dam/illumina-support/documents/documentation/system_documentation/miseq/miseq-sample-sheet-quick-ref-guide-15028392-j.pdf
    # Sample_ID only field required, though 'index' likely present for bcl2fastq
    #
    return "Sample_ID" in row


def iem_row_get_si_index(row):
    if not row_is_data_header(row):
        return -1
    return row_get_col_index(row, "index")


def iem_row_get_si2_index(row):
    if not row_is_data_header(row):
        return -1
    return row_get_col_index(row, "index2")


def row_get_si_indices(row):
    """Returns the column indices that contain fields (index, index2) which may contain 10x sample oligo sequences.

    :type row: list[string]
    :rtype: tuple[int]
    """
    if not row_is_data_header(row):
        return []
    si_indices = []

    keywords = ("index", "index2")
    for keyword in keywords:
        if keyword in row:
            si_indices.append(row.index(keyword))
        elif keyword.capitalize() in row:
            si_indices.append(row.index(keyword.capitalize()))
    return tuple(si_indices)


def rows_iem_section_generator(rows):
    """Yields groups of rows corresponding to each section of.

    an Illumina sample sheet.

    The format will be::

        {
            'section': (name of section)
            'rows': (rows except for section)
        }

    :type rows: list[list[string]]
    :rtype: generator(SampleSheetSection)
    """
    header = None
    section_rows = []
    for row in rows:
        if row_is_section_header(row):
            if header is not None:
                yield SampleSheetSection(row_get_section_name(header), section_rows)
            elif section_rows:
                # lines before header sections
                yield SampleSheetSection(None, section_rows)
            section_rows = []
            header = row
        else:
            section_rows.append(row)
    # ending section
    if header:
        yield SampleSheetSection(row_get_section_name(header), section_rows)


def get_reads_section(read_length_r1, read_length_r2):
    """Yield a Reads sample sheet section with the specified R1/R2 length.

    :rtype: SampleSheetSection
    """
    rows = [[str(read_length_r1)], [str(read_length_r2)]]
    return SampleSheetSection(SECTION_NAME_READS, rows)


def check_sample_index_collision(data_section, barcode_mismatches_param="1"):
    """Check to see if any of the sample indices would yield.

    a bcl2fastq barcode mismatch error.  This would happen if
    two sample indices on the same lane had a hamming distance of
    than barcode_mismatches_param * 2 or less.  Return a tuple (collision, err)
    that indicates whether a collision is going to occur, and if so, what
    oligos collide.

    Example (mismatches = 1)
    SI1: ATACCAGA
    SI2: ATAGCGGA (hamming distance = 2)

    A read with SI 'ATAGCAGA' has a distance 1 to both SIs, so it's ambiguous which
    sample this read would map to.

    Example (mismatches = 2)
    SI1: ATACCAGA
    SI2: AGATCTGT (hamming distance = 4)

    A read with SI 'AGACCTGA' has distance 2 to both SIs, so it's ambiguous which
    sample this read would map to.

    If the mismatch param is zero, as long there are not two records with
    the same oligo sequences on the same lane.

    For more information about the possibility of barcode mismatch permutations, see
    https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl2fastq/bcl2fastq2-v2-20-software-guide-15051736-03.pdf
    """
    if not section_is_valid_data(data_section):
        return (False, None)

    headers = data_section.rows[0]
    lane_index = -1
    if "Lane" in headers:
        lane_index = headers.index("Lane")
    if "index" in headers:
        si_index = headers.index("index")
    else:
        # what even is this if there is no SI
        return (False, None)

    si2_index = -1
    if "index2" in headers:
        si2_index = headers.index("index2")

    IndexSpec = namedtuple("IndexSpec", "lane, index, index2")
    index_specs = []
    for sample_row in data_section.rows[1:]:
        # throw out incomplete rows
        if (
            len(sample_row) <= lane_index
            or len(sample_row) <= si_index
            or (si2_index != -1 and len(sample_row) <= si2_index)
        ):
            continue
        if not sample_row:
            continue
        lane_spec = ""  # wildcard
        if lane_index != -1:
            lane_spec = sample_row[lane_index].strip()
        index_spec = sample_row[si_index].strip()
        index2_spec = ""
        if si2_index != -1:
            index2_spec = sample_row[si2_index].strip()
        index_specs.append(IndexSpec(lane_spec, index_spec, index2_spec))

    # barcode mismatches can be a combination if dual-indexed.  Also default
    # a param of "1" if not supplied.
    if barcode_mismatches_param == "":
        barcode_mismatches_param = "1"
    barcode_mismatches = [int(mismatch) for mismatch in barcode_mismatches_param.split(",")]
    if si2_index != -1 and len(barcode_mismatches) == 1:
        barcode_mismatches.append(barcode_mismatches[0])

    # now do pairwise comparison of rows to see if there are collisions
    for idx, spec1 in enumerate(index_specs):
        for spec2 in index_specs[idx + 1 :]:
            # if not on same lane (and a lane is specified; a blank lane is an inferred wildcard),
            # ignore
            if spec1.lane not in ("", spec2.lane):
                continue
            # check SI distance
            si_distance = index_hamming(spec1.index, spec2.index)
            if si_distance <= barcode_mismatches[0] * 2:
                # if dual index doesn't exist, we have a collision
                if si2_index == -1:
                    return (True, f"{spec1.index}, {spec2.index}")
                else:
                    # dual index exists; but both I7/I5 have to collide to
                    # force an ambiguity
                    si2_distance = index_hamming(spec1.index2, spec2.index2)
                    if si2_distance <= barcode_mismatches[1] * 2:
                        return (True, f"{spec1.index2}, {spec2.index2} dual index")

    return (False, None)


def get_simple_line_header(dual_index=False, specify_project=False):
    """Return the IEM header line corresponding to simple (lane/sample/index) data.

    :param specify_project: Whether to add the project column.
    :rtype: list[string]
    """
    cols = ["Lane", "Sample_ID", "Sample_Name", "index"]
    if dual_index:
        cols.append("index2")
    if specify_project:
        cols.append("Sample_Project")
    return cols


def transform_simple_line(row, fc_lane_count, project_name=None, dual_indexed=False):
    """Transform a simple (lane,sample,index) line to an IEM-compatible set of lines.

    :param row:  The input row
    :rtype: list[list[string]]
    """
    valid_row, _ = row_is_simple_data(row)
    if not valid_row:
        raise ValueError("Invalid data row: {}".format(",".join(row)))
    dual_index_row, _ = row_is_dual_index_data(row)
    dual_index_10x_si, _ = row_is_10x_dual_index(row)

    lanes = cell_expand_lanes(row[0], fc_lane_count)
    if not lanes:
        raise NoFastqDataException(
            f"Invalid data row: {row}. The flow cell has {fc_lane_count} lanes but lane {row[0]} was specified."
        )
    rows = []
    for lane in lanes:
        # blank sample name to be demuxed later
        cols = [str(lane), row[1], ""]  # lane, sample_name, sample_id
        cols.append(row[2])
        if dual_indexed:
            if dual_index_row:
                cols.append(row[3])
            elif dual_index_10x_si:
                # repeat SI in index2
                cols.append(row[2])
            # if dual_indexed by virtual of dual-index FC, leave blank
            else:
                cols.append("")
        if project_name:
            cols.append(project_name)
        rows.append(cols)
    return rows


def transform_reads_section(section, r1_read_length=None, r2_read_length=None):
    """Transform the existing Reads sample sheet section if R1/R2 are specified.

    :type section: SampleSheetSection
    :type r1_read_length: int or str
    :type r2_read_length: int or str
    :rtype: SampleSheetSection
    """
    if not section.name == SECTION_NAME_READS:
        return section

    rows = [row for row in section.rows]
    if r1_read_length is not None:
        if len(rows) > 0:
            if len(rows[0]) > 0:
                rows[0][0] = str(r1_read_length)
            else:
                rows[0].append(str(r1_read_length))
        else:
            rows.append([str(r1_read_length)])

    if r2_read_length is not None:
        if len(rows) > 1:
            if len(rows[1]) > 0:
                if len(rows[1]) > 0:
                    rows[1][0] = str(r2_read_length)
                else:
                    rows[1].append(str(r2_read_length))
        else:
            if len(rows) == 0:
                raise ValueError("Cannot have R2 override with blank R1")
            rows.append([str(r2_read_length)])

    return SampleSheetSection(SECTION_NAME_READS, rows)


def _overwrite_cell(row, idx, val, fill=""):
    """Do one of two things.

    Either:

    -- overwrite the cell at row[idx] if it exists
    -- pad the row until you can append the idx col, and then write it
    """
    while len(row) <= idx:
        row.append(fill)
    row[idx] = val
    return row


def transform_data_section(
    section,
    i7_index_length=0,
    i5_index_length=0,
    rc_sample_index2=False,
    must_single_index_flowcell=False,
    filter_dual_index=False,
    filter_single_index=False,
):
    """Find 10x sample indices within an existing sample sheet data section.

    Take an existing sample sheet data section, find 10x sample indices
    within them, and then expand the rows to change the index arguments to
    oligos.  If a Sample_Name does not exist, move the original Sample_ID into
    that slot.

    :type section: SampleSheetSection
    :param rc_sample_index2: Whether to reverse-complement the sample index2 (i5) oligo sequence.
    :param must_single_index_flowcell: Whether the original flowcell was dual indexed.
    :param filter_dual_index: Whether to emit dual-index oligos only in a dual-index config.
    :param filter_single_index: Whether to filter single-index oligos only in a single-index config.
    :rtype: SampleSheetSection
    """
    if not section_is_valid_data(section):
        # bail if the rows aren't
        logging.warning("Non-IEM data section passed to expand_section_rows_by_sample_index")
        return section

    header_row = section.rows[0]

    sample_id_idx = header_row.index("Sample_ID")
    sample_name_col_orig = False
    if "Sample_Name" in header_row:
        sample_name_idx = header_row.index("Sample_Name")
        sample_name_col_orig = True
    else:
        header_row = [col for col in header_row]
        header_row.append("Sample_Name")
        sample_name_idx = len(header_row) - 1

    si_index = iem_row_get_si_index(header_row)
    si2_index = iem_row_get_si2_index(header_row)

    # TENKIT-115 broader case (index2 found)
    if must_single_index_flowcell and si2_index != -1:
        raise SingleIndexFlowcellException(
            "Dual-index configuration found for single-index flowcell"
        )

    # TENKIT-119 if only processing single indices, splice si2
    si2_dropped = si2_index
    if filter_single_index and si2_index != -1:
        header_row = header_row[:si2_index] + header_row[si2_index + 1 :]
        si2_index = -1

    # label original sample index
    original_sample_id_idx = len(header_row)
    header_row.append("Original_Sample_ID")

    output_rows = [header_row]
    if len(section.rows) == 1:
        return SampleSheetSection(SECTION_NAME_DATA, output_rows)

    for input_row in section.rows[1:]:
        output_row = []
        # if short row, just proceed
        if len(input_row) <= si_index:
            output_rows.append(input_row)
            continue
        # found symbolic 10x dual index; replace oligo
        if input_row[si_index] in SAMPLE_DUAL_INDEX_MAP:
            if filter_single_index:
                continue
            if must_single_index_flowcell:
                raise SingleIndexFlowcellException("Dual-index oligo for single-indexed flowcell")
            i7, i5 = SAMPLE_DUAL_INDEX_MAP[input_row[si_index]]
            if rc_sample_index2:
                i5 = tk_seq.get_rev_comp(i5)
            # TENKIT-114: truncate to detected index length
            if i7_index_length and len(i7) > i7_index_length:
                i7 = i7[:i7_index_length]
            if i5_index_length and len(i5) > i5_index_length:
                i5 = i5[:i5_index_length]
            output_row = [cell for cell in input_row]
            sample_id = input_row[sample_id_idx]
            if not output_row[sample_name_idx]:
                output_row[sample_name_idx] = sample_id
            if si2_index != -1:
                _overwrite_cell(output_row, si2_index, i5)
            _overwrite_cell(output_row, si_index, i7)
            _overwrite_cell(output_row, original_sample_id_idx, input_row[sample_id_idx])
            output_rows.append(output_row)
        # found symbolic 4-oligo index; expand into multiple rows
        elif input_row[si_index] in SAMPLE_SINGLE_INDEX_MAP:
            if filter_dual_index and si2_index != -1:
                continue
            oligos = SAMPLE_SINGLE_INDEX_MAP[input_row[si_index]]
            for idx, oligo in enumerate(oligos):
                i7 = oligo
                i5 = None if si2_index == -1 else input_row[si2_index]
                if i5:
                    if rc_sample_index2:
                        i5 = tk_seq.get_rev_comp(i5)
                    if i5_index_length and len(i5) > i5_index_length:
                        i5 = i5[:i5_index_length]
                # TENKIT-119 if no dual index oligo specified alongside i7, fail
                elif i5_index_length and not filter_single_index:
                    raise DualIndexFlowcellException(
                        "Single-index 10x sample index used in dual-index flowcell"
                    )

                if i7_index_length and len(i7) > i7_index_length:
                    i7 = i7[:i7_index_length]
                output_row = [cell for cell in input_row]
                if not sample_name_col_orig:
                    _overwrite_cell(output_row, sample_name_idx, "")
                sample_id = input_row[sample_id_idx]
                output_row[sample_id_idx] = "%s_%d" % (sample_id, idx + 1)

                # if Sample_Name blank, populate with sample_id without oligo seq
                if not output_row[sample_name_idx]:
                    output_row[sample_name_idx] = sample_id
                output_row[si_index] = i7
                if si2_index != -1 and i5:
                    _overwrite_cell(output_row, si2_index, i5)

                # we dropped 'index2' from the header
                if si2_dropped != -1 and si2_index == -1:
                    output_row.pop(si2_dropped)

                # last column -- original sample id
                _overwrite_cell(output_row, original_sample_id_idx, input_row[sample_id_idx])
                output_rows.append(output_row)

        else:
            # row doesn't contain symbolic index in i7 -- just output original
            # Assume user knows what they are doing if running oligos
            # directly, and do not try to correct.

            if len(input_row) > sample_id_idx and input_row[sample_id_idx]:
                _overwrite_cell(input_row, original_sample_id_idx, input_row[sample_id_idx])
            output_rows.append(input_row)

    return SampleSheetSection(SECTION_NAME_DATA, output_rows)


def transform_samplesheet_sections(
    sections,
    r1_read_length=None,
    r2_read_length=None,
    i7_index_length=None,
    i5_index_length=None,
    rc_sample_index2=False,
    must_single_index_flowcell=False,
    filter_dual_index=False,
    filter_single_index=False,
):
    """Take a collection of SampleSheetSections and generate the set of SampleSheetSections.

    :type rows: iterable[SampleSheetSection]
    :rtype: list[SampleSheetSection]
    """
    section_dict = {section.name: section for section in sections}
    out_sections = []

    if section_dict.get(None):
        out_sections.append(section_dict[None])

    if section_dict.get(SECTION_NAME_HEADER):
        out_sections.append(section_dict[SECTION_NAME_HEADER])
    else:
        out_sections.append(section_get_default_header())

    if section_dict.get(SECTION_NAME_READS):
        out_sections.append(
            transform_reads_section(
                section_dict[SECTION_NAME_READS],
                r1_read_length=r1_read_length,
                r2_read_length=r2_read_length,
            )
        )
    elif r1_read_length is not None and r2_read_length is not None:
        out_sections.append(get_reads_section(r1_read_length, r2_read_length))

    if section_dict.get(SECTION_NAME_SETTINGS):
        out_sections.append(section_dict[SECTION_NAME_SETTINGS])

    if section_dict.get(SECTION_NAME_DATA):
        out_sections.append(
            transform_data_section(
                section_dict[SECTION_NAME_DATA],
                i7_index_length=i7_index_length,
                i5_index_length=i5_index_length,
                rc_sample_index2=rc_sample_index2,
                must_single_index_flowcell=must_single_index_flowcell,
                filter_dual_index=filter_dual_index,
                filter_single_index=filter_single_index,
            )
        )

    return out_sections


def generate_sections_from_simple_csv(
    lines,
    fc_lane_count,
    r1_read_length=None,
    r2_read_length=None,
    project_name=None,
    dual_indexed_flowcell=False,
):
    """Take lines from a simple CSV lane-sample-index layout and generate a IEM samplesheet.

    :type lines: list[list[string]]
    :type r1_read_length: int
    :type r2_read_length: int
    :type dual_indexed_flowcell: bool
    :rtype: list[SampleSheetSection]
    """
    sections = []
    sections.append(section_get_default_header())

    if r1_read_length and r2_read_length:
        sections.append(get_reads_section(r1_read_length, r2_read_length))

    dual_indexed_sheet = False
    if len(lines) > 0:
        # case 1: generate i7/i5 header if CSV contains dual-index data
        dual_indexed_sheet, _ = row_is_dual_index_data(lines[0])
    dual_indexed = dual_indexed_sheet or dual_indexed_flowcell

    section_rows = [get_simple_line_header(dual_indexed, project_name)]
    for line in lines:
        section_rows.extend(
            transform_simple_line(
                line, fc_lane_count, project_name=project_name, dual_indexed=dual_indexed
            )
        )
    sections.append(SampleSheetSection(SECTION_NAME_DATA, section_rows))
    return sections


def transform_samplesheet(
    csv_path,
    output_path,
    flowcell_lane_count=8,
    r1_read_length=None,
    r2_read_length=None,
    i7_index_length=None,
    i5_index_length=None,
    rc_sample_index2=False,
    project_name=None,
    barcode_mismatches_param="1",
    dual_indexed_flowcell=False,
    must_single_index_flowcell=False,
    filter_dual_index=False,
    filter_single_index=False,
):
    """Take an input CSV file and generate a sample sheet.

    Generates an Illumina Experiment Manager/bcl2fastq
    compatible sample sheet with 10x sample set indices replaced by their
    component oligos, and read counts embedded in the sample sheet if specified.

    FOR NOW:
    If an IEM sample sheet is supplied, the sample index, barcode length and
    project_name will not overwrite what is on the sample sheet.

    :param csv_path: The path to the CSV file to read.
    :param flowcell_lane_count: The number of lanes in the flowcell (to handle simple layout 'all' and '*' directives)
    :param r1_read_length:
    :param r2_read_length:
    :param i7_index_length:
    :param i5_index_length:
    :param rc_sample_index2: Whether the oligo sequence of the sample index2 needs to be reverse complemented.
    :param project_name: The name of the project.
    :param dual_indexed_flowcell: Whether the samplesheet had an i5 index.
    :param must_single_index_flowcell: Throw errors if 10x dual index oligos or a dual-indexed layout was detected.
    :param filter_dual_index: Whether to generate a samplesheet with dual-index sequences only.
    :param filter_single_index: Whether to generate a samplesheet with single-index sequences only.
    :return: Information about the spreadsheet in a dictionary:
       - dual_indexed: Whether the spreadsheet is dual-indexed.
    """
    csv_rows = read_csv_rows(csv_path)
    if rows_are_iem_samplesheet(csv_rows):
        sections = rows_iem_section_generator(csv_rows)
    else:
        populated_rows = [row for row in csv_rows if any(row)]
        if not rows_are_valid_csv(populated_rows):
            raise ValueError(f"Cannot figure out input type: {csv_path}")

        if row_is_simple_header(populated_rows[0]):
            populated_rows = populated_rows[1:]
        sections = generate_sections_from_simple_csv(
            populated_rows,
            flowcell_lane_count,
            r1_read_length=r1_read_length,
            r2_read_length=r2_read_length,
            project_name=project_name,
            dual_indexed_flowcell=dual_indexed_flowcell,
        )

    out_sections = transform_samplesheet_sections(
        sections,
        r1_read_length=r1_read_length,
        r2_read_length=r2_read_length,
        i7_index_length=i7_index_length,
        i5_index_length=i5_index_length,
        rc_sample_index2=rc_sample_index2,
        must_single_index_flowcell=must_single_index_flowcell,
        filter_dual_index=filter_dual_index,
        filter_single_index=filter_single_index,
    )

    output_rows = []
    for idx, section in enumerate(out_sections):
        output_rows.extend(section.to_row_array())
        if idx < len(out_sections) - 1 and not section.ends_in_blank:
            output_rows.append([])

    write_csv_rows(output_rows, output_path)

    output_info = {"dual_indexed": False}
    data_section = [section for section in out_sections if section_is_valid_data(section)]
    if len(data_section) > 0:
        headers = data_section[0].rows[0]
        output_info["dual_indexed"] = "index" in headers and "index2" in headers

        # TENKIT-106 check for index problem
        (collides, message) = check_sample_index_collision(
            data_section[0], barcode_mismatches_param
        )
        if collides:
            raise IndexAmbiguityException(message)

    return output_info
