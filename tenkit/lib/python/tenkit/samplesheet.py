#!/usr/bin/env python
#
# Copyright (c) 2016 10x Genomics, Inc. All rights reserved.
#
# Utilities for reading and generating Illumina Experiment Manager sample sheets
#

import re
import csv
import logging
import pandas as pd

from tenkit.constants import SAMPLE_INDEX_MAP
import tenkit.seq as tk_seq

NEW_SECTION_RE = re.compile(r'^\[(?P<section>\w+)\]')
CELL_RANGE_RE = re.compile(r'^(\d+)\-(\d+)$')
OLIGO_RE = re.compile(r'^[ACGNT]+$')
SECTION_NAME_HEADER = 'Header'
SECTION_NAME_READS = 'Reads'
SECTION_NAME_SETTINGS = 'Settings'
SECTION_NAME_DATA = 'Data'

I7_INDEX_COL = 'index'
I5_INDEX_COL = 'index2'


class SampleSheetSection(object):
    """
    A named section with rows of text.
    """
    def __init__(self, name, rows):
        self.name = name
        self.rows = list(rows) if rows else []

    def to_row_array(self):
        if self.name:
            output_rows = [['[%s]' % self.name]]
        else:
            output_rows = []
        output_rows.extend(self.rows)
        return output_rows

    @property
    def ends_in_blank(self):
        return self.rows and not any([val for val in self.rows[-1]])


def read_csv_rows(path):
    """
    Extract the rows from the CSV at the specified path.
    Will throw an error if the file doesn't exist.

    :type path: string
    :rtype: list[list[string]]
    """
    with open(path, 'rU') as infile:
        reader = csv.reader(infile, delimiter=',')
        rows = [row for row in reader]
        # eliminate trailing cols that have no entries (CSI-215)
        for idx, row in enumerate(rows):
            clipIndex = 0
            for col in row[::-1]:
                if not col:
                    clipIndex -= 1
                else:
                    break
            if clipIndex < 0:
                rows[idx] = rows[idx][:clipIndex]
    return rows


def write_csv_rows(rows, path):
    """
    Write CSV rows in a standard format.
    :type rows: list[list[string]]
    :type path: string
    """
    with open(path, 'w') as outfile:
        writer = csv.writer(outfile, delimiter=',', quoting=csv.QUOTE_MINIMAL)
        writer.writerows(rows)


def rows_are_valid_csv(rows):
    """
    Determine whether the rows comprise a readable simple CSV,
    with a lane number, sample and index (in that order)
    :type rows: list[list[string]]
    :rtype: bool
    """
    if not rows:
        return False

    if row_is_simple_header(rows[0]):
        data_idx = 1
    else:
        data_idx = 0

    pop_rows = [row for row in rows[data_idx:] if row]
    tuples = [row_is_simple_data(row) for row in pop_rows]
    for tup in tuples:
        if tup[1]:
            logging.warning(tup[1])

    return all([tup[0] for tup in tuples])


def row_is_simple_header(row):
    """
    Determine whether the row is a header row.  The
    three cols must be "lane","sample" and "index", in order

    :type row: list[string]
    :rtype: bool
    """
    return len(row) == 3 \
        and row[0].lower() == 'lane' \
        and row[1].lower() == 'sample' \
        and row[2].lower() == 'index'


def cell_is_valid_lane(cell):
    if cell.isdigit() and int(cell) > 0:
        return True
    if cell in ('all', '*'):
        return True
    if CELL_RANGE_RE.match(cell):
        match = CELL_RANGE_RE.match(cell)
        first = int(match.group(1))
        second = int(match.group(2))
        return first < second
    return False


def cell_expand_lanes(cell, fc_lane_count):
    """
    Given the value of the lanes in the cell and a lane
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
    elif cell in ('all', '*'):
        lanes = range(1,1+fc_lane_count)
    elif CELL_RANGE_RE.match(cell):
        match = CELL_RANGE_RE.match(cell)
        first = int(match.group(1))
        second = int(match.group(2))
        lanes = [lane for lane in range(first, 1+second) if lane <= fc_lane_count]
    return lanes


def row_is_simple_data(row):
    """
    Return whether row appears to match lane-sample-index criteria,
    and why not if there is not a match.
    :type row: list[string]
    :rtype: tuple[bool, string]
    """
    if not len(row) >= 3:
        return False, "Row has less than three columns"
    if not cell_is_valid_lane(row[0]):
        return False, "First column not a valid lane: %s" % row[0]
    if not row[1]:
        return False, "Sample name blank"
    if not (row[2] in SAMPLE_INDEX_MAP or OLIGO_RE.match(row[2])):
        return False, "Unrecognized sample index: %s" % row[2]
    return True, None


def file_is_iem_samplesheet(path):
    """
    Determine whether the specified input file is an Illumina Experiment Manager (IEM)
    sample sheet.

    :type path: string
    :rtype: bool
    """
    return rows_are_iem_samplesheet(read_csv_rows(path))


def file_is_simple_samplesheet(path):
    """
    Determine whether the specified input file is a simple CSV sample sheet
    with lanes/samples/indices.

    :type path: string
    :rtype: bool
    """
    return rows_are_valid_csv(read_csv_rows(path))


def file_get_iem_data_frame(path):
    """
    Return the IEM samplesheet data as a Pandas DataFrame,
    to perform better slicing operations.
    """
    rows = read_csv_rows(path)
    if not rows_are_iem_samplesheet(rows):
        raise ValueError("Invalid IEM samplesheet format: %s" % path)
    section_gen = rows_iem_section_generator(rows)
    for section in section_gen:
        if section_is_valid_data(section):
            # TODO this appears to be a problem if you have data columns
            # with trailing all-blank entries (see CSI-215 fix)
            df = pd.DataFrame(data=section.rows[1:], columns=section.rows[0])
            # skip tailing rows
            return df[df['Sample_ID'].notnull()]
    raise ValueError("Invalid IEM samplesheet format, no data found: %s" % path)


def rows_are_iem_samplesheet(rows):
    """
    Determine whether the rows comprise an Illumina Experiment Manager (IEM)
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
    """
    Return whether all IEM rows have a filled-in Sample_ID field.
    """
    if not row_is_data_header(rows[0]):
        return False
    sample_id_idx = rows[0].index('Sample_ID')
    for row in rows[1:]:
        if len(row) >= sample_id_idx+1:
            if not row[sample_id_idx]:
                return False
    return True


def iem_has_dual_index(path):
    rows = read_csv_rows(path)
    section_gen = rows_iem_section_generator(rows)
    for section in section_gen:
        if section_is_valid_data(section):
            header = section.rows[0]
            if 'index2' not in header:
                # we don't have an index2
                return False
            col_idx = header.index('index2')
            for row in section.rows[1:]:
                if row and not row[col_idx]:
                    # one bad apple (blank index2 col) spoils the bunch
                    return False
            return True
    return False


def section_is_valid_data(section):
    """
    Return whether the specified section contains sufficient information
    to populate bcl2fastq.
    :type rows: SampleSheetSection
    :rtype: bool
    """
    return section.name == SECTION_NAME_DATA \
           and len(section.rows) > 1 and row_is_data_header(section.rows[0])


def section_get_default_header():
    return SampleSheetSection(SECTION_NAME_HEADER, [['EMFileVersion', '4']])


def row_is_section_header(row):
    """
    Return whether or not the specified row marks a new section in the
    sample sheet (e.g., [Header] in first cell)

    :type row: list[string]
    :rtype: bool
    """
    return len(row) > 0 and NEW_SECTION_RE.match(row[0])


def row_get_section_name(row):
    """
    Return the name of the section contained in the row, if the
    row is a section header.

    :type row: list[string]
    :rtype: string
    """
    if not row_is_section_header(row):
        return None
    return NEW_SECTION_RE.match(row[0]).group('section')


def row_is_data_header(row):
    """
    Returns whether or not the row of strings is an Illumina data header line.

    :type row: list[string]
    :rtype: bool
    """
    # http://support.illumina.com/content/dam/illumina-support/documents/documentation/system_documentation/miseq/miseq-sample-sheet-quick-ref-guide-15028392-j.pdf
    # Sample_ID only field required, though 'index' likely present for bcl2fastq
    #
    return 'Sample_ID' in row


def row_get_si_indices(row):
    """
    Returns the column indices that contain fields (index, index2) which may contain 10x sample oligo sequences.

    :type row: list[string]
    :rtype: tuple[int]
    """
    if not row_is_data_header(row):
        return []
    si_indices = []

    keywords = ('index','index2')
    for keyword in keywords:
        if keyword in row:
            si_indices.append(row.index(keyword))
        elif keyword.capitalize() in row:
            si_indices.append(row.index(keyword.capitalize()))
    return tuple(si_indices)


def rows_iem_section_generator(rows):
    """
    Yields groups of rows corresponding to each section of
    an Illumina sample sheet.  The format will be
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
    """
    Yield a Reads sample sheet section with the specified R1/R2 length.
    :rtype: SampleSheetSection
    """
    rows = [[str(read_length_r1)], [str(read_length_r2)]]
    return SampleSheetSection(SECTION_NAME_READS, rows)


def get_simple_line_header(specify_project=False):
    """
    Return the IEM header line corresponding to simple (lane/sample/index) data.

    :param specify_project: Whether to add the project column.
    :rtype: list[string]
    """
    cols = ['Lane','Sample_ID','Sample_Name','index']
    if specify_project:
        cols.append('Sample_Project')
    return cols


def transform_simple_line(row, fc_lane_count, project_name=None):
    """
    Transform a simple (lane,sample,index) line to an IEM-compatible set of lines.

    :param row:  The input row
    :rtype: list[list[string]]
    """
    valid_row, err = row_is_simple_data(row)
    if not valid_row:
        raise ValueError("Invalid data row: %s" % ','.join(row))

    lanes = cell_expand_lanes(row[0], fc_lane_count)
    rows = []
    for lane in lanes:
        # blank sample name to be demuxed later
        cols = [lane,row[1],'']
        cols.append(row[2])
        if project_name:
            cols.append(project_name)
        rows.append(cols)
    return rows


def transform_reads_section(section, r1_read_length=None, r2_read_length=None):
    """
    Transform the existing Reads sample sheet section if R1/R2 are specified.

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

def _overwrite_cell(row, idx, val, fill=''):
    """
    Either:
    -- overwrite the cell at row[idx] if it exists
    -- pad the row until you can append the idx col, and then write it
    """
    while len(row) <= idx:
        row.append(fill)
    row[idx] = val
    return row


def transform_data_section(section, rc_sample_index=False):
    """
    Take an existing sample sheet data section, find 10x sample indices
    within them, and then expand the rows to change the index arguments to
    oligos.  If a Sample_Name does not exist, move the original Sample_ID into
    that slot.

    :type section: SampleSheetSection
    :param rc_sample_index: Whether to reverse-complement the sample index oligo sequence.
    :rtype: SampleSheetSection
    """
    if not section_is_valid_data(section):
        # bail if the rows aren't
        logging.warning('Non-IEM data section passed to expand_section_rows_by_sample_index')
        return section

    header_row = section.rows[0]

    sample_id_idx = header_row.index('Sample_ID')
    sample_name_col_orig = False
    if 'Sample_Name' in header_row:
        sample_name_idx = header_row.index('Sample_Name')
        sample_name_col_orig = True
    else:
        header_row = [col for col in header_row]
        header_row.append('Sample_Name')
        sample_name_idx = len(header_row)-1

    # label original sample index
    original_sample_id_idx = len(header_row)
    header_row.append('Original_Sample_ID')

    si_indices = row_get_si_indices(header_row)
    # if index columns not found, search all rows instead
    if not si_indices:
        si_indices = [idx for idx in range(len(header_row))]

    output_rows = [header_row]
    if len(section.rows) == 1:
        return SampleSheetSection(SECTION_NAME_DATA, output_rows)

    for input_row in section.rows[1:]:
        found_10x_si = False
        for col in si_indices:
            if len(input_row) <= col:
                continue
            if input_row[col] in SAMPLE_INDEX_MAP:
                found_10x_si = True
                oligos = SAMPLE_INDEX_MAP[input_row[col]]
                for idx, oligo in enumerate(oligos):
                    # flip the oligo sequence if index is RC (NextSeq I2)
                    if rc_sample_index:
                        oligo = tk_seq.get_rev_comp(oligo)
                    output_row = [cell for cell in input_row]
                    if not sample_name_col_orig:
                        _overwrite_cell(output_row, sample_name_idx, '')
                    sample_id = input_row[sample_id_idx]
                    output_row[sample_id_idx] = "%s_%d" % (sample_id, idx+1)

                    # if Sample_Name blank, populate with sample_id without oligo seq
                    if not output_row[sample_name_idx]:
                        output_row[sample_name_idx] = sample_id
                    output_row[col] = oligo

                    # last column -- original sample id
                    _overwrite_cell(output_row, original_sample_id_idx, input_row[sample_id_idx])
                    output_rows.append(output_row)

                # only replicate on first occurrence of sample index
                break
        if not found_10x_si:
            # row doesn't conform to data -- just output original
            # NOTE: reverse-complementing not applied here

            # append original sample id to end for consistency
            # iff there is something in the sample_id column
            if len(input_row) > sample_id_idx and input_row[sample_id_idx]:
                _overwrite_cell(input_row, original_sample_id_idx, input_row[sample_id_idx])
            output_rows.append(input_row)

    return SampleSheetSection(SECTION_NAME_DATA, output_rows)


def transform_samplesheet_sections(sections, r1_read_length=None, r2_read_length=None, rc_sample_index=False):
    """
    Take a collection of SampleSheetSections and generate the set of SampleSheetSections
    to write out.
    :type rows: iterable[SampleSheetSection]
    :rtype: list[SampleSheetSection]
    """
    section_dict = {section.name:section for section in sections}
    out_sections = []

    if section_dict.get(None):
        out_sections.append(section_dict[None])

    if section_dict.get(SECTION_NAME_HEADER):
        out_sections.append(section_dict[SECTION_NAME_HEADER])
    else:
        out_sections.append(section_get_default_header())

    if section_dict.get(SECTION_NAME_READS):
        out_sections.append(transform_reads_section(
            section_dict[SECTION_NAME_READS],
            r1_read_length=r1_read_length,
            r2_read_length=r2_read_length))
    elif r1_read_length is not None and r2_read_length is not None:
        out_sections.append(get_reads_section(r1_read_length, r2_read_length))

    if section_dict.get(SECTION_NAME_SETTINGS):
        out_sections.append(section_dict[SECTION_NAME_SETTINGS])

    if section_dict.get(SECTION_NAME_DATA):
        out_sections.append(transform_data_section(section_dict[SECTION_NAME_DATA], rc_sample_index=rc_sample_index))

    return out_sections


def generate_sections_from_simple_csv(
        lines,
        fc_lane_count,
        r1_read_length=None,
        r2_read_length=None,
        project_name=None):
    """
    Take lines from a simple CSV lane-sample-index layout and generate a IEM samplesheet.

    :type lines: list[list[string]]
    :type r1_read_length: int
    :type r2_read_length: int
    :rtype: list[SampleSheetSection]
    """
    sections = []
    sections.append(section_get_default_header())

    if r1_read_length and r2_read_length:
        sections.append(get_reads_section(r1_read_length, r2_read_length))

    section_rows = [get_simple_line_header(project_name)]
    for line in lines:
        section_rows.extend(transform_simple_line(
            line,
            fc_lane_count,
            project_name=project_name))
    sections.append(SampleSheetSection(SECTION_NAME_DATA, section_rows))
    return sections


def transform_samplesheet(
        csv_path, output_path,
        flowcell_lane_count=8,
        r1_read_length=None,
        r2_read_length=None,
        rc_sample_index=False,
        project_name=None):
    """
    Take an input CSV file and generate an Illumina Experiment Manager/bcl2fastq
    compatible sample sheet with 10x sample set indices replaced by their
    component oligos, and read counts embedded in the sample sheet if specified.

    FOR NOW:
    If an IEM sample sheet is supplied, the sample index, barcode length and
    project_name will not overwrite what is on the sample sheet.

    :param csv_path: The path to the CSV file to read.
    :param flowcell_lane_count: The number of lanes in the flowcell (to handle simple layout 'all' and '*' directives)
    :param r1_read_length:
    :param r2_read_length:
    :param rc_sample_index: Whether the oligo sequence of the sample index needs to be reverse complemented.
    :param project_name: The name of the project.
    :return: Information about the spreadsheet in a dictionary:
       - dual_indexed: Whether the spreadsheet is dual-indexed.
    """
    csv_rows = read_csv_rows(csv_path)
    if rows_are_iem_samplesheet(csv_rows):
        sections = rows_iem_section_generator(csv_rows)
    else:
        populated_rows = [row for row in csv_rows if any(row)]
        if not rows_are_valid_csv(populated_rows):
            raise ValueError("Cannot figure out input type: %s" % csv_path)

        if row_is_simple_header(populated_rows[0]):
            populated_rows = populated_rows[1:]
        sections = generate_sections_from_simple_csv(
            populated_rows,
            flowcell_lane_count,
            r1_read_length=r1_read_length,
            r2_read_length=r2_read_length,
            project_name=project_name)

    out_sections = transform_samplesheet_sections(
        sections, r1_read_length=r1_read_length, r2_read_length=r2_read_length, rc_sample_index=rc_sample_index)

    output_rows = []
    for idx, section in enumerate(out_sections):
        output_rows.extend(section.to_row_array())
        if idx < len(out_sections)-1 and not section.ends_in_blank:
            output_rows.append([])

    write_csv_rows(output_rows, output_path)

    output_info = {
        'dual_indexed': False
    }
    data_section = [section for section in out_sections if section_is_valid_data(section)]
    if len(data_section) > 0:
        headers = data_section[0].rows[0]
        output_info['dual_indexed'] = 'index' in headers and 'index2' in headers

    return output_info


