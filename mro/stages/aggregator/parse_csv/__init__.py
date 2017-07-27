#!/usr/bin/env python
#
# Copyright (c) 2017 10x Genomics, Inc. All rights reserved.
#

import cellranger.constants as cr_constants
import cellranger.matrix as cr_matrix
import csv
import os
import martian

__MRO__  = """
stage PARSE_CSV(
    in  path  pipestance_root,
    in  csv   aggregation_csv,
    in  bool  reanalyze,
    in  h5    matrix_h5,
    out csv   aggregation_csv,
    out map[] sample_defs,
    src py    "stages/aggregator/parse_csv",
)
"""

def main(args, outs):
    if args.reanalyze and not args.aggregation_csv:
        return # CSV not required for reanalyze
    check_paths = not args.reanalyze
    outs.sample_defs = parse_sample_sheet(args.pipestance_root, args.aggregation_csv, check_paths=check_paths)
    if args.reanalyze and args.matrix_h5:
        library_map = cr_matrix.get_gem_group_index(args.matrix_h5)
        matrix_library_ids = set([library_id for library_id, gem_group in library_map.values()])
        csv_library_ids = set([row[cr_constants.AGG_ID_FIELD] for row in outs.sample_defs])
        if matrix_library_ids != csv_library_ids:
            martian.exit("Library IDs specified in CSV (%s) do not match those contained in the input matrix (%s)"
             % (csv_library_ids, matrix_library_ids))
    copy_csv(args.aggregation_csv, outs.aggregation_csv)

def copy_csv(in_csv, out_csv):
    ''' Copy CSV, fixing formatting as needed.'''
    with open(in_csv, 'r') as reader, open(out_csv, 'w') as writer:
        writer.write(reader.read().replace('\r', '\n'))

def parse_sample_sheet(piperoot, filename, check_paths=True):
    if not os.path.exists(filename):
        martian.exit("Sample sheet does not exist: %s" % filename)

    if not os.access(filename, os.R_OK):
        martian.exit("Sample sheet is not readable, please check file permissions: %s" % filename)

    sample_defs = []
    with open(filename, 'rU') as f:
        # skip comment lines
        ff = filter(lambda row: not row.startswith('#') , f)
        reader = csv.DictReader(ff)

        # if the header row has a path in it, they forgot the header line
        fields = [field.strip() for field in reader.fieldnames]
        for field in fields:
            if os.path.exists(field):
                martian.exit("Your CSV header contains a path: %s\n  Did you forget to include a header line?" % field)

        # check required fields
        required_fields = [cr_constants.AGG_ID_FIELD, cr_constants.AGG_H5_FIELD]
        found_fields = [field.lower() for field in fields]
        for field in required_fields:
            if field not in found_fields:
                martian.exit("Your header row is missing a required field: %s." % field)

        for (i, raw_row) in enumerate(reader, start=1):
            # convert field names to lowercase and strip whitespace
            row = {k.strip().lower(): v.strip() for (k,v) in raw_row.iteritems()}

            # check required fields again (it could be missing for just this row)
            for field in required_fields:
                if field not in row:
                    martian.exit("Row %d is missing a required field: %s." % (i, field))
                if len(row[field]) == 0:
                    martian.exit("Row %d has an empty value for field: %s." % (i, field))

            if check_paths:
                # get H5 path and check that it exists
                input_h5 = expand(row[cr_constants.AGG_H5_FIELD])

                if not os.path.isabs(input_h5):
                    # assume path is relative to pipestance root
                    input_h5 = os.path.abspath(expand(os.path.join(piperoot, input_h5)))

                if not os.path.exists(input_h5):
                    martian.exit("Specified %s file does not exist: %s" % (cr_constants.AGG_H5_FIELD, input_h5))

            else:
                input_h5 = None

            row[cr_constants.AGG_H5_FIELD] = input_h5

            # actually add sample def
            sample_defs.append(row)

        # file has no sample rows
        if len(sample_defs) == 0:
            martian.exit("CSV file contains no sample data.")

    return sample_defs

def expand(p):
    return os.path.expandvars(os.path.expanduser(p))
