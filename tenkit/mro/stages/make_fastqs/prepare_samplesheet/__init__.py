#!/usr/bin/env python
#
# Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
#
import csv
import os
import shutil
import tenkit.preflight as tk_preflight
import tenkit.bcl as tk_bcl
import tenkit.lane as tk_lane
import tenkit.samplesheet as tk_sheet

__MRO__ = """
stage PREPARE_SAMPLESHEET(
    in  path    run_path,
    in  map[]   specs,
    in  string  project,
    in  string  bc_read_type,
    in  int     bc_length,
    in  string  si_read_type,
    out csv     samplesheet,
    out csv     input_samplesheet,
    out bool    dual_indexed_samplesheet,
    src py      "stages/make_fastqs/prepare_samplesheet",
)
"""


def make_csv_from_specs(specs, run_info_xml, output_dir):
    """
    Given a set of lanes-sample-indices specs, generate
    a CSV to convert into a simple layout for consumption
    by the IEM samplesheet generation code.

    :param specs: The lane-sample-index tuples (lane optional)
    :param run_info_xml: The path to the sequencer's RunInfo.xml file.  Used to determine default lanes if missing.
    :param output_dir: Where to output the simple layout CSV.
    """
    out_rows = [['lane','sample','index']]

    for spec in specs:
        lanes = spec.get('lanes', range(1, 1+tk_lane.get_flowcell_lane_count(run_info_xml)))
        sample = spec['sample']
        indices = spec['indices']

        for lane in lanes:
            for index in indices:
                out_rows.append([lane,sample,index])

    csv_path = os.path.join(output_dir, 'layout.csv')
    with open(csv_path, 'w') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerows(out_rows)
    return csv_path


def main(args, outs):
    specs = args.specs
    runinfo_path = tk_preflight.check_runinfo_xml(args.run_path)

    output_dir = os.path.dirname(outs.samplesheet)
    csv_specs = [spec for spec in specs if spec.get('csv')]
    if not csv_specs:
        csv_path = make_csv_from_specs(specs, runinfo_path, output_dir)
        outs.input_samplesheet = None
    else:
        csv_path = csv_specs[0]['csv']
        shutil.copy(csv_path, outs.input_samplesheet)

    read_info, flowcell = tk_bcl.load_run_info(runinfo_path)
    (rta_version, rc_i2_read, bcl_params) = tk_bcl.get_rta_version(args.run_path)
    read_info_by_read_type = {r['read_name']:r for r in read_info}
    r1_length = read_info_by_read_type.get('R1', {'read_length': 0})['read_length']
    r2_length = read_info_by_read_type.get('R2', {'read_length': 0})['read_length']

    rc_sample_index = (args.si_read_type == 'I2' and rc_i2_read)
    lane_count = tk_lane.get_flowcell_lane_count(runinfo_path)

    output_info = tk_sheet.transform_samplesheet(csv_path, outs.samplesheet,
        flowcell_lane_count=lane_count,
        r1_read_length=r1_length, r2_read_length=r2_length,
        rc_sample_index=rc_sample_index,
        project_name=args.project)

    outs.dual_indexed_samplesheet = output_info['dual_indexed']
