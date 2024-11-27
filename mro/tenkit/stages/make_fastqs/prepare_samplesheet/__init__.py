#!/usr/bin/env python
#
# Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
#

from __future__ import annotations

import csv
import os
import shutil

import martian

import tenkit.bcl as tk_bcl
import tenkit.lane as tk_lane
import tenkit.preflight as tk_preflight
import tenkit.samplesheet as tk_sheet

__MRO__ = """
stage PREPARE_SAMPLESHEET(
    in  path    run_path,
    in  map[]   specs,
    in  string  project,
    in  string  bc_read_type,
    in  int     bc_length,
    in  string  si_read_type,
    in  string  bcl2fastq2_args,
    in  bool    force_single_index,
    in  bool    filter_dual_index,
    in  bool    filter_single_index,
    in  bool    rc_i2_override,
    out csv     samplesheet,
    out csv     input_samplesheet,
    out bool    dual_indexed_flowcell,
    src py      "stages/make_fastqs/prepare_samplesheet",
)
"""


def get_barcode_mismatch_param(args):
    """Return the bcl2fastq2 --barcode-mismatches parameter."""
    bcl2fastq2_args = args.bcl2fastq2_args
    arglist = bcl2fastq2_args.strip().split()
    for idx, arg in enumerate(arglist):
        if arg == "--barcode-mismatches":
            return arglist[idx + 1]
    # if not present, return default value
    return "1"


def make_csv_from_specs(specs, run_info_xml, output_dir):
    """Given a set of lanes-sample-indices specs, generate.

    a CSV to convert into a simple layout for consumption
    by the IEM samplesheet generation code.

    :param specs: The lane-sample-index tuples (lane optional)
    :param run_info_xml: The path to the sequencer's RunInfo.xml file.  Used to determine default lanes if missing.
    :param output_dir: Where to output the simple layout CSV.
    """
    out_rows = [["lane", "sample", "index"]]

    for spec in specs:
        lanes = spec.get("lanes", range(1, 1 + tk_lane.get_flowcell_lane_count(run_info_xml)))
        sample = spec["sample"]
        indices = spec["indices"]

        for lane in lanes:
            for index in indices:
                out_rows.append([lane, sample, index])

    csv_path = os.path.join(output_dir, "layout.csv")
    with open(csv_path, "w") as csv_file:
        writer = csv.writer(csv_file)
        writer.writerows(out_rows)
    return csv_path


def determine_rc_i2(run_path, require_i2_detection):
    rc_i2 = tk_bcl.RTAVersionInformation(run_path).rc_i2

    if rc_i2 is None and require_i2_detection:
        # Can't determine the RC status of I2, likely due to a missing Recipe.
        # Error out here -- user will need to specify the I2 status explicitly.
        # NovaSeq I2 could be in either orientation so we don't want to proceeed with
        # without know
        msg = """
mkfastq could not determine the sequencing workflow used on this NovaSeq
flowcell. Detection of sequencing mode relies on recipe xml file in the
`Recipe/` subdirectory of the bcl run folder.
Specify `--rc-i2-override=true` for NovaSeq reagent kits v1.5  and greater, or
`--rc-i2-override=false` earlier NovaSeq reagent kits.
"""
        martian.exit(msg)

    if rc_i2 is None:
        # We can't determine the rc_i2 status- possibly because there is no i2,
        # and the caller doesn't want us to fail.
        rc_i2 = False

    return rc_i2


def main(args, outs):
    specs = args.specs
    runinfo_path = tk_preflight.check_runinfo_xml(args.run_path)

    output_dir = os.path.dirname(outs.samplesheet)
    csv_specs = [spec for spec in specs if spec.get("csv")]
    if not csv_specs:
        csv_path = make_csv_from_specs(specs, runinfo_path, output_dir)
        outs.input_samplesheet = None
    else:
        csv_path = csv_specs[0]["csv"]
        shutil.copy(csv_path, outs.input_samplesheet)

    read_info = tk_bcl.load_run_info(runinfo_path)[0]
    read_info_by_read_type = {r["read_name"]: r for r in read_info}
    r1_length = read_info_by_read_type.get("R1", {"read_length": 0})["read_length"]
    r2_length = read_info_by_read_type.get("R2", {"read_length": 0})["read_length"]
    i7_length = read_info_by_read_type.get("I1", {"read_length": 0})["read_length"]
    i5_length = 0
    has_i2 = False
    i2_info = read_info_by_read_type.get("I2")
    if i2_info is not None:
        i5_length = i2_info["read_length"]
        has_i2 = True

    if args.rc_i2_override is None:
        # This call fail the pipeline if we're on NovaSeq and we have i2 but
        # can't determine the i2 orientation.
        rc_i2_read = determine_rc_i2(args.run_path, has_i2)
    else:
        rc_i2_read = args.rc_i2_override

    martian.log_info(f"I2 Read is RC: {rc_i2_read!s}")

    lane_count = tk_lane.get_flowcell_lane_count(runinfo_path)

    barcode_mismatch_param = get_barcode_mismatch_param(args)
    try:
        tk_sheet.transform_samplesheet(
            csv_path,
            outs.samplesheet,
            flowcell_lane_count=lane_count,
            r1_read_length=r1_length,
            r2_read_length=r2_length,
            i7_index_length=i7_length,
            i5_index_length=i5_length,
            rc_sample_index2=rc_i2_read,
            project_name=args.project,
            barcode_mismatches_param=barcode_mismatch_param,
            dual_indexed_flowcell=has_i2,
            # if there is no i5 read and force-single-index is False,
            # then tell the method to fail if a dual-index config
            # is specified; noop if filter_single_index specified
            must_single_index_flowcell=(
                not args.filter_single_index and (not has_i2 and not args.force_single_index)
            ),
            filter_dual_index=args.filter_dual_index,
            filter_single_index=args.filter_single_index,
        )
    except tk_sheet.IndexAmbiguityException as e:
        martian.exit(
            f"""
The sample sheet supplied has a sample index collision.  This can happen if the
same sample index and lane were specified for multiple samples, or in certain
cases where 10x Chromium i7 Multiplex Kit and i7 Multiplex Kit N samples were
run on the same flowcell.  It is a known issue that certain sample index
combinations from these two different kits only differ by two bases-- meaning
that it is possible for the sequencer to generate index reads with sequences
that are one base away from multiple sample indices.

Please check your samplesheet to verify you do not have duplicate lane-sample
index pairs for multiple samples.  If there are no duplicates, please run
mkfastq with a --barcode-mismatches=0 argument.  (The default parameter is
--barcode-mismatches=1).  This will make bcl2fastq only accept reads that
match the sample indices exactly. The small percentage of reads that are a
single base away from multiple sample indices will be ignored.

Colliding sample index oligos: {e.message}
"""
        )
    except tk_sheet.SingleIndexFlowcellException:
        martian.exit(
            """
The samplesheet contains a sample with a dual-index (i7/i5) sample index,
but the flowcell was only run in a single-index (i7) configuration.  If you
wish to proceed, please run mkfastq with one of the following arguments:
    --force-single-index:  Use the i7 read only to demultiplex the samples,
                           even if the sample indices include an i5 index.
                           Suggested if a dual-index (e.g., SI-TT-A6) sample
                           index was mistakenly run on an i7-only flowcell.
    --filter-single-index: Only generate a samplesheet with samples identified
                           by an i7-only sample index, ignoring dual-indexed
                           samples.  Dual-indexed samples will not be
                           demultiplexed.
"""
        )
    except tk_sheet.DualIndexFlowcellException:
        martian.exit(
            """
The samplesheet contains a sample with a single-index (i7 only) sample index,
but the flowcell was run in a dual-index (i7/i5) configuration.  It is not
possible to demultiplex dual-indexed 10x samples and single-indexed 10x samples
within the same run of bcl2fastq. If you wish to proceed, please run mkfastq
with one of the following arguments:
    --filter-dual-index: Only generate a samplesheet with the samples identified
                         by i7/i5 dual-indices (e.g., SI-TT-A6), ignoring single-
                         index samples.  Single-index samples will not be
                         demultiplexed.
    --filter-single-index: Only generate a samplesheet with samples identified
                           by an i7-only sample index, ignoring dual-indexed
                           samples.  Dual-indexed samples will not be
                           demultiplexed.
NOTE: cellranger mkfastq no longer supports demultiplexing Single Cell 3' v1
flowcells. Please use cellranger 3.1 or earlier for that purpose.
"""
        )
    except tk_sheet.NoFastqDataException as e:
        martian.exit(f"No FASTQ data is available for the specified parameters:\n{e}")

    outs.dual_indexed_flowcell = has_i2
