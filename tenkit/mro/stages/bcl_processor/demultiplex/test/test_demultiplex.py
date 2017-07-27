#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#
# Test the demultiplexer

import subprocess
import os
import os.path
import glob
import tenkit.test as tk_test
from tenkit.constants import TEST_FILE_IN_DIR, TEST_FILE_OUT_DIR
from .. import *

import martian

def fastq_headers(fn):
    return [ read.header for read in FastqParser(fn).read_fastq() ]

def pick_files(pathglob):
    pathglob = os.path.join(TEST_FILE_IN_DIR, pathglob)
    return glob.glob(pathglob)

def get_headers(pathglob):
    r1_files = pick_files(pathglob)
    headers = [ x for f in r1_files for x in fastq_headers(f) ]
    headers.sort()
    return headers

def infile(path):
    return os.path.join(TEST_FILE_IN_DIR, path)

def outfile(path):
    return os.path.join(TEST_FILE_OUT_DIR, path)

in_dir = infile("demultiplex/demultiplex/Project_A/L1")
in_dir_gz = infile("demultiplex/demultiplex_gz/Project_A/L2")
in_dir_no_gz = infile("demultiplex/no_sample_index_gz/Project_A/L2")

mrs_job_name = "demultiplex_test_job"
base_out_dir = outfile("")
out_dir = os.path.join(base_out_dir, mrs_job_name,  "DEMULTIPLEX", "fork0", "files", "demultiplexed_fastq_path")
out_summary = os.path.join(base_out_dir, mrs_job_name, "DEMULTIPLEX", "fork0", "files", "demultiplex_summary.json")

class TestFunctions(tk_test.UnitTestBase):

    def setUp(self):
        martian.test_initialize(outfile(""))

        self.in_files = pick_files("demultiplex/demultiplex/*.fastq.gz")

        # Currently the demultiplexer always interleaves
        base_args = {'sample_index_error_rate': 0.15, 'interleave': 'true',
                     'rc_i2_read': 'false', 'si_read_type': 'I1'}

        self.args = base_args.copy()
        self.args['raw_fastq_path'] = infile("demultiplex/demultiplex")

        self.args_gz = base_args.copy()
        self.args_gz['raw_fastq_path'] = infile("demultiplex/demultiplex_gz")


    def write_mro(self, args):
        tpl = """
        @include "_bcl_processor_stages.mro"
        call DEMULTIPLEX(
            interleave = %(interleave)s,
            rc_i2_read = %(rc_i2_read)s,
            si_read_type = "%(si_read_type)s",
            raw_fastq_path = "%(raw_fastq_path)s",
            sample_index_error_rate = %(sample_index_error_rate)f,
        )
        """
        fn = os.path.join(base_out_dir, "test.mro")
        with open(fn, 'w') as f:
            f.write(tpl % args)

        return fn


    def run_demultiplex(self, args):
        mro_file = self.write_mro(args)
        proc = subprocess.Popen(['mrs', mro_file, mrs_job_name], cwd=base_out_dir, stdout=subprocess.PIPE)
        out, err = proc.communicate()

        if proc.returncode != 0:
            print out
            raise Exception("mrs failed on during test")

        out_path = os.path.join(out_dir, mrs_job_name, "DEMULTIPLEX", "fork0", "files", "demultiplexed_fastq_path")
        return out_path


    def get_summary(self):
        f = file(out_summary, 'r')
        summary = json.load(f)
        return summary

    def test_basic(self):
        self.run_demultiplex(self.args)

        start_headers_r1 = get_headers(in_dir + "/*_R[12]_*")
        demux_headers_r1 = get_headers(out_dir + "/read-RA_*")

        self.assertEqual(start_headers_r1, demux_headers_r1)

        summary = self.get_summary()
        self.assertTrue(summary['invalid_count'] > 0)
        self.assertTrue(all([ x > 0 for x in summary['sample_index_counts'].values()]))
        self.assertEqual(summary['num_reads'], len(start_headers_r1))


    def test_basic_gz(self):
        self.run_demultiplex(self.args_gz)

        start_headers_r1 = get_headers(in_dir_gz + "/*_R[12]_*")
        demux_headers_r1 = get_headers(out_dir + "/read-RA_*")

        self.assertEqual(start_headers_r1, demux_headers_r1)

        summary = self.get_summary()
        self.assertTrue(summary['invalid_count'] > 0)
        self.assertTrue(all([ x > 0 for x in summary['sample_index_counts'].values()]))
        self.assertEqual(summary['num_reads'], len(start_headers_r1))


    def test_interleaved(self):
        args = self.args.copy()
        args["interleave"] = 'true'
        self.run_demultiplex(args)

        start_headers_r1 = get_headers(in_dir + "/*_R[12]_*")
        demux_headers_r1 = get_headers(out_dir + "/read-RA_*")

        self.assertEqual(start_headers_r1, demux_headers_r1)

        summary = self.get_summary()
        self.assertTrue(summary['invalid_count'] > 0)
        self.assertTrue(all([ x > 0 for x in summary['sample_index_counts'].values()]))
        self.assertEqual(summary['num_reads'], len(start_headers_r1))

    def test_interleaved_gz_no_demult(self):
        args = self.args_gz.copy()
        args["interleave"] = 'true'
        args["raw_fastq_path"] = infile("demultiplex/no_sample_index_gz")
        self.run_demultiplex(args)

        start_headers_r1 = get_headers(in_dir + "/*_R[12]_*")
        demux_headers_r1 = get_headers(out_dir + "/read-RA_*")

        self.assertEqual(start_headers_r1, demux_headers_r1)

        summary = self.get_summary()
        self.assertTrue(summary['invalid_count'] > 0)
        self.assertTrue(all([ x > 0 for x in summary['sample_index_counts'].values()]))
        self.assertEqual(summary['num_reads'], len(start_headers_r1))
