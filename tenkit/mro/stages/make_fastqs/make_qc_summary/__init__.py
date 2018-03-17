#!/usr/bin/env python
#
# Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
#
from collections import defaultdict
import martian
import json
import os
import tenkit.fasta as tk_fasta
import tenkit.log_subprocess as tk_proc
import tenkit.qc as tk_qc
import tenkit.preflight as tk_preflight
import tenkit.stats as tk_stats

__MRO__ = """
stage MAKE_QC_SUMMARY(
    in  path     run_path,
    in  path     fastq_path,
    in  path     interop_path,
    in  string   bc_read_type,
    in  int      bc_start_index,
    in  int      bc_length,
    in  string   umi_read_type,
    in  int      umi_start_index,
    in  int      umi_length,
    in  bool     rc_i2_read,
    in  map      file_read_types_map,
    in  string   software_version,
    in  string   bcl2fastq_version,
    in  string   bcl2fastq_args,
    out json     qc_summary,
    out bool     completed,
    src py       "stages/make_fastqs/make_qc_summary",
) split using (
    in  string project,
    in  int    lane,
    in  string subfolder,
    in  string sample,
    in  string input_files,
)
"""
ANALYZE_RUN_PD_SUMMARY_CS_KEY_MAP = {
    'read1_q20_fraction': 'r1_q20_bases_fract',
    'read1_q30_fraction': 'r1_q30_bases_fract',
    'read2_q20_fraction': 'r2_q20_bases_fract',
    'read2_q30_fraction': 'r2_q30_bases_fract',
    'barcode_q20_fraction': 'bc_q20_bases_fract',
    'barcode_q30_fraction': 'bc_q30_bases_fract',
    'sample_index_q20_fraction': 'si_q20_bases_fract',
    'sample_index_q30_fraction': 'si_q30_bases_fract',
    'barcode_corrected_match_ratio': 'bc_on_whitelist',
    'gems': 'gems_detected',
    'total_reads': 'number_reads'
}

GEM_COUNT_THRESHOLD_ESTIMATE = 10

class ReadQCStats(object):
    def __init__(self):
        """
        Keep track of barcode sample QC metrics and infer quality rates.
        """
        self.unique_barcode_counts = defaultdict(lambda: 0)
        self.total_read_count = 0
        self.bc_total_base_count = 0
        self.bc_q30_base_count = 0
        self.exact_match_count = 0
        self.corrected_match_count = 0
        self.r1_q30_base_count = 0
        self.r1_total_base_count = 0
        self.r2_q30_base_count = 0
        self.r2_total_base_count = 0
        self.mean_qual_scores = []
        self.__locked_gem_estimate = None

    @property
    def barcode_q30_base_ratio(self):
        """
        Barcode Q30 ratio between 0 and 1.
        """
        return tk_stats.robust_divide(self.bc_q30_base_count, self.bc_total_base_count)

    @property
    def barcode_exact_match_ratio(self):
        """
        Exact barcode match ratio between 0 and 1.
        """
        return tk_stats.robust_divide(self.exact_match_count, self.total_read_count)

    @property
    def barcode_corrected_match_ratio(self):
        """
        Corrected barcode match ratio between 0 and 1.
        """
        return tk_stats.robust_divide(self.corrected_match_count, self.total_read_count)

    def update_sequence_counts(self, sequence_counts):
        """
        Add
        """
        if self.__locked_gem_estimate is not None:
            raise RuntimeError("GEM estimate has been locked; cannot add more sequences.")
        for sequence, count in sequence_counts.items():
            self.unique_barcode_counts[sequence] += count

    @property
    def read1_q30_base_ratio(self):
        """
        Barcode Q30 ratio between 0 and 1 on (non-barcode) read 1.
        """
        return tk_stats.robust_divide(self.r1_q30_base_count, self.r1_total_base_count)

    @property
    def read2_q30_base_ratio(self):
        """
        Barcode Q30 ratio between 0 and 1 on (non-barcode) read 2.
        """
        return tk_stats.robust_divide(self.r2_q30_base_count, self.r2_total_base_count)

    def add_mean_quality_score(self, mean, total_bases):
        """
        Add an average quality score to the stats
        :param average: The average quality score of the bases
        :param total_bases: The number of bases figured into the average
        """
        self.mean_qual_scores.append((mean, total_bases))

    @property
    def mean_barcode_qscore(self):
        """
        Return the average quality score of all bases tracked by this stats object.
        :return:
        """
        total_bases = float(sum([bases for qual, bases in self.mean_qual_scores]))
        if total_bases == 0:
            return 0
        return sum([qual*(float(bases)/total_bases) for qual, bases in self.mean_qual_scores])

    @property
    def total_gem_estimate(self):
        if self.__locked_gem_estimate is not None:
            return self.__locked_gem_estimate
        return len([val for val in self.unique_barcode_counts.values() if val > GEM_COUNT_THRESHOLD_ESTIMATE])

    def lock_gem_estimate(self):
        """
        Cache the GEM estimate and do not allow any more barcode updates.
        """
        self.__locked_gem_estimate = self.total_gem_estimate
        self.unique_barcode_counts = None


def shim_standard_qc_synonyms(output_dict, replace=False):
    """
    Rename or supplant all the keys generated by the
    get_illumina_sequencing_metrics method from their original names
    (dictated by the ANALYZE_RUN_PD stage) to names compatible with
    the summary output of our pipelines.

    :param output_dict: A dictionary containing metric values
    :param replace: Whether to replace the ANALYZE_RUN_PD names with standard names.  If false,
                    the new names will be added.
    """
    new_dict = dict(output_dict)
    for key, val in new_dict.items():
        if isinstance(val, dict):
            new_dict[key] = shim_standard_qc_synonyms(val, replace=replace)
    for analyze_pd_key, summary_cs_key in ANALYZE_RUN_PD_SUMMARY_CS_KEY_MAP.items():
        if new_dict.get(analyze_pd_key):
            if replace:
                new_dict[summary_cs_key] = new_dict.pop(analyze_pd_key)
            else:
                new_dict[summary_cs_key] = new_dict[analyze_pd_key]
    return new_dict


def split(args):
    chunk_dict = defaultdict(list)

    # find projects inside run dir
    dirnames = next(os.walk(args.fastq_path))[1]

    # this should be fine even if --reports-dir and --stats-dir are elsewhere
    # theoretically, you could name your sample 'Reports' or 'Stats' and
    # move the --reports-dir and --stats-dir somewhere else but at that
    # point you're just being daft
    project_dirs = [dn for dn in dirnames if dn not in ('Reports', 'Stats')]
    if not project_dirs:
        # add root folder
        project_dirs = ['.']

    for project in sorted(project_dirs):
        project_path = os.path.join(args.fastq_path, project)
        # split by detected sample (all types)
        r1_files = tk_fasta.find_input_fastq_files_bcl2fastq_demult(
            project_path, 'R1', None, None)
        r2_files = tk_fasta.find_input_fastq_files_bcl2fastq_demult(
            project_path, 'R2', None, None)
        r3_files = tk_fasta.find_input_fastq_files_bcl2fastq_demult(
            project_path, 'R3', None, None)
        r4_files = tk_fasta.find_input_fastq_files_bcl2fastq_demult(
            project_path, 'R4', None, None)
        i1_files = tk_fasta.find_input_fastq_files_bcl2fastq_demult(
            project_path, 'I1', None, None)
        i2_files = tk_fasta.find_input_fastq_files_bcl2fastq_demult(
            project_path, 'I2', None, None)

        # group files together by like sample
        all_files = r1_files + r2_files + r3_files + r4_files + i1_files + i2_files

        for path in sorted(all_files):
            file_spec = tk_fasta.IlmnFastqFile(path)

            # do not add chunk if Undetermined
            if file_spec.prefix != 'Undetermined':
                chunk_dict[(project, file_spec.lane, file_spec.prefix, file_spec.s)].append(path)

    # if all reads were undetermined or bcl2fastq didn't generate anything
    if len(chunk_dict) == 0:
        martian.exit("No FASTQs matched your sample sheet's lane, sample, and indexes. Please recheck your sample sheet.")

    chunk_defs = [{'input_files': file_list, 'project': tup[0], 'lane': tup[1], 'sample': tup[2], 'subfolder': tup[3]} for tup, file_list in sorted(chunk_dict.items())]
    return { 'chunks': chunk_defs }


def main(args, outs):
    out_base = os.path.dirname(outs.qc_summary)
    whitelist_path = tk_preflight.check_barcode_whitelist(args.barcode_whitelist)
    file_infos = [tk_fasta.IlmnFastqFile(path) for path in args.input_files]

    bc_file_type = args.file_read_types_map[args.bc_read_type]
    barcode_files = [f for f in file_infos if f.read == bc_file_type]
    outs.qc_summary = {
        'barcode': [],
        'read1': [],
        'read2': []
    }
    for idx, bf in enumerate(barcode_files):
        output_json_path = os.path.join(out_base, "output_%d_BC.json" % idx)
        subproc_args = [ 'barcodeqc', bf.filename, output_json_path, "--whitelist", whitelist_path,
                         "--bc-start-index", str(args.bc_start_index), "--bc-length",
                         str(args.bc_length)]
        if args.bc_read_type == "I2" and args.rc_i2_read:
            subproc_args.append("--rc")
        tk_proc.call(subproc_args)
        outs.qc_summary['barcode'].append(output_json_path)

    r1_files = [f for f in file_infos if args.file_read_types_map['R1'] == f.read]
    for idx, r1f in enumerate(r1_files):
        output_json_path = os.path.join(out_base, "output_%d_R1.json" % idx)
        if bc_file_type == 'R1':
            start_index = args.bc_start_index
        else:
            start_index = 0
        subproc_args = [ 'q30count', r1f.filename, output_json_path, '--read-start-index', str(start_index) ]
        tk_proc.call(subproc_args)
        outs.qc_summary['read1'].append(output_json_path)

    r2_files = [f for f in file_infos if args.file_read_types_map['R2'] == f.read]
    for idx, r2f in enumerate(r2_files):
        output_json_path = os.path.join(out_base, "output_%d_R2.json" % idx)
        if bc_file_type == 'R2':
            start_index = args.bc_start_index
        else:
            start_index = 0
        subproc_args = ['q30count', r2f.filename, output_json_path, '--read-start-index', str(start_index)]
        tk_proc.call(subproc_args)
        outs.qc_summary['read2'].append(output_json_path)


def join(args, outs, chunk_args, chunk_outs):
    outs.completed = False
    lane_sample_info = defaultdict(ReadQCStats)
    sample_info = defaultdict(ReadQCStats)

    sample_dict = defaultdict(lambda: defaultdict(list))
    for chunk_arg, chunk_out in zip(chunk_args, chunk_outs):
        sample_dict[chunk_arg.sample][chunk_arg.lane].append((chunk_arg, chunk_out))

    for sample, lane_dict in sorted(sample_dict.items()):
        sample_stats = sample_info[sample]
        for lane, chunks in sorted(lane_dict.items()):
            lane_sample_stats = lane_sample_info[(lane, sample)]
            for chunk_arg, chunk_out in chunks:
                for barcode_file in chunk_out.qc_summary['barcode']:
                    with open(barcode_file, 'r') as infile:
                        stats_dict = json.load(infile)
                        martian.log_info("%s: %d reads, %d exact, %d corrected, %d Q30, %d bases" % (
                            barcode_file,
                            stats_dict['total_read_count'],
                            stats_dict['exact_match_count'],
                            stats_dict['corrected_match_count'],
                            stats_dict['q30_base_count'],
                            stats_dict['total_base_count']
                        ))
                        for stats_info in (lane_sample_stats, sample_stats):
                            stats_info.total_read_count += stats_dict['total_read_count']
                            stats_info.exact_match_count += stats_dict['exact_match_count']
                            stats_info.corrected_match_count += stats_dict['corrected_match_count']
                            stats_info.bc_q30_base_count += stats_dict['q30_base_count']
                            stats_info.bc_total_base_count += stats_dict['total_base_count']
                            stats_info.add_mean_quality_score(stats_dict['mean_base_qscore'], stats_dict['total_base_count'])
                            stats_info.update_sequence_counts(stats_dict['validated_sequence_counts'])

                for r1_file in chunk_out.qc_summary['read1']:
                    with open(r1_file, 'r') as infile:
                        r1_stats_dict = json.load(infile)
                        martian.log_info("%s: %d Q30, %d bases" % (
                            r1_file,
                            r1_stats_dict['q30_base_count'],
                            r1_stats_dict['total_base_count'])
                        )
                        for stats_info in (lane_sample_stats, sample_stats):
                            stats_info.r1_q30_base_count += r1_stats_dict['q30_base_count']
                            stats_info.r1_total_base_count += r1_stats_dict['total_base_count']

                for r2_file in chunk_out.qc_summary['read2']:
                    with open(r2_file, 'r') as infile:
                        r2_stats_dict = json.load(infile)
                        martian.log_info("%s: %d Q30, %d bases" % (
                            r2_file,
                            r2_stats_dict['q30_base_count'],
                            r2_stats_dict['total_base_count'])
                                         )
                        for stats_info in (lane_sample_stats, sample_stats):
                            stats_info.r2_q30_base_count += r2_stats_dict['q30_base_count']
                            stats_info.r2_total_base_count += r2_stats_dict['total_base_count']
            lane_sample_stats.lock_gem_estimate()
            martian.log_info("Sample: %s, Lane: %d, GEM Estimate: %d" % (
                sample, lane, lane_sample_stats.total_gem_estimate))
        sample_stats.lock_gem_estimate()
        martian.log_info("Sample: %s, GEM Estimate: %d" % (sample, sample_stats.total_gem_estimate))

    samples_qc = {}
    lane_sample_stats = sorted(lane_sample_info.items())
    for sample, stats in sorted(sample_info.items()):
        sample_qc = {
            'all': {
                'total_reads': stats.total_read_count,
                'barcode_q30_base_ratio': stats.barcode_q30_base_ratio,
                'barcode_exact_match_ratio': stats.barcode_exact_match_ratio,
                'barcode_corrected_match_ratio': stats.barcode_corrected_match_ratio,
                'mean_barcode_qscore': stats.mean_barcode_qscore,
                'read1_q30_base_ratio': stats.read1_q30_base_ratio,
                'read2_q30_base_ratio': stats.read2_q30_base_ratio,
                'gem_count_estimate': stats.total_gem_estimate
            }
        }
        for (lane, lane_sample), lane_stats in lane_sample_stats:
            if lane_sample == sample:
                sample_qc[lane] = {
                    'total_reads': lane_stats.total_read_count,
                    'barcode_q30_base_ratio': lane_stats.barcode_q30_base_ratio,
                    'barcode_exact_match_ratio': lane_stats.barcode_exact_match_ratio,
                    'barcode_corrected_match_ratio': lane_stats.barcode_corrected_match_ratio,
                    'mean_barcode_qscore': lane_stats.mean_barcode_qscore,
                    'read1_q30_base_ratio': lane_stats.read1_q30_base_ratio,
                    'read2_q30_base_ratio': lane_stats.read2_q30_base_ratio,
                    'gem_count_estimate': lane_stats.total_gem_estimate
                }
        samples_qc[sample] = sample_qc

    try:
        # finish with sequencing run interop parsing
        sequencing_stats = tk_qc.get_illumina_sequencing_metrics(args.run_path)

        # add barcode stats
        sequencing_stats = tk_qc.infer_symbolic_read_metrics(
            sequencing_stats, "barcode", args.bc_read_type, args.bc_start_index, args.bc_length)

        # add sample index stats
        sequencing_stats = tk_qc.infer_symbolic_read_metrics(
            sequencing_stats, "sample_index", args.si_read_type, start_idx=0, length=8)

        # add UMI stats if present
        if args.umi_read_type:
            sequencing_stats = tk_qc.infer_symbolic_read_metrics(
                sequencing_stats, "umi", args.umi_read_type, args.umi_start_index, args.umi_length)
    except BaseException:
        sequencing_stats = {}

    sequencing_stats.update({'sample_qc': samples_qc})
    sequencing_stats = shim_standard_qc_synonyms(sequencing_stats, replace=True)
    sequencing_stats['10x_software_version'] = args.software_version
    sequencing_stats['bcl2fastq_version'] = args.bcl2fastq_version
    sequencing_stats['bcl2fastq_args'] = args.bcl2fastq_args

    with open(outs.qc_summary, 'w') as outfile:
        json.dump(sequencing_stats, outfile, indent=4, sort_keys=True)
    outs.completed = True
