#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# 1. Demultiplex reads using the I1 reads, if present. Initially we will detect
#    common sample indicies by looking at the reads. In  the future we will
#    accept a sample sheet
# 2. Put the output FASTQ files in a canonical location
#
# FASTQ file demultiplexer/interleaver for Illumina files.
# Takes a list of fastq files with filenames of the form:
# <prefix>_S0_L001_R1_001.fastq
# where S denotes something, L denotes lane,
# R1/R2/I1/I2 denotes Read 1, Read 2, Index Read 1, Index Read 2.
# The final 001 denotes the file number and is used to split large outputs into
# multiple files.
# If you supply multiple files that differ only by their file number, they will
# be demultiplexed
# in order and the sequences concatenated, dropping the file number index.
#
# All input fastq files must have the same <prefix> string.
#
# The tool will read an index file to determine which are the 'common' barcodes.
# Reads matching the common barcodes will be put into files labelled with the
# barcode sequence. The remaining reads will be put labelled with barcode 'X'.
#
import os
import itertools
import json
import subprocess
import numpy
import glob
import gzip
import collections
import martian
import tenkit.cache as tk_cache
import tenkit.dict_utils as tk_dict
import tenkit.seq as tk_seq
import tenkit.stats as tk_stats
from tenkit.fasta import IlmnFastqFile

__MRO__ = """
stage DEMULTIPLEX(
    in  path   raw_fastq_path,
    in  float  sample_index_error_rate,
    in  bool   interleave,
    in  bool   rc_i2_read,
    out path   demultiplexed_fastq_path,
    out json   demultiplex_summary,
    src py     "stages/bcl_processor/demultiplex",
) split using (
    in  bool   demultiplex,
    in  string common_bcs,
    in  string input_files,
    in  string read_types,
    in  int    chunk_number,
)
"""

MAX_INDICES = 1000
INVALID_SAMPLE_INDEX = 'X'

def join(args, outs, chunk_defs, chunk_outs):
    os.mkdir(outs.demultiplexed_fastq_path)

    # Move output file to final location
    for chunk_out in chunk_outs:
        for f in os.listdir(chunk_out.demultiplexed_fastq_path):
            in_file = os.path.join(chunk_out.demultiplexed_fastq_path, f)
            subprocess.call(['mv', in_file, outs.demultiplexed_fastq_path])

    # Combine result data
    r = {'num_reads':0, 'num_clusters': 0, 'invalid_count':0, 'sample_index_counts':{}}
    for chunk_out in chunk_outs:
        # We count each end of a paired-end read separately in the summary file.
        summary_counts = json.load(open(chunk_out.demultiplex_summary))
        num_clusters = sum(summary_counts.values())
        num_reads = 2 * num_clusters
        invalid_reads = summary_counts[INVALID_SAMPLE_INDEX]
        del summary_counts[INVALID_SAMPLE_INDEX]
        summary_counts = {k:2*v for (k,v) in summary_counts.iteritems()}
        r['num_clusters'] += num_clusters
        r['num_reads'] += num_reads
        r['invalid_count'] += invalid_reads
        r['sample_index_counts'] = tk_dict.add_dicts(r['sample_index_counts'], summary_counts, depth=1)
    r['invalid_frac'] = tk_stats.robust_divide(r['invalid_count'], r['num_clusters'])

    json.dump(r, open(outs.demultiplex_summary, 'w'))


def main(args, outs):
    if args.demultiplex:
        main_demultiplex_go(args, outs)
    else:
        main_demultiplex(args, outs)

class FastqRow:
    def __init__(self, header, seq, qual):
        self.header = header
        self.seq = seq
        self.qual = qual

    def write(self, f):
        f.write(self.header + "\n")
        f.write(self.seq + "\n")
        f.write("+\n")
        f.write(self.qual + "\n")

class FastqParser:
    def __init__(self, infile, rc=False):
        self.file = infile
        self.rc = rc

    def read_fastq(self):
        if self.file[-2:] == "gz":
            proc = subprocess.Popen(["gunzip", "--stdout", self.file], stdout=subprocess.PIPE)
            reader = proc.stdout
        else:
            reader = file(self.file, "r")

        while True:
            header = reader.next().strip()
            seq = reader.next().strip()
            reader.next() # incr line
            qual = reader.next().strip()

            if self.rc:
                seq = tk_seq.get_rev_comp(seq)
                qual = qual[::-1]

            yield FastqRow(header, seq, qual)

        reader.close()

class FindCommonBarcodes:
    def get_index_counts(self, fastqs, sample_size=1e6):
        #sample_per_fastq = sample_size / len(fastqs)
        index_counts = {}

        for fq in fastqs:
            n = 0
            for read in fq.read_fastq():
                index_counts[read.seq] = index_counts.get(read.seq, 0) + 1
                n += 1
                if n > sample_size:
                    break

        return index_counts

    # Look at a bunch of index reads and choose the commonly occuring barcodes.
    # return (common_barcodes, rare_barcodes)
    def pick_common_indexes(self, fastqs):
        index_counts = self.get_index_counts(fastqs)

        items_list = index_counts.items()
        items_list.sort(cmp=None, key=lambda x: x[1], reverse=True)
        total_counts = sum(v for (k,v) in items_list)

        c = 0
        i = 0
        for i in range(len(index_counts)):
            c += items_list[i][1]

            if c > 0.75 * total_counts:
                break

        # number of barcodes that account for 75% of reads
        c75 = i

        # median # of observations of barcodes accounting for the 75%
        num_obs_good_bcs = numpy.median([ count for (bc, count) in items_list[:(c75+1)] ])
        martian.log_info("Median counts of good barcodes in 1e6 reads: %s" % num_obs_good_bcs)

        min_obs_bc = max(num_obs_good_bcs / 200, 25)

        # only demultiplex a reasonable number of sample indices
        if len(items_list) > MAX_INDICES:
            min_obs_bc = max(min_obs_bc, items_list[MAX_INDICES][1])

        good_bcs = [ k for (k,v) in items_list if v > min_obs_bc ]
        noise_bcs = [ k for (k,v) in items_list if v <= min_obs_bc ]

        return (good_bcs, noise_bcs)


# Demultiplex a series of FASTQ iterators.
# The index iterator must be the first iterator
# The filename Vector{String} must all be in the same order of read type.
# Interleave map tells which output to write each of the seq_interator entries to.
def process_fastq_chunk(seq_iters, filenames, no_match_filenames, file_cache,
    _interleave_map, summary_counts, max_reads = -1):

    #out_streams = { k:[ gzip.open(x, open_file_mode) for x in v ] for (k,v) in filenames.items() }
    #no_match_out_streams = [ gzip.open(x, open_file_mode) for x in no_match_filenames ]
    valid_bcs = set(filenames.keys())

    if _interleave_map is None:
        interleave_map = range(len(seq_iters))
    else:
        interleave_map = _interleave_map

    read_iterators = itertools.izip(*seq_iters)
    n = 0

    for read_set in read_iterators:
        # Log the counts for each sample index
        bc_seq = read_set[0].seq
        if bc_seq in valid_bcs:
            summary_counts[bc_seq] += 1
        else:
            summary_counts[INVALID_SAMPLE_INDEX] += 1

        #target_streams = out_streams.get(bc_seq, no_match_out_streams)
        tfn = filenames.get(bc_seq, no_match_filenames)
        target_streams = [file_cache.get(x) for x in tfn]


        for i in range(len(read_set)):
            target_index = interleave_map[i]
            read_set[i].write(target_streams[target_index])

        n += 1
        if (n%10**5) == 0:
            martian.log_info("Reads processed %i" % n)

        if max_reads > 0 and n >= max_reads:
            break

# Demultiplex a series of FASTQ iterators.
# The index iterator must be the first iterator
# The filename Vector{String} must all be in the same order of read type.
# Interleave map tells which output to write each of the seq_interator entries to.
def process_fastq_chunk_no_demult(seq_iters, filenames, file_cache,
    _interleave_map, summary_counts, max_reads = -1):

    if _interleave_map is None:
        interleave_map = range(len(seq_iters))
    else:
        interleave_map = _interleave_map

    read_iterators = itertools.izip(*seq_iters)
    n = 0

    for read_set in read_iterators:
        # Log the counts for each sample index
        summary_counts[INVALID_SAMPLE_INDEX] += 1

        target_streams = [file_cache.get(x) for x in filenames]

        for i in range(len(read_set)):
            target_index = interleave_map[i]
            read_set[i].write(target_streams[target_index])

        n += 1
        if (n%10**5) == 0:
            martian.log_info("Reads processed %i" % n)

        if max_reads > 0 and n >= max_reads:
            break




def groupby(f, items):
    groups = collections.defaultdict(list)
    for i in items:
        groups[f(i)].append(i)
    return groups


def split(args):
    # Code supports non-interleaved mode, but we're not currently passing that argument
    #do_interleave = True

    file_glob = os.path.join(args.raw_fastq_path, "Project_*", "*", "*.fastq*")
    print file_glob
    files = glob.glob(file_glob)

    if len(files) == 0:
        martian.throw("No FASTQ files were found for this run. Perhaps there was an error in bcl2fastq, or the input data is bad?")

    file_info = [ IlmnFastqFile(x) for x in files ]

    # Some check for consistency of inputs
    #prefixes = set([x.prefix for x in file_info])

    # May need to revisit handling of multiple lanes in the future!
    # if not args.collapse_lanes and len(prefixes) > 1:
    #    martian.log_info("Observed multiple prefixes: %s" % prefixes)
    #    return 1

    file_groups = groupby(lambda x: (x.s, x.lane, x.group), file_info).items()

    # Order the demultiplex by the group filename
    file_groups.sort(key = lambda (k,files): files[0].group)

    num_files_per_group = [len(v) for (k,v) in file_groups]

    if len(set(num_files_per_group)) > 1:
        martian.throw("You are missing or have extra fastq file! Check your input files")

    read_sets = [tuple(sorted(f.read for f in grp_files)) for (grp, grp_files) in file_groups]

    if len(set(read_sets)) > 1:
        martian.throw("You don't have the same set of reads for all read groups! Check your input files!")

    # The list of read_types we are getting, eg. ["R1", "I1", "I2", "R2"]
    read_types = read_sets[0]

    index_read = args.si_read_type
    demultiplex = True

    if not (index_read in read_types):
        martian.log_info("Supplied read types: %s" % str(read_types))
        martian.log_info("Copying reads with no demultiplexing")
        demultiplex = False
        good_bcs = []

    else:
        # Set up everything we need for demultiplexing
        sort_read_types = [ index_read ]
        sort_read_types.extend(sorted([x for x in read_types if x != index_read]))
        read_types = sort_read_types

        # Figure out which barcodes are well-represented in the index file
        # We will only demultiplex these ones.
        index_files_for_calibration = [ [f for f in grp if f.read == index_read][0] for (k,grp) in file_groups]

        martian.log_info("Determining common barcodes from %d %s files" % (
            len(index_files_for_calibration), index_read))

        bcFind = FindCommonBarcodes()
        if args.rc_i2_read:
            bcFastqs = [FastqParser(f.filename, rc=(f.read == "I2")) for f in index_files_for_calibration]
        else:
            bcFastqs = [FastqParser(f.filename) for f in index_files_for_calibration]
        (good_bcs, bad_bcs) = bcFind.pick_common_indexes(bcFastqs)

        martian.log_info("Got %i common barcodes" % len(good_bcs))
        martian.log_info("Read types: %s" % str(read_types))


    chunk_defs = []
    chunk_number = 0
    chunk_len = 1

    for chunk_start in range(0, len(file_groups), chunk_len):
        grps = file_groups[chunk_start:(chunk_start+chunk_len)]

        chunk = {'demultiplex': demultiplex, 'common_bcs': good_bcs, 'read_types': read_types, 'chunk_number': chunk_number}

        chunk['input_files'] = [f.filename for (grp, file_list) in grps for f in file_list]
        chunk_defs.append(chunk)
        chunk_number += 1

    return {'chunks': chunk_defs}

def main_demultiplex_go(args, outs):
    data = {
        'common_sample_indices': args.common_bcs,
        'file_groups': [],
    }
    file_info = [IlmnFastqFile(x) for x in args.input_files]
    file_groups = groupby(lambda x: (x.s, x.lane, x.group), file_info).items()
    for (_, lane, _), input_files in file_groups:
        files = {read_type: [f for f in input_files if f.read == read_type][0].filename for read_type in args.read_types}
        data['file_groups'].append({
            'lane': lane,
            'files': files,
        })

    input_json_path = martian.make_path('godemux_input.json')
    with open(input_json_path, 'w') as f:
        json.dump(data, f)

    subproc_args = ['godemux', input_json_path, outs.demultiplexed_fastq_path,
                    outs.demultiplex_summary, '--demult-read', args.si_read_type,
                    '--chunk', str(args.chunk_number)]
    if args.rc_i2_read:
        subproc_args += ['--rci2read']
    subprocess.check_call(subproc_args)

# DEPRECATED
# This code is only here for the case where demultiplex = False
def main_demultiplex(args, outs):

    do_interleave = True
    file_info = [ IlmnFastqFile(x) for x in args.input_files ]
    file_groups = groupby(lambda x: (x.s, x.lane, x.group), file_info).items()

    demultiplex = args.demultiplex
    read_types = args.read_types
    good_bcs = args.common_bcs

    # For no interleaving:
    interleave_map = range(len(args.read_types))
    output_reads = args.read_types

    if not ("R1" in read_types) or not ("R2" in read_types):
        martian.throw("You requested interleaving, but you don't have R1 and R2 read types")

    r1_slot = read_types.index("R1")
    r2_slot = read_types.index("R2")
    interleave_map[r2_slot] = r1_slot
    output_reads = [ read_types[idx] for idx in numpy.unique(interleave_map) ]

    # Create output path
    os.mkdir(outs.demultiplexed_fastq_path)
    output_path = outs.demultiplexed_fastq_path

    # counts of each valid barcode and non-matching barcodes
    summary_counts = { bc:0 for bc in good_bcs }
    summary_counts[INVALID_SAMPLE_INDEX] = 0

    with tk_cache.FileHandleCache(open_func=gzip.open) as file_cache:
        # Iterate over the file groups
        for (k, input_files) in file_groups:
            # original path:
            # <path>/<prefix>_S0_L001_R1_001.fastq
            # new path:
            # <outpath>/read-<read_id>_si-xxxxx_lane-<lane>_chunk-<chunk>.fastq
            # input_files should have constant prefix, S, and L
            # sort input_files to match the read_types
            read_to_file_dict = { x.read:x for x in input_files }
            input_files = [ read_to_file_dict[rt] for rt in read_types ]
            output_files = [ read_to_file_dict[rt] for rt in output_reads ]

            def output_file(path, in_file, barcode):
                if do_interleave and in_file.read[0] == "R":
                    read = "RA"
                else:
                    read = in_file.read

                # Chunk over lanes to get some parallelism to speed up alignment
                f = "read-%s_si-%s_lane-%03d-chunk-%03d.fastq.gz" % (read, barcode, in_file.lane, args.chunk_number)
                return os.path.join(path, f)

            if args.rc_i2_read:
                # For NextSeq we need to RC the I2 read
                input_iters = [ FastqParser(f.filename, rc=(f.read == "I2")).read_fastq() for f in input_files ]
            else:
                input_iters = [ FastqParser(f.filename).read_fastq() for f in input_files ]

            martian.log_info("Demultiplexing from: %s" % input_files[0].filename)

            if demultiplex:
                bc_files = { bc: [output_file(output_path, f, bc) for f in output_files] for bc in good_bcs }
                err_files = [ output_file(output_path, f, "X") for f in output_files ]
                process_fastq_chunk(input_iters, bc_files, err_files, file_cache, interleave_map, summary_counts)

            else:
                out_files = [ output_file(output_path, f, 'X') for f in output_files ]
                process_fastq_chunk_no_demult(input_iters, out_files, file_cache, interleave_map, summary_counts)

        output_files = file_cache.have_opened

    # Write out the summary counts to JSON
    with open(outs.demultiplex_summary, "w") as f:
        json.dump(summary_counts, f)
