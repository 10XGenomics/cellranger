#!/usr/bin/env python
#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#
# Utilities for manipulating fasta and fastq files
#

import glob
import gzip
import numpy
import os
import re
import sys

from tenkit.constants import ILLUMINA_QUAL_OFFSET, BCL_PROCESSOR_FASTQ_MODE, ILMN_BCL2FASTQ_FASTQ_MODE

# Parse the ILMN fastq filename to get the read name, lane, read group,
# and S field. We expect a filename of the form
# <path>/<prefix>_S0_L001_R1_001.fastq
class IlmnFastqFile:
    def __init__(self, fullpath):
        fn = os.path.split(fullpath)[-1]

        dot_parts = fn.split(".")
        if dot_parts[-1] == "fastq":
            name = dot_parts[-2]
        elif len(dot_parts) > 2 and dot_parts[-2] == "fastq":
            name = dot_parts[-3]
        else:
            raise NameError("%s is not a fastq file" % fullpath)

        all_flds = name.split("_")

        flds = all_flds[-4:]
        self.prefix = "_".join(all_flds[:-4])

        self.s = flds[0][1:]
        self.lane = int(flds[1][2:])
        self.read = flds[2]
        self.group = int(flds[3])

        self.filename = fullpath

BCL_PROCESSOR_FILENAME_REGEX = re.compile(r'read-(\w\w)_si-([^_]+)_lane-(\d+)-chunk-(\d+)')

class BclProcessorFastqFile:
    """
    Parse the FASTQ name from the demux stage output, to get the
    read name, lane, sample index, and chunk.
    """
    def __init__(self, fullpath):
        filename = os.path.basename(fullpath)
        dot_parts = filename.split(".")
        if dot_parts[-1] == "fastq":
            name = dot_parts[-2]
        elif len(dot_parts) > 2 and dot_parts[-2] == "fastq":
            name = dot_parts[-3]
        else:
            raise NameError("%s is not a fastq file" % fullpath)

        name_parts = BCL_PROCESSOR_FILENAME_REGEX.match(name)
        if not name_parts:
            raise NameError("Not a demux output fastq: %s" % fullpath)
        self.read = name_parts.group(1)
        self.index = name_parts.group(2)
        self.lane = int(name_parts.group(3), 10)
        self.chunk = int(name_parts.group(4), 10)
        self.filename = fullpath

def write_read_fastq(fastq_file, name, seq, qual):
    """ Writes a single read to a fastq file
    """
    fastq_file.write('@' + name + '\n')
    fastq_file.write(seq + '\n')
    fastq_file.write('+\n')
    fastq_file.write(qual + '\n')

def write_read_fasta(fasta_file, name, seq):
    """ Writes a single read to a fastq file
    """
    fasta_file.write('>' + name + '\n')
    fasta_file.write(seq + '\n')

def write_read_pair_fastq(fastq_file, name1, seq1, qual1, name2, seq2, qual2):
    """ Writes a read-pair to a fastq file for an interleaving format
    """
    write_read_fastq(fastq_file, name1, seq1, qual1)
    write_read_fastq(fastq_file, name2, seq2, qual2)

def read_generator_fastq(fastq_file, paired_end=False):
    """ Returns an interator over a fastq file tha produces (name, seq, qual)
    If paired_end, returns both reads (assuming interleaving fastq)
    """
    line_index = 0
    for line in fastq_file:
        if line_index == 0:
            name1 = line.strip()[1:]
        elif line_index == 1:
            seq1 = line.strip()
        elif line_index == 3:
            qual1 = line.strip()
        elif line_index == 4:
            name2 = line.strip()[1:]
        elif line_index == 5:
            seq2 = line.strip()
        elif line_index == 7:
            qual2 = line.strip()

        line_index += 1
        if not(paired_end) and line_index == 4:
            line_index = 0

        if line_index == 8:
            line_index = 0

        if line_index == 0:
            if paired_end:
                yield (name1, seq1, qual1, name2, seq2, qual2)
            else:
                yield (name1, seq1, qual1)


def interleave_fastq(fastq1, fastq2, out_fastq, remove_Ns=False):
    """ Creates an interleaved fastq file from non-interleaved files
    """
    gen_1 = read_generator_fastq(fastq1)
    gen_2 = read_generator_fastq(fastq2)

    for (name1, seq1, qual1) in gen_1:
        (name2, seq2, qual2) = gen_2.next()

        if remove_Ns:
            if 'N' in seq1 or 'N' in seq2:
                continue

        write_read_pair_fastq(out_fastq, name1, seq1, qual1, name2, seq2, qual2)

def uninterleave_fastq(in_fastq, out_fastq1, out_fastq2):
    """ Takes an interleaved fastq and outputs two fastqs with the reads split
    """
    gen = read_generator_fastq(in_fastq, paired_end=True)
    for (name1, seq1, qual1, name2, seq2, qual2) in gen:
        write_read_fastq(out_fastq1, name1, seq1, qual1)
        write_read_fastq(out_fastq2, name2, seq2, qual2)

def get_qvs(qual):
    if qual is None:
        return None

    return numpy.fromstring(qual, dtype=numpy.byte) - ILLUMINA_QUAL_OFFSET

def get_bases_qual(qual, cutoff):
    if qual is None:
        return None

    qvs = numpy.fromstring(qual, dtype=numpy.byte) - ILLUMINA_QUAL_OFFSET
    return numpy.count_nonzero(qvs[qvs > cutoff])

def get_min_qual(qual):
    if qual is None or len(qual) == 0:
        return None

    return (numpy.fromstring(qual, dtype=numpy.byte) - ILLUMINA_QUAL_OFFSET).min()

def get_mean_qual(qual):
    if qual is None or len(qual) == 0:
        return None

    return (numpy.fromstring(qual, dtype=numpy.byte) - ILLUMINA_QUAL_OFFSET).mean()

def get_expected_errors(qual):
    if qual is None or len(qual) == 0:
        return None

    qvs = numpy.fromstring(qual, dtype=numpy.byte) - ILLUMINA_QUAL_OFFSET
    perr = 10.0**(-qvs/10.0)
    return perr.sum()


def find_input_fastq_files_10x_preprocess(path, read_type, sample_index, lanes, maxNs=2):
    ''' Find fastq files with a matching sample index, obeying the demultiplex filename scheme,
    with a limited number of Ns in the sample index sequence '''

    # In the case where this read wasn't taken, the glob won't match
    # anything, return no files

    # We want to pull in all sample indices that match in all non-N positions, with up to 2 Ns
    # construct the appropriate glob, then filter.
    if sample_index != '*':
        si_glob = ''.join([ "[%cN]" % bse for bse in sample_index])
    else:
        si_glob = "*"
        # Don't worry about Ns if we are just taking everything
        maxNs = 100

    if lanes is None or lanes == []:
        glb = os.path.join(path, "read-%s_si-%s_*.fastq*" % (read_type, si_glob))
        files = glob.glob(glb)
    else:
        files = []
        for lane in lanes:
            glb = os.path.join(path, "read-%s_si-%s_lane-%03d[_\-]*.fastq*" % (read_type, si_glob, int(lane)))
            files.extend(glob.glob(glb))

    good_files = []
    # filter files to remove those with > 2 Ns in the sample index
    for f in files:
        m = re.match(".*si-([A-Z]*)_", f)
        si = m.groups()[0]
        num_Ns = len([x for x in si if x == 'N'])
        if num_Ns <= maxNs:
            good_files.append(f)

    files = sorted(good_files)
    return files


def find_input_fastq_files_bcl2fastq_demult(path, read_type, sample, lanes):
    ''' Find fastq files demultiplex by bcl2fastq 2.17 -- we glob over all subdirectories beneath
        'path', for fastq files prefixed by 'sample', in the given set of 'lanes', or all
        lanes if 'lanes' is None '''

    # In the case where this read wasn't taken, the glob won't match
    # anything, return no files

    # We want to pull in all sample indices that match in all non-N positions, with up to 2 Ns
    # construct the appropriate glob, then filter.
    if sample == None:
        sample = "*"

    if lanes is None or lanes == []:
        file_pattern = "%s_*_L[0-9][0-9][0-9]_%s_[0-9][0-9][0-9].fastq*" % (sample, read_type)
        # sample sheet case (Project/Samples/fastq)
        glb = os.path.join(path, "*/%s" % file_pattern)
        files = glob.glob(glb)
        # direct folder case
        if not files:
            glb = os.path.join(path, file_pattern)
            files = glob.glob(glb)
    else:
        files = []
        for lane in lanes:
            file_pattern = "%s_*_L%03d_%s_[0-9][0-9][0-9].fastq*" % (sample, int(lane), read_type)
            glb = os.path.join(path, "*/%s" % file_pattern)
            lane_files = glob.glob(glb)
            if not lane_files:
                glb = os.path.join(path, file_pattern)
                lane_files = glob.glob(glb)
            files.extend(lane_files)

    files = sorted(files)

    # TENKIT-91 fix
    #
    # glob is limited (e.g., you can't tell it to match a non-zero amount of
    # characters that don't contain an underscore); so use the file name parser functionality
    # to see if you really matched something
    #
    # For example, with the above expression, you would match "target" and "target_not_wanted"
    # even if you only wanted "target".  One could limit the glob match to S0, but if you had a
    # valid sample with a S\d+ suffix, that may interfere as well.
    if sample != "*":
        files = [f for f in files if os.path.basename(IlmnFastqFile(f).prefix) == sample]
    return files


def find_input_file_type_with_samples(path):
    """
    Try to determine the demux type of the files at the specified path, and the unique
    sample names of the FASTQs, if applicable.

    Return either "BCL_PROCESSOR" or "ILMN_BCL2FASTQ" as the first argument.
    If ILMN_BCL2FASTQ is returned, also return an array of detected prefixes;
    the array should be empty if the detected mode is BCL_PROCESSOR mode.

    If files if neither variety are found, return None in the first parameter and
    a blank array as the sample list.

    :rtype: (str, list[str])
    """
    files = find_input_fastq_files_10x_preprocess(path, 'RA', '*', None)
    if files:
        return BCL_PROCESSOR_FASTQ_MODE, []

    files = find_input_fastq_files_bcl2fastq_demult(path, 'R1', None, None)
    if not files:
        return None, []

    ilmn_files = [IlmnFastqFile(f) for f in files]
    samples = set([os.path.basename(f.prefix) for f in ilmn_files])
    return ILMN_BCL2FASTQ_FASTQ_MODE, sorted(samples)


class AmbiguousValueError(ValueError):
    """
    Special value error signaling that the arguments need to be
    more specific.
    """
    pass


def check_fastq_types(path, fastqprefix):
    """
    Validate that the path and fastqprefix arguments are compatible and allowed.  Return
    the input_mode and sample_name that should be sent downstream.

    :param path:  The supplied input path argument.
    :param fastqprefix:  The supplied --fastqprefix sample prefix argument.  Can be a single string or collection
    :type fastqprefix: str || list[str]
    :return: The input_mode and sample_name value that should be passed into the MRO.
    :raises ValueError: If the supplied path + prefix argument are invalid.
    """
    if fastqprefix is None:
        fastqprefixes = None
        fastqprefix_outstr = 'any'
    elif isinstance(fastqprefix, basestring):
        fastqprefixes = [fastqprefix]
        fastqprefix_outstr = fastqprefix
    else:
        fastqprefixes = list(fastqprefix)
        fastqprefix_outstr = ','.join(fastqprefixes)

    demux_type, samples = find_input_file_type_with_samples(path)
    if not demux_type:
        raise ValueError("No input FASTQs were found for the requested parameters.\n\nIf your files came from bcl2fastq or mkfastq:\n - Make sure you are specifying the correct --sample(s), i.e. matching the sample sheet\n - Make sure your files follow the correct naming convention, e.g. SampleName_S1_L001_R1_001.fastq.gz (and the R2 version)\n - Make sure your --fastqs points to the correct location.\n\nRefer to the \"Specifying Input FASTQs\" page at https://support.10xgenomics.com/ for more details.\n\n")
    if demux_type == BCL_PROCESSOR_FASTQ_MODE and fastqprefix:
        raise ValueError("Cannot use --fastqprefix argument with FASTQs generated by the demux pipeline.")
    if demux_type == ILMN_BCL2FASTQ_FASTQ_MODE:
        if len(samples) > 1:
            # ambiguous fastqprefix case
            if not fastqprefix:
                samples_list = "\n".join(samples)
                raise AmbiguousValueError("The --sample argument must be specified if multiple samples were demultiplexed in a run folder.  Options:\n%s" % samples_list)
            # no overlap case
            elif not set(samples).intersection(set(fastqprefixes)):
                raise ValueError("Samples not detected among demultiplexed FASTQs: %s" % fastqprefix_outstr)
            # some overlap; legal fastqprefix case
            else:
                return ILMN_BCL2FASTQ_FASTQ_MODE, fastqprefix_outstr
        # single sample does not match fastqprefixes case
        elif fastqprefix is not None and samples[0] not in fastqprefixes:
            raise ValueError("Samples not detected among FASTQs: %s" % fastqprefix_outstr)
        # no fastqprefix case-- return the lone sample
        elif fastqprefixes is None:
            return ILMN_BCL2FASTQ_FASTQ_MODE, samples[0]
        # legal fastqprefix list case
        else:
            return ILMN_BCL2FASTQ_FASTQ_MODE, fastqprefix_outstr
    else:
        return BCL_PROCESSOR_FASTQ_MODE, "any"


def check_fastq_types_multipath(fastq_paths, fastqprefix):
    """
    Validate that at least one path is compatible with the set of fastqprefix
    arguments, and return the input mode and superset of sample names that
    should be used downstream.  Forces the input mode to be the same for each
    path; if not, an error will be returned.  Will only raise an error in this
    case, and if no path-prefix combination will select any FASTQs.

    If, like in Cell Ranger, it is legal to have chunk-level fastq modes
    instead of pipeline-level fastq modes, it is preferred to call
    check_fastq_types individually over each path to get the input mode and
    sample name(s) for each chunk.
    """
    input_modes = set([])
    sample_names = set([])
    error_messages = set([])
    for path in fastq_paths:
        try:
            input_mode, samples = check_fastq_types(path, fastqprefix)
            input_modes.add(input_mode)
            sample_names.update(samples.split(','))
        # don't handle ambiguous values if detected in individual paths
        except AmbiguousValueError, e:
            raise e
        except ValueError, e:
            sys.stderr.write("Invalid path/prefix combination: %s, %s\n" % (path, str(fastqprefix)))
            error_messages.add(e.message)

    # happens if there were no legal selections
    if len(input_modes) == 0:
        if len(error_messages) == 1:
            raise ValueError(error_messages.pop())
        else:
            raise ValueError("FASTQ selection errors:\n%s" % ('\n'.join(error_messages)))
    elif len(input_modes) > 1:
        raise ValueError("Cannot process FASTQs at same time from different demultiplexing methods.")
    # can happen if multiple paths have different samples in them
    elif not fastqprefix and len(sample_names) > 1:
        raise AmbiguousValueError("The --sample argument must be specified if multiple samples were demultiplexed across the specified run folders.  Options:\n%s" % "\n".join(sorted(sample_names)))
    else:
        return input_modes.pop(), ",".join(sorted(sample_names))

def get_run_data(fn):
    """ Parse flowcell + lane from the first FASTQ record.
    NOTE: we don't check whether there are multiple FC / lanes in this file.
    NOTE: taken from longranger/mro/stages/reads/setup_chunks
    """
    if fn[-2:] == 'gz':
        reader = gzip.open(fn)
    else:
        reader = open(fn, 'r')

    gen = read_generator_fastq(reader)

    try:
        (name, seq, qual) = gen.next()

        # Abort if this is not an Illumina-like QNAME
        match = re.search('^([^:]+):([^:]+):([^:]+):([^:]+)', name)
        if match:
            flowcell, lane = match.group(3, 4)
        else:
            flowcell, lane = '', ''

        return (flowcell, lane)
    except StopIteration:
        # empty fastq
        raise ValueError('Could not extract flowcell and lane from FASTQ file. File is empty: %s' % fn)
