#!/usr/bin/env python
#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#
# Run BWA and sort
#

import os
import tenkit.bam as tk_bam
import tenkit.fasta as tk_fasta
from exceptions import NotSupportedException
import subprocess
import log_subprocess

class Aligner(object):
    """
    Class for taking in a fastq file and running an aligner on it.
    """
    def __init__(self, in_fastq_file_name, output_name):
        self.in_fastq_file_name = in_fastq_file_name
        self.output_name = output_name

    def output_alignment(self, aligner=None, aligner_params=None, paired_end=True, genome=None, num_threads=16, read_group_header=None):
        """ Runs the given aligner on the input bam file
        If aligner is None, just copies the bam file to the new directory and then sorts and indexes it.
        """
        if aligner:
            if aligner == 'bwa':
                ref_fasta = aligner_params['ref_fasta']
                algorithm = aligner_params['algorithm']
                if paired_end:
                    bwa_align_paired(ref_fasta, self.in_fastq_file_name, self.output_name, algorithm=algorithm, num_threads=num_threads, read_group_header=read_group_header)
                else:
                    bwa_align_unpaired(ref_fasta, self.in_fastq_file_name, self.output_name, algorithm=algorithm, num_threads=num_threads, read_group_header=read_group_header)
            else:
                raise NotSupportedException('Alignment method not supported: ' + aligner)
        else:
            raise NotSupportedException('Must specify an aligner')

def bwa_read_group_arg(read_group, sample):
    ''' The -R argument to BWA must be escaped properly so that argument string contains two character \\t string, rather than actual tabs '''
    return '@RG\\tID:' + read_group + '\\tSM:' + sample

def bwa_index_ref(ref_fasta):
    """ Creates index of reference for bwa.  ref_fasta should be path to the reference fasta
    Only needs to be called once per reference.  Creates index files in the same directory as the
    reference
    """
    log_subprocess.check_call(['bwa', 'index', '-a', 'bwtsw', ref_fasta])

def bwa_align_unpaired(ref_fasta, read_fastq, out_name, algorithm='ALN', max_hits=None, read_group_header=None, num_threads=24):
    """ Runs bwa aligner on reads without using paired-information (using bam as input format).
    """
    assert(type(read_fastq) != list)

    if read_group_header is None:
        read_group_header = tk_bam.make_rg_header()

    if algorithm == 'MEM':
        # Temp file names
        sam_name = out_name + '.sam'

        sam_out_file = open(sam_name, 'w')
        log_subprocess.check_call(['bwa', 'mem', '-t', str(num_threads), '-M', '-R', read_group_header, ref_fasta, read_fastq], stdout=sam_out_file)
        sam_out_file.close()

        # Create final bam file from the sam file
        tk_bam.convert_to_bam(sam_name, out_name)

        # Remove temp files
        subprocess.check_call(['rm', sam_name])

    elif algorithm == 'ALN':
        # Temp file names
        sam_name = out_name + '.sam'
        index_name = out_name + '.sai'

        sam_out_file = open(sam_name, 'w')
        index_file = open(index_name, 'w')
        log_subprocess.check_call(['bwa', 'aln', '-t', str(num_threads), ref_fasta, read_fastq], stdout=index_file)
        index_file.close()
        if max_hits:
            log_subprocess.check_call(['bwa', 'samse', '-n', str(max_hits), ref_fasta, index_name, read_fastq], stdout=sam_out_file)
        else:
            log_subprocess.check_call(['bwa', 'samse', ref_fasta, index_name, read_fastq], stdout=sam_out_file)
        sam_out_file.close()

        # Create final bam file from the sam file
        tk_bam.convert_to_bam(sam_name, out_name)

        # Remove temp files
        subprocess.check_call(['rm', index_name])
        subprocess.check_call(['rm', sam_name])
    else:
        raise NotSupportedException('Unsupported bwa algorithm: ' + algorithm)


def bwa_align_paired(ref_fasta, read_fastq, out_name, algorithm='ALN', max_hits=None, read_group_header=None, num_threads=24):
    """Runs bwa paired-end aligner on reads using paired-end information
    Algorithm choices are currently
    MEM: Maximal Exact Matching (better for longer reads)
    ALN: Better for longer reads
    Haven't yet implemented BWA-SW
    Currently assumes the input read_fastq is in interleaved format, i.e. the reads of a pair
    are alternating.
    """
    if read_group_header is None:
        read_group_header = tk_bam.make_rg_header()

    if algorithm == 'MEM':
        devnull = open(os.devnull, 'w')
        if type(read_fastq) == list:
            assert(len(read_fastq) == 2)
            ## This restricts to primary alignments only
            out_file = open(out_name, 'w')
            ps = log_subprocess.Popen(['bwa', 'mem', '-t', str(num_threads), '-M',  '-R', read_group_header, ref_fasta, read_fastq[0], read_fastq[1]],
                                      stdout=subprocess.PIPE, stderr=devnull)
            #log_subprocess.check_call(['samtools', 'view', '-bSh', '-'], stdin=ps.stdout, stdout=out_file) # restore once bug fixed
            errors_file = open(out_name + '_ERRORS', 'w')
            log_subprocess.check_call(['samtools', 'view', '-bSh', '-'], stdin=ps.stdout, stdout=out_file, stderr=errors_file)
            out_file.close()
            errors_file.close()
        else:
            ## This restricts to primary alignments only
            out_file = open(out_name, 'w')
            ps = log_subprocess.Popen(['bwa', 'mem', '-p', '-t', str(num_threads), '-M',  '-R', read_group_header, ref_fasta, read_fastq],
                                      stdout=subprocess.PIPE, stderr=devnull)
            #log_subprocess.check_call(['samtools', 'view', '-bSh', '-'], stdin=ps.stdout, stdout=out_file) # restore once bug fixed
            errors_file = open(out_name + '_ERRORS', 'w')
            log_subprocess.check_call(['samtools', 'view', '-bSh', '-'], stdin=ps.stdout, stdout=out_file, stderr=errors_file)
            out_file.close()
            errors_file.close()
        devnull.close()

    elif algorithm == 'ALN':
        # Temp file names
        temp_fastq_name1 = out_name + '1.fastq'
        temp_fastq_name2 = out_name + '2.fastq'
        index_name_1 = out_name + '1.sai'
        index_name_2 = out_name + '2.sai'
        sam_name = out_name + '.sam'

        # Create the temp non-interleaved files
        in_fastq = open(read_fastq, 'r')
        temp_fastq1 = open(temp_fastq_name1, 'w')
        temp_fastq2 = open(temp_fastq_name2, 'w')
        tk_fasta.uninterleave_fastq(in_fastq, temp_fastq1, temp_fastq2)
        temp_fastq1.close()
        temp_fastq2.close()

        # Create the bwa index files
        index_file_1 = open(index_name_1, 'w')
        index_file_2 = open(index_name_2, 'w')
        log_subprocess.check_call(['bwa', 'aln', '-t', str(num_threads), ref_fasta, temp_fastq_name1], stdout=index_file_1)
        log_subprocess.check_call(['bwa', 'aln', '-t', str(num_threads), ref_fasta, temp_fastq_name2], stdout=index_file_2)
        index_file_1.close()
        index_file_2.close()

        # Create the sorted SAM file
        sam_out_file = open(sam_name, 'w')
        if max_hits:
            log_subprocess.check_call(['bwa', 'sampe', '-n', str(max_hits), ref_fasta, index_name_1, index_name_2, temp_fastq_name1, temp_fastq_name2], stdout=sam_out_file)
        else:
            log_subprocess.check_call(['bwa', 'sampe', ref_fasta, index_name_1, index_name_2, temp_fastq_name1, temp_fastq_name2], stdout=sam_out_file)

        sam_out_file.close()

        # Create final bam file from the sam file
        tk_bam.convert_to_bam(sam_name, out_name)

        # Clean up temporary files
        subprocess.check_call(['rm', temp_fastq_name1])
        subprocess.check_call(['rm', temp_fastq_name2])
        subprocess.check_call(['rm', index_name_1])
        subprocess.check_call(['rm', index_name_2])
        subprocess.check_call(['rm', sam_name])
    else:
        raise NotSupportedException('Unsupported bwa algorithm: ' + algorithm)
