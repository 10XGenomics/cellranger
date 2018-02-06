#!/usr/bin/env python
#
# Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
#
import os
import subprocess
import tenkit.log_subprocess as tk_subproc
import cellranger.utils as cr_utils

class Bowtie2Reference:
    """
    A class for building bowtie2 indices and aligning reads to them using bowtie2. Depends on bowtie2 and samtools.
    """

    @staticmethod
    def load_from_fasta(reference_fasta_path, index_dir=None):
        b = Bowtie2Reference()
        if not os.path.exists(reference_fasta_path):
            raise ValueError('Reference file could not be found: %s' % reference_fasta_path)

        b.reference_fasta_path = reference_fasta_path
        if index_dir is None:
            # Index will be written in the same location as reference
            b.index_path = os.path.splitext(reference_fasta_path)[0]
        else:
            b.index_path = os.path.join(index_dir, os.path.splitext(os.path.basename(reference_fasta_path))[0])
        b.indexed = False

        b._index_reference(b.index_path)
        return b

    @staticmethod
    def load_from_index(index_path):
        b = Bowtie2Reference()
        b.reference_fasta_path = None

        # Check existence and validity of index
        try:
            tk_subproc.check_call('bowtie2-inspect -n %s' % index_path, shell=True)
        except subprocess.CalledProcessError:
            raise ValueError('Bowtie2 index could not be found or was invalid')

        b.index_path = index_path
        b.indexed = True
        return b

    def __init__(self):
        self.reference_fasta_path = None
        self.index_path = None
        self.indexed = False
        pass

    def _index_reference(self, index_path, **kwargs):
        """
        Generates a bowtie2 index for the specified reference file.

        Args:
            index_path (str): path to index prefix
            **kwargs: Any additional arguments to bowtie2-build may be included. Flags may have value set to None. Values
                are not validated. Parameters with hypens in name should be defined using underscores in place of hypens.

        Notes:
            Bowtie2 generates temporary files for indexing as a side-effect.

        Examples:
            kwargs can be specified as such: myBowtie2._index_reference(index_path, large_index=None, bmax=4)

        """

        additional_arguments = cr_utils.kwargs_to_command_line_options(set(), replace_chars={'_': '-'}, **kwargs)

        command = 'bowtie2-build %s %s %s' % (additional_arguments, self.reference_fasta_path, index_path)
        tk_subproc.check_call(command, shell=True)

        self.index_path = index_path
        self.indexed = True

    def align_reads_paired(self, in_fastq_r1_fn, in_fastq_r2_fn, out_file, write_sam=False, **kwargs):
        """
        Perform paired-end alignment of reads to reference and produce BAM as output using bowtie2.

        Args:
            in_fastq_r1_fn (str): name of fastq file with R1 to align
            in_fastq_r2_fn (str): name of fastq file with R2 to align
            out_file (str):       name of BAM/SAM file to output aligned reads
            write_sam (bool):     set to True to write SAM instead of BAM
            **kwargs: Any additional arguments to bowtie2 may be included.
                      Flags may have value set to None. Values are not
                      validated except for conflicts with index and read
                      input arguments. Parameters with hypens in name should
                      be defined using underscores in place of hypens.

        Examples:
            kwargs can be specified as such: myBowtie2.align_reads_paired(f1, f2, bam, p=2, very_fast=None, N=3)
        """

        assert self.indexed

        reserved_arguments = {'x', '1', '2', 'U'}
        additional_options = cr_utils.kwargs_to_command_line_options(reserved_arguments, replace_chars={'_': '-'}, **kwargs)

        if write_sam:
            cmd = 'bowtie2 %s -x %s -1 %s -2 %s -S %s' % \
                  (additional_options, self.index_path, in_fastq_r1_fn, in_fastq_r2_fn, out_file)
        else:
            cmd = 'bowtie2 %s -x %s -1 %s -2 %s | samtools view -bS - -o %s' % \
                  (additional_options, self.index_path, in_fastq_r1_fn, in_fastq_r2_fn, out_file)

        tk_subproc.check_call(cmd, shell=True)

    def align_reads_single(self, in_fastq_fn, out_file, write_sam=False, **kwargs):
        """
        Perform single-end alignment of reads to reference and produce BAM as output using bowtie2.

        Args:
            in_fastq_fn (str): name of fastq file with reads to align
            out_file (str):    name of BAM/SAM file to output aligned reads
            write_sam (bool):  set to True to write SAM instead of BAM
            **kwargs: Any additional arguments to bowtie2 may be included.
                      Flags may have value set to None. Values are not
                      validated except for conflicts with index and read
                      input arguments. Parameters with hypens in name should
                      be defined using underscores in place of hypens.

        Examples:
            kwargs can be specified as such: myBowtie2.align_reads_paired(fastq, bam, p=2, very_fast=None, N=3)
        """

        assert self.indexed

        reserved_arguments = {'x', '1', '2', 'U'}
        additional_options = cr_utils.kwargs_to_command_line_options(reserved_arguments, replace_chars={'_': '-'}, **kwargs)

        if write_sam:
            cmd = 'bowtie2 %s -x %s -U %s -S %s' % \
                  (additional_options, self.index_path, in_fastq_fn, out_file)
        else:
            cmd = 'bowtie2 %s -x %s -U %s | samtools view -bS - -o %s' % \
                  (additional_options, self.index_path, in_fastq_fn, out_file)

        tk_subproc.check_call(cmd, shell=True)

def get_bowtie2_orientation(paired_end, strand):
    """ Determine bowtie2 orientation option given library's paired-end status and strandedness.
        Return a dict whose keys contain bowtie2 args """
    result = {}
    if paired_end:
        result['fr'] = None

    if strand == '+':
        result['norc'] = None
    elif strand == '-':
        result['nofw'] = None
    return result
