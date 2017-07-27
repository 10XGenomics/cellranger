#!/usr/bin/env python
#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#
# Functions for tabix files
#

import pysam
import subprocess
import gzip

def sort_vcf(input_filename, output_filename):
    out_file = open(output_filename, 'w')
    subprocess.check_call(['vcfsort', input_filename], stdout=out_file)
    out_file.close()

def index_vcf(filename):
    pysam.tabix_index(filename, preset='vcf', force=True)

def sort_bc_loc_tabix(input_name, output_name, temp_dir_name=None):
    """ Sorts a tabix file by barcode when it is of the form
    chrom   start   end bc_seq
    If temp_dir_name is specified, tells the sort command to use that directory to store temporary files
    """
    if not(temp_dir_name is None):
        subprocess.check_call('grep \# ' + input_name + ' > ' + output_name + ';grep -v \# ' + input_name + ' | sort -k4  -T ' + temp_dir_name + ' >> ' + output_name, shell=True)
    else:
        subprocess.check_call('grep \# ' + input_name + ' > ' + output_name + ';grep -v \# ' + input_name + ' | sort -k4  >> ' + output_name, shell=True)

def create_tabix_infile(file_name):
    return pysam.Tabixfile(file_name)

def tabix_safe_fetch(tabix_file, chrom, start, stop):
    try:
        iterator = tabix_file.fetch(chrom, start, stop)
        return iterator
    except KeyError:
        return []
    except ValueError:
        return []

def sort_unique_tabix_vcf(vcf):
    ''' Sort, uniqueify, non-destructively bgzip, and tabix a VCF.'''
    tmp = vcf.rstrip('.vcf') + '.tmp.sorted.unique.vcf'
    subprocess.check_call('cat {0} | vcfstreamsort | vcfuniq > {1}'.format(vcf, tmp), shell=True)
    subprocess.check_call('cp {0} {1}'.format(tmp, vcf), shell=True)
    subprocess.check_call('rm -f {0}'.format(tmp), shell=True)
    subprocess.check_call("bgzip -c {0} > {0}.gz".format(vcf), shell=True)
    subprocess.check_call("tabix -p vcf {0}.gz".format(vcf), shell=True)

def sort_tabix_gtf(input_gtf, output_gtf):
    ''' Tabix a GTF. This is tricky because it needs to be sorted and block-gzipped first.
    NOTE: this will mess up the conventional ordering of genes/transcripts/exons/etc.'''
    cat = "zcat" if input_gtf.endswith(".gz") else "cat"
    subprocess.check_call('{0} {1} | sort -k1,1 -k4,5n | bgzip -c > {2}'.format(cat, input_gtf, output_gtf), shell=True)
    subprocess.check_call("tabix -p gff {0}".format(output_gtf), shell=True)
