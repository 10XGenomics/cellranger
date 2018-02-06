#!/usr/bin/env python
#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#

import log_subprocess
import json
import tenkit.bam as tk_bam
from tenkit.constants import FASTA_LOCATION

def count_vars_vcf(vcf_path):
    """ Counts the number of variants in a VCF file.
    """
    in_file = open(vcf_path)
    num_records = 0
    for line in in_file:
        if line.startswith('#'):
            continue
        num_records += 1
    return num_records

def split_alt_alleles_vcf(vcf_path, out_path):
    """ Splits records with more than one ALT field into two
    """
    out_file = open(out_path, 'w')
    log_subprocess.check_call(['vcfbreakmulti', vcf_path], stdout=out_file)

def output_primitives_vcf(vcf_path, out_path):
    """ Decomposes all complex variants into SNP and indel primitives
    """
    out_file = open(out_path, 'w')
    log_subprocess.check_call(['vcfallelicprimitives', vcf_path], stdout=out_file)

def output_restrict_location_vcf(vcf_path, bed_path, genome, out_path):
    """ Outputs a vcf restricted to the locations specified in the given bed file
    """
    ref_path = FASTA_LOCATION + genome + '/' + genome + '.fa'
    out_file = open(out_path, 'w')
    print ' '.join(['vcfintersect', '-b', bed_path, '-r', ref_path, vcf_path])
    log_subprocess.check_call(['vcfintersect', '-b', bed_path, '-r', ref_path, vcf_path], stdout=out_file)
    out_file.close()

def output_intersect_vcfs(vcf1_path, vcf2_path, genome, out_path):
    """ Outputs a vcf which is the intersection of the two given vcfs.
    """
    ref_path = FASTA_LOCATION + genome + '/' + genome + '.fa'
    out_file = open(out_path, 'w')
    log_subprocess.check_call(['vcfintersect', vcf1_path, '-i', vcf2_path, '-r', ref_path], stdout=out_file)
    out_file.close()

def output_setdiff_vcfs(vcf1_path, vcf2_path, genome, out_path):
    """ Outputs a VCF file which contains variants in the first but not the second VCF file.
    """
    out_file = open(out_path, 'w')
    ref_path = FASTA_LOCATION + genome + '/' + genome + '.fa'
    log_subprocess.check_call(['vcfintersect', vcf1_path, '-i', vcf2_path, '-v', '-r', ref_path], stdout=out_file)
    out_file.close()

def output_as_tsv(vcf_path, out_path, output_gt_info=False):
    """ Outputs all of the information from the vcf file as one big tsv
    """
    out_file = open(out_path, 'w')
    if output_gt_info:
        log_subprocess.check_call(['vcf2tsv', vcf_path, '-g'], stdout=out_file)
    else:
        log_subprocess.check_call(['vcf2tsv', vcf_path], stdout=out_file)
    out_file.close()

def output_variant_depths(vcf_path, bam_path, out_path):
    """ Outputs a JSON file with chrom and pos as keys and depths as values
    corresponding to the depths of the variants in the vcf file, for the
    sequencing run in the bam file.
    The bam file needs to be sorted and index.
    """
    in_file = open(vcf_path, 'r')
    bam_file = tk_bam.create_bam_infile(bam_path)
    depths = {}
    for line in in_file:
        if line.startswith('#'):
            continue
        info = line.split('\t')
        chrom = info[0]
        pos = int(info[1])
        chrom_depths = depths.setdefault(chrom, {})
        depth = 0
        for r in bam_file.fetch(chrom, pos, pos + 1):
            depth += 1
        chrom_depths[pos] = depth

    out_file = open(out_path, 'w')
    out_file.write(json.dumps(depths))
    out_file.close()

