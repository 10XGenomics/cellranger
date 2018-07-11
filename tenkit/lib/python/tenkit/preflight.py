#!/usr/bin/env python
#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#
# Functions for preflight checks
#
import os
import re
import resource
import socket
import subprocess
import tenkit.constants as tk_constants
import tenkit.seq as tk_seq
import tenkit.reference
import martian


def is_int(s):
    try:
        int(s)
    except ValueError:
        return False
    return True


def check_file(file_type, file_path, hostname):
    if not file_path.startswith('/'):
        martian.exit("Specified %s file must be an absolute path: %s" % (file_type, file_path))
    if not os.path.exists(file_path):
        martian.exit("On machine: %s, specified %s file does not exist: %s" % (hostname, file_type, file_path))
    if os.path.isdir(file_path):
        martian.exit("On machine: %s, specified %s file is a folder: %s" % (hostname, file_type, file_path))
    if not os.access(file_path, os.R_OK):
        martian.exit("On machine: %s, specified %s file is not readable: %s" % (hostname, file_type, file_path))


def check_folder(folder_type, folder_path, hostname, permission=os.X_OK):
    if not folder_path.startswith('/'):
        martian.exit("Specified %s folder must be an absolute path: %s" % (folder_type, folder_path))
    if not os.path.exists(folder_path):
        martian.exit("On machine: %s, specified %s folder does not exist: %s" % (hostname, folder_type, folder_path))
    if not os.path.isdir(folder_path):
        martian.exit("On machine: %s, specified %s path is not a folder: %s" % (hostname, folder_type, folder_path))
    if not os.access(folder_path, permission):
        martian.exit("On machine: %s, insufficient permissions on %s folder: %s" % (hostname, folder_type, folder_path))

def check_folder_or_create(folder_type, folder_path, hostname, permission=os.X_OK):
    if not folder_path.startswith('/'):
        martian.exit("Specified %s folder must be an absolute path: %s" % (folder_type, folder_path))
    if os.path.exists(folder_path):
        if not os.path.isdir(folder_path):
            martian.exit("On machine: %s, specified %s path is not a folder: %s" % (hostname, folder_type, folder_path))
        if not os.access(folder_path, permission):
            martian.exit("On machine: %s, insufficient permissions on %s folder: %s" % (hostname, folder_type, folder_path))
    else:
        try:
            os.makedirs(folder_path)
        except:
            martian.exit("On machine: %s, could not create %s folder: %s" % (hostname, folder_type, folder_path))

def check_rta_complete(folder_path):
    """
    :return: path to valid RTAComplete.txt in folder_path
    :rtype: string
    """
    hostname = socket.gethostname()
    check_folder("sequencing run", folder_path, hostname)
    rta_complete = os.path.join(folder_path, "RTAComplete.txt")
    if not os.path.exists(rta_complete):
        martian.exit("On machine: %s, run does not appear to be complete yet.  RTAComplete.txt not found." % hostname)
    return rta_complete


def check_runinfo_xml(folder_path):
    """
    :return: path to valid RunInfo.xml in folder_path
    :rtype: string
    """
    hostname = socket.gethostname()
    check_folder("sequencing run", folder_path, hostname)
    runinfo = os.path.join(folder_path, "RunInfo.xml")
    if not os.path.exists(runinfo):
        martian.exit("On machine: %s, RunInfo.xml not found. Cannot verify run was 10X-prepped." % hostname)
    if not os.access(runinfo, os.R_OK):
        martian.exit("On machine: %s, insufficient permission to open RunInfo.xml." % hostname)
    return runinfo


def check_barcode_whitelist(whitelist_path):
    hostname = socket.gethostname()
    check_file("barcode whitelist", whitelist_path, hostname)
    return whitelist_path


def check_bed(filename, reference_path):
    fasta = tenkit.reference.open_reference(reference_path)
    nvalid = 0
    max_lines = 1000000
    with open(filename, 'r') as f:
        for idx, line in enumerate(f):
            if line.startswith('browser') or line.startswith('track') or line.startswith('-browser') or line.startswith('-track') or line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 3:
                martian.exit('Error in BED file %s\nLine %d: Too few fields. chrom, start position, and end position are required. Maybe file is not tab-delimited.' % (filename, idx))
            chrom = fields[0]
            # Check that chromosome exists in reference
            if not chrom in fasta:
                martian.exit('Error in BED file %s\nLine %d: Chromosome %s not found in the reference genome.' % (filename, idx, chrom))
            if not is_int(fields[1]):
                martian.exit('Error in BED file %s\nLine %d: Non-integer starting position %s.' % (filename, idx, fields[1]))
            if not is_int(fields[2]):
                martian.exit('Error in BED file %s\nLine %d: Non-integer ending position %s.' % (filename, idx, fields[2]))
            start = int(fields[1])
            stop = int(fields[2])
            chrom_length = len(fasta[chrom])
            if start < 0 or start > chrom_length:
                martian.exit('Error in BED file %s\nLine %d: Invalid start position %d.' % (filename, idx, start))
            if stop < 0 or stop > chrom_length:
                martian.exit('Error in BED file %s\nLine %d: Invalid stop position %d.' % (filename, idx, stop))
            nvalid += 1
            if idx >= max_lines:
                break

    if nvalid == 0:
        martian.exit('Error in BED file %s: no valid regions.' % filename)



def check_refdata(reference_path, max_contigs = None):
    hostname = socket.gethostname()

    # Determine if the reference package is a known 10X reference package
    genome = tenkit.reference.get_genome(reference_path)
    known_genome = False

    if genome is not None and tenkit.reference.is_tenx(reference_path):
        known_genome = True

        version_path = os.path.join(reference_path, "version")
        if not os.path.exists(version_path):
            return False, "Your reference does not contain the expected files, or they are not readable. Please check your reference folder on %s." % hostname

        # Known genomes get a more stringent check
        if not os.path.exists(os.path.join(reference_path, "fasta/")) or \
                not os.path.exists(os.path.join(reference_path, "genes/")):
            return False, "Your reference does not contain the expected files, or they are not readable. Please check your reference folder on %s." % hostname
    else:
        # We only require the fasta for unknown genomes
        if not os.path.exists(os.path.join(reference_path, "fasta/")):
            return False, "Your reference does not contain the expected files, or they are not readable. Please check your reference folder on %s." % hostname

    if not os.path.exists(os.path.join(reference_path, "fasta/genome.fa.flat")) or \
            not os.path.exists(os.path.join(reference_path, "fasta/genome.fa.gdx")):
        return False, "Your reference doesn't appear to be indexed. Please run the mkreference tool."

    if os.path.getmtime(os.path.join(reference_path, "fasta/genome.fa.flat")) < os.path.getmtime(os.path.join(reference_path, "fasta/genome.fa")) or \
            os.path.getmtime(os.path.join(reference_path, "fasta/genome.fa.gdx")) < os.path.getmtime(os.path.join(reference_path, "fasta/genome.fa")):
        if known_genome:
            return False, "Timestamps on the FASTA index files have changed. Please reinstall the 10X refdata tar file on %s." % hostname
        else:
            return False, "Timestamps on the FASTA index files have changed. Please reindex your reference using the mkreference tool."


    fasta = tenkit.reference.open_reference(reference_path)
    num_contigs = len(fasta.keys())

    if max_contigs is not None and num_contigs > max_contigs:
        return False, "Long Ranger supports a maximum of %d reference contigs. Your reference contains %d. Please combine small contigs into a larger contig separated by N's." % (max_contigs, num_contigs)

    max_len = max([len(v) for (k,v) in fasta.iteritems()])

    logging = "reference path %s on %s contains genome: %s." % (reference_path, hostname, str(genome))
    logging += "reference contains %d contigs. max contig length: %d." % (num_contigs, max_len)

    if max_len >= (1<<29):
        return False, "Reference contains a contig longer than 536.8Mb (2^29 bp), which is not supported due to limitations of the .bai format. Please split this contig."

    # Check for ":" in contig names -- this will break things horribly
    has_colons = any(":" in ctg_name for ctg_name in fasta.keys())
    if has_colons:
        return False, "Reference names contain colon characters: ':'. References names containing colons are not supported."

    return True, logging



def check_open_fh():
    _, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
    if 0 <= hard and hard < tk_constants.MIN_PROCESS_NOFILE:
        return False, "On machine: %s, process open file handle hard limit (%d) is less than %d. Please run 'ulimit -n %d' before restarting the pipeline." % (
            socket.gethostname(), hard, tk_constants.MIN_PROCESS_NOFILE, tk_constants.MIN_PROCESS_NOFILE)

    if not os.path.exists(tk_constants.GLOBAL_NOFILE_PATH):
        return False, "On machine: %s, %s does not exist." % (socket.gethostname(), tk_constants.GLOBAL_NOFILE_PATH)
    with open(tk_constants.GLOBAL_NOFILE_PATH) as f:
        glob_str = f.read().strip()
    if not glob_str.isdigit():
        return False, "On machine: %s, %s contains a non-integer global open file handle limit: %s." % (
            socket.gethostname(), tk_constants.GLOBAL_NOFILE_PATH, glob_str)

    glob = int(glob_str)
    if glob < tk_constants.MIN_GLOBAL_NOFILE:
        return False, "On machine: %s, global open file handle limit (%d) is less than %d. Please set the global file handle limit to %d before restarting the pipeline." % (
            socket.gethostname(), glob, tk_constants.MIN_GLOBAL_NOFILE, tk_constants.MIN_GLOBAL_NOFILE)
    return True, None

def check_is_chromium(sample_item):
    bc_in_read = sample_item.get('bc_in_read', -1)
    bc_length = sample_item.get('bc_length', -1)
    return (bc_in_read == 1 and bc_length == 16)

def check_sample_indices(sample_item, sample_index_key = "sample_indices"):
    sample_indices = sample_item[sample_index_key]
    if type(sample_indices) != list:
        return None, "Sample indices must be of type list"
    if len(sample_indices) == 0:
        return None, "Sample indices must be a non-empty list"

    new_sample_indices = []
    for sample_index in sample_indices:
        if sample_index == "any":
            return ['*'], None
        elif tk_constants.SAMPLE_INDEX_MAP.get(sample_index):
            new_sample_indices.extend(tk_constants.SAMPLE_INDEX_MAP.get(sample_index))
        elif re.match("^[%s]+$" % "".join(tk_seq.NUCS), sample_index):
            new_sample_indices.append(sample_index)
        else:
            return None, ("Sample index '%s' is not valid. Must be one of: any, SI-<number>, "
                          "SI-<plate>-<well coordinate>, 220<part number>, or "
                          "a nucleotide sequence." % sample_index)

    return new_sample_indices, None

def check_gem_groups(sample_def):
    gem_groups = [sd['gem_group'] for sd in sample_def]
    all_null = all([x is None for x in gem_groups])
    all_int = all([type(x) is int for x in gem_groups])

    # Check for self-consistent gem_group settings in the sample_def entries
    if not (all_null or all_int):
        return False, "Inconsistent gem_group tags. Please specify all gem_group tags as null, or all gem_group tags with an integer."

    # If all gem_groups are set to null, then set them all to 1
    if all_null:
        for sd in sample_def:
            sd['gem_group'] = 1

    gem_groups_sorted = sorted([sd['gem_group'] for sd in sample_def])

    # Check that numbering starts at 1
    if len(gem_groups_sorted) > 0 and gem_groups_sorted[0] != 1:
        return False, 'gem_group numbering must start at 1'

    # Check for non-contiguous gem groups
    prev = 1
    for group in gem_groups_sorted:
        if group - prev > 1:
            return False, 'gem_groups must be numbered contiguously. missing groups: %s' % range(prev + 1, group)
        prev = group

    return True, None

def check_ld_library_path():
    if os.environ.get('_TENX_LD_LIBRARY_PATH') is None:
        return False, "Environment variable $_TENX_LD_LIBRARY_PATH is not set. Please enter the 10X environment before restarting the pipeline."
    return True, None

def record_package_versions():
    results = ""
    for package in tk_constants.PACKAGE_VERSION_CMDS:
        name = package['name']
        cmd = package['cmd']

        version = "not found"
        try:
            version = subprocess.check_output(cmd, shell=True)
        except:
            pass
        results += '%s: %s\n' % (name, version)
    return results
