#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
import gzip
import martian
import os
import tenkit.constants as tk_constants
import tenkit.bam as tk_bam
import tenkit.fasta as tk_fasta
import tenkit.preflight as tk_preflight
import cellranger.chemistry as cr_chem
import cellranger.constants as cr_constants

__MRO__ = '''
stage SETUP_CHUNKS(
    in  string sample_id,
    in  map[]  sample_def,
    in  string chemistry_name,
    in  map    custom_chemistry_def,
    out map[]  chunks,
    out map    chemistry_def,
    out string barcode_whitelist,
    src py     "stages/common/setup_chunks",
)
'''

def validate_fastq_lists(filename_lists):
    """ Return true if sample indices found """
    if set(filename_lists.keys()) != set(cr_constants.FASTQ_READ_TYPES.keys()):
        martian.log_info('Read types specified: ' + ','.join(sorted(filename_lists.keys())))
        martian.log_info('Read types expected: ' + ','.join(sorted(cr_constants.FASTQ_READ_TYPES.keys())))
        martian.exit('Not all read types were specified. Exiting pipeline.')

    if len(filename_lists.values()[0]) == 0:
        return False

    if not all(len(filename_list) == len(filename_lists.values()[0]) for filename_list in filename_lists.itervalues()):
        for read_type, read_description in cr_constants.FASTQ_READ_TYPES.keys():
            martian.log_info('%s files: %s.' % (read_description, str(filename_lists[read_type])))
        martian.exit('FASTQ files differ in number. Exiting pipeline.')
    return True

def construct_chunks(filename_lists, sample_id, sample_def, reads_interleaved, chemistry):
    chunks = []

    for chunk_idx in xrange(len(filename_lists.values()[0])):
        chunk = {
            'gem_group': sample_def['gem_group'],
            'reads_interleaved': reads_interleaved,
            'read_chunks': {},
            'chemistry': chemistry,
        }

        for read_type in cr_constants.FASTQ_READ_TYPES.keys():
            filename = filename_lists[read_type][chunk_idx]
            chunk['read_chunks'][read_type] = filename

        # Build read group (@RG) string
        # Infer flowcell, lane from first fastq
        first_fastq = [fq for fq in chunk['read_chunks'].values() if fq is not None][0]
        flowcell, lane = tk_fasta.get_run_data(first_fastq)
        library_id = sample_def.get('library_id', 'MissingLibrary')
        gem_group = str(sample_def['gem_group'] or 1)
        rg_string = tk_bam.pack_rg_string(sample_id, library_id, gem_group, flowcell, lane)
        chunk['read_group'] = rg_string

        chunks.append(chunk)

    return chunks

def fill_in_missing_reads(filename_lists):
    """ Provide a list of Nones for missing reads """
    max_filenames = max(len(filename_list) for filename_list in filename_lists.itervalues())
    for read_type, filename_list in filename_lists.iteritems():
        if len(filename_list) == 0:
            filename_lists[read_type] = [None] * max_filenames

def main_bcl_processor(sample_id, sample_def, chemistry_arg, custom_chemistry_def):
    chunks = []

    sample_index_strings, msg = tk_preflight.check_sample_indices(sample_def)
    if sample_index_strings is None:
        martian.exit(msg)

    path = sample_def['read_path']
    lanes = sample_def['lanes']

    for sample_index in sample_index_strings:
        # Determine the read-type => fastq filename mapping
        try:
            chemistry_name = cr_chem.infer_sc3p_chemistry_bcl_processor(chemistry_arg, path,
                                                                        sample_index, lanes)
        except cr_chem.NoInputFastqsException:
            continue

        if chemistry_name == cr_chem.CUSTOM_CHEMISTRY_NAME:
            chemistry = custom_chemistry_def
        else:
            chemistry = cr_chem.get_chemistry(chemistry_name)

        read_type_map = cr_chem.get_read_type_map(chemistry, tk_constants.BCL_PROCESSOR_FASTQ_MODE)

        # Collect the fastq files for each read type
        filename_lists = {}
        for dest_read_type in cr_constants.FASTQ_READ_TYPES:
            src_read_type = read_type_map[dest_read_type]
            filename_lists[dest_read_type] = tk_fasta.find_input_fastq_files_10x_preprocess(
                path, src_read_type, sample_index, lanes)

        fill_in_missing_reads(filename_lists)
        if validate_fastq_lists(filename_lists):
            chunks += construct_chunks(filename_lists, sample_id, sample_def, reads_interleaved=True,
                                       chemistry=chemistry)

    return chunks

def main_ilmn_bcl2fastq(sample_id, sample_def, chemistry_arg, custom_chemistry_def):
    chunks = []

    sample_names = sample_def['sample_names']
    path = sample_def['read_path']
    lanes = sample_def['lanes']

    for sample_name in sample_names:
        # Determine the read-type => fastq filename mapping
        try:
            chemistry_name = cr_chem.infer_sc3p_chemistry_ilmn_bcl2fastq(chemistry_arg, path,
                                                                         sample_name, lanes)
        except cr_chem.NoInputFastqsException:
            continue

        if chemistry_name == cr_chem.CUSTOM_CHEMISTRY_NAME:
            chemistry = custom_chemistry_def
        else:
            chemistry = cr_chem.get_chemistry(chemistry_name)

        read_type_map = cr_chem.get_read_type_map(chemistry, tk_constants.ILMN_BCL2FASTQ_FASTQ_MODE)

        # Collect the fastq files for each read type
        filename_lists = {}
        for dest_read_type in cr_constants.FASTQ_READ_TYPES:
            src_read_type = read_type_map[dest_read_type]
            filename_lists[dest_read_type] = tk_fasta.find_input_fastq_files_bcl2fastq_demult(
                path, src_read_type, sample_name, lanes)

        fill_in_missing_reads(filename_lists)
        if validate_fastq_lists(filename_lists):
            chunks += construct_chunks(filename_lists, sample_id, sample_def, reads_interleaved=False,
                                       chemistry=chemistry)

    return chunks

def main(args, outs):
    ok, msg = tk_preflight.check_gem_groups(args.sample_def)
    if not ok:
        martian.exit(msg)

    outs.chunks = []
    for sample_def in args.sample_def:
        fastq_mode = sample_def['fastq_mode']
        chunks = []

        if fastq_mode == tk_constants.BCL_PROCESSOR_FASTQ_MODE:
            chunks = main_bcl_processor(args.sample_id, sample_def, args.chemistry_name, args.custom_chemistry_def)
        elif fastq_mode == tk_constants.ILMN_BCL2FASTQ_FASTQ_MODE:
            chunks = main_ilmn_bcl2fastq(args.sample_id, sample_def, args.chemistry_name, args.custom_chemistry_def)
        else:
            martian.throw("Unrecognized fastq_mode: %s" % fastq_mode)

        if len(chunks) == 0:
            martian.exit(cr_constants.NO_INPUT_FASTQS_MESSAGE)

        outs.chunks += chunks

    if len(outs.chunks) == 0:
        martian.exit(cr_constants.NO_INPUT_FASTQS_MESSAGE)

    check_chunk_fastqs(outs.chunks)

    check_chunk_chemistries(outs.chunks)

    # Output chemistry and barcode whitelist
    outs.chemistry_def = outs.chunks[0]['chemistry']
    outs.barcode_whitelist = cr_chem.get_barcode_whitelist(outs.chemistry_def)


def check_fastq(fastq):
    # Check if fastq is readable
    if not os.access(fastq, os.R_OK):
        martian.exit("Do not have file read permission for FASTQ file: %s" % fastq)

    # Check if fastq is gzipped
    is_gzip_fastq = True
    try:
        with gzip.open(fastq) as f:
            f.read(1)
    except:
        is_gzip_fastq = False

    if is_gzip_fastq and not fastq.endswith(cr_constants.GZIP_SUFFIX):
        martian.exit("Input FASTQ file is gzipped but filename does not have %s suffix: %s" % (fastq, cr_constants.GZIP_SUFFIX))
    if not is_gzip_fastq and fastq.endswith(cr_constants.GZIP_SUFFIX):
        martian.exit("Input FASTQ file is not gzipped but filename has %s suffix: %s" % (fastq, cr_constants.GZIP_SUFFIX))

def check_chunk_fastqs(chunks):
    for chunk in chunks:
        for key in cr_constants.FASTQ_READ_TYPES:
            fastq = chunk.get(key)
            if fastq is not None:
                check_fastq(fastq)

def check_chunk_chemistries(chunks):
    """ Ensure all samples were generated with the same chemistry. """
    unique_chemistries = set([chunk['chemistry']['name'] for chunk in chunks])
    descriptions = map(cr_chem.get_chemistry_description_from_name, list(unique_chemistries))
    if len(unique_chemistries) > 1:
        martian.exit("Found multiple chemistries: %s. Combined analysis of libraries generated with different chemistries is not supported." % ', '.join(descriptions))
