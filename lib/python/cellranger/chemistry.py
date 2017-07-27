#!/usr/bin/env python
#
# Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
#

# Determine the locations of the cell barcode, cDNA sequence, UMI, sample index.

import tenkit.constants as tk_constants
import tenkit.fasta as tk_fasta
import tenkit.stats as tk_stats
import cellranger.constants as cr_constants
import cellranger.fastq as cr_fastq
import cellranger.utils as cr_utils

class NoInputFastqsException(Exception):
    pass

# Single Cell 3' chemistries
CHEMISTRY_SC3P_V1 = {
    'name': 'SC3Pv1',
    'description': 'Single Cell 3\' v1',
    'barcode_read_type': 'I1',
    'barcode_read_offset': 0,
    'barcode_read_length': 14,
    'umi_read_type': 'R2',
    'umi_read_offset': 0,
    'umi_read_length': 10,
    'rna_read_type': 'R1',
    'rna_read_offset': 0,
    'rna_read_length': None, # To end of sequence
    'rna_read2_type': None,
    'rna_read2_offset': 0,
    'rna_read2_length': None,
    'si_read_type': 'I2',
    'si_read_offset': 0,
    'si_read_length': None,
    'strandedness': '+',
    'read_type_to_bcl_processor_filename': {
        'R1': 'RA', # Read1, Read2 interleaved
        'R2': None,
        'I1': 'I1', # Index7
        'I2': 'I2', # Index5
    },
    # Valid for the following argument to bcl2fastq:
    # --use-bases-mask=Y98,Y14,I8,Y10
    'read_type_to_bcl2fastq_filename': {
        'R1': 'R1', # Read1
        'R2': 'R3', # Read2
        'I1': 'R2', # Index7
        'I2': 'I1', # Index5
    },
    'barcode_whitelist': '737K-april-2014_rc',
    # Retain sequence after the UMI in R2
    'retain_trimmed_suffix_read': 'R2',
}

CHEMISTRY_SC3P_V2 = {
    'name': 'SC3Pv2',
    'description': 'Single Cell 3\' v2',
    'barcode_read_type': 'R1',
    'barcode_read_offset': 0,
    'barcode_read_length': 16,
    'umi_read_type': 'R1',
    'umi_read_offset': 16,
    'umi_read_length': 10,
    'rna_read_type': 'R2',
    'rna_read_offset': 0,
    'rna_read_length': None,
    'rna_read2_type': None,
    'rna_read2_offset': 0,
    'rna_read2_length': None,
    'si_read_type': 'I1',
    'si_read_offset': 0,
    'si_read_length': None,
    'strandedness': '+',
    'read_type_to_bcl_processor_filename': {
        'R1': 'RA', # Read1, Read2 interleaved
        'R2': None,
        'I1': 'I1', # Index7
        'I2': None, # Index5
    },
    # Valid for the following argument to bcl2fastq:
    # --use-bases-mask=Y26,I8,Y98
    'read_type_to_bcl2fastq_filename': {
        'R1': 'R1', # Read1
        'R2': 'R2', # Read2
        'I1': 'I1', # Index7
        'I2': None, # Index5
    },
    'barcode_whitelist': '737K-august-2016',
    # Retain sequqence after the BC-UMI in R1
    'retain_trimmed_suffix_read': 'R1',
}

SC3P_CHEMISTRIES = [
    CHEMISTRY_SC3P_V1,
    CHEMISTRY_SC3P_V2,
]

# Single Cell V(D)J
CHEMISTRY_SCVDJ = {
    'name': 'SCVDJ',
    'description': 'Single Cell V(D)J',
    'barcode_read_type': 'R1',
    'barcode_read_offset': 0,
    'barcode_read_length': 16,
    'umi_read_type': 'R1',
    'umi_read_offset': 16,
    'umi_read_length': 10,
    'rna_read_type': 'R1',
    'rna_read_offset': 26,
    'rna_read_length': None,
    'rna_read2_type': 'R2',
    'rna_read2_offset': 0,
    'rna_read2_length': None,
    'si_read_type': 'I1',
    'si_read_offset': 0,
    'si_read_length': None,
    'strandedness': '+',
    'read_type_to_bcl_processor_filename': {
        'R1': 'RA', # Read1, Read2 interleaved
        'R2': None,
        'I1': 'I1', # Index7
        'I2': None, # Index5
    },
    # Valid for the following argument to bcl2fastq:
    # --use-bases-mask=Y150,I8,Y150
    'read_type_to_bcl2fastq_filename': {
        'R1': 'R1', # Read1
        'R2': 'R2', # Read2
        'I1': 'I1', # Index7
        'I2': None, # Index5
    },
    'barcode_whitelist': '737K-august-2016',
    'retain_trimmed_suffix_read': None,
}

DEFINED_CHEMISTRIES = SC3P_CHEMISTRIES + [CHEMISTRY_SCVDJ]

# Signals the pipeline to auto-detect the chemistry
AUTO_CHEMISTRY_NAME = 'auto'

# User-defined chemistry (use the various read_type/offset/length arguments passed to the pipeline)
CUSTOM_CHEMISTRY_NAME = 'custom'

# All valid chemistry arguments to the pipeline
CHEMISTRY_ARG_NAMES = [c['name'] for c in DEFINED_CHEMISTRIES] + [AUTO_CHEMISTRY_NAME, CUSTOM_CHEMISTRY_NAME]

def check_chemistry_arg(chemistry_name):
    """ Validate a given chemistry name argument passed to the pipeline. Return (ok, error_msg) """
    if any(name == chemistry_name for name in CHEMISTRY_ARG_NAMES):
        return True, None
    return (False, 'Unrecognized chemistry name: "%s". Must be one of: %s.' % (chemistry_name, ', '.join(CHEMISTRY_ARG_NAMES)))

def check_chemistry_def(chemistry):
    """ Validate a chemistry definition dict. Return (ok, error_msg) """
    required_keys = [
        'description',
        'barcode_read_type', 'barcode_read_offset', 'barcode_read_length',
        'rna_read_type', 'rna_read_offset', 'rna_read_length',
        'rna_read2_type', 'rna_read2_offset', 'rna_read2_length',
        'umi_read_type', 'umi_read_offset', 'umi_read_length',
        'si_read_type', 'si_read_offset', 'si_read_length',
        'read_type_to_bcl_processor_filename',
        'read_type_to_bcl2fastq_filename',
        'strandedness',
    ]

    if 'name' not in chemistry:
        return (False, "Chemistry definition is missing a 'name' key")

    missing_keys = set(required_keys) - set(chemistry.keys())
    if len(missing_keys) > 0:
        return (False, "Keys missing from chemistry definition %s: %s" % (chemistry['name'], ','.join(sorted(list(missing_keys)))))

    return True, None

def check_chemistry_defs():
    """ Validate the predefined chemistry defs in this file """
    for chemistry in DEFINED_CHEMISTRIES:
        ok, msg = check_chemistry_def(chemistry)
        if not ok:
            return ok, msg

    names = [c['name'] for c in DEFINED_CHEMISTRIES]
    if len(set(names)) != len(names):
        return (False, "Chemistries share a name.")

    return True, None

def _get_read_def(chemistry, read_prefix):
    """ Retrieve a read definition from a chemistry or from a given custom chemistry """
    return cr_constants.ReadDef(chemistry['%s_type' % read_prefix],
                                chemistry['%s_offset' % read_prefix],
                                chemistry['%s_length' % read_prefix])

def get_barcode_read_def(chemistry):
    return _get_read_def(chemistry, 'barcode_read')

def get_rna_read_def(chemistry):
    return _get_read_def(chemistry, 'rna_read')

def get_rna_read2_def(chemistry):
    return _get_read_def(chemistry, 'rna_read2')

def is_paired_end(chemistry):
    return get_rna_read2_def(chemistry).read_type is not None

def get_umi_read_def(chemistry):
    return _get_read_def(chemistry, 'umi_read')

def has_umis(chemistry):
    return get_umi_read_def(chemistry).read_type is not None

def get_umi_length(chemistry):
    if has_umis(chemistry):
        return get_umi_read_def(chemistry).length
    else:
        return 0

def get_si_read_def(chemistry):
    return _get_read_def(chemistry, 'si_read')

def get_barcode_whitelist(chemistry):
    return chemistry['barcode_whitelist']

def _get_barcode_whitelist_set(chemistry):
    return set(cr_utils.load_barcode_whitelist(get_barcode_whitelist(chemistry)))

def get_read_type_map(chemistry, fastq_mode):
    """ Get the mapping of read type to fastq filename for a given chemistry. """
    if fastq_mode == tk_constants.BCL_PROCESSOR_FASTQ_MODE:
        return chemistry['read_type_to_bcl_processor_filename']
    elif fastq_mode == tk_constants.ILMN_BCL2FASTQ_FASTQ_MODE:
        return chemistry['read_type_to_bcl2fastq_filename']
    else:
        raise ValueError("Unrecognized FASTQ mode: %s" % fastq_mode)

def get_strandedness(chemistry):
    return chemistry['strandedness']

def get_description(chemistry):
    return chemistry['description']

def get_chemistry(name):
    """ Returns a chemistry definition dict for a given name """
    chemistries = filter(lambda c: c['name'] == name, DEFINED_CHEMISTRIES)
    if len(chemistries) == 0:
        raise ValueError("Could not find chemistry named %s" % name)
    if len(chemistries) > 1:
        raise ValueError("Found multiple chemistries with name %s" % name)
    return chemistries[0]

def get_chemistry_description_from_name(name):
    if name == AUTO_CHEMISTRY_NAME or name == CUSTOM_CHEMISTRY_NAME:
        return name
    else:
        return get_chemistry(name)['description']


def _compute_frac_barcodes_on_whitelist(fastqs, barcode_whitelist_set, reads_interleaved, read_def):
    """ Compute fraction of observed barcodes on the barcode whitelist """
    num_reads = 0
    barcodes_on_whitelist = 0

    for fastq in fastqs:
        barcode_reads = cr_fastq.FastqReader({read_def.read_type: fastq},
                                             read_def,
                                             reads_interleaved,
                                             None)

        for extraction in barcode_reads.in_iter:
            if num_reads == cr_constants.DETECT_CHEMISTRY_INITIAL_READS:
                break

            _, barcode, _ = extraction.read

            num_reads += 1
            if barcode in barcode_whitelist_set:
                barcodes_on_whitelist += 1

        if num_reads == cr_constants.DETECT_CHEMISTRY_INITIAL_READS:
            break

    if num_reads > 0:
        return tk_stats.robust_divide(barcodes_on_whitelist, num_reads)
    else:
        return 0.0


def _infer_sc3p_chemistry(chemistry_whitelist_frac):
    """ Infer the SC3P chemistry name from a dict of <chemistry_name>: frac_barcodes_on_whitelist. """
    best_chemistry_name, best_whitelist_frac = max(chemistry_whitelist_frac.items(), key=lambda kv: kv[1])

    if best_whitelist_frac >= cr_constants.DETECT_CHEMISTRY_MIN_FRAC_WHITELIST:
        return best_chemistry_name
    else:
        raise ValueError('Could not auto-detect Single Cell 3\' chemistry. Fraction of barcodes on whitelist was at best %0.2f%%, while we expected at least %0.2f%% for one of the chemistries.' %
                (100.0*best_whitelist_frac, 100.0*cr_constants.DETECT_CHEMISTRY_MIN_FRAC_WHITELIST))

def infer_sc3p_chemistry_bcl_processor(chemistry_arg, fastq_path, sample_index, lanes):
    """ Infer the SC3P chemistry name from the set of fastqs provided by BCL_PROCESSOR/mkfastq. """

    if chemistry_arg != AUTO_CHEMISTRY_NAME:
        return chemistry_arg

    chemistry_whitelist_fracs = {}
    num_fastqs_found = 0

    for chemistry in SC3P_CHEMISTRIES:
        barcode_read_def = get_barcode_read_def(chemistry)

        fastq_read = get_read_type_map(chemistry, tk_constants.BCL_PROCESSOR_FASTQ_MODE)[barcode_read_def.read_type]

        fastqs = tk_fasta.find_input_fastq_files_10x_preprocess(fastq_path, fastq_read,
                                                                sample_index, lanes)
        num_fastqs_found += len(fastqs)

        chemistry_whitelist_fracs[chemistry['name']] = _compute_frac_barcodes_on_whitelist(fastqs,
                                                                                          _get_barcode_whitelist_set(chemistry),
                                                                                          reads_interleaved=True,
                                                                                          read_def=barcode_read_def)

    if num_fastqs_found == 0:
        raise NoInputFastqsException()

    return _infer_sc3p_chemistry(chemistry_whitelist_fracs)

def infer_sc3p_chemistry_ilmn_bcl2fastq(chemistry_arg, fastq_path, sample_name, lanes):
    """ Infer the SC3P chemistry name from the set of fastqs provided by bcl2fastq. """

    if chemistry_arg != AUTO_CHEMISTRY_NAME:
        return chemistry_arg

    chemistry_whitelist_fracs = {}
    num_fastqs_found = 0

    for chemistry in SC3P_CHEMISTRIES:
        barcode_read_def = get_barcode_read_def(chemistry)

        fastq_read = get_read_type_map(chemistry, tk_constants.ILMN_BCL2FASTQ_FASTQ_MODE)[barcode_read_def.read_type]

        fastqs = tk_fasta.find_input_fastq_files_bcl2fastq_demult(fastq_path, fastq_read,
                                                                  sample_name, lanes)
        num_fastqs_found += len(fastqs)

        chemistry_whitelist_fracs[chemistry['name']] = _compute_frac_barcodes_on_whitelist(fastqs,
                                                                                          _get_barcode_whitelist_set(chemistry),
                                                                                          reads_interleaved=False,
                                                                                          read_def=barcode_read_def)

    if num_fastqs_found == 0:
        raise NoInputFastqsException

    return _infer_sc3p_chemistry(chemistry_whitelist_fracs)
