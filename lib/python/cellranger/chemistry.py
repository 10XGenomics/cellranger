#!/usr/bin/env python
#
# Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
#
#
# Determine the locations of the cell barcode, cDNA sequence, UMI, sample index.
#
from collections import OrderedDict
import copy
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
    'endedness': cr_constants.THREE_PRIME,
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
    'endedness': cr_constants.THREE_PRIME,
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
}

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
    # Remove BC, UMI, spacer + 2bp
    'rna_read_offset': 41,
    'rna_read_length': None,
    'rna_read2_type': 'R2',
    'rna_read2_offset': 0,
    'rna_read2_length': None,
    'si_read_type': 'I1',
    'si_read_offset': 0,
    'si_read_length': None,
    'strandedness': '+',
    'endedness': cr_constants.FIVE_PRIME,
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
}

# Single Cell V(D)J, transcript on R2 only
CHEMISTRY_SCVDJ_R2 = {
    'name': 'SCVDJ-R2',
    'description': 'Single Cell V(D)J R2-only',
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
    'strandedness': '-',
    'endedness': cr_constants.FIVE_PRIME,
    'read_type_to_bcl_processor_filename': {
        'R1': 'RA', # Read1, Read2 interleaved
        'R2': None,
        'I1': 'I1', # Index7
        'I2': None, # Index5
    },
    # Valid for the following argument to bcl2fastq:
    # --use-bases-mask=Y26,I8,Y150
    'read_type_to_bcl2fastq_filename': {
        'R1': 'R1', # Read1
        'R2': 'R2', # Read2
        'I1': 'I1', # Index7
        'I2': None, # Index5
    },
    'barcode_whitelist': '737K-august-2016',
}

# Single Cell 5' Gene Expression, paired-end
CHEMISTRY_SC5P_PE = {
    'name': 'SC5P-PE',
    'description': 'Single Cell 5\' PE',
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
    'endedness': cr_constants.FIVE_PRIME,
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
}

# Single Cell 5' Gene Expression, transcript on R2 only
CHEMISTRY_SC5P_R2 = {
    'name': 'SC5P-R2',
    'description': 'Single Cell 5\' R2-only',
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
    'strandedness': '-',
    'endedness': cr_constants.FIVE_PRIME,
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
}

# Single Cell 5' Gene Expression, transcript on R1 only
CHEMISTRY_SC5P_R1 = {
    'name': 'SC5P-R1',
    'description': 'Single Cell 5\' R1-only',
    'barcode_read_type': 'R1',
    'barcode_read_offset': 0,
    'barcode_read_length': 16,
    'umi_read_type': 'R1',
    'umi_read_offset': 16,
    'umi_read_length': 10,
    'rna_read_type': 'R1',
    # Remove BC, UMI, spacer + 2bp
    'rna_read_offset': 41,
    'rna_read_length': None,
    'rna_read2_type': None,
    'rna_read2_offset': 0,
    'rna_read2_length': None,
    'si_read_type': 'I1',
    'si_read_offset': 0,
    'si_read_length': None,
    'strandedness': '+',
    'endedness': cr_constants.FIVE_PRIME,
    # NOTE: Interleaved input is not supported!
    'read_type_to_bcl_processor_filename': {
        'R1': 'RA', # Read1, Read2 interleaved
        'R2': None,
        'I1': 'I1', # Index7
        'I2': None, # Index5
    },
    # Valid for the following argument to bcl2fastq:
    # --use-bases-mask=Y137,I8
    'read_type_to_bcl2fastq_filename': {
        'R1': 'R1', # Read1
        'R2': 'R2', # Read2
        'I1': 'I1', # Index7
        'I2': None, # Index5
    },
    'barcode_whitelist': '737K-august-2016',
}


SC3P_CHEMISTRIES = [
    CHEMISTRY_SC3P_V1,
    CHEMISTRY_SC3P_V2,
]

SC5P_CHEMISTRIES = [
    CHEMISTRY_SC5P_PE,
    CHEMISTRY_SC5P_R2,
    CHEMISTRY_SC5P_R1,
]

DEFINED_CHEMISTRIES = SC3P_CHEMISTRIES + [CHEMISTRY_SCVDJ, CHEMISTRY_SCVDJ_R2] + SC5P_CHEMISTRIES

# Aliases for human usage
CHEMISTRY_ALIASES = OrderedDict([
    ('threeprime', 'SC3P_auto'),
    ('fiveprime', 'SC5P_auto'),
])

# Signals the pipeline to auto-detect the chemistry
AUTO_CHEMISTRY_NAMES = ['auto', 'threeprime', 'fiveprime', 'SC3P_auto', 'SC5P_auto', 'SCVDJ_auto']

# User-defined chemistry (use the various read_type/offset/length arguments passed to the pipeline)
CUSTOM_CHEMISTRY_NAME = 'custom'

# All valid chemistry name arguments to the pipeline
CHEMISTRY_ARG_NAMES = CHEMISTRY_ALIASES.keys() + \
                      [c['name'] for c in DEFINED_CHEMISTRIES] + \
                      AUTO_CHEMISTRY_NAMES + \
                      [CUSTOM_CHEMISTRY_NAME]

def check_chemistry_arg(chemistry_name, allowed_chems=None):
    """ Validate a given chemistry name argument passed to the pipeline. Return (ok, error_msg) """
    check_against = allowed_chems or CHEMISTRY_ARG_NAMES

    if chemistry_name in check_against:
        return True, None

    display_list = copy.copy(check_against)
    display_list.remove(CUSTOM_CHEMISTRY_NAME)

    return (False, 'Unrecognized chemistry name: "%s". Must be one of: %s.' % (chemistry_name, ', '.join(display_list)))

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

def get_endedness(chemistry):
    return chemistry.get('endedness', cr_constants.THREE_PRIME)

def get_description(chemistry):
    return chemistry['description']

def get_chemistry(name):
    """ Returns a chemistry definition dict for a given name """

    name = CHEMISTRY_ALIASES.get(name, name)

    chemistries = filter(lambda c: c['name'] == name, DEFINED_CHEMISTRIES)
    if len(chemistries) == 0:
        raise ValueError("Could not find chemistry named %s" % name)
    if len(chemistries) > 1:
        raise ValueError("Found multiple chemistries with name %s" % name)
    return chemistries[0]

def get_chemistry_description_from_name(name):
    if name in AUTO_CHEMISTRY_NAMES or name == CUSTOM_CHEMISTRY_NAME:
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
                                             None, None)

        for read in barcode_reads.in_iter:
            if num_reads == cr_constants.DETECT_CHEMISTRY_INITIAL_READS:
                break

            _, barcode, _ = read

            num_reads += 1
            if barcode in barcode_whitelist_set:
                barcodes_on_whitelist += 1

        if num_reads == cr_constants.DETECT_CHEMISTRY_INITIAL_READS:
            break

    if num_reads > 0:
        return tk_stats.robust_divide(barcodes_on_whitelist, num_reads)
    else:
        return 0.0

def _compute_r1_length(fastqs, reads_interleaved):
    """ Infer the length of R1 """
    num_reads = 0
    r1_max_len = 0

    def get_r1_noninterleaved(read_iter):
        for _, seq, _ in read_iter:
            yield seq
    def get_r1_interleaved(read_iter):
        for _, seq, _, _, _, _ in read_iter:
            yield seq
    get_r1 = get_r1_interleaved if reads_interleaved else get_r1_noninterleaved

    for fastq in fastqs:
        with cr_utils.open_maybe_gzip(fastq, 'r') as fq_file:
            reads = tk_fasta.read_generator_fastq(fq_file, reads_interleaved)

            for r1 in get_r1(reads):
                if num_reads == cr_constants.DETECT_CHEMISTRY_INITIAL_READS:
                    break
                r1_max_len = max(len(r1), r1_max_len)
                num_reads += 1

        if num_reads == cr_constants.DETECT_CHEMISTRY_INITIAL_READS:
            break

    return r1_max_len

def infer_chemistry(chemistry_arg, fq_spec):
    """ Infers the chemistry name given a FastqSpec for a single sample index/name """
    assert fq_spec.is_single_group()

    chemistry_arg = CHEMISTRY_ALIASES.get(chemistry_arg, chemistry_arg)

    if chemistry_arg == 'SC3P_auto':
        return infer_sc3p_chemistry(fq_spec)

    elif chemistry_arg == 'SC5P_auto':
        return infer_sc5p_chemistry(fq_spec)

    elif chemistry_arg == 'SCVDJ_auto':
        return infer_scvdj_chemistry(fq_spec)

    else:
        # Must be explicit from this point on
        assert chemistry_arg not in AUTO_CHEMISTRY_NAMES
        return chemistry_arg


def infer_sc3p_chemistry(fq_spec):
    """ Infer the SC3P chemistry name.
        fq_spec (FastqSpec) for a single sample index/name """
    assert fq_spec.is_single_group()

    best_chem, best_frac = None, float('-inf')

    n_fastqs = 0

    for chemistry in SC3P_CHEMISTRIES:
        # Get the FASTQs containing the barcode for this chemistry
        barcode_read_def = get_barcode_read_def(chemistry)
        read_type = get_read_type_map(chemistry, fq_spec.fastq_mode)[barcode_read_def.read_type]

        fastqs = fq_spec.get_fastqs(read_type)
        n_fastqs += len(fastqs)

        whitelist = _get_barcode_whitelist_set(chemistry)

        wl_frac = _compute_frac_barcodes_on_whitelist(fastqs,
                                                      whitelist,
                                                      reads_interleaved=fq_spec.interleaved,
                                                      read_def=barcode_read_def)

        if wl_frac > best_frac:
            best_chem, best_frac = chemistry['name'], wl_frac

    if n_fastqs == 0:
        raise NoInputFastqsException()

    if best_frac >= cr_constants.DETECT_CHEMISTRY_MIN_FRAC_WHITELIST:
        return best_chem
    else:
        raise ValueError('Could not auto-detect Single Cell 3\' chemistry. Fraction of barcodes on whitelist was at best %0.2f%%, while we expected at least %0.2f%% for one of the chemistries.' %
                (100.0*best_frac, 100.0*cr_constants.DETECT_CHEMISTRY_MIN_FRAC_WHITELIST))


def infer_sc5p_chemistry(fq_spec):
    """ Infer the SC5P chemistry from the R1 length.
        fq_spec (FastqSpec) - for a single sample index/name """
    assert fq_spec.is_single_group()

    r1_fastqs = fq_spec.get_fastqs('R1')

    if len(r1_fastqs) == 0:
        raise NoInputFastqsException()

    r1_len = _compute_r1_length(r1_fastqs, reads_interleaved=fq_spec.interleaved)


    if fq_spec.interleaved:
        # Single-end + interleaved is unsupported.
        single_end = False
    else:
        r2_fastqs = fq_spec.get_fastqs('R2')
        single_end = len(r2_fastqs) == 0

    if single_end:
        return 'SC5P-R1'
    elif r1_len >= cr_constants.DETECT_5P_CHEMISTRY_MIN_R1_LEN_PE:
        return 'SC5P-PE'
    else:
        return 'SC5P-R2'

def infer_scvdj_chemistry(fq_spec):
    """ Infer the SCVDJ chemistry from the R1 length.
        fq_spec (FastqSpec) - for a single sample index/name """
    assert fq_spec.is_single_group()

    r1_fastqs = fq_spec.get_fastqs('R1')

    if len(r1_fastqs) == 0:
        raise NoInputFastqsException()

    r1_len = _compute_r1_length(r1_fastqs, reads_interleaved=fq_spec.interleaved)


    if fq_spec.interleaved:
        # Single-end + interleaved is unsupported.
        single_end = False
    else:
        r2_fastqs = fq_spec.get_fastqs('R2')
        single_end = len(r2_fastqs) == 0

    if single_end:
        raise ValueError("Single-end (SC5P-R1) is not supported for V(D)J.")
    elif r1_len >= cr_constants.DETECT_VDJ_CHEMISTRY_MIN_R1_LEN_PE:
        return 'SCVDJ'
    else:
        return 'SCVDJ-R2'
