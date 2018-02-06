#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#
import numpy as np
import os.path

# Chain types
TR_CHAIN_TYPE = 'TR'
IG_CHAIN_TYPE = 'IG'
VDJ_CHAIN_TYPES = [TR_CHAIN_TYPE, IG_CHAIN_TYPE]

ALL_CHAIN_TYPES = 'all'    # Do not filter on chain type
AUTO_CHAIN_TYPE = 'auto'  # Autodetect the chain type filter
CHAIN_TYPE_SPECS = VDJ_CHAIN_TYPES + [ALL_CHAIN_TYPES, AUTO_CHAIN_TYPE]

# Chains (inconsistently called "genes" in the code)
VDJ_GENES = ['TRA', 'TRB', 'TRG', 'TRD', 'IGH', 'IGK', 'IGL']
CANONICAL_VDJ_GENES = ['TRA', 'TRB', 'IGH', 'IGK', 'IGL']

CANONICAL_TR_GENES = ['TRA', 'TRB']

CHAINS_WITH_ISOTYPES = ['IGH']

# Chain pairs (inconsistently called "gene pairs" in the code)
VDJ_GENE_PAIRS = [
    ('TRA_TRB'),
    ('TRG_TRD'),
    ('IGK_IGH'),
    ('IGL_IGH'),
]
# Used to exclude cells from the per-chain median UMI counts per cell calculation
EXCLUSIVE_VDJ_GENES = [
    ['IGK', 'IGL'],
]
CANONICAL_VDJ_GENE_PAIRS = [
    ('TRA_TRB'),
    ('IGK_IGH'),
    ('IGL_IGH'),
]
CANONICAL_TR_GENE_PAIRS = [
    ('TRA_TRB'),
]

# VDJ read to contig alignment constants
VDJ_ALIGN_READ_TO_CONTIGS_MAX_INSERT_SIZE = 1000
VDJ_MISMATCH_BREAKS = range(0, 22, 2)
VDJ_DEPTH_BREAKS = range(0, 20, 1)
VDJ_SOFTCLIP_BREAKS = range(0, 52, 2)
VDJ_CONTIG_LENGTH_PERCENTILES = np.arange(0, 1.1, 0.05)

# VDJ annotation constants
VDJ_ANNOTATION_MATCH_SCORE = 2
VDJ_ANNOTATION_MISMATCH_PENALTY = 3
VDJ_ANNOTATION_GAP_OPEN_PENALTY = 5
VDJ_ANNOTATION_EXTEND_PENALTY = 1
VDJ_ANNOTATION_MIN_SCORE_RATIO = 0.75
VDJ_ANNOTATION_MIN_WORD_SIZE = 15
VDJ_PRIMER_ANNOTATION_MIN_FRACTION_MATCHED = 0.8

VDJ_5U_FEATURE_TYPES = ['5U', "5'UTR"]
VDJ_V_FEATURE_TYPES = ['V', 'L', 'V-REGION', 'L-REGION+V-REGION']
VDJ_D_FEATURE_TYPES = ['D', 'D-REGION']
VDJ_J_FEATURE_TYPES = ['J', 'J-REGION']
VDJ_C_FEATURE_TYPES = ['C', 'C-REGION']
VDJ_ORDERED_REGIONS = [VDJ_5U_FEATURE_TYPES,
                       VDJ_V_FEATURE_TYPES,
                       VDJ_D_FEATURE_TYPES,
                       VDJ_J_FEATURE_TYPES,
                       VDJ_C_FEATURE_TYPES]

# Max possible CDR3 length (in bases)
# See https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4998859/
VDJ_MAX_CDR3_LEN = 80
VDJ_MIN_CDR3_LEN = 26

# Max clonotypes for plotting frequencies
VDJ_MAX_OUTPUT_CLONOTYPES = 10

VDJ_CLONOTYPE_TYPES = ['raw', 'inferred']

VDJ_ASM_MAX_Q = 60
VDJ_QUAL_OFFSET = 33

START_CODONS = ['ATG']
STOP_CODONS = ['TAG', 'TAA', 'TGA']
CODON_TO_AA = {
    'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G',
    'TGG':'W',
    'TGT':'C', 'TGC':'C',
    'TTT':'F', 'TTC':'F',
    'TTA':'L', 'TTG':'L', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
    'ATT':'I', 'ATC':'I', 'ATA':'I',
    'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
    'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S', 'AGT':'S', 'AGC':'S',
    'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
    'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
    'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
    'TAT':'Y', 'TAC':'Y',
    'CAT':'H', 'CAC':'H',
    'CAA':'Q', 'CAG':'Q',
    'AAT':'N', 'AAC':'N',
    'AAA':'K', 'AAG':'K',
    'GAT':'D', 'GAC':'D',
    'GAA':'E', 'GAG':'E',
    'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R','AGA':'R', 'AGG':'R',
    'ATG':'M',
    # Stop codons
    'TAG':'*', 'TAA':'*', 'TGA':'*'}

CODE_PATH = os.path.dirname(os.path.abspath(__file__))
TEST_FILE_IN_DIR  = os.path.join(CODE_PATH, 'test_files', 'inputs')
TEST_FILE_OUT_DIR = os.path.join(CODE_PATH, 'test_files', 'outputs')

MEM_GB_PER_ANNOTATIONS_JSON_GB = 20

REFERENCE_FASTA_PATH = 'fasta/regions.fa'
REFERENCE_TYPE = 'V(D)J Reference'
REFERENCE_METRIC_PREFIX = 'vdj_reference_'
