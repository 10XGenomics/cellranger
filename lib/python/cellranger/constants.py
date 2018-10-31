#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
import collections
import os
from cellranger.library_constants import GENE_EXPRESSION_LIBRARY_TYPE

Interval = collections.namedtuple('Interval', ['chrom', 'start', 'end', 'length', 'strand'])
Transcript = collections.namedtuple('Transcript', ['gene', 'length', 'gc_content', 'intervals'])
Gene = collections.namedtuple('Gene', ['id', 'name', 'length', 'gc_content', 'intervals'])
Primer = collections.namedtuple('Primer', ['name', 'seq'])

ReadDef = collections.namedtuple('ReadDef', 'read_type, offset, length')

STR_DTYPE_CHAR = 'S'
UNICODE_DTYPE_CHAR = 'U'
STRING_DTYPE_CHARS = (STR_DTYPE_CHAR, UNICODE_DTYPE_CHAR)

BAM_FILE_STREAM = '-' # filename telling samtools and pysam to read stdin / write stdout

NUM_CHECK_BARCODES_FOR_ORIENTATION = 1000
REVCOMP_BARCODE_THRESHOLD = 0.75

CELLRANGER_LIB_PATH = os.path.dirname(os.path.abspath(__file__))
BARCODE_WHITELIST_PATH = os.path.join(CELLRANGER_LIB_PATH, 'barcodes')
BARCODE_WHITELIST_TRANSLATE_PATH = os.path.join(BARCODE_WHITELIST_PATH, 'translation')

NO_INPUT_FASTQS_MESSAGE = "No input FASTQs were found with the requested sample indices."

NUM_HITS_TAG = 'NH'
MULTIMAPPER_TAG = 'MM'
ANTISENSE_TAG = 'AN'
RAW_BARCODE_TAG = 'CR'
PROCESSED_BARCODE_TAG = 'CB'
RAW_BARCODE_QUAL_TAG = 'CY'
RAW_UMI_TAG = 'UR'
PROCESSED_UMI_TAG = 'UB'
UMI_QUAL_TAG = 'UY'
TRANSCRIPTS_TAG = 'TX'
GENE_IDS_TAG = 'GX'
GENE_NAMES_TAG = 'GN'
MAPPING_REGION_TAG = 'RE'
RAW_FEATURE_BARCODE_TAG = 'fr'
PROCESSED_FEATURE_BARCODE_TAG = 'fb'
FEATURE_BARCODE_QUAL_TAG = 'fq'
FEATURE_IDS_TAG = 'fx'
LIBRARY_INDEX_TAG = 'li'
EXTRA_FLAGS_TAG = 'xf'
EXTRA_FLAGS_CONF_MAPPED_TXOME   =  1
EXTRA_FLAGS_LOW_SUPPORT_UMI     =  2
EXTRA_FLAGS_GENE_DISCORDANT     =  4
EXTRA_FLAGS_UMI_COUNT           =  8
EXTRA_FLAGS_CONF_MAPPED_FEATURE = 16
GEM_GROUP_SEP = '-'

MATE_RESCUE_TAG = 'MR'

'''
Illumina read types and their internal names.
'''
FASTQ_READ_TYPES = collections.OrderedDict([
    ('R1', 'Read 1'),
    ('R2', 'Read 2'),
    ('I1', 'Index read 1'), # aka I7
    ('I2', 'Index read 2'), # aka I5
])

'''
Cell Ranger read types:
    all: all reads
    mapped: defined per region
    conf_mapped: defined per region
    conf_mapped_barcoded: read must be conf_mapped and have valid BC if BCs are used
    conf_mapped_deduped: read must be conf_mapped_barcoded, non-duplicate and have valid UMI if UMIs are used
'''
ALL_READ_TYPE = 'all'
MAPPED_READ_TYPE = 'mapped'
CONF_MAPPED_READ_TYPE = 'conf_mapped'
CONF_MAPPED_BC_READ_TYPE = 'conf_mapped_barcoded'
CONF_MAPPED_DEDUPED_READ_TYPE = 'conf_mapped_deduped_barcoded'
READ_TYPES = [
    ALL_READ_TYPE,
    MAPPED_READ_TYPE,
    CONF_MAPPED_READ_TYPE,
    CONF_MAPPED_BC_READ_TYPE,
    CONF_MAPPED_DEDUPED_READ_TYPE,
]

'''
Molecule types:
    insert: a distinct sequenced read.
    fragment: a called siPCR fragment molecule; a distinct (BC, UMI, gene, pos, strand) tuple
    cDNA: a called cDNA molecule; a distinct (BC, UMI, gene) tuple
    cDNA_candidate: a candidate cDNA molecule; a distinct (BC, UMI) tuple
'''
INSERT_MOLECULE_TYPE = 'insert'
FRAGMENT_MOLECULE_TYPE = 'fragment'
CDNA_MOLECULE_TYPE = 'cdna'
CDNA_MOLECULE_CANDIDATE_TYPE = 'cdna_candidate'
MOLECULE_TYPES = [
    INSERT_MOLECULE_TYPE,
    FRAGMENT_MOLECULE_TYPE,
    CDNA_MOLECULE_TYPE,
    CDNA_MOLECULE_CANDIDATE_TYPE,
]

'''
Barcode types:
   all_bcs: all barcodes
   filtered_bcs: cell-containing barcodes
'''
ALL_BARCODES = 'all_bcs'
FILTERED_BARCODES = 'filtered_bcs'
BARCODE_TYPES = [
    ALL_BARCODES,
    FILTERED_BARCODES,
]

'''
Regions:
    transcriptome:
      mapped: read must be mapped to sense strand of at least one transcript
      conf_mapped: read must be confidently mapped along the sense strand to a unique gene
    genome:
      mapped: read must be mapped to at least one chromosome
      conf_mapped: read must be confidently mapped to at least one chromosome
    exonic:
      mapped: read must be mapped to at least one exonic region
      conf_mapped: read must be confidently mapped to at least one exonic region
    intergenic:
      mapped: read must be mapped to at least one intergenic region
      conf_mapped: read must be confidently mapped to at least one intergenic region
    intronic:
      mapped: read must be mapped to at least one intronic region
      conf_mapped: read must be confidently mapped to at least one intronic region
'''
TRANSCRIPTOME_REGION = 'transcriptome'
GENOME_REGION = 'genome'
EXONIC_REGION = 'exonic'
INTERGENIC_REGION = 'intergenic'
INTRONIC_REGION = 'intronic'
REGIONS = [
    TRANSCRIPTOME_REGION,
    GENOME_REGION,
    EXONIC_REGION,
    INTERGENIC_REGION,
    INTRONIC_REGION,
]

REGION_TAG_MAP = {
    'E': EXONIC_REGION,
    'N': INTRONIC_REGION,
    'I': INTERGENIC_REGION,
}

'''
Duplicate types:
    cdna_pcr_uncorrected: determined by gene ID, barcode, strand, uncorrected umi
    cdna_pcr: determined by gene ID, barcode, strand, corrected umi
    si_pcr: determined by gene ID, pos, barcode, strand, corrected umi
'''
CDNA_PCR_UNCORRECTED_DUPE_TYPE = 'cdna_pcr_uncorrected'
CDNA_PCR_DUPE_TYPE = 'cdna_pcr'
SI_PCR_DUPE_TYPE = 'si_pcr'
DUPE_TYPES = [
    CDNA_PCR_UNCORRECTED_DUPE_TYPE,
    CDNA_PCR_DUPE_TYPE,
    SI_PCR_DUPE_TYPE,
]

'''
Subsampling types and target depths
'''

# Fixed depth targets in addition to quantile-based targets
SUBSAMPLE_READS_PER_CELL = [5e3, 20e3, 50e3]

RAW_SUBSAMPLE_TYPE = 'raw_rpc'
MAPPED_SUBSAMPLE_TYPE = 'conf_mapped_barcoded_filtered_bc_rpc'

ALL_SUBSAMPLE_TYPES = [
    RAW_SUBSAMPLE_TYPE,
    MAPPED_SUBSAMPLE_TYPE,
]


''' UMI properties '''
LOW_MIN_QUAL_UMI_FILTER = 'low_min_qual'
HAS_N_UMI_FILTER = 'has_n'
HOMOPOLYMER_UMI_FILTER = 'homopolymer'
PRIMER_UMI_FILTER = 'primer'
POLYT_UMI_FILTER = 'polyt_suffix'

UMI_PROPERTIES = [
    LOW_MIN_QUAL_UMI_FILTER,
    HAS_N_UMI_FILTER,
    HOMOPOLYMER_UMI_FILTER,
    PRIMER_UMI_FILTER,
    POLYT_UMI_FILTER,
]

''' Barcode properties '''
MISS_WHITELIST_BARCODE_FILTER = 'miss_whitelist'
HAS_N_BARCODE_FILTER = HAS_N_UMI_FILTER
HOMOPOLYMER_BARCODE_FILTER = HOMOPOLYMER_UMI_FILTER
LOW_MIN_QUAL_BARCODE_FILTER = LOW_MIN_QUAL_UMI_FILTER
BARCODE_MIN_QUAL_THRESHOLD = 10

BARCODE_PROPERTIES = [
    MISS_WHITELIST_BARCODE_FILTER,
    LOW_MIN_QUAL_BARCODE_FILTER,
    HAS_N_BARCODE_FILTER,
    HOMOPOLYMER_BARCODE_FILTER,
]


FORWARD_STRAND = '+'
REVERSE_STRAND = '-'
STRANDS = [FORWARD_STRAND, REVERSE_STRAND]

THREE_PRIME = 'three_prime'
FIVE_PRIME = 'five_prime'

READ_POSITION_CUTOFFS = range(0, 5100, 100)
INSERT_SIZE_CUTOFFS = range(0,1550,50)
PER_DUPE_GROUP_MAX = 100
DEFAULT_TOP_BARCODE_CUTOFF = 1000
MIN_COUNTS_PER_BARCODE = 2
MIN_COUNTS_PER_GENE = 1
ORDMAG_RECOVERED_CELLS_QUANTILE = 0.99
ORDMAG_NUM_BOOTSTRAP_SAMPLES = 100
FILTER_BARCODES_MAX_RECOVERED_CELLS_MULTIPLE = 6
DEFAULT_RECOVERED_CELLS_PER_GEM_GROUP = 3000

MINIMUM_TRIMMED_READ_LENGTH_PREFLIGHT = 16

READ_PREFIX_LENGTH = 10
TOP_N = 5
# Downsample reads when tracking top sequences to constrain mem usage
TOP_RAW_SEQUENCE_SAMPLE_RATE = 0.01
TOP_PROCESSED_SEQUENCE_SAMPLE_RATE = 0.1
HOMOPOLYMER_LENGTH = 15

DEFAULT_REPORT_TYPE = 'summary'
H5_BC_SEQUENCE_COL = 'bc_sequence'


CELLRANGER_VERSION_KEY = 'cellranger_version'

REFERENCE_METADATA_FILE = 'reference.json'
REFERENCE_STAR_PATH = 'star'
REFERENCE_FASTA_PATH = 'fasta/genome.fa'
REFERENCE_GENES_GTF_PATH = 'genes/genes.gtf'
REFERENCE_GENES_INDEX_PATH = 'pickle/genes.pickle'
REFERENCE_GENOMES_KEY = 'genomes'
REFERENCE_MEM_GB_KEY = 'mem_gb'
REFERENCE_NUM_THREADS_KEY = 'threads'

# Ref metadata keys used by GEX and VDJ
REFERENCE_FASTA_HASH_KEY = 'fasta_hash'
REFERENCE_GTF_HASH_KEY = 'gtf_hash'
REFERENCE_INPUT_FASTA_KEY = 'input_fasta_files'
REFERENCE_INPUT_GTF_KEY = 'input_gtf_files'
REFERENCE_MKREF_VERSION_KEY = 'mkref_version'
REFERENCE_VERSION_KEY = 'version'
REFERENCE_TYPE_KEY = 'type'
REFERENCE_TYPE = 'Transcriptome'
REFERENCE_METRIC_PREFIX = 'reference_'
REFERENCE_METADATA_KEYS = [
    REFERENCE_FASTA_HASH_KEY,
    REFERENCE_GTF_HASH_KEY,
    REFERENCE_INPUT_FASTA_KEY,
    REFERENCE_INPUT_GTF_KEY,
    REFERENCE_MKREF_VERSION_KEY,
    REFERENCE_VERSION_KEY,
    REFERENCE_TYPE_KEY,      # New in CR 2.0.0
    REFERENCE_GENOMES_KEY,   # New in CR 2.0.0
]

BAM_CHUNK_SIZE_GB = 0.5
MAX_BAM_CHUNKS = 256

'''
NOTE: MEM_GB_PER_THREAD is only used for a few stages where we've encountered memory oversubscription issues
on clusters without memory reservations. As of 3/15/2017, it's used by:
- RUN_PCA
- RUN_DIFFERENTIAL_EXPRESSION
- RUN_GRAPH_CLUSTERING
'''
MEM_GB_PER_THREAD = 8

COUNT_GENES_MAX_MEM_GB = 64
BYTES_PER_STR_INT_DICT_ENTRY = 180 # Empirical testing w/ str:int dicts and pympl.asizeof; len(k) = 10,14; *1.25 safety
NUM_BARCODES_PER_MEM_GB = 500000
NUM_MOLECULE_INFO_ENTRIES_PER_CHUNK = 40000000

SUPPORTED_ALIGNERS = ['star']
STAR_REQUIRED_FILES = [
    'chrLength.txt',
    'chrNameLength.txt',
    'chrName.txt',
    'chrStart.txt',
    'Genome',
    'genomeParameters.txt',
    'SA',
    'SAindex',
]
STAR_DEFAULT_HIGH_CONF_MAPQ = 255

MATRIX_REPORT_READ_TYPES = [CONF_MAPPED_BC_READ_TYPE, CONF_MAPPED_DEDUPED_READ_TYPE]
MATRIX_USE_MATRIX_FOR_READ_TYPE = [CONF_MAPPED_DEDUPED_READ_TYPE]

UMI_POLYT_SUFFIX_LENGTH = 5

NORM_MODE_MAPPED = 'mapped'
NORM_MODE_NONE = 'none'
NORM_MODES = [NORM_MODE_MAPPED, NORM_MODE_NONE]

AGG_ID_FIELD = 'library_id'
AGG_H5_FIELD = 'molecule_h5'
AGG_BATCH_FIELD = 'batch'
# to distinguish from user-defined keys
AGG_METADATA_FIELDS = [AGG_ID_FIELD, AGG_H5_FIELD, AGG_BATCH_FIELD]

MAX_INSERT_SIZE = 1000

''' Chemistry detection '''
DETECT_CHEMISTRY_MIN_FRAC_WHITELIST = 0.10
# 26bcumi + 13spacer + 11bp for gene expression
DETECT_5P_CHEMISTRY_MIN_R1_LEN_PE = 50
# 26bcumi + 13spacer + 50bp for assembler
DETECT_VDJ_CHEMISTRY_MIN_R1_LEN_PE = 90
DETECT_CHEMISTRY_INITIAL_READS = 100000


"""PIPELINE NAMES"""
PIPELINE_VDJ = 'vdj'
PIPELINE_COUNT = 'count'

PACKAGE_VERSION_CMDS = [
    {
        'name': 'mrc',
        'cmd' : 'mrc --version',
    },
    {
        'name': 'mrp',
        'cmd' : 'mrp --version',
    },
    {
        'name': 'Anaconda',
        'cmd' : 'python --version 2>&1 | cat ',
    },
    {
        'name': 'numpy',
        'cmd' : 'python -c "import numpy; print numpy.__version__"'
    },
    {
        'name': 'scipy',
        'cmd' : 'python -c "import scipy; print scipy.__version__"'
    },
    {
        'name': 'pysam',
        'cmd' : 'python -c "import pysam; print pysam.__version__"'
    },
    {
        'name': 'h5py',
        'cmd' : 'python -c "import h5py; print h5py.__version__"'
    },
    {
        'name': 'pandas',
        'cmd' : 'python -c "import pandas; print pandas.__version__"'
    },
    {
        'name': 'STAR',
        'cmd': 'STAR --version',
    },
    {
        'name': 'samtools',
        'cmd': 'samtools --version',
    }
]

# Pysam numeric codes to meaningful categories
cigar_numeric_to_category_map = {0: 'M', 1: 'I', 2: 'D', 3: 'N', 4: 'S', 5: 'H', 6: 'P', 7: '=', 8: 'X'}

CODE_PATH = os.path.dirname(os.path.abspath(__file__))
TEST_FILE_IN_DIR  = os.path.join(CODE_PATH, 'test_files', 'inputs')
TEST_FILE_OUT_DIR = os.path.join(CODE_PATH, 'test_files', 'outputs')

BARCODE_CSV_COLNAME = 'Barcode'
GENE_ID_CSV_COLNAME = 'Gene'

ALIGN_LIBRARY_TYPES = [GENE_EXPRESSION_LIBRARY_TYPE]

# Constants for barcode compatability checking
BARCODE_COMPATIBILITY_CUTOFF = 0.1
NUM_READS_TO_CHECK_BARCODE = 400000

# Field names for the `--libraries` csv input file
LIBRARIES_CSV_FIELDS = ['fastqs', 'sample', 'library_type']
