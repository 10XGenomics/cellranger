#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
import collections
import csv
import h5py
import hashlib
import itertools
import json
import numpy as np
import os
import pandas as pd
import errno
import random
import re
import shutil
import subprocess
import sys
import tables
import tenkit.bam as tk_bam
import tenkit.fasta as tk_fasta
import tenkit.log_subprocess as tk_subproc
import tenkit.seq as tk_seq
import tenkit.stats as tk_stats
import tenkit.constants as tk_constants
import cellranger.constants as cr_constants

def get_version():
    # NOTE: this makes assumptions about the directory structure
    script_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'bin')
    version_fn = os.path.join(script_dir, '..', '.version')
    if os.path.exists(version_fn):
        with open(version_fn, 'r') as f:
            output = f.read()
    else:
        output = tk_subproc.check_output(['git', 'describe', '--tags', '--always', '--dirty'], cwd=script_dir)
    return output.strip()

def fixpath(path):
    return os.path.abspath(os.path.expandvars(os.path.expanduser(path)))

def get_input_path(oldpath, is_dir=False):
    path = fixpath(oldpath)
    if not os.path.exists(path):
        sys.exit("Input file does not exist: %s" % path)
    if is_dir:
        if not os.path.isdir(path):
            sys.exit("Please provide a directory, not a file: %s" % path)
    else:
        if not os.path.isfile(path):
            sys.exit("Please provide a file, not a directory: %s" % path)
    return path

def get_input_paths(oldpaths):
    paths = []
    for oldpath in oldpaths:
        paths.append(get_input_path(oldpath))
    return paths

def get_output_path(oldpath):
    path = fixpath(oldpath)
    dirname = os.path.dirname(path)
    if not os.path.exists(dirname):
        sys.exit("Output directory does not exist: %s" % dirname)
    if not os.path.isdir(dirname):
        sys.exit("Please provide a directory, not a file: %s" % dirname)
    return path

def _load_reference_metadata_file(reference_path):
    reference_metadata_file = os.path.join(reference_path, cr_constants.REFERENCE_METADATA_FILE)
    with open(reference_metadata_file, 'r') as f:
        return json.load(f)

def get_reference_star_path(reference_path):
    return os.path.join(reference_path, cr_constants.REFERENCE_STAR_PATH)

def get_reference_genes_gtf(reference_path):
    return os.path.join(reference_path, cr_constants.REFERENCE_GENES_GTF_PATH)

def get_reference_genes_index(reference_path):
    return os.path.join(reference_path, cr_constants.REFERENCE_GENES_INDEX_PATH)

def get_reference_genome_fasta(reference_path):
    return os.path.join(reference_path, cr_constants.REFERENCE_FASTA_PATH)

def is_reference_star_indexed(fasta_path):
    for filename in cr_constants.STAR_REQUIRED_FILES:
        p = os.path.join(cr_constants.REFERENCE_STAR_PATH, filename)
        if not os.path.exists(p):
            return False
    return True

def get_fasta_hash(reference_path):
    data = _load_reference_metadata_file(reference_path)
    return data[cr_constants.REFERENCE_FASTA_HASH_KEY]

def get_gtf_hash(reference_path):
    data = _load_reference_metadata_file(reference_path)
    return data[cr_constants.REFERENCE_GTF_HASH_KEY]

def get_reference_genomes(reference_path):
    data = _load_reference_metadata_file(reference_path)
    return data[cr_constants.REFERENCE_GENOMES_KEY]

def get_reference_num_threads_request(reference_path):
    data = _load_reference_metadata_file(reference_path)
    return data[cr_constants.REFERENCE_NUM_THREADS_KEY]

def get_reference_mem_gb_request(reference_path):
    data = _load_reference_metadata_file(reference_path)
    return data[cr_constants.REFERENCE_MEM_GB_KEY]

def get_mkref_version(reference_path):
    data = _load_reference_metadata_file(reference_path)
    return data[cr_constants.REFERENCE_MKREF_VERSION_KEY]

def barcode_sort_key(read, squash_unbarcoded=False):
    formatted_bc = get_read_barcode(read)
    if squash_unbarcoded and formatted_bc is None:
        return None
    (bc, gg) = split_barcode_seq(formatted_bc)
    return gg, bc, get_read_gene_ids(read)

def pos_sort_key(read):
    return read.tid, read.pos

def cdna_pcr_dupe_func(read):
    return get_read_barcode(read), get_read_gene_ids(read)

def si_pcr_dupe_func(read):
    return cdna_pcr_dupe_func(read), read.tid, read.pos, read.mrnm, read.mpos

def format_barcode_seq(barcode, gem_group=None):
    if gem_group is not None:
        barcode += '-' + str(gem_group)
    return barcode

def format_barcode_seqs(barcode_seqs, gem_groups):
    if gem_groups is None:
        return barcode_seqs
    new_barcode_seqs = []
    unique_gem_groups = sorted(list(set(gem_groups)))
    for gg in unique_gem_groups:
        new_barcode_seqs += [format_barcode_seq(bc, gg) for bc in barcode_seqs]
    return new_barcode_seqs

def split_barcode_seq(barcode):
    if barcode is None:
        return None, None

    barcode_parts = barcode.split('-')

    barcode = barcode_parts[0]
    if len(barcode_parts) > 1:
        gem_group = int(barcode_parts[1])
    else:
        gem_group = None

    return barcode, gem_group

def is_barcode_corrected(raw_bc_seq, processed_bc_seq):
    if processed_bc_seq is None:
        return False

    bc_seq, _ = split_barcode_seq(processed_bc_seq)
    return bc_seq != raw_bc_seq

def is_umi_corrected(raw_umi_seq, processed_umi_seq):
    if processed_umi_seq is None:
        return False

    return raw_umi_seq != processed_umi_seq

def split_genes_by_genomes(genes, genomes):
    """ Returns a list of lists [genome1, genome2, ...]
    where genome1 = [gene1,gene2,...].
    Args:
      genes - list of Gene tuples
      genomes - list of genome names, e.g. ['hg19', 'mm10']
    """
    assert len(genomes) > 0

    if len(genomes) == 1:
        return [genes]

    d = collections.defaultdict(list)
    for gene in genes:
        genome = get_genome_from_str(gene.id, genomes)
        d[genome].append(gene)

    genes_per_genome = []
    for genome in genomes:
        genes_per_genome.append(d[genome])
    return genes_per_genome

def get_genome_from_str(s, genomes):
    assert len(genomes) > 0

    if s is None:
        return None

    if len(genomes) == 1:
        return genomes[0]

    for genome in genomes:
        if s.startswith(genome):
            return genome

    raise Exception('%s does not have valid associated genome' % s)

def remove_genome_from_str(s, genomes):
    assert len(genomes) > 0

    if s is None:
        return None

    if len(genomes) == 1:
        return s

    for genome in genomes:
        if s.startswith(genome):
            # Strip genome and subsequent underscore
            return s[(1+len(genome)):]

    raise Exception('%s does not have valid associated genome' % s)

def get_genome_from_read(read, chroms, genomes):
    assert len(genomes) > 0

    if len(genomes) == 1:
        return genomes[0]

    if read.is_unmapped:
        return None

    return get_genome_from_str(chroms[read.tid], genomes)

def get_read_num_hits(read):
    return _get_read_tag(read, cr_constants.NUM_HITS_TAG)

def get_read_barcode(read):
    return _get_read_tag(read, cr_constants.PROCESSED_BARCODE_TAG)

def get_read_raw_barcode(read):
    return _get_read_tag(read, cr_constants.RAW_BARCODE_TAG)

def get_read_barcode_qual(read):
    return _get_read_tag(read, cr_constants.RAW_BARCODE_QUAL_TAG)

def get_read_umi_qual(read):
    return _get_read_tag(read, cr_constants.UMI_QUAL_TAG)

def get_read_umi(read):
    return _get_read_tag(read, cr_constants.PROCESSED_UMI_TAG)

def get_read_raw_umi(read):
    return _get_read_tag(read, cr_constants.RAW_UMI_TAG)

def get_read_gene_ids(read):
    s =_get_read_tag(read, cr_constants.GENE_IDS_TAG)
    if s is None:
        return None

    return tuple(s.split(';'))

def get_read_gene_names(read):
    s =_get_read_tag(read, cr_constants.GENE_NAMES_TAG)
    if s is None:
        return None

    return tuple(s.split(';'))

def get_genes_from_transcripts(transcript_alignments, gene_index):
    read_genes = collections.OrderedDict() # TODO replace OrderedDict
    for alignment in transcript_alignments:
        gene = gene_index.get_gene_from_transcript(alignment.transcript)
        assert gene is not None
        if gene.id not in read_genes:
            read_genes[gene.id] = gene.name
    return read_genes

def get_read_transcripts_iter(read):
    s = _get_read_tag(read, cr_constants.TRANSCRIPTS_TAG)
    if s is None:
        return

    for x in s.split(';'):
        if len(x) == 0:
            continue

        parts = x.split(',')
        assert len(parts) == 3

        chrom = parts[0] if parts[0] != '' else None
        if parts[1] != '':
            strand = parts[1][0]
            pos = int(parts[1][1:])
        else:
            strand = None
            pos = None
        cigarstring = parts[2] if parts[2] != '' else None

        assert strand in cr_constants.STRANDS or strand is None

        yield chrom, strand, pos, cigarstring

def get_mapping_region(read):
    region_tag = _get_read_tag(read, 'RE')
    return cr_constants.REGION_TAG_MAP.get(region_tag, None)

def make_annotation_tags(tx, gx, gn):
    return [(cr_constants.TRANSCRIPTS_TAG, tx), (cr_constants.GENE_IDS_TAG, gx), (cr_constants.GENE_NAMES_TAG, gn)]

def iter_read_pairs(reads_iter, paired_end):
    """
    Iterate through mate pairs. In the paired end case, the number of R1:R2 alignments
    should be either N:N or 1:N / N:1 (if one end is unmapped, in which case we treat it separately).
    """
    if not paired_end:
        for read in reads_iter:
            yield read, None
    else:
        r1s, r2s = [], []
        for read in reads_iter:
            if read.is_read2:
                r2s.append(read)
            else:
                r1s.append(read)
        if len(r1s) == len(r2s):
            for r1, r2 in itertools.izip(r1s, r2s):
                yield r1, r2
        elif len(r1s) == 1:
            for r2 in r2s:
                yield None, r2
            yield r1s[0], None
        elif len(r2s) == 1:
            for r1 in r1s:
                yield r1, None
            yield None, r2s[0]
        else:
            raise Exception("Inconsistent number of R1:R2 alignments (%d:%d)" % (len(r1s), len(r2s)))

def iter_by_qname(in_genome_bam, in_trimmed_bam):
    # Iterate through multiple BAMs by qname simultaneously
    # Assume the trimmed-read-bam has every qname in the genome bam, in the same order.

    genome_bam_iter = itertools.groupby(in_genome_bam, key=lambda read: read.qname)

    if in_trimmed_bam is None:
        trimmed_bam_iter = iter(())
    else:
        trimmed_bam_iter = itertools.groupby(in_trimmed_bam, key=lambda read: read.qname)

    for (genome_qname, genome_reads), trimmed_tuple in itertools.izip_longest(genome_bam_iter,
                                                                              trimmed_bam_iter):
        trimmed_qname, trimmed_reads = trimmed_tuple or (None, [])
        genome_reads = list(genome_reads)
        trimmed_reads = list(trimmed_reads)

        assert (in_trimmed_bam is None) or trimmed_qname == genome_qname
        yield (genome_qname, genome_reads, trimmed_reads)

def set_read_tags(read, tags):
    read_tags = read.tags
    read_tags.extend(tags)
    read.tags = read_tags

def _get_read_tag(read, tag):
    try:
        r = read.opt(tag)
        if r == '':
            r = None
        return r
    except KeyError:
        return None

def get_fastq_read1(r1, r2, reads_interleaved):
    if reads_interleaved:
        name1, seq1, qual1, _, _, _ = r1
    else:
        name1, seq1, qual1 = r1
    return name1, seq1, qual1

def get_fastq_read2(r1, r2, reads_interleaved):
    if reads_interleaved:
        _, _, _, name2, seq2, qual2 = r1
    else:
        name2, seq2, qual2 = r2
    return name2, seq2, qual2


def get_primers_from_dicts(dicts):
    return [cr_constants.Primer(d['name'], d['seq']) for d in dicts] if dicts is not None else []

def get_comp_primer_name(primer):
    return primer.name + '_comp'

def get_rev_primer_name(primer):
    return primer.name + '_rev'

def get_rev_comp_primer_name(primer):
    return primer.name + '_revcomp'

def get_primer_orientation_names(primer):
    if is_homopolymer_seq(primer.seq):
        return [primer.name]
    else:
        return [primer.name,
                get_comp_primer_name(primer),
                get_rev_primer_name(primer),
                get_rev_comp_primer_name(primer)]

def get_primer_orientation_seqs(primer):
    if is_homopolymer_seq(primer.seq):
        return [primer.seq]
    else:
        return [primer.seq,
                get_comp(primer.seq),
                get_rev(primer.seq),
                tk_seq.get_rev_comp(primer.seq)]

def is_homopolymer_seq(seq):
    for nuc in tk_seq.NUCS:
        if seq == nuc * len(seq):
            return True
    return False

def is_unambiguous_nuc_seq(seq):
    return not bool(re.search(r'[^ACGT]', seq))

def is_transcriptome_read_sense(read, strandedness):
    """ - read: transcriptome-aligned read
        - strandedness: library strandedness
    """
    if strandedness is None:
        return True

    if read.is_reverse ^ read.is_read2:
        return strandedness == cr_constants.REVERSE_STRAND
    else:
        return strandedness == cr_constants.FORWARD_STRAND

def is_read_dupe_candidate(read, high_conf_mapq, use_corrected_umi=True):
    if use_corrected_umi:
        umi = get_read_umi(read)
    else:
        umi = get_read_raw_umi(read)

    return not read.is_secondary and umi and get_read_barcode(read) and \
        is_read_conf_mapped_to_transcriptome(read, high_conf_mapq)

def is_read_conf_mapped_to_transcriptome(read, high_conf_mapq):
    if read.is_unmapped:
        return False
    elif read.mapq < high_conf_mapq:
        return False
    else:
        gene_ids = get_read_gene_ids(read)
        return gene_ids is not None and len(gene_ids) == 1

def is_read_conf_mapped_to_transcriptome_barcoded(read, high_conf_mapq):
    return is_read_conf_mapped_to_transcriptome(read, high_conf_mapq) and \
        not read.is_secondary and get_read_barcode(read)

def is_read_conf_mapped_to_transcriptome_barcoded_deduped(read, high_conf_mapq):
    return is_read_conf_mapped_to_transcriptome_barcoded(read, high_conf_mapq) and \
        not read.is_duplicate

def is_read_conf_mapped_to_transcriptome_barcoded_deduped_with_umi(read, high_conf_mapq):
    return is_read_conf_mapped_to_transcriptome_barcoded_deduped(read, high_conf_mapq) and \
        get_read_umi(read)

def get_hamming_distance(kmer1, kmer2):
    kmer1_len, kmer2_len = len(kmer1), len(kmer2)

    hamming_distance = abs(kmer1_len - kmer2_len)
    for i in xrange(min(kmer1_len, kmer2_len)):
        if kmer1[i] != kmer2[i]:
            hamming_distance += 1
    return hamming_distance

def get_kmers_hamming_distance(kmers):
    min_hamming_distance = None
    for i, kmer1 in enumerate(kmers):
        for kmer2 in kmers[i+1:]:
            hamming_distance = get_hamming_distance(kmer1, kmer2)
            min_hamming_distance = min(min_hamming_distance, hamming_distance) \
                                   if min_hamming_distance is not None else hamming_distance
            if min_hamming_distance == 0:
                return min_hamming_distance
    return min_hamming_distance

def load_barcode_tsv(filename):
    barcodes = [x.strip() for x in open_maybe_gzip(filename, 'r') if not ('#' in x)]
    if len(barcodes) != len(set(barcodes)):
        raise Exception('Duplicates found in barcode whitelist: %s' % filename)
    return barcodes

def get_barcode_whitelist_path(filename):
    # Look for exact path, .txt.gz, or .txt
    if filename is None:
        return None
    elif os.path.exists(filename):
        return filename
    else:
        gz = os.path.join(cr_constants.BARCODE_WHITELIST_PATH, filename + '.txt.gz')
        if os.path.exists(gz):
            return gz

        txt = os.path.join(cr_constants.BARCODE_WHITELIST_PATH, filename + '.txt')
        return txt

def load_barcode_whitelist(filename):
    path = get_barcode_whitelist_path(filename)

    if path is None:
        return None

    if not os.path.isfile(path):
        raise NameError('Unable to find barcode whitelist: %s' % path)

    return load_barcode_tsv(path)


def load_barcode_summary(barcode_summary):
    if barcode_summary:
        with h5py.File(barcode_summary) as f:
            return list(f[cr_constants.H5_BC_SEQUENCE_COL])
    return None

def get_num_barcodes_from_barcode_summary(barcode_summary):
    if barcode_summary:
        with h5py.File(barcode_summary) as f:
            return f[cr_constants.H5_BC_SEQUENCE_COL].size

def write_h5(filename, data):
    with h5py.File(filename, 'w') as f:
        for key, value in data.iteritems():
            f[key] = value

def concatenate_h5(input_files, output_file):
    with tables.open_file(output_file, mode = 'w') as fout:
        dsets = {}
        # init datasets using the first input
        if len(input_files) > 0:
            with tables.open_file(input_files[0], mode = 'r') as fin:
                for node in fin.walk_nodes('/', 'Array'):
                    atom = tables.Atom.from_dtype(np.dtype(node.dtype))
                    dsets[node.name] = fout.create_earray('/', node.name, atom, (0,))
        # copy the data
        for input_file in input_files:
            with tables.open_file(input_file, mode = 'r') as fin:
                for (name, earray) in dsets.iteritems():
                    earray.append(fin.get_node('/', name)[:])

def compress_seq(s, bits=64):
    """ Pack a DNA sequence (no Ns!) into a 2-bit format, in a 64-bit uint """
    """ Most significant bit is set if there was an error """
    assert len(s) <= (bits/2 - 1)
    result = 0
    for nuc in s:
        if not nuc in tk_seq.NUCS_INVERSE:
            return 1L << (bits - 1)
        result = result << 2
        result = result | tk_seq.NUCS_INVERSE[nuc]
    return result

def decompress_seq(x, length, bits=64):
    x = np.uint64(x)
    assert length <= (bits/2 - 1)
    if x & (1L << (bits-1)):
        return 'N' * length
    result = bytearray(length)
    for i in xrange(length):
        result[(length-1)-i] = tk_seq.NUCS[x & np.uint64(0b11)]
        x = x >> np.uint64(2)
    return str(result)

def set_tag(read, key, old_value, new_value):
    """ Set a bam tag for a read, overwriting the previous value """
    tags = read.tags
    tags.remove((key, old_value))
    tags.append((key, new_value))
    read.tags = tags

def validate_alignment_params(aln_params):
    """ Raise exception if align params are invalid """
    keys = ['aligner', 'high_conf_mapq']
    for key in keys:
        assert key in aln_params
    assert aln_params['aligner'] in cr_constants.SUPPORTED_ALIGNERS
    assert aln_params['high_conf_mapq'] is None or aln_params['high_conf_mapq'] >= 0

def set_if_none(d,k,v):
    """ Update a dict value if it's set to None """
    if d[k] is None:
        d[k] = v

def select_alignment_params(aln_params):
    """ Select alignment parameters """
    if aln_params is None:
        return None

    result = dict(aln_params)
    if aln_params['aligner'] == 'star':
        set_if_none(result, 'high_conf_mapq', cr_constants.STAR_DEFAULT_HIGH_CONF_MAPQ)
    validate_alignment_params(result)
    return result

def build_alignment_param_metrics(align):
    """ Transform alignment param dict -> alignment metrics for output to json"""
    aln_param_metrics = {}
    for key, value in align.iteritems():
        aln_param_metrics['alignment_' + key] = value
    return aln_param_metrics

def get_high_conf_mapq(align):
    return align['high_conf_mapq']

def get_aligner(align):
    return align['aligner']

def get_thread_request_from_mem_gb(mem_gb):
    """ For systems without memory reservations, reserve multiple threads if necessary to avoid running out of memory"""
    est_threads = round(float(mem_gb) / cr_constants.MEM_GB_PER_THREAD)
    # make sure it's 1, 2, or 4
    for threads in [1, 2, 4]:
        if est_threads <= threads: return threads
    return 4

def get_mem_gb_request_from_genome_fasta(reference_path):
    in_fasta_fn = get_reference_genome_fasta(reference_path)
    genome_size_gb = float(os.path.getsize(in_fasta_fn)) / 1e9
    return max(cr_constants.MIN_MEM_GB, cr_constants.BAM_CHUNK_SIZE_GB + max(1, 2*int(genome_size_gb)))

def get_mem_gb_request_from_barcode_whitelist(barcode_whitelist_fn, gem_groups=None, use_min=True, double=False):
    barcode_whitelist = load_barcode_whitelist(barcode_whitelist_fn)

    if use_min:
        if barcode_whitelist is None:
            min_mem_gb = cr_constants.MIN_MEM_GB_NOWHITELIST
        else:
            min_mem_gb = cr_constants.MIN_MEM_GB
    else:
        min_mem_gb = 0

    if barcode_whitelist is None:
        return min_mem_gb

    if gem_groups is not None:
        num_bcs = len(barcode_whitelist)*max(gem_groups)
    else:
        num_bcs = len(barcode_whitelist)

    if double:
        return max(min_mem_gb, 2 * num_bcs / cr_constants.NUM_BARCODES_PER_MEM_GB)
    else:
        return max(min_mem_gb,     num_bcs / cr_constants.NUM_BARCODES_PER_MEM_GB)

def update_require_unique_key(dest_dict, src_dict):
    """ Update a dict w/ another dict; raise exception on duplicate keys """
    for key, value in src_dict.iteritems():
        if key in dest_dict:
            print key
        #qassert key not in dest_dict
        dest_dict[key] = value

def merge_jsons_as_dict(in_filenames):
    """ Merge a list of json files and return the result as a dictionary """
    d = {}
    for filename in in_filenames:
        if filename is None:
            continue
        with open(filename, 'r') as f:
            data = json.load(f)
            update_require_unique_key(d, data)
    return d

def get_metric_from_json(filename, key):
    with open(filename, 'r') as f:
        d = json.load(f)
    return d[key]

def format_barcode_summary_h5_key(genome, region, read_type):
    return '%s_%s_%s_reads' % (genome, region, read_type)

def downsample(rate):
    if rate is None or rate == 1.0:
        return True
    return random.random() <= rate

def numpy_groupby(values, keys):
    """ Group a collection of numpy arrays by key arrays.
        Yields (key_tuple, view_tuple) where key_tuple is the key grouped on and view_tuple is a tuple of views into the value arrays.
          values: tuple of arrays to group
          keys: tuple of sorted, numeric arrays to group by """

    if len(values) == 0:
        return
    if len(values[0]) == 0:
        return

    for key_array in keys:
        assert len(key_array) == len(keys[0])
    for value_array in values:
        assert len(value_array) == len(keys[0])

    # The indices where any of the keys differ from the previous key become group boundaries
    key_change_indices = np.logical_or.reduce(tuple(np.concatenate(([1], np.diff(key))) != 0 for key in keys))
    group_starts = np.flatnonzero(key_change_indices)
    group_ends = np.roll(group_starts, -1)
    group_ends[-1] = len(keys[0])

    for group_start, group_end in itertools.izip(group_starts, group_ends):
        yield tuple(key[group_start] for key in keys), tuple(value[group_start:group_end] for value in values)

def min_qual_below(qual, threshold):
    """ Return true if the min qual is below the threshold
         qual: a read's qual string  """
    return threshold is not None and tk_fasta.get_min_qual(qual) < threshold

def find_any_primers(seq, primer_dict):
    """ Count occurrences of primers in a query sequence
        primer_dict is of the form {name: {'seq':x, 'seq_rc':y}, ...} """
    return any(primer['seq'][0:len(seq)] in seq or \
               primer['seq_rc'][0:len(seq)] in seq \
               for primer in primer_dict.itervalues())

def is_barcode_on_whitelist(seq, whitelist):
    if whitelist is None:
        return 'N' not in seq
    else:
        return seq in whitelist

class SubprocessStream(object):
    """ Wrap a subprocess that we stream from or stream to """
    def __init__(self, *args, **kwargs):
        mode = kwargs.pop('mode', 'r')
        if mode == 'r':
            kwargs['stdout'] = subprocess.PIPE
        elif mode == 'w':
            kwargs['stdin'] = subprocess.PIPE
        else:
            raise ValueError('mode %s unsupported' % self.mode)

        kwargs['preexec_fn'] = os.setsid
        print args[0]
        sys.stdout.flush()
        self.proc = tk_subproc.Popen(*args, **kwargs)

        if mode == 'r':
            self.pipe = self.proc.stdout
        elif mode == 'w':
            self.pipe = self.proc.stdin
    def __enter__(self):
        return self
    def __iter__(self):
        return self
    def next(self):
        return self.pipe.next()
    def fileno(self):
        return self.pipe.fileno()
    def write(self, x):
        self.pipe.write(x)
    def close(self):
        self.pipe.close()
        self.proc.wait()
    def __exit__(self, tp, val, tb):
        self.close()

def open_maybe_gzip(filename, mode='r'):
    compressor = None

    if filename.endswith(cr_constants.GZIP_SUFFIX):
        compressor = 'gzip'
    elif filename.endswith(cr_constants.LZ4_SUFFIX):
        compressor = 'lz4'
    else:
        return open(filename, mode)

    assert compressor is not None

    if mode == 'r':
        return SubprocessStream([compressor, '-c', '-d', filename], mode='r')
    elif mode == 'w':
        f = open(filename, 'w')
        return SubprocessStream([compressor, '-c'], stdout=f, mode='w')
    else:
        raise ValueError("Unsupported mode for compression: %s" % mode)

def get_comp(seq):
    return str(seq).translate(tk_seq.DNA_CONVERT_TABLE)

def get_rev(seq):
    return str(seq)[::-1]

def run_command_safely(cmd, args):
    p = tk_subproc.Popen([cmd] + args, stderr=subprocess.PIPE)
    _, stderr_data = p.communicate()
    if p.returncode != 0:
        raise Exception("%s returned error code %d: %s" % (p, p.returncode, stderr_data))

class CRCalledProcessError(Exception):
    def __init__(self, msg):
        self.msg = msg
    def __str__(self):
        return self.msg

def check_completed_process(p, cmd):
    """ p   (Popen object): Subprocess
        cmd (str):          Command that was run
    """
    if p.returncode is None:
        raise CRCalledProcessError("Process did not finish: %s ." % cmd)
    elif p.returncode != 0:
        raise CRCalledProcessError("Process returned error code %d: %s ." % (p.returncode, cmd))

def mkdir(dst, allow_existing=False):
    """ Create a directory. Optionally succeed if already exists.
        Useful because transient NFS server issues may induce double creation attempts. """
    if allow_existing:
        try:
            os.mkdir(dst)
        except OSError as e:
            if e.errno == errno.EEXIST and os.path.isdir(dst):
                pass
            else:
                raise
    else:
        os.mkdir(dst)

def makedirs(dst, allow_existing=False):
    """ Create a directory recursively. Optionally succeed if already exists.
        Useful because transient NFS server issues may induce double creation attempts. """
    if allow_existing:
        try:
            os.makedirs(dst)
        except OSError as e:
            if e.errno == errno.EEXIST and os.path.isdir(dst):
                pass
            else:
                raise
    else:
        os.makedirs(dst)

def remove(f, allow_nonexisting=False):
    """ Delete a file. Allow to optionally succeed if file doesn't exist.
        Useful because transient NFS server issues may induce double deletion attempts. """
    if allow_nonexisting:
        try:
            os.remove(f)
        except OSError as e:
            if e.errno == errno.ENOENT:
                pass
            else:
                raise
    else:
        os.remove(f)

def copy(src, dst):
    """ Safely copy a file. Not platform-independent """
    run_command_safely('cp', [src, dst])

def move(src, dst):
    """ Safely move a file. Not platform-independent """
    run_command_safely('mv', [src, dst])

def copytree(src, dst, allow_existing=False):
    """ Safely recursively copy a directory. Not platform-independent """
    makedirs(dst, allow_existing=allow_existing)

    for name in os.listdir(src):
        srcname = os.path.join(src, name)
        dstname = os.path.join(dst, name)

        if os.path.isdir(srcname):
            copytree(srcname, dstname)
        else:
            copy(srcname, dstname)

def concatenate_files(out_path, in_paths, mode=''):
    with open(out_path, 'w' + mode) as out_file:
        for in_path in in_paths:
            with open(in_path, 'r' + mode) as in_file:
                shutil.copyfileobj(in_file, out_file)

def concatenate_headered_files(out_path, in_paths, mode=''):
    """ Concatenate files, taking the first line of the first file
        and skipping the first line for subsequent files.
        Asserts that all header lines are equal. """
    with open(out_path, 'w' + mode) as out_file:
        if len(in_paths) > 0:
            # Write first file
            with open(in_paths[0], 'r' + mode) as in_file:
                header = in_file.readline()
                out_file.write(header)
                shutil.copyfileobj(in_file, out_file)

        # Write remaining files
        for in_path in in_paths[1:]:
            with open(in_path, 'r' + mode) as in_file:
                this_header = in_file.readline()
                assert this_header == header
                shutil.copyfileobj(in_file, out_file)

def compute_hash_of_file(filename, block_size_bytes=2**20):
    digest = hashlib.sha1()
    with open(filename, 'rb') as f:
        for chunk in iter(lambda: f.read(block_size_bytes), b''):
            digest.update(chunk)
    return digest.hexdigest()

def flip_strand(strand):
    if strand == '+':
        return '-'
    elif strand == '-':
        return '+'
    else:
        return strand

def get_seqs(l):
    if l == 1:
        return tk_seq.NUCS
    old_seqs = get_seqs(l-1)
    new_seqs = []
    for old_seq in old_seqs:
        for base in tk_seq.NUCS:
            new_seqs.append(old_seq + base)
    return new_seqs

def load_barcode_dist(filename, barcode_whitelist, gem_group, proportions=True):
    """ Load barcode count distribution from a json file """
    # Input barcode whitelist must be an ordered type;
    # safeguard against it going out of sync with the distribution file
    assert barcode_whitelist is None or isinstance(barcode_whitelist, list)

    if not os.path.isfile(filename):
        return None

    with open(filename, 'r') as f:
        values = json.load(f)

    start = (gem_group-1)*len(barcode_whitelist)
    end = gem_group*len(barcode_whitelist)
    barcode_counts = {bc: value for bc, value in zip(barcode_whitelist, values[start:end])}
    if proportions:
        total_barcode_counts = sum(barcode_counts.values())
        barcode_dist = {bc: tk_stats.robust_divide(float(value), float(total_barcode_counts))
                        for bc, value in barcode_counts.iteritems()}
        return barcode_dist
    else:
        return barcode_counts

def get_fasta_iter(f):
    hdr = None
    seq = ''
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            if hdr:
                yield hdr, seq
            hdr = line[1:]
            seq = ''
        else:
            seq += line
    if hdr:
        yield hdr, seq

def detect_paired_end_bam(bam_filename):
    bam = tk_bam.create_bam_infile(bam_filename)
    for read in itertools.islice(bam, 100):
        if read.is_read1 or read.is_read2:
            return True
    return False
    bam.close()

def get_cigar_summary_stats(read, strand):
    """
    Get number of mismatches, insertions, deletions, ref skip, soft clip, hard clip bases from a read.
    Returns a dictionary by the element's CIGAR designation. Adds additional fields to distinguish between three and
    five prime soft-clipping for R1 and R2: R1_S_three_prime and R1_S_five_prime, etc. to account for soft-clipped local alignments.

    Args:
        read (pysam.AlignedRead): aligned read object
        strand (string): + or - to indicate library orientation (MRO argument strand, for example)

    Returns:
        dict of str,int: Key of base type to base counts for metrics. Adds additional fields to distinguish between three and
                            five prime soft-clipping: S_three_prime and S_five_prime.

    """

    statistics = {}
    cigar_tuples = read.cigar

    for i, (category,count) in enumerate(cigar_tuples):
        # Convert numeric code to category
        category = cr_constants.cigar_numeric_to_category_map[category]
        count = int(count)

        # Detect 5 prime soft-clipping
        if i == 0:

            if strand == cr_constants.REVERSE_STRAND:
                metric = 'R1_S_five_prime' if read.is_read1 else 'R2_S_three_prime'
            else:
                metric = 'R2_S_five_prime' if read.is_read1 else 'R1_S_three_prime'

            if category == 'S':
                statistics[metric] = count
            else:
                statistics[metric] = 0

        # Tally up all standard categories from BAM
        if category in statistics:
            statistics[category] += count
        else:
            statistics[category] = count

        # Detect 3 prime soft-clipping
        if i == len(cigar_tuples):
            if strand == cr_constants.REVERSE_STRAND:
                metric = 'R2_S_five_prime' if read.is_read2 else 'R1_S_three_prime'
            else:
                metric = 'R1_S_five_prime' if read.is_read2 else 'R2_S_three_prime'

            if category == 'S':
                statistics[metric] = count
            else:
                statistics[metric] = 0

    return statistics


def get_full_alignment_base_quality_scores(read):
    """
    Returns base quality scores for the full read alignment, inserting zeroes for deletions and removing
    inserted and soft-clipped bases. Therefore, only returns quality for truly aligned sequenced bases.

    Args:
        read (pysam.AlignedSegment): read to get quality scores for

    Returns:
        np.array: numpy array of quality scores

    """

    quality_scores = np.fromstring(read.qual, dtype=np.byte) - tk_constants.ILLUMINA_QUAL_OFFSET

    start_pos = 0

    for operation,length in read.cigar:
        operation = cr_constants.cigar_numeric_to_category_map[operation]

        if operation == 'D':
            quality_scores = np.insert(quality_scores, start_pos, [0] * length)
        elif operation == 'I' or operation == 'S':
            quality_scores = np.delete(quality_scores, np.s_[start_pos:start_pos + length])

        if not operation == 'I' and not operation == 'S':
            start_pos += length

    return start_pos, quality_scores


def get_unmapped_read_count_from_indexed_bam(bam_file_name):
    """
    Get number of unmapped reads from an indexed BAM file.

    Args:
        bam_file_name (str): Name of indexed BAM file.

    Returns:
        int: number of unmapped reads in the BAM

    Note:
        BAM must be indexed for lookup using samtools.

    """

    index_output = tk_subproc.check_output('samtools idxstats %s' % bam_file_name, shell=True)
    return int(index_output.strip().split('\n')[-1].split()[-1])


def kwargs_to_command_line_options(reserved_arguments=set(), sep=" ", long_prefix="--", short_prefix="-",
                                   replace_chars=dict(), **kwargs):
    """
    Convert **kwargs provided by user to a string usable command line programs as arguments.

    Args:
        reserved_arguments (set or list of str): set of arguments that this function prohibited for user use
        sep (str): separator between option/value pairs ('='for --jobmode=sge).
                   WARNING: switch args (--abc), which take no value, break if sep is not ' '
        long_prefix (str): prefix for options with more than one character ("--" for --quiet, for example)
        short_prefix (str): prefix for options with one character ("-" for -q, for example)
        replace_chars (dict): map of characters to replace in specified variable names
                            (if --align-reads is command-line option, specify align_reads with replace chars -> {'_':'-'}
        **kwargs (dict): **kawargs arguments/values to format string for

    Returns:
        str: string formatted appropriately for use as command line option. Returns no arguments provided.

    Raises:
        ValueError if user requested argument conflicts with one of the specified reserved arguments.

    """

    arguments = []
    reserved_arguments = set(reserved_arguments)

    for key, value in kwargs.iteritems():
        normalized_key = key.strip('-')

        # Replace characters for formatting
        for char, new_char in replace_chars.iteritems():
            normalized_key = normalized_key.replace(char, new_char)

        # Validate user inputs to make sure no blatant conflicts
        if normalized_key in reserved_arguments:
            raise ValueError('Specified option conflicts with reserved argument: %s. \
                             Reserved arguments are: %s' % (normalized_key, ','.join(reserved_arguments)))

        # Correctly prefix arguments
        if len(key) > 1:
            prefix = long_prefix
        else:
            prefix = short_prefix

        argument = '%s%s' % (prefix, normalized_key)
        option_value = kwargs.get(key) if kwargs.get(key) is not None else ''
        arguments.append('%s%s%s' % (argument, sep, option_value))

    return ' '.join(arguments)

def write_empty_json(filename):
    with open(filename, 'w') as f:
        json.dump({}, f)

def flatten_list(x):
    """ Flatten a hierarchy of lists """
    result = []
    for y in x:
        if isinstance(y, list):
            result.extend(flatten_list(y))
        else:
            result.append(y)
    return result

def load_barcode_csv(barcode_csv):
    """ Load a csv file of (genome,barcode) """

    bcs_per_genome = collections.defaultdict(list)
    with open(barcode_csv, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            if len(row) != 2:
                raise ValueError("Bad barcode file: %s" % barcode_csv)
            (genome, barcode) = row
            bcs_per_genome[genome].append(barcode)
    return bcs_per_genome

def load_csv_rownames(csv_file):
    rownames = np.atleast_1d(pd.read_csv(csv_file, usecols=[0]).values.squeeze())
    return rownames

def get_h5_filetype(filename):
    with tables.open_file(filename, mode = 'r') as f:
        try:
            filetype = f.get_node_attr('/', cr_constants.H5_FILETYPE_KEY)
        except AttributeError:
            filetype = None # older files lack this key
    return filetype

def save_array_h5(filename, name, arr):
    """ Save an array to the root of an h5 file """
    with tables.open_file(filename, 'w') as f:
        f.create_carray(f.root, name, obj=arr)

def load_array_h5(filename, name):
    """ Load an array from the root of an h5 file """
    with tables.open_file(filename, 'r') as f:
        return getattr(f.root, name).read()

def merge_jsons_single_level(filenames):
    """ Merge a list of toplevel-dict json files.
        Union dicts at the first level.
        Raise exception when duplicate keys encountered.
        Returns the merged dict.
    """

    merged = {}
    for filename in filenames:
        with open(filename) as f:
            d = json.load(f)

            for key, value in d.iteritems():
                if key not in merged:
                    merged[key] = value
                    continue

                # Handle pre-existing keys
                if merged[key] == value:
                    pass
                elif type(value) == dict and type(merged[key]) == dict:
                    merged[key].update(value)
                else:
                    raise ValueError("No merge strategy for key %s, value of type %s, from %s into %s" % (str(key),
                                                                                                          str(type(value)),
                                                                                                          str(value),
                                                                                                          str(merged[key])))
    return merged


def splitexts(s):
    """ Like splitext, but handle concat'd extensions like .tar.gz.
        E.g., /path/to/foo.tar.gz => ('/path/to/foo', '.tar.gz')
        Returns (str,str) - (base, extension) """
    dn = os.path.dirname(s)
    bn = os.path.basename(s)
    parts = bn.split('.')

    if len(parts) == 1:
        return (os.path.join(dn, bn), '')
    else:
        return (os.path.join(dn, parts[0]), '.' + '.'.join(parts[1:]))
