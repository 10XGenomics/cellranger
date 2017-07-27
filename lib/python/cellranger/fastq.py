#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
import collections
import itertools
import json
import numpy as np
import os
import tenkit.bam as tk_bam
import tenkit.fasta as tk_fasta
import tenkit.safe_json as tk_safe_json
import tenkit.seq as tk_seq
import tenkit.stats as tk_stats
import cellranger.constants as cr_constants
import cellranger.utils as cr_utils

# Specify which parts of the read need to be retained.
#  pre/post-slice: The prefix and postfix slices of the "discarded" regions of the read
#  bam_to_fastq: String that defines how to reconstruct the fastq from the trimmed regions and bam tags
TrimDef = collections.namedtuple('TrimDef', ['pre_slice', 'post_slice', 'bam_to_fastq'])

# The result of extracting a substring from a read.
#  read - tuple of (name, seq, qual) for extracted sequence.
#  trimmed_seq - remaining sequence
#  trimmed_qual - remaining sequence
ExtractReadResult = collections.namedtuple('ExtractReadResult', ['read',
                                                                 'trimmed_seq',
                                                                 'trimmed_qual'])

def compute_trim_defs(read_defs, destination_tags, retain_trimmed_suffix_read):
    """ Determine which portions of reads need to be retained.
        Args: read_defs - list(ReadDef)
              destination_tags - list((str,str)) - list of (seq_tag, qual_tag) or
                                                   None if the read comes from the READ/QUAL fields
              retain_trimmed_suffix_read         - The read for which to retain a trimmed suffix
        NOTE: Assumes that the read defs are contiguous - sequence in gaps between extracted seqs will be lost.
        Returns dict of {read_type (str): TrimDef} """

    assert len(read_defs) == len(destination_tags)

    # Map read defs to their destination bam tags
    dest_tags_dict = {read_def: dest_tags for read_def, dest_tags in itertools.izip(read_defs, destination_tags)}

    # Bin read defs by read type (R1, R2, etc.)
    read_types = collections.defaultdict(list)
    for read_def in read_defs:
        if read_def.read_type is None:
            continue
        read_types[read_def.read_type].append(read_def)

    trim_defs = {}

    for read_type, read_type_defs in sorted(read_types.items()):
        # Get tags that were built from this read
        read_type_defs.sort(key=lambda rd: rd.offset)
        read_type_dest_tags = [dest_tags_dict[rd] for rd in read_type_defs]

        # Construct bam_to_fastq entries
        bam_to_fastq_entries = []
        for dest_tag in read_type_dest_tags:
            if dest_tag is None:
                bam_to_fastq_entries.append('SEQ:QUAL')
            else:
                bam_to_fastq_entries.append('%s:%s' % dest_tag)

        # Special case for reads where we retain the suffix
        # Prefix trimming unimplemented
        pre_slice = slice(0,0)
        post_slice = slice(0,0)
        if retain_trimmed_suffix_read == read_type:
            bam_to_fastq_entries.append('%s:%s' % (cr_constants.TRIMMED_SEQ_TAG, cr_constants.TRIMMED_QUAL_TAG))
            post_slice = slice(read_type_defs[-1].offset + read_type_defs[-1].length, None)

        bam_to_fastq = '10x_bam_to_fastq:' + read_type + '(' + ','.join(bam_to_fastq_entries) + ')'

        trim_defs[read_type] = TrimDef(pre_slice=pre_slice,
                                       post_slice=post_slice,
                                       bam_to_fastq=bam_to_fastq)
    return trim_defs


def infer_barcode_reverse_complement(barcode_whitelist, read_iter):
    if barcode_whitelist is None:
        return False
    reg_valid_count = 0
    rc_valid_count = 0
    for name, seq, qual in itertools.islice(read_iter, cr_constants.NUM_CHECK_BARCODES_FOR_ORIENTATION):
        if seq in barcode_whitelist:
            reg_valid_count += 1
        if tk_seq.get_rev_comp(seq) in barcode_whitelist:
            rc_valid_count += 1

    frac_rc = tk_stats.robust_divide(rc_valid_count, rc_valid_count + reg_valid_count)
    return frac_rc >= cr_constants.REVCOMP_BARCODE_THRESHOLD

class BarcodeCounter:
    def __init__(self, barcode_whitelist, out_counts, gem_groups=None):
        self.barcode_counts = None
        self.barcode_index = None
        self.out_counts = out_counts
        self.barcode_seqs = cr_utils.load_barcode_whitelist(barcode_whitelist)
        if self.barcode_seqs:
            self.barcode_seqs = cr_utils.format_barcode_seqs(self.barcode_seqs, gem_groups)
            self.barcode_counts = np.zeros(len(self.barcode_seqs), dtype=np.uint32)
            self.barcode_index = {bc: index for index, bc in enumerate(self.barcode_seqs)}

    def count(self, name, seq, qual):
        if self.barcode_seqs:
            index = self.barcode_index.get(seq)
            if index is not None:
                self.barcode_counts[index] += 1

    def merge(self, gem_group, out_counts):
        if self.barcode_seqs:
            with open(out_counts, 'r') as f:
                barcode_counts = json.load(f)

            start = (gem_group-1)*len(barcode_counts)
            end = gem_group*len(barcode_counts)
            self.barcode_counts[start:end] += np.array(barcode_counts, dtype=np.uint32)

    def close(self):
        if self.barcode_seqs:
            with open(self.out_counts, 'w') as f:
                tk_safe_json.dump_numpy(list(self.barcode_counts), f)

def extract_read_maybe_paired(read_tuple, read_def, reads_interleaved, trim_def):
    """ Args: read_tuple: (name, read, qual)
        Yields: ExtractReadResult """
    if reads_interleaved and read_def.read_type == 'R1':
        name, seq, qual = cr_utils.get_fastq_read1(read_tuple, None, True)
    elif reads_interleaved and read_def.read_type == 'R2':
        name, seq, qual = cr_utils.get_fastq_read2(read_tuple, None, True)
    else:
        name, seq, qual = read_tuple

    if read_def.length is not None:
        read_slice = slice(read_def.offset, read_def.offset + read_def.length)
    else:
        read_slice = slice(read_def.offset, None)

    # Store the trimmed sequence
    if trim_def is not None:
        trimmed_seq = seq[trim_def.pre_slice] + seq[trim_def.post_slice]
        trimmed_qual = qual[trim_def.pre_slice] + qual[trim_def.post_slice]
    else:
        trimmed_seq = ''
        trimmed_qual = ''

    return ExtractReadResult((name, seq[read_slice], qual[read_slice]),
                             trimmed_seq, trimmed_qual)

def get_read_generator_fastq(fastq_open_file, read_def, reads_interleaved, trim_def):
    read_iter = tk_fasta.read_generator_fastq(fastq_open_file, paired_end=reads_interleaved and read_def.read_type in ['R1', 'R2'])
    for read_tuple in read_iter:
        yield extract_read_maybe_paired(read_tuple, read_def, reads_interleaved, trim_def)

def get_fastq_from_read_type(fastq_dict, read_def, reads_interleaved):
    if read_def.read_type == 'R2' and reads_interleaved:
        return fastq_dict.get('R1')
    else:
        return fastq_dict.get(read_def.read_type)

class FastqReader:
    # Extracts specified regions from input fastqs
    # For example, extract UMIs from the first 10 bases of the "R2" read
    def __init__(self, in_filenames, read_def, reads_interleaved, trim_defs):
        """ Args:
              in_filenames - list(str): Paths to fastq files
              read_def - ReadDef
              trim_defs - dict of {read_type (str): TrimDef} """

        self.in_fastq = None
        self.in_iter = iter([])
        self.read_def = read_def

        if in_filenames:
            in_filename = get_fastq_from_read_type(in_filenames,
                                                   read_def,
                                                   reads_interleaved)
            if in_filename:
                self.in_fastq = cr_utils.open_maybe_gzip(in_filename, 'r')

                if trim_defs is None:
                    read_trim_def = None
                else:
                    read_trim_def = trim_defs[read_def.read_type]

                self.in_iter = get_read_generator_fastq(self.in_fastq,
                                                        read_def=read_def,
                                                        reads_interleaved=reads_interleaved,
                                                        trim_def=read_trim_def)

    def close(self):
        if self.in_fastq:
            self.in_fastq.close()

class ChunkedWriter:
    # Writes sequencing read-based output, splitting into chunks as specified
    def __init__(self, base_path, max_reads_per_file):
        self.reads_per_file = 0
        self.max_reads_per_file = max_reads_per_file
        self.index = 0
        self.base_path = base_path
        self.file_paths = []
        self.curr_file = None

        cr_utils.mkdir(base_path, allow_existing=True)

    def write(self, data):
        if self.reads_per_file >= self.max_reads_per_file or self.curr_file is None:
            self.close()

            self.index += 1
            out_filename = os.path.join(self.base_path, self.generate_filename())
            self.curr_file = self.open_file(out_filename)
            self.file_paths.append(out_filename)

            self.reads_per_file = 0

        self.write_data(data)

        self.reads_per_file += 1

    def get_out_paths(self, default_len=None):
        assert self.file_paths or default_len

        if self.file_paths:
            return self.file_paths
        else:
            return [None] * default_len

    def close(self):
        if self.curr_file is not None:
            self.curr_file.close()
            self.curr_file = None

class ChunkedFastqWriter(ChunkedWriter):
    def __init__(self, *args, **kwargs):
        ChunkedWriter.__init__(self, *args, **kwargs)

    def generate_filename(self):
        return '%d.fastq' % self.index

    def open_file(self, filename):
        return open(filename, 'w')

    def write_data(self, data):
        tk_fasta.write_read_fastq(self.curr_file, *data)

class ChunkedBamWriter(ChunkedWriter):
    def __init__(self, *args, **kwargs):
        ChunkedWriter.__init__(self, *args, **kwargs)

    def generate_filename(self):
        return '%d.bam' % self.index

    def open_file(self, filename):
        # Create a dummy header to prevent samtools / pysam crashing
        return tk_bam.create_bam_outfile(filename, ['dummy'], [8])[0]

    def write_data(self, data):
        self.curr_file.write(data)

class AugmentedFastqHeader:
    """ Store 10x specific tags in fastq qname """
    TAG_SEP = '|||'
    WORD_SEP = ' '

    def __init__(self, fastq_header):
        """ Parse the fastq header """
        words = fastq_header.split(self.WORD_SEP)
        # Assume that TAG_SEP doesn't exist in the original fastq header
        fields = words[0].split(self.TAG_SEP)
        stripped_word, tag_pairs = fields[0], fields[1:]

        self.fastq_header = AugmentedFastqHeader.WORD_SEP.join([stripped_word] + words[1:])
        self.tags = zip(tag_pairs[::2], tag_pairs[1::2])

    def set_tag(self, key, value):
        for i, (k,_) in enumerate(self.tags):
            if k == key:
                self.tags[i] = (key, value)
                return
        self.tags.append((key,value))

    def get_tag(self, key):
        for k, v in self.tags:
            if k == key:
                return v
        return None

    def to_string(self):
        """ Append a TAG_SEP-delimited tag-dict onto the first word of the fastq name """
        hdr_words = self.fastq_header.split(self.WORD_SEP)
        tag_strings = [self.TAG_SEP.join(item) for item in self.tags]
        augmented_word = self.TAG_SEP.join([hdr_words[0]] + tag_strings)

        return self.WORD_SEP.join([augmented_word] + hdr_words[1:])
