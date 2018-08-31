#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
import collections
import itertools
import json
import numpy as np
import os
import re
import sys
import tenkit.constants as tk_constants
import tenkit.fasta as tk_fasta
import tenkit.safe_json as tk_safe_json
import tenkit.seq as tk_seq
import cellranger.constants as cr_constants
import cellranger.h5_constants as h5_constants
import cellranger.utils as cr_utils
import cellranger.io as cr_io

def get_bamtofastq_defs(read_defs, destination_tags):
    """ Determine which portions of reads need to be retained.
        Args: read_defs - list(ReadDef)
              destination_tags - list((str,str)) - list of (seq_tag, qual_tag) or
                                                   None if the read comes from the READ/QUAL fields
        NOTE: Assumes that the read defs are contiguous - sequence in gaps between extracted seqs will be lost.
        Returns dict of {read_type (str): str} """

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

        # Interpreted by bamtofastq
        trim_defs[read_type] = '10x_bam_to_fastq:' + read_type + '(' + ','.join(bam_to_fastq_entries) + ')'

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

    if reg_valid_count:
        return rc_valid_count >= ((rc_valid_count + reg_valid_count) *
                                  cr_constants.REVCOMP_BARCODE_THRESHOLD)
    return False

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

    def to_json(self):
        if self.barcode_seqs:
            return list(self.barcode_counts)
        else:
            return []

    def close(self):
        if self.barcode_seqs:
            with open(self.out_counts, 'w') as f:
                tk_safe_json.dump_numpy(self.to_json(), f)

    @staticmethod
    def merge_by(counter_files, keys, barcode_whitelist, gem_groups):
        """ Merge BarcodeCounters by a key.
        Args:
          counter_files (list of str): Filenames of BarcodeCounter outputs
          keys (list of str): Keys to group by
          barcode_whitelist (list of str): Same as BarcodeCounter constructor
          gem_groups (list of int): Same as BarcodeCounter constructor
        Returns:
          dict of str:dict: Keys are the group keys and dicts are serialized BarcodeCounters
        """
        distinct_keys = sorted(list(set(keys)))
        groups = {}
        for key in distinct_keys:
            groups[key] = BarcodeCounter(barcode_whitelist, None, gem_groups)
        for key, filename, gg in zip(keys, counter_files, gem_groups):
            groups[key].merge(gg, filename)
        return {key: group.to_json() for key, group in groups.iteritems()}

def extract_read_maybe_paired(read_tuple, read_def, reads_interleaved, r1_length=None, r2_length=None):
    """ Args: read_tuple: (name, read, qual)
              read_def: ReadDef object
              reads_interleaved: bool
              r1_length (int): Hard trim on 3' end of input R1
              r2_length (int): Hard trim on 3' end of input R2
        Yields: ExtractReadResult """
    if reads_interleaved and read_def.read_type == 'R1':
        name, seq, qual = cr_utils.get_fastq_read1(read_tuple, None, True)
    elif reads_interleaved and read_def.read_type == 'R2':
        name, seq, qual = cr_utils.get_fastq_read2(read_tuple, None, True)
    else:
        name, seq, qual = read_tuple

    # Apply hard trimming on input
    hard_end = sys.maxint
    if read_def.read_type == 'R1' and r1_length is not None:
        hard_end = r1_length
    elif read_def.read_type == 'R2' and r2_length is not None:
        hard_end = r2_length

    # Extract interval requested by the read def
    if read_def.length is not None:
        end = min(hard_end, read_def.offset + read_def.length)
        read_slice = slice(read_def.offset, end)
    else:
        read_slice = slice(read_def.offset, hard_end)

    return (name, seq[read_slice], qual[read_slice])

def get_read_generator_fastq(fastq_open_file, read_def, reads_interleaved, r1_length=None, r2_length=None):
    read_iter = tk_fasta.read_generator_fastq(fastq_open_file, paired_end=reads_interleaved and read_def.read_type in ['R1', 'R2'])
    for read_tuple in read_iter:
        yield extract_read_maybe_paired(read_tuple, read_def, reads_interleaved, r1_length, r2_length)

def get_feature_generator_fastq(files, extractor, interleaved, read_types, r1_length=None, r2_length=None):
    '''Extract feature barcodes from FASTQs.

    Args:
       files (list of File): FASTQ file handles for R1, R2
       extractor (FeatureExtractor): Extracts feature barcodes
       interleaved (bool): Are R1,R2 interleaved in a single file
       read_types (list of str): List of read types (e.g. R1,R2) we need to inspect
       r1_length (int): Length to hard-trim R1 to
       r2_length (int): Length to hard-trim R2 to
    Returns:
       FeatureMatchResult: Yields the feature extraction result for a read pair
'''
    assert len(files) == 2
    assert 'R1' in read_types or 'R2' in read_types

    # Apply hard trimming on input
    r1_hard_end = sys.maxint if r1_length is None else r1_length
    r2_hard_end = sys.maxint if r2_length is None else r2_length

    if interleaved:
        f = files[0]
        assert f
        # Get R1 and R2 seqs from interleaved FASTQ
        pair_iter = itertools.imap(lambda x: (x[0:3], x[3:6]),
                                   tk_fasta.read_generator_fastq(f, paired_end=True))
    else:
        r1_iter = tk_fasta.read_generator_fastq(files[0], paired_end=False) if 'R1' in read_types else iter([])
        r2_iter = tk_fasta.read_generator_fastq(files[1], paired_end=False) if 'R2' in read_types else iter([])
        pair_iter = itertools.izip_longest(r1_iter, r2_iter)

    if read_types == ['R1']:
        match_func = lambda x: extractor.extract_single_end(x[0][1][0:r1_hard_end], # seq
                                                            x[0][2][0:r1_hard_end], # qual
                                                            'R1')

    elif read_types == ['R2']:
        match_func = lambda x: extractor.extract_single_end(x[1][1][0:r2_hard_end], # seq
                                                            x[1][2][0:r2_hard_end], # qual
                                                            'R2')

    elif read_types == ['R1', 'R2']:
        match_func = lambda x: extractor.extract_paired_end(x[0][1][0:r1_hard_end], # seq
                                                            x[0][2][0:r1_hard_end], # qual
                                                            x[1][1][0:r2_hard_end], # seq
                                                            x[1][2][0:r2_hard_end]) # qual

    return itertools.imap(match_func, pair_iter)


def get_fastq_from_read_type(fastq_dict, read_def, reads_interleaved):
    ''' Use a ReadDef to determine which FASTQ files to open '''
    if read_def.read_type == 'R2' and reads_interleaved:
        return fastq_dict.get('R1')
    else:
        return fastq_dict.get(read_def.read_type)

def get_fastqs_from_feature_ref(fastq_dict, reads_interleaved, read_types):
    ''' Determine which FASTQ files to open for a FeatureExtractor'''
    assert 'R1' in read_types or 'R2' in read_types

    fastq1 = fastq_dict.get('R1') if 'R1' in read_types else None

    if reads_interleaved:
        # All reads come from the first FASTQ ("R1")
        fastq2 = fastq_dict.get('R1') if 'R2' in read_types else None
    else:
        # Read2 comes from the second FASTQ ("R2")
        fastq2 = fastq_dict.get('R2') if 'R2' in read_types else None

    return (fastq1, fastq2)

class FastqReader:
    # Extracts specified regions from input fastqs
    # For example, extract UMIs from the first 10 bases of the "R2" read
    def __init__(self, in_filenames, read_def, reads_interleaved, r1_length, r2_length):
        """ Args:
              in_filenames - Map of paths to fastq files
              read_def - ReadDef
        """

        self.in_fastq = None
        self.in_iter = iter([])
        self.read_def = read_def

        if in_filenames:
            in_filename = get_fastq_from_read_type(in_filenames,
                                                   read_def,
                                                   reads_interleaved)
            if in_filename:
                self.in_fastq = cr_io.open_maybe_gzip(in_filename, 'r')

                self.in_iter = get_read_generator_fastq(self.in_fastq,
                                                        read_def=read_def,
                                                        reads_interleaved=reads_interleaved,
                                                        r1_length=r1_length,
                                                        r2_length=r2_length)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def close(self):
        if self.in_fastq:
            self.in_fastq.close()
            self.in_fastq = None

class FastqFeatureReader:
    ''' Use a FeatureReference to extract specified feature barcodes from input fastqs
        For example, extract antibody barcodes from R2. '''

    def __init__(self, in_filenames, extractor, reads_interleaved, r1_length, r2_length):
        """ Args:
              in_filenames (dict of str -> str): Map of paths to fastq files
              feature_ref (FeatureExtractor): for extracting feature barcodes
        """

        self.in_fastqs = None
        self.in_iter = iter([])

        # Relevant read types
        read_types = extractor.get_read_types()

        if in_filenames:
            in_filenames = get_fastqs_from_feature_ref(in_filenames,
                                                       reads_interleaved,
                                                       read_types)
            if in_filenames != (None, None):
                if reads_interleaved:
                    filename = in_filenames[0] if in_filenames[0] else in_filenames[1]
                    self.in_fastqs = (cr_io.open_maybe_gzip(filename, 'r') if filename[0] else None,
                                      None)
                else:
                    self.in_fastqs = (cr_io.open_maybe_gzip(in_filenames[0], 'r') if in_filenames[0] else None,
                                      cr_io.open_maybe_gzip(in_filenames[1], 'r') if in_filenames[1] else None)

                self.in_iter = get_feature_generator_fastq(files=self.in_fastqs,
                                                           extractor=extractor,
                                                           interleaved=reads_interleaved,
                                                           read_types=read_types,
                                                           r1_length=r1_length,
                                                           r2_length=r2_length)
    def close(self):
        if self.in_fastqs:
            for f in self.in_fastqs:
                if f:
                    f.close()

class ChunkedWriter:
    # Writes sequencing read-based output, splitting into chunks as specified
    def __init__(self, base_path, max_reads_per_file):
        self.reads_per_file = 0
        self.max_reads_per_file = max_reads_per_file
        self.index = 0
        self.base_path = base_path
        self.file_paths = []
        self.curr_file = None

        cr_io.mkdir(base_path, allow_existing=True)

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
        compression = kwargs.pop('compression')

        ChunkedWriter.__init__(self, *args, **kwargs)

        if compression is None:
            self.suffix = ''
        elif compression == 'gzip':
            self.suffix = h5_constants.GZIP_SUFFIX
        elif compression == 'lz4':
            self.suffix = h5_constants.LZ4_SUFFIX
        else:
            raise ValueError('Unknown compression type: %s' % compression)
        self.compression = compression

    def generate_filename(self):
        return '%d.fastq%s' % (self.index, self.suffix or '')

    def open_file(self, filename):
        return cr_io.open_maybe_gzip(filename, 'w')

    def write_data(self, data):
        tk_fasta.write_read_fastq(self.curr_file, *data)

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

# Copied from tenkit because tenkit/preflight.py is utterly broken (imports martian)
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


class FastqSpecException(Exception):
    def __init__(self, msg):
        self.msg = msg
    def __str__(self):
        return self.msg

def _require_sample_def_key(sd, key):
    if key not in sd:
        raise FastqSpecException('Sample def is missing the key "%s."' % key)

class FastqSpec(object):
    """ All info required to find FASTQ files """
    def __init__(self, fastq_mode, read_path, lanes, sample_indices, sample_names, interleaved):
        self.fastq_mode = fastq_mode
        self.read_path = read_path
        self.lanes = lanes
        self.sample_indices = sample_indices
        self.sample_names = sample_names
        self.interleaved = interleaved

    @staticmethod
    def from_sample_def(sd):
        _require_sample_def_key(sd, 'fastq_mode')
        _require_sample_def_key(sd, 'read_path')
        _require_sample_def_key(sd, 'lanes')

        if sd['fastq_mode'] == tk_constants.ILMN_BCL2FASTQ_FASTQ_MODE:
            _require_sample_def_key(sd, 'sample_names')

            return FastqSpec(fastq_mode=sd['fastq_mode'],
                             read_path=sd['read_path'],
                             lanes=sd['lanes'],
                             sample_indices=None,
                             sample_names=sd['sample_names'],
                             interleaved=False,
            )

        elif sd['fastq_mode'] == tk_constants.BCL_PROCESSOR_FASTQ_MODE:
            _require_sample_def_key(sd, 'sample_indices')

            # Check and optionally translate sample index sets to
            # sample index sequences
            si_strings, msg = check_sample_indices(sd)
            if si_strings is None:
                raise FastqSpecException(msg)

            return FastqSpec(fastq_mode=sd['fastq_mode'],
                             read_path=sd['read_path'],
                             lanes=sd['lanes'],
                             sample_indices=si_strings,
                             sample_names=None,
                             interleaved=True,
            )

        else:
            raise FastqSpecException('Sample def contained unrecognized fastq_mode  "%s."' % sd['fastq_mode'])

    def get_group_spec_iter(self):
        """ Each group corresponds to a single sample_name or sample_index in
            the original sample_def. This method slices on that.
            Yields: (sample_index_or_name, FastQSpec) """
        if self.sample_indices is not None:
            for si in self.sample_indices:
                yield si, FastqSpec(fastq_mode=self.fastq_mode,
                                    read_path=self.read_path,
                                    lanes=self.lanes,
                                    sample_indices=[si],
                                    sample_names=None,
                                    interleaved=self.interleaved)
        elif self.sample_names is not None:
            for sn in self.sample_names:
                yield sn, FastqSpec(fastq_mode=self.fastq_mode,
                                    read_path=self.read_path,
                                    lanes=self.lanes,
                                    sample_indices=None,
                                    sample_names=[sn],
                                    interleaved=self.interleaved)
        else:
            raise ValueError("Cannote iterate over FastqSpec with no sample indices or names specified.")

    def is_single_group(self):
        """ Returns true if this spec contains a single sample index/name """
        return (self.sample_indices is not None and len(self.sample_indices) == 1) or \
            (self.sample_names is not None and len(self.sample_names) == 1)

    def get_fastqs(self, read_type):
        """ read_type (str) - One of RA,R1,R2,R3,R4,I1,I2 """

        if read_type == 'RA' and not self.interleaved:
            raise ValueError('Read type "%s" was requested but is only supported for non-interleaved FASTQs.' % read_type)

        # If interleaved, translate R1|R2 => RA
        if self.interleaved and read_type in ('R1', 'R2'):
            read_type = 'RA'

        fastqs = []

        if self.fastq_mode == tk_constants.BCL_PROCESSOR_FASTQ_MODE:
            for group, _ in self.get_group_spec_iter():
                fastqs.extend(tk_fasta.find_input_fastq_files_10x_preprocess(
                    self.read_path, read_type, group, self.lanes))

        elif self.fastq_mode == tk_constants.ILMN_BCL2FASTQ_FASTQ_MODE:
            for group, _ in self.get_group_spec_iter():
                fastqs.extend(tk_fasta.find_input_fastq_files_bcl2fastq_demult(
                    self.read_path, read_type, group, self.lanes))

        return fastqs

    def __str__(self):
        return str(self.__dict__)
