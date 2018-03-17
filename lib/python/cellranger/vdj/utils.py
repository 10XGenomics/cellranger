#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

import itertools
import json
import numpy as np
import os
import pandas as pd
import pysam
import sys
import tenkit.bam as tk_bam
import tenkit.cache as tk_cache
import cellranger.constants as cr_constants
import cellranger.fastq as cr_fastq
import cellranger.vdj.constants as vdj_constants


def parse_contig_name(name):
    """ Parse the assembler contig name """
    word = name.split(' ')[0]
    fields = word.split('_')
    bc = fields[0]
    cluster = '_'.join(fields[0:3]),
    gene = '_'.join(fields[0:4])
    isoform = '_'.join(fields[0:5])
    return bc, cluster, gene, isoform

def get_genes_in_pair(gene_pair_str):
    return gene_pair_str.split('_')

def get_recombinome_gene_name(ref_name):
    """ Get the vdj gene name for a recombinome ref name """
    for vdj_gene_name in vdj_constants.VDJ_GENES:
        if ref_name.startswith(vdj_gene_name):
            return vdj_gene_name
    return None

def load_contig_summary_table(filename):
    df = pd.read_csv(filename)
    # Index on contig id
    df.set_index('contig_id', drop=False, inplace=True, verify_integrity=True)
    return df

def get_barcode_from_contig_name(contig_name):
    return contig_name.split('_')[0]

def is_contig_filtered(contig_df, contig_id):
    """ Return true if a contig passed filter.
      contig_df (pd.DataFrame): Contains contig ids and counts. """
    if 'pass_filter' in contig_df.columns:
        return contig_df['pass_filter'][contig_id]
    else:
        # This is expected for the new contig summary coming from the ASSEMBLE stage
        sys.stderr.write("Checked if contig passes filter but contig summary did not have the column\n")
        return True

def save_contig_summary_table(data_frame, filename):
    data_frame.to_csv(filename, index=False)

def get_fastq_read_raw_barcode(fastq_read):
    return cr_fastq.AugmentedFastqHeader(fastq_read[0]).get_tag(cr_constants.RAW_BARCODE_TAG)

def get_fastq_read_barcode(fastq_read):
    return cr_fastq.AugmentedFastqHeader(fastq_read[0]).get_tag(cr_constants.PROCESSED_BARCODE_TAG)

def fastq_barcode_sort_key(fastq_read):
    """ Return barcode, qname """
    fastq_header = cr_fastq.AugmentedFastqHeader(fastq_read[0])
    bc = fastq_header.get_tag(cr_constants.PROCESSED_BARCODE_TAG)
    return bc, fastq_header.fastq_header

def load_cell_barcodes_json(filename):
    with open(filename) as f:
        return json.load(f)

def write_csv_optional(df, filename):
    """ Write a pandas dataframe to CSV if it is not None """
    if df is not None:
        write_csv(df, filename)

def write_csv(df, filename):
    """ Write a pandas dataframe to CSV in a standard way """

    # Verify that the data do not contain commas
    for col in df.select_dtypes([np.object]):
        if df[col].str.contains(',').any():
            raise ValueError("Failed write to %s: Column %s contains commas" % (filename, col))

    df.to_csv(filename, header=True, index=False, sep=',')

def write_csv_row(row, f):
    """ Write a standard CSV row to an open file """
    # Verify that the data do not contain commas
    row = map(str, row)
    for i,v in enumerate(row):
        if ',' in v:
            raise ValueError("Failed write to csv file: Column %d contains commas" % i)
    f.write('%s\n' % ','.join(row))

def format_clonotype_id(clonotype_index, inferred):
    # Takes a 0-based clonotype index and formats it into a clonotype id string
    if clonotype_index is None:
        return None
    prefix = 'inferred_clonotype' if inferred else 'clonotype'
    return '%s%d' % (prefix, 1 + clonotype_index)

def get_mem_gb_from_annotations_json(filename):
    """ Estimate mem request for loading an entire annotations json into memory """
    return int(np.ceil(vdj_constants.MEM_GB_PER_ANNOTATIONS_JSON_GB * float(os.path.getsize(filename))/1e9))

def concatenate_and_fix_bams(out_bamfile, bamfiles, drop_tags=None):
    """Concatenate bams with different references and populate tags.

    Assumes that the input bams have completely non-overlapping reference entries.
    No sorting assumed. Output reads as in the same order as input.
    Tags are stripped from the header (which is assumed to look like an AugmentedFastqHeader)
    and are stored as bam tags. Tags in drop_tags are completely dropped.

    Args:
    - drop_tags: Set of tag names to ignore. If None, than all tags will be included.
    """
    template_bam = pysam.Samfile(bamfiles[0])
    header = template_bam.header

    for bam_fn in bamfiles[1:]:
        bam = pysam.Samfile(bam_fn)
        header['SQ'].extend(bam.header['SQ'])

    refname_to_tid = { ref_head['SN']:idx for (idx, ref_head) in enumerate(header['SQ']) }

    bam_out = pysam.Samfile(out_bamfile, "wb", header=header)

    for bam_fn in bamfiles:
        bam = pysam.Samfile(bam_fn)

        for rec in bam:
            if rec.is_unmapped and rec.mate_is_unmapped:
                rec.reference_id = -1
                rec.next_reference_id = -1
            elif rec.is_unmapped and not rec.mate_is_unmapped:
                rec.next_reference_id = refname_to_tid[bam.references[rec.next_reference_id]]
                rec.reference_id = rec.next_reference_id
            elif not rec.is_unmapped and rec.mate_is_unmapped:
                rec.reference_id = refname_to_tid[bam.references[rec.reference_id]]
                rec.next_reference_id = rec.reference_id
            else:
                rec.reference_id = refname_to_tid[bam.references[rec.reference_id]]
                rec.next_reference_id = refname_to_tid[bam.references[rec.next_reference_id]]

            # Pull fields out of qname, and make tags, writing out to bam_out
            qname_split = rec.qname.split(cr_fastq.AugmentedFastqHeader.TAG_SEP)
            rec.qname = qname_split[0]

            tags = rec.tags
            for i in range(1, len(qname_split), 2):
                tag = (qname_split[i], qname_split[i+1])
                # Don't add empty tags
                if len(tag[1]) > 0 and (drop_tags is None or not tag[0] in drop_tags):
                    tags.append(tag)

            rec.tags = tags
            bam_out.write(rec)


def bam_has_seqs(filename):
    """ Return true if a bam has any reference sequences defined """
    try:
        f = tk_bam.create_bam_infile(filename)
        f.close()
        return True
    except ValueError:
        return False

""" Generator that streams items from a list of dicts [{}, {},...] """
def get_json_obj_iter(f):
        brace_count = 0
        x = ''
        in_str = False
        # Track the number of consecutive backslashes encountered preceding this char
        backslash_count = 0
        bufsize = 4096

        while True:
            buf = f.read(bufsize)
            if len(buf) == 0:
                return
            for c in buf:
                # Handle escaped quotes, possibly preceded by escaped backslashes
                # "foo\" => stay in string, 'foo"...'
                # "foo\\" => leave string, 'foo\'
                # "foo\\\" => stay in string, 'foo\"...'
                if c == '"' and (backslash_count % 2) == 0:
                    in_str = not in_str

                if c == '\\':
                    backslash_count += 1
                else:
                    backslash_count = 0

                if not in_str:
                    if c.isspace():
                        continue
                    if brace_count == 0 and c == '{':
                        brace_count = 1
                    elif brace_count == 1 and c == '}':
                        brace_count = 0
                    elif c == '{' or c == '}':
                        brace_count += 1 if c == '{' else -1
                if brace_count != 0 or len(x) > 0:
                    x += c

                if brace_count == 0 and len(x) > 0:
                    yield json.loads(x)
                    x = ''

""" Streams a list of dicts as json to a file """
class JsonDictListWriter(object):
    def __init__(self, fh):
        """ fh - file handle """
        # Track whether we need to prepend with a comma
        self.wrote_any = False
        # Start the list
        fh.write('[\n')

    def write(self, x, fh):
        """ Write a dict """
        if self.wrote_any:
            fh.write('\n,')
        json.dump(x, fh)
        self.wrote_any = True

    def finish(self, fh):
        """ Finish writing all dicts """
        # End the list
        fh.write('\n]')

class CachedJsonDictListWriters(object):
    """ Stream multiple json DictLists, using a FileHandleCache
        to limit the number of open file handles.
        This is useful when you're writing to many files from the same process
        (e.g., in a 'split' job) """
    def __init__(self, filenames):
        """ filenames list(str) - list of filenames to write to """
        self.filenames = filenames
        self.cache = tk_cache.FileHandleCache(mode='w')
        self.writers = [JsonDictListWriter(self.cache.get(fn)) for fn in filenames]

    def __enter__(self):
        return self

    def __exit__(self, e_type, e_val, e_tb):
        for writer, filename in itertools.izip(self.writers, self.filenames):
            writer.finish(self.cache.get(filename))

    def write(self, d, file_idx):
        """ d (dict) - data to write
            file_idx (int) - index of filename/writer in original given list """
        self.writers[file_idx].write(d, self.cache.get(self.filenames[file_idx]))
