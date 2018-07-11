#!/usr/bin/env python
#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#
# General utilities for manipulating nucleotide sequences
#

import string
import subprocess
import os.path
from tenkit.exceptions import NotSupportedException

NUCS = ['A', 'C', 'G', 'T']
NUCS_INVERSE = {'A':0, 'C':1, 'G':2, 'T':3}

DNA_CONVERT_TABLE = string.maketrans('ACGTacgtRYMKBDHVrymkbdhv', 'TGCAtgcaYRKMVHDByrkmvhdb')
RNA_CONVERT_TABLE = string.maketrans('ACGUacguRYMKBDHVrymkbdhv', 'UGCAugcaYRKMVHDByrkmvhdb')

IUPAC_NUC_MAP = {}
IUPAC_NUC_MAP['A'] = ['A']
IUPAC_NUC_MAP['C'] = ['C']
IUPAC_NUC_MAP['G'] = ['G']
IUPAC_NUC_MAP['T'] = ['T']
IUPAC_NUC_MAP['R'] = ['A', 'G']
IUPAC_NUC_MAP['Y'] = ['C', 'T']
IUPAC_NUC_MAP['M'] = ['C', 'A']
IUPAC_NUC_MAP['K'] = ['T', 'G']
IUPAC_NUC_MAP['W'] = ['T', 'A']
IUPAC_NUC_MAP['S'] = ['C', 'G']
IUPAC_NUC_MAP['B'] = ['C', 'T', 'G']
IUPAC_NUC_MAP['D'] = ['T', 'A', 'G']
IUPAC_NUC_MAP['H'] = ['T', 'A', 'C']
IUPAC_NUC_MAP['V'] = ['A', 'C', 'G']
IUPAC_NUC_MAP['N'] = ['T', 'A', 'C', 'G']

def get_rev_comp(seq):
    """ Reverse complement for DNA.  Included ambiguous nucleotides and retains case.
    """
    return str(seq).translate(DNA_CONVERT_TABLE)[::-1]

def get_rev_comp_rna(seq):
    """ Reverse complement for RNA.  Included ambiguous nucleotides and retains case.
    """
    return seq.translate(RNA_CONVERT_TABLE)[::-1]

def get_rec_seqs(site_seq):
    """ Code that takes an IUPAC consensus sequence and returns the full list of nucleotide
    sequences recognizable by this site seq
    """
    rec_seqs = ['']
    for n in xrange(len(site_seq)):
        nucs = IUPAC_NUC_MAP[site_seq[n]]
        new_rec_seqs = []
        for nuc in nucs:
            for rec_seq in rec_seqs:
                new_rec_seq = rec_seq + nuc
                new_rec_seqs.append(new_rec_seq)
        rec_seqs = new_rec_seqs

    return rec_seqs

def mask(seq, keep_start, keep_end):
    """Mask the sequence leaving only [keep_start, keep_end) unmasked"""
    return 'N' * keep_start + seq[keep_start:keep_end] + 'N' * (len(seq) - keep_end)

def get_cigar_map(cigar):
    """Takes a cigar (as a tuple) and returns a list giving the offsets
    for each position of a read
    """
    if cigar is None:
        return None

    cigar_map = []
    offset = 0
    for (categ, length) in cigar:
        # Aligned base
        if categ == 0:
            for j in xrange(length):
                cigar_map.append(offset)
                offset += 1
        # Insertion
        elif categ == 1:
            cigar_map.extend([None] * length)
            #for j in xrange(length):
            #    cigar_map.append(None)
        # Deletion
        elif categ == 2:
            offset += length
        # Soft-clipping
        elif categ == 4:
            cigar_map.extend([None] * length)
            #for j in xrange(length):
            #    cigar_map.append(None)
        elif categ == 5:
            pass
        else:
            raise NotSupportedException('Cigar operation not supported: ' + str(categ))

    return cigar_map

def open_maybe_gzip(filename, mode='r'):
    if filename.endswith('.gz'):
        gunzip = subprocess.Popen(['gunzip', '-c', filename],
                                  stdout=subprocess.PIPE,
                                  preexec_fn=os.setsid)
        return gunzip.stdout
    else:
        return open(filename, mode)