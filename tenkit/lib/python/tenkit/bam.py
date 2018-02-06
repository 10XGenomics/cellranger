#!/usr/bin/env python
#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#
# Utilities for manipulating bam files
#
import os
import pysam
import log_subprocess
import random
import logging
import shutil
import math
import tenkit.bio_io as tk_io
import tenkit.seq as tk_seq
import csv
import striped_smith_waterman.ssw_wrap as ssw_wrap
import bisect
import os.path

DEFAULT_RG_STRING = 'None:None:None:None:None'

def pack_rg_string(sample_id, library_id, gem_group, flowcell, lane):
    '''Make string encoding read group info, for passing around in MROs'''
    if sample_id is None:
        return DEFAULT_RG_STRING
    else:
        return ':'.join(map(str, [sample_id, library_id, gem_group, flowcell, lane]))

def make_rg_header(packed_rg_string=DEFAULT_RG_STRING):
    '''Make the RG header, matching how it's done in Lariat.'''
    result = packed_rg_string.split(':')
    if len(result) != 5:
        raise Exception("RG string must have this format - sample_id:library_id:gem_group:flowcell:lane")
    sample_id, library_id, gem_group, flowcell, lane = result
    return '@RG\\tID:{0}\\tSM:{1}\\tLB:{2}.{3}\\tPU:{0}\\tPL:ILLUMINA'.format(packed_rg_string, sample_id, library_id, gem_group)

def make_star_rg_header(packed_rg_string):
    '''Make the RG header for STAR. Returns a list of strings.'''
    return make_rg_header(packed_rg_string).split('\\t')[1:]

def make_pg_header(version, stage_name):
    ''' Generate a BAM PG header line for the stage that is manipulating the bam file '''
    pg = {
        'ID': stage_name,
        'PN': "longranger." + stage_name,
        'VN': version
    }
    return pg

def make_terminal_pg_header(version):
    ''' Generate a PG header line that denotes the entire longranger pipeline, for tracking purposes.'''
    pg = {'ID': 'longranger', 'PN': 'longranger', 'VN': version}
    return pg

def get_number_reads(bam_filename):
    return reduce(lambda x, y: x + y, [ eval('+'.join(l.rstrip('\n').split('\t')[2:]) ) for l in pysam.idxstats(bam_filename) ])

# NOTE!! pgs argument used to be called pg and take a single pg item
# When resolving to tenkit:master, clients that inject PG tags must be updated to pass a list
def create_bam_outfile(file_name, chrom_names, chrom_lengths, template=None, pgs=None, cos=None, rgs=None, replace_rg=False):
    """ Creates a bam file with given chromosome names and lengths.
    template is an existing bam file.  If it is specified, chrom_names and chrom_lengths
    are ignored. pg is dictionary specifying a 'PG' entry. ID field is required; PN/CL/PP/DS/VN fields are optional.
    rgs is a list of dicts specifiying an 'RG' entry. If replace_rg is True, the existing 'RG' entry is overwritten.
    """
    if template:
        header = template.header
        if pgs is not None:
            for pg in pgs:
                if not header.has_key('PG'):
                    header['PG'] = []
                # add in the PP field based on previous PG entry
                if len(header['PG']) > 0:
                    pp = header['PG'][-1]['ID']
                    if pp is not None:
                        pg['PP'] = pp
                header['PG'].append(pg)

        if cos is not None:
            for co in cos:
                if not header.has_key('CO'):
                    header['CO'] = []
                header['CO'].append(co)

        if rgs is not None:
            if replace_rg and header.has_key('RG') and len(rgs) > 0:
                header['RG'] = []
            for rg in rgs:
                if not header.has_key('RG'):
                    header['RG'] = []
                header['RG'].append(rg)

        bam_file = pysam.Samfile(file_name, 'wb', header=header)
        tids = {name:n for (n, name) in enumerate(template.references)}
    else:
        header = {'SQ': [{'SN': chrom_names[n], 'LN': chrom_lengths[n]} for n in xrange(len(chrom_names))]}
        bam_file = pysam.Samfile(file_name, 'wb', header=header)
        tids = {chrom_names[n]:n for n in xrange(len(chrom_names))}
    return bam_file, tids

def create_bam_infile(file_name):
    bam_file = pysam.Samfile(file_name, 'rb')
    return bam_file

def create_sam_infile(file_name):
    sam_file = pysam.Samfile(file_name, 'r')
    return sam_file

def downsample_bam(in_name, out_name, downsample_rate, restrict_chrom=None):
	""" Downsamples a bam.  Optionally also restricts the output to a single chromosome.
	"""
	f = create_bam_infile(in_name)
	g, tids = create_bam_outfile(out_name, None, None, template=f)

	if restrict_chrom is None:
		bam_iter = f
	else:
		bam_iter = f.fetch(restrict_chrom)

	should_write = {}
	for r in bam_iter:
		if should_write.has_key(r.qname):
			if should_write[r.qname]:
				g.write(r)
		else:
			if random.random() < downsample_rate:
				should_write[r.qname] = True
				g.write(r)
			else:
				should_write[r.qname] = False
	g.close()


def filter_bam(bam_iter, remove_secondary=True, min_mapq=0, remove_unmapped = False, read_filter=lambda x: True):
    for read in bam_iter:
        if remove_secondary and read.is_secondary:
            continue

        if read.mapq < min_mapq:
            continue

        if remove_unmapped and read.is_unmapped:
            continue

        if not read_filter(read):
            continue

        yield read


def write_read(bam_file, name, seq, qual, tid, pos, mapq=0, cigar=None, reverse=False, paired=False, read1=True,
               mate_tid=None, mate_pos=None, mate_reverse=False):
    """ Creates an alignedRead object and writes it to the corresponding bam file.
    If cigar is unspecified, the alignment is defaulted to a perfect alignment.
    Chrom
    """
    r = pysam.AlignedRead()
    r.qname = name
    r.seq = seq
    r.qual = qual
    r.tid = tid
    r.pos = pos
    r.mapq = mapq

    if pos == -1:
        r.is_unmapped=True

    if cigar:
        r.cigar = cigar
    else:
        r.cigar = [(0, len(seq))]

    r.is_reverse = reverse
    r.is_read1 = read1
    r.is_read2 = not(read1)

    if paired:
        r.is_paired = True

        if mate_tid:
            r.mrnm = mate_tid

        if mate_pos:
            r.mpos = mate_pos

        if mate_pos == -1:
            r.mate_is_unmapped=True

        if mate_reverse:
            r.mate_is_reverse = True
        else:
            r.mate_is_reverse = False

    bam_file.write(r)

def write_read_pair(bam_file, name, seq1, qual1, tid1, pos1, seq2, qual2, tid2, pos2, mapq1=0, mapq2=0, cigar1=None, cigar2=None,
                    reverse1=False, reverse2=True):
    """Writes both members of a read pair.  This is enforced by the two reads sharing the same name.
    """
    write_read(bam_file, name, seq1, qual1, tid1, pos1, mapq=mapq1, cigar=cigar1, reverse=reverse1, paired=True, read1=True,
               mate_tid=tid2, mate_pos=pos2, mate_reverse=reverse2)
    write_read(bam_file, name, seq2, qual2, tid2, pos2, mapq=mapq2, cigar=cigar2, reverse=reverse2, paired=True, read1=False,
              mate_tid=tid1, mate_pos=pos1, mate_reverse=reverse1)

def sort(file_name, sorted_prefix=None):
    """ Sorts and indexes the bam file given by file_name.
    """
    if sorted_prefix is None:
        sorted_prefix = file_name.replace('.bam', '') + '_sorted'

    sorted_name = sorted_prefix + ".bam"
    log_subprocess.check_call(['samtools','sort', '-o', sorted_name, file_name])

def index(file_name):
    pysam.index(str(file_name))
    if not os.path.isfile(file_name + ".bai"):
        raise RuntimeError("samtools index failed, likely bam is not sorted")

def sort_and_index(file_name, sorted_prefix=None):
    """ Sorts and indexes the bam file given by file_name.
    """
    if sorted_prefix is None:
        sorted_prefix = file_name.replace('.bam', '') + '_sorted'

    sorted_name = sorted_prefix + '.bam'
    log_subprocess.check_call(['samtools','sort', '-o', sorted_name, file_name])
    pysam.index(sorted_name)

def merge(out_file_name, input_file_names, threads = 1):
    # Note the original samtools merge call can
    # fail if the total length of the command line
    # gets too long -- use the API version instead.
    #args = ['samtools', 'merge', out_file_name]
    #args.extend(input_file_names)
    #log_subprocess.check_call(args)

    if threads > 1:
        args = ["-c", "-p", "-s", "0", "-@", str(threads)]
    else:
        args = []

    args.append(str(out_file_name))
    args.extend([str(x) for x in input_file_names])
    pysam.merge(*args)

def bam_is_empty(fn):
    if os.path.getsize(fn) > 1000000:
        return False

    bam = pysam.Samfile(fn, check_sq=False)
    try:
        bam.next()
        return False
    except StopIteration:
        return True

def concatenate(out_file_name, all_in_file_names):
    """ Concatenate a list of bam files into a final output file """

    # Filter out empty BAM files -- these cause samtools cat to generate
    # a BAM with a premature end block
    in_file_names = [f for f in all_in_file_names if not bam_is_empty(f)]

    if len(in_file_names) > 1:
        args = ['samtools', 'cat', '-o', out_file_name]
        args.extend(in_file_names)
        log_subprocess.check_call(args)
    elif len(in_file_names) == 0:
        # If all the BAMs are empty, just copy 1 over
        shutil.copy(all_in_file_names[0], out_file_name)
    else:
        shutil.copy(in_file_names[0], out_file_name)

def merge_by_name(out_file_name, in_file_names):
    """ Merge name-sorted bam files into bam file sorted by name"""
    args = ['samtools', 'merge', '-n', out_file_name]
    args.extend(in_file_names)
    log_subprocess.check_call(args)

def remove_dups(in_name, out_name):
    """ remove paired-end duplicates using samtools
    """
    log_subprocess.check_call(['samtools', 'rmdup', in_name, out_name])

def sort_by_name(file_name, sorted_prefix=None):
    """ Sorts a bam file by the read name, for paired-end
    """
    if sorted_prefix is None:
        sorted_prefix = file_name.replace('.bam', '') + '_namesorted'

    sorted_name = sorted_prefix + '.bam'
    # NOTE -- need to update our internal samtools in order to use pysam.sort
    #pysam.sort('-n', file_name, sorted_prefix)
    log_subprocess.check_call(['samtools', 'sort', '-n', file_name, sorted_prefix])

    return pysam.Samfile(sorted_name, 'rb')

def sort_by_bc(file_name, sorted_name, store_in_memory=True):
    """ Sorts a bam file by the 10X barcode (specified in the tags BC field)
    if store_in_memory is True, avoids file seeks by keeping reads in memory
    """
    in_file = create_bam_infile(file_name)
    out_file, tids = create_bam_outfile(sorted_name, None, None, template=in_file)

    if store_in_memory:
        bc_reads = {}

        for read in in_file:
            bc = tk_io.get_read_barcode(read)
            this_bc_reads = bc_reads.setdefault(bc, [])
            this_bc_reads.append(read)

        sorted_bcs = sorted(bc_reads.keys())
        for bc in sorted_bcs:
            for read in bc_reads[bc]:
                out_file.write(read)
    else:
        # Store the file offset locations (in bytes) by bc
        bc_locs = {}

        file_offset = in_file.tell()
        for read in in_file:
            bc = tk_io.get_read_barcode(read)
            this_bc_locs = bc_locs.setdefault(bc, [])
            this_bc_locs.append(file_offset)
            file_offset = in_file.tell() # has to be before next read

        sorted_bcs = sorted(bc_locs.keys())
        for bc in sorted_bcs:
            for offset in bc_locs[bc]:
                in_file.seek(offset)
                out_file.write(in_file.next())

    out_file.close()

def convert_to_bam(sam_name, bam_name):
    """ Uses samtools to create a bam_file from the samfile
    """
    sam_file = create_sam_infile(sam_name)
    bam_file, tids = create_bam_outfile(bam_name, None, None, template=sam_file)
    for read in sam_file:
        bam_file.write(read)
    sam_file.close()
    bam_file.close()

def generate_tiling_windows(input_bam, locus_size, overlap=0):
    ''' Generate a list of (chrom, start, length) loci that tile over all the references in the bam file '''

    chroms = input_bam.references
    chrom_lengths = input_bam.lengths

    loci = []
    for (chrom, length) in zip(chroms, chrom_lengths):
        start = 0
        while start + locus_size + overlap < length:
            stop = start + locus_size + overlap
            loci.append(tk_io.create_locus_info(chrom, start, stop))
            start += locus_size
        loci.append(tk_io.create_locus_info(chrom, start, length))

    return loci

import struct

def get_bsize(f, subfield_start, extra_len):
    (id1, id2, length) = struct.unpack("BBH", f.read(4))

    if id1 == 66 and id2 == 67 and length == 2:
        bsize = struct.unpack("H", f.read(2))
        return bsize[0]
    #elif f.tell() - subfield_start < extra_len:
    #    get_bsize(f, subfield_start, extra_len)
    else:
        raise Exception("BSIZE field not found -- this a valid BAM file?")

def parse_bgzf_header(f):
    cur_pos = f.tell()
    header_fmt = "BBBBIBBH"

    d = f.read(12)

    # We are at EOF when read returns an empty string
    if d == '':
        return None

    header = struct.unpack(header_fmt, d)

    # Check for a valid gzip header
    if header[0] != 31 or header[1] != 139:
        raise Exception("Not a valid gzip header")

    xlen = header[7]
    bsize = get_bsize(f, f.tell(), xlen)

    next_pos = cur_pos + bsize + 1
    f.seek(next_pos)
    return next_pos


def chunk_bam_records(input_bam, chunk_bound_key=None, chunk_size_gb=0.75,
    num_chunks=None, min_chunks=1, max_chunks=None):
    ''' Find a series BAM virtual offset ranges that cover all the reads in the BAM file.
        Each chunk will be cover at least 1 BGZF chunk, so the total number of chunks will
        be the smaller of max_chunks, and the total number of BGZF chunks in the file.
        chunk_bound_key is a function that accepts a pysam AlignedRead object and returns
        some read group identifier object.  BAM chunks returned will always start on the
        first read of a block of reads with a common group identifier.  Use this to keep
        read pairs, or blocks of reads with a common start position together when chunking.
        Set chunk_bound_key=None if there is no constraint on chunk boundaries.
        chunk_bound key can return None to signal that the current read is a valid place to
        start a chunk.
    '''
    if num_chunks is None:
        size_gb = float(os.path.getsize(input_bam.filename)) / 1e9
        num_chunks = max(min_chunks, int(math.ceil(size_gb / chunk_size_gb)))
    if max_chunks is not None:
        num_chunks = min(num_chunks, max_chunks)

    input_bam.reset()
    first_chunk = input_bam.tell() >> 16

    # Find BGZF block list
    fp = file(input_bam.filename, "rb")
    chunk_offsets = [ x for x in iter(lambda: parse_bgzf_header(fp), None) if x >= first_chunk ]

    def find_valid_virtual_offset(pos):
        ''' Scan through reads starting at given pos to find valid start site for iteration '''
        input_bam.seek(pos << 16)
        n = 0
        pos = input_bam.tell()
        last_read_key = None

        # All positions are valid chunk bounds in chunk_bound_key is None
        # unless there are no reads are left
        if chunk_bound_key is None:
            has_reads = False
            for r in input_bam:
                has_reads = True
                break
            if has_reads:
                return pos
            else:
                return None

        for r in input_bam:
            # Warn if it's taking a long time to find a good start position
            n += 1
            if n==1000:
                logging.warning("Taking a long time to find a chunk boundary. Check your chunk_boundary_test function?")

            # get the current bound_key -- a valid block begins when this key changes.
            new_key = chunk_bound_key(r)

            # The chunk_bound_key can signal that it's a valid place to split by returning None
            if new_key is None:
                return pos

            if last_read_key == None:
                last_read_key = new_key
                pos = input_bam.tell()
            elif new_key == last_read_key:
                pos = input_bam.tell()
            else:
                return pos

        # We will trim off chunks that don't have a valid
        # start position, or have no reads.
        # Any reads will get processed by the preceding chunks
        return None

    # The first read is always a valid start point
    # For subsequent chunks we iterate until the chunk_boundary_test returns True
    start_points = range(0, len(chunk_offsets), max(1, len(chunk_offsets)/num_chunks))
    start_virtual_offsets = [find_valid_virtual_offset(chunk_offsets[x]) for x in start_points[1:]]
    start_virtual_offsets = [x for x in start_virtual_offsets if x is not None]

    # Remove duplicate start sites
    start_virtual_offsets = sorted(list(set(start_virtual_offsets)))
    start_virtual_offsets.insert(0, first_chunk  << 16)

    end_virtual_offsets = start_virtual_offsets[1:]
    end_virtual_offsets.append(None)

    chunk_defs = list(zip(start_virtual_offsets, end_virtual_offsets))

    return [{ 'chunk_start': str(s), 'chunk_end': str(e) } for (s,e) in chunk_defs]

def read_bam_chunk(input_bam, chunk_def):
    ''' Iterate over the reads in input_bam contained in chunk_def.  chunk_def is a (start, end) pair
        of virtual offsets into the BAM file '''

    start, end = chunk_def

    def convert_offset(s):
        if s == 'None':
            return None
        if isinstance(s, str) or isinstance(s, unicode):
            return int(s)

        return s

    start = convert_offset(start)
    end = convert_offset(end)

    if start:
        input_bam.seek(start)
    else:
        input_bam.reset()

    while end == None or input_bam.tell() < end:
        yield input_bam.next()

def get_allele_read_info(chrom, pos, ref, alt_alleles, min_mapq_counts, min_mapq_for_mean, min_mapq_for_bc, default_indel_qual, bam, reference_pyfasta, max_reads=1000, match = 1, mismatch = -4, gap_open = -6, gap_extend = -1):
    all_alleles = [ref] + alt_alleles
    bc_qual_maps = [{} for j in xrange(len(all_alleles))]
    counts = [0 for x in all_alleles]
    diffs = [[] for x in all_alleles]
    mapq_sums = [0.0 for x in all_alleles]
    mapq_denoms = [0.0 for x in all_alleles]
    molecule_differences = [[] for x in all_alleles]
    rescued = [[] for x in all_alleles]
    num_reads = 0
    qnames = set()
    for read in bam.fetch(chrom, pos, pos + 1):
        num_reads += 1
        if read.qname in qnames:
            continue
        qnames.add(read.qname)
        if read.is_duplicate:
            continue
        if num_reads > max_reads:
            break

        is_indel_variant = False
        for allele in alt_alleles:
            if len(allele) != len(ref):
                is_indel_variant = True

        allele_index_in_read = read_contains_allele_sw(ref, all_alleles, pos, read, reference_pyfasta[chrom], match = match, mismatch = mismatch, gap_open = gap_open, gap_extend = gap_extend)
        for (allele_index, allele) in enumerate(all_alleles):
            if allele_index == allele_index_in_read:
                if dict(read.tags).get("AS") is not None and dict(read.tags).get("XS") is not None:
                    diffs[allele_index].append(float(dict(read.tags).get("AS")) - float(dict(read.tags).get("XS")))
                if dict(read.tags).get('OM') is not None:
                    if read.mapq >= 30 and dict(read.tags).get('OM') < 30:
                        rescue = 1
                    else:
                        rescue = 0
                    rescued[allele_index].append(rescue)
                if dict(read.tags).get("DM") is not None:
                    molecule_differences[allele_index].append(float(dict(read.tags).get("DM")))
                if read.mapq >= min_mapq_for_mean:
                    mapq_sums[allele_index] += read.mapq
                    mapq_denoms[allele_index] += 1
                if read.mapq >= min_mapq_counts:
                    counts[allele_index] += 1
                if read.mapq >= min_mapq_for_bc:
                    bc = tk_io.get_read_barcode(read)
                    if bc is None:
                        continue
                    cigar_map = tk_seq.get_cigar_map(read.cigar)
                    try:
                        read_offset = cigar_map.index(pos - read.pos - 1)
                    except:
                        continue
                    if allele == ref:
                        if is_indel_variant:
                            qual = str(default_indel_qual)
                        else:
                            qual = str(ord(min(read.qual[read_offset:read_offset + len(allele)])))
                        bc_quals = bc_qual_maps[allele_index].setdefault(bc, [])
                        bc_quals.append(qual)
                    # SNP
                    elif len(allele) == 1 and len(ref) == 1:
                        if is_indel_variant:
                            qual = str(default_indel_qual)
                        else:
                            qual = str(ord(read.qual[read_offset]))
                        bc_quals = bc_qual_maps[allele_index].setdefault(bc, [])
                        bc_quals.append(qual)
                    # Insert
                    elif len(allele) > len(ref) and allele.startswith(ref):
                        bc_quals = bc_qual_maps[allele_index].setdefault(bc, [])
                        bc_quals.append(str(default_indel_qual))
                    # Deletion
                    elif len(allele) < len(ref) and ref.startswith(allele):
                        bc_quals = bc_qual_maps[allele_index].setdefault(bc, [])
                        bc_quals.append(str(default_indel_qual))
                    else:
                        bc_quals = bc_qual_maps[allele_index].setdefault(bc, [])
                        bc_quals.append(str(default_indel_qual))
    bc_qual_strings = []
    for bc_qual_map in bc_qual_maps:
        bc_qual_strings.append([])
        for bc, bc_quals in bc_qual_map.iteritems():
            bc_qual_strings[-1].append(bc + '_' + '_'.join(bc_quals))


    mapq_means = [mapq_sums[i]/mapq_denoms[i] if mapq_denoms[i] > 0 else 31 for i in range(len(all_alleles))]
    return (counts, mapq_means, bc_qual_strings, molecule_differences, diffs, rescued)

class RevTuple(object):
    def __init__(self, value):
        self.value = value
    def __cmp__(self, other):
        return cmp(self.value,other[1])

def read_contains_allele_sw(ref, alleles, pos, read, reference_pyfasta, match = 1, mismatch = -1, gap_open = -6, gap_extend = -1):
    '''use a Smith-Waterman alignment to determine which allele best matches a read'''
    if read.is_unmapped or read.aend is None:
        return -1

    buf = len(ref)

    for allele in alleles:
        buf = max(buf, len(allele))
    buf += 50

    seq = read.seq
    refstart = max(0, read.pos - buf)
    refend = min(read.aend + buf, len(reference_pyfasta))

    if len(seq) < 500:
        refstart = max(0, read.pos)
        refend = min(read.aend, len(reference_pyfasta))
        readstart = 0
        readend = len(seq)
    else:
        aligned_pairs = read.aligned_pairs
        index = bisect.bisect_left(aligned_pairs, RevTuple(pos))
        startdex = max(0, index - 75)
        (readstart, refstart) = aligned_pairs[startdex]

        while startdex > 0 and (readstart == None or refstart == None):
            startdex -= 1
            (readstart, refstart) = aligned_pairs[startdex]
        while startdex < len(aligned_pairs) - 1 and (readstart == None or refstart == None):
            startdex += 1
            (readstart, refstart) = aligned_pairs[startdex]
        if readstart == None or refstart == None:
            return -1
        enddex = min(len(aligned_pairs)-1, index + 75)
        (readend, refend) = aligned_pairs[enddex]
        while enddex < len(aligned_pairs) - 1 and (readend == None or refend == None):
            enddex += 1
            (readend, refend) = aligned_pairs[enddex]
        while enddex > 0 and  (readend == None or refend == None):
            enddex -= 1
            (readend, refend) = aligned_pairs[enddex]
        if (readend == None or refend == None):
            return -1
        if enddex - startdex < 75:
            return -1

    # set up alternate reference sequences to have the alt alleles
    def setup_ref(ref, allele, pos, reference_pyfasta, refstart, refend):
        ref = reference_pyfasta[max(0,refstart-buf):(pos-1)] + allele + reference_pyfasta[(pos-1)+len(ref):min(refend+buf, len(reference_pyfasta))]
        return ref.upper()
    refs = [setup_ref(ref, allele, pos, reference_pyfasta, refstart, refend) for allele in alleles]

    best_match = -1
    best_score = -1
    best_score_equal = True
    subseq = seq[readstart:readend]

    for (index, ref) in enumerate(refs):
        # NOTE: changing gap_extend does not have any effect! gap_open controls the linear gap penalty, but it seems that there is no way to control the constant.
        aligner = ssw_wrap.Aligner(ref, match = match, mismatch = -mismatch, gap_open = -gap_open, gap_extend = -gap_extend, report_secondary=False)
        aln = aligner.align(subseq, min_score = 25, min_len = 25)
        if aln is None:
            continue

        if aln.score == best_score:
            best_score_equal = True
        if aln.score > best_score:
            best_score_equal = False
            best_score = aln.score
            best_match = index

    if best_score_equal:
        return -1

    return best_match

def pretty_print_alignment(s, r , cig, rpos):
    # this is for debug use with the striped smith waterman alignments. But can
    # be used with regular read alignments
    cig_type= []
    cig_length = []
    building_int = ''
    for c in cig:
        if c in ['0','1','2','3','4','5','6','7','8','9']:
            building_int += c
        else:
            cig_type.append(c)
            cig_length.append(int(building_int))
            building_int = ''
    sloc = 0
    rloc = rpos
    new_s = ''
    new_r = r[0:rloc]
    match_s = ''
    for i in range(rloc):
        match_s += ' '
        new_s += ' '
    matches = 0
    mismatches = 0
    insertion_starts = 0
    deletion_starts = 0
    insertion_extends = 0
    deletion_extends = 0
    for (ctype, length) in zip(cig_type,cig_length):
        if ctype == 'D':
            for i in range(length):
                new_s +=  ' '
                new_r += r[rloc]
                rloc += 1
                match_s += ' '
            deletion_starts += 1
            deletion_extends += length - 1
        elif ctype == 'I':
            for i in range(length):
                new_s += s[sloc]
                sloc += 1
                new_r += ' '
                match_s += ' '
            insertion_starts += 1
            insertion_extends += length - 1
        elif ctype == 'M':
            for i in range(length):
                if s[sloc] == r[rloc]:
                    match_s += '|'
                    matches += 1
                else:
                    match_s += '*'
                    mismatches += 1
                new_s += s[sloc]
                sloc += 1
                new_r += r[rloc]
                rloc +=1
        elif ctype == "S":
            for i in range(length):
                sloc += 1
                rloc += 1
                match_s += ' '

    new_r += r[rloc:len(r)]
    print new_r
    print match_s
    print new_s
    print str(matches)+ " matches "+str(mismatches)+ " mismatches " +str(insertion_starts)+" insertion_starts "+str(deletion_starts)+ " deletion_starts " + str(insertion_extends)+" insertion_extends "+ str(deletion_extends)+" deletion_extends"

def read_contains_allele(ref, allele, pos, read, reference_pyfasta, cigar, cigar_map):
    if allele != ref:
        (pos, ref, allele) = tk_io.normalize_variant(pos, ref, allele, reference_pyfasta)
    # reference allele
    if allele == ref:
        #cigar_map = tk_seq.get_cigar_map(read.cigar)
        try:
            read_offset = cigar_map.index(pos - read.pos - 1)
            next_read_offset = cigar_map.index(pos - read.pos)
        except:
            return False

        nucs = read.seq[read_offset:read_offset + len(allele)]

        if nucs == allele and next_read_offset == read_offset + 1:
            return True
    # SNP
    elif len(allele) == len(ref) and len(ref) == 1:
        #cigar_map = tk_seq.get_cigar_map(read.cigar)
        try:
            read_offset = cigar_map.index(pos - read.pos - 1)
            next_read_offset = cigar_map.index(pos - read.pos)
        except:
            return False

        nuc = read.seq[read_offset]

        if nuc == allele and next_read_offset == read_offset + 1:
            return True
    # Insert
    elif len(allele) > len(ref) and allele.startswith(ref):
        called_insert_seq = allele[len(ref):]
        chrom_pos = read.pos
        offset = 0
        for (categ, length) in cigar: #read.cigar:
            # Match
            if categ == 0:
                chrom_pos += length
                offset += length
            # Insertion
            elif categ == 1:
                ins_seq = read.seq[offset:offset+length]
                if chrom_pos == pos and ins_seq == called_insert_seq:
                    return True
                offset += length
            # Deletion
            elif categ == 2:
                del_length = length
                chrom_pos += length
            # Soft-clipping
            elif categ == 4:
                offset += length
            elif categ == 5:
                return False
            else:
                raise NotImplementedError('Cigar operation not supported: ' + str(categ))
    # Deletion
    elif len(allele) < len(ref) and ref.startswith(allele):
        called_del_length = len(ref) - len(allele)
        chrom_pos = read.pos
        offset = 0
        for (categ, length) in cigar: #read.cigar:
            # Match
            if categ == 0:
                chrom_pos += length
                offset += length
            # Insertion
            elif categ == 1:
                offset += length
            # Deletion
            elif categ == 2:
                del_length = length
                if chrom_pos == pos and del_length == called_del_length:
                    return True
                chrom_pos += length
            # Soft-clipping
            elif categ == 4:
                offset += length
            elif categ == 5:
                return False # skip reads w/ hard clipping
            else:
                raise NotImplementedError('Cigar operation not supported: ' + str(categ))
    else:
        # get the portion of the cigar in the allele + max(0, len(ref)-len(allele))
        read_allele_region_seq = ""
        chrom_pos = read.pos
        offset = 0
        deletions_in_allele = 0
        insertions_in_allele = 0
        positions_crossed = 0
        offset_end = -2
        base_after_allele_is_match = False # if the base after the allele isn't a match, it could be a different allele
        for (categ, length) in cigar: # read.cigar:
            # Match
            if categ == 0:
                for i in range(length):
                    if chrom_pos + i >= pos - 1 and len(read_allele_region_seq) < len(allele) and len(read.seq) > offset + i:
                        read_allele_region_seq += read.seq[offset+i]
                        positions_crossed += 1

                        if len(read_allele_region_seq) == len(allele):
                            offset_end = offset+i+1
                    # base after allele must be a match
                    if offset + i == offset_end and positions_crossed == max(len(ref), len(allele)):
                        base_after_allele_is_match = True
                        break
                chrom_pos += length
                offset += length
            # Insertion
            elif categ == 1:
                for i in range(length):
                    if chrom_pos + i >= pos - 1 and len(read_allele_region_seq) < len(allele) and len(read.seq) > offset + i:
                        read_allele_region_seq += read.seq[offset+i]
                        positions_crossed += 1
                        insertions_in_allele += 1
                        if len(read_allele_region_seq) == len(allele):
                            offset_end = offset +i + 1
                offset += length
            # Deletion
            elif categ == 2:
                for i in range(length):
                    if chrom_pos + i >= pos - 1 and len(read_allele_region_seq) < len(allele) and len(read.seq) > offset + i:
                        positions_crossed += 1
                        deletions_in_allele += 1
                del_length = length
                chrom_pos += length
            # Soft-clipping
            elif categ == 4:
                offset += length
            elif categ == 5:
                return False # skip reads w/ hard clipping
            else:
                raise NotImplementedError('Cigar operation not supported: ' + str(categ))
        if read_allele_region_seq == allele and base_after_allele_is_match and len(ref) + insertions_in_allele - deletions_in_allele == len(allele):
            return True
    return False

class BamIndex:
    def __init__(self, file_name, file_suffix, key_func, cmp_func):
        self.file_name = file_name
        self.file_suffix = file_suffix
        self.key_func = key_func
        self.cmp_func = cmp_func
        self.index_file_name = file_name + '.' + file_suffix
        self.in_bam = create_bam_infile(file_name)

    def save_key_index(self):
        f = open(self.index_file_name, 'w')
        out_index = csv.writer(f, delimiter='\t')

        self.in_bam.reset()

        last_key = None
        last_pos = self.in_bam.tell()
        for i, read in enumerate(self.in_bam):
            if i == 0:
                last_key = self.key_func(read)
                out_index.writerow((last_key, last_pos))
            else:
                key = self.key_func(read)
                cmp_value = self.cmp_func(key, last_key)
                if cmp_value < 0:
                    print read
                    print "idx: %i, key: %s, last_key: %s" % (i, str(key), str(last_key))
                    raise Exception("BAM file %s is not sorted by key" % self.file_name)
                elif cmp_value != 0:
                    out_index.writerow((key, last_pos))
                    last_key = key

            last_pos = self.in_bam.tell()

        f.close()

    def load_key_index(self):
        f = open(self.index_file_name, 'r')
        in_index = csv.reader(f, delimiter='\t')

        self.keys = []
        for row in in_index:
            if len(row) != 2:
                raise Exception("BAM index file %s has incorrect format" % self.index_file_name)
            key, pos = row
            if not pos.isdigit():
                raise Exception("BAM index file %s has incorrect format" % self.index_file_name)
            self.keys.append((key, long(pos)))

        f.close()

    def get_reads_iter_with_key(self, key):
        # Binary search through keys
        low, high = 0, len(self.keys)
        while low < high:
            mid = (low + high) / 2
            mid_key, _ = self.keys[mid]

            cmp_value = self.cmp_func(mid_key, key)
            if cmp_value < 0:
                low = mid + 1
            elif cmp_value > 0:
                high = mid
            else:
                low, high = mid, mid
                break

        _, pos = self.keys[low]
        self.in_bam.seek(pos)

        # Iterate through reads to find reads with associated key
        for read in self.in_bam:
            read_key = self.key_func(read)

            cmp_value = self.cmp_func(key, read_key)
            if cmp_value == 0:
                yield read
            elif cmp_value < 0:
                break

def qname_key_func(read):
    return read.qname

# Ported from strnum_cmp in samtools bam_sort.c
def qname_cmp_func(qname1, qname2):

    def update_char(qname, index):
        index += 1
        if index < len(qname):
            c = qname[index]
        else:
            c = ''
        return c, index

    def get_ord(c):
        if c:
            return ord(c)
        return 0

    i, j = 0, 0
    c1 = qname1[0] if qname1 else None
    c2 = qname2[0] if qname2 else None
    while c1 and c2:
        if c1.isdigit() and c2.isdigit():
            while c1 == '0':
                c1, i = update_char(qname1, i)
            while c2 == '0':
                c2, j = update_char(qname2, j)
            while c1.isdigit() and c2.isdigit() and c1 == c2:
                c1, i = update_char(qname1, i)
                c2, j = update_char(qname2, j)
            if c1.isdigit() and c2.isdigit():
                k, l = i, j
                while c1.isdigit() and c2.isdigit():
                    c1, k = update_char(qname1, k)
                    c2, l = update_char(qname2, l)
                return 1 if c1.isdigit() else (-1 if c2.isdigit() else get_ord(qname1[i]) - get_ord(qname2[j]))
            elif c1.isdigit():
                return 1
            elif c2.isdigit():
                return -1
            elif i != j:
                return 1 if i < j else -1
        else:
            if c1 != c2:
                return get_ord(c1) - get_ord(c2)
            c1, i = update_char(qname1, i)
            c2, j = update_char(qname2, j)
    return 1 if c1 else (-1 if c2 else 0)

class BamQnameIndex(BamIndex):
    def __init__(self, file_name):
        return BamIndex.__init__(self, file_name, 'qname', qname_key_func, qname_cmp_func)

    def save_qname_index(self):
        BamIndex.save_key_index(self)

    def load_qname_index(self):
        BamIndex.load_key_index(self)

    def get_reads_iter_with_qname(self, qname):
        return BamIndex.get_reads_iter_with_key(self, qname)


def bc_key_func(read):
    return tk_io.get_read_barcode(read)

def bc_cmp_func(bc1, bc2):
    if bc1 is None and bc2 is None:
        return 0
    if bc1 < bc2:
        return -1
    if bc1 > bc2:
        return 1
    return 0

class BamBCIndex(BamIndex):
    def __init__(self, file_name):
        return BamIndex.__init__(self, file_name, 'bxi', bc_key_func, bc_cmp_func)

    def save_index(self):
        BamIndex.save_key_index(self)

    def load_index(self):
        BamIndex.load_key_index(self)

    def get_reads_bc_iter(self, bc):
        return BamIndex.get_reads_iter_with_key(self, bc)


def get_insert_size(read):
    if read.tid == read.mrnm:
        return abs(read.pos - read.mpos) + len(read.seq)
    else:
        return float('Inf')
