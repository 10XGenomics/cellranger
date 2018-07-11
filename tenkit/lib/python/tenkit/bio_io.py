#!/usr/bin/env python
#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#

import re
import os
import glob
import vcf
import subprocess
from collections import OrderedDict, namedtuple
import tenkit.tabix as tk_tabix
from tenkit.regions import Regions
from tenkit.exceptions import NotSupportedException
from tenkit.constants import (PROCESSED_BARCODE_TAG, RAW_BARCODE_TAG, \
                              RAW_BARCODE_QUAL_TAG, SAMPLE_INDEX_TAG, \
                              SAMPLE_INDEX_QUAL_TAG, HAPLOTYPE_BAM_TAG, \
                              PHASE_SET_BAM_TAG)


class VariantFileReader():
    def __init__(self, file_name, file_type='VCF4.1'):
        self.file_type = file_type
        self.file_name = file_name

    def record_getter(self, restrict_type=None, fetch_chrom=None, fetch_start=None, fetch_end=None):
        """ Iterator that just returns raw records.  Useful for writing to output
        directly w/o doing additional parsing (such as when combining multiple files).
        Can be called repeatedly (i.e. re-opens the file at the instantiation of the iterator)
        """
        if not(fetch_chrom is None):
            assert(not(fetch_start is None))
            assert(not(fetch_end is None))

        if self.file_type == 'VCF4.1':
            with open(self.file_name, 'r') as in_file:

                if fetch_chrom is None:
                    record_iter = vcf.Reader(in_file)
                else:
                    try:
                        record_iter = vcf.Reader(in_file).fetch(fetch_chrom, fetch_start, end=fetch_end)
                    except KeyError:
                        # Return empty result if the chrom isn't found?
                        record_iter = []
                    except ValueError:
                        # Exception type changed to ValueError in pyvcf 0.6.8
                        # Return empty result if the chrom isn't found?
                        record_iter = []

                for record in record_iter:
                    if not(restrict_type is None):
                        var_type = record.var_type
                        if not(var_type == restrict_type):
                           continue
                    yield record


def get_var_type(ref, alt):
    """ Classifies everything besides simple SUB, Insertion or Deletion as COMPLEX
    """
    if alt is None or ref == alt:
        assert(not(ref is None))
        return 'X'
    if len(ref) == 1 and len(alt) == 1:
        return 'S'
    elif len(ref) > len(alt) and ref.startswith(alt):
        return 'D'
    elif len(ref) < len(alt) and alt.startswith(ref):
        return 'I'
    else:
        return 'C'

def normalize_variant(pos, ref, alt, reference_pyfasta):
    new_ref = ref
    new_alt = alt
    new_pos = pos
    change = True
    while change:
        old_ref = new_ref
        old_alt = new_alt
        old_pos = new_pos
        if not((new_ref == "" or new_alt == "")) and new_ref[-1] == new_alt[-1]:
            new_ref = new_ref[0:len(new_ref)-1]
            new_alt = new_alt[0:len(new_alt)-1]
        if (new_ref == "" or new_alt == "") and new_pos > 1:
            reference_base = reference_pyfasta[new_pos-2].upper()
            new_ref = reference_base + new_ref
            new_alt = reference_base + new_alt
            new_pos -= 1
        change = (new_ref != old_ref) or (new_alt != old_alt) or (new_pos != old_pos)
    if len(new_ref) == 0 or len(new_alt) == 0:
        return (old_pos, old_ref, old_alt)
    while new_ref[0] == new_alt[0] and (len(new_ref) > 1 and len(new_alt) > 1):
        new_ref = new_ref[1:]
        new_alt = new_alt[1:]
        new_pos += 1
    to_return = (new_pos, new_ref, new_alt)
    return to_return

def get_variant_iterator(vars_file, targets_file_name, restrict_locus, edge_case_filter):
    if targets_file_name is not None:
        with open(targets_file_name) as targets:
            regions = get_target_regions(targets)
    else:
        regions = None

    if restrict_locus is not None:
        chrom, start, end = get_locus_info(restrict_locus)
        # NB: add 1 to fetch_start to switch back to 1-indexing in VCF
        for record in vars_file.record_getter(fetch_chrom=chrom, fetch_start=start+1, fetch_end=end):
            if regions is None or edge_case_filter(record, regions):
                yield record
    else:
        for record in vars_file.record_getter():
            if regions is None or edge_case_filter(record, regions):
                yield record

def get_variant_iterator_inclusive(vars_file, targets_file_name, restrict_locus):
    def inclusive(record, regions):
        chrom = get_record_chrom(record)
        pos = get_record_pos(record) - 1 #put into 0 indexing for region comparison
        record_length = get_record_max_length(record)
        if chrom in regions:
            return regions[chrom].contains_point(pos) or regions[chrom].contains_point(pos+record_length)
        else:
            return False
    return get_variant_iterator(vars_file, targets_file_name, restrict_locus, inclusive)

def get_variant_iterator_pos(vars_file, targets_file_name, restrict_locus):
    def pos(record, regions):
        chrom = get_record_chrom(record)
        pos = get_record_pos(record) - 1 #put into 0 indexing for region comparison
        if chrom in regions:
            return regions[chrom].contains_point(pos)
    return get_variant_iterator(vars_file, targets_file_name, restrict_locus, pos)

def get_variant_iterator_exclusive(vars_file, targets_file_name, restrict_locus):
    def exclusive(record, regions):
        chrom = get_record_chrom(record)
        pos = get_record_pos(record) - 1 #put into 0 indexing for region comparison
        record_length = get_record_max_length(record)
        if chrom in regions:
            if regions[chrom].contains_point(pos):
                start, stop = regions[chrom].get_region_containing_point(pos)
                return pos+record_length < stop
        else:
            return False
    return get_variant_iterator(vars_file, targets_file_name, restrict_locus, exclusive)

def get_record_max_length(record):
    alts = get_record_alt_alleles(record)
    ref = get_record_ref(record)
    record_length = abs(min(get_allele_length(ref, alts[0]),0))
    if len(alts) > 1:
        record_length = abs(min(get_allele_length(ref, alts[0]), get_allele_length(ref, alts[1])))
    return record_length

class VariantFileWriter():
    def __init__(self, out_file, file_type='VCF4.1', template_file=None, template_reader=None, new_source=None,
        new_info_fields=[], new_format_fields=[], new_filters=[]):
        self.file_type = file_type
        if self.file_type == 'VCF4.1':
            if template_reader is None and template_file is not None:
                template_reader = vcf.Reader(template_file)
            elif template_reader is not None:
                pass
            else:
                metadata = OrderedDict()
                infos = OrderedDict()
                formats = OrderedDict()
                filters = OrderedDict()
                alts = OrderedDict()
                contigs = OrderedDict()
                template_reader = namedtuple('template', ['metadata', 'infos', 'formats', 'filters', 'alts', 'contigs'])
                template_reader.metadata = metadata
                template_reader.infos = infos
                template_reader.formats = formats
                template_reader.filters = filters
                template_reader.alts = alts
                template_reader.contigs = contigs

            # Add new source to metadata of header
            if not(new_source is None):
                sources = template_reader.metadata.setdefault("source", [])
                sources.append(new_source)

            # Add new info fields to header
            for info_id, info_len, info_type, info_desc, _, _ in new_info_fields:
                info_field = vcf.parser._Info(info_id, info_len, info_type, info_desc, None, None)
                template_reader.infos[info_id] = info_field

            # Add new format fields to header
            for format_id, format_len, format_type, format_desc in new_format_fields:
                format_field = vcf.parser._Format(format_id, format_len, format_type, format_desc)
                template_reader.formats[format_id] = format_field

            # Add new filters to header
            for filter_id, filter_desc in new_filters:
                filter_field = vcf.parser._Filter(filter_id, filter_desc)
                template_reader.filters[filter_id] = filter_field

            self.writer = vcf.Writer(out_file, template_reader, lineterminator='\n')
        else:
            raise NotSupportedException('File type unsupported: ' + file_type)


    def write_record(self, record):
        """ Method to write a raw record.  Useful for outputing a record directly (such as when
            combining multiple files)
        """
        if self.file_type == 'VCF4.1':
            self.writer.write_record(record)
        else:
            raise NotSupportedException('File type unsupported: ' + self.file_type)

def get_record_type(record):
    return record.var_type

def get_record_length(record):
    ''' Measure the length change of a variant '''
    ref = get_record_ref(record)
    alts = get_record_alt_alleles(record)
    alleles = [ref] + alts
    (genotype, phased) = get_record_genotype_phased(record)
    alts_in_use = [g for g in genotype if g > 0]
    if len(alts_in_use) == 0:
        return len(alts[0]) - len(ref)
    else:
        return len(alleles[alts_in_use[0]]) - len(ref)

def get_allele_length(ref, alt):
    return len(alt) - len(ref)


def is_record_bad_variant(record):
    if record.POS is None:
        return (True, "position is none, not supported")
    if record.REF is None:
        return (True, "reference allele is none, not supported")
    elif len(record.REF) > 0:
        for base in record.REF:
            if base not in {'A','C','G','T','N'}:
                return (True, "ref contains non ACGT characters, not supported")
    else:
        return (True, "ref is 0 length, not supported")
    if record.ALT is None:
        return (True, "alt allele is none, not supported")
    if record.ALT is not None and len(record.ALT) > 0:
        for alt in record.ALT:
            if alt is None:
                return (True, "alt allele is none, not supported")
            else:
                print "alt"
                print alt
                if len(str(alt)) > 0:
                    for base in str(alt):
                        if base not in {'A','C','G','T','N'}:
                            return (True, "alt contains not ACGT characters, not supported")
                else:
                    return (True, "alt is zero length, not supported")
    if record.CHROM is None:
        return (True, "chromosome is None, not supported")
    return (False, '')

def get_record_chrom(record):
    return record.CHROM

def get_record_pos(record):
    # note: this is 1-indexed
    return record.POS

def get_record_var_id(record):
    return record.ID

def get_record_filters(record):
    return record.FILTER

def check_vcf_record(record, fasta):
    ''' Returns a string containing alarms if there were any.
        Returns an empty string if there were no alarms.
        Raises ValueError on malformed input '''
    chrom = get_record_chrom(record)
    pos = get_record_pos(record) - 1
    ref = get_record_ref(record)
    alt_alleles = get_record_alt_alleles(record)

    alarms = []

    if not chrom in fasta:
        raise ValueError("VCF contains chromosome not in reference.")

    if ref != ref.upper():
        alarms.append("VCF contains reference alleles that are lower case or mixed case.")

    ref_seq = fasta[chrom][pos:(pos+len(ref))].upper()

    if ref_seq != ref:
        alarms.append("VCF contains reference alleles that do not exist at that location in reference. Could have been produced with different reference.")

    if alt_alleles is not None:
        if len(alt_alleles) == 0:
            raise ValueError("VCF contains entries with empty list of ALT alleles.")

        for allele in alt_alleles:
            if allele is None or allele == '.':
                raise ValueError("VCF contains enties with empty ALT alleles or one or more ALT alleles of '.'.")
            elif allele != allele.upper():
                alarms.append("VCF contains entries in which ALT alleles are lower case or mixed case.")
    else:
        raise ValueError("VCF contains entries with empty ALT alleles.")

    return ' '.join(alarms)

def get_record_passes_filters(record):
    ''' Record passes any VCF filtering '''
    filter_field = get_record_filters(record)
    return (filter_field is None or (filter_field is not None and (len(filter_field) == 1 and filter_field[0] == "PASS")) or filter_field == [])

def get_record_has_filter(filter, record):
    if type(record.FILTER) is list:
        for f in record.FILTER:
            if f == filter:
                return True
        return False
    else:
        return filter == record.FILTER

def get_record_ref(record):
    return record.REF

def get_record_alt_alleles(record):
    return [str(x) if not(x is None) else x for x in record.ALT]

def get_record_qual(record):
    return record.QUAL

def get_record_allele_freqs(record):
    info = record.INFO
    return info.get('AF')

def get_record_ref_allele_count(record):
    return record.INFO.get('RO')

def get_record_alt_allele_counts(record):
    return record.INFO.get('AO')

def get_record_post_homopolymer_counts(record):
    posthpc = record.INFO.get('POSTHPC')
    return posthpc if posthpc is not None else [0]

def get_record_mean_mapq_alt(record):
    mumapq = record.INFO.get('MUMAP_ALT')
    return mumapq if mumapq is not None else [-1]

def get_record_mean_mapq_ref(record):
    mumapq = record.INFO.get('MUMAP_REF')
    return mumapq if mumapq is not None else -1

def get_record_post_homopolymer_bases(record):
    posthpb = record.INFO.get('POSTHPB')
    return posthpb if posthpb is not None else ['N']

def get_record_sample_call(record, sample_num=0):
    samples = record.samples
    return samples[sample_num]

def get_record_sample_name(record, sample_num=0):
    sample_call = get_record_sample_call(record, sample_num=sample_num)
    return sample_call.sample

def get_record_genotype_phased(record, sample_num=0):
    sample_call = get_record_sample_call(record, sample_num=sample_num)
    data = sample_call.data
    gt = getattr(data, 'GT', None)
    if gt is None:
        raise Exception('GT field not present in VCF record')

    if '|' in gt:
        genotype = [int(x) if not x == '.' else x for x in gt.split('|')]
        phased = True
    else:
        genotype = [int(x) if not x == '.' else x for x in gt.split('/')]
        phased = False

    return genotype, phased

def get_record_homozygous(record):
    (geno, phased) = get_record_genotype_phased(record)
    if len(geno) < 2:
        return True
    else:
        return reduce(lambda x, y: x == y, geno)

def get_record_phase_set(record, sample_num=0):
    sample_call = get_record_sample_call(record, sample_num=sample_num)
    data = sample_call.data
    return getattr(data, 'PS', None)

def get_record_genotype_qual(record, sample_num=0):
    sample_call = get_record_sample_call(record, sample_num=sample_num)
    data = sample_call.data
    return getattr(data, 'GQ', None)

def set_record_genotype_likelihoods(record, genotype_likelihoods):
    set_data_field(record, 'GL', genotype_likelihoods)

def get_record_depth(record, sample_num=0):
    sample_call = get_record_sample_call(record, sample_num=sample_num)
    data = sample_call.data
    return getattr(data, 'DP', None)

def get_record_phase_qual(record, sample_num=0):
    sample_call = get_record_sample_call(record, sample_num=sample_num)
    data = sample_call.data
    return getattr(data, 'PQ', None)

def get_record_barcodes(record, sample_num=0):
    sample_call = get_record_sample_call(record, sample_num=sample_num)
    data = sample_call.data
    barcode_list = getattr(data, PROCESSED_BARCODE_TAG, None)

    if barcode_list is None:
        return []

    barcodes = [bcs.split(';') for bcs in barcode_list]
    for j in xrange(len(barcodes)):
        if barcodes[j] == ['']:
            barcodes[j] = []

    return barcodes

def get_sv_record_mate_chrom(record):
    assert(record.var_type == 'sv')
    alt = record.ALT[0]
    return alt.chr

def get_sv_record_mate_pos(record):
    assert(record.var_type == 'sv')
    alt = record.ALT[0]
    return alt.pos

def get_sv_record_mate_var_id(record):
    assert(record.var_type == 'sv')
    info = record.INFO
    var_id2 = info.get('MATEDID')
    return var_id2

def get_sv_record_directs(record):
    assert(record.var_type == 'sv')
    alt = record.ALT[0]
    direct1 = alt.orientation
    direct2 = alt.remoteOrientation
    return direct1, direct2

def create_translocation_records(chrom1, pos1, ref1, id1, direct1, chrom2, pos2, ref2, id2, direct2, qual, sample_name,
                                pos1_uncertain=0, pos2_uncertain=0, barcodes=None, phase_set1=None, phase1=None, phase_set2=None, phase2=None):
    alt1 = [vcf.model._Breakend(chrom2, pos2, direct1, direct2, ref1, True)]
    alt2 = [vcf.model._Breakend(chrom1, pos1, direct2, direct1, ref2, True)]

    info1 = OrderedDict()
    info1['SVTYPE'] = 'BND'
    info1['MATEID'] = id2
    info1['CIPOS'] = [0, pos1_uncertain]

    info2 = OrderedDict()
    info2['SVTYPE'] = 'BND'
    info2['MATEID'] = id1
    info2['CIPOS'] = [0, pos2_uncertain]

    format = ':'.join(['GT', 'PS', PROCESSED_BARCODE_TAG])
    if not(barcodes is None):
        barcodes_string = ';'.join(barcodes)
    else:
        barcodes_string = ''
    data_instantiator = vcf.model.make_calldata_tuple(['GT', 'PS', PROCESSED_BARCODE_TAG])
    if phase_set1 is None:
        phase_set1 = '.'
    if phase_set2 is None:
        phase_set2 = '.'

    if not(phase1 is None):
        if phase1 == 0:
            gt1 = '1|0'
        else:
            assert(phase1 == 1)
            gt1 = '0|1'
    else:
        gt1 = '0/1'

    if not(phase2 is None):
        if phase2 == 0:
            gt2 = '1|0'
        else:
            assert(phase2 == 1)
            gt2 = '0|1'
    else:
        gt2 = '0/1'

    data1 = data_instantiator(gt1, phase_set1, barcodes_string)
    data2 = data_instantiator(gt2, phase_set2, barcodes_string)

    sample_call1 = vcf.model._Call(None, sample_name, data1)
    samples1 = [sample_call1]

    sample_call2 = vcf.model._Call(None, sample_name, data2)
    samples2 = [sample_call2]

    sample_indexes = {sample_name: 0}
    record1 = vcf.model._Record(chrom1, pos1, id1, ref1, alt1, qual, [], info1, format, sample_indexes, samples=samples1)
    record2 = vcf.model._Record(chrom2, pos2, id2, ref2, alt2, qual, [], info2, format, sample_indexes, samples=samples2)

    return record1, record2

def set_record_filters(record, filters):
    record.FILTER = filters

# Emulate a named tuple, but don't require creation of a new
# type every time a field is added
class FakeNamedTuple(object):
    def __init__(self, fields, data):
        self._fields = fields
        self._data = data
        self._dict = { f:d for (f,d) in zip(fields, data)}

    def __getattr__(self, attr):
        return self._dict[attr]

    def __getitem__(self, key):
        return self._data[key]

    def _asdict(self):
        return self._dict

# Takes in a pyvcf Record object and returns the same object
# with the field for field_name replaced by the specified value
def set_data_field(record, field_name, field_val):
    assert(len(record.samples) == 1)
    new_format = record.FORMAT
    new_fields = new_format.split(':')
    if not(field_name in new_fields):
        new_fields = new_fields + [field_name]
        new_format = ':'.join(new_fields)

    sample_call = get_record_sample_call(record)
    data = sample_call.data
    data_dict = data._asdict()
    data_dict[field_name] = field_val

    new_sample_vals = []
    for field in new_fields:
        new_sample_vals.append(data_dict[field])

    # Note - the old way of passing the fields to pyVCF is memory intensive
    # because a fresh type is allocated for each call to make_calldata_tuple
    #data_instantiator = vcf.model.make_calldata_tuple(new_fields)
    #data = data_instantiator(*new_sample_vals)
    data = FakeNamedTuple(new_fields, new_sample_vals)
    sample_call.data = data
    record.samples[0] = sample_call
    record.FORMAT = new_format


def set_record_genotype_phased(record, genotype, phased):
    if phased:
        gt = str(genotype[0]) + '|' + str(genotype[1])
    else:
        gt = str(genotype[0]) + '/' + str(genotype[1])

    set_data_field(record, 'GT', gt)

def set_record_phase_set(record, phase_set):
    set_data_field(record, 'PS', phase_set)

def set_record_genotype_qual(record, genotype_qual):
    set_data_field(record, 'GQ', genotype_qual)

def set_record_phase_qual(record, phase_qual):
    set_data_field(record, 'PQ', phase_qual)

def set_record_junction_qual(record, junction_qual):
    set_data_field(record, 'JQ', junction_qual)

def set_record_barcodes(record, barcodes):
    barcode_list = [';'.join(bcs) for bcs in barcodes]

    set_data_field(record, PROCESSED_BARCODE_TAG, barcode_list)


def combine_vcfs(output_filename, input_vcf_filenames):
    tmp_filename = output_filename + ".tmp"

    for (i,fn) in enumerate(input_vcf_filenames):
        if i == 0:
            args = 'cat ' + fn
            subprocess.check_call(args + " > " + tmp_filename, shell=True)
        else:
            args = 'grep -a -v "^#" ' + fn
            ret = subprocess.call(args + " >> " + tmp_filename, shell=True)
            if ret == 2:
                raise Exception("grep call failed: " + args)

    # Sort and index the files
    tk_tabix.sort_vcf(tmp_filename, output_filename)
    tk_tabix.index_vcf(output_filename)

    try:
        os.remove(tmp_filename)
    except:
        pass


def get_locus_info(locus):
    """ Returns chrom, start and stop from locus string.
    Enforces standardization of how locus is represented.
    chrom:start_stop (start and stop should be ints or 'None')
    """
    chrom, start_stop = locus.split(':')
    if chrom == 'None':
        chrom = None

    start, stop = re.split("\.\.|-", start_stop)
    if start == 'None':
        start = None
    else:
        start = int(float(start))
    if stop == 'None':
        stop = None
    else:
        stop = int(float(stop))
    return (str(chrom), start, stop)

def create_locus_info(chrom, start, stop):
    """ Inverse to above
    """
    return str(chrom) + ':' + str(start) + '..' + str(stop)

def read_has_barcode(read):
    try:
        read.opt(RAW_BARCODE_TAG)
        return True
    except KeyError:
        return False

def get_read_barcode(read):
    '''Get the 10X barcode sequence for a read.  Returns None if no barcode attached
       or a non-whitelist sequence was observed '''
    try:
        r = read.opt(PROCESSED_BARCODE_TAG)
        if r == '':
            return None
        else:
            return r
    except KeyError:
        return None
def get_read_barcode_or_raw(read):
    try:
        r = read.opt(PROCESSED_BARCODE_TAG)
        if r == '':
            r = read.opt(RAW_BARCODE_TAG)
            return r
        else:
            return r
    except:
        r = read.opt(RAW_BARCODE_TAG)
        return r
def get_read_phase_set(read):
    try:
        r = read.opt(PHASE_SET_BAM_TAG)
        if r == '':
            return None
        return int(r)
    except KeyError:
        return None

def get_read_chimeric_alignments(read):
    for (label, value) in read.tags:
        if label == 'SA':
            if value == '':
                return None
            alignment_strs = [s.split(',') for s in value.split(';')]
            alignments = []
            for a in alignment_strs:
                if len(a) == 6:
                    # convert to 0-based
                    # return chrom, pos, strand, cigar, mapq, nm
                    alignments.append((a[0], int(a[1]) - 1, a[2], a[3], int(a[4]), int(a[5])))
            return alignments
    return None

def get_read_raw_barcode(read):
    try:
        return read.opt(RAW_BARCODE_TAG)
    except KeyError:
        return None

def get_read_barcode_qual(read):
    try:
        return read.opt(RAW_BARCODE_QUAL_TAG)
    except KeyError:
        return None

def get_read_sample_index(read):
    index = None
    for (label, value) in read.tags:
        if label == SAMPLE_INDEX_TAG:
            index = value
    return index

def get_read_sample_index_qual(read):
    index_qual = None
    for (label, value) in read.tags:
        if label == SAMPLE_INDEX_QUAL_TAG:
            index_qual = value
    return index_qual

def get_read_molecule_conf(read):
    try:
        r = read.opt('MC')
        if r == '':
            return 0.5
        return float(r)
    except KeyError:
        return 0.5

def get_read_haplotype(read):
    try:
        r = read.opt(HAPLOTYPE_BAM_TAG)
        if r == '':
            return None
        else:
            return int(r)
    except KeyError:
        return None

def get_target_regions_dict(targets_file, feature_name=None):
    """ Gets the target regions from a targets file as a chrom-indexed dictionary,
    with every entry given as a list of (start, end) tuples
    """
    targets = {}
    for line in targets_file:
        info = line.strip().split('\t')
        if line.startswith('browser') or line.startswith('track') or line.startswith('-browser') or line.startswith('-track') or line.startswith('#'):
            continue
        if len(line.strip()) == 0:
            continue
        if feature_name is not None:
            if len(info) < 4 or info[3] != feature_name:
                continue
        chrom = info[0]
        start = int(info[1])
        end = int(info[2])
        chrom_targs = targets.setdefault(chrom, [])
        chrom_targs.append((start, end))
    return targets

def get_target_regions(targets_file, feature_name=None):
    """ Gets the target regions from a targets file as a chrom-indexed dictionary,
    with every entry given as a list of Regions objects
    """
    targets_dict = get_target_regions_dict(targets_file, feature_name)
    target_regions = {}
    for (chrom, starts_ends) in targets_dict.iteritems():
        chrom_regions = Regions(regions=starts_ends)
        target_regions[chrom] = chrom_regions
    return target_regions

def get_bed_iterator(bed_file_name, locus=None):
    bed_file = open(bed_file_name)
    if not(locus is None):
        locus_chrom, locus_start, locus_stop = get_locus_info(locus)
    for line in bed_file:
        if line.startswith('browser') or line.startswith('track') or line.startswith('-browser') or line.startswith('-track') or line.startswith('#'):
            continue
        info = line.strip().split('\t')
        chrom = info[0]
        start = int(info[1])
        end = int(info[2])
        if locus is None:
            yield(chrom, start, end)
        else:
            if chrom == locus_chrom and start >= locus_start and end < locus_stop:
                yield (chrom, start, end)

def get_vc_mode(vc_precalled, vc_mode_in):

    args = vc_mode_in.split(":")
    vc_mode_a = args[0]

    # 'Normal' precalled mode -- supply a precalled and set vc_mode to null or 'disable'
    if vc_precalled is not None and (vc_mode_in is None or vc_mode_in == 'disable'):
        vc_mode = "precalled"
        precalled_file = vc_precalled
        variant_caller = None
        variant_caller_path = None

    # 'Legacy' precalled mode -- supply 'precalled:<fn>' to vc_mode. vc_precalled is ignored in this case
    elif vc_mode_a == "precalled":
        vc_mode = "precalled"
        precalled_file = args[1]
        variant_caller = None
        variant_caller_path = None

    # Modes where we call variants -- supply vc_mode as 'gatk' or 'freebayes'
    elif vc_mode_a in ["gatk", "freebayes"]:
        variant_caller = vc_mode_a

        variant_caller_path = None
        if variant_caller == "gatk":
            variant_caller_path = args[1]

        # If you also supply a precalled file, you're in precalled_plus mode,
        # otherwise you're in call mode
        if vc_precalled is None:
            vc_mode = "call"
            precalled_file = None
        else:
            vc_mode = "precalled_plus"
            precalled_file = vc_precalled

    else:
        raise Exception("Invalid variant calling config")

    print (vc_mode, variant_caller, precalled_file, variant_caller_path)
    return (vc_mode, variant_caller, precalled_file, variant_caller_path)
# get_vc_mode

################################################################################
class CopyNumberIO():
    #...........................................................................
    def __init__(self, profiles, path):
        self.profiles = profiles
        self.path = path
        if self.path[-1] != "/":
            self.path = self.path + "/"
        # if path
        self.tool_dir="/mnt/home/wei/bin/"
    # __init__

    #...........................................................................
    def write_reference(self):
        bin_size = self.profiles["BinSize"]
        reference = self.profiles["Reference"]
        reference_file = open(self.path + "Reference.txt", "w")
        reference_file.write("BinSize\t%d" % bin_size)
        for key in sorted(reference):
            reference_file.write("%s\t%f" % (key, reference[key]))
        # for key
        reference_file.close()
    # write_reference

    #...........................................................................
    def write_bedGraph(self):
        bin_size = self.profiles["BinSize"]
        self.write_reference()
        for profile_id in sorted(self.profiles["Profiles"]):
            bed_graph_file = open(self.path + profile_id + ".bedGraph", "w")
            data = self.profiles["Profiles"][profile_id]
            for chr_name in sorted(data):
                end = 0
                for value in data[chr_name]:
                    start = end + 1
                    end += bin_size
                    bed_graph_file.write("%s\t%d\t%d\t%f\n" % (chr_name, start, end, value))
                # for record
            # for chrName
            bed_graph_file.close()
        # for profile_id
    # write_bedGraph

    #...........................................................................
    def read_reference(self):
        profiles = {}
        reference = {}
        reference_file = open(self.path + "Reference.txt", "rU")
        bin_size = float("nan")
        for line in reference_file:
            line = line.strip("\n").strip("\r")
            (key, value) = line.split("\t")
            if key == "BinSize":
                bin_size = int(value)
            else:
                reference[key] = float(value)
            # if BinSize else
        # for line
        reference_file.close()
        profiles["BinSize"] = bin_size
        profiles["Reference"] = reference
        self.profiles = profiles
        return profiles
    # read_reference

    #...........................................................................
    def read_bedGraph(self):
        profiles = self.read_reference()
        #bed_graph_files = os.listdir(self.path)
        bed_graph_files = glob.glob(self.path + "*.bedGraph")
        profiles["Profiles"] = {}
        for file_name in bed_graph_files:
            profile_id = os.path.splitext(os.path.basename(file_name))[0]
            bed_graph_file = open(file_name, "r")
            profile = {}
            chr_profile = []
            previous_chr = ""
            for line in bed_graph_file:
                line = line.strip("\n").strip("\r")
                (chromosome, start, stop, count) = line.split("\t")
                if chromosome == previous_chr:
                    chr_profile.append(float(count))
                else:
                    if previous_chr != "":
                        profile[previous_chr] = chr_profile
                    # if previous_chr else
                    chr_profile = [float(count)]
                    previous_chr = chromosome
                # if chromosome else
            # for line
            profiles["Profiles"][profile_id] = profile
        # for bed_graph_file
        bed_graph_file.close();
        self.profiles = profiles
        return profiles
    # read_bedGraph

    #...........................................................................
    # See https://github.com/ENCODE-DCC/kentUtils/tree/master/src/utils/bedGraphToBigWig
    # Also see https://genome.ucsc.edu/goldenpath/help/bigWig.html
    # Also see https://genome.ucsc.edu/goldenpath/help/hubQuickStartGroups.html#multiWig
    def convert_bedGraph_to_bigWig(self, remove_bedGraph=False):
        for profile_id in sorted(self.profiles["Profiles"]):
            bedGraph_file = self.path + profile_id + ".bedGraph"
            bigWig_file = self.path + profile_id + ".bw"
            # TODO: general case - any set of chromosome names/lengths.
            # (use self.profiles["Reference"] to generate chrom.sizes)
            cmd = self.tool_dir + "bedGraphToBigWig " + bedGraph_file + " " + self.tool_dir + "hg19.chrom.sizes " + bigWig_file
            os.system(cmd)
            if remove_bedGraph:
                os.remove(outs_hp_bedGraph[hp])
            # if remove
        # for profile_id
    # convert_bedGraph_to_bigWig
# class CopyNumberIO
