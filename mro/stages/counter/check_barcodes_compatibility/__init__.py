#!/usr/bin/env python
#
# Copyright (c) 2018 10X Genomics, Inc. All rights reserved.
#
import itertools
import martian
import math
from collections import defaultdict, Counter
import cellranger.chemistry as cr_chem
import cellranger.constants as cr_constants
from cellranger.library_constants import GENE_EXPRESSION_LIBRARY_TYPE
import cellranger.utils as cr_utils
from cellranger.fastq import FastqReader
import tenkit.stats as tk_stats


__MRO__ = """
stage CHECK_BARCODES_COMPATIBILITY(
    in  map[]    chunks,
    in  string   barcode_whitelist,
    in  int      num_reads_to_check_barcode,
    in  float    barcode_compatibility_cutoff,
    out bool     barcode_compatible,
    out map      barcode_compatibility_info,
    src py       "stages/counter/check_barcodes_compatibility",
) split (
    in  map      read_chunks,
    in  map      chemistry,
    in  bool     reads_interleaved,
    out string[] sampled_barcodes,
)
"""

# throwing away the first N reads at the beginning of a fastq file
READS_BURNIN = 10000

def cosine_similarity(c1, c2):
    """ Calculate the cosine similarity of two barcode lists.

    The barcode list is represented as 'bag of barcode', which is a
    collections.Counter object (barcodes as keys and their counts as
    values. This function will calculate the cosine similarity of two
    barcode lists.

    Args:
        c1 (Counter): A Counter object of the first barcode list.
        c2 (Counter): A Counter object of the second barcode list.

    Returns:
        float: The cosine similarity of two barcode lists.
    """
    terms = set(c1).union(c2)
    dotprod = sum(c1.get(k, 0) * c2.get(k, 0) for k in terms)
    magA = math.sqrt(sum(c1.get(k, 0)**2 for k in terms))
    magB = math.sqrt(sum(c2.get(k, 0)**2 for k in terms))
    return tk_stats.robust_divide(dotprod, magA * magB)


def robust_cosine_similarity(c1, c2, robust_fraction_threshold=0.925):
    """ Calculate the robust cosine similarity of two barcode lists.
        Implements the same basic cosine distance metric, but first
        caps each count value to a threshold value such that robust_fracion_threshold
        of the total counts are greater than or equal to the threshold.
        This reduces the impact of very high count outliers.

    Args:
        c1 (Counter): A Counter object of the first barcode list.
        c2 (Counter): A Counter object of the second barcode list.

    Returns:
        float: The robust cosine similarity of two barcode lists.
    """

    thresh1 = max(1, tk_stats.NX(c1.values(), robust_fraction_threshold))
    thresh2 = max(1, tk_stats.NX(c2.values(), robust_fraction_threshold))

    terms = set(c1).union(c2)
    dotprod = sum(min(thresh1, c1.get(k, 0)) * min(thresh2, c2.get(k, 0)) for k in terms)
    magA = math.sqrt(sum(min(thresh1, c1.get(k, 0))**2 for k in terms))
    magB = math.sqrt(sum(min(thresh2, c2.get(k, 0))**2 for k in terms))
    return tk_stats.robust_divide(dotprod, magA * magB)


def split(args):
    # determine number of fastq file for each library and gem group, {gem_group : {library_type : count_of_fastq_file} }
    chunk_counts = defaultdict(lambda: defaultdict(int))
    for chunk in args.chunks:
        chunk_counts[chunk["gem_group"]][chunk["library_type"]] += 1

    single_library = True
    for gem_group in chunk_counts:
        if len(chunk_counts[gem_group]) > 1:
            single_library = False

    if single_library:
        martian.log_info('Single library in input. No need to check barcode compatibility.')
        # `[]` for the chunks will skip the main
        return {'chunks': [], 'join': {}}

    num_reads_to_check_barcode = cr_constants.NUM_READS_TO_CHECK_BARCODE if args.num_reads_to_check_barcode is None else args.num_reads_to_check_barcode
    chunks = []
    for chunk in args.chunks:
        chunk_def = chunk
        chunk_def['num_reads_per_chunk_to_check_barcode'] = int(tk_stats.robust_divide(num_reads_to_check_barcode, chunk_counts[chunk["gem_group"]][chunk["library_type"]]))
        chunks.append(chunk_def)

    return {'chunks': chunks, 'join': {'__mem_gb': 4}}

def join(args, outs, chunk_defs, chunk_outs):
    outs.barcode_compatible = True

    if chunk_outs is None or len(chunk_outs) == 0:
        return

    # aggreagate barcodes from chunk, {gem_group : {library_type : count_of_fastq_file} }
    sampled_barcodes = defaultdict(lambda: defaultdict(list))
    for chunk_def, chunk_out in zip(chunk_defs, chunk_outs):
        gem_group, lib = chunk_def.gem_group, chunk_def.library_type
        sampled_barcodes[gem_group][lib].extend(chunk_out.sampled_barcodes)

    barcodes_in_whitelist = cr_utils.load_barcode_whitelist(args.barcode_whitelist, as_set=True)
    barcode_translate_map = cr_utils.load_barcode_translate_map(args.barcode_whitelist)

    sampled_bc_counter_in_wl = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    outs.barcode_compatibility_info = {}  # record sampled barcode info

    for gem_group in sampled_barcodes:
        outs.barcode_compatibility_info[gem_group] = {}
        for lib in sampled_barcodes[gem_group]:
            sampled_bc = sampled_barcodes[gem_group][lib]
            unique_bc = set(sampled_bc)
            unique_bc_in_wl = unique_bc.intersection(barcodes_in_whitelist)

            outs.barcode_compatibility_info[gem_group][lib] = {}
            outs.barcode_compatibility_info[gem_group][lib]['num_barcodes_sampled'] = len(sampled_bc)
            outs.barcode_compatibility_info[gem_group][lib]['num_barcodes_sampled_unique'] = len(unique_bc)
            outs.barcode_compatibility_info[gem_group][lib]['num_barcodes_sampled_unique_in_whitelist'] = len(unique_bc_in_wl)

            raw_bc_counter = Counter(sampled_bc)
            bc_counter = defaultdict(int)
            for b in unique_bc_in_wl:
                if (lib != GENE_EXPRESSION_LIBRARY_TYPE) and (barcode_translate_map is not None):
                    bc_trans = barcode_translate_map.get(b, b)
                else:
                    bc_trans = b
                bc_counter[bc_trans] += raw_bc_counter[b]
            sampled_bc_counter_in_wl[gem_group][lib] = bc_counter

    barcode_compatibility_cutoff = cr_constants.BARCODE_COMPATIBILITY_CUTOFF if args.barcode_compatibility_cutoff is None else args.barcode_compatibility_cutoff

    outs.barcode_compatibility_info['pairwise_compatibility'] = {}
    # for every gem_group, check each pair of library types
    for gem_group in sampled_barcodes:
        outs.barcode_compatibility_info['pairwise_compatibility'][gem_group] = []
        for lib1, lib2 in itertools.combinations(sampled_barcodes[gem_group].keys(), 2):
            counter1, counter2 = sampled_bc_counter_in_wl[gem_group][lib1], sampled_bc_counter_in_wl[gem_group][lib2]
            overlap_size = len(set(counter1).intersection(set(counter2)))
            cosine_sim = robust_cosine_similarity(counter1, counter2)

            outs.barcode_compatibility_info['pairwise_compatibility'][gem_group].append([lib1, lib2, overlap_size, cosine_sim])
            if cosine_sim < barcode_compatibility_cutoff:
                outs.barcode_compatible = False

    # format warning/error message if incompatible
    if outs.barcode_compatible is False:
        log_msg = 'Barcodes from libraries are not compatible.'
        for gem_group in outs.barcode_compatibility_info['pairwise_compatibility']:
            for pair in outs.barcode_compatibility_info['pairwise_compatibility'][gem_group]:
                if pair[3] < barcode_compatibility_cutoff:
                    log_msg += '\n - GEM group {}: Barcodes from [{}] and [{}] have cosine similarity {:.4f}'.format(gem_group, pair[0], pair[1], pair[3])

        martian.log_info(log_msg)
        martian.exit(log_msg)
    return

def main(args, outs):
    bc_read_def = cr_chem.get_barcode_read_def(args.chemistry)
    bc_reader = FastqReader(args.read_chunks, bc_read_def, args.reads_interleaved, None, None)

    sampled_barcodes = []
    for name, seq, qual in itertools.islice(bc_reader.in_iter, READS_BURNIN, READS_BURNIN + args.num_reads_per_chunk_to_check_barcode):
        sampled_barcodes.append(seq)

    outs.sampled_barcodes = sampled_barcodes
