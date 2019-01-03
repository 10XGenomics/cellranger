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
    out map      skip_translate,
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
    outs.barcode_compatibility_info = {}  # record sampled barcode info
    outs.skip_translate = {}

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

            sampled_bc_counter_in_wl[gem_group][lib] = Counter(sampled_bc)

    barcode_compatibility_cutoff = cr_constants.BARCODE_COMPATIBILITY_CUTOFF if args.barcode_compatibility_cutoff is None else args.barcode_compatibility_cutoff

    pairwise_compatibility = {}
    exit_log_msg = "Barcodes from libraries are not compatible."

    for gem_group in sampled_barcodes:
        outs.skip_translate[gem_group] = {}
        pairwise_compatibility[gem_group] = {}
        library_types = sampled_barcodes[gem_group].keys()

        if len(library_types) < 2:
            continue

        if GENE_EXPRESSION_LIBRARY_TYPE in library_types:
            base_lib = GENE_EXPRESSION_LIBRARY_TYPE
            library_types.remove(base_lib)
            outs.skip_translate[gem_group][base_lib] = True
        else:
            # TODO: as for CR3.0, we need GEX for cell calling et al
            # at some point, we might support samples without GEX
            martian.exit("Gene expression data not found in the GEM group {}.".format(gem_group))

        base_lib_counter = sampled_bc_counter_in_wl[gem_group][base_lib]
        for lib in library_types:
            pair_key = '{}/{}'.format(base_lib, lib)
            pairwise_compatibility[gem_group][pair_key] = {}
            lib_counter = sampled_bc_counter_in_wl[gem_group][lib]

            # without translate
            overlap_size = len(set(base_lib_counter).intersection(set(lib_counter)))
            cosine_sim = robust_cosine_similarity(base_lib_counter, lib_counter)
            outs.skip_translate[gem_group][lib] = True

            # with translate
            if (lib != GENE_EXPRESSION_LIBRARY_TYPE) and (barcode_translate_map is not None):
                translated_counter = {barcode_translate_map.get(k,k) : v for (k,v) in lib_counter.iteritems()}
                overlap_size_translated = len(set(base_lib_counter).intersection(set(translated_counter)))
                cosine_sim_translated = robust_cosine_similarity(base_lib_counter, translated_counter)

                if cosine_sim_translated > cosine_sim:
                    outs.skip_translate[gem_group][lib] = False
                    overlap_size = overlap_size_translated
                    cosine_sim = cosine_sim_translated

            pairwise_compatibility[gem_group][pair_key]['overlap_size'] = overlap_size
            pairwise_compatibility[gem_group][pair_key]['cosine_similarity'] = cosine_sim
            if cosine_sim < barcode_compatibility_cutoff:
                outs.barcode_compatible = False
                exit_log_msg += '\n - GEM group {}: Barcodes from [{}] and [{}] have cosine similarity {:.4f}'.format(gem_group, base_lib, lib, cosine_sim)

    outs.barcode_compatibility_info['pairwise_compatibility'] = pairwise_compatibility
    # format warning/error message if incompatible
    if outs.barcode_compatible is False:
        martian.log_info(exit_log_msg)
        martian.exit(exit_log_msg)

    return

def main(args, outs):
    bc_read_def = cr_chem.get_barcode_read_def(args.chemistry)
    bc_reader = FastqReader(args.read_chunks, bc_read_def, args.reads_interleaved, None, None)

    sampled_barcodes = []
    for name, seq, qual in itertools.islice(bc_reader.in_iter, READS_BURNIN, READS_BURNIN + args.num_reads_per_chunk_to_check_barcode):
        sampled_barcodes.append(seq)

    outs.sampled_barcodes = sampled_barcodes
