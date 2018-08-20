#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
import cPickle
from collections import defaultdict
import json
import numpy as np
import cellranger.constants as cr_constants
import cellranger.library_constants as lib_constants
from cellranger.molecule_counter import MoleculeCounter
import cellranger.rna.library as rna_library
import cellranger.utils as cr_utils
import tenkit.safe_json as tk_safe_json
import tenkit.stats as tk_stats

__MRO__ = """
stage SUBSAMPLE_READS(
    in  h5     molecule_info,
    in  csv    filtered_barcodes,
    out json   summary,
    src py     "stages/counter/subsample_reads",
) split using (
    in  int    chunk_start,
    in  int    chunk_len,
    in  map[]  subsample_info,
    out pickle metrics,
)
"""

def get_cell_associated_barcodes(genomes, filtered_barcodes_csv):
    """ Get cell-associated barcodes by genome.
    Args:
      genomes (list of str): Genome names.
      filtered_barcodes_csv (str): Path to CSV file.
    Returns:
      dict of (str, set): Map genome to list of cell-assoc barcodes. Empty-string key is for all genomes."""
    cell_bcs = {}
    for genome in genomes:
        # Get all cell-assoc barcodes (ignoring genome) for the "" (blank) genome string
        cell_bcs[genome] = cr_utils.get_cell_associated_barcode_set(filtered_barcodes_csv,
                                                                    genome)
    # All cell-associated barcodes
    cell_bcs[''] = reduce(lambda x,y: x | y, cell_bcs.itervalues(), set())
    return cell_bcs

def split(args):
    # Get required info from the mol info
    mc = MoleculeCounter.open(args.molecule_info, 'r')

    genomes = sorted(set(f.tags.get('genome', '') for f in mc.feature_reference.feature_defs))
    cell_bcs_by_genome = get_cell_associated_barcodes(genomes, args.filtered_barcodes)

    # Get cell counts per gem group
    n_cells_per_gg = defaultdict(int)
    for bc in cell_bcs_by_genome['']:
        _, gem_group = cr_utils.split_barcode_seq(bc)
        n_cells_per_gg[gem_group] += 1

    # Assign gem group cell counts to their constituent libraries
    # TODO FIXME: Need to allow for per-library cell counts
    #   because some feature types might only have a subset of the GEX cell-assoc barcodes.
    n_cells_per_lib = np.zeros(len(mc.library_info), dtype=int)
    for lib_idx, lib in enumerate(mc.library_info):
        n_cells_per_lib[lib_idx] = n_cells_per_gg[lib['gem_group']]

    if n_cells_per_lib.sum() == 0:
        return {'chunks': []}

    library_info = mc.library_info

    raw_count_per_lib = np.array(mc.get_raw_read_pairs_per_library())
    raw_rppc_per_lib = raw_count_per_lib.astype(float) / n_cells_per_lib
    usable_count_per_lib = np.array(mc.get_usable_read_pairs_per_library())

    subsamplings = list() # track subsample info definitions

    library_types = sorted(set(lib['library_type'] for lib in library_info))
    for library_type in library_types:
        # All libraries w/ this type
        lib_indexes = np.array([i for i,lib in enumerate(library_info) if lib['library_type'] == library_type])

        # For plotting, we want a series of target depths that exist for all
        #   libraries w/ the same library type. When there's a single library
        #   per type (the common case), this is trivial - split it into deciles.
        #   But if there are multiple libraries with different depths, (e.g.,
        #   because gem-group-aggregation was used to increase cell numbers),
        #   we need to find depths that are achievable for all libraries.
        #   For now, let the lowest-depth library for a given type dictate this.
        min_raw_rppc = np.min(raw_rppc_per_lib[lib_indexes])

        # Use deciles of the raw read pairs per cell.
        deciles = np.arange(0.1, 1.1, 0.1)
        plot_targets = map(round, min_raw_rppc * deciles)

        # TODO: separate this work (internal + non)
        raw_targets = cr_constants.SUBSAMPLE_READS_PER_CELL + \
                      plot_targets

        # TODO: separate this work (internal + non)
        usable_targets = cr_constants.SUBSAMPLE_READS_PER_CELL + \
                         plot_targets

        for targets, depth_type in \
            ((raw_targets, cr_constants.RAW_SUBSAMPLE_TYPE), \
             ((usable_targets, cr_constants.MAPPED_SUBSAMPLE_TYPE)),):
            targets = sorted(list(set(map(int, targets))))
            for target_rppc in targets:
                if depth_type == cr_constants.RAW_SUBSAMPLE_TYPE:
                    # Infer the usable depth required to achieve this raw depth
                    usable_read_fracs = usable_count_per_lib.astype(float) / raw_count_per_lib
                    target_usable_counts = target_rppc * n_cells_per_lib * usable_read_fracs
                else:
                    target_usable_counts = target_rppc * n_cells_per_lib

                # Zero out libraries of the other types
                rates = np.zeros(len(library_info), dtype=float)
                rates[lib_indexes] = target_usable_counts[lib_indexes].astype(float) \
                                     / usable_count_per_lib[lib_indexes]

                # Clamp rates that are close to 1 to 1
                rates[np.absolute(rates - 1) < 1e-3] = 1

                # Zero out the libraries for which we have fewer reads than the target
                rates[rates > 1] = 0.0

                enough_data = np.any((rates > 0) & (rates <= 1))
                if not enough_data:
                    rates = np.zeros(len(rates))

                subsamplings.append({
                    'library_type': library_type,
                    'subsample_type': depth_type,
                    'target_read_pairs_per_cell': int(target_rppc),
                    'library_subsample_rates': list(map(float, rates)),
                })

    # Each chunk needs to store a piece of the mol info h5
    tgt_chunk_len = cr_constants.NUM_MOLECULE_INFO_ENTRIES_PER_CHUNK

    # Split the molecule info h5 into equi-RAM chunks
    chunks = []
    for chunk_start, chunk_len in mc.get_chunks(tgt_chunk_len, preserve_boundaries=True):
        chunks.append({
            'chunk_start': chunk_start,
            'chunk_len': chunk_len,
            'subsample_info': subsamplings,
            '__mem_gb': MoleculeCounter.estimate_mem_gb(chunk_len),
        })

    join = {
        '__mem_gb': 6,
    }

    mc.close()

    # TODO: is this really necessary w/ martian 3
    if len(chunks) == 0:
        chunks.append({
            'chunk_start': str(0),
            'chunk_len': str(0),
            'subsample_info': [],
        })

    return {'chunks': chunks, 'join': join}

def main(args, outs):
    np.random.seed(0)

    mc = MoleculeCounter.open(args.molecule_info, 'r')

    # Get cell-associated barcodes
    genomes = sorted(set(f.tags.get('genome', '') for f in mc.feature_reference.feature_defs))
    cell_bcs_by_genome = get_cell_associated_barcodes(genomes, args.filtered_barcodes)

    # Load chunk of relevant data from the mol_info
    chunk = slice(int(args.chunk_start), int(args.chunk_start) + int(args.chunk_len))
    mol_library_idx = mc.get_column_lazy('library_idx')[chunk]
    mol_read_pairs = mc.get_column_lazy('count')[chunk]
    mol_gem_group = mc.get_column_lazy('gem_group')[chunk]
    mol_barcode_idx = mc.get_column_lazy('barcode_idx')[chunk]
    mol_feature_idx = mc.get_column_lazy('feature_idx')[chunk]

    barcodes = mc.get_ref_column('barcodes')


    # Give each cell-associated barcode an integer index
    cell_bcs = sorted(list(cell_bcs_by_genome['']))
    cell_bc_to_int = {bc: i for i, bc in enumerate(cell_bcs)}


    # Give each genome an integer index
    genome_to_int = {g: i for i, g in enumerate(genomes)}
    feature_int_to_genome_int = np.array([genome_to_int[f.tags.get('genome', '')] for f in mc.feature_reference.feature_defs], dtype=int)
    mol_genome_idx = feature_int_to_genome_int[mol_feature_idx]


    # Run each subsampling task on this chunk of data
    n_tasks = len(args.subsample_info)
    n_genomes = len(genomes)
    n_cells = len(cell_bcs)

    umis_per_bc = np.zeros((n_tasks, n_genomes, n_cells))
    features_det_per_bc = np.zeros((n_tasks, n_genomes, n_cells))
    read_pairs_per_task = np.zeros((n_tasks, n_genomes))
    umis_per_task = np.zeros((n_tasks, n_genomes))

    for task_idx, task in enumerate(args.subsample_info):
        # Per-library subsampling rates
        rates_per_library = np.array(task['library_subsample_rates'], dtype=float)

        if np.count_nonzero(rates_per_library) == 0:
            continue

        mol_rate = rates_per_library[mol_library_idx]

        # Subsampled read pairs per molecule
        new_read_pairs = np.random.binomial(mol_read_pairs, mol_rate)

        # Compute tallies for each barcode
        group_keys = (mol_gem_group, mol_barcode_idx)
        group_values = (mol_feature_idx, mol_genome_idx, new_read_pairs)
        for (gg, bc_idx), (feature_idx, genome_idx, read_pairs) in \
            cr_utils.numpy_groupby(group_values, group_keys):

            barcode = cr_utils.format_barcode_seq(barcodes[bc_idx], gg)

            cell_idx = cell_bc_to_int.get(barcode)

            for this_genome_idx in xrange(len(genomes)):
                umis = np.flatnonzero((read_pairs > 0) & (genome_idx == this_genome_idx))

                # Tally UMIs and median features detected
                if barcode in cell_bcs_by_genome[genomes[this_genome_idx]]:
                    # This is a cell-associated barcode for this genome
                    umis_per_bc[task_idx, this_genome_idx, cell_idx] = len(umis)
                    features_det_per_bc[task_idx, this_genome_idx, cell_idx] = np.count_nonzero(np.bincount(feature_idx[umis]))

                # Tally numbers for duplicate fraction
                read_pairs_per_task[task_idx, this_genome_idx] += np.sum(read_pairs)
                umis_per_task[task_idx, this_genome_idx] += len(umis)

    with open(outs.metrics, 'w') as f:
        data = {
            'umis_per_bc': umis_per_bc,
            'features_det_per_bc': features_det_per_bc,
            'read_pairs': read_pairs_per_task,
            'umis': umis_per_task,
        }
        cPickle.dump(data, f, protocol = cPickle.HIGHEST_PROTOCOL)


def make_metric_name(name, library_type, genome, ss_type, ss_depth):
    lt_prefix = rna_library.get_library_type_metric_prefix(library_type)
    return '%s%s_%s_%s_%s' % (lt_prefix, genome, ss_type, ss_depth, name)

def compute_dup_frac(read_pairs, umis):
    return tk_stats.robust_divide(read_pairs - umis, read_pairs) if read_pairs > 0 else 0.0

def join(args, outs, chunk_defs, chunk_outs):
    # Merge tallies
    data = None
    for chunk in chunk_outs:
        with open(chunk.metrics) as f:
            chunk_data = cPickle.load(f)
        if data is None:
            data = chunk_data
        else:
            for k,v in data.iteritems():
                data[k] += chunk_data[k]

    # Compute metrics for each subsampling rate
    summary = {}

    with MoleculeCounter.open(args.molecule_info, 'r') as mc:
        genomes = sorted(set(f.tags.get('genome', '') for f in mc.feature_reference.feature_defs))
    cell_bcs_by_genome = get_cell_associated_barcodes(genomes, args.filtered_barcodes)

    # Give each cell-associated barcode an integer index
    cell_bcs = sorted(list(cell_bcs_by_genome['']))
    cell_bc_to_int = {bc: i for i, bc in enumerate(cell_bcs)}

    subsample_info = chunk_defs[0].subsample_info if len(chunk_defs) > 0 else []

    for i, task in enumerate(subsample_info):
        lib_type = task['library_type']
        ss_type = task['subsample_type']
        ss_depth = task['target_read_pairs_per_cell']

        if rna_library.has_genomes(lib_type):
            genome_ints = list(range(data['umis_per_bc'].shape[1]))
        else:
            genome_ints = [0]

        # Per-genome metrics
        for g in genome_ints:
            genome = genomes[g]

            # Only compute on cell-associated barcodes for this genome.
            # This only matters when there are multiple genomes present.
            cell_inds = np.array(sorted(cell_bc_to_int[bc] for bc in cell_bcs_by_genome[genome]))

            median_umis_per_cell = np.median(data['umis_per_bc'][i,g,cell_inds])
            summary[make_metric_name('subsampled_filtered_bcs_median_counts',
                                     lib_type, genome, ss_type, ss_depth)] = median_umis_per_cell

            median_features_per_cell = np.median(data['features_det_per_bc'][i,g,cell_inds])
            summary[make_metric_name('subsampled_filtered_bcs_median_unique_genes_detected',
                                     lib_type, genome, ss_type, ss_depth)] = median_features_per_cell

            dup_frac = compute_dup_frac(data['read_pairs'][i,g],  data['umis'][i,g])
            summary[make_metric_name('subsampled_duplication_frac',
                                     lib_type, genome, ss_type, ss_depth)] = dup_frac

        # Whole-dataset duplication frac
        all_read_pairs = np.sum(data['read_pairs'][i,:])
        all_umis = np.sum(data['umis'][i,:])
        dup_frac = compute_dup_frac(all_read_pairs, all_umis)
        summary[make_metric_name('subsampled_duplication_frac',
                                 lib_type, lib_constants.MULTI_REFS_PREFIX, ss_type, ss_depth)] = dup_frac


    with open(outs.summary, 'w') as f:
        json.dump(tk_safe_json.json_sanitize(summary), f, indent=4, sort_keys=True)
