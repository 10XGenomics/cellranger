#!/usr/bin/env python
#
# Copyright (c) 2018 10X Genomics, Inc. All rights reserved.
#
import cellranger.utils as cr_utils

def get_genomes_from_feature_ref(fref):
    return sorted(set(f.tags.get('genome', '') for f in fref.feature_defs))

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
