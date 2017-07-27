#!/usr/bin/env python
#
# Copyright (c) 2017 10x Genomics, Inc. All rights reserved.
#
import martian
import os

import cellranger.constants as cr_constants
import cellranger.molecule_counter as cr_mol_counter
import cellranger.utils as cr_utils

__MRO__ = """
stage AGGREGATOR_PREFLIGHT(
    in map[]  sample_defs,
    in string normalization_mode,
    src py     "stages/aggregator/aggregator_preflight",
)
"""
def main(args, outs):
    if args.normalization_mode is not None and args.normalization_mode not in cr_constants.NORM_MODES:
        martian.exit("Normalization mode must be one of: %s" % ', '.join(cr_constants.NORM_MODES))

    global_mkref_version = None
    global_fasta_hash = None
    global_gtf_hash = None
    global_chemistry = None

    # TODO make assertions about the required metrics!

    # check that all molecule files conform to spec
    files_seen = set()
    libraries_seen = set()
    for sample in args.sample_defs:
        library_id = sample[cr_constants.AGG_ID_FIELD]
        if len(library_id) == 0:
            martian.exit("Library ID cannot be empty: %s" % sample)

        if library_id in libraries_seen:
            martian.exit("Same library ID is specified on multiple rows: %s" % library_id)
        else:
            libraries_seen.add(library_id)

        mol_h5 = sample[cr_constants.AGG_H5_FIELD]

        if not os.path.exists(mol_h5):
            martian.exit("Input molecule file does not exist: %s" % mol_h5)

        if not os.access(mol_h5, os.R_OK):
            martian.exit("Input molecule file is not readable, please check file permissions: %s" % mol_h5)

        h5_filetype = cr_utils.get_h5_filetype(mol_h5)
        if h5_filetype and h5_filetype != cr_mol_counter.MOLECULE_H5_FILETYPE:
            martian.exit("Input is a %s file, but a molecule file is required" % h5_filetype)

        if mol_h5 in files_seen:
            martian.exit("Same molecule file is specified in multiple sample definitions: %s" % mol_h5)
        else:
            files_seen.add(mol_h5)

        with cr_mol_counter.MoleculeCounter.open(mol_h5, 'r') as counter:

            mol_file_version = counter.file_version
            if mol_file_version < 2:
                martian.exit("Molecule file is out of date (format version must be >= 2) %s" % mol_h5)

            mol_cr_version = counter.get_metric('cellranger_version')
            if not mol_cr_version:
                martian.exit("Molecule file was produced with old cellranger version (missing version number): %s" % mol_h5)

            mol_mkref_version = counter.get_metric('reference_mkref_version')
            if global_mkref_version is None:
                global_mkref_version = mol_mkref_version
            elif global_mkref_version != mol_mkref_version:
                martian.exit("Molecules were produced using different mkref versions (%s, %s)" % (global_mkref_version, mol_mkref_version))

            mol_fasta_hash = counter.get_metric('reference_fasta_hash')
            if global_fasta_hash is None:
                global_fasta_hash = mol_fasta_hash
            elif global_fasta_hash != mol_fasta_hash:
                martian.exit("Molecules were produced using different mkref versions (%s, %s)" % (global_fasta_hash, mol_fasta_hash))

            mol_gtf_hash = counter.get_metric('reference_gtf_hash')
            if global_gtf_hash is None:
                global_gtf_hash = mol_gtf_hash
            elif global_gtf_hash != mol_gtf_hash:
                martian.exit("Molecules were produced using different mkref versions (%s, %s)" % (global_gtf_hash, mol_gtf_hash))

            mol_chemistry = counter.get_metric('chemistry_name')
            if global_chemistry is None:
                global_chemistry = mol_chemistry
            elif global_chemistry != mol_chemistry:
                martian.exit("Molecules were produced with different chemistries (%s, %s)" % (global_chemistry, mol_chemistry))

            if counter.is_aggregated():
                martian.exit("Molecule file was aggregated from multiple samples: %s.\n Aggregated outputs cannot be re-aggregated, please pass each of the original samples instead." % mol_h5)

            if counter.nrows() == 0:
                martian.exit("Cannot aggregate file because it contains no data: %s.\n Please remove this file from the aggregation and try again." % mol_h5)

            for (gg, metrics) in counter.get_metric(cr_mol_counter.GEM_GROUPS_METRIC).iteritems():
                gg_total_reads = metrics[cr_mol_counter.GG_TOTAL_READS_METRIC]
                if gg_total_reads == 0:
                    martian.exit("Gem group %d has zero reads in file: %s\n Please re-run `cellranger count` without including this gem group." % (gg, mol_h5))
