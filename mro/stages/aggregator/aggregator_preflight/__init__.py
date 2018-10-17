#!/usr/bin/env python
#
# Copyright (c) 2017 10x Genomics, Inc. All rights reserved.
#
import martian
import cellranger.constants as cr_constants
import cellranger.molecule_counter as cr_mol_counter
from cellranger.utils import string_is_ascii

__MRO__ = """
stage AGGREGATOR_PREFLIGHT(
    in map[]  sample_defs,
    in string normalization_mode,
    src py     "stages/aggregator/aggregator_preflight",
)
"""

def incompat_msg(reason):
    return ("The datasets you are trying to aggregate were created with different {reason}s, "
            "but 'cellranger aggr' requires identical {reason}s in order to combine datasets. "
            "Please re-run 'cellranger count' with uniform {reason} in order to aggregate "
            "these data.").format(reason=reason)

def main(args, outs):
    if args.normalization_mode is not None and args.normalization_mode not in cr_constants.NORM_MODES:
        martian.exit("Normalization mode must be one of: %s" % ', '.join(cr_constants.NORM_MODES))

    global_fasta_hash = None
    global_gtf_hash = None
    global_feature_ref = None
    chemistries = set()

    # TODO make assertions about the required metrics!

    # check that all molecule files conform to spec
    libraries_seen = set()
    for sample in args.sample_defs:
        library_id = sample[cr_constants.AGG_ID_FIELD]
        if len(library_id) == 0:
            martian.exit("Library ID cannot be empty: %s" % sample)

        if not string_is_ascii(library_id):
            martian.exit("Library ID %s contains unicode characters, only ASCII is allowed." % library_id)

        if cr_constants.AGG_BATCH_FIELD in sample:
            batch_name = sample[cr_constants.AGG_BATCH_FIELD]
            if not string_is_ascii(batch_name):
                martian.exit("Batch ID %s contains unicode characters, only ASCII is allowed." % batch_name)

        if library_id in libraries_seen:
            martian.exit("Same library ID is specified on multiple rows: %s" % library_id)
        else:
            libraries_seen.add(library_id)

        mol_h5 = sample[cr_constants.AGG_H5_FIELD]

        with cr_mol_counter.MoleculeCounter.open(mol_h5, 'r') as counter:

            mol_cr_version = counter.get_metric('cellranger_version')
            if not mol_cr_version:
                martian.exit("Molecule file was produced with old cellranger version (missing version number): %s" % mol_h5)

            mol_fasta_hash = counter.get_metric('reference_fasta_hash')
            if global_fasta_hash is None:
                global_fasta_hash = mol_fasta_hash
            elif global_fasta_hash != mol_fasta_hash:
                martian.exit("{} (hashes: {} != {})".format(incompat_msg("genome reference"), global_fasta_hash, mol_fasta_hash))

            mol_gtf_hash = counter.get_metric('reference_gtf_hash')
            if global_gtf_hash is None:
                global_gtf_hash = mol_gtf_hash
            elif global_gtf_hash != mol_gtf_hash:
                martian.exit("{} (hashes: {} != {})".format(incompat_msg("annotation GTF"), global_gtf_hash, mol_gtf_hash))

            mol_feature_ref = counter.feature_reference
            if global_feature_ref is None:
                global_feature_ref = mol_feature_ref
            elif global_feature_ref != mol_feature_ref:
                martian.exit(incompat_msg("feature reference"))

            chemistry = counter.get_metric('chemistry_name')
            chemistries.add(chemistry)

            if counter.is_aggregated():
                martian.exit("Molecule file was aggregated from multiple samples: %s.\n Aggregated outputs cannot be re-aggregated, please pass each of the original samples instead." % mol_h5)

            if counter.nrows() == 0:
                martian.exit("Cannot aggregate file because it contains no data: %s.\n Please remove this file from the aggregation and try again." % mol_h5)

            for (lib_key, metrics) in counter.get_metric(cr_mol_counter.LIBRARIES_METRIC).iteritems():
                lib_total_reads = metrics[cr_mol_counter.TOTAL_READS_METRIC]
                if lib_total_reads == 0:
                    martian.exit("Library %d has zero reads in file: %s\n Please re-run `cellranger count` without including this gem group." % (lib_key, mol_h5))
