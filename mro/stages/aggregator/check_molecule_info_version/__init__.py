#!/usr/bin/env python
#
# Copyright (c) 2018 10x Genomics, Inc. All rights reserved.
#
"""
aggr preflight check + convert legacy molecule info h5 to current version
"""
import martian
import os
import datetime
import h5py
import tables

import cellranger.constants as cr_constants
import cellranger.molecule_counter as cr_mol_counter
import cellranger.io as cr_io

from cellranger.analysis.constants import CBC_MAX_NCELLS

__MRO__ = """
stage CHECK_MOLECULE_INFO_VERSION(
    in  map[] sample_defs,
    out map[] updated_sample_defs,
    src py    "stages/aggregator/check_molecule_info_version",
) split using (
    in  int   mol_h5_version,
    in  map   sample_def,
    out map   updated_sample_def,
)
"""

def split(args):
    files_seen = set()
    chunks = []

    for sample_def in args.sample_defs:

        mol_h5 = sample_def[cr_constants.AGG_H5_FIELD]

        if not os.path.exists(mol_h5):
            martian.exit("Input molecule file does not exist: %s" % mol_h5)

        if not os.access(mol_h5, os.R_OK):
            martian.exit("Input molecule file is not readable, please check file permissions: %s" % mol_h5)

        h5_filetype = cr_io.get_h5_filetype(mol_h5)
        if h5_filetype and h5_filetype != cr_mol_counter.MOLECULE_H5_FILETYPE:
            martian.exit("Input is a %s file, but a molecule file is required" % h5_filetype)

        if mol_h5 in files_seen:
            martian.exit("Same molecule file is specified in multiple sample definitions: %s" % mol_h5)
        else:
            files_seen.add(mol_h5)

        mc_fn = sample_def[cr_constants.AGG_H5_FIELD]
        mc_h5 = h5py.File(mc_fn, 'r')
        try:
            mol_h5_version = mc_h5.attrs[cr_mol_counter.FILE_VERSION_KEY]
        except KeyError:
            martian.exit("The molecule info HDF5 file (%s) was produced by an older version of Cell Ranger. Reading these files is unsupported." % mc_fn)

        if mol_h5_version == 2:
            # CRv1.3-generated files lack essential metrics introduced in CRv2 (force_cells, conf_mapped_filtered_bc_reads)
            #   hence aggregation with these files is unsupported
            mc_table = tables.open_file(mc_fn, 'r')
            mc_attrs = mc_table.get_node("/metrics")._v_attrs
            cr_major_version = int(mc_attrs["cellranger_version"].split(".", 1)[0])
            gem_metrics = dict(mc_attrs["gem_groups"])
            mc_table.close()
            if cr_major_version < 2 or not all('conf_mapped_filtered_bc_reads' in group for group in gem_metrics.itervalues()):
                martian.exit("The molecule info HDF5 file (%s) was produced by an older version of Cell Ranger. Reading these files is unsupported." % mc_fn)

            nrows = mc_h5['barcode'].shape[0]
            mem_gb = cr_mol_counter.MoleculeCounter.estimate_mem_gb(nrows, scale=4)
        else:
            mem_gb = 1

        chunks.append({
            'sample_def': sample_def,
            'mol_h5_version': mol_h5_version,
            '__mem_gb': mem_gb,
        })

    return {'chunks': chunks, 'join': {'__mem_gb': 3}}

def main(args, outs):
    outs.updated_sample_def = args.sample_def.copy()

    if args.mol_h5_version == 2:
        v2_mole_info_h5 = args.sample_def[cr_constants.AGG_H5_FIELD]
        v2_file_basename = os.path.basename(v2_mole_info_h5)
        v3_filename = '{x[0]}_v3_{x[2]}{x[1]}'.format(x=list(os.path.splitext(v2_file_basename))+[datetime.datetime.now().isoformat()])
        out_v3_mole_info_h5 = martian.make_path(v3_filename)

        cr_mol_counter.MoleculeCounter.convert_v2_to_v3(v2_mole_info_h5, out_v3_mole_info_h5)
        outs.updated_sample_def[cr_constants.AGG_H5_FIELD] = out_v3_mole_info_h5

def join(args, outs, chunk_defs, chunk_outs):
    outs.updated_sample_defs = [chunk_out.updated_sample_def for chunk_out in chunk_outs]

    if any("batch" in sample_def for sample_def in outs.updated_sample_defs):
        ncells = 0
        for sample_def in outs.updated_sample_defs:
            with cr_mol_counter.MoleculeCounter.open(sample_def[cr_constants.AGG_H5_FIELD], 'r') as mc:
                library_info = mc.get_library_info()
                ncells += sum(mc.get_num_filtered_barcodes_for_library(i) for i in xrange(len(library_info)))
        if ncells > CBC_MAX_NCELLS:
            martian.exit("You provided {:,} cells in total, but chemistry batch correction only supports up to {:,} cells.".format(
                ncells, CBC_MAX_NCELLS))
