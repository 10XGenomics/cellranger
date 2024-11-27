# Copyright (c) 2020 10x Genomics, Inc. All rights reserved.
"""Code to provide preflights for AGGR and munging of data.

From multiple versions in the preflight.
"""
from __future__ import annotations

import os

import martian

import cellranger.constants as cr_constants
import cellranger.hdf5 as cr_h5
import cellranger.molecule_counter as cr_mol_counter


def check_molecule_info_version_split(args):
    """The split method for the `CHECK_MOLECULE_INFO_VERSION` stage.

    It verifies the versions of all the input files meet requirements, and if
    necessary, schedules the main stage to update all of them.

    Args:
        args: The args to the stage

    Returns:
        dict: chunks
    """
    # pylint: disable=invalid-name
    files_seen = set()
    chunks = []

    for sample_def in args.sample_defs:
        is_spatial = args.product_type == cr_constants.SPATIAL_PRODUCT_TYPE

        mol_h5 = sample_def[cr_constants.AGG_H5_FIELD]

        if not os.path.exists(mol_h5):
            martian.exit(f"Input molecule file does not exist: {mol_h5}")

        if not os.access(mol_h5, os.R_OK):
            martian.exit(
                f"Input molecule file is not readable, please check file permissions: {mol_h5}"
            )

        if not cr_h5.is_hdf5(mol_h5):
            martian.exit(f"Input molecule file is not a valid HDF5 file: {mol_h5}")

        h5_filetype = cr_h5.get_h5_filetype(mol_h5)
        if h5_filetype and h5_filetype != cr_mol_counter.MOLECULE_H5_FILETYPE:
            martian.exit(f"Input is a {h5_filetype} file, but a molecule file is required")

        if mol_h5 in files_seen:
            martian.exit(
                f"Same molecule file is specified in multiple sample definitions: {mol_h5}"
            )
        else:
            files_seen.add(mol_h5)

        mem_gb = 1
        mc_h5, file_version = cr_mol_counter.get_h5py_file_and_version(mol_h5)
        if file_version < 3 and is_spatial:
            # Space ranger runs are all from >= version 3
            martian.exit(f"The molecule info HDF5 file ({mc_h5}) was not produced by Space Ranger.")
        elif is_spatial:
            # Verify it's a spatial whitelist
            with cr_mol_counter.MoleculeCounter.open(mol_h5, "r") as mol_info:
                if not mol_info.is_spatial_data():
                    martian.exit(
                        "The sample {} is not compatible with spaceranger aggr. Please use a matching pipeline.".format(
                            sample_def["library_id"]
                        )
                    )
        elif file_version < 2:
            martian.exit(
                "The molecule info HDF5 file (%s, version: %d) was produced by an older software version. Reading these files is unsupported."
                % (mc_h5, file_version)
            )
        elif file_version == 2:
            # CRv1.3-generated files lack essential metrics introduced in CRv2 (force_cells, conf_mapped_filtered_bc_reads)
            # hence aggregation with these files is unsupported
            metrics = cr_mol_counter.get_v2_metrics(mol_h5)
            cr_major_version = int(metrics[cr_constants.CELLRANGER_VERSION_KEY].split(".", 1)[0])
            gem_metrics = dict(metrics["gem_groups"])
            metric_in_all = all(
                "conf_mapped_filtered_bc_reads" in group for group in gem_metrics.values()
            )
            if cr_major_version < 2 or not metric_in_all:
                martian.exit(
                    f"The molecule info HDF5 file ({mol_h5}) was produced by an older software version. Reading these files is unsupported."
                )

            # pylint: disable=no-member
            nrows = mc_h5["barcode"].shape[0]
            mem_gb = cr_mol_counter.MoleculeCounter.estimate_mem_gb(nrows, scale=4)

        chunks.append(
            {
                "sample_def": sample_def,
                "mol_h5_version": file_version,
                "__mem_gb": mem_gb,
            }
        )
    return {"chunks": chunks, "join": {"__mem_gb": 3}}
