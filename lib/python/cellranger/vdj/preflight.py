#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#

from __future__ import annotations

import os
import socket

from six import ensure_binary, ensure_str

import cellranger.vdj.chain_types as chain_types
import cellranger.vdj.constants as vdj_constants
import cellranger.vdj.reference as vdj_ref
from cellranger.preflight import PreflightException
from tenkit.seq import get_rev_comp

# Inner enrichment primers
VDJ_HUMAN_TCR_INNER_PRIMERS = [b"AGTCTCTCAGCTGGTACACG", b"TCTGATGGCTCAAACACAGC"]

VDJ_HUMAN_IG_INNER_PRIMERS = [
    b"GGGAAGTTTCTGGCGGTCA",
    b"GGTGGTACCCAGTTATCAAGCAT",
    b"GTGTCCCAGGTCACCATCAC",
    b"TCCTGAGGACTGTAGGACAGC",
    b"CACGCTGCTCGTATCCGA",
    b"TAGCTGCTGGCCGC",
    b"GCGTTATCCACCTTCCACTGT",
]

VDJ_MOUSE_TCR_INNER_PRIMERS = [b"AGTCAAAGTCGGTGAACAGGCA", b"GGCCAAGCACACGAGGGTA"]

VDJ_MOUSE_IG_INNER_PRIMERS = [
    b"TACACACCAGTGTGGCCTT",
    b"CAGGCCACTGTCACACCACT",
    b"CAGGTCACATTCATCGTGCCG",
    b"GAGGCCAGCACAGTGACCT",
    b"GCAGGGAAGTTCACAGTGCT",
    b"CTGTTTGAGATCAGTTTGCCATCCT",
    b"TGCGAGGTGGCTAGGTACTTG",
    b"CCCTTGACCAGGCATCC",
    b"AGGTCACGGAGGAACCAGTTG",
    b"GGCATCCCAGTGTCACCGA",
    b"AGAAGATCCACTTCACCTTGAAC",
    b"GAAGCACACGACTGAGGCAC",
]

VDJ_HUMAN_INNER_PRIMERS = VDJ_HUMAN_TCR_INNER_PRIMERS + VDJ_HUMAN_IG_INNER_PRIMERS
VDJ_MOUSE_INNER_PRIMERS = VDJ_MOUSE_TCR_INNER_PRIMERS + VDJ_MOUSE_IG_INNER_PRIMERS
VDJ_KNOWN_INNER_PRIMERS = VDJ_HUMAN_INNER_PRIMERS + VDJ_MOUSE_INNER_PRIMERS


def check_refdata(reference_path, denovo):
    if reference_path is None and not denovo:
        raise PreflightException("Must specify --reference unless --denovo is specified.")

    if reference_path is None:
        return

    hostname = socket.gethostname()
    print(f"Checking reference_path ({reference_path}) on {hostname}...")

    required_files = [vdj_constants.REFERENCE_FASTA_PATH]

    for filename in required_files:
        p = os.path.join(reference_path, filename)

        if not os.path.isfile(p):
            raise PreflightException(
                f"Your reference does not contain the expected files, including {p}, or they are not readable. Please check your reference folder on {hostname}."
            )


def check_inner_enrichment_primers(primers_file, reference_path):
    """Check that the path is valid, contains only expected characters (ACGT) and targets C-regions."""
    # 1. Need not specify inner enrichment primers for standard human and mouse VDJ
    if primers_file is None:
        if reference_path is None:
            # If no reference is specified (in denovo mode), make sure primers are specified
            raise PreflightException(
                "You need to specify inner enrichment primers (using --inner-enrichment-primers flag) when a reference is not specified."
            )
        else:
            # Make sure that we find at least one internally designed primer that targets
            # at least one C-Region in the input reference
            for feat in vdj_ref.get_vdj_feature_iter(reference_path):
                assert isinstance(feat.region_type, bytes)
                if feat.region_type == b"C-REGION":
                    for primer in VDJ_KNOWN_INNER_PRIMERS:
                        primer_rc = get_rev_comp(primer)
                        if primer_rc in feat.sequence:
                            return

        raise PreflightException(
            f"Inner enrichment primers are required for species other than human or mouse for which primers are not provided by 10x Genomics. None of the constant regions in the reference ({reference_path}) is targeted by the known primers."
        )

    hostname = socket.gethostname()
    print(f"Checking enrichment primers ({primers_file}) on {hostname}...")

    # 2. If specified, make sure that the path exists
    if not os.path.isfile(primers_file):
        raise PreflightException(
            f"The file specifying inner enrichment primers ({primers_file}), does not exists or is not readable. Please check your path on {hostname}."
        )

    # 3. Make sure that the file is a newline separated list of ACGT sequences
    inner_primers = []
    with open(primers_file) as f:
        for i, line in enumerate(f.readlines()):
            seq = line.strip()
            if len(seq) == 0:
                raise PreflightException(
                    f"Line number {i + 1} in the inner enrichment primers file ({primers_file}) is empty. You should specify a newline separated list of primers."
                )
            for j, base in enumerate(seq):
                if base not in {"A", "C", "G", "T"}:
                    raise PreflightException(
                        f"Inner enrichment primers file ({primers_file}) contain non ACGT characters, which are not supported (Found {base} in line {i + 1}, character {j + 1}). You should specify a newline separated list of primers."
                    )
            inner_primers.append(ensure_binary(seq))

    if not inner_primers:  # Empty file
        raise PreflightException(
            f"Inner enrichment primers file ({primers_file}) contains zero entries. You should specify at least one primer"
        )

    if reference_path:
        # 4. Make sure that every primer targets at least 1 constant region in the reference.
        invalid = []
        for primer in inner_primers:
            # The inner primers are the reverse primers
            primer_rc = get_rev_comp(primer)
            found = False
            for feat in vdj_ref.get_vdj_feature_iter(reference_path):
                if feat.region_type == b"C-REGION":
                    if primer_rc in feat.sequence:
                        found = True
                        break
            if not found:
                invalid.append(primer)

        if invalid:
            invalid_str = b", ".join(invalid)
            raise PreflightException(
                f"None of the C-REGIONs in the reference {reference_path} is targeted by the following inner enrichment primer(s): {ensure_str(invalid_str)}."
            )


def check_chain(chain):
    if ensure_binary(chain) == b"TR_GD" or ensure_binary(chain) == b"TR-GD":
        raise PreflightException("Please use Cellranger multi to run Gamma/Delta TCR datasets.")
    if ensure_binary(chain) not in chain_types.CHAIN_TYPE_SPECS:
        raise PreflightException(
            "Must specify --chain as one of: {}.".format(
                ensure_str(b", ".join(chain_types.CHAIN_TYPE_SPECS))
            )
        )
