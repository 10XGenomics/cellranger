#
# Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
#
"""A file to read and ingest cell_barcodes.json files produced by PARSE_MULTI_CONFIG."""

from __future__ import annotations

import json

from cellranger.barcodes.utils import load_barcode_whitelist
from cellranger.chemistry import get_whitelist_name_from_chemistry_description


class InvalidCellBarcode(Exception):
    """Barcode is not from a given whitelist."""


class CellBarcodes(list):
    """Python equivalent of file produced by Rust code, list of barcodes as bytes."""

    def __init__(self, json_file_name):
        """Load a JSON file with a list of barcodes.

        e.g. like

        .. code-block:: python

            to_keep = [
                "AACAAAGAGCAACCAG-1",
                "AACGAAAGTGACCGAA-1",
                "AAGCGAGCAGTTGAAA-1",
                "AAGCGTTTCCCGTTCA-1",
                "AAGTGAACAAGATGGC-1",
                "AATCACGTCGTAACTG-1",
                "AATGGAAGTTAAGCAA-1",
                "AATGGAAGTTCCTTGC-1",
                "AATTCCTAGCCATTTG-1"
            ]

        Args:
            json_file_name: json path to load
        """
        super().__init__()
        with open(json_file_name) as f:
            for bc in json.load(f):
                self.append(bc.encode())

    def verify_bcs_match_whitelist(self, chemistry_description: str | None):
        """Given a chemistry description (which corresponds to a "description" in the chemistry_defs.json.

        try to load the whitelist and verify that everything in this list is on that whitelist.

        Args:
            chemistry_description: string that corresponds to a "description" in the chemistry_defs.json

        Returns:
            None or Exception
        """
        if chemistry_description is not None:
            whitelist_name = get_whitelist_name_from_chemistry_description(chemistry_description)
            if whitelist_name is not None:
                # Hard assumption here it's `[whitelist][name]` always
                whitelist = load_barcode_whitelist(whitelist_name["name"], as_set=True)
                for barcode in self:
                    stripped_bc = barcode.split(b"-")[0]
                    if stripped_bc not in whitelist:
                        raise InvalidCellBarcode(
                            f"Barcode: {barcode.decode()} is invalid for the chemistry {chemistry_description}"
                        )
