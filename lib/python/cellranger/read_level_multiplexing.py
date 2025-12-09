#
# Copyright (c) 2024 10X Genomics, Inc. All rights reserved.
#

"""Convenience functions for processing read level multiplexing data."""

from cellranger.barcodes.utils import load_probe_barcode_map
from cellranger.fast_utils import MultiGraph
from cellranger.targeted.rtl_multiplexing import get_probe_bc_whitelist


# Temporarily copied here. Need to move these functions to lib so I can import it.
def get_overhang_bc_defn(barcode_def):
    """Returns the offset and length for overhang multiplexing."""
    offset, length = None, None
    for barcode in barcode_def:
        if barcode["kind"] == "overhang":
            offset, length = (barcode["offset"], barcode["length"])
    return offset, length


def get_sample_tag_barcodes(
    multi_graph: MultiGraph, barcode_def: list[dict[str, dict[str, str]]]
) -> dict[str, list[str]]:
    """Return a mapping from sample ID to the tag barcodes associated with them."""
    probe_bc_wl_spec = get_probe_bc_whitelist(barcode_def)
    wl_map = load_probe_barcode_map(
        name=probe_bc_wl_spec["name"],
        path=probe_bc_wl_spec["translation_whitelist_path"],
        strand=probe_bc_wl_spec["strand"],
    )
    return {
        sample_id: [wl_map[tag_name] for tag_name in tag_names]
        for sample_id, tag_names in multi_graph.sample_tag_ids().items()
    }
