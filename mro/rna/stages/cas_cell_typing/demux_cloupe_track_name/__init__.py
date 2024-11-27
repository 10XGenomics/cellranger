#
# Copyright (c) 2024 10X Genomics, Inc. All rights reserved.
#
"""Demux cas track name."""

__MRO__ = """
stage DEMUX_CLOUPE_TRACK_NAME(
    in  string cas_track_name_from_user,
    in  string cas_track_name_from_ppln,
    out string cas_track_name,
    src py     "stages/cas_cell_typing/demux_cloupe_track_name",
)
"""


def main(args, outs):
    if args.cas_track_name_from_user:
        outs.cas_track_name = args.cas_track_name_from_user
    else:
        outs.cas_track_name = args.cas_track_name_from_ppln
