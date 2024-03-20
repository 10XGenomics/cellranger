# Copyright (c) 2023 10X Genomics, Inc. All rights reserved.

"""Clonotype related functions."""


def extract_clonotype_id_from_name(clonotype_name: str):
    """Extract the clonotype id from a clonotype name string."""
    return int(clonotype_name.split("clonotype")[1])
