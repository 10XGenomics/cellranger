#
# Copyright (c) 2023 10X Genomics, Inc. All rights reserved.
#
"""Hash function with no other dependencies."""
from __future__ import annotations

import hashlib


def compute_hash_of_file(filename, block_size_bytes=2**20):
    digest = hashlib.sha1()
    with open(filename, "rb") as f:
        for chunk in iter(lambda: f.read(block_size_bytes), b""):
            digest.update(chunk)
    return digest.hexdigest()
