#!/usr/bin/env python
#
# Copyright (c) 2018 10X Genomics, Inc. All rights reserved.
#
#
# Code related to the Python version of the SampleDef passed around inside the
# pipeline. Any code related to the argshim doesn't belong here.

from __future__ import annotations

from enum import Enum


# FastqMode enum from lib/rust/cr_types/src/sample_def.rs
class FastqMode(Enum):
    BCL_PROCESSOR = "BCL_PROCESSOR"
    ILMN_BCL2FASTQ = "ILMN_BCL2FASTQ"
    ONT = "ONT"
    PACBIO = "PACBIO"
    SBX = "SBX"


class SampleDef:
    """Sample definitions."""

    def __init__(
        self,
        fastq_mode: FastqMode,
        read_path: str,
        gem_group: int | None = None,
        lanes: list[int] | None = None,
        library_type: str | None = None,
        r1_length: int | None = None,
        r2_length: int | None = None,
        sample_indices: list[str] | None = None,
        sample_names: list[str] | None = None,
        subsample_rate: float | None = None,
        fastq_id: str | None = None,
    ):
        self.fastq_mode = fastq_mode
        self.gem_group = gem_group
        self.lanes = lanes
        self.library_type = library_type
        self.r1_length = r1_length
        self.r2_length = r2_length
        self.sample_indices = sample_indices
        self.sample_names = sample_names
        self.subsample_rate = subsample_rate
        self.fastq_id = fastq_id
        self.read_path = read_path

    def to_dict(self) -> dict:
        """Convert the BaseSampleDef instance to a dictionary with fastq_mode as a string."""
        return {
            "fastq_mode": (
                self.fastq_mode.value if isinstance(self.fastq_mode, FastqMode) else self.fastq_mode
            ),
            "gem_group": self.gem_group,
            "lanes": self.lanes,
            "library_type": self.library_type,
            "r1_length": self.r1_length,
            "r2_length": self.r2_length,
            "sample_indices": self.sample_indices,
            "sample_names": self.sample_names,
            "subsample_rate": self.subsample_rate,
            "fastq_id": self.fastq_id,
            "read_path": self.read_path,
        }

    @classmethod
    def from_dict(cls, data: dict):
        """Create a SampleDef instance from a dictionary."""
        data = data.copy()
        obj = cls(
            fastq_mode=FastqMode(data.pop("fastq_mode")),
            read_path=data.pop("read_path"),
            gem_group=data.pop("gem_group", None),
            lanes=data.pop("lanes", None),
            library_type=data.pop("library_type", None),
            r1_length=data.pop("r1_length", None),
            r2_length=data.pop("r2_length", None),
            sample_indices=data.pop("sample_indices", None),
            sample_names=data.pop("sample_names", None),
            subsample_rate=data.pop("subsample_rate", None),
            fastq_id=data.pop("fastq_id", None),
        )

        if data:
            unexpected_keys = ", ".join(data.keys())
            raise ValueError(f"Unexpected fields in data: {unexpected_keys}")
        else:
            return obj
