#!/usr/bin/env python3
#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#
# Utilities for manipulating fasta and fastq files
#

from __future__ import annotations

import glob
import gzip
import os
import re
import sys
from typing import IO, TYPE_CHECKING, BinaryIO, Literal, TextIO, overload

import numpy
from six import ensure_binary, ensure_str

from tenkit.constants import (
    BCL_PROCESSOR_FASTQ_MODE,
    ILLUMINA_QUAL_OFFSET,
    ILMN_BCL2FASTQ_FASTQ_MODE,
)

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator
    from io import BufferedReader


class FastqParseError(Exception):
    """Exception for malformed FASTQs."""

    def __init__(self, handle: BinaryIO | gzip.GzipFile, msg: str) -> None:
        if hasattr(handle, "name"):
            super().__init__(f"Error parsing FASTQ file {handle.name}:\n{msg}")
        else:
            # LZ4FrameFile does not have a "name" attr
            super().__init__(f"Error parsing FASTQ records from stream {handle!r}:\n{msg}")


# Parse the ILMN fastq filename to get the read name, lane, read group,
# and S field. We expect a filename of the form
# <path>/<prefix>_S0_L001_R1_001.fastq
# or if --no-lane-splitting is used
# <path>/<prefix>_S0_R1_001.fastq
class IlmnFastqFile:
    def __init__(self, fullpath: str | bytes):
        if isinstance(fullpath, str):
            fullpath = fullpath.encode()
        fn = os.path.basename(fullpath)

        dot_parts = fn.split(b".")
        if dot_parts[-1] == b"fastq":
            name = dot_parts[-2]
        elif len(dot_parts) > 2 and dot_parts[-2] == b"fastq":
            name = dot_parts[-3]
        else:
            raise NameError(f"{fullpath} is not a fastq file")

        all_flds = name.split(b"_")

        if re.findall(rb"_S\w+_L[0-9][0-9][0-9]_", name):
            flds = all_flds[-4:]
            self.prefix = b"_".join(all_flds[:-4])

            self.s = flds[0][1:]
            self.lane = int(flds[1][2:])
            self.read = flds[2]
            self.group = int(flds[3])
        else:
            flds = all_flds[-3:]
            self.prefix = b"_".join(all_flds[:-3])
            self.s = flds[0][1:]
            self.lane = None
            self.read = flds[1]
            self.group = int(flds[2])

        self.filename = fullpath


BCL_PROCESSOR_FILENAME_REGEX: re.Pattern[bytes] = re.compile(
    rb"read-(\w\w)_si-([^_]+)_lane-(\d+)-chunk-(\d+)"
)


class BclProcessorFastqFile:
    """Parse the FASTQ name from the demux stage output.

    Get the read name, lane, sample index, and chunk.
    """

    def __init__(self, fullpath: str | bytes):
        if isinstance(fullpath, str):
            fullpath = fullpath.encode()
        filename: bytes = os.path.basename(fullpath)
        dot_parts = filename.split(b".")
        if dot_parts[-1] == b"fastq":
            name = dot_parts[-2]
        elif len(dot_parts) > 2 and dot_parts[-2] == b"fastq":
            name = dot_parts[-3]
        else:
            raise NameError(f"{fullpath.decode()} is not a fastq file")

        name_parts = BCL_PROCESSOR_FILENAME_REGEX.match(name)
        if not name_parts:
            raise NameError(f"Not a demux output fastq: {fullpath.decode()}")
        self.read = name_parts.group(1)
        self.index = name_parts.group(2)
        self.lane = int(name_parts.group(3), 10)
        self.chunk = int(name_parts.group(4), 10)
        self.filename: bytes = fullpath


def write_read_fastq(
    fastq_file: BinaryIO | IO[bytes] | gzip.GzipFile, name: bytes, seq: bytes, qual: bytes
) -> None:
    """Writes a single read to a fastq file."""
    fastq_file.write(b"@")
    fastq_file.write(ensure_binary(name))
    fastq_file.write(b"\n")
    fastq_file.write(ensure_binary(seq))
    fastq_file.write(b"\n+\n")
    fastq_file.write(ensure_binary(qual))
    fastq_file.write(b"\n")


def write_read_fasta(fasta_file: TextIO, name: bytes, seq: bytes) -> None:
    """Writes a single read to a fastq file."""
    fasta_file.write(">")
    fasta_file.write(ensure_str(name))
    fasta_file.write("\n")
    fasta_file.write(ensure_str(seq))
    fasta_file.write("\n")


def write_read_pair_fastq(
    fastq_file: BinaryIO | IO[bytes],
    name1: bytes,
    seq1: bytes,
    qual1: bytes,
    name2: bytes,
    seq2: bytes,
    qual2: bytes,
):
    """Writes a read-pair to a fastq file for an interleaving format."""
    write_read_fastq(fastq_file, name1, seq1, qual1)
    write_read_fastq(fastq_file, name2, seq2, qual2)


@overload
def read_generator_fastq(
    fastq_file: BinaryIO | BufferedReader | gzip.GzipFile, paired_end: Literal[True] = ...
) -> Iterator[tuple[bytes, bytes, bytes, bytes, bytes, bytes]]: ...


@overload
def read_generator_fastq(
    fastq_file: BinaryIO | BufferedReader | gzip.GzipFile, paired_end: Literal[False] = False
) -> Iterator[tuple[bytes, bytes, bytes]]: ...


@overload
def read_generator_fastq(
    fastq_file: BinaryIO | BufferedReader | gzip.GzipFile, paired_end: bool = False
) -> (
    Iterator[tuple[bytes, bytes, bytes]] | Iterator[tuple[bytes, bytes, bytes, bytes, bytes, bytes]]
): ...


def read_generator_fastq(
    fastq_file: BinaryIO | BufferedReader | gzip.GzipFile, paired_end: bool = False
) -> (
    Iterator[tuple[bytes, bytes, bytes]] | Iterator[tuple[bytes, bytes, bytes, bytes, bytes, bytes]]
):
    """Returns an interator over a fastq file tha produces (name, seq, qual).

    If paired_end, returns both reads (assuming interleaving fastq)
    """
    line = iter(fastq_file)

    def read() -> tuple[bytes, bytes, bytes]:
        line1: bytes = ensure_binary(next(line).strip())
        if not line1.startswith(b"@"):
            raise FastqParseError(
                fastq_file,
                f"Header line = '{line1.decode()}' does not begin with '@'",
            )
        name = line1[1:]
        seq: bytes = ensure_binary(next(line).strip())
        plus: bytes = ensure_binary(next(line).strip())
        if not plus.startswith(b"+"):
            raise FastqParseError(
                fastq_file,
                f"Line 3 of FASTQ record = '{ensure_str(plus)}' does not begin with '+'",
            )
        qual = next(line).strip()

        if len(seq) != len(qual):
            raise FastqParseError(
                fastq_file,
                "FASTQ record has sequence and quality strings of unequal length:\n"
                f"seq  = {seq.decode()}\n"
                f"qual = {ensure_str(qual)}\n",
            )
        return (name, seq, qual)

    while True:
        if paired_end:
            try:
                name1, seq1, qual1 = read()
                name2, seq2, qual2 = read()
            except StopIteration:
                return
            if name1.split()[0] != name2.split()[0]:
                raise FastqParseError(
                    fastq_file,
                    "Read 1 and Read 2 header prefixes are mismatched:\n"
                    f"Read 1 header prefix = '{name1.split()[0]}'\n"
                    f"Read 2 header prefix = '{name2.split()[0]}'\n",
                )
            yield (name1, seq1, qual1, name2, seq2, qual2)
        else:
            try:
                yield read()
            except StopIteration:
                return


def uninterleave_fastq(
    in_fastq: BinaryIO,
    out_fastq1: BinaryIO | IO[bytes],
    out_fastq2: BinaryIO | IO[bytes],
):
    """Takes an interleaved fastq and outputs two fastqs with the reads split."""
    gen = read_generator_fastq(in_fastq, paired_end=True)
    for name1, seq1, qual1, name2, seq2, qual2 in gen:
        write_read_fastq(out_fastq1, name1, seq1, qual1)
        write_read_fastq(out_fastq2, name2, seq2, qual2)


@overload
def get_qvs(qual: None) -> None: ...


@overload
def get_qvs(qual: str | bytes) -> numpy.ndarray[int, numpy.dtype[numpy.byte]]: ...


def get_qvs(qual: str | bytes | None) -> numpy.ndarray[int, numpy.dtype[numpy.byte]] | None:
    if qual is None:
        return None

    return numpy.fromstring(qual, dtype=numpy.byte) - ILLUMINA_QUAL_OFFSET


@overload
def get_bases_qual(qual: None, cutoff) -> None: ...


@overload
def get_bases_qual(qual: str | bytes, cutoff: int) -> int: ...


def get_bases_qual(qual: str | bytes | None, cutoff: int) -> int | None:
    if qual is None:
        return None

    qvs = numpy.fromstring(qual, dtype=numpy.byte) - ILLUMINA_QUAL_OFFSET
    return numpy.count_nonzero(qvs[qvs > cutoff])


@overload
def get_min_qual(qual: None) -> None: ...


@overload
def get_min_qual(qual: str | bytes) -> int: ...


def get_min_qual(qual: str | bytes | None) -> int | None:
    if qual is None or len(qual) == 0:
        return None

    return numpy.fromstring(qual, dtype=numpy.byte).min() - ILLUMINA_QUAL_OFFSET


@overload
def get_expected_errors(qual: None) -> None: ...


@overload
def get_expected_errors(qual: str | bytes) -> float: ...


def get_expected_errors(qual: str | bytes | None) -> float | None:
    if qual is None or len(qual) == 0:
        return None

    qvs = numpy.fromstring(qual, dtype=numpy.byte) - ILLUMINA_QUAL_OFFSET
    perr = 10.0 ** (-qvs / 10.0)
    return perr.sum()


def find_input_fastq_files_10x_preprocess(
    path: str | bytes,
    read_type: str | bytes,
    sample_index: str | bytes,
    lanes: list[str | int] | None,
    maxNs: int = 2,
) -> list[bytes]:
    """Find fastq files with a matching sample index.

    Only finds files which obey the demultiplex filename scheme,
    with a limited number of Ns in the sample index sequence.
    """
    return _find_input_fastq_files_10x_preprocess(
        ensure_binary(path),
        ensure_binary(read_type),
        ensure_binary(sample_index),
        lanes,
        maxNs,
    )


def _find_input_fastq_files_10x_preprocess(
    path: bytes,
    read_type: bytes,
    sample_index: bytes,
    lanes: list[str | int] | None,
    maxNs: int = 2,
) -> list[bytes]:
    """Find fastq files with a matching sample index.

    This is the internal implementation with sanitized inputs.
    """
    # In the case where this read wasn't taken, the glob won't match
    # anything, return no files

    # We want to pull in all sample indices that match in all non-N positions, with up to 2 Ns
    # construct the appropriate glob, then filter.
    if sample_index != b"*":
        si_glob = b"".join(b"[%cN]" % bse for bse in sample_index)
    else:
        si_glob = b"*"
        # Don't worry about Ns if we are just taking everything
        maxNs = 100

    if not lanes:
        glb = os.path.join(glob.escape(path), b"read-%s_si-%s_*.fastq*" % (read_type, si_glob))
        files = glob.glob(glb)
    else:
        files = [
            f
            for lane in lanes
            for f in glob.iglob(
                os.path.join(
                    glob.escape(path),
                    rb"read-%s_si-%s_lane-%03d[_\-]*.fastq*" % (read_type, si_glob, int(lane)),
                )
            )
        ]

    good_files = []
    # filter files to remove those with > 2 Ns in the sample index
    for f in files:
        m = re.match(b".*si-([A-Z]*)_", f)
        assert m, f"File '{f}' doesn't match the expected naming convention"
        si: bytes = m.groups()[0]
        num_Ns = sum(1 for x in si if x == ord(b"N"))
        if num_Ns <= maxNs:
            good_files.append(f)

    files = sorted(good_files)
    return files


def find_input_fastq_files_bcl2fastq_demult(
    path: str | bytes,
    read_type: str | bytes,
    sample: str | bytes | None,
    lanes: list[str | int] | None,
) -> list[bytes]:
    """Find fastq files demultiplex by bcl2fastq 2.17.

    We glob over all subdirectories beneath
    'path', for fastq files prefixed by 'sample', in the given set of 'lanes', or all
    lanes if 'lanes' is None
    """
    return _find_input_fastq_files_bcl2fastq_demult(
        ensure_binary(path),
        ensure_binary(read_type),
        ensure_binary(sample) if sample is not None else sample,
        lanes,
    )


def _find_input_fastq_files_bcl2fastq_demult(
    path: bytes,
    read_type: bytes,
    sample: bytes | None,
    lanes: list[str | int] | None,
) -> list[bytes]:
    """Find fastq files demultiplex by bcl2fastq 2.17.

    Internal implementation with sanitized inputs.
    """

    def glob_files(path: bytes, file_pattern: bytes) -> list[bytes]:
        assert isinstance(path, bytes)
        assert isinstance(file_pattern, bytes)
        # sample sheet case (Project/Samples/fastq)
        glb = os.path.join(path, b"*", file_pattern)
        files = glob.glob(glb)
        # direct folder case
        if not files:
            glb = os.path.join(path, file_pattern)
            files = glob.glob(glb)
        return files

    # In the case where this read wasn't taken, the glob won't match
    # anything, return no files

    # We want to pull in all sample indices that match in all non-N positions, with up to 2 Ns
    # construct the appropriate glob, then filter.
    if sample is None:
        sample = b"*"

    if lanes is None or lanes == []:
        # when --no-lane-splitting is used the FASTQ files are missing the LXXX piece in the file
        # name. So don't match on the LXXX piece of the file name
        file_pattern = b"%s_*_%s_[0-9][0-9][0-9].fastq*" % (sample, read_type)
        files = glob_files(path, file_pattern)
    else:
        # if lanes are specified then we assume that --no-lane-splitting was not used and add the
        # LXXX piece to the glob pattern.
        files = []
        for lane in lanes:
            file_pattern = b"%s_*_L%03d_%s_[0-9][0-9][0-9].fastq*" % (sample, int(lane), read_type)
            files.extend(glob_files(path, file_pattern))

    files = sorted(files)

    # TENKIT-91 fix
    #
    # glob is limited (e.g., you can't tell it to match a non-zero amount of
    # characters that don't contain an underscore); so use the file name parser functionality
    # to see if you really matched something
    #
    # For example, with the above expression, you would match "target" and "target_not_wanted"
    # even if you only wanted "target".  One could limit the glob match to S0, but if you had a
    # valid sample with a S\d+ suffix, that may interfere as well.
    if sample != b"*":
        files = [f for f in files if os.path.basename(IlmnFastqFile(f).prefix) == sample]
    return files


def find_input_file_type_with_samples(path: str | bytes) -> tuple[str | None, list[bytes]]:
    """Try to determine the demux type of the files at the specified path.

    Also determine the unique sample names of the FASTQs, if applicable.

    Return either "BCL_PROCESSOR" or "ILMN_BCL2FASTQ" as the first argument.
    If ILMN_BCL2FASTQ is returned, also return an array of detected prefixes;
    the array should be empty if the detected mode is BCL_PROCESSOR mode.

    If files if neither variety are found, return None in the first parameter and
    a blank array as the sample list.

    :rtype: (str, list[str])
    """
    files = find_input_fastq_files_10x_preprocess(path, "RA", "*", None)
    if files:
        return BCL_PROCESSOR_FASTQ_MODE, []

    files = find_input_fastq_files_bcl2fastq_demult(path, "R1", None, None)
    if not files:
        return None, []

    ilmn_files = [IlmnFastqFile(f) for f in files]
    samples = {os.path.basename(f.prefix) for f in ilmn_files}
    return ILMN_BCL2FASTQ_FASTQ_MODE, sorted(samples)


class AmbiguousValueError(ValueError):
    """Value error signaling that the arguments need to be more specific."""


def _get_fastq_prefixes(
    fastqprefix: None | str | bytes | Iterable[str | bytes],
) -> tuple[list[bytes] | None, bytes]:
    if fastqprefix is None:
        fastqprefixes = None
        fastqprefix_outstr = b"any"
    elif isinstance(fastqprefix, str):
        fastqprefix = fastqprefix.encode()
        fastqprefixes = [fastqprefix]
        fastqprefix_outstr: bytes = fastqprefix
    elif isinstance(fastqprefix, bytes):
        fastqprefixes = [fastqprefix]
        fastqprefix_outstr: bytes = fastqprefix
    else:
        fastqprefixes = [ensure_binary(f) for f in fastqprefix]
        fastqprefix_outstr = b",".join(fastqprefixes)

    return fastqprefixes, fastqprefix_outstr


def check_fastq_types(
    path: str | bytes, fastqprefix: None | str | bytes | Iterable[str | bytes]
) -> tuple[str, bytes]:
    """Validate that the path and fastqprefix arguments are compatible and allowed.

    Args:
        path:  The supplied input path argument.
        fastqprefix:  The supplied --fastqprefix sample prefix argument.
            Can be a single string or collection.

    Returns:
        The input_mode and sample_name value that should be passed into the MRO.

    Raises:
        ValueError: If the supplied path + prefix argument are invalid.
    """
    fastqprefixes, fastqprefix_outstr = _get_fastq_prefixes(fastqprefix)

    demux_type, samples = find_input_file_type_with_samples(path)
    if not demux_type:
        raise ValueError(
            """No input FASTQs were found for the requested parameters.

If your files came from bcl2fastq or mkfastq:
 - Make sure you are specifying the correct --sample(s), i.e. matching the sample sheet
 - Make sure your files follow the correct naming convention, e.g. SampleName_S1_L001_R1_001.fastq.gz (and the R2 version)
 - Make sure your --fastqs points to the correct location.

Refer to the "Specifying Input FASTQs" page at https://support.10xgenomics.com/ for more details.

"""
        )
    if demux_type == BCL_PROCESSOR_FASTQ_MODE and fastqprefix:
        raise ValueError(
            "Cannot use --fastqprefix argument with FASTQs generated by the demux pipeline."
        )
    if demux_type == ILMN_BCL2FASTQ_FASTQ_MODE:
        if len(samples) > 1:
            # ambiguous fastqprefix case
            if not fastqprefix:
                samples_list = b"\n".join(samples)
                raise AmbiguousValueError(
                    f"The --sample argument must be specified if multiple samples were demultiplexed in a run folder.  Options:\n{samples_list}"
                )
            # no overlap case
            elif not set(samples).intersection(fastqprefixes or ()):
                raise ValueError(
                    f"Samples not detected among demultiplexed FASTQs: {fastqprefix_outstr.decode()}"
                )
            # some overlap; legal fastqprefix case
            else:
                return ILMN_BCL2FASTQ_FASTQ_MODE, fastqprefix_outstr
        # single sample does not match fastqprefixes case
        elif fastqprefixes is not None and samples[0] not in fastqprefixes:
            raise ValueError(f"Samples not detected among FASTQs: {fastqprefix_outstr.decode()}")
        # no fastqprefix case-- return the lone sample
        elif fastqprefixes is None:
            return ILMN_BCL2FASTQ_FASTQ_MODE, samples[0]
        # legal fastqprefix list case
        else:
            return ILMN_BCL2FASTQ_FASTQ_MODE, fastqprefix_outstr
    else:
        return BCL_PROCESSOR_FASTQ_MODE, b"any"


def check_fastq_types_multipath(
    fastq_paths: Iterable[str | bytes], fastqprefix: str | bytes
) -> tuple[str, bytes]:
    """Validate that at least one path is compatible with the set of fastqprefix arguments.

    Return the input mode and superset of sample names that
    should be used downstream.  Forces the input mode to be the same for each
    path; if not, an error will be returned.  Will only raise an error in this
    case, and if no path-prefix combination will select any FASTQs.

    If, like in Cell Ranger, it is legal to have chunk-level fastq modes
    instead of pipeline-level fastq modes, it is preferred to call
    check_fastq_types individually over each path to get the input mode and
    sample name(s) for each chunk.
    """
    input_modes = set()
    sample_names = set()
    error_messages = set()
    for path in fastq_paths:
        try:
            input_mode, samples = check_fastq_types(path, fastqprefix)
            input_modes.add(input_mode)
            sample_names.update(samples.split(b","))
        # don't handle ambiguous values if detected in individual paths
        except AmbiguousValueError as ex:
            raise ex
        except ValueError as ex:
            sys.stderr.write(
                f"Invalid path/prefix combination: {path}, {fastqprefix.decode() if isinstance(fastqprefix, bytes) else str(fastqprefix)}\n"
            )
            error_messages.add(str(ex))

    # happens if there were no legal selections
    if len(input_modes) == 0:
        if len(error_messages) == 1:
            raise ValueError(error_messages.pop())
        else:
            raise ValueError("FASTQ selection errors:\n{}".format("\n".join(error_messages)))
    elif len(input_modes) > 1:
        raise ValueError(
            "Cannot process FASTQs at same time from different demultiplexing methods."
        )
    # can happen if multiple paths have different samples in them
    elif not fastqprefix and len(sample_names) > 1:
        raise AmbiguousValueError(
            b"The --sample argument must be specified if multiple samples were demultiplexed across the specified run folders.  Options:\n%s"
            % b"\n".join(sorted(sample_names))
        )
    else:
        return input_modes.pop(), b",".join(sorted(sample_names))


def get_run_data(fn: str | bytes) -> tuple[bytes, bytes]:
    """Parse flowcell + lane from the first FASTQ record.

    NOTE: we don't check whether there are multiple FC / lanes in this file.
    NOTE: taken from longranger/mro/stages/reads/setup_chunks
    """
    if isinstance(fn, str):
        fn = fn.encode()
    with gzip.open(fn, "rb") if fn[-2:] == b"gz" else open(fn, "rb") as reader:
        gen = read_generator_fastq(reader)

        try:
            name = next(gen)[0]
        except StopIteration:
            # empty fastq
            raise ValueError(
                f"Could not extract flowcell and lane from FASTQ file. File is empty: {fn.decode()}"
            )

    # Abort if this is not an Illumina-like QNAME
    match = re.search(b"^([^:]+):([^:]+):([^:]+):([^:]+)", name)
    if match:
        flowcell, lane = match.group(3, 4)
    else:
        flowcell, lane = b"", b""

    return (flowcell, lane)
