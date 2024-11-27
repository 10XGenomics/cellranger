#!/usr/bin/env python
#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#
# Functions for preflight checks
#

from __future__ import annotations

import os
import platform
import re
import resource
import socket
import subprocess
from collections.abc import Callable, Iterable, Mapping
from ctypes import c_int, cdll
from typing import Literal, NamedTuple, TypedDict

import martian

import tenkit.reference
import tenkit.sample_index as tk_si

MIN_PROCESS_NOFILE: int = 1024
MIN_GLOBAL_NOFILE: int = 2**15
GLOBAL_NOFILE_PATH: str = "/proc/sys/fs/file-max"

OS_SUPPORT_URL = "https://support.10xgenomics.com/os-support"


def _search_file_for_re(
    filename: os.PathLike[str] | os.PathLike[bytes] | bytes | int | str, regex: str
) -> str | None:
    """Helper function for warn_deprecated_os - searches a file for a regex pattern.

    returns the line if pattern is found otherwise None
    returns None if an I/O error occurs as we expect to probe non-existent files
    """
    if os.path.exists(filename):
        searcher = re.compile(regex)
        try:
            for line in open(filename):
                if searcher.search(line):
                    return line
        except OSError:
            pass

    return None


def warn_deprecated_os() -> None:
    """Check and warn the user about a deprecated OS being used.

    Direct port of the hints here:
    https://github.com/10XDev/cellranger/blob/b8df29dde96ba65e53a8fa0e9756ee06cac753d6/bin/tenkit/common/_deprecation
    """
    for pat, name in [
        # deprecated RedHat derivatives
        (r" [56]\.", "/etc/redhat-release"),
        (r" [56]\.", "/etc/rocks-release"),
        (r" [56]\.", "/etc/os-release"),
        (r" [56]\.", "/etc/system-release"),
        # Ubuntu 13 or earlier
        (r" 1[0-3]\.", "/etc/lsb-release"),
        (r" [1-9]\.", "/etc/lsb-release"),
        # Suse 10 or 11
        (r"SUSE.* 1[01]\b", "/etc/SuSE-release"),
        # Debian 7 or earlier
        (r'PRETTY_NAME="Debian.*\b[1-7]\.', "/etc/os-release"),
    ]:
        line = _search_file_for_re(name, pat)
        if line:
            martian.alarm(
                f"""WARNING: This operating system version is unsupported or will soon be unsupported:
    {line}
Future releases of this pipeline will require a more current system.
For more information, see {OS_SUPPORT_URL}."""
            )
            return

    try:
        check_kernel_version(platform.release())
        try:
            libc_ver = os.confstr("CS_GNU_LIBC_VERSION")
        except ValueError:
            return
        check_libc_version(libc_ver)
        check_cpu_features()
    except OsVersionError as ex:
        martian.alarm(str(ex))


class OsVersionError(EnvironmentError):
    """An exception notifying that the operating system is out of date."""


def check_kernel_version(release: str) -> None:
    """Check that the kernel version is sufficiently new.

    Args:
        release (str): The kernel release, as returned by platform.release().

    Raises:
        OsVersionError: The kernel version is not sufficiently new.
    """
    if not release:
        return
    # Only check the first two numbers as custom kernels, e.g. "5.10+pk11.11", might not match
    # See CELLRANGER-6505
    version = re.match(r"\D*(\d+)\.(\d+)", release)
    if not version:
        return
    try:
        version = (int(version[1]), int(version[2]))
    except (ValueError, TypeError):  # TypeError if no regex match
        return

    if version < (3, 10):
        raise OsVersionError(
            f"""
WARNING: The kernel used by this operating system version is unsupported:
    {release}
This release requires kernel version 3.10.0 or higher.
For more information, see {OS_SUPPORT_URL}."""
        )


def check_libc_version(libc_ver: str | None) -> None:
    """Check that the libc version is sufficiently new.

    Args:
        libc_ver (str): The libc version string as returned by
                        os.confstr("CS_GNU_LIBC_VERSION")

    Raises:
        OsVersionError: The kernel version is not sufficiently new.
    """
    if not libc_ver:
        return
    libc_version = re.match(r"\D*(\d+)\.(\d+)", libc_ver)
    if not libc_version:
        return
    try:
        libc_version = (int(libc_version[1]), int(libc_version[2]))
    except ValueError:
        return
    if libc_version < (2, 17):
        raise OsVersionError(
            f"""
    WARNING: The glibc version of this operating system version is unsupported:
        {libc_ver}
    This release requires libc version 2.17 or higher.
    For more information, see {OS_SUPPORT_URL}."""
        )


def check_cpu_features() -> None:
    so_file_name = os.path.join(os.path.dirname(__file__), "cpu_support.so")
    try:
        # Load the ssw library using ctypes
        lib = cdll.LoadLibrary(so_file_name)
    except OSError as err:
        raise ImportError("Error loading cpu_support.so") from err
    cpu_supports_sse4_2 = lib.cpu_supports_sse4_2
    cpu_supports_sse4_2.restype = c_int

    if not cpu_supports_sse4_2():
        raise OsVersionError(
            f"""
WARNING: The current CPU does not support sse4.2 instructions.

For more information, see
{OS_SUPPORT_URL}."""
        )

    cpu_supports_popcnt = lib.cpu_supports_popcnt
    cpu_supports_popcnt.restype = c_int
    if not cpu_supports_popcnt():
        raise OsVersionError(
            f"""
WARNING: The current CPU does not support sse4.2 instructions.

For more information, see
{OS_SUPPORT_URL}."""
        )

    cpu_supports_avx = lib.cpu_supports_avx
    cpu_supports_avx.restype = c_int
    if not cpu_supports_avx():
        raise OsVersionError(
            f"""
WARNING: The current CPU does not support avx instructions.

Future versions of 10X Genomics software will not support CPUs older than
Intel Xeon E3 (Sandy Bridge) or AMD Opteron FX (circa 2011).

For more information, see
{OS_SUPPORT_URL}."""
        )


def is_int(s) -> bool:
    try:
        int(s)
    except ValueError:
        return False
    return True


def check_file(file_type, file_path: str, hostname) -> None:
    if not file_path.startswith("/"):
        martian.exit(f"Specified {file_type} file must be an absolute path: {file_path}")
    if not os.path.exists(file_path):
        martian.exit(
            f"On machine: {hostname}, specified {file_type} file does not exist: {file_path}"
        )
    if os.path.isdir(file_path):
        martian.exit(f"On machine: {hostname}, specified {file_type} file is a folder: {file_path}")
    if not os.access(file_path, os.R_OK):
        martian.exit(
            f"On machine: {hostname}, specified {file_type} file is not readable: {file_path}"
        )


def check_folder(folder_type, folder_path: str, hostname, permission=os.X_OK) -> None:
    if not folder_path.startswith("/"):
        martian.exit(f"Specified {folder_type} folder must be an absolute path: {folder_path}")
    if not os.path.exists(folder_path):
        martian.exit(
            f"On machine: {hostname}, specified {folder_type} folder does not exist: {folder_path}"
        )
    if not os.path.isdir(folder_path):
        martian.exit(
            f"On machine: {hostname}, specified {folder_type} path is not a folder: {folder_path}"
        )
    if not os.access(folder_path, permission):
        martian.exit(
            f"On machine: {hostname}, insufficient permissions on {folder_type} folder: {folder_path}"
        )


def check_folder_or_create(
    folder_type, folder_path: str, hostname, permission: int = os.X_OK
) -> None:
    if not folder_path.startswith("/"):
        martian.exit(f"Specified {folder_type} folder must be an absolute path: {folder_path}")
    if os.path.exists(folder_path):
        if not os.path.isdir(folder_path):
            martian.exit(
                f"On machine: {hostname}, specified {folder_type} path is not a folder: {folder_path}"
            )
        if not os.access(folder_path, permission):
            martian.exit(
                f"On machine: {hostname}, insufficient permissions on {folder_type} folder: {folder_path}"
            )
    else:
        try:
            os.makedirs(folder_path)
        except OSError:
            martian.exit(
                f"On machine: {hostname}, could not create {folder_type} folder: {folder_path}"
            )


def check_rta_complete(folder_path: str) -> str:
    """Check that the RTAComplete.txt file is present.

    :return: path to valid RTAComplete.txt in folder_path
    :rtype: string
    """
    hostname = socket.gethostname()
    check_folder("sequencing run", folder_path, hostname)
    rta_complete = os.path.join(folder_path, "RTAComplete.txt")
    if not os.path.exists(rta_complete):
        martian.exit(
            f"On machine: {hostname}, run does not appear to be complete yet.  RTAComplete.txt not found."
        )
    return rta_complete


def check_runinfo_xml(folder_path: str) -> str:
    """Checks that the RunInfo.xml file is present and readable.

    :return: path to valid RunInfo.xml in folder_path
    :rtype: string
    """
    hostname = socket.gethostname()
    check_folder("sequencing run", folder_path, hostname)
    runinfo = os.path.join(folder_path, "RunInfo.xml")
    if not os.path.exists(runinfo):
        martian.exit(
            f"On machine: {hostname}, RunInfo.xml not found. Cannot verify run was 10X-prepped."
        )
    if not os.access(runinfo, os.R_OK):
        martian.exit(f"On machine: {hostname}, insufficient permission to open RunInfo.xml.")
    return runinfo


def check_refdata(reference_path: str, max_contigs=None) -> tuple[bool, str]:
    hostname = socket.gethostname()

    # Determine if the reference package is a known 10X reference package
    genome = tenkit.reference.get_genome(reference_path)
    known_genome = False

    if genome is not None and tenkit.reference.is_tenx(reference_path):
        known_genome = True

        version_path = os.path.join(reference_path, "version")
        if not os.path.exists(version_path):
            return (
                False,
                f"Your reference does not contain the expected files, or they are not readable. Please check your reference folder on {hostname}.",
            )

        # Known genomes get a more stringent check
        if not os.path.exists(os.path.join(reference_path, "fasta/")) or not os.path.exists(
            os.path.join(reference_path, "genes/")
        ):
            return (
                False,
                f"Your reference does not contain the expected files, or they are not readable. Please check your reference folder on {hostname}.",
            )
    elif not os.path.exists(os.path.join(reference_path, "fasta/")):
        return (
            False,
            f"Your reference does not contain the expected files, or they are not readable. Please check your reference folder on {hostname}.",
        )

    if not os.path.exists(
        os.path.join(reference_path, "fasta/genome.fa.flat")
    ) or not os.path.exists(os.path.join(reference_path, "fasta/genome.fa.gdx")):
        return False, "Your reference doesn't appear to be indexed. Please run the mkref tool."

    if os.path.getmtime(os.path.join(reference_path, "fasta/genome.fa.flat")) < os.path.getmtime(
        os.path.join(reference_path, "fasta/genome.fa")
    ) or os.path.getmtime(os.path.join(reference_path, "fasta/genome.fa.gdx")) < os.path.getmtime(
        os.path.join(reference_path, "fasta/genome.fa")
    ):
        msg = """Timestamps suggest that the reference FASTA file was modified after the FASTA
index files were created. If this is the case, {}{}.
If the reference files were not deliberately modified, please run

touch {}

to fix the issue. If you do not have write permissions, please contact your
system administrator.""".format(
            (
                "please reinstall the 10X refdata\ntar file on "
                if known_genome
                else "please reindex your reference\nusing the mkref tool"
            ),
            hostname if known_genome else "",
            os.path.join(reference_path, "fasta", "genome.fa.*"),
        )
        return False, msg

    fasta = tenkit.reference.open_reference(reference_path)
    num_contigs = len(fasta)

    if max_contigs is not None and num_contigs > max_contigs:
        return (
            False,
            "Long Ranger supports a maximum of %d reference contigs. Your reference contains %d. Please combine small contigs into a larger contig separated by N's."
            % (max_contigs, num_contigs),
        )

    max_len = max(len(v) for v in fasta.values())

    logging = f"reference path {reference_path} on {hostname} contains genome: {genome!s}."
    logging += "reference contains %d contigs. max contig length: %d." % (num_contigs, max_len)

    if max_len >= (1 << 29):
        return (
            False,
            "Reference contains a contig longer than 536.8Mb (2^29 bp), which is not supported due to limitations of the .bai format. Please split this contig.",
        )

    # Check for ":" in contig names -- this will break things horribly
    has_colons = any(":" in ctg_name for ctg_name in fasta)
    if has_colons:
        return (
            False,
            "Reference names contain colon characters: ':'. References names containing colons are not supported.",
        )

    return True, logging


def check_open_fh() -> tuple[bool, str | None]:
    _, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
    if hard >= 0 and hard < MIN_PROCESS_NOFILE:
        return (
            False,
            "On machine: %s, process open file handle hard limit (%d) is less than %d. Please run 'ulimit -n %d' before restarting the pipeline."
            % (socket.gethostname(), hard, MIN_PROCESS_NOFILE, MIN_PROCESS_NOFILE),
        )

    if not os.path.exists(GLOBAL_NOFILE_PATH):
        return (
            False,
            f"On machine: {socket.gethostname()}, {GLOBAL_NOFILE_PATH} does not exist.",
        )
    with open(GLOBAL_NOFILE_PATH) as f:
        glob_str = f.read().strip()
    if not glob_str.isdigit():
        return (
            False,
            f"On machine: {socket.gethostname()}, {GLOBAL_NOFILE_PATH} contains a non-integer global open file handle limit: {glob_str}.",
        )

    glob = int(glob_str)
    if glob < MIN_GLOBAL_NOFILE:
        return (
            False,
            "On machine: %s, global open file handle limit (%d) is less than %d. Please set the global file handle limit to %d before restarting the pipeline."
            % (socket.gethostname(), glob, MIN_GLOBAL_NOFILE, MIN_GLOBAL_NOFILE),
        )
    return True, None


def check_sample_indices(
    sample_item: Mapping[str, list[str]], sample_index_key: str = "sample_indices"
) -> tuple[list | None, str | None]:
    sample_indices = sample_item[sample_index_key]
    if not sample_indices:
        return None, "Sample indices must be a non-empty list"
    if not isinstance(sample_indices, list):
        return None, "Sample indices must be of type list"

    new_sample_indices = []
    for sample_index in sample_indices:
        if sample_index == "any":
            return ["*"], None
        elif tk_si.SAMPLE_INDEX_MAP.get(sample_index):
            new_sample_indices.extend(tk_si.SAMPLE_INDEX_MAP.get(sample_index, []))
        elif index := tk_si.SAMPLE_DUAL_INDEX_MAP.get(sample_index):
            # pluck out the i7 index
            new_sample_indices.append(index[0])
        elif re.match("^[ACGTacgt]+$", sample_index):
            new_sample_indices.append(sample_index)
        else:
            return (
                None,
                (
                    f"Sample index '{sample_index}' is not valid. Must be one of: any, SI-<number>, "
                    "SI-<plate>-<well coordinate>, 220<part number>, or "
                    "a nucleotide sequence."
                ),
            )

    return new_sample_indices, None


class SampleDef(TypedDict):
    gem_group: int | None


class _SampleDefWithGemGroup(TypedDict):
    gem_group: int


def _make_some_gem_group(
    sample_def: Iterable[SampleDef], all_null: bool
) -> Iterable[_SampleDefWithGemGroup]:
    """If all gem_groups are set to null, then set them all to 1.

    This is broken out into a separate function in order to make type checkers
    work better.
    """
    if all_null:
        for sd in sample_def:
            sd["gem_group"] = 1
    return sample_def


def check_gem_groups(sample_def: Iterable[SampleDef]) -> tuple[bool, str | None]:
    gem_groups = [sd["gem_group"] for sd in sample_def]
    all_null = all(x is None for x in gem_groups)
    all_int = all(isinstance(x, int) for x in gem_groups)

    # Check for self-consistent gem_group settings in the sample_def entries
    if not (all_null or all_int):
        return (
            False,
            "Inconsistent gem_group tags. Please specify all gem_group tags as null, or all gem_group tags with an integer.",
        )

    gem_groups_sorted = sorted(sd["gem_group"] for sd in _make_some_gem_group(sample_def, all_null))

    # Check that numbering starts at 1
    if len(gem_groups_sorted) > 0 and gem_groups_sorted[0] != 1:
        return False, "gem_group numbering must start at 1"

    # Check for non-contiguous gem groups
    prev = 1
    for group in gem_groups_sorted:
        if group - prev > 1:
            return (
                False,
                f"gem_groups must be numbered contiguously. missing groups: {list(range(prev + 1, group))}",
            )
        prev = group

    # Disable multi-GEM well analysis support since it is not enabled in ATAC SLFE
    if len(set(gem_groups_sorted)) > 1:
        return (
            False,
            """Multi-GEM well analysis is not enabled in `cellranger-[atac|arc] count`. Run \
`cellranger-[atac|arc] count` on each library and then run `cellranger-[atac|arc] aggr`.""",
        )
    return True, None


def check_ld_library_path() -> tuple[Literal[False], str] | tuple[Literal[True], None]:
    if os.environ.get("_TENX_LD_LIBRARY_PATH") is None:
        return (
            False,
            "Environment variable $_TENX_LD_LIBRARY_PATH is not set. Please enter the 10X environment before restarting the pipeline.",
        )
    return True, None


class _VersionCmd(NamedTuple):
    name: str
    cmd: Callable[[], bytes]


def _call(cmd) -> Callable[[], bytes]:
    def fun():
        return subprocess.check_output(cmd)

    return fun


def _grep(
    cmd,
    search: bytes,
) -> Callable[[], bytes]:
    assert isinstance(search, bytes)

    def fun():
        output = subprocess.check_output(cmd)
        for line in output.split(b"\n"):
            if re.match(search, line):
                return line
        return b"not found"

    return fun


_PACKAGE_VERSION_CMDS: list[_VersionCmd] = [
    _VersionCmd(name="mro", cmd=_call(["mro", "--version"])),
    _VersionCmd(name="mrp", cmd=_call(["mrp", "--version"])),
    _VersionCmd(name="python", cmd=_call(["python", "--version"])),
    _VersionCmd(
        name="numpy", cmd=_call(["python", "-c", "import numpy; print(numpy.__version__)"])
    ),
    _VersionCmd(
        name="scipy", cmd=_call(["python", "-c", "import scipy; print(scipy.__version__)"])
    ),
    _VersionCmd(
        name="pysam", cmd=_call(["python", "-c", "import pysam; print(pysam.__version__)"])
    ),
    _VersionCmd(name="PyVCF", cmd=_call(["python", "-c", "import vcf; print vcf.VERSION"])),
    _VersionCmd(name="h5py", cmd=_call(["python", "-c", "import h5py; print(h5py.__version__)"])),
    _VersionCmd(
        name="pandas", cmd=_call(["python", "-c", "import pandas; print(pandas.__version__)"])
    ),
    _VersionCmd(name="bwa", cmd=_grep(["bwa"], b"^ *Version")),
    _VersionCmd(name="samtools", cmd=_call(["samtools", "--version"])),
    _VersionCmd(name="freebayes", cmd=_grep(["freebayes", "-h"], b"^version")),
]


def record_package_versions() -> str:
    results: list[str] = []
    for package in _PACKAGE_VERSION_CMDS:
        name = package.name
        cmd = package.cmd

        version = b"not found"
        try:
            version = cmd()
        except:
            pass
        results.append("{}: {}\n".format(name, version.decode(errors="replace")))
    return "".join(results)
