#!/usr/bin/env python
#
# Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
#

# pylint: disable=not-an-iterable
from __future__ import annotations

import json
import math
import pickle
from collections import OrderedDict, defaultdict
from collections.abc import Callable, Generator
from typing import Any, NamedTuple, TypeVar

import h5py
import numpy as np
from six import ensure_binary

import cellranger.barcodes.utils as bc_utils
import cellranger.h5_constants as h5_constants
import cellranger.hdf5 as cr_h5
import cellranger.rna.library as rna_library
import cellranger.utils as cr_utils

# pylint: disable=no-name-in-module,import-error
from cellranger.fast_utils import concatenate_molecule_infos
from cellranger.feature_ref import FeatureReference
from cellranger.rna.library import LIBRARY_TYPE
from cellranger.targeted.targeted_constants import TARGETING_METHOD_TL
from cellranger.wrapped_tables import tables

_T1 = TypeVar("_T1")

MOLECULE_H5_FILETYPE = "molecule"
FILE_VERSION_KEY = "file_version"

# Define set of valid possible attribute keys for all versions of molecule info HDF5 files.
# Used for validation of proper molecule info HDF5 formatting, to avoid processing other HDF5 files.
VALID_MOL_INFO_H5_ATTRIBUTE_KEYS = [
    "TITLE",
    "CLASS",
    "VERSION",
    "FILTERS",
    "PYTABLES_FORMAT_VERSION",
    h5_constants.H5_FILETYPE_KEY,
    FILE_VERSION_KEY,
]


# Version 6:
# -  adds `probe_idx` int32 dataset and `probes` group to annotate RTL reads with the probe
#    they came from.
#
# Version 5:
# - adds `umi_type` uint32 dataset to distinguish between transcriptomic and non-transcriptomic
#   umi's in intron mode. We only use 1 bit of information and the remaining 31 bits will be used
#   in the future to further annotate umis.
#
CURR_FILE_VERSION = 6

# Group whose attributes store various metadata
# This replaces the 'metrics' group used in FILE_VERSION=3
METRICS_JSON_DATASET_NAME = "metrics_json"

# Legacy from v3 of molecule_info.h5 format
# Group whose attributes store various metadata
V3_METRICS_GROUP_NAME = "metrics"

# Group that tracks which barcodes passed filters (usually means they are
# cell-associated)
BARCODE_INFO_GROUP_NAME = "barcode_info"
PROBE_GROUP_NAME = "probes"
BARCODE_IDX_COL_NAME = "barcode_idx"
FEATURE_IDX_COL_NAME = "feature_idx"
LIBRARY_IDX_COL_NAME = "library_idx"
PROBE_IDX_COL_NAME = "probe_idx"
COUNT_COL_NAME = "count"
GEM_GROUP_COL_NAME = "gem_group"
UMI_COL_NAME = "umi"
UMI_TYPE_COL_NAME = "umi_type"

BARCODE_DS_NAME = "barcodes"

# Bit flag corresponding to the Rust enum cr_types::UmiType::Txomic.
UMI_TYPE_TXOMIC = np.uint32(0x1)

MOLECULE_INFO_DEFAULT_VALUES = {UMI_TYPE_COL_NAME: UMI_TYPE_TXOMIC}

MOLECULE_INFO_COLUMNS = OrderedDict(
    [
        (GEM_GROUP_COL_NAME, np.uint16),  # Up to 65k
        (BARCODE_IDX_COL_NAME, np.uint64),
        (FEATURE_IDX_COL_NAME, np.uint32),  # Up to 4e9 features
        (LIBRARY_IDX_COL_NAME, np.uint16),  # Up to 65k
        # Available in Rust but not python yet
        # making this available should also remove the pass the in `MoleculeCounter.open` method.
        # (PROBE_IDX_COL_NAME, np.uint32),
        (UMI_COL_NAME, np.uint32),  # Up to 16-mers
        (COUNT_COL_NAME, np.uint32),  # Up to 4e9 readpairs/mol
        (UMI_TYPE_COL_NAME, np.uint32),  # Up to 32 bit flags
    ]
)

UNIMPLEMENTED_PROBE_KEYS = [PROBE_GROUP_NAME, PROBE_IDX_COL_NAME]
BARCODES = "barcodes"
MOLECULE_REF_COLUMNS = [BARCODES, "library_info"]

# Preserve contiguity of these when chunking a MoleculeCounter
CHUNK_COLUMNS = ["gem_group", BARCODE_IDX_COL_NAME]

# Common top-level metrics
GEM_GROUPS_METRIC = "gem_groups"
LIBRARIES_METRIC = "libraries"
IS_AGGREGATED_METRIC = "is_aggregated"
TARGETING_METHOD_METRIC = "targeting_method"
MOLECULE_INFO_TYPE_METRIC = "molecule_info_type"
MOLECULE_INFO_TYPE_RAW = "raw"
MOLECULE_INFO_TYPE_COUNT = "count"
ANALYSIS_PARAMETERS_METRIC = "analysis_parameters"

# Per-library metrics
TOTAL_READS_METRIC = "raw_read_pairs"
TOTAL_READS_IN_FILTERED_BARCODES_METRIC = "raw_read_pairs_in_filtered_barcodes"
DOWNSAMPLED_READS_METRIC = "downsampled_reads"
USABLE_READS_METRIC = "usable_read_pairs"
FEATURE_READS_METRIC = "feature_read_pairs"
DOWNSAMPLED_FEATURE_READS_METRIC = "downsampled_feature_read_pairs"
ON_TARGET_USABLE_READS_METRIC = "on_target_usable_read_pairs"
GG_RECOVERED_CELLS_METRIC = "recovered_cells"
GG_FORCE_CELLS_METRIC = "force_cells"
INTRON_MODE_PARAM = "include_introns"
# With CR 7.0 we changed the default to true for non-targeted assays
INTRON_MODE_HISTORIC_DEFAULT = False  # assume this value if key not present
FILTER_PROBES = "filter_probes"
FILTER_PROBES_PARAM_DEFAULT = True
NO_PROBE_FILTER = "no_probe_filter"

HDF5_COMPRESSION = "gzip"
# Number of elements per HDF5 chunk;
#   here, 1 MiB/(32 bytes per element)
HDF5_CHUNK_SIZE = 32768


# Per-barcode metadata. Sparse (not every barcode is listed)
class BarcodeInfo(NamedTuple):
    # Array-ized list of (barcode_idx, library_idx, genome_idx)
    PASS_FILTER_BARCODE_IDX = 0
    PASS_FILTER_LIBRARY_IDX = 1
    PASS_FILTER_GENOME_IDX = 2
    pass_filter: np.ndarray
    # tuples where presence indicates passing filter.
    # This is a binary 3-d sparse matrix in COO format.
    genomes: list[str]  # Genome strings for genome 0..j


BARCODE_INFO_DTYPES = {
    "pass_filter": "uint64",
    "genomes": "str",
}


def _write_barcodes(mc: MoleculeCounter, barcodes):
    """Writes the barcodes to the molecule info.

    Args:
        mc: MoleculeCounter instance
        barcodes: a numpy array of strings or iterable of strings
    """
    if isinstance(barcodes, np.ndarray) and barcodes.dtype.type is np.bytes_:
        mc.h5.create_dataset(BARCODE_DS_NAME, data=barcodes)
    elif len(barcodes) == 0:
        bc_array = np.array([], dtype="S", copy=False)
        mc.h5.create_dataset(
            BARCODE_DS_NAME,
            data=bc_array,
            compression=HDF5_COMPRESSION,
        )
    else:
        # If there are multiple barcode lengths, use the largest for the numpy dtype.
        max_barcode_len = max(len(x) for x in barcodes)
        barcode_dtype = np.dtype("S%d" % max_barcode_len)
        mc.h5.create_dataset(
            BARCODE_DS_NAME,
            data=np.fromiter(barcodes, barcode_dtype, count=len(barcodes)),
            compression=HDF5_COMPRESSION,
        )


def create_dataset(mc: MoleculeCounter, name: str):
    col_type = MOLECULE_INFO_COLUMNS[name]
    mc.columns[name] = mc.h5.create_dataset(
        name,
        (0,),
        maxshape=(None,),
        dtype=col_type,
        compression=HDF5_COMPRESSION,
        chunks=(HDF5_CHUNK_SIZE,),
    )


def get_barcode_index_to_retain(mc: MoleculeCounter, tgt_chunk_len=2000000) -> np.ndarray:
    """Get barcode indices which have nonzero counts or pass the filter."""
    barcode_info = mc.get_barcode_info()
    unique_bc_idx = barcode_info.pass_filter[:, BarcodeInfo.PASS_FILTER_BARCODE_IDX]
    for chunk_start, chunk_len in mc.get_chunks(tgt_chunk_len, preserve_boundaries=False):
        unique_bc_idx = np.union1d(
            unique_bc_idx,
            mc.get_column_lazy(BARCODE_IDX_COL_NAME)[chunk_start : chunk_start + chunk_len],
        )
    return np.array(unique_bc_idx)


## Taken from MIT licensed: https://github.com/hmallen/numpyencoder
class NumpyEncoder(json.JSONEncoder):
    """Custom encoder for numpy data types."""

    def default(self, o):
        """Encode the object.

        Args:
            o: The object to encode.

        Returns:
            A transformed version of the object to encode.
        """
        if isinstance(
            o,
            np.int_
            | np.intc
            | np.intp
            | np.int8
            | np.int16
            | np.int32
            | np.int64
            | np.uint8
            | np.uint16
            | np.uint32
            | np.uint64,
        ):
            return int(o)

        elif isinstance(o, np.float16 | np.float32 | np.float64):
            return float(o)

        elif isinstance(o, np.complex64 | np.complex128):
            return {"real": o.real, "imag": o.imag}

        elif isinstance(o, np.ndarray):
            return o.tolist()

        elif isinstance(o, np.bool_):
            return bool(o)

        elif isinstance(o, np.void):
            return None

        return json.JSONEncoder.default(self, o)


def set_file_version(mc: MoleculeCounter, version: int):
    """Set the file version.

    Args:
        mc: an moleculecounter instance
        version: The version to set
    """
    mc.file_version = version
    cr_h5.set_hdf5_attr(mc.h5, FILE_VERSION_KEY, version)


def get_v2_metrics(h5_file):
    group = tables.open_file(h5_file, "r").get_node("/metrics")
    attrset = group._v_attrs
    return {k: attrset[k] for k in attrset._f_list()}


def get_v2_library_info_and_chunks(v2_mc_in: h5py.File):
    """Method to generate library_info information for molecule info files with versions prior to this field.

     being introduced.

    Args:
        hdf_file: an open V2 or earlier molecule info file

    Returns:
        A tuple containing a list of library infos, as well as a corresponding list of chunks, one for each gem group
        with a tuple defining it as (gem_group, chunk_start, chunk_length, lib_idx)
    """

    def get_chunks_by_gem_group(gem_group_arr: np.array):
        """Return exactly one chunk per gem group."""
        # verify gem groups are sorted
        assert np.all(np.diff(gem_group_arr) >= 0)
        num_rows = gem_group_arr.shape[0]
        unique_ggs = np.unique(gem_group_arr)

        def gg_key(i):
            return gem_group_arr[i]

        chunk_iter = MoleculeCounter.get_chunks_from_partition_static(num_rows, unique_ggs, gg_key)
        chunks = []
        for lib_idx, (gg, chunk) in enumerate(zip(unique_ggs, chunk_iter)):
            chunks.append((gg, chunk[0], chunk[1], lib_idx))
        return chunks

    # Load up the gem_groups in file,
    v2_gem_groups = np.asarray(
        v2_mc_in[GEM_GROUP_COL_NAME], dtype=MOLECULE_INFO_COLUMNS[GEM_GROUP_COL_NAME]
    )
    library_info = []
    chunks = get_chunks_by_gem_group(
        v2_gem_groups
    )  # list of tuples with gem_group, start, length, lib_idx
    for gem_group, _, _, lib_idx in chunks:
        library_info.append(
            {
                GEM_GROUP_COL_NAME: int(gem_group),
                "library_id": str(lib_idx),
                LIBRARY_TYPE: rna_library.GENE_EXPRESSION_LIBRARY_TYPE,
            }
        )
    return library_info, chunks


def _get_library_info(h5f: h5py.File) -> list[dict[str, Any]]:
    """Get the library info.

    Args:
        mc: MoleculeCounter instance

    Returns:
        Dictionary of molecule infos
    """
    return json.loads(h5f["library_info"][0])


def get_library_info(mol_info_fname) -> list[dict[str, Any]]:
    """Takes a molecule info filename and loads the library info  This method.

    allows one to either read the `library_info` from new files, or generate one on the fly from
    older molecule info files without running into version compatibility errors
    Args:
        mol_info_fname: file name of a molecule info file

    Returns:
        Dictionaries of library infos
    """
    mc, version = get_h5py_file_and_version(mol_info_fname)
    if version < 3:
        return get_v2_library_info_and_chunks(mc)[0]
    else:
        return _get_library_info(mc)


def get_h5py_file_and_version(mol_info_fname, mode="r") -> tuple[h5py.File, int]:
    """Opens a molecule info h5py.File.

    Args:
        mol_info_fname: Path to a molecule info file to open

    Returns:
        A tuple with objects (h5py.File, file_version_number)
    """
    try:
        mc_h5 = h5py.File(mol_info_fname, mode=mode)
    except OSError as ex:
        raise OSError(
            f"The molecule info HDF5 file ({mol_info_fname}) is invalid. Please provide a valid HDF5 file."
        ) from ex
    if not is_valid_mol_info_h5(mc_h5):
        raise ValueError(
            f"The input molecule info HDF5 file ({mol_info_fname}) does not appear to be properly formatted. Please provide a valid file."
        )
    elif FILE_VERSION_KEY in mc_h5.attrs:
        file_version = int(mc_h5.attrs[FILE_VERSION_KEY])
    else:
        file_version = 1  # V1 doesn't have version field
    return (mc_h5, file_version)


def is_valid_mol_info_h5(mol_info):
    """Checks whether a file is a valid molecule info file.

    The molecule info HDF5 file is deemed valid if and only if all of its attribute keys exist in VALID_MOL_INFO_H5_ATTRIBUTE_KEYS

    Args:
        mol_info: Open h5py.File for the molecule info HDF5 file to be validated

    Returns:
        A boolean with True for valid files, False otherwise
    """
    return all(a in VALID_MOL_INFO_H5_ATTRIBUTE_KEYS for a in mol_info.attrs)


def barcode_whitelist_from_metrics(metrics):
    old_wl_metric_key = "chemistry_barcode_whitelist"
    if old_wl_metric_key in metrics:
        return metrics[old_wl_metric_key]
    elif isinstance(metrics["chemistry_barcode"], dict):
        return metrics["chemistry_barcode"]["whitelist"]
    else:
        barcode_def = metrics["chemistry_barcode"][0]
        GB_BARCODE_KIND = "gel_bead"
        if barcode_def["kind"] == GB_BARCODE_KIND:
            whitelist = barcode_def["whitelist"]
            if isinstance(whitelist, dict):
                return whitelist["name"]
            else:
                return whitelist
        else:
            return None


class MoleculeCounter:
    """Streams a list of tuples w/named elements to or from an h5 file."""

    # pylint: disable=too-many-public-methods

    def __init__(self):
        self.file_version = None
        self.h5: h5py.File | None = None
        self.columns: dict[str, h5py.Dataset] = OrderedDict()
        self.ref_columns: dict[str, h5py.Dataset] = OrderedDict()
        self.library_info: list[dict[str, Any]] | None = None
        self.feature_reference: FeatureReference | None = None

    def get_gb_barcode_whitelist(self):
        return barcode_whitelist_from_metrics(self.get_all_metrics())

    def get_visium_hd_slide_name(self):
        metrics = self.get_all_metrics()
        barcode_def = metrics.get("chemistry_barcode", None)
        if barcode_def and isinstance(barcode_def, list):
            whitelist = barcode_def[0]["whitelist"]
            if isinstance(whitelist, dict):
                return whitelist.get("slide", None)
        return None

    def get_gem_groups(self) -> list[int]:
        return [int(x) for x in self.get_metric(GEM_GROUPS_METRIC).keys()]

    def get_genomes(self) -> list[str]:
        return self.feature_reference.get_genomes(feature_type=rna_library.DEFAULT_LIBRARY_TYPE)

    def get_molecule_info_type(self) -> str:
        return self.get_all_metrics().get(MOLECULE_INFO_TYPE_METRIC, MOLECULE_INFO_TYPE_COUNT)

    def get_barcode_list_size(self) -> int:
        return self.get_ref_column_lazy(BARCODES).shape[0]

    def is_aggregated(self) -> bool:
        ret = self.get_metric(IS_AGGREGATED_METRIC)
        return ret if ret is not None else False

    @staticmethod
    def get_column_dtype(k) -> np.dtype:
        return np.dtype(MOLECULE_INFO_COLUMNS[k])

    @staticmethod
    def get_record_bytes() -> int:
        return sum(np.dtype(x).itemsize for x in MOLECULE_INFO_COLUMNS.values())

    @staticmethod
    def estimate_mem_gb(chunk_len: int, scale: float = 1.0, cap: bool = True) -> int:
        """Estimate memory usage in GB (not GiB) of this object given a number of records."""
        mol_entries_per_gb = int(1e9 / MoleculeCounter.get_record_bytes())
        mem_gb = math.ceil(scale * chunk_len / mol_entries_per_gb)
        if cap:
            return max(h5_constants.MIN_MEM_GB, mem_gb)
        else:
            return mem_gb

    @staticmethod
    def build_barcode_info(
        filtered_barcodes_by_genome: dict[bytes, list[bytes]], library_info, barcodes: list[bytes]
    ) -> BarcodeInfo:
        """Generate numpy arrays for per-barcode info.

        Args:
          filtered_barcodes_by_genome (dict[str,list[str]]): Keys are genomes, values are lists of filtered barcode strings.
          library_info (list[dict]): Per-library metadata.
          barcodes (list[bytes]): All barcode sequences (e.g. ['ACGT', ...]

        Returns:
          BarcodeInfo: object
        """
        # Replace a genome string with its lexicographical rank
        genome_to_idx = {g: i for i, g in enumerate(sorted(filtered_barcodes_by_genome.keys()))}

        libraries_for_gem_group = defaultdict(list)
        for lib_idx, lib in enumerate(library_info):
            libraries_for_gem_group[lib["gem_group"]].append(lib_idx)

        # Map a barcode sequence to its index into the MoleculeCounter
        #  'barcodes' array
        bc_seq_to_idx = {bc: i for i, bc in enumerate(barcodes)}

        # Populate the "pass filter" array of tuples
        pf_tuples = []
        for genome, bcs in filtered_barcodes_by_genome.items():
            genome_idx = genome_to_idx[genome]
            for bc_str in bcs:
                seq, gg = cr_utils.split_barcode_seq(bc_str)
                barcode_idx = bc_seq_to_idx[seq]

                library_inds = libraries_for_gem_group[gg]
                for library_idx in library_inds:
                    pf_tuples.append((barcode_idx, library_idx, genome_idx))

        if len(pf_tuples) > 0:
            pass_filter = np.array(pf_tuples, dtype=BARCODE_INFO_DTYPES["pass_filter"])
        else:
            pass_filter = np.zeros((0, 3), dtype=BARCODE_INFO_DTYPES["pass_filter"])

        assert pass_filter.shape[0] == len(pf_tuples)
        assert pass_filter.shape[1] == 3

        # Sort by barcode index
        pass_filter = pass_filter[np.argsort(pass_filter[:, 0]), :]

        return BarcodeInfo(
            pass_filter,
            genomes=sorted(filtered_barcodes_by_genome.keys()),
        )

    @staticmethod
    def get_filtered_barcodes(
        barcode_info: BarcodeInfo,
        library_info,
        barcodes,
        genome_idx: int | None = None,
        library_type=None,
    ):
        """Get a list of filtered barcode strings e.g. ['ACGT-1',...].

        Args:
          barcode_info (BarcodeInfo): Barcode info object.
          library_info (list of dict): Library info.
          barcodes (np.array): Barcode sequences.
          genome_idx (int): Restrict passing definition to this genome. None for no restriction.
          library_type (str): Restrict passing definition to this library type. None for no restriction.

        Returns:
          list of str
        """
        # Without restrictions, assumes passing filter in a single library or genome is sufficient
        # for a barcode to be passing filter overall.

        pass_filter = barcode_info.pass_filter

        pf_barcode_idx = pass_filter[:, BarcodeInfo.PASS_FILTER_BARCODE_IDX]
        pf_library_idx = pass_filter[:, BarcodeInfo.PASS_FILTER_LIBRARY_IDX]
        pf_genome_idx = pass_filter[:, BarcodeInfo.PASS_FILTER_GENOME_IDX]

        mask = np.ones(pass_filter.shape[0], dtype=bool)
        if genome_idx is not None:
            mask &= pf_genome_idx == genome_idx

        if library_type is not None:
            library_inds = np.array(
                [i for i, lib in enumerate(library_info) if lib[LIBRARY_TYPE] == library_type],
                dtype=MOLECULE_INFO_COLUMNS["library_idx"],
            )
            mask &= np.isin(pf_library_idx, library_inds)
        inds = np.flatnonzero(mask)

        lib_to_gg = np.array([lib["gem_group"] for lib in library_info], dtype="uint64")

        pf_gem_group = lib_to_gg[pf_library_idx[inds]]

        # Take unique, sorted barcodes (sorted by (gem_group, barcode_idx))
        gg_bcs = np.unique(np.column_stack((pf_gem_group, pf_barcode_idx[inds])), axis=0)

        # Create barcode strings
        return [
            cr_utils.format_barcode_seq(barcodes[gg_bcs[i, 1]], gg_bcs[i, 0])
            for i in range(gg_bcs.shape[0])
        ]

    @staticmethod
    def save_barcode_info(bc_info: BarcodeInfo, group: h5py.Group):
        """Save barcode info to HDF5.

        Args:
          barcode_info (BarcodeInfo): Data.
          group (h5py.Group): Output group.
        """
        group.create_dataset(
            "pass_filter",
            data=bc_info.pass_filter,
            maxshape=(None, bc_info.pass_filter.shape[1]),
            compression=HDF5_COMPRESSION,
            shuffle=True,
        )
        cr_h5.create_hdf5_string_dataset(
            group, "genomes", bc_info.genomes, compression=HDF5_COMPRESSION, shuffle=True
        )

    @staticmethod
    def load_barcode_info(group: h5py.Group) -> BarcodeInfo:
        """Load barcode info from an HDF5 group.

        Args:
          group (h5py.Group): Input group.

        Returns:
          BarcodeInfo: object
        """
        return BarcodeInfo(
            pass_filter=group["pass_filter"][:],
            genomes=cr_h5.read_hdf5_string_dataset(group["genomes"]),
        )

    def get_barcode_info(self) -> BarcodeInfo:
        return MoleculeCounter.load_barcode_info(self.h5[BARCODE_INFO_GROUP_NAME])

    @classmethod
    def open(
        cls,
        filename,
        mode: str,
        feature_ref: FeatureReference | None = None,
        barcodes=None,
        library_info=None,
        barcode_info: BarcodeInfo | None = None,
    ) -> MoleculeCounter:
        """Open a molecule info object.

        Args:
          filename (str): Filename to open or create
          mode (str): 'r' for reading, 'w' for writing.
          feature_ref (FeatureReference): Required when mode is 'w'.
          barcodes (list of str): All possible barcode sequences. Required when mode is 'w'.
          library_info (list of dict): Library metadata. Required when mode is 'w'.
          barcode_info (BarcodeInfo): Per-barcode metadata.

        Returns:
          MoleculeInfo: A new object
        """
        assert mode in ["r", "r+", "w"]

        mc = cls()

        if mode == "w":
            if feature_ref is None:
                raise ValueError(
                    "Feature reference must be specified when opening a molecule info object for writing"
                )
            if barcodes is None:
                raise ValueError(
                    "Barcodes must be specified when opening a molecule info object for writing"
                )
            if library_info is None:
                raise ValueError(
                    "Library info must be specified when opening a molecule info object for writing"
                )
            if barcode_info is None:
                raise ValueError(
                    "Barcode info must be specified when opening a molecule info object for writing"
                )

            mc.h5 = h5py.File(filename, "w")
            set_file_version(mc, CURR_FILE_VERSION)
            cr_h5.set_hdf5_attr(mc.h5, h5_constants.H5_FILETYPE_KEY, MOLECULE_H5_FILETYPE)

            # Write feature reference
            fref_group = mc.h5.create_group(h5_constants.H5_FEATURE_REF_ATTR)
            feature_ref.to_hdf5(fref_group)

            # Write barcodes
            _write_barcodes(mc, barcodes)

            # Write library info
            lib_info_json = json.dumps(library_info, indent=4, sort_keys=True)
            cr_h5.create_hdf5_string_dataset(mc.h5, "library_info", [lib_info_json])

            # Write barcode info
            g = mc.h5.create_group(BARCODE_INFO_GROUP_NAME)
            MoleculeCounter.save_barcode_info(barcode_info, g)

            # Create empty per-molecule datasets
            for name in MOLECULE_INFO_COLUMNS:
                create_dataset(mc, name)

        else:  # r or r+
            mc.h5, mc.file_version = get_h5py_file_and_version(filename, mode=mode)
            assert isinstance(mc.file_version, int)  # Help pylint's type inference.
            if mc.file_version > CURR_FILE_VERSION:
                raise ValueError(
                    "The molecule info HDF5 file (file: "
                    f"{filename}, format version {mc.file_version}) was produced "
                    "by a newer software version. Reading these files is unsupported."
                )
            if mc.file_version < 3:
                raise ValueError(
                    "The molecule info HDF5 file (file: "
                    f"{filename}, format version {mc.file_version}) was produced "
                    "by an older software version. Reading these files is unsupported."
                )
            # h5 is not actually a dict.
            for key in mc.h5.keys():  # pylint: disable=consider-using-dict-items
                if key in MOLECULE_INFO_COLUMNS:
                    mc.columns[key] = mc.h5[key]
                elif key in MOLECULE_REF_COLUMNS:
                    mc.ref_columns[key] = mc.h5[key]
                elif key == h5_constants.H5_FEATURE_REF_ATTR:
                    mc.feature_reference = FeatureReference.from_hdf5(mc.h5[key])
                elif (
                    key == V3_METRICS_GROUP_NAME
                    or key == BARCODE_INFO_GROUP_NAME
                    or key == METRICS_JSON_DATASET_NAME
                    or key in UNIMPLEMENTED_PROBE_KEYS
                ):
                    pass
                else:
                    raise AttributeError(f"Unrecognized dataset key: {key}")

            # Load library info
            mc.library_info = json.loads(cr_h5.read_hdf5_string_dataset(mc.h5["library_info"])[0])
        return mc

    def is_spatial_data(self) -> bool:
        """Using the chemistry whitelist, determine if this is Spatial (Visium) data.

        Returns:
            bool: True if the known barcode list is for Spatial data.
        """
        return (
            bc_utils.is_whitelist_spatial(self.get_gb_barcode_whitelist())
            or self.get_visium_hd_slide_name() is not None
        )

    def nrows(self) -> int:
        return self.get_column_lazy(next(iter(MOLECULE_INFO_COLUMNS))).shape[0]

    def get_chunk_key(self, idx: int) -> tuple[Any, Any]:
        return tuple(self.get_column_lazy(col)[idx] for col in CHUNK_COLUMNS)

    def set_metric(self, key, value):
        """Set a metric.

        Serialize to Pickle.
        """
        metrics = self.get_all_metrics()
        metrics[key] = value
        self.set_all_metrics(metrics)

    def get_metric(self, key: str, default: Any = None):
        """Get a metric."""
        return self.get_all_metrics().get(key, default)

    def set_all_metrics(self, metrics):
        metrics_json = json.dumps(metrics, sort_keys=True, separators=(",", ":"), cls=NumpyEncoder)
        self.h5[METRICS_JSON_DATASET_NAME] = metrics_json

    def get_all_metrics(self) -> dict[str, Any]:
        """Return a dictionary of metrics."""
        if METRICS_JSON_DATASET_NAME in self.h5:
            return json.loads(self.h5[METRICS_JSON_DATASET_NAME][()])
        elif V3_METRICS_GROUP_NAME in self.h5:
            return self._legacy_get_all_metrics()
        else:
            return {}

    def _legacy_get_all_metrics(self: MoleculeCounter) -> dict:
        return {
            k: pickle.loads(ensure_binary(v))
            for k, v in self.h5[V3_METRICS_GROUP_NAME].attrs.items()
        }

    def append_column(self, name: str, values):
        """Append an array of values to a column."""
        ds = self.columns[name]
        start = len(ds)
        end = start + len(values)
        ds.resize((end,))
        ds[start:end] = values
        self.h5.flush()

    def get_column_lazy(self, col_name: str) -> h5py.Dataset:
        """Retrieve column.

        Does not handle missing columns.

        Depending on how the file was opened,
        this may only be a file view instead of a full array.
        """
        return self.columns[col_name]

    def get_column(self, col_name: str) -> np.ndarray:
        """Load an entire column of data into memory.

        Handles missing columns.
        """
        return (
            self.get_column_lazy(col_name)[:]
            if col_name in self.columns
            else np.full(
                self.nrows(),
                MOLECULE_INFO_DEFAULT_VALUES[col_name],
                dtype=MOLECULE_INFO_COLUMNS[col_name],
            )
        )

    def set_ref_column(self, col_name: str, values):
        assert col_name in MOLECULE_REF_COLUMNS
        # pylint: disable=no-member
        self.ref_columns[col_name] = self.h5.create_carray(
            self.h5.root, col_name, obj=np.array(values)
        )

    def get_ref_column(self, col_name: str) -> np.ndarray:
        """Load a reference array into memory as a numpy array."""
        return self.get_ref_column_lazy(col_name)[:]

    def get_ref_column_lazy(self, col_name: str) -> h5py.Dataset:
        """Get a reference array as a lazy h5py Dataset."""
        return self.ref_columns[col_name]

    def get_feature_ref(self) -> FeatureReference:
        """Returns the feature reference from HDF5."""
        return FeatureReference.from_hdf5(self.h5[h5_constants.H5_FEATURE_REF_ATTR])

    def get_barcodes(self) -> np.ndarray:
        """Returns fully loaded barcode set."""
        return self.h5[BARCODE_DS_NAME][:]

    def trim_barcodes(self):
        """Changes the data in the dataset so that only barcodes with have > 0 counts and/or were pass filtered are.

        stored in the file.  Updates the barcode_idx file to indicate the new indices
        Returns:
        """
        old_barcodes = self.get_barcodes()
        old_barcode_info = self.get_barcode_info()
        new_pass_filter = old_barcode_info.pass_filter
        bc_idx_to_retain = get_barcode_index_to_retain(self)
        new_barcodes = old_barcodes[bc_idx_to_retain]
        # Make an array with indices of new barcode
        col_dtype = MOLECULE_INFO_COLUMNS[BARCODE_IDX_COL_NAME]
        new_positions = np.zeros(len(old_barcodes), dtype=col_dtype)
        for n, o in enumerate(bc_idx_to_retain):
            new_positions[o] = n

        for chnk_start, chnk_len in self.get_chunks(2000000, preserve_boundaries=False):
            future = self.h5[BARCODE_IDX_COL_NAME][chnk_start : (chnk_start + chnk_len)]
            self.h5[BARCODE_IDX_COL_NAME][chnk_start : (chnk_start + chnk_len)] = new_positions[
                future
            ]

        new_pass_filter[:, BarcodeInfo.PASS_FILTER_BARCODE_IDX] = [
            new_positions[x] for x in new_pass_filter[:, BarcodeInfo.PASS_FILTER_BARCODE_IDX]
        ]
        new_barcode_info = BarcodeInfo(
            pass_filter=new_pass_filter,
            genomes=old_barcode_info.genomes,
        )
        # write barcodes
        del self.h5[BARCODE_DS_NAME]
        _write_barcodes(self, new_barcodes)
        # Write barcode info
        del self.h5[BARCODE_INFO_GROUP_NAME]
        g = self.h5.create_group(BARCODE_INFO_GROUP_NAME)
        MoleculeCounter.save_barcode_info(new_barcode_info, g)
        return new_barcodes

    def get_column_with_indices(
        self, col_name, idxs: np.ndarray, chunk_size: int = (1 << 20)
    ) -> np.ndarray:
        """Get the column values at row indices (boolean array)."""
        col = self.get_column_lazy(col_name)
        nmol = self.nrows()
        if idxs.shape == (nmol,) and idxs.dtype == bool:
            # index boolean array
            vals = np.empty((idxs.sum(),), dtype=col.dtype)
            lwr = 0
            for begin in range(0, nmol, chunk_size):
                end = min(begin + chunk_size, nmol)
                upr = lwr + idxs[begin:end].sum()
                if lwr == upr:
                    continue
                vals[lwr:upr] = col[begin:end][idxs[begin:end]]
                lwr = upr
        elif idxs.shape[0] == 0:
            return np.empty([], dtype=col.dtype)
        else:
            # index value array
            raise NotImplementedError("index value arrays not yet supported")
        return vals

    def get_num_filtered_barcodes_for_library(self, library_idx: int) -> int:
        """Count the number of barcodes passing filter for a library.

        Args:
          library_idx (int): Index of library to count.

        Returns:
          int: Number of filtered barcodes for this library.
        """
        pass_filter = self.h5[BARCODE_INFO_GROUP_NAME]["pass_filter"][:]
        this_lib = np.flatnonzero(pass_filter[:, 1] == library_idx)
        barcode_inds = pass_filter[this_lib, 0]
        return len(np.unique(barcode_inds))

    def get_num_filtered_barcodes_for_libraries(self, library_idxs: np.ndarray) -> np.ndarray:
        """Count the number of barcodes passing filter for a range of libraries.

        This is the vectorized equivalent to

        .. code-block:: python

            np.array(
                [
                    mc.get_num_filtered_barcodes_for_library(lib_idx)
                    for lib_idx in library_idxs
                ]
            )

        Args:
          library_idxs (np.ndarray[int]): 1d array of indicies to count.

        Returns:
          np.ndarray[int]: Number of filtered barcodes for each library.
        """
        pass_filter = self.h5[BARCODE_INFO_GROUP_NAME]["pass_filter"][:]
        # pass_filter is Nx3 array of integers
        assert library_idxs.ndim == 1
        return np.count_nonzero(pass_filter[:, 1:2] == library_idxs, axis=0)

    def get_num_filtered_barcodes(self) -> int:
        """Return the number of filtered barcodes."""
        assert self.h5 is not None
        return np.unique(self.h5[BARCODE_INFO_GROUP_NAME]["pass_filter"][:, 0]).size

    def get_library_info(self) -> list[dict[str, Any]]:
        return _get_library_info(self.h5)

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        # pylint: disable=redefined-builtin
        self.close()

    def close(self):
        self.h5.close()

    def save(self):
        self.h5.close()

    @staticmethod
    def merge_barcode_infos(bc_infos: list[BarcodeInfo]) -> BarcodeInfo:
        """Merge a BarcodeInfo into another BarcodeInfo.

        Args:
          bc_infos (list of BarcodeInfo): Input BarcodeInfos.

        Returns:
          BarcodeInfo:
        """
        assert len(bc_infos) > 0
        genomes = bc_infos[0].genomes

        # Total number of barcodes with any information
        pfs = []
        for bc_info in bc_infos:
            assert bc_info.pass_filter.shape[1] == 3
            assert bc_info.genomes == genomes
            pfs.append(bc_info.pass_filter)

        new_pf = np.concatenate(pfs, axis=0)

        # Deduplicate the tuples. Unique throws an error on a zero-row array.
        if new_pf.shape[0] > 0:
            new_pf = np.unique(new_pf, axis=0)

        return BarcodeInfo(
            pass_filter=new_pf,
            genomes=genomes,
        )

    @staticmethod
    def concatenate(out_filename, in_filenames, metrics: dict | None = None) -> None:
        """Concatenate MoleculeCounter HDF5 files.

        Args:
          out_filename (str): Output HDF5 filename
          in_filenames (list of str): Input HDF5 filenames
          metrics (dict): Metrics to write
        """
        # Load reference info from first file
        first_mc = MoleculeCounter.open(in_filenames[0], "r")
        feature_ref = first_mc.get_feature_ref()
        barcodes = first_mc.get_barcodes()
        library_info = first_mc.get_library_info()

        feature_ids = [f.id for f in feature_ref.feature_defs]

        # print 'Merging barcode info'
        bc_infos = []
        for filename in in_filenames:
            with MoleculeCounter.open(filename, "r") as mc:
                bc_infos.append(mc.get_barcode_info())
        merged_bc_info = MoleculeCounter.merge_barcode_infos(bc_infos)

        # print 'Concatenating molecule info files'
        out_mc = MoleculeCounter.open(
            out_filename,
            mode="w",
            feature_ref=feature_ref,
            barcodes=barcodes,
            library_info=library_info,
            barcode_info=merged_bc_info,
        )

        # TODO: This inefficient code block is to check assumption which should always be true.
        total_rows = 0
        for filename in in_filenames:
            with MoleculeCounter.open(filename, mode="r") as in_mc:
                # Assert that these data are compatible
                assert in_mc.get_library_info() == library_info
                assert np.array_equal(in_mc.get_barcodes(), barcodes)
                fref = in_mc.get_feature_ref()
                assert [f.id for f in fref.feature_defs] == feature_ids

                # if no metrics specified, copy them from the first file
                if metrics is None:
                    metrics = in_mc.get_all_metrics()
                total_rows += in_mc.nrows()

        out_mc.set_all_metrics(metrics)
        out_mc.save()

        # Now use Rust to concatenate all the columns
        concatenate_molecule_infos(out_filename, in_filenames)

        # Slow validation check here
        with MoleculeCounter.open(out_filename, "r") as mc:
            assert total_rows == mc.nrows(), "Concatenation did not produce expected results."

    def find_last_occurrence_of_chunk_key(self, from_row: int) -> int:
        num_rows = self.nrows()
        initial_chunk_key = self.get_chunk_key(from_row)
        for i in range(from_row, num_rows):
            chunk_key = self.get_chunk_key(i)
            if chunk_key != initial_chunk_key:
                return i - 1
        return num_rows - 1

    def bisect(self, query: _T1, key_func: Callable[[int], _T1]) -> int:
        return MoleculeCounter.bisect_static(self.nrows(), query, key_func)

    @staticmethod
    def bisect_static(num_rows: int, query: _T1, key_func: Callable[[int], _T1]) -> int:
        """Performs a binary search to find the leftmost insertion point of query.

        Args:
            key_func: A function, where `key_func(i)` is the value to compare to at index i.
        """
        lo = 0
        hi = num_rows
        exists = True
        while True:
            i = (hi + lo) // 2
            curr = key_func(i)
            if curr == query:
                break
            elif hi - lo <= 1:
                # non-matching case
                exists = False
                break
            elif curr < query:
                lo = i
            else:
                hi = i

        if exists:
            # backtrack to first occurrence
            for j in range(i, -1, -1):
                curr = key_func(j)
                if curr != query:
                    return j + 1
        return 0

    def get_chunks_from_partition(
        self, values: list[_T1], key_func: Callable[[int], _T1]
    ) -> Generator[tuple[int, int], None, None]:
        return MoleculeCounter.get_chunks_from_partition_static(self.nrows(), values, key_func)

    @staticmethod
    def get_chunks_from_partition_static(
        num_rows: int, values: list[_T1], key_func: Callable[[int], _T1]
    ) -> Generator[tuple[int, int], None, None]:
        """Get chunks by partitioning on the specified values."""
        starts = [0] + [
            MoleculeCounter.bisect_static(num_rows, val, key_func) for val in values[1:]
        ]
        n = len(starts)
        for i in range(n):
            chunk_start = starts[i]
            chunk_end = starts[i + 1] if i + 1 < n else num_rows
            yield (chunk_start, chunk_end - chunk_start)

    def get_chunks(
        self, target_chunk_len: int, preserve_boundaries: bool = True
    ) -> Generator[tuple[int, int], None, None]:
        """Get chunks, optionally preserving boundaries defined by get_chunk_key().

        Yields:
            (int, int): (chunk_start, chunk_len) which are closed intervals
        """
        num_rows = self.nrows()
        chunk_start, chunk_end = 0, 0
        while chunk_end < (num_rows - 1):
            target_chunk_end = min(num_rows - 1, chunk_start + target_chunk_len - 1)
            chunk_end = (
                self.find_last_occurrence_of_chunk_key(target_chunk_end)
                if preserve_boundaries
                else target_chunk_end
            )
            chunk_len = 1 + chunk_end - chunk_start
            yield (chunk_start, chunk_len)
            chunk_start = 1 + chunk_end

    @staticmethod
    def compress_gem_group(x) -> np.uint16:
        return MOLECULE_INFO_COLUMNS[GEM_GROUP_COL_NAME](x)

    @staticmethod
    def compress_umi_seq(x: bytes, umi_bits: int) -> int:
        return cr_utils.compress_seq(x, umi_bits)

    @staticmethod
    def naive_concatenate_metric_list(metrics_list: list[dict]):
        combined_metrics: dict | None = None
        gg_metrics = {}
        lib_metrics = {}
        targeted_metrics = []
        for single_metrics in metrics_list:
            if combined_metrics is None:
                combined_metrics = {
                    k: v for k, v in single_metrics.items() if not k.startswith("target")
                }

            # Collect the targeted metrics for each input.
            targeted_metrics.append(
                {k: v for k, v in single_metrics.items() if k.startswith("target")}
            )

            # concatenate new gem groups to the metrics. if it collides with an existing
            # gem group, the old one will be overwritten.
            single_gg_metrics = single_metrics[GEM_GROUPS_METRIC]
            analysis_parameters = single_metrics.get(ANALYSIS_PARAMETERS_METRIC)
            if analysis_parameters is not None:
                for k in single_gg_metrics:
                    single_gg_metrics[k].update(analysis_parameters)
            gg_metrics.update(single_gg_metrics)
            lib_metrics.update(single_metrics[LIBRARIES_METRIC])

        assert combined_metrics is not None
        combined_metrics[GEM_GROUPS_METRIC] = gg_metrics
        combined_metrics[LIBRARIES_METRIC] = lib_metrics
        combined_metrics.pop(ANALYSIS_PARAMETERS_METRIC, None)

        # Pass through the targeting metrics if all inputs are identical.
        if all(x == targeted_metrics[0] for x in targeted_metrics):
            combined_metrics.update(targeted_metrics[0])

        return combined_metrics

    @staticmethod
    def naive_concatenate_metrics_from_h5s(mol_h5_list: list[str]):
        metrics = []
        for mol_h5 in mol_h5_list:
            with MoleculeCounter.open(mol_h5, mode="r") as counter:
                single_metrics = counter.get_all_metrics()
                metrics.append(single_metrics)
        combined_metrics = MoleculeCounter.naive_concatenate_metric_list(metrics)
        return combined_metrics

    def get_raw_read_pairs_per_library(self, all_metrics: dict | None = None) -> list[int]:
        """Get raw read pairs per library.

        Returns:
          list of int: Order is by library index

        Raises:
            ValueError if expected metrics are missing.
        """
        try:
            assert self.library_info is not None
            metric = (
                self.get_metric(LIBRARIES_METRIC)
                if all_metrics is None
                else all_metrics.get(LIBRARIES_METRIC)
            )
            return [metric[str(li)][TOTAL_READS_METRIC] for li in range(len(self.library_info))]
        except KeyError as exc:
            raise ValueError("Missing metrics in molecule counter") from exc

    def get_read_pairs_in_filtered_barcodes_per_library(
        self, all_metrics: dict | None = None
    ) -> list[int]:
        """Get read pairs in filtered barcodes per library.

        Returns:
          list of int: Order is by library index

        Raises:
            ValueError if expected metrics are missing.
        """
        try:
            assert self.library_info is not None
            metric = (
                self.get_metric(LIBRARIES_METRIC)
                if all_metrics is None
                else all_metrics.get(LIBRARIES_METRIC)
            )
            return [
                metric[str(li)][TOTAL_READS_IN_FILTERED_BARCODES_METRIC]
                for li in range(len(self.library_info))
            ]
        except KeyError:
            raise ValueError("Missing metrics in molecule counter")

    def get_transcriptomic_read_pairs_per_library(
        self, all_metrics: dict | None = None
    ) -> list[int]:
        """Get transcriptomic read pairs per library (cell + non-cell).

        Returns:
            list of int: Order is by library index

        Raises:
            ValueError if expected metrics are missing.
        """
        try:
            assert self.library_info is not None
            metric = (
                self.get_metric(LIBRARIES_METRIC)
                if all_metrics is None
                else all_metrics.get(LIBRARIES_METRIC)
            )
            return [
                metric[str(li)][FEATURE_READS_METRIC]
                for li in range(len(self.library_info))
                # TODO: should this restrict itself to GEX libraries only?
            ]
        except KeyError as exc:
            raise ValueError("Missing metrics in molecule counter") from exc

    def get_usable_read_pairs_per_library(self) -> list[int]:
        """Get usable read pairs per library.

        Returns:
          list of int: Order is by library index

        Raises:
            ValueError if expected metrics are missing.
        """
        try:
            assert self.library_info is not None
            metric = self.get_metric(LIBRARIES_METRIC)
            return [metric[str(li)][USABLE_READS_METRIC] for li in range(len(self.library_info))]
        except KeyError as exc:
            raise ValueError("Missing metrics in molecule counter") from exc

    def get_on_target_usable_read_pairs_per_library(self) -> list[int]:
        """Get usable on-target read pairs per library.

        Returns:
          list of int: Order is by library index

        Raises:
            ValueError if expected metrics are missing.
        """
        try:
            assert self.library_info is not None
            metric = self.get_metric(LIBRARIES_METRIC)
            return [
                metric[str(li)][ON_TARGET_USABLE_READS_METRIC]
                for li in range(len(self.library_info))
            ]
        except KeyError as exc:
            raise ValueError("Missing metrics in molecule counter") from exc

    def get_library_indices_by_type(self) -> dict[str, list[int]]:
        """Get indices of libraries that are GEX libraries."""
        libraries = self.get_library_info()
        library_indices_by_type = defaultdict(list)
        for lib_idx, lib in enumerate(libraries):
            library_indices_by_type[lib[LIBRARY_TYPE]].append(lib_idx)
        return library_indices_by_type

    def is_templated_ligation(self) -> bool:
        """Return True if metrics_json/targeting_method is templated_ligation."""
        return self.get_metric(TARGETING_METHOD_METRIC) == TARGETING_METHOD_TL

    @staticmethod
    def _sum_metric(mol_h5_list, metric_name, metric_type: str):
        """Combine a library- or gemgroup- level integer metric across multiple h5 files."""
        assert metric_type is LIBRARIES_METRIC or metric_type is GEM_GROUPS_METRIC
        combined = defaultdict(int)
        for mol_h5 in mol_h5_list:
            with MoleculeCounter.open(mol_h5, mode="r") as counter:
                for key, metrics in counter.get_metric(metric_type).items():
                    combined[key] += metrics[metric_name]
        return combined

    @staticmethod
    def sum_library_metric(mol_h5_list, metric_name):
        return MoleculeCounter._sum_metric(mol_h5_list, metric_name, LIBRARIES_METRIC)

    @staticmethod
    def get_total_conf_mapped_reads_in_filtered_barcodes_chunk(
        filename, filtered_bcs_set, start, length, queue
    ) -> None:
        total_mapped_reads = 0
        with MoleculeCounter.open(filename, "r", start, length) as mc:
            for barcode, gem_group, reads in zip(
                mc.get_column("barcode"), mc.get_column("gem_group"), mc.get_column("reads")
            ):
                if reads < 1:
                    continue
                if (barcode, gem_group) not in filtered_bcs_set:
                    continue
                total_mapped_reads += reads
        queue.put(total_mapped_reads)

    @staticmethod
    def is_targeted_library(library) -> bool:
        return library[
            rna_library.LIBRARY_TYPE
        ] == rna_library.GENE_EXPRESSION_LIBRARY_TYPE and rna_library.has_target_set(library)


class MergedBarcodes(list):
    """Class to hold a list of barcodes."""

    def write_to_disk(self, filename):
        """Write to disk."""
        max_barcode_len = max(len(x) for x in self)
        barcode_dtype = np.dtype("S%d" % max_barcode_len)
        with h5py.File(filename, "w") as h5:
            h5.create_dataset(
                BARCODE_DS_NAME,
                data=np.fromiter(self, barcode_dtype, count=len(self)),
                compression="gzip",
            )
            h5.close()

    @staticmethod
    def load_from_disk(filename):
        """Deserialize the data.

        :param filename:
        :return:
        """
        with h5py.File(filename, "r") as h5:
            result = MergedBarcodes()
            result.extend(h5[BARCODE_DS_NAME][:])
            h5.close()
            return result
