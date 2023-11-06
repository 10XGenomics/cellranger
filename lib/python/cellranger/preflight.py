#!/usr/bin/env python
#
# Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
#
from __future__ import annotations

import os
import socket
import sys
from typing import Any, NamedTuple

import martian

import cellranger.constants as cr_constants
import cellranger.env as cr_env
import cellranger.reference as cr_reference
import cellranger.rna.feature_ref as rna_feature_ref
import cellranger.rna.library as rna_library
import cellranger.sample_def as cr_sample_def
import cellranger.targeted.simple_utils as tgt_simple_utils
import tenkit.log_subprocess as tk_subproc
import tenkit.preflight as tk_preflight
from cellranger import csv_utils
from cellranger.fast_utils import validate_reference
from cellranger.feature_ref import FeatureDefException
from cellranger.molecule_counter import (
    LIBRARIES_METRIC,
    MOLECULE_INFO_TYPE_PERSAMPLE,
    MoleculeCounter,
)
from cellranger.reference_paths import get_reference_genomes
from cellranger.targeted.targeted_constants import (
    TARGETING_METHOD_FILE_NAMES,
    TARGETING_METHOD_TL_FILE_FORMAT,
    TargetingMethod,
)
from cellranger.targeted.targeted_spatial import SPATIAL_TARGET_DISALLOWED_PANEL_TYPES

ALLOWED_LIBRARY_TYPES = [
    rna_library.GENE_EXPRESSION_LIBRARY_TYPE,
    rna_library.CRISPR_LIBRARY_TYPE,
    rna_library.ANTIBODY_LIBRARY_TYPE,
    rna_library.ANTIGEN_LIBRARY_TYPE,
    rna_library.MULTIPLEXING_LIBRARY_TYPE,
    rna_library.CUSTOM_LIBRARY_TYPE,
]

PUBLIC_LIBRARY_TYPES = [
    rna_library.GENE_EXPRESSION_LIBRARY_TYPE,
    rna_library.CRISPR_LIBRARY_TYPE,
    rna_library.ANTIBODY_LIBRARY_TYPE,
    rna_library.CUSTOM_LIBRARY_TYPE,
]

SPATIAL_ALLOWED_LIBRARY_TYPES = [
    rna_library.GENE_EXPRESSION_LIBRARY_TYPE,
    rna_library.ANTIBODY_LIBRARY_TYPE,
]


class PreflightException(Exception):
    def __init__(self, msg):
        super().__init__(msg)
        self.msg = msg

    def __str__(self):
        return self.msg


def check(result):
    """Translate (ok,msg) style from tenkit into an exception.

    Raises:
        PreflightException: The message as an exception.
    """
    ok, msg = result
    if not ok:
        raise PreflightException(msg)


def is_int(s):
    try:
        int(s)
    except ValueError:
        return False
    return True


def check_gex_or_ab_present(sample_def):
    # Skip this check if library type is all VDJ.
    if all(x.get(rna_library.LIBRARY_TYPE) == "VDJ" for x in sample_def):
        return
    # At least one "Gene Expression" or "Antibody Capture" library is required.
    # Treat an empty library_type as GENE_EXPRESSION. This is done to maintain backward compatibility
    # with historical, v2 samples that didn't have a library_type
    if not any(
        x.get(rna_library.LIBRARY_TYPE) == rna_library.GENE_EXPRESSION_LIBRARY_TYPE
        or x.get(rna_library.LIBRARY_TYPE) == rna_library.ANTIBODY_LIBRARY_TYPE
        or x.get(rna_library.LIBRARY_TYPE) is None
        for x in sample_def
    ):
        raise PreflightException(
            "You must specify >= 1 input library with either library_type == '%s' or library_type == '%s' to run '{} count'".format(
                cr_env.product()
            )
            % (rna_library.GENE_EXPRESSION_LIBRARY_TYPE, rna_library.ANTIBODY_LIBRARY_TYPE)
        )


def check_os():
    tk_preflight.warn_deprecated_os()


def check_sample_def(sample_defs, feature_ref=None, pipeline=None, is_spatial=False):
    hostname = socket.gethostname()

    check(tk_preflight.check_gem_groups(sample_defs))
    if pipeline != cr_constants.PIPELINE_VDJ:
        check_gex_or_ab_present(sample_defs)
    # Check uniqueness of sample_def entries
    sd_entries = sorted(
        (sd.get("read_path"), sd.get("sample_names"), sd.get("sample_indices"), sd.get("lanes"))
        for sd in sample_defs
    )

    for i in range(len(sd_entries) - 1):
        if sd_entries[i] == sd_entries[i + 1]:
            msg = "Duplicated entry in the input FASTQ data. Please use a unique combination of fastq path and sample name."
            msg += "\nPath: %s" % sd_entries[i][0]
            msg += "\nNote in demux mode, a unique combination fastq path, sample indices, and lanes is required."
            raise PreflightException(msg)

    print("Checking FASTQ folder...", flush=True)
    for sample_def in sample_defs:
        read_path = sample_def["read_path"]
        if read_path.strip() == "":
            raise PreflightException("Empty fastq path specifed. Please specify an absolute path.")
        if not read_path.startswith("/"):
            raise PreflightException(
                "Specified FASTQ folder must be an absolute path: %s" % read_path
            )
        if not os.path.exists(read_path):
            raise PreflightException(
                f"On machine: {hostname}, specified FASTQ folder does not exist: {read_path}"
            )
        if os.path.isfile(read_path):
            raise PreflightException(
                "Specified FASTQ path {} is a file. The path is expected to be a directory.".format(
                    read_path
                )
            )
        if not os.access(read_path, os.X_OK):
            raise PreflightException(
                "On machine: {}, {} does not have permission to open FASTQ folder: {}".format(
                    hostname, cr_env.product(), read_path
                )
            )
        if not os.listdir(read_path):
            raise PreflightException("Specified FASTQ folder is empty: " + read_path)

        lanes = sample_def["lanes"]
        if lanes is not None:
            for lane in lanes:
                if not is_int(lane):
                    raise PreflightException("Lanes must be a comma-separated list of numbers.")

        check(tk_preflight.check_sample_indices(sample_def))

        if pipeline == cr_constants.PIPELINE_COUNT:
            if is_spatial:
                options = ", ".join("'%s'" % x for x in SPATIAL_ALLOWED_LIBRARY_TYPES)
            else:
                options = ", ".join("'%s'" % x for x in PUBLIC_LIBRARY_TYPES)
            library_type = sample_def.get(rna_library.LIBRARY_TYPE, None)

            # Check for empty library_type
            if library_type == "":
                msg = (
                    "library_type field may not be an empty string."
                    "\nThe 'library_type' field in the libraries csv"
                    " must be one of %s"
                ) % options
                raise PreflightException(msg)

            # Check for a valid library_type
            if (library_type is not None and library_type not in ALLOWED_LIBRARY_TYPES) or (
                is_spatial and library_type not in SPATIAL_ALLOWED_LIBRARY_TYPES
            ):
                print(library_type)
                msg = (
                    f"Unknown library_type: '{library_type}'."
                    "\nThe 'library_type' field in the libraries csv"
                    f" must be one of {options}"
                )
                raise PreflightException(msg)
            # Check that the library_type exists in the feature_ref
            if (
                feature_ref is not None
                and library_type is not None
                and library_type != rna_library.GENE_EXPRESSION_LIBRARY_TYPE
            ):
                if not any(x.feature_type == library_type for x in feature_ref.feature_defs):
                    msg = (
                        "You declared a library with library_type = '%s', but there are no features declared with that feature_type in the feature reference."
                        % library_type
                    )
                    msg += "\nCheck that the 'library_type' field in the libraries csv matches at least 1 entry in the 'feature_type' field in the feature reference csv"
                    raise PreflightException(msg)

        elif pipeline == cr_constants.PIPELINE_VDJ:
            # library type can be missing, or VDJ
            library_type = sample_def.get(rna_library.LIBRARY_TYPE, None)
            if library_type not in (None, rna_library.VDJ_LIBRARY_TYPE):
                msg = f"You declared a library with library_type = '{library_type}'. For the vdj pipeline, the library_type field in sample_def must be missing or '{rna_library.VDJ_LIBRARY_TYPE}'"
                raise PreflightException(msg)


STAR_REQUIRED_FILES = [
    "chrLength.txt",
    "chrNameLength.txt",
    "chrName.txt",
    "chrStart.txt",
    "Genome",
    "genomeParameters.txt",
    "SA",
    "SAindex",
]


def check_refdata(reference_path):
    hostname = socket.gethostname()

    if reference_path is None:
        raise PreflightException("Must specify a transcriptome reference path.")

    print(f"Checking reference_path ({reference_path}) on {hostname}...", flush=True)

    required_files = [cr_constants.REFERENCE_METADATA_FILE, cr_constants.REFERENCE_FASTA_PATH]
    for filename in required_files:
        p = os.path.join(reference_path, filename)
        if not os.path.isfile(p):
            raise PreflightException(
                "Your reference does not contain the expected files, or they are not readable. Please check your reference folder on %s."
                % hostname
            )

    # check for genes/genes.gtf or genes/getes.gtf.gz
    p_gtf = os.path.join(reference_path, cr_constants.REFERENCE_GENES_GTF_PATH)
    p_gtf_gz = os.path.join(reference_path, cr_constants.REFERENCE_GENES_GTF_PATH + ".gz")

    if not os.path.isfile(p_gtf) and not os.path.isfile(p_gtf_gz):
        raise PreflightException(
            "Your reference is missing gene annotations that should be present at {reference_path}/{gtf_path}[.gz], or they are not readable. Please check your reference folder on {hostname}.".format(
                reference_path=reference_path,
                gtf_path=cr_constants.REFERENCE_GENES_GTF_PATH,
                hostname=hostname,
            )
        )

    for filename in STAR_REQUIRED_FILES:
        p = os.path.join(reference_path, cr_constants.REFERENCE_STAR_PATH, filename)
        if not os.path.exists(p):
            raise PreflightException(
                "Your reference doesn't appear to be indexed. Please run the mkreference tool"
            )

    validate_reference(reference_path)


def expand_libraries_csv(csv_path):
    required_cols = {"fastqs", "sample", rna_library.LIBRARY_TYPE}
    try:
        reader = csv_utils.load_csv_filter_comments(csv_path, "libraries", required_cols)
    except csv_utils.CSVParseException as exc:
        raise PreflightException(str(exc)) from exc

    libraries = []

    for row in reader:
        print(row.keys())

        for key in row:
            if key is None or row[key] is None:
                msg = (
                    "Invalid libraries CSV file: incorrrect number of columns on line number (after excluding comment lines) %d"
                    % reader.line_num
                )
                raise PreflightException(msg)
            row[key] = row[key].strip()

        if row["sample"].strip() == "":
            raise PreflightException(
                "Empty sample field in libraries csv. Please specify an non-empty sample value for each library."
            )

        library = {
            "fastq_mode": "ILMN_BCL2FASTQ",
            "gem_group": None,
            "lanes": None,
            "read_path": row["fastqs"],
            "sample_names": [row["sample"]],
            rna_library.LIBRARY_TYPE: row[rna_library.LIBRARY_TYPE],
            "sample_indices": ["any"],
        }

        libraries.append(library)

    return libraries


def check_file_properties(path, file_desc="input"):
    """Check that the input file exists and is readable."""
    if not os.access(path, os.R_OK):
        raise PreflightException(
            f"{file_desc} file is not readable, please check file permissions: {path}"
        )

    if os.stat(path).st_size == 0:
        raise PreflightException(f"{file_desc} file {path} is empty.")


def try_load_feature_ref(transcriptome_ref_path, feature_ref_path):
    """Try creating the feature reference object.

    Raises:
        PreflightException: wrapped FeatureDefException.
    """
    try:
        feature_ref = rna_feature_ref.from_transcriptome_and_csv(
            transcriptome_ref_path, feature_ref_path
        )

        rna_feature_ref.FeatureExtractor(feature_ref)
        return feature_ref

    except FeatureDefException as ex:
        raise PreflightException(str(ex)) from ex


def check_environment():
    check(tk_preflight.check_open_fh())


class _VersionCmd(NamedTuple):
    name: str
    cmd: list[str]


_PACKAGE_VERSION_CMDS = [
    _VersionCmd(name="mro", cmd=["mro", "--version"]),
    _VersionCmd(name="mrp", cmd=["mrp", "--version"]),
    _VersionCmd(name="Anaconda", cmd=["python", "--version"]),
    _VersionCmd(name="numpy", cmd=["python", "-c", "import numpy; print(numpy.__version__)"]),
    _VersionCmd(name="scipy", cmd=["python", "-c", "import scipy; print(scipy.__version__)"]),
    _VersionCmd(name="pysam", cmd=["python", "-c", "import pysam; print(pysam.__version__)"]),
    _VersionCmd(name="h5py", cmd=["python", "-c", "import h5py; print(h5py.__version__)"]),
    _VersionCmd(name="pandas", cmd=["python", "-c", "import pandas; print(pandas.__version__)"]),
    _VersionCmd(name="STAR", cmd=["STAR", "--version"]),
    _VersionCmd(name="samtools", cmd=["samtools", "--version"]),
]


def record_package_versions(dest=sys.stdout):
    """Print the versions of various dependencies to the given destination.

    Args:
        file (file-like object): The destination to which to print the versions.
    """
    for package in _PACKAGE_VERSION_CMDS:
        name = package.name
        cmd = package.cmd

        version = tk_subproc.check_output(cmd).strip().split(b"\n", maxsplit=1)[0].decode()
        print(name, version, sep=": ", file=dest, flush=True)


def check_read_length(x):
    if x < 1:
        raise PreflightException("Specified read length must be greater than or equal to 1.")


#############################################################################################
def check_feature_preflights(sample_def, feature_ref_path):
    """If any non "Gene Expression" libraries are present then the feature-ref is required."""

    def is_not_gene_expression(sample):
        library_type = sample.get("library_type")
        if library_type is None:
            return False
        return library_type != rna_library.GENE_EXPRESSION_LIBRARY_TYPE

    if any(is_not_gene_expression(x) for x in sample_def):
        if feature_ref_path is None:
            raise PreflightException(
                "You must specify --feature-ref when using feature barcoding libraries."
            )

    def is_crispr_guide_capture(sample):
        library_type = sample.get("library_type")
        if library_type is None:
            return False
        return library_type == rna_library.CRISPR_LIBRARY_TYPE

    def is_antigen_capture(sample):
        library_type = sample.get("library_type")
        if library_type is None:
            return False
        return library_type == rna_library.ANTIGEN_LIBRARY_TYPE

    if any(is_antigen_capture(x) for x in sample_def):
        raise PreflightException(
            "Analysis of Antigen Capture libraries is unsupported in this product."
        )

    if feature_ref_path:
        try:
            feature_defs, _ = rna_feature_ref.parse_feature_def_file(
                feature_ref_path, index_offset=0
            )
        except FeatureDefException as exc:
            raise PreflightException(str(exc)) from exc

        # if CRISPR Guide Capture library is present, and there is a non targeting feature, make sure it's spelled correctly
        if any(is_crispr_guide_capture(x) for x in sample_def):
            rna_feature_ref.check_crispr_target_gene_name(feature_defs)


MINIMUM_TARGET_PANEL_GENES_PREFLIGHT = 10


def check_targeting_preflights(
    sample_def: dict[str, Any],
    reference_path: str,
    feature_reference_path: str | None,
    *,
    parse_files: bool,
    expected_targeting_method: str,
    is_spatial: bool,
):
    """_summary_.

    Args:
        sample_def (dict[str, Any]): Sample definition list
        reference_path (str): Location of the reference
        feature_reference_path (Optional[str]): Location of the feature reference path
        parse_files (bool): Flag to parse the GTF to verify that the gene is present
        expected_targeting_method (str): Type of targeting method used (Need more details)
        is_spatial (bool): Is this a spatial sample or not?

    Raises:
        PreflightException: Any Preflight specific Exception
    """
    # Get distinct, non-null target panel files
    target_panel_files = list(
        {cr_sample_def.get_target_set(sd) for sd in sample_def}.difference({None})
    )

    if len(target_panel_files) == 0:
        # Nothing to check
        return
    elif len(target_panel_files) > 1:
        # Should not be possible, except with an edited MRO file
        raise PreflightException(
            "Found multiple distinct non-null target panels in the sample_def."
        )

    target_panel = target_panel_files.pop()  # single-element set
    check_file_properties(target_panel, file_desc="The target panel or probe set csv")

    try:
        csv_utils.load_csv_filter_comments(
            filename=target_panel,
            descriptive_name=TARGETING_METHOD_FILE_NAMES.get(
                expected_targeting_method, "target panel or probe set"
            ),
            required_cols=["gene_id"],  # always required
        )
    except csv_utils.CSVParseException as err:
        raise PreflightException(str(err)) from err

    if expected_targeting_method is not None:
        check_targeting_method(expected_targeting_method, [target_panel])

    if parse_files:
        gene_index = cr_reference.NewGeneIndex.load_from_reference(reference_path)
        gene_name_to_id = {gene.name: gene.id for gene in gene_index.genes}

        try:
            target_panel_metadata, target_panel_genes, _ = tgt_simple_utils.parse_target_csv(
                target_panel,
                ref_gene_index=gene_index,
                gene_name_to_id=gene_name_to_id,
                expected_targeting_method=expected_targeting_method,
            )
        except Exception as err:
            raise PreflightException(str(err)) from err

        descriptive_name = (
            TARGETING_METHOD_FILE_NAMES.get(expected_targeting_method, expected_targeting_method)
            if expected_targeting_method is not None
            else "target panel or probe set"
        )
        if len(target_panel_genes) < MINIMUM_TARGET_PANEL_GENES_PREFLIGHT:
            raise PreflightException(
                "Ten or more genes must be specified in the {} file for compatibility with downstream analysis. Number of genes found: {}.".format(
                    descriptive_name,
                    len(target_panel_genes),
                )
            )

        # Ensure that the reference genome of the transcriptome and probe set are identical for RTL only.
        if TARGETING_METHOD_TL_FILE_FORMAT in target_panel_metadata:
            transcriptome_reference_genome = cr_reference.get_ref_name_from_genomes(
                get_reference_genomes(reference_path)
            )
            probe_set_reference_genome = target_panel_metadata["reference_genome"]
            if not transcriptome_reference_genome == probe_set_reference_genome:
                raise PreflightException(
                    "The reference genome of the transcriptome '{}' and probe set '{}' must be identical.".format(
                        transcriptome_reference_genome, probe_set_reference_genome
                    )
                )

        # Ensure that the reference transcriptome and probe set have at least one gene in common.
        if all(gene_index.gene_id_to_int(gene_id) is None for gene_id in target_panel_genes):
            raise PreflightException(
                "There are no gene IDs in common between the reference transcriptome and the {}.".format(
                    descriptive_name
                ),
            )

        panel_type = target_panel_metadata.get("panel_type", None)
        if is_spatial and panel_type in SPATIAL_TARGET_DISALLOWED_PANEL_TYPES:
            martian.alarm(
                'The provided targeted panel type "{}" is UNSUPPORTED in this product and results may be incorrect.'.format(
                    panel_type
                )
            )

        if feature_reference_path is not None:
            gene_id_to_index = {gene.id: id(gene) for gene in gene_index.genes}
            target_set_gene_indices = [
                gene_id_to_index[gene_id]
                for gene_id in target_panel_genes
                if gene_id in gene_id_to_index  # handles features that are not in the gene index
            ]
            target_sets_dict = {"": sorted(target_set_gene_indices)}
            csv_feature_defs, _ = rna_feature_ref.parse_feature_def_file(
                feature_reference_path, index_offset=len(gene_index.genes)
            )
            try:
                rna_feature_ref.check_crispr_target_gene(
                    csv_feature_defs, gene_index.genes, target_sets_dict
                )
            except FeatureDefException as e:
                raise PreflightException(str(e)) from e


def validate_targeted_compare_mol_info(
    mol_info_fn, expecting_targeted, required_metrics=[], required_library_metrics=[]
):
    if not os.path.exists(mol_info_fn):
        raise PreflightException("The molecule info file %s does not exist" % mol_info_fn)
    if not os.access(mol_info_fn, os.R_OK):
        raise PreflightException("The molecule info file %s is not readable" % mol_info_fn)
    try:
        mc = MoleculeCounter.open(mol_info_fn, "r")

        molecule_info_type = mc.get_molecule_info_type()
        if molecule_info_type == MOLECULE_INFO_TYPE_PERSAMPLE:
            raise PreflightException(
                "The multiplexed sample molecule info file %s cannot be used with targeted-compare"
                % mol_info_fn
            )

        library_info = mc.get_library_info()
        gex_lib_is_targeted = [
            rna_library.has_target_set(lib)
            for lib in library_info
            if lib["library_type"] == rna_library.GENE_EXPRESSION_LIBRARY_TYPE
        ]
        if len(gex_lib_is_targeted) == 0:
            raise PreflightException(
                "The molecule info file %s does not have any gene expression libraries."
                % mol_info_fn
            )

        library_metrics = next(iter(mc.get_all_metrics()[LIBRARIES_METRIC].values()))
        for metric in required_library_metrics:
            if metric not in library_metrics.keys():
                raise PreflightException(
                    f"The molecule info file {mol_info_fn} is too old and is missing the metric {metric}. Please rerun."
                )

        for metric_key in required_metrics:
            if mc.get_metric(metric_key) is None:
                msg = (
                    "The molecule info HDF5 file (%s) was produced by an older software version."
                    % mol_info_fn
                )
                msg += "\nReading these files is unsupported."
                raise PreflightException(msg)

        if expecting_targeted and (
            not all(gex_lib_is_targeted) or mc.get_metric("target_panel_hash") is None
        ):
            msg = (
                "The input targeted molecule info file %s does not come from a targeted analysis run. "
                % (mol_info_fn)
            )
            msg += (
                "Please rerun your targeted analysis run using '{} count --target-panel'.".format(
                    cr_env.product()
                )
            )
            raise PreflightException(msg)
        if not expecting_targeted and any(gex_lib_is_targeted):
            msg = (
                "The input parent molecule info file %s does not come from a whole transcriptome analysis run. "
                % (mol_info_fn)
            )
            msg += "Please rerun your parent analysis run using '{} count' without a target-panel.".format(
                cr_env.product()
            )
            raise PreflightException(msg)

        mc.close()

    except OSError as err:
        raise PreflightException(
            "Molecule info file %s cannot be loaded (may be corrupt)." % mol_info_fn
        ) from err
    except ValueError as err:
        raise PreflightException(str(err)) from err


def check_molecule_info_contains_library_metrics(mol_info_fn, metrics_expected=[]):
    with MoleculeCounter.open(mol_info_fn, "r") as mc:
        num_libraries = mc.get_library_info()
        if num_libraries == 0:
            raise PreflightException(
                "The molecule info file %s has no associated libraries." % mol_info_fn
            )
        # just check the any single library since it'll be the same for all
        for metric in metrics_expected:
            library_metrics = next(iter(mc.get_all_metrics()[LIBRARIES_METRIC].values()))
            if metric not in library_metrics.keys():
                raise PreflightException(
                    f"The molecule info file {mol_info_fn} is too old and is missing the metric {metric}. Please rerun."
                )


def check_target_features_same(target_sets):
    """Takes in a list of target sets."""
    # TODO: make this check the mol_info hash instead. should be able to just reuse
    # wtv function we use in aggr
    counter = 0
    global_target_set = None
    inconsistent_genes = set()
    for target_set in target_sets:
        if counter == 0:
            global_target_set = target_set
        elif target_set != global_target_set:
            inconsistent_genes.update(target_set.difference(global_target_set))
            inconsistent_genes.update(global_target_set.difference(target_set))
            raise PreflightException(
                "Target sets are not the same. Found inconsistent genes: %s" % (inconsistent_genes)
            )
        counter += 1


def check_sample_info(
    sample_def, reference_path, full_check, feature_ref_path=None, is_spatial=False
):
    feature_ref = None
    if feature_ref_path is not None:
        print("Checking feature definition file...", flush=True)
        check_file_properties(feature_ref_path, file_desc="feature reference")
        if full_check:
            # This requires loading the reference index and parsing the feature CSV,
            # and we don't want to do that e.g. during local preflight
            feature_ref = try_load_feature_ref(reference_path, feature_ref_path)

    print("Checking sample info...", flush=True)
    check_sample_def(
        sample_def, feature_ref, pipeline=cr_constants.PIPELINE_COUNT, is_spatial=is_spatial
    )


def check_targeting_method(targeting_method, target_sets):
    """Determines the targeting method specified by a given panel or target file and compares the result to the specified value."""
    if targeting_method is None:  # targeting method only passed in by CS pipeline for this check
        return

    if len(target_sets) > 1:
        raise PreflightException(f"Multiple different target sets found:\n\t{target_sets}")
    elif len(target_sets) == 0:
        raise PreflightException("A target or probe set was expected, but not found")

    expected_name = TargetingMethod.get_file_name(targeting_method)
    check_file_properties(target_sets[0], file_desc=f"The {expected_name} csv")

    try:
        target_set_metadata = tgt_simple_utils.load_target_csv_metadata(
            target_sets[0], expected_name
        )

        if len(target_set_metadata) == 0:
            raise PreflightException(
                f"The {expected_name} file {target_sets[0]} contains no metadata."
            )

        target_method_info = tgt_simple_utils.determine_targeting_method_info_from_metadata(
            target_set_metadata, targeting_method
        )

        if target_method_info is None:
            raise PreflightException(
                f"No file format version specified in the {expected_name} file header. Please check your {expected_name} file."
            )

        if target_method_info.method != targeting_method:
            raise PreflightException(
                f"The header in the {expected_name} file indicates this is a {target_method_info.file_format_tag} file. Please check your {expected_name} file."
            )

        (
            _,
            _,
            required_metadata,
            conflicting_metadata,
        ) = tgt_simple_utils.get_target_panel_or_probe_set_file_format_spec(
            targeting_method, target_method_info.file_version
        )

        tgt_simple_utils.check_target_csv_metadata(
            target_sets[0],
            TargetingMethod.get_file_name(targeting_method),
            required_metadata,
            conflicting_metadata,
        )
    except csv_utils.CSVParseException as err:
        raise PreflightException(str(err)) from err


def check_common_preflights(
    full_check,
    reference_path,
    r1_length,
    r2_length,
    recovered_cells,
    force_cells,
):
    print("Checking reference...", flush=True)
    check_refdata(reference_path)

    if r1_length is not None:
        print("Checking read 1 length...", flush=True)
        check_read_length(r1_length)
    if r2_length is not None:
        print("Checking read 2 length...", flush=True)
        check_read_length(r2_length)

    # Open file handles limit - per CELLRANGER-824, only check this on the execution machine.
    # We can tell if we're on the execution machine by looking at args.full_check
    if full_check:
        print("Checking system environment...", flush=True)
        check_environment()

    print("Checking optional arguments...", flush=True)
    recovered_cells = recovered_cells["per_gem_well"] if recovered_cells else None
    force_cells = force_cells["per_gem_well"] if force_cells else None
    if recovered_cells is not None and force_cells is not None:
        raise PreflightException(
            "Cannot specify both --force-cells and --expect-cells in the same run."
        )
