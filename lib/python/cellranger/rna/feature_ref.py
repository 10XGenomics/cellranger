#!/usr/bin/env python
#
# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
#

from __future__ import annotations

import os
import re
import string
from collections import OrderedDict
from collections.abc import Iterable
from re import Pattern
from typing import NamedTuple

from six import ensure_binary, ensure_str

import cellranger.constants
import cellranger.cr_io as cr_io
import cellranger.reference as cr_reference
import cellranger.rna.library as rna_library
import cellranger.utils as cr_utils
from cellranger import csv_utils
from cellranger.feature_ref import (
    DEFAULT_FEATURE_TAGS,
    GENOME_FEATURE_TAG,
    RESERVED_TAGS,
    FeatureDef,
    FeatureDefException,
    FeatureReference,
)
from cellranger.reference_paths import get_reference_genomes


class PatternEntry(NamedTuple):
    read: str
    regex_string: str
    regex: Pattern[bytes]
    barcode_dict: dict[bytes, FeatureDef]


# These are base FeatureDef attributes
BASE_FEATURE_FIELDS = ["id", "name", "feature_type"]

# These feature tag keys are required for feature barcode extraction
REQUIRED_TAGS = ["read", "pattern", "sequence"]


# Type conversion and validation
def convert_bool_tag(tagstr: str, tag_name: str) -> bool:
    """Validates and converts a boolean column (tag) and returns a bool.

    Args:
        tagstr (str): the value of the tag as a string
        tag_name (str): the tag name use for exception handling
    """
    ntagstr = tagstr.strip().lower()
    if ntagstr == "false":
        return False
    elif ntagstr != "true":
        raise FeatureDefException(
            f"Feature reference column/tag {tag_name} should be TRUE or FALSE, but has unexpected value: {tagstr}"
        )

    return True


TAGS_VALIDATION = {"isotype_control": convert_bool_tag}


ALLOWED_FEATURE_TYPES = [
    rna_library.CRISPR_LIBRARY_TYPE,
    rna_library.ANTIBODY_LIBRARY_TYPE,
    rna_library.ANTIGEN_LIBRARY_TYPE,
    rna_library.MULTIPLEXING_LIBRARY_TYPE,
    rna_library.CUSTOM_LIBRARY_TYPE,
]

PUBLIC_FEATURE_TYPES = [
    rna_library.CRISPR_LIBRARY_TYPE,
    rna_library.ANTIBODY_LIBRARY_TYPE,
    rna_library.CUSTOM_LIBRARY_TYPE,
]


class FeatureMatchResult(NamedTuple):
    """The result of a feature matching operation.

    0  whitelist hits - (raw, [])   where 'raw' is the longest potential sequence.
    1  whitelist hit  - (raw, [id])
    >1 whitelist hits - (raw, [ids])
    """

    barcode: bytes | None
    qual: bytes | None
    ids: list[bytes]
    indices: list[int]


ALLOWED_READS = ["R1", "R2"]
MULTI_UNICITY_CHECK = ["read", "pattern", "sequence"]


class FeatureExtractor:
    """Store a list of features (genes, antibodies, etc) and extract feature barcodes from a read sequence."""

    def __init__(self, feature_ref: FeatureReference, use_feature_types: list[str] | None = None):
        """Setup feature extraction.

        Args:
            feature_ref (FeatureReference): All feature definitions.
            use_feature_types (list of str): Feature types to extract.
        """
        self.feature_ref = feature_ref

        self.patterns = FeatureExtractor._compile(feature_ref.feature_defs, use_feature_types)

    @staticmethod
    def _compile(
        feature_defs: list[FeatureDef], use_feature_types: list[str] | None = None
    ) -> dict[tuple[str, str], PatternEntry]:
        """Prepare for feature barcode extraction.

        Args:
            feature_defs (list of FeatureDef): feature defs to compile
            use_feature_types (list of str): feature types to try to extract.
                                              None to extract all.

        Returns:
            dict of (str,str) to PatternEntry
        """
        # Store distinct patterns
        patterns: dict[tuple[str, str], PatternEntry] = {}

        for fd in feature_defs:
            # Only build patterns from searchable feature types
            if use_feature_types is not None and fd.feature_type not in use_feature_types:
                continue

            sequence = fd.tags.get("sequence")
            pattern = fd.tags.get("pattern")
            read = fd.tags.get("read")

            if not sequence:
                continue
            if not pattern:
                raise FeatureDefException(
                    f"Feature definition for {fd.id} has a sequence but no pattern specifying how to extract it."
                )
            if not read:
                raise FeatureDefException(
                    f"Feature definition for {fd.id} has a sequence but no read specifying where to extract it from."
                )

            regex_str, regex = compile_pattern(pattern, len(sequence))

            # Treat equal regex pattern strings as generating equivalent automata
            pattern_key = (read, regex_str)
            if pattern_key not in patterns:
                patterns[pattern_key] = PatternEntry(read, regex_str, regex, {})
            entry = patterns[pattern_key]

            if sequence in entry.barcode_dict:
                raise FeatureDefException(
                    f'Found two feature definitions with the same read ("{read}"), pattern ("{pattern}") and barcode sequence ("{sequence}"). This combination of values must be unique for each row in the feature definition file.'
                )

            entry.barcode_dict[sequence] = fd

        # Built from searchable features only
        return patterns

    def get_read_types(self) -> list[str]:
        """Determine which read types need to be inspected.

        Returns:
          (list of str): List of read types, e.g. ['R1']
        """
        return sorted(
            {fd.tags.get("read") for fd in self.feature_ref.feature_defs if fd.tags.get("read")}
        )

    @staticmethod
    def _filter_feature_matches(
        matches: list[tuple[bytes, bytes, FeatureDef | None, int, None]],
        any_whitelist_hits: bool,
    ) -> FeatureMatchResult:
        """Filter a list of feature matches by prioritization.

        Args:
            matches (list of (str, str, FeatureDef, int)): bc, qual, hit, min feature index, hit (corrected)
            any_whitelist_hits (bool): True if any of the hits matched the whitelist.

        Returns:
           FeatureMatchResult
        """
        if any_whitelist_hits:
            # >0 whitelist hits. Remove entries w/ no whitelist hit.
            matches = [bc_hit for bc_hit in matches if bc_hit[2]]

            # Take the longest observed sequence as the canonical sequence
            # NOTE: If there are multiple whitelist hits of equal length,
            #   we choose the first feature in the feature reference (by index).
            # However, information is not lost because the feature IDs tag
            #   contains a list of each hit's unique feature ID.
            best_hit = max(matches, key=lambda bc_hit: (len(bc_hit[0]), -bc_hit[3]))
            barcode, qual = best_hit[0:2]

            # Record ids and indices of all whitelist hits (typically 1)
            ids = [bc_hit[2].id for bc_hit in matches]
            indices = [bc_hit[2].index for bc_hit in matches]

            return FeatureMatchResult(barcode, qual, ids, indices)

        elif len(matches) > 0:
            # if there are correctable barcodes, use those instead
            if any(bc_hit[4] for bc_hit in matches):
                matches = [bc_hit for bc_hit in matches if bc_hit[4]]

            # No whitelist hits but found some sequence.
            # Whether or not we're only looking at corrected hits,
            # take the longest observed or earliest in feature reference.
            # NOTE: This may do the wrong thing if:
            #     1) There are multiple possible barcode lengths for a given pattern,
            #     2) This sequence doesn't match the whitelist at any barcode length,
            #     3) There are multiple possible corrections,
            best_hit = max(matches, key=lambda bc_hit: (len(bc_hit[0]), -bc_hit[3]))
            barcode, qual = best_hit[0:2]

            return FeatureMatchResult(barcode, qual, [], [])

        else:
            return FeatureMatchResult(None, None, [], [])

    def has_features_to_extract(self):
        """Does this class have any interesting work to do.

        Returns:
           (bool): True if this feature ref was compiled to extract any feature BCs.
        """
        return len(self.patterns) > 0

    @staticmethod
    def get_feature_type(row: dict[str, str]) -> str:
        if row["feature_type"] == rna_library.FEATURETEST_LIBRARY_TYPE:
            return rna_library.MULTIPLEXING_LIBRARY_TYPE
        return row["feature_type"]


def save_features_tsv(feature_ref, base_dir, compress):
    """Save a FeatureReference to a tsv file."""
    out_features_fn = os.path.join(ensure_binary(base_dir), b"features.tsv")
    if compress:
        out_features_fn += b".gz"

    with cr_io.open_maybe_gzip(out_features_fn, "wb") as f:
        for feature_def in feature_ref.feature_defs:
            f.write(feature_def.id)
            f.write(b"\t")
            f.write(ensure_binary(feature_def.name))
            f.write(b"\t")
            f.write(ensure_binary(feature_def.feature_type))
            f.write(b"\n")


def from_transcriptome_and_csv(
    gene_ref_path, feature_def_filename, target_features: dict[str, Iterable[int]] | None = None
) -> FeatureReference:
    """Create a FeatureReference.

    Create a FeatureReference from a transcriptome ref and a feature barcode ref.

    Args:
        gene_ref_path (str): Path to transcriptome reference. Can be None.
        feature_def_filename (str): Path to Feature Definition CSV file. Can be None.
        target_features (dictionary of list of int): Optional target set(s). Each target set
            should be a list of integers with a dictionary key representing the name. This name
            will also exist in the library_info.

    Returns:
        FeatureReference
    """
    # Load gene info
    feature_defs: list[FeatureDef] = []
    all_tag_keys = DEFAULT_FEATURE_TAGS[:]

    genomes: list[str] = get_reference_genomes(gene_ref_path)

    if gene_ref_path is not None:
        gene_index = cr_reference.NewGeneIndex.load_from_reference(gene_ref_path)

        # Stuff relevant fields of Gene tuple into FeatureDef
        for gene in gene_index.genes:
            assert isinstance(gene.id, bytes)
            genome = cr_utils.get_genome_from_str(ensure_str(gene.id), genomes)
            assert genome is not None
            fd = FeatureDef(
                index=len(feature_defs),
                id=gene.id,
                name=gene.name,
                feature_type=rna_library.GENE_EXPRESSION_LIBRARY_TYPE,
                tags={GENOME_FEATURE_TAG: genome},
            )
            feature_defs.append(fd)

    # Load feature definition file
    if feature_def_filename is not None:
        csv_feature_defs, csv_tag_keys = parse_feature_def_file(
            feature_def_filename, index_offset=len(feature_defs)
        )

        # check the CRISPR 'target_gene_id' field, if it exists
        # it needs to match a transcriptome entry
        # If target_features are present, also require that a target_gene_id
        # is contained in the target_features (i.e., Targeted Genes)
        check_crispr_target_gene(csv_feature_defs, feature_defs, target_features)

        feature_defs.extend(csv_feature_defs)
        all_tag_keys.extend(csv_tag_keys)

    return FeatureReference(feature_defs, all_tag_keys, target_features)


def check_crispr_target_gene_name(features: list[FeatureDef]):
    """Check that 'target_gene_name' and 'target_gene_id' fields specified on CRISPR.

    guides say exactly 'Non-Targeting', if applicable
    """
    for feature in features:
        if feature.feature_type == rna_library.CRISPR_LIBRARY_TYPE:
            # Check for case-insensitive "Non-Targeting" variants
            target_gene_name = feature.tags.get("target_gene_name")
            target_gene_id = feature.tags.get("target_gene_id")
            if not target_gene_name or not target_gene_id:
                continue

            reduced_gene_name = "".join(
                [char for char in target_gene_name if char not in string.punctuation]
            ).lower()
            if reduced_gene_name == "nontargeting" and target_gene_name != "Non-Targeting":
                msg = f"CRISPR: '{target_gene_name}' was specified as the target_gene_name for the feature {feature.name} in the feature reference, but it must be spelled exactly as 'Non-Targeting'."
                raise FeatureDefException(msg)
            reduced_gene_id = "".join(
                [char for char in target_gene_id if char not in string.punctuation]
            ).lower()
            if reduced_gene_id == "nontargeting" and target_gene_id != "Non-Targeting":
                msg = f"CRISPR: '{target_gene_id}' was specified as the target_gene_id for the feature {feature.name} in the feature reference, but it must spelled exactly as 'Non-Targeting'."
                raise FeatureDefException(msg)


def check_crispr_target_gene(
    features: list[FeatureDef],
    genes: list[FeatureDef],
    target_features: dict[str, Iterable[int]] | None = None,
):
    """Check that 'target_gene_id' and 'target_gene_name' fields specified on CRISPR.

    guides match a gene declared in the transcriptome. If 'target_features' are
    present, also require that each 'target_gene_id' is contained in 'target_features'
    """
    gene_id_name_map: dict[bytes, bytes] = {
        ensure_binary(g.id): ensure_binary(g.name) for g in genes
    }
    gene_id_index_map: dict[bytes, int] = {ensure_binary(g.id): id(g) for g in genes}

    special_ids: list[bytes] = [ensure_binary(x.lower()) for x in cellranger.constants.FILTER_LIST]

    for feat in features:
        target_id_str = feat.tags.get("target_gene_id")
        if target_id_str:
            target_id: bytes = ensure_binary(target_id_str)
            # Cut-out for special CRISPR guide types
            if target_id.lower() in special_ids:
                continue

            if target_id not in gene_id_name_map:
                msg = f"CRISPR: This target_gene_id ({ensure_str(target_id)}) declared for one or more guide RNAs in the feature reference does not exist in the transcriptome."
                msg += "\nPlease specify a target_gene_id that exists in the reference, or use the string 'Non-Targeting' to indicate a control guide."
                raise FeatureDefException(msg)

            if target_features is not None:
                target_index = gene_id_index_map[target_id]
                # target_features maps target set name to list of gene indices
                if not any(target_index in target_set for target_set in target_features.values()):
                    msg = (
                        f"CRISPR: {ensure_str(target_id)} was specified as the target_gene_id for a guide RNA in the feature reference, "
                        + "but this gene is not specified in the gene_id column of the target panel csv file."
                    )
                    raise FeatureDefException(msg)

            target_name = feat.tags.get("target_gene_name")
            if target_name is None or target_name == "":
                msg = f"CRISPR: No target_gene_name specified for this target_gene_id ({ensure_str(target_id)}) in the feature reference."
                raise FeatureDefException(msg)

            target_name = ensure_binary(target_name)
            if gene_id_name_map[target_id] != target_name:
                msg = f"CRISPR: You specified target_gene_id = {ensure_str(target_id)} and target_gene_name = {ensure_str(target_name)} in the feature reference.\n"
                msg += f"The transcriptome reference has gene_id = {ensure_str(target_id)} with name = {ensure_str(gene_id_name_map[target_id])}. "
                msg += "Please ensure the target_gene_name field has the correct gene name that matches the transcriptome."
                raise FeatureDefException(msg)


def validate_sequence(seq: str):
    """Validate a feature barcode sequence."""
    if len(seq) == 0:
        raise FeatureDefException("Feature sequence must be non-empty.")
    if not re.match("^[ACGTN]+$", seq):
        raise FeatureDefException(
            f'Invalid sequence: "{seq}". The only allowed characters are A, C, G, T, and N.'
        )


def compile_pattern(pattern_str, length: int):
    """Compile a feature definition pattern into a regex."""
    if "(BC)" not in pattern_str:
        raise FeatureDefException(
            f'Invalid pattern: "{pattern_str}". The pattern must contain the string "(BC)".'
        )

    # keep for reporting errors against the initial input
    input_pattern_str = pattern_str
    # We previously only supported the regex ^/$ markers for start and end,
    # but in order to be more human-friendly, we now also support 5p and 3p too
    # -- replace them here to make a valid regex
    if re.search(r"^5[Pp]?[-_]?", pattern_str):
        pattern_str = re.sub(r"^5[Pp]?[-_]?", "^", pattern_str)
    if re.search(r"[-_]?3[Pp]?$", pattern_str):
        pattern_str = re.sub(r"[-_]?3[Pp]?$", "$", pattern_str)

    check_pattern = re.sub(r"\(BC\)", "", pattern_str)
    if not re.match(r"^\^{0,1}[ACGTN]*\${0,1}$", check_pattern):
        raise FeatureDefException(
            f'Invalid pattern: "{input_pattern_str}". The pattern must optionally start with "5P", optionally end with "3P", contain exactly one instance of the string "(BC)" and otherwise contain only the characters A, C, G, T, and N.'
        )

    # Allow Ns to match anything
    regex_str = re.sub("N", ".", pattern_str)

    # Capture the feature barcode
    regex_str = re.sub(r"\(BC\)", "(.{%d,%d})" % (length, length), regex_str)

    try:
        regex = re.compile(regex_str)
    except re.error:
        raise FeatureDefException(
            f'Failed to compile the feature definition string "{regex_str}" into a regular expression: {re.error!s}.'
        )

    return regex_str, regex


def get_required_csv_columns() -> list[str]:
    """Get a list of required CSV columns."""
    return BASE_FEATURE_FIELDS + REQUIRED_TAGS


def parse_feature_def_file(filename, index_offset: int = 0) -> tuple[list[FeatureDef], list[str]]:
    """Load a CSV feature definition file.

    Args:
        filename (str): CSV filename.
        index_offset (int): Start the FeatureDef indices at this value.

    Returns:
        feature_defs (list of FeatureDef): All feature definitions
        all_tag_keys (list of str): All tag column names in original order
    """
    feature_defs: list[FeatureDef] = []
    seen_ids: set[bytes] = set()
    seen_multi_unicity: set[bytes] = set()
    required_cols = get_required_csv_columns()
    try:
        reader = csv_utils.load_csv_filter_comments(filename, "feature reference", required_cols)
    except csv_utils.CSVParseException as exc:
        raise FeatureDefException(str(exc)) from exc
    feature_index = index_offset
    all_tag_keys = [c for c in reader.fieldnames if c not in BASE_FEATURE_FIELDS and c != ""]

    # generator to rewrite reader CSVParseException as FeatureDefException
    # to avoid the world's longest try block below
    def feature_reader(reader):
        try:
            yield from reader
        except csv_utils.CSVParseException as exc:
            raise FeatureDefException(str(exc)) from exc

    for row_num, row in enumerate(feature_reader(reader)):
        # Check that there aren't extra columns, which you
        for key in row:
            if key is None:
                msg = (
                    "Your feature reference csv file contains more columns than the header on line %d. Please use a csv file with a header for each column."
                    % (row_num + 2)
                )
                msg += "\nYou might have an comma character in a field. Commas are permitted in some fields, but fields containing commas must be enclosed in quotes."
                raise FeatureDefException(msg)

        # Strip flanking whitespace from values
        # tolerate missing values
        for key in row:
            if row[key] is None:
                msg = (
                    "Your feature reference csv file is missing fields on line %d. Ensure that each row has an entry for every header column."
                    % (row_num + 2)
                )
                msg += "\nYou might have an extra comma character in the CSV header."
                raise FeatureDefException(msg)

            row[key] = row[key].strip()

        # Check for a valid library_type
        if FeatureExtractor.get_feature_type(row) not in ALLOWED_FEATURE_TYPES:
            options = " or ".join(f"'{x}'" for x in PUBLIC_FEATURE_TYPES)
            msg = (
                f"Unknown feature_type: '{FeatureExtractor.get_feature_type(row)}'."
                "\nThe 'feature_type' field in the feature reference"
                f" must be one of {options}"
            )

            raise FeatureDefException(msg)

        # Check ID uniqueness
        f_id: bytes = ensure_binary(row["id"])
        if id in seen_ids:
            raise FeatureDefException(
                'Found duplicated ID in feature reference file: "{}"'.format(row["id"])
            )
        seen_ids.add(f_id)

        multi_unicity = "-".join(map(row.get, MULTI_UNICITY_CHECK))
        if multi_unicity in seen_multi_unicity:
            read_value = row.get("read")
            pattern_value = row.get("pattern")
            sequence_value = row.get("sequence")
            raise FeatureDefException(
                f"Found two feature definitions with the same read ({read_value}), pattern ({pattern_value}) and barcode sequence ({sequence_value}). This combination of values must be unique for each row in the feature definition file."
            )
        seen_multi_unicity.add(multi_unicity)
        if "\t" in row["name"]:
            raise FeatureDefException(
                "Feature name field cannot contain tabs: '{}'".format(row["name"])
            )

        allowed_id_chars = set(string.printable) - set(string.whitespace) - set("/,'\"\\`")

        for idx, c in enumerate(row["id"]):
            if not c in allowed_id_chars:
                if c in string.whitespace:
                    raise FeatureDefException(
                        "Feature id field cannot contain whitespace: '{}'".format(row["id"])
                    )
                else:
                    msg = "Feature id field contains an illegal character at position %d: '%s'" % (
                        idx,
                        row["id"],
                    )
                    msg += "\nFeature IDs may only ASCII characters, and must not use whitespace slash, quote or comma characters."
                    raise FeatureDefException(msg)

        # Additional columns become key-value pairs
        # Maintain input order
        tag_cols = [c for c in reader.fieldnames if c not in BASE_FEATURE_FIELDS and c in row]

        tags = OrderedDict()
        for key in tag_cols:
            if key in RESERVED_TAGS:
                raise FeatureDefException(
                    f'Found invalid column name "{key}." This name cannot be used as a custom feature tag because it is reserved for internal use.'
                )
            elif key in TAGS_VALIDATION:
                row[key] = TAGS_VALIDATION[key](row[key], key)

            # Always store as strings
            tags[key] = str(row[key])

        # Validate fields
        if len(tags["sequence"]) == 0:
            raise FeatureDefException(
                "Found blank feature barcode sequence for feature id {}. The sequence column must be populated for this feature.".format(
                    row["id"]
                )
            )
        validate_sequence(tags["sequence"])

        if len(tags["pattern"]) == 0:
            raise FeatureDefException(
                "Found blank feature barcode pattern for feature id {}. The pattern column must be populated for this feature.".format(
                    row["id"]
                )
            )
        compile_pattern(tags["pattern"], len(tags["sequence"]))

        if tags["read"] not in ALLOWED_READS:
            raise FeatureDefException(
                'The feature definition file contains a read type value "{}" which is not one of the allowed read types {}.'.format(
                    tags["read"], str(ALLOWED_READS)
                )
            )

        feature_defs.append(
            FeatureDef(
                index=feature_index,
                id=f_id,
                name=row["name"],
                feature_type=FeatureExtractor.get_feature_type(row),
                tags=tags,
            )
        )

        feature_index += 1

    return feature_defs, all_tag_keys
