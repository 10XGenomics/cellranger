#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

from __future__ import annotations

import os
from typing import NamedTuple, TypedDict

from six import ensure_binary, ensure_str

import cellranger.constants as cr_constants
import cellranger.utils as cr_utils
import cellranger.vdj.chain_types as chain_types
import cellranger.vdj.constants as vdj_constants
import tenkit.safe_json as tk_safe_json
from cellranger.reference_hash import compute_hash_of_file


class VDJReferenceConstructionError(Exception):
    """Raise for customer-facing errors in VDJ reference construction."""


class VdjAnnotationFeature(NamedTuple):
    feature_id: int
    record_id: bytes
    display_name: bytes
    gene_name: bytes
    region_type: bytes
    chain_type: bytes | None
    chain: bytes | None
    isotype: bytes | None
    allele_name: bytes | None
    sequence: bytes | None


class GtfEntry(NamedTuple):
    chrom: str
    source: str
    feature: str
    start: str
    end: str
    score: str
    strand: str
    frame: str
    attributes_str: str
    attributes: dict[str, str | int]


# Fields that are stored in the ref fasta header
REF_FASTA_FIELDS = ["feature_id", "display_name"]
REF_FASTA_AUX_FIELDS = [
    "record_id",
    "gene_name",
    "region_type",
    "chain_type",
    "chain",
    "isotype",
    "allele_name",
]


ENSEMBL_FIVE_PRIME_UTR_FEATURE = "five_prime_utr"
ENSEMBL_CDS_FEATURE = "CDS"
ENSEMBL_TRANSCRIPT_FEATURE = "transcript"

ENSEMBL_VDJ_BIOTYPES = {
    "TR_C_gene",
    "TR_D_gene",
    "TR_J_gene",
    "TR_V_gene",
    "IG_C_gene",
    "IG_D_gene",
    "IG_J_gene",
    "IG_V_gene",
}


def get_vdj_reference_fasta(reference_path):
    if reference_path is not None:
        return os.path.join(reference_path, vdj_constants.REFERENCE_FASTA_PATH)
    else:
        return "/dev/null"


def infer_ensembl_vdj_feature_type(feature, biotype):
    if feature == ENSEMBL_FIVE_PRIME_UTR_FEATURE:
        return "5'UTR"

    if feature == ENSEMBL_CDS_FEATURE:
        char = biotype.split("_")[1]

        if char == "V":
            # Note: ENSEMBL "V genes" include the L region
            return "L-REGION+V-REGION"
        if char == "D":
            return "D-REGION"
        if char == "J":
            return "J-REGION"
        if char == "C":
            return "C-REGION"

    return None


def standardize_ensembl_gene_name(gene_name: str) -> str | None:
    if gene_name is None:
        return gene_name

    # capitalize
    gene_name = gene_name.upper()
    # Rename TCR-* to TR*
    if gene_name.startswith(("TCRA-", "TCRB-", "TCRG-", "TCRD-")):
        gene_name = "TR" + gene_name[3] + gene_name[5:]
    return gene_name


def infer_ensembl_vdj_chain_type(gene_name):
    """Infer e.g., TR or IG from the ensembl gene name."""
    return bytes(gene_name[0:2], "utf-8")


def infer_ensembl_vdj_chain(gene_name):
    """Infer e.g., TRA or IGH from the ensembl gene name."""
    return gene_name[0:3]


def infer_ensembl_isotype(gene_name):
    """Infer e.g., E from IGHE."""
    if len(gene_name) <= 3:
        return None
    return gene_name[3:]


def get_duplicate_feature_key(f):
    return (f.display_name, f.region_type, f.sequence, f.chain_type, f.chain, f.isotype)


def build_reference_fasta_from_fasta(
    fasta_path, reference_path, reference_name, ref_version, mkref_version
):
    """Create cellranger-compatible vdj reference files from a.

    V(D)J segment FASTA file.
    """
    seen_features: set[tuple[bytes, bytes, bytes, bytes | None, bytes | None, bytes | None]] = set()

    seen_ids: set[int] = set()
    features: list[VdjAnnotationFeature] = []

    print("Checking FASTA entries...")

    with open(fasta_path, "rb") as f:
        for header, sequence in cr_utils.get_fasta_iter(f):
            feat = parse_fasta_entry(header, sequence)

            # Enforce unique feature IDs
            if feat.feature_id in seen_ids:
                raise VDJReferenceConstructionError(
                    "Duplicate feature ID found in input FASTA: %d." % feat.feature_id
                )
            # Sanity check values
            if b" " in feat.region_type:
                raise VDJReferenceConstructionError(
                    f'Spaces not allowed in region type: "{feat.region_type}"'
                )
            if b" " in feat.gene_name:
                raise VDJReferenceConstructionError(
                    f'Spaces not allowed in gene name: "{feat.gene_name}"'
                )
            if b" " in feat.record_id:
                raise VDJReferenceConstructionError(
                    f'Spaces not allowed in record ID: "{feat.record_id}"'
                )

            key = get_duplicate_feature_key(feat)
            if key in seen_features:
                print(
                    f"Warning: Skipping duplicate entry for {feat.display_name} ({feat.region_type}, {feat.record_id})."
                )
                continue

            # Strip Ns from termini
            seq = feat.sequence
            if b"N" in seq:
                print(
                    f"Warning: Feature {(feat.display_name, feat.record_id, feat.region_type)!s} contains Ns. Stripping from the ends."
                )
                seq = seq.strip(b"N")

            if len(seq) == 0:
                print(
                    f"Warning: Feature {(feat.display_name, feat.record_id, feat.region_type)!s} is all Ns. Skipping."
                )
                continue

            # Warn on features we couldn't classify properly
            if feat.chain_type not in chain_types.VDJ_CHAIN_TYPES:
                print(
                    "Warning: Unknown chain type for: {}. Expected name to be in {}. Skipping.".format(
                        str(
                            (
                                ensure_str(feat.display_name),
                                ensure_str(feat.record_id),
                                ensure_str(feat.region_type),
                            )
                        ),
                        str(tuple(ensure_str(t) for t in chain_types.VDJ_CHAIN_TYPES)),
                    )
                )
                continue

            seen_ids.add(feat.feature_id)
            seen_features.add(key)

            # Update the sequence since we may have modified it
            feat_dict = feat._asdict()
            feat_dict.update({"sequence": seq})
            new_feat = VdjAnnotationFeature(**feat_dict)
            features.append(new_feat)
    print("...done.\n")

    if len(seen_features) == 0:
        raise VDJReferenceConstructionError(
            "An empty constant regions file was "
            "generated/detected for your custom species. "
            "Please check if there are hits to the IMGT database "
            "or run cellranger vdj in denovo mode without reference."
        )

    print("Writing sequences...")
    os.makedirs(os.path.dirname(get_vdj_reference_fasta(reference_path)))
    with open(get_vdj_reference_fasta(reference_path), "w") as out_fasta:
        for feat in features:
            out_fasta.write(convert_vdj_feature_to_fasta_entry(feat) + "\n")
    print("...done.\n")

    print("Computing hash of input FASTA file...")
    fasta_hash = compute_hash_of_file(fasta_path)
    print("...done.\n")

    print("Writing metadata JSON file into reference folder...")
    metadata = {
        cr_constants.REFERENCE_GENOMES_KEY: reference_name,
        cr_constants.REFERENCE_FASTA_HASH_KEY: fasta_hash,
        cr_constants.REFERENCE_GTF_HASH_KEY: None,
        cr_constants.REFERENCE_INPUT_FASTA_KEY: os.path.basename(fasta_path),
        cr_constants.REFERENCE_INPUT_GTF_KEY: None,
        cr_constants.REFERENCE_VERSION_KEY: ref_version,
        cr_constants.REFERENCE_MKREF_VERSION_KEY: mkref_version,
        cr_constants.REFERENCE_TYPE_KEY: vdj_constants.REFERENCE_TYPE,
    }
    with open(os.path.join(reference_path, cr_constants.REFERENCE_METADATA_FILE), "w") as json_file:
        tk_safe_json.dump_numpy(metadata, json_file, sort_keys=True, indent=4)
    print("...done.\n")


def make_display_name(gene_name: str | None, allele_name: str | None) -> str | None:
    """Make a combined gene/allele name, e.g., TRAV1-1*01."""
    if allele_name is None:
        return gene_name
    else:
        return gene_name + "*" + allele_name


def convert_vdj_feature_to_fasta_entry(feature):
    """Generate a fasta entry from a VdjAnnotationFeature."""
    fasta_fields = [getattr(feature, f) for f in REF_FASTA_FIELDS]
    aux_fields = [getattr(feature, f) for f in REF_FASTA_AUX_FIELDS]
    hdr = (
        "|".join([x.decode() if hasattr(x, "decode") else str(x) for x in fasta_fields])
        + " "
        + "|".join([x.decode() if hasattr(x, "decode") else str(x) for x in aux_fields])
    )

    return f">{hdr}\n{ensure_str(feature.sequence)}"


def parse_fasta_entry(header: bytes, sequence: bytes) -> VdjAnnotationFeature:
    """Parse a FASTA entry into a VdjAnnotationFeature object."""
    words = header.split(b" ")

    # Check the header
    if len(words) != 2:
        raise VDJReferenceConstructionError(
            f'Expected two strings separated by a space in FASTA header. Found "{header}"'
        )

    values1 = words[0].split(b"|")
    if len(values1) != len(REF_FASTA_FIELDS):
        raise VDJReferenceConstructionError(
            "First string in FASTA header (record ID) must consist of the "
            'following %d fields separated by "|": %s. Found %d values: %s'
            % (
                len(REF_FASTA_FIELDS),
                ", ".join(REF_FASTA_FIELDS),
                len(values1),
                ensure_str(b", ".join(values1)),
            )
        )

    values2 = words[1].split(b"|")
    if len(values2) != len(REF_FASTA_AUX_FIELDS):
        raise VDJReferenceConstructionError(
            "Second string in FASTA header (description) must consist of the "
            'following %d fields separated by "|": %s. Found %d values: %s'
            % (
                len(REF_FASTA_AUX_FIELDS),
                ", ".join(REF_FASTA_AUX_FIELDS),
                len(values2),
                ensure_str(b", ".join(values2)),
            )
        )

    fields = {}
    fields.update(dict(zip(REF_FASTA_FIELDS, values1)))
    fields.update(dict(zip(REF_FASTA_AUX_FIELDS, values2)))

    # Validate the feature ID
    feature_id = fields.pop("feature_id")
    try:
        feature_id = int(feature_id)
        if feature_id < 1:
            raise ValueError()
    except ValueError:
        raise VDJReferenceConstructionError(
            f'The feature ID must be an integer greater than 0. Found: "{feature_id!s}"'
        )
    fields["sequence"] = sequence

    return VdjAnnotationFeature(feature_id=feature_id, **fields)


def get_feature_id_from_aligned_ref_name(ref_name):
    """Parse an aligned ref name (i.e. from a BAM file)."""
    # There should only be 1 word, but split just in case
    words = ref_name.split(" ")
    fields = dict(zip(REF_FASTA_FIELDS, words[0].split("|")))
    return int(fields["feature_id"])


def get_vdj_feature_iter(reference_path):
    """Yield vdj features from a vdj reference fasta file."""
    if reference_path is None:
        return

    with open(get_vdj_reference_fasta(reference_path), "rb") as reference_file:
        for header, sequence in cr_utils.get_fasta_iter(reference_file):
            yield parse_fasta_entry(header, sequence)


class VdjAnnotationFeatureDict(TypedDict):
    display_name: bytes
    feature_id: int
    chain: bytes | None
    gene_name: bytes
    region_type: bytes


def convert_vdj_feature_to_dict(feature: VdjAnnotationFeature) -> VdjAnnotationFeatureDict:
    """Yield a dict."""
    return {
        "display_name": feature.display_name,
        "feature_id": feature.feature_id,
        "chain": feature.chain,
        "gene_name": feature.gene_name,
        "region_type": feature.region_type,
    }


def convert_dict_to_vdj_feature(
    d: VdjAnnotationFeatureDict, reference: VdjReference
) -> VdjAnnotationFeature:
    """Convert a dict to a VdjAnnotationFeature."""
    return reference.get_feature_by_id(d["feature_id"])


def create_dummy_feature(display_name, region_type, sequence):
    """Create a "feature" that does not correspond to an actual reference segment."""
    display_name = ensure_binary(display_name)
    region_type = ensure_binary(region_type)
    return VdjAnnotationFeature(
        feature_id=0,
        record_id=b"",
        display_name=display_name,
        gene_name=display_name,
        region_type=region_type,
        chain_type=None,
        chain=None,
        isotype=None,
        allele_name=None,
        sequence=sequence,
    )


unannotated_feature = create_dummy_feature(b"UNANNOTATED", b"UNANNOTATED", None)


class VdjReference:
    """Represents a set of V(D)J reference sequences."""

    def __init__(self, reference_path):
        self.features: dict[int, VdjAnnotationFeature] = {}

        for feature in get_vdj_feature_iter(reference_path):
            self.features[feature.feature_id] = feature

    def get_feature_by_id(self, feature_id: int) -> VdjAnnotationFeature:
        if feature_id == 0:
            return unannotated_feature
        else:
            return self.features[feature_id]
