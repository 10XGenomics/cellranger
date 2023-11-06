#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

from __future__ import annotations

import collections
import hashlib
import os
import shutil
from typing import AnyStr, NamedTuple, TypedDict

import pysam
from six import ensure_binary, ensure_str

import cellranger.constants as cr_constants
import cellranger.reference as cr_reference
import cellranger.utils as cr_utils
import cellranger.vdj.chain_types as chain_types
import cellranger.vdj.constants as vdj_constants
import tenkit.log_subprocess as tk_subproc
import tenkit.safe_json as tk_safe_json
import tenkit.seq as tk_seq


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


def get_gtf_iter(gtf_filename):
    gtf_file = open(gtf_filename)
    for line_num, line in enumerate(gtf_file):
        line = line.strip()
        if line.startswith("#"):
            continue

        row = line.split("\t")
        if len(row) != 9:
            raise ValueError(
                "Encountered malformed GTF at line %d. Expected 9 columns but found %d: %s"
                % (1 + line_num, len(row), line)
            )
        properties = cr_reference.NewGtfParser().get_properties_dict(
            row[8], line_num + 1, gtf_filename, uniquify_keys=True
        )
        yield GtfEntry(*row[:9], attributes=properties)


def get_duplicate_feature_key(f):
    return (f.display_name, f.region_type, f.sequence, f.chain_type, f.chain, f.isotype)


def _index_genome_fasta(reference_path, genome_fasta_path):
    # Note: We cannot symlink here because some filesystems in the wild
    #       do not support symlinks.
    print("Copying genome reference sequence...")
    os.makedirs(os.path.dirname(get_vdj_reference_fasta(reference_path)))
    tmp_genome_fa_path = os.path.join(reference_path, "genome.fasta")
    shutil.copy(genome_fasta_path, tmp_genome_fa_path)
    print("...done.\n")

    print("Indexing genome reference sequence...")
    tk_subproc.check_call(["samtools", "faidx", tmp_genome_fa_path])
    print("...done.\n")

    return tmp_genome_fa_path


def build_reference_fasta_from_ensembl(
    gtf_paths,
    transcripts_to_remove_path,
    genome_fasta_path,
    reference_path,
    reference_name,
    ref_version,
    mkref_version,
):
    """Create cellranger-compatible vdj reference files from a list of ENSEMBL-like GTF files.

    Input files are concatenated. No attempt to merge/reconcile information
    across them is made. Providing the files in a different order might change the
    output in cases where there are multiple entries with the same transcript id
    and the same feature type (eg. V-region).
    """
    transcripts: collections.defaultdict[tuple, list[GtfEntry]] = collections.defaultdict(list)

    if transcripts_to_remove_path:
        with open(transcripts_to_remove_path) as f:
            rm_transcripts = {line.strip() for line in f.readlines()}
    else:
        rm_transcripts = set()

    tmp_genome_fa_path = _index_genome_fasta(reference_path, genome_fasta_path)

    print("Loading genome reference sequence...")
    genome_fasta = pysam.FastaFile(tmp_genome_fa_path)
    print("...done.\n")

    print("Computing hash of genome FASTA file...")
    fasta_hash = cr_reference.compute_hash_of_file(tmp_genome_fa_path)
    print("...done.\n")

    for gtf in gtf_paths:
        print(f"Reading GTF {gtf}")

        for line_no, entry in enumerate(get_gtf_iter(gtf)):
            if not entry.feature in [ENSEMBL_FIVE_PRIME_UTR_FEATURE, ENSEMBL_CDS_FEATURE]:
                continue
            transcript_id: str = entry.attributes.get("transcript_id")
            transcript_biotype = entry.attributes.get("transcript_biotype")
            gene_biotype = entry.attributes.get("gene_biotype")
            gene_name = entry.attributes.get("gene_name")

            # Skip irrelevant biotypes
            if (
                transcript_biotype not in ENSEMBL_VDJ_BIOTYPES
                and not gene_biotype in ENSEMBL_VDJ_BIOTYPES
            ):
                continue

            # Skip blacklisted gene names
            if transcript_id in rm_transcripts:
                continue

            # Warn and skip if transcript_id missing
            if transcript_id is None:
                print("Warning: Entry on row %d has no transcript_id" % line_no)
                continue

            # Warn and skip if gene_name missing
            if gene_name is None:
                print(
                    "Warning: Transcript %s on row %d has biotype %s but no gene_name. Skipping."
                    % (transcript_id, line_no, transcript_biotype)
                )
                continue

            # Infer region type from biotype
            if transcript_biotype in ENSEMBL_VDJ_BIOTYPES:
                vdj_feature = infer_ensembl_vdj_feature_type(entry.feature, transcript_biotype)
            else:
                vdj_feature = infer_ensembl_vdj_feature_type(entry.feature, gene_biotype)

            # Warn and skip if region type could not be inferred
            if vdj_feature is None:
                print(
                    f"Warning: Transcript {transcript_id} has biotype {transcript_biotype}. Could not infer VDJ gene type. Skipping."
                )
                continue

            # Features that share a transcript_id and feature type are presumably exons
            # so keep them together.
            transcripts[(transcript_id, vdj_feature)].append(entry)

        print("...done.\n")

    print("Computing hash of genes GTF files...")
    digest = hashlib.sha1()
    # concatenate all the hashes into a string and then hash that string
    for gtf in gtf_paths:
        digest.update(cr_reference.compute_hash_of_file(gtf).encode())
    gtf_hash = digest.hexdigest()
    print("...done.\n")

    print("Fetching sequences...")
    out_fasta = open(get_vdj_reference_fasta(reference_path), "w")

    feature_id = 1
    seen_features = set()

    for (transcript_id, region_type), regions in transcripts.items():
        if not all(r.chrom == regions[0].chrom for r in regions):
            chroms: list[str] = list(sorted({r.chrom for r in regions}))
            print(
                f"Warning: Transcript {transcript_id} spans multiple contigs: {chroms!s}. Skipping."
            )
            continue

        if not all(r.strand == regions[0].strand for r in regions):
            print("Warning: Transcript %s spans multiple strands. Skipping." % transcript_id)
            continue

        chrom: str = regions[0].chrom
        strand: str = regions[0].strand
        ens_gene_name = standardize_ensembl_gene_name(regions[0].attributes["gene_name"])
        transcript_id: str = regions[0].attributes["transcript_id"]

        if chrom not in genome_fasta:
            print(
                f'Warning: Transcript {transcript_id} is on contig "{chrom}" which is not in the provided reference fasta. Skipping.'
            )
            continue

        # Build sequence
        regions.sort(key=lambda r: r.start)
        seq = b""
        for region in regions:
            # GTF coordinates are 1-based
            start, end = int(region.start) - 1, int(region.end)
            seq += ensure_binary(genome_fasta.fetch(chrom, start, end))

        # Revcomp if transcript on reverse strand
        if strand == "-":
            seq = tk_seq.get_rev_comp(seq)

        # Strip Ns from termini
        if b"N" in seq:
            print(
                "Warning: Feature %s contains Ns. Stripping from the ends."
                % str((ens_gene_name, transcript_id, region_type))
            )
            seq = seq.strip(b"N")

        if len(seq) == 0:
            print(
                "Warning: Feature %s is all Ns. Skipping."
                % str((ens_gene_name, transcript_id, region_type))
            )
            continue

        # Infer various attributes from the Ensembl gene name
        record_id = transcript_id
        gene_name = ens_gene_name
        display_name = make_display_name(gene_name=gene_name, allele_name=None)
        chain = infer_ensembl_vdj_chain(gene_name)
        chain_type = infer_ensembl_vdj_chain_type(gene_name)
        # Ensembl doesn't encode alleles
        allele_name = b"00"

        # Disallow spaces in these fields
        if " " in region_type:
            raise ValueError('Spaces not allowed in region type: "%s"' % region_type)
        if " " in gene_name:
            raise ValueError('Spaces not allowed in gene name: "%s"' % gene_name)
        if " " in record_id:
            raise ValueError('Spaces not allowed in record ID: "%s"' % record_id)

        # Warn on features we couldn't classify properly
        if chain_type not in chain_types.VDJ_CHAIN_TYPES:
            print(
                (
                    "Warning: Could not infer chain type for: %s. "
                    + "Expected the first two characters of the gene name to be in %s. Feature skipped."
                )
                % (
                    str((gene_name, record_id, region_type)),
                    str(tuple(chain_types.VDJ_CHAIN_TYPES)),
                )
            )
            continue

        if (
            region_type in vdj_constants.VDJ_C_FEATURE_TYPES
            and chain in vdj_constants.CHAINS_WITH_ISOTYPES
        ):
            isotype = infer_ensembl_isotype(ens_gene_name)
        else:
            isotype = None

        assert isinstance(chain_type, bytes)

        def _bytes_or_none(maybe_str: AnyStr | None) -> bytes | None:
            if maybe_str is None:
                return None
            return ensure_binary(maybe_str)

        gene_name = _bytes_or_none(gene_name)
        display_name = _bytes_or_none(display_name)
        region_type = _bytes_or_none(region_type)
        chain = _bytes_or_none(chain)
        isotype = _bytes_or_none(isotype)
        record_id = record_id.encode()
        feature = VdjAnnotationFeature(
            feature_id=feature_id,
            record_id=record_id,
            display_name=display_name,
            gene_name=gene_name,
            region_type=region_type,
            chain_type=chain_type,
            chain=chain,
            isotype=isotype,
            allele_name=allele_name,
            sequence=seq,
        )

        # Don't add duplicate entries
        feat_key = get_duplicate_feature_key(feature)
        if feat_key in seen_features:
            print(
                f"Warning: Skipping duplicate entry for {display_name} ({region_type}, {record_id})."
            )
            continue
        seen_features.add(feat_key)

        feature_id += 1

        out_fasta.write(convert_vdj_feature_to_fasta_entry(feature) + "\n")
    print("...done.\n")

    print("Deleting copy of genome fasta...")
    os.remove(tmp_genome_fa_path)
    os.remove(tmp_genome_fa_path + ".fai")
    print("...done.\n")

    print("Writing metadata JSON file into reference folder...")
    metadata = {
        cr_constants.REFERENCE_GENOMES_KEY: reference_name,
        cr_constants.REFERENCE_FASTA_HASH_KEY: fasta_hash,
        cr_constants.REFERENCE_GTF_HASH_KEY: gtf_hash,
        cr_constants.REFERENCE_INPUT_FASTA_KEY: os.path.basename(genome_fasta_path),
        cr_constants.REFERENCE_INPUT_GTF_KEY: ",".join(
            [os.path.basename(gtf_path) for gtf_path in gtf_paths]
        ),
        cr_constants.REFERENCE_VERSION_KEY: ref_version,
        cr_constants.REFERENCE_MKREF_VERSION_KEY: mkref_version,
        cr_constants.REFERENCE_TYPE_KEY: vdj_constants.REFERENCE_TYPE,
    }
    with open(os.path.join(reference_path, cr_constants.REFERENCE_METADATA_FILE), "w") as json_file:
        tk_safe_json.dump_numpy(metadata, json_file, sort_keys=True, indent=4)
    print("...done.\n")


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
                raise ValueError("Duplicate feature ID found in input FASTA: %d." % feat.feature_id)
            # Sanity check values
            if b" " in feat.region_type:
                raise ValueError('Spaces not allowed in region type: "%s"' % feat.region_type)
            if b" " in feat.gene_name:
                raise ValueError('Spaces not allowed in gene name: "%s"' % feat.gene_name)
            if b" " in feat.record_id:
                raise ValueError('Spaces not allowed in record ID: "%s"' % feat.record_id)

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
                    "Warning: Feature %s contains Ns. Stripping from the ends."
                    % str((feat.display_name, feat.record_id, feat.region_type))
                )
                seq = seq.strip(b"N")

            if len(seq) == 0:
                print(
                    "Warning: Feature %s is all Ns. Skipping."
                    % str((feat.display_name, feat.record_id, feat.region_type))
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

    print("Writing sequences...")
    os.makedirs(os.path.dirname(get_vdj_reference_fasta(reference_path)))
    with open(get_vdj_reference_fasta(reference_path), "w") as out_fasta:
        for feat in features:
            out_fasta.write(convert_vdj_feature_to_fasta_entry(feat) + "\n")
    print("...done.\n")

    print("Computing hash of input FASTA file...")
    fasta_hash = cr_reference.compute_hash_of_file(fasta_path)
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
        raise ValueError(
            'Expected two strings separated by a space in FASTA header. Found "%s"' % header
        )

    values1 = words[0].split(b"|")
    if len(values1) != len(REF_FASTA_FIELDS):
        raise ValueError(
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
        raise ValueError(
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
        raise ValueError(
            'The feature ID must be an integer greater than 0. Found: "%s"' % str(feature_id)
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
