#
# Copyright (c) 2023 10X Genomics, Inc. All rights reserved.
#

"""Code to make a VDJ reference."""
from __future__ import annotations

import collections
import hashlib
import os
import shutil
from typing import AnyStr

import pysam
from six import ensure_binary

from cellranger import constants as cr_constants
from cellranger import reference as cr_reference
from cellranger.reference_hash import compute_hash_of_file
from cellranger.vdj import chain_types
from cellranger.vdj import constants as vdj_constants
from cellranger.vdj.reference import (
    ENSEMBL_CDS_FEATURE,
    ENSEMBL_FIVE_PRIME_UTR_FEATURE,
    ENSEMBL_VDJ_BIOTYPES,
    GtfEntry,
    VdjAnnotationFeature,
    VDJReferenceConstructionError,
    convert_vdj_feature_to_fasta_entry,
    get_duplicate_feature_key,
    get_vdj_reference_fasta,
    infer_ensembl_isotype,
    infer_ensembl_vdj_chain,
    infer_ensembl_vdj_chain_type,
    infer_ensembl_vdj_feature_type,
    make_display_name,
    standardize_ensembl_gene_name,
)
from tenkit import log_subprocess as tk_subproc
from tenkit import safe_json as tk_safe_json
from tenkit import seq as tk_seq


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
    fasta_hash = compute_hash_of_file(tmp_genome_fa_path)
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
        digest.update(compute_hash_of_file(gtf).encode())
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
            print(f"Warning: Transcript {transcript_id} spans multiple strands. Skipping.")
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
                f"Warning: Feature {(ens_gene_name, transcript_id, region_type)!s} contains Ns. Stripping from the ends."
            )
            seq = seq.strip(b"N")

        if len(seq) == 0:
            print(
                f"Warning: Feature {(ens_gene_name, transcript_id, region_type)!s} is all Ns. Skipping."
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
            raise VDJReferenceConstructionError(
                f'Spaces not allowed in region type: "{region_type}"'
            )
        if " " in gene_name:
            raise VDJReferenceConstructionError(f'Spaces not allowed in gene name: "{gene_name}"')
        if " " in record_id:
            raise VDJReferenceConstructionError(f'Spaces not allowed in record ID: "{record_id}"')

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

    if len(seen_features) == 0:
        raise VDJReferenceConstructionError(
            "An empty constant regions file was "
            "generated/detected for your custom species. "
            "Please check if there are hits to the IMGT database "
            "or run cellranger vdj in denovo mode without reference."
        )

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


def get_gtf_iter(gtf_filename):
    gtf_file = open(gtf_filename)
    for line_num, line in enumerate(gtf_file):
        line = line.strip()
        if line.startswith("#"):
            continue

        row = line.split("\t")
        if len(row) != 9:
            raise VDJReferenceConstructionError(
                "Encountered malformed GTF at line %d. Expected 9 columns but found %d: %s"
                % (1 + line_num, len(row), line)
            )
        properties = cr_reference.NewGtfParser().get_properties_dict(
            row[8], line_num + 1, gtf_filename, uniquify_keys=True
        )
        yield GtfEntry(*row[:9], attributes=properties)
