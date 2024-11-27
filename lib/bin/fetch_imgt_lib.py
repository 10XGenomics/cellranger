#!/usr/bin/env python3

"""Downloads IMGT sequences and creates a FASTA file suitable for input into `cellranger mkvdjref --seqs`.

Creates two files:
 <prefix>-imgt-raw.fasta
 <prefix>-mkvdjref-input.fasta
Where <prefix> is the string given via the --genome arg,
  *-imgt-raw is the IMGT segments translated to FASTA,
  and *-mkvdjref-input is in the format expected by `cellranger mkvdjref --seqs`

NOTE: Writes several downloaded HTML files to the current working dir.
"""

from __future__ import annotations

import argparse
import itertools
import os
import re
import ssl
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import OrderedDict
from collections.abc import Iterable
from dataclasses import dataclass
from html.parser import HTMLParser
from io import StringIO
from textwrap import dedent
from typing import IO

# certifi is available in the conda environment
import certifi

import cellranger.vdj.constants as vdj_constants
import cellranger.vdj.reference as vdj_reference

IMGT_GENEDB_URL = "https://www.imgt.org/genedb/GENElect"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--genome", help="V(D)J reference package name, e.g., my-vdj-ref", required=True
    )
    parser.add_argument(
        "--species",
        default="Homo sapiens",
        help="IMGT species to fetch; e.g. 'Homo sapiens' (note capitalization and quotes)",
    )
    args = parser.parse_args()

    queries = make_queries(args.species)

    try:
        filenames = download_files(args.species, queries)

    except urllib.error.URLError as ex:
        print(f"Failed to download from IMGT. {ex}\n")
        print("Failed to download all files from IMGT. Exiting.")
        sys.exit(1)

    fid = 0
    # Write IMGT fasta to a file
    with (
        open(args.genome + "-imgt-raw.fasta", "w") as raw_imgt_fa,
        open(args.genome + "-mkvdjref-input.fasta", "w") as mkvdjref_fa,
    ):
        for filename in filenames:
            with open(filename) as htmlfile:
                fa_txt = _GetLastPreTag()
                fa_txt.feed(htmlfile.read())
                fa_txt = fa_txt.last_pre

            raw_imgt_fa.write(fa_txt + "\n")

            f = StringIO(str(fa_txt))
            for record in fasta_parse(f):
                fid += 1
                feature = make_feature(record, fid)
                if feature is None:
                    continue

                mkvdjref_fa.write(vdj_reference.convert_vdj_feature_to_fasta_entry(feature) + "\n")


@dataclass
class Record:
    description: str
    seq: str


def fasta_parse(f: IO) -> Iterable[Record]:
    """Parse FASTA records from an input source."""
    current_header = None
    current_seq = []
    for line in f:
        line = line.strip()
        if not line:
            continue

        if line.startswith(">"):
            if current_header:
                description = current_header[1:]
                seq = "".join(current_seq)
                current_header = line
                current_seq = []
                yield Record(description, seq)
            else:
                current_header = line
        else:
            current_seq.append(line)

    if current_header:
        description = current_header[1:]
        seq = "".join(current_seq)
        yield Record(description, seq)


def make_queries(species):
    """IMGT GENE-DB queries."""
    v_genes = ["TRAV", "TRBV", "TRDV", "TRGV", "IGHV", "IGKV", "IGLV"]
    v_label = "L-PART1 V-EXON"  # let requests turn the ' ' into a '+'
    v_query = "8.1"

    d_genes = ["TRBD", "TRDD", "IGHD"]
    d_label = None
    d_query = "7.2"

    j_genes = ["TRAJ", "TRBJ", "TRDJ", "TRGJ", "IGHJ", "IGKJ", "IGLJ"]
    j_label = None
    j_query = "7.2"

    c_genes = ["TRAC", "TRBC", "TRDC", "TRGC", "IGHC"]
    c_label = None
    c_query = "14.1"

    c_genes2 = ["IGKC", "IGLC"]
    c_label2 = None
    if species == "Mus musculus":
        c_query2 = "7.2"
    else:
        c_query2 = "14.1"

    return list(
        itertools.chain(
            itertools.product(v_genes, [v_label], [v_query]),
            itertools.product(d_genes, [d_label], [d_query]),
            itertools.product(j_genes, [j_label], [j_query]),
            itertools.product(c_genes, [c_label], [c_query]),
            itertools.product(c_genes2, [c_label2], [c_query2]),
        )
    )


def download_files(species, queries):
    """Download the sequences.

    e.g. https://www.imgt.org/genedb/GENElect?query=8.1+TRBV&species=Homo+sapiens&IMGTlabel=L-PART1+V-EXON
    """
    filenames = []
    for gene, label, number in queries:
        filename = "_".join((species.replace(" ", ""), number, gene)) + ".html"
        filenames.append(filename)
        if os.path.exists(filename):
            print(f"Already downloaded {filename}, skipping")
            continue

        # Note: IMGT is sensitive to the param order
        payload = OrderedDict(
            query=f"{number} {gene}",
            species=species,
        )

        if label:
            payload["IMGTlabel"] = label

        # Create a SSL context with default certs (from system) as well
        # as certificates from certifi (possibly more frequently updated)
        context = ssl.create_default_context()
        context.load_verify_locations(certifi.where())

        encoded_args = urllib.parse.urlencode(payload)
        used_url = f"{IMGT_GENEDB_URL}?{encoded_args}"
        try:
            print(f"Downloading {used_url} to {filename} ...")
            with urllib.request.urlopen(used_url, context=context) as r:
                if r.status != 200:
                    raise urllib.error.URLError(used_url)

                with open(filename, "wb") as f:
                    f.write(r.read())
        except urllib.error.URLError as exc:
            if isinstance(exc.reason, ssl.SSLCertVerificationError):
                print(
                    dedent(
                        """\
                    Error: Failed to verify certificate, contact your system admin.
                    Updated certificates can be installed in Ubuntu/Debian by running:
                    $ apt install -y ca-certificates
                    and in Centos/RHEL using
                    $ yum install -y ca-certificates
                    """
                    )
                )
            raise exc

        # Don't hammer the server
        time.sleep(5)
    return filenames


def make_feature(record, fid):
    """Create a VdjAnnotationFeature from a record.

    Args:
        record (str): The record to parse.
        fid (int): The feature ID.

    Returns:
        vdj_reference.VdjAnnotationFeature: The feature.
    """
    row = record.description.split("|")

    region_type = get_region_type(row[4])
    if region_type is None:
        print(f"Warning: Unrecognized IMGT region type: {row[4]}; skipping...")
        return None

    chain_type = infer_imgt_vdj_chain_type(row[1])
    chain = infer_imgt_vdj_chain(row[1])

    if (
        region_type in vdj_constants.VDJ_C_FEATURE_TYPES
        and chain in vdj_constants.CHAINS_WITH_ISOTYPES
    ):
        isotype = infer_imgt_isotype(row[1])
    else:
        isotype = None

    gene_name = re.sub("[*].*$", "", row[1])
    allele_name = infer_imgt_allele(row[1])

    return vdj_reference.VdjAnnotationFeature(
        feature_id=fid,
        record_id=row[0],
        display_name=row[1],
        gene_name=gene_name,
        region_type=region_type,
        chain_type=chain_type,
        chain=chain,
        isotype=isotype,
        allele_name=allele_name,
        sequence=str(record.seq).upper(),
    )


class _GetLastPreTag(HTMLParser):  # pylint: disable=abstract-method
    def __init__(self):
        super().__init__()
        self._last_data = None
        self.last_pre = None

    def handle_data(self, data):
        self._last_data = data

    def handle_endtag(self, tag):
        if tag == "pre":
            self.last_pre = self._last_data


# Parse the HTML files
def get_region_type(imgt_label):
    """Convert IMGT labels into CR region type strings."""
    if imgt_label == "L-PART1+V-EXON":
        return "L-REGION+V-REGION"
    elif imgt_label in ("J-REGION", "D-REGION"):
        return imgt_label
    elif (
        "EX" in imgt_label
        or "CH" in imgt_label
        or "CL" in imgt_label
        or "M" in imgt_label
        or imgt_label == "C-REGION"
    ):
        return "C-REGION"
    else:
        return None


def infer_imgt_vdj_chain_type(gene_name):
    """Infer e.g., TR or IG from the IMGT gene name."""
    return gene_name[0:2]


def infer_imgt_vdj_chain(gene_name):
    """Infer e.g., TRA or IGH from the IMGT gene name."""
    return gene_name[0:3]


def infer_imgt_isotype(gene_name):
    """Infer, e.g., E from IGHE."""
    if len(gene_name) <= 3:
        return None
    return re.sub("[*].*$", "", gene_name)[3:]


def infer_imgt_allele(gene_name):
    return re.sub("^.*[*]", "", gene_name) or None


if __name__ == "__main__":
    main()
