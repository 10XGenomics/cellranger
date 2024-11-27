#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#

from __future__ import annotations

import collections
import csv
import itertools
import json
import os
import re
import subprocess
from typing import Any, NamedTuple

import numpy as np
from six import ensure_binary, ensure_str

import cellranger.constants as cr_constants
import cellranger.cr_io as cr_io
import tenkit.seq as tk_seq


class Interval(NamedTuple):
    chrom: str | bytes  # TODO: be consistent
    start: int
    end: int
    length: int
    strand: str | None


class Gene(NamedTuple):
    id: str | bytes  # TODO: be consistent
    name: str
    length: int | None
    gc_content: float | None
    intervals: list[Interval] | None


class Transcript(NamedTuple):
    gene: Gene
    length: int | None
    gc_content: float | None
    intervals: list[Interval]


class GtfParser:
    """This is an old-style GtfParser class.

    It is used so that we can produce backwards compatible references that work
    with older releases of CR.

    For use within CR, use classes derived from `NewGtfParser`, below.

    In order for the `genes.pickle` file to unpickle in older CR builds, it must be an old-style class
    named `GtfParser`.
    """

    def __init__(self, skip_exception=False):
        # don't call this unless you're dealing with pickle backward compat issues
        if not skip_exception:
            raise Exception(
                "GtfParser is deprecated & is only present for backwards compatibility reasons. Use NewGtfParser"
            )


class GtfParseError(Exception):
    """Exception class for errors while parsing a GTF file."""

    def __init__(self, filename, msg):
        super().__init__(
            f"Error while parsing GTF file {filename}\n{msg}\n\nPlease fix your GTF and start again."
        )


class BaseReferenceError(Exception):
    """Base class for all reference construction errors."""


class GexReferenceError(BaseReferenceError):
    """Exception class for errors while constructing a GEX reference."""


class NewGtfParser:
    # Precompiled RegEx to get attributes in GTF Records
    _attribute_pattern = re.compile(r'(\S+?)\s+(".*?"|[^";\n\r]+)\s*;{0,1}')

    def gtf_reader_iter(
        self,
        filename: os.PathLike | str | bytes,
        contig_lengths: dict[str, int] | None = None,
        no_transcript_fail: bool = False,
    ):
        """Return an iteator over the rows of a GTF.

        The iterator consists of
        (row, is_comment, properties) tuples of type (str, bool, dict).
        In case an entry's properties consist of multiple key-value pairs with
        the same key, the duplicate keys are made unique. To properly translate
        the properties dict back into GTF format, use format_properties_dict.

        Args:
            contig_lengths: dictionary of {contig: length}. If specified we check
                that start/end of gtf are consistent with contig lengths
            no_transcript_fail: should we accept a gtf that has no 'transcript' row, but
                has 'exon' rows that reference a transcript? Default is False
                because of backwards compatibility with cellranger. For
                ARC we exit since this is not valid gtf.
        """
        transcript_id = "transcript_id"
        with cr_io.open_maybe_gzip(filename, "r") as f:
            reader = csv.reader(f, delimiter="\t")
            transcript_row_id = None
            for i, row in enumerate(reader):
                if len(row) == 0:
                    continue
                if row[0].startswith("#"):
                    yield row, True, None
                    continue

                if len(row) != 9:
                    raise GtfParseError(
                        filename,
                        "Invalid number of columns (%d, expect 9) in GTF line %d: %s"
                        % (len(row), i + 1, "\t".join(row)),
                    )

                start = int(row[3])
                end = int(row[4])
                if start > end:
                    raise GtfParseError(
                        filename,
                        f"Invalid GTF annotation on line {i + 1}:\n"
                        f"{row}\n"
                        f"Start position of feature = {start} > End position of feature = {end}",
                    )
                if contig_lengths:
                    if row[0] not in contig_lengths:
                        contigs = sorted(list(contig_lengths.keys()))
                        raise GtfParseError(
                            filename,
                            f"Invalid contig name encountered on GTF line {i + 1}: {row[0]}. The FASTA file "
                            f"has contigs:\n{contigs}",
                        )

                    max_len = contig_lengths[row[0]]
                    if end > max_len:
                        raise GtfParseError(
                            filename,
                            f"Invalid GTF annotation on line {i + 1}:\n"
                            f"{row}\n"
                            f"End position of feature = {end} > contig {row[0]} length = {max_len}",
                        )
                    if start < 1:
                        raise GtfParseError(
                            filename,
                            f"Invalid GTF annotation on line {i + 1}:\n"
                            f"{row}\n"
                            f"Start position of feature = {start} < 1",
                        )

                strand = row[6].encode()
                if strand not in cr_constants.STRANDS:
                    raise GtfParseError(
                        filename,
                        "Invalid strand in GTF line %d: %s" % (i + 1, "\t".join(row)),
                    )

                annotation = row[2]
                properties = self.get_properties_dict(row[8], i + 1, filename, uniquify_keys=True)
                if annotation == "transcript":
                    self.validate_transcript_id(filename, properties, i, row)
                    self.validate_gene_id(filename, properties, i, row)
                    transcript_row_id = properties[transcript_id]
                elif annotation == "exon":
                    self.validate_transcript_id(filename, properties, i, row)
                    self.validate_gene_id(filename, properties, i, row)
                    tid = properties[transcript_id]
                    if no_transcript_fail and tid != transcript_row_id:
                        raise GtfParseError(
                            filename,
                            f"Supplied GTF is invalid. This row\n{row}\n"
                            f"on line {i + 1} specifies an 'exon' annotation for a "
                            f"transcript {tid}, but there is no 'transcript' row in "
                            f"the GTF for {tid} that immediately precedes it.",
                        )
                elif annotation == "gene":
                    # NCBI for some organisms releases genes with `transcript_id "";` and
                    # that can trip up later code because genes shouldn't have transcript ids
                    # and it looks like the same transcript is on multiple chromosomes
                    tid = properties.get(transcript_id)
                    if tid is not None and tid == "":
                        del properties[transcript_id]

                yield row, False, properties

    @staticmethod
    def validate_transcript_id(filename, properties, i, row):
        if "transcript_id" not in properties:
            raise GtfParseError(
                filename,
                "Property 'transcript_id' not found in GTF line %d: %s" % (i + 1, "\t".join(row)),
            )
        if properties["transcript_id"] == "":
            raise GtfParseError(
                filename,
                "Property 'transcript_id' is empty in GTF line {:d}: {}".format(
                    i + 1, "\t".join(row)
                ),
            )
        if ";" in properties["transcript_id"]:
            raise GtfParseError(
                filename,
                "Property 'transcript_id' has invalid character ';' in GTF line %d: %s"
                % (i + 1, "\t".join(row)),
            )
        if re.search(r"\s", properties["transcript_id"]) is not None:
            raise GtfParseError(
                filename,
                "Property 'transcript_id' has invalid whitespace character in GTF line %d: %s"
                % (i + 1, "\t".join(row)),
            )

    @staticmethod
    def validate_gene_id(filename, properties: dict[str, str | int], i, row):
        if "gene_id" not in properties:
            raise GtfParseError(
                filename,
                "Property 'gene_id' not found in GTF line %d: %s" % (i + 1, "\t".join(row)),
            )
        if properties["gene_id"] == "":
            raise GtfParseError(
                filename,
                "Property 'gene_id' is empty in GTF line {:d}: {}".format(i + 1, "\t".join(row)),
            )
        if ";" in properties["gene_id"]:
            raise GtfParseError(
                filename,
                "Property 'gene_id' has invalid character ';' in GTF line %d: %s"
                % (i + 1, "\t".join(row)),
            )

    def load_gtf(
        self, in_gtf_fn: os.PathLike | str | bytes, fasta_parser: FastaParser | None = None
    ) -> tuple[dict[str, Transcript], list[Gene]]:
        transcripts: dict[str, Transcript] = {}
        gene_to_transcripts: dict[Gene, set[str]] = collections.OrderedDict()

        for row, is_comment, properties in self.gtf_reader_iter(in_gtf_fn):
            if is_comment:
                continue

            chrom, _, annotation, start, end, _, strand, _, _ = row

            if annotation != "exon":
                continue

            start = int(start) - 1
            end = int(end)
            length = abs(end - start)
            transcript_id: str = properties["transcript_id"]
            gene_id: str = properties["gene_id"]
            gene_name: str = properties.get("gene_name", gene_id)
            gene = Gene(gene_id, gene_name, None, None, None)

            if transcript_id not in transcripts:
                transcripts[transcript_id] = Transcript(gene, None, None, [])

            if gene not in gene_to_transcripts:
                gene_to_transcripts[gene] = set()

            assert transcripts[transcript_id].gene == gene
            transcripts[transcript_id].intervals.append(
                Interval(ensure_binary(chrom), start, end, length, strand)
            )
            gene_to_transcripts[gene].add(transcript_id)

        # Transcript length and GC content
        transcript_lengths: dict[str, int] = {}
        transcript_gc_contents: dict[str, float] = {}
        for transcript_id, transcript in transcripts.items():
            transcript_lengths[transcript_id] = sum(
                interval.length for interval in transcript.intervals
            )
            if fasta_parser is not None:
                transcript_gc_contents[transcript_id] = fasta_parser.get_transcript_gc_content(
                    transcript
                )

        # Gene length, GC content and start + end positions
        genes: list[Gene] = []
        for gene, transcript_ids in gene_to_transcripts.items():
            length: int = np.median(
                [transcript_lengths[transcript_id] for transcript_id in transcript_ids]
            )
            gc_content: float = (
                np.median(
                    [transcript_gc_contents[transcript_id] for transcript_id in transcript_ids]
                )
                if fasta_parser is not None
                else float("nan")
            )

            transcript_intervals = []
            for transcript_id in transcript_ids:
                transcript_intervals += transcripts[transcript_id].intervals
            transcript_intervals.sort(key=lambda interval: interval.chrom)

            intervals: list[Interval] = []
            for chrom, chrom_intervals_iter in itertools.groupby(
                transcript_intervals, lambda interval: interval.chrom
            ):
                chrom_intervals = list(chrom_intervals_iter)
                start = min(interval.start for interval in chrom_intervals)
                end = max(interval.end for interval in chrom_intervals)
                interval = Interval(chrom, start, end, end - start, None)
                intervals.append(interval)

            gene = Gene(gene.id, gene.name, length, gc_content, intervals)
            genes.append(gene)

            for transcript_id in transcript_ids:
                transcripts[transcript_id] = Transcript(
                    gene,
                    transcript_lengths[transcript_id],
                    transcript_gc_contents.get(transcript_id, float("nan")),
                    transcripts[transcript_id].intervals,
                )

        return transcripts, genes

    @staticmethod
    def get_properties_dict(
        properties_str: str, line_number: int, filename: str, uniquify_keys: bool = True
    ) -> dict[str, int | str]:
        """Parse the properties present in the 9th column of a gtf entry into a dictionary.

        If there are multiple properties with the same key, those keys will be
        made unique unless uniquify_keys is False, in which case only the final
        instance of a key will appear in the result.
        """
        if isinstance(properties_str, dict):
            return properties_str

        properties: collections.OrderedDict[str, int | str] = collections.OrderedDict()

        if uniquify_keys:
            # number of "keys" that are duplicates of ones previously seen
            repeat_num = 0

        for m in NewGtfParser._attribute_pattern.finditer(properties_str):
            # extract key-value pair
            key = m.group(1)
            value = m.group(2)
            if not key or not value:
                continue
            # parse value
            start_quote = value[0] == '"'
            end_quote = value[-1] == '"'
            if start_quote and end_quote:
                # quoted string
                value = value[1:-1]
            elif start_quote:
                raise GtfParseError(
                    filename,
                    "Error parsing GTF at line %d.  Parsed attribute began, but did not end with a quote.  Please ensure attributes that start with quotes end with them.\n Bad Attribute = %s"
                    % (line_number, value),
                )
            elif end_quote:
                raise GtfParseError(
                    filename,
                    "Error parsing GTF at line %d.  Parsed attribute ended, but did not begin with a quote.  Please ensure attributes that end with quotes end start with them.\n Bad Attribute = %s"
                    % (line_number, value),
                )
            elif value.isdigit():
                # unquoted integer
                value = int(value)
            elif '"' in key or '"' in value:
                # Leave it mostly as is (technically not against the GTF format,
                # at least as per http://mblab.wustl.edu/GTF22.html)
                # but make sure we don't have any quotes in the value, this will bjork the rust parser
                raise GtfParseError(
                    filename,
                    "Error parsing GTF at line %d.  Parsed attribute had a quote in the middle of a value.  Please ensure quotes are only used to encapsulate attribute values.\n Bad Attribute Value = %s"
                    % (line_number, value),
                )

            if '"' in key:
                raise GtfParseError(
                    filename,
                    "Error parsing GTF at line %d.  Parsed attribute had a quote in the middle of a key.  Please ensure quotes are only used to encapsulate attribute values.\n Bad Attribute Key = %s"
                    % (line_number, value),
                )
            # make key unique if necessary
            if uniquify_keys and key in properties:
                repeat_num += 1
                key = NewGtfParser._make_key_unique(key, repeat_num)

            # add a valid key-value pair to properties
            properties[key] = value
        return properties

    @staticmethod
    def format_properties_dict(properties, uniquify_keys=True):
        """Translate a properties dict into a GTF-formatted string.

        Note:
            Keys beginning with "GtfParser_" are assumed to come from
            `_make_key_unique` and will be translated. This behavior can be
            disabled by setting `uniquify_keys=False`.
        """
        properties_str = []
        for key, value in properties.items():
            # convert value to string
            if isinstance(value, int):
                value = str(value)
            else:
                value = f'"{value}"'
            # translate uniquified key if necessary
            if uniquify_keys:
                key = NewGtfParser._translate_unique_key(key)
            properties_str.append(f"{key} {value}")
        return "; ".join(properties_str) + ";"

    @staticmethod
    def _make_key_unique(key: str, num: int) -> str:
        return f"GtfParser_{key}_{num}"

    @staticmethod
    def _translate_unique_key(key: str) -> str:
        if key.startswith("GtfParser_"):
            return key[10 : key.rindex("_")]
        else:
            return key


class GtfBuilder(NewGtfParser):
    def __init__(self, in_gtf_fn, out_gtf_fn, attributes=()):
        self.in_gtf_fn = in_gtf_fn
        self.out_gtf_fn = out_gtf_fn
        self.attributes = attributes or {}

    def build_gtf(self):
        print("Writing new genes GTF file (may take 10 minutes for a 1GB input GTF file)...")
        with open(self.out_gtf_fn, "w") as f:
            writer = csv.writer(
                f, delimiter="\t", quoting=csv.QUOTE_NONE, quotechar="", lineterminator="\n"
            )
            for row, is_comment, properties in self.gtf_reader_iter(self.in_gtf_fn):
                if is_comment:
                    writer.writerow(row)
                    continue

                remove = False
                assert properties is not None
                for key, value in properties.items():
                    if key in self.attributes and value not in self.attributes[key]:
                        remove = True

                if not remove:
                    writer.writerow(row)

        print("...done\n")


# Find by relative path in lib/bin.
_GTF_TO_GENE_INDEX = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__.encode())))),
    b"bin",
    b"gtf_to_gene_index",
)


class FastaParser:
    def __init__(self, in_fasta_fn: os.PathLike | str | bytes):
        self.chroms = self.load_fasta(in_fasta_fn)

    @staticmethod
    def _get_chrom_name(fasta_header: bytes) -> bytes:
        match = re.search(rb">(\S+)", fasta_header)
        return match.groups(1)[0]

    @staticmethod
    def load_fasta(in_fasta_fn: os.PathLike | str | bytes) -> dict[bytes, bytes]:
        chroms = {}

        current_chrom, current_seq = None, None
        with open(in_fasta_fn, "rb") as f:
            for line in f:
                line = line.strip()
                if line.startswith(b">"):
                    if current_chrom is not None:
                        assert current_seq is not None
                        chroms[current_chrom] = b"".join(current_seq)
                    current_chrom = FastaParser._get_chrom_name(line)
                    current_seq = []
                else:
                    current_seq.append(line)

        if current_chrom is not None:
            assert current_seq is not None
            chroms[current_chrom] = b"".join(current_seq)

        return chroms

    def is_valid_interval(self, chrom: bytes, start: int, end: int):
        """Determine whether the half-open interval [start,end) is within the bounds of chrom."""
        assert isinstance(chrom, bytes)
        return chrom in self.chroms and start >= 0 and end <= len(self.chroms[chrom])

    def get_sequence(
        self, chrom: bytes, start: int, end: int, strand: bytes = cr_constants.FORWARD_STRAND
    ):
        """Get genomic sequence for the half-open interval [start,end)."""
        assert isinstance(chrom, bytes)
        seq = self.chroms[chrom][start:end]
        if strand == cr_constants.FORWARD_STRAND:
            return seq
        elif strand == cr_constants.REVERSE_STRAND:
            return tk_seq.get_rev_comp(seq)
        else:
            raise Exception(f"Invalid strand: {strand}")

    def get_transcript_gc_content(self, transcript_obj: Transcript) -> float:
        pattern = re.compile(b"[cCgG]")

        gc, length = 0, 0
        for interval in transcript_obj.intervals:
            assert isinstance(interval.chrom, bytes)
            if interval.chrom not in self.chroms:
                continue

            seq = self.chroms[interval.chrom][interval.start : interval.end]
            gc += len(re.findall(pattern, seq))
            length += interval.length

        if length > 0:
            return float(gc) / float(length)
        else:
            return 0


# NOTE: these stub classes are necessary to maintain backwards compatibility with old refdata (1.2 or older)
class IntervalTree:
    pass


class Region:
    pass


class IntergenicRegion(Region):
    pass


class IntronicRegion(Region):
    pass


class ExonicRegion(Region):
    pass


def map_path(module: str, name: str):
    """Choose the appropriate type for loading a given file type.

    When unpickling an object with type named `cellranger.reference.GeneIndex` in a
    pickle file, actually open it using `cellranger.reference.NewGeneIndex`.

    Otherwise dynamically load the requested module.
    """
    if module == "cellranger.reference" and name == "GeneIndex":
        return NewGeneIndex

    if module == "cellranger.reference" and name == "GtfParser":
        return NewGtfParser

    mod = __import__(module, fromlist=["__name__"])
    print(f"module: {module}, mod: {mod!s}")
    return getattr(mod, name)


class GeneIndex(GtfParser):
    def __init__(self, in_gtf_fn, in_fasta_fn, genes, transcripts, gene_ids_map):
        """Constructor.

        This should not be called directly!! Use NewGeneIndex for new usages.
        """
        self.in_gtf_fn = in_gtf_fn
        self.in_fasta_fn = in_fasta_fn
        self.transcripts = transcripts
        self.genes = genes
        self.gene_ids_map = gene_ids_map
        GtfParser.__init__(self, True)


class NewGeneIndex(NewGtfParser):
    def __init__(
        self,
        in_gtf_fn: os.PathLike | str | bytes | None = None,
        in_fasta_fn: os.PathLike | str | bytes | None = None,
    ):
        """Constructor.

        in_gtf_fn and in_fasta_fn are required, but are defaulted for backwards
        compatibility with pickled data.
        """
        if in_gtf_fn is None and in_fasta_fn is None:
            # loading from old-version pickle, which will populate the rest.
            return
        assert in_gtf_fn is not None
        assert in_fasta_fn is not None
        self.in_gtf_fn = in_gtf_fn
        self.in_fasta_fn = in_fasta_fn

        fasta_parser = FastaParser(self.in_fasta_fn)

        self.transcripts, self.genes = self.load_gtf(self.in_gtf_fn, fasta_parser=fasta_parser)
        self.gene_ids_map = {gene.id: i for i, gene in enumerate(self.genes)}
        # super(NewGeneIndex, self).__init__()

    @staticmethod
    def load_from_reference(reference_path: str | bytes):
        """Load the GeneIndex directly from the reference path.

        Directly parses the GTF and FASTA files using helper Rust tool.
        Does not require the gene.pickle file.
        """
        gene_index_fn = b"gene_index.json"

        cmd = [_GTF_TO_GENE_INDEX, ensure_binary(reference_path), gene_index_fn]
        try:
            subprocess.check_output(cmd, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as exc:
            raise RuntimeError("couldn't load reference info: " + ensure_str(exc.output)) from exc

        try:
            return NewGeneIndex.load_from_json(gene_index_fn)

        finally:
            # Delete the JSON file
            if os.path.exists(gene_index_fn):
                os.unlink(gene_index_fn)

    @staticmethod
    def load_from_json(gene_index_fn: str | bytes | os.PathLike) -> NewGeneIndex:
        """Load the GeneIndex directly from the reference path.

        Directly parses the GTF and FASTA files using helper Rust tool.
        Does not require the gene.pickle file.
        """

        def conv_gene(obj: dict[str, Any]) -> Gene:
            intervals = [conv_interval(x) for x in obj["intervals"]]
            return Gene(
                str(obj["id"]).encode(),
                str(obj["name"]),
                obj["length"],
                obj["gc_content"],
                intervals,
            )

        def conv_interval(obj: dict[str, str | int | None]) -> Interval:
            if obj["strand"]:
                s = str(obj["strand"])
            else:
                s = None
            return Interval(str(obj["chrom"]).encode(), obj["start"], obj["end"], obj["length"], s)

        def conv_transcript(obj: dict[str, Any]) -> Transcript:
            gene = conv_gene(obj["gene"])
            intervals = [conv_interval(x) for x in obj["intervals"]]
            if intervals[0].strand == "-":
                intervals.reverse()

            return Transcript(gene, obj["length"], obj["gc_content"], intervals)

        def conv_transcripts(obj: dict[str, dict[str, Any]]) -> dict[str, Transcript]:
            res = {}
            for tx_id, tx in obj.items():
                new_tx = conv_transcript(tx)
                res[str(tx_id)] = new_tx
            return res

        with open(gene_index_fn) as f:
            data = json.load(f)

        # Convert JSON into the GeneIndex namedtuples
        s = NewGeneIndex()
        s.in_gtf_fn = data["in_gtf_fn"]
        s.in_fasta_fn = data["in_fasta_fn"]
        s.transcripts = conv_transcripts(data["transcripts"])
        s.genes = [conv_gene(x) for x in data["genes"]]
        s.gene_ids_map = {gene.id: i for i, gene in enumerate(s.genes)}
        return s

    def get_transcript_length(self, transcript: str):
        assert isinstance(transcript, str)
        if transcript in self.transcripts:
            return self.transcripts[transcript].length
        return None

    def get_transcript_gc_content(self, transcript: str):
        assert isinstance(transcript, str)
        if transcript in self.transcripts:
            return self.transcripts[transcript].gc_content
        return None

    def get_gene_from_transcript(self, transcript: str):
        assert isinstance(transcript, str)
        if transcript in self.transcripts:
            return self.transcripts[transcript].gene
        return None

    def gene_id_to_int(self, gene_id: bytes) -> int | None:
        """Return the integer index of gene_id if found or None otherwise."""
        assert isinstance(gene_id, bytes)
        if gene_id in self.gene_ids_map:
            return self.gene_ids_map[gene_id]
        return None

    def get_genes(self):
        return self.genes

    def get_gene(self, gene_id: bytes) -> Gene | None:
        """Return the Gene object whose gene ID is gene_id or None otherwise."""
        assert isinstance(gene_id, bytes)
        index = self.gene_id_to_int(gene_id)
        return self.genes[index] if index is not None and index < len(self.genes) else None

    def get_gene_lengths(self) -> list[int]:
        return [gene.length for gene in self.genes]

    def get_gene_gc_contents(self):
        return [gene.gc_content for gene in self.genes]

    def get_gene_names(self):
        return [x.name for x in self.genes]

    def get_gene_ids(self):
        return [x.id for x in self.genes]

    def subset_index(self, target_genes):
        """Reduce this GeneIndex to only genes in the set target_genes."""
        self.genes = [g for g in self.genes if g.id in target_genes]
        txs = {tx_id: tx for tx_id, tx in self.transcripts.items() if tx.gene.id in target_genes}
        self.transcripts = txs
        self.gene_ids_map = {gene.id: i for i, gene in enumerate(self.genes)}
