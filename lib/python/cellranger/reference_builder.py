#!/usr/bin/env python3
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
from __future__ import annotations

import csv
import math
import os
import shutil
import subprocess
import sys
import tempfile

from six import ensure_str

import cellranger.constants as cr_constants
import tenkit.log_subprocess as tk_subproc
import tenkit.safe_json as tk_safe_json
from cellranger.reference import (
    GexReferenceError,
    GtfParseError,
    NewGtfParser,
)
from cellranger.reference_hash import compute_hash_of_file

_LIB_BIN = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "bin")

# Find by relative path in lib/bin.
_GTF_TO_GENE_INDEX = os.path.join(
    _LIB_BIN.encode(),
    b"gtf_to_gene_index",
)


class DuplicateContigNameException(Exception):
    """Used to indicate the input file has multiple contigs with the same name."""


class ReferenceBuilder(NewGtfParser):
    def __init__(
        self,
        genomes: list[str],
        in_fasta_fns,
        in_gtf_fns,
        out_dir: str | os.PathLike | None,
        ref_version,
        mkref_version,
        num_threads: int = 1,
        mem_gb: int = 4,
        prefixes_as_genomes=False,
    ):
        self.genomes = genomes
        # check that one genome is not a prefix of another. This is slightly stronger than what we
        # require, but this is reasonable to demand.
        for i, g1 in enumerate(self.genomes):
            for j, g2 in enumerate(self.genomes):
                if i != j and g1.startswith(g2):
                    raise GexReferenceError(
                        f"Supplied genome name '{g2}' is a prefix of genome name '{g1}'"
                    )

        self.in_fasta_fns = in_fasta_fns
        self.in_gtf_fns = in_gtf_fns
        self.out_dir = out_dir
        self.num_threads = num_threads
        self.mem_gb = mem_gb
        self.ref_version = ref_version
        self.mkref_version = mkref_version
        # path to final genome.fa file
        self.fasta_path = None
        # path to final genes.gtf file
        self.gtf_path: str | None = None
        self.genome_prefixes: list[str] | None = None
        self.format_genome_prefixes(prefixes_as_genomes)

    def process_fasta(self):
        """Copy FASTA to fasta/genome.fa and samtools index."""
        print("Writing genome FASTA file into reference folder...")
        self.fasta_path = os.path.join(self.out_dir, cr_constants.REFERENCE_FASTA_PATH)
        os.mkdir(os.path.dirname(self.fasta_path))
        self.write_genome_fasta(self.fasta_path)
        print("...done\n")

        print("Indexing genome FASTA file...")
        # from lib/bin
        samtools = os.path.join(
            os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "bin", "samtools"
        )
        stdout = tk_subproc.check_output(
            [samtools, "faidx", self.fasta_path], stderr=subprocess.STDOUT
        )
        output = stdout.decode()
        # We want to check if samtools told us there'd be a duplicate sequence, as this will cause sigfaults
        # downstream, see CSR-1202
        if "Ignoring duplicate sequence" in output:
            error_msg = (
                "mkref failed - your fasta file has duplicated sequence names in it."
                "Please give every contig a unique name."
                "\nsamtools faidx output is: " + output
            )
            raise DuplicateContigNameException(error_msg)
        print("...done\n")

    def process_gtf(
        self, contig_lengths: dict[str, int] | None = None, no_transcript_fail: bool = False
    ):
        """Parse and write GTF and compute gtf pickle."""
        assert self.fasta_path, "FASTA path must be provided"
        assert os.path.exists(self.fasta_path), "FASTA file not present in fasta/genome.fa"
        print("Writing genes GTF file into reference folder...")

        # write a gzip copy of the GTF.
        # note old references have a non-gzipped GTF, so cr_constants.REFERENCE_GENES_GTF_PATH
        # does not contain the ".gz" suffix. We will gzip the gtf once STAR has used it.
        assert self.out_dir is not None
        gtf_path = os.path.join(ensure_str(self.out_dir), cr_constants.REFERENCE_GENES_GTF_PATH)
        self.gtf_path = gtf_path
        os.mkdir(os.path.dirname(gtf_path))
        has_warning = self.write_genome_gtf(
            gtf_path, contig_lengths=contig_lengths, no_transcript_fail=no_transcript_fail
        )
        print("...done\n")
        return has_warning

    def gzip_gtf(self):
        # From conda, should be a sibling of the interpreter.
        bgzip = os.path.join(os.path.dirname(sys.executable), "bgzip")
        subprocess.check_call([bgzip, self.gtf_path])
        self.gtf_path = self.gtf_path + ".gz"

    def compute_metadata(self, extra_data_dict=None):
        """Compute hash of FASTA and GTF and write reference.json."""
        print("Writing genome metadata JSON file into reference folder...")
        print("Computing hash of genome FASTA file...")
        fasta_hash = compute_hash_of_file(self.fasta_path)
        print("...done\n")

        print("Computing hash of genes GTF file...")
        gtf_hash = compute_hash_of_file(self.gtf_path)
        print("...done\n")

        metadata = {
            cr_constants.REFERENCE_GENOMES_KEY: self.genomes,
            "threads": int(math.ceil(float(self.mem_gb) / 8.0)),
            cr_constants.REFERENCE_MEM_GB_KEY: self.mem_gb,
            cr_constants.REFERENCE_FASTA_HASH_KEY: fasta_hash,
            # note old references have a non-gzipped GTF, so cr_constants.REFERENCE_GENES_GTF_PATH
            # does not contain the ".gz" suffix.
            cr_constants.REFERENCE_GTF_HASH_KEY + ".gz": gtf_hash,
            cr_constants.REFERENCE_INPUT_FASTA_KEY: [
                os.path.basename(x) for x in self.in_fasta_fns
            ],
            cr_constants.REFERENCE_INPUT_GTF_KEY: [os.path.basename(x) for x in self.in_gtf_fns],
            cr_constants.REFERENCE_VERSION_KEY: self.ref_version,
            cr_constants.REFERENCE_MKREF_VERSION_KEY: self.mkref_version,
        }
        if extra_data_dict:
            metadata.update(extra_data_dict)
        new_metadata_json = os.path.join(self.out_dir, cr_constants.REFERENCE_METADATA_FILE)
        with open(new_metadata_json, "w") as f:
            tk_safe_json.dump_numpy(metadata, f, sort_keys=True, indent=4)
        print("...done\n")

    def make_star_index(self):
        """Generate STAR index for use with GEX pipelines."""
        print("Generating STAR genome index (may take over 8 core hours for a 3Gb genome)...")
        new_star_path = os.path.join(self.out_dir, cr_constants.REFERENCE_STAR_PATH)
        star = STAR(new_star_path)
        star.index_reference_with_mem_gb(
            self.fasta_path, self.gtf_path, num_threads=self.num_threads, mem_gb=self.mem_gb
        )
        print("...done.\n")

    def get_contig_lengths(self):
        """Determine contigs in reference and their lengths from the genome.fa.fai file."""
        assert self.fasta_path
        faidx = self.fasta_path + ".fai"
        assert os.path.exists(faidx)
        contig_lengths = {}
        with open(faidx) as fin:
            for line in fin:
                # if a contig name is empty the corresponding index line has leading whitespace
                fields = line.rstrip().split("\t")
                if not fields:
                    continue
                chrom = fields[0]
                length = int(fields[1])
                contig_lengths[chrom] = length
        if not contig_lengths:
            raise GexReferenceError(
                f"The samtools-constructed FASTA index file {self.fasta_path} is empty. The supplied FASTA file(s) "
                f"have no contigs: {self.in_fasta_fns}"
            )
        return contig_lengths

    def build_gex_reference(self, out_dir_exists: bool = False):
        """Construct a cellranger/spaceranger-compatible reference."""
        if not out_dir_exists:
            print(f"Creating new reference folder at {self.out_dir}")
            os.mkdir(self.out_dir)
            print("...done\n")

        self.process_fasta()

        contig_lengths = self.get_contig_lengths()
        self.process_gtf(contig_lengths, no_transcript_fail=False)

        self.validate_gtf()

        try:
            self.make_star_index()
        except subprocess.CalledProcessError as err:
            raise GexReferenceError(
                "Failed to make genome index with STAR.  This can occasionally be caused by setting the argument `memgb` too low.\n"
                f"Error was from running command '{err.cmd[0]}'\n{err}\n\n"
                "Check stdout and stderr for more information."
            ) from err

        self.gzip_gtf()

        self.compute_metadata()

    def format_genome_prefixes(self, prefixes_as_genomes: bool):
        if prefixes_as_genomes:
            self.genome_prefixes = self.genomes
            return

        if len(self.genomes) > 1:
            max_length = max(len(g) for g in self.genomes)
            self.genome_prefixes = []
            for genome in self.genomes:
                genome_prefix = genome
                if len(genome_prefix) < max_length:
                    genome_prefix += "_" * (max_length - len(genome_prefix))
                assert genome_prefix not in self.genome_prefixes
                self.genome_prefixes.append(genome_prefix)
        else:
            self.genome_prefixes = self.genomes

    def write_genome_fasta(self, out_fasta_fn: os.PathLike | str | bytes):
        # check that the first byte of fasta is >
        # helps us detect when a FASTA is gzipped for example
        for fn in self.in_fasta_fns:
            with open(fn, "rb") as fin:
                byte1 = fin.read(1)
                if byte1 == b"":
                    raise GexReferenceError(f"Input FASTA file {fn} is empty")
                if byte1 != b">":
                    raise GexReferenceError(
                        f"Input FASTA file {fn} is invalid. The first byte = {byte1!r} but it must be "
                        "'>'. Note that gzipped FASTA files cannot be processed by mkref."
                    )

        if len(self.genomes) > 1:
            with open(out_fasta_fn, "w") as f:
                for genome_prefix, in_fasta_fn in zip(self.genome_prefixes, self.in_fasta_fns):
                    with open(in_fasta_fn) as g:
                        for line in g:
                            line = line.strip()
                            if line.startswith(">"):
                                line = ">" + genome_prefix + "_" + line[1:]
                            f.write(line + "\n")
        else:
            shutil.copy(self.in_fasta_fns[0], out_fasta_fn)

    def write_genome_gtf(
        self,
        out_gtf_fn: os.PathLike | str | bytes,
        contig_lengths: dict[str, int] | None = None,
        no_transcript_fail: bool = False,
    ):
        has_warning = False
        assert self.genome_prefixes is not None
        with open(out_gtf_fn, "w") as f:
            writer = csv.writer(
                f, delimiter="\t", quoting=csv.QUOTE_NONE, quotechar="", lineterminator="\n"
            )
            for genome_prefix, in_gtf_fn in zip(self.genome_prefixes, self.in_gtf_fns):
                if len(self.genomes) > 1:
                    pfx = genome_prefix + "_"
                    len_pfx = len(pfx)

                    def prefix_func(s: str, p=pfx):
                        return p + s

                    stripped_contig_lengths = (
                        {
                            contig[len_pfx:]: length
                            for contig, length in contig_lengths.items()
                            if contig.startswith(pfx)
                        }
                        if contig_lengths
                        else None
                    )
                else:

                    def prefix_func(s: str):
                        return s

                    stripped_contig_lengths = contig_lengths

                transcript_to_chrom = {}
                cross_chrom_transcripts = set()

                num_rows = 0
                for row, is_comment, properties in self.gtf_reader_iter(
                    in_gtf_fn,
                    contig_lengths=stripped_contig_lengths,
                    no_transcript_fail=no_transcript_fail,
                ):
                    if is_comment:
                        writer.writerow(row)
                        continue

                    chrom = prefix_func(row[0])
                    row[0] = chrom

                    if "transcript_id" in properties:
                        properties["transcript_id"] = prefix_func(properties["transcript_id"])
                        curr_tx = properties["transcript_id"]
                        if curr_tx in transcript_to_chrom and transcript_to_chrom[curr_tx] != chrom:
                            # ignore recurrences of a transcript on different chromosomes - it will break the STAR index
                            cross_chrom_transcripts.add(curr_tx)
                            continue
                        transcript_to_chrom[curr_tx] = chrom
                    if "gene_id" in properties:
                        properties["gene_id"] = prefix_func(properties["gene_id"])
                    if "gene_name" in properties:
                        properties["gene_name"] = prefix_func(properties["gene_name"])

                    row[8] = self.format_properties_dict(properties, uniquify_keys=True)

                    num_rows += int(row[2] == "exon")
                    writer.writerow(row)

                if num_rows == 0:
                    raise GtfParseError(
                        in_gtf_fn, "The supplied GTF file does not contain any exon features"
                    )
                if len(cross_chrom_transcripts) > 0:
                    has_warning = True
                    print(
                        "WARNING: The following transcripts appear on multiple chromosomes in the GTF:"
                    )
                    print("\n".join(list(cross_chrom_transcripts)) + "\n")
                    print(
                        "This can indicate a problem with the reference or annotations. Only the first chromosome will be counted."
                    )
        return has_warning

    def validate_gtf(self):
        """Verifies that the GTF file of the reference can be loaded.

        Uses the Rust GTF indexing code path, and makes sure it succeeds.
        """
        with tempfile.NamedTemporaryFile(dir=self.out_dir, suffix=".json") as out_file:
            try:
                cmd = [_GTF_TO_GENE_INDEX, self.out_dir, out_file.name]
                subprocess.run(cmd, check=True, capture_output=True)

            except subprocess.CalledProcessError as exc:
                raise GexReferenceError(
                    "Error detected in GTF file: " + exc.stderr.decode()
                ) from exc


class STAR:
    def __init__(self, reference_star_path):
        self.reference_star_path = reference_star_path

    def index_reference_with_mem_gb(self, in_fasta_fn, in_gtf_fn, num_threads=1, mem_gb=None):
        # Estimate genome size based on fasta
        genome_size_b = float(os.path.getsize(in_fasta_fn))
        genome_size_gb = float(genome_size_b) / float(10**9)

        # Count number of chromsomes in fasta
        genome_num_chrs = 0
        with open(in_fasta_fn) as f:
            for line in f:
                if line.strip().startswith(">"):
                    genome_num_chrs += 1

        # Per STAR manual recommendations
        sa_index_n_bases = min(14, int(math.log(genome_size_b, 2) / 2 - 1))
        chr_bin_n_bits = min(18, int(math.log(genome_size_b / genome_num_chrs, 2)))

        if mem_gb is None:
            sa_sparse_d = None
            limit_ram = None
        else:
            # Total memory = SA memory + SA index memory + Genome memory
            # SA memory = (8 * genome size) / genomeSAsparseD
            # SA index memory = 8*(4**genomeSAindexNbases)
            # Genome memory = genome size
            sa_index_mem_gb = float(8 * (4**sa_index_n_bases)) / float(10**9)
            genome_mem_gb = genome_size_gb

            min_mem_gb = int(genome_mem_gb + sa_index_mem_gb + 3)
            if mem_gb < min_mem_gb:
                raise GexReferenceError(
                    "STAR requires at least %d GB of memory when aligning reads to your reference.\nPlease start again with --memgb=%d."
                    % (min_mem_gb, min_mem_gb)
                )

            limit_ram = mem_gb * 1024**3

            # 2 GB of buffer space
            mem_gb = max(1, mem_gb - 2)

            sa_sparse_d = float(8 * genome_size_gb) / (
                float(mem_gb) - genome_mem_gb - sa_index_mem_gb
            )
            sa_sparse_d = max(1, int(math.ceil(sa_sparse_d)))

        self.index_reference(
            in_fasta_fn,
            in_gtf_fn,
            num_threads=num_threads,
            sa_sparse_d=sa_sparse_d,
            sa_index_n_bases=sa_index_n_bases,
            chr_bin_n_bits=chr_bin_n_bits,
            limit_ram=limit_ram,
        )

    def index_reference(
        self,
        in_fasta_fn,
        in_gtf_fn,
        num_threads=1,
        sa_sparse_d=None,
        sa_index_n_bases=None,
        chr_bin_n_bits=None,
        limit_ram=None,
    ):
        if os.path.exists(self.reference_star_path):
            raise Exception(f"STAR reference path {self.reference_star_path} already exists")

        os.mkdir(self.reference_star_path)

        args = [
            os.path.join(_LIB_BIN, "STAR"),
            "--runMode",
            "genomeGenerate",
            "--genomeDir",
            self.reference_star_path,
            "--runThreadN",
            str(num_threads),
            "--genomeFastaFiles",
            in_fasta_fn,
            "--sjdbGTFfile",
            in_gtf_fn,
        ]
        if limit_ram is not None:
            args += ["--limitGenomeGenerateRAM", str(limit_ram)]
        if sa_sparse_d is not None:
            args += ["--genomeSAsparseD", str(sa_sparse_d)]
        if sa_index_n_bases is not None:
            args += ["--genomeSAindexNbases", str(sa_index_n_bases)]
        if chr_bin_n_bits is not None:
            args += ["--genomeChrBinNbits", str(chr_bin_n_bits)]

        try:
            tk_subproc.check_call(args)
        except subprocess.CalledProcessError as err:
            if err.returncode == -4:
                raise RuntimeError(
                    "mkref has failed because it is running on a computer that does not support some required instructions (AVX).  Please contact support@10xgenomics.com to learn how to work around this issue."
                )
            else:
                raise err
