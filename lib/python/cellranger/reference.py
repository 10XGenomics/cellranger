#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
import collections
import cPickle
import csv
import itertools
import json
import math
import numpy as np
import os
import subprocess
import sys
import re
import tenkit.log_subprocess as tk_subproc
import tenkit.safe_json as tk_safe_json
import tenkit.seq as tk_seq
import tenkit.bam as tk_bam
import cellranger.constants as cr_constants
import cellranger.h5_constants as h5_constants
import cellranger.io as cr_io

class GtfParser:
    GTF_ERROR_TXT = 'Please fix your GTF and start again.'

    def gtf_reader_iter(self, filename):
        with open(filename, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            for i, row in enumerate(reader):
                if len(row) == 0:
                    continue
                if row[0].startswith('#'):
                    yield row, True, None
                    continue

                if len(row) != 9:
                    sys.exit("Invalid number of columns in GTF line %d: %s\n\n%s" % (i+1, '\t'.join(row), self.GTF_ERROR_TXT))

                strand = row[6]
                if strand not in cr_constants.STRANDS:
                    sys.exit('Invalid strand in GTF line %d: %s\n\n%s' % (i+1, '\t'.join(row), self.GTF_ERROR_TXT))

                annotation = row[2]
                properties = self.get_properties_dict(row[8])
                if annotation == 'exon':
                    if 'transcript_id' not in properties:
                        sys.exit("Property 'transcript_id' not found in GTF line %d: %s\n\n%s" % (i+1, '\t'.join(row), self.GTF_ERROR_TXT))
                    if ';' in properties['transcript_id']:
                        sys.exit("Property 'transcript_id' has invalid character ';' in GTF line %d: %s\n\n%s" % (i+1, '\t'.join(row), self.GTF_ERROR_TXT))
                    if re.search(r'\s', properties['transcript_id']) is not None:
                        sys.exit("Property 'transcript_id' has invalid whitespace character in GTF line %d: %s\n\n%s" % (i+1, '\t'.join(row), self.GTF_ERROR_TXT))
                    if 'gene_id' not in properties:
                        sys.exit("Property 'gene_id' not found in GTF line %d: %s\n\n%s" % (i+1, '\t'.join(row), self.GTF_ERROR_TXT))
                    if ';' in properties['gene_id']:
                        sys.exit("Property 'gene_id' has invalid character ';' in GTF line %d: %s\n\n%s" % (i+1, '\t'.join(row), self.GTF_ERROR_TXT))

                yield row, False, properties

    def load_gtf(self, in_gtf_fn, fasta_parser=None):
        transcripts = {}
        gene_to_transcripts = collections.OrderedDict()

        for row, is_comment, properties in self.gtf_reader_iter(in_gtf_fn):
            if is_comment:
                continue

            chrom, _, annotation, start, end, _, strand, _, properties_str = row

            if annotation != "exon":
                continue

            start = int(start) - 1
            end = int(end)
            length = abs(end - start)
            transcript_id = properties['transcript_id']
            gene_id = properties['gene_id']
            gene_name = properties.get('gene_name', gene_id)
            gene = cr_constants.Gene(gene_id, gene_name, None, None, None)

            if transcript_id not in transcripts:
                transcripts[transcript_id] = cr_constants.Transcript(gene, None, None, [])

            if gene not in gene_to_transcripts:
                gene_to_transcripts[gene] = set()

            assert transcripts[transcript_id].gene == gene
            transcripts[transcript_id].intervals.append(cr_constants.Interval(chrom, start, end, length, strand))
            gene_to_transcripts[gene].add(transcript_id)

        # Transcript length and GC content
        transcript_lengths = {}
        transcript_gc_contents = {}
        for transcript_id, transcript in transcripts.iteritems():
            transcript_lengths[transcript_id] = sum([interval.length for interval in transcript.intervals])
            if fasta_parser is not None:
                transcript_gc_contents[transcript_id] = fasta_parser.get_transcript_gc_content(transcript)

        # Gene length, GC content and start + end positions
        genes = []
        for gene, transcript_ids in gene_to_transcripts.iteritems():
            length = np.median([transcript_lengths[transcript_id] for transcript_id in transcript_ids])
            gc_content = np.median([transcript_gc_contents[transcript_id] for transcript_id in transcript_ids])

            transcript_intervals = []
            for transcript_id in transcript_ids:
                transcript_intervals += transcripts[transcript_id].intervals
            transcript_intervals.sort(key=lambda interval: interval.chrom)

            intervals = []
            for chrom, chrom_intervals_iter in itertools.groupby(transcript_intervals,
                    lambda interval: interval.chrom):
                chrom_intervals = list(chrom_intervals_iter)
                start = min([interval.start for interval in chrom_intervals])
                end = max([interval.end for interval in chrom_intervals])
                interval = cr_constants.Interval(chrom, start, end, end - start, None)
                intervals.append(interval)

            gene = cr_constants.Gene(gene.id, gene.name, length, gc_content, intervals)
            genes.append(gene)

            for transcript_id in transcript_ids:
                transcripts[transcript_id] = cr_constants.Transcript(
                    gene, transcript_lengths[transcript_id],
                    transcript_gc_contents[transcript_id],
                    transcripts[transcript_id].intervals)

        return transcripts, genes

    def get_properties_dict(self, properties_str):
        if isinstance(properties_str, dict):
            return properties_str

        properties = collections.OrderedDict()
        pattern = re.compile('(\S+?)\s*"(.*?)"')
        for m in re.finditer(pattern, properties_str):
            key = m.group(1)
            value = m.group(2)
            properties[key] = value
        return properties

    def format_properties_dict(self, properties):
        properties_str = []
        for key, value in properties.iteritems():
            properties_str.append('%s "%s"' % (key, value))
        return '; '.join(properties_str)

class GtfBuilder(GtfParser):
    def __init__(self, in_gtf_fn, out_gtf_fn, attributes={}):
        self.in_gtf_fn = in_gtf_fn
        self.out_gtf_fn = out_gtf_fn
        self.attributes = attributes

    def build_gtf(self):
        print "Writing new genes GTF file (may take 10 minutes for a 1GB input GTF file)..."
        with open(self.out_gtf_fn, 'wb') as f:
            writer = csv.writer(f, delimiter='\t', quoting=csv.QUOTE_NONE, quotechar='')
            for row, is_comment, properties in self.gtf_reader_iter(self.in_gtf_fn):
                if is_comment:
                    writer.writerow(row)
                    continue

                remove = False
                for key, value in properties.iteritems():
                    if key in self.attributes and value not in self.attributes[key]:
                        remove = True

                if not remove:
                    writer.writerow(row)

        print "...done\n"

class ReferenceBuilder(GtfParser):
    def __init__(self, genomes, in_fasta_fns, in_gtf_fns, out_dir, ref_version, mkref_version, num_threads=1, mem_gb=None):
        self.genomes = genomes
        self.in_fasta_fns = in_fasta_fns
        self.in_gtf_fns = in_gtf_fns
        self.out_dir = out_dir
        self.num_threads = num_threads
        self.mem_gb = mem_gb
        self.ref_version = ref_version
        self.mkref_version = mkref_version

        self.format_genome_prefixes()

    def build_reference(self):
        print "Creating new reference folder at %s" % self.out_dir
        os.mkdir(self.out_dir)
        print "...done\n"

        print "Writing genome FASTA file into reference folder..."
        new_genome_fasta = os.path.join(self.out_dir, cr_constants.REFERENCE_FASTA_PATH)
        os.mkdir(os.path.dirname(new_genome_fasta))
        self.write_genome_fasta(new_genome_fasta)
        print "...done\n"

        print "Computing hash of genome FASTA file..."
        fasta_hash = cr_io.compute_hash_of_file(new_genome_fasta)
        print "...done\n"

        print "Indexing genome FASTA file..."
        subprocess.check_call(["samtools", "faidx", new_genome_fasta])
        print "...done\n"

        print "Writing genes GTF file into reference folder..."
        new_gene_gtf = os.path.join(self.out_dir, cr_constants.REFERENCE_GENES_GTF_PATH)
        os.mkdir(os.path.dirname(new_gene_gtf))
        self.write_genome_gtf(new_gene_gtf)
        print "...done\n"

        print "Computing hash of genes GTF file..."
        gtf_hash = cr_io.compute_hash_of_file(new_gene_gtf)
        print "...done\n"

        print "Writing genes index file into reference folder (may take over 10 minutes for a 3Gb genome)..."
        new_gene_index = os.path.join(self.out_dir, cr_constants.REFERENCE_GENES_INDEX_PATH)
        os.mkdir(os.path.dirname(new_gene_index))
        self.write_genome_gene_index(new_gene_index, new_gene_gtf, new_genome_fasta)
        print "...done\n"

        print "Writing genome metadata JSON file into reference folder..."
        metadata = {
            cr_constants.REFERENCE_GENOMES_KEY: self.genomes,
            cr_constants.REFERENCE_NUM_THREADS_KEY: int(math.ceil(float(self.mem_gb) / 8.0)),
            cr_constants.REFERENCE_MEM_GB_KEY: self.mem_gb,
            cr_constants.REFERENCE_FASTA_HASH_KEY: fasta_hash,
            cr_constants.REFERENCE_GTF_HASH_KEY: gtf_hash,
            cr_constants.REFERENCE_INPUT_FASTA_KEY: [os.path.basename(x) for x in self.in_fasta_fns],
            cr_constants.REFERENCE_INPUT_GTF_KEY: [os.path.basename(x) for x in self.in_gtf_fns],
            cr_constants.REFERENCE_VERSION_KEY: self.ref_version,
            cr_constants.REFERENCE_MKREF_VERSION_KEY: self.mkref_version,
        }
        new_metadata_json = os.path.join(self.out_dir, cr_constants.REFERENCE_METADATA_FILE)
        with open(new_metadata_json, 'w') as f:
            json.dump(tk_safe_json.json_sanitize(metadata), f, sort_keys=True, indent=4)
        print "...done\n"

        print "Generating STAR genome index (may take over 8 core hours for a 3Gb genome)..."
        new_star_path = os.path.join(self.out_dir, cr_constants.REFERENCE_STAR_PATH)
        star = STAR(new_star_path)
        star.index_reference_with_mem_gb(new_genome_fasta, new_gene_gtf,
                                         num_threads=self.num_threads,
                                         mem_gb=self.mem_gb)
        print "...done.\n"

        print ">>> Reference successfully created! <<<\n"
        print "You can now specify this reference on the command line:"
        print "cellranger --transcriptome=%s ..." % self.out_dir

    def format_genome_prefixes(self):
        if len(self.genomes) > 1:
            max_length = max([len(g) for g in self.genomes])
            self.genome_prefixes = []
            for genome in self.genomes:
                genome_prefix = genome
                if len(genome_prefix) < max_length:
                    genome_prefix += '_' * (max_length - len(genome_prefix))
                assert genome_prefix not in self.genome_prefixes
                self.genome_prefixes.append(genome_prefix)
        else:
            self.genome_prefixes = self.genomes

    def write_genome_fasta(self, out_fasta_fn):
        if len(self.genomes) > 1:
            with open(out_fasta_fn, 'w') as f:
                for genome_prefix, in_fasta_fn in itertools.izip(self.genome_prefixes, self.in_fasta_fns):
                    with open(in_fasta_fn, 'r') as g:
                        for line in g:
                            line = line.strip()
                            if line.startswith('>'):
                                line = '>' + genome_prefix + '_' + line[1:]
                            f.write(line + '\n')
        else:
            cr_io.copy(self.in_fasta_fns[0], out_fasta_fn)

    def write_genome_gtf(self, out_gtf_fn):
        with open(out_gtf_fn, 'wb') as f:
            writer = csv.writer(f, delimiter='\t', quoting=csv.QUOTE_NONE, quotechar='')
            for genome_prefix, in_gtf_fn in itertools.izip(self.genome_prefixes, self.in_gtf_fns):
                if len(self.genomes) > 1:
                    prefix_func = lambda s: '%s_%s' % (genome_prefix, s)
                else:
                    prefix_func = lambda s: s

                transcript_to_chrom = {}
                cross_chrom_transcripts = set()
                for row, is_comment, properties in self.gtf_reader_iter(in_gtf_fn):
                    if is_comment:
                        writer.writerow(row)
                        continue

                    chrom = prefix_func(row[0])
                    row[0] = chrom

                    if 'transcript_id' in properties:
                        properties['transcript_id'] = prefix_func(properties['transcript_id'])
                        curr_tx = properties['transcript_id']
                        if curr_tx in transcript_to_chrom and transcript_to_chrom[curr_tx] != chrom:
                            # ignore recurrences of a transcript on different chromosomes - it will break the STAR index
                            cross_chrom_transcripts.add(curr_tx)
                            continue
                        transcript_to_chrom[curr_tx] = chrom
                    if 'gene_id' in properties:
                        properties['gene_id'] = prefix_func(properties['gene_id'])
                    if 'gene_name' in properties:
                        properties['gene_name'] = prefix_func(properties['gene_name'])

                    row[8] = self.format_properties_dict(properties)

                    writer.writerow(row)

                if len(cross_chrom_transcripts) > 0:
                    print "WARNING: The following transcripts appear on multiple chromosomes in the GTF:"
                    print '\n'.join(list(cross_chrom_transcripts)) + '\n'
                    print "This can indicate a problem with the reference or annotations. Only the first chromosome will be counted."

    def write_genome_gene_index(self, out_pickle_fn, in_gtf_fn, in_fasta_fn):
        gene_index = GeneIndex(in_gtf_fn, in_fasta_fn)
        gene_index.save_pickle(out_pickle_fn)

class FastaParser:
    def __init__(self, in_fasta_fn):
        self.chroms = self.load_fasta(in_fasta_fn)

    @staticmethod
    def _get_chrom_name(fasta_header):
        match = re.search(r'>(\S+)', fasta_header)
        return match.groups(1)[0]

    def load_fasta(self, in_fasta_fn):
        chroms = {}

        current_chrom, current_seq = None, None
        with open(in_fasta_fn, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_chrom is not None:
                        chroms[current_chrom] = ''.join(current_seq)
                    current_chrom = FastaParser._get_chrom_name(line)
                    current_seq = []
                else:
                    current_seq.append(line)

        if current_chrom is not None:
            chroms[current_chrom] = ''.join(current_seq)

        return chroms

    def is_valid_interval(self, chrom, start, end):
        """ Determine whether the half-open interval [start,end) is within the bounds of chrom """
        return chrom in self.chroms and start >= 0 and end <= len(self.chroms[chrom])

    def get_sequence(self, chrom, start, end, strand=cr_constants.FORWARD_STRAND):
        """ Get genomic sequence for the half-open interval [start,end) """
        seq = self.chroms[chrom][start:end]
        if strand == cr_constants.FORWARD_STRAND:
            return seq
        elif strand == cr_constants.REVERSE_STRAND:
            return tk_seq.get_rev_comp(seq)
        else:
            raise Exception("Invalid strand: %s" % strand)

    def get_transcript_gc_content(self, transcript_obj):
        pattern = re.compile('[cCgG]')

        gc, length = 0, 0
        for interval in transcript_obj.intervals:
            if interval.chrom not in self.chroms:
                continue

            seq = self.chroms[interval.chrom][interval.start:interval.end]
            gc += len(re.findall(pattern, seq))
            length += interval.length

        if length > 0:
            return float(gc) / float(length)
        else:
            return 0

# NOTE: these stub classes are necessary to maintain backwards compatibility with old refdata (1.2 or older)
class IntervalTree: pass
class Region: pass
class IntergenicRegion(Region): pass
class IntronicRegion(Region): pass
class ExonicRegion(Region): pass

class GeneIndex(GtfParser):
    def __init__(self, in_gtf_fn, in_fasta_fn):
        self.in_gtf_fn = in_gtf_fn
        self.in_fasta_fn = in_fasta_fn

        fasta_parser = FastaParser(self.in_fasta_fn)

        self.transcripts, self.genes = self.load_gtf(self.in_gtf_fn, fasta_parser=fasta_parser)
        self.gene_ids_map = {gene.id: i for i, gene in enumerate(self.genes)}

    def save_pickle(self, out_pickle_fn):
        with open(out_pickle_fn, 'wb') as f:
            cPickle.dump(self, f, protocol=cPickle.HIGHEST_PROTOCOL)

    @staticmethod
    def load_pickle(in_pickle_fn):
        with open(in_pickle_fn, 'rb') as f:
            return cPickle.load(f)

    def get_transcript_length(self, transcript):
        if transcript in self.transcripts:
            return self.transcripts[transcript].length
        return None

    def get_transcript_gc_content(self, transcript):
        if transcript in self.transcripts:
            return self.transcripts[transcript].gc_content
        return None

    def get_gene_from_transcript(self, transcript):
        if transcript in self.transcripts:
            return self.transcripts[transcript].gene
        return None

    def gene_id_to_int(self, gene_id):
        if gene_id in self.gene_ids_map:
            return self.gene_ids_map[gene_id]
        return None

    def get_genes(self):
        return self.genes

    def get_gene(self, gene_id):
        return self.genes[self.gene_id_to_int(gene_id)]

    def get_gene_lengths(self):
        return [gene.length for gene in self.genes]

    def get_gene_gc_contents(self):
        return [gene.gc_content for gene in self.genes]

class STAR:
    def __init__(self, reference_star_path):
        self.reference_star_path = reference_star_path

    def index_reference_with_mem_gb(self, in_fasta_fn, in_gtf_fn, num_threads=1, mem_gb=None):
        # Estimate genome size based on fasta
        genome_size_b = float(os.path.getsize(in_fasta_fn))
        genome_size_gb = float(genome_size_b) / float(10**9)

        # Count number of chromsomes in fasta
        genome_num_chrs = 0
        with open(in_fasta_fn, 'r') as f:
            for line in f:
                if line.strip().startswith('>'):
                    genome_num_chrs += 1

        # Per STAR manual recommendations
        sa_index_n_bases = min(14, int(math.log(genome_size_b, 2)/2 - 1))
        chr_bin_n_bits = min(18, int(math.log(genome_size_b/genome_num_chrs, 2)))

        if mem_gb is None:
            sa_sparse_d = None
            limit_ram = None
        else:
            # Total memory = SA memory + SA index memory + Genome memory
            # SA memory = (8 * genome size) / genomeSAsparseD
            # SA index memory = 8*(4**genomeSAindexNbases)
            # Genome memory = genome size
            sa_index_mem_gb = float(8 * (4 ** sa_index_n_bases)) / float(10 ** 9)
            genome_mem_gb = genome_size_gb

            min_mem_gb = int(genome_mem_gb + sa_index_mem_gb + 3)
            if mem_gb < min_mem_gb:
                sys.exit("WARNING: STAR requires at least %d GB of memory when aligning reads to your reference.\nPlease start again with --memgb=%d." % (min_mem_gb, min_mem_gb))

            limit_ram = mem_gb * 1024**3

            # 2 GB of buffer space
            mem_gb = max(1, mem_gb - 2)

            sa_sparse_d = float(8 * genome_size_gb) / (float(mem_gb) - genome_mem_gb  - sa_index_mem_gb)
            sa_sparse_d = max(1, int(math.ceil(sa_sparse_d)))

        self.index_reference(in_fasta_fn, in_gtf_fn, num_threads=num_threads,
                             sa_sparse_d=sa_sparse_d,
                             sa_index_n_bases=sa_index_n_bases,
                             chr_bin_n_bits=chr_bin_n_bits,
                             limit_ram=limit_ram)

    def index_reference(self, in_fasta_fn, in_gtf_fn, num_threads=1, sa_sparse_d=None, sa_index_n_bases=None,
                        chr_bin_n_bits=None, limit_ram=None):
        if os.path.exists(self.reference_star_path):
            raise Exception('STAR reference path %s already exists' % self.reference_star_path)

        os.mkdir(self.reference_star_path)

        args = ['STAR', '--runMode', 'genomeGenerate', '--genomeDir', self.reference_star_path,
                '--runThreadN', str(num_threads), '--genomeFastaFiles', in_fasta_fn,
                '--sjdbGTFfile', in_gtf_fn]
        if limit_ram is not None:
            args += ['--limitGenomeGenerateRAM', str(limit_ram)]
        if sa_sparse_d is not None:
            args += ['--genomeSAsparseD', str(sa_sparse_d)]
        if sa_index_n_bases is not None:
            args += ['--genomeSAindexNbases', str(sa_index_n_bases)]
        if chr_bin_n_bits is not None:
            args += ['--genomeChrBinNbits', str(chr_bin_n_bits)]

        tk_subproc.check_call(args)

    def align(self, read1_fastq_fn, read2_fastq_fn,
              out_genome_bam_fn,
              threads, cwd=None,
              max_report_alignments_per_read=-1,
              read_group_tags=None):
        if cwd is None:
            cwd = os.getcwd()

        if read2_fastq_fn is None:
            read2_fastq_fn = ''

        args = [
            'STAR', '--genomeDir', self.reference_star_path,
            '--outSAMmultNmax', str(max_report_alignments_per_read),
            '--runThreadN', str(threads),
            '--readNameSeparator', 'space',
            '--outSAMunmapped', 'Within', 'KeepPairs',
            '--outSAMtype', 'SAM',
            '--outStd', 'SAM',
            '--outSAMorder', 'PairedKeepInputOrder',
        ]

        if read_group_tags is not None:
            args.append('--outSAMattrRGline')
            args.extend(read_group_tags)

        args.append('--readFilesIn')
        if read1_fastq_fn.endswith(h5_constants.GZIP_SUFFIX):
            args.append('<(gzip -c -d \'%s\')' % read1_fastq_fn)
            if read2_fastq_fn:
                args.append('<(gzip -c -d \'%s\')' % read2_fastq_fn)

        elif read1_fastq_fn.endswith(h5_constants.LZ4_SUFFIX):
            args.append('<(lz4 -c -d \'%s\')' % read1_fastq_fn)
            if read2_fastq_fn:
                args.append('<(lz4 -c -d \'%s\')' % read2_fastq_fn)

        else:
            args.append(read1_fastq_fn)
            if read2_fastq_fn:
                args.append(read2_fastq_fn)

        if out_genome_bam_fn == cr_constants.BAM_FILE_STREAM:
            # stream to pipe for downstream processing
            # NOTE: this feature is unused in the standard pipeline
            # HACK: see https://github.com/pysam-developers/pysam/issues/355
            parent_read, child_write = os.pipe()
            try:
                tk_subproc.Popen(args, stdout=child_write)
            finally:
                os.close(child_write)
            os.dup2(parent_read, sys.stdin.fileno())
            # now streaming output can be read using pysam.Samfile('-', 'r')
            # NOTE: since this does not await termination of the process, we can't reliably check the return code
        else:
            # NOTE: We'd like to pipe fastq files through a decompressor and feed those
            # streams into STAR.
            # STAR provides --readFilesCommand which will do this. But it uses a named pipe which
            # breaks on some filesystems.

            # We could also use anonymous pipes but we'd need a way to refer to them
            # on the command line and apparently not all systems support the same
            # /dev/fdN or procfs-like paths.

            # So we're forced to use the shell and process subsitution, as is recommended
            # here: https://groups.google.com/forum/#!msg/rna-star/MQdL1WxkAAw/eG6EatoOCgAJ

            # Wrap arguments in single quotes
            quoted_args = []
            for arg in args:
                if arg.startswith('<'):
                    # We want the shell to interpret this as a process substitution
                    quoted_args.append(arg)

                elif "'" in arg:
                    # We can't escape single quotes within single quoted strings.
                    # But we can concatenate different quoting mechanisms.
                    # ' => '"'"'
                    # This is relevant if the RG string contains quotes, which
                    # can happen if the user specifies such a library name.
                    arg = arg.replace("'", "'\"'\"'")
                    quoted_args.append("'%s'" % arg)

                else:
                    # Normal argument
                    quoted_args.append("'%s'" % arg)

            star_cmd = ' '.join(quoted_args)

            star = tk_subproc.Popen(star_cmd, stdout=subprocess.PIPE, cwd=cwd, shell=True, executable='bash')
            star_log = os.path.join(cwd, 'Log.out')

            with open(out_genome_bam_fn, 'w') as f:
                view_cmd = ['samtools', 'view', '-Sb', '-']
                view = tk_subproc.Popen(view_cmd, stdin=star.stdout, stdout=f, cwd=cwd)
                view.communicate()

            try:
                # Ensure that STAR process terminated so we can get a returncode
                star.communicate()
                cr_io.check_completed_process(star, args[0])

                # check samtools status
                cr_io.check_completed_process(view, ' '.join(view_cmd))

            except cr_io.CRCalledProcessError as e:
                # Give the user the path to STAR's log
                raise cr_io.CRCalledProcessError(e.msg + ' Check STAR logs for errors: %s .' % star_log)

            # check for empty BAM
            if tk_bam.bam_is_empty(out_genome_bam_fn):
                raise Exception('Aligned BAM is empty - check STAR logs for errors: %s .' % star_log )

def get_ref_name_from_genomes(genomes):
    return '_and_'.join(genomes)
