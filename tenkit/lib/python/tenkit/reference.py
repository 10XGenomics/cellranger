import os.path
import pyfasta
import re
from tenkit.regions import Regions
from collections import defaultdict
import gzip

# 10X generated references that should
# have genes, regions, and snps files
KNOWN_GENOMES = ['10X_hg19_ucsc', '10X_b37', '10X_GRCh38_no_alt_decoy']

def open_reference(reference_path):
    ''' Open a reference fasta and rename the contigs to strip any fasta comments'''
    fasta = pyfasta.Fasta(get_fasta(reference_path))

    new_fasta = {}

    for (k,v) in fasta.iteritems():
        key_prefix = k.split(" ")[0]
        new_fasta[key_prefix] = v

    return new_fasta

def get_metadata(reference_path):
    return os.path.join(reference_path, "metadata.json")

def get_fasta(reference_path):
    '''Convention for location of reference fasta in a reference path'''
    return os.path.join(reference_path, "fasta", "genome.fa")

def is_tenx(reference_path):
    return (not reference_path is None and get_genome(reference_path) in KNOWN_GENOMES)

def get_genome(reference_path):
    ''' Load the canonical name of a reference. By convention this is stored in <reference_path>/genome '''
    genome = None

    genome_name_file = os.path.join(reference_path, "genome")

    if os.path.exists(genome_name_file):
        with open(genome_name_file) as f:
            genome = f.readline().strip()

    return genome

def get_fasta_contig_index_path(reference_path):
    return os.path.join(reference_path, "fasta", "genome.fa.fai")

def get_primary_contigs(reference_path):
    return os.path.join(reference_path, 'fasta', 'primary_contigs.txt')

def load_primary_contigs(reference_path):
    '''Load set of primary contigs for variant and SV calling from reference_path.
       If now primary_contigs.txt file is specified, return all contigs. If  reference_path
       is a known 10x reference genome and has no primary_contigs.txt, filter the known bad contigs '''

    if not reference_path is None and os.path.exists(get_primary_contigs(reference_path)):
        # If we have a primary_contigs.txt file, use it
        with open(get_primary_contigs(reference_path), 'r') as f:
            primary_contigs = set([line.strip() for line in f.readlines()])

    else:
        # Default is to include all contigs
        # Otherwise implement the old contig filters
        ref = open_reference(reference_path)
        primary_contigs = set(ref.keys())

        if is_tenx(reference_path):
            primary_contigs = set(chrom for chrom in primary_contigs if not ('random' in chrom or 'U' in chrom or 'hap' in chrom or chrom == 'hs37d5'))

    return primary_contigs

def get_sex_chromosomes(reference_path):
    return os.path.join(reference_path, 'fasta', 'sex_chromosomes.tsv')

def get_centromere_regions(reference_path):
    return os.path.join(reference_path, 'regions', 'centromeres.tsv')

def load_sex_chromosomes(reference_path, chr_type):
    if not reference_path is None and os.path.exists(get_sex_chromosomes(reference_path)):
        with open(get_sex_chromosomes(reference_path),'r') as f:
            for line in f:
                tokens = line.strip().split()
                if len(tokens) > 1 and tokens[0] == chr_type:
                    return set(chrom for chrom in tokens[1:])
    else:
        return None

def load_male_chromosomes(reference_path):
    return load_sex_chromosomes(reference_path, 'male')

def load_autosomal_chromosomes(reference_path):
    return load_sex_chromosomes(reference_path, 'autosomal')



def get_sv_blacklist(reference_path):
    return os.path.join(reference_path, "regions", "sv_blacklist.bed")

def get_segdups(reference_path):
    return os.path.join(reference_path, "regions", "segdups.bedpe")

def get_loupe_genes(reference_path):
    ''' Gene / exon file consumed by Loupe '''
    return os.path.join(reference_path, "genes", "gene_annotations.gtf.gz")

def get_gene_boundaries(reference_path):
    ''' BED file of gene extents consumed by ANALYZE_SNPINDEL_CALLS '''
    return os.path.join(reference_path, "genes", "gene_annotations.gtf.gz")

def load_gene_boundaries(reference_path, protein_coding=True):
    ''' Return a dict from chromosome name to gene regions.

    If the reference path does not contain a gene annotation gtf, then
    an empty dictionary is returned.

    Args:
    - reference_path: Path to 10X reference folder
    - protein_coding: Whether to only include protein coding genes (True) or 
    all genes (False). Default is True.
    '''
    
    gene_file = get_gene_boundaries(reference_path)
    if not os.path.isfile(gene_file):
        return {}

    valid_chroms = open_reference(reference_path).keys()
            
    gene_regions = defaultdict(list)
    with gzip.open(gene_file, 'r') as gene_fn:
        for line in gene_fn:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if fields[2] == 'gene' and (not protein_coding or fields[1] == 'protein_coding'):
                chrom, start, stop = fields[0], int(fields[3]), int(fields[4])
                if not chrom in valid_chroms and chrom[0:3] != 'chr':
                    if chrom == 'MT' and 'chrM' in valid_chroms:
                        chrom = 'chrM'
                    elif 'chr' + chrom in valid_chroms:
                        chrom = 'chr' + chrom
                gene_regions[chrom].append((start, stop))
                
    for chrom, regions in gene_regions.iteritems():
        gene_regions[chrom] = Regions(regions)
    return gene_regions
                                                                        
def get_unambiguous_regions(reference_path):
    '''Calculate regions corresponding to unambiguous bases'''
    chrom_map = {}
    for chrom, seq in open_reference(reference_path).items():
        regions = [(m.start(), m.end()) for m in re.finditer('[acgtACGT]+', seq[:])]
        chrom_map[chrom] = Regions(regions=regions)
    return chrom_map
