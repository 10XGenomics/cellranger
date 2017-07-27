#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
@package pyssw
@brief Python standalone program for ssw alignment using the C library
Complete-Striped-Smith-Waterman-Library
Biopython module is require for fastq/fastq parsing
@copyright  [The MIT licence](http://opensource.org/licenses/MIT)
@author     Adrien Leger - 2014
* <adrien.leger@gmail.com>
* <adrien.leger@inserm.fr>
* <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)
"""

#~~~~~~~GLOBAL IMPORTS~~~~~~~#
# Standard library packages
import optparse
import sys
from time import time
import gzip

#~~~~~~~MAIN FUNCTION~~~~~~~#
def main (opt):

    print ("Inport subject sequence")
    # Import fasta subject
    if opt.subject.rpartition(".")[2].lower() == "gz":
        subject_handle = gzip.open(opt.subject, "r")
    else:
        subject_handle = open(opt.subject, "r")
    subject = SeqIO.read(subject_handle, "fasta")

    print ("Inport query sequences and count the number of sequences")
    # Import fasta subject
    if opt.query.rpartition(".")[2].lower() == "gz":
        nseq = count_seq(opt.query, opt.qtype, True)
        query_handle = gzip.open(opt.query, "r")
    else:
        nseq = count_seq(opt.query, opt.qtype, False)
        query_handle = open(opt.query, "r")
    query_gen = SeqIO.parse(query_handle, opt.qtype)

    print("{} contains {} sequences to align".format(opt.query, nseq))
    # Calculate a step list for the progress bar
    nseq_list = [int(nseq*i/100.0) for i in range(5,101,5)]

    print ("Initialize ssw aligner with the subject sequence")
    # Init the an Aligner object with the reference value
    ssw = Aligner(
        str(subject.seq),
        match=int(opt.match),
        mismatch=int(opt.mismatch),
        gap_open=int(opt.gap_open),
        gap_extend= int(opt.gap_extend),
        report_secondary=False,
        report_cigar=True)

    # Write the header of the SAM file
    with open("result.sam", "w") as f:
        f.write("@HD\tVN:1.0\tSO:unsorted\n")
        f.write("@SQ\tSN:{}\tLN:{}\n".format(subject.id, len(subject.seq)))
        f.write("@PG\tID:Striped-Smith-Waterman\tPN:pyssw\tVN:0.1\n")
        f.write("@CO\tScore_values = match {}, mismatch {}, gap_open {}, gap_extend {}\n".format(
            opt.match,
            opt.mismatch,
            opt.gap_open,
            opt.gap_extend))
        f.write("@CO\tFilter Options = min_score {}, min_len {}\n".format(
            opt.min_score,
            opt.min_len))

        print ("Starting alignment of queries against the subject sequence")
        start = time()
        # Align each query along the subject an write result in a SAM file
        i = 0
        for query in query_gen:

            # Find the best alignment
            if opt.reverse:
                al, orient = find_best_align (ssw, query, float(opt.min_score), int(opt.min_len))
            else:
                al, orient = ssw.align(str(query.seq), float(opt.min_score), int(opt.min_len)), True

            # If valid match found
            if al:
                f.write(sam_line(
                    qname=query.id,
                    flag=0 if orient else 16,
                    rname=subject.id,
                    pos=al.ref_begin+1,
                    cigar=al.cigar_string,
                    seq=str(query.seq),
                    qual=SeqIO.QualityIO._get_sanger_quality_str(query) if opt.qtype == "fastq" else "*",
                    tags=["AS:i:{}".format(al.score)]))

            # If no valid match found and -u flag activated (report unaligned)
            elif opt.unaligned:
                f.write(sam_line(
                    qname=query.id,
                    flag=4,
                    seq=str(query.seq),
                    qual=SeqIO.QualityIO._get_sanger_quality_str(query) if opt.qtype == "fastq" else "*"))
            # Else = match unreported

            # Progress bar
            i+=1
            if i in nseq_list:
                frac = i/float(nseq)
                t = time()-start
                print ("{} sequences \t{}% \tRemaining time = {}s".format(i, int(frac*100), round(t/frac-t, 2)))

        print ("\n{} Sequences processed in {}s".format(i, round(time()-start, 2)))

#~~~~~~~HELPER FUNCTIONS~~~~~~~#


def sam_line (qname='*', flag=4, rname='*', pos=0, mapq=0, cigar='*', rnext='*', pnext=0, tlen=0, seq='*', qual='*', tags=None):
    """
    Return a minimal sam line = by default return an undetermined sam line. Check the document
    [SAM Format Specification](http://samtools.sourceforge.net/SAM1.pdf) for a full description.
    @param qname Query template NAME
    @param flag bitwise FLAG
    @param rname Reference sequence NAME of the alignment
    @param pos 1-based leftmost mapping POSition of the first matching base
    @param mapq MAPping Quality
    @param cigar CIGAR string
    @param rnext Reference sequence name of the primary alignment of the mate
    @param pnext 1-based leftmost position of the primary alignment of the mate
    @param tlen signed observed Template LENgth
    @param seq segment SEQuence
    @param qual ASCII of base QUALity plus 33
    @param tags list of optional tags
    @return A Sam alignment line
    """
    if tags:
        return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
            qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual, " ".join(tags))
    else:
        return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
            qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual)

def find_best_align (ssw, query, min_score, min_len):

    # Align reverse and forward query
    forward_al = ssw.align(str(query.seq), min_score, min_len)
    reverse_al = ssw.align(str(query.seq.reverse_complement()), min_score, min_len)

    # Decision tree to return the best aligned sequence taking into acount the absence of result
    # by ssw_wrap in case of score filtering

    if not forward_al:
        if not reverse_al:
            return (None, None)
        else:
            return (reverse_al, False)

    else:
        if not reverse_al:
            return (forward_al, True)
        else:
            if forward_al.score >= reverse_al.score:
                return (forward_al, True)
            else:
                return (reverse_al, False)

def count_seq (filename, seq_type="fasta", gziped=False):
    """
    Count the number of sequences in a fastq or a fastq file
    @param filename Path to a valid readeable file
    @param file_type Should be either fastq or fastq. Default fasta
    @param gziped Boolean indicating if the file is gziped or not. Default False
    """
    #Standard library import
    import gzip
    from mmap import mmap

    # Verify if the file is fasta or fastq type
    assert seq_type in ["fasta", "fastq"], "The file has to be either fastq or fasta format"

    # Open the file
    if gziped:
        f = gzip.open(filename, "r")
    else:
        f = open(filename, "r")

    # FASTA Find a start line seq character ">" an increment the counter each time
    if seq_type ==  "fasta":
        nline = 0
        for line in f:
            if line[0] == ">":
                nline+=1
        f.close()
        return nline

    # FASTQ No motif to find, but 4 lines correspond to 1 sequence
    else:
        nline = 0
        for line in f:
            nline+=1
        f.close()
        return nline/4

def optparser():

    print("Parse command line options")
    # Usage and version strings
    program_name = "pyssw"
    program_version = 0.1
    version_string = "{}\t{}".format(program_name, program_version)
    usage_string = "{}.py -s subject.fasta -q fastq (or fasta) [Facultative options]".format(program_name)
    optparser = optparse.OptionParser(usage = usage_string, version = version_string)

    # Define optparser options
    hstr = "Path of the fasta file containing the subject genome sequence. Can be gziped. [REQUIRED] "
    optparser.add_option( '-s', '--subject', dest="subject", help=hstr)
    hstr = "Path of the fastq or fasta file containing the short read to be aligned. Can be gziped. [REQUIRED]"
    optparser.add_option( '-q', '--query', dest="query", help=hstr)
    hstr = "Type of the query file = fastq or fasta. [default: fastq]"
    optparser.add_option( '-t', '--qtype', dest="qtype", default="fastq", help=hstr)
    hstr = "Positive integer for weight match in genome sequence alignment. [default: 2]"
    optparser.add_option( '-m', '--match', dest="match",default=2, help=hstr)
    hstr = "Positive integer. The negative value will be used as weight mismatch in genome sequence alignment. [default: 2]"
    optparser.add_option( '-x', '--mismatch', dest="mismatch", default=2, help=hstr)
    hstr = "Positive integer. The negative value will be used as weight for the gap opening. [default: 3]"
    optparser.add_option( '-o', '--gap_open', dest="gap_open", default=3, help=hstr)
    hstr = "Positive integer. The negative value will be used as weight for the gap opening. [default: 1]"
    optparser.add_option( '-e', '--gap_extend', dest="gap_extend", default=1, help=hstr)
    hstr = "Integer. Consider alignments having a score <= as not aligned. [default: 0]"
    optparser.add_option( '-f', '--min_score', dest="min_score", default=0, help=hstr)
    hstr = "Integer. Consider alignments having a length <= as not aligned. [default: 0]"
    optparser.add_option( '-l', '--min_len', dest="min_len", default=0, help=hstr)
    hstr = "Flag. Align query in forward and reverse orientation and choose the best alignment. [Set by default]"
    optparser.add_option( '-r', '--reverse', dest="reverse", action="store_true", default=True, help=hstr)
    hstr = "Flag. Write unaligned reads in sam output [Unset by default]"
    optparser.add_option( '-u', '--unaligned', dest="unaligned", action="store_true", default=False, help=hstr)

    # Parse arg and return a dictionnary_like object of options
    opt, args = optparser.parse_args()

    if not opt.subject:
        print ("\nERROR: a subject fasta file has to be provided (-s option)\n")
        optparser.print_help()
        sys.exit()

    if not opt.query:
        print ("\nERROR: a query fasta or fastq file has to be provided (-q option)\n")
        optparser.print_help()
        sys.exit()

    return opt

#~~~~~~~TOP LEVEL INSTRUCTIONS~~~~~~~#

if __name__ == '__main__':

    # try to import Third party and local packages
    try:
        from Bio import SeqIO
    except ImportError:
        print ("ERROR: Please install Biopython package")
        sys.exit()

    try:
        from ssw_wrap import Aligner
    except ImportError:
        print ("ERROR: Please place ssw_wrap in the current directory or add its dir to python path")
        sys.exit()

    # Parse command line arguments
    opt = optparser()
    # Run the main function
    main(opt)
