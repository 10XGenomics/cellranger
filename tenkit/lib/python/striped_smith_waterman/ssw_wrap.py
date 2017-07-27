"""
@package ssw_wrap
@brief Simple python wrapper for SSW align library
To use the dynamic library libssw.so you may need to modify the LD_LIBRARY_PATH environment
variable to include the library directory (export LD_LIBRARY_PATH=$PWD) or for definitive
inclusion of the lib edit /etc/ld.so.conf and add the path or the directory containing the
library and update the cache by using /sbin/ldconfig as root
@copyright  [The MIT licence](http://opensource.org/licenses/MIT)
@author     Clement & Adrien Leger - 2014
"""

#~~~~~~~GLOBAL IMPORTS~~~~~~~#
# Standard library packages
import os
from ctypes import *

_so_file_name = os.path.dirname(os.path.abspath(__file__)) + '/' + "libssw.so"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class CAlignRes(Structure):
    """
    @class  SSWAlignRes
    @brief  ctypes Structure with s_align struct mapping returned by SSWAligner.Align func
            Correspond to the structure of the query profile
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~Ctype Structure~~~~~~~#
    _fields_ = [('score', c_uint16),
                ('score2', c_uint16),
                ('ref_begin', c_int32),
                ('ref_end', c_int32),
                ('query_begin', c_int32),
                ('query_end', c_int32),
                ('ref_end2', c_int32),
                ('cigar', POINTER(c_uint32)),
                ('cigarLen', c_int32)]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class Aligner(object):
    """
    @class  SSWAligner
    @brief Wrapper for SSW align library
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~CLASS VARIABLES~~~~~~~#

    # Dictionnary to map Nucleotide to int as expected by the SSW C library
    base_to_int = { 'A':0, 'C':1, 'G':2, 'T':3, 'N':4, 'a':0, 'c':1, 'g':2, 't':3, 'n':4}
    int_to_base = { 0:'A', 1:'C', 2:'G', 3:'T', 4:'N'}

    # Load the ssw library using ctypes
    libssw = cdll.LoadLibrary(_so_file_name)

    # Init and setup the functions pointer to map the one specified in the SSW lib
    # ssw_init method
    ssw_init = libssw.ssw_init
    ssw_init.restype = c_void_p
    ssw_init.argtypes = [POINTER(c_int8), c_int32, POINTER(c_int8), c_int32, c_int8]
    # init_destroy function
    init_destroy = libssw.init_destroy
    init_destroy.restype = None
    init_destroy.argtypes =  [c_void_p]
    # ssw_align function
    ssw_align = libssw.ssw_align
    ssw_align.restype = POINTER(CAlignRes)
    ssw_align.argtypes = [c_void_p, POINTER(c_int8), c_int32, c_uint8, c_uint8, c_uint8, c_uint16, c_int32, c_int32]
    # align_destroy function
    align_destroy = libssw.align_destroy
    align_destroy.restype = None
    align_destroy.argtypes = [POINTER(CAlignRes)]

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __repr__(self):
        msg = self.__str__()
        msg += "SCORE PARAMETERS:\n"
        msg += " Gap Weight     Open: {}     Extension: {}\n".format(-self.gap_open, -self.gap_extend)
        msg += " Align Weight   Match: {}    Mismatch: {}\n\n".format(self.match, -self.mismatch)
        msg += " Match/mismatch Score matrix\n"
        msg += " \tA\tC\tG\tT\tN\n"
        msg += " A\t{}\t{}\t{}\t{}\t{}\n".format(self.match, -self.mismatch, -self.mismatch, -self.mismatch, 0)
        msg += " C\t{}\t{}\t{}\t{}\t{}\n".format(-self.mismatch, self.match, -self.mismatch, -self.mismatch, 0)
        msg += " G\t{}\t{}\t{}\t{}\t{}\n".format(-self.mismatch, -self.mismatch, self.match, -self.mismatch, 0)
        msg += " T\t{}\t{}\t{}\t{}\t{}\n".format(-self.mismatch, -self.mismatch, -self.mismatch, self.match, 0)
        msg += " N\t{}\t{}\t{}\t{}\t{}\n\n".format(0,0,0,0,0)
        msg += "RESULT PARAMETERS:\n"
        msg += " Report cigar           {}\n".format(self.report_cigar)
        msg += " Report secondary match {}\n\n".format(self.report_secondary)
        msg += "REFERENCE SEQUENCE :\n"
        if self.ref_len <= 50:
            msg += "".join([self.int_to_base[i] for i in self.ref_seq])+"\n"
        else:
            msg += "".join([self.int_to_base[self.ref_seq[i]] for i in range(50)])+"...\n"
        msg += " Lenght :{} nucleotides\n".format(self.ref_len)
        return msg

    def __str__(self):
        return "\n<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    def __init__(self,
                ref_seq="",
                match=2,
                mismatch=2,
                gap_open=3,
                gap_extend=1,
                report_secondary=False,
                report_cigar=False):
        """
        Initialize object by creating an interface with ssw library fonctions
        A reference sequence is also assigned to the object for multiple alignment against queries
        with the align function
        @param ref_seq Reference sequence as a python string (case insensitive)
        @param match Weight for a match
        @param mismatch Absolute value of mismatch penalty
        @param gap_open Absolute value of gap open penalty
        @param gap_extend Absolute value of gap extend penalty
        @param report_secondary Report the 2nd best alignement if true
        @param report_cigar Report cigar string if true
        """

        # Store overall alignment parameters
        self.report_secondary = report_secondary
        self.report_cigar = report_cigar

        # Set gap penalties
        self.set_gap(gap_open, gap_extend)

        # Set the cost matrix
        self.set_mat(match, mismatch)

        # Set the reference sequence
        self.set_ref(ref_seq)

    #~~~~~~~SETTERS METHODS~~~~~~~#

    def set_gap(self, gap_open=3, gap_extend=1):
        """
        Store gapopen and gap extension penalties
        """
        self.gap_open = gap_open
        self.gap_extend = gap_extend


    def set_mat(self, match=2, mismatch=2):
        """
        Store match and mismatch scores then initialize a Cost matrix and fill it with match and
        mismatch values. Ambiguous base: no penalty
        """
        self.match = match
        self.mismatch = mismatch

        mat_decl = c_int8 * 25
        self.mat = mat_decl(match, -mismatch, -mismatch, -mismatch, 0,
                            -mismatch, match, -mismatch, -mismatch, 0,
                            -mismatch, -mismatch, match, -mismatch, 0,
                            -mismatch, -mismatch, -mismatch, match, 0,
                            0, 0, 0, 0, 0)

    def set_ref(self, ref_seq):
        """
        Determine the size of the ref sequence and cast it in a c type integer matrix
        """
        if ref_seq:
            self.ref_len = len(ref_seq)
            self.ref_seq = self._DNA_to_int_mat (ref_seq, self.ref_len)
        else:
            self.ref_len = 0
            self.ref_seq = ""

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def align(self, query_seq, min_score=0, min_len=0):
        """
        Perform the alignment of query against the object reference sequence
        @param query_seq Query sequence as a python string (case insensitive)
        @param min_score Minimal score of match. None will be return in case of filtering out
        @param min_len Minimal length of match. None will be return in case of filtering out
        @return A SSWAlignRes Object containing informations about the alignment.
        """
        # Determine the size of the ref sequence and cast it in a c type integer matrix
        query_len = len(query_seq)
        query_seq = self._DNA_to_int_mat (query_seq, query_len)

        # Create the query profile using the query sequence
        profile = self.ssw_init(query_seq, # Query seq in c type integers
                                c_int32(query_len), # Length of Queryseq in bites
                                self.mat, # Score matrix
                                5, # Square root of the number of elements in mat
                                2) # flag = no estimation of the best alignment score

        # Setup the mask_len parameters = distance between the optimal and suboptimal alignment
        # if < 15, the function will NOT return the suboptimal alignment information

        if query_len > 30:
            mask_len = query_len/2
        else:
            mask_len = 15

        c_result = self.ssw_align (profile, # Query profile
                                self.ref_seq, # Ref seq in c type integers
                                c_int32(self.ref_len), # Length of Refseq in bites
                                self.gap_open, # Absolute value of gap open penalty
                                self.gap_extend, # absolute value of gap extend penalty
                                1, # Bitwise FLAG for output values = return all
                                0, # Score filter = return all
                                0, # Distance filter = return all
                                mask_len) # Distance between the optimal and suboptimal alignment

        # Transform the Cstructure into a python object if score and lenght match the requirements
        if c_result and c_result.contents:
            score = c_result.contents.score
            match_len  = c_result.contents.query_end - c_result.contents.query_begin + 1
        else:
            score = -999999999999
            match_len = -10000000000

        if score >= min_score and match_len >= min_len:
            py_result = PyAlignRes(c_result, query_len, self.report_secondary, self.report_cigar)
        else:
            py_result = None

        # Free reserved space by ssw.init and ssw_init methods.
        self._init_destroy(profile)
        if c_result:
            self._align_destroy(c_result)

        # Return the object
        return py_result

    #~~~~~~~PRIVATE METHODS~~~~~~~#

    def _DNA_to_int_mat (self, seq, len_seq):
        """
        Cast a python DNA string into a Ctype int8 matrix
        """
        # Declare the matrix
        query_num_decl = c_int8 * len_seq
        query_num = query_num_decl()

        # for each letters in ATCGN transform in integers thanks to self.base_to_int
        for i in range(len_seq):
            try:
                value = self.base_to_int[seq[i]]
            # if the base is not in the canonic DNA bases assign 4 as for N
            except KeyError:
                value = 4
            finally:
                query_num[i] = value

        return query_num

    def _init_destroy(self, profile):
        """
        Free the space alocated for the matrix used by init
        """
        self.init_destroy(profile)

    def _align_destroy(self, align):
        """
        Free the space alocated for the matrix used by align
        """
        self.align_destroy(align)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class PyAlignRes(object):
    """
    @class  PyAlignRes
    @brief  Extract and verify result from a CAlignRes structure. A comprehensive python
    object is created according to user requirements (+- cigar string and secondary alignment)
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~CLASS VARIABLES~~~~~~~#

    # Load the ssw library using ctypes
    libssw = cdll.LoadLibrary(_so_file_name)

    # Init and setup the functions pointer to map the one specified in the SSW lib
    # cigar_int_to_len function
    cigar_int_to_len = libssw.cigar_int_to_len
    cigar_int_to_len.restype = c_int32
    cigar_int_to_len.argtypes = [c_int32]
    # cigar_int_to_op function
    cigar_int_to_op = libssw.cigar_int_to_op
    cigar_int_to_op.restype = c_char
    cigar_int_to_op.argtypes = [c_int32]

    #~~~~~~~FONDAMENTAL METHOD~~~~~~~#

    def __repr__(self):
        msg = self.__str__()
        msg += "OPTIMAL MATCH\n"
        msg += "Score            {}\n".format(self.score)
        msg += "Reference begin  {}\n".format(self.ref_begin)
        msg += "Reference end    {}\n".format(self.ref_end)
        msg += "Query begin      {}\n".format(self.query_begin)
        msg += "Query end        {}\n".format(self.query_end)

        if self.cigar_string:
            msg += "Cigar_string     {}\n".format(self.cigar_string)

        if self.score2:
            msg += "SUB-OPTIMAL MATCH\n"
            msg += "Score 2           {}\n".format(self.score2)
            msg += "Ref_end2          {}\n".format(self.ref_end2)

        return msg

    def __str__(self):
        return "\n<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)


    def __init__ (self, Res, query_len, report_secondary=False, report_cigar=False):
        """
        Parse CAlignRes structure and copy its values in object variables
        @param Res A CAlignRes structure
        @param query_len length of the query sequence
        @param report_secondary Report the 2nd best alignement if true
        @param report_cigar Report cigar string if true
        """
        # Parse value in the C type structure pointer
        # Minimal mandatory parameters
        self.score = Res.contents.score
        self.ref_begin = Res.contents.ref_begin
        self.ref_end = Res.contents.ref_end
        self.query_begin = Res.contents.query_begin
        self.query_end = Res.contents.query_end

        # Information for sub-optimal match if require and available
        score2 = Res.contents.score2
        if report_secondary and score2 != 0:
            self.score2 = score2
            self.ref_end2 = Res.contents.ref_end2
        else:
            self.score2 = None
            self.ref_end2 = None

        # Cigar Information if CIGAR string if require and available
        cigar_len = Res.contents.cigarLen
        if report_cigar and cigar_len > 0:
            self.cigar_string = self._cigar_string (Res.contents.cigar, cigar_len, query_len)
        else:
            self.cigar_string = None

    #~~~~~~~PRIVATE METHODS~~~~~~~#

    def _cigar_string(self, cigar, cigar_len, query_len):
        """
        Convert cigar and cigarLen into an human readable Cigar string as in SAM files
        """
        # Empty string for iterative writing of the cigar string
        cigar_string = ""

        # If the query match do not start at its first base
        # = introduce a softclip at the begining
        if self.query_begin > 0:
            op_len = self.query_begin
            op_char = "S"
            cigar_string += '{}{}'.format(op_len, op_char)

        # Iterate over the cigar (pointer to a vector of int)
        for i in range(cigar_len):
            op_len = self.cigar_int_to_len(cigar[i])
            op_char = self.cigar_int_to_op(cigar[i])
            cigar_string += '{}{}'.format(op_len, op_char)

        # If the lenght of bases aligned is shorter than the overall query length
        # = introduce a softclip at the end
        end_len = query_len - self.query_end - 1
        if  end_len != 0:
            op_len = end_len
            op_char = "S"
            cigar_string += '{}{}'.format(op_len, op_char)

        return cigar_string
