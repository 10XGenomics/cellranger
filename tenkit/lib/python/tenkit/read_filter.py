#!/usr/bin/env python
#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#

import tenkit.constants
from tenkit.bio_io import get_read_barcode

def stringent_read_filter(bam_read, require_barcode):
    ''' Test for a very high-quality read. Only reads satisfying this predicate are
        used when computing summary dup rates to avoid spurious duplicates.
        Reads must have a high MAPQ, a perfect cigar, with a fragment size somewhat
        longer than the read length '''

    return bam_read.mapq >= tenkit.constants.HIGH_CONF_MAPQ and \
            not bam_read.is_secondary and \
           ((abs(bam_read.pos - bam_read.mpos) >= tenkit.constants.MIN_MATE_OFFSET_DUP_FILTER and \
                   bam_read.tid == bam_read.rnext) or bam_read.tid != bam_read.rnext) and \
           len(bam_read.cigar) == 1 and (not require_barcode or get_read_barcode(bam_read) is not None)
