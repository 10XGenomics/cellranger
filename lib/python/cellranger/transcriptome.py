import bisect
import numpy as np
import os

import cellranger.constants as cr_constants

def make_tx_field(hit):
    if hit.pos is not None:
        fields = [hit.transcript, hit.strand + str(hit.pos), hit.cigarstring]
    else:
        fields = [hit.transcript, '', '']
    return ','.join(fields)

def make_tx_tag(hits):
    return ';'.join([make_tx_field(hit) for hit in hits])

def make_gx_tag(gene_hits):
    gx = ';'.join(gene_hits.keys())
    gn = ';'.join(gene_hits.values())
    return gx, gn

def make_region_tag(region):
    if region == cr_constants.EXONIC_REGION: return 'E'
    elif region == cr_constants.INTRONIC_REGION: return 'N'
    else: return 'I' # intergenic

def make_dummy_alignment(tx_id):
    return cr_constants.TranscriptAlignment(tx_id, None, None, None, None)

def align_to_transcriptome(read, chrom_info, transcript_info, exon_info, params):
    chrom_names, chrom_starts, chrom_bins = chrom_info
    transcript_ids, transcript_starts, transcript_ends, transcript_max_ends, transcript_strands, transcript_num_exons, transcript_exon_break_idx = transcript_info
    exon_breaks, exon_cum_lengths = exon_info

    if read.is_unmapped:
        return {}, False, None

    ref_offset = chrom_starts[read.tid]
    clipped_read_start = ref_offset + read.pos
    clipped_read_end = clipped_read_start + read.alen - 1
    tx = bisect.bisect(transcript_starts, clipped_read_end) - 1

    if tx == -1 or tx == transcript_starts.size:
        # read is at the extreme start / end of transcriptome
        return {}, False, cr_constants.INTERGENIC_REGION

    tx_hits = {}
    any_antisense = False
    any_exonic = False
    any_intronic = False
    # cycle back through overlapping transcripts
    while tx >= 0 and clipped_read_start <= max(transcript_max_ends[tx], transcript_ends[tx]):
        if clipped_read_start <= transcript_ends[tx]: # transcript overlaps the read
            alignment, stats = align_to_transcript(read, ref_offset, transcript_info, exon_info, tx, params)
            antisense, exonic, intronic = stats
            if antisense: any_antisense = True
            if exonic: any_exonic = True
            elif intronic: any_intronic = True
            if alignment is not None:
                tx_hits[alignment.transcript] = alignment
        tx -= 1

    if any_exonic: region = cr_constants.EXONIC_REGION
    elif any_intronic: region = cr_constants.INTRONIC_REGION
    else: region = cr_constants.INTERGENIC_REGION

    return tx_hits, any_antisense, region

def align_to_transcript(read, ref_offset, transcript_info, exon_info, tx, params):
    transcript_ids, transcript_starts, transcript_ends, transcript_max_ends, transcript_strands, transcript_num_exons, transcript_exon_break_idx = transcript_info
    exon_breaks, exon_cum_lengths = exon_info

    # load params
    # TODO turn these into explicit arguments?
    chemistry_negative_strand = params.get('chemistry_negative_strand', False)
    allow_paired_antisense = params.get('allow_paired_antisense', False)
    intergenic_trim_bases = params.get('intergenic_trim_bases', 0)
    intronic_trim_bases = params.get('intronic_trim_bases', 0)
    junction_trim_bases = params.get('junction_trim_bases', 0)
    region_min_overlap = params.get('region_min_overlap', 0.5)

    # lookup transcript info
    tx_start = transcript_starts[tx]
    tx_end = transcript_ends[tx]
    tx_first_exon = transcript_exon_break_idx[tx]
    tx_num_exons = transcript_num_exons[tx]
    tx_last_exon = tx_first_exon + tx_num_exons - 1 # inclusive
    star_genomic_start = ref_offset + read.pos
    star_genomic_end = star_genomic_start + read.alen - 1 # inclusive
    clipped_read_start = star_genomic_start - tx_start # convert to transcript coordinates
    clipped_read_end = star_genomic_end - tx_start

    # compute region
    left_clip, right_clip, aligned_segments = get_cigar_segments(read.cigar, clipped_read_start)
    is_exonic = is_read_exonic(aligned_segments, exon_breaks, tx_first_exon, tx_last_exon, region_min_overlap)
    is_intronic = not is_exonic and overlaps(star_genomic_start, star_genomic_end, tx_start, tx_end, 1.0)

    # discard non-exonic reads
    if not is_exonic:
        return None, (False, is_exonic, is_intronic)

    # compute strand
    tx_reverse_strand = transcript_strands[tx] == 2
    read_reverse_strand = read.is_reverse
    if read.is_paired and read.is_read2:
        read_reverse_strand = not read_reverse_strand
    if chemistry_negative_strand:
        read_reverse_strand = not read_reverse_strand
    is_antisense = tx_reverse_strand != read_reverse_strand

    stats = (is_antisense, is_exonic, is_intronic)

    # discard antisense reads
    if is_antisense and not (read.is_paired and allow_paired_antisense):
        return None, stats

    # find the beginning / ending exons
    ex_bounds = find_exons(exon_breaks, clipped_read_start, clipped_read_end, tx_first_exon, tx_last_exon,
                           intergenic_trim_bases=intergenic_trim_bases, intronic_trim_bases=intronic_trim_bases)

    # discard exon-incompatible reads
    if ex_bounds is None:
        return None, stats
    else:
        ex_start, ex_end = ex_bounds

    # compute offsets
    ex_offset = max(clipped_read_start - exon_breaks[ex_start], 0)
    tx_offset = exon_cum_lengths[ex_start / 2] + ex_offset
    tx_last_exon_length = exon_breaks[2*tx_last_exon + 1] - exon_breaks[2*tx_last_exon] + 1
    tx_length = exon_cum_lengths[tx_last_exon] + tx_last_exon_length

    # align the read to the exons
    tx_alignment = align_junctions(left_clip, right_clip, aligned_segments, exon_breaks[ex_start:ex_end+2], tolerance=junction_trim_bases)

    # discard junction-incompatible reads
    if tx_alignment is None:
        return None, stats
    else:
        tx_cigar, tx_aligned_bases = tx_alignment

    if tx_reverse_strand:
        tx_offset = tx_length - (tx_offset + tx_aligned_bases)
        tx_cigar = reversed(tx_cigar)

    tx_cigar = make_cigar_string(tx_cigar)

    # return STAR-style transcript alignment
    tx_chrom = transcript_ids[tx]
    tx_strand = cr_constants.REVERSE_STRAND if is_antisense else cr_constants.FORWARD_STRAND

    return cr_constants.TranscriptAlignment(tx_chrom, tx_strand, int(tx_offset), tx_cigar, int(tx_aligned_bases)), stats

def find_exons(exon_breaks, read_start, read_end, first_exon, last_exon, intergenic_trim_bases=0, intronic_trim_bases=0):
    # NOTE1: last_exon is inclusive, but `hi` should be exclusive, so add 1
    # NOTE2: remember that exon_breaks has the starts AND ends of each exon, so multiply exon index by 2
    # NOTE3: bisect returns the insertion point, so always subtract 1 to get the exon boundary that sits to the left of that insertion point
    ex_start = bisect.bisect(exon_breaks, read_start, lo=2*first_exon, hi=2*(last_exon+1)) - 1
    ex_end = bisect.bisect_left(exon_breaks, read_end, lo=max(ex_start,0), hi=2*(last_exon+1)) - 1
    start_overhang = 0
    end_overhang = 0

    if ex_start > 2*last_exon or ex_end < 2*first_exon:
        # read is totally out of bounds
        return None

    if ex_start % 2 == 1:
        # read overhangs exon on the left
        intergenic = ex_start < 2*first_exon
        start_overhang = exon_breaks[ex_start+1] - read_start
        trim_bases = intergenic_trim_bases if intergenic else intronic_trim_bases
        if start_overhang <= trim_bases:
            ex_start += 1
        else:
            return None

    if ex_end % 2 == 1:
        # read overhangs exon on the right
        intergenic = ex_end > 2*last_exon
        end_overhang = read_end - exon_breaks[ex_end]
        trim_bases = intergenic_trim_bases if intergenic else intronic_trim_bases
        if end_overhang <= trim_bases:
            ex_end -= 1
        else:
            return None

    return ex_start, ex_end

def is_read_exonic(aligned_segments, exon_breaks, first_exon, last_exon, min_overlap_frac=0.5):
    for read_start, read_end, _ in aligned_segments:
        ex_idx = bisect.bisect(exon_breaks, read_start, lo=2*first_exon, hi=2*(last_exon+1)) - 1
        if ex_idx > 2*last_exon:
            # read is out of bounds
            return False
        elif ex_idx % 2 == 1:
            # read starts in an intron
            exon_start = exon_breaks[ex_idx+1]
            exon_end = exon_breaks[ex_idx+2]
        else:
            # read starts in an exon
            exon_start = exon_breaks[ex_idx]
            exon_end = exon_breaks[ex_idx+1]
        segment_exonic = overlaps(read_start, read_end, exon_start, exon_end, min_overlap_frac)
        if not segment_exonic:
            return False
    return True

def overlaps(read_start, read_end, ref_start, ref_end, min_frac):
    overlap_bases = min(ref_end, read_end) - max(ref_start, read_start)
    return float(overlap_bases) / float(read_end - read_start) >= min_frac

def align_junctions(left_clip, right_clip, aligned_segments, exon_breaks, tolerance=0):
    num_exons = len(exon_breaks) / 2
    if len(aligned_segments) != num_exons:
        return None
    aligned_bases = 0
    full_cigar = left_clip
    for i in xrange(num_exons):
        read_start, read_end, read_cigar = aligned_segments[i]
        # exon_break ends are inclusive, so add 1 to make it match read_end
        exon_start, exon_end = exon_breaks[2*i], exon_breaks[2*i + 1] + 1
        aligned_bases += exon_end - exon_start
        new_cigar = read_cigar
        start_diff = exon_start - read_start
        if i == 0: # first segment
            if start_diff > 0: # overhang -> softclip
                new_cigar = mask_read_bases(new_cigar, start_diff, 4, False)
            elif start_diff < 0: # underhang -> fewer aligned bases
                aligned_bases += start_diff
        else:
            if abs(start_diff) > tolerance:
                return None
            elif start_diff > 0: # overhang -> mark as insertion
                new_cigar = mask_read_bases(new_cigar, start_diff, 1, False)
            elif start_diff < 0: # underhang -> mark as deletion
                new_cigar = mark_deleted_ref_bases(new_cigar, -start_diff, False)
        end_diff = read_end - exon_end
        if i == num_exons - 1: # last segment
            if end_diff > 0: # overhang -> softclip
                new_cigar = mask_read_bases(new_cigar, end_diff, 4, True)
            elif end_diff < 0: # underhang -> fewer aligned bases
                aligned_bases += end_diff
        else:
            if abs(end_diff) > tolerance:
                return None
            elif end_diff > 0: # overhang -> insertion
                new_cigar = mask_read_bases(new_cigar, end_diff, 1, True)
            elif end_diff < 0: # underhang -> deletion
                new_cigar = mark_deleted_ref_bases(new_cigar, -end_diff, True)
        extend_cigar(full_cigar, new_cigar)
    extend_cigar(full_cigar, right_clip)

    return full_cigar, aligned_bases

def extend_cigar(old_cigar, new_cigar):
    # extend list of cigar ops, checking the ends to see if they should be merged
    if len(new_cigar) > 0:
        if len(old_cigar) > 0 and old_cigar[-1][0] == new_cigar[0][0]:
            old_cigar[-1] = (old_cigar[-1][0], old_cigar[-1][1] + new_cigar[0][1])
            old_cigar.extend(new_cigar[1:])
        else:
            old_cigar.extend(new_cigar)

def mark_deleted_ref_bases(cigar, del_len, reverse):
    insert_point = 0 if reverse else len(cigar)
    new_cigar = cigar
    new_cigar.insert(insert_point, (2, del_len))
    return new_cigar

def mask_read_bases(cigar, mask_len, mask_op, reverse):
    # note: this assumes that refskips have been removed
    masked_cigar = []
    masked_cigar.append((mask_op, mask_len))
    consumed_bases = 0
    if reverse: cigar = cigar[::-1]
    for (op, length) in cigar:
        if consumed_bases < mask_len:
            read_bases = length if op != 2 else 0 # deletions don't consume read bases
            if consumed_bases + read_bases >= mask_len:
                retain_bases = read_bases + consumed_bases - mask_len
                masked_cigar.append((op, retain_bases))
            consumed_bases += read_bases
        else:
            masked_cigar.append((op, length))
    if reverse: masked_cigar = masked_cigar[::-1]
    return masked_cigar

def get_cigar_segments(cigar, alignment_start):
    i = 0
    j = -1

    # get clipped blocks
    left_clip = []
    right_clip = []
    if cigar[i][0] == 5: # hardclip
        left_clip.append(cigar[i])
        i += 1
    if cigar[i][0] == 4: # softclip
        left_clip.append(cigar[i])
        i += 1
    if cigar[j][0] == 5: # hardclip
        right_clip.insert(0, cigar[j])
        j -= 1
    if cigar[j][0] == 4: # softclip
        right_clip.insert(0, cigar[j])
        j -= 1

    # get aligned blocks
    aligned_segments = []
    curr_start = alignment_start
    curr_end = curr_start
    curr_cigar = []
    for (op, length) in cigar[i:len(cigar)+j+1]:
        if op == 3: # splice junction -> start a new segment
            aligned_segments.append((curr_start, curr_end, curr_cigar))
            curr_end += length
            curr_start = curr_end
            curr_cigar = []
        else: # continue the current segment
            ref_bases = length if op != 1 else 0 # insertions don't consume ref bases
            curr_end += ref_bases
            curr_cigar.append((op, length))
    aligned_segments.append((curr_start, curr_end, curr_cigar))

    return left_clip, right_clip, aligned_segments

def make_cigar_string(cigar_tuples):
    names = 'MIDNSHP=X'
    return ''.join(['%d%s' % (length, names[op]) for op, length in cigar_tuples])

def build_star_index(star_index_dir):
    # NOTE: as of v2.5.1, STAR stores all chromosomes in one array, with cumulative indexing
    chrom_name_path = os.path.join(star_index_dir, "chrName.txt")
    chrom_start_path = os.path.join(star_index_dir, "chrStart.txt")

    tx_path = os.path.join(star_index_dir, "transcriptInfo.tab")
    ex_path = os.path.join(star_index_dir, "exonInfo.tab")
    #gx_path = os.path.join(star_index_dir, "geneInfo.tab")
    #id_path = os.path.join(star_index_dir, "exonGeTrInfo.tab")

    tx_dtype = {'names': ('tx_ids', 'tx_starts', 'tx_ends', 'tx_max_ends', 'tx_strands', 'tx_num_exons', 'tx_break_idxs'), 'formats': ('object', 'u8', 'u8', 'u8', 'u1', 'u2', 'u4')}
    ex_dtype = {'names': ('ex_starts', 'ex_ends', 'ex_cum_lengths'), 'formats': ('u8', 'u8', 'u8')}
    #gx_dtype = {'names': ('gx_ids',), 'formats': ('object',)}
    #id_dtype = {'names': ('gx_idx', 'tx_idx'), 'formats': ('u4', 'u4')}

    chrom_names = np.loadtxt(chrom_name_path, dtype='object')
    chrom_starts = np.loadtxt(chrom_start_path, dtype='u8')

    tx_info = np.loadtxt(tx_path, tx_dtype, skiprows=1, unpack=True)
    ex_info = np.loadtxt(ex_path, ex_dtype, skiprows=1, unpack=True)
    #gx_info = np.loadtxt(gx_path, gx_dtype, skiprows=1, unpack=True)
    #id_info = np.loadtxt(id_path, id_dtype, skiprows=1, unpack=True, usecols=(3,4))

    # compute index range of each chrom within tx_starts
    tx_starts = tx_info[1]
    chrom_bins = [bisect.bisect_left(tx_starts, cs) for cs in chrom_starts]
    chrom_info = chrom_names, chrom_starts, chrom_bins

    # interleave exon starts/ends
    ex_starts, ex_ends, ex_cum_lengths = ex_info
    ex_breaks = np.empty((2 * ex_starts.size,), dtype=ex_starts.dtype)
    ex_breaks[0::2] = ex_starts
    ex_breaks[1::2] = ex_ends
    ex_info = ex_breaks, ex_cum_lengths

    return chrom_info, tx_info, ex_info
