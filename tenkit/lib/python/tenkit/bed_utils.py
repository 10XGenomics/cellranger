import tenkit.bio_io as tk_io
import tenkit.regions as tk_regions
import shutil

def writeOut(bed_dict, bedOut):
    with open(bedOut, "w") as fout:
        all_chroms = bed_dict.keys()
        all_chroms.sort()
        for chrom in all_chroms:
            for start, end in bed_dict[chrom]:
                fout.write("%s\t%d\t%d\n" % (chrom, start, end))

def interval_subtract(start, end, overlappings):
    final = []
    cur_start = start

    for ss, ee in overlappings:
        if ss > cur_start:
            final.append((cur_start, ss))
        cur_start = ee

    if cur_start < end:
        final.append((cur_start,end))

    return final


def merge(bed1, bed2, bedOut):
    if not bed2:
        shutil.copyfile(bed1,bedOut)
        return

    with open(bed1) as f:
        bed_dict1 = tk_io.get_target_regions(f)

    with open(bed2) as f:
        bed_dict2 = tk_io.get_target_regions(f)

    for chrom in bed_dict2:
        for start, end in bed_dict2[chrom]:
            if chrom not in bed_dict1:
                bed_dict1[chrom]=tk_regions.Regions([])
            bed_dict1[chrom].add_region((start, end))

    writeOut(bed_dict1, bedOut)


def intersect(bed1, bed2, bedOut):
    if not bed2:
        shutil.copyfile(bed1,bedOut)
        return

    with open(bed1) as f:
        bed_dict1 = tk_io.get_target_regions(f)

    with open(bed2) as f:
        bed_dict2 = tk_io.get_target_regions(f)

    all_common_chroms = [chrom for chrom in bed_dict1.keys() if chrom in bed_dict2]
    bed_dict_intersect ={}

    for chrom in all_common_chroms:
        bed_dict_intersect[chrom] = bed_dict1[chrom].intersect(bed_dict2[chrom])

    writeOut(bed_dict_intersect, bedOut)


def overlap(bed1, bed2, bedOut):
    if not bed2:
        shutil.copyfile(bed1,bedOut)
        return

    with open(bed1) as f:
        bed_dict1 = tk_io.get_target_regions(f)

    with open(bed2) as f:
        bed_dict2 = tk_io.get_target_regions(f)

    bed_dict_overlap = {}
    for chrom in bed_dict1:
        if not chrom in bed_dict_overlap:
            bed_dict_overlap[chrom] = tk_regions.Regions([])
        for start, end in bed_dict1[chrom]:
            if chrom in bed_dict2 and \
                bed_dict2[chrom].overlaps_region(start, end):
                bed_dict_overlap[chrom].add_region((start,end))

    writeOut(bed_dict_overlap, bedOut)


def no_overlap(bed1, bed2, bedOut):
    if not bed2:
        shutil.copyfile(bed1,bedOut)
        return

    with open(bed1) as f:
        bed_dict1 = tk_io.get_target_regions(f)

    with open(bed2) as f:
        bed_dict2 = tk_io.get_target_regions(f)

    bed_dict_no_overlap = {}
    for chrom in bed_dict1:
        if not chrom in bed_dict_no_overlap:
            bed_dict_no_overlap[chrom] = tk_regions.Regions([])
        for start, end in bed_dict1[chrom]:
            if not chrom in bed_dict2 or \
                not bed_dict2[chrom].overlaps_region(start, end):
                bed_dict_no_overlap[chrom].add_region((start,end))

    writeOut(bed_dict_no_overlap, bedOut)



def subtract(bed1, bed2, bedOut):
    if not bed2:
        shutil.copyfile(bed1,bedOut)
        return

    with open(bed1) as f:
        bed_dict1 = tk_io.get_target_regions(f)

    with open(bed2) as f:
        bed_dict2 = tk_io.get_target_regions(f)

    bed_dict_subtract = {}
    for chrom in bed_dict1:
        if not chrom in bed_dict_subtract:
            bed_dict_subtract[chrom] = tk_regions.Regions([])
        for start, end in bed_dict1[chrom]:
            overlappings = []
            if chrom in bed_dict2:
                overlappings = bed_dict2[chrom].overlapping_regions(start, end)
            for interval in interval_subtract(start, end, overlappings):
                bed_dict_subtract[chrom].add_region(interval)

    writeOut(bed_dict_subtract, bedOut)
