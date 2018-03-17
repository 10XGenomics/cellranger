#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#
# Given a set of FASTQs, infer 3' vs 5' chemistry.
# If a V(D)J reference is also given, try to distinguish 5' unbiased gene exprssion vs VDJ-enriched.

import martian
import json
import os
import tenkit.log_subprocess as tk_subproc
import tenkit.stats as tk_stats
import cellranger.chemistry as cr_chem
import cellranger.constants as cr_constants
import cellranger.fastq as cr_fastq
import cellranger.preflight as cr_preflight
import cellranger.utils as cr_utils
import cellranger.vdj.preflight as vdj_preflight
import cellranger.vdj.reference as vdj_ref
import enum

__MRO__ = """
stage DETECT_CHEMISTRY(
    in  string sample_id,
    in  map[]  sample_def,
    in  path   reference_path,
    in  path   vdj_reference_path,
    in  string chemistry_name_spec,
    in  string[] allowed_chems,
    out json   summary,
    out string chemistry_type,
    out txt    report,
    src py     "stages/chemistry_detector/detect_chemistry",
) split using (
)
"""

BCL_PROCESSOR_READ_TYPE = 'RA'
READ_TYPES = ['R1', 'R2']

# Number of initial reads to inspect in each FASTQ
INITIAL_READS = 10000

# Min reads to make a call whether reads are UNMAPPED
MIN_MAPPED_READS = 1000

# If the fraction of reads mapped ambiguously is more than this, we call it AMBIGUOUS_MAPPED
MAX_AMBIGUOUS_MAPPED_READ_FRAC = 0.75

# Min fraction of reads going to V(D)J ref to call as V(D)J
MIN_VDJ_READ_FRAC = 0.25

class ReadState(enum.Enum):
    SENSE_MAPPED = 0
    ANTISENSE_MAPPED = 1
    AMBIGUOUS_MAPPED = 2
    UNMAPPED = 3
    ABSENT = 4

def getReadState(sense_count, anti_count, mapped_count, total_count):

    # Check if the read is present
    if total_count == 0:
        return ReadState.ABSENT

    # Check if it is unmapped
    if mapped_count < MIN_MAPPED_READS:
        return ReadState.UNMAPPED

    ambiguous_count = mapped_count - (sense_count + anti_count)
    if (ambiguous_count > (MAX_AMBIGUOUS_MAPPED_READ_FRAC * mapped_count)) or (sense_count==anti_count):
        return ReadState.AMBIGUOUS_MAPPED

    return ReadState.SENSE_MAPPED if sense_count > anti_count else ReadState.ANTISENSE_MAPPED


def is_vdj(gex_count, vdj_count):
    return tk_stats.robust_divide(vdj_count, gex_count) >= MIN_VDJ_READ_FRAC

def run(args):
    """ Run tk_subproc.check_call and print command """
    print ' '.join(args)
    tk_subproc.check_call(args)

def map_reads(fq_spec_path, idx_path, out_path):
    run(['detect_chemistry', 'map-reads',
         idx_path, '/dev/null', out_path,
         '--initial-reads', str(INITIAL_READS),
         '--fastqs', fq_spec_path])

    # Get output
    with open(out_path) as f:
        return json.load(f)


def find_fastqs(sample_defs):
    chunks = []

    for sd in sample_defs:
        fq_spec = cr_fastq.FastqSpec.from_sample_def(sd)

        fastqs = {r: fq_spec.get_fastqs(r) for r in READ_TYPES}

        # No FASTQs found for a single sample def
        if len(fastqs) == 0:
            martian.exit(cr_constants.NO_INPUT_FASTQS_MESSAGE)

        for read_type, read_fastqs in sorted(fastqs.iteritems()):
            chunks.append({'read_type': read_type,
                           'input': read_fastqs,
                           'interleaved': fq_spec.interleaved,
            })

    # No FASTQs found at all
    if len(chunks) == 0:
        martian.exit(cr_constants.NO_INPUT_FASTQS_MESSAGE)

    return chunks

def run_preflight_checks(args):
    print "Checking sample info..."
    cr_preflight.check_sample_def(args.sample_def)

    # Don't check the GEX reference if we are just distinguishing between SCVDJ seq configs
    if 'SCVDJ' not in args.chemistry_name_spec:
        print "Checking reference..."
        cr_preflight.check_refdata(args.reference_path)

    if args.vdj_reference_path is not None:
        print "Checking V(D)J reference..."
        vdj_preflight.check_refdata(args.vdj_reference_path, True)

def prepare_transcriptome_indexes(reference_path, vdj_reference_path):
    """ Use ReadStates of R1/R2 to determine SC3Pv1 vs SC3Pv2 vs SC5P-R1 vs SC5P_auto/SCVDJ.
        Returns (chemistry_name, report, metrics)
        where report is a text report and metrics is a dict """

    ## Index the reference fasta
    fa_path = os.path.join(reference_path, cr_constants.REFERENCE_FASTA_PATH)
    new_fa_path = martian.make_path('ref.fa')

    need_index = True

    if os.path.exists(fa_path + '.fai'):
        # Look for existing .fai file (won't exist for our standard ref packages)
        martian.update_progress('Found genome FASTA index....')
        new_fa_path = fa_path
        need_index = False

    else:
        # Note: this will fail if user's fs doesn't support symlinks
        martian.update_progress('Symlinking genome FASTA...')
        os.symlink(fa_path, new_fa_path)

    if need_index:
        martian.update_progress('Indexing genome...')
        run(['samtools', 'faidx', new_fa_path])


    ## Generate a transcriptome reference from a genome ref
    martian.update_progress('Building transcriptome...')
    gtf_path = os.path.join(reference_path, cr_constants.REFERENCE_GENES_GTF_PATH)
    out_fa_path = martian.make_path('transcriptome.fa')
    # Only index the 1st encountered transcript per gene
    run(['detect_chemistry', 'get-transcripts', new_fa_path, gtf_path, out_fa_path])


    ## Build kmer index
    martian.update_progress('Building kmer index...')
    kmer_idx_path = martian.make_path('kmers.idx')

    ## Use a larger step size as the reference grows.
    ## This ensure the index size stays sane.  
    ## Should get to a step of <10 for the whole genome, which 
    ## is still 3x overlap w/ 32-mers
    fa_size = os.path.getsize(os.path.realpath(out_fa_path))
    step = fa_size / 400000000
    skip = step - 1

    index_args = ['detect_chemistry', 'index-transcripts']
    if skip > 0:
        index_args.append('--skip=%d' % skip)

    index_args.extend([out_fa_path, kmer_idx_path])
    run(index_args)

    # Build VDJ kmer index (optional)
    vdj_idx_path = None
    if vdj_reference_path is not None:
        vdj_fa_path = vdj_ref.get_vdj_reference_fasta(vdj_reference_path)
        vdj_idx_path = martian.make_path('vdj_kmers.idx')
        run(['detect_chemistry', 'index-transcripts', vdj_fa_path, vdj_idx_path])

    return (kmer_idx_path, vdj_idx_path)


def infer_sc3p_or_sc5p(chunks, kmer_idx_path, vdj_idx_path):
    """ Use ReadStates of R1/R2 to determine SC3Pv1 vs SC3Pv2 vs SC5P-R1 vs SC5P_auto/SCVDJ.
        Returns (chemistry_name, report, metrics)
        where report is a text report and metrics is a dict """

    ## Map read kmers to each strand
    martian.update_progress('Mapping reads...')

    # Prepare fastq paths
    fq_specs = {}
    for read_type in READ_TYPES:
        fq_specs[read_type] = martian.make_path('%s_in.json' % read_type)
        with open(fq_specs[read_type], 'w') as f:
            json.dump([c for c in chunks if c['read_type'] == read_type], f)

    # Map reads to gene expression reference
    metrics = {}
    for read_type in READ_TYPES:
        map_out_path = martian.make_path('%s_out.json' % read_type)
        metrics[read_type] = map_reads(fq_specs[read_type], kmer_idx_path, map_out_path)

    # Map reads to VDJ reference (optional)
    if vdj_idx_path is not None:
        for read_type in READ_TYPES:
            map_out_path = martian.make_path('vdj_%s_out.json' % read_type)
            vdj_metrics = map_reads(fq_specs[read_type], vdj_idx_path, map_out_path)

            for k,v in vdj_metrics.iteritems():
                metrics[read_type]['vdj_' + k] = v

    # Verify total read counts
    r1_total = metrics['R1']['total_reads']
    r2_total = metrics['R2']['total_reads']
    if r1_total != 0 and r2_total != 0 and r1_total != r2_total:
        martian.exit('Total read counts for R1 and R2 must be identical if both are present. There were %d R1 reads and %d R2 reads. Check that all of the FASTQ files are present.' % (r1_total, r2_total))


    ## Infer chemistry
    report = '\n'
    for read_type, m in metrics.iteritems():
        report += '%s Total Reads:     %s\n' % (read_type, str(m['total_reads']).rjust(20))
        report += '%s Sense Reads:     %s\n' % (read_type, str(m['sense_reads']).rjust(20))
        report += '%s Antisense Reads: %s\n' % (read_type, str(m['antisense_reads']).rjust(20))

    if vdj_idx_path is not None:
        for read_type, m in metrics.iteritems():
            report += '%s Sense V(D)J Reads:     %s\n' % (read_type, str(m['vdj_sense_reads']).rjust(20))
            report += '%s Antisense V(D)J Reads: %s\n' % (read_type, str(m['vdj_antisense_reads']).rjust(20))

    r1_state = getReadState(metrics['R1']['sense_reads'], metrics['R1']['antisense_reads'], metrics['R1']['mapped_reads'], metrics['R1']['total_reads'])
    r2_state = getReadState(metrics['R2']['sense_reads'], metrics['R2']['antisense_reads'], metrics['R2']['mapped_reads'], metrics['R2']['total_reads'])

    report += "\n"
    chemistry_name = None

    if (r1_state == ReadState.SENSE_MAPPED) and (r2_state == ReadState.UNMAPPED):
        chemistry_name = 'SC3Pv1'
        report += "This library is likely to be a Single Cell 3' gene expression library (v1)."

    elif r2_state == ReadState.SENSE_MAPPED:
        chemistry_name = 'SC3Pv2'
        report += "This library is likely to be a Single Cell 3' gene expression library (v2)."

    elif (r1_state == ReadState.SENSE_MAPPED) and (r2_state == ReadState.ABSENT):
        chemistry_name = 'SC5P-R1'
        report += "This library is likely to be a Single Cell 5' gene expression library (R1)."

    elif (r2_state == ReadState.ANTISENSE_MAPPED):
        r1_gex_sense_count = metrics['R1'].get('sense_reads', 0)
        r2_gex_anti_count = metrics['R2'].get('antisense_reads', 0)
        r1_vdj_sense_count = metrics['R1'].get('vdj_sense_reads', 0)
        r2_vdj_anti_count = metrics['R2'].get('vdj_antisense_reads', 0)
        if vdj_idx_path is None:
            report += "This library is likely to be a Single Cell V(D)J or Single Cell 5' gene expression library."
            chemistry_name = 'SC5P_auto'
        else:
            if (is_vdj(r1_gex_sense_count, r1_vdj_sense_count) or is_vdj(r2_gex_anti_count, r2_vdj_anti_count)):
                report += "This library is likely to be a Single Cell V(D)J library."
                chemistry_name = 'SCVDJ'

            else:
                report += "This library is likely to be a Single Cell 5' gene expression library."
                chemistry_name = 'SC5P_auto'

    else:
        report += "There was not enough information to determine the nature of the library."
        chemistry_name = None

    return chemistry_name, report, metrics


# The only reason we split is to allow VDR to kill intermediate files
def split(args):
    return {'chunks': [{'__mem_gb': 8}]}

def main(args, outs):
    # Check chemistry restrictions
    if args.allowed_chems is not None and \
       args.chemistry_name_spec not in args.allowed_chems:
        martian.exit("The chemistry name '%s' is not allowed for this pipeline. The allowed values are: %s" % (args.chemistry_name_spec, ', '.join(args.allowed_chems)))

    ## If chem explicitly specified, just check it and finish
    if args.chemistry_name_spec not in cr_chem.AUTO_CHEMISTRY_NAMES or \
       args.chemistry_name_spec == cr_chem.CUSTOM_CHEMISTRY_NAME:
        ok, msg = cr_chem.check_chemistry_arg(args.chemistry_name_spec)
        if not ok:
            martian.exit(msg)

        # Write empty json
        with open(outs.summary, 'w') as f:
            json.dump({}, f)

        outs.chemistry_type = args.chemistry_name_spec
        outs.report = None
        return

    ## Run preflight checks
    try:
        run_preflight_checks(args)
    except cr_preflight.PreflightException as e:
        martian.exit(e.msg)

    ## Find the input fastqs
    chunks = find_fastqs(args.sample_def)

    chemistry_name = args.chemistry_name_spec
    report = ''
    metrics = {}

    if args.chemistry_name_spec == 'auto':
        (txome_idx, vdj_idx) = prepare_transcriptome_indexes(args.reference_path, args.vdj_reference_path)

        auto_chemistries = {}
        for (idx, sd) in enumerate(args.sample_def):
            chunks = find_fastqs([sd])
            chemistry_name, report, metrics = infer_sc3p_or_sc5p(chunks, txome_idx, vdj_idx)
            auto_chemistries[idx] = chemistry_name

        if len(set(auto_chemistries.itervalues())) > 1:
            c = ', '.join(set(auto_chemistries.itervalues()))
            s = '\n'.join("  Sample def %d: %s" % (idx, chem) for (idx, chem) in sorted(auto_chemistries.iteritems()))
            martian.exit("Detected conflicting chemistry types (%s).\n Please run these data separately.\n%s" % (c, s))

        else:
            chemistry_name = auto_chemistries[0]


    # Further refinement:
    #   - Detect the sequencing configuration for SC5P (SC5P-PE vs SC5P-R2)
    #   - Detect the sequencing configuration for SCVDJ (SCVDJ vs SCVDJ-R2)
    #
    # The chemistry/seq-config must be consistent across all sample defs
    if chemistry_name in cr_chem.AUTO_CHEMISTRY_NAMES:
        # Map (sample_def_idx, fastq_group_name) => chemistry_name
        group_chem = {}

        for sd_idx, sd in enumerate(args.sample_def):
            fq_spec = cr_fastq.FastqSpec.from_sample_def(sd)

            # Infer chemistry for each sample index/name (aka fastq group)
            for group, group_spec in fq_spec.get_group_spec_iter():
                try:
                    group_chem[(sd_idx, group)] = cr_chem.infer_chemistry(chemistry_name, group_spec)

                except cr_chem.NoInputFastqsException:
                    # It's okay for a single sample index/name to be absent
                    continue

        if len(group_chem) == 0:
            martian.exit(cr_constants.NO_INPUT_FASTQS_MESSAGE)

        martian.log_info("Detected chemistries:")
        for (i, g) in group_chem.iteritems():
            martian.log_info("%s: %s" % (str(i), str(g)))

        # Check for multiple chemistry types
        if len(set(group_chem.itervalues())) > 1:
            c = ', '.join(map(str, sorted(list(set(group_chem.itervalues())))))
            s = ', '.join("Sample def %d/%s: %s" % (i,g,v) for ((i,g),v) in sorted(group_chem.iteritems()))
            martian.exit("Detected conflicting chemistry types (%s). Please run these data separately. %s" % (c, s))

        chemistry_name = group_chem.values()[0]

        report += "\nThe chemistry version or sequencing configuration is likely %s" % cr_chem.get_chemistry_description_from_name(chemistry_name)

    outs.chemistry_type = chemistry_name

    # Write report file
    martian.log_info(report)

    with open(outs.report, 'w') as f:
        f.write(report + "\n")

    # Write summary JSON
    metrics['chemistry'] = chemistry_name
    with open(outs.summary, 'w') as f:
        json.dump(metrics, f)

def join(args, outs, chunk_defs, chunk_outs):
    cr_utils.copy(chunk_outs[0].summary, outs.summary)
    if chunk_outs[0].report is not None:
        cr_utils.copy(chunk_outs[0].report, outs.report)
    outs.chemistry_type = chunk_outs[0].chemistry_type
