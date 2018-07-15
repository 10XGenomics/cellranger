#!/usr/bin/env python
#
# Copyright (c) 2018 10X Genomics, Inc. All rights reserved.
#

import os
import cellranger.library_constants as lib_constants
from cellranger.feature_ref import FeatureDef, FeatureReference

def save_features_bed_path(feature_ref, fullpath):
    base_dir = os.path.dirname(fullpath)
    filename = os.path.basename(fullpath)
    save_features_bed(feature_ref, base_dir, filename)

def save_features_bed(feature_ref, base_dir, filename='peaks.bed'):
    out_features_fn = os.path.join(base_dir, filename)
    with open(out_features_fn, 'w') as f:
        for feature_def in feature_ref.feature_defs:
            f.write(feature_def.name.split(':')[0] + '\t' +
                    feature_def.name.split(':')[1].split('-')[0] + '\t' +
                    feature_def.name.split(':')[1].split('-')[1] + '\n')

def save_motifs_tsv(feature_ref, base_dir):
    out_features_fn = os.path.join(base_dir, 'motifs.tsv')
    with open(out_features_fn, 'w') as f:
        for feature_def in feature_ref.feature_defs:
            f.write(feature_def.id + '\t' + feature_def.name + '\n')

# A parallel implementation "get_genome_from_str" available via utils.py
def get_genome_from_contig(s, genomes):
    assert len(genomes) > 0

    if s is None:
        return None

    if len(genomes) == 1:
        return genomes[0]

    for genome in genomes:
        if s.startswith(genome):
            return genome

    raise Exception('%s does not have valid associated genome' % s)

def from_peaks_bed(peaks_path, genomes):
    '''
    Create a FeatureReference from a peaks bed file

    Args:
        peaks_path (str): Path to peaks bed file. Can be None.
        genomes: This should be the genomes on which peaks were defined. These would be prefixes in peak names
    '''

    # Load peaks info
    feature_defs = []
    all_tag_keys = ['genome', 'derivation']

    if peaks_path:
        # Stuff relevant fields of peak tuple into FeatureDef
        peaks = None
        with open(peaks_path, 'rU') as pf:
            peaks = ["{}:{}-{}".format(*line.strip("\n").split("\t")) for line in pf]
        for peak in peaks:
            genome = get_genome_from_contig(peak, genomes)
            feature_defs.append(FeatureDef(index=len(feature_defs),
                                id=peak,
                                name=peak,
                                feature_type=lib_constants.ATACSEQ_LIBRARY_TYPE,
                                tags={'genome': genome, 'derivation': ''}))

    return FeatureReference(feature_defs, all_tag_keys)

def get_genome_from_motif(s):
    '''extracts species name from a HOCOMOCO model name format:
        <motif>_<species>.<HOCOMOCO version>.<more versions>
        http://hocomoco11.autosome.ru/human/mono?full=false'''

    if s is None:
        return None

    return s.split('.')[0].split('_')[1]


def get_name_from_motif(s):
    '''extracts motif name from a HOCOMOCO model name format:
        <motif>_<species>.<HOCOMOCO version>.<more versions>
        http://hocomoco11.autosome.ru/human/mono?full=false
        returns <motif>
    '''

    if s is None:
        return None

    return s.split('.')[0].split('_')[0]

def from_motif_list(motif_list):
    '''
    Create a FeatureReference from a motifs list.

    Args:
        motif_list (list): list of motif names. Can be None.
    '''

    # Load peaks info
    feature_defs = []
    all_tag_keys = ['genome', 'derivation']

    if motif_list:
        # Stuff into FeatureDef
        for motif in motif_list:
            genome = get_genome_from_motif(motif)
            feature_defs.append(FeatureDef(index=len(feature_defs),
                                id=motif,
                                name=get_name_from_motif(motif),
                                feature_type=lib_constants.ATACSEQ_LIBRARY_DERIVED_TYPE,
                                tags={'genome': genome, 'derivation': 'POOL'}))

    return FeatureReference(feature_defs, all_tag_keys)

def from_bed_and_motifs(peaks_path, motif_list, genomes):
    '''
    Create a FeatureReference from a bed file of peaks and then from a motifs list.

    Args:
        peaks_path (str): bed file of peaks. Can be None.
        motif_list (list): list of motif names. Can be None.
        genome: This should be the genome for which the peaks and motifs are identified.
    '''

    # Load peaks info
    feature_defs = []
    all_tag_keys = ['genome', 'derivation']

    # process peaks
    if peaks_path:
        # Stuff relevant fields of peak tuple into FeatureDef
        peaks = None
        with open(peaks_path, 'rU') as pf:
            peaks = ["{}:{}-{}".format(*line.strip("\n").split("\t")) for line in pf]
        for peak in peaks:
            genome = get_genome_from_contig(peak, genomes)
            feature_defs.append(FeatureDef(index=len(feature_defs),
                                id=peak,
                                name=peak,
                                feature_type=lib_constants.ATACSEQ_LIBRARY_TYPE,
                                tags={'genome': genome, 'derivation': ''}))

    # process motifs
    if motif_list:
        # Stuff into FeatureDef
        for motif in motif_list:
            genome = get_genome_from_motif(motif)
            feature_defs.append(FeatureDef(index=len(feature_defs),
                                id=motif,
                                name=get_name_from_motif(motif),
                                feature_type=lib_constants.ATACSEQ_LIBRARY_DERIVED_TYPE,
                                tags={'genome': genome, 'derivation': 'POOL'}))

    return FeatureReference(feature_defs, all_tag_keys)
