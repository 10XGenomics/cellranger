#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

from collections import namedtuple, OrderedDict
import csv
import itertools
import re
import os
import cellranger.io as cr_io
import cellranger.library_constants as lib_constants
import cellranger.reference as cr_reference
import cellranger.utils as cr_utils
from cellranger.feature_ref import FeatureReference, FeatureDef, FeatureDefException

# These feature tag keys are reserved for internal use.
RESERVED_TAGS = ['genome']

# These are base FeatureDef attributes
BASE_FEATURE_FIELDS = ['id', 'name', 'feature_type']

# These feature tag keys are required for feature barcode extraction
REQUIRED_TAGS = ['read', 'pattern', 'sequence']

PatternEntry = namedtuple('PatternEntry', ['read', 'regex_string', 'regex', 'barcode_dict'])

''' The result of a feature matching operation.
    0  whitelist hits - (raw, [])   where 'raw' is the longest potential sequence.
    1  whitelist hit  - (raw, [id])
    >1 whitelist hits - (raw, [ids])
'''

FeatureMatchResult = namedtuple('FeatureMatchResult', ['barcode', 'qual', 'ids', 'indices'])

ALLOWED_READS = ['R1', 'R2']

class FeatureExtractor(object):
    '''Store a list of features (genes, antibodies, etc) and extract feature barcodes from a read sequence.'''

    def __init__(self, feature_ref, use_feature_types=None):
        '''Setup feature extraction.

        Args:
            feature_ref (FeatureReference): All feature definitions.
            use_feature_types (list of str): Feature types to extract.
        '''

        self.feature_ref = feature_ref

        self.patterns = FeatureExtractor._compile(feature_ref.feature_defs,
                                                  use_feature_types)


    @classmethod
    def _compile(cls, feature_defs, use_feature_types=None):
        '''Prepare for feature barcode extraction.

        Args:
            feature_defs (list of FeatureDef): feature defs to compile
            use_feature_types (list of str): feature types to try to extract.
                                              None to extract all.

        Returns:
            dict of (str,str) to PatternEntry
        '''
        # Store distinct patterns
        patterns = {}

        for fd in feature_defs:
            # Only build patterns from searchable feature types
            if use_feature_types is not None and fd.feature_type not in use_feature_types:
                continue

            sequence = fd.tags.get('sequence')
            pattern = fd.tags.get('pattern')
            read = fd.tags.get('read')

            if not sequence:
                continue
            if not pattern:
                raise FeatureDefException('Feature definition for %s has a sequence but no pattern specifying how to extract it.' % fd.id)
            if not read:
                raise FeatureDefException('Feature definition for %s has a sequence but no read specifying where to extract it from.' % fd.id)

            regex_str, regex = compile_pattern(pattern, len(sequence))

            # Treat equal regex pattern strings as generating equivalent automata
            pattern_key = (read, regex_str)
            if pattern_key not in patterns:
                patterns[pattern_key] = PatternEntry(read, regex_str, regex, {})
            entry = patterns[pattern_key]

            if sequence in entry.barcode_dict:
                raise FeatureDefException('Found two feature definitions with the same read ("%s"), pattern ("%s") and barcode sequence ("%s"). This combination of values must be unique for each row in the feature definition file.' % (read, pattern, sequence))

            entry.barcode_dict[sequence] = fd

        # Built from searchable features only
        return patterns


    def get_read_types(self):
        '''Determine which read types need to be inspected.
        Returns:
          (list of str): List of read types, e.g. ['R1'] '''
        return sorted(list(set(fd.tags.get('read') for fd in self.feature_ref.feature_defs if fd.tags.get('read'))))

    def extract_paired_end(self, r1_seq, r1_qual, r2_seq, r2_qual):
        '''Return a FeatureMatchResult '''
        r1_any_whitelist, r1_matches = self._extract_from_seq(r1_seq, r1_qual, 'R1')
        r2_any_whitelist, r2_matches = self._extract_from_seq(r2_seq, r2_qual, 'R2')
        return self._filter_feature_matches(r1_matches + r2_matches,
                                            r1_any_whitelist or r2_any_whitelist)

    def extract_single_end(self, seq, qual, read_type):
        '''Return a FeatureMatchResult '''
        any_whitelist, matches = self._extract_from_seq(seq, qual, read_type)
        return self._filter_feature_matches(matches, any_whitelist)

    def _extract_from_seq(self, seq, qual, read_type):
        '''Extract the feature barcode from a sequence.
        Returns:
            list of (str, FeatureDef): Tuples are (barcode, feature_def)
                where barcode is the observed sequence
                and feature_def is the matching feature definition or None if there isn't one.
        '''
        matches = []
        any_whitelist_hits = False

        for ((entry_read, _), entry) in self.patterns.iteritems():
            if entry_read != read_type:
                continue

            # Extract potential barcode sequences
            regex_match = re.search(entry.regex, seq)
            if regex_match is not None:
                span = regex_match.span(1)
                bc = seq[span[0]:span[1]]
                bc_qual = qual[span[0]:span[1]]

                # Look for exact whitelist match
                hit = entry.barcode_dict.get(bc)
                if hit is not None:
                    any_whitelist_hits = True
                matches.append((bc, bc_qual, hit))

        return any_whitelist_hits, matches

    def _filter_feature_matches(self, matches, any_whitelist_hits):
        '''Filter a list of feature matches by prioritization.

        Args:
            matches (list of (str, str, FeatureDef)): bc, qual, hit.
            any_whitelist_hits (bool): True if any of the hits matched the whitelist.

        Returns:
           FeatureMatchResult
        '''

        if any_whitelist_hits:
            # >0 whitelist hits. Remove entries w/ no whitelist hit.
            matches = filter(lambda bc_hit: bc_hit[2], matches)

            # Take the longest observed sequence as the canonical sequence
            # NOTE: If there are multiple whitelist hits of equal length,
            #   we choose an arbitrary sequence as the canonical sequence.
            # However, information is not lost because the feature IDs tag
            #   contains a list of each hit's unique feature ID.
            best_hit = max(matches, key=lambda bc_hit: len(bc_hit[0]))
            barcode, qual = best_hit[0:2]

            # Record ids and indices of all whitelist hits (typically 1)
            ids = map(lambda bc_hit: bc_hit[2].id, matches)
            indices = map(lambda bc_hit: bc_hit[2].index, matches)

            return FeatureMatchResult(barcode, qual, ids, indices)

        elif len(matches) > 0:
            # No whitelist hits but found some sequence.
            # Take the longest observed sequence.
            # NOTE: This will do the wrong thing if:
            #     1) There are multiple possible barcode lengths for a given pattern,
            #     2) This sequence doesn't match the whitelist at any barcode length,
            #     3) And any but the longest observed barcode is 1-away from the whitelist.
            #   The correct thing to do is to prioritize sequences that are near-matches to the whitelist
            #   instead of just taking the longest sequence.
            best_hit = max(matches, key=lambda bc_hit: len(bc_hit[0]))
            barcode, qual = best_hit[0:2]

            return FeatureMatchResult(barcode, qual, [], [])

        else:
            return FeatureMatchResult(None, None, [], [])

    def has_features_to_extract(self):
        '''Does this class have any interesting work to do.

        Returns:
           (bool): True if this feature ref was compiled to extract any feature BCs.
        '''
        return len(self.patterns) > 0

def save_features_tsv(feature_ref, base_dir, compress):
    """Save a FeatureReference to a tsv file"""
    out_features_fn = os.path.join(base_dir, 'features.tsv')
    if compress:
        out_features_fn += '.gz'

    with cr_io.open_maybe_gzip(out_features_fn, 'w') as f:
        for feature_def in feature_ref.feature_defs:
            f.write('\t'.join((feature_def.id,
                               feature_def.name,
                               feature_def.feature_type)) + '\n')

def from_transcriptome_and_csv(gene_ref_path, feature_def_filename):
    '''Create a FeatureReference.

    Create a FeatureReference from a transcriptome ref and a feature barcode ref.

    Args:
        gene_ref_path (str): Path to transcriptome reference. Can be None.
        feature_def_filename (str): Path to Feature Definition CSV file. Can be None.
    Returns:
        FeatureReference
    '''

    # Load gene info
    feature_defs = []
    all_tag_keys = ['genome']

    genomes = cr_utils.get_reference_genomes(gene_ref_path)

    if gene_ref_path is not None:
        gene_idx_filename = cr_utils.get_reference_genes_index(gene_ref_path)
        gene_index = cr_reference.GeneIndex.load_pickle(gene_idx_filename)

        # Stuff relevant fields of Gene tuple into FeatureDef
        for gene in gene_index.genes:
            genome = cr_utils.get_genome_from_str(gene.id, genomes)
            fd = FeatureDef(index=len(feature_defs),
                            id=gene.id,
                            name=gene.name,
                            feature_type=lib_constants.GENE_EXPRESSION_LIBRARY_TYPE,
                            tags={'genome': genome,})
            feature_defs.append(fd)

    # Load feature definition file
    if feature_def_filename is not None:
        csv_feature_defs, csv_tag_keys = parse_feature_def_file(feature_def_filename,
                                                                index_offset=len(feature_defs))

        feature_defs.extend(csv_feature_defs)
        all_tag_keys.extend(csv_tag_keys)

    return FeatureReference(feature_defs, all_tag_keys)


def validate_sequence(seq):
    '''Validate a feature barcode sequence.'''
    if len(seq) == 0:
        raise FeatureDefException('Feature sequence must be non-empty.')
    if not re.match('^[ACGTN]+', seq):
        raise FeatureDefException('Invalid sequence: "%s". The only allowed characters are A, C, G, T, and N.' % seq)

def compile_pattern(pattern_str, length):
    '''Compile a feature definition pattern into a regex.'''
    if '(BC)' not in pattern_str:
        raise FeatureDefException('Invalid pattern: "%s". The pattern must contain the string "(BC)".' % pattern_str)

    check_pattern = re.sub('\(BC\)', '', pattern_str)
    if not re.match('\^{0,1}[ACGTN]*\${0,1}', check_pattern):
        raise FeatureDefException('Invalid pattern: "%s". The pattern must optionally start with "^", optionally end with "$", contain exactly one instance of the string "(BC)" and otherwise contain only the characters A, C, G, T, and N.')

    # Allow Ns to match anything
    regex_str = re.sub('N', '.', pattern_str)

    # Capture the feature barcode
    regex_str = re.sub('\(BC\)', '(.{%d,%d})' % (length, length), regex_str)

    try:
        regex = re.compile(regex_str)
    except re.error:
        raise FeatureDefException('Failed to compile the feature definition string "%s" into a regular expression: %s.' % (regex_str, str(re.error)))

    return regex_str, regex

def get_required_csv_columns():
    '''Get a list of required CSV columns. '''
    return [col for col in (BASE_FEATURE_FIELDS + REQUIRED_TAGS)]

def parse_feature_def_file(filename, index_offset=0):
    '''Load a CSV feature definition file.
    Args:
        filename (str): CSV filename.
        index_offset (int): Start the FeatureDef indices at this value.
    Returns:
        feature_defs (list of FeatureDef): All feature definitions
        all_tag_keys (list of str): All tag column names in original order
    '''

    feature_defs = []
    seen_ids = set()
    all_tag_keys = []

    with open(filename, 'rU') as f:
        # Skip comments
        rows = itertools.ifilter(lambda x: not x.startswith('#'), f)

        reader = csv.DictReader(rows)

        required_cols = get_required_csv_columns()

        feature_index = index_offset

        all_tag_keys = [c for c in reader.fieldnames if c not in BASE_FEATURE_FIELDS]

        for row in reader:
            # Check field presence
            if not set(row.keys()).issuperset(set(required_cols)):
                raise FeatureDefException('The feature definition file header must contain the following comma-delimited fields: "%s".' % ', '.join(required_cols))

            # Strip flanking whitespace from values
            for key in row.iterkeys():
                row[key] = row[key].strip()

            # Check ID uniqueness
            if row['id'] in seen_ids:
                raise FeatureDefException('Found duplicated ID in feature definition file: "%s"' % row['id'])
            seen_ids.add(row['id'])

            if "\t" in row['name']:
                raise FeatureDefException("Feature name field cannot contain tabs: '%s'" % row['name'])

            if ' ' in row['id'] or "\t" in row['id']:
                raise FeatureDefException("Feature id field cannot contain whitespace: '%s'" % row['id'])

            # Additional columns become key-value pairs
            # Maintain input order
            tag_cols = [c for c in reader.fieldnames if c not in BASE_FEATURE_FIELDS and c in row]

            tags = OrderedDict()
            for key in tag_cols:
                if key in RESERVED_TAGS:
                    raise FeatureDefException('Found invalid column name "%s." This name cannot be used as a custom feature tag because it is reserved for internal use.' % key)

                # Always store as strings
                tags[key] = str(row[key])

            # Validate fields
            if len(tags['sequence']) == 0:
                raise FeatureDefException('Found blank feature barcode sequence for feature id %s. The sequence column must be populated for this feature.' % row['id'])
            validate_sequence(tags['sequence'])

            if len(tags['pattern']) == 0:
                raise FeatureDefException('Found blank feature barcode pattern for feature id %s. The pattern column must be populated for this feature.' % row['id'])
            compile_pattern(tags['pattern'], len(tags['sequence']))

            if tags['read'] not in ALLOWED_READS:
                raise FeatureDefException('The feature definition file contains a read type value "%s" which is not one of the allowed read types %s.' % (tags['read'], str(ALLOWED_READS)))

            feature_defs.append(FeatureDef(index=feature_index,
                                           id=row['id'],
                                           name=row['name'],
                                           feature_type=row['feature_type'],
                                           tags=tags))

            feature_index += 1

    return feature_defs, all_tag_keys
