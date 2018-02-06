#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

from collections import namedtuple
import csv
import itertools
import re

FeatureDef = namedtuple('FeatureDef', ['index', 'id', 'name', 'read', 'pattern', 'sequence', 'feature_type', 'tags'])

PatternEntry = namedtuple('PatternEntry', ['read', 'regex_string', 'regex', 'barcode_dict'])

''' The result of a feature matching operation.
    0  whitelist hits - (raw, [])   where 'raw' is the longest potential sequence.
    1  whitelist hit  - (raw, [id])
    >1 whitelist hits - (raw, [ids])
'''

FeatureMatchResult = namedtuple('FeatureMatchResult', ['barcode', 'qual', 'ids', 'indices'])

NOT_COLUMNS = ['index']
OPTIONAL_COLUMNS = ['tags']

ALLOWED_READS = ['R1', 'R2']

class FeatureDefException(Exception):
    pass

class FeatureReference(object):
    ''' Store a list of features and extract them from a read sequence '''

    def __init__(self, filename, use_feature_types=None):
        ''' Load a feature definition file '''
        feature_defs = parse_feature_def_file(filename)

        # Store distinct patterns
        patterns = {}

        # Record total number of features (pre-filtering)
        self.num_total_features = len(feature_defs)

        # Filter features by type, if specified
        if use_feature_types is not None:
            feature_defs = filter(lambda x: x.feature_type in use_feature_types, feature_defs)

        for fd in feature_defs:
            regex_str, regex = compile_pattern(fd.pattern, len(fd.sequence))

            # Treat equal regex pattern strings as generating equivalent automata
            pattern_key = (fd.read, regex_str)
            if pattern_key not in patterns:
                patterns[pattern_key] = PatternEntry(fd.read, regex_str, regex, {})
            entry = patterns[pattern_key]

            if fd.sequence in entry.barcode_dict:
                raise FeatureDefException('Found two feature definitions with the same read ("%s"), pattern ("%s") and barcode sequence ("%s"). This combination of values must be unique for each row in the feature definition file.' % (fd.read, fd.pattern, fd.sequence))

            entry.barcode_dict[fd.sequence] = fd

        self.feature_defs = feature_defs
        self.patterns = patterns

    def get_read_types(self):
        ''' Determine which read types need to be inspected '''
        return sorted(list(set([fd.read for fd in self.feature_defs])))

    def extract_paired_end(self, r1_seq, r1_qual, r2_seq, r2_qual):
        ''' Return a FeatureMatchResult '''
        r1_any_whitelist, r1_matches = self._extract_from_seq(r1_seq, r1_qual, 'R1')
        r2_any_whitelist, r2_matches = self._extract_from_seq(r2_seq, r2_qual, 'R2')
        return self._filter_feature_matches(r1_matches + r2_matches,
                                            r1_any_whitelist or r2_any_whitelist)

    def extract_single_end(self, seq, qual, read_type):
        ''' Return a FeatureMatchResult '''
        any_whitelist, matches = self._extract_from_seq(seq, qual, read_type)
        return self._filter_feature_matches(matches, any_whitelist)

    def _extract_from_seq(self, seq, qual, read_type):
        ''' Extract the feature barcode from a sequence.
            Returns a list of (barcode, feature_def) tuples where barcode is the observed sequence
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
        if any_whitelist_hits:
            # >0 whitelist hits. Remove entries w/ no whitelist hit.
            matches = filter(lambda bc_hit: bc_hit[2], matches)

            # Take the longest observed sequence as the canonical sequence
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


def validate_sequence(seq):
    if len(seq) == 0:
        raise FeatureDefException('Feature sequence must be non-empty.')
    if not re.match('^[ACGT]+', seq):
        raise FeatureDefException('Invalid sequence: "%s". The only allowed characters are A, C, G, and T.' % seq)

def compile_pattern(pattern_str, length):
    ''' Compile a feature definition pattern into a regex '''
    if '(BC)' not in pattern_str:
        raise FeatureDefException('Invalid pattern: "%s". The pattern must contain the string "(BC)".' % pattern_str)

    check_pattern = re.sub('\(BC\)', '', pattern_str)
    if not re.match('\^{0,1}[ACGT]*\${0,1}', check_pattern):
        raise FeatureDefException('Invalid pattern: "%s". The pattern must optionally start with "^", optionally end with "$", contain exactly one instance of the string "(BC)" and otherwise contain only the characters A, C, G, and T')

    regex_str = re.sub('\(BC\)', '(.{%d,%d})' % (length, length), pattern_str)
    try:
        regex = re.compile(regex_str)
    except re.error:
        raise FeatureDefException('Failed to compile the feature definition string "%s" into a regular expression: %s.' % (regex_str, str(re.error)))

    return regex_str, regex

def parse_feature_def_file(filename):
    ''' Load a feature definition file '''

    feature_defs = []
    seen_ids = set()

    with open(filename, 'rU') as f:
        # Skip comments
        rows = itertools.ifilter(lambda x: not x.startswith('#'), f)

        reader = csv.DictReader(rows)

        required_cols = [col for col in FeatureDef._fields if col not in (OPTIONAL_COLUMNS + NOT_COLUMNS)]

        feature_index = 0

        for row in reader:
            # Check field presence
            if set(row.keys()) != set(required_cols):
                raise FeatureDefException('The feature definition file header must contain the following comma-delimited fields: "%s".' % ', '.join(required_cols))

            # Strip flanking whitespace from values
            for key in row.iterkeys():
                row[key] = row[key].strip()

            # Check ID uniqueness
            if row['id'] in seen_ids:
                raise FeatureDefException('Found duplicated ID in feature definition file: "%s"' % row['id'])
            seen_ids.add(row['id'])

            # Additional columns become key-value pairs
            tags = {}
            for key in set(row.keys()) - set(required_cols):
                tags[key] = row[key]
                del row[key]
            row['tags'] = tags

            # Validate fields
            validate_sequence(row['sequence'])

            compile_pattern(row['pattern'], len(row['sequence']))

            if row['read'] not in ALLOWED_READS:
                raise FeatureDefException('The feature definition file contains a read type value "%s" which is not one of the allowed read types %s.' % str(ALLOWED_READS))

            row['index'] = feature_index

            feature_defs.append(FeatureDef(**row))

            feature_index += 1

    return feature_defs
