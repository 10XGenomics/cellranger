#!/usr/bin/env python
#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#
# Takes in 2 dictionaries and forms a third dictionary,
# with values that are the sum of the two (defaulting to 0 for any absent keys)
#

def add_dicts(in_dict1, in_dict2, depth):
    out_dict = {}

    if in_dict1 is None and in_dict2 is None:
        return None

    # Keys present in first dictionary, and +/- present in second dictionary
    for (key1, value1) in in_dict1.iteritems():
        if depth == 1:
            value2 = in_dict2.get(key1, 0)
            out_dict[key1] = value1 + value2
        else:
            value2 = in_dict2.get(key1, {})
            out_dict[key1] = add_dicts(value1, value2, depth - 1)

    # Keys only present in dictionary 2
    for (key2, value2) in in_dict2.iteritems():
        if not(in_dict1.has_key(key2)):
            out_dict[key2] = value2

    return out_dict

def combine_dicts(dicts, depth):
    return reduce(lambda x,y: add_dicts(x,y,depth), dicts)

def get_key_with_max_value(dictionary):                                                                  
    if len(dictionary) == 0:
        return None
    maxbase = max(dictionary, key=dictionary.get)
    base_count = dictionary[maxbase]
    if base_count == 0:
        return None
    return max(dictionary, key=dictionary.get)
