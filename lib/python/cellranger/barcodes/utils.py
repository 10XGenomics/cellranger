#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#
# Barcode handling functions

import os

def barcode_whitelist_path(fn):
    ''' Barcode whitelist is just a text file of valid barcodes, one per line.
        Lines containing the '#' character are ignored'''

    ## first check if fn is a valid path, if so load from there
    ## second check if it is in the current code directory
    ## third grab from tenkit or raise exception
    code_path = os.path.dirname(os.path.abspath(__file__))
    if os.path.exists(fn):
        return fn
    elif os.path.exists(os.path.join(code_path, fn+".txt")):
        return os.path.join(code_path, fn + ".txt")
    else:
        return fn

def load_barcode_whitelist(fn):
    ''' Barcode whitelist is just a text file of valid barcodes, one per line.
        Lines containing the '#' character are ignored'''

    barcodes = [ x.strip() for x in open(barcode_whitelist_path(fn), 'r') if not ('#' in x) ]
    return set(barcodes)
