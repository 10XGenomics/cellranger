#!/usr/bin/env python
#
# Copyright (c) 2018 10X Genomics, Inc. All rights reserved.
#

""" RNA-specific matrix functionality """

import cellranger.rna.feature_ref as rna_feature_ref

def save_mex(matrix, base_dir, sw_version, compress=True):
    """ Save an RNA matrix in Matrix Market Exchange format
    Args:
      matrix (CountMatrix): Matrix to write.
      base_dir (str): Path to output directory.
      sw_version (str): Version of this software.
    """
    mex_metadata = {
        'software_version': sw_version,
    }
    matrix.save_mex(base_dir,
                    rna_feature_ref.save_features_tsv,
                    metadata=mex_metadata,
                    compress=compress)
