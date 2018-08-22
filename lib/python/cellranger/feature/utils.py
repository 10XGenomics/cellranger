#!/usr/bin/env python
#
# Copyright (c) 2018 10X Genomics, Inc. All rights reserved.
#
# Utils for feature-barcoding technology

import numpy as np
import os
import json
import cellranger.feature.constants as feature_constants
import tenkit.safe_json as tk_safe_json

FEATURE_JSON_KEY_MAP = {
'total_reads': 'total_reads',
'reads_per_cell': 'multi_transcriptome_total_raw_reads_per_filtered_bc',
'feature_reads_usable': 'multi_usable_reads',
'frac_valid_barcodes': 'good_bc_frac',
'frac_corrected_barcodes': 'corrected_bc_frac',
'frac_valid_umis': 'good_umi_frac',
'frac_reads_umi_corrected': 'corrected_umi_frac',
'frac_raw_feature': 'feature_bc_extracted_frac',
'frac_feature_reads': 'multi_transcriptome_conf_mapped_reads_frac',
'frac_feature_reads_usable': 'multi_transcriptome_usable_reads_frac',
'feature_reads_usable_per_cell': 'multi_usable_reads_per_filtered_bc',
'frac_reads_feature_unknown': 'unrecognized_feature_bc_frac',
'frac_reads_feature_corrected': 'corrected_feature_bc_frac',
'frac_reads_chimeras': 'low_support_umi_reads_frac',
'feature_reads_in_cells': 'multi_filtered_bcs_conf_mapped_barcoded_reads_cum_frac',
'median_umis_per_cell': 'multi_filtered_bcs_median_counts',
}

def translate_feature_metrics_json_keys(json_file_path, f_type):
    if (json_file_path is None) or not(os.path.isfile(json_file_path)):
        return {}
    return_dict = {}
    input_metrics = json.load(open(json_file_path, 'r'))

    report_prefix = feature_constants.PREFIX_FROM_FEATURE_TYPE.get(f_type, 'FEATURE') + '_'
    for key in FEATURE_JSON_KEY_MAP:
        return_dict[report_prefix + key] = input_metrics.get(report_prefix + FEATURE_JSON_KEY_MAP[key])

    return return_dict

def write_json_from_dict(input_dict, out_file_name):
    with open(out_file_name, 'w') as f:
        json.dump(tk_safe_json.json_sanitize(input_dict), f, indent=4, sort_keys=True)

def write_csv_from_dict(input_dict, out_file_name, header=None):
    with open(out_file_name, 'w') as f:
        if header is not None:
            f.write(header)
        for (key, value) in input_dict.iteritems():
            line = str(key) + ',' + str(value)  + '\n'
            f.write(line)

def get_depth_string(num_reads_per_cell):
    return str(np.round(float(num_reads_per_cell)/1000,1)) + "k"

