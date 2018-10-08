import numpy as np
import pandas as pd
import cellranger.rna.library as rna_library

HIGH_UMI_CORRECTION_THRESHOLD = 0.5
NUM_READS_THRESHOLD = 10000

def detect_aggregate_bcs(correction_data):
    """ Given a csv containing the number of reads and umi corrected reads for each barcode,
        filter out barcodes that exceed pre-defined thresholds of fraction corrected reads and fraction total reads
    """

    summary = correction_data[correction_data['library_type'] == rna_library.ANTIBODY_LIBRARY_TYPE]
    corrected_ratio = summary['umi_corrected_reads'] / summary['reads']
    total_reads = np.float(np.sum(summary['reads']))
    reads_ratio = summary['reads'] / total_reads
    correction_metrics = pd.DataFrame({'Fraction Corrected Reads': corrected_ratio, 'Fraction Total Reads': reads_ratio})
    augmented = summary.join(correction_metrics)

    ### filter out the aggregate barcodes according to set thresholds
    bcs_to_remove = {}
    for idx, row in augmented.iterrows():
        if row['Fraction Corrected Reads'] > HIGH_UMI_CORRECTION_THRESHOLD and row['reads'] > NUM_READS_THRESHOLD:
            bcs_to_remove[row['barcode']] = row['Fraction Total Reads']

    ### make a slice of the barcode summary csv containing just the aggregate barcodes
    augmented.set_index("barcode", inplace=True)
    removed_bcs_df = augmented.loc[bcs_to_remove.keys()].drop(['umis', 'candidate_dup_reads'], axis=1)
    removed_bcs_df.round = removed_bcs_df.round({'Fraction Corrected Reads': 3, 'Fraction Total Reads': 3})

    return bcs_to_remove.keys(), np.sum(bcs_to_remove.values()), removed_bcs_df
