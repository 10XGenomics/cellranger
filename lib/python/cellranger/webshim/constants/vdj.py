#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

from __future__ import annotations

import cellranger.library_constants as lib_constants
import cellranger.vdj.utils as vdj_utils
import cellranger.webshim.constants.shared as shared


def vdj_gene_pair_format(gene_pair_str):
    """Format a gene pair prefix from VdjReporter for human-usage."""
    if gene_pair_str == lib_constants.MULTI_REFS_PREFIX:
        return ""
    else:
        return "({})".format(", ".join(vdj_utils.get_genes_in_pair(gene_pair_str)))


# NOTE: Q30 metrics are currently shared with count.
#       However, RNA Read is renamed to RNA Read 1 here.
TOTAL_READ_PAIRS_METRIC = {
    "name": "VDJ_total_read_pairs",
    "display_name": "Number of Read Pairs",
    "description": "Total number of read pairs that were assigned to this library in demultiplexing.",
    "format": "integer",
}

Q30_METRICS = [
    {
        "name": "VDJ_bc_bases_with_q30_frac",
        "display_name": "Q30 Bases in Barcode",
        "description": "Fraction of cell barcode bases with Q-score >= 30, excluding very low quality/no-call (Q <= 2) bases from the denominator.",
        "format": "percent",
    },
    {
        "name": "VDJ_read_bases_with_q30_frac",
        "display_name": "Q30 Bases in RNA Read 1",
        "description": "Fraction of RNA read 1 bases with Q-score >= 30, excluding very low quality/no-call (Q <= 2) bases from the denominator.",
        "format": "percent",
    },
    {
        "name": "VDJ_read2_bases_with_q30_frac",
        "display_name": "Q30 Bases in RNA Read 2",
        "description": "Fraction of RNA read 2 bases with Q-score >= 30, excluding very low quality/no-call (Q <= 2) bases from the denominator.",
        "format": "percent",
    },
    {
        "name": "VDJ_umi_bases_with_q30_frac",
        "display_name": "Q30 Bases in UMI",
        "description": "Fraction of UMI bases with Q-score >= 30, excluding very low quality/no-call (Q <= 2) bases from the denominator.",
        "format": "percent",
    },
]

GOOD_BCS_METRIC = {
    "name": "vdj_good_bc_frac",
    "display_name": "Valid Barcodes",
    "description": "Fraction of reads with barcodes that match the whitelist after barcode correction.",
    "format": "percent",
}

NUMBER_OF_CELLS_METRIC = {
    "name": "vdj_filtered_bcs",
    "display_name": "Estimated Number of Cells",
    "description": "The number of barcodes estimated to be associated with cells that express targeted V(D)J transcripts. A barcode is cell-associated if has 2 UMIs with sufficient read support.",
    "format": "integer",
}

MEAN_READ_PAIRS_PER_CELL_METRIC = {
    "name": "vdj_total_raw_read_pairs_per_filtered_bc",
    "display_name": "Mean Read Pairs per Cell",
    "description": "Number of input read pairs divided by the estimated number of cells.",
    "format": "integer",
}

MEAN_ASSEMBLED_READ_PAIRS_PER_CELL_METRIC = {
    "name": "vdj_assemblable_read_pairs_per_filtered_bc",
    "display_name": "Mean Used Read Pairs per Cell",
    "description": "Mean number of read pairs used in assembly per cell-associated barcode. These reads must have a cell-associated barcode, map to a V(D)J gene, and have a UMI with sufficient read support.",
    "format": "integer",
}

FRAC_READS_IN_CELLS_METRIC = {
    "name": "vdj_filtered_bcs_cum_frac",
    "display_name": "Fraction Reads in Cells",
    "description": "Number of reads with cell-associated barcodes divided by the number of reads with valid barcodes.",
    "format": "percent",
}

MULTI_READS_MAPPED_TO_RECOMBINOME_METRIC = {
    "name": "multi_vdj_recombinome_mapped_reads_frac",
    "display_name": "Reads Mapped to Any V(D)J Gene",
    "description": "Fraction of reads with valid barcodes that partially or wholly map to any germline V(D)J gene segment.",
    "format": "percent",
}
READS_MAPPED_TO_RECOMBINOME_METRIC = {
    "name": "vdj_recombinome_mapped_reads_frac",
    "display_name": "Reads Mapped to %s",
    "description": "Fraction of reads with valid barcodes that map partially or wholly to a germline %s gene segment.",
    "format": "percent",
    "prefix": "canonical_vdj_genes_nomulti",
}

MEDIAN_UMIS_PER_CELL_METRIC = {
    "name": "vdj_assembly_umis_per_cell_median",
    "display_name": "Median %s UMIs per Cell",
    "description": "Median number of UMIs assigned to a %s contig per cell. For B cells, only the max of {IGK, IGL} are counted.",
    "prefix": "canonical_vdj_genes_nomulti",
    "format": "%0.1f",
}
TRA_MEDIAN_UMIS_PER_CELL_METRIC = dict(MEDIAN_UMIS_PER_CELL_METRIC, prefix_filter="TRA")
TRB_MEDIAN_UMIS_PER_CELL_METRIC = dict(MEDIAN_UMIS_PER_CELL_METRIC, prefix_filter="TRB")

MULTI_NUM_CELLS_PAIRED_VJ_SPANNING_PROD_METRIC = {
    "name": "multi_vdj_assembly_contig_pair_productive_full_len_bc_count",
    "display_name": "Number of Cells With Productive V-J Spanning Pair",
    "description": "Number of cell barcodes for which at least 1 sequence was found for each gene in a receptor pair. Both sequences must have an annotated CDR3, must be predicted to be productive, and must span the beginning of the 5' end of the V region to at least 2bp from the 3' end of the J region.",
    "format": "integer",
}
NUM_CELLS_PAIRED_VJ_SPANNING_PROD_METRIC = {
    "name": "vdj_assembly_contig_pair_productive_full_len_bc_count",
    "display_name": "Number of Cells With Productive V-J Spanning %s Pair",
    "description": "Number of cell barcodes for which at least 1 sequence was found for each gene in a receptor pair. Both sequences must have an annotated CDR3, must be predicted to be productive, and must span the beginning of the 5' end of the V region to at least 2bp from the 3' end of the J region.",
    "prefix": "canonical_vdj_gene_pairs",
    "prefix_format": vdj_gene_pair_format,
    "format": "integer",
}
# # Make this name/prefix combo static internal target purposes
# TR_NUM_CELLS_PAIRED_VJ_SPANNING_PROD_METRIC = {
#     'name': vdj_constants.CANONICAL_TR_GENE_PAIRS[0] + '_vdj_assembly_contig_pair_productive_full_len_bc_count',
#     'display_name': 'Number of Cells With Productive V-J Spanning %s Pair' % vdj_gene_pair_format(vdj_constants.CANONICAL_TR_GENE_PAIRS[0]),
#     'description': 'Number of cell barcodes for which at least 1 sequence was found for each of TRA and TRB. Both sequences must have an annotated CDR3, must be predicted to be productive, and must span the beginning of the 5\' end of the V region to at least 2bp from the the 3\' end of the J region.',
#     'format': 'integer',
# }

PAIRED_CLONOTYPE_DIVERSITY_METRIC = {
    "name": "multi_raw_vdj_paired_clonotype_diversity",
    "display_name": "Paired Clonotype Diversity",
    "description": "Effective diversity of the paired clonotypes, computed as the Inverse Simpson Index of the clonotype frequencies. A value of 1 indicates a minimally diverse sample - only one distinct clonotype was detected. A value equal to the estimated number of cells indicates a maximally diverse sample.",
    "format": "%0.2f",
}

MULTI_FRAC_CELLS_PAIRED_VJ_SPANNING_PROD_METRIC = {
    "name": "multi_vdj_assembly_contig_pair_productive_full_len_bc_frac",
    "display_name": "Cells With Productive V-J Spanning Pair",
    "description": "Fraction of cell-associated barcodes with at least one contig for each chain of the receptor pair satisfying the following: "
    + "the contig annotations span the 5' end of the V region to at least 2 bp from the 3' end of the J region of the chain, a start codon was found in the expected part of the V sequence, an in-frame CDR3 amino acid motif was found, and no stop codons were found in the aligned V-J region.",
    "format": "percent",
}
FRAC_CELLS_PAIRED_VJ_SPANNING_PROD_METRIC = {
    "name": "vdj_assembly_contig_pair_productive_full_len_bc_frac",
    "display_name": "Cells With Productive V-J Spanning %s Pair",
    "description": "Fraction of cell-associated barcodes with at least one contig for each chain of the receptor pair satisfying the following: "
    + "the contig annotations span the 5' end of the V region to at least 2 bp from the 3' end of the J region of the chain, a start codon was found in the expected part of the V sequence, an in-frame CDR3 amino acid motif was found, and no stop codons were found in the aligned V-J region.",
    "format": "percent",
    "prefix": "canonical_vdj_gene_pairs_nomulti",
    "prefix_format": vdj_gene_pair_format,
}

FRAC_CELLS_VDJ_GENE_DETECTED_METRIC = {
    "name": "vdj_assembly_contig_bc_frac",
    "display_name": "Cells With %s Contig",
    "description": "Fraction of cell-associated barcodes with at least one %s contig annotated as a full or partial V(D)J gene.",
    "format": "percent",
    "prefix": "canonical_vdj_genes_nomulti",
}

FRAC_CELLS_CDR3_DETECTED_METRIC = {
    "name": "vdj_assembly_cdr_detected_bc_frac",
    "display_name": "Cells With CDR3-annotated %s Contig",
    "description": "Fraction of cell-associated barcodes with at least one %s contig where a CDR3 was detected.",
    "format": "percent",
    "prefix": "canonical_vdj_genes_nomulti",
}

FRAC_CELLS_VJ_SPANNING_METRIC = {
    "name": "vdj_assembly_contig_full_len_bc_frac",
    "display_name": "Cells With V-J Spanning %s Contig",
    "description": "Fraction of cell-associated barcodes with at least one contig spanning the 5' end of the V region to at least 2bp from the 3' end of the J region for %s.",
    "format": "percent",
    "prefix": "canonical_vdj_genes_nomulti",
}

FRAC_CELLS_PRODUCTIVE_METRIC = {
    "name": "vdj_assembly_prod_cdr_bc_frac",
    "display_name": "Cells With Productive %s Contig",
    "description": "Fraction of cell-associated barcodes with at least one contig that spans the 5' end of the V region to at least 2bp from the 3' end of the J region for %s, has a start codon in the expected part of the V sequence, has an in-frame CDR3, and has no stop codons in the aligned V-J region.",
    "format": "percent",
    "prefix": "canonical_vdj_genes_nomulti",
}

SUMMARY_METRICS = [
    NUMBER_OF_CELLS_METRIC,
    MEAN_READ_PAIRS_PER_CELL_METRIC,
    MULTI_NUM_CELLS_PAIRED_VJ_SPANNING_PROD_METRIC,
]

SEQUENCING_METRICS = [
    TOTAL_READ_PAIRS_METRIC,
    GOOD_BCS_METRIC,
] + Q30_METRICS

ENRICHMENT_METRICS = [
    MULTI_READS_MAPPED_TO_RECOMBINOME_METRIC,
    READS_MAPPED_TO_RECOMBINOME_METRIC,
]

CELL_METRICS = [
    NUMBER_OF_CELLS_METRIC,
    MEAN_READ_PAIRS_PER_CELL_METRIC,
    MEAN_ASSEMBLED_READ_PAIRS_PER_CELL_METRIC,
    FRAC_READS_IN_CELLS_METRIC,
]

VDJ_EXPRESSION_METRICS = [
    MEDIAN_UMIS_PER_CELL_METRIC,
]

VDJ_ANNOTATION_METRICS = [
    MULTI_FRAC_CELLS_PAIRED_VJ_SPANNING_PROD_METRIC,
    FRAC_CELLS_PAIRED_VJ_SPANNING_PROD_METRIC,
    PAIRED_CLONOTYPE_DIVERSITY_METRIC,
    FRAC_CELLS_VDJ_GENE_DETECTED_METRIC,
    FRAC_CELLS_CDR3_DETECTED_METRIC,
    FRAC_CELLS_VJ_SPANNING_METRIC,
    FRAC_CELLS_PRODUCTIVE_METRIC,
]

METRICS = [
    {
        "name": "Summary",
        "metrics": SUMMARY_METRICS,
    },
    {
        "name": "Sequencing",
        "metrics": SEQUENCING_METRICS,
    },
    {
        "name": "Enrichment",
        "metrics": ENRICHMENT_METRICS,
    },
    {
        "name": "Cells",
        "metrics": CELL_METRICS,
    },
    {
        "name": "V(D)J Expression",
        "metrics": VDJ_EXPRESSION_METRICS,
    },
    {
        "name": "V(D)J Annotation",
        "metrics": VDJ_ANNOTATION_METRICS,
    },
]

CHARTS = [
    {
        "layout": {
            "title": "Barcode Rank",
            "width": 470,
            "height": 313,
            "margin": {"l": 60, "r": 0, "t": 30, "b": 40},
            "hovermode": "closest",
            "xaxis": {
                "title": "Barcodes",
                "type": "log",
            },
            "yaxis": {
                "title": "UMI counts",
                "type": "log",
            },
        },
        "data": [],
        "config": shared.CHARTS_PLOTLY_FIXED_CONFIG,
        "function": "plot_vdj_barcode_rank",
        "name": "vdj_barcode_rank",
    },
    {
        "layout": {
            "title": "Top 10 Clonotype Frequencies",
            "width": shared.DEFAULT_CHART_WIDTH,
            "height": shared.DEFAULT_CHART_HEIGHT,
            "showlegend": False,
            "xaxis": {
                "title": "Clonotype",
            },
            "yaxis": {
                "title": "Fraction of Cells",
            },
        },
        "data": [{"type": "bar"}],
        "kwargs": {
            "metric_name": "raw_vdj_clonotype_prop",
            "order_by": shared.HISTOGRAM_METRIC_ORDER_DECREASING_PROPORTION,
        },
        "config": shared.CHARTS_PLOTLY_FIXED_CONFIG,
        "function": "plot_histogram_metric",
        "name": "vdj_clonotype_histogram",
        "description": 'Top 10 clonotypes by frequency in this sample. A clonotype is defined as a unique set of CDR3 nucleotide sequences. For the full table and more details, please refer to the "clonotypes.csv" and "consensus_annotations.csv" files produced by the pipeline.',
    },
    {
        "table": {
            "width": 999,
            "height": 1000,
        },
        "title": "Top 10 Clonotype CDR3 Sequences",
        "name": "clonotype_table",
        "description": 'Top 10 clonotypes by frequency in this sample. A clonotype is defined as a unique set of CDR3 nucleotide sequences. For the full table and more details, please refer to the "clonotypes.csv" and "consensus_annotations.csv" files produced by the pipeline.',
        "function": "plot_clonotype_table",
    },
]

METRIC_ALARMS = [
    {
        "name": NUMBER_OF_CELLS_METRIC["name"],
        "format": NUMBER_OF_CELLS_METRIC["format"],
        "error": {
            "title": "No Cells Detected",
            "message": "No valid sequencing data was detected. Please check the sequencing data.",
            "test": "< 1",
        },
        "warn": {
            "title": "Low Number of Cells Detected",
            "message": "Ideal >= 10. This usually indicates poor cell quality, poor library quality, or poor sequencing quality. Application performance is likely to be affected.",
            "test": "< 10",
        },
    },
    {
        "name": MULTI_FRAC_CELLS_PAIRED_VJ_SPANNING_PROD_METRIC["name"],
        "format": MULTI_FRAC_CELLS_PAIRED_VJ_SPANNING_PROD_METRIC["format"],
        "error": {
            "title": "Low Cells With Productive V-J Spanning Pair.",
            "message": "Ideal > 30%. This can indicate poor cell quality, low yield from the RT reaction, poor specificity of the V(D)J enrichment, poor sequencing quality, or the use of an unsupported chemistry type (e.g., using Single Cell 3' for V(D)J assembly). Application performance is likely to be affected",
            "test": "< 0.2",
        },
        "warn": {
            "title": "Low Cells With Productive V-J Spanning Pair.",
            "message": "Ideal > 30%. This can indicate poor cell quality, low yield from the RT reaction, poor specificity of the V(D)J enrichment, poor sequencing quality, or the use of an unsupported chemistry type (e.g., using Single Cell 3' for V(D)J assembly). Application performance may be affected",
            "test": "< 0.3",
        },
    },
    {
        "name": MULTI_READS_MAPPED_TO_RECOMBINOME_METRIC["name"],
        "format": MULTI_READS_MAPPED_TO_RECOMBINOME_METRIC["format"],
        "error": {
            "title": "Low Fraction of Reads Mapped to Any V(D)J Gene.",
            "message": "Ideal > 60%. This can indicate poor specificity of the V(D)J enrichment, use of the wrong germline reference, or the use of an unsupported chemistry type (e.g., using Single Cell 3' for V(D)J assembly). Application performance is likely to be affected.",
            "test": "< 0.50",
        },
        "warn": {
            "title": "Low Fraction of Reads Mapped to Any V(D)J Gene.",
            "message": "Ideal > 60%. This can indicate poor specificity of the V(D)J enrichment, use of the wrong germline reference, or the use of an unsupported chemistry type (e.g., using Single Cell 3' for V(D)J assembly). Application performance may be affected.",
            "test": "< 0.60",
        },
    },
    # Note: alerts don't support programmatic prefixes, so we have to hardcode TRA/TRB here.
    {
        "name": "TRA_vdj_assembly_umis_per_cell_median",
        "format": "%0.1f",
        "warn": {
            "title": "Zero Median TRA UMIs per Cell",
            "message": "Ideal > 0. This can indicate cells with extremely low TRA expression, poor cell quality, low yield from the RT reaction, or the use of an unsupported chemistry type (e.g., using Single Cell 3' for V(D)J assembly). Application performance may be affected.",
            "test": "< 0.1",
        },
        "filters": [{"chain_type": "TR"}],
    },
    {
        "name": "TRB_vdj_assembly_umis_per_cell_median",
        "format": "%0.1f",
        "warn": {
            "title": "Zero Median TRB UMIs per Cell",
            "message": "Ideal > 0. This can indicate cells with extremely low TRB expression, poor cell quality, low yield from the RT reaction, or the use of an unsupported chemistry type (e.g., using Single Cell 3' for V(D)J assembly). Application performance may be affected.",
            "test": "< 0.1",
        },
        "filters": [{"chain_type": "TR"}],
    },
    {
        "name": "IGH_vdj_assembly_umis_per_cell_median",
        "format": "%0.1f",
        "warn": {
            "title": "Zero Median IGH UMIs per Cell",
            "message": "Ideal > 0. This can indicate cells with extremely low IGH expression, poor cell quality, low yield from the RT reaction, or the use of an unsupported chemistry type (e.g., using Single Cell 3' for V(D)J assembly). Application performance may be affected.",
            "test": "< 0.1",
        },
        "filters": [{"chain_type": "IG"}],
    },
    {
        "name": "IGK_vdj_assembly_umis_per_cell_median",
        "format": "%0.1f",
        "warn": {
            "title": "Zero Median IGK UMIs per Cell",
            "message": "Ideal > 0. This can indicate cells with extremely low IGK expression, poor cell quality, low yield from the RT reaction, or the use of an unsupported chemistry type (e.g., using Single Cell 3' for V(D)J assembly). Application performance may be affected.",
            "test": "< 0.1",
        },
        "filters": [{"chain_type": "IG"}],
    },
    {
        "name": "IGL_vdj_assembly_umis_per_cell_median",
        "format": "%0.1f",
        "warn": {
            "title": "Zero Median IGL UMIs per Cell",
            "message": "Ideal > 0. This can indicate cells with extremely low IGL expression, poor cell quality, low yield from the RT reaction, or the use of an unsupported chemistry type (e.g., using Single Cell 3' for V(D)J assembly). Application performance may be affected.",
            "test": "< 0.1",
        },
        "filters": [{"chain_type": "IG"}],
    },
    {
        "name": GOOD_BCS_METRIC["name"],
        "format": GOOD_BCS_METRIC["format"],
        "error": {
            "title": "Low Fraction Valid Barcodes",
            "message": "Ideal > 85%. This usually indicates a quality issue with the Ilumina R1 read. Application performance is likely to be affected.",
            "test": "< 0.75",
        },
        "warn": {
            "title": "Low Fraction Valid Barcodes",
            "message": "Ideal > 85%. This usually indicates a quality issue with the Ilumina R1 read. Application performance may be affected.",
            "test": "< 0.85",
        },
    },
    {
        "name": "bc_bases_with_q30_frac",
        "format": "percent",
        "error": {
            "title": "Low Barcode Q30 Fraction (Illumina R1 for Single Cell V(D)J)",
            "message": "Ideal > 90%. Application performance is likely to be affected.",
            "test": "< 0.75",
        },
        "warn": {
            "title": "Low Barcode Q30 Fraction (Illumina R1 for Single Cell V(D)J)",
            "message": "Ideal > 90%. Application performance may be affected.",
            "test": "< 0.90",
        },
    },
    {
        "name": "read_bases_with_q30_frac",
        "format": "percent",
        "error": {
            "title": "Low RNA Read Q30 Fraction (Illumina R1 for Single Cell V(D)J)",
            "message": "Ideal > 85%. Application performance is likely to be affected.",
            "test": "< 0.50",
        },
        "warn": {
            "title": "Low RNA Read Q30 Fraction (Illumina R1 for Single Cell V(D)J)",
            "message": "Ideal > 85%. Application performance may be affected.",
            "test": "< 0.85",
        },
    },
    {
        "name": "read2_bases_with_q30_frac",
        "format": "percent",
        "error": {
            "title": "Low RNA Read 2 Q30 Fraction (Illumina R2 for Single Cell V(D)J)",
            "message": "Ideal > 75%. Application performance is likely to be affected.",
            "test": "< 0.50",
        },
        "warn": {
            "title": "Low RNA Read 2 Q30 Fraction (Illumina R2 for Single Cell V(D)J)",
            "message": "Ideal > 75%. Application performance may be affected.",
            "test": "< 0.75",
        },
    },
    {
        "name": "umi_bases_with_q30_frac",
        "format": "percent",
        "error": {
            "title": "Low UMI Q30 Fraction (Illumina R1 for Single Cell V(D)J)",
            "message": "Ideal > 90%. Application performance is likely to be affected.",
            "test": "< 0.75",
        },
        "warn": {
            "title": "Low UMI Q30 Fraction (Illumina R1 for Single Cell V(D)J)",
            "message": "Ideal > 90%. Application performance may be affected.",
            "test": "< 0.90",
        },
    },
]
