#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

from __future__ import annotations

import cellranger.analysis.constants as analysis_constants
import cellranger.rna.library as rna_library
import cellranger.webshim.constants.shared as shared

REPORT_PREFIX_CRISPR = rna_library.get_library_type_metric_prefix(rna_library.CRISPR_LIBRARY_TYPE)
DISPLAY_PREFIX_CRISPR = rna_library.get_library_type_report_prefix(rna_library.CRISPR_LIBRARY_TYPE)

REPORT_PREFIX_ANTIBODY = rna_library.get_library_type_metric_prefix(
    rna_library.ANTIBODY_LIBRARY_TYPE
)
DISPLAY_PREFIX_ANTIBODY = rna_library.get_library_type_report_prefix(
    rna_library.ANTIBODY_LIBRARY_TYPE
)

REPORT_PREFIX_MULTIPLEXING = rna_library.get_library_type_metric_prefix(
    rna_library.MULTIPLEXING_LIBRARY_TYPE
)
DISPLAY_PREFIX_MULTIPLEXING = rna_library.get_library_type_report_prefix(
    rna_library.MULTIPLEXING_LIBRARY_TYPE
)

REPORT_PREFIX_CUSTOM = rna_library.get_library_type_metric_prefix(rna_library.CUSTOM_LIBRARY_TYPE)
DISPLAY_PREFIX_CUSTOM = rna_library.get_library_type_report_prefix(rna_library.CUSTOM_LIBRARY_TYPE)

MULTIPLEXING_NAME = "Tag"

PCA_PRCT_CLIP = [0.1, 99.9]
TSNE_TOTALCOUNTS_PRCT_CLIP = [5, 95]
TSNE_TOTALCOUNTS_DESCRIPTION = "UMI counts"
TSNE_GEMGROUPS_DESCRIPTION = "Gem group"
TSNE_CLUSTER_DESCRIPTION = ""
MAX_TOP_N_GENES = 50
MAX_DE_TABLE_ENTRIES = 10000
TOP_DE_GENES_MIN_MEAN = 1.0
MAX_WEBSHIM_BCS_DIM = 10000
MAX_WEBSHIM_KMEANS_K = 20
CLUSTERS_FILTER_TITLE = "Clustering Type"
PVALUE_DEEMPHASIS_CUTOFF = 0.1
SATURATION_LINE = 0.9

GEM_CALL_LABELS = [
    {
        "key": analysis_constants.GEM_CLASS_GENOME0,
        "label": "genome0",
        "color": "rgb(88,165,50)",
    },
    {
        "key": analysis_constants.GEM_CLASS_GENOME1,
        "label": "genome1",
        "color": "rgb(0,161,223)",
    },
    {
        "key": analysis_constants.GEM_CLASS_MULTIPLET,
        "label": "Multiplet",
        "color": "rgb(212,212,212)",
    },
]

TOTAL_READS_METRIC = {
    "name": "total_read_pairs",
    "display_name": "Number of Reads",
    "description": "Total number of read pairs that were assigned to this library in demultiplexing.",
    "format": "integer",
}

PRENORM_READS_METRIC = {
    "name": "pre_normalization_total_reads",
    "display_name": "Pre-Normalization Number of Reads",
    "description": "Total number of read pairs that were assigned to these libraries in demultiplexing.",
    "format": "integer",
}
POSTNORM_READS_METRIC = {
    "name": "post_normalization_total_reads",
    "display_name": "Post-Normalization Number of Reads",
    "description": "Number of read pairs after normalizing for depth among multiple libraries.",
    "format": "integer",
}

GOOD_BCS_METRIC = {
    "name": "good_bc_frac",
    "display_name": "Valid Barcodes",
    "description": "Fraction of reads with barcodes that match the whitelist after barcode correction.",
    "format": "percent",
}

GOOD_UMIS_METRIC = {
    "name": "good_umi_frac",
    "display_name": "Valid UMIs",
    "description": "Fraction of reads with valid UMIs; i.e., UMI sequences that do not contain Ns and that are not homopolymers.",
    "format": "percent",
}

GENOME_MAPPED_READS_METRIC = {
    "name": "genome_mapped_reads_frac",
    "display_name": "Reads Mapped to Genome",
    "description": "Fraction of reads that mapped to the %s genome.",
    "format": "percent",
    "prefix": "genomes",
    "hidden": "len(genomes) <= 1",
}

MULTI_GENOME_MAPPED_READS_METRIC = {
    "name": "multi_genome_mapped_reads_frac",
    "display_name": "Reads Mapped to Genome",
    "description": "Fraction of reads that mapped to the genome.",
    "format": "percent",
}

GENOME_CONF_MAPPED_READS_METRIC = {
    "name": "genome_conf_mapped_reads_frac",
    "display_name": "Reads Mapped Confidently to Genome",
    "description": "Fraction of reads that mapped uniquely to the %s genome. If a gene mapped to exonic loci from a single gene and also to non-exonic loci, it is considered uniquely mapped to one of the exonic loci.",
    "format": "percent",
    "prefix": "genomes",
    "hidden": "len(genomes) <= 1",
}

MULTI_GENOME_CONF_MAPPED_READS_METRIC = {
    "name": "multi_genome_conf_mapped_reads_frac",
    "display_name": "Reads Mapped Confidently to Genome",
    "description": "Fraction of reads that mapped uniquely to the %s genome. If a gene mapped to exonic loci from a single gene and also to non-exonic loci, it is considered uniquely mapped to one of the exonic loci.",
    "format": "percent",
}

TRANSCRIPTOME_CONF_MAPPED_READS_METRIC = {
    "name": "transcriptome_conf_mapped_reads_frac",
    "display_name": "Reads Mapped Confidently to Transcriptome",
    "description": "Fraction of reads that mapped to a unique gene in the %s transcriptome. The read must be consistent with annotated splice junctions. These reads are considered for UMI counting.",
    "format": "percent",
    "prefix": "genomes",
    "hidden": "len(genomes) <= 1",
}

MULTI_TRANSCRIPTOME_CONF_MAPPED_READS_METRIC = {
    "name": "multi_transcriptome_conf_mapped_reads_frac",
    "display_name": "Reads Mapped Confidently to Transcriptome",
    "description": "Fraction of reads that mapped to a unique gene in the transcriptome. The read must be consistent with annotated splice junctions. These reads are considered for UMI counting.",
    "format": "percent",
}

TARGETED_TRANSCRIPTOME_CONF_MAPPED_READS_METRIC = {
    "name": "multi_frac_conf_transcriptomic_reads_on_target",
    "display_name": "Reads Mapped Confidently to the Targeted Transcriptome",
    "description": "Fraction of reads that mapped to a unique and targeted gene in the transcriptome. The read must be consistent with annotated splice junctions. These reads are considered for UMI counting.",
    "format": "percent",
}

OFF_TARGET_TRANSCRIPTOME_CONF_MAPPED_READS_METRIC = {
    "name": "multi_frac_conf_transcriptomic_reads_off_target",
    "display_name": "Reads Mapped Confidently to the Non-Targeted Transcriptome",
    "description": "Fraction of reads that mapped to a unique and non-targeted gene in the transcriptome. The read must be consistent with annotated splice junctions. These reads are considered for UMI counting.",
    "format": "percent",
}

ANTISENSE_CONF_MAPPED_READS_METRIC = {
    "name": "multi_antisense_reads_frac",
    "display_name": "Reads Mapped Antisense to Gene",
    "description": "Fraction of reads confidently mapped to the transcriptome, but on the opposite strand of their annotated gene. A read is counted as antisense if it has any alignments that are consistent with an exon of a transcript but antisense to it, and has no sense alignments.",
    "format": "percent",
}

INTERGENIC_CONF_MAPPED_READS_METRIC = {
    "name": "intergenic_conf_mapped_reads_frac",
    "display_name": "Reads Mapped Confidently to Intergenic Regions",
    "description": "Fraction of reads that mapped uniquely to an intergenic region of the %s genome.",
    "format": "percent",
    "prefix": "genomes",
    "hidden": "len(genomes) <= 1",
}

MULTI_INTERGENIC_CONF_MAPPED_READS_METRIC = {
    "name": "multi_intergenic_conf_mapped_reads_frac",
    "display_name": "Reads Mapped Confidently to Intergenic Regions",
    "description": "Fraction of reads that mapped uniquely to an intergenic region of the genome.",
    "format": "percent",
}

INTRONIC_CONF_MAPPED_READS_METRIC = {
    "name": "intronic_conf_mapped_reads_frac",
    "display_name": "Reads Mapped Confidently to Intronic Regions",
    "description": "Fraction of reads that mapped uniquely to an intronic region of the %s genome.",
    "format": "percent",
    "prefix": "genomes",
    "hidden": "len(genomes) <= 1",
}

MULTI_INTRONIC_CONF_MAPPED_READS_METRIC = {
    "name": "multi_intronic_conf_mapped_reads_frac",
    "display_name": "Reads Mapped Confidently to Intronic Regions",
    "description": "Fraction of reads that mapped uniquely to an intronic region of the genome.",
    "format": "percent",
}

EXONIC_CONF_MAPPED_READS_METRIC = {
    "name": "exonic_conf_mapped_reads_frac",
    "display_name": "Reads Mapped Confidently to Exonic Regions",
    "description": "Fraction of reads that mapped uniquely to an exonic region of the %s genome.",
    "format": "percent",
    "prefix": "genomes",
    "hidden": "len(genomes) <= 1",
}

MULTI_EXONIC_CONF_MAPPED_READS_METRIC = {
    "name": "multi_exonic_conf_mapped_reads_frac",
    "display_name": "Reads Mapped Confidently to Exonic Regions",
    "description": "Fraction of reads that mapped uniquely to an exonic region of the genome.",
    "format": "percent",
}

NUMBER_OF_DETECTED_CELLS_METRIC = {
    "name": "filtered_bcs",
    "display_name": "Estimated Number of Cell Partitions",
    "description": "The number of barcodes associated with %s cell-containing partitions, estimated from the barcode UMI count distribution.",
    "format": "integer",
    "prefix": "genomes",
    "hidden": "len(genomes) <= 1",
}

MULTI_NUMBER_OF_DETECTED_CELLS_METRIC = {
    "name": "filtered_bcs_transcriptome_union",
    "display_name": "Estimated Number of Cells",
    "description": "The total number of barcodes associated with %s cell-containing partitions, estimated from the barcode count distribution.",
    "format": "integer",
}

MULTI_NUMBER_OF_ANTIBODY_DETECTED_CELLS_METRIC = {
    "name": REPORT_PREFIX_ANTIBODY + "filtered_bcs_transcriptome_union",
    "display_name": "Estimated Number of Cells",
    "description": "The total number of barcodes associated with %s cell-containing partitions, estimated from the barcode antibody count distribution.",
    "format": "integer",
}

READS_PER_DETECTED_CELL_METRIC = {
    "name": "multi_transcriptome_total_raw_reads_per_filtered_bc",
    "display_name": "Mean Reads per Cell",
    "description": "The total number of sequenced reads divided by the number of barcodes associated with cell-containing partitions.",
    "format": "integer",
}

TARGETED_READS_PER_DETECTED_CELL_METRIC = {
    "name": "total_targeted_reads_per_filtered_bc",
    "display_name": "Mean Targeted Reads per Cell",
    "description": "The total number of targeted reads in cells divided by the number of barcodes associated with cell-containing partitions.",
    "format": "integer",
}

PRENORM_READS_PER_CELL_METRIC = {
    "name": "pre_normalization_multi_transcriptome_total_raw_reads_per_filtered_bc",
    "display_name": "Pre-Normalization Mean Reads per Cell",
    "description": "The total number of sequenced reads divided by the number of barcodes associated with cell-containing partitions.",
    "format": "integer",
}
POSTNORM_READS_PER_CELL_METRIC = {
    "name": "post_normalization_multi_transcriptome_total_raw_reads_per_filtered_bc",
    "display_name": "Post-Normalization Mean Reads per Cell",
    "description": "The total number of (depth-normalized) reads divided by the number of barcodes associated with cell-containing partitions.",
    "format": "integer",
}

GENES_PER_DETECTED_CELL_METRIC = {
    "name": "filtered_bcs_median_unique_genes_detected",
    "display_name": "Median Genes per Cell",
    "description": "The median number of genes detected per %s cell-associated barcode. Detection is defined as the presence of at least 1 UMI count.",
    "format": "integer",
    "prefix": "genomes",
    "hidden": "target_set is not None",
}

TARGETED_GENES_PER_DETECTED_CELL_METRIC = {
    "name": "median_genes_per_cell_on_target",
    "display_name": "Median Targeted Genes per Cell",
    "description": "The median number of targeted genes detected per cell-associated barcode. Detection is defined as the presence of at least 1 UMI count.",
    "format": "integer",
}

COUNTS_PER_DETECTED_CELL_METRIC = {
    "name": "filtered_bcs_median_counts",
    "display_name": "Median UMI Counts per Cell",
    "description": "The median number of UMI counts per %s cell-associated barcode.",
    "format": "integer",
    "prefix": "genomes",
    "hidden": "target_set is not None",
}

TARGETED_COUNTS_PER_DETECTED_CELL_METRIC = {
    "name": "median_umis_per_cell_on_target",
    "display_name": "Median Targeted UMI Counts per Cell",
    "description": "The median number of Targeted UMI counts per cell-associated barcode.",
    "format": "integer",
}

# Deprecated
CDNA_DUPE_FRAC_METRIC = {
    "name": "multi_cdna_pcr_dupe_reads_frac",
    "display_name": "cDNA PCR Duplication",
    "description": "The fraction of candidate reads that can be attributed to the cDNA amplification process; this is the number of non-unique candidate reads divided by the total number of candidate reads. A read is a candidate if it is confidently mapped and has a valid UMI and cell barcode.",
    "format": "percent",
    "hidden": "True",
}
SEQUENCING_SATURATION_METRIC = {
    "name": "multi_cdna_pcr_dupe_reads_frac",
    "display_name": "Sequencing Saturation",
    "description": 'The fraction of reads originating from an already-observed UMI. This is a function of library complexity and sequencing depth. More specifically, this is the fraction of confidently mapped, valid cell-barcode, valid UMI reads that had a non-unique (cell-barcode, UMI, gene). This metric was called "cDNA PCR Duplication" in versions of Cell Ranger prior to 1.2.',
    "format": "percent",
    "hidden": "target_set is not None",
}
TARGETED_SEQUENCING_SATURATION_METRIC = {
    "name": "multi_cdna_pcr_dupe_reads_frac_on_target",
    "display_name": "Targeted Sequencing Saturation",
    "description": 'The fraction of targeted reads originating from an already-observed targeted UMI. This is a function of library complexity and sequencing depth. More specifically, this is the fraction of confidently mapped, valid cell-barcode, valid targeted UMI reads that had a non-unique (cell-barcode, UMI, gene). This metric was called "cDNA PCR Duplication" in versions of Cell Ranger prior to 1.2.',
    "format": "percent",
    "hidden": "target_set is None",
}

SUMMARY_GENES_PER_DETECTED_CELL_METRIC = GENES_PER_DETECTED_CELL_METRIC.copy()
SUMMARY_GENES_PER_DETECTED_CELL_METRIC["hidden"] = "len(genomes) > 1 or target_set is not None"

# These ultimately get displayed
SUMMARY_METRICS = [
    MULTI_NUMBER_OF_DETECTED_CELLS_METRIC,
    READS_PER_DETECTED_CELL_METRIC,
    POSTNORM_READS_PER_CELL_METRIC,
    SUMMARY_GENES_PER_DETECTED_CELL_METRIC,
]

MULTIPLET_RATE_METRIC = {
    "name": "filtered_bcs_inferred_multiplet_rate",
    "display_name": "Fraction GEMs with >1 Cell",
    "description": "The mean fraction of cell-associated barcodes estimated to be associated with more than one cell, calculated via bootstrap sampling and adjusting for the ratio of the two cell types.",
    "format": "percent",
}

FRAC_READS_IN_CELLS_METRIC = {
    "name": "filtered_bcs_conf_mapped_barcoded_reads_cum_frac",
    "display_name": "Fraction Reads in Cells",
    "description": "The fraction of valid-barcode, confidently-mapped-to-transcriptome reads with %s cell-associated barcodes.",
    "format": "percent",
    "prefix": "genomes",
    "hidden": "len(genomes) <= 1",
}

# Same as above, but summarized across all genomes present
MULTI_FRAC_READS_IN_CELLS_METRIC = {
    "name": "multi_filtered_bcs_conf_mapped_barcoded_reads_cum_frac",
    "display_name": "Fraction Reads in Cells",
    "description": "The fraction of valid-barcode, confidently-mapped-to-transcriptome reads with cell-associated barcodes.",
    "format": "percent",
}

MEAN_USABLE_READS_IN_CELLS_METRIC = {
    "name": "filtered_bcs_conf_mapped_barcoded_reads_per_filtered_bc",
    "display_name": "Mean Usable Reads per Cell",
    "description": "The total number of reads that have a valid barcode, valid UMI, confidently map to the transcriptome and are associated with a cell-containing GEM, divided by the number of cells.",
    "format": "integer",
    "prefix": "genomes",
}


TOTAL_GENES_DETECTED_METRIC = {
    "name": "filtered_bcs_total_unique_genes_detected",
    "display_name": "Total Genes Detected",
    "description": "The number of %s genes with at least one UMI count in any cell.",
    "format": "integer",
    "prefix": "genomes",
    "hidden": "target_set is not None",
}

TOTAL_TARGETED_GENES_DETECTED_METRIC = {
    "name": "num_genes_detected_on_target",
    "display_name": "Total Targeted Genes Detected",
    "description": "The number of targeted genes with at least one UMI count in any cell.",
    "format": "integer",
}

Q30_METRICS = [
    {
        "name": "bc_bases_with_q30_frac",
        "display_name": "Q30 Bases in Barcode",
        "description": "Fraction of cell barcode bases with Q-score >= 30, excluding very low quality/no-call (Q <= 2) bases from the denominator.",
        "format": "percent",
    },
    {
        "name": "read_bases_with_q30_frac",
        "display_name": "Q30 Bases in RNA Read",
        "description": "Fraction of RNA read bases with Q-score >= 30, excluding very low quality/no-call (Q <= 2) bases from the denominator. This is Read 1 for the Single Cell 3' v1 chemistry and Read 2 for the Single Cell 3' v2 chemistry.",
        "format": "percent",
    },
    {
        "name": "read2_bases_with_q30_frac",
        "display_name": "Q30 Bases in RNA Read 2",
        "description": "Fraction of RNA read 2 bases with Q-score >= 30, excluding very low quality/no-call (Q <= 2) bases from the denominator.",
        "format": "percent",
    },
    {
        "name": "umi_bases_with_q30_frac",
        "display_name": "Q30 Bases in UMI",
        "description": "Fraction of UMI bases with Q-score >= 30, excluding very low quality/no-call (Q <= 2) bases from the denominator.",
        "format": "percent",
    },
]

GEX_DUAL_INDEXED_Q30_METRICS = [
    {
        "name": "i1_bases_with_q30_frac",
        "display_name": "Q30 Bases in Sample Index I1",
        "description": "Fraction of I1 (or I7) read bases with Q-score >= 30, excluding very low quality/no-call (Q <= 2) bases from the denominator.",
        "format": "percent",
    },
    {
        "name": "i2_bases_with_q30_frac",
        "display_name": "Q30 Bases in Sample Index I2",
        "description": "Fraction of I2 (or I5) read bases with Q-score >= 30, excluding very low quality/no-call (Q <= 2) bases from the denominator.",
        "format": "percent",
    },
]

CRISPR_Q30_METRICS = [
    {
        "name": REPORT_PREFIX_CRISPR + "bc_bases_with_q30_frac",
        "display_name": DISPLAY_PREFIX_CRISPR + "Q30 Bases in Barcode",
        "description": "Fraction of cell barcode bases with Q-score >= 30, excluding very low quality/no-call (Q <= 2) bases from the denominator.",
        "format": "percent",
    },
    {
        "name": REPORT_PREFIX_CRISPR + "read_bases_with_q30_frac",
        "display_name": DISPLAY_PREFIX_CRISPR + "Q30 Bases in CRISPR Read",
        "description": "Fraction of CRISPR library read bases with Q-score >= 30, excluding very low quality/no-call (Q <= 2) bases from the denominator. This is Read 2 for the Single Cell 3' v3 and Single Cell 5' chemistries.",
        "format": "percent",
    },
    {
        "name": REPORT_PREFIX_CRISPR + "read2_bases_with_q30_frac",
        "display_name": DISPLAY_PREFIX_CRISPR + "Q30 Bases in CRISPR Read 2",
        "description": "Fraction of CRISPR library read 2 bases with Q-score >= 30, excluding very low quality/no-call (Q <= 2) bases from the denominator.",
        "format": "percent",
    },
    {
        "name": REPORT_PREFIX_CRISPR + "umi_bases_with_q30_frac",
        "display_name": DISPLAY_PREFIX_CRISPR + "Q30 Bases in UMI",
        "description": "Fraction of CRISPR library UMI bases with Q-score >= 30, excluding very low quality/no-call (Q <= 2) bases from the denominator.",
        "format": "percent",
    },
]

CRISPR_SEQUENCING_METRICS = [
    {
        "name": REPORT_PREFIX_CRISPR + "total_read_pairs",
        "display_name": DISPLAY_PREFIX_CRISPR + "Number of Reads",
        "description": "Total number of CRISPR library reads",
        "format": "integer",
    },
    {
        "name": REPORT_PREFIX_CRISPR + "reads_per_cell",
        "display_name": DISPLAY_PREFIX_CRISPR + "Mean Reads per Cell",
        "description": "The total number of sequenced CRISPR library reads divided by the number of barcodes associated with cell-containing partitions.",
        "format": "integer",
    },
    {
        "name": REPORT_PREFIX_CRISPR + "good_bc_frac",
        "display_name": DISPLAY_PREFIX_CRISPR + "Valid Barcodes",
        "description": "Fraction of CRISPR library reads with a barcode found in or corrected to one that is found in the whitelist",
        "format": "percent",
    },
    {
        "name": REPORT_PREFIX_CRISPR + "good_umi_frac",
        "display_name": DISPLAY_PREFIX_CRISPR + "Valid UMIs",
        "description": "Fraction of CRISPR library reads with valid UMIs; i.e., UMI sequences that do not contain Ns and that are not homopolymers.",
        "format": "percent",
    },
    {
        "name": REPORT_PREFIX_CRISPR + "multi_cdna_pcr_dupe_reads_frac",
        "display_name": DISPLAY_PREFIX_CRISPR + "Sequencing Saturation",
        "description": "The fraction of CRISPR library reads originating from an already-observed UMI. This is a function of library complexity and sequencing depth. More specifically, this is the fraction of confidently mapped, valid cell-barcode, valid UMI reads that had a non-unique (cell-barcode, UMI, CRISPR feature barcode).",
        "format": "percent",
    },
] + CRISPR_Q30_METRICS

CRISPR_APPLICATION_METRICS = [
    {
        "name": REPORT_PREFIX_CRISPR + "feature_bc_extracted_frac",
        "display_name": DISPLAY_PREFIX_CRISPR + "Fraction Reads with Putative Protospacer Sequence",
        "description": "Fraction of CRISPR library reads from which a putative protospacer sequence could be extracted",
        "format": "percent",
    },
    {
        "name": REPORT_PREFIX_CRISPR + "recognized_feature_bc_frac",
        "display_name": DISPLAY_PREFIX_CRISPR + "Fraction Guide Reads",
        "description": "Fraction of CRISPR library reads with a recognized protospacer sequence",
        "format": "percent",
    },
    {
        "name": REPORT_PREFIX_CRISPR + "frac_feature_reads_usable",
        "display_name": DISPLAY_PREFIX_CRISPR + "Fraction Guide Reads Usable",
        "description": "Fraction of CRISPR library reads with a recognized protospacer sequence, a valid UMI, and a cell-associated barcode",
        "format": "percent",
    },
    {
        "name": REPORT_PREFIX_CRISPR + "feature_reads_usable_per_cell",
        "display_name": DISPLAY_PREFIX_CRISPR + "Guide Reads Usable per Cell",
        "description": "Number of CRISPR library guide reads usable divided by the number of cell-associated barcodes",
        "format": "integer",
    },
    {
        "name": REPORT_PREFIX_CRISPR + "unrecognized_feature_bc_frac",
        "display_name": DISPLAY_PREFIX_CRISPR + "Fraction Protospacer Not Recognized",
        "description": "Among all CRISPR library reads with a putative protospacer sequence, the fraction with a protospacer sequence that was not recognized",
        "format": "percent",
    },
    {
        "name": REPORT_PREFIX_CRISPR + "feature_reads_in_cells",
        "display_name": DISPLAY_PREFIX_CRISPR + "Guide Reads in Cells",
        "description": "Among CRISPR library reads with a recognized protospacer sequence, a valid UMI, and a valid barcode, the fraction associated with cell-containing partitions",
        "format": "percent",
    },
    {
        "name": REPORT_PREFIX_CRISPR + "frac_cells_with_protospacer",
        "display_name": DISPLAY_PREFIX_CRISPR + "Cells with 1 or more protospacers detected",
        "description": "Cells with 1 or more protospacers detected",
        "format": "percent",
    },
    {
        "name": REPORT_PREFIX_CRISPR + "frac_cells_with_multiple_protospacer",
        "display_name": DISPLAY_PREFIX_CRISPR + "Cells with 2 or more protospacers detected",
        "description": "Cells with 2 or more protospacers detected",
        "format": "percent",
    },
    {
        "name": REPORT_PREFIX_CRISPR + "multi_filtered_bcs_median_counts",
        "display_name": DISPLAY_PREFIX_CRISPR
        + "Median UMIs per Cell (summed over all recognized protospacers)",
        "description": "Median UMIs per Cell (summed over all recognized protospacers)",
        "format": "integer",
    },
]


MULTIPLEXING_APPLICATION_METRICS = [
    {
        "name": REPORT_PREFIX_MULTIPLEXING + "recognized_feature_bc_frac",
        "display_name": DISPLAY_PREFIX_MULTIPLEXING + f"Fraction {MULTIPLEXING_NAME} Reads",
        "description": f"Fraction of {MULTIPLEXING_NAME} library reads with a recognized {MULTIPLEXING_NAME} sequence",
        "format": "percent",
    },
    {
        "name": REPORT_PREFIX_MULTIPLEXING + "frac_feature_reads_usable",
        "display_name": DISPLAY_PREFIX_MULTIPLEXING + f"Fraction {MULTIPLEXING_NAME} Reads Usable",
        "description": f"Fraction of {MULTIPLEXING_NAME} library reads with a recognized {MULTIPLEXING_NAME} sequence, a valid UMI, and a cell-associated barcode",
        "format": "percent",
    },
    {
        "name": REPORT_PREFIX_MULTIPLEXING + "feature_reads_usable_per_cell",
        "display_name": DISPLAY_PREFIX_MULTIPLEXING + f"{MULTIPLEXING_NAME} Reads Usable per Cell",
        "description": f"Number of {MULTIPLEXING_NAME} library guide reads usable divided by the number of cell-associated barcodes",
        "format": "integer",
    },
    {
        "name": REPORT_PREFIX_MULTIPLEXING + "unrecognized_feature_bc_frac",
        "display_name": DISPLAY_PREFIX_MULTIPLEXING
        + f"Fraction {MULTIPLEXING_NAME} Not Recognized",
        "description": f"Among all {MULTIPLEXING_NAME} library reads, the fraction with a {MULTIPLEXING_NAME} sequence that was not recognized",
        "format": "percent",
    },
    {
        "name": REPORT_PREFIX_MULTIPLEXING + "feature_reads_in_cells",
        "display_name": DISPLAY_PREFIX_MULTIPLEXING + f"{MULTIPLEXING_NAME} Reads in Cells",
        "description": f"Among {MULTIPLEXING_NAME} library reads with a recognized {MULTIPLEXING_NAME} sequence, a valid UMI, and a valid barcode, the fraction associated with cell-containing partitions",
        "format": "percent",
    },
    {
        "name": REPORT_PREFIX_MULTIPLEXING + f"frac_cells_with_{MULTIPLEXING_NAME}",
        "display_name": DISPLAY_PREFIX_MULTIPLEXING
        + f"Cells with 1 or more {MULTIPLEXING_NAME} assigned",
        "description": f"Cells with 1 or more {MULTIPLEXING_NAME} assigned",
        "format": "percent",
    },
    {
        "name": REPORT_PREFIX_MULTIPLEXING + f"frac_cells_with_multiple_{MULTIPLEXING_NAME}",
        "display_name": DISPLAY_PREFIX_MULTIPLEXING
        + f"Cells with 2 or more {MULTIPLEXING_NAME} assigned",
        "description": f"Cells with 2 or more {MULTIPLEXING_NAME} assigned",
        "format": "percent",
    },
    {
        "name": REPORT_PREFIX_MULTIPLEXING + "multi_filtered_bcs_median_counts",
        "display_name": DISPLAY_PREFIX_MULTIPLEXING
        + f"Median UMIs per Cell (summed over all recognized {MULTIPLEXING_NAME})",
        "description": f"Median UMIs per Cell (summed over all recognized {MULTIPLEXING_NAME})",
        "format": "integer",
    },
]

ANTIBODY_CELL_CALLING_METRICS = [MULTI_NUMBER_OF_ANTIBODY_DETECTED_CELLS_METRIC]

ANTIBODY_Q30_METRICS = [
    {
        "name": REPORT_PREFIX_ANTIBODY + "bc_bases_with_q30_frac",
        "display_name": DISPLAY_PREFIX_ANTIBODY + "Q30 Bases in Barcode",
        "description": "Fraction of Antibody library cell barcode bases with Q-score >= 30, excluding very low quality/no-call (Q <= 2) bases from the denominator.",
        "format": "percent",
    },
    {
        "name": REPORT_PREFIX_ANTIBODY + "read_bases_with_q30_frac",
        "display_name": DISPLAY_PREFIX_ANTIBODY + "Q30 Bases in Antibody Read",
        "description": "Fraction of Antibody library read bases with Q-score >= 30, excluding very low quality/no-call (Q <= 2) bases from the denominator. This is Read 2 for the Single Cell 3' v3 and Single Cell 5' chemistries.",
        "format": "percent",
    },
    {
        "name": REPORT_PREFIX_ANTIBODY + "read2_bases_with_q30_frac",
        "display_name": DISPLAY_PREFIX_ANTIBODY + "Q30 Bases in Antibody Read 2",
        "description": "Fraction of Antibody library read 2 bases with Q-score >= 30, excluding very low quality/no-call (Q <= 2) bases from the denominator.",
        "format": "percent",
    },
    {
        "name": REPORT_PREFIX_ANTIBODY + "umi_bases_with_q30_frac",
        "display_name": DISPLAY_PREFIX_ANTIBODY + "Q30 Bases in UMI",
        "description": "Fraction of Antibody library UMI bases with Q-score >= 30, excluding very low quality/no-call (Q <= 2) bases from the denominator.",
        "format": "percent",
    },
]


ANTIBODY_SEQUENCING_METRICS = [
    {
        "name": REPORT_PREFIX_ANTIBODY + "total_read_pairs",
        "display_name": DISPLAY_PREFIX_ANTIBODY + "Number of Reads",
        "description": "Total number of Antibody library reads",
        "format": "integer",
    },
    {
        "name": REPORT_PREFIX_ANTIBODY + "reads_per_cell",
        "display_name": DISPLAY_PREFIX_ANTIBODY + "Mean Reads per Cell",
        "description": "The total number of sequenced Antibody library reads divided by the number of barcodes associated with cell-containing partitions.",
        "format": "integer",
    },
    {
        "name": REPORT_PREFIX_ANTIBODY + "good_bc_frac",
        "display_name": DISPLAY_PREFIX_ANTIBODY + "Valid Barcodes",
        "description": "Fraction of Antibody library reads with a barcode found in or corrected to one that is found in the whitelist",
        "format": "percent",
    },
    {
        "name": REPORT_PREFIX_ANTIBODY + "good_umi_frac",
        "display_name": DISPLAY_PREFIX_ANTIBODY + "Valid UMIs",
        "description": "Fraction of Antibody library reads with valid UMIs; i.e., UMI sequences that do not contain Ns and that are not homopolymers.",
        "format": "percent",
    },
    {
        "name": REPORT_PREFIX_ANTIBODY + "multi_cdna_pcr_dupe_reads_frac",
        "display_name": DISPLAY_PREFIX_ANTIBODY + "Sequencing Saturation",
        "description": "The fraction of Antibody library reads originating from an already-observed UMI. This is a function of library complexity and sequencing depth. More specifically, this is the fraction of confidently mapped, valid cell-barcode, valid UMI reads that had a non-unique (cell-barcode, UMI, antibody feature barcode).",
        "format": "percent",
    },
] + ANTIBODY_Q30_METRICS

ANTIBODY_APPLICATION_METRICS = [
    {
        "name": REPORT_PREFIX_ANTIBODY + "recognized_feature_bc_frac",
        "display_name": DISPLAY_PREFIX_ANTIBODY + "Fraction Antibody Reads",
        "description": "Fraction of Antibody library reads that contain a recognized antibody barcode",
        "format": "percent",
    },
    {
        "name": REPORT_PREFIX_ANTIBODY + "frac_feature_reads_usable",
        "display_name": DISPLAY_PREFIX_ANTIBODY + "Fraction Antibody Reads Usable",
        "description": "Fraction of Antibody library reads that contain a recognized antibody barcode, a valid UMI, and a cell-associated barcode",
        "format": "percent",
    },
    {
        "name": REPORT_PREFIX_ANTIBODY + "feature_reads_usable_per_cell",
        "display_name": DISPLAY_PREFIX_ANTIBODY + "Antibody Reads Usable per Cell",
        "description": "Number of Antibody library reads usable divided by the number of cell-associated barcodes",
        "format": "integer",
    },
    {
        "name": REPORT_PREFIX_ANTIBODY + "reads_lost_to_aggregate_GEMs",
        "display_name": DISPLAY_PREFIX_ANTIBODY + "Fraction Antibody Reads in Aggregate Barcodes",
        "description": "Fraction of Antibody library reads that were lost after removing aggregate barcodes.",
        "format": "percent",
    },
    {
        "name": REPORT_PREFIX_ANTIBODY + "unrecognized_feature_bc_frac",
        "display_name": DISPLAY_PREFIX_ANTIBODY + "Fraction Unrecognized Antibody",
        "description": "Among all Antibody library reads, the fraction with an unrecognizable antibody barcode",
        "format": "percent",
    },
    {
        "name": REPORT_PREFIX_ANTIBODY + "feature_reads_in_cells",
        "display_name": DISPLAY_PREFIX_ANTIBODY + "Antibody Reads in Cells",
        "description": "Among Antibody library reads with a recognized antibody barcode, a valid UMI, and a valid barcode, the fraction associated with cell-containing partitions",
        "format": "percent",
    },
    {
        "name": REPORT_PREFIX_ANTIBODY + "multi_filtered_bcs_median_counts",
        "display_name": DISPLAY_PREFIX_ANTIBODY
        + "Median UMIs per Cell (summed over all recognized antibody barcodes)",
        "description": "Median UMIs per Cell (summed over all recognized antibody barcodes)",
        "format": "integer",
    },
]

MULTIPLEXING_Q30_METRICS = [
    {
        "name": REPORT_PREFIX_MULTIPLEXING + "bc_bases_with_q30_frac",
        "display_name": DISPLAY_PREFIX_MULTIPLEXING + "Q30 Bases in Barcode",
        "description": "Fraction of cell barcode bases with Q-score >= 30, excluding very low quality/no-call (Q <= 2) bases from the denominator.",
        "format": "percent",
    },
    {
        "name": REPORT_PREFIX_MULTIPLEXING + "read_bases_with_q30_frac",
        "display_name": DISPLAY_PREFIX_MULTIPLEXING + "Q30 Bases in MULTIPLEXING Read",
        "description": "Fraction of MULTIPLEXING library read bases with Q-score >= 30, excluding very low quality/no-call (Q <= 2) bases from the denominator. This is Read 2 for the Single Cell 3' v3 and Single Cell 5' chemistries.",
        "format": "percent",
    },
    {
        "name": REPORT_PREFIX_MULTIPLEXING + "read2_bases_with_q30_frac",
        "display_name": DISPLAY_PREFIX_MULTIPLEXING + "Q30 Bases in MULTIPLEXING Read 2",
        "description": "Fraction of MULTIPLEXING library read 2 bases with Q-score >= 30, excluding very low quality/no-call (Q <= 2) bases from the denominator.",
        "format": "percent",
    },
    {
        "name": REPORT_PREFIX_MULTIPLEXING + "umi_bases_with_q30_frac",
        "display_name": DISPLAY_PREFIX_MULTIPLEXING + "Q30 Bases in UMI",
        "description": "Fraction of MULTIPLEXING library UMI bases with Q-score >= 30, excluding very low quality/no-call (Q <= 2) bases from the denominator.",
        "format": "percent",
    },
]

MULTIPLEXING_SEQUENCING_METRICS = [
    {
        "name": REPORT_PREFIX_MULTIPLEXING + "total_read_pairs",
        "display_name": DISPLAY_PREFIX_MULTIPLEXING + "Number of Reads",
        "description": "Total number of MULTIPLEXING library reads",
        "format": "integer",
    },
    {
        "name": REPORT_PREFIX_MULTIPLEXING + "reads_per_cell",
        "display_name": DISPLAY_PREFIX_MULTIPLEXING + "Mean Reads per Cell",
        "description": "The total number of sequenced MULTIPLEXING library reads divided by the number of barcodes associated with cell-containing partitions.",
        "format": "integer",
    },
    {
        "name": REPORT_PREFIX_MULTIPLEXING + "good_bc_frac",
        "display_name": DISPLAY_PREFIX_MULTIPLEXING + "Valid Barcodes",
        "description": "Fraction of MULTIPLEXING library reads with a barcode found in or corrected to one that is found in the whitelist",
        "format": "percent",
    },
    {
        "name": REPORT_PREFIX_MULTIPLEXING + "multi_cdna_pcr_dupe_reads_frac",
        "display_name": DISPLAY_PREFIX_MULTIPLEXING + "Sequencing Saturation",
        "description": "The fraction of MULTIPLEXING library reads originating from an already-observed UMI. This is a function of library complexity and sequencing depth. More specifically, this is the fraction of valid cell-barcode, valid UMI, recongnized MULTIPLEXING feature barcode reads that had a non-unique (cell-barcode, UMI, MULTIPLEXING feature barcode).",
        "format": "percent",
    },
] + MULTIPLEXING_Q30_METRICS


CUSTOM_FEATURE_Q30_METRICS = [
    {
        "name": "bc_bases_with_q30_frac",
        "display_name": "%s: Q30 Bases in Barcode",
        "description": "Fraction of Custom library cell barcode bases with Q-score >= 30, excluding very low quality/no-call (Q <= 2) bases from the denominator.",
        "format": "percent",
        "prefix": "custom_features",
    },
    {
        "name": "read_bases_with_q30_frac",
        "display_name": "%s: Q30 Bases in Feature Read",
        "description": "Fraction of Custom library read bases with Q-score >= 30, excluding very low quality/no-call (Q <= 2) bases from the denominator. This is Read 2 for the Single Cell 3' v3 and Single Cell 5' chemistries.",
        "format": "percent",
        "prefix": "custom_features",
    },
    {
        "name": "read2_bases_with_q30_frac",
        "display_name": "%s: Q30 Bases in Feature Read 2",
        "description": "Fraction of Custom library read 2 bases with Q-score >= 30, excluding very low quality/no-call (Q <= 2) bases from the denominator.",
        "format": "percent",
        "prefix": "custom_features",
    },
    {
        "name": "umi_bases_with_q30_frac",
        "display_name": "%s: Q30 Bases in UMI",
        "description": "Fraction of Custom library UMI bases with Q-score >= 30, excluding very low quality/no-call (Q <= 2) bases from the denominator.",
        "format": "percent",
        "prefix": "custom_features",
    },
]


CUSTOM_FEATURE_SEQUENCING_METRICS = [
    {
        "name": "total_read_pairs",
        "display_name": "%s: Number of Reads",
        "description": "Total number of Custom library reads",
        "format": "integer",
        "prefix": "custom_features",
    },
    {
        "name": "reads_per_cell",
        "display_name": "%s: Mean Reads per Cell",
        "description": "The total number of sequenced Custom library reads divided by the number of barcodes associated with cell-containing partitions.",
        "format": "integer",
        "prefix": "custom_features",
    },
    {
        "name": "good_bc_frac",
        "display_name": "%s: Valid Barcodes",
        "description": "Fraction of Custom library reads with a barcode found in or corrected to one that is found in the whitelist",
        "format": "percent",
        "prefix": "custom_features",
    },
    {
        "name": "good_umi_frac",
        "display_name": "%s: Valid UMIs",
        "description": "Fraction of Custom library reads with valid UMIs; i.e., UMI sequences that do not contain Ns and that are not homopolymers.",
        "format": "percent",
        "prefix": "custom_features",
    },
    {
        "name": "multi_cdna_pcr_dupe_reads_frac",
        "display_name": "%s: Sequencing Saturation",
        "description": "The fraction of Custom library reads originating from an already-observed UMI. This is a function of library complexity and sequencing depth. More specifically, this is the fraction of confidently mapped, valid cell-barcode, valid UMI reads that had a non-unique (cell-barcode, UMI, feature barcode).",
        "format": "percent",
        "prefix": "custom_features",
    },
] + CUSTOM_FEATURE_Q30_METRICS

CUSTOM_FEATURE_APPLICATION_METRICS = [
    {
        "name": "recognized_feature_bc_frac",
        "display_name": "%s: Fraction Feature Reads",
        "description": "Fraction of Custom library reads that contain a recognized feature barcode",
        "format": "percent",
        "prefix": "custom_features",
    },
    {
        "name": "frac_feature_reads_usable",
        "display_name": "%s: Fraction Feature Reads Usable",
        "description": "Fraction of Custom library reads that contain a recognized feature barcode, a valid UMI, and a cell-associated barcode",
        "format": "percent",
        "prefix": "custom_features",
    },
    {
        "name": "feature_reads_usable_per_cell",
        "display_name": "%s: Feature Reads Usable per Cell",
        "description": "Number of feature Custom library reads usable divided by the number of cell-associated barcodes",
        "format": "integer",
        "prefix": "custom_features",
    },
    {
        "name": "unrecognized_feature_bc_frac",
        "display_name": "%s: Fraction Unrecognized Feature",
        "description": "Among all Custom library reads, the fraction with an unrecognizable feature barcode",
        "format": "percent",
        "prefix": "custom_features",
    },
    {
        "name": "feature_reads_in_cells",
        "display_name": "%s: Feature Reads in Cells",
        "description": "Among Custom library reads with a recognized feature barcode, a valid UMI, and a valid barcode, the fraction associated with cell-containing partitions",
        "format": "percent",
        "prefix": "custom_features",
    },
    {
        "name": "multi_filtered_bcs_median_counts",
        "display_name": "%s: Median UMIs per Cell (summed over all recognized feature barcodes)",
        "description": "Median UMIs per Cell (summed over all recognized feature barcodes)",
        "format": "integer",
        "prefix": "custom_features",
    },
]

SEQUENCING_METRICS = [
    TOTAL_READS_METRIC,
    PRENORM_READS_METRIC,
    POSTNORM_READS_METRIC,
    GOOD_BCS_METRIC,
    GOOD_UMIS_METRIC,
    CDNA_DUPE_FRAC_METRIC,
    SEQUENCING_SATURATION_METRIC,
    TARGETED_SEQUENCING_SATURATION_METRIC,
] + Q30_METRICS

MAPPING_METRICS = [
    MULTI_GENOME_MAPPED_READS_METRIC,
    GENOME_MAPPED_READS_METRIC,
    MULTI_GENOME_CONF_MAPPED_READS_METRIC,
    GENOME_CONF_MAPPED_READS_METRIC,
    MULTI_INTERGENIC_CONF_MAPPED_READS_METRIC,
    INTERGENIC_CONF_MAPPED_READS_METRIC,
    MULTI_INTRONIC_CONF_MAPPED_READS_METRIC,
    INTRONIC_CONF_MAPPED_READS_METRIC,
    MULTI_EXONIC_CONF_MAPPED_READS_METRIC,
    EXONIC_CONF_MAPPED_READS_METRIC,
    MULTI_TRANSCRIPTOME_CONF_MAPPED_READS_METRIC,
    TRANSCRIPTOME_CONF_MAPPED_READS_METRIC,
    TARGETED_TRANSCRIPTOME_CONF_MAPPED_READS_METRIC,
    OFF_TARGET_TRANSCRIPTOME_CONF_MAPPED_READS_METRIC,
    ANTISENSE_CONF_MAPPED_READS_METRIC,
]

DETECTED_CELL_METRICS = [
    MULTI_NUMBER_OF_DETECTED_CELLS_METRIC,
    NUMBER_OF_DETECTED_CELLS_METRIC,
    MULTI_FRAC_READS_IN_CELLS_METRIC,
    FRAC_READS_IN_CELLS_METRIC,
    READS_PER_DETECTED_CELL_METRIC,
    TARGETED_READS_PER_DETECTED_CELL_METRIC,
    PRENORM_READS_PER_CELL_METRIC,
    POSTNORM_READS_PER_CELL_METRIC,
    GENES_PER_DETECTED_CELL_METRIC,
    TARGETED_GENES_PER_DETECTED_CELL_METRIC,
    TOTAL_GENES_DETECTED_METRIC,
    TOTAL_TARGETED_GENES_DETECTED_METRIC,
    COUNTS_PER_DETECTED_CELL_METRIC,
    TARGETED_COUNTS_PER_DETECTED_CELL_METRIC,
]

# NOTE: these metrics only apply to barnyard samples
GEM_METRICS = [
    {
        "name": "filtered_bcs_transcriptome_union",
        "display_name": "GEMs with >0 Cells",
        "description": "The number of barcodes associated with at least one cell.",
        "format": "integer",
        "hidden": "len(genomes) <= 1",
    },
    {
        "name": "filtered_bcs_inferred_multiplets",
        "display_name": "GEMs with >1 Cell",
        "description": "The inferred number of barcodes that are associated with multiple cells.",
        "format": "integer",
    },
    MULTIPLET_RATE_METRIC,
]

# NOTE; these metrics only apply to aggregated samples
AGGREGATION_METRICS = [
    {
        "name": "frac_reads_kept",
        "display_name": "%s Fraction of Reads Kept",
        "description": "The fraction of reads kept for this input sample, after normalizing depth across samples to reduce batch effects.",
        "format": "percent",
        "prefix": "agg_batches",
    },
    {
        "name": "pre_normalization_raw_reads_per_filtered_bc",
        "display_name": "%s Pre-Normalization Total Reads per Cell",
        "description": "The mean number of total sequencing reads per cell in this input sample, prior to depth normalization.",
        "format": "integer",
        "prefix": "agg_batches",
    },
    {
        "name": "pre_normalization_cmb_reads_per_filtered_bc",
        "display_name": "%s Pre-Normalization Confidently Mapped Barcoded Reads per Cell",
        "description": "The mean number of confidently-mapped-to-transcriptome, valid-barcode reads per cell in this input sample, prior to depth normalization.",
        "format": "integer",
        "prefix": "agg_batches",
    },
]

CRISPR_AGGREGATION_METRICS = [
    {
        "name": REPORT_PREFIX_CRISPR + "frac_reads_kept",
        "display_name": "%s Fraction of Reads Kept",
        "description": "The fraction of reads kept for this input sample, after normalizing depth across samples to reduce batch effects.",
        "format": "percent",
        "prefix": "agg_batches",
    },
    {
        "name": REPORT_PREFIX_CRISPR + "pre_normalization_raw_reads_per_filtered_bc",
        "display_name": "%s Pre-Normalization Total Reads per Cell",
        "description": "The mean number of total sequencing reads per cell in this input sample, prior to depth normalization.",
        "format": "integer",
        "prefix": "agg_batches",
    },
    {
        "name": REPORT_PREFIX_CRISPR + "pre_normalization_cmb_reads_per_filtered_bc",
        "display_name": "%s Pre-Normalization Confidently Mapped Barcoded Reads per Cell",
        "description": "The mean number of confidently-mapped-to-transcriptome, valid-barcode reads per cell in this input sample, prior to depth normalization.",
        "format": "integer",
        "prefix": "agg_batches",
    },
]

ANTIBODY_AGGREGATION_METRICS = [
    {
        "name": REPORT_PREFIX_ANTIBODY + "frac_reads_kept",
        "display_name": "%s Fraction of Reads Kept",
        "description": "The fraction of reads kept for this input sample, after normalizing depth across samples to reduce batch effects.",
        "format": "percent",
        "prefix": "agg_batches",
    },
    {
        "name": REPORT_PREFIX_ANTIBODY + "pre_normalization_raw_reads_per_filtered_bc",
        "display_name": "%s Pre-Normalization Total Reads per Cell",
        "description": "The mean number of total sequencing reads per cell in this input sample, prior to depth normalization.",
        "format": "integer",
        "prefix": "agg_batches",
    },
    {
        "name": REPORT_PREFIX_ANTIBODY + "pre_normalization_cmb_reads_per_filtered_bc",
        "display_name": "%s Pre-Normalization Confidently Mapped Barcoded Reads per Cell",
        "description": "The mean number of confidently-mapped-to-transcriptome, valid-barcode reads per cell in this input sample, prior to depth normalization.",
        "format": "integer",
        "prefix": "agg_batches",
    },
]

CUSTOM_FEATURE_AGGREGATION_METRICS = [
    {
        "name": REPORT_PREFIX_CUSTOM + "_frac_reads_kept",
        "display_name": "%s Fraction of Reads Kept",
        "description": "The fraction of reads kept for this input sample, after normalizing depth across samples to reduce batch effects.",
        "format": "percent",
        "prefix": "agg_batches",
    },
    {
        "name": REPORT_PREFIX_CUSTOM + "_pre_normalization_raw_reads_per_filtered_bc",
        "display_name": "%s Pre-Normalization Total Reads per Cell",
        "description": "The mean number of total sequencing reads per cell in this input sample, prior to depth normalization.",
        "format": "integer",
        "prefix": "agg_batches",
    },
    {
        "name": REPORT_PREFIX_CUSTOM + "_pre_normalization_cmb_reads_per_filtered_bc",
        "display_name": "%s Pre-Normalization Confidently Mapped Barcoded Reads per Cell",
        "description": "The mean number of confidently-mapped-to-transcriptome, valid-barcode reads per cell in this input sample, prior to depth normalization.",
        "format": "integer",
        "prefix": "agg_batches",
    },
]

TARGETED_APPLICATION_METRICS = [
    {
        "name": "num_genes_on_target",
        "display_name": "Number of Targeted Genes",
        "description": "Number of Targeted Genes as specified in the input target panel file.",
        "format": "integer",
    },
    {
        "name": "num_genes_off_target",
        "display_name": "Number of Non-Targeted Genes",
        "description": "Number of genes in the reference genome that are not targeted.",
        "format": "integer",
    },
    {
        "name": "num_genes_quantifiable_on_target",
        "display_name": "Number of Targeted Genes >= 10 UMIs",
        "description": "Number of Targeted Genes with at least 10 UMIs in cell-associated barcodes. These genes are then considered in the calculation of per-gene enrichments.",
        "format": "integer",
    },
    {
        "name": "num_genes_quantifiable_off_target",
        "display_name": "Number of Non-Targeted Genes >= 10 UMIs",
        "description": "Number of Non-Targeted Genes with at least 10 UMIs in cell-associated barcodes. These genes are then considered in the calculation of per-gene enrichments.",
        "format": "integer",
    },
    {
        "name": "num_rpu_enriched_genes_on_target",
        "display_name": "Number of Enriched Targeted Genes",
        "description": "Genes are classified with a two-class Gaussian mixture model into two groups based on the mean Reads per UMI per gene. Enriched genes have higher mean Reads Per UMI, indicating they have been sequenced to much higher saturation.",
        "format": "integer",
    },
    {
        "name": "num_rpu_enriched_genes_off_target",
        "display_name": "Number of Enriched Non-Targeted Genes",
        "description": "Genes are classified with a two-class Gaussian mixture model into two groups based on the mean Reads per UMI per gene. Enriched genes have higher mean Reads Per UMI, indicating they have been sequenced to much higher saturation.",
        "format": "integer",
    },
    {
        "name": "mean_reads_per_umi_per_gene_on_target",
        "display_name": "Mean Reads per UMI per Targeted Gene",
        "description": "Mean number of reads per UMI per targeted gene, averaged across all genes.",
        "format": "%.2f",
    },
    {
        "name": "mean_reads_per_umi_per_gene_off_target",
        "display_name": "Mean Reads per UMI per Non-Targeted Gene",
        "description": "Mean number of reads per UMI per non-targeted gene, averaged across all genes.",
        "format": "%.2f",
    },
    {
        "name": "filtered_target_umi_reads_frac",
        "display_name": "Fraction of Reads Removed by UMI Filtering",
        "description": "Fraction of Reads confidently mapped to the Targeted Transcriptome and removed by targeted UMI filtering.",
        "format": "percent",
    },
    {
        "name": "filtered_target_umi_count_threshold",
        "display_name": "Reads per UMI threshold for UMI filtering",
        "description": "Reads per UMI threshold for UMI filtering. UMIs in targeted genes where the read support is (strictly) lower than this threshold are filtered out.",
        "format": "integer",
    },
]

# Note: these metrics only apply when batch correction is enabled
CHEMISTRY_BATCH_CORRECTION_METRICS = [
    {
        "name": "batch_effect_score_before_correction",
        "display_name": "Batch Effect Score before Correction",
        "description": "Before chemistry batch correction, the average fraction of 100 nearest neighbors belonging to the same batch, normalized by the expected value without batch effect.",
        "format": "%.2f",
    },
    {
        "name": "batch_effect_score_after_correction",
        "display_name": "Batch Effect Score after Correction",
        "description": "After chemistry batch correction, the average fraction of 100 nearest neighbors belonging to the same batch, normalized by the expected value without batch effect.",
        "format": "%.2f",
    },
]

METRICS = [
    {"name": "Summary", "metrics": SUMMARY_METRICS},
    {"name": "Sequencing", "metrics": SEQUENCING_METRICS},
    {"name": "Mapping", "metrics": MAPPING_METRICS},
    {"name": "Cells", "metrics": DETECTED_CELL_METRICS},
    {"name": "CRISPR Sequencing", "metrics": CRISPR_SEQUENCING_METRICS},
    {"name": "CRISPR Application", "metrics": CRISPR_APPLICATION_METRICS},
    {"name": "Multiplexing Sequencing", "metrics": MULTIPLEXING_SEQUENCING_METRICS},
    {"name": "Multiplexing Application", "metrics": MULTIPLEXING_APPLICATION_METRICS},
    {"name": "GEM Partitions", "metrics": GEM_METRICS},
    {"name": "Aggregation", "metrics": AGGREGATION_METRICS},
    {"name": "CRISPR Aggregation", "metrics": CRISPR_AGGREGATION_METRICS},
    {"name": "Antibody Aggregation", "metrics": ANTIBODY_AGGREGATION_METRICS},
    {"name": "Feature Aggregation", "metrics": CUSTOM_FEATURE_AGGREGATION_METRICS},
    {"name": "Antibody Cell Calling", "metrics": ANTIBODY_CELL_CALLING_METRICS},
    {"name": "Antibody Sequencing", "metrics": ANTIBODY_SEQUENCING_METRICS},
    {"name": "Antibody Application", "metrics": ANTIBODY_APPLICATION_METRICS},
    # Non- (CRISPR or Antibody) feature metrics
    {"name": "Feature Sequencing", "metrics": CUSTOM_FEATURE_SEQUENCING_METRICS},
    {"name": "Feature Application", "metrics": CUSTOM_FEATURE_APPLICATION_METRICS},
    {
        "name": "Chemistry Batch Correction",
        "metrics": CHEMISTRY_BATCH_CORRECTION_METRICS,
    },
    {"name": "Targeted Enrichment", "metrics": TARGETED_APPLICATION_METRICS},
]

# pylint: disable=line-too-long
CHARTS = [
    {
        "layout": {
            "title": "Barcode Rank",
            "width": 470,
            "height": 313,
            "margin": {"l": 60, "r": 0, "t": 30, "b": 40},
            "hovermode": "closest",
            "xaxis": {"title": "Barcodes", "type": "log"},
            "yaxis": {"title": "UMI counts", "type": "log"},
        },
        "data": [],
        "config": shared.CHARTS_PLOTLY_FIXED_CONFIG,
        "function": "plot_barcode_rank",
        "description": "The plot shows the count of filtered UMIs mapped to each barcode. Barcodes are not determined to be cell-associated strictly based on their UMI count. Instead, they could be determined based on their expression profile, or removed via Protein Aggregate Detection and Filtering and/or High Occupancy GEM Filtering. Therefore, some regions of the graph contain both cell-associated and background-associated barcodes. The color of the graph in these regions is based on the local density of barcodes that are cell-associated. Hovering over the plot displays the total number and percentage of barcodes in that region called as cells along with the number of UMI counts for those barcodes and barcode rank, ordered in descending order of UMI counts.",
        "name": "barcode_rank",
    },
    {
        "layout": {
            "title": "Cell UMI Counts",
            "width": 588,
            "height": 400,
            "margin": {"l": 65, "r": 0, "t": 30, "b": 70},
            "hovermode": "closest",
            "showlegend": True,
        },
        "data": [{"x": [], "y": [], "mode": "markers", "type": "scattergl"}],
        "config": shared.CHARTS_PLOTLY_MOVABLE_CONFIG,
        "description": "Each point represents a cell-barcode. The axes measure the total UMI counts in each cell-barcode that mapped to each transcriptome reference. The points are colored by the number of inferred cells in the GEM associated with each barcode. A multiplet represents either a GEM inferred to have encapsulated >1 cell or a barcode sequence that was shared by multiple single-cell GEMs.",
        "name": "barnyard_counts",
        "function": "plot_barnyard_barcode_counts",
    },
    {
        "layout": {
            "title": "t-SNE Projection of Cells Colored by UMI Counts",
            "width": 588,
            "height": 400,
            "margin": {"l": 70, "r": 0, "t": 30, "b": 70},
            "hovermode": "closest",
            "xaxis": {"title": "t-SNE1"},
            "yaxis": {"title": "t-SNE2"},
        },
        "data": [
            {
                "x": [],
                "y": [],
                "text": [],
                "hoverinfo": "text",
                "mode": "markers",
                "type": "scattergl",
                "marker": {
                    "color": [],
                    "size": 4,
                    "showscale": True,
                    "colorscale": "Jet",
                },
            }
        ],
        "config": shared.CHARTS_PLOTLY_MOVABLE_CONFIG,
        "name": "tsne_counts",
        "description": "Shown here are the total UMI counts for each cell-barcode. Cells with greater UMI counts likely have higher RNA content than cells with fewer UMI counts. The axes correspond to the 2-dimensional embedding produced by the t-SNE algorithm. In this space, pairs of cells that are close to each other have more similar gene expression profiles than cells that are distant from each other. The display is limited to a random subset of %d cells."
        % MAX_WEBSHIM_BCS_DIM,
        "function": "plot_tsne_totalcounts",
    },
    {
        "layout": {
            "title": "t-SNE projection of Cells Colored by Automated Clustering",
            "width": 588,
            "height": 400,
            "margin": {"l": 70, "r": 0, "t": 30, "b": 70},
            "hovermode": "closest",
            "xaxis": {"title": "t-SNE1"},
            "yaxis": {"title": "t-SNE2"},
        },
        "name": "tsne_clustering",
        "description": "These are the assignments of each cell-barcode a clusters by an automated clustering algorithm. The clustering groups together cells that have similar expression profiles. The axes correspond to the 2-dimensional embedding produced by the t-SNE algorithm. In this space, pairs of cells that are close to each other have more similar gene expression profiles than cells that are distant from each other. The display is limited to a random subset of %d cells and K-means up to K=%d. Please use Loupe Browser to view the entire dataset."
        % (MAX_WEBSHIM_BCS_DIM, MAX_WEBSHIM_KMEANS_K),
        "function": "plot_tsne",
    },
    {
        "table": {"width": 999, "height": 1000},
        "title": "Top Genes By Cluster (Log2 fold-change, p-value)",
        "name": "differential_expression",
        "description": 'The differential expression analysis seeks to find, for each cluster, genes that are more highly expressed in that cluster relative to the rest of the sample. Here a differential expression test was performed between each cluster and the rest of the sample for each gene. The Log2 fold-change (L2FC) is an estimate of the log2 ratio of expression in a cluster to that in all other cells. A value of 1.0 indicates 2-fold greater expression in the cluster of interest. The p-value is a measure of the statistical significance of the expression difference and is based on a negative binomial test. The p-value reported here has been adjusted for multiple testing via the Benjamini-Hochberg procedure. In this table you can click on a column to sort by that value. Also, in this table genes were filtered by (Mean UMI counts > %0.1f) and the top N genes by L2FC for each cluster were retained. Genes with L2FC < 0 or adjusted p-value >= %0.2f were grayed out. The number of top genes shown per cluster, N, is set to limit the number of table entries shown to %d; N=%d/K^2 where K is the number of clusters. N can range from 1 to %d. For the full table, please refer to the "differential_expression.csv" files produced by the pipeline.'
        % (
            TOP_DE_GENES_MIN_MEAN,
            PVALUE_DEEMPHASIS_CUTOFF,
            MAX_DE_TABLE_ENTRIES,
            MAX_DE_TABLE_ENTRIES,
            MAX_TOP_N_GENES,
        ),
        "function": "plot_differential_expression",
    },
    {
        "layout": {
            "title": "Sequencing Saturation",
            "width": 588,
            "height": 400,
            "margin": {"l": 70, "r": 65, "t": 30, "b": 70},
            "showlegend": False,
            "xaxis": {"title": "Mean Reads per Cell"},
            "yaxis": {"title": "Sequencing Saturation", "range": [0, 1]},
            "shapes": [
                {
                    "type": "line",
                    "x0": 0,
                    "y0": SATURATION_LINE,
                    "x1": 0,
                    "y1": SATURATION_LINE,
                    "line": {"color": "rgb(128, 128, 128)", "width": 4, "dash": "dot"},
                }
            ],
        },
        "data": [],  # data entries are built in the function
        "kwargs": {
            "metric_suffix": "subsampled_duplication_frac",
            "show_multi_genome_only": True,
        },
        "config": shared.CHARTS_PLOTLY_MOVABLE_CONFIG,
        "name": "sequencing_saturation",
        "description": "This plot shows the {} metric as a function of downsampled sequencing depth in mean reads per cell, up to the observed sequencing depth. {} is a measure of the observed library complexity, and approaches 1.0 (100%) when all converted mRNA transcripts have been sequenced. The slope of the curve near the endpoint can be interpreted as an upper bound to the benefit to be gained from increasing the sequencing depth beyond this point. The dotted line is drawn at a value reasonably approximating the saturation point.".format(
            SEQUENCING_SATURATION_METRIC["display_name"],
            SEQUENCING_SATURATION_METRIC["display_name"],
        ),
        "function": "plot_subsampled_scatterplot_metric",
    },
    {
        "layout": {
            "title": "Median Genes per Cell",
            "width": 588,
            "height": 400,
            "margin": {"l": 70, "r": 0, "t": 30, "b": 70},
            "showlegend": True,
            "xaxis": {"title": "Mean Reads per Cell"},
            "yaxis": {"title": "Median Genes per Cell"},
        },
        "data": [],  # data entries are built in the function
        "kwargs": {"metric_suffix": "subsampled_filtered_bcs_median_unique_genes_detected"},
        "kwargs_prefixes": ["references"],
        "config": shared.CHARTS_PLOTLY_MOVABLE_CONFIG,
        "name": "median_genes_per_cell",
        "description": "This plot shows the {} as a function of downsampled sequencing depth in mean reads per cell, up to the observed sequencing depth. The slope of the curve near the endpoint can be interpreted as an upper bound to the benefit to be gained from increasing the sequencing depth beyond this point.".format(
            GENES_PER_DETECTED_CELL_METRIC["display_name"]
        ),
        "function": "plot_subsampled_scatterplot_metric",
    },
]

METRIC_ALARMS = [
    {
        "name": MULTI_NUMBER_OF_DETECTED_CELLS_METRIC["name"],
        "format": MULTI_NUMBER_OF_DETECTED_CELLS_METRIC["format"],
        "error": {
            "title": "No Cells Detected",
            "message": "No valid sequencing data was detected. Please check the sequencing data.",
            "test": "< 1",
        },
        "warn": {
            "title": "Low Number of Cells Detected",
            "message": "Ideal > 100. This usually indicates poor cell handling, poor library quality, or poor sequencing quality. Application performance is likely to be affected.",
            "test": "< 100",
        },
    },
    {
        "name": MULTI_NUMBER_OF_ANTIBODY_DETECTED_CELLS_METRIC["name"],
        "format": MULTI_NUMBER_OF_ANTIBODY_DETECTED_CELLS_METRIC["format"],
        "error": {
            "title": "No Cells Detected",
            "message": "No valid sequencing data was detected. Please check the sequencing data.",
            "test": "< 1",
        },
        "warn": {
            "title": "Low Number of Cells Detected",
            "message": "Ideal > 100. This usually indicates poor cell handling, poor library quality, or poor sequencing quality. Might indicate a need for a larger antibody panel. Application performance is likely to be affected.",
            "test": "< 100",
        },
    },
    {
        "name": GOOD_BCS_METRIC["name"],
        "format": GOOD_BCS_METRIC["format"],
        "error": {
            "title": "Low Fraction Valid Barcodes",
            "message": "Ideal > 75%. This may indicate a quality issue with the R1 read. Application performance is likely to be affected.",
            "test": "< 0.5",
        },
        "warn": {
            "title": "Low Fraction Valid Barcodes",
            "message": "Ideal > 75%. This may indicate a quality issue with the R1 read. Application performance may be affected.",
            "test": "< 0.75",
        },
    },
    {
        "name": GOOD_UMIS_METRIC["name"],
        "format": GOOD_UMIS_METRIC["format"],
        "error": {
            "title": "Low Fraction Valid UMIs",
            "message": "Ideal > 75%. This may indicate a quality issue with the R1 read. Application performance is likely to be affected.",
            "test": "< 0.5",
        },
        "warn": {
            "title": "Low Fraction Valid UMIs",
            "message": "Ideal > 75%. This may indicate a quality issue with the R1 read. Application performance may be affected.",
            "test": "< 0.75",
        },
    },
    {
        "name": MULTI_TRANSCRIPTOME_CONF_MAPPED_READS_METRIC["name"],
        "format": MULTI_TRANSCRIPTOME_CONF_MAPPED_READS_METRIC["format"],
        "error": {
            "title": "Low Fraction Reads Confidently Mapped To Transcriptome",
            "message": "Ideal > 30%. This can indicate use of the wrong reference transcriptome, a reference transcriptome with overlapping genes, poor library quality, poor sequencing quality, or reads shorter than the recommended minimum. Application performance is likely to be affected.",
            "test": "< 0.20",
        },
        "warn": {
            "title": "Low Fraction Reads Confidently Mapped To Transcriptome",
            "message": "Ideal > 30%. This can indicate use of the wrong reference transcriptome, a reference transcriptome with overlapping genes, poor library quality, poor sequencing quality, or reads shorter than the recommended minimum. Application performance may be affected.",
            "test": "< 0.30",
        },
    },
    {
        "name": ANTISENSE_CONF_MAPPED_READS_METRIC["name"],
        "format": ANTISENSE_CONF_MAPPED_READS_METRIC["format"],
        "error": {
            "title": "High Fraction of Reads Mapped Antisense to Genes",
            "message": "Ideal < 10%. This can indicate use of an unsupported chemistry type (e.g. using Single Cell V(D)J for gene counting). Application performance is likely to be affected.",
            "test": "> 0.30",
        },
        "warn": {
            "title": "High Fraction of Reads Mapped Antisense to Genes",
            "message": "Ideal < 10%. This can indicate use of an unsupported chemistry type (e.g. using Single Cell V(D)J for gene counting). Application performance may be affected.",
            "test": "> 0.10",
        },
    },
    {
        "name": "bc_bases_with_q30_frac",
        "format": "percent",
        "error": {
            "title": "Low Barcode Q30 Fraction (Illumina I7 Read for Single Cell 3' v1, R1 for Single Cell 3' v2/v3 and Single Cell 5')",
            "message": "Ideal > 55%. Application performance is likely to be affected.",
            "test": "< 0.45",
        },
        "warn": {
            "title": "Low Barcode Q30 Fraction (Illumina I7 Read for Single Cell 3' v1, R1 for Single Cell 3' v2/v3 and Single Cell 5')",
            "message": "Ideal > 55%. Application performance may be affected.",
            "test": "< 0.55",
        },
    },
    {
        "name": "read_bases_with_q30_frac",
        "format": "percent",
        "error": {
            "title": "Low RNA Read Q30 Fraction (Illumina R1 for Single Cell 3' v1 and Single Cell 5' paired end, R2 for Single Cell 3' v2/v3 and Single Cell 5' R2-only)",
            "message": "Ideal > 65%. Application performance is likely to be affected.",
            "test": "< 0.55",
        },
        "warn": {
            "title": "Low RNA Read Q30 Fraction (Illumina R1 for Single Cell 3' v1 and Single Cell 5' paired end, R2 for Single Cell 3' v2/v3 and Single Cell 5' R2-only)",
            "message": "Ideal > 65%. Application performance may be affected.",
            "test": "< 0.65",
        },
    },
    {
        "name": "umi_bases_with_q30_frac",
        "format": "percent",
        "error": {
            "title": "Low UMI Q30 Fraction (Illumina R2 Read for Single Cell 3' v1, R1 for Single Cell 3' v2/v3 and Single Cell 5')",
            "message": "Ideal > 75%. Application performance is likely to be affected.",
            "test": "< 0.65",
        },
        "warn": {
            "title": "Low UMI Q30 Fraction (Illumina R2 Read for Single Cell 3' v1, R1 for Single Cell 3' v2/v3 and Single Cell 5')",
            "message": "Ideal > 75%. Application performance may be affected.",
            "test": "< 0.75",
        },
    },
    {
        "name": MULTI_FRAC_READS_IN_CELLS_METRIC["name"],
        "format": MULTI_FRAC_READS_IN_CELLS_METRIC["format"],
        "error": {
            "title": "Low Fraction Reads in Cells",
            "message": "Ideal > 70%. Application performance is likely to be affected. Many of the reads were not assigned to cell-associated barcodes. This could be caused by high levels of ambient RNA or by a significant population of cells with a low RNA content, which the algorithm did not call as cells. The latter case can be addressed by inspecting the data to determine the appropriate cell count and using --force-cells.",
            "test": "< 0.50",
        },
        "warn": {
            "title": "Low Fraction Reads in Cells",
            "message": "Ideal > 70%. Application performance may be affected. Many of the reads were not assigned to cell-associated barcodes. This could be caused by high levels of ambient RNA or by a significant population of cells with a low RNA content, which the algorithm did not call as cells. The latter case can be addressed by inspecting the data to determine the appropriate cell count and using --force-cells.",
            "test": "< 0.70",
        },
    },
    # Batch Aggregation-specific alerts
    {
        "name": "lowest_frac_reads_kept",
        "format": "percent",
        "error": {
            "title": "Low Post-Normalization Read Depth",
            "message": "Ideal > 50%. There may be large differences in sequencing depth across the input libraries. Application performance is likely to be affected.",
            "test": "< 0.25",
        },
        "warn": {
            "title": "Low Post-Normalization Read Depth",
            "message": "Ideal > 50%. There may be large differences in sequencing depth across the input libraries. Application performance may be affected.",
            "test": "< 0.50",
        },
    },
    # Antibody aggregation-specific alerts
    {
        "name": REPORT_PREFIX_ANTIBODY + "reads_lost_to_aggregate_GEMs",
        "format": "percent",
        "error": {
            "title": "High fraction of Antibody library reads coming from aggregate barcodes",
            "message": "Ideal < 5%. Detected some aggregate barcodes. These barcodes and their corresponding reads have been removed from the final matrix.",
            "test": "> 0.5",
        },
        "warn": {
            "title": "High fraction of Antibody library reads coming from aggregate barcodes",
            "message": "Ideal < 5%. Detected some aggregate barcodes. These barcodes and their corresponding reads have been removed from the final matrix.",
            "test": "> 0.05",
        },
    },
    # Targeting-specific alerts
    {
        "name": "multi_frac_conf_transcriptomic_reads_on_target",
        "format": "percent",
        "error": {
            "title": "Low Fraction Reads Confidently Mapped To Targeted Transcriptome",
            "message": "Ideal > 30%. This can indicate use of a small panel with low aggregate expression, use of the wrong target panel, use of a custom panel with high background, or inefficient targeting. Targeted performance may be affected.",
            "test": "> 0.20",
        },
        "warn": {
            "title": "Low Fraction Reads Confidently Mapped To Targeted Transcriptome",
            "message": "Ideal > 30%. This can indicate use of a small panel with low aggregate expression, use of the wrong target panel, use of a custom panel with high background, or inefficient targeting. Targeted performance may be affected.",
            "test": "> 0.30",
        },
    },
    {
        "name": "frac_on_target_genes_enriched",
        "format": "percent",
        "error": {
            "title": "Low Fraction of Targeted Genes Enriched",
            "message": "Ideal > 70%. This can indicate use of the wrong target panel, or inefficient targeting. Targeted performance may be affected.",
            "test": "> 0.70",
        },
        "warn": {
            "title": "Low Fraction of Targeted Genes Enriched",
            "message": "Ideal > 70%. This can indicate use of the wrong target panel, or inefficient targeting. Targeted performance may be affected.",
            "test": "> 0.50",
        },
    },
    {
        "name": "frac_off_target_genes_enriched",
        "format": "percent",
        "error": {
            "title": "High Fraction of Non-Targeted Genes Enriched",
            "message": "This can indicate use of the wrong target panel, or inefficient targeting. Targeted performance may be affected.",
            "test": "< 0.30",
        },
        "warn": {
            "title": "High Fraction of Non-Targeted Genes Enriched",
            "message": "This can indicate use of the wrong target panel, or inefficient targeting. Targeted performance may be affected.",
            "test": "< 0.50",
        },
    },
]

# Targeting color scheme for websummary
ON_TARGET_COLOR = "#0071D9"
OFF_TARGET_COLOR = "darkgrey"

# Targeting column name constants in feature_metrics table
TARGETING_COLNAME = "is_targeted"
READS_COLNAME = "num_reads"
UMIS_COLNAME = "num_umis"
READS_IN_CELLS_COLNAME = "num_reads_cells"
UMIS_IN_CELLS_COLNAME = "num_umis_cells"
LOG_RPU_CELLS_COLNAME = "mean_reads_per_umi_cells_log10"
RPU_COLNAME = "mean_reads_per_umi"
ENRICHMENT_COLNAME = "enriched_rpu"
FEATURE_ID_COLNAME = "feature_id"
FEATURE_NAME_COLNAME = "feature_name"
TARGET_FEATURE_METRICS_COLS = [
    FEATURE_ID_COLNAME,
    FEATURE_NAME_COLNAME,
    TARGETING_COLNAME,
    READS_IN_CELLS_COLNAME,
    UMIS_IN_CELLS_COLNAME,
    LOG_RPU_CELLS_COLNAME,
    RPU_COLNAME,
    ENRICHMENT_COLNAME,
]
