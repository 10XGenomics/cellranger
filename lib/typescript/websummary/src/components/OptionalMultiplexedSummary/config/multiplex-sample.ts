/* eslint-disable sort-keys */
import {
  confidentlyMappedReadsInCellsHelp,
  gdnaMetricsTable,
} from "./multiplex-library";

const sensitivityHelp =
  "Sensitivity of V(D)J detection measured by UMIs for each contig per cell.";

const sensitivityEnrichmentHelp =
  "Sensitivity of V(D)J detection measured by UMIs for each contig per cell and efficiency of V(D)J enrichment measured by the fraction of reads mapped to V(D)J genes.";

const MultiplexSampleConfig = {
  GEX: {
    hero: {
      title: "Cell Calling Quality",
      help: "Summary statistics about cell-associated barcodes.",
      entries: [
        {
          key: "genome",
          help: "Genome used for this analysis.",
        },
        {
          key: "total_singlets",
          help: "Number of cells called from this sample.",
        },
        {
          key: "confidently_mapped_reads_in_cells",
          help: confidentlyMappedReadsInCellsHelp,
        },
        {
          key: "median_genes_per_singlet",
          help: "The median number of genes detected per cell called from this sample. Detection is defined as the presence of at least 1 UMI count.",
        },
        {
          key: "median_umi_per_singlet",
          help: "Median number of UMIs obtained from the cells called from this sample.",
        },
        {
          key: "total_genes_detected",
          help: "The number of genes with at least one UMI count in the cells in this sample.",
        },
        {
          key: "number_of_reads_in_cells",
          help: "The total number of reads from cells called from this sample.",
        },
      ],
    },
    mapping: {
      title: "Mapping Quality (Amongst Reads From Cells Assigned To Sample)",
      help: "Mapping metrics calculated based on reads from cells assigned to this sample.",
      entries: [
        {
          key: "confidently_mapped_to_transcriptome",
          help: "Fraction of reads that mapped to a unique gene in the transcriptome. The read must be consistent with annotated splice junctions. These reads are considered for UMI counting.",
        },
        {
          key: "mapped_to_genome",
          help: "Fraction of reads that mapped to the genome.",
        },
        {
          key: "confidently_mapped_to_genome",
          help: "Fraction of reads that mapped uniquely to the genome. If a gene mapped to exonic loci from a single gene and also to non-exonic loci, it is considered uniquely mapped to one of the exonic loci.",
        },
        {
          key: "confidently_mapped_to_exonic_regions",
          help: "Fraction of reads that mapped uniquely to an exonic region of the genome.",
        },
        {
          key: "confidently_mapped_to_intronic_regions",
          help: "Fraction of reads that mapped uniquely to an intronic region of the genome.",
        },
        {
          key: "confidently_mapped_to_intergenic_regions",
          help: "Fraction of reads that mapped uniquely to an intergenic region of the genome.",
        },
        {
          key: "confidently_mapped_antisense",
          help: "Fraction of reads confidently mapped to the transcriptome, but on the opposite strand of their annotated gene. A read is counted as antisense if it has any alignments that are consistent with an exon of a transcript but antisense to it, and has no sense alignments.",
        },
        {
          key: "reads_mapped_to_probe_set",
          help: "Fraction of reads that mapped to the probe set.",
        },
        {
          key: "reads_confidently_mapped_to_probe_set",
          help: "Fraction of reads that mapped uniquely to a probe in the probe set.",
        },
        {
          key: "reads_confidently_mapped_to_filtered_probe_set",
          help: "Fraction of reads from probes that map to a unique gene. These reads are considered for UMI counting. For more information on probe filtering please visit https://www.10xgenomics.com/support",
        },
        {
          key: "reads_half_mapped_to_probe_set",
          help: "Fraction of reads that mapped to unpaired ligation products.",
        },
        {
          key: "reads_split_mapped_to_probe_set",
          help: "Fraction of reads that mapped to mispaired ligation products.",
        },
      ],
    },
    gdna: gdnaMetricsTable,
  },
  VDJ_T: {
    hero: {
      title: "Cell Calling Quality",
      help: "Summary statistics about VDJ cells called in this sample.",
      entries: [
        {
          key: "vdj_filtered_bcs",
          help: "The number of barcodes estimated to be associated with T cells.",
        },
        {
          key: "multi_vdj_assembly_contig_pair_productive_full_len_bc_frac",
          help: "Number of cell barcodes for which at least 1 full-length productive sequence was found for each chain of the (TRA, TRB) receptor pair.",
        },
        {
          key: "TRA_TRB_vdj_assembly_contig_pair_productive_full_len_bc_frac",
          help: "Fraction of cell-associated barcodes with at least one productive contig for each chain of the (TRA, TRB) receptor pair.",
        },
        {
          key: "TRA_vdj_assembly_prod_cdr_bc_frac",
          help: "Fraction of cell-associated barcodes with at least one contig that spans the 5' end of the V region to the 3' end of the J region for TRA, has a start codon in the expected part of the V sequence, has an in-frame CDR3, and has no stop codons in the aligned V-J region.",
        },
        {
          key: "TRB_vdj_assembly_prod_cdr_bc_frac",
          help: "Fraction of cell-associated barcodes with at least one contig that spans the 5' end of the V region to the 3' end of the J region for TRB, has a start codon in the expected part of the V sequence, has an in-frame CDR3, and has no stop codons in the aligned V-J region.",
        },
        {
          key: "vdj_filtered_bcs_cum_frac",
          help: "Number of reads with cell-associated barcodes divided by the number of reads with valid barcodes.",
        },
      ],
    },
    sensitivity: {
      title: "VDJ Sensitivity",
      help: sensitivityHelp,
      entries: [
        {
          key: "TRA_vdj_assembly_umis_per_cell_median",
          help: "Median number of UMIs assigned to a TRA contig per cell.",
        },
        {
          key: "TRB_vdj_assembly_umis_per_cell_median",
          help: "Median number of UMIs assigned to a TRB contig per cell.",
        },
      ],
    },
    sensitivityEnrichment: {
      title: "VDJ Sensitivity and Enrichment",
      help: sensitivityEnrichmentHelp,
      entries: [
        {
          key: "TRA_vdj_assembly_umis_per_cell_median",
          help: "Median number of UMIs assigned to a TRA contig per cell.",
        },
        {
          key: "TRB_vdj_assembly_umis_per_cell_median",
          help: "Median number of UMIs assigned to a TRB contig per cell.",
        },
        {
          key: "multi_vdj_recombinome_mapped_reads_frac",
          help: "Fraction of reads with valid barcodes that partially or wholly map to any germline V(D)J gene segment.",
        },
        {
          key: "TRA_vdj_recombinome_mapped_reads_frac",
          help: "Fraction of reads with valid barcodes that map partially or wholly to a germline TRA gene segment.",
        },
        {
          key: "TRB_vdj_recombinome_mapped_reads_frac",
          help: "Fraction of reads with valid barcodes that map partially or wholly to a germline TRB gene segment.",
        },
      ],
    },
    pairedClonotypeDiversity: {
      title: "Clonotype Metrics",
      entries: [
        {
          key: "multi_raw_vdj_paired_clonotype_diversity",
          help: "Effective diversity of the paired clonotypes, computed as the Inverse Simpson Index of the clonotype frequencies. A value of 1 indicates a minimally diverse sample - only one distinct clonotype was detected. A value equal to the estimated number of cells indicates a maximally diverse sample.",
        },
      ],
    },
  },
  VDJ_T_GD: {
    hero: {
      title: "Cell Calling Quality",
      help: "Summary statistics about VDJ cells called in this sample.",
      entries: [
        {
          key: "vdj_filtered_bcs",
          help: "The number of barcodes estimated to be associated with T cells.",
        },
        {
          key: "multi_vdj_assembly_contig_pair_productive_full_len_bc_frac",
          help: "Fraction of cell-associated barcodes with at least one productive contig for each chain of the receptor pair. A productive contig satisfies the following conditions: the contig annotations span the 5' end of the V region to the 3' end of the J region of the chain, a start codon was found in the expected part of the V sequence, an in-frame CDR3 amino acid motif was found, and no stop codons were found in the aligned V-J region.",
        },
        {
          key: "TRG_TRD_vdj_assembly_contig_pair_productive_full_len_bc_frac",
          help: "Fraction of cell-associated barcodes with at least one productive contig for each chain of the (TRG, TRD) receptor pair.",
        },
        {
          key: "TRG_vdj_assembly_prod_cdr_bc_frac",
          help: "Fraction of cell-associated barcodes with at least one contig that spans the 5' end of the V region to the 3' end of the J region for TRG, has a start codon in the expected part of the V sequence, has an in-frame CDR3, and has no stop codons in the aligned V-J region.",
        },
        {
          key: "TRD_vdj_assembly_prod_cdr_bc_frac",
          help: "Fraction of cell-associated barcodes with at least one contig that spans the 5' end of the V region to the 3' end of the J region for TRD, has a start codon in the expected part of the V sequence, has an in-frame CDR3, and has no stop codons in the aligned V-J region.",
        },
        {
          key: "vdj_filtered_bcs_cum_frac",
          help: "Number of reads with cell-associated barcodes divided by the number of reads with valid barcodes.",
        },
      ],
    },
    sensitivity: {
      title: "VDJ Sensitivity",
      help: sensitivityHelp,
      entries: [
        {
          key: "TRG_vdj_assembly_umis_per_cell_median",
          help: "Median number of UMIs assigned to a TRG contig per cell.",
        },
        {
          key: "TRD_vdj_assembly_umis_per_cell_median",
          help: "Median number of UMIs assigned to a TRD contig per cell.",
        },
      ],
    },
    sensitivityEnrichment: {
      title: "VDJ Sensitivity and Enrichment",
      help: sensitivityEnrichmentHelp,
      entries: [
        {
          key: "TRG_vdj_assembly_umis_per_cell_median",
          help: "Median number of UMIs assigned to a TRG contig per cell.",
        },
        {
          key: "TRD_vdj_assembly_umis_per_cell_median",
          help: "Median number of UMIs assigned to a TRD contig per cell.",
        },
        {
          key: "multi_vdj_recombinome_mapped_reads_frac",
          help: "Fraction of reads with valid barcodes that partially or wholly map to any germline V(D)J gene segment.",
        },
        {
          key: "TRG_vdj_recombinome_mapped_reads_frac",
          help: "Fraction of reads with valid barcodes that map partially or wholly to a germline TRG gene segment.",
        },
        {
          key: "TRD_vdj_recombinome_mapped_reads_frac",
          help: "Fraction of reads with valid barcodes that map partially or wholly to a germline TRD gene segment.",
        },
      ],
    },
    pairedClonotypeDiversity: {
      title: "Clonotype Metrics",
      entries: [
        {
          key: "multi_raw_vdj_paired_clonotype_diversity",
          help: "Effective diversity of the paired clonotypes, computed as the Inverse Simpson Index of the clonotype frequencies. A value of 1 indicates a minimally diverse sample - only one distinct clonotype was detected. A value equal to the estimated number of cells indicates a maximally diverse sample.",
        },
      ],
    },
  },
  VDJ_B: {
    hero: {
      title: "Cell Calling Quality",
      help: "Summary statistics about VDJ cells called in this sample.",
      entries: [
        {
          key: "vdj_filtered_bcs",
          help: "The number of barcodes estimated to be associated with B cells.",
        },
        {
          key: "multi_vdj_assembly_contig_pair_productive_full_len_bc_frac",
          help: "Fraction of cell-associated barcodes with at least one productive contig for each chain of the receptor pair. A productive contig satisfies the following conditions: the contig annotations span the 5' end of the V region to the 3' end of the J region of the chain, a start codon was found in the expected part of the V sequence, an in-frame CDR3 amino acid motif was found, and no stop codons were found in the aligned V-J region.",
        },
        {
          key: "IGK_IGH_vdj_assembly_contig_pair_productive_full_len_bc_frac",
          help: "Fraction of cell-associated barcodes with at least one productive contig for each chain of the (IGK, IGH) receptor pair.",
        },
        {
          key: "IGL_IGH_vdj_assembly_contig_pair_productive_full_len_bc_frac",
          help: "Fraction of cell-associated barcodes with at least one productive contig for each chain of the (IGL, IGH) receptor pair.",
        },
        {
          key: "IGH_vdj_assembly_prod_cdr_bc_frac",
          help: "Fraction of cell-associated barcodes with at least one contig that spans the 5' end of the V region to the 3' end of the J region for IGH, has a start codon in the expected part of the V sequence, has an in-frame CDR3, and has no stop codons in the aligned V-J region.",
        },
        {
          key: "IGK_vdj_assembly_prod_cdr_bc_frac",
          help: "Fraction of cell-associated barcodes with at least one contig that spans the 5' end of the V region to the 3' end of the J region for IGK, has a start codon in the expected part of the V sequence, has an in-frame CDR3, and has no stop codons in the aligned V-J region.",
        },
        {
          key: "IGL_vdj_assembly_prod_cdr_bc_frac",
          help: "Fraction of cell-associated barcodes with at least one contig that spans the 5' end of the V region to the 3' end of the J region for IGL, has a start codon in the expected part of the V sequence, has an in-frame CDR3, and has no stop codons in the aligned V-J region.",
        },
        {
          key: "vdj_filtered_bcs_cum_frac",
          help: "Number of reads with cell-associated barcodes divided by the number of reads with valid barcodes.",
        },
      ],
    },
    sensitivity: {
      title: "VDJ Sensitivity",
      help: sensitivityHelp,
      entries: [
        {
          key: "IGH_vdj_assembly_umis_per_cell_median",
          help: "Median number of UMIs assigned to a IGH contig per cell.",
        },
        {
          key: "IGK_vdj_assembly_umis_per_cell_median",
          help: "Median number of UMIs assigned to a IGK contig per cell.",
        },
        {
          key: "IGL_vdj_assembly_umis_per_cell_median",
          help: "Median number of UMIs assigned to a IGL contig per cell.",
        },
      ],
    },
    sensitivityEnrichment: {
      title: "VDJ Sensitivity and Enrichment",
      help: sensitivityEnrichmentHelp,
      entries: [
        {
          key: "IGH_vdj_assembly_umis_per_cell_median",
          help: "Median number of UMIs assigned to a IGH contig per cell.",
        },
        {
          key: "IGK_vdj_assembly_umis_per_cell_median",
          help: "Median number of UMIs assigned to a IGK contig per cell.",
        },
        {
          key: "IGL_vdj_assembly_umis_per_cell_median",
          help: "Median number of UMIs assigned to a IGL contig per cell.",
        },
        {
          key: "multi_vdj_recombinome_mapped_reads_frac",
          help: "Fraction of reads with valid barcodes that partially or wholly map to any germline V(D)J gene segment.",
        },
        {
          key: "IGH_vdj_recombinome_mapped_reads_frac",
          help: "Fraction of reads with valid barcodes that map partially or wholly to a germline IGH gene segment.",
        },
        {
          key: "IGK_vdj_recombinome_mapped_reads_frac",
          help: "Fraction of reads with valid barcodes that map partially or wholly to a germline IGK gene segment.",
        },
        {
          key: "IGL_vdj_recombinome_mapped_reads_frac",
          help: "Fraction of reads with valid barcodes that map partially or wholly to a germline IGL gene segment.",
        },
      ],
    },
    pairedClonotypeDiversity: {
      title: "Clonotype Metrics",
      entries: [
        {
          key: "multi_raw_vdj_paired_clonotype_diversity",
          help: "Effective diversity of the paired clonotypes, computed as the Inverse Simpson Index of the clonotype frequencies. A value of 1 indicates a minimally diverse sample - only one distinct clonotype was detected. A value equal to the estimated number of cells indicates a maximally diverse sample.",
        },
      ],
    },
  },
  AB: {
    hero: {
      title: "Antibody Expression",
      help: "Summary statistics about the detection and association with cells for antibody reads in this sample.",
      entries: [
        {
          key: "total_singlets",
          help: "Number of cells called from this sample. Cell calling is based on gene expression data when present.",
        },
        {
          key: "fraction_antibody_reads",
          help: "Fraction of read pairs that contain a recognized antibody Feature Barcode.",
        },
        {
          key: "reads_in_cells",
          help: "The fraction of valid-barcode, valid-UMI, recognized antibody Feature Barcode reads with cell-associated barcodes.",
        },
        {
          key: "median_umis_per_singlet",
          help: "Median number of UMIs obtained from cells called from this sample.",
        },
        {
          key: "antibody_reads_usable_per_cell",
          help: "Mean number of usable reads (valid UMI, recognized antibody Feature Barcode) sequenced from cells called from this sample.",
        },
        {
          key: "fraction_reads_in_aggregate_barcodes",
          help: "Fraction of read pairs with valid barcodes that were removed because they are aggregates out of all reads with valid barcodes that are assigned to this sample (not just reads from cells).",
        },
        {
          key: "number_of_reads_in_cells",
          help: "The total number of reads from cells associated with this sample.",
        },
      ],
    },
  },
  CRISPR: {
    hero: {
      title: "Guide Expression",
      help: "Summary statistics about the detection and association with cells for CRISPR reads in this sample.",
      entries: [
        {
          key: "total_singlets",
          help: "Number of cells called from this sample.",
        },
        {
          key: "reads_in_cells",
          help: "The fraction of valid-barcode, valid-UMI, recognized guide Feature Barcode reads with cell-associated barcodes.",
        },
        {
          key: "median_umis_per_singlet",
          help: "Median number of UMIs obtained from the cells called from this sample.",
        },
        {
          key: "cells_with_one_or_more_protospacers_detected",
          help: "Fraction of cells with one or more protospacers detected. In the multiplexing case, only cell-associated barcodes assigned exactly one CMO are included in this calculation.",
        },
        {
          key: "cells_with_two_or_more_protospacers_detected",
          help: "Fraction of cells with two or more protospacers detected. In the multiplexing case, only cell-associated barcodes assigned exactly one CMO are included in this calculation.",
        },
        {
          key: "fraction_reads_with_putative_protospacer",
          help: "Fraction of CRISPR library reads from which a putative protospacer sequence could be extracted.",
        },
        {
          key: "fraction_guide_reads",
          help: "Fraction of CRISPR library reads with a recognized protospacer sequence.",
        },
        {
          key: "guide_reads_usable_per_cell",
          help: "Mean number of usable reads (valid UMI, recognized protospacer sequence) sequenced from the cells called from this sample.",
        },
        {
          key: "number_of_reads_in_cells",
          help: "The total number of reads from cells associated with this sample.",
        },
      ],
    },
  },
  AG: {
    hero: {
      title: "Antigen Expression",
      help: "Summary statistics about the detection and association with cells for antigen reads in this sample.",
      entries: [
        {
          key: "feature_type",
          help: "The feature type used for computing the metrics.",
        },
        {
          key: "total_singlets",
          help: "Number of cells called from this sample from the respective feature type (gene expression or VDJ).",
        },
        {
          key: "median_umis_per_singlet",
          help: "Median number of antigen UMIs obtained from cells called from this sample.",
        },
        {
          key: "antigen_reads_usable_per_cell",
          help: "Mean number of usable reads (valid UMI, recognized antigen-barcode) sequenced from cells called from this sample.",
        },
      ],
    },
  },
  Custom: {
    hero: {
      title: "Feature Expression",
      help: "Summary statistics about the detection and association with cells for custom feature reads in this sample.",
      entries: [
        {
          key: "total_singlets",
          help: "Number of cells called from this sample.",
        },
        {
          key: "median_umis_per_singlet",
          help: "Median number of UMIs obtained from the cells called from this sample.",
        },
        {
          key: "feature_reads_usable_per_cell",
          help: "Mean number of usable reads (valid UMI, recognized feature-barcode sequence) sequenced from the cells called from this sample.",
        },
      ],
    },
  },
};
// TODO
// Add `as const`, which then enables us to do the following:
// type ExtractAllKeys<T> =
//   T extends { key: infer K extends string } ? K :
//   T extends (infer U)[] ? ExtractAllKeys<U> :
//   T extends object ? ExtractAllKeys<T[keyof T]> :
//   never;
//
// ExtractAllKeys<typeof MultiplexSampleConfig.GEX.hero.entries>
// | ExtractAllKeys<typeof MultiplexSampleConfig.GEX.mapping.entries>
//
// To give all valid keys per config, which we can then use to
// get compile time checks on how we use key literals.

export default MultiplexSampleConfig;
