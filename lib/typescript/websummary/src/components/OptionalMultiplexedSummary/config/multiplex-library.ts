/* eslint-disable sort-keys */
import { MetricEntryConfig, MetricTableConfig } from "./types";

export const q30BarcodeHelp = [
  {
    msg: "Fraction of cell barcode bases with Q-score >= 30, excluding very low quality/no-call (Q <= 2) bases from the denominator. The cell barcode is the combination of the GEM barcode and probe barcode.",
    isRtl: true,
  },
  {
    msg: "Fraction of cell barcode bases with Q-score >= 30, excluding very low quality/no-call (Q <= 2) bases from the denominator.",
    isRtl: false,
  },
];

export const q30RnaRead1Help = [
  {
    msg: "Fraction of probe read bases (R2) with Q-score >= 30, excluding very low quality/no-call (Q <= 2) bases from the denominator.",
    isRtl: true,
  },
  {
    msg: "Fraction of RNA Read bases with Q-score >= 30, excluding very low quality/no-call (Q <= 2) bases from the denominator.",
    isRtl: false,
  },
];

export const validBarcodesHelp = [
  {
    msg: "Fraction of reads with barcodes that are present in the whitelist after barcode correction. The cell barcode is the combination of the GEM barcode and probe barcode and both parts must be valid following correction.",
    isRtl: true,
  },
  {
    msg: "Fraction of reads with barcodes that are present in the whitelist after barcode correction.",
    isRtl: false,
  },
];

export const cellCallingQualityHelp =
  "Summary statistics about cell-associated barcodes.";

const featureBarcodingLibraryQualityTableEntries: MetricEntryConfig[] = [
  {
    key: "mean_reads_per_cell_associated_partition",
    help: "The total number of sequenced read pairs divided by the number of cell-associated barcodes.",
  },
  {
    key: "sequencing_saturation",
    help: "The fraction of reads originating from an already-observed UMI. This fraction is a function of library complexity and sequencing depth. It is calculated as the ratio of reads with non-unique (cell barcode, UMI, Feature Barcode) combinations to the total number of valid reads (those with recognized cell barcode, UMI, and Feature Barcode).",
  },
  {
    key: "valid_barcodes",
    help: validBarcodesHelp,
  },
  {
    key: "valid_gem_barcodes",
    help: "Fraction of reads with GEM barcodes that are present in the whitelist after barcode correction.",
  },
  {
    key: "valid_probe_barcodes",
    help: "Fraction of reads with probe barcodes that are present in the whitelist after barcode correction.",
  },
  {
    key: "valid_umis",
    help: "Fraction of reads with valid UMI sequences; i.e. UMI sequences that do not contain Ns and that are not homopolymers.",
  },
];

const featureBarcodingSequencingMetricsTableEntries: MetricEntryConfig[] = [
  {
    key: "fastq_id",
    help: "Unique identifier for each sequencing run.",
  },
  {
    key: "number_of_reads",
    help: "Total number of read pairs sequenced during this run.",
  },
  {
    key: "q30_barcode",
    help: q30BarcodeHelp,
  },
  {
    key: "q30_gem_barcode",
    help: "Fraction of GEM barcode bases with Q-score >= 30, excluding very low quality/no-call (Q <= 2) bases from the denominator.",
  },
  {
    key: "q30_probe_barcode",
    help: "Fraction of probe barcode bases (or antibody multiplexing barcode bases for Flex with Antibody Feature Barcode) with Q-score >= 30, excluding very low quality/no-call (Q <= 2) bases from the denominator.",
  },
  {
    key: "q30_umi",
    help: "Fraction of UMI bases with Q-score >= 30, excluding very low quality/no-call (Q <= 2) bases from the denominator.",
  },
  {
    key: "q30_read1",
    help: q30RnaRead1Help,
  },
  {
    key: "q30_read2",
    help: "Fraction of RNA Read 2 bases with Q-score >= 30, excluding very low quality/no-call (Q <= 2) bases from the denominator.",
  },
];

export const confidentlyMappedReadsInCellsHelp = [
  {
    msg: "The fraction of valid-barcode, valid-UMI, confidently-mapped-to-the-whole-probe-set reads with cell-associated barcodes.",
    isRtl: true,
  },
  {
    msg: "The fraction of valid-barcode, valid-UMI, confidently-mapped-to-transcriptome reads with cell-associated barcodes.",
    isRtl: false,
  },
];

const probeBarcodeMetricsTable: MetricTableConfig = {
  title: "Metrics per probe barcode",
  help: "Distribution of cells and UMIs across probe barcodes in this library.",
  entries: [
    {
      key: "probe_barcode_id",
      help: "The identifier of this probe barcode.",
    },
    {
      key: "sample_id",
      help: "The identifier of the sample associated with this probe barcode.",
    },
    {
      key: "umi_per_probe_barcode",
      help: "Number and fraction of UMIs for this probe barcode amongst all UMIs for that library type in the raw feature-barcode matrix.",
    },
    {
      key: "cells_per_probe_barcode",
      help: "Number and fraction of cells per probe barcode amongst all cells detected in this GEM well. Cell calling is based on gene expression data when present.",
    },
  ],
  orderBy: ["sample_id", "probe_barcode_id"],
};

export const gdnaMetricsTable: MetricTableConfig = {
  title: "UMIs from Genomic DNA",
  help: "Estimation of Genomic DNA level in this sample.",
  entries: [
    {
      key: "estimated_gdna_content",
      help: "The estimated fraction of filtered UMIs derived from genomic DNA based on the discordance between probes targeting exon-junction-spanning regions and non-exon-junction-spanning regions.",
    },
    {
      key: "estimated_gdna_unspliced_threshold",
      help: "The estimated number of UMIs derived from genomic DNA for each probe targeting non-exon-junction-spanning regions. A probe not spanning an exon junction with a total UMI count below this value has a high likelihood of its UMIs being derived primarily from hybridization to genomic DNA rather than the mRNA. For details, please visit http://10xgen.com/gdna-estimation",
    },
  ],
};

const LibraryShared = {
  VDJ: {
    libraryQuality: {
      title: "VDJ Library Quality",
      help: "Metrics related to overall library quality.",
      entries: [
        {
          key: "vdj_total_raw_read_pairs_per_filtered_bc",
          help: "The total number of sequenced read pairs divided by the number of cell-associated barcodes.",
        },
        {
          key: "vdj_assemblable_read_pairs_per_filtered_bc",
          help: "Mean number of read pairs used in assembly per cell-associated barcode. These reads must have a cell-associated barcode, map to a V(D)J gene, and have a UMI with sufficient read support.",
        },
        {
          key: "vdj_good_bc_frac",
          help: "Fraction of reads with barcodes that are present in the whitelist after barcode correction.",
        },
      ],
    },
    sequencing: {
      title: "VDJ Sequencing Run Quality",
      help: "Metrics per sequencing run.",
      entries: [
        {
          key: "fastq_id",
          help: "Unique identifier for each sequencing run.",
        },
        {
          key: "number_of_reads",
          help: "Total number of read pairs sequenced during this run.",
        },
        {
          key: "q30_barcode",
          help: q30BarcodeHelp,
        },
        {
          key: "q30_gem_barcode",
          help: "Fraction of GEM barcode bases with Q-score >= 30, excluding very low quality/no-call (Q <= 2) bases from the denominator.",
        },
        {
          key: "q30_probe_barcode",
          help: "Fraction of GEM barcode bases with Q-score >= 30, excluding very low quality/no-call (Q <= 2) bases from the denominator.",
        },
        {
          key: "q30_umi",
          help: "Fraction of UMI bases with Q-score >= 30, excluding very low quality/no-call (Q <= 2) bases from the denominator.",
        },
        {
          key: "q30_read1",
          help: q30RnaRead1Help,
        },
        {
          key: "q30_read2",
          help: "Fraction of RNA Read 2 bases with Q-score >= 30, excluding very low quality/no-call (Q <= 2) bases from the denominator.",
        },
      ],
      orderBy: ["fastq_id"],
    },
    perOcmBarcode: {
      title: "Metrics per OCM Barcode ID",
      help: "Distribution of VDJ cells across OCM barcodes in this library.",
      entries: [
        {
          key: "ocm_barcode_id",
          help: "The identifier of this OCM barcode.",
        },
        {
          key: "sample_id",
          help: "The identifier of the sample associated with this OCM barcode.",
        },
        {
          key: "vdj_cells_per_tag",
          help: "Number and fraction of VDJ cells per OCM Barcode ID amongst all VDJ cells detected in this GEM well.",
        },
      ],
      orderBy: ["ocm_barcode_id", "sample_id"],
    },
    perHashtag: {
      title: "Metrics per Hashtag ID",
      help: "Distribution of VDJ cells across Hashtags in this library.",
      entries: [
        {
          key: "hashtag_id",
          help: "The identifier of this Hashtag.",
        },
        {
          key: "sample_id",
          help: "The identifier of the sample associated with this Hashtag.",
        },
        {
          key: "vdj_cells_per_tag",
          help: "Number and fraction of VDJ cells per Hashtag ID amongst all VDJ cells detected in this GEM well.",
        },
      ],
      orderBy: ["sample_id", "hashtag_id"],
    },
  },
};

const MultiplexLibraryConfig = {
  GEX: {
    hero: {
      title: "Cell Calling Quality",
      help: cellCallingQualityHelp,
      entries: [
        {
          key: "cell_associated_partitions",
          help: "The number of barcodes identified by the cell-calling algorithm as containing a cell. Barcodes removed by Protein Aggregate Detection and Filtering or High Occupancy GEM Filtering are not counted.",
        },
        {
          key: "reads_in_cell_associated_partitions",
          help: confidentlyMappedReadsInCellsHelp,
        },
        {
          // Note: this metric's value is inverted during transformation so it
          // represents the fraction of cells not in high-occupancy GEMs.
          key: "rtl_multiplexing_fraction_cells_in_high_occupancy_gems",
          help: "Fraction of cell-associated barcodes from initial cell calls that remain after high occupancy GEM filtering. Cell calling is performed and all barcodes associated with any GEMs that have significantly higher probe barcodes per GEM than we would expect from optimal chip loading are removed to mitigate higher than expected barcode collision rates.",
        },
      ],
    },
    mapping: {
      title: "Mapping Quality",
      help: "Mapping metrics calculated based on all reads in the library.",
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
    libraryQuality: {
      title: "Library Quality",
      help: "Metrics related to overall library quality.",
      entries: [
        {
          key: "mean_reads_per_cell_associated_partition",
          help: "The total number of sequenced read pairs divided by the number of cell-associated barcodes.",
        },
        {
          key: "sequencing_saturation",
          help: "The fraction of reads originating from an already-observed UMI. This fraction is a function of library complexity and sequencing depth. It is calculated as the ratio of reads with non-unique (cell barcode, UMI, gene) combinations to the total number of valid reads (those with recognized cell barcode, UMI, and genes).",
        },
        {
          key: "valid_barcodes",
          help: validBarcodesHelp,
        },
        {
          key: "valid_gem_barcodes",
          help: "Fraction of reads with GEM barcodes that are present in the whitelist after barcode correction.",
        },
        {
          key: "valid_probe_barcodes",
          help: "Fraction of reads with probe barcodes that are present in the whitelist after barcode correction.",
        },
        {
          key: "valid_umis",
          help: "Fraction of reads with valid UMI sequences; i.e. UMI sequences that do not contain Ns and that are not homopolymers.",
        },
      ],
    },
    sequencing: {
      title: "Sequencing Run Quality",
      help: "Metrics per sequencing run.",
      entries: [
        {
          key: "fastq_id",
          help: "Unique identifier for each sequencing run.",
        },
        {
          key: "number_of_reads",
          help: "Total number of read pairs sequenced during this run.",
        },
        {
          key: "q30_barcode",
          help: q30BarcodeHelp,
        },
        {
          key: "q30_gem_barcode",
          help: "Fraction of GEM barcode bases with Q-score >= 30, excluding very low quality/no-call (Q <= 2) bases from the denominator.",
        },
        {
          key: "q30_probe_barcode",
          help: "Fraction of probe barcode bases (or antibody multiplexing barcode bases for Flex with Antibody Feature Barcode) with Q-score >= 30, excluding very low quality/no-call (Q <= 2) bases from the denominator.",
        },
        {
          key: "q30_umi",
          help: "Fraction of UMI bases with Q-score >= 30, excluding very low quality/no-call (Q <= 2) bases from the denominator.",
        },
        {
          key: "q30_read1",
          help: q30RnaRead1Help,
        },
        {
          key: "q30_read2",
          help: "Fraction of RNA Read 2 bases with Q-score >= 30, excluding very low quality/no-call (Q <= 2) bases from the denominator.",
        },
      ],
      orderBy: ["fastq_id"],
    },
    perRtlProbeBarcode: probeBarcodeMetricsTable,
    perOcmBarcode: {
      title: "Metrics per OCM Barcode",
      help: "Distribution of cells and UMIs across OCM barcodes in this library.",
      entries: [
        {
          key: "ocm_barcode_id",
          help: "The identifier of this OCM barcode.",
        },
        {
          key: "sample_id",
          help: "The identifier of the sample associated with this OCM barcode.",
        },
        {
          key: "umi_per_ocm_barcode",
          help: "Number and fraction of UMIs for this overhang amongst all UMIs for that library type in the raw feature-barcode matrix.",
        },
        {
          key: "cells_per_ocm_barcode",
          help: "Number and fraction of cells per OCM Barcode ID amongst all cells detected in this GEM well. Cell calling is based on gene expression data when present.",
        },
      ],
      orderBy: ["ocm_barcode_id", "sample_id"],
    },
    gdna: gdnaMetricsTable,
  },
  AB: {
    hero: {
      title: "Antibody Expression",
      help: "Summary statistics about the detection and association with cells for antibody reads in the library.",
      entries: [
        {
          key: "cell_associated_partitions",
          help: "The number of barcodes identified by the cell-calling algorithm as containing a cell. Barcodes removed by Protein Aggregate Detection and Filtering or High Occupancy GEM Filtering are not counted.",
        },
        {
          key: "fraction_antibody_reads",
          help: "Fraction of read pairs that contain a recognized antibody Feature Barcode.",
        },
        {
          key: "reads_in_cell_associated_partitions",
          help: "The fraction of valid-barcode, valid-UMI, recognized antibody Feature Barcode reads with cell-associated barcodes.",
        },
        {
          key: "fraction_antibody_reads_usable",
          help: "Fraction of read pairs that contain a recognized antibody Feature Barcode, a valid UMI, and a cell-associated barcode.",
        },
        {
          key: "fraction_reads_in_aggregate_barcodes",
          help: "Fraction of read pairs with valid barcodes that were removed because they are aggregates.",
        },
      ],
    },
    libraryQuality: {
      title: "Antibody Library Quality",
      help: "Metrics related to overall library quality.",
      entries: featureBarcodingLibraryQualityTableEntries,
    },
    sequencing: {
      title: "Antibody Sequencing Run Quality",
      help: "Metrics per sequencing run.",
      entries: featureBarcodingSequencingMetricsTableEntries,
      orderBy: ["fastq_id"],
    },
    rtl_probe_barcode_metrics: probeBarcodeMetricsTable,
  },
  CRISPR: {
    hero: {
      title: "CRISPR Guide Expression",
      help: "Summary statistics about the detection and association with cells for CRISPR reads in this library.",
      entries: [
        {
          key: "cell_associated_partitions",
          help: "The number of barcodes identified by the cell-calling algorithm as containing a cell. Barcodes removed by Protein Aggregate Detection and Filtering or High Occupancy GEM Filtering are not counted.",
        },
        {
          key: "fraction_guide_reads",
          help: "Fraction of CRISPR library reads with a recognized protospacer sequence.",
        },
        {
          key: "reads_in_cell_associated_partitions",
          help: "Among CRISPR library reads with a recognized protospacer sequence, a valid UMI, and a valid barcode, the fraction with cell-associated barcodes.",
        },
        {
          key: "mean_reads_per_cell_associated_partition",
          help: "The total number of sequenced read pairs divided by the number of cell-associated barcodes.",
        },
        {
          key: "fraction_guide_reads_usable",
          help: "Fraction of CRISPR library reads with a recognized protospacer sequence, a valid UMI, and a cell-associated barcode.",
        },
        {
          key: "fraction_reads_with_putative_protospacer",
          help: "Fraction of CRISPR library reads from which a putative protospacer sequence could be extracted.",
        },
        {
          key: "fraction_protospacer_not_recognized",
          help: "Among all CRISPR library reads with a putative protospacer sequence, the fraction with a protospacer sequence that did not match any specified in the Feature Reference CSV file provided to Cell Ranger.",
        },
      ],
    },
    libraryQuality: {
      title: "CRISPR Library Quality",
      help: "Metrics related to overall library quality.",
      entries: featureBarcodingLibraryQualityTableEntries,
    },
    sequencing: {
      title: "CRISPR Sequencing Run Quality",
      help: "Metrics per sequencing run.",
      entries: featureBarcodingSequencingMetricsTableEntries,
      orderBy: ["fastq_id"],
    },
    rtl_probe_barcode_metrics: probeBarcodeMetricsTable,
  },
  VDJ_T: {
    enrichment: {
      title: "Enrichment",
      help: "Summary statistics about efficiency of V(D)J enrichment measured by the fraction of reads mapped to V(D)J genes.",
      entries: [
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
    cellCalling: {
      title: "Cell Calling Quality",
      help: cellCallingQualityHelp,
      entries: [
        {
          key: "vdj_filtered_bcs",
          help: "The number of barcodes estimated to be associated with cells that express targeted V(D)J transcripts.",
        },
        {
          key: "multi_vdj_assembly_contig_pair_productive_full_len_bc_frac",
          help: "Fraction of cell-associated barcodes with at least one productive contig for each chain of the receptor pair. A productive contig satisfies the following conditions: the contig annotations span the 5' end of the V region to the 3' end of the J region of the chain, a start codon was found in the expected part of the V sequence, an in-frame CDR3 amino acid motif was found, and no stop codons were found in the aligned V-J region.",
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
    ...LibraryShared.VDJ,
  },
  VDJ_T_GD: {
    enrichment: {
      title: "Enrichment",
      help: "Summary statistics about efficiency of V(D)J enrichment measured by the fraction of reads mapped to V(D)J genes.",
      entries: [
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
    cellCalling: {
      title: "Cell Calling Quality",
      help: cellCallingQualityHelp,
      entries: [
        {
          key: "vdj_filtered_bcs",
          help: "The number of barcodes estimated to be associated with cells that express targeted V(D)J transcripts.",
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
    ...LibraryShared.VDJ,
  },
  VDJ_B: {
    enrichment: {
      title: "Enrichment",
      help: "Summary statistics about efficiency of V(D)J enrichment measured by the fraction of reads mapped to V(D)J genes.",
      entries: [
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
    cellCalling: {
      title: "Cell Calling Quality",
      help: cellCallingQualityHelp,
      entries: [
        {
          key: "vdj_filtered_bcs",
          help: "The number of barcodes estimated to be associated with cells that express targeted V(D)J transcripts.",
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
    ...LibraryShared.VDJ,
  },
  AG: {
    antigen_reads_table: {
      title: "Antigen Reads",
      entries: [
        {
          key: "fraction_antigen_reads",
          help: "Fraction of read pairs that contain a recognized antigen-barcode.",
        },
        {
          key: "reads_in_cell_associated_partitions",
          help: "The fraction of valid-barcode, valid-UMI, recognized antigen-barcode reads with cell-associated barcodes.",
        },
        {
          key: "fraction_antigen_reads_usable",
          help: "Fraction of read pairs that contain a recognized antigen-barcode, a valid UMI, and a cell-associated barcode.",
        },
        {
          key: "fraction_reads_in_aggregate_barcodes",
          help: "Fraction of read pairs with valid barcodes that were removed because they are aggregates.",
        },
      ],
    },
    libraryQuality: {
      title: "Antigen Library Quality",
      help: "Metrics related to overall library quality.",
      entries: featureBarcodingLibraryQualityTableEntries,
    },
    sequencing: {
      title: "Antigen Sequencing Run Quality",
      help: "Metrics per sequencing run.",
      entries: featureBarcodingSequencingMetricsTableEntries,
      orderBy: ["fastq_id"],
    },
  },
  CMO: {
    hero: {
      title: "Multiplexing Quality",
      entries: [
        {
          key: "cell_associated_partitions",
          help: "Cell-associated barcodes include cell barcodes that are either (1) assigned to a sample, (2) identified as multiplets, or (3) not assigned any Hashtags or CMOs. For a more detailed explanation, please see Technical Note CG000475 on https://www.10xgenomics.com/support",
        },
        {
          key: "samples_assigned_at_least_one_singlet",
          help: "Number of samples to which at least one cell was assigned. Only cell-associated barcodes assigned exactly one CMO were assigned to a sample.",
        },
        {
          key: "singlet_capture_ratio",
          help: "Ratio between the number of singlets (i.e. cell-associated barcodes assigned exactly one CMO) obtained and the number of singlets expected in this experiment.",
        },
        {
          key: "singlets_assigned_to_a_sample",
          help: "Number and fraction of cells assigned to a sample amongst all cells detected in this GEM well. Only cell-associated barcodes assigned exactly one CMO were assigned to a sample.",
        },
        {
          key: "cell_associated_partitions_identified_as_multiplets",
          help: "Cell-associated barcodes that were assigned more than one CMO and hence determined to be multiplets.",
        },
        {
          key: "cell_associated_partitions_not_assigned_any_tags",
          help: "Cell-associated barcodes that either (i) did not have enough CMO molecules above background or (ii) could not be confidently assigned to a singlet or multiplet state.",
        },
        {
          key: "median_cmo_umis_per_singlet",
          help: "Median number of CMO UMIs captured per cell-associated barcode assigned exactly one CMO.",
        },
        {
          key: "fraction_cmo_reads",
          help: "Fraction of reads that contain a recognized CMO sequence.",
        },
        {
          key: "reads_in_cell_associated_partitions",
          help: "The fraction of valid-barcode, valid-UMI, recognized multiplexing-barcode reads with cell-associated barcodes.",
        },
        {
          key: "fraction_cmo_reads_usable",
          help: "Fraction of read pairs that contain a recognized CMO sequence, a valid UMI, and a cell-associated barcode.",
        },
        {
          key: "fraction_reads_from_multiplets",
          help: "Amongst all sequenced read pairs, fraction with a cell-barcode identified as a multiplet.",
        },
      ],
    },
    perCmoBarcode: {
      title: "Metrics per CMO",
      entries: [
        {
          key: "gem_well_cmo",
          help: "Metrics in this table are provided for each CMO.",
        },
        {
          key: "sample_id",
          help: "The identifier of the sample associated with this CMO.",
        },
        {
          key: "cmo_signal_to_background_ratio",
          help: "Computed as the difference between labeled and unlabelled mean CMO counts (log scale) divided by the variance.",
        },
        {
          key: "singlets_assigned_to_cmo",
          help: "Fraction of cells assigned this particular CMO (and only this CMO) amongst all cells detected in this GEM well.",
        },
        {
          key: "cmo_reads_in_cell_associated_partitions",
          help: "Amongst all reads with a valid barcode, valid UMI, and this particular CMO sequence, fraction arising from cell-containing partitions.",
        },
      ],
      orderBy: ["sample_id", "gem_well_cmo"],
    },
    libraryQuality: {
      title: "CMO Library Quality",
      help: "Metrics related to overall library quality.",
      entries: [
        {
          key: "mean_reads_per_cell_associated_partition",
          help: "The total number of sequenced read pairs divided by the number of cell-associated barcodes.",
        },
        {
          key: "valid_barcodes",
          help: validBarcodesHelp,
        },
        {
          key: "valid_umis",
          help: "Fraction of reads with valid UMI sequences; i.e. UMI sequences that do not contain Ns and that are not homopolymers.",
        },
      ],
    },
    sequencing: {
      title: "CMO Sequencing Run Quality",
      help: "Metrics per sequencing run.",
      entries: [
        {
          key: "fastq_id",
          help: "Unique identifier for each sequencing run.",
        },
        {
          key: "number_of_reads",
          help: "Total number of read pairs sequenced during this run.",
        },
        {
          key: "q30_barcode",
          help: q30BarcodeHelp,
        },
        {
          key: "q30_umi",
          help: "Fraction of UMI bases with Q-score >= 30, excluding very low quality/no-call (Q <= 2) bases from the denominator.",
        },
        {
          key: "q30_read1",
          help: q30RnaRead1Help,
        },
        {
          key: "q30_read2",
          help: "Fraction of RNA Read 2 bases with Q-score >= 30, excluding very low quality/no-call (Q <= 2) bases from the denominator.",
        },
      ],
      orderBy: ["fastq_id"],
    },
  },
  Hashtag: {
    hero: {
      title: "Multiplexing Quality",
      help: "Summary statistics about sample multiplexing and cell assignments.",
      entries: [
        {
          key: "cell_associated_partitions",
          help: "Number of cell-associated barcodes called as containing one or more cells.",
        },
        {
          key: "samples_assigned_at_least_one_singlet",
          help: "Number of samples to which at least one cell was assigned. Only cell-associated barcodes assigned exactly one Hashtag were assigned to a sample.",
        },
        {
          key: "singlets_assigned_to_a_sample",
          help: "Number and fraction of cells assigned to a sample amongst all cells detected in this GEM well. Only cell-associated barcodes assigned exactly one Hashtag were assigned to a sample.",
        },
        {
          key: "cell_associated_partitions_identified_as_multiplets",
          help: "Cell-associated barcodes that were assigned more than one Hashtag and hence determined to be multiplets.",
        },
        {
          key: "cell_associated_partitions_not_assigned_any_tags",
          help: "Cell-associated barcodes that either (i) did not have enough Hashtag molecules above background or (ii) could not be confidently assigned to a singlet or multiplet state.",
        },
        {
          key: "median_hashtag_umis_per_singlet",
          help: "Median number of Hashtag UMIs captured per cell-associated barcode assigned exactly one Hashtag.",
        },
        {
          key: "hashtag_singlet_capture_ratio",
          help: "Ratio between the number of singlets (i.e. cell-associated barcodes assigned exactly one Hashtag) obtained and the number of singlets expected in this experiment according to Poisson statistics.",
        },
      ],
    },
    perHashtag: {
      title: "Metrics per Hashtag",
      help: "Cell assignments and signal quality metrics across each hashtag.",
      entries: [
        {
          key: "gem_well_hashtag",
          help: "Metrics in this table are provided for each Hashtag.",
        },
        {
          key: "sample_id",
          help: "The identifier of the sample associated with this hashtag.",
        },
        {
          key: "hashtag_reads_in_cell_associated_partitions",
          help: "Amongst all reads with a valid barcode, valid UMI for this particular Hashtag, fraction arising from cell-containing partitions.",
        },
        {
          key: "singlets_assigned_to_hashtag",
          help: "Fraction of cells assigned this particular Hashtag (and only this Hashtag) amongst all cells detected in this GEM well.",
        },
        {
          key: "hashtag_signal_to_background_ratio",
          help: "Computed as the difference between labeled and unlabelled mean Hashtag counts (log scale) divided by the variance.",
        },
      ],
      orderBy: ["sample_id", "gem_well_hashtag"],
    },
  },
  Custom: {
    hero: {
      title: "Feature Expression",
      entries: [
        {
          key: "fraction_feature_reads",
          help: "Fraction of reads that contain a recognized feature-barcode sequence.",
        },
        {
          key: "reads_in_cell_associated_partitions",
          help: "The fraction of valid-barcode, valid-UMI, recognized feature-barcode reads with cell-associated barcodes.",
        },
        {
          key: "fraction_feature_reads_usable",
          help: "Fraction of read pairs that contain a recognized feature-barcode, a valid UMI, and a cell-associated barcode",
        },
      ],
    },
    libraryQuality: {
      title: "Library Quality",
      help: "Metrics related to overall library quality.",
      entries: featureBarcodingLibraryQualityTableEntries,
    },
    sequencing: {
      title: "Sequencing Run Quality",
      help: "Metrics per sequencing run.",
      entries: featureBarcodingSequencingMetricsTableEntries,
      orderBy: ["fastq_id"],
    },
  },
};

export default MultiplexLibraryConfig;
