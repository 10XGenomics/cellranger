/* eslint-disable sort-keys */
import MultiplexLibraryConfig, {
  cellCallingQualityHelp,
  gdnaMetricsTable,
} from "./multiplex-library";
import MultiplexSampleConfig from "./multiplex-sample";
import { MetricEntryConfig, MetricTableConfig } from "./types";

const selectConfigs = (
  table: MetricTableConfig,
  category: "Library" | "Cells",
  ...keys: string[]
): MetricEntryConfig[] =>
  keys.map((key) => {
    const cfg = table.entries.find((cfg) => cfg.key === key);
    if (cfg === undefined) {
      throw new Error(`Config key not found: ${key}`);
    }
    return { ...cfg, category: category };
  });

const copyConfigs = (
  table: MetricTableConfig,
  category: "Library" | "Cells",
): MetricEntryConfig[] =>
  table.entries.map((cfg) => ({ ...cfg, category: category }));

const libraryQuality: MetricTableConfig = {
  title: "Library Quality",
  help: MultiplexLibraryConfig.GEX.libraryQuality.help,
  entries: selectConfigs(
    MultiplexLibraryConfig.GEX.libraryQuality,
    "Library",
    "mean_reads_per_cell_associated_partition",
    "sequencing_saturation",
    "valid_barcodes",
    "valid_umis",
  ),
};

const sequencingQuality: MetricTableConfig = {
  title: "Sequencing Run Quality",
  help: MultiplexLibraryConfig.GEX.sequencing.help,
  entries: selectConfigs(
    MultiplexLibraryConfig.GEX.sequencing,
    "Library",
    "fastq_id",
    "number_of_reads",
    "q30_barcode",
    "q30_umi",
    "q30_read1",
    "q30_read2",
  ),
  orderBy: ["fastq_id"],
};

const vdjLibraryQuality: MetricTableConfig = {
  title: "VDJ Library Quality",
  help: MultiplexLibraryConfig.VDJ_B.libraryQuality.help,
  entries: selectConfigs(
    MultiplexLibraryConfig.VDJ_B.libraryQuality,
    "Library",
    "vdj_total_raw_read_pairs_per_filtered_bc",
    "vdj_assemblable_read_pairs_per_filtered_bc",
    "vdj_good_bc_frac",
  ),
};

const vdjSequencingQuality: MetricTableConfig = {
  title: "VDJ Sequencing Run Quality",
  help: MultiplexLibraryConfig.VDJ_B.sequencing.help,
  entries: selectConfigs(
    MultiplexLibraryConfig.VDJ_B.sequencing,
    "Library",
    "fastq_id",
    "number_of_reads",
    "q30_barcode",
    "q30_umi",
    "q30_read1",
    "q30_read2",
  ),
};

const vdjClonotypeMetrics: MetricTableConfig = {
  title: "Clonotype Metrics",
  entries: copyConfigs(
    MultiplexSampleConfig.VDJ_T.pairedClonotypeDiversity,
    "Cells",
  ),
};

export const SinglePlexConfig = {
  GEX: {
    cellCalling: {
      title: "Cell Calling Quality",
      help: MultiplexLibraryConfig.GEX.hero.help,
      entries: [
        ...selectConfigs(
          MultiplexLibraryConfig.GEX.hero,
          "Library",
          "cell_associated_partitions",
          "reads_in_cell_associated_partitions",
        ),
        ...selectConfigs(
          MultiplexSampleConfig.GEX.hero,
          "Cells",
          "median_genes_per_singlet",
          "median_umi_per_singlet",
          "total_genes_detected",
        ),
      ],
    },
    cellCallingBarnyard: {
      title: "Cell Calling Quality",
      help: MultiplexLibraryConfig.GEX.hero.help,
      entries: [
        ...selectConfigs(
          MultiplexLibraryConfig.GEX.hero,
          "Library",
          "cell_associated_partitions",
          "reads_in_cell_associated_partitions",
        ),
        {
          key: "filtered_bcs_observed_all",
          help: "The number of GEMs with at least one cell.",
        },
        {
          key: "filtered_bcs_inferred_multiplets",
          help: "The number of GEMs with at least two cells.",
        },
        {
          key: "filtered_bcs_inferred_multiplet_rate",
          help: "The fraction of GEMs with two or more cells.",
        },
      ],
    },
    barnyardPerGenome: {
      title: "Cells per genome",
      entries: [
        ...selectConfigs(
          MultiplexSampleConfig.GEX.hero,
          "Cells",
          "genome",
          "total_singlets",
          "confidently_mapped_reads_in_cells",
          "median_genes_per_singlet",
          "median_umi_per_singlet",
          "total_genes_detected",
        ),
        {
          key: "mean_reads_per_cell",
          help: "Mean number of read pairs sequenced from the cells called from this sample.",
          category: "Cells" as const,
        },
      ],
      orderBy: ["total_singlets"],
    },
    mapping: {
      title: "Mapping Quality",
      help: MultiplexLibraryConfig.GEX.mapping.help,
      entries: copyConfigs(MultiplexLibraryConfig.GEX.mapping, "Library"),
    },
    libraryQuality: libraryQuality,
    sequencing: sequencingQuality,
    gdna_metrics: gdnaMetricsTable,
  },
  VDJ_B: {
    cellCalling: {
      title: "Cell Calling Quality",
      help: cellCallingQualityHelp,
      entries: [
        ...selectConfigs(
          MultiplexLibraryConfig.VDJ_B.cellCalling,
          "Library",
          "vdj_filtered_bcs",
        ),
        ...selectConfigs(
          MultiplexSampleConfig.VDJ_B.hero,
          "Cells",
          "multi_vdj_assembly_contig_pair_productive_full_len_bc_frac",
          "IGK_IGH_vdj_assembly_contig_pair_productive_full_len_bc_frac",
          "IGL_IGH_vdj_assembly_contig_pair_productive_full_len_bc_frac",
          "IGH_vdj_assembly_prod_cdr_bc_frac",
          "IGK_vdj_assembly_prod_cdr_bc_frac",
          "IGL_vdj_assembly_prod_cdr_bc_frac",
        ),
        ...selectConfigs(
          MultiplexLibraryConfig.VDJ_B.cellCalling,
          "Library",
          "vdj_filtered_bcs_cum_frac",
        ),
      ],
    },
    enrichmentSensitivity: {
      title: "VDJ Sensitivity and Enrichment",
      help: MultiplexSampleConfig.VDJ_B.sensitivityEnrichment.help,
      entries: [
        ...copyConfigs(MultiplexSampleConfig.VDJ_B.sensitivity, "Cells"),
        ...copyConfigs(MultiplexLibraryConfig.VDJ_B.enrichment, "Library"),
      ],
    },
    libraryQuality: vdjLibraryQuality,
    sequencingQuality: vdjSequencingQuality,
    clonotype: vdjClonotypeMetrics,
  },
  VDJ_T: {
    cellCalling: {
      title: "Cell Calling Quality",
      help: cellCallingQualityHelp,
      entries: [
        ...selectConfigs(
          MultiplexLibraryConfig.VDJ_T.cellCalling,
          "Library",
          "vdj_filtered_bcs",
        ),
        {
          key: "multi_vdj_assembly_contig_pair_productive_full_len_bc_count",
          help: "Number of cell barcodes for which at least 1 full-length productive sequence was found for each chain of the (TRA, TRB) receptor pair.",
        },
        ...selectConfigs(
          MultiplexSampleConfig.VDJ_T.hero,
          "Cells",
          "TRA_TRB_vdj_assembly_contig_pair_productive_full_len_bc_frac",
          "TRA_vdj_assembly_prod_cdr_bc_frac",
          "TRB_vdj_assembly_prod_cdr_bc_frac",
        ),
        ...selectConfigs(
          MultiplexLibraryConfig.VDJ_T.cellCalling,
          "Library",
          "vdj_filtered_bcs_cum_frac",
        ),
      ],
    },
    enrichmentSensitivity: {
      title: "VDJ Enrichment and Sensitivity",
      help: MultiplexSampleConfig.VDJ_T.sensitivityEnrichment.help,
      entries: [
        ...copyConfigs(MultiplexSampleConfig.VDJ_T.sensitivity, "Cells"),
        ...copyConfigs(MultiplexLibraryConfig.VDJ_T.enrichment, "Library"),
      ],
    },
    libraryQuality: vdjLibraryQuality,
    sequencingQuality: vdjSequencingQuality,
    clonotype: vdjClonotypeMetrics,
  },
  VDJ_T_GD: {
    cellCalling: {
      title: "Cell Calling Quality",
      help: cellCallingQualityHelp,
      entries: [
        ...selectConfigs(
          MultiplexLibraryConfig.VDJ_T_GD.cellCalling,
          "Library",
          "vdj_filtered_bcs",
        ),
        {
          key: "multi_vdj_assembly_contig_pair_productive_full_len_bc_count",
          help: "Number of cell barcodes for which at least 1 full-length productive sequence was found for each chain of the (TRG, TRD) receptor pair.",
        },
        ...selectConfigs(
          MultiplexSampleConfig.VDJ_T_GD.hero,
          "Cells",
          "TRG_TRD_vdj_assembly_contig_pair_productive_full_len_bc_frac",
          "TRG_vdj_assembly_prod_cdr_bc_frac",
          "TRD_vdj_assembly_prod_cdr_bc_frac",
        ),
        ...selectConfigs(
          MultiplexLibraryConfig.VDJ_T_GD.cellCalling,
          "Library",
          "vdj_filtered_bcs_cum_frac",
        ),
      ],
    },
    enrichmentSensitivity: {
      title: "VDJ Sensitivity and Enrichment",
      help: MultiplexSampleConfig.VDJ_T_GD.sensitivityEnrichment.help,
      entries: [
        ...copyConfigs(MultiplexSampleConfig.VDJ_T_GD.sensitivity, "Cells"),
        ...copyConfigs(MultiplexLibraryConfig.VDJ_T_GD.enrichment, "Library"),
      ],
    },
    libraryQuality: vdjLibraryQuality,
    sequencingQuality: vdjSequencingQuality,
    clonotype: vdjClonotypeMetrics,
  },
  AB: {
    hero: {
      title: "Antibody Expression",
      help: MultiplexLibraryConfig.AB.hero.help,
      entries: [
        ...selectConfigs(
          MultiplexLibraryConfig.AB.hero,
          "Library",
          "cell_associated_partitions",
          "fraction_antibody_reads",
          "reads_in_cell_associated_partitions",
          "fraction_antibody_reads_usable",
        ),
        ...selectConfigs(
          MultiplexSampleConfig.AB.hero,
          "Cells",
          "antibody_reads_usable_per_cell",
        ),
        ...selectConfigs(
          MultiplexLibraryConfig.AB.hero,
          "Library",
          "fraction_reads_in_aggregate_barcodes",
        ),
        ...selectConfigs(
          MultiplexSampleConfig.AB.hero,
          "Cells",
          "median_umis_per_singlet",
        ),
      ],
    },
    libraryQuality: libraryQuality,
    sequencingQuality: sequencingQuality,
  },
  AG: {
    hero: {
      title: "Antigen Reads",
      entries: [
        ...copyConfigs(
          MultiplexLibraryConfig.AG.antigen_reads_table,
          "Library",
        ),
        {
          key: "fraction_unknown_antigen",
          help: "Fraction of read pairs with an unrecognized antigen-barcode sequence.",
          category: "Library",
        } as MetricEntryConfig,
      ],
    },
    expression: {
      title: "Antigen Expression",
      entries: copyConfigs(MultiplexSampleConfig.AG.hero, "Cells"),
    },
    libraryQuality: libraryQuality,
    sequencingQuality: sequencingQuality,
  },
  CRISPR: {
    hero: {
      title: "CRISPR Guides Expression",
      help: MultiplexLibraryConfig.CRISPR.hero.help,
      entries: [
        ...selectConfigs(
          MultiplexLibraryConfig.CRISPR.hero,
          "Library",
          "fraction_guide_reads",
          "reads_in_cell_associated_partitions",
          "fraction_guide_reads_usable",
          "fraction_reads_with_putative_protospacer",
          "fraction_protospacer_not_recognized",
        ),
        ...selectConfigs(
          MultiplexSampleConfig.CRISPR.hero,
          "Cells",
          "cells_with_one_or_more_protospacers_detected",
          "cells_with_two_or_more_protospacers_detected",
          "median_umis_per_singlet",
        ),
      ],
    },
    libraryQuality: libraryQuality,
    sequencingQuality: sequencingQuality,
  },
  Custom: {
    hero: {
      title: "Feature Expression",
      entries: [
        ...copyConfigs(MultiplexLibraryConfig.Custom.hero, "Library"),
        ...selectConfigs(
          MultiplexSampleConfig.Custom.hero,
          "Cells",
          "feature_reads_usable_per_cell",
          "median_umis_per_singlet",
        ),
      ],
    },
    libraryQuality: libraryQuality,
    sequencingQuality: sequencingQuality,
  },
};
