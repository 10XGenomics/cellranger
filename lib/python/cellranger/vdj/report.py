#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

from __future__ import annotations

from collections import defaultdict

import numpy as np
import pandas as pd
from six import ensure_str

import cellranger.constants as cr_constants
import cellranger.library_constants as lib_constants
import cellranger.report as cr_report  # pylint: disable=no-name-in-module
import cellranger.segment as cr_segment
import cellranger.utils as cr_utils
import cellranger.vdj.annotations as vdj_annot
import cellranger.vdj.constants as vdj_constants
import cellranger.vdj.reference as vdj_reference
import cellranger.vdj.utils as vdj_utils

VDJ_MISMATCH_BREAKS = list(range(0, 22, 2))
VDJ_DEPTH_BREAKS = list(range(0, 20, 1))
VDJ_SOFTCLIP_BREAKS = list(range(0, 52, 2))
VDJ_CONTIG_LENGTH_PERCENTILES = np.arange(0, 1.1, 0.05)

from collections.abc import Collection, Iterable
from typing import Any


def vdj_metrics_dict(prefix=""):
    """Construct the metrics dictionary for VDJ metrics."""
    # This is a method, rather than a global, so that it doesn't have to
    # construct the giant dictionary whenever the module is imported.

    # shared sub-parts of VDJ_METRICS, to reduce redundancy and memory usage.
    prefixes = "prefixes"
    prefix_list = [
        "vdj_genes",
        "gem_groups",
        "vdj_gene_pairs",
        "canonical_vdj_gene_pairs",
        "vdj_clonotype_types",
        "primer_names_plus_rc",
    ]
    vdj_genes = prefix_list[:1]
    vdj_genes_prefix = {prefixes: vdj_genes}
    count_gem_groups = (cr_report.CountMetric, {prefixes: prefix_list[1:2]})
    vdj_gene_pairs_prefix = {prefixes: prefix_list[2:3]}
    vdj_canonical_and_clonotype_prefix = {prefixes: prefix_list[3:5]}
    vdj_clonotype_types = prefix_list[4:5]
    primer_names_prefix = {prefixes: prefix_list[5:6]}
    args = "args"
    contig_length_hist = (
        cr_report.HistogramMetric,
        {args: [VDJ_CONTIG_LENGTH_PERCENTILES]},
    )
    depth_break_hist = (cr_report.HistogramMetric, {args: [VDJ_DEPTH_BREAKS]})
    mismatch_break_hist = (cr_report.HistogramMetric, {args: [VDJ_MISMATCH_BREAKS]})
    softclip_break_hist = (cr_report.HistogramMetric, {args: [VDJ_SOFTCLIP_BREAKS]})
    empty_mean_metric = (cr_report.MeanMetric, {})
    empty_median_metric = (cr_report.MedianMetric, {})
    empty_percent_metric = (cr_report.PercentMetric, {})
    percent_genes_metric = (cr_report.PercentMetric, vdj_genes_prefix)
    median_genes_metric = (cr_report.MedianMetric, vdj_genes_prefix)
    topn_clonotype_metric = (
        cr_report.TopNMetric,
        {args: [vdj_constants.VDJ_MAX_OUTPUT_CLONOTYPES], prefixes: vdj_clonotype_types},
    )
    percent_clonotype_metric = (cr_report.PercentMetric, {prefixes: vdj_clonotype_types})
    percent_gene_pairs_metric = (
        cr_report.PercentMetric,
        vdj_gene_pairs_prefix,
    )
    percent_canonical_clono_metric = (cr_report.PercentMetric, vdj_canonical_and_clonotype_prefix)
    diversity_canonical_clono_metric = (
        cr_report.EffectiveDiversityMetric,
        vdj_canonical_and_clonotype_prefix,
    )

    ret_dict = {
        "barcode_reads": (
            cr_report.DictionaryMetric,
            {
                "kwargs": {"report_type": "barcodes"},
                prefixes: ["references", "regions", "read_types"],
            },
        ),
        # Dummy metrics to add canonical VDJ gene prefixes
        "vdj_dummy_metric": (
            cr_report.CountMetric,
            {
                prefixes: [
                    "canonical_vdj_genes",
                    "canonical_vdj_gene_pairs",
                    "canonical_vdj_genes_nomulti",
                    "canonical_vdj_gene_pairs_nomulti",
                ]
            },
        ),
        # VDJ general metrics
        "vdj_filtered_bcs": (cr_report.CountMetric, {}),
        "vdj_filtered_bcs_cum_frac": empty_percent_metric,
        "vdj_filter_bcs_rpu_threshold": count_gem_groups,
        "vdj_filter_bcs_umi_threshold": count_gem_groups,
        "vdj_filter_bcs_confidence": count_gem_groups,
        "vdj_corrected_bc_frac": empty_percent_metric,
        "vdj_corrected_umi_frac": empty_percent_metric,
        "vdj_good_bc_frac": empty_percent_metric,
        "vdj_total_raw_reads_per_filtered_bc": empty_percent_metric,
        "vdj_total_raw_read_pairs_per_filtered_bc": empty_percent_metric,
        "vdj_assemblable_read_pairs_per_filtered_bc": empty_percent_metric,
        "vdj_filtered_bcs_relative_difference_from_recovered_cells": empty_percent_metric,
        "vdj_sequencing_efficiency": empty_percent_metric,
        # VDJ read-filter metrics
        "vdj_recombinome_mapped_reads_frac": percent_genes_metric,
        "vdj_recombinome_mapped_read_pairs_frac": percent_genes_metric,
        "vdj_recombinome_antisense_read_pairs_frac": percent_genes_metric,
        "vdj_recombinome_mapped_read1_frac": percent_genes_metric,
        "vdj_recombinome_mapped_read2_frac": percent_genes_metric,
        "vdj_recombinome_chimeric_read_pairs_frac": percent_genes_metric,
        "vdj_recombinome_total_umis_per_cell_distribution": (
            cr_report.NumpyHistogramMetric,
            {args: [100]},
        ),
        "vdj_recombinome_total_umis_per_cell_median": empty_median_metric,
        "vdj_recombinome_umis_per_cell_distribution": (
            cr_report.NumpyHistogramMetric,
            {args: [100], prefixes: vdj_genes},
        ),
        "vdj_recombinome_umis_per_cell_median": median_genes_metric,
        # VDJ UMI filtering metrics
        "vdj_recombinome_readpairs_per_umi_distribution": (
            cr_report.DictionaryMetric,
            vdj_genes_prefix,
        ),
        "vdj_recombinome_low_support_reads_frac": percent_genes_metric,
        "vdj_recombinome_readpairs_per_umi_n50": (
            cr_report.CountMetric,
            {prefixes: prefix_list[:2]},
        ),
        # VDJ primer trimming metrics
        "vdj_trim_read1_frac": (cr_report.PercentMetric, primer_names_prefix),
        "vdj_trim_read2_frac": (cr_report.PercentMetric, primer_names_prefix),
        # VDJ assembly metrics
        "vdj_assembly_assembled_bcs_frac": empty_percent_metric,
        "vdj_assembly_umis_per_cell_distribution": (
            cr_report.NumpyHistogramMetric,
            {args: [100], prefixes: vdj_genes},
        ),
        "vdj_assembly_umis_per_cell_median": median_genes_metric,
        "vdj_assembly_mean_unassembled_umis_frac": empty_mean_metric,
        "vdj_assembly_contigs_per_bc": (
            cr_report.NumpyHistogramMetric,
            {args: [5], prefixes: vdj_genes},
        ),
        "vdj_assembly_contig_bc_frac": percent_genes_metric,
        "vdj_assembly_contig_pair_detected_bc_frac": percent_gene_pairs_metric,
        "vdj_assembly_pairing_efficiency": percent_gene_pairs_metric,
        "vdj_assembly_contig_full_len_bc_frac": percent_genes_metric,
        "vdj_assembly_contig_full_len_q40_bc_frac": percent_genes_metric,
        "vdj_assembly_contig_pair_full_len_bc_frac": percent_gene_pairs_metric,
        "vdj_assembly_contig_pair_full_len_q40_bc_frac": percent_gene_pairs_metric,
        "vdj_assembly_contig_frac": percent_genes_metric,
        "vdj_assembly_chimeric_contig_frac": empty_percent_metric,
        "vdj_assembly_unannotated_contig_frac": empty_percent_metric,
        "vdj_assembly_v_detected_contig_frac": percent_genes_metric,
        "vdj_assembly_vd_detected_contig_frac": percent_genes_metric,
        "vdj_assembly_vj_detected_contig_frac": percent_genes_metric,
        "vdj_assembly_full_len_contig_frac": percent_genes_metric,
        "vdj_assembly_median_contig_length": median_genes_metric,
        "vdj_assembly_contig_pair_productive_full_len_bc_count": (
            cr_report.CountMetric,
            vdj_gene_pairs_prefix,
        ),
        "vdj_assembly_contig_pair_productive_full_len_bc_frac": percent_gene_pairs_metric,
        "vdj_assembly_gt0prodcdr_contig_pair_productive_full_len_bc_frac": percent_gene_pairs_metric,
        "vdj_assembly_cdr_detected_bc_frac": percent_genes_metric,
        "vdj_assembly_prod_cdr_bc_frac": percent_genes_metric,
        "vdj_assembly_gt1cdr_cdrPosBc_frac": percent_genes_metric,
        "vdj_assembly_gt1prodCdr_cdrPosBc_frac": percent_genes_metric,
        "vdj_assembly_mult_cdrs_cdrPosBc_frac": percent_genes_metric,
        # VDJ contig filtering metrics
        "vdj_high_conf_prod_contig_frac": percent_genes_metric,
        # VDJ clonotype inference metrics
        "vdj_clonotype_count": (cr_report.CountMetric, vdj_canonical_and_clonotype_prefix),
        "vdj_clonotype_freq": topn_clonotype_metric,
        "vdj_clonotype_prop": topn_clonotype_metric,
        "vdj_clonotype_diversity": diversity_canonical_clono_metric,
        "vdj_paired_clonotype_diversity": diversity_canonical_clono_metric,
        "vdj_unassigned_clonotype_bc_frac": percent_clonotype_metric,
        "vdj_paired_clonotype_frac": percent_canonical_clono_metric,
        "vdj_paired_clonotype_bc_frac": percent_canonical_clono_metric,
        "cdrs_per_bc_histogram": (
            cr_report.NumpyHistogramMetric,
            {args: [10], prefixes: vdj_clonotype_types},
        ),
        "major_clonotype_bc_frac": percent_clonotype_metric,
        "vdj_clonotype_gt1_v_annotations_contig_frac": percent_clonotype_metric,
        "vdj_clonotype_gt1_j_annotations_contig_frac": percent_clonotype_metric,
        "vdj_clonotype_consensus_wrong_cdr_contig_frac": percent_clonotype_metric,
        # ALIGN_READS_TO_CONTIGS_METRICS
        # Read/contig agreement
        "vdj_contig_mapping_frac": empty_percent_metric,
        "vdj_contig_chimera_frac": empty_percent_metric,
        "vdj_contig_improper_pair_frac": empty_percent_metric,
        "vdj_contig_antisense_pair_frac": empty_percent_metric,
        "vdj_contig_single_read_mapping_frac": empty_percent_metric,
        "vdj_contig_mismatch_frac": empty_percent_metric,
        "vdj_contig_insertion_frac": empty_percent_metric,
        "vdj_contig_deletion_frac": empty_percent_metric,
        "vdj_contig_mismatch_histogram": mismatch_break_hist,
        "vdj_contig_insertion_histogram": mismatch_break_hist,
        "vdj_contig_deletion_histogram": mismatch_break_hist,
        # Coverage
        "vdj_contig_min_depth_histogram": depth_break_hist,
        "vdj_contig_min_q40_depth_histogram": depth_break_hist,
        "vdj_contig_mean_depth": empty_mean_metric,
        "vdj_contig_mean_q40_depth": empty_mean_metric,
        "vdj_contig_mean_depth_histogram": depth_break_hist,
        "vdj_contig_mean_q40_depth_histogram": depth_break_hist,
        # Overall coverage (heuristic for how coverage along contig piles up)
        "vdj_contig_depth_at_contig_percentiles": contig_length_hist,
        "vdj_contig_q40_depth_at_contig_percentiles": contig_length_hist,
        # Soft clipping
        "vdj_contig_soft_clipped_bases_r1_three_prime_mean": empty_mean_metric,
        "vdj_contig_soft_clipped_bases_r1_five_prime_mean": empty_mean_metric,
        "vdj_contig_soft_clipped_bases_r1_three_prime_histogram": softclip_break_hist,
        "vdj_contig_soft_clipped_bases_r1_five_prime_histogram": softclip_break_hist,
        "vdj_contig_soft_clipped_bases_r1_three_prime_frac": empty_percent_metric,
        "vdj_contig_soft_clipped_bases_r1_five_prime_frac": empty_percent_metric,
        "vdj_contig_soft_clipped_bases_r2_three_prime_mean": empty_mean_metric,
        "vdj_contig_soft_clipped_bases_r2_five_prime_mean": empty_mean_metric,
        "vdj_contig_soft_clipped_bases_r2_three_prime_histogram": softclip_break_hist,
        "vdj_contig_soft_clipped_bases_r2_five_prime_histogram": softclip_break_hist,
        "vdj_contig_soft_clipped_bases_r2_three_prime_frac": empty_percent_metric,
        "vdj_contig_soft_clipped_bases_r2_five_prime_frac": empty_percent_metric,
        "vdj_contig_median_insert_size": empty_median_metric,
        "vdj_contig_insert_size_histogram": (
            cr_report.HistogramMetric,
            {args: [cr_constants.INSERT_SIZE_CUTOFFS]},
        ),
        "vdj_contig_iqr_insert_size": (cr_report.IQRMetric, {}),
    }
    return {prefix + k: v for k, v in ret_dict.items()}


def is_barcode_corrected(raw_bc_seq: bytes, processed_bc_seq: bytes | None) -> bool:
    if processed_bc_seq is None:
        return False

    bc_seq, _ = cr_utils.split_barcode_seq(processed_bc_seq)
    return bc_seq != raw_bc_seq


class VdjReporter(cr_report.Reporter):
    def __init__(self, prefix="", **kwargs):
        kwargs.update({"metrics_dict": vdj_metrics_dict(prefix)})

        self.vdj_genes = [lib_constants.MULTI_REFS_PREFIX] + vdj_constants.VDJ_GENES
        self.canonical_vdj_genes = [
            lib_constants.MULTI_REFS_PREFIX
        ] + vdj_constants.CANONICAL_VDJ_GENES
        self.canonical_vdj_genes_nomulti = vdj_constants.CANONICAL_VDJ_GENES
        self.vdj_gene_pairs = [lib_constants.MULTI_REFS_PREFIX] + vdj_constants.VDJ_GENE_PAIRS
        self.canonical_vdj_gene_pairs = [
            lib_constants.MULTI_REFS_PREFIX
        ] + vdj_constants.CANONICAL_VDJ_GENE_PAIRS
        self.canonical_vdj_gene_pairs_nomulti = vdj_constants.CANONICAL_VDJ_GENE_PAIRS
        self.vdj_clonotype_types = vdj_constants.VDJ_CLONOTYPE_TYPES

        if "vdj_reference_path" in kwargs:
            self.vdj_reference_path = kwargs["vdj_reference_path"]
            self.vdj_reference = vdj_reference.VdjReference(kwargs["vdj_reference_path"])
            kwargs.pop("vdj_reference_path")

        self.primer_names_plus_rc = []
        if "primers" in kwargs:
            self.primer_names_plus_rc = [primer.name for primer in kwargs["primers"]] + [
                primer.name + "_rc" for primer in kwargs["primers"]
            ]

        cr_report.Reporter.__init__(self, **kwargs)

    def vdj_filter_barcodes_cb(
        self, cell_barcodes, barcodes, counts, total_read_pairs, recovered_cells
    ):
        self._get_metric_attr("vdj_filtered_bcs").set_value(len(cell_barcodes))
        cell_barcodes = set(cell_barcodes)
        cell_read_pairs = 0
        barcoded_read_pairs = 0
        for bc, count in zip(barcodes, counts):
            if bc in cell_barcodes:
                cell_read_pairs += count
            barcoded_read_pairs += count

        self._get_metric_attr("vdj_filtered_bcs_cum_frac").set_value(
            cell_read_pairs, barcoded_read_pairs
        )
        self._get_metric_attr("vdj_total_raw_read_pairs_per_filtered_bc").set_value(
            total_read_pairs, len(cell_barcodes)
        )

        self._get_metric_attr(
            "vdj_filtered_bcs_relative_difference_from_recovered_cells"
        ).set_value(len(cell_barcodes) - recovered_cells, recovered_cells)

    def vdj_assembly_cb(
        self,
        umi_summary_df: pd.DataFrame | None,
        annotation_dicts: Iterable[vdj_annot.AnnotationDict],
        reference: vdj_reference.VdjReference,
        prefix="",
    ):
        # build the chain:num_umis dict using contig_annotations_json
        numis_by_chain: dict[str, int] = {key: 0 for key in self.vdj_genes}
        for annotation_dict in annotation_dicts:
            annot = vdj_annot.AnnotatedContig.from_dict(annotation_dict, reference)
            if prefix == "":
                if annot.productive and annot.high_confidence and annot.is_single_chain():
                    chain = annot.contig_chains().pop()
                    numis_by_chain[ensure_str(chain)] += annot.umi_count
                    numis_by_chain[lib_constants.MULTI_REFS_PREFIX] += annot.umi_count
            elif annot.high_confidence and annot.is_single_chain():
                chain = annot.contig_chains().pop()
                numis_by_chain[ensure_str(chain)] += annot.umi_count
                numis_by_chain[lib_constants.MULTI_REFS_PREFIX] += annot.umi_count

        def exclusive_count_mismatch(chain):
            assert isinstance(chain, str)
            for chains in vdj_constants.EXCLUSIVE_VDJ_GENES:  # type: list[str]
                if chain in chains:
                    return numis_by_chain[chain] != max(numis_by_chain[c] for c in chains)
            return False

        for chain in self.vdj_genes:  # type: str
            assert isinstance(chain, str)
            # Don't count this cell if this chain is exceeded in UMIs by another chain
            #   in its first exclusive set. E.g., if IGK is present and IGL is absent
            #   then don't include the 0 in IGL's median.
            if exclusive_count_mismatch(chain):
                continue

            self._get_metric_attr(prefix + "vdj_assembly_umis_per_cell_median", chain).add(
                numis_by_chain[chain]
            )
            self._get_metric_attr(prefix + "vdj_assembly_umis_per_cell_distribution", chain).add(
                numis_by_chain[chain]
            )

        if umi_summary_df is None:
            frac_unassigned = 0.0
        else:
            good_umi_df = umi_summary_df[umi_summary_df["good_umi"] & (umi_summary_df["reads"] > 2)]
            frac_unassigned = np.mean(pd.isnull(good_umi_df.contigs))

        self._get_metric_attr(prefix + "vdj_assembly_mean_unassembled_umis_frac").add(
            frac_unassigned
        )

    def vdj_barcode_contig_cb(
        self,
        barcode,
        contigs: Collection[tuple[Any, bytes]],
        annotation_dicts: Collection[vdj_annot.AnnotationDict],
        reference: vdj_reference.VdjReference,
        prefix="",
    ):
        bc_fields = {"barcode": barcode}
        assert len(contigs) == len(annotation_dicts)

        # Count whether the assembly was succesfully generated
        self._get_metric_attr(prefix + "vdj_assembly_assembled_bcs_frac").add(
            1, filter=len(contigs) > 0
        )

        contig_type_counts = {}

        # Compute contig metrics
        full_len_genes = set()
        q40_genes = set()
        # Count productive (and full length) contigs by gene
        productive_contigs = defaultdict(int)
        # Distinct CDR sequences per gene
        gene_cdrs = defaultdict(set)
        gene_prod_cdrs = defaultdict(set)
        cdr_counts = {gene: defaultdict(int) for gene in self.vdj_genes}
        for gene in self.vdj_genes:
            assert isinstance(gene, str)

        for (_, contig_seq), annotation_dict in zip(contigs, annotation_dicts):
            annotation = vdj_annot.AnnotatedContig.from_dict(annotation_dict, reference)
            assert isinstance(contig_seq, bytes)
            if not annotation.high_confidence:
                continue

            # Assign a chain to this contig. Note: gene here is actually chain.
            v_hits = annotation.get_region_hits(vdj_constants.VDJ_V_FEATURE_TYPES)
            j_hits = annotation.get_region_hits(vdj_constants.VDJ_J_FEATURE_TYPES)
            contig_gene = annotation.get_single_chain()

            # First count unannotated vs annotated
            if contig_gene is None:
                # Unannotated
                self._get_metric_attr(prefix + "vdj_assembly_unannotated_contig_frac").add(
                    1, filter=True
                )
            else:
                contig_gene = ensure_str(contig_gene)
                assert isinstance(contig_gene, str)
                self._get_metric_attr(prefix + "vdj_assembly_unannotated_contig_frac").add(
                    1, filter=False
                )

            # Next count chimeras - overwrites contig_gene for use below
            chimeric_contig_frac = self._get_metric_attr(
                prefix + "vdj_assembly_chimeric_contig_frac"
            )
            if contig_gene == "Multi":
                # Chimeric
                chimeric_contig_frac.add(1, filter=True)
                contig_gene = None
            else:
                chimeric_contig_frac.add(1, filter=False)

            contig_type_counts[contig_gene] = 1 + contig_type_counts.get(contig_gene, 0)

            # Gene-specific per-contig metrics
            for gene in self.vdj_genes:
                is_gene = gene in (contig_gene, lib_constants.MULTI_REFS_PREFIX)

                # Contig length
                if is_gene:
                    self._get_metric_attr(prefix + "vdj_assembly_median_contig_length", gene).add(
                        len(contig_seq)
                    )

                # Contig VDJ detection
                self._get_metric_attr(prefix + "vdj_assembly_v_detected_contig_frac", gene).add(
                    1, filter=is_gene and len(v_hits) > 0
                )

                self._get_metric_attr(prefix + "vdj_assembly_vj_detected_contig_frac", gene).add(
                    1, filter=is_gene and len(v_hits) > 0 and len(j_hits) > 0
                )

                is_full_len_contig = annotation.has_full_length_vj_hit()

                self._get_metric_attr(prefix + "vdj_assembly_full_len_contig_frac", gene).add(
                    1, filter=is_gene and is_full_len_contig
                )

                if is_gene and gene != lib_constants.MULTI_REFS_PREFIX:
                    if is_full_len_contig:
                        full_len_genes.add(gene)
                        if annotation.productive:
                            productive_contigs[gene] += 1
                    if annotation.has_cdr():
                        gene_cdrs[gene].add(annotation.get_cdr_seq())
                        cdr_counts[gene][annotation.get_cdr_seq()] += 1
                        if annotation.productive:
                            gene_prod_cdrs[gene].add(annotation.get_cdr_seq())
                    quals = annotation.get_vj_quals()
                    if not quals is None and len(quals) > 0 and np.all(quals > 40):
                        q40_genes.add(gene)

        for gene in self.vdj_genes:
            if gene == lib_constants.MULTI_REFS_PREFIX:
                bc_gene_count = sum(contig_type_counts.get(g, 0) for g in self.vdj_genes)
                gene_is_full_len = len(full_len_genes) > 0
                gene_is_highQ = len(q40_genes) > 0
            else:
                bc_gene_count = contig_type_counts.get(gene, 0)
                gene_is_full_len = gene in full_len_genes
                gene_is_highQ = gene in q40_genes

            bc_fields[gene + "_count"] = bc_gene_count
            bc_fields[gene + "_full_len_count"] = sum(1 for g in full_len_genes if g == gene)
            bc_fields[gene + "_q40_count"] = sum(1 for g in full_len_genes if g == gene)

            self._get_metric_attr(prefix + "vdj_assembly_contigs_per_bc", gene).add(bc_gene_count)

            self._get_metric_attr(prefix + "vdj_assembly_contig_bc_frac", gene).add(
                1, filter=bc_gene_count > 0
            )

            self._get_metric_attr(prefix + "vdj_assembly_contig_full_len_bc_frac", gene).add(
                1, filter=gene_is_full_len
            )

            self._get_metric_attr(prefix + "vdj_assembly_contig_full_len_q40_bc_frac", gene).add(
                1, filter=gene_is_highQ
            )

            if gene != lib_constants.MULTI_REFS_PREFIX:
                self._get_metric_attr(prefix + "vdj_assembly_cdr_detected_bc_frac", gene).add(
                    1, filter=len(gene_cdrs[gene]) > 0
                )

                self._get_metric_attr(prefix + "vdj_assembly_prod_cdr_bc_frac", gene).add(
                    1, filter=len(gene_prod_cdrs[gene]) > 0
                )

                bc_fields[gene + "_cdr_count"] = len(gene_cdrs[gene])
                bc_fields[gene + "_prod_cdr_count"] = len(gene_prod_cdrs[gene])

                if len(gene_cdrs[gene]) > 0:
                    self._get_metric_attr(prefix + "vdj_assembly_gt1cdr_cdrPosBc_frac", gene).add(
                        1, filter=len(gene_cdrs[gene]) > 1
                    )

                    self._get_metric_attr(
                        prefix + "vdj_assembly_gt1prodCdr_cdrPosBc_frac", gene
                    ).add(1, filter=len(gene_prod_cdrs[gene]) > 1)

                    # Chech for multiple occurrences of the same CDR
                    mult_cdrs = np.any(np.array([cdr_counts[gene][g] for g in gene_cdrs[gene]]) > 1)
                    self._get_metric_attr(
                        prefix + "vdj_assembly_mult_cdrs_cdrPosBc_frac", gene
                    ).add(1, filter=mult_cdrs)

        pairs_detected = set()
        full_len_pairs = set()
        q40_pairs = set()
        productive_pairs = set()
        for gene_pair in self.vdj_gene_pairs:
            if all(
                contig_type_counts.get(gene, 0) > 0
                for gene in vdj_utils.get_genes_in_pair(gene_pair)
            ):
                pairs_detected.add(gene_pair)
            if gene_pair != lib_constants.MULTI_REFS_PREFIX and all(
                (gene in full_len_genes) for gene in vdj_utils.get_genes_in_pair(gene_pair)
            ):
                full_len_pairs.add(gene_pair)
            if gene_pair != lib_constants.MULTI_REFS_PREFIX and all(
                (productive_contigs[gene] > 0) for gene in vdj_utils.get_genes_in_pair(gene_pair)
            ):
                productive_pairs.add(gene_pair)
            if gene_pair != lib_constants.MULTI_REFS_PREFIX and all(
                gene in q40_genes for gene in vdj_utils.get_genes_in_pair(gene_pair)
            ):
                q40_pairs.add(gene_pair)

        for gene_pair in self.vdj_gene_pairs:
            self._get_metric_attr(
                prefix + "vdj_assembly_contig_pair_detected_bc_frac", gene_pair
            ).add(
                1,
                filter=gene_pair in pairs_detected
                or (gene_pair == lib_constants.MULTI_REFS_PREFIX and len(pairs_detected) > 0),
            )

            self._get_metric_attr(
                prefix + "vdj_assembly_contig_pair_full_len_bc_frac", gene_pair
            ).add(
                1,
                filter=gene_pair in full_len_pairs
                or (gene_pair == lib_constants.MULTI_REFS_PREFIX and len(full_len_pairs) > 0),
            )

            if gene_pair != lib_constants.MULTI_REFS_PREFIX:
                any_detected = any(
                    gene in full_len_genes for gene in vdj_utils.get_genes_in_pair(gene_pair)
                )
                # Fraction of cells with pair over cells with any gene of the pair
                self._get_metric_attr(prefix + "vdj_assembly_pairing_efficiency", gene_pair).add(
                    int(any_detected), filter=gene_pair in full_len_pairs
                )

            self._get_metric_attr(
                prefix + "vdj_assembly_contig_pair_full_len_q40_bc_frac", gene_pair
            ).add(
                1,
                filter=gene_pair in q40_pairs
                or (gene_pair == lib_constants.MULTI_REFS_PREFIX and len(q40_pairs) > 0),
            )

            prod_pair = gene_pair in productive_pairs or (
                gene_pair == lib_constants.MULTI_REFS_PREFIX and len(productive_pairs) > 0
            )
            if prod_pair:
                self._get_metric_attr(
                    prefix + "vdj_assembly_contig_pair_productive_full_len_bc_count", gene_pair
                ).add(1)

            self._get_metric_attr(
                prefix + "vdj_assembly_contig_pair_productive_full_len_bc_frac", gene_pair
            ).add(1, filter=prod_pair)

            # Frac cells paired, conditional on the presence of a single productive contig
            if len(productive_contigs) > 0:
                self._get_metric_attr(
                    prefix + "vdj_assembly_gt0prodcdr_contig_pair_productive_full_len_bc_frac",
                    gene_pair,
                ).add(1, filter=prod_pair)

        return bc_fields

    def contig_mapping_frac_statistics_cb(self, contig_read_iter, contig_length, strand):
        """Calculate basic mapping statistics on how well VDJ reads map back to contigs.

        Used in ALIGN_READS_TO_CONTIGS_PD

        Args:
            contig_read_iter: iterator of reads in aligned to a single contig
            strand (string): + or - to indicate library orientation (MRO argument strand, for example)
        """
        # Get all the metrics from reporter
        mapping_frac: cr_report.PercentMetric = self._get_metric_attr("vdj_contig_mapping_frac")
        chimera_frac: cr_report.PercentMetric = self._get_metric_attr("vdj_contig_chimera_frac")
        read_pair_improper_pair_frac: cr_report.PercentMetric = self._get_metric_attr(
            "vdj_contig_improper_pair_frac"
        )
        antisense_pair_frac: cr_report.PercentMetric = self._get_metric_attr(
            "vdj_contig_antisense_pair_frac"
        )
        single_read_mapped_frac: cr_report.PercentMetric = self._get_metric_attr(
            "vdj_contig_single_read_mapping_frac"
        )

        mismatch_histogram: cr_report.HistogramMetric = self._get_metric_attr(
            "vdj_contig_mismatch_histogram"
        )
        insertion_histogram: cr_report.HistogramMetric = self._get_metric_attr(
            "vdj_contig_insertion_histogram"
        )
        deletion_histogram: cr_report.HistogramMetric = self._get_metric_attr(
            "vdj_contig_deletion_histogram"
        )
        mismatch_frac: cr_report.PercentMetric = self._get_metric_attr("vdj_contig_mismatch_frac")
        insertion_frac: cr_report.PercentMetric = self._get_metric_attr("vdj_contig_insertion_frac")
        deletion_frac: cr_report.PercentMetric = self._get_metric_attr("vdj_contig_deletion_frac")

        min_depth_histogram: cr_report.HistogramMetric = self._get_metric_attr(
            "vdj_contig_min_depth_histogram"
        )
        min_q40_depth_histogram: cr_report.HistogramMetric = self._get_metric_attr(
            "vdj_contig_min_q40_depth_histogram"
        )
        mean_depth: cr_report.MeanMetric = self._get_metric_attr("vdj_contig_mean_depth")
        mean_q40_depth: cr_report.MeanMetric = self._get_metric_attr("vdj_contig_mean_q40_depth")
        mean_depth_histogram: cr_report.HistogramMetric = self._get_metric_attr(
            "vdj_contig_mean_depth_histogram"
        )
        mean_q40_depth_histogram: cr_report.HistogramMetric = self._get_metric_attr(
            "vdj_contig_mean_q40_depth_histogram"
        )

        depth_at_contig_percentiles: cr_report.HistogramMetric = self._get_metric_attr(
            "vdj_contig_depth_at_contig_percentiles"
        )
        q40_depth_at_contig_percentiles: cr_report.HistogramMetric = self._get_metric_attr(
            "vdj_contig_q40_depth_at_contig_percentiles"
        )

        soft_clipped_r1_three_prime_histogram: cr_report.HistogramMetric = self._get_metric_attr(
            "vdj_contig_soft_clipped_bases_r1_three_prime_histogram"
        )
        soft_clipped_r1_five_prime_histogram: cr_report.HistogramMetric = self._get_metric_attr(
            "vdj_contig_soft_clipped_bases_r1_five_prime_histogram"
        )
        soft_clipped_r1_three_prime_mean: cr_report.MeanMetric = self._get_metric_attr(
            "vdj_contig_soft_clipped_bases_r1_three_prime_mean"
        )
        soft_clipped_r1_five_prime_mean: cr_report.MeanMetric = self._get_metric_attr(
            "vdj_contig_soft_clipped_bases_r1_five_prime_mean"
        )
        soft_clipped_r1_three_prime_frac: cr_report.PercentMetric = self._get_metric_attr(
            "vdj_contig_soft_clipped_bases_r1_three_prime_frac"
        )
        soft_clipped_r1_five_prime_frac: cr_report.PercentMetric = self._get_metric_attr(
            "vdj_contig_soft_clipped_bases_r1_five_prime_frac"
        )
        soft_clipped_r2_three_prime_histogram: cr_report.HistogramMetric = self._get_metric_attr(
            "vdj_contig_soft_clipped_bases_r2_three_prime_histogram"
        )
        soft_clipped_r2_five_prime_histogram: cr_report.HistogramMetric = self._get_metric_attr(
            "vdj_contig_soft_clipped_bases_r2_five_prime_histogram"
        )
        soft_clipped_r2_three_prime_mean: cr_report.MeanMetric = self._get_metric_attr(
            "vdj_contig_soft_clipped_bases_r2_three_prime_mean"
        )
        soft_clipped_r2_five_prime_mean: cr_report.MeanMetric = self._get_metric_attr(
            "vdj_contig_soft_clipped_bases_r2_five_prime_mean"
        )
        soft_clipped_r2_three_prime_frac: cr_report.PercentMetric = self._get_metric_attr(
            "vdj_contig_soft_clipped_bases_r2_three_prime_frac"
        )
        soft_clipped_r2_five_prime_frac: cr_report.PercentMetric = self._get_metric_attr(
            "vdj_contig_soft_clipped_bases_r2_five_prime_frac"
        )
        median_insert_size: cr_report.MedianMetric = self._get_metric_attr(
            "vdj_contig_median_insert_size"
        )
        insert_size_histogram: cr_report.HistogramMetric = self._get_metric_attr(
            "vdj_contig_insert_size_histogram"
        )
        iqr_insert_size: cr_report.IQRMetric = self._get_metric_attr("vdj_contig_iqr_insert_size")

        # Calculate metrics over all reads
        contig_coverage = None
        contig_q40_coverage = None

        for read in contig_read_iter:
            # Initialize contig coverage to all zeroes across the contig
            if contig_coverage is None:
                contig_coverage = np.array([0] * contig_length)
                contig_q40_coverage = np.array([0] * contig_length)
            assert contig_q40_coverage is not None

            # Mapping metrics
            mapping_frac.add(1, filter=not read.is_unmapped)

            # Mismatch metrics and aggregation of coverage metrics (only mapped reads)
            if not read.is_unmapped:
                # Base category counts
                alignment_stats = cr_segment.get_cigar_summary_stats(read, strand)
                insertions = alignment_stats.get("I", 0)
                deletions = alignment_stats.get("D", 0)
                mismatches = dict(read.tags).get("NM") - insertions - deletions

                mismatch_histogram.add(mismatches)
                insertion_histogram.add(insertions)
                deletion_histogram.add(deletions)
                mismatch_frac.add(1, filter=mismatches > 0)
                insertion_frac.add(1, filter=insertions > 0)
                deletion_frac.add(1, filter=deletions > 0)

                # Insert size (only compute once per read-pair)
                if read.is_read1 and not read.mate_is_unmapped:
                    insert_size = np.abs(read.tlen)
                    if insert_size <= cr_constants.MAX_INSERT_SIZE:
                        median_insert_size.add(int(insert_size))
                        insert_size_histogram.add(int(insert_size))
                        iqr_insert_size.add(int(insert_size))

                # Soft clipping
                soft_clipped_r1_three_prime_mean.add(alignment_stats.get("R1_S_three_prime", 0))
                soft_clipped_r1_five_prime_mean.add(alignment_stats.get("R1_S_five_prime", 0))
                soft_clipped_r2_three_prime_mean.add(alignment_stats.get("R2_S_three_prime", 0))
                soft_clipped_r2_five_prime_mean.add(alignment_stats.get("R2_S_five_prime", 0))

                soft_clipped_r1_three_prime_histogram.add(
                    alignment_stats.get("R1_S_three_prime", 0)
                )
                soft_clipped_r1_five_prime_histogram.add(alignment_stats.get("R1_S_five_prime", 0))
                soft_clipped_r2_three_prime_histogram.add(
                    alignment_stats.get("R2_S_three_prime", 0)
                )
                soft_clipped_r2_five_prime_histogram.add(alignment_stats.get("R2_S_five_prime", 0))

                soft_clipped_r1_three_prime_frac.add(
                    1, alignment_stats.get("R1_S_three_prime", 0) > 0
                )
                soft_clipped_r1_five_prime_frac.add(
                    1, alignment_stats.get("R1_S_five_prime", 0) > 0
                )
                soft_clipped_r2_three_prime_frac.add(
                    1, alignment_stats.get("R2_S_three_prime", 0) > 0
                )
                soft_clipped_r2_five_prime_frac.add(
                    1, alignment_stats.get("R2_S_five_prime", 0) > 0
                )

                # Aggregate coverage metrics (track over all reads in a contig)
                (
                    alignment_length,
                    quality_scores,
                ) = cr_segment.get_full_alignment_base_quality_scores(read)
                q40_coverage = (quality_scores > 40).astype(int)
                all_coverage = (quality_scores > 0).astype(int)

                read_start = read.pos

                contig_q40_coverage[read_start : read_start + alignment_length] += q40_coverage
                contig_coverage[read_start : read_start + alignment_length] += all_coverage

            # Read pair metrics (only consider one read in pair so do not double count mate)
            if read.is_read1:
                chimera_frac.add(
                    1,
                    filter=read.rname != read.rnext
                    and not (read.is_unmapped or read.mate_is_unmapped),
                )
                read_pair_improper_pair_frac.add(1, filter=read.is_reverse == read.mate_is_reverse)
                antisense_pair_frac.add(
                    1,
                    filter=(
                        read.is_read1
                        and not read.is_reverse
                        and strand == cr_constants.REVERSE_STRAND
                    ),
                )
                single_read_mapped_frac.add(1, read.is_unmapped != read.mate_is_unmapped)

        if contig_coverage is None:
            return
        assert contig_q40_coverage is not None

        # Coverage (perform on aggregated values over contigs)
        min_q40_depth_histogram.add(contig_q40_coverage.min())
        min_depth_histogram.add(contig_coverage.min())

        mean_q40_depth.add(contig_q40_coverage.mean())
        mean_depth.add(contig_coverage.mean())

        mean_q40_depth_histogram.add(contig_q40_coverage.mean())
        mean_depth_histogram.add(contig_coverage.mean())

        # General pattern of coverage across length of contigs
        # Split each contig into the same number of bins. For each bin, compute the
        # average coverage, and sum the per-bin average coverage across contigs.
        for depth, q40_depth, percentile in zip(
            np.array_split(contig_coverage, len(VDJ_CONTIG_LENGTH_PERCENTILES)),
            np.array_split(contig_q40_coverage, len(VDJ_CONTIG_LENGTH_PERCENTILES)),
            VDJ_CONTIG_LENGTH_PERCENTILES,
        ):
            depth_at_contig_percentiles.add(percentile, int(depth.mean()))
            q40_depth_at_contig_percentiles.add(percentile, int(q40_depth.mean()))

    def vdj_barcode_cb(self, raw_bc, processed_bc):
        self._get_metric_attr("vdj_corrected_bc_frac").add(
            1, filter=is_barcode_corrected(raw_bc, processed_bc)
        )
        self._get_metric_attr("vdj_good_bc_frac").add(1, filter=processed_bc)

        if processed_bc:
            self._get_metric_attr("barcode_reads").add(processed_bc)
