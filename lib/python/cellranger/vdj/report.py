#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

from collections import defaultdict
import itertools
import numpy as np
import pandas as pd
import cellranger.constants as cr_constants
import cellranger.report as cr_report
import cellranger.utils as cr_utils
import cellranger.vdj.annotations as vdj_annot
import cellranger.vdj.constants as vdj_constants
import cellranger.library_constants as lib_constants
import cellranger.vdj.reference as vdj_reference
import cellranger.vdj.utils as vdj_utils


VDJ_METRICS = {
    'barcode_reads': (cr_report.DictionaryMetric, {'kwargs': {'report_type': 'barcodes'}, 'prefixes': ['references', 'regions', 'read_types']}),

    # Dummy metrics to add canonical VDJ gene prefixes
    'vdj_dummy_metric': (cr_report.CountMetric, {'prefixes': ['canonical_vdj_genes', 'canonical_vdj_gene_pairs', 'canonical_vdj_genes_nomulti', 'canonical_vdj_gene_pairs_nomulti']}),

    # VDJ general metrics
    'vdj_filtered_bcs': (cr_report.CountMetric, {}),
    'vdj_filtered_bcs_cum_frac': (cr_report.PercentMetric, {}),
    'vdj_filter_bcs_rpu_threshold': (cr_report.CountMetric, {'prefixes': ['gem_groups']}),
    'vdj_filter_bcs_umi_threshold': (cr_report.CountMetric, {'prefixes': ['gem_groups']}),
    'vdj_filter_bcs_confidence': (cr_report.CountMetric, {'prefixes': ['gem_groups']}),
    'vdj_corrected_bc_frac': (cr_report.PercentMetric, {}),
    'vdj_corrected_umi_frac': (cr_report.PercentMetric, {}),
    'vdj_good_bc_frac': (cr_report.PercentMetric, {}),
    'vdj_total_raw_reads_per_filtered_bc': (cr_report.PercentMetric, {}),
    'vdj_total_raw_read_pairs_per_filtered_bc': (cr_report.PercentMetric, {}),
    'vdj_assemblable_read_pairs_per_filtered_bc':(cr_report.PercentMetric, {}),
    'vdj_filtered_bcs_relative_difference_from_recovered_cells': (cr_report.PercentMetric, {}),
    'vdj_sequencing_efficiency': (cr_report.PercentMetric, {}),

    # VDJ read-filter metrics
    'vdj_recombinome_mapped_reads_frac': (cr_report.PercentMetric, {'prefixes': ['vdj_genes']}),
    'vdj_recombinome_mapped_read_pairs_frac': (cr_report.PercentMetric, {'prefixes': ['vdj_genes']}),
    'vdj_recombinome_antisense_read_pairs_frac': (cr_report.PercentMetric, {'prefixes': ['vdj_genes']}),
    'vdj_recombinome_mapped_read1_frac': (cr_report.PercentMetric, {'prefixes': ['vdj_genes']}),
    'vdj_recombinome_mapped_read2_frac': (cr_report.PercentMetric, {'prefixes': ['vdj_genes']}),
    'vdj_recombinome_chimeric_read_pairs_frac':(cr_report.PercentMetric, {'prefixes': ['vdj_genes']}),
    'vdj_recombinome_total_umis_per_cell_distribution': (cr_report.NumpyHistogramMetric, {'args': [100]}),
    'vdj_recombinome_total_umis_per_cell_median': (cr_report.MedianMetric, {}),
    'vdj_recombinome_umis_per_cell_distribution': (cr_report.NumpyHistogramMetric, {'args': [100], 'prefixes': ['vdj_genes']}),
    'vdj_recombinome_umis_per_cell_median': (cr_report.MedianMetric, {'prefixes': ['vdj_genes']}),

    # VDJ UMI filtering metrics
    'vdj_recombinome_readpairs_per_umi_distribution': (cr_report.DictionaryMetric, {'prefixes': ['vdj_genes']}),
    'vdj_recombinome_low_support_reads_frac': (cr_report.PercentMetric, {'prefixes': ['vdj_genes']}),
    'vdj_recombinome_readpairs_per_umi_n50': (cr_report.CountMetric, {'prefixes': ['vdj_genes', 'gem_groups']}),

    # VDJ primer trimming metrics
    'vdj_trim_read1_frac': (cr_report.PercentMetric, {'prefixes': ['primer_names_plus_rc']}),
    'vdj_trim_read2_frac': (cr_report.PercentMetric, {'prefixes': ['primer_names_plus_rc']}),

    # VDJ assembly metrics
    'vdj_assembly_assembled_bcs_frac': (cr_report.PercentMetric, {}),
    'vdj_assembly_umis_per_cell_distribution': (cr_report.NumpyHistogramMetric, {'args': [100], 'prefixes': ['vdj_genes']}),
    'vdj_assembly_umis_per_cell_median': (cr_report.MedianMetric, {'prefixes': ['vdj_genes']}),
    'vdj_assembly_mean_unassembled_umis_frac': (cr_report.MeanMetric, {}),
    'vdj_assembly_contigs_per_bc': (cr_report.NumpyHistogramMetric, {'args': [5], 'prefixes': ['vdj_genes']}),
    'vdj_assembly_contig_bc_frac': (cr_report.PercentMetric, {'prefixes': ['vdj_genes']}),
    'vdj_assembly_contig_pair_detected_bc_frac': (cr_report.PercentMetric, {'prefixes': ['vdj_gene_pairs']}),
    'vdj_assembly_pairing_efficiency': (cr_report.PercentMetric, {'prefixes': ['vdj_gene_pairs']}),
    'vdj_assembly_contig_full_len_bc_frac': (cr_report.PercentMetric, {'prefixes': ['vdj_genes']}),
    'vdj_assembly_contig_full_len_q40_bc_frac': (cr_report.PercentMetric, {'prefixes': ['vdj_genes']}),
    'vdj_assembly_contig_pair_full_len_bc_frac': (cr_report.PercentMetric, {'prefixes': ['vdj_gene_pairs']}),
    'vdj_assembly_contig_pair_full_len_q40_bc_frac': (cr_report.PercentMetric, {'prefixes': ['vdj_gene_pairs']}),

    'vdj_assembly_contig_frac': (cr_report.PercentMetric, {'prefixes': ['vdj_genes']}),
    'vdj_assembly_chimeric_contig_frac': (cr_report.PercentMetric, {}),
    'vdj_assembly_unannotated_contig_frac': (cr_report.PercentMetric, {}),
    'vdj_assembly_v_detected_contig_frac': (cr_report.PercentMetric, {'prefixes': ['vdj_genes']}),
    'vdj_assembly_vd_detected_contig_frac': (cr_report.PercentMetric, {'prefixes': ['vdj_genes']}),
    'vdj_assembly_vj_detected_contig_frac': (cr_report.PercentMetric, {'prefixes': ['vdj_genes']}),
    'vdj_assembly_full_len_contig_frac': (cr_report.PercentMetric, {'prefixes': ['vdj_genes']}),

    'vdj_assembly_median_contig_length': (cr_report.MedianMetric, {'prefixes': ['vdj_genes']}),

    'vdj_assembly_contig_pair_productive_full_len_bc_count': (cr_report.CountMetric, {'prefixes': ['vdj_gene_pairs']}),
    'vdj_assembly_contig_pair_productive_full_len_bc_frac': (cr_report.PercentMetric, {'prefixes': ['vdj_gene_pairs']}),
    'vdj_assembly_gt0prodcdr_contig_pair_productive_full_len_bc_frac': (cr_report.PercentMetric, {'prefixes': ['vdj_gene_pairs']}),
    'vdj_assembly_cdr_detected_bc_frac': (cr_report.PercentMetric, {'prefixes': ['vdj_genes']}),
    'vdj_assembly_prod_cdr_bc_frac': (cr_report.PercentMetric, {'prefixes': ['vdj_genes']}),
    'vdj_assembly_gt1cdr_cdrPosBc_frac': (cr_report.PercentMetric, {'prefixes': ['vdj_genes']}),
    'vdj_assembly_gt1prodCdr_cdrPosBc_frac': (cr_report.PercentMetric, {'prefixes': ['vdj_genes']}),
    'vdj_assembly_mult_cdrs_cdrPosBc_frac': (cr_report.PercentMetric, {'prefixes': ['vdj_genes']}),

    # VDJ contig filtering metrics
    'vdj_high_conf_prod_contig_frac': (cr_report.PercentMetric, {'prefixes': ['vdj_genes']}),

    # VDJ clonotype inference metrics
    'vdj_clonotype_count': (cr_report.CountMetric, {'prefixes': ['canonical_vdj_gene_pairs', 'vdj_clonotype_types']}),
    'vdj_clonotype_freq': (cr_report.TopNMetric, {'args': [vdj_constants.VDJ_MAX_OUTPUT_CLONOTYPES],
                                                  'prefixes': ['vdj_clonotype_types']}),
    'vdj_clonotype_prop': (cr_report.TopNMetric, {'args': [vdj_constants.VDJ_MAX_OUTPUT_CLONOTYPES],
                                                  'prefixes': ['vdj_clonotype_types']}),
    'vdj_clonotype_diversity': (cr_report.EffectiveDiversityMetric, {'prefixes': ['canonical_vdj_gene_pairs', 'vdj_clonotype_types']}),
    'vdj_paired_clonotype_diversity': (cr_report.EffectiveDiversityMetric, {'prefixes': ['canonical_vdj_gene_pairs', 'vdj_clonotype_types']}),
    'vdj_unassigned_clonotype_bc_frac': (cr_report.PercentMetric, {'prefixes': ['vdj_clonotype_types']}),
    'vdj_paired_clonotype_frac': (cr_report.PercentMetric, {'prefixes': ['canonical_vdj_gene_pairs', 'vdj_clonotype_types']}),
    'vdj_paired_clonotype_bc_frac': (cr_report.PercentMetric, {'prefixes': ['canonical_vdj_gene_pairs', 'vdj_clonotype_types']}),

    'cdrs_per_bc_histogram': (cr_report.NumpyHistogramMetric, {'args': [10], 'prefixes': ['vdj_clonotype_types']}),
    'major_clonotype_bc_frac': (cr_report.PercentMetric, {'prefixes': ['vdj_clonotype_types']}),
    'vdj_clonotype_gt1_v_annotations_contig_frac': (cr_report.PercentMetric, {'prefixes': ['vdj_clonotype_types']}),
    'vdj_clonotype_gt1_j_annotations_contig_frac': (cr_report.PercentMetric, {'prefixes': ['vdj_clonotype_types']}),
    'vdj_clonotype_consensus_wrong_cdr_contig_frac': (cr_report.PercentMetric, {'prefixes': ['vdj_clonotype_types']}),

    # ALIGN_READS_TO_CONTIGS_METRICS
    # Read/contig agreement
    'vdj_contig_mapping_frac': (cr_report.PercentMetric, {}),
    'vdj_contig_chimera_frac': (cr_report.PercentMetric, {}),
    'vdj_contig_improper_pair_frac': (cr_report.PercentMetric, {}),
    'vdj_contig_antisense_pair_frac': (cr_report.PercentMetric, {}),
    'vdj_contig_single_read_mapping_frac': (cr_report.PercentMetric, {}),
    'vdj_contig_mismatch_frac': (cr_report.PercentMetric, {}),
    'vdj_contig_insertion_frac': (cr_report.PercentMetric, {}),
    'vdj_contig_deletion_frac': (cr_report.PercentMetric, {}),
    'vdj_contig_mismatch_histogram': (cr_report.HistogramMetric, {'args': [vdj_constants.VDJ_MISMATCH_BREAKS]}),
    'vdj_contig_insertion_histogram': (cr_report.HistogramMetric, {'args': [vdj_constants.VDJ_MISMATCH_BREAKS]}),
    'vdj_contig_deletion_histogram': (cr_report.HistogramMetric, {'args': [vdj_constants.VDJ_MISMATCH_BREAKS]}),

    # Coverage
    'vdj_contig_min_depth_histogram': (cr_report.HistogramMetric, {'args': [vdj_constants.VDJ_DEPTH_BREAKS]}),
    'vdj_contig_min_q40_depth_histogram': (cr_report.HistogramMetric, {'args': [vdj_constants.VDJ_DEPTH_BREAKS]}),
    'vdj_contig_mean_depth': (cr_report.MeanMetric, {}),
    'vdj_contig_mean_q40_depth': (cr_report.MeanMetric, {}),
    'vdj_contig_mean_depth_histogram': (cr_report.HistogramMetric, {'args': [vdj_constants.VDJ_DEPTH_BREAKS]}),
    'vdj_contig_mean_q40_depth_histogram': (cr_report.HistogramMetric, {'args': [vdj_constants.VDJ_DEPTH_BREAKS]}),

    # Overall coverage (heuristic for how coverage along contig piles up)
    'vdj_contig_depth_at_contig_percentiles': (cr_report.HistogramMetric, {'args': [vdj_constants.VDJ_CONTIG_LENGTH_PERCENTILES]}),
    'vdj_contig_q40_depth_at_contig_percentiles': (cr_report.HistogramMetric, {'args': [vdj_constants.VDJ_CONTIG_LENGTH_PERCENTILES]}),

    # Soft clipping
    'vdj_contig_soft_clipped_bases_r1_three_prime_mean': (cr_report.MeanMetric, {}),
    'vdj_contig_soft_clipped_bases_r1_five_prime_mean': (cr_report.MeanMetric, {}),
    'vdj_contig_soft_clipped_bases_r1_three_prime_histogram': (cr_report.HistogramMetric, {'args': [vdj_constants.VDJ_SOFTCLIP_BREAKS]}),
    'vdj_contig_soft_clipped_bases_r1_five_prime_histogram': (cr_report.HistogramMetric, {'args': [vdj_constants.VDJ_SOFTCLIP_BREAKS]}),
    'vdj_contig_soft_clipped_bases_r1_three_prime_frac': (cr_report.PercentMetric, {}),
    'vdj_contig_soft_clipped_bases_r1_five_prime_frac': (cr_report.PercentMetric, {}),

    'vdj_contig_soft_clipped_bases_r2_three_prime_mean': (cr_report.MeanMetric, {}),
    'vdj_contig_soft_clipped_bases_r2_five_prime_mean': (cr_report.MeanMetric, {}),
    'vdj_contig_soft_clipped_bases_r2_three_prime_histogram': (cr_report.HistogramMetric, {'args': [vdj_constants.VDJ_SOFTCLIP_BREAKS]}),
    'vdj_contig_soft_clipped_bases_r2_five_prime_histogram': (cr_report.HistogramMetric, {'args': [vdj_constants.VDJ_SOFTCLIP_BREAKS]}),
    'vdj_contig_soft_clipped_bases_r2_three_prime_frac': (cr_report.PercentMetric, {}),
    'vdj_contig_soft_clipped_bases_r2_five_prime_frac': (cr_report.PercentMetric, {}),

    'vdj_contig_median_insert_size':    (cr_report.MedianMetric, {}),
    'vdj_contig_insert_size_histogram': (cr_report.HistogramMetric, {'args': [cr_constants.INSERT_SIZE_CUTOFFS]}),
    'vdj_contig_iqr_insert_size':       (cr_report.IQRMetric, {}),

}

class VdjReporter(cr_report.Reporter):
    def __init__(self, **kwargs):
        kwargs.update({'metrics_dict': VDJ_METRICS})

        self.vdj_genes = [lib_constants.MULTI_REFS_PREFIX] + vdj_constants.VDJ_GENES
        self.canonical_vdj_genes = [lib_constants.MULTI_REFS_PREFIX] + vdj_constants.CANONICAL_VDJ_GENES
        self.canonical_vdj_genes_nomulti = vdj_constants.CANONICAL_VDJ_GENES
        self.vdj_gene_pairs = [lib_constants.MULTI_REFS_PREFIX] + vdj_constants.VDJ_GENE_PAIRS
        self.canonical_vdj_gene_pairs = [lib_constants.MULTI_REFS_PREFIX] + vdj_constants.CANONICAL_VDJ_GENE_PAIRS
        self.canonical_vdj_gene_pairs_nomulti = vdj_constants.CANONICAL_VDJ_GENE_PAIRS
        self.vdj_clonotype_types = vdj_constants.VDJ_CLONOTYPE_TYPES

        if 'vdj_reference_path' in kwargs:
            self.vdj_reference_path = kwargs['vdj_reference_path']
            self.vdj_reference = vdj_reference.VdjReference(kwargs['vdj_reference_path'])
            kwargs.pop('vdj_reference_path')

        self.primer_names_plus_rc = []
        if 'primers' in kwargs:
            self.primer_names_plus_rc = [primer.name for primer in kwargs['primers']] + \
                                        [primer.name + '_rc' for primer in kwargs['primers']]

        cr_report.Reporter.__init__(self, **kwargs)

    def vdj_filter_barcodes_cb(self, cell_barcodes, barcodes, counts,
                               total_read_pairs,
                               recovered_cells):
        self._get_metric_attr('vdj_filtered_bcs').set_value(len(cell_barcodes))
        cell_barcodes = set(cell_barcodes)
        cell_read_pairs = 0
        barcoded_read_pairs = 0
        for bc, count in itertools.izip(barcodes, counts):
            if bc in cell_barcodes:
                cell_read_pairs += count
            barcoded_read_pairs += count

        self._get_metric_attr('vdj_filtered_bcs_cum_frac').set_value(cell_read_pairs, barcoded_read_pairs)
        self._get_metric_attr('vdj_total_raw_read_pairs_per_filtered_bc').set_value(total_read_pairs, len(cell_barcodes))

        self._get_metric_attr('vdj_filtered_bcs_relative_difference_from_recovered_cells').set_value(len(cell_barcodes) - recovered_cells, recovered_cells)

    def vdj_recombinome_bam_cb(self, read1, read2, bam, strand):
        """ Report on a single read-to-germline-segment alignment.
            Returns (chain1, chain2); chain annotations for each read """
        paired_end = read2 is not None

        if paired_end:
            assert read1.is_secondary == read2.is_secondary
            if read1.is_secondary or read2.is_secondary:
                return (None, None)
        else:
            if read1.is_secondary:
                return (None, None)

        read1_mapped = not read1.is_unmapped
        read2_mapped = paired_end and not read2.is_unmapped
        pair_mapped = (read1_mapped and read2_mapped) or (not paired_end) and read1_mapped

        if read1_mapped:
            ref_name = bam.references[read1.tid]
            feature_id = vdj_reference.get_feature_id_from_aligned_ref_name(ref_name)
            read1_chain = self.vdj_reference.get_feature_by_id(feature_id).chain
        else:
            read1_chain = None

        if read2_mapped:
            ref_name = bam.references[read2.tid]
            feature_id = vdj_reference.get_feature_id_from_aligned_ref_name(ref_name)
            read2_chain = self.vdj_reference.get_feature_by_id(feature_id).chain
        else:
            read2_chain = None

        for chain in vdj_constants.VDJ_GENES + [lib_constants.MULTI_REFS_PREFIX]:
            if chain == lib_constants.MULTI_REFS_PREFIX:
                read1_this_chain = True
                read2_this_chain = True
            else:
                read1_this_chain = read1_chain == chain
                read2_this_chain = read2_chain == chain

            mapped_reads_frac = self._get_metric_attr('vdj_recombinome_mapped_reads_frac', chain)
            mapped_read1_frac = self._get_metric_attr('vdj_recombinome_mapped_read1_frac', chain)
            mapped_read2_frac = self._get_metric_attr('vdj_recombinome_mapped_read2_frac', chain)
            mapped_read_pairs_frac = self._get_metric_attr('vdj_recombinome_mapped_read_pairs_frac', chain)
            chimeric_read_pairs_frac = self._get_metric_attr('vdj_recombinome_chimeric_read_pairs_frac', chain)
            antisense_read_pairs_frac = self._get_metric_attr('vdj_recombinome_antisense_read_pairs_frac', chain)

            mapped_reads_frac.add(1, filter=read1_mapped and read1_this_chain)
            mapped_reads_frac.add(1, filter=read2_mapped and read2_this_chain)
            mapped_read1_frac.add(1, filter=read1_mapped and read1_this_chain)
            mapped_read2_frac.add(1, filter=read2_mapped and read2_this_chain)

            is_pair_mapped = pair_mapped and read1_this_chain and read2_this_chain \
                             or (not paired_end) and pair_mapped and read1_this_chain
            mapped_read_pairs_frac.add(1, filter=is_pair_mapped)

            chimeric_read_pairs_frac.add(1, filter=paired_end \
                                         and pair_mapped \
                                         and read1_chain != read2_chain \
                                         and (read1_this_chain or read2_this_chain))

            is_antisense = is_pair_mapped \
                           and ((strand == '+' and read1.is_reverse and (not read2.is_reverse)) \
                                or (strand == '-' and (not read1.is_reverse) and read2.is_reverse)) \
                           or (not paired_end) and ((strand == '+' and read1.is_reverse) \
                                                    or (strand == '-' and not read1.is_reverse))
            antisense_read_pairs_frac.add(1, filter=is_antisense)

        return (read1_chain, read2_chain)

    def vdj_assembly_cb(self, summary_df, umi_summary_df, annotation_dicts, reference):

        # From contig name to its chain (eg. TRA, TRB)
        contig_chains = {}
        for annotation_dict in annotation_dicts:
            annotation = vdj_annot.AnnotatedContig.from_dict(annotation_dict, reference)
            if not annotation.high_confidence:
                continue

            if annotation.is_single_chain():
                contig_chain = list(annotation.contig_chains())[0]
                contig_chains[annotation.contig_name] = contig_chain

        good_umis = set() if umi_summary_df is None else \
                    set([str(int(s)) for s in umi_summary_df[umi_summary_df['good_umi']]['umi_id']])

        # From chain to the umis assigned to contigs of that chain
        chain_umis = defaultdict(list)
        if not summary_df is None:
            for _, row in summary_df.iterrows():
                if row.contig_name in contig_chains:
                    # If the contig was assigned to a single chain...
                    chain = contig_chains[row.contig_name]
                    # ...get all umis for that contig
                    umis = [u for u in row.umi_list.split(',') if u in good_umis]
                    chain_umis[chain].extend(umis)
                    chain_umis[lib_constants.MULTI_REFS_PREFIX].extend(umis)

        numis_by_chain = {c: len(set(chain_umis.get(c, []))) for c in self.vdj_genes}

        for chain in self.vdj_genes:
            # Don't count this cell if this chain is exceeded in UMIs by another chain
            #   in its first exclusive set. E.g., if IGK is present and IGL is absent
            #   then don't include the 0 in IGL's median.
            exclusive_sets = filter(lambda s: chain in s, vdj_constants.EXCLUSIVE_VDJ_GENES)
            if len(exclusive_sets) > 0:
                counts = [numis_by_chain[c] for c in exclusive_sets[0]]

                if numis_by_chain[chain] != max(counts):
                    continue

            metric = self._get_metric_attr('vdj_assembly_umis_per_cell_median', chain)
            metric.add(numis_by_chain[chain])
            metric = self._get_metric_attr('vdj_assembly_umis_per_cell_distribution', chain)
            metric.add(numis_by_chain[chain])

        if umi_summary_df is None:
            frac_unassigned = 0.0
        else:
            good_umi_df = umi_summary_df[umi_summary_df['good_umi'] & (umi_summary_df['reads'] > 2)]
            frac_unassigned = np.mean(pd.isnull(good_umi_df.contigs))

        metric = self._get_metric_attr('vdj_assembly_mean_unassembled_umis_frac')
        metric.add(frac_unassigned)


    def vdj_barcode_contig_cb(self, barcode, contigs, annotation_dicts, reference):
        bc_fields = {}
        bc_fields['barcode'] = barcode
        assert len(contigs) == len(annotation_dicts)

        annotations = []
        for (_, contig_seq), annotation_dict in itertools.izip(contigs, annotation_dicts):
            annotations.append(vdj_annot.AnnotatedContig.from_dict(annotation_dict, reference))

        # Count whether the assembly was succesfully generated
        self._get_metric_attr('vdj_assembly_assembled_bcs_frac').add(1, filter=len(contigs) > 0)

        contig_type_counts = {}

        # Compute contig metrics
        full_len_genes = set()
        q40_genes = set()
        # Count productive (and full length) contigs by gene
        productive_contigs = defaultdict(int)
        # Distinct CDR sequences per gene
        gene_cdrs = defaultdict(set)
        gene_prod_cdrs = defaultdict(set)
        cdr_counts = {}
        for gene in self.vdj_genes:
            cdr_counts[gene] = defaultdict(int)

        for (contig_name, contig_seq), annotation in itertools.izip(contigs, annotations):
            if not annotation.high_confidence:
                continue

            # Assign a chain to this contig. Note: gene here is actually chain.
            v_hits = annotation.get_region_hits(vdj_constants.VDJ_V_FEATURE_TYPES)
            j_hits = annotation.get_region_hits(vdj_constants.VDJ_J_FEATURE_TYPES)
            contig_gene = annotation.get_single_chain()

            # First count unannotated vs annotated
            if contig_gene is None:
                # Unannotated
                self._get_metric_attr('vdj_assembly_unannotated_contig_frac').add(1, filter=True)
            else:
                self._get_metric_attr('vdj_assembly_unannotated_contig_frac').add(1, filter=False)

            # Next count chimeras - overwrites contig_gene for use below
            chimeric_contig_frac = self._get_metric_attr('vdj_assembly_chimeric_contig_frac')
            if contig_gene == 'Multi':
                # Chimeric
                chimeric_contig_frac.add(1, filter=True)
                contig_gene =  None
            else:
                chimeric_contig_frac.add(1, filter=False)


            contig_type_counts[contig_gene] = 1 + contig_type_counts.get(contig_gene, 0)

            # Gene-specific per-contig metrics
            _, _, assembly_gene, _ = vdj_utils.parse_contig_name(contig_name)
            for gene in self.vdj_genes:
                is_gene = gene == contig_gene or gene == lib_constants.MULTI_REFS_PREFIX

                # Contig length
                if is_gene:
                    contigs_median_len = self._get_metric_attr('vdj_assembly_median_contig_length', gene)
                    contigs_median_len.add(len(contig_seq))

                # Contig VDJ detection
                v_contig_frac = self._get_metric_attr('vdj_assembly_v_detected_contig_frac', gene)
                v_contig_frac.add(1, filter=is_gene and len(v_hits) > 0)

                vj_contig_frac = self._get_metric_attr('vdj_assembly_vj_detected_contig_frac', gene)
                vj_contig_frac.add(1, filter=is_gene and len(v_hits) > 0 and len(j_hits) > 0)

                is_full_len_contig = annotation.has_full_length_vj_hit()

                full_len_contig_frac = self._get_metric_attr('vdj_assembly_full_len_contig_frac', gene)
                full_len_contig_frac.add(1, filter=is_gene and is_full_len_contig)

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

            contigs_per_bc = self._get_metric_attr('vdj_assembly_contigs_per_bc', gene)
            contigs_per_bc.add(bc_gene_count)

            contig_bc_frac = self._get_metric_attr('vdj_assembly_contig_bc_frac', gene)
            contig_bc_frac.add(1, filter=bc_gene_count > 0)

            full_len_bc_frac = self._get_metric_attr('vdj_assembly_contig_full_len_bc_frac', gene)
            full_len_bc_frac.add(1, filter=gene_is_full_len)

            highQ_frac = self._get_metric_attr('vdj_assembly_contig_full_len_q40_bc_frac', gene)
            highQ_frac.add(1, filter=gene_is_highQ)

            if gene != lib_constants.MULTI_REFS_PREFIX:
                metric = self._get_metric_attr('vdj_assembly_cdr_detected_bc_frac', gene)
                metric.add(1, filter=len(gene_cdrs[gene]) > 0)

                metric = self._get_metric_attr('vdj_assembly_prod_cdr_bc_frac', gene)
                metric.add(1, filter=len(gene_prod_cdrs[gene]) > 0)

                bc_fields[gene + "_cdr_count"] = len(gene_cdrs[gene])
                bc_fields[gene + "_prod_cdr_count"] = len(gene_prod_cdrs[gene])

                if len(gene_cdrs[gene]) > 0:
                    metric = self._get_metric_attr('vdj_assembly_gt1cdr_cdrPosBc_frac', gene)
                    metric.add(1, filter=len(gene_cdrs[gene]) > 1)

                    metric = self._get_metric_attr('vdj_assembly_gt1prodCdr_cdrPosBc_frac', gene)
                    metric.add(1, filter=len(gene_prod_cdrs[gene]) > 1)

                    # Chech for multiple occurrences of the same CDR
                    mult_cdrs = np.any(np.array([cdr_counts[gene][g] for g in gene_cdrs[gene]]) > 1)
                    metric = self._get_metric_attr('vdj_assembly_mult_cdrs_cdrPosBc_frac', gene)
                    metric.add(1, filter=mult_cdrs)

        pairs_detected = set()
        full_len_pairs = set()
        q40_pairs = set()
        productive_pairs = set()
        for gene_pair in self.vdj_gene_pairs:
            if all(contig_type_counts.get(gene, 0) > 0 for gene in vdj_utils.get_genes_in_pair(gene_pair)):
                pairs_detected.add(gene_pair)
            if gene_pair != lib_constants.MULTI_REFS_PREFIX and \
               all((gene in full_len_genes) for gene in vdj_utils.get_genes_in_pair(gene_pair)):
                full_len_pairs.add(gene_pair)
            if gene_pair != lib_constants.MULTI_REFS_PREFIX and \
               all((productive_contigs[gene] > 0) for gene in vdj_utils.get_genes_in_pair(gene_pair)):
                productive_pairs.add(gene_pair)
            if gene_pair != lib_constants.MULTI_REFS_PREFIX and \
               all(gene in q40_genes for gene in vdj_utils.get_genes_in_pair(gene_pair)):
                q40_pairs.add(gene_pair)

        for gene_pair in self.vdj_gene_pairs:
            pair_detected_bc_frac = self._get_metric_attr('vdj_assembly_contig_pair_detected_bc_frac', gene_pair)
            pair_detected_bc_frac.add(1, filter=gene_pair in pairs_detected or (gene_pair == lib_constants.MULTI_REFS_PREFIX and len(pairs_detected) > 0))

            pair_full_len_bc_frac = self._get_metric_attr('vdj_assembly_contig_pair_full_len_bc_frac', gene_pair)
            pair_full_len_bc_frac.add(1, filter=gene_pair in full_len_pairs or (gene_pair == lib_constants.MULTI_REFS_PREFIX and len(full_len_pairs) > 0))

            if gene_pair != lib_constants.MULTI_REFS_PREFIX:
                efficiency = self._get_metric_attr('vdj_assembly_pairing_efficiency', gene_pair)
                any_detected = any(gene in full_len_genes for gene in vdj_utils.get_genes_in_pair(gene_pair))
                # Fraction of cells with pair over cells with any gene of the pair
                efficiency.add(int(any_detected), filter=gene_pair in full_len_pairs)

            pair_full_len_q40_bc_frac = self._get_metric_attr('vdj_assembly_contig_pair_full_len_q40_bc_frac', gene_pair)
            pair_full_len_q40_bc_frac.add(1, filter=gene_pair in q40_pairs or (gene_pair == lib_constants.MULTI_REFS_PREFIX and len(q40_pairs) > 0))

            prod_pair = gene_pair in productive_pairs or (gene_pair == lib_constants.MULTI_REFS_PREFIX and len(productive_pairs) > 0)
            if prod_pair:
                prod_pair_full_len_bc_count = self._get_metric_attr('vdj_assembly_contig_pair_productive_full_len_bc_count', gene_pair)
                prod_pair_full_len_bc_count.add(1)

            prod_pair_full_len_bc_frac = self._get_metric_attr('vdj_assembly_contig_pair_productive_full_len_bc_frac', gene_pair)
            prod_pair_full_len_bc_frac.add(1, filter=prod_pair)

            # Frac cells paired, conditional on the presence of a single productive contig
            if len(productive_contigs) > 0:
                hasprod_prod_pair_full_len_bc_frac = self._get_metric_attr('vdj_assembly_gt0prodcdr_contig_pair_productive_full_len_bc_frac', gene_pair)
                hasprod_prod_pair_full_len_bc_frac.add(1, filter=prod_pair)

        return bc_fields

    def contig_mapping_frac_statistics_cb(self, contig_read_iter, contig_length, strand):
        """
        Calculate basic mapping statistics on how well VDJ reads map back to the contigs they were assembled from.
        Used in ALIGN_READS_TO_CONTIGS_PD

        Args:
            contig_read_iter: iterator of reads in aligned to a single contig
            strand (string): + or - to indicate library orientation (MRO argument strand, for example)

        """

        # Get all the metrics from reporter
        mapping_frac = self._get_metric_attr('vdj_contig_mapping_frac')
        chimera_frac = self._get_metric_attr('vdj_contig_chimera_frac')
        read_pair_improper_pair_frac = self._get_metric_attr('vdj_contig_improper_pair_frac')
        antisense_pair_frac = self._get_metric_attr('vdj_contig_antisense_pair_frac')
        single_read_mapped_frac = self._get_metric_attr('vdj_contig_single_read_mapping_frac')

        mismatch_histogram = self._get_metric_attr('vdj_contig_mismatch_histogram')
        insertion_histogram = self._get_metric_attr('vdj_contig_insertion_histogram')
        deletion_histogram = self._get_metric_attr('vdj_contig_deletion_histogram')
        mismatch_frac = self._get_metric_attr('vdj_contig_mismatch_frac')
        insertion_frac = self._get_metric_attr('vdj_contig_insertion_frac')
        deletion_frac = self._get_metric_attr('vdj_contig_deletion_frac')

        min_depth_histogram = self._get_metric_attr('vdj_contig_min_depth_histogram')
        min_q40_depth_histogram = self._get_metric_attr('vdj_contig_min_q40_depth_histogram')
        mean_depth = self._get_metric_attr('vdj_contig_mean_depth')
        mean_q40_depth = self._get_metric_attr('vdj_contig_mean_q40_depth')
        mean_depth_histogram = self._get_metric_attr('vdj_contig_mean_depth_histogram')
        mean_q40_depth_histogram = self._get_metric_attr('vdj_contig_mean_q40_depth_histogram')

        depth_at_contig_percentiles = self._get_metric_attr('vdj_contig_depth_at_contig_percentiles')
        q40_depth_at_contig_percentiles = self._get_metric_attr('vdj_contig_q40_depth_at_contig_percentiles')

        soft_clipped_r1_three_prime_histogram = self._get_metric_attr('vdj_contig_soft_clipped_bases_r1_three_prime_histogram')
        soft_clipped_r1_five_prime_histogram = self._get_metric_attr('vdj_contig_soft_clipped_bases_r1_five_prime_histogram')
        soft_clipped_r1_three_prime_mean = self._get_metric_attr('vdj_contig_soft_clipped_bases_r1_three_prime_mean')
        soft_clipped_r1_five_prime_mean = self._get_metric_attr('vdj_contig_soft_clipped_bases_r1_five_prime_mean')
        soft_clipped_r1_three_prime_frac = self._get_metric_attr('vdj_contig_soft_clipped_bases_r1_three_prime_frac')
        soft_clipped_r1_five_prime_frac = self._get_metric_attr('vdj_contig_soft_clipped_bases_r1_five_prime_frac')
        soft_clipped_r2_three_prime_histogram = self._get_metric_attr('vdj_contig_soft_clipped_bases_r2_three_prime_histogram')
        soft_clipped_r2_five_prime_histogram = self._get_metric_attr('vdj_contig_soft_clipped_bases_r2_five_prime_histogram')
        soft_clipped_r2_three_prime_mean = self._get_metric_attr('vdj_contig_soft_clipped_bases_r2_three_prime_mean')
        soft_clipped_r2_five_prime_mean = self._get_metric_attr('vdj_contig_soft_clipped_bases_r2_five_prime_mean')
        soft_clipped_r2_three_prime_frac = self._get_metric_attr('vdj_contig_soft_clipped_bases_r2_three_prime_frac')
        soft_clipped_r2_five_prime_frac = self._get_metric_attr('vdj_contig_soft_clipped_bases_r2_five_prime_frac')
        median_insert_size = self._get_metric_attr('vdj_contig_median_insert_size')
        insert_size_histogram = self._get_metric_attr('vdj_contig_insert_size_histogram')
        iqr_insert_size = self._get_metric_attr('vdj_contig_iqr_insert_size')

        # Calculate metrics over all reads
        contig_coverage = None
        contig_q40_coverage = None

        for read in contig_read_iter:
            # Initialize contig coverage to all zeroes across the contig
            if contig_coverage is None:
                contig_coverage = np.array([0] * contig_length)
                contig_q40_coverage = np.array([0] * contig_length)

            # Mapping metrics
            mapping_frac.add(1, filter=not read.is_unmapped)

            # Mismatch metrics and aggregation of coverage metrics (only mapped reads)
            if not read.is_unmapped:
                # Base category counts
                alignment_stats = cr_utils.get_cigar_summary_stats(read, strand)
                insertions = alignment_stats.get('I', 0)
                deletions = alignment_stats.get('D', 0)
                mismatches = dict(read.tags).get('NM') - insertions - deletions

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
                        median_insert_size.add(insert_size)
                        insert_size_histogram.add(insert_size)
                        iqr_insert_size.add(insert_size)

                # Soft clipping
                soft_clipped_r1_three_prime_mean.add(alignment_stats.get('R1_S_three_prime', 0))
                soft_clipped_r1_five_prime_mean.add(alignment_stats.get('R1_S_five_prime', 0))
                soft_clipped_r2_three_prime_mean.add(alignment_stats.get('R2_S_three_prime', 0))
                soft_clipped_r2_five_prime_mean.add(alignment_stats.get('R2_S_five_prime', 0))

                soft_clipped_r1_three_prime_histogram.add(alignment_stats.get('R1_S_three_prime', 0))
                soft_clipped_r1_five_prime_histogram.add(alignment_stats.get('R1_S_five_prime', 0))
                soft_clipped_r2_three_prime_histogram.add(alignment_stats.get('R2_S_three_prime', 0))
                soft_clipped_r2_five_prime_histogram.add(alignment_stats.get('R2_S_five_prime', 0))

                soft_clipped_r1_three_prime_frac.add(1, alignment_stats.get('R1_S_three_prime', 0) > 0)
                soft_clipped_r1_five_prime_frac.add(1, alignment_stats.get('R1_S_five_prime', 0) > 0)
                soft_clipped_r2_three_prime_frac.add(1, alignment_stats.get('R2_S_three_prime', 0) > 0)
                soft_clipped_r2_five_prime_frac.add(1, alignment_stats.get('R2_S_five_prime', 0) > 0)

                # Aggregate coverage metrics (track over all reads in a contig)
                alignment_length, quality_scores = cr_utils.get_full_alignment_base_quality_scores(read)
                q40_coverage = (quality_scores > 40).astype(int)
                all_coverage = (quality_scores > 0).astype(int)

                read_start = read.pos

                contig_q40_coverage[read_start:read_start + alignment_length] += q40_coverage
                contig_coverage[read_start:read_start + alignment_length] += all_coverage

            # Read pair metrics (only consider one read in pair so do not double count mate)
            if read.is_read1:
                chimera_frac.add(1, filter=read.rname != read.rnext and not (read.is_unmapped or read.mate_is_unmapped))
                read_pair_improper_pair_frac.add(1, filter=read.is_reverse == read.mate_is_reverse)
                antisense_pair_frac.add(1, filter=(read.is_read1 and not read.is_reverse and strand == cr_constants.REVERSE_STRAND))
                single_read_mapped_frac.add(1, read.is_unmapped != read.mate_is_unmapped)

        if contig_coverage is None:
            return

        # Coverage (perform on aggregated values over contigs)
        min_q40_depth_histogram.add(contig_q40_coverage.min())
        min_depth_histogram.add(contig_coverage.min())

        mean_q40_depth.add(contig_q40_coverage.mean())
        mean_depth.add(contig_coverage.mean())

        mean_q40_depth_histogram.add(contig_q40_coverage.mean())
        mean_depth_histogram.add(contig_coverage.mean())

        # General pattern of coverage across length of contigs
        # Split each contig into the same number of bins. For each bin, compute the
        # average coverage.
        depth_percentile_binned = [int(x.mean()) for x in
                                   list(np.array_split(contig_coverage, len(vdj_constants.VDJ_CONTIG_LENGTH_PERCENTILES)))]
        q40_depth_percentile_binned = [int(x.mean()) for x in
                                       list(np.array_split(contig_q40_coverage, len(vdj_constants.VDJ_CONTIG_LENGTH_PERCENTILES)))]

        # Now sum the per-bin average coverage across contigs.
        for i, (depth, q40_depth) in enumerate(zip(depth_percentile_binned, q40_depth_percentile_binned)):
            depth_at_contig_percentiles.add(vdj_constants.VDJ_CONTIG_LENGTH_PERCENTILES[i], depth)
            q40_depth_at_contig_percentiles.add(vdj_constants.VDJ_CONTIG_LENGTH_PERCENTILES[i], q40_depth)


    def vdj_barcode_cb(self, raw_bc, processed_bc):
        self._get_metric_attr('vdj_corrected_bc_frac').add(1, filter=cr_utils.is_barcode_corrected(raw_bc, processed_bc))
        self._get_metric_attr('vdj_good_bc_frac').add(1, filter=processed_bc)

        if processed_bc:
            self._get_metric_attr('barcode_reads').add(processed_bc)
