#!/usr/bin/env python3
#
# Copyright (c) 2025 10X Genomics, Inc. All rights reserved.
#
"""Martian Python types."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, TypedDict

if TYPE_CHECKING:
    from cellranger.mro_types.filetypes import (
        bam,
        bam_bai,
        bam_csi,
        bed,
        csv,
        fa,
        fasta,
        fasta_fai,
        fastq,
        geojson,
        h5,
        json,
        npy,
        path,
        pb,
        tiff,
        tsv,
        vcf,
        vloupe_,
    )


class WhitelistSpec(TypedDict):
    """Different ways to specify a whitelist of sequences (typically barcode sequence)."""

    name: str
    translation: bool
    strand: str
    translation_whitelist_path: str
    slide: str
    part: str


class UmiReadComponent(TypedDict):
    """UMI read components."""

    read_type: str
    offset: int
    length: int  # The length of the UMI. At most this number of bases will be extracted for use as a UMI.
    min_length: int  # If a shorter UMI can be used, add it here. None indicates that the full


class RnaReadComponent(TypedDict):
    """RNA read components."""

    read_type: str
    offset: int
    length: int
    min_length: int


class BarcodeReadComponent(TypedDict):
    """Barcode read components."""

    read_type: str
    kind: str
    offset: int
    length: int
    whitelist: WhitelistSpec


class ChemistryDef(TypedDict):
    """Chemistry definition."""

    name: str
    description: str
    endedness: str
    strandedness: str
    barcode: list[BarcodeReadComponent]
    umi: list[UmiReadComponent]
    rna: RnaReadComponent
    rna2: RnaReadComponent
    barcode_extraction: dict


class TranscriptomeFileInfo(TypedDict):
    """Transcriptome files information extracted from <reference_path>/reference.json."""

    fasta_hash: str
    gtf_hash: str
    gtf_hash_gz: str
    input_fasta_files: list[str]
    input_gtf_files: list[str]
    mem_gb: float
    mkref_version: str
    reference_path: path


class ReferenceInfo(TypedDict):
    """Reference genome/transcriptome information."""

    genomes: list[str]
    version: str
    transcriptome_info: TranscriptomeFileInfo | None


class ProbeBCDef(TypedDict):
    """Probe barcode definition."""

    id: str
    sequence: list[str]
    offset: int
    length: int


class SubsamplingDef(TypedDict):
    """The dictionary of metadata returned by make_subsamplings."""

    library_type: str
    subsample_type: str
    target_read_pairs_per_cell: int
    library_subsample_rates: list[float]


class CellCallingParam(TypedDict):
    """Cell calling parameter."""

    per_gem_well: float
    per_sample: dict[str, float]


class CellCalling(TypedDict):
    """Carries options to customize the cell calling mode."""

    recovered_cells: CellCallingParam
    force_cells: CellCallingParam
    emptydrops_minimum_umis: CellCallingParam
    global_minimum_umis: CellCallingParam
    max_mito_percent: CellCallingParam
    cell_barcodes: json
    override_mode: str
    override_library_types: list[str]
    disable_ab_aggregate_detection: bool
    disable_high_occupancy_gem_detection: bool


# mro/rna/_rtl_ab_comparison_stages_pd.mro
class BetaDistribution(TypedDict):
    """Beta distribution parameters."""

    alpha: float
    beta: float


class StudentTParams(TypedDict):
    """Student's t-distribution parameters."""

    df: float
    loc: float
    scale: float


class GenNormalParams(TypedDict):
    """Generalized normal distribution parameters."""

    loc: float
    scale: float
    beta: float


class LogNormalParams(TypedDict):
    """Log-normal distribution parameters."""

    shape: float
    scale: float
    loc: float


class GammaParams(TypedDict):
    """Gamma distribution parameters."""

    alpha: float
    beta: float
    loc: float


class GeneticDemuxParams(TypedDict):
    """Genetic demultiplexing parameters."""

    num_clusters: int
    candidate_snps: vcf
    known_genotypes: KnownGenotypes
    min_ref_umi: int
    min_alt_umi: int
    min_ref_cells: int
    min_alt_cells: int
    restarts: int


class KnownGenotypes(TypedDict):
    """Known genotypes specification."""

    vcf: vcf
    sample_names: list[str]


# mro/rna/stages/multi/
class FlexVdjParams(TypedDict):
    """Flex VDJ scheme parameters."""

    flex_vdj_scheme: str | None


class BarcodeAssignments(TypedDict):
    """Forced sample barcode assignments."""

    sample_barcodes: json | None
    cells_per_tag: json | None
    non_singlet_barcodes: json | None


class SampleBamFile(TypedDict):
    """Per-sample BAM file and its index."""

    sample: str
    bam_file: bam | None
    bai_index_file: bam_bai | None
    csi_index_file: bam_csi | None


class SampleMetrics(TypedDict):
    """Per-sample metrics files."""

    sample: str
    summary: json
    per_barcode_metrics: csv | None


class SampleMoleculeInfo(TypedDict):
    """Per-sample molecule info."""

    sample: str
    h5_file: h5
    summary: json


class SampleMatrices(TypedDict):
    """Per-sample count matrices."""

    sample: str
    filtered_matrix_h5: h5
    filtered_matrix_mex: path
    raw_matrix_h5: h5
    raw_probe_bc_matrix: h5 | None
    raw_matrix_mex: path
    filtered_barcodes: csv
    aggregate_barcodes: csv | None
    per_probe_metrics: csv | None


class SampleSlfeOuts(TypedDict):
    """Per-sample SLFE outputs."""

    sample: str
    bam_file: bam | None
    bai_index_file: bam_bai | None
    csi_index_file: bam_csi | None
    metrics_summary: json
    per_barcode_metrics: csv | None
    molecule_info: h5
    filtered_matrix_h5: h5
    filtered_matrix_mex: path
    raw_matrix_h5: h5
    raw_probe_bc_matrix: h5 | None
    raw_matrix_mex: path
    filtered_barcodes: csv
    aggregate_barcodes: csv | None
    per_probe_metrics: csv | None


# mro/rna/_sc_rna_aggregator_stages_pd.mro
class MetaAnalysisDef(TypedDict):
    """Meta-analysis definition for PD stages."""

    sample_id: str
    sample_desc: str
    pipeline_name: str
    pipeline_path: path
    pipeline_version: str
    genome: dict
    target_set: str
    sample_defs: list[dict]  # TODO: be more specific with the type
    sequencing_libraries: dict  # TODO: be more specific with the type


# mro/rna/_cr_aggr_stages.mro
class VdjAggrCsvLibrary(TypedDict):
    """VDJ aggregation CSV library definition."""

    library_id: str
    vdj_contig_info: path
    donor: str
    origin: str
    meta: dict[str, str]


class VdjAggrInput(TypedDict):
    """VDJ aggregation input definition."""

    libraries: list[VdjAggrCsvLibrary]


class FilterSwitch(TypedDict):
    """Enclone filter switches."""

    asm_shared_contig: bool
    enclone_shared_contig: bool
    enclone_multiplet: bool
    enclone_umi: bool
    enclone_cross: bool


# mro/rna/stages/common
class SampleDef(TypedDict):
    """Sample definition."""

    fastq_mode: str
    gem_group: int
    lanes: list[int]
    library_type: str
    r1_length: int
    r2_length: int
    read_path: path
    sample_indices: list[str]
    sample_names: list[str]
    subsample_rate: float
    fastq_id: str


class VdjInputsCS(TypedDict):
    """VDJ inputs for customer-facing stages."""

    sample_def: list[SampleDef]
    denovo: bool
    inner_enrichment_primers: path
    chain_type: str
    physical_library_id: str
    r1_length: int
    r2_length: int
    chemistry_spec: str
    max_reads_per_barcode: int


class V1PatternFixArgs(TypedDict):
    """V1 pattern fix arguments."""

    v1_filtered_fbm: h5
    v1_pattern_type: int


class SegmentationInputs(TypedDict):
    """Segmentation input parameters."""

    user_provided_segmentations: geojson
    square_barcode_to_cell_map: csv
    instance_mask_tiff: tiff
    instance_mask_npy: npy
    max_nucleus_diameter_px: int
    barcode_assignment_distance_micron: int


# mro/rna/_sc_antibody_analyses_pd_stages.mro
class VdjCellsPD(TypedDict):
    """VDJ cells PD struct."""

    vdj_t: None | json
    vdj_t_gd: None | json
    vdj_b: None | json


class WhiteList(TypedDict):
    """Whitelist specification."""

    kind: str
    name: str


class CaptureSeq(TypedDict):
    """Capture sequence definition."""

    length: int
    offset: int
    read_type: str
    sequence: str


class ProbeBarcode(TypedDict):
    """Probe barcode definition."""

    length: int
    offset: int
    read_type: str
    whitelist: WhiteList


class Target(TypedDict):
    """Target definition."""

    length: int
    offset: int
    read_type: str
    type: str
    whitelist: WhiteList


class Umi(TypedDict):
    """UMI definition."""

    length: int
    offset: int
    read_type: str


class ChemDef(TypedDict):
    """Chemistry definition."""

    name: str
    umi: list[Umi]
    target: list[Target] | None
    probe_barcode: list[ProbeBarcode]
    capture_seq: list[CaptureSeq]


class ReducedSampleDef(TypedDict):
    """Reduced sample definition."""

    fastq_mode: str
    lanes: list[int]
    read_path: path
    sample_indices: list[str]


# mro/rna/_cr_lib_stages.mro
class VdjInputs(TypedDict):
    """VDJ inputs struct (internal, used by SPLIT_VDJ_INPUTS and related stages)."""

    sample_def: list[SampleDef]
    chemistry_spec: str
    custom_chemistry_def: ChemistryDef
    primers: list[dict[str, Any]]
    subsample_rate: float
    initial_reads: int
    primer_initial_reads: int
    special_genomic_regions: list[str]
    denovo: bool
    r1_length: int
    r2_length: int
    inner_enrichment_primers: path
    chain_type: str
    physical_library_id: str
    max_reads_per_barcode: int


# mro/rna/_sc_multi_defs.mro
class VdjOutputsCS(TypedDict):
    """VDJ customer-facing outputs struct."""

    metrics_summary_csv: None | csv
    clonotypes: None | csv
    filtered_contig_fasta: None | fasta
    filtered_contig_fastq: None | fastq
    filtered_contig_annotations_csv: None | csv
    all_contig_fasta: None | fasta
    all_contig_fasta_fai: None | fasta_fai
    all_contig_fastq: None | fastq
    all_contig_bam: None | bam
    all_contig_bam_bai: None | bam_bai
    all_contig_annotations_json: None | json
    all_contig_annotations_bed: None | bed
    all_contig_annotations_csv: None | csv
    cell_barcodes: None | json
    consensus_fasta: None | fasta
    consensus_fasta_fai: None | fasta_fai
    consensus_bam: None | bam
    consensus_bam_bai: None | bam_bai
    consensus_annotations_csv: None | csv
    concat_ref_fasta: None | fasta
    concat_ref_fasta_fai: None | fasta_fai
    concat_ref_bam: None | bam
    concat_ref_bam_bai: None | bam_bai
    vloupe: None | vloupe_
    airr_rearrangement: None | tsv
    donor_regions: None | fa
    vdj_contig_info: None | pb


# mro/rna/_sc_rna_meta_analysis_stages_pd.mro
class TargetedMetadata(TypedDict):
    """Targeted metadata for targeted meta-analysis."""

    sample_id: str
    sample_idx: int
    matched_control_id: str
    matched_control_idx: int
    is_matched_control: bool
    matched_targeted_ids: list[str]
    target_sets: dict[str, csv]
    target_set_complements: dict[str, csv]
    target_set_csvs: dict[str, csv]
    target_set_metadata_csvs: dict[str, csv]
    orig_filtered_barcodes_csv: csv
    matched_filtered_barcodes_csv: csv
    num_cells_to_subset: int
