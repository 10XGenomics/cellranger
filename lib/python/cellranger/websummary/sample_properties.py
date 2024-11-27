# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
"""These define information about the sample required to generate a web summary."""
# pylint: disable=too-few-public-methods,missing-docstring,too-many-arguments


from __future__ import annotations

import json
from dataclasses import dataclass

import martian
from typing_extensions import Self

from cellranger.chemistry import (
    CHEMISTRY_SC3P_LT,
    ChemistryDef,
    ChemistryDefs,
    get_primary_chemistry_def,
)
from cellranger.version import get_version

CYTASSIST_RUN_NAME_KEY = "cytassistRunName"
CYTASSIST_SERIAL_KEY = "cytassistInstrumentSerial"
CYTASSIST_SOFTWARE_VERSION = "cytassistInstrumentSoftwareVersion"
CYTASSIST_NOT_FOUND_SENTINEL = "Not Found"


class SampleProperties:
    """Sample properties."""

    sample_id: str
    sample_desc: str
    version_from_git: bool

    def __init__(self, sample_id, sample_desc, version_from_git=False):
        self.sample_id = sample_id
        self.sample_desc = sample_desc
        if not version_from_git:
            self.version = martian.get_pipelines_version()
        else:
            self.version = get_version()


@dataclass
class CytassistRunProperties:
    cytassist_run_name: str
    cytassist_instrument_serial: str
    cytassist_instrument_software_version: str

    @classmethod
    def from_json(cls, json_path: str | bytes) -> Self:
        with open(json_path) as f:
            metadata = json.load(f)

        return cls(
            cytassist_run_name=(
                metadata.get(CYTASSIST_RUN_NAME_KEY)
                if metadata.get(CYTASSIST_RUN_NAME_KEY)
                else CYTASSIST_NOT_FOUND_SENTINEL
            ),
            cytassist_instrument_serial=(
                metadata.get(CYTASSIST_SERIAL_KEY)
                if metadata.get(CYTASSIST_SERIAL_KEY)
                else CYTASSIST_NOT_FOUND_SENTINEL
            ),
            cytassist_instrument_software_version=(
                metadata.get(CYTASSIST_SOFTWARE_VERSION)
                if metadata.get(CYTASSIST_SOFTWARE_VERSION)
                else CYTASSIST_NOT_FOUND_SENTINEL
            ),
        )


# pylint: disable=too-many-instance-attributes
class CountSampleProperties(SampleProperties):
    """Sample properties for Count, Aggr, Reanalyze, Spatial web summaries."""

    genomes: list[str]
    is_spatial: bool
    target_set: str | None
    target_panel_summary: str | None
    feature_ref_path: str | None
    include_introns: bool
    reorientation_mode: str | None
    filter_probes: bool | None
    aligner: str | None
    redundant_loupe_alignment: bool
    loupe_alignment_file: str | None
    v1_pattern_fix: dict | None
    default_layout: bool | None
    override_id: bool | None
    slide_id_mismatch: bool | None
    is_visium_hd: bool | None
    cmdline: str | None
    itk_error_string: str | None

    # pylint: disable=too-many-locals
    def __init__(
        self,
        sample_id,
        sample_desc,
        genomes,
        version_from_git=False,
        is_spatial=False,
        target_set=None,
        target_panel_summary=None,
        feature_ref_path=None,
        include_introns=False,
        reorientation_mode=None,
        filter_probes=None,
        aligner=None,
        redundant_loupe_alignment=False,
        loupe_alignment_file=None,
        v1_pattern_fix=None,
        default_layout=False,
        override_id=False,
        slide_id_mismatch=False,
        is_visium_hd=False,
        cmdline=None,
        itk_error_string=None,
    ):
        super().__init__(sample_id, sample_desc, version_from_git=version_from_git)
        self.genomes = genomes
        self.is_spatial = is_spatial
        self.target_set = target_set
        self.target_panel_summary = target_panel_summary
        self.feature_ref_path = feature_ref_path
        self.include_introns = include_introns
        self.reorientation_mode = reorientation_mode
        self.filter_probes = filter_probes
        self.aligner = aligner
        self.redundant_loupe_alignment = redundant_loupe_alignment
        self.loupe_alignment_file = loupe_alignment_file
        self.v1_pattern_fix = v1_pattern_fix
        self.default_layout = default_layout
        self.override_id = override_id
        self.slide_id_mismatch = slide_id_mismatch
        self.is_visium_hd = is_visium_hd
        self.cmdline = cmdline
        self.itk_error_string = itk_error_string

    @property
    def is_targeted(self):
        return self.target_set is not None

    @property
    def is_lt(self):
        return False


class ExtendedCountSampleProperties(CountSampleProperties):
    """Properties for a count run."""

    reference_path: str
    chemistry_defs: ChemistryDefs
    cytassist_run_properties: CytassistRunProperties | None

    # pylint: disable=too-many-locals
    def __init__(
        self,
        sample_id,
        sample_desc,
        genomes,
        reference_path,
        chemistry_defs,
        target_set=None,
        target_panel_summary=None,
        feature_ref_path=None,
        version_from_git=False,
        is_spatial=False,
        include_introns=False,
        reorientation_mode=None,
        filter_probes=None,
        disable_ab_aggregate_detection=False,
        aligner=None,
        redundant_loupe_alignment=False,
        loupe_alignment_file=None,
        v1_pattern_fix=None,
        default_layout=False,
        override_id=False,
        slide_id_mismatch=False,
        is_visium_hd=False,
        cmdline=None,
        cytassist_run_metrics=None,
        itk_error_string=None,
    ):
        super().__init__(
            sample_id,
            sample_desc,
            genomes,
            version_from_git=version_from_git,
            is_spatial=is_spatial,
            target_set=target_set,
            target_panel_summary=target_panel_summary,
            feature_ref_path=feature_ref_path,
            include_introns=include_introns,
            reorientation_mode=reorientation_mode,
            filter_probes=filter_probes,
            aligner=aligner,
            redundant_loupe_alignment=redundant_loupe_alignment,
            loupe_alignment_file=loupe_alignment_file,
            v1_pattern_fix=v1_pattern_fix,
            default_layout=default_layout,
            override_id=override_id,
            slide_id_mismatch=slide_id_mismatch,
            is_visium_hd=is_visium_hd,
            cmdline=cmdline,
            itk_error_string=itk_error_string,
        )
        self.reference_path = reference_path
        self.chemistry_defs = chemistry_defs
        self.disable_ab_aggregate_detection = disable_ab_aggregate_detection
        if cytassist_run_metrics:
            self.cytassist_run_properties = CytassistRunProperties.from_json(cytassist_run_metrics)
        else:
            self.cytassist_run_properties = None

    def chemistry_description(self) -> str:
        """Return the chemistry description of the primary chemistry."""
        return get_primary_chemistry_def(self.chemistry_defs)["description"]

    @property
    def is_lt(self) -> bool:
        """Return whether any chemistry is low throughput."""
        return CHEMISTRY_SC3P_LT in self.chemistry_defs.values()


class AggrCountSampleProperties(CountSampleProperties):
    """Properties from an Aggr Run."""

    def __init__(
        self,
        sample_id,
        sample_desc,
        genomes,
        agg_batches,
        is_spatial,
        target_set=None,
        target_panel_summary=None,
        feature_ref_path=None,
        version_from_git=False,
    ):
        super().__init__(
            sample_id,
            sample_desc,
            genomes,
            version_from_git=version_from_git,
            is_spatial=is_spatial,
            target_set=target_set,
            target_panel_summary=target_panel_summary,
            feature_ref_path=feature_ref_path,
        )
        self.agg_batches = agg_batches


class VdjSampleProperties(SampleProperties):
    chemistry_def: ChemistryDef
    chain_type: str

    def __init__(self, sample_id, sample_desc, chemistry_def, chain_type, version_from_git=False):
        super().__init__(sample_id, sample_desc, version_from_git=version_from_git)
        self.chemistry_def = chemistry_def
        self.chain_type = chain_type

    def chemistry_description(self) -> str:
        """Return the chemistry description."""
        return self.chemistry_def["description"]

    @property
    def is_spatial(self):
        return False


class SampleDataPaths:  # pylint: disable=too-many-instance-attributes
    def __init__(
        self,
        summary_path=None,
        barcode_summary_path=None,
        analysis_path=None,
        filtered_barcodes_path=None,
        feature_metrics_path=None,
        antibody_histograms_path=None,
        antibody_treemap_path=None,
        raw_normalized_heatmap_path=None,
        isotype_scatter_path=None,
        gex_fbc_correlation_heatmap_path=None,
        antigen_histograms_path=None,
        antigen_treemap_path=None,
        vdj_clonotype_summary_path=None,
        vdj_barcode_support_path=None,
        vdj_cell_barcodes_path=None,
    ):
        assert filtered_barcodes_path is None or vdj_cell_barcodes_path is None
        self.summary_path = summary_path
        self.barcode_summary_path = barcode_summary_path
        self.analysis_path = analysis_path
        self.filtered_barcodes_path = filtered_barcodes_path
        self.feature_metrics_path = feature_metrics_path
        self.antibody_histograms_path = antibody_histograms_path
        self.antibody_treemap_path = antibody_treemap_path
        self.raw_normalized_heatmap_path = raw_normalized_heatmap_path
        self.isotype_scatter_path = isotype_scatter_path
        self.gex_fbc_correlation_heatmap_path = gex_fbc_correlation_heatmap_path
        self.antigen_histograms_path = antigen_histograms_path
        self.antigen_treemap_path = antigen_treemap_path
        self.vdj_clonotype_summary_path = vdj_clonotype_summary_path
        self.vdj_barcode_support_path = vdj_barcode_support_path
        self.vdj_cell_barcodes_path = vdj_cell_barcodes_path


class CellTypeDataPaths:  # pylint: disable=too-many-instance-attributes
    def __init__(self, summary_path=None, cell_type_bar_chart=None, cell_type_umap_plot=None):
        self.summary_path = summary_path
        self.cell_type_bar_chart = cell_type_bar_chart
        self.cell_type_umap_plot = cell_type_umap_plot
