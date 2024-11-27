#!/usr/bin/env python
#
# Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
#

"""Here, you could find everything shown in the summary tab."""


from __future__ import annotations

import os.path
from typing import TYPE_CHECKING

import cellranger.constants as cr_constants
import cellranger.molecule_counter as cr_mc
import cellranger.rna.library as rna_library
import cellranger.vdj.constants as vdj_constants
import cellranger.webshim.common as cr_webshim
import cellranger.webshim.constants.shared as shared_constants
import cellranger.websummary.plotly_tools as pltly
from cellranger.feature.throughputs import (
    HT_THROUGHPUT,
    INCONSISTENT_THROUGHPUT_METRIC,
    LT_THROUGHPUT,
    MT_THROUGHPUT,
    THROUGHPUT_INFERRED_METRIC,
)
from cellranger.h5_constants import H5_CHEMISTRY_DESC_KEY
from cellranger.reference_paths import get_ref_name_from_genomes
from cellranger.targeted.targeted_constants import TARGETING_METHOD_HC, TARGETING_METHOD_TL
from cellranger.webshim.data import FILTERED_BCS_TRANSCRIPTOME_UNION
from cellranger.websummary.isotypes import ALT_LINK_HELP_TEXT
from cellranger.websummary.metrics import INFO_THRESHOLD, WARNING_THRESHOLD
from cellranger.websummary.sample_properties import (
    AggrCountSampleProperties,
    CountSampleProperties,
    ExtendedCountSampleProperties,
    SampleProperties,
    VdjSampleProperties,
)

if TYPE_CHECKING:
    from cellranger.webshim.data import SampleData

ALARMS = shared_constants.ALARMS

_ANTIBODY_reads_lost_to_aggregate_GEMs = "ANTIBODY_reads_lost_to_aggregate_GEMs"
_ANTIGEN_reads_lost_to_aggregate_GEMs = "ANTIGEN_reads_lost_to_aggregate_GEMs"
ARC_V1_MULTIOME_DESCRIPTION = "Single Cell Multiome ATAC + Gene Expression v1"

# Feature barcoding internel name <-> display name
FEATURE_BARCODINGS = [
    rna_library.CRISPR_METRIC_PREFIX,
    rna_library.ANTIBODY_METRIC_PREFIX,
    rna_library.ANTIGEN_METRIC_PREFIX,
    rna_library.CUSTOM_METRIC_PREFIX,
]

FB_DISPLAY_NAME = {
    rna_library.CRISPR_METRIC_PREFIX: "CRISPR",
    rna_library.ANTIBODY_METRIC_PREFIX: "Antibody",
    rna_library.ANTIGEN_METRIC_PREFIX: "Antigen",
    rna_library.CUSTOM_METRIC_PREFIX: "Custom Feature",
}

RANK_PLOT_HELP = [
    [
        "Barcode Rank Plot",
        [
            "The plot shows the count of filtered UMIs mapped to each barcode. Barcodes are not determined to be cell-associated strictly based on their UMI count. Instead, they could be determined based on their expression profile, or removed via Protein Aggregate Detection and Filtering and/or High Occupancy GEM Filtering. Therefore, some regions of the graph contain both cell-associated and background-associated barcodes. The color of the graph in these regions is based on the local density of barcodes that are cell-associated. Hovering over the plot displays the total number and percentage of barcodes in that region called as cells along with the number of UMI counts for those barcodes and barcode rank, ordered in descending order of UMI counts."
        ],
    ]
]

ZOOM_IMAGE_HELP = [
    [
        "Tissue Detection and Fiducial Alignment",
        [
            (
                "Shows the tissue image in gray tones with an overlay of the aligned fiducial frame (open blue circles) and the capture area spots (gray circles). "
                "For the latter, the circles filled in red denote selected tissue-associated spots and the remaining open gray circles denote unselected spots. "
                "Hover mouse cursor over the image to magnify the view. "
                "Confirm fiducial frame aligns well with fiducial spots, e.g. the corner shapes match, and confirm selection of tissue-covered spots. "
                "If the result shows poor fiducial alignment or tissue detection, consider sharing the image with support@10xgenomics.com so we can improve the algorithm. "
                "Otherwise, perform manual alignment and spot selection with Loupe Browser."
            )
        ],
    ],
]

REGISTRATION_IMAGE_HELP = [
    "CytAssist Image Alignment",
    [
        (
            "Shows the CytAssist image aligned onto the microscope image. "
            "Click-drag the opacity slider to blend the two images to confirm good alignment, i.e. tissue boundaries and features should remain fixed in place. "
            "For QC purposes, fluorescence microscopy images are inverted to have a light background. "
            "If alignment is poor, rerun with Loupe Browser manual alignment.",
        )
    ],
]


COMMON_CELL_CALLING_METRIC_KEYS = [
    "filtered_bcs",
    "filtered_bcs_conf_mapped_barcoded_reads_cum_frac",
    "multi_transcriptome_total_raw_reads_per_filtered_bc",
    "filtered_reads_per_filtered_bc",
]

CELL_CALLING_METRIC_KEYS = COMMON_CELL_CALLING_METRIC_KEYS + [
    "filtered_bcs_median_counts",
    "filtered_bcs_median_unique_genes_detected",
    "filtered_bcs_total_unique_genes_detected",
]

TARGETED_CELL_CALLING_METRIC_KEYS_HC = COMMON_CELL_CALLING_METRIC_KEYS + [
    "total_targeted_reads_per_filtered_bc",
    "median_umis_per_cell_on_target",
    "median_genes_per_cell_on_target",
    "num_genes_detected_on_target",
]

TARGETED_CELL_CALLING_METRIC_KEYS_TL = COMMON_CELL_CALLING_METRIC_KEYS + [
    "median_umis_per_cell_on_target",
    "median_genes_per_cell_on_target",
    "num_genes_detected_on_target",
]

ANTIBODY_filtered_bcs_transcriptome_union = "ANTIBODY_filtered_bcs_transcriptome_union"
ANTIBODY_CELL_CALLING_METRIC_KEYS = [
    ANTIBODY_filtered_bcs_transcriptome_union,
    "ANTIBODY_multi_transcriptome_total_raw_reads_per_filtered_bc",
]

CELL_CALLING_ALARM_KEYS = ["filtered_bcs_conf_mapped_barcoded_reads_cum_frac"]
TARGETED_CELL_CALLING_ALARM_KEYS_HC = [
    "filtered_bcs_conf_mapped_barcoded_reads_cum_frac",
    FILTERED_BCS_TRANSCRIPTOME_UNION,
    "filtered_reads_per_filtered_bc",
]

# metric keys for sequencing (GEX and feature barcoding)
SEQUENCING_METRIC_KEYS = [
    "total_read_pairs",  # CR Number of reads
    "sequenced_reads_count",  # SR Number of reads.
    "unprocessed_read_pairs",
    "good_bc_frac",
    "good_umi_frac",
    "multi_cdna_pcr_dupe_reads_frac",
    "bc_bases_with_q30_frac",
    "read_bases_with_q30_frac",
    "read2_bases_with_q30_frac",
    "umi_bases_with_q30_frac",
    "square_008um.bins_under_tissue_frac",  # HD specific
]

TARGETED_SEQUENCING_METRIC_KEYS = [
    "total_read_pairs",  # CR Number of reads
    "sequenced_reads_count",
    "unprocessed_read_pairs",
    "subsampled_frac",
    "good_bc_frac",
    "good_umi_frac",
    "multi_cdna_pcr_dupe_reads_frac_on_target",
    "bc_bases_with_q30_frac",
    "read_bases_with_q30_frac",
    "read2_bases_with_q30_frac",
    "umi_bases_with_q30_frac",
    "square_008um.bins_under_tissue_frac",  # HD specific
]

SEQUENCING_ALARM_KEYS = [
    "good_bc_frac",
    "good_umi_frac",
    "bc_bases_with_q30_frac",
    "read_bases_with_q30_frac",
    "umi_bases_with_q30_frac",
    "fraction_bc_outside_image",
]

AGGREGATION_METRIC_KEYS = [
    "frac_reads_kept",
    "pre_normalization_raw_reads_per_filtered_bc",
    "pre_normalization_cmb_reads_per_filtered_bc",
    "pre_normalization_targeted_cmb_reads_per_filtered_bc",
]

# metric keys for feature barcoding application
FB_APP_METRIC_KEYS = {
    "CRISPR": [
        "CRISPR_feature_bc_extracted_frac",
        "CRISPR_recognized_feature_bc_frac",
        "CRISPR_frac_feature_reads_usable",
        "CRISPR_feature_reads_usable_per_cell",
        "CRISPR_unrecognized_feature_bc_frac",
        "CRISPR_feature_reads_in_cells",
        "CRISPR_frac_cells_with_protospacer",
        "CRISPR_frac_cells_with_multiple_protospacer",
        "CRISPR_multi_filtered_bcs_median_counts",
    ],
    "ANTIBODY": [
        "ANTIBODY_recognized_feature_bc_frac",
        "ANTIBODY_frac_feature_reads_usable",
        "ANTIBODY_feature_reads_usable_per_cell",
        _ANTIBODY_reads_lost_to_aggregate_GEMs,
        "ANTIBODY_unrecognized_feature_bc_frac",
        "ANTIBODY_feature_reads_in_cells",
        "ANTIBODY_multi_filtered_bcs_median_counts",
    ],
    "ANTIGEN": [
        "ANTIGEN_recognized_feature_bc_frac",
        "ANTIGEN_frac_feature_reads_usable",
        "ANTIGEN_feature_reads_usable_per_cell",
        _ANTIGEN_reads_lost_to_aggregate_GEMs,
        "ANTIGEN_unrecognized_feature_bc_frac",
        "ANTIGEN_feature_reads_in_cells",
        "ANTIGEN_multi_filtered_bcs_median_counts",
    ],
    "Custom": [
        "Custom_recognized_feature_bc_frac",
        "Custom_frac_feature_reads_usable",
        "Custom_feature_reads_usable_per_cell",
        "Custom_unrecognized_feature_bc_frac",
        "Custom_feature_reads_in_cells",
        "Custom_multi_filtered_bcs_median_counts",
    ],
}

# targeted hero metric keys for targeted GEX -- hacky way to override hero metrics
# dict keys are the targeted hero metrics to display, values are the original GEX
# hero metrics, which are also the keys used to display the values in the html template

HERO_METRIC_MAPPING = {
    TARGETING_METHOD_TL: {
        FILTERED_BCS_TRANSCRIPTOME_UNION: FILTERED_BCS_TRANSCRIPTOME_UNION,
        "multi_transcriptome_total_raw_reads_per_filtered_bc": "multi_transcriptome_total_raw_reads_per_filtered_bc",
        "median_genes_per_cell_on_target": "filtered_bcs_median_unique_genes_detected",
    },
    TARGETING_METHOD_HC: {
        FILTERED_BCS_TRANSCRIPTOME_UNION: FILTERED_BCS_TRANSCRIPTOME_UNION,
        "multi_transcriptome_total_raw_reads_per_filtered_bc": "multi_transcriptome_total_raw_reads_per_filtered_bc",
        "multi_transcriptome_targeted_conf_mapped_reads_frac": "filtered_bcs_median_unique_genes_detected",
    },
}


MAPPING_KEYS = [
    "genome_mapped_reads_frac",
    "genome_conf_mapped_reads_frac",
    "intergenic_conf_mapped_reads_frac",
    "intronic_conf_mapped_reads_frac",
    "exonic_conf_mapped_reads_frac",
    "transcriptome_conf_mapped_reads_frac",
    "multi_transcriptome_targeted_conf_mapped_reads_frac",
    "antisense_reads_frac",
]

TEMP_LIG_MAPPING_KEYS = [
    "genome_mapped_reads_frac",
    "genome_conf_mapped_reads_frac",
    "multi_transcriptome_targeted_conf_mapped_reads_frac",
    "multi_transcriptome_half_mapped_reads_frac",
    "multi_transcriptome_split_mapped_reads_frac",
]

MAPPING_ALARMS = [
    "transcriptome_conf_mapped_reads_frac",
    "antisense_reads_frac",
]

TARGETED_MAPPING_ALARMS = [
    "multi_transcriptome_targeted_conf_mapped_reads_frac",
    "targeted_unsupported_panel",
] + MAPPING_ALARMS

TEMP_LIG_MAPPING_ALARMS = [
    "genome_conf_mapped_reads_frac",
    "multi_transcriptome_targeted_conf_mapped_reads_frac",
    "multi_transcriptome_half_mapped_reads_frac",
    "multi_transcriptome_split_mapped_reads_frac",
]
HD_TEMP_LIG_MAPPING_ALARMS = TEMP_LIG_MAPPING_ALARMS + ["estimated_gdna_content"]


def get_empty_rank_plot():
    """Generates a template for the barcode rank plot.

    The template can be used by components making different types of rank plots.
    config/layout are pre-filled, the consumer is responsible for adding the data.

    Returns:
        dict: data for a plotly plot
    """
    return {
        "config": pltly.PLOT_CONFIG,
        "layout": {
            "title": "Barcode Rank Plot",
            "xaxis": {
                "title": "Barcodes",
                "type": "log",
                "showline": True,
                "zeroline": False,
                "fixedrange": False,
            },
            "yaxis": {
                "title": "UMI counts",
                "type": "log",
                "showline": True,
                "zeroline": False,
                "fixedrange": False,
            },
            "font": pltly.DEFAULT_WEB_FONT,
            "hovermode": "closest",
        },
        "data": [],
    }


def add_data(websummary_json, alarm_list, input_data):
    """Adds data to global dictionary."""
    if input_data is None:
        return

    if ALARMS in input_data:
        alarm_list.extend(input_data[ALARMS])
        del input_data[ALARMS]
    websummary_json.update(input_data)
    return


def hero_metrics(metadata, sample_data: SampleData, species_list):
    if sample_data is None or sample_data.summary is None:
        return None

    def _generate_data_with_alarms(data, alarm_keys):
        alarms = metadata.gen_metric_list(sample_data.summary, alarm_keys, species_list)
        new_alarms = [metric.alarm_dict for metric in alarms if metric.alarm_dict]
        if new_alarms:
            data[ALARMS] = new_alarms
        return data

    # For FB-only web summaries, only use cell counts and total antibody reads
    if sample_data.is_antibody_only():
        data = {}
        for key in ANTIBODY_CELL_CALLING_METRIC_KEYS:
            metrics = metadata.gen_metric_list(sample_data.summary, [key], species_list)
            for metric in metrics:
                # remove the prefix, so the key name matches with the input json expects
                new_key = metric.key.replace("ANTIBODY_", "")
                data[new_key] = metric.gen_metric_dict()

        alarm_keys = [ANTIBODY_filtered_bcs_transcriptome_union]
        data_with_alarms = _generate_data_with_alarms(data, alarm_keys)

        return data_with_alarms

    if sample_data.targeting_method is not None:
        data = {}
        for key, new_key in HERO_METRIC_MAPPING[sample_data.targeting_method].items():
            metrics = metadata.gen_metric_list(sample_data.summary, [key], species_list)
            for metric in metrics:
                # translate the key, so the key name matches with the input json expects
                data[new_key] = metric.gen_metric_dict()

        data_with_alarms = _generate_data_with_alarms(data, [])
        return data_with_alarms

    data = {}
    for key in [
        FILTERED_BCS_TRANSCRIPTOME_UNION,
        "multi_transcriptome_total_raw_reads_per_filtered_bc",
    ]:
        metrics = metadata.gen_metric_list(sample_data.summary, [key], species_list)
        for metric in metrics:
            data[metric.key] = metric.gen_metric_dict()

    is_barnyard = len(species_list) > 1
    if not is_barnyard:
        median_unique_genes = "filtered_bcs_median_unique_genes_detected"
        metrics = metadata.gen_metric_list(sample_data.summary, [median_unique_genes], species_list)
        for metric in metrics:
            data[median_unique_genes] = metric.gen_metric_dict()

    alarm_keys = [FILTERED_BCS_TRANSCRIPTOME_UNION]
    data_with_alarms = _generate_data_with_alarms(data, alarm_keys)

    return data_with_alarms


def _get_chemistry_description(
    sample_data: SampleData,
    sample_properties: SampleProperties | None = None,
    sample_def: dict | None = None,
) -> str:
    """This function returns the chemistry description of the given sample data.

    If sample_def is provided (Aggr), the function uses the molecule_h5 file to retrieve the chemistry description.
    If sample_def is not provided, the function uses the sample_data and sample_properties to determine the chemistry description.
    If the sample is a spatial sample, the function also appends the assay type to the chemistry description.

    Args:
        sample_data (required): The sample data for which the chemistry description is to be returned
        sample_properties (optional): Additional properties of the sample that may be needed to determine the chemistry description.
        sample_def (optional): A dictionary containing the molecule_h5 file and other information needed to retrieve the chemistry description from the file.

    Returns:
        chemistry (string): Chemistry description
    """
    if sample_def:
        # Used by aggr
        mol_info = cr_mc.MoleculeCounter.open(sample_def["molecule_h5"], "r")
        chemistry_description = mol_info.get_metric("chemistry_description")
        is_spatial = mol_info.is_spatial_data()
    else:
        assert sample_data.summary is not None
        assert isinstance(sample_properties, CountSampleProperties | VdjSampleProperties)
        chemistry_description = (
            sample_data.summary.get("chemistry_description")
            or sample_properties.chemistry_description()
        )
        is_spatial = sample_properties.is_spatial

    if not is_spatial:
        return chemistry_description

    if sample_data.targeting_method == TARGETING_METHOD_TL:
        # covers both fresh/fixed frozen and FFPE products
        # Determine the assay from the chemistry.
        assay_version = {
            "Visium V1 Slide": "v1",
            "Visium V2 Slide": "v1",
            "Visium V3 Slide": "v2",
            "Visium V4 Slide": "v2",
            "Visium V5 Slide": "v2",
            "Visium HD v1": "v1",
        }.get(chemistry_description)
        assay_base = {"Visium HD v1": "HD probe-based"}.get(chemistry_description, "Probe-based")
        assay = f"{assay_base} {assay_version}" if assay_version else assay_base
    elif sample_data.targeting_method == TARGETING_METHOD_HC:
        assay = "Targeted"
    else:
        assay = "3'"
    sanitised_chemistry_description = {"Visium HD v1": "Visium HD H1 Slide"}.get(
        chemistry_description, chemistry_description
    )
    return f"{sanitised_chemistry_description} - {assay}"


def _add_throughput_to_chemistry(chemistry_desc: str, sample_data: SampleData) -> str:
    throughput_inferred = sample_data.summary.get(THROUGHPUT_INFERRED_METRIC)

    chemistry_with_throughput = chemistry_desc
    if throughput_inferred == HT_THROUGHPUT and chemistry_desc[-2:] not in [
        LT_THROUGHPUT,
        MT_THROUGHPUT,
        HT_THROUGHPUT,
    ]:
        chemistry_with_throughput = f"{chemistry_desc} {HT_THROUGHPUT}"

    if (
        chemistry_desc.endswith(HT_THROUGHPUT)
        and throughput_inferred
        and throughput_inferred != HT_THROUGHPUT
    ):
        sample_data.summary[INCONSISTENT_THROUGHPUT_METRIC] = HT_THROUGHPUT
    return chemistry_with_throughput


# pylint: disable=too-many-branches
def pipeline_info_table(
    sample_data: SampleData,
    sample_properties: SampleProperties,
    pipeline: str,
    metadata=None,
    species_list: list[str] | None = None,
    sample_defs: list[dict] | None = None,
):
    """Generates a table of general pipeline information.

    Args:
        sample_data:
        sample_properties:
        pipeline:
        metadata:
        species_list:
        sample_defs: Map of sample_defs passed in from PARSE_CSV
    """
    assert isinstance(sample_properties, SampleProperties)

    if sample_data is None or sample_data.summary is None:
        return None

    alarms = []

    # Get the chemistry description for the sample(s)
    if isinstance(sample_properties, AggrCountSampleProperties) and sample_defs:
        chemistry_row = []
        for sample_def in sample_defs:
            chemistry = _get_chemistry_description(sample_data=sample_data, sample_def=sample_def)
            chemistry_row.append([f"Chemistry ({sample_def['library_id']})", chemistry])
    elif not sample_properties.is_spatial and isinstance(
        sample_properties, ExtendedCountSampleProperties
    ):
        # We may have different chemistries per library; need to collect them.
        # If chemistry for every library is the same, only produce one line.
        # Do not do this for spatial.
        chem_desc_unique = set(
            chem_def["description"] for chem_def in sample_properties.chemistry_defs.values()
        )
        if len(chem_desc_unique) == 1:
            chemistry_row = [
                ["Chemistry", _add_throughput_to_chemistry(chem_desc_unique.pop(), sample_data)]
            ]
        else:
            chemistry_row = [
                [
                    f"Chemistry ({lib_type})",
                    _add_throughput_to_chemistry(chem_def["description"], sample_data),
                ]
                for (lib_type, chem_def) in sample_properties.chemistry_defs.items()
            ]
            # Ensure stable ordering by alphabetized library type.
            chemistry_row.sort()
    else:
        chemistry = _get_chemistry_description(
            sample_data=sample_data, sample_properties=sample_properties
        )
        chemistry_row = [["Chemistry", _add_throughput_to_chemistry(chemistry, sample_data)]]

    if pipeline == shared_constants.PIPELINE_AGGR:
        run_identifier = "Run"
    else:
        run_identifier = "Sample"
    rows = [
        [" ".join([run_identifier, "ID"]), sample_properties.sample_id],
        [" ".join([run_identifier, "Description"]), sample_properties.sample_desc],
    ]
    rows.extend(chemistry_row)

    if sample_data.summary.get("spatial_slide_info", None) is not None:
        rows.append(["Slide Serial Number", sample_data.summary["spatial_slide_info"]])

    if pipeline in shared_constants.PIPELINE_COUNT and not sample_properties.is_spatial:
        rows.append(["Include introns", str(sample_properties.include_introns)])

    if isinstance(sample_properties, ExtendedCountSampleProperties):
        if sample_properties.target_set:
            if sample_data.targeting_method == TARGETING_METHOD_TL:
                rows.append(["Probe Set Name", sample_properties.target_set])
            else:
                rows.append(["Target Panel Name", sample_properties.target_set])
            rows.append(["Number of Genes Targeted", sample_data.summary["num_genes_on_target"]])
        if (
            hasattr(sample_properties, "feature_ref_path")
            and sample_properties.feature_ref_path is not None
        ):
            rows.append(["Feature Reference", os.path.basename(sample_properties.feature_ref_path)])
        if sample_properties.disable_ab_aggregate_detection:
            rows.append(["Antibody Aggregate Filtering", "Disabled"])

    # Find references in the summary
    if isinstance(sample_properties, AggrCountSampleProperties) and not sample_properties.genomes:
        rows.append(
            [
                cr_constants.REFERENCE_TYPE,
                "Not applicable for aggr with feature barcoding-only samples",
            ]
        )
    elif pipeline in [shared_constants.PIPELINE_AGGR, shared_constants.PIPELINE_REANALYZE]:
        genomes = sample_properties.genomes
        if genomes is not None:
            rows.append(
                [
                    cr_constants.REFERENCE_TYPE,
                    get_ref_name_from_genomes(genomes),
                ]
            )
    else:
        reference_metric_prefixes = [
            cr_constants.REFERENCE_METRIC_PREFIX,
            vdj_constants.REFERENCE_METRIC_PREFIX,
        ]
        # Find all references in the summary
        for prefix in reference_metric_prefixes:
            ref_type_key = f"{prefix}{cr_constants.REFERENCE_TYPE_KEY}"
            if ref_type_key in sample_data.summary:
                ref_type = sample_data.summary[ref_type_key]

                ref_version_key = f"{prefix}{cr_constants.REFERENCE_VERSION_KEY}"
                if ref_version_key in sample_data.summary:
                    ref_version = f"-{sample_data.summary.get(ref_version_key)}"
                else:
                    ref_version = ""

                ref_name_key = f"{prefix}{cr_constants.REFERENCE_GENOMES_KEY}"
                if ref_name_key in sample_data.summary:
                    ref_name = sample_data.summary.get(ref_name_key)
                    if isinstance(ref_name, list):
                        ref_name = get_ref_name_from_genomes(ref_name)

                    rows.append([ref_type, f"{ref_name}{ref_version}"])

    # add pipeline version
    rows.append(["Pipeline Version", sample_properties.version])

    # add reorientation mode description to Sample table but not for aggr
    if (
        hasattr(sample_properties, "is_spatial")
        and sample_properties.is_spatial
        and pipeline not in shared_constants.PIPELINE_AGGR
    ):
        if sample_properties.reorientation_mode in ("rotation", "rotation+mirror"):
            rows.append(["Image Reorientation", "On"])
        else:
            rows.append(["Image Reorientation", "Off"])
        if sample_properties.loupe_alignment_file is None:
            rows.append(["Loupe Manual Alignment", "None"])
        else:
            rows.append(
                [
                    "Loupe Manual Alignment",
                    os.path.basename(sample_properties.loupe_alignment_file),
                ]
            )

    # add filter_probes mode description to the Sample table
    if sample_data.targeting_method == TARGETING_METHOD_TL:
        rows.append(["Filter Probes", "Off" if sample_properties.filter_probes is False else "On"])

    # add cytassist metadata
    if (
        isinstance(sample_properties, ExtendedCountSampleProperties)
        and sample_properties.is_spatial
        and sample_properties.cytassist_run_properties
    ):
        rows.extend(
            [
                [
                    "CytAssist Run Name",
                    sample_properties.cytassist_run_properties.cytassist_run_name,
                ],
                [
                    "CytAssist Instrument Serial ID",
                    sample_properties.cytassist_run_properties.cytassist_instrument_serial,
                ],
                [
                    "CytAssist Instrument Software Version",
                    sample_properties.cytassist_run_properties.cytassist_instrument_software_version,
                ],
            ]
        )

    pipeline_info = {
        "header": ["Sample"],
        "rows": rows,
    }
    to_return = {"pipeline_info_table": pipeline_info}

    # We want to alarm if ARC is used in GEX chemistry.
    if metadata is not None and species_list is not None:
        alarms = metadata.gen_metric_list(
            sample_data.summary, [H5_CHEMISTRY_DESC_KEY, "inconsistent_throughput"], species_list
        )
        new_alarms = [metric.alarm_dict for metric in alarms if metric.alarm_dict]
        if new_alarms:
            to_return[ALARMS] = new_alarms
    if (
        _get_chemistry_description(sample_data=sample_data, sample_properties=sample_properties)
        == ARC_V1_MULTIOME_DESCRIPTION
    ):
        to_return.setdefault(ALARMS, []).append(
            {
                "formatted_value": None,
                "title": "Unsupported workflow used",
                "message": "Multiome Gene Expression only analysis is not a supported workflow. Results may vary.",
                "level": WARNING_THRESHOLD,
            }
        )

    # If --aligner was used to force the aligner warn the user
    if (
        isinstance(sample_properties, CountSampleProperties)
        and sample_properties.aligner
        and sample_properties.is_spatial
    ):
        if ALARMS not in to_return:
            to_return[ALARMS] = []
        aligner_info_alarm = {
            "formatted_value": None,
            "title": "Force Aligner Used",
            "message": f"The --aligner option was used to set the sequencing read aligner to {sample_properties.aligner}. Incorrect usage of this option will lead to unusable data and low mapping metrics. Please contact support@10xgenomics.com for any further questions.",
            "level": INFO_THRESHOLD,
        }
        to_return[ALARMS].extend([aligner_info_alarm])

    # pattern downsampling
    if (
        isinstance(sample_properties, CountSampleProperties)
        and sample_properties.v1_pattern_fix
        and sample_properties.is_spatial
    ):
        if ALARMS not in to_return:
            to_return[ALARMS] = []
        downsample_info_alarm = {
            "formatted_value": None,
            "title": "UMI Downsampling Performed",
            "message": "The --v1-filtered-fbm option was used to downsample affected spots.",
            "level": INFO_THRESHOLD,
        }
        to_return[ALARMS].extend([downsample_info_alarm])

    # If default layout was used for Visium HD
    if (
        isinstance(sample_properties, CountSampleProperties)
        and sample_properties.default_layout
        and sample_properties.is_visium_hd
    ):
        if ALARMS not in to_return:
            to_return[ALARMS] = []
        default_layout_alarm = {
            "formatted_value": None,
            "title": "Default Layout Used",
            "message": "Visium HD slide ID was not provided. Analysis was performed using a generic slide layout, and will result in incorrect assignment of barcodes under tissue and misalignment of microscope image to UMI data. For best results, re-run with correct slide ID.",
            "level": WARNING_THRESHOLD,
        }
        to_return[ALARMS].extend([default_layout_alarm])

    # If Slide ID override was used for Visium HD
    if (
        isinstance(sample_properties, CountSampleProperties)
        and sample_properties.override_id
        and sample_properties.is_visium_hd
        and sample_properties.loupe_alignment_file is None
    ):
        if ALARMS not in to_return:
            to_return[ALARMS] = []
        override_id_alarm = {
            "formatted_value": None,
            "title": "Override ID Used",
            "message": "Visium HD slide ID entered on the CytAssist instrument was overridden by the inputs provided to the Space Ranger pipeline. Confirm that the slide ID matches the slide used to generate the dataset.",
            "level": WARNING_THRESHOLD,
        }
        to_return[ALARMS].extend([override_id_alarm])

    # If Slide ID in Loupe does not match slide ID from cytassist metadata
    if (
        isinstance(sample_properties, CountSampleProperties)
        and sample_properties.slide_id_mismatch  # TODO: FIX THIS
        and sample_properties.is_visium_hd
    ):
        if ALARMS not in to_return:
            to_return[ALARMS] = []
        slide_id_mismatch_alarm = {
            "formatted_value": None,
            "title": "Slide ID Mismatch",
            "message": "Visium HD slide ID entered on the CytAssist instrument was overridden by the slide ID provided to the Loupe Manual Aligner. Confirm that the slide ID matches the slide used to generate the dataset.",
            "level": WARNING_THRESHOLD,
        }
        to_return[ALARMS].extend([slide_id_mismatch_alarm])

    if (
        isinstance(sample_properties, CountSampleProperties)
        and sample_properties.itk_error_string
        and sample_properties.is_visium_hd
    ):
        if ALARMS not in to_return:
            to_return[ALARMS] = []
        itk_error_alarm = {
            "formatted_value": None,
            "title": "Tissue Registration",
            "message": "CytAssist to microscope image registration metrics indicate a possible poor alignment. "
            "Review your registration carefully. If automated registration performs poorly, perform manual alignment via "
            "<a href='https://www.10xgenomics.com/support/software/loupe-browser/latest' target='_blank' title='Loupe Browser' rel='noopener noreferrer'>Loupe Browser</a>.",
            "level": WARNING_THRESHOLD,
        }
        to_return[ALARMS].extend([itk_error_alarm])

    # If Isotype Normalization was used
    if sample_data.summary.get("ANTIBODY_isotype_normalized") is not None:
        rows.append(["Isotype Normalization", sample_data.summary["ANTIBODY_isotype_normalized"]])
        if sample_data.summary.get("ANTIBODY_isotype_normalized") == "Off":
            isotype_norm_info_alarm = {
                "formatted_value": None,
                "title": "Isotype Normalization",
                "message": f"Feature reference CSV does not contain isotype controls. Space Ranger did not perform {ALT_LINK_HELP_TEXT}.",
                "level": INFO_THRESHOLD,
            }
            if ALARMS not in to_return:
                to_return[ALARMS] = []
                to_return[ALARMS].extend([isotype_norm_info_alarm])
    # If loupe alignment file contains redundant information: image registration info when only Cytassist image provided
    to_return = _display_loupe_warning(sample_properties, to_return)

    return to_return


def create_table_with_alarms(
    table_key, title, metric_keys, alarm_keys, metadata, sample_data, species_list
):
    """Sequencing info for GEX."""
    if sample_data is None or sample_data.summary is None:
        return None

    data_dict = {}
    # Account for imaging block not having metrics
    if len(metric_keys):
        metrics = metadata.gen_metric_list(sample_data.summary, metric_keys, species_list)
        if metrics:
            # Not all requested metrics will appear, and we only generate help text for those
            # that show up in the table
            observed_keys = {x.parent_metric_info.name for x in metrics}
            filtered_keys = [x for x in metric_keys if x in observed_keys]
            data_dict["help"] = {
                "title": title,
                "data": metadata.gen_metric_helptext(filtered_keys),
            }
            data_dict["table"] = {
                "rows": [[metric.name, metric.value_string] for metric in metrics]
            }
    else:
        data_dict["help"] = {"title": title, "data": []}
        data_dict["table"] = {}
    if not data_dict:
        return None

    result = {table_key: data_dict}

    # Alerts.
    if alarm_keys:
        alarms = metadata.gen_metric_list(sample_data.summary, alarm_keys, species_list)
        # If a metric is from a barnyard and the cumulative version of the metric should be tested
        # we do not check the non-cumulative metrics
        alarms = [
            x
            for x in alarms
            if not (
                x.is_barnyard and x.parent_metric_info.include_cumulative and not x.is_cumulative
            )
        ]
        new_alarms = [metric.alarm_dict for metric in alarms if metric.alarm_dict]

        if new_alarms:
            result[ALARMS] = new_alarms

    return result


def sequencing_table(metadata, sample_data, species_list, is_targeted=False, is_hd=False):
    """Sequencing info for GEX."""
    metric_keys = SEQUENCING_METRIC_KEYS if not is_targeted else TARGETED_SEQUENCING_METRIC_KEYS
    if is_hd:
        metric_keys.append("filtered_bcs_conf_mapped_barcoded_reads_cum_frac")
    return create_table_with_alarms(
        "sequencing",
        "Sequencing",
        metric_keys,
        SEQUENCING_ALARM_KEYS,
        metadata,
        sample_data,
        species_list,
    )


def feature_barcode_sequencing_table(metadata, sample_data, species_list, feature_barcode):
    metric_keys = [f"{feature_barcode}_{i}" for i in SEQUENCING_METRIC_KEYS]
    alarm_keys = [f"{feature_barcode}_{i}" for i in SEQUENCING_ALARM_KEYS]

    return create_table_with_alarms(
        f"{feature_barcode.upper()}_sequencing",
        f"{FB_DISPLAY_NAME[feature_barcode]} Sequencing",
        metric_keys,
        alarm_keys,
        metadata,
        sample_data,
        species_list,
    )


def feature_barcode_application_table(metadata, sample_data, species_list, feature_barcode):
    """Feature barcoding application metric."""
    return create_table_with_alarms(
        f"{feature_barcode.upper()}_application",
        f"{FB_DISPLAY_NAME[feature_barcode]} Application",
        FB_APP_METRIC_KEYS[feature_barcode],
        FB_APP_METRIC_KEYS[feature_barcode],
        metadata,
        sample_data,
        species_list,
    )


def mapping_table(metadata, sample_data, species_list):
    """Mapping info table."""
    if sample_data.is_targeted():
        # Post library targeting
        if sample_data.targeting_method == TARGETING_METHOD_HC:
            alarm_keys = TARGETED_MAPPING_ALARMS
        # template ligation targeting
        elif sample_data.targeting_method == TARGETING_METHOD_TL and not sample_data.is_visium_hd:
            alarm_keys = TEMP_LIG_MAPPING_ALARMS
        elif sample_data.targeting_method == TARGETING_METHOD_TL and sample_data.is_visium_hd:
            alarm_keys = HD_TEMP_LIG_MAPPING_ALARMS
    else:
        alarm_keys = MAPPING_ALARMS

    return create_table_with_alarms(
        "mapping",
        "Mapping",
        (
            MAPPING_KEYS
            if "targeting_method" not in sample_data.summary
            or sample_data.targeting_method == TARGETING_METHOD_HC
            else TEMP_LIG_MAPPING_KEYS
        ),
        alarm_keys,
        metadata,
        sample_data,
        species_list,
    )


def cell_or_spot_calling_table(
    metadata,
    sample_data,
    sample_properties,
    species_list,
    metric_keys,
    alarm_keys,
):
    """Cell calling data (table and plot)."""
    # TODO: Barnyard not currently in spatial
    is_barnyard = len(species_list) > 1
    if is_barnyard:
        metric_keys.insert(0, FILTERED_BCS_TRANSCRIPTOME_UNION)

    if sample_properties.is_spatial and not sample_properties.is_visium_hd:
        tbl_name = "Spots"
    elif sample_properties.is_spatial and sample_properties.is_visium_hd:
        tbl_name = "Squares"
    else:
        tbl_name = "Cells"

    table_dict = create_table_with_alarms(
        "cells", tbl_name, metric_keys, alarm_keys, metadata, sample_data, species_list
    )
    if table_dict is None:
        return None

    # add image
    if sample_properties.is_spatial:
        return table_dict
    else:
        # the data we are interested in is in data_dict["cells"]
        data_dict = table_dict["cells"]
        to_return = {}
        # Be sure to bubble up alarms
        if ALARMS in table_dict:
            to_return[ALARMS] = table_dict[ALARMS]
        chart = get_empty_rank_plot()
        knee_plot = cr_webshim.plot_barcode_rank(chart, sample_properties, sample_data)
        if knee_plot:
            data_dict["barcode_knee_plot"] = knee_plot
            data_dict["help"]["data"] = data_dict["help"]["data"] + RANK_PLOT_HELP
        to_return["cells"] = data_dict
        return to_return


def summary_image_table(
    metadata,
    sample_data,
    species_list,
    metric_keys,
    alarm_keys,
    zoom_images,
    regist_images=None,
):
    """Image block for spatial samples."""
    tbl_name = "Image"
    # Keeping this core function just in case we want to throw image based alarms in the future
    table_dict = create_table_with_alarms(
        "image", tbl_name, metric_keys, alarm_keys, metadata, sample_data, species_list
    )
    if table_dict is None:
        return None
    # add tissue detection image
    table_dict["image"]["zoom_images"] = zoom_images
    # add tissue detection help
    table_dict["image"]["help"]["data"] = ZOOM_IMAGE_HELP
    if regist_images:
        # add cytassist tissue registration image
        table_dict["image"]["regist_images"] = regist_images
        # add cytassist tissue registration help
        table_dict["image"]["help"]["data"].append(REGISTRATION_IMAGE_HELP)
    return table_dict


def batch_correction_table(metadata, sample_data, species_list):
    metric_keys = [
        "batch_effect_score_before_correction",
        "batch_effect_score_after_correction",
    ]

    return create_table_with_alarms(
        "batch_correction",
        "Chemistry Batch Correction",
        metric_keys,
        None,
        metadata,
        sample_data,
        species_list,
    )


def normalization_summary_table(metadata, sample_data, sample_properties):
    """Report normalization metrics in aggr.

    The trick here is to use the batch prefix as the
    a species/genome prefix, and define these metric as species_specific.

    TODO: the above trick doesn't generate good web summary if it's a barnyard aggr sample.
    """
    if not isinstance(sample_properties, AggrCountSampleProperties):
        return None
    metric_keys = [
        "pre_normalization_total_reads",
        "post_normalization_total_reads",
        "pre_normalization_multi_transcriptome_total_raw_reads_per_filtered_bc",
        "post_normalization_multi_transcriptome_total_raw_reads_per_filtered_bc",
    ]

    alarm_keys = []
    batches = sample_properties.agg_batches

    table = create_table_with_alarms(
        "normalization_summary",
        "Normalization Summary",
        metric_keys,
        alarm_keys,
        metadata,
        sample_data,
        batches,
    )
    # reorganize the table (from long into wide)
    # Intended output:
    #    \t Total Number of Reads \t Mean Reads per Cell
    # Pre-Normalization \t XXXXX \t YYYYY
    # Post-Normalization \t ZZZZ \t RRRRR
    if (
        table is not None
        and "normalization_summary" in table
        and "table" in table["normalization_summary"]
        and "rows" in table["normalization_summary"]["table"]
    ):
        PRE_NORM = "Pre-Normalization"
        POST_NORM = "Post-Normalization"
        data = {}
        for row_title, row_value in table["normalization_summary"]["table"]["rows"]:
            if PRE_NORM in row_title:
                metric_name = row_title.split(PRE_NORM + " ")[1]
                if PRE_NORM not in data:
                    data[PRE_NORM] = {}
                if metric_name not in data[PRE_NORM]:
                    data[PRE_NORM][metric_name] = [row_value]
            if POST_NORM in row_title:
                metric_name = row_title.split(POST_NORM + " ")[1]
                if POST_NORM not in data:
                    data[POST_NORM] = {}
                if metric_name not in data[POST_NORM]:
                    data[POST_NORM][metric_name] = [row_value]
        metrics = data[PRE_NORM].keys()
        new_rows = []
        for row_name in [PRE_NORM, POST_NORM]:
            row_content = [row_name]
            for metric in metrics:
                row_content.append(data[row_name][metric][0])
            new_rows.append(row_content)

        table["normalization_summary"]["table"]["header"] = [" "] + list(metrics)
        table["normalization_summary"]["table"]["rows"] = new_rows

    return table


# reorganize the table (from long into wide)
# Intended output:
# Sample ID \t Fraction of Reads Kept \t Pre-Normalization Total Reads per Cell \t Pre-Normalization Confidently Mapped Barcoded Reads per Cell
# samp1 \t XXXXX \t YYYYY \t AAAAA
# samp2 \t ZZZZ \t RRRRR \t BBBBB
def split_table_by_samples(table, table_name):
    """Reorganize the table (from long into wide), one row for each sample id."""
    if (
        table is not None
        and table_name in table
        and "table" in table[table_name]
        and "rows" in table[table_name]["table"]
    ):
        data = {}
        all_metrics = []
        all_samples = []
        for row in table[table_name]["table"]["rows"]:
            metric_name_with_sample = row[0]
            metric_value = row[1]
            if "(" in metric_name_with_sample:
                metric_name = metric_name_with_sample.split(" (")[0]
                sample_name = metric_name_with_sample.split(" (")[1][0:-1]
            else:
                # in reanalyze
                metric_name = metric_name_with_sample
                sample_name = "NA"
            if sample_name not in all_samples:
                all_samples.append(sample_name)
            if sample_name not in data:
                data[sample_name] = {}
            if metric_name not in data[sample_name]:
                data[sample_name][metric_name] = metric_value
            if metric_name not in all_metrics:
                all_metrics.append(metric_name)

        new_rows = []
        for row_name in all_samples:
            row_content = [row_name]
            for metric_name in all_metrics:
                row_content.append(data[row_name].get(metric_name, float("NaN")))
            new_rows.append(row_content)

        table[table_name]["table"]["header"] = ["Sample ID"] + all_metrics
        table[table_name]["table"]["rows"] = new_rows


def aggregation_by_sample_table(metadata, sample_data, sample_properties):
    """Report per-sample normalization metrics in aggr."""
    if not isinstance(sample_properties, AggrCountSampleProperties):
        return None

    alarm_keys = [
        "lowest_frac_reads_kept",
    ]
    batches = sample_properties.agg_batches
    table = create_table_with_alarms(
        "aggregation_by_sample",
        "Gene Expression Aggregation by Sample",
        AGGREGATION_METRIC_KEYS,
        alarm_keys,
        metadata,
        sample_data,
        batches,
    )
    split_table_by_samples(table, "aggregation_by_sample")
    return table


def feature_barcode_aggregation_table(metadata, sample_data, sample_properties, feature_barcode):
    if not isinstance(sample_properties, AggrCountSampleProperties):
        return None
    metric_keys = [f"{feature_barcode}_{i}" for i in AGGREGATION_METRIC_KEYS]
    batches = sample_properties.agg_batches
    table = create_table_with_alarms(
        f"{feature_barcode}_aggregation",
        f"{FB_DISPLAY_NAME[feature_barcode]} Aggregation by Sample",
        metric_keys,
        None,
        metadata,
        sample_data,
        batches,
    )
    split_table_by_samples(table, f"{feature_barcode}_aggregation")
    return table


def _display_loupe_warning(sample_properties, to_return):
    """Display warning to websummary.

    if loupe file contains redundant information
    """
    if (
        isinstance(sample_properties, CountSampleProperties)
        and sample_properties.redundant_loupe_alignment
    ):
        if ALARMS not in to_return:
            to_return[ALARMS] = []
        loupe_info_alarm = {
            "formatted_value": None,
            "title": "Unused Loupe alignment info",
            "message": "Loupe alignment file contains unused tissue registration info because the microscope image used to produce the alignment file was not input to Space Ranger.\
                 If this was intentional, then you can proceed with analyzing your data. Please contact support@10xgenomics.com for any further questions.",
            "level": INFO_THRESHOLD,
        }
        to_return[ALARMS].extend([loupe_info_alarm])
    return to_return
