#
# Copyright (c) 2025 10X Genomics, Inc. All rights reserved.
#

"""Structify cell annotation outs."""

import os
from dataclasses import asdict, dataclass

import martian

from cellranger.cr_io import hard_link, recursive_hard_link_dict

__MRO__ = """
struct CellTypeWebSummaryBundle(
    string          cell_annotation_model_family    "Cell annotation model family used.",
    json            cell_type_interactive_bar_chart "Interactive bar chart of cell type composition.",
    json            cell_types_box_plot             "Box plot of cell type UMI counts.",
    json            cell_types_umap_plot            "UMAP plot of cell types.",
    json            diffexp                         "Differential expression table.",
    json            cell_annotation_metrics         "Cell annotation metrics.",
    CellTypePdPlots cell_type_pd_plots              "Cell annotation PD plots.",
    json            cell_annotation_metadata        "Cell annotation metadata.",
    CellTypeInfo    cell_type_info                  "Cell type information.",
    json            cell_type_spatial_plot          "Spatial plot of cell types.",
    bool            aggregate_and_annotate          "Whether metacell aggregation was used before annotation.",
)


struct CellTypes(
    csv     cell_types,
    json.gz cell_annotation_results,
    csv     cell_annotation_differential_expression,
)

struct CellAnnotationLoupeInputs(
    string model_name,
    string cell_annotation_cloupe_name,
    json   cloupe_cas_types,
)

struct CellTypesAllMethods(
    map<CellTypes> cell_types,
)

struct WebSummaryCellTypeInputs(
    map<CellTypeWebSummaryBundle> cell_type_websummary_bundles,
)

stage STRUCTIFY_CELL_ANNOTATION_OUTS(
    in  string[]                      model_families,
    in  CellTypeWebSummaryBundle[]    cell_annotation_websummary_bundles,
    in  CellTypes[]                   all_cell_types,
    in  json[]                        cloupe_cas_types,
    in  string[]                      cell_annotation_cloupe_names,
    in  json[]                        cell_annotation_analysis_metrics_in,
    out map<json>                     cell_annotation_analysis_metrics,
    out map<json>                     cell_annotation_metrics,
    out CellAnnotationLoupeInputs[]   loupe_file_inputs,
    out WebSummaryCellTypeInputs      websummary_inputs,
    out WebSummaryCellTypeInputs      count_websummary_inputs,
    out CellTypesAllMethods           cell_types_all_methods,
    src py                            "stages/cas_cell_typing/structify_cell_annotation_outs",
)
"""


@dataclass
class CellAnnotationLoupeInputs:
    """Class corresponding to martian struct."""

    model_name: str | None = None
    cell_annotation_cloupe_name: str | None = None
    cloupe_cas_types: str | None = None


WEBSUMMARY_BUNDLE_NOT_PATHS = [
    "cell_annotation_model_family",
    "cell_type_info",
    "aggregate_and_annotate",
]


def main(args, outs):  # pylint: disable=too-many-locals
    assert (
        len(args.model_families)
        == len(args.cell_annotation_websummary_bundles)
        == len(args.all_cell_types)
        == len(args.cloupe_cas_types)
        == len(args.cell_annotation_cloupe_names)
        == len(args.cell_annotation_analysis_metrics_in)
    ), "Input vectors to stage are not of same length"

    cell_annotation_metrics_list = []
    cell_type_websummary_bundles_outs = {}
    for model_name, cell_annotation_websummary_bundle in zip(
        args.model_families, args.cell_annotation_websummary_bundles
    ):
        cell_type_websummary_bundle_dict = {}
        for key, flpath in cell_annotation_websummary_bundle.items():
            if key not in WEBSUMMARY_BUNDLE_NOT_PATHS and isinstance(flpath, dict):
                nested_dict = {}
                for nested_key, nested_flpath in flpath.items():
                    if nested_flpath:
                        nested_flname = os.path.basename(nested_flpath)
                        tmp_file_path = martian.make_path(
                            f"{model_name}_{key}_{nested_key}_{nested_flname}"
                        ).decode()
                        hard_link(nested_flpath, tmp_file_path)
                        nested_dict[nested_key] = tmp_file_path
                    else:
                        nested_dict[nested_key] = None
                cell_type_websummary_bundle_dict[key] = nested_dict
            elif key not in WEBSUMMARY_BUNDLE_NOT_PATHS:
                if flpath:
                    flname = os.path.basename(flpath)
                    tmp_file_path = martian.make_path(f"{model_name}_{key}_{flname}").decode()
                    hard_link(flpath, tmp_file_path)
                    cell_type_websummary_bundle_dict[key] = tmp_file_path
                else:
                    cell_type_websummary_bundle_dict[key] = None
            else:
                cell_type_websummary_bundle_dict[key] = flpath
        cell_annotation_metrics_list.append(
            cell_type_websummary_bundle_dict.get("cell_annotation_metrics")
        )
        if cell_type_websummary_bundle_dict:
            cell_type_websummary_bundles_outs[model_name] = cell_type_websummary_bundle_dict.copy()
    if cell_type_websummary_bundles_outs:
        outs.websummary_inputs = {
            "cell_type_websummary_bundles": cell_type_websummary_bundles_outs.copy()
        }
    else:
        outs.websummary_inputs = None

    tmp_cell_types_dict = {}
    for model_name, cell_types in zip(args.model_families, args.all_cell_types):
        if cell_types:
            tmp_cell_types_dict[model_name] = cell_types
    if tmp_cell_types_dict:
        outs.cell_types_all_methods = {"cell_types": recursive_hard_link_dict(tmp_cell_types_dict)}
    else:
        outs.cell_types_all_methods = None

    outs.loupe_file_inputs = []
    for model_name, cell_annotation_cloupe_name, cloupe_cas_types in zip(
        args.model_families, args.cell_annotation_cloupe_names, args.cloupe_cas_types
    ):
        tmp_file_path = None
        if cloupe_cas_types:
            tmp_file_path = martian.make_path(f"{model_name}_cloupe_cas_types.json").decode()
            hard_link(cloupe_cas_types, tmp_file_path)
        outs.loupe_file_inputs.append(
            asdict(
                CellAnnotationLoupeInputs(
                    model_name=model_name,
                    cell_annotation_cloupe_name=cell_annotation_cloupe_name,
                    cloupe_cas_types=tmp_file_path,
                )
            )
        )

    outs.cell_annotation_metrics = {}
    outs.cell_annotation_analysis_metrics = {}
    for model_name, cell_annotation_metrics, cell_annotation_analysis_metrics in zip(
        args.model_families,
        cell_annotation_metrics_list,
        args.cell_annotation_analysis_metrics_in,
    ):
        if cell_annotation_metrics:
            tmp_cell_annotation_metrics = martian.make_path(
                f"{model_name}_cell_annotation_metrics.json"
            ).decode()
            hard_link(cell_annotation_metrics, tmp_cell_annotation_metrics)
            outs.cell_annotation_metrics[model_name] = tmp_cell_annotation_metrics
        else:
            outs.cell_annotation_metrics[model_name] = None

        if cell_annotation_analysis_metrics:
            tmp_cell_annotation_analysis_metrics = martian.make_path(
                f"{model_name}_cell_annotation_analysis_metrics.json"
            ).decode()
            hard_link(cell_annotation_analysis_metrics, tmp_cell_annotation_analysis_metrics)
            outs.cell_annotation_analysis_metrics[model_name] = tmp_cell_annotation_analysis_metrics
        else:
            outs.cell_annotation_analysis_metrics[model_name] = None
