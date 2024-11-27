# Copyright (c) 2018 10X Genomics, Inc. All rights reserved.
"""Tools for managing metrics with metadata attached."""
from __future__ import annotations

# Copied from cellranger-atac (commit 5ccda578ed71e24289c6ea9bb6dec4d9de5d11dc)
import csv
import os
from collections import OrderedDict

import pandas as pd
from six import ensure_str

from cellranger.rna.library import add_multi_prefix, add_species_prefix
from cellranger.vdj.chain_types import load_chains_from_chain_type
from cellranger.webshim.data import SampleData
from tenkit.safe_json import desanitize_value

VALID_THRESHOLD = "VALID"
WARNING_THRESHOLD = "WARN"
ERROR_THRESHOLD = "ERROR"
INFO_THRESHOLD = "INFO"

# We enforce that spaceranger, VDJ and cellranger have the same fields in their metrics files
# to avoid errors if shared code expects a field but it is missing.  Space Ranger
# and VDJ (8/6/2019) do not currently support barnyard experiments, but we still add include_cumulative and
# is_species_specific to it.
required_cols = {
    "summary",
    "category",
    "full_name",
    "alert_warn_detail",
    "alert_error_detail",
    "acceptable",
    "targeted",
    "alert_warn_detail_cs",
    "alert_error_detail_cs",
    "acceptable_cs",
    "targeted_cs",
    "evaluate_type",
    "format_type",
    "is_species_specific",
    "include_cumulative",
    "alert_warn_name",
    "alert_error_name",
    "help_description",
}


def get_maybe_nested_key(data: dict, key: str):
    """Given a key which may be nested using the dot notation, return the value.

    This function will return None if the key is not found. The maximum depth of nesting
    is expected to be very small (~2) so we use recursion.
    """
    if key in data:
        return data[key]
    elif "." in key:
        left, right = key.split(".", 1)
        if left in data:
            return get_maybe_nested_key(data[left], right)

    return None


def _validate_expected_keys_present(keys):
    missing_keys = required_cols - keys
    assert len(missing_keys) == 0, "Expected keys in Metrics CSV file not found."


def convert_value(func):
    """Decorator to convert from strings "NaN", "Inf", "-Inf" back to floats.

    This is used to enable comparisons with numeric types.

    Args:
        func: A function with value, target as args
    """

    def _convert_and_test(value, target):
        value = desanitize_value(value)
        return func(value, target)

    return _convert_and_test


def load_metric_data(metric_csv_file):
    """Loads the metric metadata as a dictionary by parsing the input csv file.

    Entries can point to another entry. For example 'alert_warn_detail_cs' field
    could just say 'alert_warn_detail', meaning we want to use the content in
    'alert_warn_detail' for 'alert_warn_detail_cs' as well. Only a single level of
    redirection is supported. i.e You should not have A -> B -> C, instead use A->C
    and B->C. No circular references too (i.e A->B and B->A). The code does not
    explicitly check for these, but these simple rules should be followed when building
    the metrics csv.
    """
    data = pd.read_csv(
        ensure_str(metric_csv_file), index_col=0, comment="#", float_precision="high"
    )

    # Code to resolve single level of redirections
    headers = set(data.keys().tolist())
    _validate_expected_keys_present(headers)
    for i, row in data.iterrows():
        for k, v in row.items():
            if v in headers:
                data.at[i, k] = row[v]
    # Return the data as a dictionary
    metric_data = OrderedDict()
    for key in data.index.values:
        if key in metric_data:
            raise OSError("Metrics file cannot contain duplicated keys")
        metric_data[key] = data.loc[key]

        # Makes empty values None instead of NaN where appropriate
        # if do in context of entire dataframe vs. one key at a time can still get NaNs
        metric_data[key] = metric_data[key].where(pd.notnull(metric_data[key]), None)
    return metric_data


def output_metrics_csv_from_annotations(
    metric_annotations, sample_data, out_fname, species_list=None
):
    """Produces a metrics csv file from.

    :param metric_annotations: An instance of MetricAnnotations that specifies which metrics will
    be output in the summary file based on the boolean `summary` column in the associated CSV file.
    :param sample_data: An instance of SampleData containing the necessary metrics we're looking for
    :param out_fname: Output csv file name
    :param species_list: A list of strings with each species (typically from the genomes of the h5)
    :return: None
    """
    assert isinstance(metric_annotations, MetricAnnotations)
    assert isinstance(sample_data, SampleData)
    metric_annotations.output_metrics_csv(sample_data, out_fname, species_list)


class MetricAnnotations:
    def __init__(self, metric_file=None, intron_mode_alerts=False):
        """Load in metric information from the associated csv file."""
        metric_file = "metrics.csv" if metric_file is None else metric_file
        file_path = os.path.join(os.path.dirname(__file__), metric_file)
        self.metric_data = load_metric_data(file_path)

        if intron_mode_alerts:
            file_path = os.path.join(os.path.dirname(__file__), "intron_mode_metrics.csv")
            self._override_metric_settings(file_path)
        self.source_file = file_path

    def gen_metric(
        self, key, value, species=None, is_barnyard=False, debug=True, is_cumulative=False
    ):
        """Returns a single metric object for the given key and value.

         Alerts are raised when metric falls outside
        normal range, which depends on debug status.

        if it's barnyard sample, add species suffix to name and alert_name.
        """
        metric_info = self.metric_data[key]
        name = metric_info.full_name
        alert_name_suffix = ""
        if species:
            key = f"{species}_{key}"
            name += f" ({species})" if is_barnyard else ""
            alert_name_suffix += f" ({species})" if is_barnyard else ""

        # alarm ranges are dependent on debug, which indicates internal or customer-facing.
        alert_name_map = {
            WARNING_THRESHOLD: (
                metric_info.alert_warn_name + alert_name_suffix
                if metric_info.alert_warn_name
                else None
            ),
            ERROR_THRESHOLD: (
                metric_info.alert_error_name + alert_name_suffix
                if metric_info.alert_error_name
                else None
            ),
        }
        # This appears under detail when an alert is raised
        alert_detail_map = {
            WARNING_THRESHOLD: (
                metric_info.alert_warn_detail if debug else metric_info.alert_warn_detail_cs
            ),
            ERROR_THRESHOLD: (
                metric_info.alert_error_detail if debug else metric_info.alert_error_detail_cs
            ),
        }
        acceptable = metric_info.acceptable if debug else metric_info.acceptable_cs
        targeted = metric_info.targeted if debug else metric_info.targeted_cs

        return Metric(
            key,
            name,
            value,
            metric_info,
            alert_detail_map=alert_detail_map,
            acceptable=acceptable,
            targeted=targeted,
            evaluate_type=metric_info.evaluate_type,
            format_type=metric_info.format_type,
            category=metric_info.category,
            alert_name_map=alert_name_map,
            is_barnyard=is_barnyard,
            is_cumulative=is_cumulative,
        )

    def compile_summary_metrics(self, value_dict, species_list=None):
        """Processes a metrics dictionary and select summary metrics based on 2nd column in.

        metrics.csv for keys provided or all registered keys in metrics.csv

        Note: The ATAC version of this returns a dictionary with key/value pairs, but this
        version will return a list of tuples with formatted-name/value
        """
        keys = (key for key, value in self.metric_data.items() if value.summary)
        mets = self.gen_metric_list(value_dict, keys, species_list)
        return mets

    def gen_metric_helptext(self, keys):
        """Processes a metrics dictionary and generates helptext for keys if present in metrics.csv."""
        output = []
        for key in keys:
            if key in self.metric_data:
                metric_info = self.metric_data[key]
                if metric_info.help_description is not None:
                    full_name = metric_info.full_name
                    output.append([full_name, [metric_info.help_description]])
            else:
                print(f"{key} not found in registered metrics")
        return output

    def gen_metric_list(self, value_dict, keys, species_list, debug=True):
        """Returns a list of metric objects for the provided keys, using the value dictionary to get values for.

        each metric.  When metrics are species-specific, a list of species is required.  Alerts are raised when
        metrics fall outside normal ranges, which depend on debug status.
        """
        output = []
        is_barnyard = len(species_list) > 1
        is_antibody_only = len(species_list) == 0

        for key in keys:
            # A metric might be requested but not in the CSV for this product, we skip it.
            # e.g. `read2_bases_with_q30_frac` exists in the sequencing table for 5', but
            # not for spatial.
            if key not in self.metric_data:
                continue
            metric_info = self.metric_data[key]
            if metric_info.is_species_specific:
                if (species_list is None or not species_list) and not is_antibody_only:
                    raise ValueError("Must provide a species list for species-specific metrics")

                if is_barnyard and metric_info.include_cumulative:
                    subkey = add_multi_prefix(key)
                    value = get_maybe_nested_key(value_dict, subkey)
                    if value is not None:
                        output.append(
                            self.gen_metric(
                                key,
                                value,
                                debug=debug,
                                is_barnyard=is_barnyard,
                                is_cumulative=True,
                            )
                        )

                for species in species_list:
                    subkey = add_species_prefix(species, key)
                    value = get_maybe_nested_key(value_dict, subkey)
                    if value is not None:
                        output.append(
                            self.gen_metric(
                                key,
                                value,
                                species,
                                is_barnyard,
                                debug,
                                is_cumulative=False,
                            )
                        )
            else:
                value = get_maybe_nested_key(value_dict, key)
                if value is not None:
                    output.append(self.gen_metric(key, value, debug=debug))

        return output

    def output_metrics_csv(self, sample_data, filename, species_list=None):
        assert isinstance(sample_data, SampleData)
        mets = self.compile_summary_metrics(sample_data.summary, species_list=species_list)
        header = [x.name for x in mets]
        values = [x.value for x in mets]
        with open(filename, "w") as f:
            writer = csv.writer(f, lineterminator="\n")
            writer.writerows([header, values])

    def _override_metric_settings(self, filename):
        """In order to enable metric settings to inherit from others, we can load up a second file with.

        settings that add or extend the previous ones.

        Args:
            filename: Filename of metrics.csv with override settings

        Returns:
            None
        """
        new_metric_data = load_metric_data(filename)
        self.metric_data.update(new_metric_data)


class SpatialMetricAnnotations(MetricAnnotations):
    def __init__(self, intron_mode_alerts=False):
        super().__init__("spatial_metrics.csv", intron_mode_alerts=intron_mode_alerts)


class SpatialAggrMetricAnnotations(SpatialMetricAnnotations):
    def __init__(self):
        super().__init__()
        file_path = os.path.join(os.path.dirname(__file__), "spatial_aggr_metrics.csv")
        self._override_metric_settings(file_path)


class SpatialTargetedMetricAnnotations(SpatialMetricAnnotations):
    def __init__(self):
        super().__init__()
        file_path = os.path.join(os.path.dirname(__file__), "spatial_targeted_metrics.csv")
        self._override_metric_settings(file_path)


class SpatialTargetedAggrMetricAnnotations(SpatialAggrMetricAnnotations):
    def __init__(self):
        super().__init__()
        file_path = os.path.join(os.path.dirname(__file__), "spatial_targeted_aggr_metrics.csv")
        self._override_metric_settings(file_path)


class SpatialTemplateLigationAggrMetricAnnotations(SpatialAggrMetricAnnotations):
    def __init__(self):
        super().__init__()
        file_path = os.path.join(
            os.path.dirname(__file__), "spatial_template_ligation_aggr_metrics.csv"
        )
        self._override_metric_settings(file_path)


class SpatialTemplateLigationMetricAnnotations(SpatialMetricAnnotations):
    def __init__(self):
        super().__init__()
        file_path = os.path.join(os.path.dirname(__file__), "spatial_template_ligation_metrics.csv")
        self._override_metric_settings(file_path)


class SpatialHDTemplateLigationMetricAnnotations(SpatialTemplateLigationMetricAnnotations):
    def __init__(self):
        super().__init__()
        file_path = os.path.join(
            os.path.dirname(__file__), "spatial_hd_template_ligation_metrics.csv"
        )
        self._override_metric_settings(file_path)


class TargetedAggrMetricAnnotations(MetricAnnotations):
    def __init__(self):
        super().__init__()
        file_path = os.path.join(os.path.dirname(__file__), "targeted_aggr_metrics.csv")
        self._override_metric_settings(file_path)


class TemplateLigationAggrMetricAnnotations(MetricAnnotations):
    def __init__(self):
        super().__init__()
        file_path = os.path.join(os.path.dirname(__file__), "template_ligation_aggr_metrics.csv")
        self._override_metric_settings(file_path)


class TargetedMetricAnnotations(MetricAnnotations):
    def __init__(self):
        super().__init__()
        file_path = os.path.join(os.path.dirname(__file__), "targeted_metrics.csv")
        self._override_metric_settings(file_path)


class LTMetricAnnotations(MetricAnnotations):
    def __init__(self, intron_mode_alerts=False):
        super().__init__(intron_mode_alerts=intron_mode_alerts)
        file_path = os.path.join(os.path.dirname(__file__), "lt_metrics.csv")
        self._override_metric_settings(file_path)


class Metric:
    """Contains metadata about a single metric along with methods for evaluating its value with respect to that.

    metadata.
    """

    def __init__(
        self,
        key,
        name,
        value,
        parent_metric_info,
        alert_detail_map=None,
        acceptable=None,
        targeted=None,
        evaluate_type=None,
        format_type=None,
        category=None,
        alert_name_map=None,
        is_barnyard=False,
        is_cumulative=False,
    ):
        """:param key: Metric identifier.

        :param name: Full Metric Name
        :param value: The value of the metric
        :param parent_metric_info: Information from the row in the metrics.csv relevant to this metric.
        :param alert_detail_map: Details text for both warning and error alerts
        :param acceptable: Value to determines the error level
        :param targeted: Value to determine the warning level
        :param evaluate_type: lt/gt/range/exists/is_true/is_false/None
        :param format_type: int/percentage/flt/float
        :param category: A category name this metric belongs to (e.g. Mapping)
        :param alert_name_map: Alert text for both warning and error alerts
        :param is_barnyard: True if this metric was generated from a barnyard experiment
        :param is_cumulative: True if this metric is accumulated over multiple genomes
        """
        self.key = key
        self.name = name
        self.value = value
        self.alert_detail_map = alert_detail_map
        self.is_barnyard = is_barnyard
        self.is_cumulative = is_cumulative
        self.parent_metric_info = parent_metric_info
        self.acceptable = self._cast_threhold(acceptable)
        self.targeted = self._cast_threhold(targeted)
        self.alert_name_map = alert_name_map

        function_map = {
            None: lambda x, y: True,
            "lt": self._less_than,
            "gt": self._greater_than,
            "range": self._in_range,
            "exists": self._exists,
            "is_true": self._is_true,
            "is_false": self._is_false,
            "is_equal": self._is_equal,
            "is_not_equal": self._is_not_equal,
        }
        if evaluate_type not in function_map:
            raise ValueError(f"Unknown evaluation type: {evaluate_type}")
        self.evaluate_type = evaluate_type
        self.evaluation_function = function_map[evaluate_type]

        if format_type is None:
            format_type = "flat"
        if format_type not in ["flat", "float", "percentage", "int", "scientific"]:
            raise ValueError(f"Unknown format type: {format_type}")
        self.format_type = format_type

        if category is None:
            category = "General"
        self.category = category

    @staticmethod
    def _cast_threhold(t):
        try:
            t_str = str(t)
            if "-" in t_str and (not t_str.startswith("-")):
                res = [float(i.strip()) for i in t.split("-")]
            else:
                res = float(t)
        except (ValueError, TypeError):
            res = t
        return res

    def gen_metric_dict(self):
        threshold = self.threshold_type
        translate_dict = {
            VALID_THRESHOLD: "pass",
            WARNING_THRESHOLD: "warn",
            ERROR_THRESHOLD: "error",
        }
        return {
            "threshold": translate_dict[threshold],
            "metric": self.value_string,
            "name": self.name,
        }

    @property
    def alarm_dict(self):
        threshold = self.threshold_type
        if threshold == VALID_THRESHOLD:
            return {}
        return {
            "raw_value": self.value,
            "formatted_value": self.value_string,
            "raised": True,
            "parent": self.key,
            "title": self.alert_name_map[threshold],
            "message": self.alert_detail_map[threshold],
            "level": threshold,
            "test": "",
            "id": self.key,
        }

    # Evaluation methods
    @staticmethod
    @convert_value
    def _less_than(value, target):
        return value < target

    @staticmethod
    @convert_value
    def _greater_than(value, target):
        return value > target

    @staticmethod
    def _in_range(value, target):
        return (value > target[0]) and (value < target[1])

    @staticmethod
    def _exists(value, target):
        return value is not None

    @staticmethod
    def _is_true(value, target):
        return value is True

    @staticmethod
    def _is_equal(value, target):
        return value == target

    @staticmethod
    def _is_not_equal(value, target):
        return value != target

    @staticmethod
    def _is_false(value, target):
        return value is False

    @property
    def threshold_type(self):
        if self.targeted is None and self.acceptable is None:
            return VALID_THRESHOLD
        elif self.targeted is None:
            # Only acceptable - error if we don't meet them.
            if self.evaluation_function(self.value, self.acceptable):
                return VALID_THRESHOLD
            return ERROR_THRESHOLD
        elif self.acceptable is None:
            # Only targets - warn if we don't meet them.
            if self.evaluation_function(self.value, self.targeted):
                return VALID_THRESHOLD
            return WARNING_THRESHOLD
        elif self.evaluation_function(self.value, self.targeted):
            return VALID_THRESHOLD
        elif self.evaluation_function(self.value, self.acceptable):
            return WARNING_THRESHOLD
        else:
            return ERROR_THRESHOLD

    @property
    def color(self):
        if self.targeted is None and self.acceptable is None:
            # Meaningless, return grey hex code
            return "BEBEBE"
        threshold = self.threshold_type
        if threshold == VALID_THRESHOLD:
            return "B4FFB4"
        elif threshold == WARNING_THRESHOLD:
            return "FFFFB4"
        else:
            return "FFB4B4"

    @property
    def acceptable_string(self):
        return self._format_target_value(self.acceptable)

    @property
    def targeted_string(self):
        return self._format_target_value(self.targeted)

    @property
    def value_string(self):
        return self._format_value(self.value)

    def _format_target_value(self, value):
        if value is None:
            return ""

        if self.evaluate_type == "exists":
            return "Exists"
        if self.evaluate_type == "is_true":
            return "True"
        if self.evaluate_type == "is_false":
            return "False"

        if self.evaluate_type == "lt":
            return f"< {self._format_value(value)}"
        if self.evaluate_type == "gt":
            return f"> {self._format_value(value)}"
        if self.evaluate_type == "range":
            return f"{self._format_value(value[0])} - {self._format_value(value[1])}"
        raise AssertionError("unreachable - invalid evaluate_type")

    def _format_value(self, value):
        if value is None or value == "NaN":
            return "None"

        if self.format_type == "flat":
            return f"{value}"
        elif self.format_type == "float":
            return f"{value:,.2f}"
        elif self.format_type == "percentage":
            return f"{value:.1%}"
        elif self.format_type == "int":
            return f"{value:,.0f}"
        elif self.format_type == "scientific":
            return f"{value:.1e}"
        raise AssertionError("unreachable - invalid format_type")


class VDJMetricAnnotations:
    def __init__(self):
        """Load metric information from the associated csv file."""
        ############ TODO: CHANGE THIS ###########################
        file_path = os.path.join(os.path.dirname(__file__), "vdj_metrics.csv")
        self.metric_data = load_metric_data(file_path)

    def gen_metric(self, key, value, debug=False, chain=None):
        """Returns a single metric object for the given key and value.

         Alerts are raised when metric falls outside
        normal range, which depends on debug status.
        """
        metric_info = self.metric_data[key]
        name = metric_info.full_name

        key = key.format(chain=chain)
        name = name.format(chain=chain)

        alert_name_map = {
            WARNING_THRESHOLD: (
                metric_info.alert_warn_name.format(chain=chain)
                if metric_info.alert_warn_name
                else None
            ),
            ERROR_THRESHOLD: (
                metric_info.alert_error_name.format(chain=chain)
                if metric_info.alert_error_name
                else None
            ),
        }

        # alarm ranges are dependent on debug, which indicates internal or customer-facing.
        # This appears under detail when an alert is raised
        alert_warn_detail = (
            metric_info.alert_warn_detail if debug else metric_info.alert_warn_detail_cs
        )
        alert_error_detail = (
            metric_info.alert_error_detail if debug else metric_info.alert_error_detail_cs
        )
        alert_detail_map = {
            WARNING_THRESHOLD: alert_warn_detail.format(chain=chain) if alert_warn_detail else None,
            ERROR_THRESHOLD: alert_error_detail.format(chain=chain) if alert_error_detail else None,
        }
        acceptable = metric_info.acceptable if debug else metric_info.acceptable_cs
        targeted = metric_info.targeted if debug else metric_info.targeted_cs

        return Metric(
            key,
            name,
            value,
            metric_info,
            alert_detail_map=alert_detail_map,
            acceptable=acceptable,
            targeted=targeted,
            evaluate_type=metric_info.evaluate_type,
            format_type=metric_info.format_type,
            category=metric_info.category,
            alert_name_map=alert_name_map,
        )

    def gen_metric_list(self, value_dict, keys, debug=False, chain_type=None):
        """Returns a list of metric objects for the provided keys, using the value dictionary to get values for.

        each metric. Alerts are raised when metrics fall outside normal ranges, which depend on debug status.
        """
        output = []
        for key in keys:
            metric_info = self.metric_data[key]
            if metric_info.is_chain_specific:
                assert chain_type is not None, f"Got no chain type for {key}"
                for chain in load_chains_from_chain_type(chain_type):
                    full_key = key.format(chain=chain)
                    if full_key in value_dict:
                        output.append(
                            self.gen_metric(key, value_dict[full_key], debug=debug, chain=chain)
                        )
                    else:
                        print(f"{full_key} not found in metrics")
            elif key in value_dict:
                output.append(self.gen_metric(key, value_dict[key], debug=debug))
            else:
                print(f"{key} not found in metrics")
        return output

    def gen_metric_helptext(self, keys, chain_type=None):
        """Processes a metrics dictionary and generates helptext for keys if present in the metrics csv."""
        output = []
        for key in keys:
            if key in self.metric_data:
                metric_info = self.metric_data[key]
                if metric_info.is_chain_specific:
                    assert chain_type is not None, f"Got no chain type for {key}"
                    for chain in load_chains_from_chain_type(chain_type):
                        full_name = metric_info.full_name.format(chain=chain)
                        output += [[full_name, [metric_info.help_description.format(chain=chain)]]]
                elif metric_info.help_description is not None:
                    full_name = metric_info.full_name
                    output += [[full_name, [metric_info.help_description]]]
            else:
                print(f"{key} not found in registered metrics")
        return output
