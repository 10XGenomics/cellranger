# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
"""Python versions of the React Components available in the web summary.

The Web Summary data is serialized as JSON and for simple objects we just pass
dictionaries around, but also use typed classes for some aspects of the
Components.
"""


from __future__ import annotations

import json
from dataclasses import dataclass

import numpy as np

import cellranger.websummary.plotly_tools as pltly
from cellranger.websummary.numeric_converters import _round_float, round_floats_in_list


class React10X:
    # pylint: disable=too-few-public-methods
    """A base class used to indicate data is related to React Components used in the websummary."""

    def __init__(self):
        pass

    def to_dict_for_json(self):
        """Convert this object to a dictionary suitable for json serialization.

        We currently serialize our Web Summary data with json.dumps
        which only understands a few built-in types.  To support serialization of the web summary
        data, types can implement this method to provide custom serialization of their values.

        Caution:
            Edits to this returned value will be reflected in the original instance
        """
        return vars(self)


class ReactComponentEncoder(json.JSONEncoder):
    # pylint: disable=too-few-public-methods,arguments-differ
    """An Encoder class to enable json.dumps to interpret React10X objects."""

    # pylint: disable=method-hidden
    def default(self, o):
        if isinstance(o, React10X):
            return o.to_dict_for_json()

        if isinstance(o, bytes):
            return o.decode()

        if isinstance(o, np.integer):
            return int(o)

        # we have an unhandled type -- this call generates an exception
        return json.JSONEncoder.default(self, o)


class ResourceCollection(React10X):
    """A class to store resources shared across multiple React components."""

    # Sets the maximum number of unique items we can store
    _STRING_SIZE = 3
    # Key to append to resouce keys
    RESOURCE_PREFIX = b"_resources"

    def __init__(self):
        super().__init__()
        self.data = {}

    def add_shared_resource(self, item):
        """Generates a unique id to refer to a shared resource, form is sequential AAB, AAC, AAD, etc.

        Args:
            item: an item to add to the dictionary

        Returns:
            A string of length _STRING_SIZE that is the key for the added value.
        """
        cur_items = len(self) + 1
        ascii_lowend = ord("A")
        ascii_highend = ord("Z") + 1
        val_range = ascii_highend - ascii_lowend
        max_val = val_range**ResourceCollection._STRING_SIZE - 1
        if cur_items > max_val:
            raise ValueError("Limit on total number of shared resources is exceeded.")
        vals = [
            ascii_lowend + (cur_items // (val_range**i)) % val_range
            for i in range(ResourceCollection._STRING_SIZE - 1, -1, -1)
        ]
        key = bytes(bytearray(vals))
        self.data[key.decode()] = item
        return (ResourceCollection.RESOURCE_PREFIX + b"_" + key).decode()

    def __len__(self):
        """Used to determine if any data is present."""
        return len(self.data)

    def to_dict_for_json(self):
        return self.data


@dataclass
class SampleData(React10X):
    """Holds values used to populate the title and header of a websummary."""

    id: str  # pylint: disable=invalid-name
    description: str
    pipeline: str


class CompiledWebSummaryData(React10X):
    """A class that enables optimized storage in the web summary.

    Before serializing
    we look for shared data in the object by object `id` for complex objects or by
    equality for a set of some primitive arrays (e.g. lists of floats).
    """

    _RESOURCES_STRING = "resources"

    def __init__(self):
        super().__init__()
        # make pylint happy
        self.resources = None
        setattr(self, CompiledWebSummaryData._RESOURCES_STRING, ResourceCollection())
        # make sure making pylint happy didn't break anything
        assert self.resources is not None

    def compile(self):
        """Recursively examine the object and look for shared arrays it can combine."""
        # my_new_self = copy.deepcopy(self)
        # TODO: Implement

    def to_dict_for_json(self):
        """Convert to a dictionary suitable for json serialization.

        Same as base class, but we convert the python publicaly accessible `resources` to the
        javascript specified `_resources` string that uses an underscore to avoid namespace conflicts
        """
        data = vars(self)
        if CompiledWebSummaryData._RESOURCES_STRING in data:
            values = data[CompiledWebSummaryData._RESOURCES_STRING]
            del data[CompiledWebSummaryData._RESOURCES_STRING]
        else:
            values = []

        if len(values) > 0:
            data[ResourceCollection.RESOURCE_PREFIX.decode()] = values
        return data


class BiplotsWithThresholding(React10X):  # pylint: disable=too-few-public-methods
    """A wrapper for the corresponding React Component."""

    def __init__(
        self,
        vectors,
        vector_names,
        show_transparency_slider,
        show_thresholding_slider,
        threshold_values,
        threshold_slider_values,
        threshold_slider_title,
        unfiltered_name,
        filtered_name,
        percent_filtered_name,
        excluded_categories,
        assignments,
        assignment_colors,
    ):  # pylint: disable=too-many-arguments
        """Construct an instance.

        Args:
            vectors: [list of lists of float values],
            vector_names:[list of names for each vector in vectors]
            show_transparency_slider: boolean to show slider or not
            threshold_values: values to threshold points on
            threshold_slider_values: [min, max, step] for slider
            threshold_slider_title: string to display over threshold slider
            unfiltered_name: string to label current unfiltered count with
            filtered_name: string to label current filtered count with
            percent_filtered_name: string to label percentage filtered
            excluded_categories: list of assignments not to include in filtering counts or
                           once filter bar is moved to the right at all
            assignments: [list of values to label the points with],
            assignment_colors: {dictionary of assignment_label:colors}
        """
        super().__init__()
        self.vectors = vectors
        self.vector_names = vector_names
        self.show_transparency_slider = show_transparency_slider
        self.show_thresholding_slider = show_thresholding_slider
        self.threshold_values = threshold_values
        self.threshold_slider_values = threshold_slider_values
        self.threshold_slider_title = threshold_slider_title
        self.unfiltered_name = unfiltered_name
        self.filtered_name = filtered_name
        self.percent_filtered_name = percent_filtered_name
        self.excluded_categories = excluded_categories
        self.assignments = assignments
        self.assignment_colors = assignment_colors


class ClusteringData(React10X):  # pylint: disable=too-few-public-methods
    """POD class for data on a clustering."""

    def __init__(self, key, clustering, data):
        super().__init__()
        self.key = key
        self.name = clustering.description
        self.data = data


class HeaderWithHelp(React10X):  # pylint: disable=too-few-public-methods
    """POD class for equivalent React component."""

    def __init__(self, title, help_text=None):
        """Create a new HeaderWithHelp Component.

        Args:
            title (str): The title text
            help_text (str): Optional string of plain text or html that will
                appear when the question mark is clicked.
        """
        super().__init__()
        self.title = title
        if help_text is not None:
            # pylint: disable=invalid-name
            self.helpText = help_text


class PlotWithHeader(React10X):  # pylint: disable=too-few-public-methods,invalid-name
    """POD class for equivalent React component."""

    def __init__(self, title, plot, help_text=None):
        """Make a plot with a websummary style component.

        Args:
            title: Title to display
            plot: plotly plot data
            help_text: Optional help text
        """
        super().__init__()
        self.help = HeaderWithHelp(title, help_text)
        self.plot = plot


class Clusterings(React10X):  # pylint: disable=too-few-public-methods,invalid-name
    """Stores similar data for a collection of different possible clusterings.

    Examples of such include differential expression tables or plots.
    """

    def __init__(self):
        super().__init__()
        self.clusterings = []

    def add_clustering(self, data):
        assert isinstance(data, ClusteringData)
        self.clusterings.append(data)


class SharedCoordinatePlotCollection(Clusterings):  # pylint: disable=too-few-public-methods
    """Class used to store a collection of plots that share layout/config coordinates.

    Used for serialization to the WebSummary
    where this type is interpreted and expanded in the JavaScript code.
    """

    def __init__(self, x, y, config, layout, plt_type, marker):
        """Construct an instance.

        Args:
            x: Vector of doubles
            y: Vector of doubles
            config: dict with plotly config elements
            layout: dict with plotly layout elements
            plt_type: either `scatter` or `scattergl`
            marker: The marker to use in the scatter plot JSON, e.g. {"opacity": 0.9, "size": 4}
        """
        super().__init__()
        self.x = round_floats_in_list(x)
        self.y = round_floats_in_list(y)
        self.config = config
        self.layout = layout
        self.type = plt_type
        self.marker = marker

    def add_clustering(self, key, clustering):  # pylint: disable=arguments-differ
        new_clustering = ClusteringData(key, clustering, self._clustering_to_clusters(clustering))
        super().add_clustering(new_clustering)

    def _clustering_to_clusters(self, clustering):
        clustering_labels = clustering.clusters
        num_cells = len(clustering_labels)
        data = []
        for i in range(max(clustering_labels)):
            index = i + 1
            name = f"Cluster {index}"
            indices = np.where(clustering_labels == index)[0]
            prop = len(indices) * 1.0 / num_cells
            new_cluster = {
                "name": name,
                "indices": list(indices),
                "type": self.type,
                "mode": "markers",
                "marker": self.marker,
                "text": f"{name}: {prop:.1%}",
            }
            data.append(new_cluster)
        return data

    def to_dict_for_json(self):
        """We drop the unusesd variables."""
        del self.type
        del self.marker
        return vars(self)


class ClusteringSelector(React10X):
    """Mimics the data requirements and model of the React Component in ClusteringSelector.js."""

    def __init__(self, plot_help_txt, table_help_txt):
        """Initialize the data for a ClusteringSelector Component, which displays a top row of plots.

        one each on the left and right side, and a bottom row with a table.  Both the plots and the
        table have a header above them.  The plot on the right side can either be a single plot
        that does not change as different clusterings are selected (as in Cell Ranger) or a list of
        plots equal in size to the number of clusters used (as in Space Ranger).

        Args:
            plot_help_txt: dict of data needed for the React component `DynamicHelptext` that
                will be displayed above the plots
            table_help_txt: dict of data needed for the React Component `HeaderWithHelp` that
                will be displayed above the table.
        """
        super().__init__()
        self._left_plots = None
        self._right_plots = None
        self._tables = None
        self.table_help_txt = table_help_txt
        self.plot_help_txt = plot_help_txt

    def _validate_same_clusters(self):
        """Verifies that array clusters match.

        The plots/tables should all be arrays with the same
        cluster represented at each position in them.
        """
        to_check = [
            x
            for x in [self._left_plots, self._right_plots, self._tables]
            if x and isinstance(x, Clusterings)
        ]
        if len(to_check) > 1:
            baseline = to_check[0].clusterings
            for alt in to_check[1:]:
                if any(
                    x.key != y.key or x.name != y.name for x, y in zip(baseline, alt.clusterings)
                ):
                    raise ValueError(
                        "The clusterings for the plots and tables of a "
                        "clustering selector must be the same and in order."
                    )

    @property
    def right_plots(self):
        return self._right_plots

    @right_plots.setter
    def right_plots(self, value):
        """The right plots must be a SharedCoordinatePlotCollection."""
        assert isinstance(value, SharedCoordinatePlotCollection)
        self._right_plots = value
        self._validate_same_clusters()

    @property
    def left_plots(self):
        return self._left_plots

    @left_plots.setter
    def left_plots(self, value):
        """Generate left plots.

        Args:
            value: Either a single plot, or a list of plots.
        """
        assert isinstance(value, dict | SharedCoordinatePlotCollection)
        self._left_plots = value
        self._validate_same_clusters()

    @property
    def tables(self):
        return self._tables

    @tables.setter
    def tables(self, value):
        assert isinstance(value, Clusterings)
        self._tables = value
        self._validate_same_clusters()

    def to_dict_for_json(self):
        l_plots = self._left_plots
        if isinstance(l_plots, SharedCoordinatePlotCollection):
            l_plots = l_plots.to_dict_for_json()
        return {
            "left_plots": l_plots,
            "right_plots": self._right_plots,
            "plot_help": self.plot_help_txt,
            "table_help": self.table_help_txt,
            "tables": self._tables,
        }


class DifferentialExpressionTableValue(React10X):
    # pylint: disable=too-few-public-methods
    """Holds the data used to render a cell in the DifExp table."""

    def __init__(self, pval, log2_fold_change, greyed_out):
        super().__init__()
        assert isinstance(pval, float)
        assert isinstance(log2_fold_change, float)
        assert isinstance(greyed_out, bool)
        self.pval = _round_float(pval)
        self.log2_fold_change = _round_float(log2_fold_change)
        self.greyed_out = greyed_out

    def to_dict_for_json(self):
        oval = {"p": self.pval, "l": self.log2_fold_change}
        # Only output this if True
        # And set it as 1 to use less space than the boolean in JSON
        # And be interpreted correctly by the web summary code (copied below)
        # const gcval = get(row.original, gc);
        # const col = gcval ? '#DDD' : default_col;
        if self.greyed_out:
            oval["g"] = 1
        return oval


class BarnyardPanel(React10X):
    # pylint: disable=too-few-public-methods
    """Data for a barnyard plot on analysis tab."""

    def __init__(self, plot, gems):
        super().__init__()
        assert plot is not None
        assert gems is not None
        self.plot = plot
        self.gems = gems

    def to_dict_for_json(self):
        return {"plot": self.plot, "gems": self.gems}


def _change_all_plots(d, plt_func):
    """Recursive function to update all elements in a web summary.

    Args:
        d: the data to iterate/recurse through
        plt_func: A function that takes a dictionary or object representing plot data and
            modifies it
    """
    # pylint: disable=invalid-name
    if isinstance(d, dict) and all(key in d for key in ["data", "layout", "config"]):
        plt_func(d)
    elif isinstance(d, dict):
        for v in d.values():
            _change_all_plots(v, plt_func)
    elif isinstance(d, list):
        for v in d:
            _change_all_plots(v, plt_func)
    elif isinstance(d, SharedCoordinatePlotCollection):
        plt_func(d)
    elif isinstance(d, ClusteringSelector):
        _change_all_plots(d.left_plots, plt_func)
        _change_all_plots(d.right_plots, plt_func)


class WebSummaryData(React10X):
    """Class to store the web summary data as an object instead of as pure JSON.

    This class is a builder and will be invalid after any call to to_dict_for_json
    as this method may modify or change the data to fit conventions of the
    web summary
    """

    _SUMMARY_TAB = "summary_tab"
    _ANALYSIS_TAB = "analysis_tab"
    _ANTIBODY_TAB = "antibody_tab"
    _ANTIGEN_TAB = "antigen_tab"
    _CELL_TYPES_TAB = "cell_types_tab"
    # Only used by space ranger
    _ALIGNMENT_TAB = "alignment_tab"
    _ALARMS = "alarms"
    _SPATIAL_TAB = "spatial_tab"
    _CLUSTERING_SELECTOR = "clustering_selector"
    _ANTIBODY_CLUSTERING_SELECTOR = "antibody_clustering_selector"
    _VDJ_TAB = "vdj_analysis_tab"
    _DIAGNOSTICS = "diagnostics"

    def __init__(self, sample_properties, command, pipeline):
        super().__init__()
        self._data = {
            WebSummaryData._ALARMS: {WebSummaryData._ALARMS: []},
            "sample": {
                "id": sample_properties.sample_id,
                "description": sample_properties.sample_desc,
                "command": command,
                "subcommand": pipeline,
            },
            WebSummaryData._SUMMARY_TAB: {},
            WebSummaryData._ANALYSIS_TAB: {},
            WebSummaryData._ANTIBODY_TAB: {},
            WebSummaryData._CELL_TYPES_TAB: {},
            WebSummaryData._ANTIGEN_TAB: {},
            WebSummaryData._ALIGNMENT_TAB: {},
            WebSummaryData._SPATIAL_TAB: {},
            WebSummaryData._VDJ_TAB: {},
            WebSummaryData._DIAGNOSTICS: {},
            ResourceCollection.RESOURCE_PREFIX: ResourceCollection(),
        }
        self._clustering_selector = None
        self._antibody_clustering_selector = None

    def _clear_if_empty(self, value):
        if not self._data[value]:
            del self._data[value]

    def _set_all_plots_font(self, font_name=pltly.DEFAULT_WEB_FONT):
        def format_fn(x):
            return pltly.add_font_to_layout(x, font=font_name)

        _change_all_plots(self._data, format_fn)

    def _share_plot_images(self) -> None:
        """Scans through all the plots, and de-duplicates shared images.

        Modifications are in place.
        """
        images_str = "images"
        source_str = "source"

        def try_get_source(plt):
            layout = pltly.get_layout_from_chart(plt)
            image = layout.get(images_str)
            if image is not None and len(image) == 1:
                # Assume we only ever have one image
                source = image[0].get(source_str)
                if source is not None:
                    return source
            return None

        images = {}

        def counting_closure(plt):
            source = try_get_source(plt)
            if source is not None:
                if source in images:
                    images[source] += 1
                else:
                    images[source] = 1

        _change_all_plots(self._data, counting_closure)

        replacements = {}
        for key, value in images.items():
            if value > 1:
                resource_id = self.shared_resources.add_shared_resource(key)
                replacements[key] = resource_id

        if len(replacements) > 0:

            def replacement_closure(plt):
                source = try_get_source(plt)
                if source is not None:
                    layout = pltly.get_layout_from_chart(plt)
                    layout[images_str][0][source_str] = replacements[source]

            _change_all_plots(self._data, replacement_closure)

    # pylint: disable=invalid-name
    def to_dict_for_json(self):
        """Should only be called once at the end."""
        self._check_valid()
        self._format_clustering_selector()
        self._format_antibody_clustering_selector()
        self._set_all_plots_font()
        self._share_plot_images()

        to_filter = [
            WebSummaryData._ALIGNMENT_TAB,
            WebSummaryData._ANALYSIS_TAB,
            WebSummaryData._ANTIBODY_TAB,
            WebSummaryData._CELL_TYPES_TAB,
            WebSummaryData._ANTIGEN_TAB,
            WebSummaryData._SPATIAL_TAB,
            WebSummaryData._VDJ_TAB,
            WebSummaryData._DIAGNOSTICS,
            ResourceCollection.RESOURCE_PREFIX,
        ]
        for k in to_filter:
            self._clear_if_empty(k)

        to_return = {"summary": self._data}
        if ResourceCollection.RESOURCE_PREFIX in self._data:
            to_return[ResourceCollection.RESOURCE_PREFIX.decode()] = self._data[
                ResourceCollection.RESOURCE_PREFIX
            ]
            del self._data[ResourceCollection.RESOURCE_PREFIX]
        # Invalidate future calls by deleting data leading to errors
        del self._data
        return to_return

    def _check_valid(self):
        """Make sure we are not calling class after formatting and altering it."""
        if not hasattr(self, "_data"):
            raise KeyError(
                "The web summary was accessed after the data was returned and "
                "the instance invalidated."
            )

    def _format_clustering_selector(self):
        """If a clustering selector is present, prepare it for serialization to JSON.

        Also, invalidate the reference.
        """
        if self._clustering_selector:
            self.analysis_tab[WebSummaryData._CLUSTERING_SELECTOR] = self._clustering_selector
            self._clustering_selector = None

    def _format_antibody_clustering_selector(self):
        """If a clustering selector is present, prepare it for serialization to JSON.

        Also, invalidate the reference.
        """
        if self._antibody_clustering_selector:
            self.antibody_tab[WebSummaryData._ANTIBODY_CLUSTERING_SELECTOR] = (
                self._antibody_clustering_selector
            )
            self._antibody_clustering_selector = None

    def _get_safe(self, key):
        self._check_valid()
        return self._data.get(key)

    def _set_tab(self, key, value):
        self._check_valid()
        self._data[key] = value

    @property
    def summary_tab(self):
        return self._get_safe(WebSummaryData._SUMMARY_TAB)

    @property
    def alarms(self):
        self._check_valid()
        return self._data[WebSummaryData._ALARMS][WebSummaryData._ALARMS]

    @property
    def analysis_tab(self):
        return self._get_safe(WebSummaryData._ANALYSIS_TAB)

    @property
    def antibody_tab(self):
        return self._get_safe(WebSummaryData._ANTIBODY_TAB)

    @property
    def antigen_tab(self):
        return self._get_safe(WebSummaryData._ANTIGEN_TAB)

    @property
    def cell_types_tab(self):
        return self._get_safe(WebSummaryData._CELL_TYPES_TAB)

    @property
    def diagnostics(self):
        return self._get_safe(WebSummaryData._DIAGNOSTICS)

    @diagnostics.setter
    def diagnostics(self, value):
        if value is not None:
            self._set_tab(WebSummaryData._DIAGNOSTICS, value)

    @property
    def shared_resources(self):
        return self._get_safe(ResourceCollection.RESOURCE_PREFIX)

    @property
    def clustering_selector(self):
        self._check_valid()
        return self._clustering_selector

    @clustering_selector.setter
    def clustering_selector(self, value):
        if value is None:  # We can continually set this to None
            self._clustering_selector = value
        else:
            self._check_valid()
            assert isinstance(value, ClusteringSelector)
            self._clustering_selector = value

    @property
    def antibody_clustering_selector(self):
        self._check_valid()
        return self._antibody_clustering_selector

    @antibody_clustering_selector.setter
    def antibody_clustering_selector(self, value):
        if value is None:  # We can continually set this to None
            self._antibody_clustering_selector = value
        else:
            self._check_valid()
            assert isinstance(value, ClusteringSelector)
            self._antibody_clustering_selector = value

    @property
    def alignment_tab(self):
        return self._get_safe(WebSummaryData._ALIGNMENT_TAB)

    @alignment_tab.setter
    def alignment_tab(self, value):
        self._set_tab(WebSummaryData._ALIGNMENT_TAB, value)

    @property
    def vdj_tab(self):
        return self._get_safe(WebSummaryData._VDJ_TAB)

    @vdj_tab.setter
    def vdj_tab(self, value):
        self._set_tab(WebSummaryData._VDJ_TAB, value)

    @property
    def spatial_tab(self):
        self._check_valid()
        return self._data[WebSummaryData._SPATIAL_TAB]
