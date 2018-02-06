#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
import collections
import copy
import csv
import itertools
import json
import math
import random
import re
import numpy as np
import sys
import tenkit.safe_json as tk_safe_json
import cellranger.analysis.clustering as cr_clustering
import cellranger.analysis.singlegenome as cr_sg_analysis
import cellranger.analysis.multigenome as cr_mg_analysis
import cellranger.constants as cr_constants
import cellranger.reference as cr_reference
import cellranger.report as cr_report
import cellranger.utils as cr_utils
import cellranger.vdj.constants as vdj_constants
import cellranger.vdj.report as vdj_report
import cellranger.webshim.constants.shared as shared_constants
import cellranger.webshim.constants.gex as ws_gex_constants
import cellranger.webshim.constants.vdj as ws_vdj_constants
from cellranger.webshim.data import SampleData
import cellranger.webshim.template as template

def max_norm(values):
    min_value = min(values)
    values = [v - min_value for v in values]

    max_value = max(values)
    if max_value > 0:
        return [float(v) / float(max_value) for v in values]
    return values

def add_prefix(prefix, name):
    if '%s' in name:
        return name % prefix
    else:
        return '_'.join([prefix, name]) if prefix else name

def format_name(display_name, prefix, prefixes, prefix_format_func=None):
    if prefix is not None and prefix_format_func is not None:
        prefix = prefix_format_func(prefix)

    # Default multi -> '' if no format func given
    if prefix_format_func is None and prefix == cr_constants.MULTI_REFS_PREFIX:
        prefix = ''

    if len(prefixes) > 1 or '%s' in display_name:
        display_name = add_prefix(prefix, display_name)

    # Replace underscores w/ spaces
    display_name = display_name.replace('_', ' ')

    # Collapse whitespace
    display_name = re.sub('\s+', ' ', display_name)

    return display_name

def format_description(description, prefix, prefixes, prefix_format_func=None):
    if prefix is not None and prefix_format_func is not None:
        prefix = prefix_format_func(prefix)

    # Default multi -> '' if no format func given
    if prefix_format_func is None and prefix == cr_constants.MULTI_REFS_PREFIX:
        prefix = ''

    if '%s' in description:
        # Escape stray percents
        description = re.sub('%([^s])', '%%\\1', description)

        # Only add the prefix if there are multiple possibilities
        s = str(prefix) if len(prefixes) > 1 and prefix else ''
        description = description % s

    # Collapse whitespace
    description = re.sub('\s+', ' ', description)

    return description

def format_value(value, format_type):
    if format_type == 'string':
        return str(value)

    value = value if not math.isnan(float(value)) else 0

    if format_type == 'percent':
        return '{:.1%}'.format(float(value))
    elif format_type == 'integer':
        return '{:,d}'.format(int(value))
    elif format_type[0] == '%':
        return format_type % float(value)
    raise Exception('Invalid format type: %s' % format_type)

def lookup_name(data, name):
    name_parts = name.split('/', 1)
    data = data.get(name_parts[0])
    if len(name_parts) == 1 or data is None:
        return data
    else:
        return lookup_name(data, name_parts[1])

def add_alarm(value, formatted_value, level, alarm_dict, alarms):
    # The numeric value may have become a string via JSON ('NaN', 'Infinity', '-Infinity')
    # Note that by default python's json encoder outputs a bare NaN symbol which is even worse
    test = 'float("%s") %s' % (str(value), alarm_dict['test'])
    raised = eval(test, {}, {})
    if raised:
        title = alarm_dict['title']
        message = alarm_dict['message']
        alarms[title] = {
            'title': title,
            'message': message,
            'value': formatted_value,
            'level': level,
        }
    return raised

def add_alarms(data, alarms_dict, alarms):
    full_name = alarms_dict['name']    # Prefixes on alarm metrics not implemented
    value = lookup_name(data, full_name)
    if value is None:
        return

    formatted_value = format_value(value, alarms_dict['format'])
    error_dict = alarms_dict.get(shared_constants.ALARM_ERROR)
    if error_dict is not None and add_alarm(value, formatted_value, shared_constants.ALARM_ERROR, error_dict, alarms):
        return
    warn_dict = alarms_dict.get(shared_constants.ALARM_WARN)
    if warn_dict is not None:
        add_alarm(value, formatted_value, shared_constants.ALARM_WARN, warn_dict, alarms)

def add_table_rows(data, name, metric_dict, rows, style_func, target_func, prefixes=None):
    """ Args:
          data (dict) - Dict containing summary metrics
          name (str) - Metric name (key in data dict)
          metric_dict (dict) - Metric display definition dict
    """
    format_type = metric_dict['format']
    display_name = metric_dict.get('display_name', name)
    description = metric_dict.get('description', display_name)

    if not prefixes:
        prefixes = [None]

    values = []
    for prefix in prefixes:
        # Find the metric and construct a table row
        full_name = add_prefix(prefix, name)
        value = lookup_name(data, full_name)

        if value is not None:
            prefix_format_func = metric_dict.get('prefix_format')
            formatted_name = format_name(display_name, prefix, prefixes,
                                         prefix_format_func=prefix_format_func)
            formatted_description = format_description(description, prefix, prefixes,
                                                       prefix_format_func=prefix_format_func)
            formatted_value = format_value(value, format_type)

            style = style_func(metric_dict, value)

            rows.append([
                {
                    'v': formatted_name,
                    'f': formatted_description,
                    's': style,
                },
                {
                    'v': target_func(metric_dict),
                    's': style,
                },
                {
                    'v': formatted_value,
                    's': style,
                    'r': value,
                },
            ])
            values.append(value)

def is_metric_hidden(sample_properties, metric_dict):
    if 'hidden' in metric_dict:
        return eval(metric_dict['hidden'], sample_properties, sample_properties)
    return False

def build_tables(sample_properties, table_dicts, alarm_table_dicts, sample_data,
                 style_func=lambda *args: "",
                 target_func=lambda *args: "",
                 metric_headers=[],
                 all_prefixes={}):

    tables = []
    for table_dict in table_dicts:
        table_name = table_dict['name']
        table_metrics = table_dict['metrics']

        rows = []
        for metric_dict in table_metrics:
            name = metric_dict['name']
            prefix = metric_dict.get('prefix')

            if is_metric_hidden(sample_properties, metric_dict):
                continue

            if prefix is not None:
                # Search sample properties for prefix values
                prefixes = sample_properties.get(prefix)

                # Otherwise search the given prefix values
                if not prefixes:
                    prefixes = all_prefixes.get(prefix)

                # Apply the prefix filter for this metric
                if 'prefix_filter' in metric_dict:
                    if prefixes is not None:
                        prefixes = filter(lambda x: x in metric_dict['prefix_filter'], prefixes)
                    else:
                        sys.stderr.write('Warning: no metric prefix values to filter for metric %s. The prefix name %s is probably wrong.\n' % (name, prefix))

                add_table_rows(sample_data.summary, name, metric_dict, rows, style_func, target_func, prefixes=prefixes)
            else:
                add_table_rows(sample_data.summary, name, metric_dict, rows, style_func, target_func)

        if len(rows) > 0:
            tables.append({
                'name': table_name,
                'rows': rows,
                'headers': metric_headers,
            })

    # Build alarm table
    alarms = collections.OrderedDict()
    for alarm_dict in alarm_table_dicts:
        add_alarms(sample_data.summary, alarm_dict, alarms)

    return tables, alarms.values()

def convert_numpy_array_to_line_chart(array, ntype):
    array = np.sort(array)[::-1]

    rows = []
    previous_count = None
    for (index,), count in np.ndenumerate(array):
        if index == 0 or index == len(array)-1:
            rows.append([index, ntype(count)])
        elif previous_count != count:
            previous_index = rows[-1][0]
            if previous_index != index - 1:
                rows.append([index - 1, ntype(previous_count)])
            rows.append([index, ntype(count)])
        previous_count = count
    return rows

def _plot_barcode_rank(chart, counts, num_cells):
    """ Generate a generic barcode rank plot """
    rows = convert_numpy_array_to_line_chart(counts, int)

    for i, row in enumerate(rows):
        index, count = row[0], row[1]
        if index < num_cells:
            series_list = [chart['data'][0]]
        elif index == num_cells:
            # Connect the two lines
            series_list = [chart['data'][0], chart['data'][1]]
        else:
            series_list = [chart['data'][1]]

        for series in series_list:
            series['x'].append(index)
            series['y'].append(count)

    # Handle case where there is no background
    bg_series = chart['data'][1]
    if len(bg_series['x']) == 0:
        bg_series['x'].append(0)
        bg_series['y'].append(0)

    return chart

def plot_barcode_rank(chart, sample_properties, sample_data):
    """ Generate the RNA counter barcode rank plot """
    if sample_properties.get('genomes') is None or sample_data.barcode_summary is None:
        return None

    if len(sample_properties['genomes']) == 0:
        return None

    # UMI counts per BC across all genomes present
    if len(sample_properties['genomes']) > 1:
        genome = cr_constants.MULTI_REFS_PREFIX
    else:
        genome = sample_properties['genomes'][0]

    key = cr_utils.format_barcode_summary_h5_key(genome,
                                                 cr_constants.TRANSCRIPTOME_REGION,
                                                 cr_constants.CONF_MAPPED_DEDUPED_READ_TYPE)

    if key in sample_data.barcode_summary:
        counts_per_bc = sample_data.barcode_summary[key][:]
        return _plot_barcode_rank(chart, counts_per_bc, sample_data.num_cells)
    else:
        # Not guaranteed to exist, depending on pipeline
        pass


def plot_vdj_barcode_rank(chart, sample_properties, sample_data):
    """ Generate the VDJ barcode rank plot """
    bc_support = sample_data.vdj_barcode_support
    if bc_support is None:
        return

    summary_data = sample_data.summary
    num_cells = summary_data.get('vdj_filtered_bcs')
    if num_cells is None:
        return

    support = np.sort(bc_support['count'].as_matrix())[::-1]

    return _plot_barcode_rank(chart, support, num_cells)


def plot_clonotype_table(chart, sample_properties, sample_data):
    if sample_data.vdj_clonotype_summary is None:
        return None

    clonotypes = sample_data.vdj_clonotype_summary.iloc[0:10]

    # This column used to be called 'cdr3s'; allow the webshim to work on older data
    cdr3_aa_col = 'cdr3s_aa'
    if cdr3_aa_col not in clonotypes:
        cdr3_aa_col = 'cdr3s'

    col_defs = collections.OrderedDict([
        ('clonotype_id', {'label': 'Clonotype ID',
                          'format': 'string',
                          'title': 'Clonotype ID',
                          'style': 'text-align: left'}),
        (cdr3_aa_col,     {'label': 'CDR3s',
                           'format': 'string',
                           'title': 'CDR3s in clonotype',
                           'style': 'text-align: left'}),
        ('frequency',    {'label': 'Frequency',
                          'format': 'integer',
                          'title': 'Number of cells with clonotype',
                          'style': 'text-align: right'}),
        ('proportion',   {'label': 'Proportion',
                          'format': '%0.4f',
                          'title': 'Fraction of cell with clonotype',
                          'style': 'text-align: right'}),
    ])

    cols = []
    for name, col_def in col_defs.iteritems():
        if name not in clonotypes:
            raise ValueError('Column not found in clonotype summary: %s' % name)
        cols.append({
            'label': col_defs[name]['label'],
            'title': col_defs[name]['title'],
        })

    rows = []
    for _, cl_row in clonotypes.iterrows():
        row = []
        for col_name, col_def in col_defs.iteritems():
            value = cl_row[col_name]
            formatted_value = format_value(value, col_def['format'])

            # Make the CDR3 list bit more readable
            formatted_value = formatted_value.replace(';', '; ')

            row.append({
                'v': tk_safe_json.json_sanitize(value),
                'f': formatted_value,
                's': col_def['style'],
            })
        rows.append(row)

    chart['table'].update({'rows': rows, 'cols': cols})

    return chart


def plot_histogram_metric(chart, sample_properties, sample_data, **kwargs):
    """ Plot a HistogramMetric from the summary json """
    summary_data = sample_data.summary
    items = summary_data.get(kwargs['metric_name'], {}).items()
    if len(items) < 1:
        return None

    ordering = kwargs.get('order_by', shared_constants.HISTOGRAM_METRIC_DEFAULT_ORDERING)

    if ordering == shared_constants.HISTOGRAM_METRIC_ORDER_INTEGER_BIN:
        items.sort(key=lambda x: convert_to_int_gracefully(x[0]))

    elif ordering == shared_constants.HISTOGRAM_METRIC_ORDER_DECREASING_FREQUENCY:
        items.sort(key=lambda x: -convert_to_int_gracefully(x[1]))

    elif ordering == shared_constants.HISTOGRAM_METRIC_ORDER_DECREASING_PROPORTION:
        items.sort(key=lambda x: -convert_to_float_gracefully(x[1]))

    x, y = zip(*items)
    chart['data'][0].update({'x': x, 'y': y})

    return chart


def plot_barnyard_barcode_counts(chart, sample_properties, sample_data):
    analysis = sample_data.analysis
    if not isinstance(analysis, cr_mg_analysis.MultiGenomeAnalysis):
        return None

    chart['data'] = []
    for label_info in ws_gex_constants.GEM_CALL_LABELS:
        name = label_info['label']
        name = name.replace('genome0', analysis.result['genome0'])
        name = name.replace('genome1', analysis.result['genome1'])
        chart['data'].append({
            'x': [],
            'y': [],
            'name': name,
            'mode': 'markers',
            'marker': {
                'color': label_info['color'],
            },
        })

    call_to_series = {v['key']: i for i, v in enumerate(ws_gex_constants.GEM_CALL_LABELS)}
    for count0, count1, call in itertools.izip(analysis.result['count0'],
                                               analysis.result['count1'],
                                               analysis.result['call']):
        series = chart['data'][call_to_series[call]]
        series['x'].append(int(count0))
        series['y'].append(int(count1))

    chart['layout']['xaxis'] = {
        'title': '%s UMI counts' % analysis.result['genome0'],
        'rangemode': 'tozero',
        'autorange': True,
    }
    chart['layout']['yaxis'] = {
        'title': '%s UMI counts' % analysis.result['genome1'],
        'rangemode': 'tozero',
        'autorange': True,
    }

    return chart

def plot_preprocess(analysis, sample_properties):
    if analysis is None:
        return None
    if isinstance(analysis, cr_sg_analysis.SingleGenomeAnalysis) and analysis.is_zero_matrix():
        return None
    if isinstance(analysis, cr_mg_analysis.MultiGenomeAnalysis):
        return analysis

    # Limit the number of K-means clusterings displayed to limit the HTML filesize
    new_clusterings = {}
    for key, clu in analysis.clusterings.iteritems():
        if not (clu.clustering_type == cr_clustering.CLUSTER_TYPE_KMEANS and \
                clu.num_clusters > ws_gex_constants.MAX_WEBSHIM_KMEANS_K):
            new_clusterings[key] = clu
    analysis.clusterings = new_clusterings

    return analysis

def load_sample_data(sample_properties, sample_data_paths):
    return SampleData(sample_properties, sample_data_paths, plot_preprocess)

def plot_tsne(chart, sample_properties, sample_data):
    """ Plot cells in t-SNE space, colored by clustering label """
    analysis = sample_data.analysis
    if not analysis or len(sample_properties['genomes']) > 1 or not isinstance(analysis, cr_sg_analysis.SingleGenomeAnalysis):
        return None

    args = [analysis.get_tsne().transformed_tsne_matrix,
            ws_gex_constants.TSNE_CLUSTER_DESCRIPTION,
            1, 2]
    return clustering_plot_func(chart, sample_properties, sample_data, plot_dimensions, args)

def plot_dimensions(chart, transformed_matrix, description, pc1, pc2,
                    clip=None, clustering=None, diff_expr=None, values=None,
                    original_cluster_sizes=None):
    """ Plot cells in a 2-d space, colored by clustering label """
    assert (values is not None and clustering is None) or \
        (clustering is not None and values is None and original_cluster_sizes is not None)

    if clustering is not None:
        values = clustering.clusters

    if original_cluster_sizes is None:
        value_freqs = np.bincount(values)[1:].tolist()
    else:
        # Cluster size array starts w/ cluster1 at index0; make it equivalent to np.bincount on arbitrary values
        value_freqs = [0] + original_cluster_sizes.tolist()

    n, m = transformed_matrix.shape
    if m < max(pc1, pc2):
        return None

    max_value = values.max()

    chart['data'] = [None] * max_value
    for value in np.unique(values):
        chart['data'][value - 1] = {
            'x': [],
            'y': [],
            'name': '%s %s - %s cells' % (
                description, format_value(value, 'integer'), format_value(value_freqs[value], 'integer'),
            ),
            'hoverinfo': 'name',
            'mode': 'markers',
            'type': 'scattergl',
            'marker': {
                'size': 4,
            },
        }

    for i in xrange(n):
        r1 = shared_constants.DATA_VALUE_FORMAT % transformed_matrix[i, pc1-1]
        r2 = shared_constants.DATA_VALUE_FORMAT % transformed_matrix[i, pc2-1]
        value = values[i]

        series = chart['data'][value - 1]
        series['x'].append(r1)
        series['y'].append(r2)

    if clip is not None:
        xmin, xmax = np.percentile(transformed_matrix[:, pc1-1], clip)
        ymin, ymax = np.percentile(transformed_matrix[:, pc2-1], clip)
        chart['layout']['xaxis'] = {
            'range': [xmin, xmax],
        }
        chart['layout']['yaxis'] = {
            'range': [ymin, ymax],
        }

    chart['config'] = shared_constants.CHARTS_PLOTLY_MOVABLE_CONFIG
    return chart

def plot_dimensions_color(chart, transformed_matrix, values, description, vmin, vmax, pc1, pc2):
    _, m = transformed_matrix.shape
    if m < max(pc1, pc2):
        return None

    series = chart['data'][0]

    index_order = range(transformed_matrix.shape[0])
    random.shuffle(index_order)
    for i in index_order:
        r1 = shared_constants.DATA_VALUE_FORMAT % transformed_matrix[i, pc1-1]
        r2 = shared_constants.DATA_VALUE_FORMAT % transformed_matrix[i, pc2-1]
        value = int(values[i])
        text = '%s: %s' % (description, format_value(value, 'integer'))

        if value < vmin:
            value = vmin
        if value > vmax:
            value = vmax

        series['x'].append(r1)
        series['y'].append(r2)
        series['marker']['color'].append(value)
        series['text'].append(text)

    return chart

def plot_tsne_totalcounts(chart, sample_properties, sample_data):
    """ Plot cells colored by total counts """
    analysis = sample_data.analysis
    if not analysis or len(sample_properties['genomes']) > 1:
        return None

    reads_per_bc = analysis.matrix.get_reads_per_bc()
    vmin, vmax = np.percentile(reads_per_bc, ws_gex_constants.TSNE_TOTALCOUNTS_PRCT_CLIP)

    return plot_dimensions_color(chart, analysis.get_tsne().transformed_tsne_matrix,
                                 reads_per_bc,
                                 ws_gex_constants.TSNE_TOTALCOUNTS_DESCRIPTION,
                                 vmin, vmax,
                                 1, 2)

def _plot_differential_expression(chart, analysis, clustering=None, diff_expr=None, original_cluster_sizes=None):
    n_clusters = clustering.clusters.max()

    # Get the union of top DE genes
    top_genes = set()

    # Limit the number of entries in the DE table
    n_genes = int(np.floor(float(ws_gex_constants.MAX_DE_TABLE_ENTRIES) / (n_clusters**2)))
    if n_genes < 1:
        n_genes = 1
    elif n_genes > ws_gex_constants.MAX_TOP_N_GENES:
        n_genes = ws_gex_constants.MAX_TOP_N_GENES

    cols = [
        {'type': 'string', 'label': 'Gene ID'},
        {'type': 'string', 'label': 'Gene name'},
    ]

    for i in xrange(n_clusters):
        # Filter genes by mean count and sort by log2 fold-change, descending
        means =   diff_expr.data[:, 0+3*i]
        log2fcs = diff_expr.data[:, 1+3*i]

        keep_indices = np.flatnonzero(means >= ws_gex_constants.TOP_DE_GENES_MIN_MEAN)
        top_gene_indices = keep_indices[log2fcs[keep_indices].argsort()[::-1]][:n_genes]

        for j in top_gene_indices:
            top_genes.add(analysis.matrix.int_to_gene_id(j))

        cols.append({'type': 'number',
                     'label': 'L2FC',
                     'title': 'Log2 fold-change in cluster %d vs other cells' % (i+1)})
        cols.append({'type': 'number',
                     'label': 'p-value',
                     'title': 'Adjusted p-value of differential expression in cluster %d' % (i+1)})

    rows = []
    for gene_id in top_genes:
        i = analysis.matrix.gene_id_to_int(gene_id)
        gene_name = analysis.matrix.gene_id_to_name(gene_id)

        row = [gene_id, gene_name]
        for j in xrange(n_clusters):
            log2fc =      diff_expr.data[i, 1+(3*j)]
            adj_p_value = diff_expr.data[i, 2+(3*j)]

            if log2fc <= 0 or adj_p_value >= ws_gex_constants.PVALUE_DEEMPHASIS_CUTOFF:
                style = '#DDD'
            else:
                style = '#000'

            row.append({'v': tk_safe_json.json_sanitize(log2fc),
                        'f': format_value(log2fc, '%.2f'),
                        's': style})
            row.append({'v': tk_safe_json.json_sanitize(adj_p_value),
                        'f': format_value(adj_p_value, '%.0e'),
                        's': style})


        rows.append(row)

    # Sort by log2fc, descending, in first cluster
    if n_clusters > 0:
        rows = sorted(rows, key=lambda row: row[2]['v'], reverse=True)

    chart['table'].update({'rows': rows, 'cols': cols})
    return chart

def plot_differential_expression(chart, sample_properties, sample_data):
    return clustering_plot_func(chart, sample_properties, sample_data, _plot_differential_expression, [sample_data.analysis])

def clustering_plot_func(chart, sample_properties, sample_data, plot_func, args=[], kwargs={}):
    if len(sample_properties['genomes']) > 1:
        return None

    analysis = sample_data.analysis
    if analysis is None:
        return None

    new_charts = []
    for clustering_key, clustering in analysis.clusterings.iteritems():
        kwargs['clustering'] = clustering
        kwargs['original_cluster_sizes'] = sample_data.original_cluster_sizes[clustering_key]
        kwargs['diff_expr'] = analysis.differential_expression[clustering_key]
        new_chart = plot_func(copy.deepcopy(chart), *args, **kwargs)
        if new_chart is not None:
            new_chart['filters'] = {ws_gex_constants.CLUSTERS_FILTER_TITLE: clustering.description}
            new_charts.append(new_chart)
    return new_charts

def make_chart_filters(sample_properties, analysis):
    if analysis is None or len(sample_properties['genomes']) > 1 or not isinstance(analysis, cr_sg_analysis.SingleGenomeAnalysis):
        return {}

    filter_values = map(lambda x: x.description, cr_clustering.sort_clusterings(analysis.clusterings.values()))

    return {
        ws_gex_constants.CLUSTERS_FILTER_TITLE: {
            'values': filter_values,
            'selected': filter_values[0],
        },
    }

def plot_subsampled_scatterplot_metric(chart, sample_properties, sample_data, **kwargs):
    """
    Modifies chart data entry to add traces for subsampled metrics with specified suffixes.
    Metrics must take the form <reference>_<subsample_type>_<subsample_depth>_<metric_suffix>, where metric suffix is specified as a kwarg

    KWARGS:
        metric_suffix (str): suffix for the subsampled metric given a metric of the form <reference>_<subsample_type>_<subsample_depth>_<metric_suffix>

    """
    summary_data = sample_data.summary or {}

    subsample_type = kwargs.get('subsample_type', 'raw_rpc')
    ref_prefix = kwargs.get('ref_prefix', '(.+)_')

    # Regular expression to match <reference>_<subsample_type>_<subsample_depth>_<metric_suffix> pattern
    metric_pattern = '^%s(%s)_([0-9]+)_%s' % (ref_prefix,
                                              subsample_type,
                                              kwargs.get('metric_suffix', None))

    metric_search_results = [re.search(metric_pattern, key) for key in summary_data.keys()]

    # Only return chart when multiple references are present for multi-genome only metrics
    if kwargs.get('multi_genome_only', False):
        references = set([result.group(1) for result in metric_search_results if result is not None])
        if len(references) <= 2:
            return

    # Find relevant metrics and extract data
    points = []

    for key, search_result in zip(summary_data.keys(), metric_search_results):
        # Only count corresponding keys and non-zero values
        if search_result and summary_data[key] > 0:
            genome = search_result.group(1)
            subsample_type = search_result.group(2).replace('_', ' ')

            if kwargs.get('show_multi_genome_only', False) and (genome != 'multi'):
                continue
            subsample_depth = int(search_result.group(3))

            if genome == 'multi':
                trace_label = ''
            else:
                trace_label = genome

            points.append((subsample_depth,
                           float(summary_data[key]),
                           trace_label,
                           genome,
                           subsample_type,
                       ))

    # Don't produce chart if there's no data
    if len(points) == 0:
        return None

    # Sort extracted points by subsample depth (makes line plot connect properly)
    sorted_points = sorted(points, key=lambda x: (x[2], x[0]))

    # Build up a set of traces for <reference>_<subsample_type> pairs
    traces = {}

    for point in sorted_points:
        depth, value, trace_label, genome, subsample_type = point
        trace_key = (genome, subsample_type)

        if trace_key not in traces:
            # add (0,0)
            traces[trace_key] = {
                'x': [0, depth],
                'y': [0, value],
                'name': trace_label,
                'mode': 'lines',
                'line': {
                    'width': 3,
                },
            }
        else:
            traces[trace_key]['x'].append(depth)
            traces[trace_key]['y'].append(value)

    # Set extents of saturation line
    if chart.get('layout') and chart['layout'].get('shapes') and chart['layout']['shapes'][0]['type'] == 'line':
        chart['layout']['shapes'][0]['x1'] = max([max(trace['x']) for trace in traces.itervalues()])

    # Add new data into chart
    chart['data'] = [v for k,v in sorted(traces.items())]
    return chart

def build_charts(sample_properties, chart_dicts, sample_data, module=None):
    modules = [module, globals()] if module else [globals()]

    filters = make_chart_filters(sample_properties, sample_data.analysis)

    charts = []
    for chart_dict in chart_dicts:
        chart_dict = copy.deepcopy(chart_dict)
        function = chart_dict.pop('function')
        for module in modules:
            f = module.get(function)
            if f is not None:
                break
        kwargs = chart_dict.pop('kwargs', {})

        new_chart_obj = f(chart_dict, sample_properties, sample_data, **kwargs)
        if new_chart_obj is None:
            continue

        new_charts = new_chart_obj if isinstance(new_chart_obj, list) else [new_chart_obj]
        charts.extend(new_charts)

    return charts, filters

def filter_vdj_prefixes(all_prefixes, sample_properties):
    """ Only get subset of metric prefix values """
    chain_filter = sample_properties.get('chain_type')
    if chain_filter is None:
        return all_prefixes

    # NOTE: Assumes chains (TRA, TRB, etc) are prefixed with the chain_type (TR or IG)
    result = {}
    for key, values in all_prefixes.iteritems():
        # Only filter any prefix that is a candidate for filtering
        # (i.e., contains some values prefixed by the selected chain type)
        if values is not None and \
           any(v.startswith(chain_filter) for v in values):
            result[key] = [v for v in values if v.startswith(chain_filter) or v == cr_constants.MULTI_REFS_PREFIX]
        else:
            result[key] = values
    return result


def filter_vdj_alarms(all_alarms, sample_properties):
    """ Only get subset of metric alarms """
    chain_filter = sample_properties.get('chain_type')
    if chain_filter is None:
        return all_alarms

    result = []
    for alarm in all_alarms:
        # No filters specified; don't filter
        if 'filters' not in alarm:
            result.append(alarm)
            continue

        for f in alarm['filters']:
            if f.get('chain_type') == chain_filter:
                result.append(alarm)

    return result


def get_constants_for_pipeline(pipeline, sample_properties):
    """ Get the appropriate metrics/alarms/charts for a pipeline """
    if pipeline == shared_constants.PIPELINE_VDJ:
        metrics, alarms, charts = ws_vdj_constants.METRICS, ws_vdj_constants.METRIC_ALARMS, ws_vdj_constants.CHARTS

        metric_prefixes = filter_vdj_prefixes(vdj_report.VdjReporter().get_all_prefixes(),
                                              sample_properties)

        alarms = filter_vdj_alarms(alarms, sample_properties)

    else:
        metrics, alarms, charts = ws_gex_constants.METRICS, ws_gex_constants.METRIC_ALARMS, ws_gex_constants.CHARTS

        metric_prefixes = cr_report.Reporter().get_all_prefixes()

    return metrics, alarms, charts, metric_prefixes


def build_web_summary_json(sample_properties, sample_data, pipeline):
    """ sample_properties - dict
        sample_data - *SampleData class
        pipeline - string """
    view = copy.deepcopy(sample_properties)

    metrics, alarms, charts, all_prefixes = get_constants_for_pipeline(pipeline, sample_properties)

    tables, alarms = build_tables(sample_properties, metrics, alarms, sample_data, all_prefixes=all_prefixes)
    if tables:
        view['tables'] = tables
    if alarms:
        view['alarms'] = alarms

    charts, filters = build_charts(sample_properties, charts,
                                   sample_data=sample_data)
    if charts:
        view['charts'] = charts
    if filters:
        view['filters'] = filters

    # Selected metrics that the web summary template needs
    info = build_info_dict(sample_properties, sample_data, pipeline)
    if info:
        view['info'] = info

    return view


def build_info_dict(sample_properties, sample_data, pipeline):
    """ Add miscellaneous metrics required by the web summary template """
    if sample_data.summary is None:
        return None
    info = {}

    info['chemistry_description'] = sample_data.summary.get('chemistry_description')

    info['references'] = []

    if pipeline in [shared_constants.PIPELINE_AGGR, shared_constants.PIPELINE_REANALYZE]:
        genomes = sample_properties.get('genomes')
        if genomes is not None:
            info['references'].append({
                'type': cr_constants.REFERENCE_TYPE,
                'name': cr_reference.get_ref_name_from_genomes(genomes)
            })

    else:
        # Find all references in the summary
        reference_metric_prefixes = [cr_constants.REFERENCE_METRIC_PREFIX,
                                     vdj_constants.REFERENCE_METRIC_PREFIX]

        for prefix in reference_metric_prefixes:
            type_metric = '%s%s' % (prefix, cr_constants.REFERENCE_TYPE_KEY)

            if type_metric not in sample_data.summary:
                continue

            ref_type = sample_data.summary[type_metric]

            name_metric = '%s%s' % (prefix, cr_constants.REFERENCE_GENOMES_KEY)
            if name_metric not in sample_data.summary:
                raise ValueError("Reference metadata metric %s not found in metrics summary." % name_metric)
            ref_name = sample_data.summary.get(name_metric)

            if isinstance(ref_name, list):
                ref_name = cr_reference.get_ref_name_from_genomes(ref_name)

            info['references'].append({
                'type': ref_type,
                'name': ref_name
            })

    return info


def build_web_summary_html(filename, sample_properties, sample_data, pipeline,
                           template_dir=None, alerts_output_filename=None):
    view = build_web_summary_json(sample_properties, sample_data, pipeline)

    if not view:
        return

    with open(filename, 'w') as f:
        f.write(template.convert_webshim_json_to_html(view, pipeline, template_dir=template_dir))

    # Write raised alerts to a json file
    if alerts_output_filename is not None:
        with open(alerts_output_filename, 'w') as f:
            json.dump(tk_safe_json.json_sanitize(view.get('alarms', [])), f, indent=4, sort_keys=True)


def build_metrics_summary_csv(filename, sample_properties, sample_data, pipeline):
    metrics, alarms, charts, all_prefixes = get_constants_for_pipeline(pipeline, sample_properties)

    tables, _ = build_tables(sample_properties, metrics, alarms, sample_data, all_prefixes=all_prefixes)
    if not tables:
        sys.stderr.write("No metrics tables were generated, skipping CSV generation.\n")
        return

    csv_metrics = collections.OrderedDict()
    for table in tables:
	if not table:
		continue
        for metric, _, value in table['rows']:
            if type(metric) == dict:
                metric = metric['v']
            if type(value) == dict:
                value = value['v']
            if metric not in csv_metrics:
                csv_metrics[metric] = value

    with open(filename, 'wb') as f:
        writer = csv.writer(f)
        writer.writerow(csv_metrics.keys())
        writer.writerow(csv_metrics.values())


def convert_to_int_gracefully(s):
    try:
        return int(s)
    except ValueError:
        return sys.maxint

def convert_to_float_gracefully(s):
    try:
        return float(s)
    except ValueError:
        return sys.float_info.max
