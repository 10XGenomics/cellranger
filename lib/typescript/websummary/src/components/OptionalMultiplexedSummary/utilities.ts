import uniqWith from "lodash/uniqWith";

import expandPlots from "../../shared/websummary/components/CellRanger/PreProcess";
import { MetricEntryConfig } from "./config/types";
import {
  Fraction,
  type JsonMetricSummary,
  type LibraryType,
  type LibraryWebSummary,
  metricIsFraction,
  metricIsNumber,
  metricIsString,
  MetricValueType,
  type SampleWebSummary,
} from "./types";

/**
 * Some plots have data arranged in a way that can't be rendered by the components by default.
 * This preprocesses that data to get it ready for those components. This isn't idempotent,
 * if it's called more than once, there will be errors.
 * @param {object} data the top level data passed to this component.
 */
export const preprocess = (data: SampleWebSummary) => {
  // TODO
  // Keypaths may need to be updated once data is finalized.
  // This is the only plot that won't require reassignment.
  const plot = data.gex_tab?.content.clustering_and_diffexp_plots?.right_plots;
  if (plot) {
    expandPlots(plot);
  }
  const plot2 =
    data.antibody_tab?.content.clustering_and_diffexp_plots?.right_plots;
  if (plot2) {
    expandPlots(plot2);
  }
  const plot3 =
    data.antigen_tab?.content.clustering_and_diffexp_plots?.right_plots;
  if (plot3) {
    expandPlots(plot3);
  }
};

/**
 * Given a category, library type, and a list of metrics, returns the metrics
 * that are associated with the library type.
 */
export const getMetricsForLibraryType = (
  libraryType: LibraryType,
  data: JsonMetricSummary[],
) => {
  const filteredData = data.filter(
    (metric) => metric.library_type === libraryType,
  );
  return uniqWith(
    filteredData,
    (a, b) => a.key === b.key && a.grouping_key === b.grouping_key,
  );
};

/**
 * Formatting utility functions for metrics.
 */

const truncateFloat = (value: number, precision: number) => {
  const multiplier = Math.pow(10, precision || 0);
  return Math.round(value * multiplier) / multiplier;
};

function getToLocaleString(value: number, precision: number = 0) {
  return value.toLocaleString("en-US", {
    maximumFractionDigits: precision,
    minimumFractionDigits: precision,
  });
}

export function getFloatFromFraction(value: Fraction) {
  return value.numerator / value.denominator;
}

function getPercentageFromFloat(value: number, precision: number = 0) {
  return `${truncateFloat(value * 100, precision)}%`;
}

export function formatNumercValue(metric: {
  config: {
    type: MetricValueType;
  };
  value: number;
}) {
  switch (metric.config.type) {
    case "usize":
      return getToLocaleString(metric.value);
    case "FloatAsInt":
      return getToLocaleString(metric.value, 0);
    case "f64":
      return getToLocaleString(metric.value, 1);
    case "Percent":
      return getPercentageFromFloat(metric.value, 1);
  }
}

/**
 * Given a metric, returns a formatted value based on the metric's type.
 */
export function getFormattedValue(metric: JsonMetricSummary) {
  if (metricIsNumber(metric)) {
    return formatNumercValue(metric);
  }

  if (metricIsFraction(metric)) {
    const numeratorStr = getToLocaleString(metric.value.numerator);
    const percentage =
      metric.value.denominator > 0
        ? getPercentageFromFloat(getFloatFromFraction(metric.value), 1)
        : "---%";
    return `${numeratorStr} (${percentage})`;
  }

  if (metricIsString(metric)) {
    return metric.value;
  }

  return "-";
}

/**
 * Given a single key and a list of metrics, return the first metric that matches.
 *
 * Return undefined if there is no match.
 */
export const getSingleMetric = (
  key: MetricEntryConfig,
  data: JsonMetricSummary[],
) => {
  const matches = getRowMetrics([key], data);
  return matches.length > 0 ? matches[0] : undefined;
};

/**
 * Given a list of keys and a list of metrics, returns the metrics that are
 * associated with the keys.
 */
export const getRowMetrics = (
  keys: MetricEntryConfig[],
  data: JsonMetricSummary[],
) => {
  return keys
    .map((entry) => {
      return data.find(
        (metric) =>
          (entry.category === undefined ||
            entry.category === metric.category) &&
          metric.key === entry.key,
      );
    })
    .filter(Boolean) as JsonMetricSummary[];
};

export const hasAlerts = (data: LibraryWebSummary | SampleWebSummary) => {
  return Object.keys(data).some((key: keyof typeof data) => {
    const value = data[key];
    if (value && typeof value !== "string" && "alerts" in value) {
      return value.alerts.length > 0;
    }
    return false;
  });
};
