import groupBy from "lodash/groupBy";

import {
  ConditionalMetricHelp,
  MetricEntryConfig,
  MetricTableConfig,
} from "../config/types";
import {
  type ExperimentalDesign,
  type JsonMetricSummary,
  metricIsEmpty,
  metricIsFraction,
  metricIsNumber,
  metricIsString,
} from "../types";
import { getFloatFromFraction } from "../utilities";

export const getRenderedHelpText = (text: string) => {
  const renderHelpText = () => {
    const urlRegex = /(https?:\/\/[^\s]+)/g;
    return text
      .split(urlRegex)
      .map((part) => {
        if (urlRegex.test(part)) {
          return `<a href="${part}" class="tooltip-link" target="_blank" rel="noopener noreferrer">${part}</a>`;
        }
        return part;
      })
      .join("");
  };
  return `<span>${renderHelpText()}</span>`;
};

/**
 * Get help text for the provided metric, possibly using an override from the table config.
 */
export const getMetricHelpText = (
  metric: Pick<JsonMetricSummary, "config" | "key">,
  config: Pick<MetricTableConfig, "entries" | "help">,
  experimentalDesign: Pick<ExperimentalDesign, "is_rtl">,
): string => {
  // FIXME: we should get the metric view config when we filter the metrics, rather
  // than searching for it here.
  const configuredHelp = config.entries.find(
    (entry) => entry.key === metric.key,
  )?.help;

  if (configuredHelp) {
    if (typeof configuredHelp === "string") {
      return configuredHelp;
    }
    const condHelp = selectConditionalHelp(configuredHelp, experimentalDesign);
    if (condHelp !== null) {
      return condHelp;
    }
  }
  // No help text override, pass through.
  return metric.config.help || "";
};

/**
 * Check each conditional metric help option against the experimental design.
 * Return the first help text that matches, or null if none match (configurations
 * should ensure this never happens).
 */
const selectConditionalHelp = (
  options: ConditionalMetricHelp[],
  experimentalDesign: Pick<ExperimentalDesign, "is_rtl">,
): string | null => {
  for (const opt of options) {
    if (optMatchesDesign(opt, experimentalDesign)) {
      return opt.msg;
    }
  }
  return null;
};

/**
 * Return true if the conditions for the provided help text match the experimental design.
 */
const optMatchesDesign = (
  opt: ConditionalMetricHelp,
  experimentalDesign: Pick<ExperimentalDesign, "is_rtl">,
): boolean => {
  let matches = true;
  if (opt.isRtl !== undefined) {
    matches &&= opt.isRtl === experimentalDesign.is_rtl;
  }
  return matches;
};

function compareMetrics(
  a: JsonMetricSummary,
  b: JsonMetricSummary,
  direction: "asc" | "desc",
) {
  // values should always sort before null
  if (!metricIsEmpty(a) && metricIsEmpty(b)) {
    return -1;
  } else if (metricIsEmpty(a) && !metricIsEmpty(b)) {
    return 1;
  }

  if (metricIsString(a) && metricIsString(b)) {
    return direction === "asc"
      ? a.value.localeCompare(b.value)
      : b.value.localeCompare(a.value);
  }

  if (metricIsNumber(a) && metricIsNumber(b)) {
    return direction === "asc" ? a.value - b.value : b.value - a.value;
  }

  if (metricIsFraction(a) && metricIsFraction(b)) {
    const fa = getFloatFromFraction(a.value);
    const fb = getFloatFromFraction(b.value);
    return direction === "asc" ? fa - fb : fb - fa;
  }

  console.warn("Attempting to compare metrics of different types", a, b);
  return 0;
}

/**
 * Sorts table row by `orderBy` and `direction`.
 */
export function sortMultiRowMetrics(
  tableRows: JsonMetricSummary[][],
  orderBy: string[],
  direction: "asc" | "desc",
) {
  const sorted = [...tableRows];
  return sorted.sort((a, b) => {
    for (const field of orderBy) {
      const aMetric = a.find((metric) => metric.key === field);
      const bMetric = b.find((metric) => metric.key === field);
      if (aMetric !== undefined && bMetric !== undefined) {
        return compareMetrics(aMetric, bMetric, direction);
      }
    }
    return 0;
  });
}

/**
 * Takes metrics and groups them together based on grouping keys,
 * while also filtering out non-table metrics.
 */
export const processMultiRowMetrics = (
  metrics: JsonMetricSummary[],
  config: Pick<MetricTableConfig, "entries">,
): JsonMetricSummary[][] => {
  const { entries } = config;

  // use the first entry in table to collect grouping keys
  const groupingKeys = metrics
    .filter((metric) => metric.key === entries[0].key)
    .map((metric) => metric?.grouping_key);

  // filter out non-table metrics
  const filteredMetrics = metrics.filter(
    (metric) =>
      entries.some((entry) => entry.key === metric.key) &&
      groupingKeys.includes(metric.grouping_key),
  );

  // group metrics into table rows by grouping key
  const groupedRows = groupBy(filteredMetrics, "grouping_key");

  // populate table structure, insert null placeholders where necessary
  const tableRows = cleanTableStructure(entries, Object.values(groupedRows));
  return tableRows;
};

/**
 * Given rows of table metrics and the table configuration, clean the rows.
 *
 * Ensure that every row has the same number of entries - insert placeholders
 * into rows that are missing a value for a particular metric for which at least
 * one other row has a value.
 */
export const cleanTableStructure = (
  entries: MetricEntryConfig[],
  tableRows: JsonMetricSummary[][],
): JsonMetricSummary[][] => {
  // identify metrics that do not have values in any rows
  // also get the metric configs for generating placeholder entries for missing metrics
  const missingMetrics = new Set(entries.map((entry) => entry.key));
  const exemplarMetrics = new Map();

  for (const row of tableRows) {
    if (missingMetrics.size === 0) {
      break;
    }
    for (const metric of row) {
      missingMetrics.delete(metric.key);
      exemplarMetrics.set(metric.key, metric);
    }
  }

  // trim down to only metrics that we have at least one instance of
  const filteredEntries = entries.filter(
    (entry) => !missingMetrics.has(entry.key),
  );

  const createPlaceholder = (key: string): JsonMetricSummary => {
    const exemplar = exemplarMetrics.get(key);
    if (!exemplar) {
      // This should never happen.
      console.error("Missing exemplar metric for key", key);
      return {
        alerts: [],
        category: "Library",
        config: {
          alerts: [],
          header: "",
          type: "String",
        },
        key,
        library_type: "Gene Expression",
        value: null,
      };
    }
    return placeholderMetric(exemplar);
  };

  return tableRows.map((row) =>
    filteredEntries.map(
      (entry) =>
        row.find((metric) => metric.key === entry.key) ||
        createPlaceholder(entry.key),
    ),
  );
};

export const placeholderMetric = (
  template: JsonMetricSummary,
): JsonMetricSummary => ({
  ...template,
  alerts: [],
  config: {
    ...template.config,
    type: "String",
  },
  grouping_key: null,
  value: null,
});
