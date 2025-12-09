import { useMemo } from "react";

import { MetricTableConfig } from "../../config/types";
import { getMetricHelpText } from "../../MetricsTable/utils";
import {
  type ExperimentalDesign,
  type JsonMetricSummary,
  metricIsEmpty,
  metricIsFraction,
  metricIsNumber,
} from "../../types";
import { getFloatFromFraction } from "../../utilities";
import { getAlertThreshold } from "../alertUtils";
import {
  groupMultiSampleMetricsByKey,
  type MetricWithBarcodes,
} from "../hooks";
import { isPlateBased, parsePlateBarcodeId } from "./plateUtils";
import { SampleMetricsBarcodes } from "./utils";

type DescribedValue = { description?: string; value: number };
type Threshold = ReturnType<typeof getAlertThreshold>;
/**
 * Groups a plot's values with a title for the plot, e.g.
 * Plate A and its metrics.
 * In the non plate case, there is just one group with an empty
 * title for the entire run.
 */
type Group = { title: string; values: DescribedValue[] };

function getValue(metric: JsonMetricSummary) {
  if (metricIsFraction(metric)) {
    return getFloatFromFraction(metric.value);
  }

  if (!metricIsNumber(metric)) {
    return null;
  }

  return metric.value;
}

function isNumericOrEmpty(
  metricsWithBarcodes: ReadonlyArray<MetricWithBarcodes>,
) {
  return metricsWithBarcodes.every(
    ({ metric }) =>
      metricIsFraction(metric) ||
      metricIsEmpty(metric) ||
      metricIsNumber(metric),
  );
}

/**
 * For each unique metric type across samples, maps a plate to that metric's values
 * and then outputs data that can be used in a beeswarm plot.
 * E.g. the current metric we're working on is 'total singlets' in the caller,
 * iterate over all samples' 'total_singlets' values, group them by plate via
 * the barcode so that: {A: [1, 2, 3], B: [4, 5, 6], etc.} which can then be assigned
 * back to 'total_singlets' in the caller.
 */
function metricsGroupsByPlate(
  barcodesWithMetrics: ReadonlyArray<MetricWithBarcodes>,
) {
  // For each unique metric type across samples, maps a plate to its metric values.
  type PlateToValues = Record<string, DescribedValue[]>;
  const plateWithValues = barcodesWithMetrics.reduce<PlateToValues>(
    (plateMap, sampleMetric) => {
      const { barcodes, metric, sampleId } = sampleMetric;
      const value = getValue(metric);
      if (value === null) {
        // e.g. If it's a string or empty metric.
        return plateMap;
      }

      barcodes.forEach((barcodeId) => {
        const description = [
          "Sample ID:<br>",
          `${sampleId}<br><br>`,
          "Value:<br>",
          `${value}`,
        ].join("");

        const chunks = parsePlateBarcodeId(barcodeId);
        if (!chunks) {
          return;
        }
        const { plate } = chunks;
        const id = `Plate ${plate}`;

        plateMap[id] = [...(plateMap[id] ?? []), { description, value }];
      });
      return plateMap;
    },
    {},
  );

  const groups = Object.keys(plateWithValues).reduce<Group[]>(
    (groups, entity) => [
      ...groups,
      { title: entity, values: plateWithValues[entity] },
    ],
    [],
  );
  return groups;
}

/**
 * For each unique metric type across samples, maps a barcode to that metric's values
 * and then outputs data that can be used in a beeswarm plot.
 * Should be used when we're not dealing with plate-based data.
 */
function createSingleMetricGroup(
  barcodesWithMetrics: ReadonlyArray<MetricWithBarcodes>,
) {
  // For each unique metric type across samples, maps the entire run to its metric values.
  type RunValues = DescribedValue[];

  const runWithValues = barcodesWithMetrics.reduce<RunValues>(
    (runMetrics, sampleMetric) => {
      const { metric, sampleId } = sampleMetric;
      const value = getValue(metric);
      if (value === null) {
        // e.g. If it's a string or empty metric.
        return runMetrics;
      }

      const description = [
        "Sample ID:<br>",
        `${sampleId}<br><br>`,
        "Value:<br>",
        `${value}`,
      ].join("");

      return [...runMetrics, { description, value }];
    },
    [],
  );

  return [{ title: "", values: runWithValues }];
}

/**
 * Produce metrics that can be visualized in a beeswarm plot.
 * Configs are used to setup the help text for each metric.
 * Metric keys determine for which metrics we generate data.
 */
export const useBeeswarmMetricsData = (
  sampleMetricsBarcodes: SampleMetricsBarcodes,
  configs: MetricTableConfig[],
  experimentalDesign: Pick<ExperimentalDesign, "is_rtl" | "include_introns">,
  metricKeys: string[],
) => {
  return useMemo(() => {
    const data: Array<{
      groups: Array<Group>;
      helpText: string;
      key: string;
      title: string;
      typedThreshold: Threshold;
    }> = [];

    metricKeys.forEach((key) => {
      const barcodesMetrics = groupMultiSampleMetricsByKey(
        sampleMetricsBarcodes,
        key,
      );
      if (barcodesMetrics.length === 0) {
        console.warn(`no metrics produced for key ${key}`);
        return;
      }

      // The config should be consistent across metrics of the same type.
      const config = barcodesMetrics[0].metric.config;

      // Note:
      // Similarly to usePlateMetricsData, we can't plot non-numeric metrics,
      // so we should skip them early on. Unlike usePlateMetricsData, we have
      // to do a direct check.
      if (!isNumericOrEmpty(barcodesMetrics)) {
        return;
      }

      let helpText = "";
      for (const tableConfig of configs) {
        helpText = getMetricHelpText({ config, key }, tableConfig, {
          is_rtl: experimentalDesign.is_rtl,
        });
        if (helpText) {
          break;
        }
      }

      const plateBased = isPlateBased(
        experimentalDesign,
        barcodesMetrics[0].barcodes,
      );
      const groups = plateBased
        ? metricsGroupsByPlate(barcodesMetrics)
        : createSingleMetricGroup(barcodesMetrics);

      data.push({
        groups,
        helpText,
        key,
        title: config.header,
        typedThreshold: getAlertThreshold(config.alerts, experimentalDesign),
      });
    });

    return data;
  }, [sampleMetricsBarcodes, configs, experimentalDesign, metricKeys]);
};
