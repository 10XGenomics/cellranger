import { interpolateMagma } from "d3-scale-chromatic";
import { useMemo } from "react";

import { MetricTableConfig } from "../../config/types";
import { getMetricHelpText } from "../../MetricsTable/utils";
import {
  ExperimentalDesign,
  JsonMetricSummary,
  MetricConfig,
  metricIsFraction,
  metricIsNumber,
} from "../../types";
import { formatNumercValue, getFormattedValue } from "../../utilities";
import { getAlertThreshold } from "../alertUtils";
import {
  groupMultiSampleMetricsByKey,
  type MetricWithBarcodes,
} from "../hooks";
import { ColorScale, PlateDataMap } from "../types";
import { defaultPlateData, parsePlateBarcodeId } from "./plateUtils";
import { SampleMetricsBarcodes } from "./utils";

type Bounds = NonNullable<ReturnType<typeof getBounds>>;
type Threshold = ReturnType<typeof getAlertThreshold>;
type ColorInterpolator = (x: number) => string;

function getBounds(
  config: MetricConfig,
  sampleMetrics: { metric: JsonMetricSummary }[],
) {
  let min = Infinity;
  let max = -Infinity;

  switch (config.type) {
    case "usize":
    case "FloatAsInt":
    case "f64": {
      // For these, we need to iterate through the values to find the min/max and
      // map the colors to them.
      sampleMetrics.forEach(({ metric }) => {
        if (metricIsNumber(metric)) {
          min = Math.min(min, metric.value);
          max = Math.max(max, metric.value);
        }
      });
      break;
    }
    case "Percent": {
      // For this, it's already a normalized percent, so we know min is 0, max is 1.
      min = 0;
      max = 1;
      break;
    }
    case "CountAndPercent": {
      // For this, we have a numerator/denominator where the the quotient makes
      // the percent, but we need to determine the denominator as a our max,
      // and the min is the smallest numerator.
      sampleMetrics.forEach(({ metric }) => {
        if (metricIsFraction(metric)) {
          min = Math.min(min, metric.value.numerator);
          max = Math.max(max, metric.value.denominator);
        }
      });
      break;
    }

    // Anything else is a non-numeric value and we can't handle (e.g. 'string')
    default:
      break;
  }

  if (!Number.isFinite(min) || !Number.isFinite(max)) {
    return null;
  }

  // If every value is the same, then there is no max or min, so we should
  // put some sane placeholders in.
  min = min === max ? 0 : min;
  max = min === max ? 1 : max;

  return { max, min, type: config.type };
}

function getValue(metric: JsonMetricSummary) {
  if (metricIsFraction(metric)) {
    return metric.value.numerator;
  }

  if (!metricIsNumber(metric)) {
    return null;
  }

  return metric.value;
}

/**
 * Produces one group of PlateWell visualizations for one metric key, e.g.
 * "total_singlets": {plate1: colors, plate2: colors, etc.}
 */
function metricsPlates(
  bounds: Bounds,
  colorInterpolator: ColorInterpolator,
  metricsWithBarcodes: ReadonlyArray<MetricWithBarcodes>,
) {
  const platesWellColors = metricsWithBarcodes.reduce<PlateDataMap>(
    (plateMap, sampleMetric) => {
      const { barcodes, metric, sampleId } = sampleMetric;
      barcodes.forEach((barcodeId) => {
        const chunks = parsePlateBarcodeId(barcodeId);
        if (!chunks) {
          return plateMap;
        }

        const { plate, column, row } = chunks;
        const well = `${row}${column}` as const;

        const value = getValue(metric);
        if (value === null) {
          // e.g. If it's a string or empty metric.
          return plateMap;
        }

        if (!plateMap[plate]) {
          plateMap[plate] = defaultPlateData();
        }

        const formattedValue = getFormattedValue(metric);

        const normalizedValue =
          (value - bounds.min) / (bounds.max - bounds.min);
        const color = colorInterpolator(normalizedValue);
        const tooltip = [
          `Sample barcode ID: ${barcodeId}<br>`,
          `Sample ID: ${sampleId}<br><br>`,
          `Value: ${formattedValue}<br>`,
        ].join("");

        plateMap[plate][well] = { color, tooltip };
      });
      return plateMap;
    },
    {},
  );

  return platesWellColors;
}

function colorScale(bounds: Bounds, colorInterpolator: ColorInterpolator) {
  const chunkCount = 32;
  const midpointIndex = Math.round((chunkCount - 1) * 0.5);
  const scale: ColorScale = [];

  function getLabel(i: number) {
    if (i === 0) {
      return bounds.min;
    } else if (i === midpointIndex) {
      return (bounds.min + bounds.max) * 0.5;
    } else if (i === chunkCount - 1) {
      return bounds.max;
    }
    return undefined;
  }

  function formatLabel(numericLabel: number | undefined) {
    if (numericLabel === undefined) {
      return undefined;
    }
    const formatted = formatNumercValue({
      config: bounds,
      value: numericLabel,
    });
    return formatted ?? numericLabel.toFixed(1);
  }

  for (let i = 0; i < chunkCount; i++) {
    // Normalized value from 0 to 1 because this is what
    // the color interpolator expects.
    const normalized = i / (chunkCount - 1);
    const label = formatLabel(getLabel(i));
    scale.push({
      color: colorInterpolator(normalized),
      label,
    });
  }

  return scale;
}

/**
 * Produce metrics that can be visualized in a plate setting (similar to usePlateAlerts),
 * where a numeric value of the metric then maps to some color value from a 2D
 * color scale.
 * Configs are used to setup the help text for each metric.
 */
export const usePlateMetricsData = (
  sampleMetricsBarcodes: SampleMetricsBarcodes,
  configs: MetricTableConfig[],
  experimentalDesign: Pick<ExperimentalDesign, "is_rtl" | "include_introns">,
  metricKeys: string[],
) => {
  return useMemo(() => {
    const data: Array<{
      bounds: Bounds;
      colorScale: ColorScale;
      helpText: string;
      key: string;
      plateDataMap: PlateDataMap;
      title: string;
      typedThreshold: Threshold;
    }> = [];

    const colorInterpolator = interpolateMagma;

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
      // We also want to filter out non-numeric metrics here.
      // getBounds is basically doing double duty here, because when
      // it returns null it's likely that it's because the metric is
      // non-numeric.
      const bounds = getBounds(config, barcodesMetrics);
      if (!bounds) {
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

      data.push({
        bounds,
        colorScale: colorScale(bounds, colorInterpolator),
        helpText,
        key,
        plateDataMap: metricsPlates(bounds, colorInterpolator, barcodesMetrics),
        title: config.header,
        typedThreshold: getAlertThreshold(
          config.alerts,
          experimentalDesign,
          bounds,
        ),
      });
    });

    return data;
  }, [sampleMetricsBarcodes, configs, experimentalDesign, metricKeys]);
};
