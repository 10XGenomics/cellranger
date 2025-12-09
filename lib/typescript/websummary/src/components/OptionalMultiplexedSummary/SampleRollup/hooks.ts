import { useMemo } from "react";

import {
  JsonMetricSummary,
  LibraryType,
  MultiWebSummarySampleData,
  SampleWebSummary,
} from "../types";
import { getMetricsForLibraryType, hasAlerts } from "../utilities";
import { getAlertSeverity, getAlertSeverityFromAlerts } from "../utils";
import { SampleMetricsBarcodes } from "./SecitonWithVisualizations/utils";

/**
 * Return the rolled-up alert severity for a specific tab.
 */
const getAlertSeverityForTab = (
  sampleData: MultiWebSummarySampleData[],
  tabKey: keyof Omit<
    SampleWebSummary,
    "id" | "description" | "multiplexing_barcode_ids"
  >,
) => {
  return getAlertSeverity(sampleData, (data) => {
    const alerts = data.data[tabKey]?.alerts || [];
    return getAlertSeverityFromAlerts(alerts);
  });
};

export const useSampleHasAlerts = (data: MultiWebSummarySampleData[]) =>
  useMemo(() => data.map((data) => hasAlerts(data.data)), [data]);

export const useLibraryAlertSeverity = (data: MultiWebSummarySampleData[]) =>
  useMemo(
    /* eslint-disable sort-keys */
    () => ({
      gex: getAlertSeverityForTab(data, "gex_tab"),
      vdjT: getAlertSeverityForTab(data, "vdj_t_tab"),
      vdjTGD: getAlertSeverityForTab(data, "vdj_t_gd_tab"),
      vdjB: getAlertSeverityForTab(data, "vdj_b_tab"),
      antibody: getAlertSeverityForTab(data, "antibody_tab"),
      antigen: getAlertSeverityForTab(data, "antigen_tab"),
      crispr: getAlertSeverityForTab(data, "crispr_tab"),
      customFeature: getAlertSeverityForTab(data, "custom_feature_tab"),
      /* eslint-enable sort-keys */
    }),
    [data],
  );

export const useFilteredPerSampleData = (
  data: MultiWebSummarySampleData[],
  filter: string,
  sampleHasAlerts: boolean[],
  showOnlyWithAlerts: boolean,
) =>
  useMemo(() => {
    const filtered: MultiWebSummarySampleData[] = [];

    data.forEach((d, i) => {
      const hasAlerts = sampleHasAlerts[i];
      const failsFilter =
        filter !== "" && !d.data.id.toLowerCase().includes(filter);
      const failsAlerts = showOnlyWithAlerts && !hasAlerts;
      if (failsFilter || failsAlerts) {
        return;
      }

      filtered.push(d);
    });

    return filtered;
  }, [data, filter, sampleHasAlerts, showOnlyWithAlerts]);

/**
 * Take existing metrics and filter them down to a specific library
 * type, e.g. GEX, VDJ T, AB, etc.
 */
export const useMetricsPerSamplePerLibraryType = (
  metricsPerSample: JsonMetricSummary[][],
  libraryType: LibraryType,
) =>
  useMemo(
    () => metricsPerSample.map((m) => getMetricsForLibraryType(libraryType, m)),
    [libraryType, metricsPerSample],
  );

export type MetricWithBarcodes = {
  barcodes: string[];
  metric: JsonMetricSummary;
  sampleId: string;
};

/*
 * For a metric's key, e.g. "singlets_assigned_to_other_samples",
 * associate each metric back to that key via its sample.
 */
export function groupMultiSampleMetricsByKey(
  sampleMetricsBarcodes: SampleMetricsBarcodes,
  metricKey: string,
) {
  const barcodesMetrics: Array<MetricWithBarcodes> = [];
  sampleMetricsBarcodes.forEach((metricsBarcodes) => {
    const {
      sampleId,
      barcodes: barcodeIds,
      metrics: metricPerSample,
    } = metricsBarcodes;
    metricPerSample.forEach((metric) => {
      if (metric.key !== metricKey) {
        return;
      }

      barcodesMetrics.push({
        barcodes: barcodeIds,
        metric,
        // However, for a given key, e.g. 'total_singlets', the sampleId
        // must be kept with its metric and barcodes.
        sampleId,
      });
    });
  });
  return barcodesMetrics;
}
