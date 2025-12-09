import {
  AlertConfig,
  AlertLevel,
  ExperimentalDesign,
  JsonMetricSummary,
  MetricConfig,
  SampleWebSummary,
} from "./types";

// Helper function to create metric summaries for testing
export function makeMetric(
  key: string,
  grouping_key: number,
  value: string | number | { numerator: number; denominator: number },
): JsonMetricSummary {
  return {
    alerts: [],
    category: "Test Category",
    config: {
      alerts: [],
      header: "Test Header",
      type: typeof value === "string" ? "String" : "f64",
    },
    grouping_key,
    key,
    library_type: "Gene Expression",
    value,
  };
}

export function makeConfigMetric(
  config: Partial<MetricConfig>,
  ...args: Parameters<typeof makeMetric>
) {
  const m = makeMetric(...args);
  return {
    ...m,
    config: {
      ...m.config,
      ...config,
    },
  };
}

export function makeAlert(level: AlertLevel) {
  return {
    formatted_value: "Test",
    level,
    message: "Test",
    title: "Test Alert",
  };
}

export function makeAlertConfig(alertConfig: Partial<AlertConfig>) {
  return {
    conditions: {
      include_introns: null,
      is_rtl: null,
    },
    detail: "Test Detail",
    error_threshold: null,
    error_title: null,
    if_metric_is: null,
    warn_threshold: null,
    warn_title: null,
    ...alertConfig,
  };
}

export function makeMultiWebsummaryData(
  data: Partial<SampleWebSummary> = {},
  metrics?: JsonMetricSummary[],
) {
  return {
    data: {
      description: "test description",
      gex_tab: {
        alerts: [],
        content: {},
      },
      id: "test id",
      multiplexing_barcode_ids: ["A-A1", "B-H12"],
      ...data,
    },
    metrics: metrics ?? [makeMetric("metric_1", 1, 5)],
  };
}

export function makeExperimentalDesign(
  design: Partial<ExperimentalDesign> = {},
) {
  return {
    csv: "",
    include_introns: false,
    is_barnyard: false,
    is_rtl: false,
    multiplexing_method: null,
    ...design,
  };
}
